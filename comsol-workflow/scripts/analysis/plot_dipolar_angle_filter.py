"""Angle-limited complex EM reconstruction from the saved Fig.1 spectrum."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import c, epsilon_0, mu_0

from scripts.analysis.analyze_dipolar_singularity import WORK, CONFIG, digest, plt, save
from comsol_workflow.figure_output import save_figure_formats
from comsol_workflow.output_paths import table_path

ANGLES = (0.5, 1., 2., 2.5, 3., 3.5, 5., 10., 90.)
STEM = 'Fig.10_angle_filtered_surface_EM'
Z0 = np.sqrt(mu_0 / epsilon_0)


def cone_fields(q, electric, k0, theta):
    """Hard polar-angle low pass of complex E; upward H from Maxwell's relation."""
    if not np.isfinite(theta) or not 0 < theta <= 90:
        raise ValueError('Polar half-angle must be in (0, 90] degrees')
    q = np.asarray(q)
    if electric.shape != (len(q), len(q), 3) or not np.isfinite(electric).all():
        raise ValueError('Expected finite complex vector spectrum on the square q grid')
    qx, qy = np.meshgrid(q, q)
    radius2 = qx*qx + qy*qy
    mask = (radius2 < k0*k0) & (radius2 <= (k0*np.sin(np.deg2rad(theta)))**2)
    n = np.stack([qx, qy, np.sqrt(np.maximum(k0*k0-radius2, 0))], -1) / k0
    e = np.where(mask[..., None], electric, 0)
    h = np.cross(n, e) / Z0
    return e, h, mask


def energy(electric, magnetic):
    return (epsilon_0*np.sum(abs(electric)**2, axis=-1)
            + mu_0*np.sum(abs(magnetic)**2, axis=-1)) / 4


def run():
    from scripts.analysis.analyze_dipolar_complete import padded_transform_roi

    data_dir, logs = WORK/'01_results', WORK/'80_logs'
    source = data_dir/'Fig1_full_BZ_angular_spectrum_2x.npz'
    identity = logs/'finite_model_identity.json'
    config = CONFIG/'config.json'
    sources = {str(p): digest(p) for p in (source, identity, config)}
    paths = [WORK/'10_overview'/f'{STEM}.png', WORK/'11_pdf'/f'{STEM}.pdf',
             data_dir/f'{STEM}.npz', table_path(logs, f'{STEM}.csv'), logs/f'{STEM}.json',
             WORK/'12_reports'/f'{STEM}.md']
    if any(p.exists() for p in paths):
        raise FileExistsError('Angle-filter outputs already exist; do not silently overwrite')
    with np.load(source) as d:
        q, angular = d['q_m_inv'], d['E_up_surface_angular']
        z_surface, z_input = float(d['z_surface_m']), float(d['z_input_m'])
        source_fft_size = int(d['fft_size'])
        convention = str(d['time_convention'])
    assert convention == 'exp(+i omega t)'
    frequency = json.loads(identity.read_text(encoding='utf-8'))['frequency_THz']['re']*1e12
    k0 = 2*np.pi*frequency/c
    dq = float(np.diff(q).mean())
    np.testing.assert_allclose(q, (np.arange(len(q))-len(q)//2)*dq, atol=1e-7, rtol=1e-12)
    # The same q bins can be inverted on a coarser x grid without interpolation.
    # This spacing also resolves the highest intensity frequency, 2*k0.
    inverse_size = 2001
    dx = 2*np.pi/(inverse_size*dq)
    assert dx < np.pi/(2*k0) and inverse_size > len(q)
    half = int(np.ceil(100e-6/dx))
    points = 2*half+1
    x = np.arange(-half, half+1)*dx
    qx, qy = np.meshgrid(q, q)
    kz = np.sqrt(np.maximum(k0*k0-qx*qx-qy*qy, 0))
    propagating = qx*qx+qy*qy < k0*k0
    power_weights = dq*dq/(2*Z0*(2*np.pi)**2)*kz/k0
    full_power = float(np.sum(power_weights[propagating]*np.sum(abs(angular[propagating])**2, axis=-1)))
    e_all, h_all, masks, rows, inverse_errors = [], [], [], [], []
    for theta in ANGLES:
        ae, ah, mask = cone_fields(q, angular, k0, theta)
        ids = np.flatnonzero(mask.any(axis=0))
        lo, hi = int(ids[0]), int(ids[-1])+1
        reconstructed = []
        for field in (ae, ah):
            reconstructed.append(np.stack([padded_transform_roi(field[lo:hi, lo:hi, j], dx,
                inverse_size, points, inverse_transform=True) for j in range(3)], -1))
        e, h = reconstructed
        # Independent direct Fourier sums check the phase, sign, origin and scale.
        for iy, ix in ((half, half), (half+7, half-11)):
            phase = np.exp(-1j*(qx[mask]*x[ix]+qy[mask]*x[iy]))*dq*dq/(2*np.pi)**2
            for field, actual in ((ae, e), (ah, h)):
                expected = phase @ field[mask]
                error = float(np.linalg.norm(actual[iy, ix]-expected)/max(np.max(np.linalg.norm(actual, axis=-1)), 1e-300))
                assert error < 1e-10
                inverse_errors.append(error)
        retained = float(np.sum(power_weights[mask]*np.sum(abs(ae[mask])**2, axis=-1)))
        density = energy(e, h)
        rows.append(dict(theta_max_deg=theta, NA=float(np.sin(np.deg2rad(theta))),
                         retained_bins=int(mask.sum()), retained_power_W=retained,
                         retained_power_fraction=retained/full_power,
                         surface_peak_J_m3=float(density.max()),
                         surface_center_J_m3=float(density[half, half])))
        e_all.append(e); h_all.append(h); masks.append(mask)
        print(f'theta <= {theta:g} deg: power {retained/full_power:.3%}; inverse complete', flush=True)
    np.testing.assert_array_equal(masks[-1], propagating)
    assert np.all(np.diff([r['retained_power_fraction'] for r in rows]) >= 0)
    np.testing.assert_allclose(rows[-1]['retained_power_fraction'], 1., atol=1e-14)
    e_all, h_all, masks = np.array(e_all), np.array(h_all), np.array(masks)
    near = energy(e_all, h_all)
    near_scale = float(near.max())
    center = len(q)//2
    far_scale = float(energy(angular[center, center], np.cross([0, 0, 1], angular[center, center])/Z0))
    full_h = np.cross(np.stack([qx/k0, qy/k0, kz/k0], -1), np.where(propagating[..., None], angular, 0))/Z0
    far = energy(np.where(propagating[..., None], angular, 0), full_h)/far_scale
    np.savez_compressed(paths[2], theta_max_deg=ANGLES, q_m_inv=q, masks=masks,
                        x_um=x*1e6, E_surface_V_m=e_all, H_surface_A_m=h_all,
                        u_EM_J_m3=near, shared_near_scale_J_m3=near_scale,
                        angular_EM_center_scale=far_scale, frequency_Hz=frequency,
                        z_surface_m=z_surface, inverse_fft_size=inverse_size,
                        source_fft_size=source_fft_size, source_name=source.name)
    with paths[3].open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    cfg = json.loads(config.read_text(encoding='utf-8'))
    boundary = np.asarray(cfg['inner_boundary']); boundary = np.vstack([boundary, boundary[0]])
    footprint = np.asarray(cfg['finite_boundary']); footprint = np.vstack([footprint, footprint[0]])
    zoom = float(np.max(abs(footprint)))
    selection = np.flatnonzero(abs(q/k0) <= np.sin(np.deg2rad(10))+dq/k0)
    qa = q[selection]
    edges = np.rad2deg(np.arcsin(np.r_[qa-dq/2, qa[-1]+dq/2]/k0))
    plt.rcParams.update({'font.size': 10, 'xtick.direction': 'in', 'ytick.direction': 'in', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, len(ANGLES), figsize=(3.15*len(ANGLES)+.4, 6.8), layout='constrained')
    zoom_selection = np.flatnonzero(abs(x*1e6) <= zoom)
    zoom_peaks = [float(a[np.ix_(zoom_selection, zoom_selection)].max()) for a in near]
    zoom_bars = []
    phi = np.linspace(0, 2*np.pi, 721)
    extent = [x[0]*1e6-dx*5e5, x[-1]*1e6+dx*5e5]*2
    for i, (theta, row) in enumerate(zip(ANGLES, rows)):
        top = axes[0, i]
        far_im = top.pcolormesh(edges, edges, (far*masks[i])[np.ix_(selection, selection)],
                               cmap='magma', vmin=0, vmax=1, shading='flat', rasterized=True)
        title = f'{theta:g} deg' if theta < 90 else '90 deg (NA=1)'
        top.set(title=rf'$\theta\leq$ {title}'+'\n'+f'Power retained: {row["retained_power_fraction"]:.1%}',
                xlim=(-10, 10), ylim=(-10, 10), aspect='equal',
                xlabel=r'$\theta_x$ (deg)', ylabel=r'$\theta_y$ (deg)', xticks=[-10, 0, 10], yticks=[-10, 0, 10])
        radius = np.sin(np.deg2rad(theta))
        top.plot(np.rad2deg(np.arcsin(radius*np.cos(phi))),
                 np.rad2deg(np.arcsin(radius*np.sin(phi))), color='white', lw=.8,
                 ls='-', label='mask boundary')
        if theta > 10:
            top.text(.5, .94, 'Mask boundary outside view', transform=top.transAxes,
                     ha='center', va='top', color='.8', fontsize=8)
        ax = axes[1, i]
        near_im = ax.imshow(near[i], origin='lower', extent=extent, cmap='magma',
                            norm=plt.Normalize(0, zoom_peaks[i]), interpolation='nearest')
        ax.plot(footprint[:, 0], footprint[:, 1], color='.65', lw=.5, ls='--')
        ax.plot(boundary[:, 0], boundary[:, 1], color='white', lw=.6)
        ax.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), aspect='equal',
               xlabel=r'$x$ ($\mu$m)', ylabel=r'$y$ ($\mu$m)',
               xticks=[-zoom, 0, zoom], yticks=[-zoom, 0, zoom])
        ax.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%.0f'))
        ax.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%.0f'))
        if i:
            top.set_ylabel('')
            ax.set_ylabel('')
        bar = fig.colorbar(near_im, ax=ax, pad=.025, fraction=.045,
                           ticks=np.linspace(0, zoom_peaks[i], 3), format='%.3g',
                           label=r'$u_{EM}$ (J m$^{-3}$)')
        zoom_bars.append((bar, [ax]))
    far_bar = fig.colorbar(far_im, ax=list(axes[0]), pad=.015, fraction=.025, label=r'$\mathcal{U}_{EM}/\mathcal{U}_{EM}(0)$')
    far_bar.set_ticks(np.linspace(0, 1, 6))
    for bar in [far_bar, *(bar for bar, _ in zoom_bars)]:
        bar.ax.tick_params(direction='in', labelsize=8)
        bar.ax.yaxis.label.set_rotation(270); bar.ax.yaxis.labelpad=14
        bar.ax.yaxis.label.set_size(9)
    fig.suptitle(f'Angular low-pass filtering | Surface EM field at z = {z_surface*1e9:g} nm (air side)\n'
                 'Rows: filtered angular spectrum / structure-scale surface field', fontsize=15)
    fig.supxlabel('Complex E and H are filtered before inverse Fourier transformation. White: mask boundary (top) / cavity boundary (bottom); dashed: slab footprint.\n'
                  'Near-field panels: independent physical scales (J/m^3). All scales are linear; evanescent fields are excluded.', fontsize=10)
    save(fig, STEM, dpi=250)
    # Fixed square panels keep the many individual colorbars from expanding gaps.
    fig.set_layout_engine('none')
    width, height = fig.get_size_inches()
    for r in range(2):
        for i in range(len(ANGLES)):
            axes[r, i].set_position([(.6+3.15*i)/width, (.8+2.7*(1-r))/height,
                                     2.1/width, 2.1/height])
    for bar, panel in [(far_bar, axes[0, -1]), *((bar, panels[0]) for bar, panels in zoom_bars)]:
        box = panel.get_position()
        bar.ax.set_box_aspect(None)
        bar.ax.set_position([box.x1+.06/width, box.y0, .065/width, box.height])
        bar.ax.tick_params(labelsize=8)
        bar.ax.yaxis.labelpad=14
    fig.canvas.draw()
    save_figure_formats(fig, WORK, STEM, dpi=250, bbox_inches='tight')

    metadata = dict(source=str(source), source_sha256=sources, theta_max_deg=list(ANGLES),
                    frequency_Hz=frequency, z_surface_m=z_surface, z_input_m=z_input,
                    cutoff='hypot(kx,ky) <= k0*sin(theta_max), AND hypot(kx,ky) < k0',
                    theta_definition='polar half-angle from +z in air, not a square theta_x/theta_y crop',
                    field='upward propagating complex E/H, evaluated on slab top AIR side; not total dielectric near field',
                    source_already_backpropagated=True, extra_propagation_applied=False,
                    source_fft_size=source_fft_size, inverse_fft_size=inverse_size,
                    q_step_m_inv=dq, center_angle_step_deg=float(np.rad2deg(np.arcsin(dq/k0))),
                    x_step_um=dx*1e6, near_view_um=[-zoom, zoom], zoom_view_um=[-zoom, zoom],
                    stored_near_range_um=[float(x[0]*1e6), float(x[-1]*1e6)],
                    panel_rows=['filtered_angular_spectrum', 'structure_scale_surface_field'],
                    mask_outline='solid; exact circular k-space cutoff mapped to theta_x/theta_y',
                    mask_outline_outside_view_deg=[a for a in ANGLES if a > 10],
                    inverse_direct_sum_max_relative_error=max(inverse_errors),
                    full_NA_mask_exact=True, retained_power_monotonic=True,
                    shared_near_scale_J_m3=near_scale,
                    near_color_scale='linear', far_color_scale='linear', coordinate_scale='linear',
                    near_shared_scale_rows=[], zoom_color_scale='linear', zoom_color_units='J/m^3',
                    zoom_color_limits=[[0, peak] for peak in zoom_peaks],
                    zoom_normalization='none; independent physical color limits within the displayed window',
                    source_far_EM_center=far_scale, far_color_limits=[0, 1],
                    power_formula='sum |A_E|^2*(kz/k0)*dq^2 / (2*Z0*(2*pi)^2)',
                    rows=rows, COMSOL_started=False)
    for p, before in sources.items():
        assert digest(Path(p)) == before, p
    paths[4].write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    table = '\n'.join(f'| {r["theta_max_deg"]:g}° | {r["NA"]:.6f} | {r["retained_power_fraction"]:.3%} | {r["surface_peak_J_m3"]:.6g} |' for r in rows)
    paths[5].write_text(f'''# Fig.10：角度滤波后的表面 EM 近场

[PNG](../10_overview/{STEM}.png) · [PDF](../11_pdf/{STEM}.pdf)

读取 Fig.1 对应的复数向上电场角谱 {source.name}，不是对 EM 强度图取平方根或逆变换。输入来自 z={z_input*1e6:g} μm 的空气 E/H，已带相位回溯至板上表面 z={z_surface*1e9:g} nm；本次不重复传播，不启动 COMSOL。

θ 是相对 +z 法线的圆锥半角。保留 kx²+ky²≤k0²sin²θ 且 kx²+ky²<k0²，使用硬截止；各分量复相位不变。它是横向空间频率低通，90°表示完整传播光锥，排除 kz=0 的掠射边界。

取 A_H=(k/k0)×A_E/Z0；对过滤后的三个 E 和三个 H 分量分别做物理逆傅里叶变换：
E(rho)=sum A_E(q) exp(-i q·rho) dq²/(2π)²，H同理。沿用项目正号分析/负号重构约定；数值实现使用与该约定配对的变换，不能根据 FFT/iFFT 函数名称猜测符号。

最后计算空气侧时间平均能量密度 u_EM=(epsilon0|E|²+mu0|H|²)/4。该场不含倏逝分量，不是介质内部或完整表面储能场。Ex、Ey、Ez、Hx、Hy、Hz的复数数组及完整 EM 数值保存在同名NPZ。

图仅保留两行，并收紧面板间距。第一行为过滤后的 EM 能量角谱，固定显示±10°；90°参考仍保留整个传播光锥，未受该视窗裁剪。第二行为板上表面空气侧的结构范围近场放大图，各图独立使用0到当前视窗内峰值的线性色标，显示实际能量密度 u_EM（J/m³）；不同图中相同颜色不代表相同能量密度。第一行按原角谱中心值归一化，共用0–1线性色标。所有坐标轴均为线性。原±100 μm宽视野一行已移除，原始复场和完整视野数据仍保留在NPZ中。 第一行的实线表示圆形k空间掩膜在角度坐标中的实际边界；90°边界位于±10°视窗之外，已在对应面板注明。

| θ截止 | NA | 保留向上功率 | 表面峰值 J/m³ |
| --- | --- | --- | --- |
{table}

功率比例使用 kz/k0 加权的复电场角谱积分，与未过滤90°结果比较，不把图像像素或局域能量简单相加称作辐射功率。

角谱步长为{dq:.8g} m^-1，中心角度步长约{metadata['center_angle_step_deg']:.4f}°；0.5°孔径受有限角谱网格限制，未伪造更高分辨率。硬截止可产生旁瓣；较小截止角通常导致空间展宽。有限空气面测量窗口和传播场假设仍适用。

逆变换保留原 q 网格，采用{inverse_size}点计算网格，在x/y方向按{dx*1e6:.6g} μm取值；该间距满足最高辐射强度空间频率的Nyquist条件，未对图片插值或平滑。直接复数 Fourier 求和校验最大相对误差{max(inverse_errors):.3g}；90°掩膜恢复全部传播点，保留功率随θ单调增加。原角谱和来源文件SHA256未改变。

离线复现入口：scripts.analysis.plot_dipolar_angle_filter。已有同名输出时拒绝静默覆盖。原Fig.1–9保持不变。
''', encoding='utf-8')
    print(json.dumps({'figure': str(paths[0]), 'rows': rows, 'checks': max(inverse_errors)}, indent=2))


if __name__ == '__main__':
    run()
