"""Offline audit: periodic Gamma radiation, finite bright center, and 19-k hypothesis.

Run from the repository root: python -B -m scripts.analysis.plot_dipolar_gamma_mechanism
Uses saved mesh=5 fields only; never starts COMSOL or modifies source results.
"""
from __future__ import annotations

from comsol_workflow.output_paths import table_links
from comsol_workflow.output_paths import output_path

import argparse
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from comsol_workflow.farfield_fft import load_air_field
from scripts.analysis.visualize_dipolar_broadening_2d import OUTPUT, SOURCE, colorbar
from comsol_workflow.dipolar_lattice_analysis import load_data, window
from scripts.analysis.plot_dipolar_volume_radiation import CAVITY_SUPPORT, align_bar, index_map, save
from comsol_workflow.dipolar_radiation import pairs, planar_K0_to_far
from scripts.run_main.run_dipolar_periodic_reference import OUT
from comsol_workflow.result_io import digest


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def quadrant_integrals(air):
    """QI, QII, QIII, QIV; share samples on either axis between adjacent regions."""
    x, y = air['x_axis'], air['y_axis']
    dx, dy = np.diff(x).mean(), np.diff(y).mean()
    np.testing.assert_allclose(np.diff(x), dx, atol=1e-12)
    np.testing.assert_allclose(np.diff(y), dy, atol=1e-12)
    np.testing.assert_allclose(x, -x[::-1], atol=1e-12)
    np.testing.assert_allclose(y, -y[::-1], atol=1e-12)
    wx = [(1 - np.sign(x))/2, (1 + np.sign(x))/2]
    wy = [(1 - np.sign(y))/2, (1 + np.sign(y))/2]
    values = np.array([[np.sum(air[c]*wy[iy][:, None]*wx[ix][None, :])*dx*dy*1e-12
                        for c in ('Ex', 'Ey')] for iy, ix in [(1, 1), (1, 0), (0, 0), (0, 1)]])
    total = np.array([air[c].sum()*dx*dy*1e-12 for c in ('Ex', 'Ey')])
    np.testing.assert_allclose(values.sum(axis=0), total, rtol=1e-12, atol=abs(total[1])*1e-13)
    return values, total


def k_axes(ax, limit, letter='K'):
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect='equal',
           xlabel=rf'${letter}_x/(2\pi/a)$', ylabel=rf'${letter}_y/(2\pi/a)$')
    ax.set_xticks([-limit, 0, limit]); ax.set_yticks([-limit, 0, limit])


def complex_axes(ax, limits):
    ax.axhline(0, color='.8', lw=.6); ax.axvline(0, color='.8', lw=.6)
    ax.set(xlim=limits, ylim=limits, aspect='equal',
           xlabel=r'Re$(F/R)$', ylabel=r'Im$(F/R)$')


def arrow(ax, start, end, color, width=1.4):
    ax.annotate('', (end.real, end.imag), (start.real, start.imag),
                arrowprops={'arrowstyle': '->', 'color': color, 'lw': width})


def run(support_only=False):
    d = load_data(SOURCE)
    summary_path = output_path(OUTPUT, 's5_volume_summary.json')
    summary = read_json(summary_path)
    assert summary['basis_support'] == CAVITY_SUPPORT
    meta = read_json(output_path(OUTPUT, 'electric_19_summary.json'))
    manifest = read_json(OUT/'99_config/manifest.json')
    assert summary['mesh'] == manifest['mesh'] == 5
    assert summary['zeta'] == manifest['geometry']['zeta'] == 1.156
    assert summary['metadata']['native_air_error'] < summary['metadata']['conjugate_air_error']
    assert digest(SOURCE) == manifest['source_Hz_sha256']
    from pathlib import Path
    air_path = Path(meta['source'])
    assert digest(air_path) == meta['sha256']
    source_paths = [summary_path, SOURCE, air_path, output_path(OUTPUT, 's5_volume_radiation.npz'),
                    output_path(OUTPUT, 's5_cavity_support_control.npz'), OUT/'99_config/manifest.json']
    with np.load(output_path(OUTPUT, 's5_volume_radiation.npz')) as z:
        data = {key: z[key] for key in ['A', 'T0', 'K0_terms', 'k_um_inv', 'direct_air_F',
                'full_volume_F', 'direct_air_map_F', 'K_axis_over_2pi_a']}
    with np.load(output_path(OUTPUT, 's5_cavity_support_control.npz')) as z:
        cavity_t, cavity_terms = z['T0'], z['terms']
    with np.load(output_path(OUTPUT, 's5_finite_volume_o3.npz')) as z:
        omega, z_air = complex(z['omega_rad_s']), float(z['z_air_m'])
    np.testing.assert_allclose(data['k_um_inv'], d['k19'], atol=1e-12, rtol=0)
    np.testing.assert_allclose(data['K0_terms'], data['A'][:, :, None]*data['T0'], rtol=1e-12)
    c = []
    for n in range(19):
        path = OUT/f'01_results/k{n:02d}_reference.npz'
        assert digest(path) == read_json(OUT/f'01_results/k{n:02d}_audit.json')['data_sha256']
        source_paths.append(path)
        with np.load(path) as z:
            np.testing.assert_allclose(z['k_um_inv'], d['k19'][n], atol=1e-12, rtol=0)
            band = int(z['gamma_band_to_mode'][1])
            c.append(z['air_mean_E_at_z0_over_Hz_rms'][1, 1, band])
    c = np.array(c)
    ref = data['direct_air_F'][1]
    cube = data['K0_terms']/ref
    non_gamma = cube[1:, 1]
    np.testing.assert_allclose(data['T0'], cavity_t, atol=1e-24, rtol=1e-12)
    np.testing.assert_allclose(data['K0_terms'], cavity_terms, atol=1e-12, rtol=1e-12)
    assert abs(non_gamma).max() < 1e-12
    gamma_cavity = cube[0, 1]
    gamma_cladding = np.zeros(2, complex)  # Basis support is zero outside cavity.
    groups = np.array([cube[0, 1], non_gamma.sum(axis=0),
                       np.delete(cube, 1, axis=1).sum(axis=(0, 1)),
                       data['full_volume_F']/ref - cube.sum(axis=(0, 1))])
    np.testing.assert_allclose(groups.sum(axis=0), data['full_volume_F']/ref, atol=1e-12)
    air = load_air_field(air_path)
    quadrants_c, air_c = quadrant_integrals(air)
    np.testing.assert_allclose(planar_K0_to_far(air_c, z_air, omega), data['direct_air_F'],
                               rtol=1e-10, atol=abs(ref)*1e-12)
    quadrant_f = planar_K0_to_far(quadrants_c, z_air, omega)/ref
    s0 = window(-d['k19'], d)
    brute_s0 = np.exp(-1j*d['cell_centers']@d['k19'].T).mean(axis=0)
    np.testing.assert_allclose(s0, brute_s0, atol=1e-13, rtol=0)
    assert np.max(abs(s0[1:])) < 1e-12
    np.testing.assert_allclose(quadrant_f[:, 1], np.full(4, .25), atol=1e-12)
    np.testing.assert_allclose(quadrant_f[:, 0].sum(), 0, atol=1e-12)
    axis = data['K_axis_over_2pi_a']; middle = int(np.argmin(abs(axis)))
    actual = data['direct_air_map_F']/ref
    np.testing.assert_allclose(actual[middle, middle], data['direct_air_F']/ref, atol=1e-12)
    power = abs(c)**2; periodic_scale = power.max()
    xy = d['k19']/d['scale']; limit = float(axis[-1])
    xx, yy = np.meshgrid(axis, axis)
    kernel = window(np.stack([xx, yy], axis=-1)*d['scale'] - d['k19'][1], d)
    np.testing.assert_allclose(kernel[middle, middle], s0[1], atol=1e-13)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'xtick.direction': 'in', 'ytick.direction': 'in', 'pdf.fonttype': 42})

    gamma_total_ratio = float(power[0].sum()/power.sum(axis=1).max())
    if not support_only:
        fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.6), layout='constrained')
        index_map(axes[0], xy); axes[0].set_title('a  Original 19 k points')
        k_axes(axes[0], limit, 'k')
        for pol, name in enumerate(['Ex', 'Ey']):
            ax = axes[pol+1]
            im = ax.scatter(*xy.T, c=power[:, pol]/periodic_scale, s=155, cmap='magma',
                            norm=LogNorm(1e-10, 1), edgecolors='.65', linewidths=.4)
            ax.annotate(r'$\Gamma$', (0, 0), xytext=(9, 8), textcoords='offset points')
            ax.set_title(f'{chr(98+pol)}  Periodic py: {name} radiation'); k_axes(ax, limit, 'k')
        bar = fig.colorbar(im, ax=axes[1:], pad=.02,
                          label=r'$|c_\alpha(k)|^2 / \max_{n,\beta}|c_\beta(k_n)|^2$')
        fig.suptitle(r'Periodic py radiation | $\zeta=1.156$, mesh=5', fontsize=14)
        fig.supxlabel(f'Gamma / largest sampled total squared amplitude = {gamma_total_ratio:.3g}  |  '
                      'All eigenfields use the same Hz-RMS convention; no pointwise intensity normalization.', fontsize=9)
        align_bar(fig, bar, axes[1:]); save(fig, '19_s5_periodic_py_radiation')

        fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.8), layout='constrained')
        for pol, name in enumerate(['Ex', 'Ey']):
            ax = axes[pol]
            im = ax.imshow(abs(actual[:, :, pol])**2, origin='lower', extent=[-limit, limit]*2,
                           cmap='magma', norm=LogNorm(1e-8, 1))
            ax.set_title(f'{chr(97+pol)}  Actual finite cavity: {name}'); k_axes(ax, limit)
        bar = fig.colorbar(im, ax=axes[:2], pad=.02, label=r'$|F_\alpha^{air}(K)|^2 / |R|^2$')
        ax = axes[2]; ax.set_title('c  Four spatial quadrants at K=0')
        complex_axes(ax, (-.12, 1.12)); ax.set_ylim(-.62, .62)
        colors = ['#b44c42', '#28618a']
        for pol in range(2):
            ends = np.r_[0j, np.cumsum(quadrant_f[:, pol])]
            for n, (start, end) in enumerate(zip(ends[:-1], ends[1:])):
                arrow(ax, start, end, colors[pol], 1.8)
                if pol == 1:
                    ax.text((start.real+end.real)/2, -.075, f'Q{n+1}', ha='center', fontsize=8)
            ax.scatter(ends[-1].real, ends[-1].imag, color=colors[pol], marker='*', s=70,
                       label=['Ex sum = 0', 'Ey sum = 1'][pol], zorder=4)
        inset = ax.inset_axes([.16, .58, .40, .35])
        ends = np.r_[0j, np.cumsum(quadrant_f[:, 0])]
        for start, end in zip(ends[:-1], ends[1:]):
            arrow(inset, start, end, colors[0], 1.2)
        inset.axhline(0, color='.8', lw=.5); inset.axvline(0, color='.8', lw=.5)
        inset.set(xlim=(-.065, .012), ylim=(-.017, .035), title='Ex detail')
        inset.set_aspect('equal'); inset.tick_params(labelsize=7)
        ax.legend(frameon=False, fontsize=9, loc='lower right')
        fig.suptitle(r'Finite bright center: complex sums of the saved air field | $R=F_y^{air}(0)$', fontsize=14)
        fig.supxlabel('Q1: x>0,y>0; Q2: x<0,y>0; Q3: x<0,y<0; Q4: x>0,y<0.  '
                      'Spatial symmetry reconstruction; not an attribution to internal k points.', fontsize=9)
        align_bar(fig, bar, axes[:2]); save(fig, '20_s5_finite_center_symmetry')

    fig, axes = plt.subplots(2, 3, figsize=(14.4, 9.4), layout='constrained',
                             gridspec_kw={'wspace': .18})
    ax = axes[0, 0]
    im = ax.imshow(abs(kernel)**2, origin='lower', extent=[-limit, limit]*2,
                   cmap='magma', norm=LogNorm(1e-8, 1))
    ax.plot(0, 0, '+', color='#55bac3', ms=7)
    ax.annotate('K=0: exact null', (0, 0), xytext=(7, -13), textcoords='offset points',
                color='#55bac3', fontsize=8)
    ax.plot(*xy[1], '+', color='white', ms=6)
    ax.set_title('a  Cavity-window envelope, k index 1'); k_axes(ax, limit)
    colorbar(fig, ax, im, r'$|S(K-k_1)|^2$')
    for pol, name in enumerate(['Ex', 'Ey']):
        ax = axes[0, pol+1]; complex_axes(ax, (-.01, .01))
        ax.scatter(0, 0, s=25, color='#28618a')
        ax.set_title(f'{chr(98+pol)}  Non-Gamma py: {name} at K=0')
        ax.text(.04, .95, 'Exact cavity-window zero', transform=ax.transAxes, va='top', fontsize=9)
        ax.text(.04, .85, f'Max numerical |d_n/R| = {abs(non_gamma[:, pol]).max():.2g}',
                transform=ax.transAxes, va='top', fontsize=9)
        ax.text(.04, .08, 'Roundoff has no physical interference phase', transform=ax.transAxes, fontsize=8)
    ax = axes[1, 0]
    ax.semilogy(np.arange(19), abs(s0), 'o', ms=4, color='#28618a')
    ax.set(title='d  Original-window zeros at K=0',
           xlabel='Original k index', ylabel=r'$|S(-k_n)|$', xticks=[0, 3, 6, 9, 12, 15, 18],
           ylim=(1e-17, 10))
    ax.text(.05, .62, 'Original 1141-cell cavity window\n'
            r'$S(-k_n)=\delta_{n0}$'+'\nOff-Gamma values below are roundoff.', transform=ax.transAxes, fontsize=9)
    ax = axes[1, 1]
    magnitudes = abs(np.array([gamma_cavity[1], cube[:, :, 1].sum(), 1+0j]))
    ax.bar(np.arange(3), magnitudes, width=.55, color=['#28618a', '#7c6c9d', '#b44c42'])
    ax.set(yscale='log', ylim=(1e-4, 3), title='e  Cavity-only Ey center amplitudes',
           ylabel=r'$|F_y(0)/R|$', xticks=np.arange(3), xticklabels=['Gamma py', 'Cavity 19', 'Direct air'])
    ax.tick_params(axis='x', labelsize=8)
    for n, value in enumerate(magnitudes): ax.text(n, value*1.25, f'{value:.3g}', ha='center', fontsize=9)
    ax.text(.04, .96, 'Reference basis is zero in cladding', transform=ax.transAxes, va='top', fontsize=8)
    ax = axes[1, 2]
    values = np.r_[groups[:, 1], data['full_volume_F'][1]/ref, 1+0j]
    x = np.arange(len(values))
    ax.bar(x-.17, values.real, .34, color='#28618a', label='Real')
    ax.bar(x+.17, values.imag, .34, color='#b44c42', label='Imaginary')
    ax.axhline(0, color='.7', lw=.6)
    ax.set(title='f  Ey center: complete complex budget', ylabel=r'$F_y(0)/R$',
           xticks=x, xticklabels=['Gamma py', 'Off-Gamma py', 'Other bands', 'Field remainder', 'Full volume', 'Direct air'])
    ax.tick_params(axis='x', labelrotation=33, labelsize=8)
    for label in ax.get_xticklabels(): label.set_ha('right')
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Does the original 19-k expansion explain the bright center?', fontsize=15)
    fig.supxlabel(f'Cavity-only reference; cavity electric L2 residual {summary["electric_L2_residual_inside_cavity"]:.1%}.  '
                  f'Cavity 19 Ey amplitude: {abs(cube[:, :, 1].sum()):.3%} of actual center.  '
                  f'Volume/air discrepancy: {summary["volume_vs_air_relative_error"]:.1%}.', fontsize=9)
    save(fig, '21_s5_around_gamma_test')

    metrics = {
        'mesh': summary['mesh'], 'zeta': summary['zeta'], 'point_count': 19,
        'periodic_gamma_c_over_Hz_rms_ohm': pairs(c[0]),
        'periodic_gamma_total_squared_amplitude_over_max19': gamma_total_ratio,
        'periodic_gamma_Ey_squared_amplitude_over_max19_Ey': float(power[0, 1]/power[:, 1].max()),
        'nonGamma_py_K0_over_actual_Ey': pairs(non_gamma.sum(axis=0)),
        'nonGamma_py_coherent_over_incoherent': [None, None],
        'basis_support': CAVITY_SUPPORT,
        'nonGamma_origin_status': 'analytic zero; numerical residue has no physical coherence ratio',
        'nonGamma_py_Ey_amplitude_over_actual': float(abs(non_gamma[:, 1].sum())),
        'quadrant_amplitudes_over_actual_Ey': pairs(quadrant_f),
        'quadrant_order': ['x+,y+', 'x-,y+', 'x-,y-', 'x+,y-'],
        'Gamma_py_cavity_part_same_coefficient': pairs(gamma_cavity),
        'Gamma_py_cladding_part_same_coefficient': pairs(gamma_cladding),
        'cavity_support_19_Ey_amplitude_over_actual': float(abs(cavity_terms[:, :, 1].sum()/ref)),
        'window_nonGamma_max_abs_K0': float(abs(s0[1:]).max()),
        'electric_L2_residual_inside_cavity': summary['electric_L2_residual_inside_cavity'],
        'volume_vs_air_relative_error': summary['volume_vs_air_relative_error'],
        'direct_Ey_peak_K_over_2pi_a': [float(axis[v]) for v in np.unravel_index(abs(actual[:, :, 1]).argmax(), actual.shape[:2])[::-1]],
        'checks': 'source identities, native phase, A*T, complete budget, quadrants/full air, window/direct sum, center pixel passed',
        'source_sha256': {str(path): digest(path) for path in source_paths},
        'interpretation': 'periodic weak Gamma and finite bright normal emission coexist; original 19-k envelope attribution is not established',
    }
    (output_path(OUTPUT, 's5_gamma_mechanism_metrics.json')).write_text(json.dumps(metrics, indent=2), encoding='utf-8')
    with (output_path(OUTPUT, 's5_gamma_mechanism_points.csv')).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index', 'kx_over_2pi_a', 'ky_over_2pi_a', 'abs_S_at_K0',
                         'periodic_cx_re_ohm', 'periodic_cx_im_ohm', 'periodic_cy_re_ohm', 'periodic_cy_im_ohm',
                         'finite_py_Ex_over_R_re', 'finite_py_Ex_over_R_im', 'finite_py_Ey_over_R_re', 'finite_py_Ey_over_R_im'])
        for n in range(19):
            writer.writerow([n, *xy[n], abs(s0[n]), *np.stack([c[n].real, c[n].imag], axis=-1).ravel(),
                             *np.stack([cube[n, 1].real, cube[n, 1].imag], axis=-1).ravel()])
    report = table_links(f'''# cavity-only 周期参考与有限结构中心辐射

固定 ζ={summary['zeta']}、mesh={summary['mesh']}、mode19，保留原19个 k 点、每点四条带。py 使用用户命名，Cartesian Ex/Ey 不交换。仅后处理现有解，未启动 COMSOL。

**当前周期参考基底只覆盖 cavity cells，在 cladding 和外围留白区域恒为零。** 原全结构延拓分析已被取代，历史文件保存在 `../90_history/s5_superseded_allcell_20260921.zip`。

## 计算范围

1. 实际有限结构的复电场仍来自原始 COMSOL 解，包含实际 cavity/cladding 相互作用；没有重新求解移除 cladding 的结构。
2. 以 cavity 周期模式构造基底：`psi_nb(R+rho)=E_cell,nb(rho)*exp(-i k_n.R)`，仅当 R 属于 cavity 时非零。
3. 仅在 cavity 内最小化 `integral epsilon_r |E_finite - sum A_nb psi_nb|² dV`，同时使用 Ex/Ey/Ez 求复系数 A。它们已经重新投影，不沿用旧的全结构系数。
4. 仅在 cavity 材料体积内积分 `T_alpha,nb(0)=k0²/(4*pi)*integral (epsilon_r-1) psi_alpha,nb exp(+ik0 z) dV`。所有基函数响应使用有限模式同一个复频率；非零观察 K 保留横向投影和传播相位。
5. `d_alpha,nb(K)=A_nb T_alpha,nb(K)`。每个 k 的图是四条带复振幅相干相加后的模方；总图再对19个 k 相干相加。所有可比较图使用同一个 R=F_y^air(0)，不逐点归一化。

cladding 不参与 A 的拟合，也不参与 cavity 基底的 T 积分。残差定义为完整实际材料体积分减去 cavity 19点的贡献，包含 cavity 内未重构场、全部实际 cladding 场和外围留白区；不能把残差唯一称为 cladding 辐射。

## 更新后的结果

- cavity 内三维电场相对 L2 残差：{summary['electric_L2_residual_inside_cavity']:.6%}。这与旧全结构误差的积分区域不同。
- 19点 Ey 中心复振幅/R：{cube[:, :, 1].sum():.10g}，其模为实际中心幅度的 **{abs(cube[:, :, 1].sum()):.6%}**。这是幅度比，不是能量百分比。
- Γ-py 的 Ey 中心复振幅/R：{gamma_cavity[1]:.10g}；其 cladding 基底贡献按支持范围定义为零。
- 原 cavity 窗口 `S(q)=mean_R exp(i q.R)` 对原商格点满足 `S(-k_n)=delta_(n,0)`。18个非Γ的 K=0 项解析为零；窗口数值残余最大 {abs(s0[1:]).max():.3g}，未强制改写复数据为零。
- 非Γ-py的 Ey 合成数值残余/R为 {non_gamma[:, 1].sum():.5g}，不解释其相位，也不计算零点的相干/非相干比。
- 完整材料体积分与独立空气面中心复幅差仍为 {summary['volume_vs_air_relative_error']:.6%}；限制支持范围没有消除这个独立对照差异。

## 图19–21的含义

图19的周期辐射系数采用 `c_alpha(k)=exp(+ikz z)*mean_cell[E_alpha exp(+ik.rho)]/Hz_rms`。Γ 的总平方幅度/19点最大值为 {gamma_total_ratio:.8g}，|cx|={abs(c[0,0]):.8g} Ω、|cy|={abs(c[0,1]):.8g} Ω。固定mesh支持弱辐射，但不证明严格数学零点。

图20从实际空气面直接积分 `C_alpha(K)=integral E_alpha(rho,z_air)exp(+iK.rho)dA`，再以 `F=i kz/(2*pi)*exp(+ikz z_air)*C` 转为球面振幅，未使用周期基底。Ey中心仍是当前观察范围内的峰值。四个镜像象限各贡献R/4，Ex相消；这是原对称性重构的实空间检查，不是非Γ内部态归因。

图21a/d展示 cavity 窗口和中心正交零，b/c明确显示非Γ中心零点，e比较 cavity Γ-py、cavity 19点总和与实际中心幅度，f保留完整复数预算。图12/13继续使用index map首格和Ex/Ey共同色标，但现在只计算 cavity cells。

因此，周期 Γ 弱辐射与实际有限结构中心亮点共存的事实保留；当前 cavity-only 的19点分解仍未定量解释亮点，不能把非Γ中心的数值舍入残余当作相干增强。

## 产物与复现

- 图12–18及图21的PNG/PDF已更新；独立周期图19和实际空气面图20的物理数据没有变化。
- 当前数据：`../01_results/s5_volume_radiation.npz`、`../80_logs/s5_volume_summary.json`、`../80_logs/s5_volume_K0_terms.csv`；支持范围已显式记为cavity-only。
- 中心控制结果：`../80_logs/s5_cavity_support_control.json` / `../01_results/s5_cavity_support_control.npz`；逐点审计：`../80_logs/s5_gamma_mechanism_points.csv`；指标与源哈希：`../80_logs/s5_gamma_mechanism_metrics.json`。
- 历史延拓分析只在 `../90_history/s5_superseded_allcell_20260921.zip` 中作为已取代记录保留。原模型、有限场及周期参考数据未修改。
''')
    (output_path(OUTPUT, 's5_gamma_mechanism_report.md')).write_text(report, encoding='utf-8')
    (output_path(OUTPUT, 's5_mesh5_report.md')).write_text(
        '# 固定 mesh=5：当前采用 cavity-only 周期参考\n\n'
        'cavity 周期场的投影和辐射积分仅覆盖 cavity cells，在 cladding 内为零。'
        '当前完整说明及图12–18、图21的计算依据见 '
        '[cavity-only 计算说明](s5_gamma_mechanism_report.md)。\n\n'
        f'cavity 内场重构误差 {summary["electric_L2_residual_inside_cavity"]:.4%}；'
        f'19点 Ey 中心合成幅度为实际值的 {abs(cube[:, :, 1].sum()):.4%}。'
        '原19点非Γ包络在中心为正交零；当前分解不能解释完整亮点。\n\n'
        '已取代的全结构延拓数据与说明保存在 `../90_history/s5_superseded_allcell_20260921.zip`。'
        '原始COMSOL模型及复场数据未改变。\n', encoding='utf-8')
    validation_path = output_path(OUTPUT, 's5_validation_summary.json')
    validation = read_json(validation_path)
    validation['basis_support'] = CAVITY_SUPPORT
    validation['status'] = 'fixed_mesh5_cavity_only_workflow_complete_mechanism_not_established'
    validation['scientific_gates'].update({
        'full_vector_internal_projection': f'complete; cavity-only electric L2 residual {summary["electric_L2_residual_inside_cavity"]:.4%}',
        'independent_finite_radiation_response': 'complete; cavity-supported reference, original FEM quadrature',
        'complete_farfield_closure': f'not established; volume/air discrepancy {summary["volume_vs_air_relative_error"]:.4%}; cavity19/air complex-amplitude error {summary["basis_vs_air_relative_error"]:.4%}'})
    validation_path.write_text(json.dumps(validation, indent=2), encoding='utf-8')

    print(json.dumps({key: value for key, value in metrics.items() if key != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-only', action='store_true', help='Update figure 21 and reports; preserve independent figures 19/20')
    run(support_only=parser.parse_args().support_only)
