"""Offline Fig.8 interference and Fig.9 directional cell budgets; no COMSOL."""
from __future__ import annotations

from comsol_workflow.output_paths import table_links, table_path
import argparse
import json

import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter
from scipy.constants import epsilon_0, mu_0

from comsol_workflow.hex_lattice_utils import unit_cell_corners
from scripts.analysis import analyze_dipolar_complete as a
from scripts.analysis.analyze_dipolar_cell_radiation import budgets
from scripts.analysis.analyze_dipolar_singularity import digest
from scripts.analysis.analyze_dipolar_spatial_filter import NAMES

FIG8 = "Fig.8_sector_interference"
FIG9 = "Fig.9_unitcell_sector_contributions"
REPORT = "sector_interference_report.md"
HEADING = "## Fig.8–9：上下与左右分组的相干干涉"
GROUPS = ["Top + bottom", "Left/right four"]
GROUP_IDS = [[1, 4], [0, 2, 3, 5]]


def em_energy(e, h):
    return (epsilon_0 * np.sum(abs(e)**2, axis=-1)
            + mu_0 * np.sum(abs(h)**2, axis=-1)) / 4


def interference_terms(ea, ha, eb, hb, denominator):
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("The common full-field center denominator must be positive")
    ia, ib = em_energy(ea, ha) / denominator, em_energy(eb, hb) / denominator
    ke = epsilon_0 * np.real(ea * eb.conj()) / (2 * denominator)
    kh = mu_0 * np.real(ha * hb.conj()) / (2 * denominator)
    cross = (ke + kh).sum(axis=-1)
    total = em_energy(ea + eb, ha + hb) / denominator
    np.testing.assert_allclose(ia + ib + cross, total, rtol=1e-12, atol=1e-14)
    return dict(I_A=ia, I_B=ib, self_sum=ia + ib, cross=cross, total=total,
                C_A=ia + cross / 2, C_B=ib + cross / 2,
                cross_E_components=ke, cross_H_components=kh)


def cell_sector_weights(xy, origin):
    """Whole-cell directional attribution, with equal weights on shared rays."""
    local = np.asarray(xy) - origin
    angle = np.arctan2(local[:, 1], local[:, 0])
    directions = np.deg2rad(np.arange(30, 360, 60))
    distance = abs(np.angle(np.exp(1j * (angle[None] - directions[:, None]))))
    weights = (distance <= np.pi / 6 + 1e-12).astype(float)
    weights[:, np.linalg.norm(local, axis=1) < 1e-12] = 1
    weights /= weights.sum(axis=0)
    np.testing.assert_allclose(weights.sum(axis=0), 1, rtol=0, atol=1e-15)
    return weights


def self_check():
    e = np.array([1., 0., 0.], dtype=complex)
    n = np.array([0., 0., 1.])
    h = np.cross(n, e) / a.Z0
    u = float(em_energy(e, h))
    for factor, expected in [(1, (2., 4.)), (-1, (-2., 0.)), (1j, (0., 2.))]:
        result = interference_terms(e, h, factor * e, factor * h, u)
        np.testing.assert_allclose([result['cross'], result['total']], expected, atol=1e-14)
    orthogonal = np.array([0., 1., 0.], dtype=complex)
    result = interference_terms(e, h, orthogonal, np.cross(n, orthogonal) / a.Z0, u)
    np.testing.assert_allclose([result['cross'], result['total']], [0., 2.], atol=1e-14)
    xy = np.array([[0., 0.], [1., 0.], [0., 1.], [0., -1.], [.3, .8]])
    weights = cell_sector_weights(xy, np.zeros(2))
    np.testing.assert_allclose(weights[:, 0], 1 / 6)
    np.testing.assert_allclose(weights[:, 1], [.5, 0, 0, 0, 0, .5])
    np.testing.assert_allclose(weights[:, 2], [0, 1, 0, 0, 0, 0])
    np.testing.assert_allclose(cell_sector_weights(-xy, np.zeros(2)), np.roll(weights, 3, axis=0))
    fields = np.array([[1+2j, .3j], [-.4j, 2.], [3., -1j], [1j, -.7], [.5, 1j]])
    total, _, signed = budgets(fields)
    np.testing.assert_allclose((weights @ fields).sum(axis=0), total)
    np.testing.assert_allclose((weights @ signed).sum(), signed.sum())
    print("Synthetic coherent/antiphase/orthogonal and sector-weight checks passed.", flush=True)


def plot_interference(q, maps, center):
    half = len(q) // 2
    dq = float(q[1] - q[0])
    edges = np.rad2deg(np.arcsin(np.r_[q - dq / 2, q[-1] + dq / 2] / a.K0))
    angle = np.rad2deg(np.arcsin(q / a.K0))
    specs = [("I_A", "a  Top + bottom"), ("I_B", "b  Left/right four"),
             ("self_sum", "c  Self terms: A + B"),
             ("cross", "d  Interference between A and B"),
             ("total", "e  Coherent sum: A + B + interference")]
    fig, axes = a.plt.subplots(2, 3, figsize=(17.2, 11.2), layout="constrained")
    for ax, (key, title) in zip(axes.flat, specs):
        signed = key == "cross"
        norm = Normalize(-1 if signed else 0, 1)
        assert maps[key].min() >= norm.vmin - 1e-12 and maps[key].max() <= 1 + 1e-12
        artist = ax.pcolormesh(edges, edges, maps[key], cmap="RdBu_r" if signed else "magma",
                               norm=norm, shading="flat", rasterized=True)
        ax.set(title=title, xlabel=r"$\theta_x$ (°)", ylabel=r"$\theta_y$ (°)",
               xlim=(-10, 10), ylim=(-10, 10), aspect="equal",
               xticks=[-10, -5, 0, 5, 10], yticks=[-10, -5, 0, 5, 10])
        label = (r"$K_{AB}/\mathcal{U}_{\rm full}(0)$ (1)" if signed
                 else r"$\mathcal{U}/\mathcal{U}_{\rm full}(0)$ (1)")
        bar = a.colorbar(fig, ax, artist, label)
        bar.set_ticks(np.linspace(norm.vmin, norm.vmax, 5))
        ax.text(.035, .045, f"Center: {center[key]:+.2%}" if signed else f"Center: {center[key]:.2%}",
                transform=ax.transAxes, fontsize=10, color="black" if signed else "white")
        assert ax.get_aspect() == 1.0 and artist.get_clim() == (norm.vmin, 1.)
    ax = axes[1, 2]
    incoherent, cross, total = (maps[k][half] for k in ["self_sum", "cross", "total"])
    ax.fill_between(angle, incoherent, total, where=cross >= 0, color="#d6604d", alpha=.22)
    ax.fill_between(angle, incoherent, total, where=cross < 0, color="#4393c3", alpha=.22)
    ax.plot(angle, total, color="#202020", lw=1.6, label="Coherent total")
    ax.plot(angle, incoherent, color="#737373", lw=1.4, ls="--", label="Self terms")
    ax.plot(angle, cross, color="#b2182b", lw=1.4, label="Interference")
    ax.axhline(0, color="#aaaaaa", lw=.6)
    ax.set(title=r"f  Central cut: $\theta_y=0^\circ$", xlabel=r"$\theta_x$ (°)",
           ylabel=r"Normalized EM angular spectrum (1)", xlim=(-10, 10),
           ylim=(min(-.05, float(np.floor(cross.min() * 20) / 20 - .02)), 1.05),
           xticks=[-10, -5, 0, 5, 10], box_aspect=1)
    assert ax.get_ylim()[0] < cross.min() and ax.get_ylim()[1] > total.max()
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.text(.98, .97, f"Center\n{center['self_sum']:.4f} + {center['cross']:.4f} = 1", ha="right", va="top",
            transform=ax.transAxes, fontsize=10)
    fig.suptitle("Far field intensity | Top + bottom (A) and Left/right four (B)\n"
                 "Common normalization: full-field center = 1", fontsize=14)
    fig.supxlabel(f"At center: signed contribution of A = {center['C_A']:.2%}; B = {center['C_B']:.2%}. "
                  "Each receives half the cross term.\nRed interference increases intensity; blue interference decreases it.", fontsize=10)
    assert len(axes.flat) == 6
    a.save(fig, FIG8, dpi=a.FIGURE_DPI)


def plot_cell_groups(data, weights, signed, cfg, checks, pad_fraction):
    corners = data['centers_um'][:, None, :] + unit_cell_corners(float(cfg['a']))[None, :, :]
    regions = data['regions'][:-1]
    limits = checks['figure']['signed_color_limits']
    bound = checks['figure']['xy_limit_um']
    norm = Normalize(*limits)
    boundary = np.asarray(cfg['inner_boundary'])
    outer = np.asarray(cfg['finite_boundary'])
    origin = (boundary.max(axis=0) + boundary.min(axis=0)) / 2
    maps = np.empty((2, 3, len(signed)))
    fig, axes = a.plt.subplots(2, 3, figsize=(19.4, 11.4), layout="constrained")
    for row, group in enumerate(GROUPS):
        for col, region in enumerate(['cavity', 'cladding', 'cavity + cladding']):
            ax = axes[row, col]
            selected = np.ones(len(signed), bool) if col == 2 else regions == region
            w = weights[row] * selected
            values = maps[row, col] = w * signed
            active = w > 0
            ax.add_collection(PolyCollection(corners, facecolors='#e9e9e9', edgecolors='none', rasterized=True))
            cells = PolyCollection(corners[active], array=values[active], cmap='RdBu_r', norm=norm,
                                   edgecolors='none', rasterized=True)
            ax.add_collection(cells)
            for polygon, width in [(boundary, .85), (outer, .6)]:
                closed = np.vstack([polygon, polygon[0]])
                ax.plot(closed[:, 0], closed[:, 1], color='#444444', lw=width)
            for theta in np.deg2rad(np.arange(0, 360, 60)):
                end = origin + 2 * bound * np.array([np.cos(theta), np.sin(theta)])
                ax.plot([origin[0], end[0]], [origin[1], end[1]], color='#838383', lw=.55, ls='--')
            ax.set(title=f"{'abcdef'[3 * row + col]}  {group} | {region}",
                   xlabel=r"$x$ (µm)", ylabel=r"$y$ (µm)", xlim=(-bound, bound), ylim=(-bound, bound), aspect='equal')
            ax.text(.035, .045, f"Signed sum: {values.sum():+.2%}", transform=ax.transAxes, fontsize=11,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.8, pad=2))
            bar = a.colorbar(fig, ax, cells, r"$w_{gj}C_j/I_{\rm source}(0)$ (1)")
            bar.formatter = ScalarFormatter(useMathText=True)
            bar.formatter.set_powerlimits((0, 0))
            bar.update_ticks()
            assert ax.get_aspect() == 1.0 and cells.get_clim() == tuple(limits)
    np.testing.assert_allclose(maps[:, 2], maps[:, 0] + maps[:, 1], rtol=0, atol=1e-18)
    np.testing.assert_allclose(maps[:, 2].sum(axis=0), signed, rtol=1e-14, atol=1e-18)
    fig.suptitle("Directional grouping of Fig.7e | Signed coherent contribution at the far-field center\n"
                 "Rows: two upper/lower sectors and four side sectors | Common color scale from Fig.7e", fontsize=14)
    fig.supxlabel(f"Pad, kept separate without directional assignment: {pad_fraction:+.2%}. "
                  "Cells + pad = 100%. Gray cells are outside the selected group.\n"
                  "Whole cells grouped by center; shared-ray cells have half weight. This is a cell-integral map.", fontsize=10)
    a.save(fig, FIG9, dpi=a.FIGURE_DPI)
    return maps


def run():
    self_check()
    previous = [p for folder in [a.O, a.WORK / '11_pdf']
                for i in range(1, 8) for p in folder.glob(f'Fig.{i}_*')]
    previous += list((a.WORK / '00_model').rglob('*.mph'))
    snapshot = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in previous}
    spatial_path = a.D / 'spatial_filter_six_sectors.npz'
    cell_path = a.D / 'unitcell_center_radiation.npz'
    spatial_checks = json.loads((a.L / 'spatial_filter_checks.json').read_text(encoding='utf8'))
    cell_checks = json.loads((a.L / 'unitcell_center_radiation_checks.json').read_text(encoding='utf8'))
    cfg = json.loads((a.CONFIG / 'config.json').read_text(encoding='utf8'))
    with np.load(spatial_path) as d:
        q, e, h = d['q_m_inv'], d['merged_E_up_angular'], d['merged_H_up_angular']
        gram, old_far = d['center_gram'], d['merged_far_EM_normalized']
        denominator = float(d['far_denominator'])
    assert np.isfinite(e).all() and np.isfinite(h).all()
    amplitude_errors = {}
    for label, fields in [('E', e), ('H', h)]:
        scale = float(abs(fields[0]).max())
        assert scale > 0
        # Normalize before comparison so symmetry-forced zeros use the field's physical scale.
        np.testing.assert_allclose((fields[1] + fields[2]) / scale, fields[0] / scale,
                                   rtol=1e-11, atol=1e-13)
        amplitude_errors[label] = float(np.linalg.norm(fields[1] + fields[2] - fields[0]) / np.linalg.norm(fields[0]))
    half = len(q) // 2
    assert abs(q[half]) < 1e-12
    np.testing.assert_allclose(em_energy(e[0, half, half], h[0, half, half]), denominator, rtol=1e-12, atol=0)
    terms = interference_terms(e[1], h[1], e[2], h[2], denominator)
    for key, saved in [('total', old_far[0]), ('I_A', old_far[1]), ('I_B', old_far[2])]:
        np.testing.assert_allclose(terms[key], saved, rtol=1e-11, atol=1e-14)
    center = {k: float(v[half, half]) for k, v in terms.items() if v.ndim == 2}
    for key, ids in zip(['C_A', 'C_B'], GROUP_IDS):
        np.testing.assert_allclose(center[key], gram.real[ids].sum(), rtol=1e-12, atol=1e-14)
    center.update(Top_self=float(gram[1, 1].real), Bottom_self=float(gram[4, 4].real),
                  Top_Bottom_cross=float(2 * gram[1, 4].real),
                  Top_Bottom_correlation_real=float((gram[1, 4] / np.sqrt(gram[1, 1].real * gram[4, 4].real)).real))
    np.testing.assert_allclose(center['Top_self'] + center['Bottom_self'] + center['Top_Bottom_cross'], center['I_A'], rtol=1e-11)
    a.plt.rcParams.update({'font.size': 10, 'xtick.direction': 'in', 'ytick.direction': 'in', 'pdf.fonttype': 42})
    plot_interference(q, terms, center)
    np.savez_compressed(a.D / (FIG8 + '.npz'), q_m_inv=q, **terms,
                        far_denominator=denominator, normalization='same-window full EM angular spectrum at center')
    pd.DataFrame([{'quantity': k, 'normalized_center_value': v} for k, v in center.items()]).to_csv(table_path(a.L, FIG8 + '_center.csv'), index=False)

    data = a.load(cell_path.name)
    assert 'geometric cell clipping' in str(data['integration_method'])
    assert data['regions'][-1] == 'pad' and len(data['centers_um']) + 1 == len(data['F_center_V'])
    np.testing.assert_allclose(data['frequency_Hz'], a.FREQ, rtol=1e-13)
    f = data['F_center_V']
    ft, own, signed_w = budgets(f)
    np.testing.assert_allclose(ft, data['F_total_V'], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(signed_w, data['signed_W_sr'], rtol=1e-12, atol=1e-16)
    source_center = float(data['source_center_W_sr'])
    np.testing.assert_allclose(np.sum(abs(ft)**2) / (2 * a.Z0), source_center, rtol=1e-12)
    signed = signed_w[:-1].sum(axis=1) / source_center
    origin = np.asarray(spatial_checks['origin_um'])
    boundary = np.asarray(cfg['inner_boundary'])
    np.testing.assert_allclose(origin, (boundary.max(axis=0) + boundary.min(axis=0)) / 2, atol=1e-12)
    weights = cell_sector_weights(data['centers_um'], origin)
    np.testing.assert_allclose(cell_sector_weights(2 * origin - data['centers_um'], origin), np.roll(weights, 3, axis=0), atol=1e-14)
    merged_weights = np.array([weights[ids].sum(axis=0) for ids in GROUP_IDS])
    np.testing.assert_allclose((weights @ f[:-1]).sum(axis=0) + f[-1], ft, rtol=1e-12, atol=1e-12)
    pad_fraction = float(signed_w[-1].sum() / source_center)
    np.testing.assert_allclose(signed.sum() + pad_fraction, 1, rtol=1e-12)
    cell_maps = plot_cell_groups(data, merged_weights, signed, cfg, cell_checks, pad_fraction)
    rows = []
    for kind, labels, allocation in [('six sectors', NAMES, weights), ('two groups', GROUPS, merged_weights)]:
        for label, w in zip(labels, allocation):
            for region in ['cavity', 'cladding']:
                wr = w * (data['regions'][:-1] == region)
                fg = wr @ f[:-1]
                contribution = float(np.real(fg * ft.conj()).sum() / (2 * a.Z0))
                np.testing.assert_allclose(contribution / source_center, wr @ signed, rtol=1e-12, atol=1e-14)
                rows.append(dict(partition=kind, group=label, region=region, weighted_cell_count=float(wr.sum()),
                                 Fx_re_V=float(fg[0].real), Fx_im_V=float(fg[0].imag),
                                 Fy_re_V=float(fg[1].real), Fy_im_V=float(fg[1].imag),
                                 group_self_W_sr=float(np.sum(abs(fg)**2) / (2 * a.Z0)),
                                 signed_W_sr=contribution, signed_fraction=contribution / source_center))
    for region in ['cavity', 'cladding', 'pad', 'all']:
        keep = np.ones(len(f), bool) if region == 'all' else data['regions'] == region
        fg = f[keep].sum(axis=0)
        contribution = float(np.real(fg * ft.conj()).sum() / (2 * a.Z0))
        rows.append(dict(partition='complete region', group='all directions', region=region,
                         weighted_cell_count=int(keep[:-1].sum()), Fx_re_V=float(fg[0].real), Fx_im_V=float(fg[0].imag),
                         Fy_re_V=float(fg[1].real), Fy_im_V=float(fg[1].imag),
                         group_self_W_sr=float(np.sum(abs(fg)**2) / (2 * a.Z0)),
                         signed_W_sr=contribution, signed_fraction=contribution / source_center))
    pd.DataFrame(rows).to_csv(table_path(a.L, FIG9 + '_groups.csv'), index=False)
    np.savez_compressed(a.D / (FIG9 + '.npz'), centers_um=data['centers_um'], regions=data['regions'][:-1],
                        sector_names=NAMES, sector_weights=weights, group_names=GROUPS, group_weights=merged_weights,
                        F_six_sectors_V=weights @ f[:-1], F_two_groups_V=merged_weights @ f[:-1],
                        signed_cell_fraction=signed, grouped_cell_fraction=cell_maps,
                        F_pad_V=f[-1], pad_signed_fraction=pad_fraction, F_total_V=ft,
                        source_center_W_sr=source_center, origin_um=origin,
                        grouping='whole-cell centers; shared-ray weights; pad unassigned')
    write_report(center, cell_maps.sum(axis=-1), pad_fraction, spatial_checks, cell_checks, source_center)
    sources = [spatial_path, cell_path, a.L / 'spatial_filter_checks.json',
               a.L / 'unitcell_center_radiation_checks.json', a.CONFIG / 'config.json']
    unchanged = all((p.stat().st_size, p.stat().st_mtime_ns) == snapshot[str(p)] for p in previous)
    assert unchanged
    checks = dict(COMSOL_started=False, inputs={str(p): dict(sha256=digest(p), bytes=p.stat().st_size) for p in sources},
                  center=center, coherent_amplitude_relative_errors=amplitude_errors, angle_grid_size=len(q), existing_fft_size=spatial_checks['fft_size'],
                  interference_identity_max_absolute_error=float(abs(terms['self_sum'] + terms['cross'] - terms['total']).max()),
                  source_signed_sum_including_pad=float(signed.sum() + pad_fraction), pad_fraction=pad_fraction,
                  cell_weight_sum_max_error=float(abs(weights.sum(axis=0) - 1).max()),
                  shared_ray_cells=int(np.sum((weights > 0).sum(axis=0) == 2)),
                  origin_cells=int(np.sum((weights > 0).sum(axis=0) == 6)),
                  cell_group_signed_fractions=cell_maps.sum(axis=-1).tolist(),
                  old_figures_and_models_unchanged=unchanged, preserved_file_count=len(previous),
                  source_to_air_complex_error=cell_checks['cell_sum_vs_air_complex_relative_error'],
                  source_to_air_intensity_error=cell_checks['source_vs_air_intensity_relative_error'],
                  independent_cell_convergence=False, whole_cell_grouping_not_sector_volume_clipping=True,
                  figures=[FIG8, FIG9], dpi=a.FIGURE_DPI)
    a.write_json('sector_interference_checks.json', checks)
    print(json.dumps({k: v for k, v in checks.items() if k != 'inputs'}, ensure_ascii=False, indent=2), flush=True)


def write_report(center, grouped, pad, spatial_checks, cell_checks, source_center):
    table = '\n'.join(f'| {name} | {grouped[i,0]:+.6%} | {grouped[i,1]:+.6%} | {grouped[i,2]:+.6%} |'
                      for i, name in enumerate(GROUPS))
    report = table_links(rf"""
{HEADING}

本节对应用户关于“两组之间的相长干涉如何计算并可视化”和“把图 e 按同样方向分组画出”的两条批注。
新增 [Fig.8_sector_interference.png](../10_overview/Fig.8_sector_interference.png)、[Fig.9_unitcell_sector_contributions.png](../10_overview/Fig.9_unitcell_sector_contributions.png)，均有../11_pdf/ 中的同名 PDF、380 dpi；原 Fig.1–7 与 MPH 保留。
本轮只读现有缓存，没有启动 COMSOL、重新求解或修改原始复场。没有增加目录。

### Fig.8：从复角谱计算两组之间的干涉

**数据起点与前序操作。** `01_results/spatial_filter_six_sectors.npz` 的 merged_E_up_angular、merged_H_up_angular、q_m_inv、far_denominator 与 center_gram。
这些数据来自 Fig.5：实际有限模式空气复 E/H → 分离向上波 → NA=1 筛选 → 保留传播相位回溯到 z=100 nm → 在有限表面视窗切成六扇区 → 对切向 E/H 分别 FFT → 线性向上波投影。
空间分区作用于复场，不对 near field intensity 或 EM 强度图做 FFT。原 FFT 每轴补零到 {spatial_checks['fft_size']}；本图使用同一角度网格，无插值增密或新的 FFT。
A=Top+Bottom，B=Left/right four。先相加各扇区复振幅取得 E_A/H_A、E_B/H_B；从输入缓存核对 E_A+E_B=E_full 与 H_A+H_B=H_full。
同一相位基准和 z=100 nm 参考面保持不变，不为各分区选择新相位。

**能量、自项和干涉。** 定义 EM 角谱量 U(E,H)=(epsilon0 sum|E|²+mu0 sum|H|²)/4，并以同一有限表面窗口的 D=U_full(0)>0 作为公共分母。
电、磁场均保留全部三分量；下式中点积只包含相同 Cartesian 分量的交叉项，不把正交 Ex 与 Ey 通道交叉相乘。

\[
I_A=\frac{{\epsilon_0|\mathbf E_A|^2+\mu_0|\mathbf H_A|^2}}{{4D}},\quad
I_B=\frac{{\epsilon_0|\mathbf E_B|^2+\mu_0|\mathbf H_B|^2}}{{4D}},
\]
\[
\kappa_{{AB}}=\frac{{\epsilon_0\operatorname{{Re}}(\mathbf E_A\cdot\mathbf E_B^*)
+\mu_0\operatorname{{Re}}(\mathbf H_A\cdot\mathbf H_B^*)}}{{2D}},\qquad
I_{{\rm full}}=I_A+I_B+\kappa_{{AB}}.
\]

直接从复场乘积计算干涉，再用“总强度−自项和”逐像素独立核对代数恒等式；不从图的颜色或相位零度的实部箭头推断干涉。
正 kappa 表示相长、负 kappa 表示相消，零表示净交叉项为零。负干涉不是负能量。
本次中心为红色相长区，而 theta_y=0 截线的两侧存在蓝色相消区：两组干涉同时增强中心、压低部分旁侧角度的强度。
I_A+I_B 是删除两组间交叉项后的数学对照，不代表当前相干本征模式真的变成两个非相干光源。
本图沿用 Fig.5 的无量纲 EM 角谱定义，不把它称为 W/sr 或积分辐射功率。角度坐标 theta_x=asin(kx/k0)、theta_y=asin(ky/k0)，±10°，没有额外角度 Jacobian。

**逐面板说明。**

| 面板 | 从输入得到的运算和数据 | 物理意义与用途 |
| --- | --- | --- |
| a | 取 A 的复 E/H，求 I_A | 只保留上下两区表面场所得的远场自项；已包含 Top 与 Bottom 之间的干涉 |
| b | 取 B 的复 E/H，求 I_B | 只保留左右四区表面场所得的远场自项；已包含四区内部的干涉 |
| c | 逐角度相加 I_A+I_B | 给出不包含 A–B 交叉项的强度基准，用于比较相干叠加带来的变化 |
| d | 按复场点积直接计算 kappa_AB | 直接显示哪些出射角被两组干涉增强或削弱；红为增强、蓝为削弱 |
| e | 对 E_A+E_B、H_A+H_B 求完整 EM 强度 | 实际同窗完整远场，必须逐像素等于 a+b+d |
| f | 取 theta_y=0 原网格截线，对比 c、d、e | 填色显示总强度相对自项和的增减，直接量化中心亮点中的干涉部分 |

a/b/c/e 全部使用 magma 和线性色标 0–1；d 使用 RdBu_r 和 −1–1，所有图均使用同一个 D。
负色标用于交叉项，不用于实际强度。没有各面板独立取峰值，也没有 Fig.5 中 0–0.04 与 0–1 的混合显示上限。
所有中心数字都来自缓存的严格 kx=ky=0 网格点，而非邻域平均、最大像素或整个亮斑的角积分。

**中心数值与图 e 的贡献口径。** 中心 I_A={center['I_A']:.12g}、I_B={center['I_B']:.12g}、kappa_AB={center['cross']:.12g}，三者相加为 1。
按照两两干涉平均分配的同一约定，A 的归一化贡献 I_A+kappa_AB/2={center['C_A']:.9%}，B 的归一化贡献为 {center['C_B']:.9%}。
所以 A 的合并自项 {center['I_A']:.6%} 不能直接称为 A 对完整中心强度的贡献比例；与 B 的交叉项必须纳入。
Top 和 Bottom 的各自中心自项均约 {center['Top_self']:.12g}，两者之间交叉项为 +{center['Top_Bottom_cross']:.12g}，归一化相干内积实部为 {center['Top_Bottom_correlation_real']:.12g}。
上下相加近乎完全相长。本图 d 的交叉项是“上下组与左右组之间”的干涉，与 Top–Bottom 两个单独扇区之间的干涉不是同一项。

**输出和核对。** `01_results/{FIG8}.npz` 保存所有二维自项、交叉项、相干总量、两组有符号贡献以及电、磁各分量的交叉项。
`80_logs/{FIG8}_center.csv` 保存中心预算。复振幅线性闭合、完整缓存强度、中心 Gram 矩阵贡献三条路径同时核对。
这仍是有限窗口表面孔径分解。相对 Fig.1 直接角谱，继承的中心复振幅差为 {spatial_checks['full_vs_Fig1_center_complex_relative_error']:.6%}，±10°复场差为 {spatial_checks['full_vs_Fig1_complex_relative_error']:.6%}；并非本轮新做的窗口收敛。

### Fig.9：对 Fig.7 e 的完整晶胞做同方向归组

**数据起点。** `01_results/unitcell_center_radiation.npz` 的 centers_um、regions、F_center_V、F_total_V、signed_W_sr，以及现有检查 JSON 和保存 config.json 的真实晶胞几何。
F_j 已由实际 finite slab 的全厚度极化源积分获得，经过 quarter 矢量镜像和跨胞四面体几何裁切修正；本轮不重新做体积分，不使用旧 pointbin 数据。
每胞保持既有复振幅、相位和有符号贡献；没有逐胞平滑、拟合或归一到自身储能。

\[
\mathbf F_j(0)=\frac{{k_0^2}}{{4\pi}}\int_{{V_j}}(\epsilon_r-1)(E_x,E_y)e^{{+ik_0z}}dV,
\quad I_0=\frac{{|\mathbf F_{{\rm total}}|^2}}{{2Z_0}},\quad
C_j=\frac{{\operatorname{{Re}}(\mathbf F_j\cdot\mathbf F_{{\rm total}}^*)}}{{2Z_0}}.
\]

图9所有数值除以同一个 I0={source_center:.12g} W/sr，包含完整 cavity、cladding 和 pad，分母严格沿用 Fig.7 e。
这是严格正上方的 dP/dOmega，不是全 NA 积分功率；负值或超过 100% 来自有符号干涉分配，不是独立能量份额或移除材料后重求解的损耗变化。

**方向权重。** 从 Fig.5 的同一原点出发，按晶胞中心的方位角归入中心方向为 30°、90°、150°、210°、270°、330° 的六个 60° 扇区。
边界射线上的完整晶胞在两区各分 1/2；原点晶胞每区 1/6。对每胞恒有 sum_g w_gj=1。
A 将 Top 和 Bottom 权重相加，B 将其余四区权重相加，然后分别筛选 cavity、cladding。
本规则是完整晶胞的方向归组约定。即使一个非轴心胞跨过射线，也按其中心归组；没有重新裁切体源到扇区。共享射线晶胞的半权不表示已求得两半体积各自的真实复场。
因此本图精确重现此前对逐胞表的方向统计，不冒充严格扇区体积分。

\[
\mathbf F_{{g,r}}=\sum_{{j\in r}}w_{{gj}}\mathbf F_j,\qquad
C_{{g,r}}=\sum_{{j\in r}}w_{{gj}}C_j
=\frac{{\operatorname{{Re}}(\mathbf F_{{g,r}}\cdot\mathbf F_{{\rm total}}^*)}}{{2Z_0}}.
\]

先加复振幅再投影于同一个总场，与直接加权原 C_j 的结果逐组核对。这里使用全结构总场作为相干参照，不能改用组内 F_g 计算自项替代 C_g。

**逐面板说明。** 每个六边形显示的是 w_gj C_j/I0，沿用图 e 的 RdBu_r 色标与 xy 范围。红表示对总中心场的正贡献，蓝表示负贡献。
灰色晶胞不属于当前材料/方向组合，不是负值或零辐射。实线是实际 cavity–cladding 和外侧晶胞区域边界；虚线给出六扇区方向。
图中没有把 pad 涂成晶胞，也没有将胞积分图标为 z=0 局部场分布。

| 面板 | 选择与操作 | 数据意义和用途 |
| --- | --- | --- |
| a | A 方向权重 × cavity 的原 C_j | 显示上下方向 cavity 的中心增强/抵消位置 |
| b | A 方向权重 × cladding 的原 C_j | 显示上下方向 cladding 的中心贡献位置 |
| c | 逐胞 a+b | 上下方向完整晶胞的净贡献；未混入未分方向的 pad |
| d | B 方向权重 × cavity 的原 C_j | 显示左右四区 cavity 的中心贡献 |
| e | B 方向权重 × cladding 的原 C_j | 显示左右四区 cladding 的中心贡献 |
| f | 逐胞 d+e | 左右四区完整晶胞的净贡献 |

每幅图左下角为对应地图数值之和。第一行与第二行逐胞相加恢复 Fig.7 e，而列间有 c=a+b、f=d+e。

| 方向组 | cavity | cladding | 完整晶胞净和，不含 pad |
| --- | ---: | ---: | ---: |
{table}

pad 的整体有符号贡献为 {pad:+.6%}，单独保留；本轮没有足够的 pad 方向子积分，因此不人为均摊到 A/B。
全体晶胞净和加 pad 恢复 100%。每个组的复 F、独立组自项和有符号贡献保存在 `80_logs/{FIG9}_groups.csv`；CSV 同时保存六个单扇区与两组合并。
`01_results/{FIG9}.npz` 保存逐胞六区/两组权重、地图、复组振幅、pad 和公共分母，原逐胞 F_j 可追踪回输入文件。

**两图如何对应及限制。** 两图使用相同方向和同一种交叉项平均分配约定，但表面场区域与材料源区域不一一对应。
传播、NA 筛选及回溯是非局域算子：先把材料源分区再传播，通常不等于先传播完整源再把表面场分区。
Fig.8 的 36%/64% 类表面贡献不应直接要求与 Fig.9 的源贡献相等；归一化基准也分别来自表面重建和体源重建。
体源到独立空气场的中心复振幅差仍为 {cell_checks['cell_sum_vs_air_complex_relative_error']:.6%}，强度差为 {cell_checks['source_vs_air_intensity_relative_error']:.6%}，单胞弱辐射尚未独立收敛。
本图只把现有数据中的方向分布、相长/相消与归一化口径展示清楚，不把分组代数闭合当作新的独立物理验证。

**复现与验证。** `.venv\Scripts\python.exe -B -m scripts.analysis.analyze_dipolar_sector_interference` 只运行本离线补图；加 `--check` 只做同相、反相、相差90°、正交通道及扇区边界/镜像合成检查。
输入 SHA256、实际中心预算、逐胞与分组闭合误差、旧图及 MPH 保留检查见 `80_logs/sector_interference_checks.json`；图面合同写入既有 figure_manifest.json。
""")
    (a.R / REPORT).write_text(report, encoding='utf8')
    path = a.R / 'report.md'
    before = path.read_text(encoding='utf8').split('\n' + HEADING)[0]
    path.write_text(before.rstrip() + '\n' + report, encoding='utf8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    if parser.parse_args().check:
        self_check()
    else:
        run()
