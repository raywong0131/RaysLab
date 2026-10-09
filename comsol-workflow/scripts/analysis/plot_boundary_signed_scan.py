"""Reintegrate exported Hz, validate tracking, and plot signed boundary indicators."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.run_main import run_boundary_analysis as s4
from comsol_workflow.boundary_integrals import signed_indicator, complex_columns


def reciprocal_segments(x, y):
    """Insert a display-only break between finite reciprocal branches of opposite sign."""
    xx, yy = [], []
    for j, (a, b) in enumerate(zip(x, y)):
        if j and (not np.isfinite(b*y[j-1]) or b*y[j-1] < 0):
            xx.append(np.nan)
            yy.append(np.nan)
        xx.append(a)
        yy.append(b)
    return xx, yy


def analyze(root):
    from scripts.run_sweep.run_boundary_surface_signed_scan import validate_fields
    root = Path(root)
    manifest = s4.json_read(root/'99_config/config.json')
    rows = []
    previous = None
    for case in manifest['cases']:
        source = Path(case.get('source_directory', case['directory']))
        metadata = source / ('execution.json' if case['zeta'] == 1 else 'export.json')
        audit = s4.json_read(metadata)
        assert audit['status'] == 'complete'
        data, totals = validate_fields(source, case['zeta'])
        h, w, xy = data['surface_Hz_raw'], data['area_weights_m2'], data['area_xy_m']
        if previous is None:
            overlap = 1.
        else:
            old_xy, old_h, old_w = previous
            np.testing.assert_allclose(xy, old_xy, atol=1e-20, rtol=0)
            np.testing.assert_allclose(w, old_w, atol=1e-30, rtol=1e-12)
            overlap = abs(np.dot(w, old_h.conj()*h))/np.sqrt(np.dot(w, abs(old_h)**2)*np.dot(w, abs(h)**2))
            if overlap < .9:
                raise ValueError(f"Mode tracking overlap below .9 at {case['zeta']}: {overlap}")
        previous = (xy.copy(), h.copy(), w.copy())
        indicator = signed_indicator(totals['QAx'], totals['QBx'])
        np.testing.assert_allclose(indicator['s']**2+indicator['q_perp']**2, indicator['r']**2, rtol=1e-12, atol=1e-25)
        kind = audit.get('field_kind', 'single_eigenmode')
        row = {'zeta': case['zeta'], 'eta': .96, 'mesh': 5, 'z_nm': 100, 'edge_order': 128,
            **indicator, 'f_A_raw': totals['QAx'], 'f_B_raw': totals['QBx'],
            'surface_overlap_previous': float(overlap), 'field_kind': kind,
            'Q': audit.get('Q', np.nan), 'source_directory': str(source),
            'surface_odd_x_error': audit['surface_parity']['odd_x_error'],
            'surface_even_y_error': audit['surface_parity']['even_y_error']}
        rows.append(complex_columns(row))
    table = pd.DataFrame(rows)
    assert table.zeta.tolist() == manifest['zetas']
    table.to_csv(root/'s_zeta.csv', index=False)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), layout='constrained')
    for col, (lo, hi, label) in enumerate(((1, 1.2, 'Full scan'), (1.15, 1.16, 'Near 1.156'))):
        part = table[(table.zeta >= lo) & (table.zeta <= hi)]
        axes[0, col].plot(part.zeta, part.s, 'o-', color='#28749a', ms=4, lw=1)
        xx, yy = reciprocal_segments(part.zeta.to_numpy(), part.inv_s.to_numpy())
        axes[1, col].plot(xx, yy, 'o-', color='#bd543f', ms=4, lw=1)
        axes[0, col].set_title(f'{chr(97+col)}  {label}: signed projection')
        axes[1, col].set_title(f'{chr(99+col)}  {label}: reciprocal')
        for ax in axes[:, col]:
            ax.axhline(0, color='.75', lw=.6)
            ax.axvline(1.156, color='.55', lw=.8, ls=':')
            ax.set(xlim=(lo-.015*(hi-lo), hi+.015*(hi-lo)), xlabel=r'$\zeta$ (1)')
        axes[0, col].set_ylabel(r'$s$ (1)')
        axes[1, col].set_ylabel(r'$1/s$ (1; symlog)')
        axes[1, col].set_yscale('symlog', linthresh=10)
        axes[1, col].set_ylim(min(-10, table.inv_s.min()*1.8), max(10, table.inv_s.max()*1.8))
    fig.suptitle(r'$\eta=0.96$ | $z=100$ nm | $p_y$ at $\Gamma$ | mesh 5 | 128 points per edge', fontsize=12)
    fig.supxlabel('Dots: computed data. Dotted line: nominated position 1.156. Reciprocal lines break at sign changes.', fontsize=9)
    fig.savefig(root/'s_and_inverse.png', dpi=200)
    fig.savefig(root/'s_and_inverse.pdf')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.8, 3.8), layout='constrained')
    ax.plot(table.zeta, table.q_perp, 'o-', ms=3, lw=1, label=r'$q_\perp$: orthogonal projection')
    ax.plot(table.zeta, table.r, 's-', ms=3, lw=1, label=r'$r=|f_A+f_B|/(|f_A|+|f_B|)$')
    ax.axhline(0, color='.75', lw=.6)
    ax.axvline(1.156, color='.55', lw=.8, ls=':')
    ax.set(xlabel=r'$\zeta$ (1)', ylabel='Dimensionless (1)', title='Complex cancellation diagnostics')
    ax.legend(frameon=False)
    fig.savefig(root/'cancellation_diagnostics.png', dpi=200)
    fig.savefig(root/'cancellation_diagnostics.pdf')
    plt.close(fig)
    brackets = [[float(a), float(b)] for a, b, sa, sb in zip(table.zeta[:-1], table.zeta[1:], table.s[:-1], table.s[1:]) if sa*sb < 0]
    selected = table[table.zeta == 1.156].iloc[0]
    s4.write_json(root/'analysis_summary.json', {'status': 'complete', 'point_count': len(table),
        'sign_change_brackets': brackets, 'minimum_adjacent_overlap': float(table.surface_overlap_previous.min()),
        'zeta1_156': {k: float(selected[k]) for k in ('s', 'inv_s', 'q_perp', 'r')},
        'mixed_state_zetas': table.loc[table.field_kind != 'single_eigenmode', 'zeta'].tolist(),
        'normal': 'dielectric to air', 'phase_convention': 'exp(-i omega t)',
        'reciprocal_at_exact_zero': 'undefined; NaN in CSV; no epsilon or clipping',
        'raw_integral_units': 'A for physical Hz in A/m; eigenmode amplitude is arbitrary'})
    report = f'''# ζ 扫描结果

固定 η=0.96、b=245 nm、a=820 nm、板厚200 nm、折射率3.3、Γ点用户py模式，mesh=5。
取z=100 nm板层侧的Hz复数场，对六个三角孔的18条边分别作128点高斯积分；法向由介质指向空气。
A包含孔1、4，B包含孔2、3、5、6；各组完整求和，不重复乘孔数。

s = Re[(fA+fB)conj(fA)] / [|fA|(|fA|+|fB|)]。
q采用相同表达式的虚部，r=|fA+fB|/(|fA|+|fB|)，验证了r²=s²+q²。
s、1/s、q、r均无量纲；全局复相位及场幅度缩放不影响结果。原始fA、fB的单位为A，但本征模幅度任意。

共{len(table)}个ζ点，{manifest['new_solve_count']}个新Γ求解，{len(table)-manifest['new_solve_count']}个已有同设置导出通过几何、积分点、网格记录和来源哈希核验后复用。
全区间步长0.01，1.150至1.160加密至0.001。积分CSV从原始128点场重新积分生成。
相邻表面场重叠最低{table.surface_overlap_previous.min():.8f}。

ζ=1.156时：s={selected.s:.12g}，1/s={selected.inv_s:.12g}，q={selected.q_perp:.12g}，r={selected.r:.12g}。
实际采样点的s变号区间：{brackets}。它们仅是离散变号区间，不是已收敛的精确零点。
s=0只表示A方向的投影为零；还需q=0才能称复积分完全抵消。没有强制1.156为零。

ζ={table.loc[table.field_kind != 'single_eigenmode', 'zeta'].tolist()}使用近简并p子空间的常数复系数组合辨认py对称性，并非严格单一本征模；原始两模频率与Q保留在来源记录，不赋予组合场一个Q。
mesh=5未做网格收敛验证，细微变号或倒数大值不能直接认定为物理发散。
z方向选择板层侧并不消除孔边界横向场迹的采样敏感性。

图中只有真实点和相邻点直线，无拟合、平滑、epsilon或人为零点平移。1/s使用保留正负的对称对数坐标（|1/s|≤10区间线性），跨s变号处断线；恰好s=0时CSV为空值，表示未定义。
主图为s_and_inverse.png/PDF（上排s，下排1/s，左全区间，右局部）；cancellation_diagnostics.png/PDF给出q和r。
数据为s_zeta.csv，所有案例及来源哈希见99_config/config.json。
'''
    (root/'report.md').write_text(report, encoding='utf-8')
    print(table[['zeta', 's', 'inv_s', 'q_perp', 'r', 'surface_overlap_previous']].to_string(index=False), flush=True)
    return table


if __name__ == '__main__':
    import sys
    analyze(Path(sys.argv[1]))
