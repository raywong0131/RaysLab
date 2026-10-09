"""Actual finite-slab cells: coherent radiation into the strict normal direction.

Offline only. Full-thickness native E/H exports, no periodic replacement fields.
Run with --check for synthetic checks, --plot-only to redraw saved cell integrals.
"""
from __future__ import annotations

from comsol_workflow.output_paths import table_links, table_path
import argparse
import json

import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.path import Path as PolygonPath
from scipy.spatial import cKDTree

from comsol_workflow.hex_lattice_utils import unit_cell_corners
from scripts.analysis import analyze_dipolar_complete as a
from scripts.analysis.analyze_dipolar_singularity import digest

STEM = "unitcell_center_radiation"
FIGURE = "Fig.7_" + STEM
HEADING = "## Fig.7：实际 finite slab 的逐晶胞中心辐射"


def grouped_sum(ids, values, count):
    """Sum vector quantities without constructing a cells-by-points matrix."""
    values = np.asarray(values)
    result = np.empty((count, values.shape[1]), dtype=values.dtype)
    for p, value in enumerate(values.T):
        result[:, p] = np.bincount(ids, weights=value.real, minlength=count)
        if np.iscomplexobj(values):
            result[:, p] += 1j * np.bincount(ids, weights=value.imag, minlength=count)
    return result


def reflected_center_integrals(xyz_um, electric, weights, epsilon_r, ids,
                               mirrors, k0):
    """Exact z retardation and the saved mode's x/y/z electric-field parities.

    mirrors maps each full-cell ID to its reflected cell ID, with pad last.
    Inputs occupy x,y,z > 0; points on symmetry planes are not volume nodes.
    """
    factor = k0**2 / (4 * np.pi)
    source = electric[:, :2] * (2 * weights * (epsilon_r - 1)
                                * np.cos(k0 * xyz_um[:, 2] * 1e-6))[:, None]
    partial = grouped_sum(ids, source, mirrors.shape[1]) * factor
    result = np.zeros_like(partial)
    for mapping, signs in zip(mirrors, [(1, 1), (-1, 1), (-1, 1), (1, 1)]):
        # A reflected cell on an axis can map to itself: add, never overwrite.
        result[mapping] += partial * signs
    return result


def budgets(amplitudes):
    """Self terms and symmetric allocation of interference; not unique losses."""
    total = amplitudes.sum(axis=0)
    own = abs(amplitudes)**2 / (2 * a.Z0)
    signed = np.real(amplitudes * total.conj()) / (2 * a.Z0)
    return total, own, signed


def self_check():
    # All quadrants, cells crossing mirror planes, and a separate pad bin.
    centers = np.array([(x, y) for y in [-1., 0., 1.] for x in [-1., 0., 1.]])
    tree = cKDTree(centers)
    mirrors = np.array([np.r_[tree.query(centers * (sx, sy))[1], len(centers)]
                        for sx, sy in [(1, 1), (-1, 1), (1, -1), (-1, -1)]])
    rng = np.random.default_rng(31)
    xyz = rng.uniform([.01, .01, .003], [1.2, 1.2, .1], (37, 3))
    electric = rng.normal(size=(37, 3)) + 1j * rng.normal(size=(37, 3))
    weights = rng.uniform(.1, 1., 37) * 1e-20
    eps = np.full(37, 10.89)
    ids = tree.query(xyz[:, :2])[1]
    ids[-2:] = len(centers)
    k0 = 4.1e6
    computed = reflected_center_integrals(xyz, electric, weights, eps, ids, mirrors, k0)
    expected = np.zeros_like(computed)
    for m, (sx, sy) in enumerate([(1, 1), (-1, 1), (1, -1), (-1, -1)]):
        for sz in [-1, 1]:
            source = electric[:, :2] * (sx * sy, 1)
            source *= (weights * (eps - 1) * np.exp(1j * k0 * sz * xyz[:, 2] * 1e-6))[:, None]
            expected += grouped_sum(mirrors[m, ids], source, len(centers) + 1) * k0**2 / (4 * np.pi)
    np.testing.assert_allclose(computed, expected, rtol=2e-13, atol=1e-20)
    np.testing.assert_allclose(computed.sum(0)[0], 0., atol=1e-20)
    signals = np.array([[1., 2j], [-2., 1j], [.2j, -4j]])
    total, own, signed = budgets(signals)
    np.testing.assert_allclose(signed.sum(0), abs(total)**2 / (2 * a.Z0))
    cross = sum(np.real(signals[i] * signals[j].conj()) / a.Z0
                for i in range(len(signals)) for j in range(i))
    np.testing.assert_allclose(own.sum(0) + cross, signed.sum(0))
    assert np.any(signed < 0)
    print("PASS: full-volume direct sum, mirror-plane cells, vector parity, interference allocation.", flush=True)


def compute():
    cfg_path = a.CONFIG / "config.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf8"))
    points = cfg["bulk_points"] + cfg["cladding_points"]
    nc = len(cfg["bulk_points"])
    count = len(points)
    centers = np.array([[p["x"], p["y"]] for p in points])
    tree = cKDTree(centers)
    poly = PolygonPath(unit_cell_corners(float(cfg["a"])))
    mirrors = []
    for sx, sy in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
        distances, mapping = tree.query(centers * (sx, sy))
        assert distances.max() < 1e-10 and len(np.unique(mapping)) == count
        mirrors.append(np.r_[mapping, count])
    mirrors = np.array(mirrors)
    amplitudes = np.zeros((count + 1, 2), complex)
    energy = np.zeros((count + 1, 6))
    source_l1 = np.zeros((count + 1, 2))
    volume = np.zeros(count + 1)
    node_count = np.zeros(count + 1, dtype=np.int64)
    independent = np.zeros(2, complex)
    sources = [cfg_path, a.L / "finite_model_identity.json",
               a.L / "finite_boundary_conditions.json"]
    near_boundary = 0
    total_points = 0
    for region in ["cavity", "cladding", "pad"]:
        path = a.D / f"finite_native_{region}_EH.npz"
        sources.append(path)
        data = a.load(path.name)
        xyz, w = data["xyz_um"], data["weights_octant_m3"]
        electric, magnetic, eps = data["E_V_m"], data["H_A_m"], data["epsilon_r"]
        assert np.all(np.isfinite(electric)) and np.all(np.isfinite(magnetic))
        assert np.all(w > 0) and np.all(xyz > 0)
        assert np.array_equal(eps > 1, data["material_mask"])
        _, ids = tree.query(xyz[:, :2])
        local = xyz[:, :2] - centers[ids]
        inside = poly.contains_points(local, radius=1e-10)
        if region == "pad":
            assert not np.any(inside)
            ids[:] = count
        else:
            assert np.all(inside)
            assert np.all(ids < nc) if region == "cavity" else np.all(ids >= nc)
            # A count only: this does not estimate an integration error.
            near_boundary += np.count_nonzero(~poly.contains_points(local, radius=-.002))
        amplitudes += reflected_center_integrals(xyz, electric, w, eps, ids, mirrors, a.K0)
        density = np.column_stack([a.epsilon_0 * eps[:, None] * abs(electric)**2,
                                   a.mu_0 * abs(magnetic)**2]) / 4
        partial_energy = grouped_sum(ids, density * (2 * w[:, None]), count + 1)
        partial_l1 = grouped_sum(ids, abs(electric[:, :2]) * (2 * w * (eps - 1))[:, None], count + 1)
        partial_volume = np.bincount(ids, weights=2 * w, minlength=count + 1)
        partial_count = np.bincount(ids, minlength=count + 1)
        for mapping in mirrors:
            energy[mapping] += partial_energy
            source_l1[mapping] += partial_l1
            volume[mapping] += partial_volume
            node_count[mapping] += 2 * partial_count
        independent[1] += np.sum(electric[:, 1] * w * (eps - 1)
                                  * np.cos(a.K0 * xyz[:, 2] * 1e-6)) * 8 * a.K0**2 / (4 * np.pi)
        total_points += len(w)
        print(f"{region}: integrated {len(w):,} native volume nodes", flush=True)

    assert np.all(node_count > 0)
    clipped_path=a.D/"unitcell_clipped_integrals.npz"
    if not clipped_path.exists():
        raise FileNotFoundError("Run scripts.analysis.audit_dipolar_cell_quadrature --all-cells to build cell-clipped integrals first")
    clipped=a.load(clipped_path.name)
    assert str(clipped["config_sha256"])==digest(cfg_path)
    np.testing.assert_allclose(clipped["cell_centers_um"],centers,rtol=0,atol=1e-12)
    np.testing.assert_allclose(clipped["baseline_F_total_V"],amplitudes.sum(0),rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose(clipped["baseline_U_components_J"],energy.sum(0),rtol=1e-11)
    amplitudes=clipped["F_center_V"]
    energy=clipped["energy_components_J"]
    source_l1=clipped["source_L1_V_m2"]
    volume=clipped["volume_m3"]
    exact_cell_volume=float(a.load("cell_quadrature_audit_all.npz")["exact_cell_volume_m3"])
    np.testing.assert_allclose(volume[:-1],exact_cell_volume,rtol=1e-7,atol=0)
    sources.extend([clipped_path,a.L/"cell_quadrature_audit_all.json"])
    clip_audit=json.loads((a.L/"cell_quadrature_audit_all.json").read_text(encoding="utf8"))
    total, own, signed = budgets(amplitudes)
    intensity = float(np.sum(abs(total)**2) / (2 * a.Z0))
    uref = a.native()["U_ref_J"]
    cancel = np.divide(abs(amplitudes) * 4 * np.pi / a.K0**2, source_l1,
                       out=np.zeros_like(source_l1), where=source_l1 > 0)
    assert np.all(cancel <= 1 + 1e-12)
    source_files = [a.D / "full_mesh_radiation.npz", a.D / "direct_radiative_fields.npz",
                    a.D / "mesh_energy_regions.npz", a.L / "native_full_integrals.json"]
    sources.extend(source_files)
    mesh = a.load("full_mesh_radiation.npz")["F_V"]
    air = a.load("direct_radiative_fields.npz")["F_V"]
    mid = len(mesh) // 2
    mesh_center, air_center = mesh[mid, mid, :2], air[mid, mid, :2]
    native = np.array(a.native()["full_material_center_F_V_pairs"])
    native_center = native[:, 0] + 1j * native[:, 1]
    regions = np.array(["cavity"] * nc + ["cladding"] * (count - nc) + ["pad"])
    shells = np.array([p["shell"] for p in points] + [-1])
    aggregates = []
    shell_rows = []
    for label, mask in [(r, regions == r) for r in ["cavity", "cladding", "pad"]] + [("all", np.ones(count + 1, bool))]:
        field = amplitudes[mask].sum(0)
        aggregates.append(dict(region=label,cell_count=int(mask[:-1].sum()),
            U_EM_J=float(energy[mask].sum()), U_EM_over_Uref=float(energy[mask].sum() / uref),
            Fx_re_V=float(field[0].real), Fx_im_V=float(field[0].imag),
            Fy_re_V=float(field[1].real), Fy_im_V=float(field[1].imag),
            group_self_W_sr=float(np.sum(abs(field)**2) / (2 * a.Z0)),
            sum_unit_self_W_sr=float(own[mask].sum()),
            signed_contribution_W_sr=float(signed[mask].sum()),
            signed_fraction_of_source_center=float(signed[mask].sum() / intensity)))
    for shell in np.unique(shells):
        mask = shells == shell
        f = amplitudes[mask].sum(0)
        shell_rows.append(dict(shell=int(shell),region="pad" if shell < 0 else str(regions[mask][0]),
            count=int(mask.sum()),U_EM_J=float(energy[mask].sum()),
            Fy_re_V=float(f[1].real),Fy_im_V=float(f[1].imag),
            sum_self_W_sr=float(own[mask].sum()),
            signed_contribution_W_sr=float(signed[mask].sum()),
            signed_fraction_of_source_center=float(signed[mask].sum() / intensity)))
    regional_energy = np.array([energy[regions == r].sum(0) for r in ["cavity", "cladding", "pad"]])
    expected_energy = a.load("mesh_energy_regions.npz")["energy_components_J"]
    checks = dict(COMSOL_started=False, angular_scope="strict normal direction, kx=ky=0 only",
        source_model="actual induced polarization in homogeneous air Green function; real part of eigenfrequency",
        frequency_Hz=a.FREQ,cavity_cells=nc,cladding_cells=count-nc,pad_is_unit_cell=False,
        units={"F_center":"V","source_L1":"V*m^2","energy":"J","intensity":"W/sr","length":"um"},
        integration_method="exact hexagonal cell clipping; affine source/density reconstructed from the same four native Gauss samples on cut tetrahedra",
        cell_volume_relative_error=clip_audit["clipped_max_cell_volume_error"],
        pointbin_volume_relative_error=clip_audit["old_max_cell_volume_error"],
        original_octant_nodes=total_points,original_nodes_within_about_1nm_of_cell_boundary=int(near_boundary),
        volume_m3=float(volume.sum()), U_ref_J=uref,U_mesh_J=float(energy.sum()),
        stored_energy_relative_error_to_COMSOL=float(abs(energy.sum() / uref - 1)),
        regional_energy_relative_error=float(np.linalg.norm(regional_energy-expected_energy) / np.linalg.norm(expected_energy)),
        cell_sum_vs_unpartitioned_relative_error=float(np.linalg.norm(total-independent) / np.linalg.norm(independent)),
        cell_sum_vs_existing_mesh_relative_error=float(np.linalg.norm(total-mesh_center) / np.linalg.norm(mesh_center)),
        cell_sum_vs_COMSOL_native_relative_error=float(np.linalg.norm(total-native_center) / np.linalg.norm(native_center)),
        cell_sum_vs_air_complex_relative_error=float(np.linalg.norm(total-air_center) / np.linalg.norm(air_center)),
        source_center_W_sr=intensity,air_center_W_sr=float(np.sum(abs(air_center)**2) / (2*a.Z0)),
        source_vs_air_intensity_relative_error=float(abs(np.sum(abs(total)**2) / np.sum(abs(air_center)**2) - 1)),
        sum_signed_contribution_relative_error=float(abs(signed.sum() / intensity - 1)),
        sum_unit_self_W_sr=float(own.sum()),pair_interference_W_sr=float(intensity-own.sum()),
        cell_quadrature_independently_converged=False,
        note="Point-bin boundary error corrected by conservative geometric clipping. Affine reconstruction on cut tetrahedra still requires independent source convergence.",
        region_budget=aggregates,
        inputs={str(p):{"sha256":digest(p),"bytes":p.stat().st_size} for p in sources})
    checks["independent_center_closure_passed"] = (checks["cell_sum_vs_air_complex_relative_error"] <= .02
                                                  and checks["source_vs_air_intensity_relative_error"] <= .05)
    assert checks["cell_sum_vs_unpartitioned_relative_error"] < 1e-11
    assert checks["cell_sum_vs_existing_mesh_relative_error"] < 1e-9
    assert checks["cell_volume_relative_error"] < 1e-7
    np.testing.assert_allclose(energy.sum(0),expected_energy.sum(0),rtol=1e-11)
    assert checks["sum_signed_contribution_relative_error"] < 1e-12
    rows = []
    phase_ref = np.angle(total[1])
    for j in range(count + 1):
        p = points[j] if j < count else {}
        row = dict(cell_id=j if j < count else "pad",region=regions[j],i=p.get("i"),j=p.get("j"),
                   shell=p.get("shell"),x_um=p.get("x"),y_um=p.get("y"),
                   volume_m3=volume[j],original_pointbin_node_count=node_count[j],
                   U_EM_J=energy[j].sum(),U_EM_over_Uref=energy[j].sum()/uref)
        for q, channel in enumerate(["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]):
            row[f"U_{channel}_J"] = energy[j,q]
        for q, channel in enumerate(["Ex", "Ey"]):
            f = amplitudes[j,q]
            row.update({f"F_{channel}_re_V":f.real,f"F_{channel}_im_V":f.imag,
                        f"F_{channel}_phase_relative_to_total_Ey_deg":(float(np.rad2deg(np.angle(f*np.exp(-1j*phase_ref)))) if abs(f)>0 else None),
                        f"self_{channel}_W_sr":own[j,q],f"signed_{channel}_W_sr":signed[j,q],
                        f"cancellation_{channel}":cancel[j,q]})
        row.update(self_total_W_sr=own[j].sum(),signed_total_W_sr=signed[j].sum(),
                   signed_fraction_of_source_center=signed[j].sum()/intensity,
                   self_per_cell_U_s_inv_sr=own[j].sum()/energy[j].sum())
        rows.append(row)
    pd.DataFrame(rows).to_csv(table_path(a.L, f"{STEM}_cells.csv"),index=False)
    pd.DataFrame(aggregates).to_csv(table_path(a.L, f"{STEM}_regions.csv"),index=False)
    pd.DataFrame(shell_rows).to_csv(table_path(a.L, f"{STEM}_shells.csv"),index=False)
    np.savez_compressed(a.D / f"{STEM}.npz",centers_um=centers,regions=regions,shells=shells,
        F_center_V=amplitudes,energy_components_J=energy,source_L1_V_m2=source_l1,
        cell_cancellation=cancel,self_W_sr=own,signed_W_sr=signed,volume_m3=volume,
        original_pointbin_node_count=node_count,F_total_V=total,F_air_center_V=air_center,
        U_ref_J=uref,source_center_W_sr=intensity,frequency_Hz=a.FREQ,
        angular_scope="kx=ky=0; upward only",integration_method="geometric cell clipping; affine interpolation on cut FEM tetrahedra",normalization="source center includes cells and pad")
    a.write_json(f"{STEM}_checks.json",checks)
    print(json.dumps({k:v for k,v in checks.items() if k != "inputs"},ensure_ascii=False,indent=2),flush=True)


def plot_report():
    data = a.load(f"{STEM}.npz")
    checks = json.loads((a.L / f"{STEM}_checks.json").read_text(encoding="utf8"))
    cfg = json.loads((a.CONFIG / "config.json").read_text(encoding="utf8"))
    if "cell_volume_relative_error" not in checks or "integration_method" not in data:
        raise RuntimeError("Recompute the cell budget with validated clipped integrals before plotting")
    centers = data["centers_um"]
    corners = centers[:,None,:] + unit_cell_corners(float(cfg["a"]))[None,:,:]
    u = data["energy_components_J"][:-1].sum(1) / float(data["U_ref_J"])
    self_terms = data["self_W_sr"][:-1] / float(data["source_center_W_sr"])
    signed = data["signed_W_sr"][:-1].sum(1) / float(data["source_center_W_sr"])
    phase = np.angle(data["F_center_V"][:-1,1] * data["F_total_V"][1].conj())
    boundary = np.asarray(cfg["inner_boundary"])
    boundary = np.vstack([boundary,boundary[0]])
    boundary_full = np.asarray(cfg["finite_boundary"])
    bound = float(np.ceil(abs(boundary_full).max()))
    top_u = 10.**np.ceil(np.log10(u.max()))
    top_self = 10.**np.ceil(np.log10(self_terms.max()))
    top_signed = float(abs(signed).max())
    specs = [
        (u,"a  Cell EM stored energy",r"$U_j/U_{\mathrm{ref}}$ (1)","magma",a.LogNorm(top_u*1e-5,top_u)),
        (self_terms[:,0],"b  Ex: cell self radiation",r"$I_{j,x}/I_{\mathrm{source}}(0)$ (1)","magma",a.LogNorm(top_self*1e-5,top_self)),
        (self_terms[:,1],"c  Ey: cell self radiation",r"$I_{j,y}/I_{\mathrm{source}}(0)$ (1)","magma",a.LogNorm(top_self*1e-5,top_self)),
        (phase,"d  Ey: phase relative to total",r"$\arg(F_{j,y}F_{\mathrm{total},y}^{*})$ (rad)","twilight_shifted",Normalize(-np.pi,np.pi)),
        (signed,"e  Signed coherent contribution",r"$C_j/I_{\mathrm{source}}(0)$ (1)","RdBu_r",Normalize(-top_signed,top_signed)),
        (data["cell_cancellation"][:-1,1],"f  Ey: intracell cancellation",r"$|\int_j\chi E_y e^{ik_0z}dV|/\int_j|\chi E_y|dV$ (1)","magma",a.LogNorm(1e-5,1)),
    ]
    a.plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    fig,axes = a.plt.subplots(2,3,figsize=(17,11.5),layout="constrained")
    import matplotlib.patheffects as pe
    for ax,(values,title,label,cmap,norm) in zip(axes.flat,specs):
        display_values = np.maximum(values,norm.vmin) if isinstance(norm,a.LogNorm) else values
        cells = PolyCollection(corners,array=display_values,cmap=cmap,norm=norm,
                               edgecolors="none",rasterized=True)
        ax.add_collection(cells)
        ax.plot(boundary[:,0],boundary[:,1],color="white",lw=.7,
                path_effects=[pe.Stroke(linewidth=1.3,foreground="#444444"),pe.Normal()])
        ax.set(title=title,xlabel=r"$x$ (µm)",ylabel=r"$y$ (µm)",
               xlim=(-bound,bound),ylim=(-bound,bound),aspect="equal")
        ax.set_facecolor("#ececec")
        bar = a.colorbar(fig,ax,cells,label)
        if isinstance(norm,a.LogNorm):
            lo,hi = int(round(np.log10(norm.vmin))),int(round(np.log10(norm.vmax)))
            bar.set_ticks(10.**np.arange(lo,hi+1))
        if cmap == "twilight_shifted":
            bar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            bar.set_ticklabels([r"$-\pi$",r"$-\pi/2$","0",r"$\pi/2$",r"$\pi$"])
        assert ax.get_aspect() == 1.0
    fig.suptitle("Actual finite-slab cells | Radiation into the strict normal direction\n"
                 f"{checks['cavity_cells']} cavity + {checks['cladding_cells']} cladding cells | Full-thickness polarization-source integrals",fontsize=14)
    fig.supxlabel("Each hexagon is a cell integral, not a local field pixel. White outline: cavity-cladding interface.\n"
                  "Pad retained in the budget. Source-to-air closure remains unverified.",fontsize=10)
    a.save(fig,FIGURE,dpi=a.FIGURE_DPI)
    table = "\n".join(f"| {r['region']} | {r['U_EM_over_Uref']:.6%} | {r['Fy_re_V']:.6g} {r['Fy_im_V']:+.6g}i | {r['group_self_W_sr']:.6g} | {r['signed_contribution_W_sr']:.6g} | {r['signed_fraction_of_source_center']:.6%} |"
                      for r in checks["region_budget"])
    report = table_links(rf"""
{HEADING}

文件：[Fig.7_unitcell_center_radiation.png](../10_overview/Fig.7_unitcell_center_radiation.png)，../11_pdf/ 中的同名 PDF，380 dpi；2×3 六边形晶胞统计图。
本节只研究严格法向 kx=ky=0 的中心辐射，不把中心数值称为全 NA 积分功率。本轮没有启动 COMSOL。

**问题与辐射单元。** 使用真实 finite mode19（用户 py 基模）在 slab 中的场，每个实际 cavity/cladding 晶胞是一组诱导极化源。
1141 个 cavity 胞与 3540 个 cladding 胞全部保留；外侧无孔 pad 单独保留为补充区域，不冒充晶胞。
单胞源来自完整结构的自洽场，因此已经受到相邻晶胞及界面的影响。分离一个源积分，不是取出这个晶胞重新求解；不能解释为隔离晶胞的独立损耗率。

**原始数据与几何。** 读取 01_results/finite_native_cavity_EH.npz、finite_native_cladding_EH.npz、finite_native_pad_EH.npz。
它们保存原生四面体求积节点 xyz_um、weights_octant_m3、E_V_m、H_A_m、epsilon_r 和 material_mask；正权重四点规则为二次精度。
本次处理原始八分之一薄板的 {checks['original_octant_nodes']:,} 个体节点，覆盖整个厚度的正半部；200 nm 全厚度通过模型的 z 奇偶恢复。
z=0 场图不能代替这个体积分。本轮不读取周期场或 Bloch 拟合系数，也不使用 NA 滤波后的近场作为辐射源。
晶胞中心、i/j、shell、cavity/cladding 身份及 a 从保存 config.json 读取。当前几何没有 cladding 位移；六边形取项目 unit_cell_corners(a)。
先用实际六边形棱柱裁切跨胞四面体，对其四个原生 Gauss 样本确定的仿射源/能量密度进行几何精确积分；完整落在单胞内的四面体保留原规则。
这个修正改变的是实际积分域及权重，原 FEM 场样本和相位保持不变。pad 是同一仿射体场在晶胞并集外的守恒余量，未对独立空气场做幅值校准。
旧“按求积点归属分配整份权重”结果保存在 pointbin 快照；其 cavity 散斑已确认主要来自跨胞四面体的权重误分配。
源文件完整路径、SHA256、大小、频率及边界条件记录见 ../80_logs/{STEM}_checks.json。

**从实际场到辐射源。** 非磁性介质采用空气为参考背景，\(\chi=\epsilon_r-1\)，COMSOL 时间约定为 \(e^{{+i\omega t}}\)：

\[
\mathbf P(\mathbf r)=\epsilon_0\chi(\mathbf r)\mathbf E(\mathbf r),\qquad
\mathbf J_p(\mathbf r)=+i\omega\mathbf P(\mathbf r).
\]

孔内空气的 chi=0，自动没有极化源；储能计算仍包含孔内 E/H。没有将 u_EM、|E|² 或 z=0 场图拿来做辐射 Fourier 变换。
频率使用保存本征频率的实部 {a.FREQ/1e12:.12g} THz，与现有空气/体源计算一致；这仍是实频 Green 算子作用于复本征场的近似。

**一般方向与本轮中心积分。** 若远处场写为 \(\mathbf E(R\hat{{\mathbf n}})=\mathbf F(\hat{{\mathbf n}})e^{{-ik_0R}}/R\)，
在均匀空气 Green 模型中，第 j 个体单元的远场系数为

\[
\mathbf F_j(\hat{{\mathbf n}})=\frac{{k_0^2}}{{4\pi}}
(\mathbf I-\hat{{\mathbf n}}\hat{{\mathbf n}}^T)
\int_{{V_j}}\chi\mathbf E(\mathbf r)e^{{+ik_0\hat{{\mathbf n}}\cdot\mathbf r}}\,dV.
\]

本轮仅求 \(\hat{{\mathbf n}}=\hat{{\mathbf z}}\)，因此 Fz=0，Ex/Ey 是两个正交通道，保留准确的厚度传播相位：

\[
F_{{j,\alpha}}(0)=\frac{{k_0^2}}{{4\pi}}
\int_{{V_j}}\chi E_\alpha e^{{+ik_0z}}\,dV,\quad\alpha=x,y.
\]

本模式 Ex/Ey 关于 z 为偶，正半厚度贡献乘 2cos(k0 z)；关于 x、y，Ex 均为奇，Ey 均为偶。
先积分正象限内的部分晶胞，再将其映射到实际镜像晶胞；跨轴晶胞累加镜像部分，不能每胞一律乘八。
z 相位没有近似成 1，胞内场没有替换成常数或点偶极子。法向没有横向传播相位差，但原始复场的胞间相位全部保留。
全域 Ex 的相消部分由保存的 quarter 边界条件及镜像构造保证，不能作为独立的全结构对称性验证。
本项目的 x-axis PEC 指 y=0，y-axis PMC 指 x=0，不能把轴名称误读成法向坐标。跨 x=0 或 y=0 的完整胞 Ex 源积分为零；Ey 为偶，不会因此产生暗十字。
H 是轴矢量：x 镜像的 H 符号为 (+,-,-)，y 镜像为 (+,-,+)，z 镜像为 (-,-,+)。所有操作只乘实数 ±1，不取共轭，不为各胞重新选相位。
实际 z=0 原始截面中，部分轴上导出迹值并非严格满足理论奇分量零值；本轮体积分节点均严格离轴，没有读取这些迹值。
此项作为额外的边界采样/有限元迹值待核查事项记录在 ../80_logs/cell_quarter_symmetry_audit.json，不将图像离轴镜像完全一致误称为独立物理验证。

**完整电磁能量与强度。** 每胞完整储能为

\[
U_j=\frac14\int_{{V_j}}\left(\epsilon_0\epsilon_r|\mathbf E|^2+\mu_0|\mathbf H|^2\right)dV.
\]

六分量分别保存。远场 H 系数为 \(\hat{{\mathbf z}}\times\mathbf F_j/Z_0\)，空气中电、磁能相等；
完整 EM 能量系数 \(\lim R^2u_{{EM}}=\epsilon_0|\mathbf F|^2/2\)，乘 c 得到方向辐射功率密度：
\(I=dP/d\Omega=|\mathbf F|^2/(2Z_0)\)，单位 W/sr。
这已包含完整 EM 能流，不能再叠加一份“磁功率”。它与前图的 EM 角谱权重定义不同，当前 CSV 的 W/sr 明确是 dP/dΩ。
本征场振幅任意，因此同时保存 U_ref 归一化量及 I_j/U_j；后者的单位为 s^-1 sr^-1，是单胞源自项对储能的比值，不是独立腔 Q。

**怎样定义每胞贡献。** 原始可追踪量为复振幅 F_j。每胞自项 \(I_j=|F_{{j,x}}|^2/(2Z_0)+|F_{{j,y}}|^2/(2Z_0)\) 非负，但不能相加得到总强度。
令 \(\mathbf F_T=\sum_j\mathbf F_j\)，包含 pad，并采用对两两干涉各分配一半的对称约定：

\[
C_j=\frac{{\operatorname{{Re}}(\mathbf F_j\cdot\mathbf F_T^*)}}{{2Z_0}},\quad
\sum_jC_j=I_T,\quad
I_T=\sum_j I_j+\sum_{{j<l}}\frac{{\operatorname{{Re}}(\mathbf F_j\cdot\mathbf F_l^*)}}{{Z_0}}.
\]

Cj>0 表示该胞在这个相干态中增强中心，Cj<0 表示抵消。负值不是负能量；它也不是移除实际材料后重新求解的功率变化。
这种贡献划分依赖明确的干涉分配约定，不声称存在唯一、正定、可加的单胞辐射损耗。Ex/Ey 分别核算，没有把正交通道相互干涉。
逐胞复振幅保存在 NPZ/CSV，可直接恢复任意一对的交叉项，避免输出无必要的数千万行配对表。

**六个面板的数据、操作与用途。** 所有面板中，一个六边形代表一个体积分值；不是点图，也不是局部场像素或 z=0 场分布。

| 面板 | 原始数据及操作 | 输出与归一化 | 要回答的问题 |
| --- | --- | --- | --- |
| a | 全厚度实际 E/H、epsilon_r 与体权重，积分六分量能量 | Uj/U_ref；U_ref={float(data['U_ref_J']):.12g} J；log [{top_u*1e-5:g},{top_u:g}] | 储能在哪些胞，是否等同于辐射贡献的位置 |
| b | 每胞实际 Ex 极化源，保留 z 相位积分后取模平方/(2Z0) | Ijx/I_source(0)；与 c 同一 log [{top_self*1e-5:g},{top_self:g}] | 单胞可有 Ex 辐射，但全结构是否通过干涉抵消 |
| c | 与 b 相同，对 Ey 处理 | Ijy/I_source(0)，共用 b 分母和色标 | 每胞 Ey 源的法向辐射自项有多大 |
| d | 每胞复 Fy 乘总 Fy 的共轭再取 argument | 相位 rad，[-pi,pi]；循环色图 | 胞间振幅为什么增强或相消；极弱振幅的相位须结合 c 阅读 |
| e | 复 F_j 与总 F_T 的内积取实部/(2Z0) | Cj/I_source(0)，线性 ±{top_signed:.9g}，RdBu_r | 哪些胞增强、哪些胞抵消中心；胞加 pad 后和为 1 |
| f | abs[integral chi Ey exp(ik0z)dV] / integral abs(chi Ey)dV | 无量纲 [0,1]，log [1e-5,1] | 每胞内部向正上方辐射抵消的程度；越小抵消越充分 |

所有 log 色条不超过五个数量级；低于色标的值只在显示上使用底色，原数据不截断。白线是实际 cavity/cladding 分界。
图中共用的辐射分母为当前体源模型总中心 {checks['source_center_W_sr']:.12g} W/sr，包含 pad，未把每胞独立归一到 1。
pad 没有画成晶胞，其贡献完整列入下面的预算。胞内相消指标不是 Q/BIC 判据；零附近还会受到求积误差影响。

**当前区域预算。** 分组自项是先相加该区域所有胞的复振幅再平方；与“逐胞自项之和”不同。C 的分母及相位参照均包含 pad。

| 区域 | EM 储能/U_ref | 合并 Fy (V) | 分组自项 (W/sr) | 有符号贡献 C (W/sr) | C/I_source |
| --- | ---: | --- | ---: | ---: | ---: |
{table}

这些百分数是有符号相干贡献，可小于 0 或超过 100%，不是光子来源概率。pad 自项仍是整圈 pad 的自项，不能当作一个普通晶胞的自项。
本结果让每个实际胞的内部相消、胞间相位与中心辐射建立联系，不预设贡献一定只集中于 interface。

**验证与适用范围。** 逐胞加 pad 相对不分区体积分的复振幅误差 {checks['cell_sum_vs_unpartitioned_relative_error']:.3e}；
相对既有原生 mesh 中心误差 {checks['cell_sum_vs_existing_mesh_relative_error']:.3e}；相对 COMSOL intorder=6 原生积分误差 {checks['cell_sum_vs_COMSOL_native_relative_error']:.6%}。
三域六分量储能相对旧点分配的变化 {checks['regional_energy_relative_error']:.3e}（积分域修正，不是全能量误差）；全 EM 储能相对 COMSOL 原生积分差 {checks['stored_energy_relative_error_to_COMSOL']:.6%}。
有符号贡献和相对总强度误差 {checks['sum_signed_contribution_relative_error']:.3e}。合成源测试独立比较显式八象限体积积分、跨轴晶胞、矢量奇偶及干涉预算。

**独立空气参考仍有差异。** 空气 E/H 提取的实际中心为 {checks['air_center_W_sr']:.12g} W/sr；
本体源模型相对它的中心复振幅误差 {checks['cell_sum_vs_air_complex_relative_error']:.6%}，强度差 {checks['source_vs_air_intensity_relative_error']:.6%}。
未通过原 2%/5% 独立闭合标准；本轮仍给出逐胞定量结果，并将其标为体源模型中的贡献，不能直接宣称是已验证的真实中心辐射百分比。
没有整体校准幅值、重新对齐相位或删除 pad 来强行闭合。

**已修正的数值问题及剩余精度限制。** 旧点归属规则使完全相同的 cavity 晶胞体积波动最高约 0.727%，强相消的辐射残差被这一误差放大。
当前几何裁切后，全体晶胞体积最大相对误差为 {checks['cell_volume_relative_error']:.3e}；48 个原始正象限胞的对照中，代表相邻两胞的 Ey 强度比由 12844.6 降至 1.051。
这个对照保持场、mesh、相位和镜像规则不变，只修正跨胞积分。跨轴单胞的 Ex 相消仍保留，不能将其与随机胞间散斑混为一谈。
几何裁切精确，但跨胞四面体内采用四点确定的仿射源/密度；这不是原 COMSOL 高阶场的精确再评价，也不是 mesh/每胞弱辐射的独立收敛证明。
总积分守恒是此重构的代数性质，不能当成独立的物理闭合。后续若需要更精确的单胞绝对值，应在裁切子单元上直接评价保存解并验证求积阶数收敛。
本轮没有新增 COMSOL 取点。诊断数据见 cell_quadrature_audit_patch/all 的 NPZ、CSV、JSON；原结果按 pointbin 命名保留。

**输出与复现。** ../01_results/{STEM}.npz 保存 4681 胞及 pad 的 Ex/Ey 复 F、六分量 U、相消指标、原始积分量和参照。
../80_logs/{STEM}_cells.csv 每行一个胞，含实际 i/j、shell、位置、区域、储能、相位、自项及有符号贡献；最后一行为 pad。
同目录 {STEM}_regions.csv 和 {STEM}_shells.csv 保存相干分组预算；{STEM}_checks.json 保存身份哈希和全部核对。
先执行 `.venv\Scripts\python.exe -B -m scripts.analysis.audit_dipolar_cell_quadrature --all-cells` 构建几何裁切缓存，
再执行 `.venv\Scripts\python.exe -B -m scripts.analysis.analyze_dipolar_cell_radiation` 完成预算与绘图；`--check` 只跑合成检查，`--plot-only` 从缓存绘图。
本轮增加 Fig.7 及数值表，保留 Fig.1–6 和所有 MPH；没有启动 COMSOL，也未新建目录。
""")
    audit_path=a.L/"cell_quadrature_diagnosis.json"
    if audit_path.exists():
        audit=json.loads(audit_path.read_text(encoding="utf8"))
        report+="\n**保持同一实际场的局部对照。** 下列所有晶胞都在原始正象限内，远离镜像轴。\n\n"
        report+="| 通道、相邻胞 | 原点分配强度比 | 裁切仿射积分强度比 |\n| --- | ---: | ---: |\n"
        for row in audit["adjacent_cell_examples"]:
            report+=f"| Ey，{row['first_cell']} / {row['second_cell']} | {row['pointbin_Ey_intensity_ratio']:.6g} | {row['clipped_affine_Ey_intensity_ratio']:.6g} |\n"
        row=audit["Ex_adjacent_example"]
        report+=f"| Ex，{row['first_cell']} / {row['second_cell']} | {row['pointbin_intensity_ratio']:.6g} | {row['clipped_affine_intensity_ratio']:.6g} |\n"
        report+=f"\n第一对晶胞的完整 E/H 能量内积归一重叠为 {audit['full_EM_shape_overlap_first_pair']:.9%}；胞内场形状相近，原积分却产生四个数量级的强度比。\n"
        report+="旧 cavity 的 Ey 胞内相消指标中位数约 0.606%，与旧分胞体积误差的量级相近，因此边界权重误差足以主导微弱辐射残差。这个尺度比较不是严格误差上界，直接裁切对照才是定位依据。\n"
        report+="\n**裁切积分的具体算法。** 以原四面体顶点的重心坐标 lambda 表示位置，四个 Gauss 节点的重心坐标组成矩阵 B。\n"
        report+="对每个被积量 s（复极化源分量、六分量能量密度及源模长），用 s_vertex=B^(-1)s_Gauss 定义单元内仿射重构。\n"
        report+="原四面体与单胞六条竖直边界半空间逐面求交；将所得凸多面体各面与内部点组成四面体，准确求得 m_j=integral_cell_intersection lambda dV。\n"
        report+="胞内积分为 m_j^T B^(-1)s_Gauss；对完全在单胞内的原单元仍用原四点规则。完整电磁储能使用同一裁切方法，不只修正电场源。\n"
        report+="仿射多项式的积分与分区加和是精确的；原有限元高阶场与其仿射近似之间的差异仍需独立取点/阶数检查。没有跨胞拟合、空间平滑或强制连续化。\n"
    (a.R / f"{STEM}_report.md").write_text(report,encoding="utf8")
    path = a.R / "report.md"
    text = path.read_text(encoding="utf8").split("\n" + HEADING)[0]
    sector_report=a.R / "sector_interference_report.md"
    appendix="\n" + sector_report.read_text(encoding="utf8") if sector_report.exists() else ""
    path.write_text(text.rstrip() + "\n" + report + appendix,encoding="utf8")
    checks["figure"] = dict(name=FIGURE,dpi=a.FIGURE_DPI,panels=6,
        quantity="hexagonal cell integrals, not local field samples",xy_limit_um=bound,
        U_color_limits=[top_u*1e-5,top_u],self_color_limits=[top_self*1e-5,top_self],
        signed_color_limits=[-top_signed,top_signed],phase_color_limits=[-np.pi,np.pi])
    a.write_json(f"{STEM}_checks.json",checks)
    print("Fig.7 PNG/PDF and detailed per-cell report saved.",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--plot-only",action="store_true")
    args = parser.parse_args()
    if args.check:
        self_check()
    else:
        if not args.plot_only:
            self_check()
            compute()
        plot_report()
