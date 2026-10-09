"""Offline mode19 analysis. No COMSOL imports or solver calls.

Run from repository root: python -B -m scripts.analysis.analyze_dipolar_singularity
Only writes classified folders in the designated dipolar-analysis workspace.
"""
from __future__ import annotations

from comsol_workflow.output_paths import table_path
from comsol_workflow.dipolar_radiation import forward, inverse, backprop, energy_group_gram, budget_rows

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from matplotlib.path import Path as PolygonPath
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
from scipy.constants import epsilon_0, mu_0, c
from scipy.fft import fft2, ifft2, fftshift, ifftshift
from scipy.spatial import ConvexHull

from comsol_workflow.farfield_fft import load_air_field, symmetric_pad
from comsol_workflow.hex_lattice_utils import unit_cell_corners
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
from comsol_workflow.figure_output import save_figure_formats

ROOT = Path(__file__).resolve().parents[2]
INPUTS = load_dipolar_inputs()
SERIES = INPUTS.series_dir
MODE = INPUTS.mode_dir
WORK = INPUTS.output_dir
REPORTS = WORK / "12_reports"
REF = INPUTS.reference_dir
CONFIG = INPUTS.source_config_dir
Z0 = np.sqrt(mu_0 / epsilon_0)
SOURCES = {}
METRICS = {}


def digest(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def register(path):
    path = Path(path).resolve()
    SOURCES[str(path)] = {"sha256": digest(path), "bytes": path.stat().st_size}
    return path


def read_json(path):
    return json.loads(register(path).read_text(encoding="utf-8"))


def write_json(name, data):
    (WORK / "80_logs" / name).write_text(json.dumps(data, indent=2, ensure_ascii=False,
                                               allow_nan=False) + "\n", encoding="utf-8")


def log(message):
    print(message, flush=True)
    with (WORK / "80_logs/analysis.log").open("a", encoding="utf-8") as f:
        f.write(message + "\n")












def save(fig, name, dpi=190):
    fig.get_layout_engine().set(w_pad=.20,h_pad=.07,wspace=.12,hspace=.06)
    files = [path.relative_to(WORK).as_posix() for path in
             save_figure_formats(fig, WORK, name, dpi=dpi, bbox_inches="tight")]
    manifest_path=WORK/"80_logs/figure_manifest.json"
    metadata=json.loads(manifest_path.read_text(encoding="utf8")) if manifest_path.exists() else {}
    metadata[name]={"dpi":dpi,"files":files,
        "panels":[{"title":ax.get_title(),"xlabel":ax.get_xlabel(),"ylabel":ax.get_ylabel(),
        "xlim":list(ax.get_xlim()),"ylim":list(ax.get_ylim()),
        "data_colors":[{"cmap":artist.get_cmap().name,"clim":list(artist.get_clim()),
                        "scale":type(artist.norm).__name__}
                       for artist in [*ax.images,*ax.collections] if artist.get_array() is not None]}
        for ax in fig.axes]}
    write_json("figure_manifest.json",metadata)
    plt.close(fig)


def colorbar(fig, ax, artist, label):
    cax = make_axes_locatable(ax).append_axes("right", size="4.5%", pad=.08)
    bar = fig.colorbar(artist, cax=cax)
    bar.set_label(label, rotation=270, labelpad=18)
    bar.ax.tick_params(direction="in",labelsize=8)
    bar.ax.yaxis.label.set_size(9)
    return bar


def maps(name, arrays, titles, axis, label, suptitle, logscale=False, limit=None, overlay=None):
    count = len(arrays)
    columns = 3 if count in (3, 6) else 2
    rows = int(np.ceil(count/columns))
    fig, axes = plt.subplots(rows, columns, figsize=(5*columns, 4.6*rows), squeeze=False,
                             layout="constrained")
    maximum = max(float(np.max(a)) for a in arrays)
    norm = LogNorm(maximum*1e-6, maximum) if logscale else plt.Normalize(0, maximum)
    extent = [axis[0], axis[-1]]*2
    for ax, data, title in zip(axes.flat, arrays, titles):
        image = ax.imshow(data, origin="lower", extent=extent, cmap="magma", norm=norm)
        ax.set(title=title, aspect="equal")
        is_k = "farfield" in name or "cached" in name
        ax.set_xlabel(r"$q_x/(2\pi/a)$" if is_k else r"$x$ ($\mu$m)")
        ax.set_ylabel(r"$q_y/(2\pi/a)$" if is_k else r"$y$ ($\mu$m)")
        if is_k:
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
        if limit:
            ax.set_xlim(-limit, limit); ax.set_ylim(-limit, limit)
        if overlay is not None:
            ax.plot(*overlay.T, color="white", lw=.65, alpha=.65)
        colorbar(fig, ax, image, label)
    for ax in list(axes.flat)[count:]:
        ax.set_visible(False)
    fig.suptitle(suptitle, fontsize=13)
    save(fig, name)


def analyze_air(cfg, eig, geometry):
    air = load_air_field(register(MODE / "11_simulation_exports/E_air.parquet"))
    N = int(cfg["fft_size"]); middle = N//2
    dx = float(np.diff(air["x_axis"]).mean())*1e-6
    np.testing.assert_allclose(np.diff(air["x_axis"])*1e-6, dx, rtol=1e-10)
    np.testing.assert_allclose(np.diff(air["y_axis"])*1e-6, dx, rtol=1e-10)
    frequency = float(eig["re"])*1e12; k0 = 2*np.pi*frequency/c
    # Saved geometry explicitly records thickness, and its upper-half export is centered at z=0.
    height_text = geometry["simulation_common"]["slab_height"].replace(" ", "")
    if not height_text.endswith("[nm]"):
        raise ValueError("Expected explicit saved slab_height in nm")
    z_surface = float(height_text[:-4])*1e-9/2
    air_meta = read_json(CONFIG / "air_field_metadata.json")
    z_air = float(air_meta["air_plane_z_um"])*1e-6
    d = z_air-z_surface
    assert d > 0 and air_meta["air_plane_expression"] == "z_slab_top + H_air"
    METRICS.update(z_surface_um=z_surface*1e6, z_air_um=z_air*1e6, backprop_distance_um=d*1e6)
    q = 2*np.pi*fftshift(np.fft.fftfreq(N, dx))
    qx, qy = np.meshgrid(q, q)
    qr2 = qx*qx+qy*qy; mask = qr2 < k0*k0
    kz = np.sqrt(np.maximum(k0*k0-qr2, 0))
    n = np.stack((qx/k0, qy/k0, kz/k0), axis=-1)
    spectrum = np.stack([forward(air[a], dx, N) for a in ("Ex", "Ey", "Ez")], axis=-1)
    spectrum *= mask[..., None]
    center = spectrum[middle, middle].copy()
    direct = np.array([air[a].sum()*dx*dx for a in ("Ex", "Ey", "Ez")])
    METRICS["fft_center_vs_direct_relative_error"] = float(np.linalg.norm(center-direct)/np.linalg.norm(direct))
    assert METRICS["fft_center_vs_direct_relative_error"] < 1e-10
    longitudinal = np.sum(n*spectrum, axis=-1)
    projected = spectrum - n*longitudinal[..., None]
    METRICS["air_longitudinal_squared_norm_fraction"] = float(np.sum(abs(longitudinal)**2)/np.sum(abs(spectrum)**2))
    projected *= mask[..., None]
    phase = np.exp(1j*kz*d)
    near = np.stack([inverse(spectrum[..., j]*phase, dx) for j in range(3)], -1)
    # P1 raw component transform remains intact; this separate transverse field supports Poynting checks.
    et = np.stack([inverse(projected[..., j]*phase, dx) for j in range(3)], -1)
    h_spectrum = np.cross(n, projected)/Z0
    ht = np.stack([inverse(h_spectrum[..., j]*phase, dx) for j in range(3)], -1)
    full_axis_um = (np.arange(N)-middle)*dx*1e6
    ex_energy = epsilon_0/4*abs(near)**2
    ue_t = epsilon_0/4*np.sum(abs(et)**2, axis=-1)
    uh = mu_0/4*np.sum(abs(ht)**2, axis=-1)
    sz = .5*np.cross(et, ht.conj())[..., 2].real
    dq = float(q[1]-q[0])
    spectral_power = float(np.sum((kz/k0)*np.sum(abs(projected)**2, axis=-1))*dq*dq/(2*Z0*(2*np.pi)**2))
    plane_power = float(sz.sum()*dx*dx)
    METRICS["outgoing_air_power_W_native_eigenmode_scale"] = spectral_power
    METRICS["poynting_parseval_relative_error"] = abs(plane_power-spectral_power)/spectral_power
    assert METRICS["poynting_parseval_relative_error"] < 1e-10
    k_integral = np.sum(abs(spectrum)**2)*dq*dq/(2*np.pi)**2
    METRICS["electric_parseval_relative_error"] = float(abs(np.sum(abs(near)**2)*dx*dx-k_integral)/k_integral)
    back_center = near.sum(axis=(0, 1))*dx*dx
    expected_center = center*np.exp(1j*k0*d)
    METRICS["surface_integral_center_relative_error"] = float(np.linalg.norm(back_center-expected_center)/np.linalg.norm(expected_center))
    assert METRICS["surface_integral_center_relative_error"] < 1e-10
    mask09 = qr2 < (.9*k0)**2
    METRICS["NA09_fraction_of_NA1_electric_spectral_norm"] = float(np.sum(abs(spectrum[mask09])**2)/np.sum(abs(spectrum)**2))
    METRICS["NA09_fraction_of_NA1_radiative_power"] = float(np.sum((kz/k0)[mask09]*np.sum(abs(projected[mask09])**2, axis=-1))*dq*dq/(2*Z0*(2*np.pi)**2)/spectral_power)
    surface09 = np.stack([backprop(spectrum[..., j], kz, d, mask09, dx) for j in range(2)], -1)
    old09 = np.stack([backprop(spectrum[..., j], kz, z_air, mask09, dx) for j in range(2)], -1)
    METRICS["surface_E_xy_change_NA09_to_NA1_L2"] = float(np.linalg.norm(near[..., :2]-surface09)/np.linalg.norm(near[..., :2]))
    METRICS["surface_vs_midplane_NA09_E_xy_L2"] = float(np.linalg.norm(surface09-old09)/np.linalg.norm(surface09))
    # Spatial split refers to the filtered field, not independently localized volume radiation sources.
    centers = np.array([[p["x"], p["y"]] for p in geometry["bulk_points"]])
    corners = unit_cell_corners(float(geometry["a"]))
    hull = ConvexHull((centers[:, None, :]+corners[None, :, :]).reshape(-1, 2))
    xx, yy = np.meshgrid(full_axis_um, full_axis_um)
    distance = np.full(xx.shape, -np.inf)
    for nx, ny, b in hull.equations:
        distance = np.maximum(distance, nx*xx+ny*yy+b)
    a = float(geometry["a"])
    regions = np.where(distance < -a, 0, np.where(distance <= a, 1, 2))
    labels = ["cavity_interior", "interface_envelope_band", "exterior_and_filter_tails"]
    integrals = np.array([near[regions==i, :2].sum(axis=0)*dx*dx for i in range(3)])
    np.testing.assert_allclose(integrals.sum(axis=0), expected_center[:2], atol=np.linalg.norm(center)*1e-10)
    factor = (k0/(2*np.pi))**2/(2*Z0)
    pd.DataFrame(budget_rows(labels, integrals, factor, "filtered-surface spatial regions; W/sr")).to_csv(table_path(WORK / '80_logs', 'surface_center_budget.csv'), index=False)
    fractions = []
    for i, name in enumerate(labels):
        fractions.append(dict(region=name, Ex_electric_fraction=float(ex_energy[..., 0][regions==i].sum()/ex_energy[..., 0].sum()),
                              Ey_electric_fraction=float(ex_energy[..., 1][regions==i].sum()/ex_energy[..., 1].sum()),
                              Ey_center_amplitude_re=float((integrals[i,1]/expected_center[1]).real),
                              Ey_center_amplitude_im=float((integrals[i,1]/expected_center[1]).imag)))
    pd.DataFrame(fractions).to_csv(table_path(WORK / '80_logs', 'surface_region_energy.csv'),index=False)
    METRICS["surface_region_fractions"] = fractions
    f0 = 1j*k0/(2*np.pi)*np.exp(1j*k0*z_air)*center
    METRICS["center_Ex_over_Ey_power"] = float(abs(f0[0])**2/abs(f0[1])**2)
    METRICS["center_dP_dOmega_W_sr"] = (abs(f0[:2])**2/(2*Z0)).tolist()
    # Exact mirror statistics reflect the quarter-simulation reconstruction, not a separate symmetry solve.
    METRICS["surface_mirror_checks"] = {
        "Ex_odd_x_relative_L2":float(np.linalg.norm(near[...,0]+near[:,::-1,0])/np.linalg.norm(near[...,0])),
        "Ex_odd_y_relative_L2":float(np.linalg.norm(near[...,0]+near[::-1,:,0])/np.linalg.norm(near[...,0])),
        "Ey_even_x_relative_L2":float(np.linalg.norm(near[...,1]-near[:,::-1,1])/np.linalg.norm(near[...,1])),
        "Ey_even_y_relative_L2":float(np.linalg.norm(near[...,1]-near[::-1,:,1])/np.linalg.norm(near[...,1]))}
    crop = np.flatnonzero(abs(full_axis_um) <= max(abs(air["x_axis"]))+dx*1e6*.1)
    sl = slice(crop[0], crop[-1]+1)
    near_axis = full_axis_um[sl]
    np.savez_compressed(WORK/"01_results/radiative_fields.npz", q_axis_m_inv=q, lightcone_flat_indices=np.flatnonzero(mask),
                        E_angular_lc_V_m=spectrum[mask], E_transverse_angular_lc_V_m=projected[mask],
                        x_full_um=full_axis_um, x_roi_um=near_axis, E_surface_roi_V_m=near[sl,sl],
                        E_transverse_surface_roi_V_m=et[sl,sl], H_outgoing_surface_roi_A_m=ht[sl,sl],
                        region_center_integrals_V_m=integrals, region_names=labels,
                        z_surface_m=z_surface, z_air_m=z_air, dx_m=dx, frequency_Hz=frequency)
    boundary = hull.points[hull.vertices]
    boundary = np.vstack([boundary,boundary[0]])
    maps("Fig.3_radiative_nearfield", [ex_energy[sl,sl,j] for j in (0,1)]+[ex_energy[sl,sl].sum(axis=-1)],
         ["Ex: radiative electric energy", "Ey: radiative electric energy", "Total electric energy"],
         near_axis, r"$u_e$ (J m$^{-3}$)", "Surface radiative field | NA=1 | z=0.10 um | common native eigenmode scale",overlay=boundary)
    maps("Fig.4_outgoing_EM_nearfield",[ue_t[sl,sl],uh[sl,sl],(ue_t+uh)[sl,sl]],
         ["Electric: transverse outgoing field","Magnetic: reconstructed outgoing field","Complete outgoing EM energy"],
         near_axis,r"$u$ (J m$^{-3}$)","Air-side outgoing-wave reconstruction | independent slab H is unavailable",overlay=boundary)
    scale = 2*np.pi/(a*1e-6)
    far_axis = q/scale
    sf = np.flatnonzero(abs(far_axis) <= .10); sf=slice(sf[0],sf[-1]+1)
    F = projected[ sf,sf ] * (1j*kz[sf,sf]/(2*np.pi)*np.exp(1j*kz[sf,sf]*z_air))[...,None]
    far_power = abs(F)**2/(2*Z0)
    maps("Fig.5_farfield_channels",[far_power[...,0],far_power[...,1],far_power.sum(axis=-1)],
         ["Ex Cartesian contribution","Ey Cartesian contribution","Full vector radiation"],
         far_axis[sf],r"$dP/d\Omega$ (W sr$^{-1}$)",
         "Far field near normal emission | NA=1 | common native eigenmode scale",logscale=True)
    write_json("analysis_config.json",dict(mode_idx=19,user_mode="py",frequency_Hz=frequency,Q=float(eig["q"]),
               NA=1,fft_size=N,source_grid=air["Ex"].shape[0],z_air_um=z_air*1e6,z_surface_um=z_surface*1e6,
               backprop_distance_um=d*1e6,time_convention="exp(+i omega t)",forward_sign="+i q.r",
               nearfield="independent raw Ex/Ey/Ez filtering; transverse projection is separately saved for outgoing EM energy",
               normalization="native common eigenmode scale; full stored U_ref unavailable",
               workspace_layout={"01_results":"raw arrays","80_logs":"CSV/JSON/logs","10_overview":"PNG previews","11_pdf":"PDF figures","12_reports":"Markdown reports/specs"}))
    del spectrum, projected, near, et, ht, surface09, old09
    log("Air transform, surface backpropagation, energy/flux checks and figures complete.")
    return f0



def periodic_z0():
    """Recompute z=0 images from native fields, not prior derived results."""
    manifest=read_json(REF/"99_config/manifest.json")
    geometry_path=register(Path(manifest["geometry_source"]))
    assert digest(geometry_path)==manifest["geometry_sha256"]
    geometry=json.loads(geometry_path.read_text(encoding="utf8"))
    register(ROOT/"scripts/run_main/run_dipolar_periodic_reference.py")  # provenance only, never imported
    with np.load(register(REF/"01_results/k00_reference.npz")) as d:
        assert np.allclose(d["k_um_inv"],0)
        points=d["rho_points_um"]
        mode=int(d["gamma_band_to_mode"][1])
        e=d["E_xy_center_native"][mode]
    inside_holes=np.zeros(len(points),bool)
    for hole in geometry["holes"]:
        inside_holes |= PolygonPath(np.asarray(hole)).contains_points(points)
    cfg=read_json(CONFIG/"config.json")
    eps=np.where(inside_holes,1.,float(cfg["simulation_common"]["refractive_index"])**2)
    ue=epsilon_0/4*eps[:,None]*abs(e)**2
    xy=np.round(points,12)
    x=np.unique(xy[:,0]);y=np.unique(xy[:,1])
    for axis in (x,y):
        np.testing.assert_allclose(np.diff(axis),np.diff(axis)[0],rtol=1e-8)
    image=np.full((len(y),len(x),2),np.nan)
    image[np.searchsorted(y,xy[:,1]),np.searchsorted(x,xy[:,0])]=ue
    maximum=float(np.nanmax(image.sum(-1)))
    np.savez_compressed(WORK/"01_results/periodic_z0_fields.npz",x_um=x,y_um=y,z_nm=0,
                        E_xy_native=e,rho_points_um=points,epsilon_r=eps,
                        electric_energy_density_grid=image,display_reference_max=maximum)
    from pyarrow.parquet import read_schema
    columns=read_schema(register(MODE/"11_simulation_exports/Hz_center.parquet")).names
    METRICS["slab_z0_available_columns"]=columns
    METRICS["finite_slab_z0_E_status"]="missing; no interpolation from nonzero-z quadrature"
    METRICS["periodic_z0_status"]="native Ex/Ey, material-weighted image; independent reference amplitude"
    return x,y,image,maximum


def report():
    m=METRICS
    content=f"""# mode19：当前可视化框架与数据缺口

本版按最新要求从已保存的原始场重新计算；没有启动 COMSOL，没有读取旧 S5 辐射重建、
旧投影系数或旧误差预算作为本轮结果。输出按 01_results/ 数组、80_logs/ 日志、10_overview/ PNG、11_pdf/ PDF、12_reports/ Markdown 分类；00_model/ 保留模型。
图片使用 Fig.1_ 等命名；数字文件夹前缀不再用于图片。

## 当前图组

| 图 | 内容与数据状态 |
| --- | --- |
| [Fig.1 总览](../10_overview/Fig.1_analysis_framework.png) | 周期场 → 有限薄板 → 可辐射表面场 → 远场；Ex/Ey/完整电磁量三行 |
| [Fig.2 薄板截面](../10_overview/Fig.2_finite_slab_z0.png) | 有限腔 z=0 nm；缺少原始 Ex/Ey 和完整 E/H，明确 NO DATA 占位 |
| [Fig.3 可辐射近场](../10_overview/Fig.3_radiative_nearfield.png) | NA=1 回溯至上表面空气侧，Ex/Ey/总电能密度 |
| [Fig.4 出射电磁场](../10_overview/Fig.4_outgoing_EM_nearfield.png) | 横向投影后按向上平面波关系重建 H，展示电、磁和总能量密度 |
| [Fig.5 远场](../10_overview/Fig.5_farfield_channels.png) | 从原始 E_air 重新求角谱，展示中心邻域 Ex/Ey/完整矢量 dP/dΩ |
| [Fig.6 贡献链占位](../10_overview/Fig.6_component_chain.png) | 六类来源逐阶段的完整 EM、Ex/Ey 表面场及中心预算；数据缺失时保留空面板 |

周期 Γ-py 面板已改为原始 z=0 nm Ex/Ey 的二维场图；材料加权量为 ε|Eα|²/4，
采用 Ex/Ey 共用的参考峰值归一化。该周期参考具有独立振幅，不能与有限模式直接按颜色比较能量。
截面采样映射到原始规则网格，采样范围外保持空白，不用散点大小表示能量。

有限薄板中心 parquet 字段为 {m['slab_z0_available_columns']}，只包含 Hz。
上一版点图是**每个晶胞体积分的电能**，不是 z=0 的实际场分布，已经撤下。
体求积节点不含 z=0，不把它们平移、外推或插值成该平面的实际电场。
所有 finite slab 面板现在指定 z=0 nm；缺数面板没有数值场，也不以 Hz 替代 Ex/Ey。

内部截面取 z=0 nm；NA=1 可辐射近场仍回溯到真实上表面 z={m['z_surface_um']*1000:.6g} nm 的空气侧。
不能把空气传播算子延伸到介质内部，并将其结果称为实际 slab 场。

## 上一版 03、04、05、06 的含义和处理

- **原 03_farfield_channels**：Ex、Ey 的 Cartesian 辐射贡献及全矢量 dP/dΩ。
  横纵轴是空气出射面内波矢 q，颜色是每单位立体角的辐射功率；不是内部 Bloch k 权重。
  它直接关联中心亮点，保留为本版 Fig.5，并从原始空气复场重新计算。
- **原 04_cached_chain_gap**：旧 19-k 四带 cavity 重建、旧直接远场和两者复场差的平方幅度对比。
  这是旧字典的误差诊断，不能作为新的机制解释；本版移除，未使用其缓存或数值结论。
- **原 04_component_chain**：六类来源的贡献链框架，不是已完成的因果分解。
  本版 Fig.6 仅保留待填的数据位置，删除旧投影数值旁注。
- **原 05_center_budget**：旧中心辐射的自身项、干涉、体源和空气参考的柱状比较；本版移除。
- **原 06_closure_symmetry**：变换一致性、镜像对称性和旧独立闭合误差的柱状比较；本版移除。
  当前重新计算的实现检查保留在 80_logs/validation.json，不单独绘柱状图。

## 本轮计算与检查

时间约定 exp(+iωt)，空间 Fourier 分析使用正号，回溯乘 exp(+ikz d)。
空气面 z={m['z_air_um']:.6g} μm，回溯距离 {m['backprop_distance_um']:.6g} μm。
Ex/Ey/Ez 均以复场分别变换，保留相位。Ex/Ey 在法向是两个正交通道，离轴总量还须包含 Ez。

- FFT 中心与直接空间复积分相对差：{m['fft_center_vs_direct_relative_error']:.3e}。
- 表面复积分检查：{m['surface_integral_center_relative_error']:.3e}。
- 电场 Parseval 检查：{m['electric_parseval_relative_error']:.3e}。
- 出射 H 重建后的角谱/穿面功率检查：{m['poynting_parseval_relative_error']:.3e}。
- 原始空气角谱纵向平方范数占比：{m['air_longitudinal_squared_norm_fraction']:.6%}。
- 中心 Ex/Ey 功率比：{m['center_Ex_over_Ey_power']:.3e}，Ex 处于数值零量级。

Fig.3 保留原始分量滤波场；Fig.4 使用独立保存的横向投影结果和推导出的向上出射 H。
推导的空气 H 不是薄板磁场，不能填充完整薄板储能缺口。
本征模振幅任意，有限场图中的 J、W 采用同一原生数值尺度；没有伪造完整 U_ref 归一化。
一致性检查不等于独立 Maxwell 因果闭合，当前只确认中心由 Ey 主导，尚未完成来源归因。

## 需要补充的数据或计算（本轮不执行）

1. 从已保存 mode19 解补导出 z=0 nm 的完整 E/H 与坐标、材料 ε/μ、域掩膜，
   覆盖 cavity、cladding 和外缘；保存复振幅和统一相位。取得后绘制真实二维场分布。
2. 完整体储能及周期分解需要相同体求积点上的有限/周期 E/H、材料和权重。
   二维截面图本身不能给出整个薄板的完整储能。
3. 六类来源在同一材料与评价频率下重新计算完整 NA=1 复辐射角谱，
   再逐项回溯近场，构建自身项与干涉预算；不用旧 19-k 辐射缓存填图。
4. 需要独立空气 E/H 或额外高度/孔径导出，检查向上/向下波与有限窗口误差。
5. 先补导出已有解；只有保存解或精度不足时，再列出新求解范围。当前未启动 COMSOL。
"""
    report_path=REPORTS/"report.md"
    if report_path.exists():
        previous=report_path.read_text(encoding="utf8")
        begin="<!-- BEGIN FIGURE METHODS -->"
        end="<!-- END FIGURE METHODS -->"
        if begin in previous:
            if previous.count(begin)!=1 or previous.count(end)!=1 or previous.index(end)<previous.index(begin):
                raise ValueError("Malformed detailed figure-methods block; report preserved")
            details=previous[previous.index(begin):previous.index(end)+len(end)]
            content=content.replace("## 上一版",details+"\n\n## 上一版",1)
    report_path.write_text(content,encoding="utf8")
    write_json("validation.json",{"offline_execution":"complete","COMSOL_started":False,
        "legacy_derived_results_reused":False,"mechanism":"not_closed",
        "finite_slab_z0":"missing raw E/H; placeholders","metrics":m})



def self_check():
    n=33;dx=.2; axis=(np.arange(n)-n//2)*dx
    xx,yy=np.meshgrid(axis,axis);k=2*np.pi/(n*dx)
    field=np.exp(-1j*(3*k*xx+2*k*yy))
    a=forward(field,dx,n)
    assert np.unravel_index(abs(a).argmax(),a.shape)==(n//2+2,n//2+3)
    np.testing.assert_allclose(inverse(a,dx),field,atol=2e-14)
    kz=np.full((n,n),4.); mask=np.ones((n,n),bool)
    evolved=a*np.exp(-1j*kz*.7)
    np.testing.assert_allclose(backprop(evolved,kz,.7,mask,dx),field,atol=2e-14)
    rng=np.random.default_rng(4);b=rng.normal(size=(20,3))+1j*rng.normal(size=(20,3))
    f=rng.normal(size=20)+1j*rng.normal(size=20);coef=np.diag(np.array([.3+1j,.2,-.7j]))
    g=energy_group_gram(b.conj().T@b,b.conj().T@f,float(np.vdot(f,f).real),coef)
    real_fields=np.column_stack([b@x for x in coef]+[f-b@coef.sum(axis=0)])
    np.testing.assert_allclose(g,real_fields.conj().T@real_fields,atol=1e-12)
    rows=budget_rows(["one","two"],np.array([[1,1],[1,-1]],complex),1,"test")
    totals=[r["value"] for r in rows if r["kind"]=="total"]
    assert totals==[4.,0.]
    print("PASS: physical Fourier sign, inverse scale, backprop phase, nonorthogonal energy/residual, coherent/cancelling budgets")


def main():
    for name in ("01_results","80_logs","10_overview","11_pdf","12_reports"):
        (WORK/name).mkdir(parents=True,exist_ok=True)
    assert not [p for p in WORK.iterdir() if p.is_file()], "Workspace root must contain directories only"
    (WORK/"80_logs/analysis.log").write_text("",encoding="utf-8")
    plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    self_check()
    write_json("figure_manifest.json",{})
    log("Offline execution started. No solver imports/calls. Existing sources read-only.")
    cfg=read_json(MODE/"12_farfield_FFT/farfield_summary.json")
    eig=pd.read_csv(register(MODE/"10_overview/eigenfrequency.csv")).iloc[0]
    assert int(eig["mode_idx"])==19
    geometry=read_json(CONFIG/"config.json")
    register(ROOT/"docs/spec/plan-execute_20260921_dipolar_singularity_analysis.md")
    register(Path(__file__))
    analyze_air(cfg,eig,geometry)
    framework()
    log("Checking source fingerprints and classified output layout.")
    for path,data in SOURCES.items():
        assert digest(path)==data["sha256"],path
    write_json("source_manifest.json",{"sources":SOURCES,"COMSOL_started":False,
               "input_policy":"read only; arrays under 01_results, tables/logs under 80_logs, PNG under 10_overview, PDF under 11_pdf, Markdown under 12_reports"})
    report()
    log("Offline deliverables finished. Complete slab EM and causal radiation closure remain explicitly pending.")
    print(json.dumps(METRICS,indent=2),flush=True)



def placeholder(ax, reason, title=None):
    """Missing data are blank, never a numerical zero or simulated field."""
    ax.set_facecolor("#f5f5f5")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#d7d7d7")
    if title:
        ax.set_title(title, fontsize=10)
    ax.text(.5,.62,"NO DATA",ha="center",va="center",transform=ax.transAxes,
            fontsize=14,color="#777777",weight="bold")
    ax.text(.5,.35,reason,ha="center",va="center",transform=ax.transAxes,
            fontsize=9,color="#666666",linespacing=1.6)


def framework():
    """Assemble available fields and explicit gaps; no solver or new field fit."""
    with np.load(WORK/"01_results/radiative_fields.npz") as data:
        near=data["E_surface_roi_V_m"]
        transverse=data["E_transverse_surface_roi_V_m"]
        magnetic=data["H_outgoing_surface_roi_A_m"]
        axis=data["x_roi_um"]
        q=data["q_axis_m_inv"]
        indices=data["lightcone_flat_indices"]
        angular=data["E_transverse_angular_lc_V_m"]
        frequency=float(data["frequency_Hz"])
        z_surface=float(data["z_surface_m"])
    ref_x,ref_y,ref_energy,ref_max=periodic_z0()
    raw_u=epsilon_0/4*abs(near)**2
    u_em=epsilon_0/4*np.sum(abs(transverse)**2,-1)+mu_0/4*np.sum(abs(magnetic)**2,-1)
    k0=2*np.pi*frequency/c
    iy,ix=np.unravel_index(indices,(len(q),len(q)))
    kz2=np.maximum(k0*k0-q[ix]**2-q[iy]**2,0)
    lc_power=abs(angular)**2*kz2[:,None]/(8*np.pi**2*Z0)
    select=np.flatnonzero(abs(q/k0)<1.035);lo,hi=select[0],select[-1]+1
    far=np.full((hi-lo,hi-lo,3),np.nan)
    far[iy-lo,ix-lo]=lc_power
    fig,axes=plt.subplots(3,4,figsize=(22,13),layout="constrained")
    fig.get_layout_engine().set(w_pad=.16,h_pad=.07,wspace=.10,hspace=.06)
    nearmax=float(max(raw_u[...,:2].max(),u_em.max()))
    pmax=float(lc_power.sum(axis=1).max())
    for row,pol in enumerate(("Ex","Ey")):
        ax=axes[row,0]
        im=ax.imshow(ref_energy[...,row]/ref_max,origin="lower",
                     extent=[ref_x[0],ref_x[-1],ref_y[0],ref_y[-1]],cmap="magma",vmin=0,vmax=1)
        ax.set(title=f"Gamma-py: {pol}, z=0 nm",aspect="equal",
               xlabel=r"$x$ ($\mu$m)",ylabel=r"$y$ ($\mu$m)")
        colorbar(fig,ax,im,r"$u_{e,\alpha}/\max(u_{e,x}+u_{e,y})$")
        placeholder(axes[row,1],f"Raw {pol}(x,y,z=0) unavailable\nNo cell-integral or off-plane substitute",
                    f"Finite slab: {pol}, z=0 nm")
    placeholder(axes[2,0],"Periodic Hx / Hy / Hz unavailable\nFull EM normalization pending","Periodic reference: full EM")
    placeholder(axes[2,1],"Raw full E/H at z=0 unavailable\nFull energy-density map pending","Finite slab: full EM, z=0 nm")
    for row in range(3):
        ax=axes[row,2]
        values=raw_u[...,row] if row<2 else u_em
        im=ax.imshow(values,origin="lower",extent=[axis[0],axis[-1]]*2,
                     cmap="magma",vmin=0,vmax=nearmax)
        title=f"{('Ex','Ey')[row]} electric energy" if row<2 else "Outgoing EM energy (inferred H)"
        ax.set(title=title,aspect="equal",xlabel=r"$x$ ($\mu$m)",ylabel=r"$y$ ($\mu$m)")
        colorbar(fig,ax,im,r"Energy density (J m$^{-3}$)")
        ax=axes[row,3]
        values=far[...,row] if row<2 else far.sum(-1)
        im=ax.imshow(values,origin="lower",extent=[q[lo]/k0,q[hi-1]/k0]*2,
                     cmap="magma",norm=LogNorm(pmax*1e-6,pmax))
        ax.add_patch(Circle((0,0),1,fill=False,color=".55",lw=.7))
        ax.plot(0,0,"+",color="#b9d7de",ms=4,mew=.6)
        title=f"{('Ex','Ey')[row]} Cartesian contribution" if row<2 else "Full-vector radiation"
        ax.set(title=title,aspect="equal",xlabel=r"$q_x/k_0$",ylabel=r"$q_y/k_0$")
        colorbar(fig,ax,im,r"$dP/d\Omega$ (W sr$^{-1}$)")
    fig.suptitle(
        "Mode19 / py | Periodic reference  ->  finite slab  ->  radiative surface  ->  far field\n"
        f"Surface z={z_surface*1e6:.2f} um; NA=1 | Blank panels = missing data | Arrows show analysis stages, not proven attribution",
        fontsize=14)
    fig.supxlabel(
        "Slab sections: z=0 nm. Radiative surface: z=100 nm (air side). Periodic panels: independent reference shape scale.\n"
        "Air EM panels assume upward transverse waves. Ex/Ey are exhaustive orthogonal channels only at normal emission.",
        fontsize=10)
    save(fig,"Fig.1_analysis_framework")

    fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout="constrained")
    for ax,component in zip(axes,("Ex","Ey","Full EM")):
        placeholder(ax,"Raw z=0 E/H export required\nOnly Hz is available at this plane",
                    f"Finite slab: {component}, z=0 nm")
    fig.suptitle("Actual finite-slab cross section | no off-plane interpolation or per-cell integral replacement")
    save(fig,"Fig.2_finite_slab_z0")

    # The traceable six-source chain keeps its slots even when H / NA=1 responses are missing.
    names=["Gamma-py","Around-Gamma py","Other cavity bands",
           "Cavity reconstruction residual","Actual cladding","Outer pad / other material"]
    titles=["Internal full EM energy","Ex radiative surface (NA=1)",
            "Ey radiative surface (NA=1)","Center Ex / Ey power budget"]
    fig,axes=plt.subplots(6,4,figsize=(17,15),layout="constrained")
    for row,name in enumerate(names):
        for col,title in enumerate(titles):
            reason="Matched full E/H data unavailable" if col==0 else (
                "Matched full-NA complex response\nnot available" if col<3 else
                "Matched source-to-air chain\nnot yet closed")
            placeholder(axes[row,col],reason,title if row==0 else None)
        axes[row,0].set_ylabel(name,fontsize=10)
    fig.suptitle("Contribution chain | explicit data gaps\nFull EM and matched NA=1 responses are pending; no legacy reconstruction reused",fontsize=14)
    fig.supxlabel("Self terms and interference must be tracked together. A remainder defined from the direct field is not independent source attribution.",fontsize=10)
    save(fig,"Fig.6_component_chain")

    write_json("framework_status.json",{
        "entry":"10_overview/Fig.1_analysis_framework.png",
        "missing_chain":"10_overview/Fig.6_component_chain.png",
        "slab_z_nm":0,"slab_Ex_Ey_status":"missing raw z=0 fields; blank placeholders",
        "periodic_z_nm":0,"surface_z_nm":z_surface*1e9,
        "legacy_derived_results_reused":False,"COMSOL_started":False})
    print("Visualization framework and explicit missing-data panels saved.",flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--self-check",action="store_true")
    parser.add_argument("--framework-only",action="store_true",
                        help="Use saved offline arrays; preserve NO DATA slots.")
    args=parser.parse_args()
    if args.self_check:
        self_check()
    elif args.framework_only:
        METRICS.update(json.loads((WORK/"80_logs/validation.json").read_text(encoding="utf8"))["metrics"])
        plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
        framework()
        report()
    else:
        main()
