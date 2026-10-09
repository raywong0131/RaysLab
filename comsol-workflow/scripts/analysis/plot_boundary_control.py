"""Plot the saved C6v p-pair analysis; optional read-only Ey field extraction."""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import matplotlib as mpl

from comsol_workflow.boundary_plotting import _complex, _load_area_field, build_edge_figures
from scripts.analysis.decompose_boundary_c6v import field_pair
from scripts.analysis.audit_s31_boundary import native_integral
from comsol_workflow.boundary_integrals import EPS0
from scripts.run_main import run_boundary_analysis as s4


def load_control(session):
    session = Path(session).resolve()
    source = session / "results/symmetry_decomposition_v1"
    audit = s4.json_read(source / "audit.json")
    case = Path(audit["source_case"]).resolve()
    if not case.is_relative_to(session) or audit["status"] != "complete_pending_scientific_review":
        raise ValueError("Expected a completed decomposition inside this session")
    if max(audit["validation_parity"].values()) >= audit["parity_threshold"]:
        raise ValueError("Saved p-pair symmetry check did not pass")
    for name, sha in audit["source_sha256"].items():
        if s4.digest(case / name) != sha:
            raise ValueError(f"Source case changed: {name}")
    manifest = s4.json_read(case.parent / "99_config/s4_manifest.json")
    geometry = s4.json_read(case / "99_config/realized_geometry.json")
    modes = pd.read_csv(s4.data_directory(case) / "eigenfrequencies.csv")
    indices = audit["mode_indices"]
    if [r["mode_idx"] for r in audit["coefficients"]] != indices:
        raise ValueError("Saved coefficient ordering disagrees with mode ordering")
    coefficients = np.array([complex(r["coefficient_re"], r["coefficient_im"]) for r in audit["coefficients"]])
    samples = pd.read_csv(source / "area_comparison.csv").rename(columns={"samples": "area_samples"})
    fine = samples[samples.z_nm == 0].sort_values("refinement").tail(1).copy()
    outer = np.asarray(geometry["outer"]) * 1e-6
    holes = [np.asarray(h) * 1e-6 for h in geometry["holes"]]
    fine["area_m2"] = abs(s4.polygon_area(outer))
    for axis in ("x", "y"):
        for part in ("re", "im"):
            fine[f"E{axis}_integral_{part}"] = fine[f"E{axis}_area_{part}"] * fine.area_m2
    fine["normalization_factor_re"], fine["normalization_factor_im"] = 1., 0.
    frequency = modes.loc[indices, "re"]
    meta = {"mode_idx": -1, "mesh": manifest["cases"][0]["mesh"],
            "frequency_label": f"p-pair {frequency.min():.5f} / {frequency.max():.5f} THz"}
    return SimpleNamespace(session=session, source=source, audit=audit, case=case, manifest=manifest,
                           modes=modes, indices=indices, coefficients=coefficients, samples=samples,
                           fine=fine, outer=outer, holes=holes, meta=meta)


def extract_field(data):
    """Read only the two saved fields, using the previously accepted coefficients."""
    logs = data.session / "80_logs"
    path = logs / "area_field_z0.npz"
    if path.exists() or (logs / "field_extraction.json").exists():
        raise FileExistsError("Field cache already exists; use plotting without --extract-field")
    settings = data.manifest["s4"]
    refinement = int(data.fine.refinement.iloc[0])
    xy, weights, _ = s4.area_rule(data.outer, data.holes, settings["area_order"],
                                  settings["area_subdivisions"] * refinement)
    record = {"eigensolves": 0, "mesh_runs": 0, "model_saved": False,
              "source_case": str(data.case), "mode_indices": data.indices,
              "coefficients": data.audit["coefficients"], "quadrature_samples": len(xy),
              "source_audit_sha256": s4.digest(data.source / "audit.json"), "status": "starting"}
    record["preflight"] = s4.runtime_preflight()
    logs.mkdir(parents=True, exist_ok=True)
    client = model = None
    try:
        import mph
        print("Loading saved C6v model; no mesh, solve or model save", flush=True)
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics checkout failed")
        model = client.load(str(data.case / "00_model/s4_gamma.mph"))
        params = model.java.param()
        if any(abs(float(params.evaluate(k))) > 1e-12 for k in ("kx", "ky")):
            raise ValueError("Source solution is not Gamma")
        if not np.isclose(float(params.evaluate("a")), data.manifest["constants"]["A"] * 1e-6, rtol=1e-12, atol=0):
            raise ValueError("Source lattice mismatch")
        freq = np.asarray(model.java.result().numerical("gev1").getReal()).T
        if not np.allclose(freq[:, :2], data.modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Source frequencies changed")
        sampler = s4.Sampler(SimpleNamespace(model=model))

        def plane(mode, points, z, expressions):
            return np.einsum("ncm,m->nc", field_pair(sampler, data.indices, points, z, expressions), data.coefficients)

        values = plane(-1, xy, 0., ["ewfd.Ex", "ewfd.Ey"])
        for i, axis in enumerate(("x", "y")):
            expected = _complex(data.fine, f"E{axis}_integral")[0]
            if not np.isclose(weights @ values[:, i], expected, rtol=1e-9,
                              atol=max(weights @ abs(values[:, i]) * 1e-10, 1e-30)):
                raise ValueError("Read-back field disagrees with saved decomposition integral")
        s4.export_area_field(logs, SimpleNamespace(plane=plane), -1, data.outer,
                             xy, weights, values, 1., refinement)
        _load_area_field(data.session, data.fine, data.meta, refinement)
        record.update(status="complete", field_sha256=s4.digest(path))
        print(f"Field cache verified: {path}", flush=True)
    except Exception as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all(s4.digest(data.case / p) == sha for p, sha in data.audit["source_sha256"].items())
        s4.write_json(logs / "field_extraction.json", record)
    if not record["source_unchanged"]:
        raise RuntimeError("Source model/case changed")


def gradient_average(integrals, coefficients, omegas, area):
    """Convert each weighted Hz derivative using its own complex frequency."""
    integrals, coefficients, omegas = (np.asarray(v, complex) for v in (integrals, coefficients, omegas))
    if (integrals.ndim != 1 or integrals.shape != coefficients.shape or integrals.shape != omegas.shape
            or not all(np.isfinite(v).all() for v in (integrals, coefficients, omegas))
            or np.any(omegas == 0) or not np.isfinite(area) or area <= 0):
        raise ValueError("Invalid gradient conversion inputs")
    return -1j * np.sum(coefficients * integrals / omegas) / (EPS0 * area)


def extract_gradient(data):
    logs = data.session / "80_logs"
    paths = [logs / name for name in ("hz_gradient_integrals.csv", "hz_gradient_components.csv", "gradient_extraction.json")]
    if any(path.exists() for path in paths):
        raise FileExistsError("Gradient extraction already exists; reuse the saved data")
    record = {"status": "starting", "source_audit_sha256": s4.digest(data.source / "audit.json"),
              "eigensolves": 0, "mesh_runs": 0, "model_saved": False, "z_nm": 0,
              "derivative": "d(laginterp(2,ewfd.Hz),x)", "preflight": s4.runtime_preflight()}
    client = model = None
    totals, components = [], []
    try:
        import mph
        print("Reading saved model: weighted Hz x-derivative; no solve", flush=True)
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics checkout failed")
        model = client.load(str(data.case / "00_model/s4_gamma.mph"))
        if any(abs(float(model.java.param().evaluate(k))) > 1e-12 for k in ("kx", "ky")):
            raise ValueError("Saved solution is not Gamma")
        freq = np.asarray(model.java.result().numerical("gev1").getReal()).T
        if not np.allclose(freq[:, :2], data.modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Frequency identity mismatch")
        geometry = s4.json_read(data.case / "99_config/realized_geometry.json")
        eps_air, eps_slab = [float(geometry["simulation_config"][key])**2
                             for key in ("air_refractive_index", "slab_refractive_index")]
        area = abs(s4.polygon_area(data.outer))
        omegas = 2*np.pi*1e12*(data.modes.loc[data.indices, "re"].to_numpy()
                              - 1j*data.modes.loc[data.indices, "im"].to_numpy())

        def retain(method, level, values, count):
            average = gradient_average(values, data.coefficients, omegas, area)
            totals.append({"method": method, "level": level, "z_nm": 0., "samples": count,
                           "area_m2": area, "Ey_gradient": average})
            for mode, value, coefficient, omega in zip(data.indices, values, data.coefficients, omegas):
                components.append({"method": method, "level": level, "mode_idx": mode,
                                   "weighted_dxHz_integral_A": value, "coefficient": coefficient,
                                   "omega_rad_s": omega, "Ey_gradient_contribution": -1j*coefficient*value/(omega*EPS0*area)})
            print(f"{method} level {level}: <Ey>_gradient={average} V/m", flush=True)

        sampler = s4.Sampler(SimpleNamespace(model=model))
        settings = data.manifest["s4"]
        for refinement in (1, 2):
            xy, weights, air = s4.area_rule(data.outer, data.holes, settings["area_order"],
                                           settings["area_subdivisions"] * refinement)
            values = field_pair(sampler, data.indices, xy, 0., ["d(laginterp(2,ewfd.Hz),x)/(1[A/m^2])"])[:, 0]
            retain("python", refinement, (weights * np.where(air, 1/eps_air, 1/eps_slab)) @ values, len(xy))
        comp = model.java.component("comp1")
        bottom = list(comp.selection("sel_bottom").entities(2))
        domains = list(comp.selection("sel_slab_layer").entities(3))
        dielectric = list(comp.selection("sel_slab_dom").entities(3))
        updown = np.asarray(comp.geom("geom1").getUpDown())
        for order in (4, 8, 12):
            values = []
            for mode in data.indices:
                total, selected_area = 0j, 0j
                for domain in domains:
                    faces = [face for face in bottom if domain in updown[:, face-1]]
                    if len(faces) != 1:
                        raise ValueError("Ambiguous bottom material face")
                    epsilon = eps_slab if domain in dielectric else eps_air
                    value = native_integral(model, "s4_plot_dxHz", "IntSurface", faces,
                        ["1", f"side({domain},d(laginterp(2,ewfd.Hz),x))/{epsilon}"],
                        ["m^2", "A"], mode, order)
                    selected_area += value[0]
                    total += value[1]
                if not np.isclose(selected_area, area, rtol=1e-8, atol=0):
                    raise ValueError("Native gradient area mismatch")
                values.append(total)
            retain("native", order, values, None)
        s4.write_csv(paths[0], totals)
        s4.write_csv(paths[1], components)
        record.update(status="complete", output_sha256={path.name: s4.digest(path) for path in paths[:2]})
    except Exception as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all(s4.digest(data.case / p) == sha for p, sha in data.audit["source_sha256"].items())
        s4.write_json(paths[2], record)
    if not record["source_unchanged"]:
        raise RuntimeError("Source files changed")


def load_gradient(data):
    logs = data.session / "80_logs"
    record = s4.json_read(logs / "gradient_extraction.json")
    if (record["status"] != "complete" or not record["source_unchanged"]
            or record["source_audit_sha256"] != s4.digest(data.source / "audit.json")
            or not all(s4.digest(logs / name) == sha for name, sha in record["output_sha256"].items())):
        raise ValueError("Gradient cache identity mismatch")
    frame = pd.read_csv(logs / "hz_gradient_integrals.csv")
    if list(zip(frame.method, frame.level)) != [("python", 1), ("python", 2), ("native", 4), ("native", 8), ("native", 12)]:
        raise ValueError("Unexpected gradient quadrature levels")
    _complex(frame, "Ey_gradient")
    return frame


def control_figures(data):
    cache = s4.json_read(data.session / "80_logs/field_extraction.json")
    if (cache["status"] != "complete" or cache["source_audit_sha256"] != s4.digest(data.source / "audit.json")
            or cache["field_sha256"] != s4.digest(data.session / "80_logs/area_field_z0.npz")):
        raise ValueError("Field cache identity changed")
    fine = data.fine
    field = _load_area_field(data.session, fine, data.meta, int(fine.refinement.iloc[0]))
    reference = _complex(fine, "Ey_integral")[0]
    if abs(reference) <= 32 * np.finfo(float).eps * (field["weights_m2"] @ abs(field["electric_vm"][:, 1])):
        raise ValueError("Ey area integral is numerically unresolved")
    edges = pd.read_csv(data.source / "edge_integrals.csv")
    edges = edges[(edges.z_nm == 0) & (edges.refinement == fine.refinement.iloc[0])].sort_values(["hole", "edge"])
    expected = [(h, e) for h in range(1, 7) for e in range(1, 4)]
    if list(zip(edges.hole, edges.edge)) != expected:
        raise ValueError("Expected 18 unique edges at the finest z=0 section")
    geometric_edges = s4.hole_edges(data.holes)
    geometry = pd.DataFrame([{"hole": e.hole, "edge": e.edge, "group": e.group,
        "x0_m": e.start[0], "y0_m": e.start[1], "x1_m": e.end[0], "y1_m": e.end[1],
        "nx_dielectric_to_air": e.normal[0], "ny_dielectric_to_air": e.normal[1]} for e in geometric_edges])
    if not np.array_equal(edges.group, geometry.group):
        raise ValueError("Group labels disagree with geometry")
    boundary = {}
    for component, electric in (("x", "y"), ("y", "x")):
        # Ly/Lx already include EACH mode's complex Maxwell coefficient.
        terms = _complex(edges, "L" + electric)
        total = _complex(fine, "L" + electric)[0]
        if abs(terms.sum() - total) > 1e-10 * max(np.sum(abs(terms)), 1e-30):
            raise ValueError("Edge sum disagrees with area balance")
        boundary[component] = terms
    groups = pd.read_csv(data.source / "group_totals.csv")
    groups = groups[(groups.z_nm == 0) & (groups.refinement == fine.refinement.iloc[0])]
    if len(groups) != 1:
        raise ValueError("Expected one matching group total")
    for group in ("A", "B"):
        actual = _complex(edges[edges.group == group], "Ly").sum()
        if not np.isclose(actual, _complex(groups, group + "_group_Ly")[0], rtol=1e-9, atol=1e-12):
            raise ValueError("Group sum disagrees with edges")
    phase = np.conj(reference) / abs(reference)
    displayed_field = {**field, "grid_electric_vm": field["grid_electric_vm"] * phase}
    displayed_boundary = {key: value * phase for key, value in boundary.items()}
    figures = build_edge_figures(geometry, edges.assign(z_over_H=0.), fine, data.meta, displayed_field, displayed_boundary, reference,
                                dimensional=True)
    figures.pop("s4_edges_py_Ex_reference").clear()
    gradient = load_gradient(data)
    direct = _complex(fine, "Ey_area")[0]
    derivative = _complex(gradient[(gradient.method == "python") & (gradient.level == fine.refinement.iloc[0])], "Ey_gradient")[0]
    formulas = [r"$\langle E_y\rangle=\frac{1}{A_{\mathrm{cell}}}\int_{\mathrm{cell}}E_y(\mathbf{r})\,d^2r$",
                r"$\langle E_y\rangle=-\frac{i}{\epsilon_0 A_{\mathrm{cell}}}\sum_m\frac{a_m}{\omega_m}"
                r"\int_{\mathrm{cell}}\frac{1}{\epsilon}\partial_xH_{z,m}\,d^2r$",
                r"$\langle E_y\rangle=\frac{i\Delta(1/\epsilon)}{\epsilon_0 A_{\mathrm{cell}}}"
                r"\sum_m\frac{a_m}{\omega_m}\oint_S H_{z,m}n_x\,dl$"]
    fig = figures["s4_edges_py"]
    fig._suptitle.set_text(fig._suptitle.get_text() + " | common phase")
    for ax, formula, raw_value in zip(fig.axes, formulas, (direct, derivative, boundary["x"].sum())):
        value = raw_value * phase
        value_label = f"{value.real:.5g}" if abs(value.imag) <= 1e-12*max(abs(value), 1e-30) else f"{value.real:.5g}{value.imag:+.5g}i"
        ax.text(.5, -.29, formula + "\n" + rf"$={value_label}\;\mathrm{{V/m}}$",
                transform=ax.transAxes, ha="center", va="top", fontsize=10)
    normalized = edges[["hole", "edge", "group"]].copy()
    normalized["py_contribution_re"] = (boundary["x"] / direct).real
    normalized["py_contribution_im"] = (boundary["x"] / direct).imag
    return figures, normalized


def write_analysis(data, normalized):
    logs = data.session / "80_logs"
    reference = _complex(data.fine, "Ey_area")[0]
    gradient = load_gradient(data)
    derivative = _complex(gradient[(gradient.method == "python") & (gradient.level == data.fine.refinement.iloc[0])], "Ey_gradient")[0]
    gradient_native = _complex(gradient[(gradient.method == "native") & (gradient.level == 12)], "Ey_gradient")[0]
    gradient_coarse = _complex(gradient[(gradient.method == "python") & (gradient.level == 1)], "Ey_gradient")[0]
    boundary = _complex(data.fine, "Ly")[0]
    phase = np.conj(reference) / abs(reference)
    comparison = [{"method": name, "z_nm": 0., "Ey_average": value,
                   "display_Ey_average": value*phase, "phase_rotation": phase, "unit": "V/m"}
                  for name, value in (("direct_Ey_area", reference), ("weighted_dxHz_area", derivative), ("hole_boundary_sum", boundary))]
    s4.write_csv(logs / "method_comparison.csv", comparison)
    native = pd.read_csv(data.source / "native_balance.csv")
    native_fine = native[native.order == native.order.max()]
    native_reference = _complex(native_fine, "Ey_area")[0]
    rows = []
    for method, frame, divisor in (("python", data.samples, reference), ("native", native, native_reference)):
        for _, row in frame.iterrows():
            terms = {key: complex(row[key + "_re"], row[key + "_im"]) / divisor
                     for key in ("Ey_area", "Ly", "Oy", "Zy", "predicted_Ey")}
            rows.append({"method": method, "z_nm": row.get("z_nm", 0.),
                         "level": row.get("refinement", row.get("order")),
                         "trace": row.get("trace", "point_interpolation"),
                         "reference_Ey_area": divisor, **terms,
                         "boundary_residual": abs(terms["Ly"] - terms["Ey_area"]),
                         "full_residual": abs(terms["predicted_Ey"] - terms["Ey_area"])})
    s4.write_csv(logs / "normalized_integrals.csv", rows)
    total = complex(normalized.py_contribution_re.sum(), normalized.py_contribution_im.sum())
    expected = _complex(data.fine, "Ly")[0] / reference
    if not np.isclose(total, expected, rtol=1e-9, atol=1e-12):
        raise ValueError("Reported edge sum changed")
    full = _complex(data.fine, "predicted_Ey")[0] / reference
    levels = data.samples[data.samples.z_nm == 0].sort_values("refinement")
    area_change = abs(_complex(levels, "Ey_area")[-1] - _complex(levels, "Ey_area")[-2]) / abs(reference)
    native_error = abs(reference / native_reference - 1)
    lines = ["# C6v 对照：面积与边界积分的图示分析", "",
             "主图顺序：Ey二维场、边界几何、18边有符号贡献；只显示z=0 nm。",
             "主图位于 `../10_overview/s4_edges_py.png`，场、边贡献及三个计算结果均使用V/m，不作面积归一化。",
             "保留源数据的实际幅度，三种计算共同旋转相位，使直接Ey面积平均为正实数；不除以面积平均的大小。",
             f"统一旋转角为 {-np.degrees(np.angle(reference)):.8g}°；旋转系数的模严格为1，原始相对相位和相对误差不变。",
             "色彩/柱高显示这个公共相位下的实部；图中公式值使用同一旋转，原始值和显示值均保存在method_comparison.csv。",
             "先前的归一化诊断图 s31_diagnostic_py 保留原样，独立作为相对误差参考。", "",
             "| 主图 | 计算方法 | 原始值（V/m） | 图中公共相位值（V/m） |", "|---|---|---:|---:|",
             f"| 1 | 直接Ey面积积分/A_cell | {reference:.9g} | {reference*phase:.9g} |",
             f"| 2 | 加权Hz横向导数面积积分及Maxwell系数 | {derivative:.9g} | {derivative*phase:.9g} |",
             f"| 3 | 当前18条孔边界贡献加和 | {boundary:.9g} | {boundary*phase:.9g} |", "",
             "图2显式使用d(laginterp(2,ewfd.Hz),x)重建导数；并非用图3或Ey减厚度项反推。",
             f"图2原生12阶结果为 {gradient_native:.9g} V/m；Python细档与其相对差 {abs(derivative-gradient_native)/abs(gradient_native):.3%}。",
             f"图2两档Python变化/细档结果为 {abs(derivative-gradient_coarse)/abs(derivative):.3%}，导数积分仍有采样敏感性。",
             "三维模型的直接Ey还包含厚度导数项，因此后两式是待比较的二维关系，不能预设三值相同。", "",
             f"- 直接面积平均：{reference:.12g} V/m；面积积分：{_complex(data.fine, 'Ey_integral')[0]:.12g} V·m。",
             f"- 18边和/面积：**{total:.8g}**；与1的差的模为 **{abs(total-1):.6g}（{abs(total-1):.3%}）**。",
             f"- 加入外边界与厚度项：(B+O+Z)/A = {full:.8g}，完整复数残差 **{abs(full-1):.3%}**。",
             f"- 两档面积积分变化/细档面积：{area_change:.3%}；细档面积与原生面积差：{native_error:.3%}。", "",
             "这些数据不支持仅保留孔边界项的等式在此三维结果上成立。厚度项不可忽略；",
             "加入厚度和外边界后残差减小，但不能据此宣称已完成有限元网格收敛或公式验证。", "",
             "## 原生积分核验", "",
             f"主图不作归一化。以下相对误差分析单独以原生面积平均 {native_reference:.12g} V/m 为基准。",
             "默认继承、空气侧、介质侧都是已有结果；不能因为某一种残差较小就选择它作为真值。", "",
             "| 边界迹 | B/A | O/A | Z/A | 完整复数残差 |",
             "|---|---:|---:|---:|---:|"]
    for _, row in native_fine.iterrows():
        ratios = [complex(row[key + "_re"], row[key + "_im"]) / native_reference for key in ("Ly", "Oy", "Zy")]
        lines.append(f"| {row.trace} | {ratios[0]:.7g} | {ratios[1]:.7g} | {ratios[2]:.7g} | {abs(sum(ratios)-1):.3%} |")
    fd = pd.read_csv(data.source / "native_thickness_fd.csv").sort_values("step_nm").head(1)
    z_fd = _complex(fd, "Zy_fd")[0] / native_reference
    z_native = _complex(native_fine, "Zy")[0] / native_reference
    fd_errors = abs((_complex(native_fine, "Ly") + _complex(native_fine, "Oy")) / native_reference + z_fd - 1)
    lines += ["", "上表厚度项采用逐单元Lagrange重建；独立截面差分提供额外交叉检查。",
              f"最小差分步长 {fd.step_nm.iloc[0]:g} nm 的 Z/A = {z_fd:.7g}；与重建Z的相对差 {abs(z_fd-z_native)/abs(z_native):.3%}。",
              f"使用该差分Z后，各边界迹的完整残差范围为 {fd_errors.min():.3%}–{fd_errors.max():.3%}。"]
    lines += ["", "## 对称态与误差解释", "",
              f"两模式频率差 {data.audit['frequency_split_ghz']:.6g} GHz；组合态并非精确单频本征态。",
              "采用已核验的常数复系数；每个模式用自身复频率计算边界贡献，再相加。",
              "场采样重现原34816点面积积分；没有重新拟合相位、逐点强制对称或把边界和强制设为1。",
              "积分阶数检查衡量当前离散场上的数值积分，不能替代网格收敛；裸边界差也不能全部称作数值误差。",
              "z=50/100 nm在原始数据和 normalized_integrals.csv 中仅作参考。",
              "Hz导数的两档Python和4/8/12阶原生核验见hz_gradient_integrals.csv；各模式复频率分量见hz_gradient_components.csv。", "",
              "## 来源与复现", "",
              "源积分：`../results/symmetry_decomposition_v1/`；既有结果全部保持原始字节。",
              "新增 area_field_z0.npz 为读取保存模型得到的组合态Ex/Ey积分样本及显示网格，mode_idx=-1表示组合态。",
              "field_extraction.json记录复系数、源审计哈希及缓存哈希；normalized_edge_contributions.csv为18边完整复数贡献。",
              "重新绘图无需COMSOL：", "", "```powershell",
              f'uv run python -B -m scripts.analysis.plot_boundary_control "{data.session}"', "```", ""]
    (logs / "plot_analysis.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session")
    parser.add_argument("--extract-field", action="store_true", help="Read saved MPH to create missing field cache; never solve")
    parser.add_argument("--extract-gradient", action="store_true", help="Read saved MPH to integrate weighted Hz x-derivatives; never solve")
    args = parser.parse_args(argv)
    data = load_control(args.session)
    before = {p: s4.digest(p) for p in (data.session / "results").rglob("*") if p.is_file()}
    if args.extract_field:
        extract_field(data)
    if args.extract_gradient:
        extract_gradient(data)
    figures, normalized = control_figures(data)
    overview = data.session / "10_overview"
    overview.mkdir(exist_ok=True)
    try:
        for name, fig in figures.items():
            for extension in ("png", "svg"):
                path = overview / f"{name}.{extension}"
                with mpl.rc_context({"svg.fonttype": "none"}):
                    fig.savefig(path, dpi=300, facecolor="white")
                print(path, flush=True)
        normalized.to_csv(data.session / "80_logs/normalized_edge_contributions.csv", index=False)
        dimensional_edges = normalized.rename(columns={"py_contribution_re": "Ey_contribution_re", "py_contribution_im": "Ey_contribution_im"}).copy()
        raw = (normalized.py_contribution_re.to_numpy() + 1j*normalized.py_contribution_im.to_numpy()) * _complex(data.fine, "Ey_area")[0]
        dimensional_edges["Ey_contribution_re"], dimensional_edges["Ey_contribution_im"] = raw.real, raw.imag
        reference = _complex(data.fine, "Ey_area")[0]
        display = raw * np.conj(reference) / abs(reference)
        dimensional_edges["display_Ey_contribution_re"], dimensional_edges["display_Ey_contribution_im"] = display.real, display.imag
        dimensional_edges["unit"] = "V/m"
        dimensional_edges.to_csv(data.session / "80_logs/edge_contributions_vm.csv", index=False)
        write_analysis(data, normalized)
    finally:
        for fig in figures.values():
            fig.clear()
    if not all(s4.digest(p) == sha for p, sha in before.items()):
        raise RuntimeError("Existing source results changed")
    s4.write_json(data.session / "80_logs/plot_verification.json", {
        "source_files_unchanged": len(before), "source_sha256": {p.relative_to(data.session).as_posix(): sha for p, sha in before.items()},
        "new_eigensolves": 0, "field_cache_reused": not args.extract_field, "gradient_cache_reused": not args.extract_gradient,
        "source_analysis": "constant complex p-pair combination; each mode retains its own frequency"})


if __name__ == "__main__":
    main()
