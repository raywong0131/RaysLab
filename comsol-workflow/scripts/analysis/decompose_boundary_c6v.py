"""Read-only C6v p-pair decomposition; never mesh, solve, or save an MPH."""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
import traceback

import numpy as np
import pandas as pd

from scripts.run_main import run_boundary_analysis as s4
from scripts.analysis.audit_s31_boundary import native_integral
from comsol_workflow.boundary_integrals import EPS0


def weighted_parity(h, hx, hy, weights):
    root = np.sqrt(np.asarray(weights, float))
    return s4.parity_errors(root*h, root*hx, root*hy)


def symmetry_coefficients(h, hx, hy, weights):
    """Two constant complex coefficients; no pointwise symmetrization."""
    h, hx, hy = (np.asarray(v, complex) for v in (h, hx, hy))
    w = np.asarray(weights, float)
    if h.ndim != 2 or h.shape[1] != 2 or hx.shape != h.shape or hy.shape != h.shape:
        raise ValueError("Expected two fields at identical physical sample points")
    if w.shape != (len(h),) or np.any(w <= 0) or not all(np.isfinite(v).all() for v in (h, hx, hy, w)):
        raise ValueError("Invalid field samples or weights")
    w = w / w.sum()
    scales = np.sqrt(w @ abs(h)**2)
    if np.any(scales == 0):
        raise ValueError("Rank-deficient p subspace")
    v, vx, vy = (a/scales for a in (h, hx, hy))
    gram = v.conj().T @ (w[:, None]*v)
    d, u = np.linalg.eigh(gram)
    if d[0] <= 1e-10*d[-1]:
        raise ValueError("Rank-deficient p subspace")
    whiten = u @ np.diag(d**-0.5)
    dx, dy = vx+v, vy-v
    cost = dx.conj().T @ (w[:, None]*dx) + dy.conj().T @ (w[:, None]*dy)
    reduced = whiten.conj().T @ cost @ whiten
    _, vectors = np.linalg.eigh((reduced+reduced.conj().T)/2)
    return (whiten @ vectors[:, 0])/scales


def scaled_terms(qx, qy, ox, oy, zx, zy, omegas, area, eps_air, eps_slab, coefficients):
    """Preserve the original complex frequency of EACH component."""
    keys = ("Lx", "Ly", "Ox", "Oy", "Zx", "Zy", "predicted_Ex", "predicted_Ey")
    total = dict.fromkeys(keys, 0j)
    for j, (omega, alpha) in enumerate(zip(omegas, coefficients)):
        part = s4.maxwell_terms(qx[j], qy[j], ox[j], oy[j], zx[j], zy[j],
                                omega, area, eps_air, eps_slab)
        for key in keys:
            total[key] += alpha*part[key]
    return total


def field_pair(sampler, indices, xy, z, expressions):
    return np.stack([sampler.plane(i, xy, z, expressions) for i in indices], axis=-1)


def python_integrals(sampler, indices, coefficients, omegas, outer, holes, settings, eps, output):
    area = abs(s4.polygon_area(outer))
    edges, outer_edges = s4.hole_edges(holes), s4.polygon_edges(outer)
    all_edges, groups, balances = [], [], []
    exprs = ["ewfd.Ex", "ewfd.Ey", "d(laginterp(2,ewfd.Hx),z)/(1[A/m^2])",
             "d(laginterp(2,ewfd.Hy),z)/(1[A/m^2])"]
    for z_nm in settings["z_nm"]:
        for refinement in (1, 2):
            xy, w, air = s4.area_rule(outer, holes, settings["area_order"], settings["area_subdivisions"]*refinement)
            values = field_pair(sampler, indices, xy, z_nm*1e-9, exprs)
            inv_eps = np.where(air, 1/eps[0], 1/eps[1])
            e = np.einsum("ncm,m->nc", values[:, :2], coefficients)
            a = w @ e / area
            dz = np.einsum("n,ncm->cm", w*inv_eps, values[:, 2:])
            hole_values, outer_values = [], []
            for edge in [*edges, *outer_edges]:
                points, weights = s4.line_rule(edge, settings["edge_order"]*refinement)
                hz = field_pair(sampler, indices, points, z_nm*1e-9, ["ewfd.Hz"])[:, 0]
                integrals = weights @ hz
                (hole_values if edge.hole else outer_values).append(integrals)
                if edge.hole:
                    p = 1j/(omegas*EPS0*area)*(1/eps[0]-1/eps[1])
                    weighted = np.dot(coefficients*p, integrals)
                    all_edges.append({"z_nm": z_nm, "refinement": refinement, "hole": edge.hole,
                                      "edge": edge.edge, "group": edge.group, "length_m": edge.length,
                                      "nx": edge.normal[0], "ny": edge.normal[1],
                                      "integral_hz": np.dot(coefficients, integrals),
                                      "Lx": -edge.normal[1]*weighted, "Ly": edge.normal[0]*weighted})
            hole_values, outer_values = np.asarray(hole_values), np.asarray(outer_values)
            normals, normal_outer = np.array([e.normal for e in edges]), np.array([e.normal for e in outer_edges])
            q, o = normals.T @ hole_values, normal_outer.T @ outer_values / eps[1]
            terms = scaled_terms(*q, *o, *dz, omegas, area, *eps, coefficients)
            chosen = [r for r in all_edges if r["z_nm"] == z_nm and r["refinement"] == refinement]
            ga = sum(r["Ly"] for r in chosen if r["group"] == "A")
            gb = sum(r["Ly"] for r in chosen if r["group"] == "B")
            groups.append({"z_nm": z_nm, "refinement": refinement, "A_group_Ly": ga,
                           "B_group_Ly": gb, "sum_Ly": ga+gb,
                           "cancellation_ratio": abs(ga+gb)/(abs(ga)+abs(gb))})
            balances.append({"z_nm": z_nm, "refinement": refinement, "samples": len(xy),
                             "Ex_area": a[0], "Ey_area": a[1], "mean_abs_Ey": w @ abs(e[:, 1])/area,
                             **terms})
            s4.write_csv(output / "edge_integrals.csv", all_edges)
            s4.write_csv(output / "group_totals.csv", groups)
            s4.write_csv(output / "area_comparison.csv", balances)
            print("Python quadrature", z_nm, refinement, s4.complex_columns(balances[-1]), flush=True)
    return balances


def run(case_dir, output):
    case_dir, output = Path(case_dir).resolve(), Path(output).resolve()
    if s4.SESSION_ROOT is None or not output.is_relative_to(s4.SESSION_ROOT.resolve()):
        raise ValueError("Explicit S4 session root is required; output must stay inside it")
    if output == case_dir or output.is_relative_to(case_dir):
        raise ValueError("Do not write into the original solved case")
    output.mkdir(parents=True, exist_ok=False)
    source = case_dir / "00_model/s4_gamma.mph"
    before = {p.relative_to(case_dir).as_posix(): s4.digest(p) for p in case_dir.rglob("*") if p.is_file()}
    record = {"source_case": str(case_dir), "source_sha256": before, "eigensolves": 0,
              "method": "constant complex p-subspace combination, Hz odd-x/even-y",
              "single_eigenstate": False, "status": "starting", "script_sha256": s4.digest(__file__),
              "source_files": s4.source_hashes()}
    s4.write_json(output / "audit.json", record)
    client = model = None
    try:
        manifest = s4.json_read(case_dir.parent / "99_config/s4_manifest.json")
        geometry = s4.json_read(case_dir / "99_config/realized_geometry.json")
        case, = [c for c in manifest["cases"] if c["case_id"] == case_dir.name]
        if s4.digest(case_dir.parent/case["snapshot"]) != case["snapshot_sha256"]:
            raise ValueError("Snapshot identity mismatch")
        if case["zeta"] != 1:
            raise ValueError("This explicit decomposition is restricted to the C6v control")
        modes = pd.read_csv(s4.data_directory(case_dir)/"eigenfrequencies.csv")
        indices = list(modes.index[(modes.p_weight >= manifest["s4"]["p_weight_minimum"]) & modes.is_valid])
        if len(indices) != 2:
            raise ValueError("Expected exactly two independently identified p candidates")
        omegas = 2*np.pi*1e12*(modes.loc[indices, "re"].to_numpy()-1j*modes.loc[indices, "im"].to_numpy())
        record["preflight"] = s4.runtime_preflight()
        import mph
        print("Loading existing C6v Gamma model; no solve", flush=True)
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics checkout failed")
        model = client.load(str(source))
        p = model.java.param()
        if any(abs(float(p.evaluate(k))) > 1e-12 for k in ("kx", "ky")):
            raise ValueError("Saved solution is not Gamma")
        a = float(p.evaluate("a"))
        if not np.isclose(a, manifest["constants"]["A"]*1e-6, rtol=1e-12, atol=0):
            raise ValueError("Lattice identity mismatch")
        freq = np.asarray(model.java.result().numerical("gev1").getReal()).T
        if not np.allclose(freq[:, :2], modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Frequency identity mismatch")
        sampler = s4.Sampler(SimpleNamespace(model=model))
        outer = np.asarray(geometry["outer"])*1e-6
        holes = [np.asarray(h)*1e-6 for h in geometry["holes"]]
        axis = np.linspace(-.3*a, .3*a, 17)
        xx, yy = np.meshgrid(axis, axis)
        xy = np.column_stack([xx.ravel(), yy.ravel()])
        fields = [field_pair(sampler, indices, xy*reflection, 0, ["ewfd.Hz"])[:, 0]
                  for reflection in ([1, 1], [-1, 1], [1, -1])]
        coefficients = symmetry_coefficients(*fields, np.ones(len(xy)))
        training = weighted_parity(*(h @ coefficients for h in fields), np.ones(len(xy)))
        # Independent material-conforming points over the entire unit cell.
        coords, weights, _ = s4.area_rule(outer, holes, 8, 2)
        validation_fields = [field_pair(sampler, indices, coords*reflection, 0, ["ewfd.Hz"])[:, 0]
                             for reflection in ([1, 1], [-1, 1], [1, -1])]
        coefficients *= s4.phase_and_norm(validation_fields[0] @ coefficients, weights, -coords[:, 0])
        validation = weighted_parity(*(h @ coefficients for h in validation_fields), weights)
        combined = validation_fields[0] @ coefficients
        freq_bar = np.dot(weights*np.conj(combined), validation_fields[0] @ (coefficients*omegas))/np.dot(weights, abs(combined)**2)
        freq_spread = np.sqrt(weights @ abs(validation_fields[0] @ (coefficients*(omegas-freq_bar)))**2 /
                              (weights @ abs(combined)**2))/abs(freq_bar)
        record.update(mode_indices=[int(i) for i in indices], training_parity=training,
                      validation_parity=validation, validation_samples=len(coords),
                      parity_threshold=manifest["s4"]["parity_tolerance"],
                      frequency_split_ghz=float(abs(omegas[0].real-omegas[1].real)/(2*np.pi*1e9)),
                      relative_frequency_spread=float(freq_spread),
                      coefficients=[s4.complex_columns({"mode_idx": int(i), "coefficient": v}) for i, v in zip(indices, coefficients)],
                      frequency_rule="Maxwell conversion uses each original omega; no common-frequency substitution")
        s4.write_json(output / "audit.json", record)
        print("Symmetry decomposition", record["training_parity"], record["validation_parity"], flush=True)
        if max(validation.values()) >= record["parity_threshold"]:
            record["status"] = "symmetry_not_resolved"
            return
        eps = [float(geometry["simulation_config"][k])**2 for k in ("air_refractive_index", "slab_refractive_index")]
        python_integrals(sampler, indices, coefficients, omegas, outer, holes, manifest["s4"], eps, output)
        air_rows = []
        points, w, _ = s4.area_rule(outer, holes, 8, 4)
        area = abs(s4.polygon_area(outer))
        for z in (.6*float(p.evaluate("H_air")), .75*float(p.evaluate("H_air"))):
            e = field_pair(sampler, indices, points, z, ["ewfd.Ex", "ewfd.Ey"])
            integrals = np.einsum("n,ncm->cm", w, e)/area
            direct = integrals @ coefficients
            propagated = integrals @ (coefficients*np.exp(-1j*omegas/299792458*np.sqrt(eps[0])*z))
            air_rows.append({"z_m": z, "cx": direct[0], "cy": direct[1],
                             "cx_at_z0": propagated[0], "cy_at_z0": propagated[1]})
        s4.write_csv(output/"air_radiation.csv", air_rows)
        native_checks(model, sampler, indices, coefficients, omegas, outer, holes, eps, output)
        record["status"] = "complete_pending_scientific_review"
    except Exception as exc:
        record.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all((case_dir/p).is_file() and s4.digest(case_dir/p) == sha for p, sha in before.items())
        s4.write_json(output / "audit.json", record)
        if record["status"] == "complete_pending_scientific_review":
            write_report(output, record)
        print(str(output), flush=True)


def native_checks(model, sampler, indices, coefficients, omegas, outer, holes, eps, output):
    """Native integration of each component, followed by the SAME linear sum."""
    from jpype.types import JInt
    area = abs(s4.polygon_area(outer))
    comp = model.java.component("comp1")
    bottom = list(comp.selection("sel_bottom").entities(2))
    domains = list(comp.selection("sel_slab_layer").entities(3))
    diel_ids = list(comp.selection("sel_slab_dom").entities(3))
    if len(diel_ids) != 1:
        raise ValueError("Expected one connected slab domain")
    updown = np.asarray(comp.geom("geom1").getUpDown())

    def face(domain):
        matches = [b for b in bottom if domain in updown[:, b-1]]
        if len(matches) != 1:
            raise ValueError("Ambiguous bottom material face")
        return matches[0]

    native_area, native_z, edge_rows = [], [], []
    for order in (4, 8, 12):
        for mode in indices:
            v = native_integral(model, "s4cs_area", "IntSurface", bottom,
                                ["1", "ewfd.Ex", "ewfd.Ey"], ["m^2", "V*m", "V*m"], mode, order)
            if not np.isclose(v[0].real, area, rtol=1e-8, atol=0):
                raise ValueError("Native area mismatch")
            native_area.append({"mode_idx": int(mode), "order": order, "Ex": v[1]/area, "Ey": v[2]/area})
            total = np.zeros(3, complex)
            for domain in domains:
                epsilon = eps[1] if domain in diel_ids else eps[0]
                exprs = [f"side({domain},{e})/{epsilon}" for e in
                         ("ewfd.Hx", "d(laginterp(2,ewfd.Hx),z)", "d(laginterp(2,ewfd.Hy),z)")]
                total += native_integral(model, "s4cs_z", "IntSurface", [face(domain)], exprs,
                                         ["A*m", "A", "A"], mode, order)
            native_z.append({"mode_idx": int(mode), "order": order, "Hx0": total[0], "dzHx": total[1], "dzHy": total[2]})
        print("Native surface order", order, flush=True)
    s4.write_csv(output/"native_area_components.csv", native_area)
    s4.write_csv(output/"native_thickness_components.csv", native_z)
    edges = [*s4.hole_edges(holes), *s4.polygon_edges(outer)]
    for index, edge in enumerate(edges):
        selection = comp.selection().create(f"s4cs_edge{index}", "Ball")
        selection.set("entitydim", JInt(1))
        selection.set("condition", "intersects")
        midpoint = (edge.start+edge.end)/(2e-6)
        for key, v in zip(("posx", "posy", "posz"), [*midpoint, 0]):
            selection.set(key, float(v))
        selection.set("r", 1e-6)
        entities = list(selection.entities(1))
        if len(entities) != 1:
            raise ValueError("Ambiguous native line entity")
        exprs = ["1", "ewfd.Hz", f"side({face(diel_ids[0])},ewfd.Hz)"]
        if edge.hole:
            air_ids = list(comp.selection(f"sel_hole_dom_{edge.hole-1}").entities(3))
            if len(air_ids) != 1:
                raise ValueError("Ambiguous hole domain")
            exprs.append(f"side({face(air_ids[0])},ewfd.Hz)")
        for order in (4, 8, 12):
            for mode in indices:
                v = native_integral(model, "s4cs_line", "IntLine", entities, exprs,
                                    ["m", *(["A"]*(len(exprs)-1))], mode, order)
                if not np.isclose(v[0].real, edge.length, rtol=1e-8, atol=0):
                    raise ValueError("Native edge length mismatch")
                edge_rows.append({"mode_idx": int(mode), "order": order, "hole": edge.hole,
                                  "edge": edge.edge, "group": edge.group, "nx": edge.normal[0], "ny": edge.normal[1],
                                  "inherited": v[1], "dielectric": v[2], "air": v[3] if edge.hole else v[2]})
    s4.write_csv(output/"native_edge_components.csv", edge_rows)
    balances = []
    for order in (4, 8, 12):
        ar = [next(r for r in native_area if r["mode_idx"] == m and r["order"] == order) for m in indices]
        dr = [next(r for r in native_z if r["mode_idx"] == m and r["order"] == order) for m in indices]
        for trace in ("inherited", "air", "dielectric"):
            q, o = [], []
            for mode in indices:
                lines = [r for r in edge_rows if r["mode_idx"] == mode and r["order"] == order]
                q.append([sum(r[trace]*r[n] for r in lines if r["hole"]) for n in ("nx", "ny")])
                o.append([sum(r[trace]*r[n] for r in lines if not r["hole"])/eps[1] for n in ("nx", "ny")])
            dz = np.array([[r["dzHx"], r["dzHy"]] for r in dr]).T
            terms = scaled_terms(*np.asarray(q).T, *np.asarray(o).T, *dz, omegas, area, *eps, coefficients)
            ax = np.dot(coefficients, [r["Ex"] for r in ar])
            ay = np.dot(coefficients, [r["Ey"] for r in ar])
            balances.append({"order": order, "trace": trace, "Ex_area": ax, "Ey_area": ay, **terms,
                             "boundary_error_relative_Ey": abs(terms["Ly"]-ay)/abs(ay),
                             "full_error_relative_Ey": abs(terms["predicted_Ey"]-ay)/abs(ay)})
    s4.write_csv(output/"native_balance.csv", balances)
    # Independent positive-z section differences. Native z=0 uses explicit material sides.
    cut = model.java.result().dataset().create("s4cs_cut", "CutPlane")
    cut.set("data", "dset1")
    cut.set("quickplane", "xy")
    fd_rows = []
    dielectric_area = area-sum(abs(s4.polygon_area(h)) for h in holes)
    for mode, alpha, omega in zip(indices, coefficients, omegas):
        values = {0: next(r["Hx0"] for r in native_z if r["mode_idx"] == mode and r["order"] == 12)}
        for z_nm in (.5, 1, 2, 4):
            cut.set("quickz", f"{z_nm}[nm]")
            v = native_integral(model, "s4cs_cutint", "IntSurface", None,
                                ["1", f"ewfd.Hx/if(dom=={diel_ids[0]},{eps[1]},{eps[0]})", f"if(dom=={diel_ids[0]},1,0)"],
                                ["m^2", "A*m", "m^2"], mode, 12, "s4cs_cut")
            if not np.allclose(v[[0, 2]].real, [area, dielectric_area], rtol=1e-8, atol=0):
                raise ValueError("Cut-plane material-area mismatch")
            values[z_nm] = v[1]
        for h in (.5, 1, 2):
            dz = (-3*values[0]+4*values[h]-values[2*h])/(2*h*1e-9)
            fd_rows.append({"mode_idx": int(mode), "step_nm": h, "Zy_contribution": alpha*1j/(omega*EPS0*area)*dz})
    s4.write_csv(output/"native_thickness_fd.csv",
                 [{"step_nm": h, "Zy_fd": sum(r["Zy_contribution"] for r in fd_rows if r["step_nm"] == h)} for h in (.5, 1, 2)])
    print("Native trace and thickness diagnostics complete", flush=True)


def write_report(output, record):
    def value(row, name):
        return complex(row[name+"_re"], row[name+"_im"])

    native = pd.read_csv(output/"native_balance.csv")
    sample = pd.read_csv(output/"area_comparison.csv")
    finest = native[native.order == 12]
    area = value(finest.iloc[0], "Ey_area")
    lines = ["# S4 C6v 对照：对称性分解与积分核验", "",
             f"源案例：`{record['source_case']}`。无新网格、无本征求解、无 MPH 保存。",
             "", "## 对称性分解", "",
             f"两个 p 模式的频率差为 {record['frequency_split_ghz']:.8g} GHz。",
             f"训练点误差：{record['training_parity']}。独立全晶胞复核：{record['validation_parity']}。",
             f"原阈值保留为 {record['parity_threshold']}；未逐点强制镜像。",
             f"组合态的相对频率离散残差为 {record['relative_frequency_spread']:.6g}。",
             "这不是一个新的精确单频本征态。每个原始模式用自身复频率换算 Maxwell 积分，随后用同一复系数线性相加。所有电、磁场采用同一个 z=0 Hz RMS 归一化及公共相位。",
             "", "## z=0 的直接比较", "",
             f"原生面积平均 Ey = {area:.12g} V/m。",
             "B 为孔边界项，O 为晶胞外边界项，Z 为厚度导数项；符号约定与既有 S4 程序一致。",
             "", "| 边界迹 | B/A | O/A | Z/A | 完整相对残差 |", "|---|---:|---:|---:|---:|"]
    for _, r in finest.iterrows():
        lines.append(f"| {r['trace']} | {value(r, 'Ly')/area:.7g} | {value(r, 'Oy')/area:.7g} | {value(r, 'Zy')/area:.7g} | {r.full_error_relative_Ey:.5g} |")
    lines += ["", "## Python 积分检查", "", "| 节点数 | z (nm) | Ey/A_native(z=0) |", "|---:|---:|---:|"]
    for _, r in sample.iterrows():
        lines.append(f"| {int(r.samples)} | {r.z_nm:g} | {value(r, 'Ey_area')/area:.9g} |")
    lines += ["", "逐边及 A/B 分组见 `edge_integrals.csv`、`group_totals.csv`；三种原生边界迹原始分量见 `native_edge_components.csv`。",
              "空气平面辐射保存在 `air_radiation.csv`，与晶胞内截面平均分开。厚度截面差分见 `native_thickness_fd.csv`。",
              "", "## 证据边界", "",
              "本次只处理一个网格的已有解；积分阶数稳定不等于有限元网格收敛，不依据某种边界迹更接近目标就选其作为物理真值。",
              "边界项与面积项是否一致、三维厚度项是否可忽略，必须分别看表中完整复数结果。对照并未将剩余数值误差强制归零。",
              f"原案例的 {len(record['source_sha256'])} 个文件 SHA-256 前后一致：{record['source_unchanged']}。",
              "", "方法引用：[COMSOL 6.3 的 side 与 laginterp 定义](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_definitions.21.036.html)。", ""]
    (output/"report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(args.case_dir, args.output)


if __name__ == "__main__":
    main()
