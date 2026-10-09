"""Read saved S4 fields and audit Hz derivatives; never mesh, solve or save MPH."""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np

from comsol_workflow.boundary_integrals import EPS0
from scripts.run_main import run_boundary_analysis as s4

METHODS = ["direct_dx", "lag1", "lag2", "lag3", "fd_2nm", "fd_1nm", "fd_0.5nm"]
STEPS_NM = (2., 1., .5)


def contains(points, polygon):
    """Convex polygon membership, independent of winding, SI coordinates."""
    points, polygon = np.asarray(points), np.asarray(polygon)
    edge = np.roll(polygon, -1, axis=0) - polygon
    rel = points[:, None, :] - polygon[None, :, :]
    cross = edge[None, :, 0] * rel[:, :, 1] - edge[None, :, 1] * rel[:, :, 0]
    return np.all(cross >= -1e-28, axis=1) | np.all(cross <= 1e-28, axis=1)


def distance(points, polygons):
    result = np.full(len(points), np.inf)
    for polygon in polygons:
        for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
            edge = b - a
            t = np.clip((points-a) @ edge / (edge @ edge), 0, 1)
            result = np.minimum(result, np.linalg.norm(points-a-t[:, None]*edge, axis=1))
    return result


def material_ids(points, holes):
    ids = np.zeros(len(points), int)
    for n, hole in enumerate(holes, 1):
        ids[contains(points, hole)] = n
    return ids


def fd_mask(points, outer, holes, step):
    plus, minus = points + [step, 0.], points - [step, 0.]
    base = material_ids(points, holes)
    return (contains(points, outer) & contains(plus, outer) & contains(minus, outer)
            & (distance(points, [outer, *holes]) > step)
            & (material_ids(plus, holes) == base) & (material_ids(minus, holes) == base))


def central_difference(plus, minus, step):
    if step <= 0:
        raise ValueError("Finite-difference step must be positive")
    return (plus-minus)/(2*step)


def convert_derivatives(gradients, omegas, epsilon):
    """Input shape (mode, method, point); convert before coherent combination."""
    return -1j * gradients / (np.asarray(omegas)[:, None, None]*EPS0*epsilon[None, None, :])


def coherent_pair(values, coefficients):
    return np.einsum("m,m...->...", coefficients, values)


def metrics(reference, prediction, weights, mask, area):
    w, ref, pred = weights[mask], reference[mask], prediction[mask]
    if not len(w) or w.sum() <= 0:
        return None
    residual = pred-ref
    denom = np.sum(w*abs(ref)**2)
    return {"samples": len(w), "area_fraction": float(w.sum()/area),
            "Ey_mean_fullcell": w@ref/area, "Hz_Ey_mean_fullcell": w@pred/area,
            "residual_mean_fullcell": w@residual/area,
            "Ey_mean_abs_fullcell": float(w@abs(ref)/area),
            "relative_L2": float(np.sqrt(np.sum(w*abs(residual)**2)/denom)) if denom > 0 else None}


def extra_points(outer, holes):
    points, labels = [], []
    # Separate horizontal cuts through A and B holes, avoiding exact interfaces.
    xs = np.linspace(outer[:, 0].min(), outer[:, 0].max(), 501)
    for label, y in ((1, 0.), (2, .22e-6)):
        xy = np.column_stack([xs, np.full(len(xs), y)])
        keep = contains(xy, outer) & (distance(xy, [outer, *holes]) > .2e-9)
        points.extend(xy[keep]); labels.extend([label]*int(keep.sum()))
    probe_rows = []
    for hole_id, hole in enumerate(holes, 1):
        for edge_id, (a, b) in enumerate(zip(hole, np.roll(hole, -1, axis=0)), 1):
            edge = b-a
            normal = np.array([-edge[1], edge[0]])/np.linalg.norm(edge)
            # Orient normal from air to dielectric for unambiguous probe-side labels.
            midpoint = (a+b)/2
            if contains((midpoint+normal*.1e-9)[None, :], hole)[0]:
                normal *= -1
            for fraction in (.25, .5, .75):
                for offset_nm in (.5, 1., 2., 5., 10.):
                    for side in (-1, 1):
                        xy = a+fraction*edge+side*offset_nm*1e-9*normal
                        if not contains(xy[None, :], outer)[0]:
                            continue
                        probe_rows.append({"extra_index": len(points), "hole": hole_id,
                                           "edge": edge_id, "fraction": fraction,
                                           "offset_nm": offset_nm, "side": "air" if side < 0 else "dielectric"})
                        points.append(xy); labels.append(3)
    return np.asarray(points), np.asarray(labels), probe_rows


def source_hashes(data):
    paths = [data.case / p for p in data.audit["source_sha256"]]
    paths += [data.source / "audit.json", data.session / "80_logs/area_field_z0.npz",
              s4.ROOT / "scripts/parameter.json"]
    return {str(p): s4.digest(p) for p in paths}


def extract(data, output):
    logs = output / "80_logs"
    audit_path, cache_path = logs / "audit.json", logs / "local_fields.npz"
    if audit_path.exists() or cache_path.exists():
        raise FileExistsError("Existing extraction is immutable; use --report-only or a new task directory")
    cache = np.load(data.session / "80_logs/area_field_z0.npz")
    area_xy, weights = cache["xy_m"], cache["weights_m2"]
    if not np.allclose(cache["outer_m"], data.outer, rtol=1e-12, atol=1e-20):
        raise ValueError("Saved quadrature geometry differs")
    geometry = s4.json_read(data.case / "99_config/realized_geometry.json")
    eps_air, eps_diel = [float(geometry["simulation_config"][k])**2
                         for k in ("air_refractive_index", "slab_refractive_index")]
    extra, labels, probes = extra_points(data.outer, data.holes)
    xy = np.concatenate([area_xy, extra])
    ids = material_ids(xy, data.holes)
    epsilon = np.where(ids > 0, eps_air, eps_diel)
    omegas = 2*np.pi*1e12*(data.modes.loc[data.indices, "re"].to_numpy()
                          - 1j*data.modes.loc[data.indices, "im"].to_numpy())
    before = source_hashes(data)
    record = {"status": "starting", "source_case": str(data.case), "source_sha256": before,
              "eigensolves": 0, "mesh_runs": 0, "model_saved": False,
              "z_nm": 0., "quadrature_samples": len(area_xy), "extra_samples": len(extra),
              "mode_indices": data.indices, "coefficients": data.audit["coefficients"],
              "methods": METHODS, "steps_nm": STEPS_NM,
              "source_manifest": data.manifest, "source_geometry": geometry,
              "script_sha256": s4.digest(__file__),
              "phase_convention": "exp(-i omega t); all native phasors and frequencies conjugated",
              "normalization": "raw modes preserved; py uses existing fixed coefficients"}
    s4.write_json(audit_path, record)
    client = model = None
    try:
        record["preflight"] = s4.runtime_preflight()
        s4.write_json(audit_path, record)
        import mph
        print("Starting COMSOL 6.3 for saved-field evaluation only", flush=True)
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics license unavailable")
        model = client.load(str(data.case / "00_model/s4_gamma.mph"))
        params = model.java.param()
        identity = {k: float(params.evaluate(k)) for k in ("a", "H", "kx", "ky")}
        if abs(identity["kx"]) > 1e-12 or abs(identity["ky"]) > 1e-12:
            raise ValueError("Saved model is not at Gamma")
        if not np.isclose(identity["a"], data.manifest["constants"]["A"]*1e-6, rtol=1e-12, atol=0):
            raise ValueError("Lattice identity mismatch")
        if not np.isclose(identity["H"], 200e-9, rtol=1e-12, atol=0):
            raise ValueError("Expected the approved 200 nm slab")
        freq = np.asarray(model.java.result().numerical("gev1").getReal()).T
        if not np.allclose(freq[:, :2], data.modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Saved eigenfrequency identity mismatch")
        record["loaded_parameters_SI"] = identity
        record["loaded_frequencies_THz"] = freq[:, :2].tolist()
        record["status"] = "extracting"
        s4.write_json(audit_path, record)
        sampler = s4.Sampler(SimpleNamespace(model=model))
        fields, derivatives = [], []
        valid = np.ones((len(METHODS), len(xy)), bool)
        for n, step_nm in enumerate(STEPS_NM, 4):
            valid[n] = fd_mask(xy, data.outer, data.holes, step_nm*1e-9)
        expressions = ["ewfd.Ex", "ewfd.Ey", "ewfd.Hz", "d(ewfd.Hz,x)/(1[A/m^2])"]
        expressions += [f"d(laginterp({p},ewfd.Hz),x)/(1[A/m^2])" for p in (1, 2, 3)]
        for mode in data.indices:
            print(f"Mode {mode}: E, H and four derivative evaluations at {len(xy)} points", flush=True)
            values = sampler.plane(mode, xy, 0., expressions)
            fields.append(values[:, :3])
            grad = np.full((len(METHODS), len(xy)), np.nan+1j*np.nan)
            grad[:4] = values[:, 3:].T
            for n, step_nm in enumerate(STEPS_NM, 4):
                h, mask = step_nm*1e-9, valid[n]
                paired = np.concatenate([xy[mask]+[h, 0.], xy[mask]-[h, 0.]])
                print(f"Mode {mode}: centered difference h={step_nm:g} nm, {mask.sum()} valid points", flush=True)
                values_h = sampler.plane(mode, paired, 0., ["ewfd.Hz"])[:, 0]
                plus, minus = np.split(values_h, 2)
                grad[n, mask] = central_difference(plus, minus, h)
            derivatives.append(grad)
        fields, derivatives = np.asarray(fields), np.asarray(derivatives)
        combined = coherent_pair(fields, data.coefficients)
        old = cache["electric_vm"]
        relative = float(np.linalg.norm(combined[:len(area_xy), :2]-old)/np.linalg.norm(old))
        record["electric_cache_relative_L2"] = relative
        if relative > 1e-9:
            raise ValueError(f"Direct E does not reproduce saved cache: {relative}")
        np.savez_compressed(cache_path, xy_m=xy, weights_m2=weights, fields=fields,
                            derivatives=derivatives, valid=valid, epsilon=epsilon, material_id=ids,
                            interface_distance_m=distance(xy, data.holes),
                            extra_labels=labels, outer_m=data.outer, holes_m=np.asarray(data.holes),
                            omegas=omegas, coefficients=data.coefficients,
                            mode_indices=data.indices, methods=np.asarray(METHODS))
        s4.write_csv(logs / "probe_locations.csv", probes)
        record.update(status="complete_pending_scientific_review", cache_sha256=s4.digest(cache_path))
        print(f"Cache verified against old E data; relative L2={relative:.3g}", flush=True)
    except Exception as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        try:
            if model is not None:
                client.remove(model)
        finally:
            record["source_unchanged"] = all(s4.digest(p) == sha for p, sha in before.items())
            s4.write_json(audit_path, record)
    if not record["source_unchanged"]:
        raise RuntimeError("Source files changed during validation")


def make_report(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    logs, figures = output / "80_logs", output / "10_overview"
    audit = s4.json_read(logs / "audit.json")
    if audit["status"] != "complete_pending_scientific_review" or not audit["source_unchanged"]:
        raise ValueError("Cannot report an incomplete extraction")
    if s4.digest(logs / "local_fields.npz") != audit["cache_sha256"]:
        raise ValueError("Field cache changed")
    z = np.load(logs / "local_fields.npz")
    weights, fields, deriv = z["weights_m2"], z["fields"], z["derivatives"]
    n, area = len(weights), float(weights.sum())
    converted = convert_derivatives(deriv, z["omegas"], z["epsilon"])
    py_e = coherent_pair(fields, z["coefficients"])[:, 1]
    py_h = coherent_pair(converted, z["coefficients"])
    mode_names = [f"raw_mode{int(i)}" for i in z["mode_indices"]] + ["py"]
    refs = [*fields[:, :, 1], py_e]
    preds = [*converted, py_h]
    air = z["material_id"][:n] > 0
    near = z["interface_distance_m"][:n] <= 5e-9
    regions = {"all": np.ones(n, bool), "air_bulk": air & ~near,
               "dielectric_bulk": ~air & ~near, "interface_air": near & air,
               "interface_dielectric": near & ~air}
    rows, sensitivity, common_rows = [], [], []
    common_fd = np.all(z["valid"][4:, :n], axis=0)
    for name, ref, pred in zip(mode_names, refs, preds):
        for j, method in enumerate(METHODS):
            for region, mask in regions.items():
                m = metrics(ref[:n], pred[j, :n], weights, mask & z["valid"][j, :n], area)
                if m is not None:
                    total_error = np.sum(weights[z["valid"][j, :n]] * abs((pred[j, :n]-ref[:n])[z["valid"][j, :n]])**2)
                    selected = mask & z["valid"][j, :n]
                    error_share = np.sum(weights[selected]*abs((pred[j, :n]-ref[:n])[selected])**2)/total_error if total_error > 0 else 0.
                    rows.append({"mode": name, "method": method, "region": region,
                                 "fraction_of_total_squared_error": float(error_share), **m})
            common = z["valid"][j, :n] & z["valid"][2, :n]
            m = metrics(pred[2, :n], pred[j, :n], weights, common, area)
            if m is not None:
                sensitivity.append({"mode": name, "method": method, "reference": "lag2", **m})
            common_result = metrics(ref[:n], pred[j, :n], weights, common_fd, area)
            common_difference = metrics(pred[2, :n], pred[j, :n], weights, common_fd, area)
            common_rows.append({"mode": name, "method": method,
                                "relative_L2_difference_from_lag2": common_difference["relative_L2"],
                                **common_result})
    s4.write_csv(logs / "method_region_comparison.csv", rows)
    s4.write_csv(logs / "derivative_sensitivity.csv", sensitivity)
    s4.write_csv(logs / "common_mask_comparison.csv", common_rows)
    phase = np.exp(-1j*np.angle(weights@py_e[:n]))
    probe = pd.read_csv(logs / "probe_locations.csv")
    probe_i = n+probe.extra_index.to_numpy(int)
    for j, method in enumerate(METHODS):
        probe[method+"_re_Vm"] = (phase*py_h[j, probe_i]).real
        probe[method+"_im_Vm"] = (phase*py_h[j, probe_i]).imag
    probe["Ey_re_Vm"] = (phase*py_e[probe_i]).real
    probe["Ey_im_Vm"] = (phase*py_e[probe_i]).imag
    probe.to_csv(logs / "interface_probe_values.csv", index=False)
    plt.rcParams.update({"font.size": 10, "xtick.direction": "in", "ytick.direction": "in",
                         "figure.dpi": 140, "axes.spines.top": True, "axes.spines.right": True})
    figures.mkdir(parents=True, exist_ok=True)
    xy = z["xy_m"]
    maps = [phase*py_e[:n], phase*py_h[2, :n], phase*(py_h[2, :n]-py_e[:n])]
    limit = max(float(np.max(abs(v.real))) for v in maps)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), layout="constrained")
    for ax, values, title in zip(axes, maps, ("Direct $E_y$", "$H_z$ derivative: lag2", "lag2 minus direct")):
        artist = ax.scatter(xy[:n, 0]*1e6, xy[:n, 1]*1e6, c=values.real, s=1.1,
                            cmap="RdBu_r", vmin=-limit, vmax=limit, rasterized=True)
        for poly in [z["outer_m"], *z["holes_m"]]:
            closed = np.vstack([poly, poly[0]])*1e6
            ax.plot(*closed.T, color="0.2", lw=.45)
        ax.set(aspect="equal", xlabel="$x$ (µm)", ylabel="$y$ (µm)", title=title)
    fig.colorbar(artist, ax=axes, label="Common-phase real field (V/m)", shrink=.85)
    fig.suptitle("$p_y$, Γ, z=0; quadrature samples, common phase and field scale")
    fig.savefig(figures / "local_field_comparison_py.png"); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
    colors = {"direct": "black", "lag1": "#4477aa", "lag2": "#ee7733", "lag3": "#228833", "fd_0.5nm": "#aa3377"}
    for col, label in enumerate((1, 2)):
        ids = n+np.flatnonzero(z["extra_labels"] == label)
        xx = xy[ids, 0]*1e6
        series = [("direct", py_e[ids])]+[(METHODS[j], py_h[j, ids]) for j in (1, 2, 3, 6)]
        for row, component in enumerate((np.real, np.imag)):
            ax = axes[row, col]
            for name, values in series:
                vv = component(phase*values).copy()
                # Break plotted curves at material transitions; never bridge a discontinuity.
                transition = np.r_[False, np.diff(z["material_id"][ids]) != 0]
                vv[transition] = np.nan
                ax.plot(xx, vv, label=name, color=colors[name], lw=1, alpha=.85)
            ax.set(xlabel="$x$ (µm)", ylabel=("Re" if row == 0 else "Im")+" field (V/m)",
                   title=f"y={xy[ids[0],1]*1e6:.2f} µm")
            ax.axhline(0, color=".7", lw=.5)
    axes[0, 0].legend(ncol=2, fontsize=8)
    fig.suptitle("$p_y$: identical points, phase and amplitude for every method")
    fig.savefig(figures / "local_field_cuts_py.png"); plt.close(fig)

    # The raw modes and py have different amplitudes; these dimensionless errors are comparable.
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), layout="constrained")
    for name in mode_names:
        rr = [r for r in rows if r["mode"] == name and r["region"] == "all" and r["method"] != "direct_dx"]
        axes[0].plot([r["method"] for r in rr], [r["relative_L2"] for r in rr], "o-", label=name)
    axes[0].set(ylabel="Weighted relative L2 residual", title="Residual against direct $E_y$")
    axes[0].tick_params(axis="x", rotation=30); axes[0].legend(fontsize=8)
    for method in ("lag1", "lag2", "lag3", "fd_0.5nm"):
        rr = [r for r in rows if r["mode"] == "py" and r["method"] == method and r["region"] != "all"]
        axes[1].plot([r["region"] for r in rr], [r["relative_L2"] for r in rr], "o-", label=method)
    axes[1].set(ylabel="Weighted relative L2 residual", title="$p_y$: location of residual")
    axes[1].tick_params(axis="x", rotation=25); axes[1].legend(fontsize=8)
    fig.savefig(figures / "derivative_sensitivity.png"); plt.close(fig)

    text = ["# S4 电场与磁场导数逐点验证", "",
            "只读保存的三维 Γ 模型，z=0；未划分网格、未求解、未保存 MPH。",
            "原始两模式各按其复频率转换，py 使用既有系数；所有源文件散列保持不变。",
            f"直接电场复现旧缓存的相对 L2 误差：{audit['electric_cache_relative_L2']:.6g}。",
            "有限差分排除跨材料及晶胞边界的点；其均值是有效子集对完整晶胞均值的贡献，不等于全晶胞积分。",
            "分区均值均以完整晶胞面积归一化；界面区为距离孔边不超过 5 nm 的两侧。",
            "", "## py：整体比较", "",
            "下表复均值使用同一相位，使完整晶胞直接 Ey 均值为正实数；原始复值见 CSV。", "",
            "| 方法 | 有效面积比例 | 直接 Ey (V/m) | 导数重建 Ey (V/m) | 相对 L2 残差 |",
            "|---|---:|---:|---:|---:|"]
    for r in rows:
        if r["mode"] == "py" and r["region"] == "all":
            text.append(f"| {r['method']} | {r['area_fraction']:.6f} | {phase*r['Ey_mean_fullcell']:.6g} | {phase*r['Hz_Ey_mean_fullcell']:.6g} | {r['relative_L2']:.6g} |")
    text += ["", "## 统一差分取点范围", "",
             "以下所有方法使用三个差分掩膜的交集，直接 Ey、材料权重和面积分节点完全相同。", "",
             "| 方法 | 相对 Ey 的 L2 残差 | 相对 lag2 的 L2 差异 |",
             "|---|---:|---:|"]
    for r in common_rows:
        if r["mode"] == "py":
            text.append(f"| {r['method']} | {r['relative_L2']:.6g} | {r['relative_L2_difference_from_lag2']:.6g} |")
    text += ["", "## 分区误差来源", "",
             "| 分区 | 面积比例 | 占 lag2 总平方误差比例 |",
             "|---|---:|---:|"]
    for r in rows:
        if r["mode"] == "py" and r["method"] == "lag2" and r["region"] != "all":
            text.append(f"| {r['region']} | {r['area_fraction']:.6g} | {r['fraction_of_total_squared_error']:.6g} |")
    text += ["", "## 判读限制", "",
             "重建阶数和差分步长的变化用于检查数值敏感性，不能替代网格收敛。有限差分也可能跨越有限元内部单元面。",
             "未把三维厚度项加入目标表达式；本报告只检验直接 Ey 与 Hz 的面内导数，不据此独立归因剩余差异。",
             "空间场形相似并不能自动保证较小的净积分具有相同相对精度；必须同时检查逐点残差与完整复积分。", ""]
    text += ["空间图只绘制实际面积分节点，空白缝隙不代表零场。截线图在材料切换处断开，避免连接不连续的取值。", ""]
    (logs / "report.md").write_text("\n".join(text), encoding="utf-8")
    s4.write_json(logs / "report_verification.json", {
        "status": "complete", "plot_phase_re": float(phase.real), "plot_phase_im": float(phase.imag),
        "report_script_sha256": s4.digest(__file__),
        "phase_changes_amplitude": False, "field_per_method_normalization": False,
        "figures": [p.name for p in sorted(figures.glob("*.png"))],
        "region_masks_partition": bool(np.all(sum(regions[k].astype(int) for k in regions if k != "all") == 1)),
        "area_m2": area, "quadrature_samples": n})
    print(f"Report and diagnostic figures: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_session", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    allowed = (s4.ROOT / "scripts/.out/unit_cell_2D").resolve()
    if output.parent != allowed or output == args.source_session.resolve():
        raise ValueError("Output must be a new independent task directly under unit_cell_2D")
    if s4.SESSION_ROOT != output:
        raise ValueError("COMSOL_WORKFLOW_S4_SESSION_ROOT must match this task directory")
    if not args.report_only:
        from scripts.analysis.plot_boundary_control import load_control
        data = load_control(args.source_session)
        extract(data, output)
    make_report(output)


if __name__ == "__main__":
    main()
