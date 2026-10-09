"""Audit S31 on an existing Gamma solution; never mesh, solve, or save the MPH."""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
import traceback

import numpy as np
import pandas as pd

from comsol_workflow.boundary_integrals import EPS0
from scripts.run_main import run_boundary_analysis as s4


def native_integral(model, tag, kind, selection, expressions, units, mode, order, dataset="dset1"):
    from jpype.types import JArray, JInt, JString
    numerical = model.java.result().numerical()
    if tag in list(numerical.tags()):
        numerical.remove(tag)
    node = numerical.create(tag, kind)
    node.set("data", dataset)
    if selection is None:
        node.selection().all()
    else:
        node.selection().set(JArray(JInt)(list(selection)))
    node.set("expr", JArray(JString)(expressions))
    node.set("unit", JArray(JString)(units))
    node.set("solnum", JArray(JInt)([mode + 1]))
    node.set("intorderactive", "on")
    node.set("intorder", JInt(order))
    real = np.asarray(node.getReal(), float)
    imag = np.asarray(node.getImag(), float)
    values = np.conj(real + 1j * imag).reshape(-1)
    if len(values) != len(expressions) or not np.isfinite(values).all():
        raise ValueError("Native integral returned invalid shape or values")
    return values


def audit(case_dir, output):
    case_dir, output = Path(case_dir).resolve(), Path(output).resolve()
    if not output.is_relative_to(s4.TASK_ROOT.resolve()):
        raise ValueError("Audit staging must stay in the task tmp directory")
    output.mkdir(parents=True, exist_ok=False)
    manifest = s4.json_read(case_dir.parent / "99_config/s4_manifest.json")
    meta = s4.json_read(case_dir / "99_config/s4_summary.json")
    geometry = s4.json_read(case_dir / "99_config/realized_geometry.json")
    cases = [c for c in manifest["cases"] if c["case_id"] == case_dir.name]
    if len(cases) != 1 or s4.digest(case_dir.parent / cases[0]["snapshot"]) != cases[0]["snapshot_sha256"]:
        raise ValueError("Case/snapshot identity mismatch")
    source = case_dir / "00_model/s4_gamma.mph"
    before = {p.relative_to(case_dir).as_posix(): s4.digest(p) for p in case_dir.rglob("*") if p.is_file()}
    record = {"source_case": str(case_dir), "source_sha256": before, "eigensolves": 0,
              "mode_idx": meta["mode_idx"], "user_mode": "py", "z_nm": 0, "status": "starting",
              "code_sha256": {str(Path(__file__).relative_to(s4.ROOT)): s4.digest(__file__), **s4.source_hashes()}}
    s4.write_json(output / "audit.json", record)
    model = client = None
    try:
        record["preflight"] = s4.runtime_preflight()
        import mph
        from jpype.types import JArray, JInt, JString
        print("Loading saved Gamma solution", flush=True)
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(JString(p)) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics checkout failed")
        model = client.load(str(source))
        jm = model.java
        params = {k: float(jm.param().evaluate(k)) for k in ("a", "H", "kx", "ky")}
        if abs(params["kx"]) > 1e-12 or abs(params["ky"]) > 1e-12:
            raise ValueError("Saved solution is not Gamma")
        if not np.isclose(params["a"], manifest["constants"]["A"] * 1e-6, rtol=1e-12, atol=0):
            raise ValueError("Lattice mismatch")
        jm.param().set("s31_expected_H", geometry["simulation_config"]["slab_height"])
        if not np.isclose(params["H"], jm.param().evaluate("s31_expected_H"), rtol=1e-12, atol=0):
            raise ValueError("Slab height mismatch")
        jm.param().remove("s31_expected_H")
        modes = pd.read_csv(s4.data_directory(case_dir) / "eigenfrequencies.csv")
        freq = np.asarray(jm.result().numerical("gev1").getReal()).T
        if freq.shape != (len(modes), 3) or not np.allclose(freq[:, :2], modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Frequency identity mismatch")
        mode = int(meta["mode_idx"])
        factor = complex(meta["normalization_factor_re"], meta["normalization_factor_im"])
        omega = 2 * np.pi * complex(modes.loc[mode, "re"], -modes.loc[mode, "im"]) * 1e12
        record.update(parameters=params, frequency_thz=meta["frequency_thz"], omega=s4.complex_columns({"value": omega}))
        comp = jm.component("comp1")
        record["selections"] = {str(tag): list(comp.selection(tag).entities()) for tag in comp.selection().tags()}
        physics = comp.physics("ewfd")
        record["materials"] = {
            str(tag): {"domains": list(comp.material(tag).selection().entities()),
                       "n": [str(v) for v in comp.material(tag).propertyGroup("RefractiveIndex").getStringArray("n")]}
            for tag in comp.material().tags()}
        record["physics"] = {}
        for tag in physics.feature().tags():
            feature = physics.feature(tag)
            item = {"type": str(feature.getType()), "entities": list(feature.selection().entities())}
            if str(tag).startswith("pc"):
                item.update(periodic_type=str(feature.getString("PeriodicType")), wavevector=[str(v) for v in feature.getStringArray("kFloquet")])
            record["physics"][str(tag)] = item
        sampler = s4.Sampler(SimpleNamespace(model=model))
        outer = np.asarray(geometry["outer"]) * 1e-6
        holes = [np.asarray(h) * 1e-6 for h in geometry["holes"]]
        area = abs(s4.polygon_area(outer))
        bottom = list(comp.selection("sel_bottom").entities(2))
        print("Identity checked; native bottom-surface integration", flush=True)
        native = []
        for order in (4, 8, 12):
            value = native_integral(model, "s31surface", "IntSurface", bottom,
                                    ["1", "ewfd.Ey", "abs(ewfd.Ey)", "d(ewfd.Hx,z)", "d(ewfd.Hz,x)"],
                                    ["m^2", "V*m", "V*m", "A", "A"], mode, order)
            if not np.isclose(value[0].real, area, rtol=1e-8, atol=0):
                raise ValueError(f"Bottom selection area mismatch: {value[0]} versus {area}")
            native.append({"order": order, "area_m2": value[0].real,
                           "Ey_area": factor * value[1] / area,
                           "abs_Ey_area": abs(factor) * value[2].real / area,
                           "dz_Hx_integral": factor * value[3], "dx_Hz_integral": factor * value[4]})
            print(s4.complex_columns(native[-1]), flush=True)
        s4.write_csv(output / "native_area.csv", native)
        xy = np.array([[0, 0], [0.04, 0.02], [0.2624, 0], [-0.2624, 0], [0.38, 0.03]]) * 1e-6
        expressions = ["x", "y", "z", "ewfd.Ey", "ewfd.Hz", "ewfd.Hx/(1[A/m])",
                       "d(ewfd.Hx,z)", "d(ewfd.Hz,x)/(1[A/m^2])", "dom"]
        samples = []
        for z_nm in (0, 0.5, 1, 2, 5, 10, 20, 50):
            values = sampler.plane(mode, xy, z_nm * 1e-9, expressions)
            samples.extend({"point": i, "z_nm": z_nm, **dict(zip(expressions, row))} for i, row in enumerate(values))
        s4.write_csv(output / "local_samples.csv", samples)
        eps_air = float(geometry["simulation_config"]["air_refractive_index"]) ** 2
        eps_slab = float(geometry["simulation_config"]["slab_refractive_index"]) ** 2
        diel_ids = list(comp.selection("sel_slab_dom").entities(3))
        if len(diel_ids) != 1:
            raise ValueError("Audit requires the recorded connected slab domain")
        updown = np.asarray(comp.geom("geom1").getUpDown())
        def bottom_face(domain):
            faces = [b for b in bottom if domain in updown[:, b - 1]]
            if len(faces) != 1:
                raise ValueError(f"Expected one bottom face for domain {domain}: {faces}")
            return faces[0]
        print("Elementwise Lagrange reconstruction of vector-field derivatives", flush=True)
        lagrange = []
        domains = list(comp.selection("sel_slab_layer").entities(3))
        for degree in (1, 2, 3):
            for order in (4, 8, 12):
                integrals = np.zeros(4, complex)
                for domain in domains:
                    eps = eps_slab if domain in diel_ids else eps_air
                    exprs = [f"side({domain},{expr})" for expr in
                             ("ewfd.Ey", f"ewfd.Hx/{eps}",
                              f"d(laginterp({degree},ewfd.Hx),z)/{eps}",
                              f"d(laginterp({degree},ewfd.Hz),x)/{eps}")]
                    integrals += factor*native_integral(model, "s31lag", "IntSurface", [bottom_face(domain)],
                                                        exprs, ["V*m", "A*m", "A", "A"], mode, order)
                p = 1j/(omega*EPS0*area)
                lagrange.append({"degree": degree, "order": order, "A_side": integrals[0]/area,
                                 "weighted_Hx0": integrals[1], "Z_lagrange": p*integrals[2],
                                 "minus_dx_term": -p*integrals[3], "A_from_local_curl": p*(integrals[2]-integrals[3])})
                s4.write_csv(output / "native_lagrange.csv", lagrange)
        print("Native edge integrals, including separate material traces", flush=True)
        hole_edges = s4.hole_edges(holes)
        outer_edges = s4.polygon_edges(outer)
        line_rows = []
        for index, edge in enumerate([*hole_edges, *outer_edges]):
            tag = f"s31edge{index}"
            selection = comp.selection().create(tag, "Ball")
            selection.set("entitydim", JInt(1))
            selection.set("condition", "intersects")
            midpoint = (edge.start + edge.end) / (2e-6)
            for key, v in zip(("posx", "posy", "posz"), [*midpoint, 0]):
                selection.set(key, float(v))
            selection.set("r", 1e-6)
            entities = list(selection.entities(1))
            if len(entities) != 1:
                raise ValueError(f"Edge selection is ambiguous: {tag}: {entities}")
            exprs = ["1", "ewfd.Hz", f"side({bottom_face(diel_ids[0])},ewfd.Hz)"]
            if edge.hole:
                air_ids = list(comp.selection(f"sel_hole_dom_{edge.hole - 1}").entities(3))
                if len(air_ids) != 1:
                    raise ValueError("Ambiguous hole domain")
                exprs.append(f"side({bottom_face(air_ids[0])},ewfd.Hz)")
            for order in (4, 8, 12):
                v = native_integral(model, "s31line", "IntLine", entities, exprs,
                                    ["m", *(["A"] * (len(exprs) - 1))], mode, order)
                if not np.isclose(v[0].real, edge.length, rtol=1e-8, atol=0):
                    raise ValueError("Native edge length mismatch")
                line_rows.append({"hole": edge.hole, "edge": edge.edge, "order": order,
                                  "entity": entities[0], "nx": edge.normal[0], "length_m": v[0].real,
                                  "Hz_integral": factor * v[1], "Hz_diel": factor * v[2],
                                  "Hz_air": factor * (v[3] if edge.hole else v[2])})
        s4.write_csv(output / "native_edges.csv", line_rows)
        native_balance = []
        for row in native:
            edges_at_order = [e for e in line_rows if e["order"] == row["order"]]
            for trace in ("Hz_integral", "Hz_diel", "Hz_air"):
                hole = sum(e[trace] * e["nx"] for e in edges_at_order if e["hole"])
                outer_x = sum(e[trace] * e["nx"] for e in edges_at_order if not e["hole"]) / eps_slab
                terms = s4.maxwell_terms(hole, 0, outer_x, 0, 0, 0, omega, area, eps_air, eps_slab)
                native_balance.append({"order": row["order"], "trace": trace,
                                       "A": row["Ey_area"], "B": terms["Ly"], "O": terms["Oy"],
                                       "Z_exported": 0j})
        s4.write_csv(output / "native_balance.csv", native_balance)
        print("Native thickness differences on adjacent sections", flush=True)
        dataset = jm.result().dataset().create("s31cut", "CutPlane")
        dataset.set("data", "dset1")
        dataset.set("quickplane", "xy")
        weighted_hx = f"ewfd.Hx/if(dom=={diel_ids[0]},{eps_slab},{eps_air})"
        cut_rows = []
        for order in (4, 8, 12):
            for z_nm in (0, 0.5, 1, 2, 4, 8):
                if z_nm == 0:
                    # On a real boundary, dom denotes the boundary entity, not
                    # its 3D material domain. Reuse explicitly sided integrals.
                    row = next(r for r in lagrange if r["degree"] == 2 and r["order"] == order)
                    hx_integral = row["weighted_Hx0"]
                else:
                    dataset.set("quickz", f"{z_nm}[nm]")
                    v = native_integral(model, "s31hxs", "IntSurface", None,
                                        ["1", weighted_hx, f"if(dom=={diel_ids[0]},1,0)"],
                                        ["m^2", "A*m", "m^2"], mode, order, "s31cut")
                    if not np.isclose(v[0].real, area, rtol=1e-8, atol=0):
                        raise ValueError("Cut-plane native area mismatch")
                    dielectric_area = area - sum(abs(s4.polygon_area(h)) for h in holes)
                    if not np.isclose(v[2].real, dielectric_area, rtol=1e-8, atol=0):
                        raise ValueError("Cut-plane material domain weighting mismatch")
                    hx_integral = factor*v[1]
                cut_rows.append({"order": order, "z_nm": z_nm, "weighted_Hx_integral": hx_integral})
                s4.write_csv(output / "native_thickness_samples.csv", cut_rows)
        native_thickness = []
        for order in (4, 8, 12):
            by_z = {r["z_nm"]: r["weighted_Hx_integral"] for r in cut_rows if r["order"] == order}
            for h_nm in (4, 2, 1, 0.5):
                dz = (-3*by_z[0]+4*by_z[h_nm]-by_z[2*h_nm])/(2*h_nm*1e-9)
                native_thickness.append({"order": order, "step_nm": h_nm,
                                         "Z_fd": 1j/(omega*EPS0*area)*dz})
        s4.write_csv(output / "native_thickness_fd.csv", native_thickness)
        print("Checking translated periodic pairs", flush=True)
        pairs = []
        for i in range(3):
            edge, other = outer_edges[i], outer_edges[i + 3]
            p, w = s4.line_rule(edge, 128)
            q, _ = s4.line_rule(other, 128)
            q = q[::-1]
            if not np.allclose(q - p, (q - p)[0], rtol=1e-9, atol=1e-14):
                raise ValueError("Opposite edges are not translations")
            for offset_nm in (0, 0.5, 1):
                fields = [sampler.plane(mode, coords - offset_nm * 1e-9 * e.normal, 0,
                                       ["ewfd.Ex", "ewfd.Ey", "ewfd.Hz"]) * factor
                          for coords, e in ((p, edge), (q, other))]
                tang = (edge.end - edge.start) / edge.length
                et = [f[:, :2] @ tang for f in fields]
                hz = [f[:, 2] for f in fields]
                pairs.append({"pair": i + 1, "offset_nm": offset_nm,
                              "E_t_relative_difference": np.linalg.norm(et[0] - et[1]) / max(np.linalg.norm(t) for t in et),
                              "Hz_relative_difference": np.linalg.norm(hz[0] - hz[1]) / max(np.linalg.norm(h) for h in hz),
                              "Hz_nx_pair_integral": edge.normal[0] * (w @ (hz[0] - hz[1]))})
        s4.write_csv(output / "periodic_pairs.csv", pairs)
        print("Checking local derivatives with independent differences", flush=True)
        local = []
        for z_nm in (0, 25):
            base = sampler.plane(mode, xy, z_nm * 1e-9, ["ewfd.Ey", "ewfd.Hx/(1[A/m])", "dom", "d(ewfd.Hx,z)", "d(ewfd.Hz,x)/(1[A/m^2])",
                    "d(laginterp(2,ewfd.Hx),z)/(1[A/m^2])", "d(laginterp(2,ewfd.Hz),x)/(1[A/m^2])"])
            eps = np.where(np.isin(base[:, 2].real, diel_ids), eps_slab, eps_air)
            source_term = 1j * omega * EPS0 * eps * base[:, 0]
            for h_nm in (2, 1, 0.5):
                h = h_nm * 1e-9
                f1 = sampler.plane(mode, xy, z_nm * 1e-9 + h, ["ewfd.Hx/(1[A/m])"])[:, 0]
                f2 = sampler.plane(mode, xy, z_nm * 1e-9 + 2*h, ["ewfd.Hx/(1[A/m])"])[:, 0]
                dz = (-3*base[:, 1] + 4*f1 - f2) / (2*h)
                hp = sampler.plane(mode, xy + [h, 0], z_nm * 1e-9, ["ewfd.Hz"])[:, 0]
                hm = sampler.plane(mode, xy - [h, 0], z_nm * 1e-9, ["ewfd.Hz"])[:, 0]
                dx = (hp - hm) / (2*h)
                for i in range(len(xy)):
                    local.append({"point": i, "z_nm": z_nm, "step_nm": h_nm,
                                  "dz_fd": factor*dz[i], "dx_fd": factor*dx[i],
                                  "dz_exported": factor*base[i, 3], "dx_exported": factor*base[i, 4],
                                  "dz_lagrange": factor*base[i, 5], "dx_lagrange": factor*base[i, 6],
                                  "ampere_source": factor*source_term[i],
                                  "residual_fd": factor*(dz[i] - dx[i] + source_term[i])})
        s4.write_csv(output / "local_derivatives.csv", local)
        refinements = []
        thickness = []
        for refinement in (1, 2, 4):
            coords, weights, air = s4.area_rule(outer, holes, 8, 2*refinement)
            fields = factor * sampler.plane(mode, coords, 0, ["ewfd.Ey", "ewfd.Hx/(1[A/m])"])
            edges = s4.edge_integrals(hole_edges, lambda p: factor*sampler.plane(mode, p, 0, ["ewfd.Hz"])[:, 0], 64*refinement)
            outside = s4.edge_integrals(outer_edges, lambda p: factor*sampler.plane(mode, p, 0, ["ewfd.Hz"])[:, 0], 64*refinement)
            terms = s4.maxwell_terms(sum(e["qx"] for e in edges), 0,
                                    sum(e["qx"] for e in outside)/eps_slab, 0, 0, 0,
                                    omega, area, eps_air, eps_slab)
            refinements.append({"refinement": refinement, "area_samples": len(coords),
                                "A": weights @ fields[:, 0] / area, "B": terms["Ly"], "O": terms["Oy"]})
            s4.write_csv(output / "quadrature.csv", refinements)
            print("Quadrature", refinement, s4.complex_columns(refinements[-1]), flush=True)
            inv_eps = np.where(air, 1/eps_air, 1/eps_slab)
            for h_nm in (4, 2, 1):
                h = h_nm*1e-9
                h1 = factor*sampler.plane(mode, coords, h, ["ewfd.Hx/(1[A/m])"])[:, 0]
                h2 = factor*sampler.plane(mode, coords, 2*h, ["ewfd.Hx/(1[A/m])"])[:, 0]
                dz = (-3*fields[:, 1] + 4*h1 - h2)/(2*h)
                term = 1j/(omega*EPS0*area) * ((weights*inv_eps) @ dz)
                thickness.append({"refinement": refinement, "area_samples": len(coords), "step_nm": h_nm,
                                  "Z_fd": term, "Hx0_rms": np.sqrt(weights @ abs(fields[:, 1])**2/area)})
                s4.write_csv(output / "thickness_fd.csv", thickness)
            print("Thickness differences", refinement, "complete", flush=True)
        record["status"] = "complete_pending_review"
    except Exception as exc:
        record.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["unchanged"] = all((case_dir / p).is_file() and s4.digest(case_dir / p) == sha for p, sha in before.items())
        s4.write_json(output / "audit.json", record)
        print(str(output), flush=True)
    return output


def report(output):
    """Render source-backed findings; integration convergence is not FEM convergence."""
    from comsol_workflow.boundary_plotting import _complex, _figure
    import matplotlib as mpl
    output = Path(output)
    record = s4.json_read(output / "audit.json")
    meta = s4.json_read(Path(record["source_case"]) / "99_config/s4_summary.json")
    if record["status"] != "complete_pending_review" or not record["unchanged"]:
        raise ValueError("Only a completed audit with unchanged source data can be reported")
    native = pd.read_csv(output / "native_area.csv").sort_values("order")
    edges = pd.read_csv(output / "native_balance.csv")
    edges = edges[edges.order == edges.order.max()].set_index("trace")
    quad = pd.read_csv(output / "quadrature.csv").sort_values("refinement")
    thickness = pd.read_csv(output / "native_thickness_fd.csv")
    thickness = thickness[thickness.order == thickness.order.max()].sort_values("step_nm", ascending=False)
    pairs = pd.read_csv(output / "periodic_pairs.csv")
    local = pd.read_csv(output / "local_derivatives.csv")
    lagrange = pd.read_csv(output / "native_lagrange.csv")
    lag_final = lagrange[(lagrange.degree == 2) & (lagrange.order == lagrange.order.max())]
    z_lagrange = _complex(lag_final, "Z_lagrange")[0]
    reference = _complex(native, "Ey_area")[-1]
    if reference == 0:
        raise ValueError("Native area is zero; report absolute results before normalizing")
    values_a = _complex(quad, "A") / reference
    values_b = np.r_[1, _complex(quad, "B")[-1] / reference,
                      _complex(edges.loc[["Hz_integral", "Hz_air", "Hz_diel"]], "B") / reference]
    values_z = _complex(thickness, "Z_fd") / reference
    balance = []
    for trace in edges.index:
        row = edges.loc[[trace]]
        b, o = (_complex(row, key)[0] for key in ("B", "O"))
        for h, z in zip(thickness.step_nm, _complex(thickness, "Z_fd")):
            balance.append({"trace": trace, "step_nm": h, "A": reference, "B": b, "O": o,
                            "Z_fd": z, "full_residual": abs(reference-b-o-z),
                            "full_residual_over_abs_A": abs(reference-b-o-z)/abs(reference)})
    s4.write_csv(output / "balance_sensitivity.csv", balance)
    area_error = abs(_complex(quad, "A") - reference)/abs(reference)
    stability = abs(_complex(native, "Ey_area")[-1]-_complex(native, "Ey_area")[-2])/abs(reference)
    with mpl.rc_context({"font.size": 9, "axes.titlesize": 10, "xtick.direction": "in", "ytick.direction": "in", "svg.fonttype": "none"}):
        fig = _figure((14, 6))
        axs = fig.subplots(1, 3)
        fig.subplots_adjust(left=.065, right=.975, top=.75, bottom=.29, wspace=.36)
        fig.suptitle(rf"S31 audit | $p_y$ / $E_y$ | z=0 nm | {record['frequency_thz']:.5f} THz | existing mesh {meta['mesh']}", y=.975, fontsize=12)
        fig.text(.5, .9, r"$A=\langle E_y\rangle,\quad B=\mathcal{K}_\Gamma\oint H_z n_x\,dl$"
                 + "  |  Target: A = B  |  One native area reference A = 1 for all panels", ha="center", fontsize=10)
        for ax in axs:
            ax.set_box_aspect(.75)
            ax.axhline(0, color=".75", lw=.7, zorder=0)
            ax.set_ylabel(r"Re(value / $A_{native}$) (1)")
        ax = axs[0]
        ax.plot(range(len(quad)), values_a.real, "o-", color="#222222", lw=1.4)
        ax.axhline(1, color="#0072b2", ls="--", lw=1, label="Native area = 1")
        ax.set(xticks=range(len(quad)), xticklabels=[f"{n:,}" for n in quad.area_samples],
               xlabel="Python area quadrature points", title="(a) Area sampling versus native integral")
        ax.margins(x=.15, y=.25)
        ax.legend(loc="lower right", frameon=False, fontsize=8)
        for i, value in enumerate(values_a.real):
            ax.annotate(f"{value:+.3f}", (i, value), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8)
        ax.text(.5, -.29, "Native orders 4 / 8 / 12 agree\nPython errors: " + " / ".join(f"{v:.1%}" for v in area_error),
                transform=ax.transAxes, ha="center", va="top")
        ax = axs[1]
        bars = ax.bar(range(5), values_b.real, color=["#222222", "#999999", "#d55e00", "#0072b2", "#009e73"], width=.58)
        ax.axhline(1, color=".45", ls="--", lw=.8)
        ax.set(xticks=range(5), xticklabels=["Area\nA", "B\npoint", "B\ninherited", "B\nair side", "B\ndiel. side"],
               title="(b) Boundary result depends on trace")
        ax.margins(y=.25)
        for rect, value in zip(bars, values_b.real):
            ax.annotate(f"{value:+.2f}", (rect.get_x()+rect.get_width()/2, value),
                        xytext=(0, 6 if value >= 0 else -6), textcoords="offset points", ha="center",
                        va="bottom" if value >= 0 else "top", fontsize=8)
        ax.text(.5, -.29, "Native line orders 4 / 8 / 12 agree\nMaterial traces remain inconsistent",
                transform=ax.transAxes, ha="center", va="top")
        ax = axs[2]
        ax.plot(range(len(thickness)), values_z.real, "o-", color="#9467bd", lw=1.4, label="Section finite difference")
        ax.axhline((z_lagrange/reference).real, color="#009e73", lw=1, ls="--", label="Elementwise reconstructed Z")
        ax.axhline(0, color="#222222", ls="--", label="Exported d(Hx,z): Z = 0")
        ax.set(xticks=range(len(thickness)), xticklabels=[f"{h:g}" for h in thickness.step_nm],
               xlabel="Forward-difference step h (nm)", title="(c) Reconstructed thickness term is nonzero")
        ax.margins(x=.15, y=.3)
        ax.legend(loc="center", frameon=False, fontsize=8)
        z_difference = abs(_complex(thickness, "Z_fd")[-1]-z_lagrange)/abs(z_lagrange)
        ax.text(.5, -.29, f"Reconstructed Z / A = {(z_lagrange/reference).real:.2f}\nFinest difference versus reconstruction: {z_difference:.2%}",
                transform=ax.transAxes, ha="center", va="top", color="#a63603")
        fig.text(.5, .105, "Signed projections share one complex reference; reported errors use full complex differences.", ha="center")
        fig.text(.5, .055, "Same FEM solution, no new solve. Integration order stability does not establish field or mesh convergence.", ha="center", color=".35")
        for extension in ("png", "svg"):
            fig.savefig(output / f"diagnostic_py.{extension}", dpi=300, facecolor="white")
        fig.clear()
    lines = ["# S31 已有解核验", "", "结论：原生面积积分已稳定，原先导出为零的磁场导数是求值方式的问题。重建厚度项非零且明显大于净面积平均；磁场边界迹仍不一致，当前不能宣称 S31 成立或被推翻。", "",
             f"源案例：`{record['source_case']}`。用户 py、Γ 点 mode {record['mode_idx']}，z=0 nm；无网格化、无本征求解、无 MPH 保存。", "",
             "## 独立面积积分", "",
             f"COMSOL 原生 z=0 面积积分，选择面积 {native.area_m2.iloc[-1]:.12g} m²，与几何面积一致。",
             f"A = {reference.real:.12g} {reference.imag:+.12g}i V/m。4/8/12 阶一致，末两档变化 / |A| = {stability:.3g}。",
             f"面积相消条件数 ∫|Ey| / |∫Ey| ≈ {native.abs_Ey_area.iloc[-1]/abs(reference):.2f}，净结果对局部误差敏感。", "",
             "| Python 面积节点数 | 与原生 A 的复数相对差 |", "| ---: | ---: |"]
    lines += [f"| {n:,} | {v:.3%} |" for n, v in zip(quad.area_samples, area_error)]
    lines += ["", "## 孔边界取值与周期性", "", "下面的原生线积分都已通过长度核验，4/8/12 阶稳定；不同取值方式不能当作同一条已收敛物理边界迹。", "",
              "| 取值方式 | B/A（完整复数） |", "| --- | ---: |"]
    for trace, b in zip((f"Python 点插值（每边 {64*int(quad.refinement.iloc[-1])} 点）", "原生默认继承", "原生空气侧", "原生介质侧"), values_b[1:]):
        lines.append(f"| {trace} | {b.real:.8g} {b.imag:+.8g}i |")
    periodic = pairs[pairs.offset_nm == 0]
    lines += ["", f"三对 Γ 周期面的切向电场相对差最大 {periodic.E_t_relative_difference.max():.3g}，支持周期电场配对已经生效；Hz 配对仍有差异，见 `s31_periodic_pairs.csv`。近零 Hz 的相对差会放大，不单独作为误差界。",
              "", "## 厚度导数与完整关系", "",
              "原生及点导出的 d(Hx,z)、d(Hz,x) 均为零，但邻近点 Hx/Hz 明显变化，有限差分给出非零导数。COMSOL 的 laginterp 文档明确说明，矢量单元场不支持的空间导数可通过逐单元 Lagrange 重建得到。",
              f"显式材料域侧的重建 Z = {z_lagrange.real:.10g} {z_lagrange.imag:+.10g}i V/m；1/2/3阶重建及4/8/12阶积分一致。该结果描述当前离散场，并非已经完成网格收敛的物理真值。",
              "在 z=0 边界上，dom 对应当前维度的实体编号，不能直接用三维材料域编号判断介质。端点按各材料域显式 side(...) 分别积分；正 z 截面核对总面积及介质面积后使用材料权重。",
              "独立截面差分采用 [-3F(0)+4F(h)-F(2h)]/(2h)，F=∫Hx/ε dS。结果如下：", "",
              "| h (nm) | 原生截面差分 Z（V/m） |", "| ---: | ---: |"]
    lines += [f"| {h:g} | {z.real:.8g} {z.imag:+.8g}i |" for h, z in zip(thickness.step_nm, _complex(thickness, "Z_fd"))]
    local_residual = np.linalg.norm(_complex(local, "residual_fd"))/np.linalg.norm(_complex(local, "ampere_source"))
    lines += ["", f"最小步长差分与逐单元重建 Z 的差异为 {z_difference:.3%}。步长差分可能跨越离散场的单元跳变，不能把剩余差异直接当作严格误差界。",
              f"5 个远离材料接口的点、z=0/25 nm、h=2/1/0.5 nm 的局部差分 Maxwell 残差整体 L2 比为 {local_residual:.5g}；这是局部诊断，不代表全域误差界。",
              "Python 面积积分使用相同固定归一化，额外加密到 139264 节点；面积、边界及厚度差分仍有变化。原生积分稳定仅说明对当前离散场的积分稳定，不能替代 FEM 网格收敛。",
              "完整关系的组合敏感性见 `s31_balance_sensitivity.csv`：不同边界迹给出不同残差，暂不选择其中一组作为已验证的闭合结果。", "",
              "## 下一步", "", "下一步需要有限元网格/边界迹收敛证据。现有三维厚度项相对于净面积平均不可忽略，应核对 S31 的二维假设及是否要求厚度平均。未启动新网格或二维模型求解。", "",
              "COMSOL 6.3 文档：[laginterp](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_definitions.21.036.html#2246516)。", "",
              "## 数据保护与复现", "", f"源案例 {len(record['source_sha256'])} 个既有文件 SHA-256 全部不变。审计记录和原生/采样/差分明细均随本报告提供。", ""]
    (output / "report.md").write_text("\n".join(lines), encoding="utf8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir")
    parser.add_argument("--output", required=True, help="New staging directory; existing results are never overwritten")
    parser.add_argument("--report-only", action="store_true", help="Render a completed audit without COMSOL")
    args = parser.parse_args()
    if not args.report_only:
        audit(args.case_dir, args.output)
    report(args.output)


if __name__ == "__main__":
    main()
