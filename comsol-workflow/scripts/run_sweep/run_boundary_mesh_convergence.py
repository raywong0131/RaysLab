"""Remesh the saved S4 Gamma model, preserving physics and original data."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from types import SimpleNamespace
import time
import traceback

import numpy as np
import pandas as pd

from scripts.run_main import run_boundary_analysis as s4
from scripts.analysis import decompose_boundary_c6v as dec
from scripts.analysis.audit_s31_boundary import native_integral
from scripts.analysis.validate_boundary_local_fields import convert_derivatives, metrics, distance
from comsol_workflow.boundary_integrals import EPS0


def subspace_scores(reference, candidates, weights):
    """Weighted projection onto a reference p subspace; independent of pair rotation."""
    w = np.asarray(weights)/np.sum(weights)
    r, v = np.asarray(reference), np.asarray(candidates)
    gram = r.conj().T @ (w[:, None]*r)
    eigenvalues = np.linalg.eigvalsh(gram)
    if eigenvalues[0] <= 1e-10*eigenvalues[-1]:
        raise ValueError("Reference p subspace is rank deficient")
    overlaps = r.conj().T @ (w[:, None]*v)
    projected = np.sum(overlaps.conj()*np.linalg.solve(gram, overlaps), axis=0).real
    norms = w @ abs(v)**2
    if np.any(norms <= 0):
        raise ValueError("Zero candidate field")
    return np.clip(projected/norms, 0., 1.)


def complex_change(previous, current):
    difference = abs(current-previous)
    scale = max(abs(previous), abs(current))
    return difference, difference/scale if scale > 0 else 0.


def snapshot_for_mesh(source, mesh):
    import copy
    result = copy.deepcopy(source)
    result["mesh_size"] = int(mesh)
    return result


def mesh_stats(model):
    mesh = model.java.component("comp1").mesh("mesh1")
    kinds = [str(k) for k in mesh.getTypes()]
    result = {"auto_size": int(mesh.autoMeshSize()), "automatic": bool(mesh.isAutomatic()),
              "vertices": int(mesh.getNumVertex()),
              "elements": {k: int(mesh.getNumElem(k)) for k in kinds}}
    if "tet" in kinds:
        ids = np.asarray(mesh.getElemEntity("tet"), int)
        slab = list(model.java.component("comp1").selection("sel_slab_layer").entities(3))
        result["slab_layer_tetrahedra"] = int(np.isin(ids, slab).sum())
    return result


def explicit_frequencies(model):
    """Dimensionless expressions prevent automatic damping-unit changes after a solve."""
    node = model.java.result().numerical().create("s4mc_freq", "EvalGlobal")
    node.set("data", "dset1")
    node.set("expr", ["ewfd.omega/(2*pi*1[THz])", "ewfd.damp/(2*pi*1[THz])", "ewfd.Qfactor"])
    node.set("unit", ["1", "1", "1"])
    freq = np.asarray(node.getReal()).T[:, :3]
    # Qfactor is unsigned; do not change the original damping sign.
    if not np.allclose(np.abs(freq[:, 0]/(2*freq[:, 1])), freq[:, 2], rtol=1e-7, atol=1e-10):
        raise ValueError("Frequency/damping units inconsistent with Q")
    return freq


def native_gradient(model, indices, coefficients, omegas, eps, area, output):
    comp = model.java.component("comp1")
    bottom = list(comp.selection("sel_bottom").entities(2))
    domains = list(comp.selection("sel_slab_layer").entities(3))
    dielectric = list(comp.selection("sel_slab_dom").entities(3))
    updown = np.asarray(comp.geom("geom1").getUpDown())
    rows = []
    for order in (8, 12):
        totals = np.zeros(len(indices), complex)
        for j, mode in enumerate(indices):
            selected_area = 0j
            for domain in domains:
                faces = [b for b in bottom if domain in updown[:, b-1]]
                if len(faces) != 1:
                    raise ValueError("Ambiguous bottom material face")
                epsilon = eps[1] if domain in dielectric else eps[0]
                value = native_integral(model, "s4mc_gradient", "IntSurface", faces,
                    ["1", f"side({domain},d(laginterp(2,ewfd.Hz),x))/{epsilon}"],
                    ["m^2", "A"], mode, order)
                selected_area += value[0]
                totals[j] += value[1]
            if not np.isclose(selected_area, area, rtol=1e-8, atol=0):
                raise ValueError("Native material areas do not sum to the cell")
        ey = -1j*np.sum(coefficients*totals/omegas)/(EPS0*area)
        rows.append({"order": order, "Ey_gradient": ey})
    s4.write_csv(output / "native_gradient.csv", rows)
    return rows[-1]["Ey_gradient"]


def analyze(model, data, level, output, reference=None):
    logs = output / "80_logs"
    logs.mkdir(parents=True)
    sampler = s4.Sampler(SimpleNamespace(model=model))
    outer, holes = data.outer, data.holes
    area = abs(s4.polygon_area(outer))
    qc, wc, _ = s4.area_rule(outer, holes, 8, 2)
    qf, wf, air = s4.area_rule(outer, holes, 8, 4)
    freq = explicit_frequencies(model)
    if not np.isfinite(freq).all():
        raise ValueError("Nonfinite frequencies")
    candidate_h = np.column_stack([sampler.plane(i, qc, 0., ["ewfd.Hz"])[:, 0]
                                   for i in range(len(freq))])
    if reference is None:
        if not np.allclose(freq[:, :2], data.modes[["re", "im"]], rtol=1e-10, atol=1e-12):
            raise ValueError("Baseline eigenfrequency identity mismatch")
        indices = data.indices
        reference_basis = candidate_h[:, indices]
    else:
        reference_basis = reference["basis"]
        scores = subspace_scores(reference_basis, candidate_h, wc)
        indices = sorted(np.argsort(scores)[-2:].tolist(), key=lambda i: freq[i, 0])
        if min(scores[indices]) < .8:
            raise ValueError(f"p-pair identity unresolved: projection scores={scores}")
        if min(scores[indices])-max(np.delete(scores, indices)) < .3:
            raise ValueError("p subspace ambiguous among candidate modes")
    scores = subspace_scores(reference_basis, candidate_h, wc)
    s4.write_csv(logs / "eigenfrequencies.csv", [
        {"mode_idx": i, "re_THz": f[0], "im_THz": f[1], "q": f[2],
         "reference_p_projection": scores[i], "selected": i in indices}
        for i, f in enumerate(freq)])
    omegas = 2*np.pi*1e12*(freq[indices, 0]-1j*freq[indices, 1])
    axis = np.linspace(-.3*float(model.java.param().evaluate("a")),
                        .3*float(model.java.param().evaluate("a")), 17)
    xx, yy = np.meshgrid(axis, axis)
    qt = np.column_stack([xx.ravel(), yy.ravel()])
    training_fields = [dec.field_pair(sampler, indices, qt*mirror, 0., ["ewfd.Hz"])[:, 0]
                       for mirror in ([1, 1], [-1, 1], [1, -1])]
    coefficients = dec.symmetry_coefficients(*training_fields, np.ones(len(qt)))
    validation = [candidate_h[:, indices]] + [
        dec.field_pair(sampler, indices, qc*mirror, 0., ["ewfd.Hz"])[:, 0]
        for mirror in ([-1, 1], [1, -1])]
    anchor = -qc[:, 0] if reference is None else reference["combined"]
    coefficients *= s4.phase_and_norm(validation[0]@coefficients, wc, anchor)
    parity = dec.weighted_parity(*(f@coefficients for f in validation), wc)
    if max(parity.values()) >= data.manifest["s4"]["parity_tolerance"]:
        raise ValueError(f"Independent symmetry check failed: {parity}")
    combined_coarse = validation[0]@coefficients
    overlap = 1. if reference is None else abs(wc@(np.conj(reference["combined"])*combined_coarse))/np.sqrt(
        (wc@abs(reference["combined"])**2)*(wc@abs(combined_coarse)**2))
    if overlap < .9:
        raise ValueError(f"Cross-mesh py field overlap is too small: {overlap}")
    eps = [float(data.geometry["simulation_config"][key])**2
           for key in ("air_refractive_index", "slab_refractive_index")]
    fields = dec.field_pair(sampler, indices, qf, 0.,
        ["ewfd.Ex", "ewfd.Ey", "ewfd.Hz", "d(laginterp(2,ewfd.Hz),x)/(1[A/m^2])"])
    epsilon = np.where(air, eps[0], eps[1])
    combined = np.einsum("ncm,m->nc", fields, coefficients)
    gradients = fields[:, 3, :].T[:, None, :]
    reconstructed = convert_derivatives(gradients, omegas, epsilon)[:, 0, :]
    ey_from_h = coefficients @ reconstructed
    area_ey, gradient_ey = wf@combined[:, 1]/area, wf@ey_from_h/area
    local_rows = []
    interface = distance(qf, holes) <= 5e-9
    masks = {"all": np.ones(len(qf), bool), "air_bulk": air & ~interface,
             "dielectric_bulk": ~air & ~interface,
             "interface_air": air & interface, "interface_dielectric": ~air & interface}
    for mode, ey, eh in [("py", combined[:, 1], ey_from_h)] + [
        (f"raw_mode{indices[j]}", fields[:, 1, j], reconstructed[j]) for j in range(2)]:
        for region, mask in masks.items():
            local_rows.append({"mode": mode, "region": region,
                               **metrics(ey, eh, wf, mask, area)})
    s4.write_csv(logs / "local_residual.csv", local_rows)
    np.savez_compressed(logs / "fields_z0.npz", xy_m=qf, weights_m2=wf, raw_fields=fields,
                        coefficients=coefficients, omegas=omegas, mode_indices=indices,
                        combined_fields=combined[:, :3], reconstructed_Ey=ey_from_h,
                        air=air, outer_m=outer, holes_m=holes)
    # Reuse the exact existing edge/area implementation and its two quadrature levels.
    settings = {**data.manifest["s4"], "z_nm": [0.]}
    balances = dec.python_integrals(sampler, indices, coefficients, omegas, outer, holes,
                                   settings, eps, logs)
    if not np.isclose(balances[-1]["Ey_area"], area_ey, rtol=1e-9, atol=1e-10):
        raise ValueError("Independent field and integral extraction disagree")
    # Native checks retain all three edge traces; no trace chosen to force agreement.
    dec.native_checks(model, sampler, indices, coefficients, omegas, outer, holes, eps, logs)
    native_g = native_gradient(model, indices, coefficients, omegas, eps, area, logs)
    native = pd.read_csv(logs / "native_balance.csv")
    nr = native[(native.order == 12) & (native.trace == "inherited")].iloc[0]
    cv = lambda row, key: complex(row[key+"_re"], row[key+"_im"])
    row = {"mesh": level, "frequency_p_low_THz": freq[indices[0], 0],
           "frequency_p_high_THz": freq[indices[1], 0],
           "p_split_GHz": abs(freq[indices[0], 0]-freq[indices[1], 0])*1000,
           "field_overlap_baseline": float(overlap), **parity,
           "Ey_area": area_ey, "Ey_gradient": gradient_ey,
           "Ey_boundary": balances[-1]["Ly"],
           "Ey_area_native": cv(nr, "Ey_area"), "Ey_gradient_native": native_g,
           "Ey_boundary_native": cv(nr, "Ly"),
           "Ey_thickness_native": cv(nr, "Zy"), "Ey_outer_native": cv(nr, "Oy"),
           "local_relative_L2": local_rows[0]["relative_L2"],
           "area_quadrature_relative": complex_change(balances[0]["Ey_area"], area_ey)[1],
           "boundary_quadrature_relative": complex_change(balances[0]["Ly"], balances[-1]["Ly"])[1],
           "gradient_python_native_relative": complex_change(gradient_ey, native_g)[1],
           "Hz_rms_fine": float(np.sqrt(wf@abs(combined[:, 2])**2/area)),
           "coefficients": [s4.complex_columns({"mode_idx": int(i), "value": v}) for i, v in zip(indices, coefficients)]}
    s4.write_json(logs / "summary.json", s4.complex_columns(row))
    return row, {"basis": reference_basis, "combined": combined_coarse} if reference is None else reference


def run(source, output, levels):
    from scripts.analysis.plot_boundary_control import load_control
    output = Path(output).resolve()
    if s4.SESSION_ROOT != output or output.parent != (s4.ROOT / "scripts/.out/unit_cell_2D").resolve():
        raise ValueError("Require a dedicated direct-child S4 session directory")
    execution = output / "99_config/execution.json"
    if execution.exists():
        raise FileExistsError("Execution already attempted; do not overwrite results")
    if levels != [4, 3, 2]:
        raise ValueError("Approved scope is exactly three new meshes: 4, 3, 2")
    data = load_control(source)
    data.geometry = s4.json_read(data.case / "99_config/realized_geometry.json")
    c, = data.manifest["cases"]
    if (c["b0_nm"], c["eta"], c["zeta"], c["mesh"]) != (245., .95, 1., 5):
        raise ValueError("Unexpected baseline geometry or mesh")
    source_mph = data.case / "00_model/s4_gamma.mph"
    before = {str(data.case / p): sha for p, sha in data.audit["source_sha256"].items()}
    before[str(data.source / "audit.json")] = s4.digest(data.source / "audit.json")
    before[str(s4.ROOT / "scripts/parameter.json")] = s4.digest(s4.ROOT / "scripts/parameter.json")
    for p, sha in before.items():
        if s4.digest(p) != sha:
            raise ValueError("Baseline source hash mismatch")
    record = {"status": "preflight", "pid": os.getpid(), "levels": [5, *levels],
              "source_case": str(data.case), "source_sha256": before,
              "source_manifest": data.manifest, "source_geometry": data.geometry,
              "new_solves": 0, "completed_levels": [], "script_sha256": s4.digest(__file__),
              "user_mode": "py", "internal_component": "px",
              "preflight": s4.runtime_preflight()}
    s4.write_json(execution, record)
    client = model = None
    reference, rows, previous_stats = None, [], None
    try:
        import mph
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("COMSOL/Wave Optics license checkout failed")
        for level in [5, *levels]:
            case = output / f"mesh{level}"
            if case.exists():
                raise FileExistsError(f"Refusing to overwrite {case}")
            (case / "00_model").mkdir(parents=True)
            (case / "99_config").mkdir()
            snapshot = snapshot_for_mesh(s4.json_read(data.case / "99_config/parameter.json"), level)
            s4.write_json(case / "99_config/parameter.json", snapshot)
            s4.load_shared_parameters(case / "99_config/parameter.json")
            record.update(status="loading", current_mesh=level)
            s4.write_json(execution, record)
            print(f"Mesh {level}: loading the unchanged source MPH", flush=True)
            model = client.load(str(source_mph))
            params = model.java.param()
            if any(abs(float(params.evaluate(k))) > 1e-12 for k in ("kx", "ky")):
                raise ValueError("Source is not Gamma")
            if not np.allclose([float(params.evaluate("a")), float(params.evaluate("H"))],
                               [820e-9, 200e-9], rtol=1e-12, atol=0):
                raise ValueError("Geometry identity mismatch")
            mesh = model.java.component("comp1").mesh("mesh1")
            if int(mesh.autoMeshSize()) != 5 or not mesh.isAutomatic():
                raise ValueError("Expected an unchanged physics-controlled mesh-5 source")
            client.java.showProgress(str(case / "00_model/comsol_progress.log"))
            start = time.monotonic()
            if level != 5:
                record["status"] = "meshing"
                s4.write_json(execution, record)
                print(f"Mesh {level}: remeshing, same geometry and physics", flush=True)
                model.java.sol("sol1").clearSolutionData()
                mesh.autoMeshSize(level)
                mesh.run()
            stats = mesh_stats(model)
            if previous_stats and stats["elements"].get("tet", 0) <= previous_stats["elements"].get("tet", 0):
                raise ValueError("Nominal refinement did not increase the actual tetrahedron count")
            stats["mesh_seconds"] = time.monotonic()-start
            s4.write_json(case / "99_config/mesh.json", stats)
            print(f"Mesh {level} built: {stats}", flush=True)
            if level != 5:
                record["status"] = "solving"
                s4.write_json(execution, record)
                start = time.monotonic()
                print(f"Mesh {level}: Gamma eigensolve", flush=True)
                model.java.sol("sol1").runAll()
                stats["solve_seconds"] = time.monotonic()-start
                record["new_solves"] += 1
                model.save(str(case / "00_model/s4_gamma.mph"))
                stats["mph_sha256"] = s4.digest(case / "00_model/s4_gamma.mph")
                print(f"Mesh {level}: solved and new checkpoint saved", flush=True)
            else:
                stats.update(solve_seconds=0., source_mph=str(source_mph))
            s4.write_json(case / "99_config/mesh.json", stats)
            record["status"] = "analyzing"
            s4.write_json(execution, record)
            row, reference = analyze(model, data, level, case, reference)
            row.update(tetrahedra=stats["elements"]["tet"], vertices=stats["vertices"],
                       slab_tetrahedra=stats["slab_layer_tetrahedra"])
            rows.append(row)
            s4.write_csv(output / "80_logs/convergence.csv", [{k: v for k, v in r.items() if k != "coefficients"} for r in rows])
            record["completed_levels"].append(level)
            s4.write_json(execution, record)
            previous_stats = stats
            client.remove(model)
            model = None
            print(f"Mesh {level}: extraction complete", flush=True)
        record["status"] = "complete_pending_scientific_review"
    except Exception as exc:
        record.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all(s4.digest(p) == h for p, h in before.items())
        s4.write_json(execution, record)
    if not record["source_unchanged"]:
        raise RuntimeError("Baseline or shared parameters changed")
    report(output)


def reprocess(source, output):
    """Re-evaluate existing checkpoints after correcting only exported frequency units."""
    from scripts.analysis.plot_boundary_control import load_control
    output = Path(output).resolve()
    if output != s4.SESSION_ROOT:
        raise ValueError("Session root mismatch")
    original = s4.json_read(output / "99_config/execution.json")
    if original["completed_levels"] != [5, 4, 3, 2] or not original["source_unchanged"]:
        raise ValueError("Require a complete, unchanged mesh series")
    target = output / "corrected_analysis"
    target.mkdir(exist_ok=False)
    data = load_control(source)
    data.geometry = s4.json_read(data.case / "99_config/realized_geometry.json")
    source_paths = [data.case / "00_model/s4_gamma.mph"] + [output / f"mesh{v}/00_model/s4_gamma.mph" for v in (4, 3, 2)]
    before = {str(p): s4.digest(p) for p in source_paths}
    record = {"status": "starting", "eigensolves": 0, "source_sha256": before,
              "reason": "Explicit THz scaling for damping; previous new-mesh integrals invalid",
              "script_sha256": s4.digest(__file__), "completed_levels": [],
              "preflight": s4.runtime_preflight()}
    s4.write_json(output / "99_config/reanalysis.json", record)
    client = model = None
    reference, rows = None, []
    try:
        import mph
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("License unavailable")
        for level, source_mph in zip((5, 4, 3, 2), source_paths):
            print(f"Mesh {level}: explicit-unit reanalysis, no new solve", flush=True)
            model = client.load(str(source_mph))
            raw_units = [str(u) for u in model.java.result().numerical("gev1").getStringArray("unit")]
            stats = mesh_stats(model)
            case = target / f"mesh{level}"
            row, reference = analyze(model, data, level, case, reference)
            row.update(tetrahedra=stats["elements"]["tet"], vertices=stats["vertices"],
                       slab_tetrahedra=stats["slab_layer_tetrahedra"])
            rows.append(row)
            s4.write_json(case / "80_logs/unit_audit.json", {"original_gev1_units": raw_units,
                "frequency_method": "omega/(2*pi*1[THz]), damp/(2*pi*1[THz]); unit=1",
                "source_mph": str(source_mph), "source_sha256": before[str(source_mph)]})
            s4.write_csv(target / "80_logs/convergence.csv", [{k:v for k,v in r.items() if k != "coefficients"} for r in rows])
            client.remove(model)
            model = None
            record["completed_levels"].append(level)
            s4.write_json(output / "99_config/reanalysis.json", record)
        record["status"] = "complete_pending_scientific_review"
    except Exception as exc:
        record.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all(s4.digest(p) == h for p,h in before.items())
        s4.write_json(output / "99_config/reanalysis.json", record)
    report(target)
    original.update(active_analysis="corrected_analysis", original_analysis_valid=False,
                    original_analysis_issue="COMSOL changed the damping display unit after remeshing; see reanalysis.json")
    s4.write_json(output / "99_config/execution.json", original)
    (output / "80_logs/READ_FIRST.md").write_text(
        "# Corrected analysis\n\nUse ../corrected_analysis/80_logs/report.md and convergence.csv.\n"
        "Original mesh4/3/2 analysis used a damping column in Hz as THz and is invalid.\n"
        "The solved MPH checkpoints are unchanged; corrected analysis explicitly scales both frequency components.\n", encoding="utf-8")


def scale_mesh_sizes(mesh, factor):
    """Preserve mesh topology/periodic pairs while scaling both global and local sizes."""
    if not 0 < factor < 1:
        raise ValueError("Refinement factor must be between zero and one")
    mesh.automatic(False)
    changes = []

    def walk(parent, prefix):
        for tag in parent.feature().tags():
            feature = parent.feature(str(tag))
            if str(feature.getType()) in ("Size", "MeshSizeDefault"):
                feature.set("custom", "on")
                for key in ("hmax", "hmin"):
                    old = str(feature.getString(key))
                    new = f"({old})*{factor:.12g}"
                    feature.set(key, new)
                    changes.append({"feature": "/".join(prefix+[str(tag)]), "property": key,
                                    "old": old, "new": str(feature.getString(key))})
            walk(feature, prefix+[str(tag)])

    walk(mesh, [])
    if len(changes) < 4:
        raise ValueError("Expected global and dielectric size constraints")
    return changes


def run_sizes(source, output):
    """True h-refinement; automatic grades alone leave this slab's size fixed."""
    from scripts.analysis.plot_boundary_control import load_control
    output = Path(output).resolve()
    if s4.SESSION_ROOT != output:
        raise ValueError("Session root mismatch")
    target = output / "size_refinement"
    target.mkdir(exist_ok=False)
    data = load_control(source)
    data.geometry = s4.json_read(data.case / "99_config/realized_geometry.json")
    source_mph = data.case / "00_model/s4_gamma.mph"
    before = {str(data.case / p): sha for p,sha in data.audit["source_sha256"].items()}
    before[str(s4.ROOT / "scripts/parameter.json")] = s4.digest(s4.ROOT / "scripts/parameter.json")
    if not all(s4.digest(p) == h for p,h in before.items()):
        raise ValueError("Source identity mismatch")
    record = {"status": "starting", "source_sha256": before, "factors": [1., .5, .25, .125],
              "new_solves": 0, "completed_factors": [], "preflight": s4.runtime_preflight(),
              "script_sha256": s4.digest(__file__), "frequency_units": "explicit THz and Q consistency check"}
    status_path = target / "99_config/execution.json"
    s4.write_json(status_path, record)
    rows, reference, previous = [], None, None
    client = model = None
    try:
        import mph
        client = mph.start(version="6.3")
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL", "WAVEOPTICS")):
            raise RuntimeError("License unavailable")
        for factor in record["factors"]:
            label = f"h{factor:g}"
            case = target / label
            (case / "00_model").mkdir(parents=True)
            model = client.load(str(source_mph))
            params = model.java.param()
            if any(abs(float(params.evaluate(k))) > 1e-12 for k in ("kx", "ky")):
                raise ValueError("Source must be at Gamma")
            client.java.showProgress(str(case / "00_model/comsol_progress.log"))
            record.update(current_factor=factor, status="loading")
            s4.write_json(status_path, record)
            mesh = model.java.component("comp1").mesh("mesh1")
            changes = []
            if factor != 1:
                record["status"] = "meshing"
                s4.write_json(status_path, record)
                model.java.sol("sol1").clearSolutionData()
                changes = scale_mesh_sizes(mesh, factor)
                print(f"h/h0={factor}: actual size changes {changes}", flush=True)
                mesh.clearMesh()
                mesh.run()
            stats = mesh_stats(model)
            stats.update(h_over_h0=factor, changed_constraints=changes)
            s4.write_json(case / "99_config/mesh.json", stats)
            print(f"h/h0={factor}: built {stats['elements']['tet']} tetrahedra; slab={stats['slab_layer_tetrahedra']}", flush=True)
            if previous and stats["slab_layer_tetrahedra"] <= previous["slab_layer_tetrahedra"]*1.25:
                raise ValueError("Slab mesh did not actually refine")
            if factor != 1:
                record["status"] = "solving"
                s4.write_json(status_path, record)
                start = time.monotonic()
                model.java.sol("sol1").runAll()
                stats["solve_seconds"] = time.monotonic()-start
                model.save(str(case / "00_model/s4_gamma.mph"))
                stats["mph_sha256"] = s4.digest(case / "00_model/s4_gamma.mph")
                record["new_solves"] += 1
                s4.write_json(case / "99_config/mesh.json", stats)
                print(f"h/h0={factor}: solved, checkpoint saved", flush=True)
            record["status"] = "analyzing"
            s4.write_json(status_path, record)
            row, reference = analyze(model, data, factor, case, reference)
            row.update(tetrahedra=stats["elements"]["tet"], vertices=stats["vertices"],
                       slab_tetrahedra=stats["slab_layer_tetrahedra"], mesh_control="size_scale")
            rows.append(row)
            s4.write_csv(target / "80_logs/convergence.csv", [{k:v for k,v in r.items() if k != "coefficients"} for r in rows])
            record["completed_factors"].append(factor)
            s4.write_json(status_path, record)
            previous = stats
            client.remove(model)
            model = None
            print(f"h/h0={factor}: analysis complete", flush=True)
        record["status"] = "complete_pending_scientific_review"
    except Exception as exc:
        record.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        record["source_unchanged"] = all(s4.digest(p) == h for p,h in before.items())
        s4.write_json(status_path, record)
    report(target)
    original = s4.json_read(output / "99_config/execution.json")
    original.update(active_analysis="size_refinement", original_analysis_valid=False,
                    original_analysis_issue="Unchanged slab mesh and damping display-unit error; use size_refinement")
    s4.write_json(output / "99_config/execution.json", original)


def report(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output = Path(output)
    df = pd.read_csv(output / "80_logs/convergence.csv")
    values = {key: df[key+"_re"].to_numpy()+1j*df[key+"_im"].to_numpy() for key in
              ("Ey_area", "Ey_gradient", "Ey_boundary", "Ey_area_native", "Ey_gradient_native", "Ey_boundary_native")}
    phase = np.exp(-1j*np.angle(values["Ey_area_native"][0]))
    changes = []
    for i in range(1, len(df)):
        for key, val in values.items():
            absolute, relative = complex_change(val[i-1], val[i])
            changes.append({"coarse_mesh": float(df.mesh.iloc[i-1]), "fine_mesh": float(df.mesh.iloc[i]),
                            "quantity": key, "absolute_change_Vm": absolute, "relative_change": relative})
    s4.write_csv(output / "80_logs/adjacent_changes.csv", changes)
    size_control = "mesh_control" in df and (df.mesh_control == "size_scale").all()
    label = "h/h0" if size_control else "mesh"
    lines = ["# S4 Γ 点网格收敛检查", "", "几何和三维模型不变，仅加密网格。" +
             ("h/h0 为相对原始最大/最小单元尺寸的缩放比例，不是自动网格等级。" if size_control else "自动等级试算。"),
             "各档采用相同 Hz RMS 归一化，相位与 mesh 5 的 py 场对齐，不分别归一化电场积分。",
             "主要网格比较使用原生 12 阶积分，并用 8 阶检验积分阶数稳定性。三种边界迹均保留，不挑选最接近 Ey 的结果。",
             "三种积分保留完整复数；下表共同相位使基准的原生直接 Ey 平均为正实数。", "",
             f"| {label} | 四面体数 | 原生直接 Ey (V/m) | 原生 Hz 导数 (V/m) | 原生孔边界和 (V/m) | 局部 L2 残差 |",
             "|---:|---:|---:|---:|---:|---:|"]
    for i, r in df.iterrows():
        lines.append(f"| {r.mesh:g} | {int(r.tetrahedra)} | {phase*values['Ey_area_native'][i]:.7g} | {phase*values['Ey_gradient_native'][i]:.7g} | {phase*values['Ey_boundary_native'][i]:.7g} | {r.local_relative_L2:.6g} |")
    if size_control:
        checks = []
        for i, row in df.iterrows():
            logs = output / f"h{row.mesh:g}" / "80_logs"
            balance = pd.read_csv(logs / "native_balance.csv")
            gradient = pd.read_csv(logs / "native_gradient.csv")
            def z(record, key):
                return complex(record[key+"_re"], record[key+"_im"])
            b8 = balance[(balance.order == 8) & (balance.trace == "inherited")].iloc[0]
            b12 = balance[(balance.order == 12) & (balance.trace == "inherited")].iloc[0]
            g8 = z(gradient[gradient.order == 8].iloc[0], "Ey_gradient")
            g12 = z(gradient[gradient.order == 12].iloc[0], "Ey_gradient")
            air = z(balance[(balance.order == 12) & (balance.trace == "air")].iloc[0], "Ly")
            dielectric = z(balance[(balance.order == 12) & (balance.trace == "dielectric")].iloc[0], "Ly")
            boundary = z(b12, "Ly")
            ey = z(b12, "Ey_area")
            order_change = max(complex_change(z(b8,k),z(b12,k))[1] for k in ("Ey_area", "Ly"))
            checks.append({"h_over_h0": float(row.mesh),
                           "native_8_to_12_max_relative": max(order_change, complex_change(g8,g12)[1]),
                           "air_dielectric_trace_spread_over_boundary": abs(air-dielectric)/max(abs(boundary),1e-30),
                           "gradient_boundary_difference_over_gradient": abs(g12-boundary)/max(abs(g12),1e-30),
                           "Ey_gradient_difference_over_Ey": abs(ey-g12)/max(abs(ey),1e-30)})
        s4.write_csv(output / "80_logs/integration_checks.csv", checks)
        lines += ["", "## 原生求积与边界迹检查", "",
                  "表中均为无量纲比值；边界迹差以继承迹边界和为尺度，导数–边界差以导数积分为尺度，电场–导数差以直接电场积分为尺度。",
                  "积分阶数稳定仅检验当前有限元场的求积，不替代网格收敛。", "",
                  "| h/h0 | 8→12 阶最大相对变化 | 空气/介质边界迹差 | 导数–边界差 | 电场–导数差 |",
                  "|---:|---:|---:|---:|---:|"]
        for check in checks:
            lines.append("| " + " | ".join(f"{v:.7g}" for v in check.values()) + " |")
    lines += ["", "## 固定节点的 Python 求积对照", "",
              "该求积不贴合有限元内部单元边界，不能只凭节点总数判断误差可忽略。", "",
              f"| {label} | Ey 与原生积分差 | Hz 导数与原生积分差 | 边界和与原生积分差 |",
              "|---:|---:|---:|---:|"]
    for i,r in df.iterrows():
        differences = [complex_change(values[k][i],values[k+"_native"][i])[1]
                       for k in ("Ey_area", "Ey_gradient", "Ey_boundary")]
        lines.append(f"| {r.mesh:g} | {differences[0]:.7g} | {differences[1]:.7g} | {differences[2]:.7g} |")
    lines += ["", "## 最细两档的完整复差", "",
              "| 量 | 绝对变化 (V/m) | 相对变化 |", "|---|---:|---:|"]
    for r in changes:
        if r["fine_mesh"] == float(df.mesh.iloc[-1]):
            lines.append(f"| {r['quantity']} | {r['absolute_change_Vm']:.7g} | {r['relative_change']:.7g} |")
    frequency_change = max(abs(df[k].iloc[-1]-df[k].iloc[-2])/abs(df[k].iloc[-1])
                           for k in ("frequency_p_low_THz", "frequency_p_high_THz"))
    final_changes = [r["relative_change"] for r in changes if r["fine_mesh"] == float(df.mesh.iloc[-1])]
    integral_stable = max(final_changes) < .01
    native_changes = [r["relative_change"] for r in changes if r["fine_mesh"] == float(df.mesh.iloc[-1]) and r["quantity"].endswith("_native")]
    lines += ["", f"两个 p 模末两档的最大相对频率变化：{frequency_change:.7g}。",
              f"频率变化满足预设 1e-4 判据：{frequency_change < 1e-4}。",
              f"三项原生复积分末两档变化均小于 1%：{max(native_changes)<.01}。",
              f"六项 Python/原生复积分末两档变化均小于 1%：{integral_stable}。",
              "是否已收敛还须核对求积变化与逐边迹敏感性；两档相近不是严格误差界。",
              "各 mesh 的 native_balance.csv 保留继承、空气、介质三种边界迹；未选择最接近电场的那一种作为真值。",
              "双模组合不是新的精确单频本征态；各模态频率在 Maxwell 换算中分开使用。",
              "本轮未改变三维近似，也未通过调整幅值或相位强制两个表达相等。", ""]
    gap = abs(values["Ey_area_native"][-1]-values["Ey_gradient_native"][-1])/abs(values["Ey_area_native"][-1])
    magnetic_gap = abs(values["Ey_gradient_native"][-1]-values["Ey_boundary_native"][-1])/abs(values["Ey_gradient_native"][-1])
    lines += ["## 本轮判断", "",
              f"最细网格的电场–导数积分相对复差为 {gap:.4%}；导数–孔边界和的相对差为 {magnetic_gap:.4%}。",
              "网格加密改善了两种磁场积分的一致性，但未消除它们与直接电场积分的大幅差异。",
              "该结果不支持将大幅差异主要归因于粗网格；它也不能单独判定差异的物理来源。",
              "本轮未全部达到预设收敛判据，不将结果表述为所有量均已严格收敛，也不自动追加更细网格或修改 SI。", ""]
    (output / "80_logs/report.md").write_text("\n".join(lines), encoding="utf-8")
    plt.rcParams.update({"font.size": 10, "xtick.direction": "in", "ytick.direction": "in"})
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")
    x = df.tetrahedra.to_numpy()
    for key, label, color in (("Ey_area_native", "Direct $E_y$", "black"), ("Ey_gradient_native", "$H_z$ derivative", "#4477aa"),
                              ("Ey_boundary_native", "Hole boundary sum", "#ee7733")):
        v = phase*values[key]
        axes[0].plot(x, v.real, "o-", label=label, color=color)
        axes[1].plot(x, v.imag, "o-", label=label, color=color)
    axes[0].set(ylabel="Re mean field (V/m)")
    axes[1].set(ylabel="Im mean field (V/m)")
    axes[2].plot(x, df.local_relative_L2, "o-", color="#228833")
    axes[2].set(ylabel="Relative L2 residual: derivative vs $E_y$")
    for ax in axes:
        ax.set(xscale="log", xlabel="Number of tetrahedra")
        ax.tick_params(which="both", direction="in")
    axes[0].legend(fontsize=8)
    fig.suptitle("$p_y$ at Γ: fixed geometry; native integrals for mesh comparison")
    (output / "10_overview").mkdir(exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(output / f"10_overview/mesh_convergence_py.{ext}", dpi=170)
    plt.close(fig)
    print(f"Convergence report: {output / '80_logs/report.md'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_session", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--reprocess", action="store_true")
    parser.add_argument("--size-refinement", action="store_true")
    args = parser.parse_args()
    if args.size_refinement:
        run_sizes(args.source_session, args.output)
    elif args.reprocess:
        reprocess(args.source_session, args.output)
    elif args.report_only:
        report(args.output)
    else:
        run(args.source_session, args.output, [4, 3, 2])


if __name__ == "__main__":
    main()
