"""S4 preparation and single-Gamma-point COMSOL execution.

Use from repository root: python -m scripts.run_main.run_boundary_analysis --help.
Preparation never imports mph, builds a mesh, or changes shared parameters.
"""
from __future__ import annotations

from comsol_workflow.result_io import digest, read_json as json_read, write_json, write_csv

import argparse
import ast
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
from types import SimpleNamespace
import uuid

import numpy as np

from comsol_workflow.output_paths import data_directory

from comsol_workflow.boundary_integrals import area_rule, complex_columns, edge_integrals, group_totals, hole_edges, line_rule, maxwell_terms, parity_errors, phase_and_norm, polygon_area, polygon_edges, length_only_prediction, segment_minimum
from scripts.run_main.parameter_config import (
    load_shared_parameters, parameter_path_from_environment, UNIT_CELL_2D_OUTPUT_ROOT,
)

ROOT = Path(__file__).resolve().parents[2]
VERSION = "s4-boundary-v3"


def configured_session_root():
    """Opt in to one task container; never redirect to an arbitrary directory."""
    value = os.environ.get("COMSOL_WORKFLOW_S4_SESSION_ROOT")
    if not value:
        return None
    directory = Path(value).resolve()
    if directory.parent != UNIT_CELL_2D_OUTPUT_ROOT.resolve():
        raise ValueError("S4 session root must be a direct child of scripts/.out/unit_cell_2D")
    return directory


SESSION_ROOT = configured_session_root()
TASK_ROOT = SESSION_ROOT / "tmp" if SESSION_ROOT is not None else ROOT / "tmp"










def source_constants():
    """Read existing literal defaults without importing the solver module."""
    wanted = {"A", "D", "SLAB_HEIGHT", "REFRACTIVE_INDEX", "MODE_TYPE", "WAVELENGTH"}
    source = ROOT / "scripts/run_main/run_band_pair.py"
    result = {}
    for node in ast.parse(source.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in wanted:
                    result[target.id] = ast.literal_eval(node.value)
    if result.keys() != wanted:
        raise ValueError("Band constants changed: inspect the existing builder before execution")
    return result


def source_hashes():
    paths = [Path(__file__), ROOT / "comsol_workflow/boundary_integrals.py",
             ROOT / "comsol_workflow/boundary_plotting.py",
             ROOT / "comsol_workflow/output_paths.py",
             ROOT / "comsol_workflow/result_io.py",
             ROOT / "scripts/run_main/run_band_pair.py",
             ROOT / "scripts/run_main/parameter_config.py",
             ROOT / "comsol_workflow/geometry_utils.py",
             ROOT / "comsol_workflow/basis_utils.py",
             ROOT / "comsol_workflow/simulation_spatial/hexagon_unit_cell.py"]
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in paths}


def prepare(args):
    source = Path(args.parameter_file or parameter_path_from_environment()).resolve()
    parameters = load_shared_parameters(source)
    settings = parameters.unit_cell_2d
    run_id = args.run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:6]
    if not re.fullmatch(r"[a-zA-Z0-9_-]{1,96}", run_id):
        raise ValueError("Invalid run ID")
    directory = TASK_ROOT / ("s4_" + run_id)
    if directory.exists():
        raise FileExistsError(f"Refusing to reuse output directory: {directory}")
    eta = settings.cell.eta if args.eta is None else args.eta
    zetas = [1.156] if args.zeta is None else args.zeta
    if len(set(zetas)) != len(zetas) or not zetas:
        raise ValueError("Use distinct zeta values")
    if len(zetas) > args.max_cases:
        raise ValueError("Case count exceeds explicit --max-cases budget")
    meshes = args.mesh if args.mesh is not None else [parameters.mesh_size]
    if len(set(meshes)) != len(meshes) or len(meshes) * len(zetas) > args.max_cases:
        raise ValueError("Mesh/case count exceeds the case budget")
    constants = source_constants()
    cases, snapshots = [], []
    for mesh in meshes:
        for zeta in zetas:
            snapshot = json_read(source)
            snapshot["mesh_size"] = mesh
            snapshot["unit_cell_2d"].update(eta=eta, zeta=zeta)
            # Snapshot validation uses the shared loader, not an independent schema.
            case_id = f"case{len(cases):03d}_eta{eta:g}_zeta{zeta:g}_mesh{mesh}"
            cases.append({"case_id": case_id, "b0_nm": settings.cell.b0_nm,
                          "eta": eta, "zeta": zeta, "mesh": mesh,
                          "snapshot": f"{case_id}/99_config/parameter.json"})
            snapshots.append(snapshot)
    directory.mkdir(parents=True)
    (directory / "99_config").mkdir()
    (directory / "99_config/parameter_source.json").write_bytes(source.read_bytes())
    for case, snapshot in zip(cases, snapshots):
        path = directory / case["snapshot"]
        write_json(path, snapshot)
        load_shared_parameters(path)
        case["snapshot_sha256"] = digest(path)
    manifest = {
        "version": VERSION, "run_id": run_id, "status": "prepared_not_solved",
        "source_parameter": str(source), "source_parameter_sha256": digest(source),
        "source_files": source_hashes(), "constants": constants, "cases": cases,
        "output_directory": str(result_directory(run_id, cases)),
        "source_zeta": settings.cell.zeta, "configured_target_bands": list(settings.target_bands),
        "execution": {"k_over_G": [0, 0], "solve_count": len(cases),
                      "eigenmode_pair_count": settings.eigenmode_pair_count,
                      "expected_eigenvalue_count": 2 * settings.eigenmode_pair_count,
                      "shift": "c_const/1.55[um]", "full_k_scan": False},
        "s4": {"edge_order": args.edge_order, "area_order": args.area_order,
               "area_subdivisions": args.area_subdivisions,
               "z_nm": args.z_nm, "primary_z_nm": 0.0,
               "surface_side": "slab_layer_domains_at_exact_z",
               "trace_offset_over_a": args.trace_offset_over_a,
               "parity_tolerance": 0.05, "identity_tolerance": 1e-3,
               "p_weight_minimum": 0.5, "mode_frequency_tolerance_thz": settings.frequency_tolerance_thz,
               "minimum_mode_overlap": settings.minimum_overlap,
               "minimum_ambiguity_gap": settings.minimum_ambiguity_gap,
               "coordinate_frame": "existing COMSOL x/y; no axis swap",
               "target": "manuscript py: Hz odd under x reflection, even under y reflection",
               "internal_component_expected": "px; verify rather than relabel",
               "phasor_export": "exp(-i*omega*t); complex E,H,frequency conjugated from COMSOL",
               "normal": "dielectric_to_air", "normalization": "cell_Hz_RMS_at_z0"},
        "limitations": ["3D half-slab model, not exact 2D TE",
                        "No automatic degenerate-subspace rotation",
                        "No cache is considered a Gamma solution without recorded k identity",
                        "No automatic root search or full parameter-space scan"],
    }
    write_json(directory / "99_config/s4_manifest.json", manifest)
    write_json(directory / "99_config/prepared_identity.json",
               {"manifest_sha256": digest(directory / "99_config/s4_manifest.json")})
    assert "mph" not in sys.modules, "prepare imported COMSOL"
    print(json.dumps({"prepared": str(directory), "cases": cases,
                      "execution": manifest["execution"], "constants": constants,
                      "output": manifest["output_directory"], "s4": manifest["s4"],
                      "configured_target_bands": manifest["configured_target_bands"],
                      "source_zeta": settings.cell.zeta}, indent=2))
    return directory


def result_directory(run_id, cases):
    first = cases[0]
    label = f"b{first['b0_nm']:g}_eta{first['eta']:g}_zeta{first['zeta']:g}_mesh{first['mesh']}"
    output_root = SESSION_ROOT / "results" if SESSION_ROOT is not None else UNIT_CELL_2D_OUTPUT_ROOT
    return (output_root / f"S4_{label}_Gamma_{len(cases)}case_{run_id}").resolve()


def prepared_manifest(directory):
    directory = Path(directory).resolve()
    if directory.parent != TASK_ROOT.resolve() or not directory.name.startswith("s4_"):
        raise ValueError(f"Preparation must be a task directory under {TASK_ROOT}/s4_<run_id>")
    path = directory / "99_config/s4_manifest.json"
    if digest(path) != json_read(directory / "99_config/prepared_identity.json")["manifest_sha256"]:
        raise RuntimeError("Manifest changed after preparation")
    manifest = json_read(path)
    if manifest["version"] != VERSION or directory.name != "s4_" + manifest["run_id"]:
        raise ValueError("Preparation identity does not match this entry")
    if manifest["source_files"] != source_hashes():
        raise RuntimeError("Code changed after preparation; prepare a new run")
    for case in manifest["cases"]:
        snapshot = (directory / case["snapshot"]).resolve()
        if not snapshot.is_relative_to(directory):
            raise ValueError("Snapshot must stay inside the preparation directory")
        if digest(snapshot) != case["snapshot_sha256"]:
            raise RuntimeError("Parameter snapshot changed after preparation")
        load_shared_parameters(snapshot)
    expected = result_directory(manifest["run_id"], manifest["cases"])
    if Path(manifest["output_directory"]).resolve() != expected:
        raise ValueError("Formal output identity does not match prepared parameters")
    return manifest


def runtime_preflight(*, allow_concurrent=False):
    if os.name != "nt":
        raise RuntimeError("Solver entry is configured for the Windows COMSOL host")
    bin_dir = Path("C:/Program Files/COMSOL/COMSOL63/Multiphysics/bin/win64")
    if not (bin_dir / "comsol.exe").is_file():
        raise RuntimeError("COMSOL 6.3 executable was not found")
    os.environ["PATH"] = str(bin_dir) + os.pathsep + os.environ.get("PATH", "")
    command = "$ErrorActionPreference='Stop'; Get-CimInstance Win32_Process | Where-Object {$_.Name -match '^(python|comsol|comsolmphserver|comsolserver).*\\.exe$'} | Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Compress"
    proc = subprocess.run(["powershell", "-NoProfile", "-Command", command],
                          capture_output=True, text=True, check=True)
    active = json.loads(proc.stdout) if proc.stdout.strip() else []
    active = [active] if isinstance(active, dict) else active
    if any(p["Name"].lower().startswith("comsol") for p in active) and not allow_concurrent:
        raise RuntimeError("Existing COMSOL processes detected; do not share resources implicitly: " + str(active))
    license_dir = bin_dir.parents[1] / "license"
    license_file = license_dir / "license.dat"
    license_status = "Local license file; product availability requires COMSOL API checkout"
    if re.search(r"^SERVER\s", license_file.read_text(errors="replace"), re.MULTILINE):
        license_check = subprocess.run(
            [str(license_dir / "win64/lmutil.exe"), "lmstat", "-c", str(license_file), "-a"],
            capture_output=True, text=True, timeout=30)
        license_status = license_check.stdout + license_check.stderr
        if license_check.returncode or "license server UP" not in license_status or "LMCOMSOL: UP" not in license_status:
            raise RuntimeError("COMSOL license service is unavailable: " + license_status)
    return {"comsol_executable": shutil.which("comsol"), "processes_before_start": active,
            "concurrent_explicitly_allowed": bool(allow_concurrent),
            "license_status": license_status,
            "license": "file/service checked; product checkout checked before model build"}


class Sampler:
    """Arbitrary 3D points, explicit real/imaginary data, SI result convention."""
    def __init__(self, runner):
        from jpype.types import JArray, JDouble, JInt, JString
        self.JArray, self.JDouble, self.JInt, self.JString = JArray, JDouble, JInt, JString
        self.model = runner.model.java
        self.node = self.model.result().numerical().create("s4interp", "Interp")
        self.node.set("data", "dset1")
        self.node.set("ext", 0.0)
        self.node.set("matherr", "on")
        self.node.set("coorderr", "on")
        self.node.set("recover", "off")
        self.slab_top = float(self.model.param().evaluate("H")) / 2
        self.slab_domains = self.JArray(self.JInt, 1)(
            list(self.model.component("comp1").selection("sel_slab_layer").entities(3)))

    def __call__(self, mode, xyz_m, expressions):
        xyz_m = np.asarray(xyz_m, float)
        result = []
        # COMSOL's automatic display unit can be um even for x,y,z. Divide by
        # an explicit SI unit before extraction so every returned value has the
        # stated SI numerical magnitude, independently of the model display.
        units = {"x": "m", "y": "m", "z": "m", "ewfd.Ex": "V/m", "ewfd.Ey": "V/m",
                 "ewfd.Hz": "A/m", "d(ewfd.Hx,z)": "A/m^2", "d(ewfd.Hy,z)": "A/m^2"}
        evaluated = [f"({e})/(1[{units[e]}])" if e in units else e for e in expressions]
        self.node.set("expr", self.JArray(self.JString, 1)(evaluated))
        self.node.set("solnum", self.JArray(self.JInt, 1)([int(mode) + 1]))
        for chunk in np.array_split(xyz_m, max(1, int(np.ceil(len(xyz_m) / 4096)))):
            # geom1 lengthUnit is um; evaluated field variables use SI units.
            self.node.setInterpolationCoordinates(self.JArray(self.JDouble, 2)((chunk.T / 1e-6).tolist()))
            real = np.asarray(self.node.getData(), float)
            imag = np.asarray(self.node.getImagData(), float) if self.node.isComplex() else np.zeros_like(real)
            data = (real + 1j * imag).reshape(len(expressions), -1).T
            if data.shape != (len(chunk), len(expressions)) or not np.isfinite(data).all():
                raise ValueError("COMSOL interpolation returned missing/nonfinite field values")
            result.append(np.conj(data))
        return np.concatenate(result)

    def plane(self, mode, xy, z, expressions):
        # Exact bottom/top planes use the slab-layer domain trace, including
        # hole air. Do not displace a requested interface plane into a bulk layer.
        if 0 <= z <= self.slab_top + 1e-16:
            self.node.selection().set(self.slab_domains)
        else:
            self.node.selection().all()
        return self(mode, np.column_stack([xy, np.full(len(xy), z)]), expressions)


def select_mode(runner, sampler, modes, a, z, output, settings, previous=None):
    # Same fixed physical points in every case and each mirror. Never compare
    # arbitrary COMSOL cut-plane point order across independent meshes.
    axis = np.linspace(-0.3 * a, 0.3 * a, 17)
    xx, yy = np.meshgrid(axis, axis)
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    rows, fields, eligible = [], {}, []
    for index, m in modes.iterrows():
        if not bool(m["is_valid"]):
            continue
        h = sampler.plane(index, xy, z, ["ewfd.Hz"])[:, 0]
        hx = sampler.plane(index, xy * [-1, 1], z, ["ewfd.Hz"])[:, 0]
        hy = sampler.plane(index, xy * [1, -1], z, ["ewfd.Hz"])[:, 0]
        errors = parity_errors(h, hx, hy)
        overlap = (float(abs(np.vdot(previous["field"], h)) /
                         (np.linalg.norm(previous["field"]) * np.linalg.norm(h)))
                   if previous is not None else None)
        row = {"mode_idx": int(index), "frequency_thz": float(m["re"]),
               "p_weight": float(m["p_weight"]), "px_weight_internal": float(m["px_weight"]),
               "py_weight_internal": float(m["py_weight"]), **errors, "overlap_previous": overlap}
        rows.append(row)
        fields[int(index)] = h
        if (m["p_weight"] >= settings["p_weight_minimum"]
                and max(errors.values()) < settings["parity_tolerance"]):
            if previous is None or (overlap >= settings["minimum_mode_overlap"]
                    and abs(m["re"] - previous["frequency_thz"]) <= settings["mode_frequency_tolerance_thz"]):
                eligible.append(row)
    write_csv(output / "mode_identification.csv", rows)
    if len(eligible) > 1 and previous is not None:
        eligible.sort(key=lambda r: r["overlap_previous"], reverse=True)
        if eligible[0]["overlap_previous"] - eligible[1]["overlap_previous"] >= settings["minimum_ambiguity_gap"]:
            eligible = eligible[:1]
    if len(eligible) != 1:
        raise RuntimeError("Target mode missing, mixed or ambiguous. Inspect mode_identification.csv; no forced relabeling")
    chosen = eligible[0]
    idx = chosen["mode_idx"]
    # Stable spatial reference, independent of a possibly vanishing Ey average.
    factor = phase_and_norm(fields[idx], np.ones(len(xy)),
                            previous["field"] if previous is not None else -xy[:, 0])
    state = {"field": fields[idx] * factor, "frequency_thz": chosen["frequency_thz"]}
    return idx, chosen, state, factor / abs(factor)


def section_coordinates(settings, height):
    values = np.asarray(settings["z_nm"], float)
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
            or values[0] != 0 or np.any(np.diff(values) <= 0)
            or np.any(values > height * 0.5e9 + 1e-6)):
        raise ValueError("Sections must start at z=0 and increase within the closed half-slab")
    return [(float(nm), float(nm * 1e-9), float(nm * 1e-9 / height)) for nm in values]


def export_area_field(output, sampler, mode, outer, xy, weights, values, factor, refinement):
    """Save actual z=0 quadrature fields and a separate display raster, in SI units."""
    from matplotlib.path import Path as PolygonPath
    lower, upper = np.min(outer, axis=0), np.max(outer, axis=0)
    spacing = (upper[0] - lower[0]) / 201
    nx, ny = 201, int(np.ceil((upper[1] - lower[1]) / spacing))
    x = np.linspace(lower[0], upper[0], nx, endpoint=False) + (upper[0] - lower[0]) / (2 * nx)
    y = np.linspace(lower[1], upper[1], ny, endpoint=False) + (upper[1] - lower[1]) / (2 * ny)
    xx, yy = np.meshgrid(x, y)
    grid_xy = np.column_stack([xx.ravel(), yy.ravel()])
    inside = PolygonPath(outer).contains_points(grid_xy)
    grid_e = np.full((len(grid_xy), 2), np.nan + 1j * np.nan, complex)
    grid_e[inside] = factor * sampler.plane(mode, grid_xy[inside], 0.0, ["ewfd.Ex", "ewfd.Ey"])
    path = Path(output) / "area_field_z0.npz"
    np.savez_compressed(path, xy_m=xy, weights_m2=weights, electric_vm=values[:, :2],
                        grid_x_m=x, grid_y_m=y, grid_electric_vm=grid_e.reshape(ny, nx, 2),
                        inside=inside.reshape(ny, nx), outer_m=outer, z_nm=0.0,
                        refinement=refinement, mode_idx=mode, normalization_factor=factor)
    return path


def evaluate_case(runner, sampler, mode, omega, outer, holes, settings, phase, output):
    a = float(runner.model.java.param().evaluate("a"))
    height = float(runner.model.java.param().evaluate("H"))
    sections = section_coordinates(settings, height)
    eps_air = float(runner.config.air_refractive_index) ** 2
    eps_slab = float(runner.config.slab_refractive_index) ** 2
    area = abs(polygon_area(outer))
    edges = hole_edges(holes)
    geometry_rows = [{"hole": e.hole, "edge": e.edge, "group": e.group,
                      "x0_m": e.start[0], "y0_m": e.start[1], "x1_m": e.end[0],
                      "y1_m": e.end[1], "length_m": e.length,
                      "nx_dielectric_to_air": e.normal[0], "ny_dielectric_to_air": e.normal[1]}
                     for e in edges]
    write_csv(output / "edge_geometry.csv", geometry_rows)
    order, divisions = settings["area_order"], settings["area_subdivisions"]
    xy0, w0, _ = area_rule(outer, holes, order, divisions)
    z0 = 0.0
    hz0 = sampler.plane(mode, xy0, z0, ["ewfd.Hz"])[:, 0]
    rms = np.sqrt(np.dot(w0, abs(hz0) ** 2) / area)
    if not np.isfinite(rms) or rms <= 0:
        raise ValueError("Zero/nonfinite Hz RMS at z=0")
    factor = phase / rms
    calibration = sampler.plane(mode, xy0[:5], z0, ["x", "y", "z"])
    if not np.allclose(calibration.real, np.column_stack([xy0[:5], np.full(5, z0)]), rtol=1e-8, atol=1e-14):
        raise RuntimeError("COMSOL interpolation coordinate units do not match SI calibration")
    # Vector-element magnetic fields need elementwise Lagrange reconstruction
    # for spatial derivatives; direct d(H,z) can silently evaluate to zero.
    expressions = ["ewfd.Ex", "ewfd.Ey", "d(laginterp(2,ewfd.Hx),z)/(1[A/m^2])",
                   "d(laginterp(2,ewfd.Hy),z)/(1[A/m^2])"]
    all_edges, totals, comparisons, traces = [], [], [], []
    for z_nm, z, fraction in sections:
        section = {"z_nm": z_nm, "z_over_H": fraction,
                   "section_role": "display" if z_nm == 0 else "reference",
                   "surface_side": "slab_layer_domains_at_exact_z"}
        for refinement in (1, 2):
            edge_order = settings["edge_order"] * refinement
            sample = lambda xy: factor * sampler.plane(mode, xy, z, ["ewfd.Hz"])[:, 0]
            rows = edge_integrals(edges, sample, edge_order)
            group = group_totals(rows)
            for row in rows:
                all_edges.append({**section, "refinement": refinement, **row})
            xy, weights, air = area_rule(outer, holes, order, divisions * refinement)
            inv_eps = np.where(air, 1 / eps_air, 1 / eps_slab)
            values = factor * sampler.plane(mode, xy, z, expressions)
            if z_nm == 0 and refinement == 2:
                export_area_field(output, sampler, mode, outer, xy, weights, values, factor, refinement)
            ex, ey = weights @ values[:, :2] / area
            dx, dy = (weights * inv_eps) @ values[:, 2:]
            outer_rows = edge_integrals(polygon_edges(outer), sample, edge_order)
            ox = sum(r["qx"] for r in outer_rows) / eps_slab
            oy = sum(r["qy"] for r in outer_rows) / eps_slab
            terms = maxwell_terms(group["Qx"], group["Qy"], ox, oy, dx, dy,
                                  omega, area, eps_air, eps_slab)
            scale = max(abs(terms["prefactor"] * terms["jump"]) * group["Dq"]
                        + abs(terms["Zx"]) + abs(terms["Zy"]), 1e-300)
            comparisons.append({**section, "refinement": refinement,
                                "derivative_evaluation": "elementwise_laginterp2",
                                "area_m2": area, "area_quadrature_m2": sum(weights),
                                "area_samples": len(xy), "Hz_RMS_raw": rms,
                                "normalization_factor": factor, "Ex_area": ex, "Ey_area": ey,
                                "Ex_integral": ex * area, "Ey_integral": ey * area,
                                **terms, "amplitude_scale": scale,
                                "Ex_identity_error": abs(ex - terms["predicted_Ex"]) / scale,
                                "Ey_identity_error": abs(ey - terms["predicted_Ey"]) / scale,
                                "outer_error": (abs(terms["Ox"]) + abs(terms["Oy"])) / scale})
            totals.append({**section, "refinement": refinement, **group})
        # Both material traces, two offsets, no corner points. Record actual
        # domain IDs without relying on fixed COMSOL entity numbering.
        for e in edges:
            xy, _ = line_rule(e, min(settings["edge_order"], 24))
            for shrink in (1, 0.5):
                offset = a * settings["trace_offset_over_a"] * shrink
                air_values = sampler.plane(mode, xy + offset * e.normal, z, ["ewfd.Hz", "dom"])
                diel_values = sampler.plane(mode, xy - offset * e.normal, z, ["ewfd.Hz", "dom"])
                trace_scale = max(np.linalg.norm(air_values[:, 0]), np.linalg.norm(diel_values[:, 0]), 1e-300)
                traces.append({"hole": e.hole, "edge": e.edge, **section,
                               "offset_m": offset,
                               "Hz_trace_relative_difference": np.linalg.norm(air_values[:, 0] - diel_values[:, 0]) / trace_scale,
                               "air_domain_ids": json.dumps(np.unique(air_values[:, 1].real).tolist()),
                               "dielectric_domain_ids": json.dumps(np.unique(diel_values[:, 1].real).tolist()),
                               "distinct_domain_each_point": bool(np.all(air_values[:, 1] != diel_values[:, 1]))})
    # Independent G=0 outgoing field in homogeneous air, two heights.
    xy, weights, _ = area_rule(outer, [], order, divisions * 2)
    hair = float(runner.model.java.param().evaluate("H_air"))
    air_rows = []
    for z in (0.6 * hair, 0.75 * hair):
        val = factor * sampler.plane(mode, xy, z, ["ewfd.Ex", "ewfd.Ey"])
        c = weights @ val / area
        # SI e^-iωt: upward outgoing wave exp(+ikz*z).
        c0 = c * np.exp(-1j * omega / 299792458 * np.sqrt(eps_air) * z)
        air_rows.append({"z_m": z, "cx": c[0], "cy": c[1], "cx_at_z0": c0[0], "cy_at_z0": c0[1]})
    write_csv(output / "edge_integrals.csv", all_edges)
    write_csv(output / "group_totals.csv", totals)
    write_csv(output / "area_comparison.csv", comparisons)
    write_csv(output / "interface_traces.csv", traces)
    write_csv(output / "air_radiation.csv", air_rows)
    fine = [r for r in comparisons if r["refinement"] == 2]
    max_error = max(max(r["Ex_identity_error"], r["Ey_identity_error"]) for r in fine)
    coarse = [r for r in comparisons if r["refinement"] == 1]
    quadrature_change = max(max(abs(f["Ex_area"] - c["Ex_area"]), abs(f["Ey_area"] - c["Ey_area"]))
                            / f["amplitude_scale"] for f, c in zip(fine, coarse))
    air_scale = max(np.linalg.norm([air_rows[0]["cx_at_z0"], air_rows[0]["cy_at_z0"]]), 1e-300)
    air_height_error = float(np.linalg.norm([air_rows[1]["cx_at_z0"] - air_rows[0]["cx_at_z0"],
                                            air_rows[1]["cy_at_z0"] - air_rows[0]["cy_at_z0"]]) / air_scale)
    return {"identity_max_error": max_error, "area_quadrature_change": quadrature_change,
            "derivative_evaluation": "elementwise_laginterp2",
            "sections_nm": settings["z_nm"], "primary_z_nm": 0.0, "normalization_z_nm": 0.0,
            "surface_side": "slab_layer_domains_at_exact_z",
            "normalization_factor_re": float(factor.real), "normalization_factor_im": float(factor.imag),
            "all_trace_domains_distinct": all(r["distinct_domain_each_point"] for r in traces),
            "air_height_relative_error": air_height_error,
            "scientific_acceptance": "pending",
            "identity_gate": bool(max_error < settings["identity_tolerance"]),
            "limitations": ["Mesh convergence and physical zero location not established by one case",
                            "Air-height relative error is ill-conditioned near a radiation zero",
                            "Frozen-field and length-only approximations require a separate comparison"]}


def frozen_predictions(sampler, mode, holes, settings, factor, height, targets, current_id):
    """Evaluate all requested geometries in the first solved field per mesh.

    The source field is read at actual moved coordinates, not scaled edge means.
    A separate length-only column keeps those original means unchanged.
    """
    rows = []
    for z_nm, z, fraction in section_coordinates(settings, height):
        sample = lambda xy: factor * sampler.plane(mode, xy, z, ["ewfd.Hz"])[:, 0]
        reference = edge_integrals(hole_edges(holes), sample, settings["edge_order"] * 2)
        for case_id, target_holes in targets:
            edges = hole_edges(target_holes)
            for name, values in (
                ("length_only", length_only_prediction(reference, edges)),
                ("frozen_spatial_Hz", edge_integrals(edges, sample, settings["edge_order"] * 2)),
            ):
                rows.append({"source_case": current_id, "target_case": case_id,
                             "z_nm": z_nm, "z_over_H": fraction, "approximation": name, **group_totals(values)})
    return rows


def run(args):
    preparation = Path(args.run_dir).resolve()
    manifest = prepared_manifest(preparation)
    directory = Path(manifest["output_directory"])
    # Do not silently reuse partially solved models or overwrite previous output.
    if directory.exists():
        raise FileExistsError("Execution already attempted; preserve results and prepare a new run")
    try:
        preflight = runtime_preflight()
    except Exception as exc:
        write_json(preparation / "99_config/preflight.json", {"status": "failed", "error": str(exc)})
        raise
    write_json(preparation / "99_config/preflight.json", preflight)
    print(json.dumps({"output": str(directory), "execution": manifest["execution"],
                      "cases": manifest["cases"], "mode": manifest["s4"]["target"]}, indent=2), flush=True)
    from scripts.run_main import run_band_pair as band
    from comsol_workflow.geometry_utils import create_hexagon_design
    from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationConfig
    import mph
    from jpype import JClass
    previous, rows, reference_meshes, predictions = None, [], set(), []
    try:
        mph.start(version="6.3")
        model_util = JClass("com.comsol.model.util.ModelUtil")
        if not model_util.checkoutLicense("COMSOL", "WAVEOPTICS"):
            raise RuntimeError("COMSOL/WAVEOPTICS license checkout failed")
        preflight["product_checkout"] = "COMSOL and WAVEOPTICS passed"
        write_json(preparation / "99_config/preflight.json", preflight)
        directory.mkdir(parents=True)
        for relative in ["99_config/s4_manifest.json", "99_config/parameter_source.json",
                         *[c["snapshot"] for c in manifest["cases"]]]:
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(preparation / relative, target)
        write_json(directory / "99_config/preflight.json", preflight)
        write_json(directory / "99_config/execution.json", {"status": "starting", "pid": os.getpid()})
        for c in manifest["cases"]:
            case_dir = directory / c["case_id"]
            output = case_dir / "01_results"
            output.mkdir()
            model_dir = case_dir / "00_model"
            model_dir.mkdir()
            params = load_shared_parameters(directory / c["snapshot"])
            hp = band.get_hole_params(params.unit_cell_2d.cell.to_fourier_params("s4_boundary"))
            outer, holes, info = create_hexagon_design(band.A, hp)
            if min(v["min_dist"] for v in info) < band.D:
                raise ValueError("Geometry violates existing minimum-clearance contract")
            config = SimulationConfig(slab_height=band.SLAB_HEIGHT, slab_refractive_index=band.REFRACTIVE_INDEX,
                                      wavelength=band.WAVELENGTH, mode_type=band.MODE_TYPE,
                                      eigenmode_count=params.unit_cell_2d.eigenmode_pair_count,
                                      eigenfrequency_shift="c_const/1.55[um]", mesh_auto_size=params.mesh_size)
            write_json(case_dir / "99_config/realized_geometry.json",
                       {"length_unit": "um", "outer": np.asarray(outer).tolist(),
                        "holes": np.asarray(holes).tolist(), "hole_params": hp.tolist(),
                        "simulation_config": asdict(config), "k_over_G": [0, 0]})
            print("Solving Gamma: " + c["case_id"], flush=True)
            with band.ReusableSimulationRun(config) as runner:
                model_util.showProgress(str(model_dir / "comsol_progress.log"))
                runner.build_and_run(band.A, [np.asarray(h).tolist() for h in holes],
                                     {"kx": 0.0, "ky": 0.0}, mesh_auto_size=params.mesh_size)
                runner.model.save(str(model_dir / "s4_gamma.mph"))
                print("Gamma solved; checkpoint saved; extracting fields", flush=True)
                modes = band.persist_k_point_solution(runner, output)
                sampler = Sampler(runner)
                height = float(runner.model.java.param().evaluate("H"))
                mode, chosen, previous, phase = select_mode(runner, sampler, modes, band.A * 1e-6,
                    0.0, output, manifest["s4"], previous)
                omega = 2 * np.pi * complex(modes.loc[mode, "re"], -modes.loc[mode, "im"]) * 1e12
                summary = evaluate_case(runner, sampler, mode, omega, np.asarray(outer) * 1e-6,
                    [np.asarray(h) * 1e-6 for h in holes], manifest["s4"], phase, output)
                if c["mesh"] not in reference_meshes:
                    targets = []
                    for target in manifest["cases"]:
                        if target["mesh"] != c["mesh"]:
                            continue
                        tp = load_shared_parameters(directory / target["snapshot"])
                        thp = band.get_hole_params(tp.unit_cell_2d.cell.to_fourier_params("s4_boundary"))
                        _, tholes, tinfo = create_hexagon_design(band.A, thp)
                        if min(v["min_dist"] for v in tinfo) < band.D:
                            raise ValueError("Prediction geometry violates minimum clearance")
                        targets.append((target["case_id"], [np.asarray(h) * 1e-6 for h in tholes]))
                    factor = complex(summary["normalization_factor_re"], summary["normalization_factor_im"])
                    predictions.extend(frozen_predictions(sampler, mode,
                        [np.asarray(h) * 1e-6 for h in holes], manifest["s4"], factor,
                        height, targets, c["case_id"]))
                    write_csv(directory / "80_logs/frozen_field_predictions.csv", predictions)
                    reference_meshes.add(c["mesh"])
                row = {**c, **chosen, **summary, "status": "computed_pending_review"}
                write_json(case_dir / "99_config/s4_summary.json", row)
                from comsol_workflow.boundary_plotting import plot_s4_case
                plot_s4_case(case_dir)
                rows.append(row)
                write_json(directory / "99_config/execution.json", {"status": "running", "pid": os.getpid(), "cases": rows})
        write_json(directory / "99_config/execution.json", {"status": "complete_pending_scientific_review", "cases": rows})
        report(directory)
    except Exception as exc:
        write_json((directory if directory.exists() else preparation) / "99_config/execution.json", {"status": "failed", "cases": rows,
                   "error": str(exc), "traceback": traceback.format_exc()})
        raise
    finally:
        if digest(manifest["source_parameter"]) != manifest["source_parameter_sha256"]:
            print("WARNING: shared parameters changed since preparation; this run used immutable snapshots", file=sys.stderr)


def report(directory):
    directory = Path(directory).resolve()
    manifest = json_read(directory / "99_config/s4_manifest.json")
    logs = directory / "80_logs"
    logs.mkdir(parents=True, exist_ok=True)
    summaries = []
    for c in manifest["cases"]:
        path = directory / c["case_id"] / "99_config/s4_summary.json"
        if path.exists():
            summaries.append(json_read(path))
    write_csv(logs / "case_table.csv", [{k: v for k, v in s.items() if not isinstance(v, (list, dict))} for s in summaries])
    # Keep physical air amplitude, corrected section Ey and bare boundary sum
    # distinct. None is substituted for another when estimating a minimum.
    series = {}
    for c in manifest["cases"]:
        result = data_directory(directory / c["case_id"])
        for filename, quantity, filter_key in (
            ("group_totals.csv", "Qx", "refinement"),
            ("area_comparison.csv", "Ey_area", "refinement"),
            ("air_radiation.csv", "cy_at_z0", None),
        ):
            path = result / filename
            if not path.exists():
                continue
            with path.open(encoding="utf-8", newline="") as f:
                for r in csv.DictReader(f):
                    if filter_key and int(r[filter_key]) != 2:
                        continue
                    section = r.get("z_over_H", r.get("z_m"))
                    key = (c["eta"], c["mesh"], quantity, section)
                    v = complex(float(r[quantity + "_re"]), float(r[quantity + "_im"]))
                    series.setdefault(key, []).append((c["zeta"], v))
    estimates = []
    for (eta, mesh, quantity, section), values in series.items():
        values.sort()
        for (z0, v0), (z1, v1) in zip(values[:-1], values[1:]):
            estimate = segment_minimum(z0, z1, v0, v1)
            if estimate is not None:
                estimates.append({"eta": eta, "mesh": mesh, "quantity": quantity,
                                  "section": section, "zeta_left": z0, "zeta_right": z1,
                                  "scientific_acceptance": "requires_solved_point_and_convergence",
                                  **estimate})
    write_csv(logs / "zero_estimates.csv", estimates)
    text = ["# S4 execution report", "", f"Run: {manifest['run_id']}", "",
            f"Computed cases: {len(summaries)} / {len(manifest['cases'])}.", "",
            "Scientific acceptance remains pending. A completed solve is not proof of a radiation zero.", "",
            "Inspect mode_identification.csv, edge_integrals.csv, group_totals.csv, area_comparison.csv,",
            "interface_traces.csv and air_radiation.csv in each case's 80_logs (legacy: 01_results).", "",
            "The 3D correction terms and the outgoing air-plane amplitude are distinct outputs.", "",
            "frozen_field_predictions.csv compares fixed edge means and the fixed spatial Hz field.", "",
            "When available, zero_estimates.csv contains complex linear-interpolation diagnostics, not solved zeros.", "",
            "Mesh convergence and A/B parameter trends still require scientific review."]
    (logs / "acceptance.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    print(str(logs / "acceptance.md"))


def check_model(args):
    """Read an existing solution to verify the backend; never solve or save it."""
    source = Path(args.model).resolve()
    if source.suffix.lower() != ".mph" or not source.is_file():
        raise ValueError("Supply an existing MPH checkpoint")
    output = TASK_ROOT / ("s4_backend_check_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:6])
    output.mkdir(parents=True)
    source_sha = digest(source)
    audit = {"source_mph": str(source), "source_sha256": source_sha,
             "eigensolves": 0, "purpose": "API validation only, not S4 scientific evidence"}
    model, client = None, None
    try:
        audit["preflight"] = runtime_preflight()
        import mph
        client = mph.start(version="6.3")
        model = client.load(str(source))
        a = float(model.java.param().evaluate("a"))
        height = float(model.java.param().evaluate("H"))
        sampler = Sampler(SimpleNamespace(model=model))
        xyz = np.array([[.15 * a, .1 * a, .2 * height],
                        [-.15 * a, .1 * a, .2 * height], [.2 * a, -.1 * a, .3 * height]])
        expressions = ["x", "y", "z", "ewfd.Ex", "ewfd.Ey", "ewfd.Hz",
                       "d(ewfd.Hx,z)", "d(ewfd.Hy,z)", "dom"]
        values = sampler(0, xyz, expressions)
        audit["requested_coordinates_m"] = xyz.tolist()
        audit["returned_coordinates_m"] = values[:, :3].real.tolist()
        if not np.allclose(values[:, :3].real, xyz, rtol=1e-8, atol=1e-14):
            raise RuntimeError("Interpolation-coordinate SI calibration failed")
        audit.update(status="passed", coordinate_calibration="passed",
                     kx_per_m=float(model.java.param().evaluate("kx")),
                     ky_per_m=float(model.java.param().evaluate("ky")),
                     field_expressions=expressions, samples_finite=True,
                     imaginary_Hz_max=float(np.max(abs(values[:, 5].imag))))
        write_csv(output / "api_samples.csv", [dict(zip(expressions, row)) for row in values])
    except Exception as exc:
        audit.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        if model is not None:
            client.remove(model)
        audit["source_unchanged"] = source_sha == digest(source)
        write_json(output / "backend_check.json", audit)
        print(str(output), flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare", help="write snapshots; do not import/start COMSOL")
    prep.add_argument("--parameter-file")
    prep.add_argument("--run-id")
    prep.add_argument("--eta", type=float)
    prep.add_argument("--zeta", type=float, nargs="+")
    prep.add_argument("--mesh", type=int, nargs="+")
    prep.add_argument("--max-cases", type=int, default=1)
    prep.add_argument("--edge-order", type=int, default=64)
    prep.add_argument("--area-order", type=int, default=8)
    prep.add_argument("--area-subdivisions", type=int, default=2)
    prep.add_argument("--z-nm", type=float, nargs="+", default=[0.0, 50.0, 100.0])
    prep.add_argument("--trace-offset-over-a", type=float, default=1e-5)
    execute = sub.add_parser("run", help="explicitly start Gamma eigensolves for a prepared run")
    execute.add_argument("run_dir")
    rep = sub.add_parser("report", help="read completed outputs; do not start COMSOL")
    rep.add_argument("run_dir")
    check = sub.add_parser("check-model", help="read an existing MPH to verify APIs; no eigensolve or save")
    check.add_argument("model")
    args = p.parse_args(argv)
    if args.command == "prepare":
        if (args.edge_order < 4 or args.area_order < 2 or args.area_subdivisions < 1
                or not args.z_nm or args.z_nm[0] != 0 or not np.isfinite(args.z_nm).all()
                or any(b <= a for a, b in zip(args.z_nm, args.z_nm[1:]))
                or not 0 < args.trace_offset_over_a < 1e-3):
            p.error("Invalid quadrature/section/trace settings; sections must start at z=0 and increase")
        prepare(args)
    elif args.command == "run":
        run(args)
    elif args.command == "report":
        report(args.run_dir)
    else:
        check_model(args)


if __name__ == "__main__":
    main()
