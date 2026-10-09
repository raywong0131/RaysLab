#!/usr/bin/env python3
"""Scan c(k) around Gamma for the cavity unit cell and run analysis by default."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.band_connector import (  # noqa: E402
    assign_mode_candidates,
    field_overlap,
    load_hz_center,
)
from comsol_workflow.geometry_utils import (  # noqa: E402
    create_hexagon_design,
    visualize_hexagon_design,
)
from comsol_workflow.simulation_spatial.hexagon_unit_cell import (  # noqa: E402
    SimulationConfig,
)
from comsol_workflow.unit_cell_2d_analysis import (  # noqa: E402
    POLARIZATION_DEFINITION_VERSION,
    SYMMETRY_GAUGE_VERSION,
    SYMMETRY_RENDER_VERSION,
    formal_mode_label,
    run_analysis,
    unit_cell_2d_series_name,
)
from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import (  # noqa: E402
    UNIT_CELL_2D_OUTPUT_ROOT,
    load_shared_parameters,
    parameter_path_from_environment,
    unit_cell_2d_metadata,
)


RESULTS_DIRNAME = "01_results"
OVERVIEW_DIRNAME = "10_overview"
MODEL_DIRNAME = "00_model"
LOGS_DIRNAME = "80_logs"
CONFIG_DIRNAME = "99_config"
POINT_METADATA_FILENAME = "point_metadata.json"
MODEL_FILENAME = "unit_cell_2D.mph"
K_STRING_PRECISION = 6
SCAN_DEFINITION_VERSION = "unit-cell-2d-cartesian-quarter-or-full-v2"
TRACKING_DEFINITION_VERSION = "gamma-anchor-multipred-hungarian-repair-v1"
GEOMETRY_DEFINITION_VERSION = "hexagon-unit-cell-2d-fourier-v2"
FIELD_EXTRACTION_VERSION = "air-plane-batched-ex-ey-planar-l2-v1"


def make_grid_points(
    q_axis_over_g: Sequence[float],
    *,
    sampling_domain: str = "full",
) -> list[dict[str, object]]:
    if sampling_domain not in {"quarter", "full"}:
        raise ValueError("sampling_domain must be 'quarter' or 'full'")
    axis = tuple(float(value) for value in q_axis_over_g)
    points = []
    for qy in axis:
        for qx in axis:
            qx = 0.0 if abs(qx) < 0.5 * 10 ** (-K_STRING_PRECISION) else qx
            qy = 0.0 if abs(qy) < 0.5 * 10 ** (-K_STRING_PRECISION) else qy
            points.append(
                {
                    "direction": "unit_cell_2d",
                    "kx": qx,
                    "ky": qy,
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "kx_str": f"{qx:.{K_STRING_PRECISION}f}",
                    "ky_str": f"{qy:.{K_STRING_PRECISION}f}",
                    "k_norm": math.hypot(qx, qy),
                    "is_comsol_solved": True,
                    "sampling_domain": sampling_domain,
                }
            )
    points.sort(
        key=lambda point: (
            max(abs(float(point["kx"])), abs(float(point["ky"]))),
            float(point["k_norm"]),
            float(point["ky"]),
            float(point["kx"]),
        )
    )
    for solve_order, point in enumerate(points):
        point["solve_order"] = solve_order
    return points


def resolved_series_dir(parameters) -> Path:
    settings = parameters.unit_cell_2d
    name = unit_cell_2d_series_name(
        b0_nm=settings.cell.b0_nm,
        eta=settings.cell.eta,
        zeta=settings.cell.zeta,
        mesh_size=parameters.mesh_size,
        q_max_over_g=settings.q_max_over_g,
        q_points_per_axis=settings.q_points_per_axis,
        target_bands=settings.target_bands,
        is_quarter=settings.is_quarter,
    )
    return UNIT_CELL_2D_OUTPUT_ROOT / name


def point_key(point: dict[str, object]) -> tuple[float, float]:
    return float(point["kx"]), float(point["ky"])


def _point_dir(results_dir: Path, point: dict[str, object]) -> Path:
    return band.k_point_dir(results_dir, point)


def _save_runner_model(runner, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    runner.model.save(str(path))


def _cache_identity(parameters) -> dict[str, object]:
    return {
        "workflow": "run_unit_cell_2d",
        "mesh_size": parameters.mesh_size,
        "eigenmode_pair_count": parameters.unit_cell_2d.eigenmode_pair_count,
        "expected_eigenvalue_count": (
            2 * parameters.unit_cell_2d.eigenmode_pair_count
        ),
        "eigenfrequency_shift": "c_const/1.55[um]",
        "comsol_target_version": "6.3",
        "geometry_definition_version": GEOMETRY_DEFINITION_VERSION,
        "scan_definition_version": SCAN_DEFINITION_VERSION,
        "tracking_definition_version": TRACKING_DEFINITION_VERSION,
        "field_extraction_version": FIELD_EXTRACTION_VERSION,
        "G_definition": "4*pi/sqrt(3)/a",
        "k_string_precision": K_STRING_PRECISION,
        "isQuarter": parameters.unit_cell_2d.is_quarter,
        "sampling_domain": parameters.unit_cell_2d.sampling_domain,
        "solved_q_axis_over_G": list(parameters.unit_cell_2d.q_axis_over_g),
        "plot_q_axis_over_G": list(parameters.unit_cell_2d.plot_q_axis_over_g),
        "solved_point_count": parameters.unit_cell_2d.q_points_per_axis**2,
        "plot_coordinate_count": len(
            parameters.unit_cell_2d.plot_q_axis_over_g
        ) ** 2,
        "symmetry_group": "C2v",
        "symmetry_render_version": SYMMETRY_RENDER_VERSION,
        "symmetry_gauge_version": SYMMETRY_GAUGE_VERSION,
        "unit_cell_2d": unit_cell_2d_metadata(parameters.unit_cell_2d),
        "polarization_definition_version": POLARIZATION_DEFINITION_VERSION,
    }


def _load_cached_rows(
    point_dir: Path,
    identity: dict[str, object],
    target_bands: Sequence[str],
) -> list[dict[str, object]] | None:
    try:
        payload = json.loads(
            (point_dir / POINT_METADATA_FILENAME).read_text(encoding="utf-8")
        )
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None
    rows = payload.get("rows")
    if payload.get("cache_identity") != identity or not isinstance(rows, list):
        return None
    if {row.get("band_label") for row in rows} != set(target_bands):
        return None
    if not (point_dir / "eigenfrequencies.csv").is_file():
        return None
    for row in rows:
        if row.get("match_status") not in {"gamma", "matched", "repaired"}:
            continue
        try:
            mode_idx = int(row["mode_idx"])
        except (KeyError, TypeError, ValueError, OverflowError):
            return None
        field_path = point_dir / "eigenmodes" / f"{mode_idx:02d}_Hz_center.parquet"
        if not field_path.is_file():
            return None
    return rows


def _inward_predecessors(
    point: dict[str, object],
    states: dict[tuple[float, float], dict[str, object]],
    limit: int = 3,
) -> list[dict[str, object]]:
    radius = float(point["k_norm"])
    qx, qy = point_key(point)
    eligible = []
    for (px, py), state in states.items():
        predecessor_radius = math.hypot(px, py)
        if predecessor_radius >= radius - 1e-15:
            continue
        eligible.append((math.hypot(qx - px, qy - py), state))
    eligible.sort(key=lambda item: item[0])
    return [state for _distance, state in eligible[:limit]]


def _nearest_predecessors(
    point: dict[str, object],
    states: dict[tuple[float, float], dict[str, object]],
    limit: int = 4,
) -> list[dict[str, object]]:
    """Return deterministic neighbours for the single post-scan repair pass."""
    qx, qy = point_key(point)
    eligible = [
        (math.hypot(qx - px, qy - py), py, px, state)
        for (px, py), state in states.items()
        if not (math.isclose(px, qx, abs_tol=1e-15)
                and math.isclose(py, qy, abs_tol=1e-15))
    ]
    eligible.sort(key=lambda item: (item[0], item[1], item[2]))
    return [state for _distance, _py, _px, state in eligible[:limit]]


def _complex_columns(prefix: str, value: complex) -> dict[str, float]:
    return {
        f"{prefix}_re": float(complex(value).real),
        f"{prefix}_im": float(complex(value).imag),
    }


def _invalid_row(
    point: dict[str, object],
    band_label: str,
    status: str,
    *,
    evidence: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        **point,
        "band_label": band_label,
        "formal_band_label": formal_mode_label(band_label),
        "mode_idx": None,
        "match_status": status,
        "frequency_thz": None,
        "q_factor": None,
        "selected_score": None,
        "runner_up_score": None,
        "ambiguity_gap": None,
        "frequency_delta_thz": None,
        "candidate_count": 0,
        "predecessors": "[]",
        "candidate_scores": "{}",
        "cx_raw_re": None,
        "cx_raw_im": None,
        "cy_raw_re": None,
        "cy_raw_im": None,
        "normalization_kind": "planar_l2",
        "normalization_value": None,
        "cx_norm_re": None,
        "cx_norm_im": None,
        "cy_norm_re": None,
        "cy_norm_im": None,
        "field_plane": "air",
        "polarization_definition_version": POLARIZATION_DEFINITION_VERSION,
        **(evidence or {}),
    }


def select_point_modes(
    point: dict[str, object],
    modes: pd.DataFrame,
    point_dir: Path,
    target_bands: Sequence[str],
    states_by_band: dict[str, dict[tuple[float, float], dict[str, object]]],
    settings,
    *,
    predecessor_strategy: str = "inward",
) -> dict[str, dict[str, object]]:
    valid_modes = modes[modes["is_valid"].astype(bool)]
    candidate_indices = [int(index) for index in valid_modes.index]
    fields = {
        mode_idx: load_hz_center(point_dir, mode_idx)
        for mode_idx in candidate_indices
    }
    labels = []
    predecessor_sets = {}
    previous_frequencies = []
    overlap_rows = []
    for label in target_bands:
        if predecessor_strategy == "inward":
            predecessors = _inward_predecessors(point, states_by_band[label])
        elif predecessor_strategy == "nearest":
            predecessors = _nearest_predecessors(point, states_by_band[label])
        else:
            raise ValueError(f"Unknown predecessor strategy: {predecessor_strategy}")
        if not predecessors:
            continue
        labels.append(label)
        predecessor_sets[label] = predecessors
        previous_frequencies.append(
            float(np.mean([item["frequency_thz"] for item in predecessors]))
        )
        overlap_rows.append(
            [
                float(
                    np.mean(
                        [
                            field_overlap(item["field"], fields[mode_idx])
                            for item in predecessors
                        ]
                    )
                )
                for mode_idx in candidate_indices
            ]
        )
    selected: dict[str, dict[str, object]] = {}
    if not labels or not candidate_indices:
        return selected
    overlaps = np.asarray(overlap_rows, dtype=float)
    candidate_frequencies = [
        float(modes.loc[index, "re"]) for index in candidate_indices
    ]
    assignments = assign_mode_candidates(
        previous_frequencies,
        candidate_frequencies,
        overlaps,
        settings.frequency_tolerance_thz,
        settings.minimum_overlap,
    )
    for row_index, (label, assignment) in enumerate(zip(labels, assignments)):
        scores = overlaps[row_index]
        if not assignment["matched"]:
            continue
        candidate_position = int(assignment["candidate_index"])
        mode_idx = candidate_indices[candidate_position]
        selected_score = float(scores[candidate_position])
        frequency_differences = np.abs(
            np.asarray(candidate_frequencies, dtype=float)
            - previous_frequencies[row_index]
        )
        alternative_mask = (
            (frequency_differences <= settings.frequency_tolerance_thz)
            & (scores >= settings.minimum_overlap)
        )
        alternative_mask[candidate_position] = False
        alternatives = scores[alternative_mask]
        runner_up = float(np.max(alternatives)) if alternatives.size else 0.0
        selected[label] = {
            "mode_idx": mode_idx,
            "field": fields[mode_idx],
            "frequency_thz": float(modes.loc[mode_idx, "re"]),
            "q_factor": float(modes.loc[mode_idx, "q"]),
            "selected_score": selected_score,
            "runner_up_score": runner_up,
            "ambiguity_gap": selected_score - runner_up,
            "frequency_delta_thz": float(assignment["frequency_difference"]),
            "candidate_count": len(candidate_indices),
            "predecessors": [
                [item["qx_over_G"], item["qy_over_G"]]
                for item in predecessor_sets[label]
            ],
            "candidate_scores": {
                str(index): float(score)
                for index, score in zip(candidate_indices, scores)
            },
        }
    return selected


def _match_passes_ambiguity(match: dict[str, object], settings) -> bool:
    return float(match["ambiguity_gap"]) >= settings.minimum_ambiguity_gap


def _matched_row(
    point: dict[str, object],
    band_label: str,
    match: dict[str, object],
    polarization: dict[str, object],
    status: str,
) -> dict[str, object]:
    return {
        **point,
        "band_label": band_label,
        "formal_band_label": formal_mode_label(band_label),
        "mode_idx": int(match["mode_idx"]),
        "match_status": status,
        "frequency_thz": float(match["frequency_thz"]),
        "q_factor": float(match["q_factor"]),
        "selected_score": float(match["selected_score"]),
        "runner_up_score": float(match["runner_up_score"]),
        "ambiguity_gap": float(match["ambiguity_gap"]),
        "frequency_delta_thz": float(match["frequency_delta_thz"]),
        "candidate_count": int(match["candidate_count"]),
        "predecessors": json.dumps(match["predecessors"]),
        "candidate_scores": json.dumps(match["candidate_scores"]),
        **_complex_columns("cx_raw", polarization["cx_raw"]),
        **_complex_columns("cy_raw", polarization["cy_raw"]),
        "normalization_kind": polarization["normalization_kind"],
        "normalization_value": float(polarization["normalization_value"]),
        **_complex_columns("cx_norm", polarization["cx_norm"]),
        **_complex_columns("cy_norm", polarization["cy_norm"]),
        "field_plane": polarization["field_plane"],
        "polarization_definition_version": POLARIZATION_DEFINITION_VERSION,
    }


def _restore_cached_states(
    rows: Sequence[dict[str, object]],
    point_dir: Path,
    states_by_band: dict[str, dict[tuple[float, float], dict[str, object]]],
) -> None:
    for row in rows:
        if row["match_status"] not in {"gamma", "matched", "repaired"}:
            continue
        mode_idx = int(row["mode_idx"])
        state = dict(row)
        state["field"] = load_hz_center(point_dir, mode_idx)
        states_by_band[str(row["band_label"])][
            (float(row["qx_over_G"]), float(row["qy_over_G"]))
        ] = state


def _prepare_geometry(parameters, model_dir: Path):
    unit_cell_params = parameters.unit_cell_2d.cell.to_fourier_params(
        "unit_cell_2d_p_bic"
    )
    hole_params = band.get_hole_params(unit_cell_params)
    hexagon, holes, info = create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            "Cavity unit cell violates clearance: "
            f"min_dist={minimum_clearance:g}, d={band.D:g}"
        )
    visualize_hexagon_design(
        hexagon,
        holes,
        info,
        filename=model_dir / "geometry.png",
        annotate=False,
    )
    return unit_cell_params, hole_params, holes, minimum_clearance


def run_unit_cell_2d(
    output_dir: Path | None = None,
    *,
    skip_analysis: bool = False,
) -> dict[str, object]:
    parameters = load_shared_parameters(parameter_path_from_environment())
    settings = parameters.unit_cell_2d
    if settings.normalization_kind != "planar_l2":
        raise ValueError("Only planar_l2 normalization is currently implemented")
    series_dir = resolved_series_dir(parameters) if output_dir is None else Path(output_dir)
    results_dir = series_dir / RESULTS_DIRNAME
    overview_dir = series_dir / OVERVIEW_DIRNAME
    model_dir = series_dir / MODEL_DIRNAME
    logs_dir = series_dir / LOGS_DIRNAME
    config_dir = series_dir / CONFIG_DIRNAME
    for directory in (model_dir, results_dir, overview_dir, logs_dir, config_dir):
        directory.mkdir(parents=True, exist_ok=True)

    identity = _cache_identity(parameters)
    config_path = config_dir / "config.json"
    if config_path.exists():
        existing = json.loads(config_path.read_text(encoding="utf-8"))
        if existing.get("cache_identity") != identity:
            raise RuntimeError(
                f"Existing unit-cell 2D series has a different identity: {series_dir}"
            )
    points = make_grid_points(
        settings.q_axis_over_g,
        sampling_domain=settings.sampling_domain,
    )
    unit_cell_params, hole_params, holes, minimum_clearance = _prepare_geometry(
        parameters, model_dir
    )
    simulation_config = SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift="c_const/1.55[um]",
        eigenmode_count=settings.eigenmode_pair_count,
        mesh_auto_size=parameters.mesh_size,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )
    band.write_json_atomic(
        config_path,
        {
            "workflow": "run_unit_cell_2d",
            "cache_identity": identity,
            "parameter_path": str(parameters.source_path),
            "geometry_source": "scripts/parameter.json:unit_cell_2d",
            "unit_cell_2d_geometry": settings.cell.to_dict(),
            "unit_cell_2d_fourier_params": unit_cell_params,
            "hole_params": hole_params,
            "minimum_clearance": minimum_clearance,
            "lattice_constant_um": band.A,
            "G_per_um": 4.0 * math.pi / math.sqrt(3.0) / band.A,
            "simulation": band.point_solver_identity(simulation_config),
            "points": points,
            "automatic_analysis": {
                "configured": settings.analysis_enabled,
                "skip_analysis": bool(skip_analysis),
            },
        },
    )

    print(f"Unit-cell 2D output: {series_dir}")
    print(f"Geometry source: {parameters.source_path} -> unit_cell_2d")
    print(
        f"Solved grid ({settings.sampling_domain}): "
        f"{len(settings.q_axis_over_g)}x{len(settings.q_axis_over_g)} "
        f"({len(points)} real COMSOL points), "
        f"plot={len(settings.plot_q_axis_over_g)}x"
        f"{len(settings.plot_q_axis_over_g)}, "
        f"target={','.join(formal_mode_label(item) for item in settings.target_bands)}"
    )
    print(
        f"Mesh={parameters.mesh_size}, "
        f"eigenvalue_pairs={settings.eigenmode_pair_count}, "
        f"expected_eigenvalues={2 * settings.eigenmode_pair_count}, "
        "shift=c_const/1.55[um]"
    )

    rows: list[dict[str, object]] = []
    states_by_band = {label: {} for label in settings.target_bands}
    runner = None
    cache_hits = 0
    solve_count = 0
    repair_attempt_count = 0
    repair_success_count = 0
    repair_solve_count = 0
    hole_vertices = [hole.tolist() for hole in holes]
    try:
        for point_index, point in enumerate(points, start=1):
            point_dir = _point_dir(results_dir, point)
            cached_rows = _load_cached_rows(
                point_dir, identity, settings.target_bands
            )
            if cached_rows is not None:
                cache_hits += 1
                rows.extend(cached_rows)
                _restore_cached_states(cached_rows, point_dir, states_by_band)
                print(
                    f"[{point_index}/{len(points)}] "
                    f"kx={point['kx_str']} ky={point['ky_str']} (cached)"
                )
                continue
            print(
                f"[{point_index}/{len(points)}] "
                f"kx={point['kx_str']} ky={point['ky_str']} (run)"
            )
            if runner is None:
                runner = band.ReusableSimulationRun(simulation_config)
            runner.build_and_run(
                band.A,
                hole_vertices,
                {"kx": float(point["kx"]), "ky": float(point["ky"])},
                mesh_auto_size=parameters.mesh_size,
            )
            solve_count += 1
            modes = band.persist_k_point_solution(runner, point_dir)
            point_rows: list[dict[str, object]] = []
            if math.isclose(float(point["k_norm"]), 0.0, abs_tol=1e-15):
                gamma_modes = band.select_gamma_bands(modes)
                selected = {}
                for label in settings.target_bands:
                    if label not in gamma_modes:
                        raise ValueError(f"Gamma selection has no target band {label!r}")
                    mode_idx = int(gamma_modes[label])
                    selected[label] = {
                        "mode_idx": mode_idx,
                        "field": load_hz_center(point_dir, mode_idx),
                        "frequency_thz": float(modes.loc[mode_idx, "re"]),
                        "q_factor": float(modes.loc[mode_idx, "q"]),
                        "selected_score": 1.0,
                        "runner_up_score": 0.0,
                        "ambiguity_gap": 1.0,
                        "frequency_delta_thz": 0.0,
                        "candidate_count": int(modes["is_valid"].astype(bool).sum()),
                        "predecessors": [],
                        "candidate_scores": {str(mode_idx): 1.0},
                    }
                status = "gamma"
            else:
                selected = select_point_modes(
                    point,
                    modes,
                    point_dir,
                    settings.target_bands,
                    states_by_band,
                    settings,
                )
                status = "matched"
            for label in settings.target_bands:
                match = selected.get(label)
                if match is None:
                    point_rows.append(
                        _invalid_row(point, label, "missing")
                    )
                    continue
                if (
                    status != "gamma"
                    and not _match_passes_ambiguity(match, settings)
                ):
                    point_rows.append(
                        _invalid_row(
                            point,
                            label,
                            "ambiguous",
                            evidence={
                                "selected_score": match["selected_score"],
                                "runner_up_score": match["runner_up_score"],
                                "ambiguity_gap": match["ambiguity_gap"],
                                "candidate_count": match["candidate_count"],
                                "predecessors": json.dumps(match["predecessors"]),
                                "candidate_scores": json.dumps(
                                    match["candidate_scores"]
                                ),
                            },
                        )
                    )
                    continue
                polarization = runner.compute_polarization_data(
                    int(match["mode_idx"]), settings.field_plane
                )
                row = _matched_row(point, label, match, polarization, status)
                point_rows.append(row)
                state = dict(row)
                state["field"] = match["field"]
                states_by_band[label][point_key(point)] = state
            band.write_json_atomic(
                point_dir / POINT_METADATA_FILENAME,
                {"cache_identity": identity, "rows": point_rows},
            )
            rows.extend(point_rows)

        # One deterministic repair pass. It may use completed neighbours from
        # any direction, but it never interpolates c(k): every accepted repair
        # re-solves that k point before extracting its polarization.
        for point in points:
            point_rows = [
                row
                for row in rows
                if math.isclose(
                    float(row["qx_over_G"]), float(point["qx_over_G"]),
                    abs_tol=1e-15,
                )
                and math.isclose(
                    float(row["qy_over_G"]), float(point["qy_over_G"]),
                    abs_tol=1e-15,
                )
            ]
            invalid_labels = [
                str(row["band_label"])
                for row in point_rows
                if row["match_status"] in {"missing", "ambiguous"}
            ]
            if not invalid_labels:
                continue
            repair_attempt_count += len(invalid_labels)
            point_dir = _point_dir(results_dir, point)
            if runner is None:
                runner = band.ReusableSimulationRun(simulation_config)
            runner.build_and_run(
                band.A,
                hole_vertices,
                {"kx": float(point["kx"]), "ky": float(point["ky"])},
                mesh_auto_size=parameters.mesh_size,
            )
            solve_count += 1
            repair_solve_count += 1
            # Refresh both the candidate table and Hz fields after re-solving;
            # near a degeneracy, solnum ordering from the earlier solve is not
            # assumed to remain stable.
            modes = band.persist_k_point_solution(runner, point_dir)
            selected = select_point_modes(
                point,
                modes,
                point_dir,
                invalid_labels,
                states_by_band,
                settings,
                predecessor_strategy="nearest",
            )
            accepted = {
                label: match
                for label, match in selected.items()
                if _match_passes_ambiguity(match, settings)
            }
            if not accepted:
                continue
            for label, match in accepted.items():
                polarization = runner.compute_polarization_data(
                    int(match["mode_idx"]), settings.field_plane
                )
                replacement = _matched_row(
                    point, label, match, polarization, "repaired"
                )
                existing = next(
                    row for row in point_rows if row["band_label"] == label
                )
                existing.clear()
                existing.update(replacement)
                state = dict(replacement)
                state["field"] = match["field"]
                states_by_band[label][point_key(point)] = state
                repair_success_count += 1
            band.write_json_atomic(
                point_dir / POINT_METADATA_FILENAME,
                {"cache_identity": identity, "rows": point_rows},
            )
        if runner is not None:
            _save_runner_model(runner, model_dir / MODEL_FILENAME)
    finally:
        if runner is not None:
            runner.clear()

    grid = pd.DataFrame(rows).sort_values(
        ["band_label", "qy_over_G", "qx_over_G"]
    )
    band.write_csv_atomic(logs_dir / "polarization_grid.csv", grid)
    grid.to_parquet(results_dir / "polarization_grid.parquet", index=False)
    tracking_columns = [
        "qx_over_G",
        "qy_over_G",
        "band_label",
        "mode_idx",
        "match_status",
        "frequency_thz",
        "q_factor",
        "selected_score",
        "runner_up_score",
        "ambiguity_gap",
        "frequency_delta_thz",
        "candidate_count",
        "predecessors",
        "candidate_scores",
    ]
    band.write_csv_atomic(
        logs_dir / "band_tracking_2d.csv", grid[tracking_columns]
    )
    summary = {
        "status": "solved",
        "series_dir": str(series_dir),
        "point_count": len(points),
        "isQuarter": settings.is_quarter,
        "sampling_domain": settings.sampling_domain,
        "solved_q_axis_over_G": list(settings.q_axis_over_g),
        "plot_q_axis_over_G": list(settings.plot_q_axis_over_g),
        "solved_point_count": len(points),
        "plot_coordinate_count": len(settings.plot_q_axis_over_g) ** 2,
        "plot_coordinate_count_is_data_point_count": False,
        "row_count": len(grid),
        "cache_hits": cache_hits,
        "solve_count": solve_count,
        "repair_attempt_count": repair_attempt_count,
        "repair_success_count": repair_success_count,
        "repair_solve_count": repair_solve_count,
        "match_status_counts": grid["match_status"].value_counts().to_dict(),
        "analysis_requested": not skip_analysis,
        "analysis_completed": False,
    }
    summary_path = config_dir / "run_summary.json"
    band.write_json_atomic(summary_path, summary)
    if not skip_analysis:
        try:
            analysis_summary = run_analysis(
                series_dir,
                angle_methods=settings.angle_methods,
                mask_relative_threshold=settings.mask_relative_threshold,
            )
        except Exception:
            summary["status"] = "analysis_failed"
            band.write_json_atomic(summary_path, summary)
            raise
        summary["status"] = "complete"
        summary["analysis_completed"] = True
        summary["analysis_summary"] = str(
            series_dir / CONFIG_DIRNAME / "winding_summary.json"
        )
        summary["analysis_output_count"] = len(analysis_summary["outputs"])
        band.write_json_atomic(summary_path, summary)
    print(f"Unit-cell 2D workflow status: {summary['status']}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run cavity unit-cell 2D k scan and automatic analysis."
    )
    parser.add_argument(
        "--skip-analysis",
        action="store_true",
        help="Skip the default offline analysis stage for diagnostics only.",
    )
    return parser


def main(argv: list[str] | None = None) -> dict[str, object]:
    args = build_parser().parse_args(argv)
    return run_unit_cell_2d(skip_analysis=args.skip_analysis)


if __name__ == "__main__":
    main()
