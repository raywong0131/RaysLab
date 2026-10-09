"""Expensive COMSOL stage for data-driven unit-cell valley refinement.

This module is deliberately imported only after ``--prepare-only`` has exited.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Sequence

import matplotlib.tri as mtri
import numpy as np
import pandas as pd

from comsol_workflow.band_connector import load_hz_center
from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationConfig
from comsol_workflow.unit_cell_2d_valley_refinement import (
    K_PRECISION,
    SOURCE_BAND_LABEL,
    load_source_scan,
    merge_polarization_points,
    plot_comsol_refined_b4,
    plot_convergence_repair_sampling_plan,
)
from scripts.run_main import run_band_pair as band
from scripts.run_main import run_unit_cell_2d as unit2d


LOCK_FILENAME = "valley_refinement.lock"
MODEL_FILENAME = "unit_cell_2D_valley_refinement.mph"
REFINED_TABLE_FILENAME = "refined_polarization_points.parquet"


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    band.write_json_atomic(path, payload)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _plan_hash(plan: dict[str, object]) -> str:
    stable = {
        "source_hashes": plan["source_hashes"],
        "sampling_parameters": plan["sampling_parameters"],
        "round_1_points": [
            point for point in plan["points"] if int(point["round"]) == 1
        ],
    }
    payload = json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _point_from_row(row) -> dict[str, object]:
    qx = float(row.qx_over_G)
    qy = float(row.qy_over_G)
    return {
        "direction": "unit_cell_2d_valley_refinement",
        "kx": qx,
        "ky": qy,
        "qx_over_G": qx,
        "qy_over_G": qy,
        "kx_str": f"{qx:.{K_PRECISION}f}",
        "ky_str": f"{qy:.{K_PRECISION}f}",
        "k_norm": math.hypot(qx, qy),
        "solve_order": int(row.solve_order),
        "plan_source": str(row.plan_source),
        "source_id": int(row.source_id),
        "round": int(row.round),
        "recenter_index": int(getattr(row, "recenter_index", 0)),
        "stencil_index": int(row.stencil_index),
        "stencil_axis": str(getattr(row, "stencil_axis", "junction")),
        "axis_index": int(getattr(row, "axis_index", row.stencil_index)),
    }


def _canonical_coordinate(source, qx: float, qy: float) -> tuple[float, float]:
    method = getattr(source, "canonical_coordinate", None)
    if callable(method):
        return method(qx, qy)
    if int(getattr(source, "is_quarter", 0)):
        return abs(float(qx)), abs(float(qy))
    return float(qx), float(qy)


def _refined_point_dir(series_dir: Path, point: dict[str, object]) -> Path:
    return (
        series_dir
        / "01_results"
        / "refined_points"
        / f"kx={point['kx_str']}_ky={point['ky_str']}"
    )


def _cache_identity(parameters, source_hashes, plan_hash: str) -> dict[str, object]:
    return {
        "workflow": "unit_cell_2d_valley_refinement",
        "source_hashes": source_hashes,
        "sampling_plan_hash": plan_hash,
        "mesh_size": parameters.mesh_size,
        "eigenmode_pair_count": parameters.unit_cell_2d.eigenmode_pair_count,
        "expected_eigenvalue_count": 4,
        "eigenfrequency_shift": "c_const/1.55[um]",
        "geometry_definition_version": unit2d.GEOMETRY_DEFINITION_VERSION,
        "tracking_definition_version": unit2d.TRACKING_DEFINITION_VERSION,
        "field_extraction_version": unit2d.FIELD_EXTRACTION_VERSION,
        "normalization_kind": "planar_l2",
        "field_plane": "air",
        "isQuarter": parameters.unit_cell_2d.is_quarter,
        "sampling_domain": parameters.unit_cell_2d.sampling_domain,
        "symmetry_render_version": unit2d.SYMMETRY_RENDER_VERSION,
        "symmetry_gauge_version": unit2d.SYMMETRY_GAUGE_VERSION,
    }


def _detect_comsol_processes() -> list[str]:
    if os.name != "nt":
        return []
    try:
        result = subprocess.run(
            ["tasklist", "/FO", "CSV", "/NH"],
            capture_output=True,
            text=True,
            check=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    return [
        line
        for line in result.stdout.splitlines()
        if any(token in line.lower() for token in ("comsol", "mphserver"))
    ]


def _acquire_lock(series_dir: Path) -> Path:
    path = series_dir / "99_config" / LOCK_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError(f"Valley-refinement lock already exists: {path}") from exc
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump({"pid": os.getpid()}, handle)
    return path


def _restore_source_states(source) -> dict[str, dict[tuple[float, float], dict[str, object]]]:
    states = {SOURCE_BAND_LABEL: {}}
    for row in source.frame.itertuples(index=False):
        point_dir = (
            source.series_dir
            / "01_results"
            / "k_points"
            / f"kx={float(row.qx_over_G):.{K_PRECISION}f}_ky={float(row.qy_over_G):.{K_PRECISION}f}"
        )
        mode_idx = int(row.mode_idx)
        field_path = point_dir / "eigenmodes" / f"{mode_idx:02d}_Hz_center.parquet"
        if not field_path.is_file():
            raise FileNotFoundError(f"Missing coarse Hz predecessor: {field_path}")
        state = row._asdict()
        state["field"] = load_hz_center(point_dir, mode_idx)
        states[SOURCE_BAND_LABEL][
            (float(row.qx_over_G), float(row.qy_over_G))
        ] = state
    return states


def _load_cached_row(
    point_dir: Path, identity: dict[str, object]
) -> dict[str, object] | None:
    path = point_dir / unit2d.POINT_METADATA_FILENAME
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    rows = payload.get("rows")
    if payload.get("cache_identity") != identity or not isinstance(rows, list) or len(rows) != 1:
        return None
    row = rows[0]
    if row.get("match_status") not in {"matched", "repaired"}:
        return None
    mode_idx = int(row["mode_idx"])
    if not (point_dir / "eigenmodes" / f"{mode_idx:02d}_Hz_center.parquet").is_file():
        return None
    return row


def _solve_point(
    runner,
    point: dict[str, object],
    point_dir: Path,
    holes: Sequence[Sequence[Sequence[float]]],
    parameters,
    states_by_band,
    identity,
) -> tuple[dict[str, object], int]:
    settings = parameters.unit_cell_2d
    point = dict(point)
    point["is_comsol_solved"] = True
    point["sampling_domain"] = settings.sampling_domain
    solve_count = 0
    final_row = None
    for attempt, strategy in enumerate(("nearest", "nearest"), start=1):
        runner.build_and_run(
            band.A,
            holes,
            {"kx": float(point["kx"]), "ky": float(point["ky"])},
            mesh_auto_size=parameters.mesh_size,
        )
        solve_count += 1
        modes = band.persist_k_point_solution(runner, point_dir)
        selected = unit2d.select_point_modes(
            point,
            modes,
            point_dir,
            [SOURCE_BAND_LABEL],
            states_by_band,
            settings,
            predecessor_strategy=strategy,
        )
        match = selected.get(SOURCE_BAND_LABEL)
        if match is None or not unit2d._match_passes_ambiguity(match, settings):
            continue
        polarization = runner.compute_polarization_data(
            int(match["mode_idx"]), settings.field_plane
        )
        status = "matched" if attempt == 1 else "repaired"
        final_row = unit2d._matched_row(
            point, SOURCE_BAND_LABEL, match, polarization, status
        )
        final_row.update(
            {
                "plan_source": point["plan_source"],
                "source_id": point["source_id"],
                "round": point["round"],
                "recenter_index": point["recenter_index"],
                "stencil_index": point["stencil_index"],
                "stencil_axis": point["stencil_axis"],
                "axis_index": point["axis_index"],
            }
        )
        state = dict(final_row)
        state["field"] = match["field"]
        states_by_band[SOURCE_BAND_LABEL][unit2d.point_key(point)] = state
        break
    if final_row is None:
        final_row = unit2d._invalid_row(
            point, SOURCE_BAND_LABEL, "ambiguous_after_repair"
        )
        final_row.update(
            {
                "plan_source": point["plan_source"],
                "source_id": point["source_id"],
                "round": point["round"],
                "recenter_index": point["recenter_index"],
                "stencil_index": point["stencil_index"],
                "stencil_axis": point["stencil_axis"],
                "axis_index": point["axis_index"],
            }
        )
    _atomic_json(
        point_dir / unit2d.POINT_METADATA_FILENAME,
        {"cache_identity": identity, "rows": [final_row]},
    )
    return final_row, solve_count


def _next_round_points(
    completed: pd.DataFrame,
    previous_plan: pd.DataFrame,
    source,
    settings,
    round_number: int,
) -> pd.DataFrame:
    valid = completed.loc[completed.match_status.isin(["matched", "repaired"])].copy()
    valid["abs_cy_raw"] = np.hypot(valid.cy_raw_re, valid.cy_raw_im)
    existing = {
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(source.frame.qx_over_G, source.frame.qy_over_G)
    }
    existing.update(
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(previous_plan.qx_over_G, previous_plan.qy_over_G)
    )
    records = []
    ratio = float(settings.refinement_ratio) ** (round_number - 1)
    for (plan_source, source_id), group in valid.groupby(["plan_source", "source_id"]):
        if plan_source == "valley_cross_stencil":
            # Tangent samples provide longitudinal display coverage only.  If
            # they participate in the minimum search, an anchor is pulled
            # along the valley toward whichever tangent endpoint is lower.
            normal_group = group.loc[group.stencil_axis == "normal"]
            if normal_group.empty:
                continue
            center = normal_group.loc[normal_group.abs_cy_raw.idxmin()]
            previous_group = previous_plan.loc[
                (previous_plan.plan_source == plan_source)
                & (previous_plan.source_id == source_id)
            ]
            if previous_group.empty:
                continue
            normal = previous_group[["normal_x", "normal_y"]].dropna().iloc[0].to_numpy(dtype=float)
            tangent = previous_group[["tangent_x", "tangent_y"]].dropna().iloc[0].to_numpy(dtype=float)
            normal_offsets = np.linspace(
                -float(settings.normal_half_width_over_g) / ratio,
                float(settings.normal_half_width_over_g) / ratio,
                int(settings.normal_points_per_round),
            )
            tangent_offsets = np.linspace(
                -float(settings.tangent_half_width_over_g) / ratio,
                float(settings.tangent_half_width_over_g) / ratio,
                int(settings.tangent_points_per_round),
            )
            offsets_xy = [
                ("normal", index, offset * normal[0], offset * normal[1])
                for index, offset in enumerate(normal_offsets)
            ] + [
                ("tangent", index, offset * tangent[0], offset * tangent[1])
                for index, offset in enumerate(tangent_offsets)
            ]
        elif plan_source == "junction_grid":
            center = group.loc[group.abs_cy_raw.idxmin()]
            offsets = np.linspace(
                -float(settings.junction_half_width_over_g) / ratio,
                float(settings.junction_half_width_over_g) / ratio,
                int(settings.junction_points_per_axis),
            )
            junction_pairs = [(dx, dy) for dy in offsets for dx in offsets]
            offsets_xy = [
                ("junction", index, dx, dy)
                for index, (dx, dy) in enumerate(junction_pairs)
            ]
            normal = np.asarray([np.nan, np.nan])
            tangent = np.asarray([np.nan, np.nan])
        else:
            # Segment gap-fill points improve longitudinal real-node coverage
            # but do not define an adaptive minimum-search stencil.
            continue
        for stencil_index, (stencil_axis, axis_index, dx, dy) in enumerate(offsets_xy):
            qx, qy = _canonical_coordinate(
                source,
                float(center.qx_over_G + dx),
                float(center.qy_over_G + dy),
            )
            if abs(qx) > source.q_max + 1e-12 or abs(qy) > source.q_max + 1e-12:
                continue
            key = (f"{qx:.{K_PRECISION}f}", f"{qy:.{K_PRECISION}f}")
            if key in existing:
                continue
            existing.add(key)
            records.append(
                {
                    "plan_source": plan_source,
                    "source_id": int(source_id),
                    "round": round_number,
                    "recenter_index": 0,
                    "stencil_index": stencil_index,
                    "stencil_axis": stencil_axis,
                    "axis_index": axis_index,
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "center_qx_over_G": float(center.qx_over_G),
                    "center_qy_over_G": float(center.qy_over_G),
                    "offset_over_G": math.hypot(dx, dy),
                    "normal_x": float(normal[0]),
                    "normal_y": float(normal[1]),
                    "tangent_x": float(tangent[0]),
                    "tangent_y": float(tangent[1]),
                    "status": "planned",
                    "kx_str": key[0],
                    "ky_str": key[1],
                }
            )
    start = len(previous_plan)
    records.sort(key=lambda item: (math.hypot(item["qx_over_G"], item["qy_over_G"]), item["qy_over_G"], item["qx_over_G"]))
    for offset, record in enumerate(records):
        record["point_id"] = start + offset
        record["solve_order"] = start + offset
    return pd.DataFrame(records)


def _logical_stencil_center(
    completed: pd.DataFrame,
    source_frame: pd.DataFrame,
    plan_group: pd.DataFrame,
) -> pd.Series | None:
    """Return a solved center even when coordinate dedup omitted this batch row."""
    if plan_group.empty:
        return None
    qx = float(plan_group.iloc[0].center_qx_over_G)
    qy = float(plan_group.iloc[0].center_qy_over_G)
    for frame in (completed, source_frame):
        matches = frame.loc[
            np.isclose(
                frame.qx_over_G.astype(float), qx, atol=0.5 * 10**-K_PRECISION
            )
            & np.isclose(
                frame.qy_over_G.astype(float), qy, atol=0.5 * 10**-K_PRECISION
            )
        ]
        if matches.empty:
            continue
        row = matches.iloc[0].copy()
        row["abs_cy_raw"] = math.hypot(
            float(row.cy_raw_re), float(row.cy_raw_im)
        )
        return row
    return None


def _edge_recenter_points(
    completed: pd.DataFrame,
    previous_plan: pd.DataFrame,
    source,
    settings,
    round_number: int,
    recenter_index: int,
) -> pd.DataFrame:
    """Append windows only when a complete logical stencil has an edge minimum."""
    valid_completed = completed.loc[
        completed.match_status.isin(["matched", "repaired"])
    ].copy()
    valid_completed["abs_cy_raw"] = np.hypot(
        valid_completed.cy_raw_re, valid_completed.cy_raw_im
    )
    current = completed.loc[
        (completed["round"] == round_number)
        & (completed["recenter_index"] == recenter_index - 1)
        & completed.match_status.isin(["matched", "repaired"])
    ].copy()
    if current.empty:
        return pd.DataFrame()
    current["abs_cy_raw"] = np.hypot(current.cy_raw_re, current.cy_raw_im)
    existing = {
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(source.frame.qx_over_G, source.frame.qy_over_G)
    }
    existing.update(
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(previous_plan.qx_over_G, previous_plan.qy_over_G)
    )
    records = []
    ratio = float(settings.refinement_ratio) ** (round_number - 1)
    for (plan_source, source_id), group in current.groupby(["plan_source", "source_id"]):
        plan_group = previous_plan.loc[
            (previous_plan.plan_source == plan_source)
            & (previous_plan.source_id == source_id)
            & (previous_plan["round"] == round_number)
            & (previous_plan.recenter_index == recenter_index - 1)
        ]
        if plan_group.empty:
            continue
        logical_center = _logical_stencil_center(
            valid_completed, source.frame, plan_group
        )
        if plan_source == "valley_cross_stencil":
            search = group.loc[group.stencil_axis == "normal"].copy()
            if search.empty:
                continue
            search["_logical_center"] = False
            if logical_center is not None:
                center_key = (
                    f"{float(logical_center.qx_over_G):.{K_PRECISION}f}",
                    f"{float(logical_center.qy_over_G):.{K_PRECISION}f}",
                )
                search_keys = {
                    (
                        f"{float(qx):.{K_PRECISION}f}",
                        f"{float(qy):.{K_PRECISION}f}",
                    )
                    for qx, qy in zip(search.qx_over_G, search.qy_over_G)
                }
                if center_key not in search_keys:
                    center = logical_center.to_dict()
                    center.update(
                        {
                            "stencil_axis": "normal",
                            "axis_index": int(settings.normal_points_per_round) // 2,
                            "_logical_center": True,
                        }
                    )
                    search = pd.concat(
                        [search, pd.DataFrame([center])], ignore_index=True
                    )
            minimum = search.loc[search.abs_cy_raw.idxmin()]
            if bool(minimum.get("_logical_center", False)):
                continue
            if int(minimum.axis_index) not in {
                0,
                int(settings.normal_points_per_round) - 1,
            }:
                continue
            normal = plan_group[["normal_x", "normal_y"]].dropna().iloc[0].to_numpy(dtype=float)
            tangent = plan_group[["tangent_x", "tangent_y"]].dropna().iloc[0].to_numpy(dtype=float)
            normal_offsets = np.linspace(
                -float(settings.normal_half_width_over_g) / ratio,
                float(settings.normal_half_width_over_g) / ratio,
                int(settings.normal_points_per_round),
            )
            tangent_offsets = np.linspace(
                -float(settings.tangent_half_width_over_g) / ratio,
                float(settings.tangent_half_width_over_g) / ratio,
                int(settings.tangent_points_per_round),
            )
            offsets_xy = [
                ("normal", index, offset * normal[0], offset * normal[1])
                for index, offset in enumerate(normal_offsets)
            ] + [
                ("tangent", index, offset * tangent[0], offset * tangent[1])
                for index, offset in enumerate(tangent_offsets)
            ]
        elif plan_source == "junction_grid":
            width = int(settings.junction_points_per_axis)
            search = group.copy()
            search["_logical_center"] = False
            if logical_center is not None:
                center_key = (
                    f"{float(logical_center.qx_over_G):.{K_PRECISION}f}",
                    f"{float(logical_center.qy_over_G):.{K_PRECISION}f}",
                )
                search_keys = {
                    (
                        f"{float(qx):.{K_PRECISION}f}",
                        f"{float(qy):.{K_PRECISION}f}",
                    )
                    for qx, qy in zip(search.qx_over_G, search.qy_over_G)
                }
                if center_key not in search_keys:
                    center = logical_center.to_dict()
                    center.update(
                        {
                            "stencil_index": width * width // 2,
                            "_logical_center": True,
                        }
                    )
                    search = pd.concat(
                        [search, pd.DataFrame([center])], ignore_index=True
                    )
            minimum = search.loc[search.abs_cy_raw.idxmin()]
            if bool(minimum.get("_logical_center", False)):
                continue
            index = int(minimum.stencil_index)
            x_index, y_index = index % width, index // width
            if x_index not in {0, width - 1} and y_index not in {0, width - 1}:
                continue
            offsets = np.linspace(
                -float(settings.junction_half_width_over_g) / ratio,
                float(settings.junction_half_width_over_g) / ratio,
                width,
            )
            junction_pairs = [(dx, dy) for dy in offsets for dx in offsets]
            offsets_xy = [
                ("junction", index, dx, dy)
                for index, (dx, dy) in enumerate(junction_pairs)
            ]
            normal = np.asarray([np.nan, np.nan])
            tangent = np.asarray([np.nan, np.nan])
        else:
            continue
        for stencil_index, (stencil_axis, axis_index, dx, dy) in enumerate(offsets_xy):
            qx, qy = _canonical_coordinate(
                source,
                float(minimum.qx_over_G + dx),
                float(minimum.qy_over_G + dy),
            )
            if abs(qx) > source.q_max + 1e-12 or abs(qy) > source.q_max + 1e-12:
                continue
            key = (f"{qx:.{K_PRECISION}f}", f"{qy:.{K_PRECISION}f}")
            if key in existing:
                continue
            existing.add(key)
            records.append(
                {
                    "plan_source": plan_source,
                    "source_id": int(source_id),
                    "round": round_number,
                    "recenter_index": recenter_index,
                    "stencil_index": stencil_index,
                    "stencil_axis": stencil_axis,
                    "axis_index": axis_index,
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "center_qx_over_G": float(minimum.qx_over_G),
                    "center_qy_over_G": float(minimum.qy_over_G),
                    "offset_over_G": math.hypot(dx, dy),
                    "normal_x": float(normal[0]),
                    "normal_y": float(normal[1]),
                    "tangent_x": float(tangent[0]),
                    "tangent_y": float(tangent[1]),
                    "status": "planned_recenter",
                    "kx_str": key[0],
                    "ky_str": key[1],
                }
            )
    start = len(previous_plan)
    records.sort(key=lambda item: (math.hypot(item["qx_over_G"], item["qy_over_G"]), item["qy_over_G"], item["qx_over_G"]))
    for offset, record in enumerate(records):
        record["point_id"] = start + offset
        record["solve_order"] = start + offset
    return pd.DataFrame(records)


def build_convergence_repair_targets(
    completed: pd.DataFrame,
    previous_plan: pd.DataFrame,
    settings,
) -> pd.DataFrame:
    """Describe every logical window, including those needing no new coordinate."""
    valid = completed.loc[
        completed.match_status.isin(["matched", "repaired"])
        & (completed["round"] == 1)
        & (completed.recenter_index == 0)
    ].copy()
    valid["abs_cy_raw"] = np.hypot(valid.cy_raw_re, valid.cy_raw_im)
    records = []
    repair_recenter_index = int(previous_plan.recenter_index.max()) + 1
    for (plan_source, source_id), group in valid.groupby(
        ["plan_source", "source_id"]
    ):
        plan_group = previous_plan.loc[
            (previous_plan.plan_source == plan_source)
            & (previous_plan.source_id == source_id)
            & (previous_plan["round"] == 1)
            & (previous_plan.recenter_index == 0)
        ]
        if plan_group.empty:
            continue
        if plan_source == "valley_cross_stencil":
            normal_group = group.loc[group.stencil_axis == "normal"]
            if normal_group.empty:
                continue
            center = normal_group.loc[normal_group.abs_cy_raw.idxmin()]
            edge = int(center.axis_index) in {
                int(normal_group.axis_index.min()),
                int(normal_group.axis_index.max()),
            }
            ratio = 1.0 if edge else float(settings.refinement_ratio)
            repair_kind = "round1_normal_recenter" if edge else "round2_normal_refine"
            normal = (
                plan_group[["normal_x", "normal_y"]]
                .dropna()
                .iloc[0]
                .to_numpy(dtype=float)
            )
            tangent = (
                plan_group[["tangent_x", "tangent_y"]]
                .dropna()
                .iloc[0]
                .to_numpy(dtype=float)
            )
        elif plan_source == "junction_grid":
            center = group.loc[group.abs_cy_raw.idxmin()]
            width = int(settings.junction_points_per_axis)
            index = int(center.stencil_index)
            x_index, y_index = index % width, index // width
            edge = x_index in {0, width - 1} or y_index in {0, width - 1}
            ratio = 1.0 if edge else float(settings.refinement_ratio)
            repair_kind = "round1_junction_recenter" if edge else "round2_junction_refine"
            normal = np.asarray([np.nan, np.nan])
            tangent = np.asarray([np.nan, np.nan])
        else:
            continue
        records.append(
            {
                "plan_source": str(plan_source),
                "source_id": int(source_id),
                "round": 2,
                "recenter_index": repair_recenter_index,
                "repair_kind": repair_kind,
                "center_qx_over_G": float(center.qx_over_G),
                "center_qy_over_G": float(center.qy_over_G),
                "normal_x": float(normal[0]),
                "normal_y": float(normal[1]),
                "tangent_x": float(tangent[0]),
                "tangent_y": float(tangent[1]),
                "refinement_ratio_used": ratio,
            }
        )
    columns = [
        "plan_source",
        "source_id",
        "round",
        "recenter_index",
        "repair_kind",
        "center_qx_over_G",
        "center_qy_over_G",
        "normal_x",
        "normal_y",
        "tangent_x",
        "tangent_y",
        "refinement_ratio_used",
    ]
    return pd.DataFrame(records, columns=columns).sort_values(
        ["plan_source", "source_id"]
    ).reset_index(drop=True)


def build_convergence_repair_plan(
    completed: pd.DataFrame,
    previous_plan: pd.DataFrame,
    source,
    settings,
) -> pd.DataFrame:
    """Build one review-only normal-first correction pass from solved nodes."""
    targets = build_convergence_repair_targets(completed, previous_plan, settings)
    existing = {
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(source.frame.qx_over_G, source.frame.qy_over_G)
    }
    existing.update(
        (f"{float(qx):.{K_PRECISION}f}", f"{float(qy):.{K_PRECISION}f}")
        for qx, qy in zip(previous_plan.qx_over_G, previous_plan.qy_over_G)
    )
    records: list[dict[str, object]] = []
    for target in targets.itertuples(index=False):
        center_qx = float(target.center_qx_over_G)
        center_qy = float(target.center_qy_over_G)
        ratio = float(target.refinement_ratio_used)
        normal = np.asarray([target.normal_x, target.normal_y], dtype=float)
        tangent = np.asarray([target.tangent_x, target.tangent_y], dtype=float)
        if target.plan_source == "valley_cross_stencil":
            normal_offsets = np.linspace(
                -float(settings.normal_half_width_over_g) / ratio,
                float(settings.normal_half_width_over_g) / ratio,
                int(settings.normal_points_per_round),
            )
            tangent_offsets = np.linspace(
                -float(settings.tangent_half_width_over_g) / ratio,
                float(settings.tangent_half_width_over_g) / ratio,
                int(settings.tangent_points_per_round),
            )
            offsets_xy = [
                ("normal", index, offset * normal[0], offset * normal[1])
                for index, offset in enumerate(normal_offsets)
            ] + [
                ("tangent", index, offset * tangent[0], offset * tangent[1])
                for index, offset in enumerate(tangent_offsets)
            ]
        else:
            offsets = np.linspace(
                -float(settings.junction_half_width_over_g) / ratio,
                float(settings.junction_half_width_over_g) / ratio,
                int(settings.junction_points_per_axis),
            )
            offsets_xy = [
                ("junction", index, dx, dy)
                for index, (dx, dy) in enumerate(
                    (dx, dy) for dy in offsets for dx in offsets
                )
            ]
        for stencil_index, (stencil_axis, axis_index, dx, dy) in enumerate(
            offsets_xy
        ):
            qx, qy = _canonical_coordinate(
                source,
                center_qx + dx, center_qy + dy
            )
            if abs(qx) > source.q_max + 1e-12 or abs(qy) > source.q_max + 1e-12:
                continue
            key = (f"{qx:.{K_PRECISION}f}", f"{qy:.{K_PRECISION}f}")
            if key in existing:
                continue
            existing.add(key)
            records.append(
                {
                    "plan_source": target.plan_source,
                    "source_id": int(target.source_id),
                    "round": 2,
                    "recenter_index": int(target.recenter_index),
                    "repair_kind": target.repair_kind,
                    "stencil_index": stencil_index,
                    "stencil_axis": stencil_axis,
                    "axis_index": axis_index,
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "center_qx_over_G": center_qx,
                    "center_qy_over_G": center_qy,
                    "offset_over_G": math.hypot(dx, dy),
                    "normal_x": float(normal[0]),
                    "normal_y": float(normal[1]),
                    "tangent_x": float(tangent[0]),
                    "tangent_y": float(tangent[1]),
                    "status": "proposed_convergence_repair",
                    "kx_str": key[0],
                    "ky_str": key[1],
                }
            )
    records.sort(
        key=lambda item: (
            math.hypot(item["qx_over_G"], item["qy_over_G"]),
            item["qy_over_G"],
            item["qx_over_G"],
        )
    )
    start = len(previous_plan)
    for offset, record in enumerate(records):
        record["point_id"] = start + offset
        record["solve_order"] = start + offset
    result = pd.DataFrame(records)
    projected = len(previous_plan) + len(result)
    if projected > int(settings.maximum_new_points):
        raise RuntimeError(
            "Convergence repair would exceed maximum_new_points: "
            f"{projected} > {settings.maximum_new_points}"
        )
    return result


def prepare_convergence_repair(series_dir: Path, parameters) -> dict[str, object]:
    """Create a correction CSV/PNG without changing the approved solve plan."""
    series_dir = Path(series_dir).resolve()
    plan = json.loads(
        (series_dir / "99_config" / "sampling_plan.json").read_text(
            encoding="utf-8"
        )
    )
    source = load_source_scan(plan["source_series"])
    completed = pd.read_parquet(
        series_dir / "01_results" / REFINED_TABLE_FILENAME
    )
    previous_plan = pd.read_csv(series_dir / "80_logs" / "refinement_points.csv")
    targets = build_convergence_repair_targets(
        completed,
        previous_plan,
        parameters.unit_cell_2d.valley_refinement,
    )
    repair = build_convergence_repair_plan(
        completed,
        previous_plan,
        source,
        parameters.unit_cell_2d.valley_refinement,
    )
    csv_path = series_dir / "80_logs" / "convergence_repair_points.csv"
    repair.to_csv(csv_path, index=False)
    targets_path = series_dir / "80_logs" / "convergence_repair_targets.csv"
    targets.to_csv(targets_path, index=False)
    figure_path = (
        series_dir / "00_model" / "valley_convergence_repair_sampling_plan.png"
    )
    plot_convergence_repair_sampling_plan(
        source, completed, repair, figure_path
    )
    summary = {
        "status": "prepared_for_review",
        "source_series": str(source.series_dir),
        "existing_refined_point_count": len(previous_plan),
        "proposed_new_point_count": len(repair),
        "audited_target_count": len(targets),
        "projected_refined_point_count": len(previous_plan) + len(repair),
        "maximum_new_points": int(
            parameters.unit_cell_2d.valley_refinement.maximum_new_points
        ),
        "repair_points": str(csv_path),
        "repair_points_sha256": _sha256(csv_path),
        "repair_targets": str(targets_path),
        "repair_targets_sha256": _sha256(targets_path),
        "sampling_plan_figure": str(figure_path),
        "comsol_started": False,
    }
    _atomic_json(
        series_dir / "99_config" / "convergence_repair_plan.json", summary
    )
    return summary


def _coordinate_key(qx: object, qy: object) -> tuple[str, str]:
    return (
        f"{float(qx):.{K_PRECISION}f}",
        f"{float(qy):.{K_PRECISION}f}",
    )


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.stem + ".tmp" + path.suffix)
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _combine_approved_repair_plan(
    previous_plan: pd.DataFrame,
    repair: pd.DataFrame,
    source,
    settings,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate an immutable approved CSV and append only genuinely new nodes."""
    required = {
        "plan_source",
        "source_id",
        "round",
        "recenter_index",
        "stencil_index",
        "stencil_axis",
        "axis_index",
        "qx_over_G",
        "qy_over_G",
        "center_qx_over_G",
        "center_qy_over_G",
        "normal_x",
        "normal_y",
        "tangent_x",
        "tangent_y",
        "point_id",
        "solve_order",
    }
    missing = sorted(required.difference(repair.columns))
    if missing:
        raise ValueError(f"Approved repair CSV is missing columns: {missing}")
    repair = repair.copy()
    if int(getattr(source, "is_quarter", 0)):
        canonical = [
            _canonical_coordinate(source, qx, qy)
            for qx, qy in zip(repair.qx_over_G, repair.qy_over_G)
        ]
        repair["qx_over_G"] = [item[0] for item in canonical]
        repair["qy_over_G"] = [item[1] for item in canonical]
    repair["kx_str"] = [
        _coordinate_key(qx, qy)[0]
        for qx, qy in zip(repair.qx_over_G, repair.qy_over_G)
    ]
    repair["ky_str"] = [
        _coordinate_key(qx, qy)[1]
        for qx, qy in zip(repair.qx_over_G, repair.qy_over_G)
    ]
    repair_keys = list(zip(repair.kx_str, repair.ky_str))
    if len(set(repair_keys)) != len(repair_keys):
        raise ValueError("Approved repair CSV contains duplicate coordinates")
    source_keys = {
        _coordinate_key(qx, qy)
        for qx, qy in zip(source.frame.qx_over_G, source.frame.qy_over_G)
    }
    if source_keys.intersection(repair_keys):
        raise ValueError("Approved repair CSV overlaps the coarse COMSOL grid")
    previous = previous_plan.copy()
    previous["kx_str"] = [
        _coordinate_key(qx, qy)[0]
        for qx, qy in zip(previous.qx_over_G, previous.qy_over_G)
    ]
    previous["ky_str"] = [
        _coordinate_key(qx, qy)[1]
        for qx, qy in zip(previous.qx_over_G, previous.qy_over_G)
    ]
    previous_keys = set(zip(previous.kx_str, previous.ky_str))
    new_mask = [key not in previous_keys for key in repair_keys]
    already_appended = repair.loc[[not value for value in new_mask]]
    if not already_appended.empty:
        previous_by_key = {
            key: row
            for key, (_, row) in zip(
                zip(previous.kx_str, previous.ky_str), previous.iterrows()
            )
        }
        for row in already_appended.itertuples(index=False):
            existing = previous_by_key[(row.kx_str, row.ky_str)]
            if (
                str(existing.plan_source) != str(row.plan_source)
                or int(existing.source_id) != int(row.source_id)
            ):
                raise ValueError(
                    "An approved repair coordinate already exists with different "
                    f"provenance: {(row.kx_str, row.ky_str)}"
                )
    new_points = repair.loc[new_mask].copy()
    combined = pd.concat([previous, new_points], ignore_index=True, sort=False)
    combined_keys = list(zip(combined.kx_str, combined.ky_str))
    if len(set(combined_keys)) != len(combined_keys):
        raise ValueError("Combined refinement plan is not coordinate-unique")
    if len(combined) > int(settings.maximum_new_points):
        raise RuntimeError(
            "Approved repair would exceed maximum_new_points: "
            f"{len(combined)} > {settings.maximum_new_points}"
        )
    return combined, new_points


def audit_refinement_convergence(
    source,
    completed: pd.DataFrame,
    repair_plan: pd.DataFrame,
    settings,
) -> dict[str, object]:
    """Audit repaired logical stencils and the real-node triangulation."""
    valid = completed.loc[
        completed.match_status.isin(["matched", "repaired"])
    ].copy()
    valid["abs_cy_raw"] = np.hypot(valid.cy_raw_re, valid.cy_raw_im)
    node_frame = pd.concat([source.frame, valid], ignore_index=True, sort=False)
    node_frame["coordinate_key"] = [
        ":".join(_coordinate_key(qx, qy))
        for qx, qy in zip(node_frame.qx_over_G, node_frame.qy_over_G)
    ]
    node_frame["abs_cy_raw"] = np.hypot(
        node_frame.cy_raw_re, node_frame.cy_raw_im
    )
    node_by_key = {
        str(row.coordinate_key): row for row in node_frame.itertuples(index=False)
    }
    group_rows: list[dict[str, object]] = []
    for (plan_source, source_id), group in repair_plan.groupby(
        ["plan_source", "source_id"]
    ):
        center_qx = float(group.center_qx_over_G.iloc[0])
        center_qy = float(group.center_qy_over_G.iloc[0])
        repair_kind = str(group.repair_kind.iloc[0])
        ratio = (
            1.0
            if "recenter" in repair_kind
            else float(settings.refinement_ratio)
        )
        expected: list[tuple[str, str]] = []
        edge = True
        minimum_key = None
        if plan_source == "valley_cross_stencil":
            direction_rows = group[["normal_x", "normal_y"]].dropna()
            if direction_rows.empty:
                group_rows.append(
                    {
                        "plan_source": plan_source,
                        "source_id": int(source_id),
                        "repair_kind": repair_kind,
                        "complete": False,
                        "edge_minimum": True,
                        "reason": "missing normal direction",
                    }
                )
                continue
            normal = direction_rows.iloc[0].to_numpy(dtype=float)
            offsets = np.linspace(
                -float(settings.normal_half_width_over_g) / ratio,
                float(settings.normal_half_width_over_g) / ratio,
                int(settings.normal_points_per_round),
            )
            expected = [
                _coordinate_key(
                    *_canonical_coordinate(
                        source,
                        center_qx + offset * normal[0],
                        center_qy + offset * normal[1],
                    )
                )
                for offset in offsets
            ]
        elif plan_source == "junction_grid":
            offsets = np.linspace(
                -float(settings.junction_half_width_over_g) / ratio,
                float(settings.junction_half_width_over_g) / ratio,
                int(settings.junction_points_per_axis),
            )
            expected = [
                _coordinate_key(
                    *_canonical_coordinate(source, center_qx + dx, center_qy + dy)
                )
                for dy in offsets
                for dx in offsets
            ]
        expected = [
            key
            for key in expected
            if abs(float(key[0])) <= source.q_max + 10**-K_PRECISION
            and abs(float(key[1])) <= source.q_max + 10**-K_PRECISION
        ]
        found = [node_by_key.get(":".join(key)) for key in expected]
        complete = bool(expected) and all(row is not None for row in found)
        if complete:
            magnitudes = np.asarray(
                [float(row.abs_cy_raw) for row in found], dtype=float
            )
            minimum_index = int(np.argmin(magnitudes))
            minimum_key = expected[minimum_index]
            if plan_source == "valley_cross_stencil":
                edge = minimum_index in {0, len(expected) - 1}
            else:
                width = int(settings.junction_points_per_axis)
                x_index = minimum_index % width
                y_index = minimum_index // width
                edge = x_index in {0, width - 1} or y_index in {0, width - 1}
        group_rows.append(
            {
                "plan_source": str(plan_source),
                "source_id": int(source_id),
                "repair_kind": repair_kind,
                "expected_node_count": len(expected),
                "found_node_count": sum(row is not None for row in found),
                "complete": complete,
                "edge_minimum": bool(edge),
                "minimum_coordinate": (
                    list(minimum_key) if minimum_key is not None else None
                ),
            }
        )
    merged = merge_polarization_points(source.frame, valid)
    triangulation_ok = True
    triangulation_error = None
    triangle_count = 0
    minimum_triangle_area = None
    maximum_edge_length = None
    long_edge_limit = float(source.coarse_step) * math.sqrt(2.0) * 1.001
    try:
        x = merged.qx_over_G.to_numpy(dtype=float)
        y = merged.qy_over_G.to_numpy(dtype=float)
        triangulation = mtri.Triangulation(x, y)
        triangles = triangulation.triangles
        triangle_count = int(len(triangles))
        p0 = np.column_stack([x[triangles[:, 0]], y[triangles[:, 0]]])
        p1 = np.column_stack([x[triangles[:, 1]], y[triangles[:, 1]]])
        p2 = np.column_stack([x[triangles[:, 2]], y[triangles[:, 2]]])
        areas = 0.5 * np.abs(
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        )
        edges = np.concatenate(
            [
                np.linalg.norm(p1 - p0, axis=1),
                np.linalg.norm(p2 - p1, axis=1),
                np.linalg.norm(p0 - p2, axis=1),
            ]
        )
        minimum_triangle_area = float(np.min(areas))
        maximum_edge_length = float(np.max(edges))
        triangulation_ok = bool(
            minimum_triangle_area > 1e-16
            and maximum_edge_length <= long_edge_limit
        )
    except (RuntimeError, ValueError) as exc:
        triangulation_ok = False
        triangulation_error = str(exc)
    incomplete_groups = sum(not row["complete"] for row in group_rows)
    edge_groups = sum(row["edge_minimum"] for row in group_rows)
    invalid_modes = int(len(completed) - len(valid))
    point_limit_ok = len(completed) <= int(settings.maximum_new_points)
    converged = bool(
        group_rows
        and incomplete_groups == 0
        and edge_groups == 0
        and invalid_modes == 0
        and point_limit_ok
        and triangulation_ok
    )
    return {
        "status": "converged" if converged else "incomplete",
        "converged": converged,
        "refined_point_count": int(len(completed)),
        "valid_mode_count": int(len(valid)),
        "invalid_mode_count": invalid_modes,
        "audited_group_count": len(group_rows),
        "incomplete_stencil_count": incomplete_groups,
        "edge_minimum_count": edge_groups,
        "point_limit": int(settings.maximum_new_points),
        "point_limit_ok": point_limit_ok,
        "triangulation": {
            "ok": triangulation_ok,
            "triangle_count": triangle_count,
            "minimum_triangle_area": minimum_triangle_area,
            "maximum_edge_length": maximum_edge_length,
            "maximum_allowed_edge_length": long_edge_limit,
            "error": triangulation_error,
        },
        "groups": group_rows,
    }


def solve_approved_convergence_repair(
    series_dir: Path,
    parameters,
    *,
    resume: bool,
    allow_concurrent_comsol: bool,
) -> dict[str, object]:
    """Solve only the immutable, separately reviewed convergence-repair CSV."""
    series_dir = Path(series_dir).resolve()
    repair_summary_path = (
        series_dir / "99_config" / "convergence_repair_plan.json"
    )
    repair_summary = json.loads(repair_summary_path.read_text(encoding="utf-8"))
    repair_path = Path(repair_summary["repair_points"]).resolve()
    expected_hash = str(repair_summary.get("repair_points_sha256", ""))
    if not expected_hash or _sha256(repair_path) != expected_hash:
        raise RuntimeError("Approved convergence-repair CSV hash mismatch")
    targets_path = Path(repair_summary["repair_targets"]).resolve()
    expected_targets_hash = str(
        repair_summary.get("repair_targets_sha256", "")
    )
    if (
        not expected_targets_hash
        or _sha256(targets_path) != expected_targets_hash
    ):
        raise RuntimeError("Approved convergence-repair target hash mismatch")
    if repair_summary.get("status") not in {
        "prepared_for_review",
        "solving",
        "incomplete",
    }:
        raise RuntimeError(
            "Convergence-repair plan is not in an executable state: "
            f"{repair_summary.get('status')}"
        )
    plan_path = series_dir / "99_config" / "sampling_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    source = load_source_scan(plan["source_series"])
    if source.source_hashes != plan["source_hashes"]:
        raise RuntimeError("Source hashes changed after sampling-plan approval")
    detected = _detect_comsol_processes()
    if detected and not allow_concurrent_comsol:
        raise RuntimeError(
            "Another COMSOL process is visible; rerun only after it exits or pass "
            "--allow-concurrent-comsol explicitly"
        )
    points_path = series_dir / "80_logs" / "refinement_points.csv"
    previous_plan = pd.read_csv(points_path)
    repair = pd.read_csv(repair_path)
    repair_targets = pd.read_csv(targets_path)
    points, new_points = _combine_approved_repair_plan(
        previous_plan,
        repair,
        source,
        parameters.unit_cell_2d.valley_refinement,
    )
    refined_path = series_dir / "01_results" / REFINED_TABLE_FILENAME
    completed = pd.read_parquet(refined_path)
    completed_keys = {
        _coordinate_key(qx, qy)
        for qx, qy in zip(completed.qx_over_G, completed.qy_over_G)
    }
    if not resume and completed_keys.intersection(
        _coordinate_key(qx, qy)
        for qx, qy in zip(repair.qx_over_G, repair.qy_over_G)
    ):
        raise RuntimeError("Repair results already exist but --no-resume was requested")
    if not new_points.empty:
        _atomic_csv(points, points_path)
        plan["points"] = (
            points.astype(object)
            .where(pd.notna(points), None)
            .to_dict(orient="records")
        )
        plan["sampling_plan_version"] = "convergence-repair-1"
        _atomic_json(plan_path, plan)
    lock = _acquire_lock(series_dir)
    repair_summary.update(
        {
            "status": "solving",
            "approved_csv_verified": True,
            "approved_csv_sha256": expected_hash,
            "appended_new_point_count": int(len(new_points)),
            "comsol_started": False,
            "pid": os.getpid(),
        }
    )
    _atomic_json(repair_summary_path, repair_summary)
    simulation_config = SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift="c_const/1.55[um]",
        eigenmode_count=2,
        mesh_auto_size=parameters.mesh_size,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )
    _unit_params, _hole_params, holes, _clearance = unit2d._prepare_geometry(
        parameters, series_dir / "00_model"
    )
    hole_vertices = [hole.tolist() for hole in holes]
    states = _restore_source_states(source)
    for row in completed.itertuples(index=False):
        point = {
            "kx_str": _coordinate_key(row.qx_over_G, row.qy_over_G)[0],
            "ky_str": _coordinate_key(row.qx_over_G, row.qy_over_G)[1],
        }
        unit2d._restore_cached_states(
            [row._asdict()], _refined_point_dir(series_dir, point), states
        )
    rows = completed.to_dict(orient="records")
    runner = None
    cache_hits = 0
    solve_count = 0
    unresolved = 0
    identity = _cache_identity(parameters, source.source_hashes, _plan_hash(plan))
    try:
        for row in repair.itertuples(index=False):
            point = _point_from_row(row)
            key = (point["kx_str"], point["ky_str"])
            if key in completed_keys:
                cache_hits += 1
                continue
            point_dir = _refined_point_dir(series_dir, point)
            cached = _load_cached_row(point_dir, identity) if resume else None
            if cached is not None:
                result = cached
                cache_hits += 1
                unit2d._restore_cached_states([cached], point_dir, states)
            else:
                if runner is None:
                    runner = band.ReusableSimulationRun(simulation_config)
                    repair_summary["comsol_started"] = True
                    _atomic_json(repair_summary_path, repair_summary)
                result, point_solves = _solve_point(
                    runner,
                    point,
                    point_dir,
                    hole_vertices,
                    parameters,
                    states,
                    identity,
                )
                solve_count += point_solves
            rows.append(result)
            completed_keys.add(key)
            if result["match_status"] == "ambiguous_after_repair":
                unresolved += 1
            completed = pd.DataFrame(rows)
            _atomic_parquet(completed, refined_path)
        if runner is not None:
            unit2d._save_runner_model(
                runner, series_dir / "00_model" / MODEL_FILENAME
            )
    except Exception as exc:
        repair_summary.update(
            {"status": "incomplete", "error": f"{type(exc).__name__}: {exc}"}
        )
        _atomic_json(repair_summary_path, repair_summary)
        raise
    finally:
        if runner is not None:
            runner.clear()
        lock.unlink(missing_ok=True)
    completed = pd.DataFrame(rows)
    completed.to_csv(
        series_dir / "80_logs" / "band_tracking_refined.csv", index=False
    )
    convergence = audit_refinement_convergence(
        source,
        completed,
        repair_targets,
        parameters.unit_cell_2d.valley_refinement,
    )
    convergence["repair_points_sha256"] = expected_hash
    convergence["repair_targets_sha256"] = expected_targets_hash
    convergence_path = series_dir / "99_config" / "convergence_summary.json"
    _atomic_json(convergence_path, convergence)
    repair_summary.update(
        {
            "status": "complete" if convergence["converged"] else "incomplete",
            "cache_hits": cache_hits,
            "solve_count": solve_count,
            "unresolved_mode_count": unresolved,
            "convergence_summary": str(convergence_path),
        }
    )
    repair_summary.pop("error", None)
    _atomic_json(repair_summary_path, repair_summary)
    run_summary_path = series_dir / "99_config" / "run_summary.json"
    run_summary = json.loads(run_summary_path.read_text(encoding="utf-8"))
    run_summary.update(
        {
            "status": "solved_pending_plot" if convergence["converged"] else "incomplete",
            "planned_point_count": len(points),
            "valid_point_count": int(
                completed.match_status.isin(["matched", "repaired"]).sum()
            ),
            "unresolved_point_count": unresolved,
            "repair_cache_hits": cache_hits,
            "repair_solve_count": solve_count,
            "convergence_summary": str(convergence_path),
        }
    )
    _atomic_json(run_summary_path, run_summary)
    return repair_summary


def solve_valley_refinement(
    series_dir: Path,
    parameters,
    *,
    resume: bool,
    allow_concurrent_comsol: bool,
) -> dict[str, object]:
    """Solve all planned points, appending later rounds from real minima only."""
    series_dir = Path(series_dir).resolve()
    plan_path = series_dir / "99_config" / "sampling_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    source = load_source_scan(plan["source_series"])
    if source.source_hashes != plan["source_hashes"]:
        raise RuntimeError("Source hashes changed after sampling-plan approval")
    detected = _detect_comsol_processes()
    if detected and not allow_concurrent_comsol:
        raise RuntimeError(
            "Another COMSOL process is visible; rerun only after it exits or pass "
            "--allow-concurrent-comsol explicitly"
        )
    lock = _acquire_lock(series_dir)
    summary_path = series_dir / "99_config" / "run_summary.json"
    points_path = series_dir / "80_logs" / "refinement_points.csv"
    points = pd.read_csv(points_path)
    _atomic_json(summary_path, {"status": "solving", "output_dir": str(series_dir), "pid": os.getpid()})
    simulation_config = SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift="c_const/1.55[um]",
        eigenmode_count=2,
        mesh_auto_size=parameters.mesh_size,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )
    _unit_params, _hole_params, holes, _clearance = unit2d._prepare_geometry(
        parameters, series_dir / "00_model"
    )
    hole_vertices = [hole.tolist() for hole in holes]
    states = _restore_source_states(source)
    rows: list[dict[str, object]] = []
    runner = None
    cache_hits = 0
    solve_count = 0
    unresolved = 0

    def solve_batch(batch: pd.DataFrame, identity: dict[str, object]) -> None:
        nonlocal runner, cache_hits, solve_count, unresolved
        for row in batch.itertuples(index=False):
            point = _point_from_row(row)
            point_dir = _refined_point_dir(series_dir, point)
            cached = _load_cached_row(point_dir, identity) if resume else None
            if cached is not None:
                cache_hits += 1
                rows.append(cached)
                unit2d._restore_cached_states([cached], point_dir, states)
                continue
            if runner is None:
                runner = band.ReusableSimulationRun(simulation_config)
            solved, point_solves = _solve_point(
                runner, point, point_dir, hole_vertices, parameters, states, identity
            )
            rows.append(solved)
            solve_count += point_solves
            if solved["match_status"] == "ambiguous_after_repair":
                unresolved += 1

    try:
        for round_number in range(1, int(parameters.unit_cell_2d.valley_refinement.refinement_rounds) + 1):
            round_points = points.loc[points["round"] == round_number]
            if round_number > 1 and round_points.empty:
                completed = pd.DataFrame(rows)
                appended = _next_round_points(
                    completed, points, source,
                    parameters.unit_cell_2d.valley_refinement,
                    round_number,
                )
                if not appended.empty:
                    points = pd.concat([points, appended], ignore_index=True)
                    if len(points) > parameters.unit_cell_2d.valley_refinement.maximum_new_points:
                        raise RuntimeError("Adaptive sampling exceeded maximum_new_points")
                    points.to_csv(points_path, index=False)
                    plan["points"] = points.to_dict(orient="records")
                    plan["sampling_plan_version"] = round_number
                    _atomic_json(plan_path, plan)
                round_points = appended
            identity = _cache_identity(parameters, source.source_hashes, _plan_hash(plan))
            solve_batch(round_points, identity)
            for recenter_index in range(
                1,
                int(parameters.unit_cell_2d.valley_refinement.max_edge_recenters) + 1,
            ):
                recentered = _edge_recenter_points(
                    pd.DataFrame(rows),
                    points,
                    source,
                    parameters.unit_cell_2d.valley_refinement,
                    round_number,
                    recenter_index,
                )
                if recentered.empty:
                    break
                points = pd.concat([points, recentered], ignore_index=True)
                if len(points) > parameters.unit_cell_2d.valley_refinement.maximum_new_points:
                    raise RuntimeError("Adaptive recentering exceeded maximum_new_points")
                points.to_csv(points_path, index=False)
                plan["points"] = points.to_dict(orient="records")
                plan["sampling_plan_version"] = (
                    f"round-{round_number}-recenter-{recenter_index}"
                )
                _atomic_json(plan_path, plan)
                solve_batch(recentered, identity)
            pd.DataFrame(rows).to_parquet(
                series_dir / "01_results" / REFINED_TABLE_FILENAME, index=False
            )
        if runner is not None:
            unit2d._save_runner_model(runner, series_dir / "00_model" / MODEL_FILENAME)
    finally:
        if runner is not None:
            runner.clear()
        lock.unlink(missing_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(series_dir / "80_logs" / "band_tracking_refined.csv", index=False)
    status = "solved" if unresolved == 0 else "incomplete"
    summary = {
        "status": status,
        "output_dir": str(series_dir),
        "planned_point_count": len(points),
        "valid_point_count": int(frame.match_status.isin(["matched", "repaired"]).sum()),
        "unresolved_point_count": unresolved,
        "cache_hits": cache_hits,
        "solve_count": solve_count,
        "comsol_started": runner is not None,
    }
    _atomic_json(summary_path, summary)
    if unresolved:
        raise RuntimeError(f"{unresolved} refined points remain ambiguous; final B4 is blocked")
    return summary


def plot_completed_refinement(series_dir: Path, parameters) -> Path:
    series_dir = Path(series_dir).resolve()
    plan = json.loads((series_dir / "99_config" / "sampling_plan.json").read_text(encoding="utf-8"))
    source = load_source_scan(plan["source_series"])
    refined_path = series_dir / "01_results" / REFINED_TABLE_FILENAME
    if not refined_path.is_file():
        raise FileNotFoundError(f"No completed refined COMSOL table: {refined_path}")
    refined = pd.read_parquet(refined_path)
    invalid = refined.loc[~refined.match_status.isin(["matched", "repaired"])]
    if not invalid.empty:
        raise RuntimeError("Refined point table contains unresolved modes; refusing final B4")
    merged = merge_polarization_points(source.frame, refined)
    merged_path = series_dir / "01_results" / "merged_polarization_points.parquet"
    merged.to_parquet(merged_path, index=False)
    output = series_dir / "10_overview" / "B4_cy_magnitude_log_comsol_refined.png"
    plot_comsol_refined_b4(
        merged,
        output,
        q_max=source.q_max,
        display_points_per_axis=parameters.unit_cell_2d.valley_refinement.display_points_per_axis,
        is_quarter=source.is_quarter,
    )
    summary_path = series_dir / "99_config" / "run_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    convergence_path = series_dir / "99_config" / "convergence_summary.json"
    if convergence_path.is_file():
        convergence = json.loads(convergence_path.read_text(encoding="utf-8"))
    else:
        convergence = {
            "status": "missing",
            "converged": False,
            "reason": "No convergence audit was produced by the solve stage",
        }
        _atomic_json(convergence_path, convergence)
    converged = bool(convergence.get("converged", False))
    summary.update(
        {
            "status": "complete" if converged else "incomplete",
            "merged_points": len(merged),
            "convergence_summary": str(convergence_path),
        }
    )
    if converged:
        summary["final_figure"] = str(output)
        summary.pop("provisional_figure", None)
    else:
        summary["provisional_figure"] = str(output)
        summary.pop("final_figure", None)
    _atomic_json(summary_path, summary)
    return output
