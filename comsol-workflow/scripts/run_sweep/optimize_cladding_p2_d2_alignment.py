#!/usr/bin/env python3
"""Quickly align cladding p2@Gamma to a fixed cavity with a d1/d2 soft target."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import (  # noqa: E402
    ACTIVE_PARAMETERS,
    CellParameters,
    UNIT_CELL_OUTPUT_ROOT,
)
from scripts.run_sweep import cladding_bandgap_alignment as shared_scan  # noqa: E402


FIXED_CAVITY = CellParameters(b0_nm=245.0, eta=0.960, zeta=1.156)
START_CLADDING = CellParameters(b0_nm=242.0, eta=0.980, zeta=0.930)
MESH_AUTO_SIZE = 5
EIGENMODE_COUNT = ACTIVE_PARAMETERS.unit_cell_eigenmode_count
EIGENFREQUENCY_SHIFT = "c_const/1.55[um]"
P2_ALIGNMENT_TOLERANCE_THZ = 0.02

ZETA_BOUNDS = (0.890, 0.970)
ZETA_COARSE_STEP = 0.004
ZETA_FINE_STEP = 0.001
ETA_BOUNDS = (0.940, 0.999)
ETA_COARSE_STEP = 0.004
ETA_FINE_STEP = 0.001
B0_BOUNDS_NM = (238.0, 246.0)
B0_COARSE_STEP_NM = 0.5
B0_FINE_STEP_NM = 0.1
MAX_DIRECTION_STEPS = 10

GAMMA_POINT = shared_scan.GAMMA_POINT
OUT_DIR = UNIT_CELL_OUTPUT_ROOT / (
    "cladding_p2_d2_alignment_"
    f"{FIXED_CAVITY.label('cav')}_"
    f"start-{START_CLADDING.label('clad')}_mesh{MESH_AUTO_SIZE}"
)


def make_cell_params(cell: CellParameters, name: str) -> dict[str, float | str]:
    """Use the shared compact-cell conversion used by the primary workflows."""
    return cell.to_fourier_params(name)


def simulation_config() -> band.SimulationConfig:
    """Return the task-fixed mesh5 unit-cell solver settings."""
    return band.SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift=EIGENFREQUENCY_SHIFT,
        eigenmode_count=EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )


def gamma_point_is_complete(
    root: Path,
    config: band.SimulationConfig | None = None,
) -> bool:
    """Check the task-local Gamma table and fields for every analyzed valid mode."""
    point_dir = band.k_point_dir(Path(root), GAMMA_POINT)
    return band.k_point_is_complete(point_dir, config or simulation_config())


def run_gamma_point(
    root: Path,
    holes: list[list[list[float]]],
) -> pd.DataFrame:
    """Solve Gamma and analyze only the frequency window used for p/d selection."""
    point_dir = band.k_point_dir(Path(root), GAMMA_POINT)
    csv_path = point_dir / "eigenfrequencies.csv"
    config = simulation_config()
    if gamma_point_is_complete(root, config):
        return pd.read_csv(csv_path)

    eigenmodes_dir = point_dir / "eigenmodes"
    eigenmodes_dir.mkdir(parents=True, exist_ok=True)
    with band.SimulationRun(config) as sim_run:
        sim_run.build_and_run(
            band.A,
            holes,
            {"kx": 0.0, "ky": 0.0},
        )
        frame = pd.DataFrame(
            sim_run.get_eigenfrequencies(),
            columns=["re", "im", "q"],
        )
        rows = []
        for mode_idx, mode in frame.iterrows():
            in_window = bool(
                band.FREQUENCY_MIN_THZ
                <= float(mode["re"])
                <= band.FREQUENCY_MAX_THZ
            )
            row: dict[str, object] = {
                "in_analysis_window": in_window,
                "is_valid": False,
            }
            if in_window:
                valid, peaks = band._mode_is_valid(sim_run, int(mode_idx))
                row.update(peaks)
                row["is_valid"] = bool(valid)
                if valid:
                    coordinates_re, hz_re = sim_run.get_2d_fields(
                        int(mode_idx), "ewfd.Hz", "center"
                    )
                    coordinates_im, hz_im = sim_run.get_2d_fields(
                        int(mode_idx), "ewfd.Hz*(-i)", "center"
                    )
                    hz_im_aligned = band.interpolate_field(
                        coordinates_im,
                        hz_im,
                        coordinates_re,
                    )
                    pd.DataFrame(
                        {
                            "x": coordinates_re[:, 0],
                            "y": coordinates_re[:, 1],
                            "re": np.real(hz_re),
                            "im": np.real(hz_im_aligned),
                        }
                    ).to_parquet(
                        eigenmodes_dir / f"{int(mode_idx):02d}_Hz_center.parquet"
                    )
                    row.update(band.load_mode_composition(point_dir, int(mode_idx)))
            rows.append(row)

    analysis = pd.DataFrame(rows, index=frame.index)
    frame = pd.concat([frame, analysis], axis=1)
    for column in (
        "s_weight",
        "p_weight",
        "d_weight",
        "f_weight",
        "px_weight",
        "py_weight",
        "dx_weight",
        "dy_weight",
        "dominant_subspace",
        "dominant_mode",
    ):
        if column not in frame:
            frame[column] = np.nan
    band.write_csv_atomic(csv_path, frame)
    band.write_json_atomic(
        point_dir / band.POINT_SOLVER_CONFIG_FILENAME,
        band.point_solver_identity(config),
    )
    return frame


def run_gamma_case(
    root: Path,
    cell: CellParameters,
    name: str,
) -> tuple[pd.DataFrame, dict[str, int], dict[str, object]]:
    """Run or reuse one Gamma point and identify its two p and two d modes."""
    params = make_cell_params(cell, name)
    hole_params = band.get_hole_params(params)
    _hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            "Unit-cell geometry violates clearance: "
            f"min_dist={minimum_clearance:g}, d={band.D:g}"
        )
    frame = run_gamma_point(Path(root), [hole.tolist() for hole in holes])
    selected = band.select_gamma_bands(frame)
    chosen = frame.loc[[selected[label] for label in ("p1", "p2", "d1", "d2")]]
    if not chosen["is_valid"].astype(bool).all():
        raise RuntimeError("Selected Gamma p/d modes include an invalid mode")
    return frame, selected, {
        "params": params,
        "hole_params": hole_params,
        "triangles": holes,
        "minimum_clearance": minimum_clearance,
        "selected_gamma_modes": selected,
    }


def _selected_frequency(
    frame: pd.DataFrame,
    selected: dict[str, int],
    label: str,
) -> float:
    value = float(frame.loc[int(selected[label]), "re"])
    if not math.isfinite(value):
        raise ValueError(f"Selected {label}@Gamma frequency is not finite")
    return value


def load_cavity_reference() -> dict[str, object]:
    """Compute or reuse the fixed cavity p2/d1 Gamma reference."""
    root = OUT_DIR / "cavity_reference"
    frame, selected, metadata = run_gamma_case(
        root,
        FIXED_CAVITY,
        "fixed_cavity_p2_d1_reference",
    )
    reference = {
        "cavity": FIXED_CAVITY.to_dict(),
        "mesh_auto_size": MESH_AUTO_SIZE,
        "eigenmode_count": EIGENMODE_COUNT,
        "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
        "cavity_p2_frequency_thz": _selected_frequency(frame, selected, "p2"),
        "cavity_d1_frequency_thz": _selected_frequency(frame, selected, "d1"),
        **metadata,
    }
    band.write_json_atomic(root / "reference.json", reference)
    return reference


def candidate_key(
    b0_nm: float,
    eta: float,
    zeta: float,
) -> tuple[float, float, float]:
    return round(float(b0_nm), 1), round(float(eta), 3), round(float(zeta), 3)


def candidate_dir(b0_nm: float, eta: float, zeta: float) -> Path:
    b0_nm, eta, zeta = candidate_key(b0_nm, eta, zeta)
    return OUT_DIR / "candidates" / f"clad({b0_nm:.1f}-{eta:.3f}-{zeta:.3f})"


def summarize_candidate(
    cell: CellParameters,
    stage: str,
    frame: pd.DataFrame,
    selected: dict[str, int],
    reference: dict[str, object],
) -> dict[str, object]:
    """Measure the primary p2 alignment and secondary d1/d2 soft error."""
    p2_frequency = _selected_frequency(frame, selected, "p2")
    d2_frequency = _selected_frequency(frame, selected, "d2")
    p2_signed_error = p2_frequency - float(reference["cavity_p2_frequency_thz"])
    d_soft_signed_error = d2_frequency - float(reference["cavity_d1_frequency_thz"])
    p2_error = abs(p2_signed_error)
    d_soft_error = abs(d_soft_signed_error)
    aligned = p2_error <= P2_ALIGNMENT_TOLERANCE_THZ + 1.0e-12
    return {
        **cell.to_dict(),
        "stage": str(stage),
        "status": "p2_aligned" if aligned else "p2_unaligned",
        "feasible": True,
        "p2_aligned": bool(aligned),
        "cavity_p2_frequency_thz": float(reference["cavity_p2_frequency_thz"]),
        "cladding_p2_frequency_thz": p2_frequency,
        "p2_signed_error_thz": p2_signed_error,
        "p2_error_thz": p2_error,
        "cavity_d1_frequency_thz": float(reference["cavity_d1_frequency_thz"]),
        "cladding_d2_frequency_thz": d2_frequency,
        "d1_d2_signed_error_thz": d_soft_signed_error,
        "d1_d2_soft_error_thz": d_soft_error,
        "p2_mode_idx": int(selected["p2"]),
        "d2_mode_idx": int(selected["d2"]),
    }


def failed_candidate(
    cell: CellParameters,
    stage: str,
    reference: dict[str, object],
    error: Exception,
) -> dict[str, object]:
    return {
        **cell.to_dict(),
        "stage": str(stage),
        "status": "candidate_error",
        "feasible": False,
        "p2_aligned": False,
        "cavity_p2_frequency_thz": float(reference["cavity_p2_frequency_thz"]),
        "cladding_p2_frequency_thz": np.nan,
        "p2_signed_error_thz": np.nan,
        "p2_error_thz": np.nan,
        "cavity_d1_frequency_thz": float(reference["cavity_d1_frequency_thz"]),
        "cladding_d2_frequency_thz": np.nan,
        "d1_d2_signed_error_thz": np.nan,
        "d1_d2_soft_error_thz": np.nan,
        "p2_mode_idx": np.nan,
        "d2_mode_idx": np.nan,
        "error": f"{type(error).__name__}: {error}",
    }


def rank_candidates(rows: pd.DataFrame) -> pd.DataFrame:
    """Rank p2 first, using d1/d2 only after p2 enters tolerance."""
    if rows.empty:
        return rows.copy()
    ranked = rows.copy()
    feasible = ranked["feasible"].fillna(False).astype(bool)
    aligned = feasible & ranked["p2_aligned"].fillna(False).astype(bool)
    p_error = pd.to_numeric(ranked["p2_error_thz"], errors="coerce").fillna(np.inf)
    d_error = pd.to_numeric(
        ranked["d1_d2_soft_error_thz"], errors="coerce"
    ).fillna(np.inf)
    ranked["_tier"] = np.select([aligned, feasible], [0, 1], default=2)
    ranked["_primary"] = np.where(aligned, d_error, p_error)
    ranked["_secondary"] = np.where(aligned, p_error, d_error)
    ranked = ranked.sort_values(
        ["_tier", "_primary", "_secondary", "evaluation_order"],
        kind="mergesort",
    ).drop(columns=["_tier", "_primary", "_secondary"])
    ranked = ranked.reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def row_score(row: dict[str, object]) -> tuple[float, float, float, float]:
    """Return the same lexicographic objective used in the persisted ranking."""
    if not bool(row.get("feasible", False)):
        return math.inf, math.inf, math.inf, float(row.get("evaluation_order", math.inf))
    p_error = float(row["p2_error_thz"])
    d_error = float(row["d1_d2_soft_error_thz"])
    if bool(row.get("p2_aligned", False)):
        return 0.0, d_error, p_error, float(row.get("evaluation_order", math.inf))
    return 1.0, p_error, d_error, float(row.get("evaluation_order", math.inf))


def _replace_coordinate(
    cell: CellParameters,
    coordinate: str,
    value: float,
) -> CellParameters:
    values = cell.to_dict()
    values[coordinate] = float(value)
    return CellParameters(**values)


def scan_coordinate(
    current: CellParameters,
    coordinate: str,
    coarse_step: float,
    fine_step: float,
    bounds: tuple[float, float],
    stage: str,
    evaluate: Callable[[CellParameters, str], dict[str, object]],
) -> CellParameters:
    """Probe both directions, walk while improving, then refine one fine step."""
    current_row = evaluate(current, stage)
    origin_value = float(getattr(current, coordinate))
    probes: list[tuple[int, CellParameters, dict[str, object]]] = []
    for direction in (-1, 1):
        value = origin_value + direction * coarse_step
        if bounds[0] - 1.0e-12 <= value <= bounds[1] + 1.0e-12:
            candidate = _replace_coordinate(current, coordinate, value)
            probes.append((direction, candidate, evaluate(candidate, stage)))

    best_direction = 0
    best_cell = current
    best_row = current_row
    for direction, candidate, row in probes:
        if row_score(row) < row_score(best_row):
            best_direction = direction
            best_cell = candidate
            best_row = row

    if best_direction:
        for _step_index in range(MAX_DIRECTION_STEPS - 1):
            value = float(getattr(best_cell, coordinate)) + best_direction * coarse_step
            if not (bounds[0] - 1.0e-12 <= value <= bounds[1] + 1.0e-12):
                break
            candidate = _replace_coordinate(best_cell, coordinate, value)
            row = evaluate(candidate, stage)
            if row_score(row) < row_score(best_row):
                best_cell, best_row = candidate, row
            else:
                break

    fine_candidates: list[
        tuple[int, CellParameters, dict[str, object]]
    ] = [(0, best_cell, best_row)]
    center = float(getattr(best_cell, coordinate))
    for direction in (-1, 1):
        value = center + direction * fine_step
        if bounds[0] - 1.0e-12 <= value <= bounds[1] + 1.0e-12:
            candidate = _replace_coordinate(best_cell, coordinate, value)
            fine_candidates.append(
                (direction, candidate, evaluate(candidate, f"{stage}_fine"))
            )
    fine_direction, fine_best_cell, fine_best_row = min(
        fine_candidates,
        key=lambda item: row_score(item[2]),
    )
    if fine_direction:
        for _step_index in range(MAX_DIRECTION_STEPS - 1):
            value = (
                float(getattr(fine_best_cell, coordinate))
                + fine_direction * fine_step
            )
            if not (bounds[0] - 1.0e-12 <= value <= bounds[1] + 1.0e-12):
                break
            candidate = _replace_coordinate(fine_best_cell, coordinate, value)
            row = evaluate(candidate, f"{stage}_fine")
            if row_score(row) < row_score(fine_best_row):
                fine_best_cell, fine_best_row = candidate, row
            else:
                break
    return fine_best_cell


def run_search(
    reference: dict[str, object],
    evaluator: Callable[[CellParameters, str, dict[str, object]], dict[str, object]] | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run the requested zeta -> eta -> b0 coordinate-priority search."""
    records: dict[tuple[float, float, float], dict[str, object]] = {}

    def persist() -> pd.DataFrame:
        ranked = rank_candidates(pd.DataFrame(records.values()))
        band.write_csv_atomic(OUT_DIR / "scan_summary.csv", ranked)
        return ranked

    def evaluate(cell: CellParameters, stage: str) -> dict[str, object]:
        key = candidate_key(cell.b0_nm, cell.eta, cell.zeta)
        if key not in records:
            if evaluator is None:
                row = evaluate_candidate(cell, stage, reference)
            else:
                row = evaluator(cell, stage, reference)
            row = dict(row)
            row["evaluation_order"] = len(records) + 1
            records[key] = row
            persist()
        return records[key]

    current = START_CLADDING
    start = evaluate(current, "start")
    if not bool(start.get("feasible", False)):
        raise RuntimeError(f"Starting cladding candidate failed: {start.get('error', start['status'])}")
    current = scan_coordinate(
        current,
        "zeta",
        ZETA_COARSE_STEP,
        ZETA_FINE_STEP,
        ZETA_BOUNDS,
        "zeta",
        evaluate,
    )
    current = scan_coordinate(
        current,
        "eta",
        ETA_COARSE_STEP,
        ETA_FINE_STEP,
        ETA_BOUNDS,
        "eta",
        evaluate,
    )
    current = scan_coordinate(
        current,
        "b0_nm",
        B0_COARSE_STEP_NM,
        B0_FINE_STEP_NM,
        B0_BOUNDS_NM,
        "b0",
        evaluate,
    )

    recheck = [(current, evaluate(current, "zeta_recheck"))]
    for offset in (-ZETA_FINE_STEP, ZETA_FINE_STEP):
        value = current.zeta + offset
        if ZETA_BOUNDS[0] <= value <= ZETA_BOUNDS[1]:
            candidate = _replace_coordinate(current, "zeta", value)
            recheck.append((candidate, evaluate(candidate, "zeta_recheck")))
    current = min(recheck, key=lambda item: row_score(item[1]))[0]

    ranked = persist()
    best = dict(evaluate(current, "final"))
    return ranked, best


def evaluate_candidate(
    cell: CellParameters,
    stage: str,
    reference: dict[str, object],
) -> dict[str, object]:
    """Run or load a cladding candidate and persist its restart metadata."""
    b0_nm, eta, zeta = candidate_key(cell.b0_nm, cell.eta, cell.zeta)
    cell = CellParameters(b0_nm=b0_nm, eta=eta, zeta=zeta)
    root = candidate_dir(b0_nm, eta, zeta)
    summary_path = root / "candidate_summary.json"
    if summary_path.is_file():
        try:
            cached = json.loads(summary_path.read_text(encoding="utf-8"))
            if (
                math.isclose(
                    float(cached["cavity_p2_frequency_thz"]),
                    float(reference["cavity_p2_frequency_thz"]),
                    abs_tol=1.0e-12,
                )
                and math.isclose(
                    float(cached["cavity_d1_frequency_thz"]),
                    float(reference["cavity_d1_frequency_thz"]),
                    abs_tol=1.0e-12,
                )
                and int(cached.get("eigenmode_count", -1)) == EIGENMODE_COUNT
            ):
                print(f"[{stage}] {cell.label('clad')} (cached)", flush=True)
                return cached
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            pass

    print(f"[{stage}] {cell.label('clad')} Gamma", flush=True)
    root.mkdir(parents=True, exist_ok=True)
    try:
        frame, selected, metadata = run_gamma_case(
            root,
            cell,
            "cladding_p2_d2_alignment",
        )
        summary = summarize_candidate(cell, stage, frame, selected, reference)
        error = ""
    except Exception as exc:  # keep completed candidates restartable
        metadata = {}
        summary = failed_candidate(cell, stage, reference, exc)
        error = summary["error"]
    summary["eigenmode_count"] = EIGENMODE_COUNT
    band.write_json_atomic(
        root / "config.json",
        {
            "cladding": cell.to_dict(),
            "stage": stage,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "eigenmode_count": EIGENMODE_COUNT,
            "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
            "reference": reference,
            "summary": summary,
            "error": error,
            **metadata,
        },
    )
    band.write_json_atomic(summary_path, summary)
    return summary


def scan_config() -> dict[str, object]:
    return {
        "workflow": "optimize_cladding_p2_d2_alignment",
        "fixed_cavity": FIXED_CAVITY.to_dict(),
        "start_cladding": START_CLADDING.to_dict(),
        "parameter_priority": ["zeta", "eta", "b0_nm"],
        "primary_target": "cladding p2@Gamma -> cavity p2@Gamma",
        "soft_target": "cladding d2@Gamma -> cavity d1@Gamma",
        "p2_alignment_tolerance_thz": P2_ALIGNMENT_TOLERANCE_THZ,
        "mesh_auto_size": MESH_AUTO_SIZE,
        "eigenmode_count": EIGENMODE_COUNT,
        "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
        "ranges": {
            "zeta": list(ZETA_BOUNDS),
            "eta": list(ETA_BOUNDS),
            "b0_nm": list(B0_BOUNDS_NM),
        },
        "coarse_steps": {
            "zeta": ZETA_COARSE_STEP,
            "eta": ETA_COARSE_STEP,
            "b0_nm": B0_COARSE_STEP_NM,
        },
        "fine_steps": {
            "zeta": ZETA_FINE_STEP,
            "eta": ETA_FINE_STEP,
            "b0_nm": B0_FINE_STEP_NM,
        },
        "max_direction_steps": MAX_DIRECTION_STEPS,
        "output_dir": OUT_DIR,
    }


def preflight() -> None:
    config = scan_config()
    print("Cladding p2/d2 quick-alignment preflight")
    print(f"  fixed cavity: {FIXED_CAVITY.label('cav')}")
    print(f"  start cladding: {START_CLADDING.label('clad')}")
    print("  targets: p2 -> p2 (primary); cavity d1 -> cladding d2 (soft)")
    print("  parameter priority: zeta -> eta -> b0")
    print(
        f"  solver: mesh{MESH_AUTO_SIZE}, modes={EIGENMODE_COUNT}, "
        f"shift={EIGENFREQUENCY_SHIFT}"
    )
    print(f"  ranges: {config['ranges']}")
    print(f"  output: {OUT_DIR}")
    print(
        "  first probes after start: "
        f"zeta={START_CLADDING.zeta - ZETA_COARSE_STEP:.3f}, "
        f"{START_CLADDING.zeta + ZETA_COARSE_STEP:.3f}"
    )


def main(argv: list[str] | None = None) -> dict[str, object] | None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print and persist the preflight configuration without starting COMSOL",
    )
    args = parser.parse_args(argv)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    preflight()
    band.write_json_atomic(OUT_DIR / "scan_config.json", scan_config())
    if args.dry_run:
        return None

    reference = load_cavity_reference()
    ranked, best = run_search(reference)
    result = {
        "fixed_cavity": FIXED_CAVITY.to_dict(),
        "cladding": {
            "b0_nm": float(best["b0_nm"]),
            "eta": float(best["eta"]),
            "zeta": float(best["zeta"]),
        },
        "summary": best,
        "candidate_count": int(len(ranked)),
        "reference": reference,
    }
    band.write_json_atomic(OUT_DIR / "best_parameters.json", result)
    print(
        "Best cladding: "
        f"clad({float(best['b0_nm']):.1f}-{float(best['eta']):.3f}-"
        f"{float(best['zeta']):.3f}); "
        f"p2 error={float(best['p2_error_thz']):.6g} THz; "
        f"d1/d2 soft error={float(best['d1_d2_soft_error_thz']):.6g} THz",
        flush=True,
    )
    print(f"Results: {OUT_DIR}", flush=True)
    return result


if __name__ == "__main__":
    main()
