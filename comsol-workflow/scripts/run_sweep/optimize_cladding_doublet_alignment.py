#!/usr/bin/env python3
"""Align C6 cladding p/d doublets to fixed cavity Gamma modes."""

from __future__ import annotations

import json
import math
import shutil
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


FIXED_CAVITY = CellParameters(b0_nm=245.0, eta=0.960, zeta=1.156)
START_CLADDING = CellParameters(b0_nm=240.2, eta=0.967, zeta=1.0)
CLADDING_ZETA = 1.0
MESH_AUTO_SIZE = 5
B0_MIN_NM = 235.0
B0_MAX_NM = 255.0
ETA_MIN = 0.900
ETA_MAX = 0.999
B0_RESOLUTION_NM = 0.1
ETA_RESOLUTION = 0.001
ALIGNMENT_TOLERANCE_THZ = 0.02
DEGENERACY_TOLERANCE_THZ = 0.02
MAX_NEWTON_ITERATIONS = 5
JACOBIAN_B0_STEP_NM = 1.0
JACOBIAN_ETA_STEP = 0.010
MAX_B0_UPDATE_NM = 3.0
MAX_ETA_UPDATE = 0.020
LOCAL_B0_OFFSETS_NM = (-0.2, -0.1, 0.0, 0.1, 0.2)
LOCAL_ETA_OFFSETS = (-0.002, -0.001, 0.0, 0.001, 0.002)
RUN_FULL_BAND_VALIDATION = False

GAMMA_POINT = {
    "direction": "gamma_m",
    "kx": 0.0,
    "ky": 0.0,
    "kx_str": "0.000000",
    "ky_str": "0.000000",
    "k_norm": 0.0,
    "path_fraction": 0.0,
}

OUT_DIR = UNIT_CELL_OUTPUT_ROOT / (
    "cladding_doublet_alignment_to_cavity_p1_d1_"
    f"{FIXED_CAVITY.label('cav')}_zeta1_mesh{MESH_AUTO_SIZE}"
)
LEGACY_D2_CACHE_DIR = UNIT_CELL_OUTPUT_ROOT / (
    "cladding_doublet_alignment_"
    f"{FIXED_CAVITY.label('cav')}_zeta1_mesh{MESH_AUTO_SIZE}"
)


def quantize_b0(value: float) -> float:
    bounded = min(max(float(value), B0_MIN_NM), B0_MAX_NM)
    return round(bounded / B0_RESOLUTION_NM) * B0_RESOLUTION_NM


def quantize_eta(value: float) -> float:
    bounded = min(max(float(value), ETA_MIN), ETA_MAX)
    return round(bounded / ETA_RESOLUTION) * ETA_RESOLUTION


def candidate_key(b0_nm: float, eta: float) -> tuple[float, float]:
    return round(quantize_b0(b0_nm), 1), round(quantize_eta(eta), 3)


def make_cladding_params(b0_nm: float, eta: float) -> dict[str, float | str]:
    """Return the exact-C6 compact-to-Fourier cladding parameters."""
    return {
        "name": "cladding_c6_doublet_alignment",
        "r_f0": float(eta),
        "b_square_f0": (float(b0_nm) / (band.B_0 * 1000.0)) ** 2,
        "b_square_f3": 0.0,
    }


def simulation_config() -> band.SimulationConfig:
    return band.SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift=band.EIGENFREQUENCY_SHIFT,
        eigenmode_count=band.EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )


def load_cavity_targets() -> dict[str, object]:
    """Run or reuse the fixed cavity Gamma point and return p1/d1 targets."""
    root = OUT_DIR / "cavity_reference"
    root.mkdir(parents=True, exist_ok=True)
    params = FIXED_CAVITY.to_fourier_params("cavity_p_bic")
    hole_params = band.get_hole_params(params)
    hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            f"Cavity clearance {minimum_clearance:g} is below d={band.D:g}"
        )
    geometry_path = root / "geometry.png"
    if not geometry_path.is_file():
        band.visualize_hexagon_design(
            hexagon,
            holes,
            info,
            filename=geometry_path,
            annotate=False,
        )
    config = simulation_config()
    frame = band.run_k_point(
        root,
        [hole.tolist() for hole in holes],
        GAMMA_POINT,
        config,
    )
    selected = band.select_gamma_bands(frame)
    targets = {
        "cavity_p1_frequency_thz": float(frame.loc[int(selected["p1"]), "re"]),
        "cavity_d1_frequency_thz": float(frame.loc[int(selected["d1"]), "re"]),
        "source_unit_cell_dir": root,
        "source_config": root / "config.json",
        "source_selected_bands": band.k_point_dir(root, GAMMA_POINT)
        / "eigenfrequencies.csv",
    }
    band.write_json_atomic(
        root / "config.json",
        {
            "fixed_cavity": FIXED_CAVITY.to_dict(),
            "params": params,
            "hole_params": hole_params,
            "triangles": holes,
            "minimum_clearance": minimum_clearance,
            "selected_gamma_modes": selected,
            "simulation": band.point_solver_identity(config),
            "targets": targets,
        },
    )
    return targets


def candidate_dir(b0_nm: float, eta: float) -> Path:
    b0_nm, eta = candidate_key(b0_nm, eta)
    return OUT_DIR / "candidates" / f"b0_{b0_nm:.1f}_eta_{eta:.3f}"


def legacy_candidate_dir(b0_nm: float, eta: float) -> Path:
    """Return the old d2-target cache location for the same geometry."""
    b0_nm, eta = candidate_key(b0_nm, eta)
    return LEGACY_D2_CACHE_DIR / "candidates" / f"b0_{b0_nm:.1f}_eta_{eta:.3f}"


def prepare_geometry(
    b0_nm: float,
    eta: float,
) -> tuple[list[list[list[float]]], dict[str, object]]:
    params = make_cladding_params(b0_nm, eta)
    hole_params = band.get_hole_params(params)
    hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    del hexagon
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            f"Cladding clearance {minimum_clearance:g} is below d={band.D:g}"
        )
    side_spread = float(np.ptp(hole_params[:, 2]))
    if side_spread > 1.0e-12:
        raise RuntimeError(
            f"zeta=1 did not produce equal triangle sides; spread={side_spread:g} um"
        )
    return [hole.tolist() for hole in holes], {
        "params": params,
        "hole_params": hole_params,
        "triangles": holes,
        "minimum_clearance": minimum_clearance,
        "triangle_side_spread_um": side_spread,
    }


def summarize_candidate(
    b0_nm: float,
    eta: float,
    gamma_modes: pd.DataFrame,
    targets: dict[str, object],
    stage: str,
) -> dict[str, object]:
    """Measure doublet means, numerical splittings, and both target errors."""
    b0_nm, eta = candidate_key(b0_nm, eta)
    result: dict[str, object] = {
        "stage": str(stage),
        "b0_nm": b0_nm,
        "eta": eta,
        "zeta": CLADDING_ZETA,
        "status": "invalid_modes",
        "feasible": False,
        "accepted": False,
        "cavity_p1_frequency_thz": float(targets["cavity_p1_frequency_thz"]),
        "cavity_d1_frequency_thz": float(targets["cavity_d1_frequency_thz"]),
        "cladding_p_doublet_frequency_thz": np.nan,
        "cladding_d_doublet_frequency_thz": np.nan,
        "p_signed_error_thz": np.nan,
        "d_signed_error_thz": np.nan,
        "p_error_thz": np.nan,
        "d_error_thz": np.nan,
        "maximum_alignment_error_thz": np.nan,
        "rms_alignment_error_thz": np.nan,
        "p_splitting_thz": np.nan,
        "d_splitting_thz": np.nan,
    }
    try:
        selected = band.select_gamma_bands(gamma_modes)
        p_indices = [int(selected[label]) for label in ("p1", "p2")]
        d_indices = [int(selected[label]) for label in ("d1", "d2")]
        chosen = gamma_modes.loc[p_indices + d_indices]
        if not chosen["is_valid"].astype(bool).all():
            result["status"] = "invalid_selected_mode"
            return result
        p_values = pd.to_numeric(
            gamma_modes.loc[p_indices, "re"], errors="coerce"
        ).to_numpy(dtype=float)
        d_values = pd.to_numeric(
            gamma_modes.loc[d_indices, "re"], errors="coerce"
        ).to_numpy(dtype=float)
        if not np.isfinite(np.concatenate([p_values, d_values])).all():
            result["status"] = "nonfinite_frequency"
            return result
    except (KeyError, TypeError, ValueError):
        return result

    p_frequency = float(np.mean(p_values))
    d_frequency = float(np.mean(d_values))
    p_split = float(np.ptp(p_values))
    d_split = float(np.ptp(d_values))
    p_signed = p_frequency - float(targets["cavity_p1_frequency_thz"])
    d_signed = d_frequency - float(targets["cavity_d1_frequency_thz"])
    p_error = abs(p_signed)
    d_error = abs(d_signed)
    maximum_error = max(p_error, d_error)
    rms_error = math.sqrt(0.5 * (p_signed**2 + d_signed**2))
    degeneracy_ok = (
        p_split <= DEGENERACY_TOLERANCE_THZ
        and d_split <= DEGENERACY_TOLERANCE_THZ
    )
    accepted = (
        degeneracy_ok
        and p_error <= ALIGNMENT_TOLERANCE_THZ
        and d_error <= ALIGNMENT_TOLERANCE_THZ
    )
    result.update(
        {
            "status": "accepted" if accepted else "ok",
            "feasible": True,
            "accepted": bool(accepted),
            "cladding_p_doublet_frequency_thz": p_frequency,
            "cladding_d_doublet_frequency_thz": d_frequency,
            "p_signed_error_thz": p_signed,
            "d_signed_error_thz": d_signed,
            "p_error_thz": p_error,
            "d_error_thz": d_error,
            "maximum_alignment_error_thz": maximum_error,
            "rms_alignment_error_thz": rms_error,
            "p_splitting_thz": p_split,
            "d_splitting_thz": d_split,
        }
    )
    return result


def rank_candidates(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return rows.copy()
    ranked = rows.copy()
    feasible = ranked["feasible"].fillna(False).astype(bool)
    accepted = ranked["accepted"].fillna(False).astype(bool)
    degeneracy_error = np.maximum(
        pd.to_numeric(ranked["p_splitting_thz"], errors="coerce"),
        pd.to_numeric(ranked["d_splitting_thz"], errors="coerce"),
    )
    ranked["_tier"] = np.select([accepted, feasible], [0, 1], default=2)
    ranked["_maximum"] = pd.to_numeric(
        ranked["maximum_alignment_error_thz"], errors="coerce"
    ).fillna(np.inf)
    ranked["_rms"] = pd.to_numeric(
        ranked["rms_alignment_error_thz"], errors="coerce"
    ).fillna(np.inf)
    ranked["_degeneracy"] = pd.Series(degeneracy_error, index=ranked.index).fillna(
        np.inf
    )
    ranked = ranked.sort_values(
        ["_tier", "_maximum", "_rms", "_degeneracy", "b0_nm", "eta"],
        kind="mergesort",
    ).drop(columns=["_tier", "_maximum", "_rms", "_degeneracy"])
    ranked = ranked.reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def evaluate_candidate(
    b0_nm: float,
    eta: float,
    targets: dict[str, object],
    stage: str,
) -> dict[str, object]:
    b0_nm, eta = candidate_key(b0_nm, eta)
    root = candidate_dir(b0_nm, eta)
    summary_path = root / "candidate_summary.json"
    if summary_path.is_file():
        cached = json.loads(summary_path.read_text(encoding="utf-8"))
        try:
            targets_match = math.isclose(
                float(cached["cavity_p1_frequency_thz"]),
                float(targets["cavity_p1_frequency_thz"]),
                abs_tol=1.0e-12,
            ) and math.isclose(
                float(cached["cavity_d1_frequency_thz"]),
                float(targets["cavity_d1_frequency_thz"]),
                abs_tol=1.0e-12,
            )
        except (KeyError, TypeError, ValueError):
            targets_match = False
        if targets_match and int(cached.get("eigenmode_count", -1)) == band.EIGENMODE_COUNT:
            cached["stage"] = str(stage)
            print(f"[{stage}] b0={b0_nm:.1f} nm eta={eta:.3f} (cached)")
            return cached

    root.mkdir(parents=True, exist_ok=True)
    holes, geometry = prepare_geometry(b0_nm, eta)
    point_dir = band.k_point_dir(root, GAMMA_POINT)
    legacy_root = legacy_candidate_dir(b0_nm, eta)
    legacy_point_dir = band.k_point_dir(legacy_root, GAMMA_POINT)
    config = simulation_config()
    if band.k_point_is_complete(point_dir, config):
        source_root = root
        state = "cached"
    elif band.k_point_is_complete(legacy_point_dir, config):
        source_root = legacy_root
        state = "legacy-cache"
    else:
        source_root = root
        state = "run"
    print(f"[{stage}] b0={b0_nm:.1f} nm eta={eta:.3f} Gamma ({state})")
    gamma_modes = band.run_k_point(source_root, holes, GAMMA_POINT, config)
    summary = summarize_candidate(b0_nm, eta, gamma_modes, targets, stage)
    summary["eigenmode_count"] = band.EIGENMODE_COUNT
    band.write_json_atomic(
        root / "config.json",
        {
            "b0_nm": b0_nm,
            "eta": eta,
            "zeta": CLADDING_ZETA,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "eigenmode_count": band.EIGENMODE_COUNT,
            "targets": targets,
            **geometry,
        },
    )
    band.write_json_atomic(summary_path, summary)
    return summary


def _finite_difference_probe(
    value: float,
    step: float,
    lower: float,
    upper: float,
) -> float:
    if value + step <= upper + 1.0e-12:
        return value + step
    if value - step >= lower - 1.0e-12:
        return value - step
    raise RuntimeError("Optimization range is smaller than the derivative step")


def adaptive_search(
    targets: dict[str, object],
    evaluator: Callable[[float, float, dict[str, object], str], dict[str, object]] = evaluate_candidate,
) -> tuple[dict[str, object], pd.DataFrame]:
    """Use a bounded local Jacobian, then an exact-resolution neighborhood."""
    records: dict[tuple[float, float], dict[str, object]] = {}

    def evaluate(b0_nm: float, eta: float, stage: str) -> dict[str, object]:
        key = candidate_key(b0_nm, eta)
        row = evaluator(*key, targets, stage)
        records[key] = row
        band.write_csv_atomic(
            OUT_DIR / "scan_summary.csv",
            rank_candidates(pd.DataFrame(records.values())),
        )
        return row

    current = candidate_key(
        START_CLADDING.b0_nm,
        min(START_CLADDING.eta, ETA_MAX),
    )
    for iteration in range(MAX_NEWTON_ITERATIONS):
        base = evaluate(*current, f"adaptive_{iteration}")
        if bool(base.get("accepted", False)):
            break
        if not bool(base.get("feasible", False)):
            raise RuntimeError(f"Adaptive base candidate is invalid: {base['status']}")

        b_probe = quantize_b0(
            _finite_difference_probe(
                current[0], JACOBIAN_B0_STEP_NM, B0_MIN_NM, B0_MAX_NM
            )
        )
        eta_probe = quantize_eta(
            _finite_difference_probe(
                current[1], JACOBIAN_ETA_STEP, ETA_MIN, ETA_MAX
            )
        )
        b_row = evaluate(b_probe, current[1], f"jacobian_b0_{iteration}")
        eta_row = evaluate(current[0], eta_probe, f"jacobian_eta_{iteration}")
        rows = (base, b_row, eta_row)
        if not all(bool(row.get("feasible", False)) for row in rows):
            raise RuntimeError("A Jacobian probe produced invalid selected modes")

        residual = np.array(
            [base["p_signed_error_thz"], base["d_signed_error_thz"]],
            dtype=float,
        )
        jacobian = np.column_stack(
            [
                (
                    np.array(
                        [b_row["p_signed_error_thz"], b_row["d_signed_error_thz"]],
                        dtype=float,
                    )
                    - residual
                )
                / (b_probe - current[0]),
                (
                    np.array(
                        [
                            eta_row["p_signed_error_thz"],
                            eta_row["d_signed_error_thz"],
                        ],
                        dtype=float,
                    )
                    - residual
                )
                / (eta_probe - current[1]),
            ]
        )
        update, *_ = np.linalg.lstsq(jacobian, -residual, rcond=None)
        update[0] = float(np.clip(update[0], -MAX_B0_UPDATE_NM, MAX_B0_UPDATE_NM))
        update[1] = float(np.clip(update[1], -MAX_ETA_UPDATE, MAX_ETA_UPDATE))
        proposed = candidate_key(current[0] + update[0], current[1] + update[1])
        if proposed == current:
            break
        current = proposed

    preliminary = rank_candidates(pd.DataFrame(records.values()))
    if preliminary.empty or not bool(preliminary.iloc[0]["feasible"]):
        raise RuntimeError("Adaptive search found no feasible cladding candidate")
    center_b0 = float(preliminary.iloc[0]["b0_nm"])
    center_eta = float(preliminary.iloc[0]["eta"])
    for b_offset in LOCAL_B0_OFFSETS_NM:
        for eta_offset in LOCAL_ETA_OFFSETS:
            b0_nm, eta = candidate_key(center_b0 + b_offset, center_eta + eta_offset)
            if (b0_nm, eta) not in records:
                evaluate(b0_nm, eta, "local_refinement")

    ranked = rank_candidates(pd.DataFrame(records.values()))
    best = ranked.iloc[0].drop(labels="rank").to_dict()
    return best, ranked


def run_full_band_validation(best: dict[str, object]) -> dict[str, object]:
    output_dir = OUT_DIR / (
        "full_band_validation_"
        f"clad({float(best['b0_nm']):.1f}-{float(best['eta']):.3f}-1)"
    )
    return band.run_band_workflow(
        output_dir,
        FIXED_CAVITY.to_fourier_params("cavity_p_bic"),
        make_cladding_params(float(best["b0_nm"]), float(best["eta"])),
    )


def preflight(targets: dict[str, object]) -> None:
    print("Cladding C6 doublet alignment preflight")
    print(f"  cavity: {FIXED_CAVITY.to_dict()} (fixed)")
    print(f"  cladding start: {START_CLADDING.to_dict()}")
    print(
        "  target: cladding p1/p2@Gamma mean -> "
        f"cavity p1@Gamma={float(targets['cavity_p1_frequency_thz']):.12g} THz"
    )
    print(
        "  target: cladding d1/d2@Gamma mean -> "
        f"cavity d1@Gamma={float(targets['cavity_d1_frequency_thz']):.12g} THz"
    )
    print(
        f"  cladding: b0=[{B0_MIN_NM:g}, {B0_MAX_NM:g}] nm, "
        f"eta=[{ETA_MIN:g}, {ETA_MAX:g}], zeta=1 exactly"
    )
    print(
        f"  precision: b0={B0_RESOLUTION_NM:g} nm, eta={ETA_RESOLUTION:g}; "
        f"tolerance={ALIGNMENT_TOLERANCE_THZ:g} THz"
    )
    print(
        f"  mesh: auto size {MESH_AUTO_SIZE}; "
        f"solver shift={band.EIGENFREQUENCY_SHIFT}; modes={band.EIGENMODE_COUNT}"
    )
    print(f"  cache/output: {OUT_DIR}")
    validation = "enabled" if RUN_FULL_BAND_VALIDATION else "disabled"
    print(f"  full-band validation: {validation}")


def main() -> None:
    targets = load_cavity_targets()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    preflight(targets)
    band.write_json_atomic(
        OUT_DIR / "scan_config.json",
        {
            "workflow": "optimize_cladding_doublet_alignment",
            "shared_parameters": ACTIVE_PARAMETERS.to_metadata(),
            "fixed_cavity": FIXED_CAVITY.to_dict(),
            "start_cladding": START_CLADDING.to_dict(),
            "fixed_cladding_zeta": CLADDING_ZETA,
            "b0_range_nm": [B0_MIN_NM, B0_MAX_NM],
            "eta_range": [ETA_MIN, ETA_MAX],
            "precision": {
                "b0_nm": B0_RESOLUTION_NM,
                "eta": ETA_RESOLUTION,
            },
            "alignment_tolerance_thz": ALIGNMENT_TOLERANCE_THZ,
            "degeneracy_tolerance_thz": DEGENERACY_TOLERANCE_THZ,
            "targets": targets,
        },
    )
    best, ranked = adaptive_search(targets)
    band.write_csv_atomic(OUT_DIR / "scan_summary.csv", ranked)
    band.write_json_atomic(
        OUT_DIR / "best_parameters.json",
        {
            "cavity": FIXED_CAVITY.to_dict(),
            "cladding": {
                "b0_nm": float(best["b0_nm"]),
                "eta": float(best["eta"]),
                "zeta": CLADDING_ZETA,
            },
            "summary": best,
            "targets": targets,
        },
    )
    if not bool(best.get("accepted", False)):
        raise RuntimeError(
            "No candidate met both alignment and degeneracy tolerances; "
            f"closest result is recorded in {OUT_DIR / 'best_parameters.json'}"
        )

    if RUN_FULL_BAND_VALIDATION:
        run_full_band_validation(best)
    candidates = OUT_DIR / "candidates"
    if candidates.is_dir():
        shutil.rmtree(candidates)
    print(f"Best cladding parameters: b0={best['b0_nm']:.1f} nm, eta={best['eta']:.3f}, zeta=1")
    print(f"Results: {OUT_DIR}")


if __name__ == "__main__":
    main()
