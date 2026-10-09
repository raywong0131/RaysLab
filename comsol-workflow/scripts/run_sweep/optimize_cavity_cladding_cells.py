#!/usr/bin/env python3
"""Optimize cavity user-py Q at Gamma, then align cladding user-px at Gamma."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import UNIT_CELL_OUTPUT_ROOT  # noqa: E402


MESH_AUTO_SIZE = 9
FOURIER_REFERENCE_B0_NM = 230.0

CAVITY_B0_NM = 245.0
CAVITY_ETA_MIN = 0.950
CAVITY_ETA_MAX = 0.970
CAVITY_ZETA_MIN = 1.146
CAVITY_ZETA_MAX = 1.166
CAVITY_COARSE_STEP = 0.005
CAVITY_FINE_STEP = 0.001
CAVITY_FINE_HALF_WIDTH = 0.001
CAVITY_REFERENCE_FREQUENCY_THZ = 198.4493726243485
CAVITY_FREQUENCY_TOLERANCE_THZ = 0.1

CLADDING_B0_MIN_NM = 239.0
CLADDING_B0_MAX_NM = 245.0
CLADDING_B0_COARSE_STEP_NM = 1.0
CLADDING_B0_FINE_STEP_NM = 0.1
CLADDING_ZETA_MIN = 0.910
CLADDING_ZETA_MAX = 0.950
CLADDING_ZETA_COARSE_STEP = 0.005
CLADDING_ZETA_FINE_STEP = 0.001
CLADDING_ETA_MIN = 0.970
CLADDING_ETA_MAX = 0.990
CLADDING_ETA_COARSE_STEP = 0.005
CLADDING_ETA_FINE_STEP = 0.001
CLADDING_INITIAL_B0_NM = 242.0
CLADDING_INITIAL_ETA = 0.980
CLADDING_INITIAL_ZETA = 0.930
CLADDING_ALIGNMENT_TOLERANCE_THZ = 0.05

MIN_P_WEIGHT = 0.90
MIN_USER_COMPONENT_SHARE = 0.80
MIN_GAMMA_GAP_THZ = 0.0
K_STRING_PRECISION = 6
SCREEN_K_MAGNITUDES = (0.0, 0.01, 0.03, 0.06, 0.10)
FINAL_K_MAGNITUDES = tuple(
    sorted(
        {
            *(round(value, 6) for value in np.arange(0.0, 0.200001, 0.005)),
            0.002,
            0.004,
            0.006,
        }
    )
)
FINAL_CAVITY_CANDIDATE_COUNT = 3

# The user accepted this densely validated candidate as the closest practical
# cavity result.  Reuse its cached data and continue directly with cladding.
USE_ACCEPTED_CAVITY_CANDIDATE = True
ACCEPTED_CAVITY_ETA = 0.960
ACCEPTED_CAVITY_ZETA = 1.155

# Keep the root short enough for deeply nested Windows k-point/parquet paths.
# The complete parameter ranges are recorded in scan_config.json.
OUT_DIR = UNIT_CELL_OUTPUT_ROOT / "cell_pair_opt_mesh9"


def decimal_grid(start: float, stop: float, step: float, digits: int) -> tuple[float, ...]:
    """Return an inclusive decimal grid without floating-point endpoint drift."""
    scale = 10**digits
    lo = round(float(start) * scale)
    hi = round(float(stop) * scale)
    stride = round(float(step) * scale)
    if stride <= 0:
        raise ValueError("grid step must be positive")
    return tuple(value / scale for value in range(lo, hi + 1, stride))


def make_cell_params(
    name: str,
    b0_nm: float,
    eta: float,
    zeta: float,
) -> dict[str, float | str]:
    """Translate user compact geometry parameters into the existing Fourier input."""
    return {
        "name": str(name),
        "r_f0": float(eta),
        "b_square_f0": (float(b0_nm) / FOURIER_REFERENCE_B0_NM) ** 2,
        "b_square_f3": (float(zeta) ** 2 - 1.0) / 2.0,
    }


def make_cavity_params(eta: float, zeta: float) -> dict[str, float | str]:
    return make_cell_params("cavity_py_gamma_optimization", CAVITY_B0_NM, eta, zeta)


def make_cladding_params(
    b0_nm: float,
    eta: float,
    zeta: float,
) -> dict[str, float | str]:
    return make_cell_params("cladding_px_alignment", b0_nm, eta, zeta)


def internal_component_for_user_mode(user_mode: str) -> str:
    """Map user-visible p-mode names onto stored project-internal weight fields."""
    mapping = {"px": "py_weight", "py": "px_weight"}
    try:
        return mapping[str(user_mode)]
    except KeyError as error:
        raise ValueError(f"Unknown user p mode: {user_mode}") from error


def identify_user_p_mode(
    gamma_modes: pd.DataFrame,
    selected_modes: dict[str, int],
    user_mode: str,
) -> tuple[str, int, pd.Series, float]:
    """Select the p band with the largest share of the requested user mode."""
    component = internal_component_for_user_mode(user_mode)
    candidates = []
    for label in ("p1", "p2"):
        mode_idx = int(selected_modes[label])
        row = gamma_modes.loc[mode_idx]
        px = float(row["px_weight"])
        py = float(row["py_weight"])
        total = px + py
        share = float(row[component]) / total if math.isfinite(total) and total > 0 else np.nan
        candidates.append((share, label, mode_idx, row))
    finite = [item for item in candidates if math.isfinite(item[0])]
    if not finite:
        raise ValueError(f"Selected p modes have no finite {user_mode} composition")
    share, label, mode_idx, row = max(finite, key=lambda item: item[0])
    return label, mode_idx, row, float(share)


def make_direction_points(
    direction: str,
    magnitudes: Sequence[float],
) -> list[dict[str, float | str]]:
    """Build one Gamma branch using the unit-cell workflow point schema."""
    endpoints = {"gamma_m": 0.5, "gamma_k": 1.0 / math.sqrt(3.0)}
    if direction not in endpoints:
        raise ValueError(f"Unknown direction: {direction}")
    endpoint = endpoints[direction]
    points = []
    for raw in sorted({round(float(value), K_STRING_PRECISION) for value in magnitudes}):
        if not 0.0 <= raw <= endpoint:
            raise ValueError(f"k magnitude {raw:g} is outside {direction}")
        kx = raw if direction == "gamma_k" else 0.0
        ky = raw if direction == "gamma_m" else 0.0
        points.append(
            {
                "direction": direction,
                "kx": kx,
                "ky": ky,
                "kx_str": f"{kx:.{K_STRING_PRECISION}f}",
                "ky_str": f"{ky:.{K_STRING_PRECISION}f}",
                "k_norm": raw,
                "path_fraction": raw / endpoint,
                "path_coordinate": (-1.0 if direction == "gamma_m" else 1.0)
                * raw,
            }
        )
    return points


def branch_points(
    magnitudes: Sequence[float],
) -> dict[str, list[dict[str, float | str]]]:
    return {
        direction: make_direction_points(direction, magnitudes)
        for direction in ("gamma_m", "gamma_k")
    }


def unique_points(
    branches: dict[str, list[dict[str, float | str]]],
) -> list[dict[str, float | str]]:
    unique = {}
    for direction in ("gamma_m", "gamma_k"):
        for point in branches[direction]:
            unique.setdefault((point["kx_str"], point["ky_str"]), point)
    return list(unique.values())


def simulation_config() -> band.SimulationConfig:
    """Use the main unit-cell physics with an explicit mesh9 setting."""
    return band.SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift=band.EIGENFREQUENCY_SHIFT,
        eigenmode_count=band.EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )


def prepare_geometry(
    case_dir: Path,
    params: dict[str, float | str],
) -> tuple[list[list[list[float]]], dict[str, object]]:
    """Resolve, validate, and persist one compact unit-cell geometry."""
    case_dir.mkdir(parents=True, exist_ok=True)
    hole_params = band.get_hole_params(params)
    hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            f"Unit-cell clearance {minimum_clearance:g} is below d={band.D:g}"
        )
    geometry_path = case_dir / "geometry.png"
    if not geometry_path.exists():
        band.visualize_hexagon_design(
            hexagon,
            holes,
            info,
            filename=geometry_path,
            annotate=False,
        )
    return [hole.tolist() for hole in holes], {
        "hole_params": hole_params,
        "triangles": holes,
        "minimum_clearance": minimum_clearance,
    }


def run_target_band_case(
    case_dir: Path,
    params: dict[str, float | str],
    user_mode: str,
    magnitudes: Sequence[float],
    stage: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run/resume k points and track the p band matching one user-visible mode."""
    case_dir = Path(case_dir)
    holes, geometry = prepare_geometry(case_dir, params)
    branches = branch_points(magnitudes)
    points = unique_points(branches)
    config = simulation_config()
    for point_index, point in enumerate(points, start=1):
        point_dir = band.k_point_dir(case_dir, point)
        state = "cached" if band.k_point_is_complete(point_dir, config) else "run"
        print(
            f"[{stage}] {point_index}/{len(points)} "
            f"kx={point['kx_str']} ky={point['ky_str']} ({state})"
        )
        band.run_k_point(case_dir, holes, point, config)

    gamma_dir = band.k_point_dir(case_dir, branches["gamma_m"][0])
    gamma_modes = pd.read_csv(gamma_dir / "eigenfrequencies.csv")
    selected_modes = band.select_gamma_bands(gamma_modes)
    target_label, target_mode_idx, target_gamma, user_share = identify_user_p_mode(
        gamma_modes, selected_modes, user_mode
    )
    tracked = []
    for direction in ("gamma_m", "gamma_k"):
        frame = band.track_selected_branch(case_dir, branches[direction], selected_modes)
        tracked.append(frame[frame["band_label"] == target_label].copy())
    target = pd.concat(tracked, ignore_index=True)
    target.insert(0, "stage", str(stage))
    target.insert(1, "user_mode", str(user_mode))
    target.insert(2, "user_mode_share_gamma", user_share)
    band.write_csv_atomic(case_dir / f"selected_user_{user_mode}.csv", target)
    metadata = {
        "stage": stage,
        "user_mode": user_mode,
        "target_band_label": target_label,
        "target_mode_idx": int(target_mode_idx),
        "target_gamma_frequency_thz": float(target_gamma["re"]),
        "target_gamma_q": float(target_gamma["q"]),
        "target_user_mode_share": user_share,
        "selected_gamma_modes": selected_modes,
        "branches": branches,
        **geometry,
    }
    return target, metadata


def summarize_cavity_candidate(
    eta: float,
    zeta: float,
    target: pd.DataFrame,
    metadata: dict[str, object],
    stage: str,
) -> dict[str, object]:
    """Apply cavity frequency/composition constraints and Q-at-Gamma objective."""
    result = {
        "stage": stage,
        "eta": round(float(eta), 3),
        "zeta": round(float(zeta), 3),
        "status": "invalid",
        "feasible": False,
        "peak_at_gamma": False,
        "gamma_frequency_thz": float(metadata["target_gamma_frequency_thz"]),
        "frequency_drift_thz": float(metadata["target_gamma_frequency_thz"])
        - CAVITY_REFERENCE_FREQUENCY_THZ,
        "gamma_q": float(metadata["target_gamma_q"]),
        "gamma_p_weight": np.nan,
        "user_py_share": float(metadata["target_user_mode_share"]),
        "peak_excess_dex": np.nan,
        "peak_k_signed": np.nan,
        "peak_direction": "",
        "target_band_label": metadata["target_band_label"],
    }
    if (target["match_status"] == "unmatched").any():
        result["status"] = "unmatched_target"
        return result
    gamma = target[np.isclose(pd.to_numeric(target["k_norm"], errors="coerce"), 0.0)]
    if gamma.empty:
        result["status"] = "missing_gamma"
        return result
    gamma_row = gamma.iloc[0]
    result["gamma_p_weight"] = float(gamma_row["p_weight"])
    if not bool(gamma_row["is_valid"]):
        result["status"] = "invalid_gamma"
        return result
    if abs(float(result["frequency_drift_thz"])) > CAVITY_FREQUENCY_TOLERANCE_THZ:
        result["status"] = "frequency_drift_outside_tolerance"
        return result
    if float(result["gamma_p_weight"]) < MIN_P_WEIGHT:
        result["status"] = "p_weight_below_minimum"
        return result
    if float(result["user_py_share"]) < MIN_USER_COMPONENT_SHARE:
        result["status"] = "user_py_share_below_minimum"
        return result
    q = pd.to_numeric(target["q"], errors="coerce")
    if not np.all(np.isfinite(q) & (q > 0.0)):
        result["status"] = "invalid_q"
        return result
    nonzero = target[~np.isclose(pd.to_numeric(target["k_norm"]), 0.0)].copy()
    if nonzero.empty:
        result["status"] = "missing_nonzero_k"
        return result
    nonzero["peak_excess_dex"] = np.log10(nonzero["q"].astype(float)) - math.log10(
        float(result["gamma_q"])
    )
    peak = nonzero.loc[nonzero["peak_excess_dex"].idxmax()]
    peak_excess = float(peak["peak_excess_dex"])
    peak_at_gamma = peak_excess <= 1.0e-9
    result.update(
        {
            "status": "ok",
            "feasible": True,
            "peak_at_gamma": bool(peak_at_gamma),
            "peak_excess_dex": peak_excess,
            "peak_k_signed": (
                0.0
                if peak_at_gamma
                else (-1.0 if peak["direction"] == "gamma_m" else 1.0)
                * float(peak["k_norm"])
            ),
            "peak_direction": "gamma" if peak_at_gamma else str(peak["direction"]),
        }
    )
    return result


def rank_cavity_candidates(rows: pd.DataFrame) -> pd.DataFrame:
    """Prefer valid Gamma maxima, then larger Gamma Q and smaller frequency drift."""
    ranked = rows.copy()
    if ranked.empty:
        return ranked.assign(rank=pd.Series(dtype=int))
    feasible = ranked["feasible"].fillna(False).astype(bool)
    at_gamma = ranked["peak_at_gamma"].fillna(False).astype(bool)
    ranked["_tier"] = np.where(feasible & at_gamma, 0, np.where(feasible, 1, 2))
    ranked["_primary"] = np.where(
        ranked["_tier"] == 0,
        -pd.to_numeric(ranked["gamma_q"], errors="coerce").fillna(-np.inf),
        pd.to_numeric(ranked["peak_excess_dex"], errors="coerce").fillna(np.inf),
    )
    ranked["_drift"] = pd.to_numeric(
        ranked["frequency_drift_thz"], errors="coerce"
    ).abs().fillna(np.inf)
    ranked = ranked.sort_values(
        ["_tier", "_primary", "_drift", "eta", "zeta"], kind="mergesort"
    ).drop(columns=["_tier", "_primary", "_drift"])
    ranked = ranked.reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def cavity_candidate_dir(eta: float, zeta: float) -> Path:
    return (
        OUT_DIR
        / "cavity"
        / "candidates"
        / f"eta_{float(eta):.3f}_zeta_{float(zeta):.3f}"
    )


def evaluate_cavity(
    eta: float,
    zeta: float,
    stage: str,
    magnitudes: Sequence[float],
) -> tuple[pd.DataFrame, dict[str, object]]:
    eta = round(float(eta), 3)
    zeta = round(float(zeta), 3)
    root = cavity_candidate_dir(eta, zeta)
    params = make_cavity_params(eta, zeta)
    target, metadata = run_target_band_case(root, params, "py", magnitudes, stage)
    summary = summarize_cavity_candidate(eta, zeta, target, metadata, stage)
    band.write_json_atomic(
        root / "config.json",
        {
            "eta": eta,
            "zeta": zeta,
            "params": params,
            "eigenmode_count": band.EIGENMODE_COUNT,
            "summary": summary,
            **metadata,
        },
    )
    return target, summary


def persist_cavity_summary(records: dict[tuple[float, float], dict[str, object]]) -> pd.DataFrame:
    ranked = rank_cavity_candidates(pd.DataFrame(records.values()))
    band.write_csv_atomic(OUT_DIR / "cavity" / "scan_summary.csv", ranked)
    return ranked


def optimize_cavity() -> tuple[dict[str, object], pd.DataFrame]:
    """Run coarse 2D, local 0.001 refinement, and dense validation."""
    records: dict[tuple[float, float], dict[str, object]] = {}
    target_rows: dict[tuple[float, float], pd.DataFrame] = {}

    def evaluate(eta: float, zeta: float, stage: str, magnitudes=SCREEN_K_MAGNITUDES):
        key = (round(float(eta), 3), round(float(zeta), 3))
        target_rows[key], records[key] = evaluate_cavity(*key, stage, magnitudes)
        persist_cavity_summary(records)
        return records[key]

    for eta in decimal_grid(CAVITY_ETA_MIN, CAVITY_ETA_MAX, CAVITY_COARSE_STEP, 3):
        for zeta in decimal_grid(
            CAVITY_ZETA_MIN, CAVITY_ZETA_MAX, CAVITY_COARSE_STEP, 3
        ):
            evaluate(eta, zeta, "coarse")

    coarse = persist_cavity_summary(records)
    feasible = coarse[coarse["feasible"].astype(bool)]
    if feasible.empty:
        raise RuntimeError("No coarse cavity candidate satisfies the frequency constraints")
    center = feasible.iloc[0]
    fine_etas = decimal_grid(
        max(CAVITY_ETA_MIN, float(center["eta"]) - CAVITY_FINE_HALF_WIDTH),
        min(CAVITY_ETA_MAX, float(center["eta"]) + CAVITY_FINE_HALF_WIDTH),
        CAVITY_FINE_STEP,
        3,
    )
    fine_zetas = decimal_grid(
        max(CAVITY_ZETA_MIN, float(center["zeta"]) - CAVITY_FINE_HALF_WIDTH),
        min(CAVITY_ZETA_MAX, float(center["zeta"]) + CAVITY_FINE_HALF_WIDTH),
        CAVITY_FINE_STEP,
        3,
    )
    for eta in fine_etas:
        for zeta in fine_zetas:
            key = (round(eta, 3), round(zeta, 3))
            if key not in records:
                evaluate(*key, "fine")

    fine_ranked = persist_cavity_summary(records)
    dense_candidates = fine_ranked[fine_ranked["feasible"].astype(bool)].head(
        FINAL_CAVITY_CANDIDATE_COUNT
    )
    dense_rows = []
    for row in dense_candidates.itertuples(index=False):
        dense_rows.append(evaluate(row.eta, row.zeta, "final", FINAL_K_MAGNITUDES))
    final_ranked = rank_cavity_candidates(pd.DataFrame(dense_rows))
    if final_ranked.empty or not bool(final_ranked.iloc[0]["peak_at_gamma"]):
        raise RuntimeError("No densely validated cavity candidate has its Q maximum at Gamma")
    best = final_ranked.iloc[0].drop(labels="rank").to_dict()
    best_key = (round(float(best["eta"]), 3), round(float(best["zeta"]), 3))
    best_dir = OUT_DIR / "cavity" / "best"
    best_dir.mkdir(parents=True, exist_ok=True)
    band.write_csv_atomic(best_dir / "selected_user_py.csv", target_rows[best_key])
    band.write_json_atomic(
        best_dir / "best_config.json",
        {
            "summary": best,
            "b0_nm": CAVITY_B0_NM,
            "params": make_cavity_params(*best_key),
            "reference_frequency_thz": CAVITY_REFERENCE_FREQUENCY_THZ,
            "frequency_tolerance_thz": CAVITY_FREQUENCY_TOLERANCE_THZ,
            "mesh_auto_size": MESH_AUTO_SIZE,
        },
    )
    band.write_csv_atomic(best_dir / "dense_validation_summary.csv", final_ranked)
    return best, target_rows[best_key]


def load_accepted_cavity_candidate() -> tuple[dict[str, object], pd.DataFrame]:
    """Promote the user-accepted cached cavity candidate without rerunning it."""
    eta = round(float(ACCEPTED_CAVITY_ETA), 3)
    zeta = round(float(ACCEPTED_CAVITY_ZETA), 3)
    candidate_dir = cavity_candidate_dir(eta, zeta)
    config_path = candidate_dir / "config.json"
    target_path = candidate_dir / "selected_user_py.csv"
    if not config_path.is_file() or not target_path.is_file():
        raise FileNotFoundError(
            "Accepted cavity candidate is incomplete: "
            f"expected {config_path} and {target_path}"
        )

    payload = json.loads(config_path.read_text(encoding="utf-8"))
    if int(payload.get("eigenmode_count", -1)) != band.EIGENMODE_COUNT:
        raise RuntimeError("Accepted cavity mode count does not match parameter.json")
    summary = dict(payload["summary"])
    if not math.isclose(float(summary["eta"]), eta, abs_tol=5e-7):
        raise RuntimeError("Accepted cavity eta does not match its cached summary")
    if not math.isclose(float(summary["zeta"]), zeta, abs_tol=5e-7):
        raise RuntimeError("Accepted cavity zeta does not match its cached summary")
    if not bool(summary.get("feasible", False)):
        raise RuntimeError("Accepted cavity candidate is not frequency-feasible")
    if abs(float(summary["frequency_drift_thz"])) > CAVITY_FREQUENCY_TOLERANCE_THZ:
        raise RuntimeError("Accepted cavity candidate exceeds the frequency tolerance")

    summary["selection_status"] = "user_accepted_closest_to_gamma"
    target_rows = pd.read_csv(target_path)
    best_dir = OUT_DIR / "cavity" / "best"
    best_dir.mkdir(parents=True, exist_ok=True)
    band.write_csv_atomic(best_dir / "selected_user_py.csv", target_rows)
    band.write_json_atomic(
        best_dir / "best_config.json",
        {
            "summary": summary,
            "b0_nm": CAVITY_B0_NM,
            "params": make_cavity_params(eta, zeta),
            "reference_frequency_thz": CAVITY_REFERENCE_FREQUENCY_THZ,
            "frequency_tolerance_thz": CAVITY_FREQUENCY_TOLERANCE_THZ,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "source_candidate_dir": candidate_dir,
        },
    )
    return summary, target_rows


GAMMA_POINT = make_direction_points("gamma_m", (0.0,))[0]


def cladding_candidate_dir(b0_nm: float, eta: float, zeta: float) -> Path:
    return (
        OUT_DIR
        / "cladding"
        / "candidates"
        / f"b0_{float(b0_nm):.1f}_eta_{float(eta):.3f}_zeta_{float(zeta):.3f}"
    )


def evaluate_cladding(
    b0_nm: float,
    eta: float,
    zeta: float,
    target_frequency_thz: float,
    stage: str,
) -> dict[str, object]:
    """Evaluate user-facing cladding px at Gamma against the cavity py target."""
    b0_nm = round(float(b0_nm), 1)
    eta = round(float(eta), 3)
    zeta = round(float(zeta), 3)
    root = cladding_candidate_dir(b0_nm, eta, zeta)
    params = make_cladding_params(b0_nm, eta, zeta)
    holes, geometry = prepare_geometry(root, params)
    point_dir = band.k_point_dir(root, GAMMA_POINT)
    config = simulation_config()
    state = "cached" if band.k_point_is_complete(point_dir, config) else "run"
    print(
        f"[{stage}] cladding b0={b0_nm:.1f} eta={eta:.3f} zeta={zeta:.3f} "
        f"Gamma ({state})"
    )
    frame = band.run_k_point(root, holes, GAMMA_POINT, config)
    selected = band.select_gamma_bands(frame)
    label, mode_idx, mode, share = identify_user_p_mode(frame, selected, "px")
    p_indices = [int(selected[name]) for name in ("p1", "p2")]
    d_indices = [int(selected[name]) for name in ("d1", "d2")]
    lower_edge = max(float(frame.loc[index, "re"]) for index in p_indices)
    upper_edge = min(float(frame.loc[index, "re"]) for index in d_indices)
    gap = upper_edge - lower_edge
    frequency = float(mode["re"])
    alignment_error = abs(frequency - float(target_frequency_thz))
    midgap = 0.5 * (lower_edge + upper_edge)
    feasible = (
        bool(mode["is_valid"])
        and float(mode["p_weight"]) >= MIN_P_WEIGHT
        and share >= MIN_USER_COMPONENT_SHARE
        and gap > MIN_GAMMA_GAP_THZ
    )
    summary = {
        "stage": stage,
        "b0_nm": b0_nm,
        "eta": eta,
        "zeta": zeta,
        "status": "ok" if feasible else "invalid",
        "feasible": bool(feasible),
        "within_alignment_tolerance": bool(
            feasible and alignment_error <= CLADDING_ALIGNMENT_TOLERANCE_THZ
        ),
        "cavity_py_frequency_thz": float(target_frequency_thz),
        "cladding_px_frequency_thz": frequency,
        "signed_alignment_error_thz": frequency - float(target_frequency_thz),
        "alignment_error_thz": alignment_error,
        "cladding_px_share": share,
        "cladding_px_p_weight": float(mode["p_weight"]),
        "target_band_label": label,
        "target_mode_idx": int(mode_idx),
        "gap_lower_edge_thz": lower_edge,
        "gap_upper_edge_thz": upper_edge,
        "gamma_gap_thz": gap,
        "midgap_thz": midgap,
        "midgap_error_thz": abs(midgap - float(target_frequency_thz)),
        "selected_gamma_modes": selected,
    }
    band.write_json_atomic(
        root / "config.json",
        {
            "params": params,
            "summary": summary,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "eigenmode_count": band.EIGENMODE_COUNT,
            **geometry,
        },
    )
    return summary


def rank_cladding_candidates(rows: pd.DataFrame) -> pd.DataFrame:
    """Rank valid candidates by alignment, then midgap placement and perturbation."""
    ranked = rows.copy()
    if ranked.empty:
        return ranked.assign(rank=pd.Series(dtype=int))
    feasible = ranked["feasible"].fillna(False).astype(bool)
    aligned = feasible & ranked["within_alignment_tolerance"].fillna(False).astype(bool)
    ranked["_tier"] = np.select([aligned, feasible], [0, 1], default=2)
    ranked["_alignment"] = pd.to_numeric(
        ranked["alignment_error_thz"], errors="coerce"
    ).fillna(np.inf)
    ranked["_midgap"] = pd.to_numeric(
        ranked["midgap_error_thz"], errors="coerce"
    ).fillna(np.inf)
    ranked["_distance"] = (
        (pd.to_numeric(ranked["b0_nm"]) - CLADDING_INITIAL_B0_NM).abs()
        / (CLADDING_B0_MAX_NM - CLADDING_B0_MIN_NM)
        + (pd.to_numeric(ranked["zeta"]) - CLADDING_INITIAL_ZETA).abs()
        / (CLADDING_ZETA_MAX - CLADDING_ZETA_MIN)
        + (pd.to_numeric(ranked["eta"]) - CLADDING_INITIAL_ETA).abs()
        / (CLADDING_ETA_MAX - CLADDING_ETA_MIN)
    )
    ranked = ranked.sort_values(
        ["_tier", "_alignment", "_midgap", "_distance", "b0_nm", "zeta", "eta"],
        kind="mergesort",
    ).drop(columns=["_tier", "_alignment", "_midgap", "_distance"])
    ranked = ranked.reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def optimize_cladding(target_frequency_thz: float) -> dict[str, object]:
    """Tune b0 first, zeta second, and eta only when alignment still fails."""
    records: dict[tuple[float, float, float], dict[str, object]] = {}

    def evaluate(b0_nm: float, eta: float, zeta: float, stage: str):
        key = (round(float(b0_nm), 1), round(float(eta), 3), round(float(zeta), 3))
        if key not in records:
            records[key] = evaluate_cladding(*key, target_frequency_thz, stage)
            persist()
        return records[key]

    def persist() -> pd.DataFrame:
        ranked = rank_cladding_candidates(pd.DataFrame(records.values()))
        band.write_csv_atomic(OUT_DIR / "cladding" / "scan_summary.csv", ranked)
        return ranked

    for b0_nm in decimal_grid(
        CLADDING_B0_MIN_NM,
        CLADDING_B0_MAX_NM,
        CLADDING_B0_COARSE_STEP_NM,
        1,
    ):
        evaluate(b0_nm, CLADDING_INITIAL_ETA, CLADDING_INITIAL_ZETA, "b0_coarse")
    best = persist().iloc[0]
    for b0_nm in decimal_grid(
        max(CLADDING_B0_MIN_NM, float(best["b0_nm"]) - 1.0),
        min(CLADDING_B0_MAX_NM, float(best["b0_nm"]) + 1.0),
        CLADDING_B0_FINE_STEP_NM,
        1,
    ):
        evaluate(b0_nm, CLADDING_INITIAL_ETA, CLADDING_INITIAL_ZETA, "b0_fine")
    best = persist().iloc[0]

    for zeta in decimal_grid(
        CLADDING_ZETA_MIN,
        CLADDING_ZETA_MAX,
        CLADDING_ZETA_COARSE_STEP,
        3,
    ):
        evaluate(float(best["b0_nm"]), CLADDING_INITIAL_ETA, zeta, "zeta_coarse")
    best = persist().iloc[0]
    for zeta in decimal_grid(
        max(CLADDING_ZETA_MIN, float(best["zeta"]) - 0.005),
        min(CLADDING_ZETA_MAX, float(best["zeta"]) + 0.005),
        CLADDING_ZETA_FINE_STEP,
        3,
    ):
        evaluate(float(best["b0_nm"]), CLADDING_INITIAL_ETA, zeta, "zeta_fine")
    best = persist().iloc[0]

    for b0_nm in decimal_grid(
        max(CLADDING_B0_MIN_NM, float(best["b0_nm"]) - 0.5),
        min(CLADDING_B0_MAX_NM, float(best["b0_nm"]) + 0.5),
        CLADDING_B0_FINE_STEP_NM,
        1,
    ):
        evaluate(b0_nm, CLADDING_INITIAL_ETA, float(best["zeta"]), "b0_retune")
    best = persist().iloc[0]

    if not bool(best["within_alignment_tolerance"]):
        for eta in decimal_grid(
            CLADDING_ETA_MIN,
            CLADDING_ETA_MAX,
            CLADDING_ETA_COARSE_STEP,
            3,
        ):
            evaluate(float(best["b0_nm"]), eta, float(best["zeta"]), "eta_coarse")
        best = persist().iloc[0]
        for eta in decimal_grid(
            max(CLADDING_ETA_MIN, float(best["eta"]) - 0.002),
            min(CLADDING_ETA_MAX, float(best["eta"]) + 0.002),
            CLADDING_ETA_FINE_STEP,
            3,
        ):
            evaluate(float(best["b0_nm"]), eta, float(best["zeta"]), "eta_fine")
        best = persist().iloc[0]
        for b0_nm in decimal_grid(
            max(CLADDING_B0_MIN_NM, float(best["b0_nm"]) - 0.5),
            min(CLADDING_B0_MAX_NM, float(best["b0_nm"]) + 0.5),
            CLADDING_B0_FINE_STEP_NM,
            1,
        ):
            evaluate(b0_nm, float(best["eta"]), float(best["zeta"]), "eta_b0_retune")
        best = persist().iloc[0]

    if not bool(best["within_alignment_tolerance"]):
        raise RuntimeError(
            "No cladding candidate reaches the configured Gamma alignment tolerance"
        )
    result = best.drop(labels="rank").to_dict()
    best_dir = OUT_DIR / "cladding" / "best"
    best_dir.mkdir(parents=True, exist_ok=True)
    band.write_json_atomic(
        best_dir / "best_config.json",
        {
            "summary": result,
            "params": make_cladding_params(
                result["b0_nm"], result["eta"], result["zeta"]
            ),
            "mesh_auto_size": MESH_AUTO_SIZE,
        },
    )
    return result


def validate_best_pair(
    cavity_best: dict[str, object],
    cladding_best: dict[str, object],
) -> dict[str, object]:
    """Run the complete shared-path band workflow for the optimized cell pair."""
    root = OUT_DIR / "best_pair_full_band"
    cavity_params = make_cavity_params(cavity_best["eta"], cavity_best["zeta"])
    cladding_params = make_cladding_params(
        cladding_best["b0_nm"], cladding_best["eta"], cladding_best["zeta"]
    )
    result = band.run_band_workflow(root, cavity_params, cladding_params)
    cavity_gamma = result["cavity"][
        np.isclose(pd.to_numeric(result["cavity"]["k_norm"]), 0.0)
    ].drop_duplicates("band_label").set_index("mode_idx", drop=False)
    cladding_gamma = result["cladding"][
        np.isclose(pd.to_numeric(result["cladding"]["k_norm"]), 0.0)
    ].drop_duplicates("band_label").set_index("mode_idx", drop=False)
    cavity_selected = {
        row.band_label: int(row.mode_idx) for row in cavity_gamma.itertuples()
    }
    cladding_selected = {
        row.band_label: int(row.mode_idx) for row in cladding_gamma.itertuples()
    }
    _, _, cavity_mode, cavity_share = identify_user_p_mode(
        cavity_gamma, cavity_selected, "py"
    )
    _, _, cladding_mode, cladding_share = identify_user_p_mode(
        cladding_gamma, cladding_selected, "px"
    )
    validation = {
        "cavity_py_frequency_thz": float(cavity_mode["re"]),
        "cavity_py_q": float(cavity_mode["q"]),
        "cavity_py_share": cavity_share,
        "cladding_px_frequency_thz": float(cladding_mode["re"]),
        "cladding_px_share": cladding_share,
        "alignment_error_thz": abs(
            float(cladding_mode["re"]) - float(cavity_mode["re"])
        ),
        "mesh_auto_size": MESH_AUTO_SIZE,
        "output_dir": root,
    }
    config_dir = root / "99_config"
    config_dir.mkdir(parents=True, exist_ok=True)
    band.write_json_atomic(config_dir / "pair_validation.json", validation)
    return validation


def main() -> dict[str, object]:
    """Run the restartable two-stage optimization and full-band validation."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    band.write_json_atomic(
        OUT_DIR / "scan_config.json",
        {
            "mesh_auto_size": MESH_AUTO_SIZE,
            "cavity": {
                "b0_nm": CAVITY_B0_NM,
                "eta_range": [CAVITY_ETA_MIN, CAVITY_ETA_MAX],
                "zeta_range": [CAVITY_ZETA_MIN, CAVITY_ZETA_MAX],
                "coarse_step": CAVITY_COARSE_STEP,
                "fine_step": CAVITY_FINE_STEP,
                "reference_frequency_thz": CAVITY_REFERENCE_FREQUENCY_THZ,
                "frequency_tolerance_thz": CAVITY_FREQUENCY_TOLERANCE_THZ,
                "screen_k_magnitudes": SCREEN_K_MAGNITUDES,
                "final_k_magnitudes": FINAL_K_MAGNITUDES,
                "target_user_mode": "py",
                "color": "red",
                "use_accepted_candidate": USE_ACCEPTED_CAVITY_CANDIDATE,
                "accepted_eta": ACCEPTED_CAVITY_ETA,
                "accepted_zeta": ACCEPTED_CAVITY_ZETA,
            },
            "cladding": {
                "initial": {
                    "b0_nm": CLADDING_INITIAL_B0_NM,
                    "eta": CLADDING_INITIAL_ETA,
                    "zeta": CLADDING_INITIAL_ZETA,
                },
                "b0_range_nm": [CLADDING_B0_MIN_NM, CLADDING_B0_MAX_NM],
                "eta_range": [CLADDING_ETA_MIN, CLADDING_ETA_MAX],
                "zeta_range": [CLADDING_ZETA_MIN, CLADDING_ZETA_MAX],
                "alignment_tolerance_thz": CLADDING_ALIGNMENT_TOLERANCE_THZ,
                "target_user_mode": "px",
                "color": "blue",
            },
        },
    )
    if USE_ACCEPTED_CAVITY_CANDIDATE:
        cavity_best, _ = load_accepted_cavity_candidate()
    else:
        cavity_best, _ = optimize_cavity()
    cladding_best = optimize_cladding(float(cavity_best["gamma_frequency_thz"]))
    validation = validate_best_pair(cavity_best, cladding_best)
    result = {
        "cavity_best": cavity_best,
        "cladding_best": cladding_best,
        "validation": validation,
    }
    band.write_json_atomic(OUT_DIR / "optimization_result.json", result)
    return result


if __name__ == "__main__":
    main()
