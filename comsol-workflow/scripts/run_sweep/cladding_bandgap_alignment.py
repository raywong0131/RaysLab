#!/usr/bin/env python3
"""Align project-internal p_y (user-facing p_x) and validate the p-d bandgap."""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import UNIT_CELL_OUTPUT_ROOT  # noqa: E402


B0_MIN_NM = 230.0
B0_MAX_NM = 240.0
B0_COARSE_STEP_NM = 1.0
B0_FINE_STEP_NM = 0.1
ZETA_MIN = 0.830
ZETA_MAX = 0.850
ZETA_COARSE_STEP = 0.005
ZETA_FINE_STEP = 0.001
ETA_MIN = 0.950
ETA_MAX = 0.970
ETA_COARSE_STEP = 0.005
ETA_FINE_STEP = 0.001
FINE_HALF_WIDTH = 0.005
INITIAL_B0_NM = 235.0
INITIAL_ETA = 0.960
INITIAL_ZETA = 0.840
ALIGNMENT_TOLERANCE_THZ = 0.05
MIN_GAMMA_GAP_THZ = 2.0
FOURIER_REFERENCE_B0_NM = 230.0
OUT_DIR = UNIT_CELL_OUTPUT_ROOT / "cladding_bandgap_alignment"
GAMMA_POINT = {
    "direction": "gamma_m",
    "kx": 0.0,
    "ky": 0.0,
    "kx_str": "0.000000",
    "ky_str": "0.000000",
    "k_norm": 0.0,
    "path_fraction": 0.0,
}


def decimal_grid(
    start: float,
    stop: float,
    step: float,
    digits: int,
) -> tuple[float, ...]:
    """Return an inclusive decimal grid without floating endpoint drift."""
    scale = 10**digits
    start_integer = round(float(start) * scale)
    stop_integer = round(float(stop) * scale)
    step_integer = round(float(step) * scale)
    return tuple(value / scale for value in range(start_integer, stop_integer + 1, step_integer))


def coarse_b0_values() -> tuple[float, ...]:
    return decimal_grid(B0_MIN_NM, B0_MAX_NM, B0_COARSE_STEP_NM, 1)


def coarse_zeta_values() -> tuple[float, ...]:
    return decimal_grid(ZETA_MIN, ZETA_MAX, ZETA_COARSE_STEP, 3)


def fine_zeta_values(center: float) -> tuple[float, ...]:
    low = max(ZETA_MIN, round(float(center) - FINE_HALF_WIDTH, 3))
    high = min(ZETA_MAX, round(float(center) + FINE_HALF_WIDTH, 3))
    return decimal_grid(low, high, ZETA_FINE_STEP, 3)


def coarse_eta_values() -> tuple[float, ...]:
    return decimal_grid(ETA_MIN, ETA_MAX, ETA_COARSE_STEP, 3)


def fine_eta_values(center: float) -> tuple[float, ...]:
    low = max(ETA_MIN, round(float(center) - FINE_HALF_WIDTH, 3))
    high = min(ETA_MAX, round(float(center) + FINE_HALF_WIDTH, 3))
    return decimal_grid(low, high, ETA_FINE_STEP, 3)


def make_cladding_params(
    b0_nm: float,
    eta: float,
    zeta: float,
) -> dict[str, float | str]:
    """Translate the user-facing cladding geometry parameters."""
    return {
        "name": "cladding_bandgap_alignment",
        "r_f0": float(eta),
        "b_square_f0": (float(b0_nm) / FOURIER_REFERENCE_B0_NM) ** 2,
        "b_square_f3": (float(zeta) ** 2 - 1.0) / 2.0,
    }


def _component_share(first: float, second: float) -> float:
    total = float(first) + float(second)
    if not math.isfinite(total) or total <= 0.0:
        return float("nan")
    return float(second) / total


def identify_py_mode(
    gamma_modes: pd.DataFrame,
    selected_modes: dict[str, int],
) -> tuple[int, pd.Series]:
    """Return project-internal p_y, which is user-facing p_x."""
    candidates = []
    for label in ("p1", "p2"):
        mode_idx = int(selected_modes[label])
        row = gamma_modes.loc[mode_idx]
        share = _component_share(row["px_weight"], row["py_weight"])
        candidates.append((share, mode_idx, row))
    finite = [candidate for candidate in candidates if math.isfinite(candidate[0])]
    if not finite:
        raise ValueError("Selected p modes have no finite px/py composition")
    _share, mode_idx, row = max(finite, key=lambda candidate: candidate[0])
    return mode_idx, row


def summarize_candidate(
    b0_nm: float,
    eta: float,
    zeta: float,
    stage: str,
    gamma_modes: pd.DataFrame,
    selected_modes: dict[str, int],
    cavity_py_frequency_thz: float,
) -> dict[str, object]:
    """Calculate alignment and frequency-defined Gamma p-d gap metrics."""
    result = {
        "b0_nm": round(float(b0_nm), 1),
        "eta": round(float(eta), 3),
        "zeta": round(float(zeta), 3),
        "stage": str(stage),
        "status": "invalid_modes",
        "feasible": False,
        "within_alignment_tolerance": False,
        "cavity_py_frequency_thz": float(cavity_py_frequency_thz),
        "cladding_py_mode_idx": np.nan,
        "cladding_py_frequency_thz": np.nan,
        "cladding_py_share": np.nan,
        "signed_alignment_error_thz": np.nan,
        "alignment_error_thz": np.nan,
        "gap_lower_edge_thz": np.nan,
        "gap_upper_edge_thz": np.nan,
        "gap_width_thz": np.nan,
        "midgap_thz": np.nan,
        "midgap_error_thz": np.nan,
    }
    try:
        p_indices = [int(selected_modes[label]) for label in ("p1", "p2")]
        d_indices = [int(selected_modes[label]) for label in ("d1", "d2")]
        chosen_indices = p_indices + d_indices
        rows = gamma_modes.loc[chosen_indices]
        if not rows["is_valid"].astype(bool).all():
            return result
        frequencies = pd.to_numeric(rows["re"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(frequencies).all():
            result["status"] = "nonfinite_frequency"
            return result
        if any(gamma_modes.loc[index, "dominant_subspace"] != "p" for index in p_indices):
            result["status"] = "invalid_p_subspace"
            return result
        if any(gamma_modes.loc[index, "dominant_subspace"] != "d" for index in d_indices):
            result["status"] = "invalid_d_subspace"
            return result
        py_mode_idx, py_row = identify_py_mode(gamma_modes, selected_modes)
    except (KeyError, TypeError, ValueError):
        return result

    py_frequency = float(py_row["re"])
    py_share = _component_share(py_row["px_weight"], py_row["py_weight"])
    lower_edge = max(float(gamma_modes.loc[index, "re"]) for index in p_indices)
    upper_edge = min(float(gamma_modes.loc[index, "re"]) for index in d_indices)
    gap_width = upper_edge - lower_edge
    midgap = 0.5 * (lower_edge + upper_edge)
    signed_alignment = py_frequency - float(cavity_py_frequency_thz)
    alignment_error = abs(signed_alignment)
    midgap_error = abs(midgap - float(cavity_py_frequency_thz))
    within_tolerance = alignment_error <= ALIGNMENT_TOLERANCE_THZ + 1e-12
    feasible = gap_width >= MIN_GAMMA_GAP_THZ - 1e-12
    result.update(
        {
            "status": (
                "ok_aligned"
                if feasible and within_tolerance
                else "ok_unaligned"
                if feasible
                else "gap_below_minimum"
            ),
            "feasible": bool(feasible),
            "within_alignment_tolerance": bool(feasible and within_tolerance),
            "cladding_py_mode_idx": int(py_mode_idx),
            "cladding_py_frequency_thz": py_frequency,
            "cladding_py_share": py_share,
            "signed_alignment_error_thz": signed_alignment,
            "alignment_error_thz": alignment_error,
            "gap_lower_edge_thz": lower_edge,
            "gap_upper_edge_thz": upper_edge,
            "gap_width_thz": gap_width,
            "midgap_thz": midgap,
            "midgap_error_thz": midgap_error,
        }
    )
    return result


def rank_candidates(summary: pd.DataFrame) -> pd.DataFrame:
    """Rank candidates by alignment tier, midgap error, and parameter priority."""
    if summary.empty:
        return summary.copy()
    ranked = summary.copy()
    feasible = ranked["feasible"].fillna(False).astype(bool)
    aligned = feasible & ranked["within_alignment_tolerance"].fillna(False).astype(bool)
    ranked["ranking_tier"] = np.select([aligned, feasible], [0, 1], default=2)
    ranked["ranking_primary"] = np.where(
        aligned,
        pd.to_numeric(ranked["midgap_error_thz"], errors="coerce"),
        pd.to_numeric(ranked["alignment_error_thz"], errors="coerce"),
    )
    ranked["ranking_secondary"] = np.where(
        aligned,
        pd.to_numeric(ranked["alignment_error_thz"], errors="coerce"),
        pd.to_numeric(ranked["midgap_error_thz"], errors="coerce"),
    )
    ranked[["ranking_primary", "ranking_secondary"]] = ranked[
        ["ranking_primary", "ranking_secondary"]
    ].fillna(np.inf)
    ranked["stage_order"] = ranked["stage"].map(
        {
            "b0": 0,
            "zeta_coarse": 1,
            "zeta_fine": 2,
            "eta_coarse": 3,
            "eta_fine": 4,
        }
    ).fillna(99)
    ranked["parameter_distance"] = (
        (ranked["b0_nm"] - INITIAL_B0_NM).abs() / (B0_MAX_NM - B0_MIN_NM)
        + (ranked["zeta"] - INITIAL_ZETA).abs() / (ZETA_MAX - ZETA_MIN)
        + (ranked["eta"] - INITIAL_ETA).abs() / (ETA_MAX - ETA_MIN)
    )
    ranked = ranked.sort_values(
        [
            "ranking_tier",
            "ranking_primary",
            "ranking_secondary",
            "stage_order",
            "parameter_distance",
            "b0_nm",
            "zeta",
            "eta",
        ],
        kind="mergesort",
    ).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1, dtype=int)
    return ranked


def scan_root() -> Path:
    """Return the stable output root for the fixed cavity reference."""
    return Path(OUT_DIR) / "cav_235_0.960_1.156"


def candidate_dir(b0_nm: float, eta: float, zeta: float) -> Path:
    return (
        scan_root()
        / "candidates"
        / (
            f"clad_b0_{float(b0_nm):.1f}_eta_{float(eta):.3f}_"
            f"zeta_{float(zeta):.3f}"
        )
    )


def _simulation_config():
    return band.SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift=band.EIGENFREQUENCY_SHIFT,
        eigenmode_count=band.EIGENMODE_COUNT,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )


def run_gamma_case(
    case_dir: Path,
    params: dict[str, float | str],
) -> tuple[pd.DataFrame, dict[str, int], dict[str, object]]:
    """Run or load one Gamma point and return its selected p/d modes."""
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    hole_params = band.get_hole_params(params)
    hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            "Unit-cell geometry violates clearance: "
            f"min_dist={minimum_clearance:g}, d={band.D:g}"
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
    frame = band.run_k_point(
        case_dir,
        [hole.tolist() for hole in holes],
        GAMMA_POINT,
        _simulation_config(),
    )
    selected_modes = band.select_gamma_bands(frame)
    metadata = {
        "hole_params": hole_params,
        "triangles": holes,
        "minimum_clearance": minimum_clearance,
    }
    return frame, selected_modes, metadata


def load_cavity_reference() -> dict[str, object]:
    """Load the fixed cavity Gamma point and identify user-facing p_x."""
    root = scan_root() / "cavity_reference"
    frame, selected_modes, metadata = run_gamma_case(root, band.BULK_PARAMS)
    mode_idx, mode = identify_py_mode(frame, selected_modes)
    reference = {
        "py_mode_idx": int(mode_idx),
        "py_frequency_thz": float(mode["re"]),
        "py_share": _component_share(mode["px_weight"], mode["py_weight"]),
        "selected_gamma_modes": selected_modes,
        **metadata,
    }
    band.write_json_atomic(
        root / "config.json",
        {"params": band.BULK_PARAMS, "reference": reference},
    )
    return reference


def evaluate_candidate(
    b0_nm: float,
    eta: float,
    zeta: float,
    stage: str,
    cavity_reference: dict[str, object],
) -> dict[str, object]:
    """Run or reuse one cladding Gamma candidate and persist its metrics."""
    b0_nm = round(float(b0_nm), 1)
    eta = round(float(eta), 3)
    zeta = round(float(zeta), 3)
    root = candidate_dir(b0_nm, eta, zeta)
    params = make_cladding_params(b0_nm, eta, zeta)
    try:
        frame, selected_modes, metadata = run_gamma_case(root, params)
        summary = summarize_candidate(
            b0_nm,
            eta,
            zeta,
            stage,
            frame,
            selected_modes,
            float(cavity_reference["py_frequency_thz"]),
        )
        error = ""
    except Exception as exc:  # keep other candidates restartable after one failure
        selected_modes = {}
        metadata = {}
        error = f"{type(exc).__name__}: {exc}"
        summary = {
            "b0_nm": b0_nm,
            "eta": eta,
            "zeta": zeta,
            "stage": str(stage),
            "status": "candidate_error",
            "feasible": False,
            "within_alignment_tolerance": False,
            "cavity_py_frequency_thz": float(cavity_reference["py_frequency_thz"]),
            "cladding_py_mode_idx": np.nan,
            "cladding_py_frequency_thz": np.nan,
            "cladding_py_share": np.nan,
            "signed_alignment_error_thz": np.nan,
            "alignment_error_thz": np.nan,
            "gap_lower_edge_thz": np.nan,
            "gap_upper_edge_thz": np.nan,
            "gap_width_thz": np.nan,
            "midgap_thz": np.nan,
            "midgap_error_thz": np.nan,
        }
    band.write_json_atomic(
        root / "config.json",
        {
            "params": params,
            "stage": stage,
            "selected_gamma_modes": selected_modes,
            "summary": summary,
            "error": error,
            **metadata,
        },
    )
    return summary


def plot_optimization_trace(summary: pd.DataFrame, output_path: Path) -> None:
    """Plot the three search metrics against their controlling parameters."""
    ordered = summary.sort_values("evaluation_order")
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    panels = (
        (
            "b0_nm",
            "alignment_error_thz",
            r"Cladding $b_0$ (nm)",
            r"$|\Delta f_{p_x}|$ (THz)",
            r"$p_x$ alignment",
        ),
        (
            "zeta",
            "midgap_error_thz",
            r"Cladding $\zeta$",
            "Midgap error (THz)",
            "Midgap offset",
        ),
        (
            "eta",
            "gap_width_thz",
            r"Cladding $\eta$",
            r"$\Gamma$ p-d gap (THz)",
            r"$\Gamma$ gap width",
        ),
    )
    stages = list(dict.fromkeys(ordered["stage"].astype(str)))
    for axis, (x_column, y_column, x_label, y_label, title) in zip(axes, panels):
        for stage in stages:
            rows = ordered[ordered["stage"].astype(str) == stage]
            axis.scatter(
                pd.to_numeric(rows[x_column], errors="coerce"),
                pd.to_numeric(rows[y_column], errors="coerce"),
                s=15,
                alpha=0.75,
                label=stage,
            )
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
    axes[0].axhline(
        ALIGNMENT_TOLERANCE_THZ,
        color="0.35",
        linestyle="--",
        linewidth=0.9,
        label="alignment tolerance",
    )
    axes[2].axhline(
        MIN_GAMMA_GAP_THZ,
        color="0.35",
        linestyle="--",
        linewidth=0.9,
        label="minimum gap",
    )
    legend_items = {}
    for axis in axes:
        handles, labels = axis.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            legend_items.setdefault(label, handle)
    figure.legend(
        legend_items.values(),
        legend_items.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=min(5, len(legend_items)),
        fontsize=7,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=250)
    plt.close(figure)


def persist_summary(
    records: dict[tuple[float, float, float], dict[str, object]],
    root: Path,
    plot: bool = True,
) -> pd.DataFrame:
    """Persist the current deterministic candidate ranking."""
    ranked = rank_candidates(pd.DataFrame(records.values()))
    band.write_csv_atomic(Path(root) / "scan_summary.csv", ranked)
    if plot:
        plot_optimization_trace(ranked, Path(root) / "optimization_trace.png")
    return ranked


def retune_b0(
    eta: float,
    zeta: float,
    stage: str,
    evaluate,
    seed_b0: float,
) -> dict[str, object]:
    """Tune b0 on its discrete grid, using signed-error bracketing when possible."""
    seed = min(B0_MAX_NM, max(B0_MIN_NM, round(float(seed_b0), 1)))
    rows = [evaluate(seed, eta, zeta, stage)]

    def valid_rows():
        return sorted(
            [row for row in rows if bool(row["feasible"])],
            key=lambda row: float(row["b0_nm"]),
        )

    def find_crossings():
        crossings = []
        ordered = valid_rows()
        for left, right in zip(ordered[:-1], ordered[1:]):
            left_error = float(left["signed_alignment_error_thz"])
            right_error = float(right["signed_alignment_error_thz"])
            if math.isfinite(left_error) and math.isfinite(right_error):
                if left_error * right_error <= 0.0:
                    crossings.append((left, right))
        return crossings

    if not bool(rows[0].get("within_alignment_tolerance", False)):
        for distance in range(1, round(B0_MAX_NM - B0_MIN_NM) + 1):
            values = {
                round(seed - distance * B0_COARSE_STEP_NM, 1),
                round(seed + distance * B0_COARSE_STEP_NM, 1),
            }
            for value in sorted(values):
                if B0_MIN_NM <= value <= B0_MAX_NM:
                    rows.append(evaluate(value, eta, zeta, stage))
            if find_crossings() or any(
                bool(row.get("within_alignment_tolerance", False)) for row in rows
            ):
                break

    crossings = find_crossings()

    if crossings:
        left, right = min(
            crossings,
            key=lambda pair: min(
                float(pair[0]["alignment_error_thz"]),
                float(pair[1]["alignment_error_thz"]),
            ),
        )
        lo = round(float(left["b0_nm"]) * 10)
        hi = round(float(right["b0_nm"]) * 10)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            row = evaluate(mid / 10.0, eta, zeta, stage)
            rows.append(row)
            if (
                float(left["signed_alignment_error_thz"])
                * float(row["signed_alignment_error_thz"])
                <= 0.0
            ):
                right, hi = row, mid
            else:
                left, lo = row, mid
        rows.extend(
            [
                evaluate(lo / 10.0, eta, zeta, stage),
                evaluate(hi / 10.0, eta, zeta, stage),
            ]
        )
    elif rows:
        best = rank_candidates(pd.DataFrame(rows)).iloc[0]
        low = max(B0_MIN_NM, float(best["b0_nm"]) - B0_FINE_STEP_NM)
        high = min(B0_MAX_NM, float(best["b0_nm"]) + B0_FINE_STEP_NM)
        rows.extend(
            evaluate(value, eta, zeta, stage)
            for value in decimal_grid(low, high, B0_FINE_STEP_NM, 1)
        )
    else:
        raise RuntimeError("b0 search produced no candidates")
    return rank_candidates(pd.DataFrame(rows)).iloc[0].to_dict()


def optimize_zeta(
    eta: float,
    evaluate,
    seed_b0: float,
    coarse_stage: str,
    fine_stage: str,
) -> dict[str, object]:
    """Optimize zeta while retuning b0 for every zeta value."""
    coarse_best = [
        retune_b0(eta, zeta, coarse_stage, evaluate, seed_b0)
        for zeta in coarse_zeta_values()
    ]
    best_coarse = rank_candidates(pd.DataFrame(coarse_best)).iloc[0].to_dict()
    fine_best = [
        retune_b0(
            eta,
            zeta,
            fine_stage,
            evaluate,
            float(best_coarse["b0_nm"]),
        )
        for zeta in fine_zeta_values(float(best_coarse["zeta"]))
    ]
    return rank_candidates(pd.DataFrame(coarse_best + fine_best)).iloc[0].to_dict()


def run_search(
    evaluate=None,
    cavity_reference: dict[str, object] | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run the restartable b0, zeta, and conditional eta search."""
    root = scan_root()
    root.mkdir(parents=True, exist_ok=True)
    reference = load_cavity_reference() if cavity_reference is None else cavity_reference
    evaluator = evaluate_candidate if evaluate is None else evaluate
    records: dict[tuple[float, float, float], dict[str, object]] = {}

    band.write_json_atomic(
        root / "scan_config.json",
        {
            "cavity_reference": reference,
            "b0_range_nm": [B0_MIN_NM, B0_MAX_NM],
            "b0_steps_nm": [B0_COARSE_STEP_NM, B0_FINE_STEP_NM],
            "zeta_range": [ZETA_MIN, ZETA_MAX],
            "zeta_steps": [ZETA_COARSE_STEP, ZETA_FINE_STEP],
            "eta_range": [ETA_MIN, ETA_MAX],
            "eta_steps": [ETA_COARSE_STEP, ETA_FINE_STEP],
            "fine_half_width": FINE_HALF_WIDTH,
            "alignment_tolerance_thz": ALIGNMENT_TOLERANCE_THZ,
            "minimum_gamma_gap_thz": MIN_GAMMA_GAP_THZ,
        },
    )

    def cached_evaluate(b0_nm, eta, zeta, stage):
        key = (round(float(b0_nm), 1), round(float(eta), 3), round(float(zeta), 3))
        if key not in records:
            row = dict(evaluator(*key, stage, reference))
            row["evaluation_order"] = len(records) + 1
            records[key] = row
        return records[key]

    stage_one = retune_b0(
        INITIAL_ETA,
        INITIAL_ZETA,
        "b0",
        cached_evaluate,
        INITIAL_B0_NM,
    )
    persist_summary(records, root)
    stage_two = optimize_zeta(
        INITIAL_ETA,
        cached_evaluate,
        float(stage_one["b0_nm"]),
        "zeta_coarse",
        "zeta_fine",
    )
    ranked = persist_summary(records, root)
    aligned = ranked[
        ranked["feasible"].astype(bool)
        & ranked["within_alignment_tolerance"].astype(bool)
    ]

    if aligned.empty:
        eta_coarse_best = []
        for eta in coarse_eta_values():
            eta_coarse_best.append(
                optimize_zeta(
                    eta,
                    cached_evaluate,
                    float(stage_two["b0_nm"]),
                    "eta_coarse",
                    "eta_coarse",
                )
            )
            persist_summary(records, root, plot=False)
        best_eta_coarse = rank_candidates(pd.DataFrame(eta_coarse_best)).iloc[0]
        for eta in fine_eta_values(float(best_eta_coarse["eta"])):
            optimize_zeta(
                eta,
                cached_evaluate,
                float(best_eta_coarse["b0_nm"]),
                "eta_fine",
                "eta_fine",
            )
            persist_summary(records, root, plot=False)
        ranked = persist_summary(records, root)

    feasible = ranked[ranked["feasible"].astype(bool)]
    if feasible.empty:
        raise RuntimeError("No cladding candidate preserves a valid Gamma p-d gap")
    best = ranked.iloc[0].to_dict()
    return ranked, best


def _gamma_mode_table(
    selected_bands: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Reduce the duplicated directional Gamma rows to one row per band label."""
    matched = selected_bands[selected_bands["match_status"] != "unmatched"].copy()
    gamma = matched[
        np.isclose(
            pd.to_numeric(matched["path_coordinate"], errors="coerce"),
            0.0,
        )
    ]
    gamma = gamma.drop_duplicates("band_label", keep="first").set_index("band_label")
    labels = ("p1", "p2", "d1", "d2")
    missing = [label for label in labels if label not in gamma.index]
    if missing:
        raise RuntimeError(f"Full-band Gamma rows are missing labels: {missing}")
    table = gamma.loc[list(labels)].reset_index(drop=True)
    return table, {label: index for index, label in enumerate(labels)}


def validate_full_band(
    best: dict[str, object],
    cavity_reference: dict[str, object],
) -> dict[str, object]:
    """Run the selected full band and validate Gamma and path-wide p-d gaps."""
    best_dir = scan_root() / "best"
    full_band_dir = best_dir / "full_band"
    b0_nm = float(best["b0_nm"])
    eta = float(best["eta"])
    zeta = float(best["zeta"])
    cladding_params = make_cladding_params(b0_nm, eta, zeta)
    band.write_json_atomic(
        best_dir / "best_config.json",
        {
            "summary": best,
            "cavity_reference": cavity_reference,
            "cladding_params": cladding_params,
        },
    )
    band.write_csv_atomic(best_dir / "best_gamma_summary.csv", pd.DataFrame([best]))

    result = band.run_band_workflow(
        full_band_dir,
        band.BULK_PARAMS,
        cladding_params,
    )
    cavity = result["cavity"]
    cladding = result["cladding"]
    if (cavity["match_status"] == "unmatched").any() or (
        cladding["match_status"] == "unmatched"
    ).any():
        raise RuntimeError("Full-band validation contains unmatched modes")

    cavity_gamma, cavity_selected = _gamma_mode_table(cavity)
    cavity_py_idx, cavity_py = identify_py_mode(cavity_gamma, cavity_selected)
    cladding_gamma, cladding_selected = _gamma_mode_table(cladding)
    gamma_summary = summarize_candidate(
        b0_nm,
        eta,
        zeta,
        "full_band",
        cladding_gamma,
        cladding_selected,
        float(cavity_py["re"]),
    )
    gamma_summary["cavity_py_mode_idx"] = int(cavity_py_idx)
    gamma_summary["cavity_reference_difference_thz"] = abs(
        float(cavity_py["re"]) - float(cavity_reference["py_frequency_thz"])
    )

    matched = cladding[cladding["match_status"] != "unmatched"]
    p_frequencies = pd.to_numeric(
        matched.loc[matched["band_label"].str.startswith("p"), "re"],
        errors="coerce",
    )
    d_frequencies = pd.to_numeric(
        matched.loc[matched["band_label"].str.startswith("d"), "re"],
        errors="coerce",
    )
    full_path_gap = float(d_frequencies.min() - p_frequencies.max())
    if not math.isfinite(full_path_gap) or full_path_gap <= 0.0:
        raise RuntimeError("Full-path cladding p-d gap closed")

    cavity_points = cavity[["direction", "path_coordinate"]].drop_duplicates()
    cladding_points = cladding[["direction", "path_coordinate"]].drop_duplicates()
    validation = {
        "gamma_summary": gamma_summary,
        "full_path_gap_thz": full_path_gap,
        "cavity_sample_count": int(len(cavity_points)),
        "cladding_sample_count": int(len(cladding_points)),
        "sampling_counts_match": bool(len(cavity_points) == len(cladding_points)),
        "full_band_dir": full_band_dir,
    }
    band.write_json_atomic(best_dir / "validation.json", validation)
    return validation


def main() -> dict[str, object]:
    """Run the Gamma search and the final full-band validation."""
    cavity_reference = load_cavity_reference()
    summary, best = run_search(cavity_reference=cavity_reference)
    validation = validate_full_band(best, cavity_reference)
    return {
        "best": best,
        "validation": validation,
        "candidate_count": int(len(summary)),
    }


if __name__ == "__main__":
    main()
