#!/usr/bin/env python3
"""Scan cavity-cell zeta for a p2 Q maximum at Gamma."""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import UNIT_CELL_OUTPUT_ROOT  # noqa: E402


B0_NM = 235.0
ETA = 0.96
FOURIER_REFERENCE_B0_NM = 230.0
ZETA_MIN = 1.100
ZETA_MAX = 1.200
COARSE_ZETA_STEP = 0.010
FINE_ZETA_STEP = 0.001
FINE_HALF_WIDTH = 0.010
K_STRING_PRECISION = 6
CLADDING_GAP_MIN_THZ = 192.478192022
CLADDING_GAP_MAX_THZ = 197.269008408
MIN_GAMMA_P_WEIGHT = 0.9
OUT_DIR = UNIT_CELL_OUTPUT_ROOT / "unit_cell_p2_gamma_scan"
SCREEN_K_MAGNITUDES = (
    0.000,
    0.002,
    0.004,
    0.006,
    0.010,
    0.015,
    0.020,
    0.030,
    0.050,
    0.075,
    0.100,
    0.125,
    0.150,
    0.175,
    0.200,
)
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
FINAL_CANDIDATE_COUNT = 3


def make_cavity_params(zeta: float) -> dict[str, float | str]:
    """Translate ``(b0, eta, zeta)`` into the existing Fourier parameters."""
    zeta = float(zeta)
    if not ZETA_MIN <= zeta <= ZETA_MAX:
        raise ValueError(f"zeta must be in [{ZETA_MIN:.3f}, {ZETA_MAX:.3f}]")
    return {
        "name": "cavity_p_bic_zeta_scan",
        "r_f0": ETA,
        "b_square_f0": (B0_NM / FOURIER_REFERENCE_B0_NM) ** 2,
        "b_square_f3": (zeta**2 - 1.0) / 2.0,
    }


def coarse_zeta_values() -> tuple[float, ...]:
    """Return the inclusive three-decimal coarse scan grid."""
    start = round(ZETA_MIN * 1000)
    stop = round(ZETA_MAX * 1000)
    step = round(COARSE_ZETA_STEP * 1000)
    return tuple(value / 1000.0 for value in range(start, stop + 1, step))


def fine_zeta_values(center: float) -> tuple[float, ...]:
    """Return an inclusive, globally clipped three-decimal refinement grid."""
    center_milli = round(float(center) * 1000)
    low = max(round(ZETA_MIN * 1000), center_milli - round(FINE_HALF_WIDTH * 1000))
    high = min(round(ZETA_MAX * 1000), center_milli + round(FINE_HALF_WIDTH * 1000))
    step = round(FINE_ZETA_STEP * 1000)
    return tuple(value / 1000.0 for value in range(low, high + 1, step))


def make_direction_points(
    direction: str,
    magnitudes: Sequence[float],
) -> list[dict[str, float | str]]:
    """Build the existing unit-cell k-point schema on one Gamma branch."""
    endpoints = {"gamma_m": 0.5, "gamma_k": 1.0 / math.sqrt(3.0)}
    if direction not in endpoints:
        raise ValueError(f"Unknown direction: {direction}")
    endpoint = endpoints[direction]
    points = []
    for raw_magnitude in magnitudes:
        magnitude = round(float(raw_magnitude), K_STRING_PRECISION)
        if not 0.0 <= magnitude <= endpoint:
            raise ValueError(f"k magnitude {magnitude:g} is outside {direction}")
        kx = magnitude if direction == "gamma_k" else 0.0
        ky = magnitude if direction == "gamma_m" else 0.0
        sign = -1.0 if direction == "gamma_m" else 1.0
        points.append(
            {
                "direction": direction,
                "kx": kx,
                "ky": ky,
                "kx_str": f"{kx:.{K_STRING_PRECISION}f}",
                "ky_str": f"{ky:.{K_STRING_PRECISION}f}",
                "k_norm": magnitude,
                "path_fraction": magnitude / endpoint,
                "path_coordinate": sign * magnitude,
            }
        )
    return points


def _candidate_base(zeta: float, stage: str) -> dict[str, object]:
    return {
        "zeta": round(float(zeta), 3),
        "stage": str(stage),
        "status": "missing_p2",
        "feasible": False,
        "peak_at_gamma": False,
        "gamma_frequency_thz": np.nan,
        "gamma_q": np.nan,
        "gamma_p_weight": np.nan,
        "peak_excess_dex": np.nan,
        "peak_k_signed": np.nan,
        "peak_direction": "",
        "sampled_point_count": 0,
    }


def summarize_p2_candidate(
    zeta: float,
    p2_rows: pd.DataFrame,
    stage: str,
) -> dict[str, object]:
    """Apply hard constraints and summarize the sampled p2 Q maximum."""
    result = _candidate_base(zeta, stage)
    required = {
        "direction",
        "band_label",
        "k_norm",
        "q",
        "re",
        "p_weight",
        "is_valid",
        "match_status",
    }
    if not required.issubset(p2_rows.columns):
        return result
    rows = p2_rows[p2_rows["band_label"] == "p2"].copy()
    result["sampled_point_count"] = int(len(rows))
    if rows.empty:
        return result
    if (rows["match_status"] == "unmatched").any():
        result["status"] = "unmatched_p2"
        return result

    gamma = rows[np.isclose(pd.to_numeric(rows["k_norm"], errors="coerce"), 0.0)]
    if gamma.empty:
        result["status"] = "missing_gamma_p2"
        return result
    gamma_row = gamma.iloc[0]
    result.update(
        {
            "gamma_frequency_thz": float(gamma_row["re"]),
            "gamma_q": float(gamma_row["q"]),
            "gamma_p_weight": float(gamma_row["p_weight"]),
        }
    )
    if not bool(gamma_row["is_valid"]):
        result["status"] = "invalid_gamma"
        return result
    if not CLADDING_GAP_MIN_THZ <= result["gamma_frequency_thz"] <= CLADDING_GAP_MAX_THZ:
        result["status"] = "gamma_frequency_outside_gap"
        return result
    if result["gamma_p_weight"] < MIN_GAMMA_P_WEIGHT:
        result["status"] = "gamma_p_weight_below_min"
        return result

    q_values = pd.to_numeric(rows["q"], errors="coerce").to_numpy(dtype=float)
    if not np.all(np.isfinite(q_values) & (q_values > 0.0)):
        result["status"] = "invalid_q"
        return result
    gamma_q = float(result["gamma_q"])
    nonzero = rows[~np.isclose(pd.to_numeric(rows["k_norm"]), 0.0)].copy()
    if nonzero.empty:
        result["status"] = "missing_nonzero_k"
        return result
    nonzero["peak_excess_dex"] = np.log10(nonzero["q"].astype(float)) - math.log10(
        gamma_q
    )
    peak_row = nonzero.loc[nonzero["peak_excess_dex"].idxmax()]
    peak_excess = float(peak_row["peak_excess_dex"])
    peak_at_gamma = peak_excess <= 0.0
    sign = -1.0 if peak_row["direction"] == "gamma_m" else 1.0
    result.update(
        {
            "status": "ok",
            "feasible": True,
            "peak_at_gamma": bool(peak_at_gamma),
            "peak_excess_dex": peak_excess,
            "peak_k_signed": 0.0 if peak_at_gamma else sign * float(peak_row["k_norm"]),
            "peak_direction": "gamma" if peak_at_gamma else str(peak_row["direction"]),
        }
    )
    return result


def rank_candidates(summary: pd.DataFrame) -> pd.DataFrame:
    """Rank feasible Gamma peaks first, then nearest off-Gamma candidates."""
    if summary.empty:
        return summary.assign(rank=pd.Series(dtype=int))
    ranked = summary.copy()
    feasible = ranked["feasible"].fillna(False).astype(bool)
    at_gamma = ranked["peak_at_gamma"].fillna(False).astype(bool)
    ranked["_category"] = np.where(feasible & at_gamma, 0, np.where(feasible, 1, 2))
    ranked["_primary"] = np.where(
        ranked["_category"] == 0,
        -pd.to_numeric(ranked["gamma_q"], errors="coerce").fillna(-np.inf),
        pd.to_numeric(ranked["peak_excess_dex"], errors="coerce").fillna(np.inf),
    )
    ranked["_secondary"] = -pd.to_numeric(
        ranked["gamma_q"], errors="coerce"
    ).fillna(-np.inf)
    ranked = ranked.sort_values(
        ["_category", "_primary", "_secondary", "zeta"], kind="mergesort"
    ).drop(columns=["_category", "_primary", "_secondary"])
    ranked = ranked.reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def scan_root() -> Path:
    """Return the stable root for the fixed b0/eta scan."""
    return Path(OUT_DIR) / f"b0_{B0_NM:g}_eta_{ETA:.3f}"


def candidate_dir(zeta: float) -> Path:
    return scan_root() / f"zeta_{float(zeta):.3f}"


def _unique_points(
    branches: dict[str, list[dict[str, float | str]]],
) -> list[dict[str, float | str]]:
    unique = {}
    for direction in ("gamma_m", "gamma_k"):
        for point in branches[direction]:
            unique.setdefault((point["kx_str"], point["ky_str"]), point)
    return list(unique.values())


def _simulation_config():
    return band.SimulationConfig(
        slab_height=band.SLAB_HEIGHT,
        slab_refractive_index=band.REFRACTIVE_INDEX,
        eigenfrequency_shift=band.EIGENFREQUENCY_SHIFT,
        eigenmode_count=band.EIGENMODE_COUNT,
        mode_type=band.MODE_TYPE,
        wavelength=band.WAVELENGTH,
    )


def run_candidate(
    zeta: float,
    stage: str,
    magnitudes: Sequence[float] = SCREEN_K_MAGNITUDES,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run or resume one zeta candidate and return tracked p2 rows plus summary."""
    params = make_cavity_params(zeta)
    root = candidate_dir(zeta)
    root.mkdir(parents=True, exist_ok=True)
    hole_params = band.get_hole_params(params)
    hexagon, holes, info = band.create_hexagon_design(band.A, hole_params)
    minimum_clearance = float(min(item["min_dist"] for item in info))
    if minimum_clearance < band.D:
        raise ValueError(
            f"zeta={zeta:.3f} violates clearance: "
            f"min_dist={minimum_clearance:g}, d={band.D:g}"
        )
    geometry_path = root / "geometry.png"
    if not geometry_path.exists():
        band.visualize_hexagon_design(
            hexagon,
            holes,
            info,
            filename=geometry_path,
            annotate=False,
        )

    branches = {
        direction: make_direction_points(direction, magnitudes)
        for direction in ("gamma_m", "gamma_k")
    }
    hole_vertices = [hole.tolist() for hole in holes]
    config = _simulation_config()
    unique_points = _unique_points(branches)
    for point_index, point in enumerate(unique_points, start=1):
        point_root = band.k_point_dir(root, point)
        state = "cached" if band.k_point_is_complete(point_root, config) else "run"
        print(
            f"[zeta={zeta:.3f} {stage}] {point_index}/{len(unique_points)} "
            f"kx={point['kx_str']} ky={point['ky_str']} ({state})"
        )
        band.run_k_point(root, hole_vertices, point, config)

    gamma_dir = band.k_point_dir(root, branches["gamma_m"][0])
    selected_modes = band.select_gamma_bands(
        pd.read_csv(gamma_dir / "eigenfrequencies.csv")
    )
    tracked_p2 = []
    for direction in ("gamma_m", "gamma_k"):
        tracked = band.track_selected_branch(root, branches[direction], selected_modes)
        tracked_p2.append(tracked[tracked["band_label"] == "p2"].copy())
    p2_rows = pd.concat(tracked_p2, ignore_index=True)
    p2_rows.insert(0, "zeta", round(float(zeta), 3))
    p2_rows.insert(1, "stage", str(stage))
    summary = summarize_p2_candidate(zeta, p2_rows, stage)
    band.write_csv_atomic(root / "selected_p2.csv", p2_rows)
    band.write_json_atomic(
        root / "config.json",
        {
            "zeta": round(float(zeta), 3),
            "stage": stage,
            "b0_nm": B0_NM,
            "eta": ETA,
            "params": params,
            "hole_params": hole_params,
            "triangles": holes,
            "minimum_clearance": minimum_clearance,
            "selected_gamma_modes": selected_modes,
            "branches": branches,
            "summary": summary,
        },
    )
    return p2_rows, summary


def plot_scan_summary(summary: pd.DataFrame, output_path: Path) -> None:
    """Plot the search metrics against zeta with all axes explicitly labeled."""
    ranked = summary.sort_values("zeta")
    feasible = ranked["feasible"].fillna(False).astype(bool)
    figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
    panels = (
        ("peak_k_signed", "Signed peak |k|/G", "Q-peak location"),
        ("peak_excess_dex", "Peak excess log10(Q/QΓ)", "Off-Γ Q excess"),
        ("gamma_frequency_thz", "p2 Γ frequency (THz)", "Γ frequency"),
        ("gamma_q", "p2 Q at Γ", "Γ Q"),
    )
    for axis, (column, ylabel, title) in zip(axes.flat, panels):
        axis.plot(
            ranked.loc[feasible, "zeta"],
            ranked.loc[feasible, column],
            color="#2563eb",
            marker="o",
            linewidth=1.2,
        )
        if (~feasible).any():
            axis.scatter(
                ranked.loc[~feasible, "zeta"],
                ranked.loc[~feasible, column],
                color="#9ca3af",
                marker="x",
                label="infeasible",
            )
        axis.set_xlabel("ζ")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
    axes[0, 0].axhline(0.0, color="0.4", linestyle="--", linewidth=0.8)
    axes[0, 1].axhline(0.0, color="0.4", linestyle="--", linewidth=0.8)
    axes[1, 0].axhspan(
        CLADDING_GAP_MIN_THZ,
        CLADDING_GAP_MAX_THZ,
        color="#bfdbfe",
        alpha=0.35,
        label="cladding gap",
    )
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].set_yscale("log")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def plot_best_q(p2_rows: pd.DataFrame, output_path: Path) -> None:
    """Plot the densely validated p2 Q-k curve on the signed branch axis."""
    figure, axis = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    styles = {"gamma_m": ("Γ–M", "-"), "gamma_k": ("Γ–K", "--")}
    for direction, (label, linestyle) in styles.items():
        rows = p2_rows[
            (p2_rows["direction"] == direction)
            & (p2_rows["match_status"] != "unmatched")
            & (pd.to_numeric(p2_rows["k_norm"], errors="coerce") <= 0.2)
        ].sort_values("k_norm")
        k = pd.to_numeric(rows["k_norm"], errors="coerce").to_numpy(dtype=float)
        q = pd.to_numeric(rows["q"], errors="coerce").to_numpy(dtype=float)
        usable = np.isfinite(k) & np.isfinite(q) & (q > 0.0)
        sign = -1.0 if direction == "gamma_m" else 1.0
        axis.plot(
            sign * k[usable],
            np.log10(q[usable]),
            color="#dc2626",
            linestyle=linestyle,
            marker="o",
            markersize=3,
            linewidth=1.2,
            label=label,
        )
    axis.axvline(0.0, color="0.6", linestyle="--", linewidth=0.8)
    axis.set_xlim(-0.2, 0.2)
    axis.set_xlabel("Signed |k|/G  (Γ–M < 0, Γ–K > 0)")
    axis.set_ylabel("log10(Q)")
    axis.set_title("Best cavity p2 Q–k")
    axis.grid(True, alpha=0.25)
    axis.legend()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def _persist_summary(records: dict[float, dict[str, object]], root: Path) -> pd.DataFrame:
    ranked = rank_candidates(pd.DataFrame(records.values()))
    band.write_csv_atomic(root / "scan_summary.csv", ranked)
    plot_scan_summary(ranked, root / "zeta_scan.png")
    return ranked


def main() -> dict[str, object]:
    """Run the restartable coarse, fine, and dense validation stages."""
    root = scan_root()
    root.mkdir(parents=True, exist_ok=True)
    band.write_json_atomic(
        root / "scan_config.json",
        {
            "b0_nm": B0_NM,
            "eta": ETA,
            "zeta_range": [ZETA_MIN, ZETA_MAX],
            "coarse_zeta_step": COARSE_ZETA_STEP,
            "fine_zeta_step": FINE_ZETA_STEP,
            "fine_half_width": FINE_HALF_WIDTH,
            "screen_k_magnitudes": SCREEN_K_MAGNITUDES,
            "final_k_magnitudes": FINAL_K_MAGNITUDES,
            "cladding_gap_thz": [CLADDING_GAP_MIN_THZ, CLADDING_GAP_MAX_THZ],
            "minimum_gamma_p_weight": MIN_GAMMA_P_WEIGHT,
        },
    )

    records: dict[float, dict[str, object]] = {}
    rows_by_zeta: dict[float, pd.DataFrame] = {}
    for zeta in coarse_zeta_values():
        key = round(float(zeta), 3)
        rows_by_zeta[key], records[key] = run_candidate(
            key, "coarse", SCREEN_K_MAGNITUDES
        )
        _persist_summary(records, root)

    coarse_ranked = rank_candidates(pd.DataFrame(records.values()))
    feasible_coarse = coarse_ranked[coarse_ranked["feasible"].astype(bool)]
    if feasible_coarse.empty:
        _persist_summary(records, root)
        raise RuntimeError(
            "No coarse zeta candidate satisfies the Gamma frequency and p-weight constraints"
        )
    center = float(feasible_coarse.iloc[0]["zeta"])

    fine_grid = fine_zeta_values(center)
    for zeta in fine_grid:
        key = round(float(zeta), 3)
        if key in records:
            continue
        rows_by_zeta[key], records[key] = run_candidate(
            key, "fine", SCREEN_K_MAGNITUDES
        )
        _persist_summary(records, root)

    fine_keys = {round(float(value), 3) for value in fine_grid}
    fine_ranked = rank_candidates(
        pd.DataFrame([records[key] for key in fine_keys if key in records])
    )
    final_candidates = fine_ranked[fine_ranked["feasible"].astype(bool)].head(
        FINAL_CANDIDATE_COUNT
    )
    if final_candidates.empty:
        _persist_summary(records, root)
        raise RuntimeError("No feasible fine zeta candidate is available for dense validation")

    final_records = []
    for zeta in final_candidates["zeta"].to_numpy(dtype=float):
        key = round(float(zeta), 3)
        rows_by_zeta[key], final_summary = run_candidate(
            key, "final", FINAL_K_MAGNITUDES
        )
        records[key] = final_summary
        final_records.append(final_summary)
        _persist_summary(records, root)
    validated = rank_candidates(pd.DataFrame(final_records))
    best_summary = validated.iloc[0].drop(labels="rank").to_dict()
    best_zeta = round(float(best_summary["zeta"]), 3)
    best_rows = rows_by_zeta[best_zeta]

    best_dir = root / "best"
    best_dir.mkdir(parents=True, exist_ok=True)
    band.write_csv_atomic(best_dir / "selected_p2.csv", best_rows)
    band.write_json_atomic(
        best_dir / "best_config.json",
        {
            "summary": best_summary,
            "b0_nm": B0_NM,
            "eta": ETA,
            "zeta": best_zeta,
            "params": make_cavity_params(best_zeta),
            "final_k_magnitudes": FINAL_K_MAGNITUDES,
        },
    )
    plot_best_q(best_rows, best_dir / "cavity_p2_q_vs_k.png")
    _persist_summary(records, root)
    return best_summary


if __name__ == "__main__":
    main()
