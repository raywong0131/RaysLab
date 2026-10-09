#!/usr/bin/env python3
"""Analyze cavity Hz flatness and cladding-shell trends for selected finite modes."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.hex_lattice_utils import unit_cell_corners  # noqa: E402
from comsol_workflow.lattice_fourier_postprocess import (  # noqa: E402
    CENTER_NORMALIZED_INTENSITY_CMAP,
    CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL,
    CENTER_NORMALIZED_INTENSITY_TICKS,
    CENTER_NORMALIZED_INTENSITY_TICK_LABELS,
    CENTER_NORMALIZED_INTENSITY_VCENTER,
    CENTER_NORMALIZED_INTENSITY_VMAX,
    CENTER_NORMALIZED_INTENSITY_VMIN,
    INTENSITY_MAPS_NORMH_FILENAME,
)


SUMMARY_FILENAME = "finite_flatness_trend.csv"
SHELL_FILENAME = "finite_flatness_shells.csv"
CELL_FILENAME = "finite_flatness_cells.csv"
CANONICAL_FULL_CELL_MAP_FILENAME = "cavity_cladding_intensity_maps.png"
COMPATIBILITY_FULL_CELL_MAP_FILENAME = INTENSITY_MAPS_NORMH_FILENAME
LINEAR_CENTER_MAP_CMAP = CENTER_NORMALIZED_INTENSITY_CMAP
LINEAR_CENTER_MAP_VMIN = CENTER_NORMALIZED_INTENSITY_VMIN
LINEAR_CENTER_MAP_VCENTER = CENTER_NORMALIZED_INTENSITY_VCENTER
LINEAR_CENTER_MAP_VMAX = CENTER_NORMALIZED_INTENSITY_VMAX
LINEAR_CENTER_MAP_TICKS = CENTER_NORMALIZED_INTENSITY_TICKS
LINEAR_CENTER_MAP_TICK_LABELS = CENTER_NORMALIZED_INTENSITY_TICK_LABELS
LINEAR_CENTER_MAP_TITLE = (
    "Cavity And Cladding Unit-Cell Mean Hz Intensity\n"
    r"Linear scale; colorbar midpoint marks $I_R/I_0=1$"
)
LINEAR_CENTER_MAP_COLORBAR_LABEL = CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL
REQUIRED_TREND_COLUMNS = {
    "shift_factor",
    "folder",
    "mode_idx",
    "gamma_k_weight_fraction",
}


def lattice_shell(indices: np.ndarray) -> np.ndarray:
    indices = np.asarray(indices, dtype=int)
    return np.maximum.reduce(
        [
            np.abs(indices[:, 0]),
            np.abs(indices[:, 1]),
            np.abs(indices[:, 0] + indices[:, 1]),
        ]
    )


def normalized_intensity_metrics(intensities: np.ndarray) -> dict[str, float]:
    intensities = np.asarray(intensities, dtype=float)
    if intensities.ndim != 1 or intensities.size == 0:
        raise ValueError("intensities must be a non-empty one-dimensional array")
    if not np.isfinite(intensities).all() or np.any(intensities < 0.0):
        raise ValueError("intensities must be finite and non-negative")
    mean = float(np.mean(intensities))
    if mean <= 0.0:
        raise ValueError("mean intensity must be positive")
    normalized = intensities / mean
    cv = float(np.std(normalized))
    return {
        "mean": mean,
        "cv": cv,
        "flatness": float(1.0 / (1.0 + cv * cv)),
        "dmax": float(np.max(np.abs(normalized - 1.0))),
    }


def recover_cavity_cell_data(npz_path: Path) -> dict[str, np.ndarray | float]:
    with np.load(npz_path) as data:
        f_components = np.asarray(data["F"], dtype=complex)
        fields_by_cyclic_index = np.fft.fft(f_components, axis=0, norm="ortho")
        intensity_by_cyclic_index = np.mean(
            np.abs(fields_by_cyclic_index) ** 2,
            axis=1,
        )
        cyclic_indices = np.asarray(data["cyclic_indices"], dtype=int)
        cell_intensity = intensity_by_cyclic_index[cyclic_indices]
        weight_fraction = np.asarray(data["weight_fraction"], dtype=float)
        return {
            "mode_idx": int(data["mode_idx"]),
            "cell_indices": np.asarray(data["cell_indices"], dtype=int),
            "cell_centers": np.asarray(data["cell_centers"], dtype=float),
            "rho_points": np.asarray(data["rho_points"], dtype=float),
            "cell_intensity": cell_intensity,
            "gamma_k_weight_fraction": float(weight_fraction[0]),
            "parseval_error": float(data["parseval_error"]),
            "reconstruction_error": float(data["reconstruction_error"]),
        }


def deduplicate_complex_field(field_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    required = {"x", "y", "re", "im"}
    missing = required - set(field_df.columns)
    if missing:
        raise ValueError(f"field table missing columns: {sorted(missing)}")
    points = np.round(field_df[["x", "y"]].to_numpy(dtype=float), decimals=12)
    values = field_df["re"].to_numpy(dtype=float) + 1j * field_df["im"].to_numpy(
        dtype=float
    )
    unique_points, inverse, counts = np.unique(
        points,
        axis=0,
        return_inverse=True,
        return_counts=True,
    )
    summed_real = np.bincount(inverse, weights=np.real(values))
    summed_imag = np.bincount(inverse, weights=np.imag(values))
    unique_values = (summed_real + 1j * summed_imag) / counts
    return unique_points, unique_values


class ComplexFieldInterpolator:
    def __init__(self, points: np.ndarray, values: np.ndarray):
        triangulation = Delaunay(np.asarray(points, dtype=float))
        self.linear = LinearNDInterpolator(
            triangulation,
            np.asarray(values, dtype=complex),
            fill_value=np.nan + 1j * np.nan,
        )
        self.nearest = NearestNDInterpolator(points, values)

    def __call__(self, query_points: np.ndarray) -> np.ndarray:
        values = np.asarray(self.linear(query_points), dtype=complex)
        missing = ~np.isfinite(np.real(values)) | ~np.isfinite(np.imag(values))
        if np.any(missing):
            values[missing] = self.nearest(np.asarray(query_points)[missing])
        return values


def interpolate_cell_intensities(
    interpolator: ComplexFieldInterpolator,
    centers: np.ndarray,
    rho_points: np.ndarray,
    *,
    chunk_size: int = 64,
) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    rho_points = np.asarray(rho_points, dtype=float)
    intensity = np.empty(len(centers), dtype=float)
    for start in range(0, len(centers), chunk_size):
        stop = min(start + chunk_size, len(centers))
        query = centers[start:stop, None, :] + rho_points[None, :, :]
        values = interpolator(query.reshape(-1, 2)).reshape(
            stop - start,
            len(rho_points),
        )
        intensity[start:stop] = np.mean(np.abs(values) ** 2, axis=1)
    return intensity


def shell_summary(
    *,
    shift_factor: float,
    mode_idx: int,
    shells: np.ndarray,
    intensities: np.ndarray,
    cavity_mean: float,
    bulk_radius: int,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for shell in sorted(np.unique(shells)):
        values = np.asarray(intensities[shells == shell], dtype=float)
        mean = float(np.mean(values))
        relative = values / cavity_mean
        rows.append(
            {
                "shift_factor": shift_factor,
                "mode_idx": mode_idx,
                "region": "cavity" if shell <= bulk_radius else "cladding",
                "shell": int(shell),
                "relative_shell": int(shell - bulk_radius),
                "n_cells": int(len(values)),
                "mean_intensity": mean,
                "mean_over_cavity": float(mean / cavity_mean),
                "min_over_cavity": float(np.min(relative)),
                "p95_over_cavity": float(np.quantile(relative, 0.95)),
                "max_over_cavity": float(np.max(relative)),
                "shell_cv": float(np.std(values) / mean) if mean > 0.0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def cladding_guardrails(
    *,
    cladding_intensity: np.ndarray,
    cladding_shells: np.ndarray,
    cavity_mean: float,
    cavity_outer_mean_ratio: float,
    bulk_radius: int,
) -> dict[str, float]:
    shell_df = shell_summary(
        shift_factor=0.0,
        mode_idx=0,
        shells=cladding_shells,
        intensities=cladding_intensity,
        cavity_mean=cavity_mean,
        bulk_radius=bulk_radius,
    ).sort_values("shell")
    profile = shell_df["mean_over_cavity"].to_numpy(dtype=float)
    relative = np.asarray(cladding_intensity, dtype=float) / cavity_mean
    upward = float(np.maximum(np.diff(profile), 0.0).sum())
    transition_profile = np.concatenate([[cavity_outer_mean_ratio], profile])
    return {
        "cladding_hotspot_max": float(np.max(relative)),
        "cladding_hotspot_p95": float(np.quantile(relative, 0.95)),
        "first_cladding_mean_ratio": float(profile[0]),
        "first_cladding_over_cavity_outer": float(
            profile[0] / cavity_outer_mean_ratio
        ),
        "cladding_upward_variation": upward,
        "boundary_plus_cladding_upward_variation": float(
            np.maximum(np.diff(transition_profile), 0.0).sum()
        ),
        "cladding_max_shell_cv": float(shell_df["shell_cv"].max()),
    }


def _case_paths(series_root: Path, trend_row) -> dict[str, Path]:
    case_dir = series_root / str(trend_row.folder)
    mode_dir = case_dir / "01_results" / f"mode{int(trend_row.mode_idx)}"
    return {
        "case": case_dir,
        "config": case_dir / "99_config" / "config.json",
        "npz": mode_dir / "13_lattice_fourier_Hz" / "finite_lattice_fourier_hz.npz",
        "field": mode_dir / "11_simulation_exports" / "Hz_center.parquet",
    }


def analyze_selected_mode(series_root: Path, trend_row) -> tuple[
    dict[str, object], pd.DataFrame, pd.DataFrame
]:
    paths = _case_paths(series_root, trend_row)
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing selected-mode inputs: {missing}")

    config = json.loads(paths["config"].read_text(encoding="utf-8"))
    period = float(config["a"])
    bulk_radius = int(config["bulk_radius"])
    cladding_layers = int(config["cladding_layers"])
    cavity = recover_cavity_cell_data(paths["npz"])
    if int(cavity["mode_idx"]) != int(trend_row.mode_idx):
        raise ValueError(f"{paths['npz']} mode_idx disagrees with selected trend mode")
    gamma_npz = float(cavity["gamma_k_weight_fraction"])
    if not math.isclose(
        gamma_npz,
        float(trend_row.gamma_k_weight_fraction),
        rel_tol=0.0,
        abs_tol=2e-12,
    ):
        raise ValueError(
            f"{paths['npz']} Gamma weight {gamma_npz} disagrees with trend "
            f"{trend_row.gamma_k_weight_fraction}"
        )

    cavity_indices = np.asarray(cavity["cell_indices"], dtype=int)
    cavity_centers = np.asarray(cavity["cell_centers"], dtype=float)
    cavity_shells = lattice_shell(cavity_indices)
    cavity_intensity = np.asarray(cavity["cell_intensity"], dtype=float)
    cavity_metrics = normalized_intensity_metrics(cavity_intensity)
    cavity_mean = cavity_metrics["mean"]
    cavity_outer_mean_ratio = float(
        np.mean(cavity_intensity[cavity_shells == bulk_radius]) / cavity_mean
    )

    field_df = pd.read_parquet(paths["field"], columns=["x", "y", "re", "im"])
    field_points, field_values = deduplicate_complex_field(field_df)
    interpolator = ComplexFieldInterpolator(field_points, field_values)

    interpolated_cavity_intensity = interpolate_cell_intensities(
        interpolator,
        cavity_centers,
        np.asarray(cavity["rho_points"], dtype=float),
    )
    cavity_interp_normalized = interpolated_cavity_intensity / np.mean(
        interpolated_cavity_intensity
    )
    cavity_npz_normalized = cavity_intensity / cavity_mean
    interpolation_rmse = float(
        np.sqrt(np.mean((cavity_interp_normalized - cavity_npz_normalized) ** 2))
    )

    cladding_points = config["cladding_points"]
    cladding_centers = np.array(
        [[point["x"], point["y"]] for point in cladding_points],
        dtype=float,
    )
    cladding_indices = np.array(
        [[point["i"], point["j"]] for point in cladding_points],
        dtype=int,
    )
    cladding_shells = np.array(
        [point.get("shell", lattice_shell(cladding_indices)[idx]) for idx, point in enumerate(cladding_points)],
        dtype=int,
    )
    expected_shells = np.arange(
        bulk_radius + 1,
        bulk_radius + cladding_layers + 1,
        dtype=int,
    )
    if not np.array_equal(np.unique(cladding_shells), expected_shells):
        raise ValueError(
            f"{paths['config']} has cladding shells {np.unique(cladding_shells)}, "
            f"expected {expected_shells}"
        )
    cladding_intensity = interpolate_cell_intensities(
        interpolator,
        cladding_centers,
        np.asarray(cavity["rho_points"], dtype=float),
    )
    guardrails = cladding_guardrails(
        cladding_intensity=cladding_intensity,
        cladding_shells=cladding_shells,
        cavity_mean=cavity_mean,
        cavity_outer_mean_ratio=cavity_outer_mean_ratio,
        bulk_radius=bulk_radius,
    )

    shift_factor = float(trend_row.shift_factor)
    mode_idx = int(trend_row.mode_idx)
    summary = trend_row._asdict()
    summary.update(
        {
            "period": period,
            "bulk_radius": bulk_radius,
            "cladding_layers": cladding_layers,
            "n_cavity_cells": int(len(cavity_intensity)),
            "n_cladding_cells": int(len(cladding_intensity)),
            "n_rho_points": int(len(cavity["rho_points"])),
            "gamma_k_weight_fraction_npz": gamma_npz,
            "cavity_intensity_cv": cavity_metrics["cv"],
            "cavity_intensity_flatness": cavity_metrics["flatness"],
            "cavity_intensity_dmax": cavity_metrics["dmax"],
            "cavity_outer_mean_ratio": cavity_outer_mean_ratio,
            "cavity_interpolation_normalized_rmse": interpolation_rmse,
            "parseval_error": float(cavity["parseval_error"]),
            "reconstruction_error": float(cavity["reconstruction_error"]),
            **guardrails,
        }
    )

    shell_df = pd.concat(
        [
            shell_summary(
                shift_factor=shift_factor,
                mode_idx=mode_idx,
                shells=cavity_shells,
                intensities=cavity_intensity,
                cavity_mean=cavity_mean,
                bulk_radius=bulk_radius,
            ),
            shell_summary(
                shift_factor=shift_factor,
                mode_idx=mode_idx,
                shells=cladding_shells,
                intensities=cladding_intensity,
                cavity_mean=cavity_mean,
                bulk_radius=bulk_radius,
            ),
        ],
        ignore_index=True,
    )

    cavity_cells = pd.DataFrame(
        {
            "shift_factor": shift_factor,
            "mode_idx": mode_idx,
            "region": "cavity",
            "shell": cavity_shells,
            "relative_shell": cavity_shells - bulk_radius,
            "i": cavity_indices[:, 0],
            "j": cavity_indices[:, 1],
            "x": cavity_centers[:, 0],
            "y": cavity_centers[:, 1],
            "intensity": cavity_intensity,
            "intensity_over_cavity": cavity_intensity / cavity_mean,
        }
    )
    cladding_cells = pd.DataFrame(
        {
            "shift_factor": shift_factor,
            "mode_idx": mode_idx,
            "region": "cladding",
            "shell": cladding_shells,
            "relative_shell": cladding_shells - bulk_radius,
            "i": cladding_indices[:, 0],
            "j": cladding_indices[:, 1],
            "x": cladding_centers[:, 0],
            "y": cladding_centers[:, 1],
            "intensity": cladding_intensity,
            "intensity_over_cavity": cladding_intensity / cavity_mean,
        }
    )
    return summary, shell_df, pd.concat([cavity_cells, cladding_cells], ignore_index=True)


def load_selected_trend(path: Path) -> pd.DataFrame:
    selected = pd.read_csv(path)
    missing = REQUIRED_TREND_COLUMNS - set(selected.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    selected = selected.copy()
    for column in ("shift_factor", "mode_idx", "gamma_k_weight_fraction"):
        selected[column] = pd.to_numeric(selected[column], errors="raise")
    if selected["shift_factor"].duplicated().any():
        raise ValueError(f"{path} contains duplicate shift_factor rows")
    return selected.sort_values("shift_factor")


def _shift_axis_label(summary: pd.DataFrame) -> str:
    if "shift_kind" in summary.columns:
        kinds = summary["shift_kind"].dropna().astype(str)
        if not kinds.empty and kinds.eq("directional").all():
            return "CLADDING X SHIFT / A"
    return "CLADDING_INWARD_SHIFT / A"


def _set_shift_axis(ax, shifts: np.ndarray, *, label: str) -> None:
    ax.set_xticks(shifts)
    ax.set_xlim(float(np.min(shifts) - 0.005), float(np.max(shifts) + 0.005))
    ax.set_xlabel(label)
    ax.grid(True, alpha=0.22, linestyle="--", linewidth=0.7)


def save_cavity_trend(summary: pd.DataFrame, path: Path) -> None:
    shifts = summary["shift_factor"].to_numpy(dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(
        shifts,
        summary["gamma_k_weight_fraction"],
        color="#2563eb",
        marker="o",
        linewidth=1.8,
        label=r"$P_\Gamma$",
    )
    axes[0].set_ylabel(r"Gamma weight $P_\Gamma$")
    axes[0].set_title("Finite-Cavity Hz Flatness Versus Cladding Shift")
    axes[0].legend(loc="best")
    axes[0].grid(True, alpha=0.22, linestyle="--", linewidth=0.7)

    axes[1].plot(
        shifts,
        100.0 * summary["cavity_intensity_cv"],
        color="#d97706",
        marker="s",
        linewidth=1.7,
        label="cavity intensity CV",
    )
    axes[1].plot(
        shifts,
        100.0 * summary["cavity_intensity_dmax"],
        color="#6b7280",
        marker="^",
        linewidth=1.5,
        linestyle="--",
        label="maximum cell deviation",
    )
    axes[1].set_ylabel("Intensity variation (%)")
    axes[1].legend(loc="best")
    _set_shift_axis(axes[1], shifts, label=_shift_axis_label(summary))
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_cladding_trend(summary: pd.DataFrame, path: Path) -> None:
    shifts = summary["shift_factor"].to_numpy(dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(
        shifts,
        summary["cladding_hotspot_max"],
        color="#c2410c",
        marker="o",
        linewidth=1.8,
        label="maximum cladding cell",
    )
    axes[0].plot(
        shifts,
        summary["cladding_hotspot_p95"],
        color="#f59e0b",
        marker="s",
        linewidth=1.6,
        label="cladding cell P95",
    )
    axes[0].axhline(1.0, color="#374151", linewidth=1.0, linestyle=":")
    axes[0].set_ylabel("Intensity / cavity mean")
    axes[0].set_title("Cladding Hz Guardrails Versus Cladding Shift")
    axes[0].legend(loc="best")
    axes[0].grid(True, alpha=0.22, linestyle="--", linewidth=0.7)

    axes[1].plot(
        shifts,
        summary["boundary_plus_cladding_upward_variation"],
        color="#2563eb",
        marker="o",
        linewidth=1.7,
        label="boundary + outward upward variation",
    )
    axes[1].plot(
        shifts,
        summary["cladding_max_shell_cv"],
        color="#6b7280",
        marker="^",
        linewidth=1.5,
        linestyle="--",
        label="maximum within-shell CV",
    )
    axes[1].set_ylabel("Dimensionless variation")
    axes[1].legend(loc="best")
    _set_shift_axis(axes[1], shifts, label=_shift_axis_label(summary))
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_shell_profiles(
    shells: pd.DataFrame,
    path: Path,
    *,
    shift_axis_label: str = "CLADDING_INWARD_SHIFT / A",
) -> None:
    cladding = shells.loc[shells["region"].eq("cladding")].copy()
    shifts = np.sort(cladding["shift_factor"].unique())
    cmap = plt.get_cmap("viridis")
    norm = Normalize(float(shifts.min()), float(shifts.max()))
    fig, ax = plt.subplots(figsize=(8.2, 5.5), constrained_layout=True)
    for shift in shifts:
        group = cladding.loc[cladding["shift_factor"].eq(shift)].sort_values(
            "relative_shell"
        )
        ax.plot(
            group["relative_shell"],
            group["mean_over_cavity"],
            color=cmap(norm(float(shift))),
            marker="o",
            markersize=3.8,
            linewidth=1.4,
        )
    ax.set_yscale("log")
    ax.set_xlabel("Cladding shell (1 = adjacent to cavity)")
    ax.set_ylabel("Shell mean intensity / cavity mean")
    ax.set_title("Cladding-Shell Mean Hz Intensity Profiles")
    ax.set_xticks(sorted(cladding["relative_shell"].unique()))
    ax.grid(True, alpha=0.22, linestyle="--", linewidth=0.7)
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        label=shift_axis_label,
    )
    colorbar.set_ticks(shifts)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _cell_patches(group: pd.DataFrame, period: float) -> list[Polygon]:
    base = unit_cell_corners(period)
    return [
        Polygon(base + np.array([row.x, row.y]), closed=True)
        for row in group.itertuples(index=False)
    ]


def _make_shift_panel_axes(count: int):
    """Create a compact panel grid with exactly ``count`` active axes."""
    if count <= 0:
        raise ValueError("at least one shift is required for cell-map plotting")
    max_columns = 5
    nrows = math.ceil(count / max_columns)
    ncols = math.ceil(count / nrows)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.5 * ncols, 3.6 * nrows),
        constrained_layout=True,
        squeeze=False,
    )
    all_axes = list(axes.ravel())
    active_axes = all_axes[:count]
    for ax in all_axes[count:]:
        ax.set_visible(False)
    return fig, active_axes


def save_cell_maps(cells: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    shifts = summary["shift_factor"].to_numpy(dtype=float)
    cavity_values = cells.loc[
        cells["region"].eq("cavity"), "intensity_over_cavity"
    ].to_numpy(dtype=float)
    log2_values = np.log2(np.maximum(cavity_values, np.finfo(float).tiny))
    bound = max(1.0, float(np.quantile(np.abs(log2_values), 0.99)))
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
    cmap = plt.get_cmap("coolwarm")
    fig, active_axes = _make_shift_panel_axes(len(shifts))
    for ax, shift in zip(active_axes, shifts, strict=True):
        group = cells.loc[
            cells["shift_factor"].eq(shift) & cells["region"].eq("cavity")
        ]
        period = float(summary.loc[summary["shift_factor"].eq(shift), "period"].iloc[0])
        collection = PatchCollection(
            _cell_patches(group, period),
            cmap=cmap,
            norm=norm,
            edgecolor="none",
        )
        collection.set_array(
            np.log2(group["intensity_over_cavity"].to_numpy(dtype=float))
        )
        ax.add_collection(collection)
        ax.autoscale_view()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"shift={shift:.3f}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Cavity Unit-Cell Mean Hz Intensity", fontsize=14)
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=active_axes,
        label=r"$\log_2(I_R / \overline{I}_{cavity})$",
        shrink=0.86,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_full_cell_maps(cells: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    shifts = summary["shift_factor"].to_numpy(dtype=float)
    log2_values = np.log2(
        np.maximum(
            cells["intensity_over_cavity"].to_numpy(dtype=float),
            np.finfo(float).tiny,
        )
    )
    vmin = min(-1.0, float(np.quantile(log2_values, 0.01)))
    vmax = max(1.0, float(np.quantile(log2_values, 0.99)))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")
    fig, active_axes = _make_shift_panel_axes(len(shifts))
    for ax, shift in zip(active_axes, shifts, strict=True):
        group = cells.loc[cells["shift_factor"].eq(shift)]
        period = float(
            summary.loc[summary["shift_factor"].eq(shift), "period"].iloc[0]
        )
        collection = PatchCollection(
            _cell_patches(group, period),
            cmap=cmap,
            norm=norm,
            edgecolor="none",
        )
        collection.set_array(
            np.log2(group["intensity_over_cavity"].to_numpy(dtype=float))
        )
        ax.add_collection(collection)
        ax.autoscale_view()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"shift={shift:.3f}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Cavity And Cladding Unit-Cell Mean Hz Intensity", fontsize=14)
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=active_axes,
        label=r"$\log_2(I_R / \overline{I}_{cavity})$",
        shrink=0.86,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)


def normalize_cells_by_cavity_center(cells: pd.DataFrame) -> pd.DataFrame:
    required = {"shift_factor", "region", "i", "j", "intensity"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"cell table missing required columns: {sorted(missing)}")
    center_rows = cells.loc[
        cells["region"].eq("cavity") & cells["i"].eq(0) & cells["j"].eq(0),
        ["shift_factor", "intensity"],
    ].rename(columns={"intensity": "center_cell_intensity"})
    if center_rows["shift_factor"].duplicated().any():
        raise ValueError("cell table has multiple cavity center cells for one shift")
    expected_shifts = set(cells["shift_factor"].unique())
    if set(center_rows["shift_factor"]) != expected_shifts:
        raise ValueError("cell table is missing a cavity center cell for one or more shifts")
    if (center_rows["center_cell_intensity"] <= 0.0).any():
        raise ValueError("cavity center-cell intensity must be positive")
    normalized = cells.merge(center_rows, on="shift_factor", validate="many_to_one")
    normalized["intensity_over_center"] = (
        normalized["intensity"] / normalized["center_cell_intensity"]
    )
    return normalized


def save_full_cell_maps_linear_center_normalized(
    cells: pd.DataFrame,
    summary: pd.DataFrame,
    path: Path,
    *,
    cmap_name: str = LINEAR_CENTER_MAP_CMAP,
) -> None:
    cells = normalize_cells_by_cavity_center(cells)
    shifts = summary["shift_factor"].to_numpy(dtype=float)
    # Keep the requested 0--2 linear scale, but make the physical reference
    # I_R / I_0 = 1 the visual midpoint.  The colormap is configurable so
    # palettes can be compared without changing the quantitative mapping.
    norm = TwoSlopeNorm(
        vmin=LINEAR_CENTER_MAP_VMIN,
        vcenter=LINEAR_CENTER_MAP_VCENTER,
        vmax=LINEAR_CENTER_MAP_VMAX,
    )
    cmap = plt.get_cmap(cmap_name)
    fig, active_axes = _make_shift_panel_axes(len(shifts))
    for ax, shift in zip(active_axes, shifts, strict=True):
        group = cells.loc[cells["shift_factor"].eq(shift)]
        period = float(
            summary.loc[summary["shift_factor"].eq(shift), "period"].iloc[0]
        )
        collection = PatchCollection(
            _cell_patches(group, period),
            cmap=cmap,
            norm=norm,
            edgecolor="none",
        )
        collection.set_array(
            np.clip(
                group["intensity_over_center"].to_numpy(dtype=float),
                LINEAR_CENTER_MAP_VMIN,
                LINEAR_CENTER_MAP_VMAX,
            )
        )
        ax.add_collection(collection)

        ax.autoscale_view()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"shift={shift:.3f}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(LINEAR_CENTER_MAP_TITLE, fontsize=14)
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=active_axes,
        label=LINEAR_CENTER_MAP_COLORBAR_LABEL,
        shrink=0.86,
    )
    colorbar.set_ticks(LINEAR_CENTER_MAP_TICKS)
    colorbar.set_ticklabels(LINEAR_CENTER_MAP_TICK_LABELS)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_analysis(
    selected_modes_csv: Path,
    *,
    series_root: Path | None = None,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    selected_modes_csv = Path(selected_modes_csv)
    selected = load_selected_trend(selected_modes_csv)
    series_root = (
        Path(series_root)
        if series_root is not None
        else selected_modes_csv.parent.parent
    )
    output_dir = (
        Path(output_dir)
        if output_dir is not None
        else selected_modes_csv.parent / "flatness_trend"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries: list[dict[str, object]] = []
    shell_frames: list[pd.DataFrame] = []
    cell_frames: list[pd.DataFrame] = []
    for trend_row in selected.itertuples(index=False):
        print(
            f"Analyzing shift={float(trend_row.shift_factor):.3f}, "
            f"mode={int(trend_row.mode_idx)}",
            flush=True,
        )
        summary, shell_df, cell_df = analyze_selected_mode(series_root, trend_row)
        summaries.append(summary)
        shell_frames.append(shell_df)
        cell_frames.append(cell_df)

    summary_df = pd.DataFrame(summaries).sort_values("shift_factor")
    shell_df = pd.concat(shell_frames, ignore_index=True).sort_values(
        ["shift_factor", "shell"]
    )
    cell_df = pd.concat(cell_frames, ignore_index=True).sort_values(
        ["shift_factor", "region", "shell", "i", "j"]
    )

    summary_path = output_dir / SUMMARY_FILENAME
    shell_path = output_dir / SHELL_FILENAME
    cell_path = output_dir / CELL_FILENAME
    summary_df.to_csv(summary_path, index=False)
    shell_df.to_csv(shell_path, index=False)
    cell_df.to_csv(cell_path, index=False)

    cavity_trend_path = output_dir / "cavity_flatness_trend.png"
    cladding_trend_path = output_dir / "cladding_guardrail_trend.png"
    shell_profile_path = output_dir / "cladding_shell_profiles.png"
    cell_map_path = output_dir / "cavity_intensity_maps.png"
    full_cell_map_path = output_dir / CANONICAL_FULL_CELL_MAP_FILENAME
    intensity_maps_normh_path = output_dir / COMPATIBILITY_FULL_CELL_MAP_FILENAME
    save_cavity_trend(summary_df, cavity_trend_path)
    save_cladding_trend(summary_df, cladding_trend_path)
    save_shell_profiles(
        shell_df,
        shell_profile_path,
        shift_axis_label=_shift_axis_label(summary_df),
    )
    save_cell_maps(cell_df, summary_df, cell_map_path)
    save_full_cell_maps_linear_center_normalized(
        cell_df,
        summary_df,
        full_cell_map_path,
    )
    shutil.copyfile(full_cell_map_path, intensity_maps_normh_path)
    return {
        "summary": summary_path,
        "shells": shell_path,
        "cells": cell_path,
        "cavity_trend": cavity_trend_path,
        "cladding_trend": cladding_trend_path,
        "shell_profiles": shell_profile_path,
        "cell_maps": cell_map_path,
        "full_cell_maps": full_cell_map_path,
        "linear_center_cell_maps": intensity_maps_normh_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-modes-csv", type=Path, required=True)
    parser.add_argument("--series-root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        outputs = run_analysis(
            args.selected_modes_csv,
            series_root=args.series_root,
            output_dir=args.output_dir,
        )
    except Exception as exc:
        print(f"Could not generate finite flatness trend: {exc}", file=sys.stderr)
        return 1
    for label, path in outputs.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
