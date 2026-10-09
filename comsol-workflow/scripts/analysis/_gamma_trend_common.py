"""Shared helpers for Gamma-weight trend analysis scripts."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_frequency_thz(expression: str) -> float:
    match = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        expression,
    )
    if match is None:
        raise ValueError(f"Could not parse frequency from {expression!r}")
    return float(match.group(0))


def select_closest_mode_by_weight(
    mode_scores: pd.DataFrame,
    *,
    target_frequency: float,
    source: str,
    weight_column: str,
    weight_threshold: float,
    additional_numeric_columns: Iterable[str] = (),
) -> pd.Series:
    required = {"frequency", weight_column, *additional_numeric_columns}
    missing = required - set(mode_scores.columns)
    if missing:
        raise ValueError(f"{source} missing required columns: {sorted(missing)}")

    df = mode_scores.copy()
    for column in sorted(required):
        df[column] = pd.to_numeric(df[column], errors="coerce")
    numeric_values = df[list(required)].to_numpy(dtype=float)
    finite_rows = np.isfinite(numeric_values).all(axis=1)
    candidates = df.loc[
        finite_rows & (df[weight_column] > float(weight_threshold))
    ].copy()
    if candidates.empty:
        raise ValueError(
            f"{source} has no mode with valid frequency and "
            f"{weight_column} > {weight_threshold:g}"
        )

    candidates["frequency_distance_to_target"] = (
        candidates["frequency"] - float(target_frequency)
    ).abs()
    sort_columns = ["frequency_distance_to_target"]
    if "mode_idx" in candidates.columns:
        candidates["mode_idx"] = pd.to_numeric(
            candidates["mode_idx"], errors="coerce"
        )
        sort_columns.append("mode_idx")
    candidates = candidates.sort_values(sort_columns, ascending=True)
    return candidates.iloc[0]


def load_trend_csv(
    path: Path,
    *,
    required_plot_columns: Iterable[str],
    target_frequency: float | None = None,
    allowed_shift_factors: Iterable[float] | None = None,
    require_unique_target: bool = False,
) -> tuple[pd.DataFrame, float]:
    path = Path(path)
    selected_modes = pd.read_csv(path)
    required = set(required_plot_columns)
    missing = required - set(selected_modes.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")

    selected_modes = selected_modes.copy()
    for column in sorted(required):
        selected_modes[column] = pd.to_numeric(
            selected_modes[column], errors="coerce"
        )
    invalid_rows = selected_modes[list(required)].isna().any(axis=1)
    if invalid_rows.any():
        raise ValueError(f"{path} has non-numeric values in required plot columns")

    if target_frequency is None:
        if "target_frequency" not in selected_modes.columns:
            raise ValueError(
                f"{path} missing target_frequency column; pass --target-frequency"
            )
        target_values = pd.to_numeric(
            selected_modes["target_frequency"], errors="coerce"
        ).dropna()
        if target_values.empty:
            raise ValueError(f"{path} has no numeric target_frequency values")
        unique_targets = target_values.unique()
        if require_unique_target and len(unique_targets) != 1:
            raise ValueError(f"{path} has inconsistent target_frequency values")
        target_frequency = float(unique_targets[0])
    else:
        target_frequency = float(target_frequency)

    if allowed_shift_factors is not None:
        selected_modes = selected_modes.loc[
            selected_modes["shift_factor"].isin(list(allowed_shift_factors))
        ].copy()
        if selected_modes.empty:
            raise ValueError(f"{path} has no rows within configured shift factors")

    selected_modes = selected_modes.sort_values("shift_factor")
    return selected_modes, target_frequency


def save_weight_frequency_plot(
    selected_modes: pd.DataFrame,
    path: Path,
    *,
    target_frequency: float,
    weight_column: str,
    title: str,
    target_label: str,
) -> None:
    shift_factors = sorted(selected_modes["shift_factor"].unique())
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    weight_line, = ax.plot(
        selected_modes["shift_factor"],
        selected_modes[weight_column],
        marker="o",
        linewidth=1.6,
        markersize=5.5,
        color="#2563eb",
        label=weight_column,
    )
    ax.set_xlabel("CLADDING_INWARD_SHIFT / A")
    ax.set_ylabel(weight_column)
    ax.tick_params(axis="y", labelcolor="#2563eb")
    ax.set_title(title)
    ax.set_xticks(shift_factors)
    ax.set_xlim(min(shift_factors) - 0.005, max(shift_factors) + 0.005)
    ax.grid(True, alpha=0.28, linestyle="--", linewidth=0.7)

    freq_ax = ax.twinx()
    freq_line, = freq_ax.plot(
        selected_modes["shift_factor"],
        selected_modes["frequency"],
        marker="s",
        linewidth=1.5,
        markersize=5.0,
        color="#dc2626",
        label="frequency",
    )
    target_line = freq_ax.axhline(
        target_frequency,
        color="#dc2626",
        linestyle="--",
        linewidth=1.1,
        alpha=0.8,
        label=target_label,
    )
    freq_ax.set_ylabel("frequency (THz)")
    freq_ax.tick_params(axis="y", labelcolor="#dc2626")

    lines = [weight_line, freq_line, target_line]
    ax.legend(
        lines,
        [line.get_label() for line in lines],
        loc="upper left",
        bbox_to_anchor=(0.02, 0.82),
        borderaxespad=0.0,
        fontsize=8,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)
