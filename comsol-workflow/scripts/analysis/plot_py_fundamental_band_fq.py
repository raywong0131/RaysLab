#!/usr/bin/env python3
"""Plot cavity/cladding bands beside finite-cavity Q and mark the user p_y mode."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.transforms import offset_copy

from scripts.run_main import run_band_pair as band


FREQUENCY_LIMITS = (190.0, 210.0)
FREQUENCY_TICKS = (190.0, 195.0, 200.0, 205.0, 210.0)
K_LIMIT = 0.3
OUTPUT_NAME = "p_y_fundamental_band_f_Q.png"
HIGHLIGHT_COLOR = "#ef2b16"
MODE_COLOR = "#2f93c5"
BAND_BOX_ASPECT = 1.62
QUALITY_BOX_ASPECT = 1.73
FONT_SIZE = 14.0
TITLE_FONT_SIZE = 16.0
LEGEND_FONT_SIZE = 12.0


def _read_csv(path: Path, required: set[str]) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    return frame


def _selected_py_mode(finite_run: Path, merged: pd.DataFrame) -> pd.Series:
    paths = list(
        finite_run.glob(
            "shift*/symmetry_1_xPEC_yPMC/10_overview/analysis_mode_selection.csv"
        )
    )
    if len(paths) != 1:
        raise ValueError(f"Expected one p_y selection table, found {len(paths)}")
    selection = _read_csv(
        paths[0],
        {
            "mode_idx",
            "frequency",
            "q",
            "target_internal_mode",
            "fundamental_eligible",
            "analysis_selected",
        },
    )
    chosen = selection[
        selection["analysis_selected"].astype(str).str.lower().eq("true")
    ]
    if len(chosen) != 1:
        raise ValueError(f"Expected one selected p_y fundamental, found {len(chosen)}")
    mode = chosen.iloc[0]
    if str(mode["target_internal_mode"]) != "px":
        raise ValueError("User p_y must map to the project-internal px component")
    if str(mode["fundamental_eligible"]).lower() != "true":
        raise ValueError("Selected p_y mode is not fundamental_eligible")

    matches = merged[
        (pd.to_numeric(merged["symmetry_id"], errors="coerce") == 1)
        & (pd.to_numeric(merged["mode_idx"], errors="coerce") == int(mode["mode_idx"]))
        & np.isclose(
            pd.to_numeric(merged["re"], errors="coerce"),
            float(mode["frequency"]),
        )
    ]
    if len(matches) != 1 or not np.isclose(float(matches.iloc[0]["q"]), float(mode["q"])):
        raise ValueError("Selected p_y frequency/Q does not match the merged table")
    return mode


def _band_character(rows: pd.DataFrame, label: str) -> tuple[np.ndarray, object]:
    if label.startswith("p"):
        values = band.p_x_fraction(
            pd.to_numeric(rows["px_weight"], errors="coerce"),
            pd.to_numeric(rows["py_weight"], errors="coerce"),
        )
        color_map = LinearSegmentedColormap.from_list(
            "p_character",
            [band.P_Y_COLOR, "white", band.P_X_COLOR],
        )
        return values, color_map
    first = pd.to_numeric(rows["dy_weight"], errors="coerce").to_numpy(float)
    second = pd.to_numeric(rows["dx_weight"], errors="coerce").to_numpy(float)
    total = first + second
    values = np.full(total.shape, 0.5, dtype=float)
    np.divide(second, total, out=values, where=total > 0.0)
    return np.clip(values, 0.0, 1.0), LinearSegmentedColormap.from_list(
        "d_character",
        ["#b0b0b0", "black"],
    )


def _plot_bands(axis, frame: pd.DataFrame, title: str) -> None:
    for label in ("p1", "p2", "d1", "d2"):
        rows = frame[frame["band_label"] == label].copy()
        direction = rows["direction"].astype(str)
        k_norm = pd.to_numeric(rows["k_norm"], errors="coerce").to_numpy(float)
        rows["plot_x"] = np.where(direction.eq("gamma_m"), -k_norm, k_norm)
        rows = rows.sort_values("plot_x").drop_duplicates("plot_x", keep="first")
        x = rows["plot_x"].to_numpy(float)
        y = pd.to_numeric(rows["re"], errors="coerce").to_numpy(float)
        matched = rows["match_status"].astype(str).ne("unmatched").to_numpy()
        character, color_map = _band_character(rows, label)
        points = np.column_stack([x, y])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        weights = 0.5 * (character[:-1] + character[1:])
        usable = matched[:-1] & matched[1:]
        usable &= np.all(np.isfinite(segments), axis=(1, 2)) & np.isfinite(weights)
        collection = LineCollection(
            segments[usable],
            cmap=color_map,
            norm=Normalize(0.0, 1.0),
            linewidth=2.0,
        )
        collection.set_array(weights[usable])
        axis.add_collection(collection)

    axis.set_xlim(-K_LIMIT, K_LIMIT)
    axis.set_ylim(*FREQUENCY_LIMITS)
    axis.set_yticks(FREQUENCY_TICKS)
    axis.set_xticks(
        [-K_LIMIT, 0.0, K_LIMIT],
        labels=[r"$0.3M\Gamma$", r"$\Gamma$", r"$0.3\Gamma K$"],
    )
    axis.get_xticklabels()[0].set_ha("left")
    axis.get_xticklabels()[-1].set_ha("right")
    axis.set_title(title, fontsize=TITLE_FONT_SIZE, pad=9)


def plot(unit_cell_overview: Path, finite_run: Path, output_path: Path) -> dict[str, object]:
    unit_cell_overview = unit_cell_overview.resolve(strict=True)
    finite_run = finite_run.resolve(strict=True)
    finite_overview = finite_run / "overview"
    band_columns = {
        "direction",
        "band_label",
        "k_norm",
        "re",
        "match_status",
        "px_weight",
        "py_weight",
        "dx_weight",
        "dy_weight",
    }
    cavity = _read_csv(unit_cell_overview / "cavity" / "selected_bands.csv", band_columns)
    cladding = _read_csv(
        unit_cell_overview / "cladding" / "selected_bands.csv",
        band_columns,
    )
    merged = _read_csv(
        finite_overview / "eigenfrequencies.csv",
        {"re", "q", "mode_idx", "symmetry_id"},
    )
    selected = _selected_py_mode(finite_run, merged)
    frequency = pd.to_numeric(merged["re"], errors="coerce").to_numpy(float)
    quality_factor = pd.to_numeric(merged["q"], errors="coerce").to_numpy(float)
    finite = np.isfinite(frequency) & np.isfinite(quality_factor)
    if int(finite.sum()) < 20:
        raise ValueError("Too few finite f-Q observations for a scatter plot")

    py_frequency = float(selected["frequency"])
    py_q = float(selected["q"])
    q_upper = max(1000.0, math.ceil(float(np.max(quality_factor[finite])) / 1000.0) * 1000.0)
    figure = plt.figure(figsize=(12.0, 6.47))
    grid = figure.add_gridspec(
        1,
        3,
        width_ratios=(1.065, 1.055, 1.0),
        left=0.08,
        right=0.98,
        bottom=0.13,
        top=0.84,
        wspace=0.02,
    )
    axes = [figure.add_subplot(grid[0, index]) for index in range(3)]
    axes[0].set_box_aspect(BAND_BOX_ASPECT)
    axes[1].set_box_aspect(BAND_BOX_ASPECT)
    axes[2].set_box_aspect(QUALITY_BOX_ASPECT)
    for axis in axes:
        axis.set_xscale("linear")
        axis.set_yscale("linear")
    _plot_bands(axes[0], cavity, "Cavity")
    _plot_bands(axes[1], cladding, "Cladding")
    axes[0].set_ylabel("Frequency (THz)", fontsize=FONT_SIZE, labelpad=12)
    axes[1].tick_params(labelleft=False)

    q_axis = axes[2]
    highlight_transform = offset_copy(
        q_axis.transData,
        fig=figure,
        x=-1.0,
        y=1.0,
        units="points",
    )
    q_axis.scatter(
        quality_factor[finite],
        frequency[finite],
        s=58,
        color="#075378",
        edgecolors="#06384f",
        linewidths=0.3,
        alpha=0.96,
        zorder=3,
    )
    q_axis.scatter(
        quality_factor[finite],
        frequency[finite],
        s=43,
        color=MODE_COLOR,
        edgecolors="none",
        alpha=0.98,
        zorder=4,
    )
    q_axis.scatter(
        quality_factor[finite],
        frequency[finite],
        s=2,
        color="#e5f8ff",
        edgecolors="none",
        alpha=0.7,
        transform=highlight_transform,
        zorder=5,
    )
    q_axis.scatter(
        [py_q],
        [py_frequency],
        s=110,
        color="#9f1e12",
        edgecolors="none",
        zorder=6,
    )
    q_axis.scatter(
        [py_q],
        [py_frequency],
        s=82,
        color=HIGHLIGHT_COLOR,
        edgecolors="none",
        zorder=7,
    )
    q_axis.scatter(
        [py_q],
        [py_frequency],
        s=3.5,
        color="#ffd5ca",
        edgecolors="none",
        alpha=0.75,
        transform=highlight_transform,
        zorder=8,
    )
    q_axis.text(
        0.96,
        0.94,
        r"$p_y$ fundamental mode",
        transform=q_axis.transAxes,
        ha="right",
        va="top",
        color=HIGHLIGHT_COLOR,
        fontsize=FONT_SIZE,
    )
    q_axis.set_xlim(0.0, q_upper)
    q_axis.set_ylim(*FREQUENCY_LIMITS)
    q_axis.set_xticks([0.0, 4000.0, 8000.0])
    q_axis.set_yticks(FREQUENCY_TICKS)
    q_axis.get_xticklabels()[0].set_ha("left")
    q_axis.get_xticklabels()[-1].set_ha("right")
    q_axis.set_title("Quality factor", fontsize=TITLE_FONT_SIZE, pad=9)
    q_axis.tick_params(labelleft=False)

    for axis in axes:
        axis.tick_params(
            direction="in",
            which="both",
            top=True,
            right=True,
            labelsize=FONT_SIZE,
        )
        axis.tick_params(axis="x", pad=6)
        axis.patch.set_alpha(0.0)
        for spine in axis.spines.values():
            spine.set_linewidth(1.0)
    figure.canvas.draw()
    line_y = figure.transFigure.inverted().transform(
        axes[0].transData.transform((0.0, py_frequency))
    )[1]
    figure.add_artist(
        Line2D(
            [axes[0].get_position().x0, axes[2].get_position().x1],
            [line_y, line_y],
            transform=figure.transFigure,
            color=HIGHLIGHT_COLOR,
            linestyle="--",
            linewidth=1.2,
            zorder=-0.5,
            clip_on=False,
        )
    )
    figure.legend(
        handles=[
            Line2D([], [], color=band.P_X_COLOR, lw=2.0, label=r"$p_y$"),
            Line2D([], [], color=band.P_Y_COLOR, lw=2.0, label=r"$p_x$"),
            Line2D([], [], color="#b0b0b0", lw=2.0, label=r"$d_{xy}$"),
            Line2D([], [], color="black", lw=2.0, label=r"$d_{x^2+y^2}$"),
        ],
        loc="upper right",
        bbox_to_anchor=(0.98, 0.985),
        ncol=4,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        handlelength=2.2,
        columnspacing=1.3,
    )

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.stem}.tmp.png")
    figure.savefig(temporary, dpi=300, facecolor="white")
    plt.close(figure)
    os.replace(temporary, output_path)
    return {
        "output_path": str(output_path),
        "plotted_mode_count": int(finite.sum()),
        "frequency_limits_thz": list(FREQUENCY_LIMITS),
        "q_limits": [0.0, q_upper],
        "x_scale": "linear",
        "y_scale": "linear",
        "py_frequency_thz": py_frequency,
        "py_q": py_q,
        "py_source_mode_idx": int(selected["mode_idx"]),
        "py_internal_component": "px",
        "py_user_component": "py",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("unit_cell_overview", type=Path)
    parser.add_argument("finite_run", type=Path)
    parser.add_argument("--output-name", default=OUTPUT_NAME)
    args = parser.parse_args()
    output = args.finite_run / "overview" / args.output_name
    print(json.dumps(plot(args.unit_cell_overview, args.finite_run, output), indent=2))


if __name__ == "__main__":
    main()
