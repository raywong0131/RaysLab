#!/usr/bin/env python3
"""Generate and check common bulk/side/corner finite lattice patches."""

from __future__ import annotations

from tests._paths import output_root
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]

from comsol_workflow.finite_patch_lattice import (
    BasisName,
    BulkSideCornerConfig,
    COLUMNS,
    centered_axis,
    generate_lattice,
    lattice_unit_cell_corners,
    oblique_box_corners,
    part_intervals,
    unit_cell_edges,
    unit_cell_spacing_for_part,
    uv_to_xy,
)
from comsol_workflow.polygon_utils import polygon_area


OUT_DIR = output_root() / "test_bulk_side_corner_lattice_generation"
DPI = 220


@dataclass(frozen=True)
class PreviewCase:
    name: str
    basis: BasisName
    config: BulkSideCornerConfig


CASES = [
    PreviewCase(
        "square_small_gap",
        "orthogonal",
        BulkSideCornerConfig(a=0.82, b=1.12, gap=0.08, bulk_periods=5, boundary_periods=3),
    ),
    PreviewCase(
        "square_large_gap",
        "orthogonal",
        BulkSideCornerConfig(a=0.82, b=1.12, gap=0.58, bulk_periods=5, boundary_periods=3),
    ),
    PreviewCase(
        "hex_small_gap",
        "hex",
        BulkSideCornerConfig(a=0.82, b=1.12, gap=0.08, bulk_periods=5, boundary_periods=3),
    ),
    PreviewCase(
        "hex_large_gap",
        "hex",
        BulkSideCornerConfig(a=0.82, b=1.12, gap=0.58, bulk_periods=5, boundary_periods=3),
    ),
]

COLORS = {
    "bulk": "#2563eb",
    "side": "#d97706",
    "corner": "#9333ea",
}


def reset_output_dir() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.glob("*"):
        if path.is_file():
            path.unlink()


def case_stem(case: PreviewCase) -> str:
    config = case.config
    raw = (
        f"{case.name}_a_{config.a:g}_b_{config.b:g}_gap_{config.gap:g}_"
        f"bulk_{config.bulk_periods}_boundary_{config.boundary_periods}"
    )
    return raw.replace(".", "p")


def expected_counts(config: BulkSideCornerConfig) -> dict[str, int]:
    n_bulk = config.bulk_periods**2
    n_side = 4 * config.boundary_periods * config.bulk_periods
    n_corner = 4 * config.boundary_periods**2
    return {
        "bulk": n_bulk,
        "side": n_side,
        "corner": n_corner,
        "total": n_bulk + n_side + n_corner,
    }


def check_lattice(case: PreviewCase, points: pd.DataFrame) -> dict[str, int | float | str]:
    config = case.config
    counts = expected_counts(config)
    assert list(points.columns) == COLUMNS
    assert len(points) == counts["total"]
    assert points["region"].value_counts().to_dict() == {
        "bulk": counts["bulk"],
        "side": counts["side"],
        "corner": counts["corner"],
    }
    assert not points.duplicated(["part", "u", "v"]).any()

    x_expected, y_expected = uv_to_xy(points["u"].to_numpy(), points["v"].to_numpy(), case.basis)
    assert np.allclose(points["x"].to_numpy(), x_expected)
    assert np.allclose(points["y"].to_numpy(), y_expected)

    _, centered = centered_axis(config.bulk_periods, config.a)
    bulk = points.loc[points["part"].eq("bulk")]
    assert np.allclose(np.sort(bulk["u"].unique()), centered)
    assert np.allclose(np.sort(bulk["v"].unique()), centered)

    for part in points["part"].unique():
        spacing_u, spacing_v = unit_cell_spacing_for_part(str(part), config)
        cell = lattice_unit_cell_corners(spacing_u, spacing_v, case.basis)
        assert len(cell) >= 4
        assert polygon_area(cell) > 0.0

    inner_edge, outer_edge = unit_cell_edges(config)
    assert outer_edge > inner_edge

    return {
        "case": case.name,
        "basis": case.basis,
        "a": config.a,
        "b": config.b,
        "gap": config.gap,
        "bulk_periods": config.bulk_periods,
        "boundary_periods": config.boundary_periods,
        "bulk_count": counts["bulk"],
        "side_count": counts["side"],
        "corner_count": counts["corner"],
        "total_count": counts["total"],
        "inner_edge": inner_edge,
        "outer_edge": outer_edge,
    }


def draw_closed(ax, vertices: np.ndarray, **kwargs) -> None:
    closed = np.vstack([vertices, vertices[0]])
    ax.plot(closed[:, 0], closed[:, 1], **kwargs)


def draw_part_boxes(ax, case: PreviewCase) -> None:
    for region, part, u0, u1, v0, v1 in part_intervals(case.config):
        color = COLORS[region]
        polygon = oblique_box_corners(u0, u1, v0, v1, case.basis)
        ax.add_patch(Polygon(
            polygon,
            closed=True,
            facecolor=color,
            edgecolor=color,
            linewidth=0.8,
            alpha=0.08,
            zorder=1,
        ))
        center_x, center_y = uv_to_xy((u0 + u1) / 2.0, (v0 + v1) / 2.0, case.basis)
        ax.text(float(center_x), float(center_y), part, color=color, ha="center", va="center", fontsize=7, zorder=7)


def draw_lattice_cells(ax, case: PreviewCase, points: pd.DataFrame) -> None:
    config = case.config
    cell_cache = {
        part: lattice_unit_cell_corners(*unit_cell_spacing_for_part(str(part), config), case.basis)
        for part in points["part"].unique()
    }
    for row in points.itertuples(index=False):
        color = COLORS[str(row.region)]
        corners = cell_cache[str(row.part)] + np.array([row.x, row.y], dtype=float)
        ax.add_patch(Polygon(
            corners,
            closed=True,
            facecolor=color,
            edgecolor=color,
            linewidth=0.35,
            alpha=0.15 if row.region != "bulk" else 0.22,
            zorder=2,
        ))


def annotate_axes(ax, case: PreviewCase, points: pd.DataFrame) -> None:
    config = case.config
    _, coords = centered_axis(config.bulk_periods, config.a)
    half_span = (config.bulk_periods - 1) * config.a / 2.0
    center_u = float(coords[int(np.argmin(np.abs(coords)))])

    bulk_x, bulk_y = uv_to_xy(center_u, half_span, case.basis)
    side_x, side_y = uv_to_xy(center_u, half_span + config.gap, case.basis)
    ax.plot([bulk_x, side_x], [bulk_y, side_y], color="#dc2626", linewidth=1.2, zorder=8)
    ax.scatter([bulk_x, side_x], [bulk_y, side_y], s=8, c="#dc2626", zorder=9)
    ax.annotate(
        f"gap = {config.gap:g}",
        xy=((float(bulk_x) + float(side_x)) / 2.0, (float(bulk_y) + float(side_y)) / 2.0),
        xytext=(8, 16),
        textcoords="offset points",
        color="#b91c1c",
        fontsize=9,
    )

    bulk = points.loc[points["region"].eq("bulk"), ["x", "y"]].to_numpy(dtype=float)
    origin_index = int(np.argmin(np.linalg.norm(bulk, axis=1)))
    origin = bulk[origin_index]
    distances = np.linalg.norm(bulk - origin, axis=1)
    distances[origin_index] = np.inf
    neighbor = bulk[int(np.argmin(distances))]
    ax.plot([origin[0], neighbor[0]], [origin[1], neighbor[1]], color="#1d4ed8", linewidth=1.2, zorder=8)
    ax.scatter([origin[0], neighbor[0]], [origin[1], neighbor[1]], s=8, c="#1d4ed8", zorder=9)
    ax.text((origin[0] + neighbor[0]) / 2.0, (origin[1] + neighbor[1]) / 2.0, f"a = {config.a:g}", color="#1d4ed8", fontsize=9)

    top_v = half_span + config.gap
    b_u = center_u + config.a
    b_start_x, b_start_y = uv_to_xy(b_u, top_v, case.basis)
    b_end_x, b_end_y = uv_to_xy(b_u, top_v + config.b, case.basis)
    ax.plot([b_start_x, b_end_x], [b_start_y, b_end_y], color="#b45309", linewidth=1.2, zorder=8)
    ax.scatter([b_start_x, b_end_x], [b_start_y, b_end_y], s=8, c="#b45309", zorder=9)
    ax.text(float(b_end_x), float(b_end_y), f"b = {config.b:g}", color="#92400e", fontsize=9)

    counts = expected_counts(config)
    info = (
        f"basis: {case.basis}\n"
        f"a: {config.a:g}\n"
        f"b: {config.b:g}\n"
        f"gap: {config.gap:g}\n"
        f"bulk: {counts['bulk']}\n"
        f"side: {counts['side']}\n"
        f"corner: {counts['corner']}"
    )
    ax.text(
        0.02,
        0.02,
        info,
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        ha="left",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#999999", "alpha": 0.92},
        zorder=10,
    )


def plot_lattice(path: Path, case: PreviewCase, points: pd.DataFrame, *, annotated: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 8.6))
    if annotated:
        draw_part_boxes(ax, case)
    draw_lattice_cells(ax, case, points)

    for region, color, label in [
        ("corner", COLORS["corner"], "corner(b,b)"),
        ("side", COLORS["side"], "side(a,b)"),
        ("bulk", COLORS["bulk"], "bulk(a,a)"),
    ]:
        subset = points.loc[points["region"].eq(region), ["x", "y"]].to_numpy(dtype=float)
        ax.scatter(subset[:, 0], subset[:, 1], s=8, c=color, label=label, linewidths=0, zorder=4)

    inner_edge, outer_edge = unit_cell_edges(case.config)
    inner = oblique_box_corners(-inner_edge, inner_edge, -inner_edge, inner_edge, case.basis)
    outer = oblique_box_corners(-outer_edge, outer_edge, -outer_edge, outer_edge, case.basis)
    draw_closed(ax, outer, color="#d97706", linewidth=1.2, label="outer cell boundary", zorder=6)
    draw_closed(ax, inner, color="#2563eb", linewidth=1.2, label="bulk cell boundary", zorder=7)

    if annotated:
        annotate_axes(ax, case, points)

    xy = np.vstack([points[["x", "y"]].to_numpy(dtype=float), inner, outer])
    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), case.config.a)
    pad = 0.08 * span + max(case.config.a, case.config.b)
    ax.set_xlim(float(xy[:, 0].min() - pad), float(xy[:, 0].max() + pad))
    ax.set_ylim(float(xy[:, 1].min() - pad), float(xy[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"{case.name}: common bulk/side/corner lattice", pad=10)
    ax.grid(True, alpha=0.16)
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main() -> None:
    reset_output_dir()
    summaries = []
    for case in CASES:
        points = generate_lattice(case.config, case.basis)
        summary = check_lattice(case, points)
        summaries.append(summary)

        stem = case_stem(case)
        points.to_csv(OUT_DIR / f"{stem}_centers.csv", index=False)
        plot_lattice(OUT_DIR / f"{stem}.png", case, points, annotated=False)
        plot_lattice(OUT_DIR / f"{stem}_annotated.png", case, points, annotated=True)
        print(
            f"{case.name}: total={summary['total_count']}, "
            f"bulk={summary['bulk_count']}, side={summary['side_count']}, "
            f"corner={summary['corner_count']}"
        )

    pd.DataFrame(summaries).to_csv(OUT_DIR / "summary.csv", index=False)
    print(f"wrote common lattice previews to {OUT_DIR}")


if __name__ == "__main__":
    main()
