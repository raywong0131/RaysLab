#!/usr/bin/env python3
"""Preview bulk-only regular-hex lattice edges and gap offsets."""

from __future__ import annotations

from tests._paths import output_root
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

from comsol_workflow.hex_lattice_utils import (
    cell_polygon as _cell_polygon,
    index_shell as index_shell_value,
    lattice_layer_continuous_shell as _lattice_layer_continuous_shell,
    normal_gap_to_shell,
    regular_hex_outline,
    shell_indices,
    shell_value,
    unit_cell_corners as _unit_cell_corners,
    uv_to_xy,
)
from comsol_workflow.polygon_utils import (
    cell_has_outside_region_part,
    exterior_boundary_segments,
    offset_polygon,
    ordered_boundary_loop,
    polygons_intersect,
)


OUT_DIR = output_root() / "test_bulk_cladding_lattice_generation"
DPI = 220
A_BULK = 0.82
BULK_RADIUS = 3
BACKGROUND_EXTRA_PERIODS = 4
CLADDING_LAYERS = 4
SMALL_GAP = 0.08
LARGE_GAP = 0.65
COLUMNS = ["region", "is_bulk", "i", "j", "u", "v", "x", "y"]
CLADDING_COLUMNS = ["region", "shell", "layer", "i", "j", "u", "v", "x", "y"]


@dataclass(frozen=True)
class PreviewCase:
    name: str
    gap: float


@dataclass(frozen=True)
class BoundaryMethod:
    slug: str
    title: str
    base_label: str
    expanded_label: str
    base_boundary: np.ndarray
    expanded_boundary: np.ndarray
    expanded_color: str


@dataclass(frozen=True)
class CladdingParameter:
    name: str
    a_cladding: float


@dataclass(frozen=True)
class CladdingSelection:
    first_shell: int
    retained_centers: pd.DataFrame


CASES = [
    PreviewCase("small_gap", gap=SMALL_GAP),
    PreviewCase("large_gap", gap=LARGE_GAP),
]

CLADDING_PARAMETERS = [
    CladdingParameter("bulk_large_period", a_cladding=0.62),
    CladdingParameter("similar_period", a_cladding=0.82),
    CladdingParameter("cladding_large_period", a_cladding=1.12),
]


def unit_cell_corners(period: float = A_BULK) -> np.ndarray:
    return _unit_cell_corners(period)


def cell_polygon(center_x: float, center_y: float, period: float = A_BULK) -> np.ndarray:
    return _cell_polygon(center_x, center_y, period)


def shell_edges(gap: float) -> dict[str, float]:
    bulk_center_shell = BULK_RADIUS * A_BULK
    bulk_envelope_edge_shell = bulk_center_shell + 0.5 * A_BULK
    envelope_gap_shell = bulk_envelope_edge_shell + normal_gap_to_shell(gap)
    background_shell = envelope_gap_shell + BACKGROUND_EXTRA_PERIODS * A_BULK
    return {
        "bulk_center": bulk_center_shell,
        "bulk_envelope_edge": bulk_envelope_edge_shell,
        "envelope_gap": envelope_gap_shell,
        "background": background_shell,
    }


def generate_bulk_lattice(case: PreviewCase) -> pd.DataFrame:
    edges = shell_edges(case.gap)
    index_radius = int(np.ceil(edges["background"] / A_BULK)) + 3
    axis = np.arange(-index_radius, index_radius + 1, dtype=int)
    i, j = np.meshgrid(axis, axis, indexing="xy")
    u = i.astype(float) * A_BULK
    v = j.astype(float) * A_BULK
    center_shell = shell_value(u, v)
    keep = center_shell <= edges["background"] + 1e-10
    is_bulk = center_shell[keep] <= edges["bulk_center"] + 1e-10
    x, y = uv_to_xy(u[keep], v[keep])
    return pd.DataFrame({
        "region": np.where(is_bulk, "bulk", "background"),
        "is_bulk": is_bulk,
        "i": i[keep].astype(float).ravel(),
        "j": j[keep].astype(float).ravel(),
        "u": u[keep].ravel(),
        "v": v[keep].ravel(),
        "x": np.asarray(x, dtype=float).ravel(),
        "y": np.asarray(y, dtype=float).ravel(),
    }, columns=COLUMNS)


def bulk_cell_polygons(points: pd.DataFrame) -> list[np.ndarray]:
    bulk = points.loc[points["is_bulk"]]
    return [cell_polygon(float(row.x), float(row.y)) for row in bulk.itertuples(index=False)]


def cladding_selection(period: float, expanded_boundary: np.ndarray, layer_count: int = CLADDING_LAYERS) -> CladdingSelection:
    max_radius = float(np.max(np.linalg.norm(expanded_boundary, axis=1)))
    max_shell = int(np.ceil(max_radius / period)) + layer_count + 6
    first_shell = None
    retained_rows = []

    for shell in range(max_shell + 1):
        shell_all_intersect = True
        for i, j in shell_indices(shell):
            x, y = uv_to_xy(i * period, j * period)
            polygon = cell_polygon(float(x), float(y), period)
            if not polygons_intersect(polygon, expanded_boundary):
                shell_all_intersect = False
            if cell_has_outside_region_part(polygon, expanded_boundary):
                retained_rows.append({
                    "region": "cladding",
                    "shell": float(shell),
                    "layer": np.nan,
                    "i": float(i),
                    "j": float(j),
                    "u": float(i * period),
                    "v": float(j * period),
                    "x": float(x),
                    "y": float(y),
                })
        if shell_all_intersect:
            first_shell = shell

    if first_shell is None:
        raise RuntimeError("Could not find a cladding shell intersecting the expanded boundary.")

    last_shell = first_shell + layer_count - 1
    retained_centers = pd.DataFrame(retained_rows, columns=CLADDING_COLUMNS)
    retained_centers = retained_centers.loc[retained_centers["shell"] <= last_shell].copy()
    retained_centers["layer"] = retained_centers["shell"] - first_shell + 1
    return CladdingSelection(first_shell=first_shell, retained_centers=retained_centers)


def generate_cladding_lattice(period: float, first_shell: int, layer_count: int = CLADDING_LAYERS) -> pd.DataFrame:
    last_shell = first_shell + layer_count - 1
    axis = np.arange(-last_shell, last_shell + 1, dtype=int)
    i, j = np.meshgrid(axis, axis, indexing="xy")
    shell = index_shell_value(i, j)
    keep = (shell >= first_shell) & (shell <= last_shell)
    u = i.astype(float) * period
    v = j.astype(float) * period
    x, y = uv_to_xy(u[keep], v[keep])
    return pd.DataFrame({
        "region": "cladding",
        "shell": shell[keep].astype(float).ravel(),
        "layer": (shell[keep] - first_shell + 1).astype(float).ravel(),
        "i": i[keep].astype(float).ravel(),
        "j": j[keep].astype(float).ravel(),
        "u": u[keep].ravel(),
        "v": v[keep].ravel(),
        "x": np.asarray(x, dtype=float).ravel(),
        "y": np.asarray(y, dtype=float).ravel(),
    }, columns=CLADDING_COLUMNS)


def nearest_spacing(points: pd.DataFrame) -> float:
    coords = points.loc[points["is_bulk"], ["x", "y"]].to_numpy(dtype=float)
    deltas = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    distances[distances <= 1e-10] = np.inf
    return float(np.min(distances))


def check_points(
    case: PreviewCase,
    points: pd.DataFrame,
    boundary: np.ndarray,
    gap_boundary: np.ndarray,
    lattice_shell: np.ndarray,
) -> dict[str, float | int | str]:
    assert list(points.columns) == COLUMNS
    assert set(points["region"]) == {"bulk", "background"}
    assert not points.duplicated(["x", "y"]).any()
    assert np.isclose(nearest_spacing(points), A_BULK)
    assert len(boundary) > 0
    assert len(gap_boundary) == len(boundary)
    assert len(lattice_shell) > 0

    return {
        "case": case.name,
        "a_bulk": A_BULK,
        "gap": case.gap,
        "bulk_radius": BULK_RADIUS,
        "bulk_count": int(points["is_bulk"].sum()),
        "background_count": int((~points["is_bulk"]).sum()),
        "bulk_nearest_spacing": nearest_spacing(points),
        "true_boundary_vertices": len(boundary),
        "lattice_shell_vertices": len(lattice_shell),
    }


def draw_closed_line(ax, vertices: np.ndarray, **kwargs) -> None:
    closed = np.vstack([vertices, vertices[0]])
    ax.plot(closed[:, 0], closed[:, 1], **kwargs)


def draw_lattice_cells(ax, points: pd.DataFrame) -> None:
    background_segments = []
    bulk_patches = []
    for point in points.itertuples(index=False):
        polygon = cell_polygon(float(point.x), float(point.y))
        if bool(point.is_bulk):
            bulk_patches.append(Polygon(
                polygon,
                closed=True,
                facecolor="#dbeafe",
                edgecolor="#60a5fa",
                linewidth=0.45,
                alpha=0.88,
                zorder=2,
            ))
        else:
            for idx, start in enumerate(polygon):
                background_segments.append([start, polygon[(idx + 1) % len(polygon)]])

    if background_segments:
        ax.add_collection(LineCollection(background_segments, colors="#d1d5db", linewidths=0.28, alpha=0.55, zorder=1))
    for patch in bulk_patches:
        ax.add_patch(patch)


def draw_inner_bulk_cells(ax, points: pd.DataFrame) -> None:
    bulk = points.loc[points["is_bulk"]]
    for point in bulk.itertuples(index=False):
        ax.add_patch(Polygon(
            cell_polygon(float(point.x), float(point.y)),
            closed=True,
            facecolor="#dbeafe",
            edgecolor="#60a5fa",
            linewidth=0.45,
            alpha=0.9,
            zorder=6,
        ))
    coords = bulk[["x", "y"]].to_numpy(dtype=float)
    ax.scatter(coords[:, 0], coords[:, 1], s=8, c="#2563eb", label="inner bulk centers", linewidths=0, alpha=0.9, zorder=7)


def draw_cladding_cells(ax, cladding: pd.DataFrame, period: float) -> None:
    patches = []
    for point in cladding.itertuples(index=False):
        patches.append(Polygon(
            cell_polygon(float(point.x), float(point.y), period),
            closed=True,
            facecolor="#dcfce7",
            edgecolor="#22c55e",
            linewidth=0.42,
            alpha=0.72,
            zorder=1,
        ))
    for patch in patches:
        ax.add_patch(patch)

    coords = cladding[["x", "y"]].to_numpy(dtype=float)
    ax.scatter(coords[:, 0], coords[:, 1], s=5, c="#16a34a", label="cladding centers", linewidths=0, alpha=0.65, zorder=2)


def boundary_methods(
    case: PreviewCase,
    true_boundary: np.ndarray,
    true_gap_boundary: np.ndarray,
    lattice_shell: np.ndarray,
) -> list[BoundaryMethod]:
    edges = shell_edges(case.gap)
    return [
        BoundaryMethod(
            slug="envelope_gap",
            title="envelope offset",
            base_label="envelope bulk edge",
            expanded_label="envelope edge + gap",
            base_boundary=regular_hex_outline(edges["bulk_envelope_edge"]),
            expanded_boundary=regular_hex_outline(edges["envelope_gap"]),
            expanded_color="#ef4444",
        ),
        BoundaryMethod(
            slug="miter_offset",
            title="miter offset",
            base_label="true cell-union edge",
            expanded_label="true edge + gap",
            base_boundary=true_boundary,
            expanded_boundary=true_gap_boundary,
            expanded_color="#a855f7",
        ),
        BoundaryMethod(
            slug="lattice_layer_shell",
            title="lattice-layer continuous shell",
            base_label="true cell-union edge",
            expanded_label="lattice-layer continuous shell",
            base_boundary=true_boundary,
            expanded_boundary=lattice_shell,
            expanded_color="#f97316",
        ),
    ]


def setup_cladding_axes(
    ax,
    case: PreviewCase,
    points: pd.DataFrame,
    method: BoundaryMethod,
    cladding_parameter: CladdingParameter,
    *,
    annotated: bool,
) -> tuple[dict[str, float | int | str], pd.DataFrame]:
    selection = cladding_selection(cladding_parameter.a_cladding, method.expanded_boundary)
    first_shell = selection.first_shell
    cladding = generate_cladding_lattice(cladding_parameter.a_cladding, first_shell)
    retained_centers = selection.retained_centers.copy()

    draw_cladding_cells(ax, cladding, cladding_parameter.a_cladding)
    ax.add_patch(Polygon(
        method.expanded_boundary,
        closed=True,
        facecolor="white",
        edgecolor="none",
        linewidth=0.0,
        zorder=4,
    ))
    draw_inner_bulk_cells(ax, points)

    draw_closed_line(
        ax,
        method.base_boundary,
        color="#111827",
        linestyle="-",
        linewidth=1.35,
        label=method.base_label,
        zorder=8,
    )
    draw_closed_line(
        ax,
        method.expanded_boundary,
        color=method.expanded_color,
        linestyle="-",
        linewidth=1.65,
        label=method.expanded_label,
        zorder=9,
    )

    cladding_xy = cladding[["x", "y"]].to_numpy(dtype=float)
    bulk_xy = points.loc[points["is_bulk"], ["x", "y"]].to_numpy(dtype=float)
    xy = np.vstack([cladding_xy, bulk_xy, method.base_boundary, method.expanded_boundary])
    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), A_BULK)
    pad = 0.06 * span + max(A_BULK, cladding_parameter.a_cladding)
    ax.set_xlim(float(xy[:, 0].min() - pad), float(xy[:, 0].max() + pad))
    ax.set_ylim(float(xy[:, 1].min() - pad), float(xy[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.16)
    ax.set_title(f"{case.name}: {cladding_parameter.name}, {method.title}", pad=10)

    summary = {
        "case": case.name,
        "method": method.slug,
        "parameter": cladding_parameter.name,
        "a_bulk": A_BULK,
        "a_cladding": cladding_parameter.a_cladding,
        "gap": case.gap,
        "first_cladding_shell": first_shell,
        "cladding_layers": CLADDING_LAYERS,
        "cladding_cells": int(len(cladding)),
        "retained_center_count": int(len(retained_centers)),
    }
    context = {
        "case": case.name,
        "method": method.slug,
        "parameter": cladding_parameter.name,
        "a_bulk": A_BULK,
        "a_cladding": cladding_parameter.a_cladding,
        "gap": case.gap,
        "first_cladding_shell": first_shell,
        "cladding_layers": CLADDING_LAYERS,
    }
    for key, value in context.items():
        retained_centers[key] = value
    retained_centers = retained_centers[[*context.keys(), *CLADDING_COLUMNS]]

    if annotated:
        info = (
            f"a_bulk = {A_BULK:g}\n"
            f"a_cladding = {cladding_parameter.a_cladding:g}\n"
            f"gap = {case.gap:g}\n"
            f"first shell = {first_shell}\n"
            f"cladding layers = {CLADDING_LAYERS}\n"
            f"cladding cells = {len(cladding)}"
        )
        ax.text(
            0.02,
            0.02,
            info,
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
            ha="left",
            bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#999999", "alpha": 0.92},
            zorder=10,
        )
    return summary, retained_centers


def setup_method_axes(
    ax,
    case: PreviewCase,
    points: pd.DataFrame,
    method: BoundaryMethod,
    *,
    annotated: bool,
) -> None:
    draw_lattice_cells(ax, points)

    background = points.loc[~points["is_bulk"], ["x", "y"]].to_numpy(dtype=float)
    bulk = points.loc[points["is_bulk"], ["x", "y"]].to_numpy(dtype=float)
    ax.scatter(background[:, 0], background[:, 1], s=3, c="#9ca3af", label="background bulk-period lattice", linewidths=0, alpha=0.45, zorder=3)
    ax.scatter(bulk[:, 0], bulk[:, 1], s=8, c="#2563eb", label="highlighted bulk centers", linewidths=0, alpha=0.9, zorder=4)

    draw_closed_line(
        ax,
        method.base_boundary,
        color="#111827",
        linestyle="-",
        linewidth=1.35,
        label=method.base_label,
        zorder=7,
    )
    draw_closed_line(
        ax,
        method.expanded_boundary,
        color=method.expanded_color,
        linestyle="-",
        linewidth=1.65,
        label=method.expanded_label,
        zorder=7,
    )

    xy = np.vstack([points[["x", "y"]].to_numpy(dtype=float), method.base_boundary, method.expanded_boundary])
    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), A_BULK)
    pad = 0.06 * span + A_BULK
    ax.set_xlim(float(xy[:, 0].min() - pad), float(xy[:, 0].max() + pad))
    ax.set_ylim(float(xy[:, 1].min() - pad), float(xy[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.16)
    ax.set_title(f"{case.name}: {method.title}", pad=10)

    if annotated:
        info = (
            f"a_bulk = {A_BULK:g}\n"
            f"gap = {case.gap:g}\n"
            f"bulk centers = {int(len(bulk))}\n"
            f"background centers = {int(len(background))}\n"
            f"base vertices = {len(method.base_boundary)}\n"
            f"expanded vertices = {len(method.expanded_boundary)}"
        )
        ax.text(
            0.02,
            0.02,
            info,
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
            ha="left",
            bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#999999", "alpha": 0.92},
            zorder=10,
        )


def plot_method_single(
    path: Path,
    case: PreviewCase,
    points: pd.DataFrame,
    method: BoundaryMethod,
    *,
    annotated: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    setup_method_axes(ax, case, points, method, annotated=annotated)
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def plot_method_comparison(
    path: Path,
    generated: list[tuple[PreviewCase, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]],
    method_slug: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(generated), figsize=(7.0 * len(generated), 7.2), constrained_layout=True)
    method_title = ""
    for ax, (case, points, true_boundary, true_gap_boundary, lattice_shell) in zip(np.atleast_1d(axes), generated):
        methods = {method.slug: method for method in boundary_methods(case, true_boundary, true_gap_boundary, lattice_shell)}
        method = methods[method_slug]
        method_title = method.title
        setup_method_axes(ax, case, points, method, annotated=True)
    fig.suptitle(f"bulk-only {method_title}", fontsize=14)
    axes[0].legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def plot_cladding_comparison(
    path: Path,
    generated: list[tuple[PreviewCase, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]],
    cladding_parameter: CladdingParameter,
    method_slug: str,
) -> tuple[list[dict[str, float | int | str]], list[pd.DataFrame]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(generated), figsize=(7.0 * len(generated), 7.2), constrained_layout=True)
    summaries = []
    retained_centers = []
    method_title = ""
    for ax, (case, points, true_boundary, true_gap_boundary, lattice_shell) in zip(np.atleast_1d(axes), generated):
        methods = {method.slug: method for method in boundary_methods(case, true_boundary, true_gap_boundary, lattice_shell)}
        method = methods[method_slug]
        method_title = method.title
        summary, centers = setup_cladding_axes(ax, case, points, method, cladding_parameter, annotated=True)
        summaries.append(summary)
        retained_centers.append(centers)
    fig.suptitle(f"cladding {cladding_parameter.name}: {method_title}", fontsize=14)
    axes[0].legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return summaries, retained_centers


def case_stem(case: PreviewCase) -> str:
    return f"{case.name}_a_bulk_{A_BULK:g}_gap_{case.gap:g}".replace(".", "p")


def reset_output_dir() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.glob("*"):
        if path.is_file():
            path.unlink()


def main() -> None:
    reset_output_dir()
    summaries = []
    generated = []

    for case in CASES:
        points = generate_bulk_lattice(case)
        polygons = bulk_cell_polygons(points)
        true_boundary_segments = exterior_boundary_segments(polygons)
        true_boundary = ordered_boundary_loop(true_boundary_segments)
        true_gap_boundary = offset_polygon(true_boundary, case.gap)
        lattice_shell = _lattice_layer_continuous_shell(BULK_RADIUS, A_BULK, case.gap)

        summaries.append(check_points(case, points, true_boundary, true_gap_boundary, lattice_shell))
        generated.append((case, points, true_boundary, true_gap_boundary, lattice_shell))

        stem = case_stem(case)
        points.to_csv(OUT_DIR / f"{stem}_centers.csv", index=False)
        for method in boundary_methods(case, true_boundary, true_gap_boundary, lattice_shell):
            plot_method_single(
                OUT_DIR / f"{stem}_{method.slug}.png",
                case,
                points,
                method,
                annotated=True,
            )

        print(
            f"{case.name}: "
            f"bulk={summaries[-1]['bulk_count']}, "
            f"background={summaries[-1]['background_count']}, "
            f"a_bulk={A_BULK:g}, gap={case.gap:g}, "
            f"true_edge_vertices={summaries[-1]['true_boundary_vertices']}, "
            f"lattice_shell_vertices={summaries[-1]['lattice_shell_vertices']}"
        )

    for method_slug in ["envelope_gap", "miter_offset", "lattice_layer_shell"]:
        plot_method_comparison(OUT_DIR / f"{method_slug}_comparison.png", generated, method_slug)
    cladding_summaries = []
    cladding_retained_centers = []
    for cladding_parameter in CLADDING_PARAMETERS:
        for method_slug in ["envelope_gap", "miter_offset", "lattice_layer_shell"]:
            path = OUT_DIR / f"cladding_{cladding_parameter.name}_{method_slug}_comparison.png"
            method_summaries, method_retained_centers = plot_cladding_comparison(path, generated, cladding_parameter, method_slug)
            cladding_summaries.extend(method_summaries)
            cladding_retained_centers.extend(method_retained_centers)
    pd.DataFrame(summaries).to_csv(OUT_DIR / "summary.csv", index=False)
    pd.DataFrame(cladding_summaries).to_csv(OUT_DIR / "cladding_summary.csv", index=False)
    if cladding_retained_centers:
        pd.concat(cladding_retained_centers, ignore_index=True).to_csv(OUT_DIR / "cladding_retained_centers.csv", index=False)
    print(f"\nWrote lattice edge previews to {OUT_DIR}")


if __name__ == "__main__":
    main()
