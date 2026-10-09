#!/usr/bin/env python3
"""Preview finite cavity hole polygons clipped by a miter-offset boundary."""

from __future__ import annotations

from tests._paths import output_root
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis
from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.hex_lattice_utils import (
    LatticePoint,
    boundary_from_lattice_cells,
    bulk_points as _bulk_points,
    cell_polygon as _cell_polygon,
    lattice_points as _lattice_points,
    retained_cladding_points as _retained_cladding_points,
    unit_cell_corners as _unit_cell_corners,
)
from comsol_workflow.polygon_utils import (
    clip_polygon_to_outside_region,
    offset_polygon,
    validate_clipped_hole,
)


OUT_DIR = output_root() / "test_cavity_hole_polygon_clipping"
DPI = 220

A = 0.82
R_0 = A / 3.0
B_0 = 0.23
BULK_RADIUS = 3
GEOMETRIC_GAPS = [0.05, 0.20, 0.35, 0.50, 0.65]
CLADDING_LAYERS = 4

BULK_PARAMS = {
    "name": "cavity_p_bic",
    "r_f0": 0.97,
    "b_square_f0": 1.0,
    "b_square_f3": 0.123,
}
CLADDING_PARAMS = {
    "name": "bulk_gap_centered",
    "r_f0": 0.95,
    "b_square_f0": 0.95,
}


def get_hole_params(params: dict[str, float]) -> np.ndarray:
    r_f0 = params["r_f0"]
    b_square_f0 = params["b_square_f0"]

    r_f = np.zeros(6, dtype=float)
    theta_f = np.zeros(6, dtype=float)
    b_square_f = np.zeros(6, dtype=float)
    phi_f = np.zeros(6, dtype=float)

    r_f[0] = 1.0
    b_square_f[0] = 1.0

    for name, value in params.items():
        if name in {"name", "r_f0", "b_square_f0"}:
            continue
        if name.startswith("r_f"):
            r_f[int(name.removeprefix("r_f"))] = float(value)
        elif name.startswith("theta_f"):
            theta_f[int(name.removeprefix("theta_f"))] = float(value)
        elif name.startswith("b_square_f"):
            b_square_f[int(name.removeprefix("b_square_f"))] = float(value)
        elif name.startswith("phi_f"):
            phi_f[int(name.removeprefix("phi_f"))] = float(value)

    r_f = R_0 * r_f0 * r_f
    b_square_f = B_0 * B_0 * b_square_f0 * b_square_f

    fourier_min_values = get_fourier_basis_min_nonzero_values(6)
    r = fourier_to_standard_basis(r_f / fourier_min_values, 6)
    theta = fourier_to_standard_basis(theta_f / fourier_min_values, 6)
    b_square = fourier_to_standard_basis(b_square_f / fourier_min_values, 6)
    b = np.sqrt(np.maximum(b_square, 0.0))
    phi = fourier_to_standard_basis(phi_f / fourier_min_values, 6)
    return np.stack([r, theta, b, phi], axis=1)


def create_unit_cell_holes(params: dict[str, float]) -> list[np.ndarray]:
    _, holes, _ = create_hexagon_design(A, get_hole_params(params))
    return list(holes)


def unit_cell_corners(period: float = A) -> np.ndarray:
    return _unit_cell_corners(period)


def cell_polygon(center_x: float, center_y: float, period: float = A) -> np.ndarray:
    return _cell_polygon(center_x, center_y, period)


def lattice_points(max_shell: int, period: float = A) -> list[LatticePoint]:
    return _lattice_points(max_shell, period)


def bulk_points() -> list[LatticePoint]:
    return _bulk_points(BULK_RADIUS, A)


def bulk_boundary(points: list[LatticePoint]) -> np.ndarray:
    return boundary_from_lattice_cells(points, A)


def retained_cladding_points(boundary: np.ndarray) -> list[LatticePoint]:
    _, retained = _retained_cladding_points(A, boundary, CLADDING_LAYERS)
    return retained


def translate_polygon(polygon: np.ndarray, x: float, y: float) -> np.ndarray:
    return np.asarray(polygon, dtype=float) + np.array([x, y], dtype=float)


def build_simulation_holes(gap: float) -> tuple[list[np.ndarray], dict[str, object]]:
    bulk = bulk_points()
    inner_boundary = bulk_boundary(bulk)
    expanded_boundary = offset_polygon(inner_boundary, gap)
    cladding = retained_cladding_points(expanded_boundary)

    bulk_unit_holes = create_unit_cell_holes(BULK_PARAMS)
    cladding_unit_holes = create_unit_cell_holes(CLADDING_PARAMS)

    holes: list[np.ndarray] = []
    original_holes: list[np.ndarray] = []
    bulk_hole_count = 0
    cladding_hole_count = 0
    clipped_hole_count = 0
    removed_cladding_holes = 0

    for point in bulk:
        for hole in bulk_unit_holes:
            translated = translate_polygon(hole, point.x, point.y)
            original_holes.append(translated)
            holes.append(translated)
            bulk_hole_count += 1

    for point in cladding:
        for hole in cladding_unit_holes:
            translated = translate_polygon(hole, point.x, point.y)
            original_holes.append(translated)
            clipped = clip_polygon_to_outside_region(translated, expanded_boundary)
            validate_clipped_hole(translated, clipped, expanded_boundary)
            if not clipped:
                removed_cladding_holes += 1
                continue
            was_clipped = len(clipped) != 1
            if not was_clipped:
                piece = clipped[0]
                was_clipped = piece.shape != translated.shape or not np.allclose(piece, translated, atol=1e-9)
            if was_clipped:
                clipped_hole_count += 1
            for piece in clipped:
                holes.append(piece)
                cladding_hole_count += 1

    metadata = {
        "a": A,
        "bulk_radius": BULK_RADIUS,
        "geometric_gap": gap,
        "cladding_layers": CLADDING_LAYERS,
        "bulk_params": BULK_PARAMS,
        "cladding_params": CLADDING_PARAMS,
        "bulk_cells": len(bulk),
        "cladding_candidate_cells": len(cladding),
        "bulk_holes": bulk_hole_count,
        "cladding_holes": cladding_hole_count,
        "clipped_cladding_holes": clipped_hole_count,
        "removed_cladding_holes": removed_cladding_holes,
        "total_hole_polygons": len(holes),
        "inner_boundary": inner_boundary,
        "expanded_boundary": expanded_boundary,
        "bulk_points": bulk,
        "cladding_points": cladding,
        "original_holes": original_holes,
    }
    return holes, metadata


def closed_line(vertices: np.ndarray) -> np.ndarray:
    return np.vstack([vertices, vertices[0]])


def add_hole_patches(
    ax,
    holes: list[np.ndarray],
    *,
    facecolor: str,
    edgecolor: str,
    linewidth: float,
    alpha: float,
    zorder: int,
    label: str | None = None,
) -> None:
    for idx, hole in enumerate(holes):
        ax.add_patch(Polygon(
            hole,
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
            label=label if idx == 0 else "_nolegend_",
        ))


def limits_from_polygons(polygons: list[np.ndarray], pad_scale: float = 0.08) -> tuple[float, float, float, float]:
    pts = np.vstack(polygons)
    span = max(float(np.ptp(pts[:, 0])), float(np.ptp(pts[:, 1])), A)
    pad = pad_scale * span
    return (
        float(pts[:, 0].min() - pad),
        float(pts[:, 0].max() + pad),
        float(pts[:, 1].min() - pad),
        float(pts[:, 1].max() + pad),
    )


def apply_limits(ax, limits: tuple[float, float, float, float]) -> None:
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.16)


def set_limits_from_polygons(ax, polygons: list[np.ndarray], pad_scale: float = 0.08) -> None:
    apply_limits(ax, limits_from_polygons(polygons, pad_scale=pad_scale))


def draw_geometry_panel(ax, holes: list[np.ndarray], metadata: dict[str, object], limits: tuple[float, float, float, float]) -> None:
    add_hole_patches(ax, holes, facecolor="#111827", edgecolor="#111827", linewidth=0.35, alpha=0.9, zorder=3)
    apply_limits(ax, limits)
    ax.set_title(f"gap = {metadata['geometric_gap']:g} um")


def draw_geometry(path: Path, holes: list[np.ndarray], metadata: dict[str, object]) -> None:
    limits = limits_from_polygons([*holes, metadata["expanded_boundary"]])
    fig, ax = plt.subplots(figsize=(8.0, 8.0), constrained_layout=True)
    draw_geometry_panel(ax, holes, metadata, limits)
    fig.suptitle("simulation-compatible clipped hole polygons", fontsize=14)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def draw_unit_cells(ax, points: list[LatticePoint], *, facecolor: str, edgecolor: str, alpha: float, zorder: int, label: str) -> None:
    for point in points:
        ax.add_patch(Polygon(cell_polygon(point.x, point.y), closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=0.35, alpha=alpha, zorder=zorder))
    coords = np.array([[point.x, point.y] for point in points], dtype=float)
    ax.scatter(coords[:, 0], coords[:, 1], s=5, c=edgecolor, linewidths=0, alpha=0.7, label=label, zorder=zorder + 1)


def draw_annotated_panel(ax, holes: list[np.ndarray], metadata: dict[str, object], limits: tuple[float, float, float, float]) -> None:
    bulk_points_list = metadata["bulk_points"]
    cladding_points_list = metadata["cladding_points"]
    inner_boundary = metadata["inner_boundary"]
    expanded_boundary = metadata["expanded_boundary"]
    original_holes = metadata["original_holes"]

    draw_unit_cells(ax, bulk_points_list, facecolor="#dbeafe", edgecolor="#3b82f6", alpha=0.55, zorder=1, label="bulk p-BIC cells")
    draw_unit_cells(ax, cladding_points_list, facecolor="#dcfce7", edgecolor="#22c55e", alpha=0.38, zorder=0, label="outside bandgap cells")
    ax.add_patch(Polygon(expanded_boundary, closed=True, facecolor="white", edgecolor="none", linewidth=0.0, zorder=2))
    draw_unit_cells(ax, bulk_points_list, facecolor="#dbeafe", edgecolor="#3b82f6", alpha=0.55, zorder=3, label="_nolegend_")

    add_hole_patches(
        ax,
        original_holes,
        facecolor="#fca5a5",
        edgecolor="#ef4444",
        linewidth=0.25,
        alpha=0.22,
        zorder=4,
        label="uncut holes",
    )
    add_hole_patches(
        ax,
        holes,
        facecolor="#991b1b",
        edgecolor="#7f1d1d",
        linewidth=0.35,
        alpha=0.88,
        zorder=5,
        label="clipped holes",
    )

    inner_closed = closed_line(inner_boundary)
    expanded_closed = closed_line(expanded_boundary)
    ax.plot(inner_closed[:, 0], inner_closed[:, 1], color="#111827", linewidth=0.75, label="bulk cell-union edge", zorder=6)
    ax.plot(expanded_closed[:, 0], expanded_closed[:, 1], color="#a855f7", linewidth=0.9, label="miter edge + gap", zorder=7)

    apply_limits(ax, limits)
    ax.set_title(f"gap = {metadata['geometric_gap']:g} um")
    info = (
        f"a = {A:g} um\n"
        f"bulk = {BULK_PARAMS['name']}\n"
        f"outside = {CLADDING_PARAMS['name']}\n"
        f"gap = {metadata['geometric_gap']:g} um\n"
        f"bulk cells = {metadata['bulk_cells']}\n"
        f"outside cells = {metadata['cladding_candidate_cells']}\n"
        f"hole polygons = {metadata['total_hole_polygons']}\n"
        f"clipped outside holes = {metadata['clipped_cladding_holes']}"
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
 

def draw_annotated_geometry(path: Path, holes: list[np.ndarray], metadata: dict[str, object]) -> None:
    limits = limits_from_polygons([*metadata["original_holes"], *holes, metadata["expanded_boundary"]])
    fig, ax = plt.subplots(figsize=(8.6, 8.6), constrained_layout=True)
    draw_annotated_panel(ax, holes, metadata, limits)
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.suptitle("annotated clipped cavity hole polygons", fontsize=14)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def gap_stem(gap: float) -> str:
    return f"gap_{gap:g}".replace(".", "p")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.glob("*.png"):
        path.unlink()

    written_paths = []
    for gap in GEOMETRIC_GAPS:
        holes, metadata = build_simulation_holes(gap)
        stem = gap_stem(gap)
        geometry_path = OUT_DIR / f"{stem}_geometry.png"
        annotated_path = OUT_DIR / f"{stem}_geometry_annotated.png"
        draw_geometry(geometry_path, holes, metadata)
        draw_annotated_geometry(annotated_path, holes, metadata)
        written_paths.extend([geometry_path, annotated_path])

    print(f"Wrote {len(written_paths)} gap-specific figures to {OUT_DIR}")
    for path in written_paths:
        print(f"  {path.name}")


if __name__ == "__main__":
    main()
