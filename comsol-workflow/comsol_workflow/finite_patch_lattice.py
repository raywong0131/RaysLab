from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd


BasisName = Literal["orthogonal", "hex"]
COLUMNS = ["region", "part", "layer", "i", "j", "u", "v", "x", "y"]


@dataclass(frozen=True)
class BulkSideCornerConfig:
    a: float
    b: float
    gap: float
    bulk_periods: int
    boundary_periods: int

    def __post_init__(self) -> None:
        if self.a <= 0.0:
            raise ValueError("a must be positive")
        if self.b <= 0.0:
            raise ValueError("b must be positive")
        if self.gap < 0.0:
            raise ValueError("gap must be non-negative")
        if self.bulk_periods <= 0:
            raise ValueError("bulk_periods must be positive")
        if self.boundary_periods <= 0:
            raise ValueError("boundary_periods must be positive")


def primitive_vectors(basis: BasisName = "orthogonal") -> tuple[np.ndarray, np.ndarray]:
    if basis == "orthogonal":
        return np.array([1.0, 0.0]), np.array([0.0, 1.0])
    if basis == "hex":
        return np.array([1.0, 0.0]), np.array([0.5, np.sqrt(3.0) / 2.0])
    raise ValueError("basis must be 'orthogonal' or 'hex'")


def uv_to_xy(u, v, basis: BasisName = "orthogonal") -> tuple[np.ndarray, np.ndarray]:
    u_arr = np.asarray(u, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    e_u, e_v = primitive_vectors(basis)
    x = u_arr * e_u[0] + v_arr * e_v[0]
    y = u_arr * e_u[1] + v_arr * e_v[1]
    return x, y


def centered_axis(periods: int, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    if periods <= 0:
        raise ValueError("periods must be positive")
    indices = np.arange(periods, dtype=float) - (periods - 1) / 2.0
    return indices, indices * spacing


def boundary_offsets(config: BulkSideCornerConfig) -> np.ndarray:
    half_span = (config.bulk_periods - 1) * config.a / 2.0
    layers = np.arange(config.boundary_periods, dtype=float)
    return half_span + config.gap + layers * config.b


def unit_cell_edges(config: BulkSideCornerConfig) -> tuple[float, float]:
    inner_edge = config.bulk_periods * config.a / 2.0
    outer_edge = float(boundary_offsets(config)[-1] + config.b / 2.0)
    return inner_edge, outer_edge


def part_intervals(config: BulkSideCornerConfig) -> list[tuple[str, str, float, float, float, float]]:
    half_span = (config.bulk_periods - 1) * config.a / 2.0
    inner = half_span + config.gap
    outer = inner + (config.boundary_periods - 1) * config.b
    return [
        ("bulk", "bulk", -half_span, half_span, -half_span, half_span),
        ("side", "plus_v", -half_span, half_span, inner, outer),
        ("side", "minus_v", -half_span, half_span, -outer, -inner),
        ("side", "plus_u", inner, outer, -half_span, half_span),
        ("side", "minus_u", -outer, -inner, -half_span, half_span),
        ("corner", "corner_plus_u_plus_v", inner, outer, inner, outer),
        ("corner", "corner_minus_u_plus_v", -outer, -inner, inner, outer),
        ("corner", "corner_plus_u_minus_v", inner, outer, -outer, -inner),
        ("corner", "corner_minus_u_minus_v", -outer, -inner, -outer, -inner),
    ]


def oblique_box_corners(
    u0: float,
    u1: float,
    v0: float,
    v1: float,
    basis: BasisName = "orthogonal",
) -> np.ndarray:
    u = np.array([u0, u1, u1, u0], dtype=float)
    v = np.array([v0, v0, v1, v1], dtype=float)
    x, y = uv_to_xy(u, v, basis)
    return np.column_stack([x, y])


def _frame(
    region: str,
    part: str,
    layer,
    i,
    j,
    u,
    v,
    basis: BasisName,
) -> pd.DataFrame:
    x, y = uv_to_xy(u, v, basis)
    return pd.DataFrame({
        "region": region,
        "part": part,
        "layer": np.asarray(layer, dtype=int).ravel(),
        "i": np.asarray(i, dtype=float).ravel(),
        "j": np.asarray(j, dtype=float).ravel(),
        "u": np.asarray(u, dtype=float).ravel(),
        "v": np.asarray(v, dtype=float).ravel(),
        "x": np.asarray(x, dtype=float).ravel(),
        "y": np.asarray(y, dtype=float).ravel(),
    }, columns=COLUMNS)


def bulk_points(config: BulkSideCornerConfig, basis: BasisName = "orthogonal") -> pd.DataFrame:
    axis, coords = centered_axis(config.bulk_periods, config.a)
    i, j = np.meshgrid(axis, axis, indexing="xy")
    u, v = np.meshgrid(coords, coords, indexing="xy")
    layer = np.zeros(u.size, dtype=int)
    return _frame("bulk", "bulk", layer, i, j, u, v, basis)


def side_strip_points(
    config: BulkSideCornerConfig,
    part: str,
    sign: float,
    fixed_axis: Literal["u", "v"],
    basis: BasisName = "orthogonal",
) -> pd.DataFrame:
    axis, coords = centered_axis(config.bulk_periods, config.a)
    offsets = boundary_offsets(config)
    layers = np.arange(config.boundary_periods, dtype=int)

    layer_grid = np.tile(layers[:, None], (1, config.bulk_periods))
    axis_grid = np.tile(axis, (config.boundary_periods, 1))
    coord_grid = np.tile(coords, (config.boundary_periods, 1))
    offset_grid = np.tile(offsets[:, None], (1, config.bulk_periods))

    if fixed_axis == "u":
        i = sign * layer_grid
        j = axis_grid
        u = sign * offset_grid
        v = coord_grid
    elif fixed_axis == "v":
        i = axis_grid
        j = sign * layer_grid
        u = coord_grid
        v = sign * offset_grid
    else:
        raise ValueError("fixed_axis must be 'u' or 'v'")

    return _frame("side", part, layer_grid, i, j, u, v, basis)


def side_strips(config: BulkSideCornerConfig, basis: BasisName = "orthogonal") -> pd.DataFrame:
    return pd.concat([
        side_strip_points(config, "plus_v", 1.0, "v", basis),
        side_strip_points(config, "minus_v", -1.0, "v", basis),
        side_strip_points(config, "plus_u", 1.0, "u", basis),
        side_strip_points(config, "minus_u", -1.0, "u", basis),
    ], ignore_index=True)


def corner_block_points(
    config: BulkSideCornerConfig,
    part: str,
    sign_u: float,
    sign_v: float,
    basis: BasisName = "orthogonal",
) -> pd.DataFrame:
    offsets = boundary_offsets(config)
    layers = np.arange(config.boundary_periods, dtype=int)
    u_layer, v_layer = np.meshgrid(layers, layers, indexing="ij")
    u_offset, v_offset = np.meshgrid(offsets, offsets, indexing="ij")
    layer = np.maximum(u_layer, v_layer)
    return _frame(
        "corner",
        part,
        layer,
        sign_u * u_layer,
        sign_v * v_layer,
        sign_u * u_offset,
        sign_v * v_offset,
        basis,
    )


def corner_blocks(config: BulkSideCornerConfig, basis: BasisName = "orthogonal") -> pd.DataFrame:
    return pd.concat([
        corner_block_points(config, "corner_plus_u_plus_v", 1.0, 1.0, basis),
        corner_block_points(config, "corner_minus_u_plus_v", -1.0, 1.0, basis),
        corner_block_points(config, "corner_plus_u_minus_v", 1.0, -1.0, basis),
        corner_block_points(config, "corner_minus_u_minus_v", -1.0, -1.0, basis),
    ], ignore_index=True)


def generate_lattice(config: BulkSideCornerConfig, basis: BasisName = "orthogonal") -> pd.DataFrame:
    return pd.concat(
        [bulk_points(config, basis), side_strips(config, basis), corner_blocks(config, basis)],
        ignore_index=True,
    )


def unit_cell_spacing_for_part(part: str, config: BulkSideCornerConfig) -> tuple[float, float]:
    if part == "bulk":
        return config.a, config.a
    if part in {"plus_v", "minus_v"}:
        return config.a, config.b
    if part in {"plus_u", "minus_u"}:
        return config.b, config.a
    if part.startswith("corner_"):
        return config.b, config.b
    raise ValueError(f"unknown lattice part: {part}")


def _dedupe_polygon_vertices(vertices: np.ndarray, tolerance: float = 1e-10) -> np.ndarray:
    if len(vertices) == 0:
        return vertices
    kept = [vertices[0]]
    for vertex in vertices[1:]:
        if np.linalg.norm(vertex - kept[-1]) > tolerance:
            kept.append(vertex)
    if len(kept) > 1 and np.linalg.norm(kept[0] - kept[-1]) <= tolerance:
        kept.pop()
    return np.asarray(kept, dtype=float)


def _clip_polygon_with_halfplane(polygon: np.ndarray, normal: np.ndarray, offset: float) -> np.ndarray:
    tolerance = 1e-12
    if len(polygon) == 0:
        return polygon
    clipped: list[np.ndarray] = []
    previous = polygon[-1]
    previous_inside = float(previous @ normal) <= offset + tolerance
    for current in polygon:
        current_inside = float(current @ normal) <= offset + tolerance
        if current_inside != previous_inside:
            direction = current - previous
            denominator = float(direction @ normal)
            if abs(denominator) > tolerance:
                fraction = (offset - float(previous @ normal)) / denominator
                clipped.append(previous + fraction * direction)
        if current_inside:
            clipped.append(current)
        previous = current
        previous_inside = current_inside
    return _dedupe_polygon_vertices(np.asarray(clipped, dtype=float))


def lattice_unit_cell_corners(
    spacing_u: float,
    spacing_v: float,
    basis: BasisName = "orthogonal",
) -> np.ndarray:
    if spacing_u <= 0.0 or spacing_v <= 0.0:
        raise ValueError("spacing_u and spacing_v must be positive")

    e_u, e_v = primitive_vectors(basis)
    lattice_vectors = [spacing_u * e_u, spacing_v * e_v]
    size = 4.0 * max(spacing_u, spacing_v)
    polygon = np.array([
        [-size, -size],
        [size, -size],
        [size, size],
        [-size, size],
    ], dtype=float)

    search_radius = max(3, int(np.ceil(2.0 * max(spacing_u, spacing_v) / min(spacing_u, spacing_v))) + 2)
    for i in range(-search_radius, search_radius + 1):
        for j in range(-search_radius, search_radius + 1):
            if i == 0 and j == 0:
                continue
            vector = i * lattice_vectors[0] + j * lattice_vectors[1]
            norm_sq = float(vector @ vector)
            if norm_sq <= 0.0:
                continue
            polygon = _clip_polygon_with_halfplane(polygon, vector, norm_sq / 2.0)

    center = polygon.mean(axis=0)
    angles = np.arctan2(polygon[:, 1] - center[1], polygon[:, 0] - center[0])
    return polygon[np.argsort(angles)]
