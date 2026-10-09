from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .polygon_utils import (
    exterior_boundary_segments,
    ordered_boundary_loop,
    polygon_fully_inside_region,
    polygons_intersect,
)


@dataclass(frozen=True)
class LatticePoint:
    i: int
    j: int
    shell: int
    x: float
    y: float


def uv_to_xy(u, v) -> tuple[np.ndarray, np.ndarray]:
    u_arr = np.asarray(u, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    return u_arr + 0.5 * v_arr, np.sqrt(3.0) * v_arr / 2.0


def shell_value(u, v) -> np.ndarray:
    u_arr = np.asarray(u, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    return np.maximum.reduce([np.abs(u_arr), np.abs(v_arr), np.abs(u_arr + v_arr)])


def index_shell(i, j) -> np.ndarray:
    i_arr = np.asarray(i, dtype=int)
    j_arr = np.asarray(j, dtype=int)
    return np.maximum.reduce([np.abs(i_arr), np.abs(j_arr), np.abs(i_arr + j_arr)])


def normal_gap_to_shell(gap: float) -> float:
    return 2.0 * gap / np.sqrt(3.0)


def regular_hex_outline(shell: float) -> np.ndarray:
    u = np.array([shell, 0.0, -shell, -shell, 0.0, shell], dtype=float)
    v = np.array([0.0, shell, shell, 0.0, -shell, -shell], dtype=float)
    x, y = uv_to_xy(u, v)
    return np.column_stack([x, y])


def unit_cell_corners(period: float) -> np.ndarray:
    radius = period / np.sqrt(3.0)
    angles = np.pi / 6.0 + np.arange(6, dtype=float) * np.pi / 3.0
    return np.column_stack([radius * np.cos(angles), radius * np.sin(angles)])


def cell_polygon(center_x: float, center_y: float, period: float) -> np.ndarray:
    return unit_cell_corners(period) + np.array([center_x, center_y], dtype=float)


def lattice_points(max_shell: int, period: float, min_shell: int = 0) -> list[LatticePoint]:
    points = []
    for i in range(-max_shell, max_shell + 1):
        for j in range(-max_shell, max_shell + 1):
            shell = int(index_shell(i, j))
            if shell < min_shell or shell > max_shell:
                continue
            x, y = uv_to_xy(i * period, j * period)
            points.append(LatticePoint(i=i, j=j, shell=shell, x=float(x), y=float(y)))
    return points


def shell_indices(shell: int) -> list[tuple[int, int]]:
    if shell == 0:
        return [(0, 0)]
    indices = []
    for i in range(-shell, shell + 1):
        for j in range(-shell, shell + 1):
            if int(index_shell(i, j)) == shell:
                indices.append((i, j))
    return indices


def boundary_from_lattice_cells(points: list[LatticePoint], period: float) -> np.ndarray:
    cells = [cell_polygon(point.x, point.y, period) for point in points]
    return ordered_boundary_loop(exterior_boundary_segments(cells))


def bulk_points(radius: int, period: float) -> list[LatticePoint]:
    return lattice_points(radius, period)


def bulk_boundary(radius: int, period: float) -> np.ndarray:
    return boundary_from_lattice_cells(bulk_points(radius, period), period)


def shell_cells_all_intersect_boundary_region(period: float, shell: int, boundary: np.ndarray) -> bool:
    for i, j in shell_indices(shell):
        x, y = uv_to_xy(i * period, j * period)
        if not polygons_intersect(cell_polygon(float(x), float(y), period), boundary):
            return False
    return True


def first_cladding_shell(period: float, boundary: np.ndarray, layer_count: int = 0, extra_shells: int = 6) -> int:
    max_radius = float(np.max(np.linalg.norm(boundary, axis=1)))
    max_shell = int(np.ceil(max_radius / period)) + layer_count + extra_shells
    first_shell = None
    for shell in range(max_shell + 1):
        if shell_cells_all_intersect_boundary_region(period, shell, boundary):
            first_shell = shell
    if first_shell is None:
        raise RuntimeError("Could not find a cladding shell intersecting the boundary.")
    return first_shell


def retained_cladding_points(period: float, boundary: np.ndarray, layer_count: int) -> tuple[int, list[LatticePoint]]:
    first_shell = first_cladding_shell(period, boundary, layer_count=layer_count)
    last_shell = first_shell + layer_count - 1
    candidates = lattice_points(last_shell, period)
    retained = []
    for point in candidates:
        polygon = cell_polygon(point.x, point.y, period)
        if not polygon_fully_inside_region(polygon, boundary):
            retained.append(point)
    return first_shell, retained


def boundary_for_extra_layers(radius: int, period: float, extra_layers: int) -> np.ndarray:
    return boundary_from_lattice_cells(lattice_points(radius + extra_layers, period), period)


def ray_boundary_intersection(direction: np.ndarray, boundary: np.ndarray) -> np.ndarray:
    best_distance = -np.inf
    best_point = None
    for idx, start in enumerate(boundary):
        end = boundary[(idx + 1) % len(boundary)]
        edge = end - start
        matrix = np.column_stack([direction, -edge])
        det = float(np.linalg.det(matrix))
        if abs(det) < 1e-12:
            continue
        ray_distance, edge_fraction = np.linalg.solve(matrix, start)
        if ray_distance >= -1e-10 and -1e-10 <= edge_fraction <= 1.0 + 1e-10 and ray_distance > best_distance:
            best_distance = float(ray_distance)
            best_point = ray_distance * direction
    if best_point is None:
        return np.zeros(2, dtype=float)
    return np.asarray(best_point, dtype=float)


def lattice_layer_continuous_shell(radius: int, period: float, gap: float) -> np.ndarray:
    layer_position = gap / period
    inner_layer = int(np.floor(layer_position + 1e-12))
    fraction = float(layer_position - inner_layer)
    if fraction < 1e-12:
        return boundary_for_extra_layers(radius, period, inner_layer)

    inner_boundary = boundary_for_extra_layers(radius, period, inner_layer)
    outer_boundary = boundary_for_extra_layers(radius, period, inner_layer + 1)
    shell = []
    for outer_point in outer_boundary:
        norm = float(np.linalg.norm(outer_point))
        if norm < 1e-12:
            shell.append(outer_point)
            continue
        direction = outer_point / norm
        inner_point = ray_boundary_intersection(direction, inner_boundary)
        shell.append(inner_point + fraction * (outer_point - inner_point))
    return np.asarray(shell, dtype=float)
