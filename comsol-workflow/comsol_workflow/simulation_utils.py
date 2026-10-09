from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable, Literal, Sequence

import numpy as np
from scipy.optimize import linprog

import mph
from jpype import JClass
from jpype.types import JArray, JDouble, JInt, JString

from .energy_recovery import integrate_field, interpolate_field


SimulationMode = Literal["unit_cell", "finite_size", "finite_quarter", "strip_1d"]
BoundaryShape = Literal[
    "hexagon",
    "hexagon_boundary",
    "quarter_hexagon_boundary",
    "quarter_rectangle",
    "square",
    "rectangle",
]


def _isotropic_tensor(value: str) -> list[str]:
    return [value, "0", "0", "0", value, "0", "0", "0", value]


def _polygon_signed_area(vertices: np.ndarray) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _turn_cross_z(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    ab = b - a
    bc = c - b
    return float(ab[0] * bc[1] - ab[1] * bc[0])


def _polygon_inward_halfspaces(
    vertices: np.ndarray,
    idx: int,
    tolerance: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    normals = []
    offsets = []
    for vertex_idx, start in enumerate(vertices):
        end = vertices[(vertex_idx + 1) % vertices.shape[0]]
        edge = end - start
        length = float(np.linalg.norm(edge))
        if length <= tolerance:
            raise ValueError(
                f"LayerSpec.holes[{idx}] must not contain zero-length edges; "
                f"edge starts at vertex {vertex_idx}"
            )
        normal = np.array([-edge[1], edge[0]], dtype=float) / length
        normals.append(normal)
        offsets.append(float(np.dot(normal, start)))
    return np.asarray(normals, dtype=float), np.asarray(offsets, dtype=float)


def _point_min_boundary_distance(point: np.ndarray, normals: np.ndarray, offsets: np.ndarray) -> float:
    return float(np.min(normals @ point - offsets))


def _max_inscribed_ball(
    normals: np.ndarray,
    offsets: np.ndarray,
    idx: int,
    tolerance: float = 1e-12,
) -> tuple[np.ndarray, float]:
    # Maximize r subject to each edge half-space containing a radius-r ball.
    a_ub = np.column_stack([-normals, np.ones(normals.shape[0])])
    b_ub = -offsets
    result = linprog(
        c=np.array([0.0, 0.0, -1.0]),
        A_ub=a_ub,
        b_ub=b_ub,
        bounds=[(None, None), (None, None), (0.0, None)],
        method="highs",
    )
    if not result.success:
        raise ValueError(f"LayerSpec.holes[{idx}] failed to compute an interior selection seed: {result.message}")

    center = np.asarray(result.x[:2], dtype=float)
    radius = min(float(result.x[2]), _point_min_boundary_distance(center, normals, offsets))
    if radius <= tolerance:
        raise ValueError(
            f"LayerSpec.holes[{idx}] has no positive-radius interior ball for COMSOL selection"
        )
    return center, radius


def _hole_selection_seed(
    vertices: np.ndarray,
    idx: int,
    tolerance: float = 1e-12,
) -> tuple[np.ndarray, float]:
    normals, offsets = _polygon_inward_halfspaces(vertices, idx, tolerance=tolerance)
    return _max_inscribed_ball(normals, offsets, idx, tolerance=tolerance)


def _validate_hole_polygon(vertices: np.ndarray, idx: int, tolerance: float = 1e-12) -> None:
    if not np.all(np.isfinite(vertices)):
        raise ValueError(f"LayerSpec.holes[{idx}] must contain only finite coordinates")

    area = _polygon_signed_area(vertices)
    if area <= tolerance:
        raise ValueError(
            f"LayerSpec.holes[{idx}] must be counter-clockwise and non-degenerate; "
            f"signed area is {area:g}"
        )

    _polygon_inward_halfspaces(vertices, idx, tolerance=tolerance)

    for vertex_idx in range(vertices.shape[0]):
        cross = _turn_cross_z(
            vertices[vertex_idx],
            vertices[(vertex_idx + 1) % vertices.shape[0]],
            vertices[(vertex_idx + 2) % vertices.shape[0]],
        )
        if cross < -tolerance:
            raise ValueError(
                f"LayerSpec.holes[{idx}] must be convex; "
                f"found clockwise turn near vertex {vertex_idx + 1} with cross product {cross:g}"
            )


@dataclass(frozen=True)
class LayerSpec:
    height: str
    holes: list[np.ndarray]
    refractive_index: str
    extinction_coefficient: str = "0"
    label: str | None = None
    hole_selection_centers: list[np.ndarray] = field(init=False, repr=False)
    hole_selection_max_radii: list[float] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.height, str) or not self.height.strip():
            raise ValueError("LayerSpec.height must be a non-empty COMSOL expression string")
        if not isinstance(self.refractive_index, str) or not self.refractive_index.strip():
            raise ValueError("LayerSpec.refractive_index must be a non-empty string")
        if not isinstance(self.extinction_coefficient, str) or not self.extinction_coefficient.strip():
            raise ValueError("LayerSpec.extinction_coefficient must be a non-empty string")
        if self.label is not None and (not isinstance(self.label, str) or not self.label.strip()):
            raise ValueError("LayerSpec.label must be None or a non-empty string")

        if self.holes is None or isinstance(self.holes, (str, bytes)):
            raise TypeError("LayerSpec.holes must be a sequence of 2D polygon arrays")

        normalized_holes = []
        hole_selection_centers = []
        hole_selection_max_radii = []
        for idx, hole in enumerate(self.holes):
            array = np.asarray(hole, dtype=float)
            if array.ndim != 2 or array.shape[1] != 2:
                raise ValueError(f"LayerSpec.holes[{idx}] must have shape (N, 2)")
            if array.shape[0] < 3:
                raise ValueError(f"LayerSpec.holes[{idx}] must have at least 3 vertices")
            _validate_hole_polygon(array, idx)
            center, max_radius = _hole_selection_seed(array, idx)
            normalized_holes.append(array)
            hole_selection_centers.append(center)
            hole_selection_max_radii.append(max_radius)

        object.__setattr__(self, "holes", normalized_holes)
        object.__setattr__(self, "hole_selection_centers", hole_selection_centers)
        object.__setattr__(self, "hole_selection_max_radii", hole_selection_max_radii)


@dataclass(frozen=True)
class BoundarySpec:
    shape: BoundaryShape
    a: float
    b: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.shape, str):
            raise TypeError("BoundarySpec.shape must be a string")

        shape = self.shape.lower()
        if shape not in {
            "hexagon",
            "hexagon_boundary",
            "quarter_hexagon_boundary",
            "quarter_rectangle",
            "square",
            "rectangle",
        }:
            raise ValueError(
                "BoundarySpec.shape must be 'hexagon', 'hexagon_boundary', "
                "'quarter_hexagon_boundary', 'quarter_rectangle', 'square', "
                "or 'rectangle'"
            )
        object.__setattr__(self, "shape", shape)

        a = float(self.a)
        if a <= 0:
            raise ValueError("BoundarySpec.a must be positive")
        object.__setattr__(self, "a", a)

        if shape in {"hexagon", "hexagon_boundary", "quarter_hexagon_boundary", "square"}:
            if self.b is not None:
                raise ValueError("BoundarySpec.b must be None for hexagon, hexagon_boundary, and square boundaries")
            return

        if self.b is None:
            raise ValueError("BoundarySpec.b is required for rectangle boundaries")
        b = float(self.b)
        if b <= 0:
            raise ValueError("BoundarySpec.b must be positive")
        object.__setattr__(self, "b", b)


@dataclass(frozen=True)
class _LayerStackEntry:
    index: int
    height_param: str
    model_height_expr: str
    z_lo_param: str
    z_hi_param: str
    z_mid_param: str
    z_lo_expr: str
    z_hi_expr: str
    z_mid_expr: str


@dataclass(frozen=True)
class _LayerStack:
    entries: list[_LayerStackEntry]
    z_slab_top_param: str
    z_slab_top_expr: str


@dataclass(frozen=True)
class _BoundaryDefinition:
    footprint_label: str
    footprint_table: list[list[str]]
    ga_expr: str
    gb_expr: str | None = None


@dataclass(frozen=True)
class SimulationConfig:
    wavelength: str = "1550 [nm]"
    air_height: str = "lda0"
    air_layer_top_factor: float = 1.0
    pml_layer_top_factor: float = 1.5
    air_refractive_index: str = "1"
    air_extinction_coefficient: str = "0"
    eigenfrequency_shift: str = "c_const/1.55[um]"
    eigenmode_count: int = 8
    air_cutplane_z: str = "0.75*H_air"
    selection_tolerance: str = "1 [nm]"
    mode_type: str = "TE"
    mesh_auto_size: int = 9
    symmetry_x_boundary: str | None = None
    symmetry_y_boundary: str | None = None

    @staticmethod
    def _layer_distance_expr(factor: float) -> str:
        if factor == 1.0:
            return "H_air"
        return f"H_air*{factor:g}"

    @property
    def air_layer_distance(self) -> str:
        return self._layer_distance_expr(self.air_layer_top_factor)

    @property
    def pml_layer_distance(self) -> str:
        return self._layer_distance_expr(self.pml_layer_top_factor)


def _sum_expr(expressions: Sequence[str]) -> str:
    if not expressions:
        return "0 [um]"
    return " + ".join(expressions)


def _make_layer_stack(layers: Sequence[LayerSpec]) -> _LayerStack:
    model_heights = []
    entries = []
    for idx, _layer in enumerate(layers):
        height_param = f"H_layer_{idx}"
        model_height_expr = f"{height_param}/2" if idx == 0 else height_param
        z_lo_expr = _sum_expr(model_heights)
        z_hi_expr = _sum_expr([*model_heights, model_height_expr])
        z_mid_expr = f"({z_lo_expr} + {z_hi_expr})/2"
        entries.append(_LayerStackEntry(
            index=idx,
            height_param=height_param,
            model_height_expr=model_height_expr,
            z_lo_param=f"z_layer_{idx}_lo",
            z_hi_param=f"z_layer_{idx}_hi",
            z_mid_param=f"z_layer_{idx}_mid",
            z_lo_expr=z_lo_expr,
            z_hi_expr=z_hi_expr,
            z_mid_expr=z_mid_expr,
        ))
        model_heights.append(model_height_expr)
    return _LayerStack(entries=entries, z_slab_top_param="z_slab_top", z_slab_top_expr=_sum_expr(model_heights))


def _normalize_layers(layers: Sequence[LayerSpec]) -> list[LayerSpec]:
    if layers is None or isinstance(layers, (str, bytes)):
        raise TypeError("layers must be a non-empty sequence of LayerSpec objects")

    normalized = list(layers)
    if not normalized:
        raise ValueError("layers must contain at least one LayerSpec")
    for idx, layer in enumerate(normalized):
        if not isinstance(layer, LayerSpec):
            raise TypeError(
                f"layers[{idx}] must be a LayerSpec. "
                "The old build_and_run(a, holes, k=None) API is no longer supported."
            )
    return normalized


def _normalize_simulation_mode(simulation_mode: str) -> SimulationMode:
    if not isinstance(simulation_mode, str):
        raise TypeError(
            "simulation_mode must be 'unit_cell', 'finite_size', "
            "'finite_quarter', or 'strip_1d'"
        )
    mode = simulation_mode.lower()
    if mode not in {"unit_cell", "finite_size", "finite_quarter", "strip_1d"}:
        raise ValueError(
            "simulation_mode must be 'unit_cell', 'finite_size', "
            "'finite_quarter', or 'strip_1d'"
        )
    return mode


def _normalize_k(k: dict[str, float] | None, simulation_mode: SimulationMode) -> tuple[dict[str, float], bool]:
    if simulation_mode in {"finite_size", "finite_quarter"}:
        if k is not None:
            raise ValueError("finite_size simulations do not accept k/Floquet input")
        return {"kx": 0.0, "ky": 0.0}, False

    if simulation_mode == "strip_1d":
        if k is None:
            return {"kx": 0.0, "ky": 0.0}, False
        if not isinstance(k, dict) or set(k) != {"kx"}:
            raise TypeError("strip_1d k must be None or a dict with only 'kx'")
        return {"kx": float(k["kx"]), "ky": 0.0}, True

    if k is None:
        return {"kx": 0.0, "ky": 0.0}, False
    if not isinstance(k, dict) or "kx" not in k or "ky" not in k:
        raise TypeError("k must be None or a dict with 'kx' and 'ky'")
    return {"kx": float(k["kx"]), "ky": float(k["ky"])}, True


def _validate_build_inputs(
    boundary: BoundarySpec,
    layers: Sequence[LayerSpec],
    simulation_mode: str,
    k: dict[str, float] | None,
) -> tuple[BoundarySpec, list[LayerSpec], SimulationMode, dict[str, float], bool]:
    if not isinstance(boundary, BoundarySpec):
        raise TypeError(
            "boundary must be a BoundarySpec. "
            "The old build_and_run(a, holes, k=None) API is no longer supported."
        )
    normalized_layers = _normalize_layers(layers)
    normalized_mode = _normalize_simulation_mode(simulation_mode)
    if normalized_mode == "strip_1d" and boundary.shape != "rectangle":
        raise ValueError("strip_1d simulations require a rectangle BoundarySpec")
    if normalized_mode == "finite_quarter" and boundary.shape not in {
        "quarter_hexagon_boundary",
        "quarter_rectangle",
    }:
        raise ValueError(
            "finite_quarter simulations require a quarter_hexagon_boundary or "
            "quarter_rectangle BoundarySpec"
        )
    k_values, use_floquet = _normalize_k(k, normalized_mode)
    return boundary, normalized_layers, normalized_mode, k_values, use_floquet


def _boundary_definition(boundary: BoundarySpec) -> _BoundaryDefinition:
    if boundary.shape == "hexagon":
        return _BoundaryDefinition(
            footprint_label="Hex",
            footprint_table=[
                ["a/sqrt(3)*cos(pi/6)", "a/sqrt(3)*sin(pi/6)"],
                ["0", "a/sqrt(3)"],
                ["-a/sqrt(3)*cos(pi/6)", "a/sqrt(3)*sin(pi/6)"],
                ["-a/sqrt(3)*cos(pi/6)", "-a/sqrt(3)*sin(pi/6)"],
                ["0", "-a/sqrt(3)"],
                ["a/sqrt(3)*cos(pi/6)", "-a/sqrt(3)*sin(pi/6)"],
            ],
            ga_expr="4*pi/sqrt(3)/a",
        )
    if boundary.shape == "hexagon_boundary":
        return _BoundaryDefinition(
            footprint_label="HexagonBoundary",
            footprint_table=[
                ["a/2", "0"],
                ["a/4", "sqrt(3)*a/4"],
                ["-a/4", "sqrt(3)*a/4"],
                ["-a/2", "0"],
                ["-a/4", "-sqrt(3)*a/4"],
                ["a/4", "-sqrt(3)*a/4"],
            ],
            ga_expr="2*pi/a",
        )
    if boundary.shape == "quarter_hexagon_boundary":
        return _BoundaryDefinition(
            footprint_label="QuarterHexagonBoundary",
            footprint_table=[
                ["0", "0"],
                ["a/2", "0"],
                ["a/4", "sqrt(3)*a/4"],
                ["0", "sqrt(3)*a/4"],
            ],
            ga_expr="2*pi/a",
        )
    if boundary.shape == "quarter_rectangle":
        return _BoundaryDefinition(
            footprint_label="QuarterRectangle",
            footprint_table=[
                ["0", "0"],
                ["a/2", "0"],
                ["a/2", "b/2"],
                ["0", "b/2"],
            ],
            ga_expr="2*pi/a",
            gb_expr="2*pi/b",
        )
    if boundary.shape == "square":
        return _BoundaryDefinition(
            footprint_label="Square",
            footprint_table=[
                ["a/2", "a/2"],
                ["-a/2", "a/2"],
                ["-a/2", "-a/2"],
                ["a/2", "-a/2"],
            ],
            ga_expr="2*pi/a",
        )
    return _BoundaryDefinition(
        footprint_label="Rectangle",
        footprint_table=[
            ["a/2", "b/2"],
            ["-a/2", "b/2"],
            ["-a/2", "-b/2"],
            ["a/2", "-b/2"],
        ],
        ga_expr="2*pi/a",
        gb_expr="2*pi/b",
    )


def _wall_seed_points(boundary: BoundarySpec, z: float) -> dict[str, list[tuple[float, float, float]]]:
    if boundary.shape == "hexagon":
        r_in = boundary.a / 2.0
        return {
            "pc1": [
                (r_in, 0.0, z),
                (-r_in, 0.0, z),
            ],
            "pc2": [
                (r_in / 2.0, boundary.a * float(np.sqrt(3)) / 4.0, z),
                (-r_in / 2.0, -boundary.a * float(np.sqrt(3)) / 4.0, z),
            ],
            "pc3": [
                (-r_in / 2.0, boundary.a * float(np.sqrt(3)) / 4.0, z),
                (r_in / 2.0, -boundary.a * float(np.sqrt(3)) / 4.0, z),
            ],
        }
    if boundary.shape == "square":
        half_width = boundary.a / 2.0
        return {
            "pc1": [
                (half_width, 0.0, z),
                (-half_width, 0.0, z),
            ],
            "pc2": [
                (0.0, half_width, z),
                (0.0, -half_width, z),
            ],
        }
    if boundary.shape == "hexagon_boundary":
        root3 = float(np.sqrt(3.0))
        return {
            "side_hex_ne_sw": [
                (3.0 * boundary.a / 8.0, root3 * boundary.a / 8.0, z),
                (-3.0 * boundary.a / 8.0, -root3 * boundary.a / 8.0, z),
            ],
            "side_hex_n_s": [
                (0.0, root3 * boundary.a / 4.0, z),
                (0.0, -root3 * boundary.a / 4.0, z),
            ],
            "side_hex_nw_se": [
                (-3.0 * boundary.a / 8.0, root3 * boundary.a / 8.0, z),
                (3.0 * boundary.a / 8.0, -root3 * boundary.a / 8.0, z),
            ],
        }
    if boundary.shape == "quarter_hexagon_boundary":
        root3 = float(np.sqrt(3.0))
        return {
            "symmetry_x_axis": [
                (boundary.a / 4.0, 0.0, z),
            ],
            "symmetry_y_axis": [
                (0.0, root3 * boundary.a / 8.0, z),
            ],
            "outer_hex_ne": [
                (3.0 * boundary.a / 8.0, root3 * boundary.a / 8.0, z),
            ],
            "outer_hex_n": [
                (boundary.a / 8.0, root3 * boundary.a / 4.0, z),
            ],
        }
    if boundary.shape == "quarter_rectangle":
        return {
            "symmetry_x_axis": [
                (boundary.a / 4.0, 0.0, z),
            ],
            "symmetry_y_axis": [
                (0.0, float(boundary.b) / 4.0, z),
            ],
            "outer_rect_x": [
                (boundary.a / 2.0, float(boundary.b) / 4.0, z),
            ],
            "outer_rect_y": [
                (boundary.a / 4.0, float(boundary.b) / 2.0, z),
            ],
        }

    half_x = boundary.a / 2.0
    half_y = float(boundary.b) / 2.0
    return {
        "side_x": [
            (half_x, 0.0, z),
            (-half_x, 0.0, z),
        ],
        "side_y": [
            (0.0, half_y, z),
            (0.0, -half_y, z),
        ],
    }


def _boundary_large_extent(boundary: BoundarySpec) -> float:
    if boundary.shape in {"rectangle", "quarter_rectangle"}:
        return max(boundary.a, float(boundary.b)) * 10.0
    return boundary.a * 10.0


def _boundary_face_seed_xy(boundary: BoundarySpec) -> tuple[float, float]:
    if boundary.shape == "quarter_hexagon_boundary":
        return boundary.a / 8.0, float(np.sqrt(3.0)) * boundary.a / 16.0
    if boundary.shape == "quarter_rectangle":
        return boundary.a / 4.0, float(boundary.b) / 4.0
    return 0.0, 0.0


def normalize_symmetry_boundary(condition: str) -> str:
    normalized = str(condition).upper()
    if normalized not in {"PEC", "PMC"}:
        raise ValueError(f"symmetry boundary must be 'PEC' or 'PMC', got {condition!r}")
    return normalized


def hz_reflection_sign(condition: str) -> float:
    """Return the TE-like Hz parity associated with a mirror boundary."""
    return 1.0 if normalize_symmetry_boundary(condition) == "PEC" else -1.0


def electric_reflection_matrix(axis: str, condition: str) -> np.ndarray:
    """Return the polar-vector E transform for reflection across an axis.

    ``axis='x'`` reflects ``(x, y) -> (x, -y)`` and therefore describes the
    boundary lying along the x axis.  ``axis='y'`` reflects
    ``(x, y) -> (-x, y)``.
    """
    normalized = normalize_symmetry_boundary(condition)
    if axis == "x":
        coordinate_reflection = np.diag([1.0, -1.0, 1.0])
    elif axis == "y":
        coordinate_reflection = np.diag([-1.0, 1.0, 1.0])
    else:
        raise ValueError("axis must be 'x' or 'y'")
    symmetry_sign = -1.0 if normalized == "PEC" else 1.0
    return symmetry_sign * coordinate_reflection


def reconstruct_quarter_scalar_field(
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    x_reflection_sign: float,
    y_reflection_sign: float,
    tolerance: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Mirror a first-quadrant scalar field into all four quadrants."""
    coordinates = np.asarray(coordinates, dtype=float)
    values = np.asarray(values)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("coordinates must have shape (N, 2)")
    if values.ndim != 1 or len(values) != len(coordinates):
        raise ValueError("values must have shape (N,)")
    if not np.isfinite(x_reflection_sign) or not np.isfinite(y_reflection_sign):
        raise ValueError("reflection signs must be finite")

    scale = max(float(np.max(np.abs(coordinates), initial=0.0)), 1.0)
    tol = scale * 1e-10 if tolerance is None else float(tolerance)
    if tol < 0.0 or not np.isfinite(tol):
        raise ValueError("tolerance must be finite and non-negative")
    if np.any(coordinates[:, 0] < -tol) or np.any(coordinates[:, 1] < -tol):
        raise ValueError("quarter-field coordinates must lie in x >= 0, y >= 0")

    quarter = coordinates.copy()
    quarter[np.abs(quarter) <= tol] = 0.0
    quarter_values = values.copy()
    if y_reflection_sign < 0.0:
        quarter_values[quarter[:, 0] == 0.0] = 0.0
    if x_reflection_sign < 0.0:
        quarter_values[quarter[:, 1] == 0.0] = 0.0

    coordinate_parts = [quarter]
    value_parts = [quarter_values]
    reflect_y = quarter[:, 0] > tol
    if np.any(reflect_y):
        reflected = quarter[reflect_y].copy()
        reflected[:, 0] *= -1.0
        coordinate_parts.append(reflected)
        value_parts.append(y_reflection_sign * quarter_values[reflect_y])
    reflect_x = quarter[:, 1] > tol
    if np.any(reflect_x):
        reflected = quarter[reflect_x].copy()
        reflected[:, 1] *= -1.0
        coordinate_parts.append(reflected)
        value_parts.append(x_reflection_sign * quarter_values[reflect_x])
    reflect_both = reflect_x & reflect_y
    if np.any(reflect_both):
        reflected = -quarter[reflect_both]
        coordinate_parts.append(reflected)
        value_parts.append(
            x_reflection_sign * y_reflection_sign * quarter_values[reflect_both]
        )
    return np.vstack(coordinate_parts), np.concatenate(value_parts)


def reconstruct_quarter_electric_grid(
    quarter_values: np.ndarray,
    quarter_domain_mask: np.ndarray,
    *,
    x_boundary: str,
    y_boundary: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Mirror a regular first-quadrant E grid into a centered odd full grid."""
    quarter_values = np.asarray(quarter_values, dtype=complex)
    quarter_domain_mask = np.asarray(quarter_domain_mask, dtype=bool)
    if quarter_values.ndim != 3 or quarter_values.shape[0] != 3:
        raise ValueError("quarter_values must have shape (3, N, N)")
    if quarter_values.shape[1] != quarter_values.shape[2]:
        raise ValueError("quarter_values must use a square grid")
    if quarter_domain_mask.shape != quarter_values.shape[1:]:
        raise ValueError("quarter_domain_mask shape must match the field grid")

    half_size = quarter_values.shape[1]
    full_size = 2 * half_size - 1
    midpoint = half_size - 1
    full_values = np.zeros((3, full_size, full_size), dtype=complex)
    full_mask = np.zeros((full_size, full_size), dtype=bool)
    full_values[:, midpoint:, midpoint:] = quarter_values
    full_mask[midpoint:, midpoint:] = quarter_domain_mask

    mirror_x_axis = electric_reflection_matrix("x", x_boundary)
    mirror_y_axis = electric_reflection_matrix("y", y_boundary)

    def transform(matrix: np.ndarray, field: np.ndarray) -> np.ndarray:
        return np.einsum("ab,bij->aij", matrix, field)

    upper_left_source = quarter_values[:, :, 1:][:, :, ::-1]
    full_values[:, midpoint:, :midpoint] = transform(
        mirror_y_axis,
        upper_left_source,
    )
    full_mask[midpoint:, :midpoint] = quarter_domain_mask[:, 1:][:, ::-1]

    lower_right_source = quarter_values[:, 1:, :][:, ::-1, :]
    full_values[:, :midpoint, midpoint:] = transform(
        mirror_x_axis,
        lower_right_source,
    )
    full_mask[:midpoint, midpoint:] = quarter_domain_mask[1:, :][::-1, :]

    lower_left_source = quarter_values[:, 1:, 1:][:, ::-1, ::-1]
    full_values[:, :midpoint, :midpoint] = transform(
        mirror_x_axis @ mirror_y_axis,
        lower_left_source,
    )
    full_mask[:midpoint, :midpoint] = quarter_domain_mask[1:, 1:][::-1, ::-1]
    full_values[:, ~full_mask] = 0.0
    return full_values, full_mask


def _set_refractive_index_material(material, refractive_index: str, extinction_coefficient: str) -> None:
    material.propertyGroup("RefractiveIndex").set("n", "")
    material.propertyGroup("RefractiveIndex").set("ki", "")
    material.propertyGroup("RefractiveIndex").set(
        "n",
        JArray(JString, 1)(_isotropic_tensor(refractive_index)),
    )
    material.propertyGroup("RefractiveIndex").set(
        "ki",
        JArray(JString, 1)(_isotropic_tensor(extinction_coefficient)),
    )


class SimulationRun:
    def __init__(
        self,
        config: SimulationConfig | None = None,
        progress_log_path: str | os.PathLike[str] | None = None,
    ):
        self.config = config or SimulationConfig()
        self.progress_log_path = os.fspath(progress_log_path) if progress_log_path is not None else None
        self.client = mph.start()
        self.model = self.client.create("Model")
        self._plotting_initialized = False

    def __del__(self):
        self.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear()
        return False

    def clear(self):
        if getattr(self, "model", None) is not None:
            self.client.remove(self.model)
            self.model = None

    def build_geometry(
        self,
        boundary: BoundarySpec,
        layers: Sequence[LayerSpec],
        simulation_mode: SimulationMode,
        k: dict[str, float] | None = None,
    ) -> None:
        """Build a multilayer COMSOL model without solving it."""
        boundary, layers, simulation_mode, k_values, use_floquet = _validate_build_inputs(
            boundary,
            layers,
            simulation_mode,
            k,
        )

        jmodel = self.model.java
        config = self.config
        boundary_def = _boundary_definition(boundary)
        layer_stack = _make_layer_stack(layers)

        jmodel.param().set("a", f"{boundary.a} [um]")
        if boundary.shape in {"rectangle", "quarter_rectangle"}:
            jmodel.param().set("b", f"{boundary.b} [um]")
        jmodel.param().set("H_air", config.air_height)
        jmodel.param().set("lda0", config.wavelength)
        jmodel.param().set("selection_tol", config.selection_tolerance)
        jmodel.param().set("Ga", boundary_def.ga_expr)
        if boundary_def.gb_expr is not None:
            jmodel.param().set("Gb", boundary_def.gb_expr)
        jmodel.param().set("kx", f"{k_values['kx']}*Ga")
        ky_scale = "Gb" if boundary_def.gb_expr is not None else "Ga"
        jmodel.param().set("ky", f"{k_values['ky']}*{ky_scale}")

        for layer, entry in zip(layers, layer_stack.entries):
            jmodel.param().set(entry.height_param, layer.height)
            jmodel.param().set(entry.z_lo_param, entry.z_lo_expr)
            jmodel.param().set(entry.z_hi_param, entry.z_hi_expr)
            jmodel.param().set(entry.z_mid_param, entry.z_mid_expr)
        jmodel.param().set(layer_stack.z_slab_top_param, layer_stack.z_slab_top_expr)

        selection_tol = float(jmodel.param().evaluate("selection_tol", "um"))
        for layer, entry in zip(layers, layer_stack.entries):
            for hole_idx, max_radius in enumerate(layer.hole_selection_max_radii):
                if max_radius <= 2.0 * selection_tol:
                    raise ValueError(
                        f"LayerSpec.holes[{hole_idx}] in layer {entry.index} cannot fit a "
                        f"COMSOL selection ball with radius {selection_tol:g} [um] and 2x margin; "
                        f"maximum interior radius is {max_radius:g} [um]"
                    )

        jmodel.component().create("comp1", True)
        jmodel.component("comp1").geom().create("geom1", 3)
        jmodel.component("comp1").mesh().create("mesh1")
        jmodel.component("comp1").geom("geom1").lengthUnit("um")

        polygon_counter = 1

        def next_polygon_tag() -> str:
            nonlocal polygon_counter
            tag = f"pol{polygon_counter}"
            polygon_counter += 1
            return tag

        def create_workplane_with_polygons(
            workplane_tag: str,
            z_expr: str,
            footprint_label: str,
            holes: Sequence[np.ndarray],
        ) -> None:
            jmodel.component("comp1").geom("geom1").create(workplane_tag, "WorkPlane")
            jmodel.component("comp1").geom("geom1").feature(workplane_tag).set("quickz", z_expr)
            jmodel.component("comp1").geom("geom1").feature(workplane_tag).set("unite", True)

            footprint_tag = next_polygon_tag()
            wp_geom = jmodel.component("comp1").geom("geom1").feature(workplane_tag).geom()
            wp_geom.create(footprint_tag, "Polygon")
            wp_geom.feature(footprint_tag).label(footprint_label)
            wp_geom.feature(footprint_tag).set("source", "table")
            wp_geom.feature(footprint_tag).set("table", JArray(JString, 2)(boundary_def.footprint_table))

            for hole_idx, hole in enumerate(holes):
                hole_tag = next_polygon_tag()
                hole_table = [[f"{x} [um]", f"{y} [um]"] for x, y in hole]
                wp_geom.create(hole_tag, "Polygon")
                wp_geom.feature(hole_tag).label(f"Hole_{hole_idx + 1}")
                wp_geom.feature(hole_tag).set("source", "table")
                wp_geom.feature(hole_tag).set("table", JArray(JString, 2)(hole_table))

        for layer, entry in zip(layers, layer_stack.entries):
            workplane_tag = f"wp{entry.index + 1}"
            create_workplane_with_polygons(
                workplane_tag,
                entry.z_lo_param,
                f"{boundary_def.footprint_label} layer {entry.index + 1}",
                layer.holes,
            )
            extrude_tag = f"ext{entry.index + 1}"
            jmodel.component("comp1").geom("geom1").create(extrude_tag, "Extrude")
            jmodel.component("comp1").geom("geom1").feature(extrude_tag).setIndex("distance", entry.model_height_expr, 0)
            jmodel.component("comp1").geom("geom1").feature(extrude_tag).selection("input").set(workplane_tag)

        top_workplane_tag = f"wp{len(layers) + 1}"
        create_workplane_with_polygons(
            top_workplane_tag,
            layer_stack.z_slab_top_param,
            f"{boundary_def.footprint_label} air",
            [],
        )
        top_extrude_tag = f"ext{len(layers) + 1}"
        jmodel.component("comp1").geom("geom1").create(top_extrude_tag, "Extrude")
        jmodel.component("comp1").geom("geom1").feature(top_extrude_tag).set(
            "distance",
            JArray(JString, 1)([config.air_layer_distance, config.pml_layer_distance]),
        )
        jmodel.component("comp1").geom("geom1").feature(top_extrude_tag).set(
            "scale",
            JArray(JDouble, 2)([[1.0, 1.0], [1.0, 1.0]]),
        )
        jmodel.component("comp1").geom("geom1").feature(top_extrude_tag).set(
            "displ",
            JArray(JDouble, 2)([[0.0, 0.0], [0.0, 0.0]]),
        )
        jmodel.component("comp1").geom("geom1").feature(top_extrude_tag).set(
            "twist",
            JArray(JInt, 1)([0, 0]),
        )
        jmodel.component("comp1").geom("geom1").feature(top_extrude_tag).selection("input").set(top_workplane_tag)

        jmodel.component("comp1").geom("geom1").run()

        H_air_val = float(jmodel.param().evaluate("H_air", "um"))
        z_slab_top = float(jmodel.param().evaluate(layer_stack.z_slab_top_param, "um"))
        selection_tol = float(jmodel.param().evaluate("selection_tol", "um"))

        z_pml_lo = z_slab_top + H_air_val * config.air_layer_top_factor
        z_top = z_slab_top + H_air_val * config.pml_layer_top_factor
        z_air_mid = z_slab_top + H_air_val * config.air_layer_top_factor / 2.0
        large = _boundary_large_extent(boundary)

        def box_sel(name, label, dim, xmn, xmx, ymn, ymx, zmn, zmx, cond="inside") -> None:
            jmodel.component("comp1").selection().create(name, "Box")
            selection = jmodel.component("comp1").selection(name)
            selection.set("entitydim", JInt(dim))
            selection.set("xmin", xmn)
            selection.set("xmax", xmx)
            selection.set("ymin", ymn)
            selection.set("ymax", ymx)
            selection.set("zmin", zmn)
            selection.set("zmax", zmx)
            selection.set("condition", cond)
            if label:
                selection.label(label)

        def exterior_ball_sel(name, label, x, y, z) -> None:
            jmodel.component("comp1").selection().create(name, "Ball")
            selection = jmodel.component("comp1").selection(name)
            selection.set("entitydim", JInt(2))
            selection.set("inputent", "selections")
            selection.set("input", JArray(JString, 1)(["sel_exterior_boundaries"]))
            selection.set("condition", "intersects")
            selection.set("groupcontang", "on")
            selection.set("posx", x)
            selection.set("posy", y)
            selection.set("posz", z)
            selection.set("r", selection_tol)
            if label:
                selection.label(label)

        jmodel.component("comp1").selection().create("sel_exterior_boundaries")
        selection = jmodel.component("comp1").selection("sel_exterior_boundaries")
        selection.geom("geom1", JInt(3), JInt(2), JArray(JString, 1)(["exterior"]))
        selection.all()
        selection.label("exterior_boundaries")

        if boundary.shape in {"quarter_hexagon_boundary", "quarter_rectangle"}:
            box_sel(
                "sel_bottom",
                "bottom",
                2,
                -large,
                large,
                -large,
                large,
                -selection_tol,
                selection_tol,
                "inside",
            )
            box_sel(
                "sel_top_face",
                "top",
                2,
                -large,
                large,
                -large,
                large,
                z_top - selection_tol,
                z_top + selection_tol,
                "inside",
            )
        else:
            face_seed_x, face_seed_y = _boundary_face_seed_xy(boundary)
            exterior_ball_sel("sel_bottom", "bottom", face_seed_x, face_seed_y, 0.0)
            exterior_ball_sel("sel_top_face", "top", face_seed_x, face_seed_y, z_top)

        for layer, entry in zip(layers, layer_stack.entries):
            z_lo = float(jmodel.param().evaluate(entry.z_lo_param, "um"))
            z_hi = float(jmodel.param().evaluate(entry.z_hi_param, "um"))
            z_mid = float(jmodel.param().evaluate(entry.z_mid_param, "um"))

            extent_name = f"sel_layer_{entry.index}_extent"
            box_sel(
                extent_name,
                f"layer_{entry.index}_extent",
                3,
                -large,
                large,
                -large,
                large,
                z_lo - selection_tol,
                z_hi + selection_tol,
            )

            hole_domain_names = []
            for hole_idx, hole_center in enumerate(layer.hole_selection_centers):
                hole_sel_name = f"sel_layer_{entry.index}_hole_dom_{hole_idx}"
                jmodel.component("comp1").selection().create(hole_sel_name, "Ball")
                hole_selection = jmodel.component("comp1").selection(hole_sel_name)
                hole_selection.set("entitydim", JInt(3))
                hole_selection.set("condition", "intersects")
                hole_selection.set("posx", float(hole_center[0]))
                hole_selection.set("posy", float(hole_center[1]))
                hole_selection.set("posz", z_mid)
                hole_selection.set("r", selection_tol)
                hole_domain_names.append(hole_sel_name)

            material_selection_name = f"sel_layer_{entry.index}_mat_dom"
            if hole_domain_names:
                hole_union_name = f"sel_layer_{entry.index}_hole_doms"
                jmodel.component("comp1").selection().create(hole_union_name, "Union")
                hole_union = jmodel.component("comp1").selection(hole_union_name)
                hole_union.set("entitydim", JInt(3))
                hole_union.set("input", JArray(JString, 1)(hole_domain_names))
                hole_union.label(f"layer_{entry.index}_hole_doms")

                jmodel.component("comp1").selection().create(material_selection_name, "Difference")
                material_selection = jmodel.component("comp1").selection(material_selection_name)
                material_selection.set("entitydim", JInt(3))
                material_selection.set("add", JArray(JString, 1)([extent_name]))
                material_selection.set("subtract", JArray(JString, 1)([hole_union_name]))
                material_selection.label(f"layer_{entry.index}_mat_dom")
            else:
                jmodel.component("comp1").selection().create(material_selection_name, "Union")
                material_selection = jmodel.component("comp1").selection(material_selection_name)
                material_selection.set("entitydim", JInt(3))
                material_selection.set("input", JArray(JString, 1)([extent_name]))
                material_selection.label(f"layer_{entry.index}_mat_dom")

        box_sel(
            "sel_top_pml_dom",
            "top_pml_dom",
            3,
            -large,
            large,
            -large,
            large,
            z_pml_lo - selection_tol,
            z_top + selection_tol,
        )
        jmodel.component("comp1").selection().create("sel_pml_dom", "Union")
        selection = jmodel.component("comp1").selection("sel_pml_dom")
        selection.set("entitydim", JInt(3))
        selection.set("input", JArray(JString, 1)(["sel_top_pml_dom"]))
        selection.label("pml_dom")

        jmodel.component("comp1").selection().create("sel_all_dom")
        selection = jmodel.component("comp1").selection("sel_all_dom")
        selection.geom("geom1", JInt(3))
        selection.all()
        selection.label("all_domains")

        jmodel.component("comp1").selection().create("sel_physical_dom", "Difference")
        selection = jmodel.component("comp1").selection("sel_physical_dom")
        selection.set("entitydim", JInt(3))
        selection.set("add", JArray(JString, 1)(["sel_all_dom"]))
        selection.set("subtract", JArray(JString, 1)(["sel_pml_dom"]))
        selection.label("physical_domains")

        wall_selection_names = []
        wall_seeds = _wall_seed_points(boundary, z_air_mid)
        for wall_key, seeds in wall_seeds.items():
            sub_names = []
            for side_idx, (x, y, z) in enumerate(seeds):
                selection_name = f"sel_{wall_key}_w{side_idx}"
                exterior_ball_sel(selection_name, None, x, y, z)
                sub_names.append(selection_name)
            union_name = f"sel_{wall_key}"
            jmodel.component("comp1").selection().create(union_name, "Union")
            wall_selection = jmodel.component("comp1").selection(union_name)
            wall_selection.set("entitydim", JInt(2))
            wall_selection.label(f"{wall_key}_walls")
            wall_selection.set("input", JArray(JString, 1)(sub_names))
            wall_selection_names.append(wall_key)

        jmodel.component("comp1").material().create("mat1", "Common")
        air_material = jmodel.component("comp1").material("mat1")
        air_material.label("Air")
        air_material.propertyGroup().create("RefractiveIndex", "Refractive Index")
        _set_refractive_index_material(
            air_material,
            config.air_refractive_index,
            config.air_extinction_coefficient,
        )

        for layer, entry in zip(layers, layer_stack.entries):
            material_tag = f"mat{entry.index + 2}"
            jmodel.component("comp1").material().create(material_tag, "Common")
            material = jmodel.component("comp1").material(material_tag)
            material.label(layer.label or f"Layer {entry.index}")
            material.selection().named(f"sel_layer_{entry.index}_mat_dom")
            material.propertyGroup().create("RefractiveIndex", "Refractive Index")
            _set_refractive_index_material(
                material,
                layer.refractive_index,
                layer.extinction_coefficient,
            )

        jmodel.component("comp1").coordSystem().create("pml1", "PML")
        jmodel.component("comp1").coordSystem("pml1").selection().named("sel_pml_dom")

        jmodel.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1")
        jmodel.component("comp1").physics("ewfd").create("sctr1", "Scattering", 2)
        jmodel.component("comp1").physics("ewfd").feature("sctr1").selection().named("sel_top_face")

        mode_type = config.mode_type.upper()
        if mode_type == "TM":
            bottom_feature = "pec_bottom"
            bottom_condition = "PerfectElectricConductor"
        elif mode_type == "TE":
            bottom_feature = "pmc_bottom"
            bottom_condition = "PerfectMagneticConductor"
        else:
            raise ValueError("mode_type must be 'TE' or 'TM'")

        jmodel.component("comp1").physics("ewfd").create(bottom_feature, bottom_condition, 2)
        jmodel.component("comp1").physics("ewfd").feature(bottom_feature).selection().named("sel_bottom")

        if simulation_mode == "unit_cell":
            for idx, wall_key in enumerate(wall_selection_names, start=1):
                feature_tag = f"pc{idx}"
                jmodel.component("comp1").physics("ewfd").create(feature_tag, "PeriodicCondition", 2)
                jmodel.component("comp1").physics("ewfd").feature(feature_tag).selection().named(f"sel_{wall_key}")
                if use_floquet:
                    jmodel.component("comp1").physics("ewfd").feature(feature_tag).set("PeriodicType", "Floquet")
                    jmodel.component("comp1").physics("ewfd").feature(feature_tag).set(
                        "kFloquet",
                        JArray(JString, 1)(["kx", "ky", "0"]),
                    )
        elif simulation_mode == "strip_1d":
            jmodel.component("comp1").physics("ewfd").create("sctr2", "Scattering", 2)
            jmodel.component("comp1").physics("ewfd").feature("sctr2").selection().named("sel_side_y")

            jmodel.component("comp1").physics("ewfd").create("pc1", "PeriodicCondition", 2)
            jmodel.component("comp1").physics("ewfd").feature("pc1").selection().named("sel_side_x")
            if use_floquet:
                jmodel.component("comp1").physics("ewfd").feature("pc1").set("PeriodicType", "Floquet")
                jmodel.component("comp1").physics("ewfd").feature("pc1").set(
                    "kFloquet",
                    JArray(JString, 1)(["kx", "0", "0"]),
                )
        elif simulation_mode == "finite_quarter":
            symmetry_conditions = {
                "symmetry_x_axis": config.symmetry_x_boundary,
                "symmetry_y_axis": config.symmetry_y_boundary,
            }
            for axis_name, condition in symmetry_conditions.items():
                normalized_condition = "" if condition is None else str(condition).upper()
                if normalized_condition not in {"PEC", "PMC"}:
                    raise ValueError(
                        f"{axis_name} boundary condition must be 'PEC' or 'PMC', "
                        f"got {condition!r}"
                    )
                for candidate in ("PEC", "PMC"):
                    feature_tag = f"bc_{axis_name}_{candidate.lower()}"
                    feature_type = (
                        "PerfectElectricConductor"
                        if candidate == "PEC"
                        else "PerfectMagneticConductor"
                    )
                    jmodel.component("comp1").physics("ewfd").create(
                        feature_tag,
                        feature_type,
                        2,
                    )
                    feature = jmodel.component("comp1").physics("ewfd").feature(
                        feature_tag
                    )
                    feature.selection().named(f"sel_{axis_name}")
                    feature.active(candidate == normalized_condition)

            outer_wall_keys = (
                ("outer_rect_x", "outer_rect_y")
                if boundary.shape == "quarter_rectangle"
                else ("outer_hex_ne", "outer_hex_n")
            )
            for idx, wall_key in enumerate(
                outer_wall_keys,
                start=2,
            ):
                feature_tag = f"sctr{idx}"
                jmodel.component("comp1").physics("ewfd").create(
                    feature_tag,
                    "Scattering",
                    2,
                )
                jmodel.component("comp1").physics("ewfd").feature(
                    feature_tag
                ).selection().named(f"sel_{wall_key}")
        else:
            for idx, wall_key in enumerate(wall_selection_names, start=2):
                feature_tag = f"sctr{idx}"
                jmodel.component("comp1").physics("ewfd").create(feature_tag, "Scattering", 2)
                jmodel.component("comp1").physics("ewfd").feature(feature_tag).selection().named(f"sel_{wall_key}")

        jmodel.study().create("std1")
        jmodel.study("std1").create("eig", "Eigenfrequency")
        jmodel.study("std1").feature("eig").set("shift", config.eigenfrequency_shift)
        jmodel.study("std1").feature("eig").set("neigsactive", True)
        jmodel.study("std1").feature("eig").set("neigs", JInt(config.eigenmode_count))

    def run_simulation(
        self,
        mesh_save_path: str | os.PathLike[str] | None = None,
        mesh_builder: Callable[[object, SimulationConfig], None] | None = None,
    ) -> None:
        """Mesh, solve, and initialize result datasets for the current model."""
        jmodel = self.model.java
        config = self.config
        mesh_save_path_str = os.fspath(mesh_save_path) if mesh_save_path is not None else None

        model_util = None
        progress_log_path = getattr(self, "progress_log_path", None)
        if progress_log_path is not None:
            progress_log_dir = os.path.dirname(progress_log_path)
            if progress_log_dir:
                os.makedirs(progress_log_dir, exist_ok=True)
            model_util = JClass("com.comsol.model.util.ModelUtil")
            model_util.showProgress(progress_log_path)

        try:
            if mesh_builder is None:
                jmodel.component("comp1").mesh("mesh1").autoMeshSize(config.mesh_auto_size)
            else:
                mesh_builder(jmodel, config)
            jmodel.component("comp1").mesh("mesh1").run()
            if mesh_save_path_str is not None:
                self.model.save(mesh_save_path_str)
            jmodel.study("std1").createAutoSequences("all")
            jmodel.sol("sol1").runAll()
        finally:
            if model_util is not None:
                model_util.showProgress(False)

        jmodel.result().numerical().create("gev1", "EvalGlobal")
        jmodel.result().numerical("gev1").label("Eigenfrequencies (ewfd)")
        jmodel.result().numerical("gev1").set("data", "dset1")
        jmodel.result().numerical("gev1").set("expr", ["ewfd.omega/2/pi", "ewfd.damp/2/pi", "ewfd.Qfactor"])
        jmodel.result().numerical("gev1").set("unit", ["THz", "THz", "1"])

        jmodel.result().table().create("tbl1", "Table")
        jmodel.result().numerical("gev1").set("table", "tbl1")
        jmodel.result().numerical("gev1").run()
        jmodel.result().numerical("gev1").setResult()

        jmodel.result().dataset().create("cpl1", "CutPlane")
        jmodel.result().dataset("cpl1").set("quickplane", "xy")
        jmodel.result().dataset("cpl1").set("quickz", "0")

        jmodel.result().dataset().create("cpl2", "CutPlane")
        jmodel.result().dataset("cpl2").set("quickplane", "xy")
        jmodel.result().dataset("cpl2").set("quickz", f"z_slab_top + {config.air_cutplane_z}")

        jmodel.result().dataset().create("cpl3", "CutPlane")
        jmodel.result().dataset("cpl3").set("quickplane", "yz")

        jmodel.result().dataset().create("cpl4", "CutPlane")
        jmodel.result().dataset("cpl4").set("quickplane", "xz")

        self.plane_datasets = {
            "center": "cpl1",
            "air": "cpl2",
            "yz": "cpl3",
            "xz": "cpl4",
        }
        jmodel.result().numerical().create("int1", "Interp")
        jmodel.result().numerical().create("int_grid", "Interp")
        jmodel.result().numerical().create("max3d", "MaxVolume")

    def set_finite_quarter_symmetry(
        self,
        x_boundary: str,
        y_boundary: str,
    ) -> None:
        """Activate one PEC/PMC feature per finite-quarter symmetry axis."""
        conditions = {
            "symmetry_x_axis": normalize_symmetry_boundary(x_boundary),
            "symmetry_y_axis": normalize_symmetry_boundary(y_boundary),
        }
        physics = self.model.java.component("comp1").physics("ewfd")
        for axis_name, active_condition in conditions.items():
            for candidate in ("PEC", "PMC"):
                physics.feature(
                    f"bc_{axis_name}_{candidate.lower()}"
                ).active(candidate == active_condition)

    def prepare_reused_finite_quarter_model(self) -> dict[str, object]:
        """Upgrade a solved legacy quarter model for mesh-preserving reruns."""
        jmodel = self.model.java
        required = {
            "study": (jmodel.study(), ("std1",)),
            "solution": (jmodel.sol(), ("sol1",)),
            "dataset": (
                jmodel.result().dataset(),
                ("dset1", "cpl1", "cpl2", "cpl3", "cpl4"),
            ),
            "numerical": (
                jmodel.result().numerical(),
                ("gev1", "int1", "int_grid", "max3d"),
            ),
        }
        for kind, (container, tags) in required.items():
            missing = [tag for tag in tags if not container.hasTag(tag)]
            if missing:
                raise ValueError(
                    "Cannot reuse finite-quarter model: missing "
                    f"{kind} tag(s) {missing!r}"
                )

        physics = jmodel.component("comp1").physics("ewfd")
        upgraded: list[str] = []
        disabled_legacy: list[str] = []
        for axis_name in ("symmetry_x_axis", "symmetry_y_axis"):
            legacy_tag = f"bc_{axis_name}"
            if physics.feature().hasTag(legacy_tag):
                physics.feature(legacy_tag).active(False)
                disabled_legacy.append(legacy_tag)
            for condition, feature_type in (
                ("PEC", "PerfectElectricConductor"),
                ("PMC", "PerfectMagneticConductor"),
            ):
                feature_tag = f"bc_{axis_name}_{condition.lower()}"
                if not physics.feature().hasTag(feature_tag):
                    physics.create(feature_tag, feature_type, 2)
                    upgraded.append(feature_tag)
                feature = physics.feature(feature_tag)
                feature.selection().named(f"sel_{axis_name}")
                feature.active(False)

        study = jmodel.study("std1").feature("eig")
        study.set("shift", self.config.eigenfrequency_shift)
        study.set("neigsactive", True)
        study.set("neigs", JInt(self.config.eigenmode_count))
        solution = jmodel.sol("sol1")
        if not solution.feature().hasTag("e1"):
            raise ValueError(
                "Cannot reuse finite-quarter model: sol1 is missing "
                "the automatic eigenvalue solver feature 'e1'"
            )
        eigen_solver = solution.feature("e1")
        eigen_solver.set("shift", self.config.eigenfrequency_shift)
        eigen_solver.set("neigs", JInt(self.config.eigenmode_count))
        return {
            "disabled_legacy_boundary_features": disabled_legacy,
            "created_boundary_features": upgraded,
            "eigenfrequency_shift": self.config.eigenfrequency_shift,
            "eigenmode_count": self.config.eigenmode_count,
        }

    def rerun_eigenfrequency_solution(self) -> None:
        """Recompute sol1 without rebuilding geometry or mesh."""
        solution = self.model.java.sol("sol1")
        model_util = None
        progress_log_path = getattr(self, "progress_log_path", None)
        if progress_log_path is not None:
            progress_log_dir = os.path.dirname(progress_log_path)
            if progress_log_dir:
                os.makedirs(progress_log_dir, exist_ok=True)
            model_util = JClass("com.comsol.model.util.ModelUtil")
            model_util.showProgress(progress_log_path)
        try:
            solution.clearSolutionData()
            solution.runAll()
        finally:
            if model_util is not None:
                model_util.showProgress(False)

    def copy_solution(
        self,
        solution_tag: str,
        dataset_tag: str,
        *,
        label: str | None = None,
    ) -> dict[str, str]:
        """Copy the solved working sequence and expose it through a dataset."""
        jmodel = self.model.java
        if jmodel.sol().hasTag(solution_tag):
            raise ValueError(f"COMSOL solution tag already exists: {solution_tag}")
        if jmodel.result().dataset().hasTag(dataset_tag):
            raise ValueError(f"COMSOL dataset tag already exists: {dataset_tag}")
        jmodel.sol("sol1").copySolution(solution_tag)
        copied = jmodel.sol(solution_tag)
        if label is not None:
            copied.label(label)
        jmodel.result().dataset().create(dataset_tag, "Solution")
        dataset = jmodel.result().dataset(dataset_tag)
        dataset.set("solution", solution_tag)
        if label is not None:
            dataset.label(label)
        return {"solution_tag": solution_tag, "dataset_tag": dataset_tag}

    def build_and_run(
        self,
        boundary: BoundarySpec,
        layers: Sequence[LayerSpec],
        simulation_mode: SimulationMode = "unit_cell",
        k: dict[str, float] | None = None,
        mesh_save_path: str | os.PathLike[str] | None = None,
        mesh_builder: Callable[[object, SimulationConfig], None] | None = None,
    ) -> None:
        """Build and run a multilayer COMSOL eigenfrequency simulation."""
        self.build_geometry(boundary, layers, simulation_mode, k)
        self.run_simulation(mesh_save_path=mesh_save_path, mesh_builder=mesh_builder)

    def setup_plotting(self):
        """Lazily initialize COMSOL plot groups and export objects."""
        if self._plotting_initialized:
            return
        jmodel = self.model.java

        jmodel.result().create("pg1", "PlotGroup2D")
        jmodel.result("pg1").label("2D Field (ewfd)")
        jmodel.result("pg1").run()
        jmodel.result("pg1").create("surf1", "Surface")
        jmodel.result("pg1").feature("surf1").set("expr", "ewfd.normH")
        jmodel.result("pg1").run()
        jmodel.result().export().create("plot1", "pg1", "surf1", "Plot")
        jmodel.result().export().create("img1", "pg1", "Image")
        jmodel.result().export("img1").set("target", "file")

        jmodel.result().create("pg2", "PlotGroup3D")
        jmodel.result("pg2").run()
        jmodel.result("pg2").label("3D Plot")
        jmodel.result("pg2").create("mslc1", "Multislice")
        jmodel.result("pg2").feature("mslc1").set("multiplanexmethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("xcoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("multiplaneymethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("ycoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("multiplanezmethod", "coord")
        jmodel.result("pg2").feature("mslc1").set("zcoord", JInt(0))
        jmodel.result("pg2").feature("mslc1").set("expr", "ewfd.normE")
        jmodel.result("pg2").run()
        jmodel.result().export().create("img2", "pg2", "Image")
        jmodel.result().export("img2").set("target", "file")

        self._plotting_initialized = True

    def export_eigenfrequencies(self, save_path):
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java
        jmodel.result().table("tbl1").save(os.path.join(save_path, "eigenfrequencies.txt"))

    def get_eigenfrequencies(self, dataset: str = "dset1"):
        jmodel = self.model.java
        jmodel.result().numerical("gev1").set("data", dataset)
        jmodel.result().numerical("gev1").run()
        data = np.asarray(jmodel.result().numerical("gev1").getReal()).T
        return data

    def export_2d_fields(
        self,
        eigenmode_idx,
        expr,
        expr_name,
        save_path,
        plane="center",
        export_data=True,
        export_image=True,
        image_filename=None,
        color_table=None,
        color_scale_mode=None,
        show_colorbar=None,
    ):
        self.setup_plotting()
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java
        jmodel.result("pg1").run()
        jmodel.result("pg1").set("data", self.plane_datasets[plane])
        jmodel.result("pg1").set("looplevel", JArray(JInt, 1)([eigenmode_idx + 1]))
        surface = jmodel.result("pg1").feature("surf1")
        surface.set("expr", expr)
        if color_table is not None:
            surface.set("colortable", color_table)
        if color_scale_mode is not None:
            surface.set("colorscalemode", color_scale_mode)
        if show_colorbar is not None:
            surface.set("colorlegend", "on" if show_colorbar else "off")
        jmodel.result("pg1").run()
        if export_data:
            jmodel.result().export("plot1").set(
                "filename",
                os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_{plane}_2d.txt"),
            )
            jmodel.result().export("plot1").run()
        if export_image:
            filename = image_filename or f"{eigenmode_idx:02d}_{expr_name}_{plane}_2d.png"
            from .field_plotting import save_hz_simulation_plot

            if save_hz_simulation_plot(
                self, eigenmode_idx, expr, os.path.join(save_path, filename), plane,
            ):
                return
            image_export = jmodel.result().export("img1")
            if show_colorbar is not None:
                image_export.set("options2d", "on")
                image_export.set("legend2d", "on" if show_colorbar else "off")
            image_export.set(
                "pngfilename",
                os.path.join(save_path, filename),
            )
            image_export.run()

    def get_2d_fields(self, eigenmode_idx, expr, plane="center"):
        jmodel = self.model.java
        jmodel.result().numerical("int1").set("data", self.plane_datasets[plane])
        jmodel.result().numerical("int1").set("expr", JArray(JString, 1)([expr]))
        jmodel.result().numerical("int1").set("solnum", JArray(JInt, 1)([eigenmode_idx + 1]))
        jmodel.result().numerical("int1").run()

        coords = np.asarray(jmodel.result().numerical("int1").getCoordinates()).T
        values = np.asarray(jmodel.result().numerical("int1").getData()).reshape(-1)
        return coords, values

    def get_fields_at_coordinates(self, eigenmode_idx, expressions, coordinates, dataset="dset1"):
        expressions = list(expressions)
        if not expressions:
            raise ValueError("expressions must not be empty")
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("coordinates must have shape (N, 3)")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates must contain only finite values")

        jmodel = self.model.java
        try:
            interpolation = jmodel.result().numerical("int_grid")
        except Exception:
            jmodel.result().numerical().create("int_grid", "Interp")
            interpolation = jmodel.result().numerical("int_grid")

        interpolation.set("data", dataset)
        interpolation.set("expr", JArray(JString, 1)(expressions))
        interpolation.set("solnum", JArray(JInt, 1)([eigenmode_idx + 1]))
        interpolation.setInterpolationCoordinates(
            JArray(JDouble, 2)(coordinates.T.tolist())
        )
        interpolation.run()

        expected_shape = (len(expressions), len(coordinates))

        def solution_slice(raw, name):
            values = np.asarray(raw, dtype=float)
            if values.ndim == 3:
                if values.shape[1] != 1:
                    raise ValueError(f"Expected one solution axis in {name}, got {values.shape}")
                values = values[:, 0, :]
            if values.shape != expected_shape:
                raise ValueError(
                    f"Expected {name} interpolation shape {expected_shape}, got {values.shape}"
                )
            return values

        real = solution_slice(interpolation.getData(), "real")
        imag = solution_slice(interpolation.getImagData(), "imaginary")
        return real + 1j * imag

    def get_3d_field_max_value(self, eigenmode_idx, expr):
        jmodel = self.model.java
        max_feature = jmodel.result().numerical("max3d")
        position_scale = self._maxvolume_position_scale_to_um()

        max_feature.set("data", "dset1")
        max_feature.set("expr", JArray(JString, 1)([expr]))
        max_feature.set("includepos", "on")
        max_feature.set("innerinput", "manual")
        max_feature.set("solnum", JArray(JInt, 1)([eigenmode_idx + 1]))
        max_feature.selection().named("sel_physical_dom")

        result = np.asarray(max_feature.computeResult(), dtype=float).reshape(-1)
        if result.size < 4:
            raise ValueError(f"Expected max value plus 3D position for {expr}, got {result}")

        peak = {
            "value": float(result[0]),
            "x": float(result[1]) * position_scale,
            "y": float(result[2]) * position_scale,
            "z": float(result[3]) * position_scale,
        }

        z_slab_top = float(jmodel.param().evaluate("z_slab_top", "um"))
        H_air = float(jmodel.param().evaluate("H_air", "um"))
        tolerance = float(jmodel.param().evaluate("selection_tol", "um"))
        z_top = z_slab_top + H_air * self.config.pml_layer_top_factor
        if not -tolerance <= peak["z"] <= z_top + tolerance:
            raise ValueError(
                f"MaxVolume position for {expr} is outside model z range: "
                f"z={peak['z']}, expected {-tolerance} <= z <= {z_top + tolerance}. "
                f"Raw result: {result}"
            )
        return peak

    def _maxvolume_position_scale_to_um(self) -> float:
        if hasattr(self, "_maxvolume_position_scale"):
            return self._maxvolume_position_scale

        jmodel = self.model.java
        probe_tag = "max3d_unit_probe"
        jmodel.result().numerical().create(probe_tag, "MaxVolume")
        try:
            max_feature = jmodel.result().numerical(probe_tag)
            max_feature.set("data", "dset1")
            max_feature.set("expr", JArray(JString, 1)(["x"]))
            max_feature.set("unit", JArray(JString, 1)(["um"]))
            max_feature.set("includepos", "on")
            max_feature.set("innerinput", "manual")
            max_feature.set("solnum", JArray(JInt, 1)([1]))
            max_feature.selection().named("sel_physical_dom")
            result = np.asarray(max_feature.computeResult(), dtype=float).reshape(-1)
        finally:
            jmodel.result().numerical().remove(probe_tag)

        self._maxvolume_position_scale = abs(float(result[0]) / float(result[1]))

        return self._maxvolume_position_scale

    def compute_polarization(self, eigenmode_idx):
        jmodel = self.model.java

        kx = float(jmodel.param().evaluate("kx", "um^-1"))
        ky = float(jmodel.param().evaluate("ky", "um^-1"))

        pts_ex_re, ex_re = self.get_2d_fields(eigenmode_idx, "ewfd.Ex", "air")
        pts_ex_im, ex_im = self.get_2d_fields(eigenmode_idx, "ewfd.Ex*(-i)", "air")
        ex_comp = ex_re + 1j * interpolate_field(pts_ex_im, ex_im, pts_ex_re)

        pts_ey_re, ey_re = self.get_2d_fields(eigenmode_idx, "ewfd.Ey", "air")
        pts_ey_im, ey_im = self.get_2d_fields(eigenmode_idx, "ewfd.Ey*(-i)", "air")
        ey_comp = ey_re + 1j * interpolate_field(pts_ey_im, ey_im, pts_ey_re)

        phase_ex = np.exp(1j * (kx * pts_ex_re[:, 0] + ky * pts_ex_re[:, 1]))
        phase_ey = np.exp(1j * (kx * pts_ey_re[:, 0] + ky * pts_ey_re[:, 1]))

        area_x = float(np.real(integrate_field(pts_ex_re, np.ones(len(pts_ex_re)))))
        area_y = float(np.real(integrate_field(pts_ey_re, np.ones(len(pts_ey_re)))))

        cx = integrate_field(pts_ex_re, ex_comp * phase_ex) / np.sqrt(area_x)
        cy = integrate_field(pts_ey_re, ey_comp * phase_ey) / np.sqrt(area_y)

        return complex(cx), complex(cy)

    def export_3d_fields(self, eigenmode_idx, expr, expr_name, save_path):
        self.setup_plotting()
        os.makedirs(save_path, exist_ok=True)
        jmodel = self.model.java
        jmodel.result("pg2").run()
        jmodel.result("pg2").set("data", "dset1")
        jmodel.result("pg2").set("looplevel", JArray(JInt, 1)([eigenmode_idx + 1]))
        jmodel.result("pg2").feature("mslc1").set("expr", expr)
        jmodel.result("pg2").run()
        jmodel.result().export("img2").set(
            "pngfilename",
            os.path.join(save_path, f"{eigenmode_idx:02d}_{expr_name}_3d.png"),
        )
        jmodel.result().export("img2").run()
