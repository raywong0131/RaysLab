from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .cladding_shift_profile import (
    ELLIPSE_SHIFT_GEOMETRY,
    CladdingShiftGeometry,
    CladdingShiftProfile,
    cladding_shift_group,
    ellipse_inward_shift,
)
from .hex_lattice_utils import (
    LatticePoint,
    boundary_from_lattice_cells,
    cell_polygon,
    index_shell,
    lattice_points,
)


FINITE_GEOMETRIES = frozenset({"hex", "square"})
SQUARE_GROUP_HOLE_INDICES = {1: (2, 3, 4), 2: (0, 1, 5)}
_EPSILON = 1e-10


@dataclass(frozen=True)
class FiniteGeometrySpec:
    shape: str
    period: float
    cavity_layers: int
    cladding_layers: int
    empty_buffer_periods: float = 0.5

    def __post_init__(self) -> None:
        if self.shape not in FINITE_GEOMETRIES:
            raise ValueError("finite geometry shape must be 'hex' or 'square'")
        if not math.isfinite(self.period) or self.period <= 0.0:
            raise ValueError("finite geometry period must be positive and finite")
        if type(self.cavity_layers) is not int or self.cavity_layers <= 0:
            raise ValueError("cavity_layers must be a positive integer")
        if type(self.cladding_layers) is not int or self.cladding_layers <= 0:
            raise ValueError("cladding_layers must be a positive integer")
        if (
            not math.isfinite(self.empty_buffer_periods)
            or self.empty_buffer_periods < 0.0
        ):
            raise ValueError("empty_buffer_periods must be finite and non-negative")


@dataclass(frozen=True)
class FiniteFootprint:
    shape: str
    width: float
    height: float
    vertices: np.ndarray
    cavity_half_width: float | None = None
    cavity_half_height: float | None = None
    structured_half_width: float | None = None
    structured_half_height: float | None = None

    @property
    def boundary_shape(self) -> str:
        return "hexagon_boundary" if self.shape == "hex" else "rectangle"

    def contains_xy(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x_arr = np.asarray(x, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        tolerance = max(self.width, self.height, 1.0) * 1e-12
        if self.shape == "square":
            return (
                (np.abs(x_arr) <= self.width / 2.0 + tolerance)
                & (np.abs(y_arr) <= self.height / 2.0 + tolerance)
            )
        root3 = math.sqrt(3.0)
        return (
            (np.abs(y_arr) <= root3 * self.width / 4.0 + tolerance)
            & (
                root3 * np.abs(x_arr) + np.abs(y_arr)
                <= root3 * self.width / 2.0 + tolerance
            )
        )


@dataclass(frozen=True)
class FragmentPlacement:
    m: int
    n: int
    point: LatticePoint
    group_id: int
    region: str
    support: tuple[float, float, float, float]
    layer: int
    outside_x: bool
    outside_y: bool
    shift: np.ndarray
    polygons: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class FiniteGeometryPlan:
    spec: FiniteGeometrySpec
    placements: tuple[FragmentPlacement, ...]
    records: tuple[dict[str, object], ...]
    cavity_analysis_cells: tuple[LatticePoint, ...]
    cladding_parent_cells: tuple[LatticePoint, ...]
    footprint: FiniteFootprint
    validity_boundary: np.ndarray
    metadata: dict[str, object]

    @property
    def holes(self) -> tuple[np.ndarray, ...]:
        return tuple(np.asarray(record["polygon"], dtype=float) for record in self.records)


def square_fragment_cell_polygon(
    point: LatticePoint,
    group_id: int,
    period: float,
) -> np.ndarray:
    """Return the actual left/right half of a triangular-lattice hex cell."""
    if group_id not in SQUARE_GROUP_HOLE_INDICES:
        raise ValueError("square fragment group_id must be 1 or 2")
    polygon = cell_polygon(point.x, point.y, period)
    indices = (1, 2, 3, 4) if group_id == 1 else (4, 5, 0, 1)
    return polygon[list(indices)]


def square_cell_plot_records(
    plan: FiniteGeometryPlan,
) -> tuple[dict[str, object], ...]:
    """Merge compatible square fragments into true hex-cell plot records."""
    if plan.spec.shape != "square":
        raise ValueError("square cell plot records require a square geometry plan")
    by_parent: dict[tuple[int, int], list[FragmentPlacement]] = {}
    for placement in plan.placements:
        by_parent.setdefault((placement.m, placement.n), []).append(placement)

    records: list[dict[str, object]] = []
    for placements in by_parent.values():
        ordered = sorted(placements, key=lambda item: item.group_id)
        merge_parent = (
            len(ordered) == 2
            and ordered[0].region == ordered[1].region
            and ordered[0].layer == ordered[1].layer
            and np.allclose(ordered[0].shift, ordered[1].shift, rtol=0.0, atol=1e-12)
        )
        if merge_parent:
            first = ordered[0]
            records.append(
                {
                    "region": first.region,
                    "point": first.point,
                    "group_id": 0,
                    "layer": first.layer,
                    "outside_x": any(item.outside_x for item in ordered),
                    "outside_y": any(item.outside_y for item in ordered),
                    "modulation_shift": np.asarray(first.shift, dtype=float),
                    "polygon": cell_polygon(
                        first.point.x, first.point.y, plan.spec.period
                    ),
                }
            )
            continue
        for placement in ordered:
            records.append(
                {
                    "region": placement.region,
                    "point": placement.point,
                    "group_id": placement.group_id,
                    "layer": placement.layer,
                    "outside_x": placement.outside_x,
                    "outside_y": placement.outside_y,
                    "modulation_shift": np.asarray(placement.shift, dtype=float),
                    "polygon": square_fragment_cell_polygon(
                        placement.point,
                        placement.group_id,
                        plan.spec.period,
                    ),
                }
            )
    return tuple(records)


def _translate_polygon(polygon: np.ndarray, translation: Sequence[float]) -> np.ndarray:
    return np.asarray(polygon, dtype=float) + np.asarray(translation, dtype=float)


def _hex_boundary_vertices(a: float) -> np.ndarray:
    root3 = math.sqrt(3.0)
    return np.array(
        [
            [a / 2.0, 0.0],
            [a / 4.0, root3 * a / 4.0],
            [-a / 4.0, root3 * a / 4.0],
            [-a / 2.0, 0.0],
            [-a / 4.0, -root3 * a / 4.0],
            [a / 4.0, -root3 * a / 4.0],
        ],
        dtype=float,
    )


def _hex_boundary_size(points: Sequence[LatticePoint], period: float, pad: float) -> float:
    vertices = np.vstack([cell_polygon(point.x, point.y, period) for point in points])
    root3 = math.sqrt(3.0)
    x = vertices[:, 0]
    y = vertices[:, 1]
    minimum = max(
        float(4.0 * np.max(np.abs(y)) / root3),
        float(2.0 * np.max(np.abs(root3 * x + y)) / root3),
        float(2.0 * np.max(np.abs(root3 * x - y)) / root3),
    )
    return minimum + 2.0 * pad


def _ideal_bulk_hex_vertices(radius: int, period: float) -> np.ndarray:
    indices = np.array(
        [[radius, 0], [0, radius], [-radius, radius], [-radius, 0], [0, -radius], [radius, -radius]],
        dtype=float,
    )
    return np.column_stack(
        [
            period * (indices[:, 0] + 0.5 * indices[:, 1]),
            period * math.sqrt(3.0) * indices[:, 1] / 2.0,
        ]
    )


def _normal_fan_feature(point: LatticePoint, radius: int, period: float) -> tuple[str, int]:
    p = np.array([point.x, point.y], dtype=float)
    vertices = _ideal_bulk_hex_vertices(radius, period)
    best_distance = math.inf
    best_kind = "side"
    best_idx = 0
    for idx, start in enumerate(vertices):
        end = vertices[(idx + 1) % len(vertices)]
        edge = end - start
        t = float(np.dot(p - start, edge) / np.dot(edge, edge))
        t_clamped = float(np.clip(t, 0.0, 1.0))
        distance = float(np.sum((p - (start + t_clamped * edge)) ** 2))
        if distance < best_distance:
            best_distance = distance
            if 1e-10 < t_clamped < 1.0 - 1e-10:
                best_kind, best_idx = "side", idx
            else:
                best_kind = "corner"
                best_idx = idx if t_clamped <= 1e-10 else (idx + 1) % 6
    return best_kind, best_idx


def _normal_fan_direction(kind: str, idx: int, radius: int, period: float) -> np.ndarray:
    vertices = _ideal_bulk_hex_vertices(radius, period)

    def side_normal(side_idx: int) -> np.ndarray:
        edge = vertices[(side_idx + 1) % 6] - vertices[side_idx]
        normal = np.array([-edge[1], edge[0]], dtype=float)
        return normal / np.linalg.norm(normal)

    if kind == "side":
        return side_normal(idx)
    direction = side_normal((idx - 1) % 6) + side_normal(idx)
    return direction / np.linalg.norm(direction)


def _hex_shift(
    point: LatticePoint,
    spec: FiniteGeometrySpec,
    geometry: CladdingShiftGeometry,
    profile: CladdingShiftProfile,
    x_factor: float,
    y_factor: float,
) -> tuple[int, np.ndarray]:
    layer = point.shell - spec.cavity_layers + 1
    envelope = profile.envelope(layer)
    if geometry.kind == ELLIPSE_SHIFT_GEOMETRY:
        return layer, np.asarray(
            ellipse_inward_shift(
                point.x,
                point.y,
                spec.period,
                x_factor,
                y_factor,
                envelope=envelope,
            ),
            dtype=float,
        )
    cavity_reference_shell = max(spec.cavity_layers - 1, 1)
    kind, idx = _normal_fan_feature(point, cavity_reference_shell, spec.period)
    factor = x_factor if cladding_shift_group(kind, idx) == "x" else y_factor
    magnitude = envelope * factor * spec.period
    if kind == "corner":
        magnitude *= 2.0 / math.sqrt(3.0)
    return layer, magnitude * _normal_fan_direction(
        kind, idx, cavity_reference_shell, spec.period
    )


def _build_hex_plan(
    spec: FiniteGeometrySpec,
    cavity_holes: Sequence[np.ndarray],
    cladding_holes: Sequence[np.ndarray],
    geometry: CladdingShiftGeometry,
    profile: CladdingShiftProfile,
    x_factor: float,
    y_factor: float,
) -> FiniteGeometryPlan:
    cavity_points = lattice_points(spec.cavity_layers - 1, spec.period)
    cladding_points = lattice_points(
        spec.cavity_layers + spec.cladding_layers - 1,
        spec.period,
        min_shell=spec.cavity_layers,
    )
    placements: list[FragmentPlacement] = []
    records: list[dict[str, object]] = []
    for region, points, templates in (
        ("bulk", cavity_points, cavity_holes),
        ("cladding", cladding_points, cladding_holes),
    ):
        for point in points:
            if region == "bulk":
                layer, shift = 0, np.zeros(2, dtype=float)
            else:
                layer, shift = _hex_shift(
                    point, spec, geometry, profile, x_factor, y_factor
                )
            polygons = tuple(
                _translate_polygon(hole, (point.x + shift[0], point.y + shift[1]))
                for hole in templates
            )
            support_polygon = cell_polygon(point.x, point.y, spec.period)
            support = (
                float(np.min(support_polygon[:, 0])),
                float(np.max(support_polygon[:, 0])),
                float(np.min(support_polygon[:, 1])),
                float(np.max(support_polygon[:, 1])),
            )
            placement = FragmentPlacement(
                m=point.i,
                n=point.j,
                point=point,
                group_id=0,
                region=region,
                support=support,
                layer=layer,
                outside_x=region == "cladding",
                outside_y=region == "cladding",
                shift=shift,
                polygons=polygons,
            )
            placements.append(placement)
            records.extend(
                {
                    "region": region,
                    "point": point,
                    "group_id": 0,
                    "layer": layer,
                    "parent_layer": point.shell + 1,
                    "cladding_shift_layer": (
                        layer if region == "cladding" else None
                    ),
                    "polygon": polygon,
                    "modulation_shift": shift,
                }
                for polygon in polygons
            )
    all_points = [*cavity_points, *cladding_points]
    boundary_a = _hex_boundary_size(
        all_points,
        spec.period,
        spec.empty_buffer_periods * spec.period,
    )
    footprint = FiniteFootprint(
        shape="hex",
        width=boundary_a,
        height=math.sqrt(3.0) * boundary_a / 2.0,
        vertices=_hex_boundary_vertices(boundary_a),
    )
    validity_radius = spec.cavity_layers - 1 + spec.cladding_layers // 2
    validity_boundary = boundary_from_lattice_cells(
        lattice_points(validity_radius, spec.period), spec.period
    )
    metadata = {
        "finite_geometry": "hex",
        "finite_geometry_schema": "hex_full_cells_v1",
        "finite_geometry_shift_formula": (
            "hex_radial_ellipse_v1"
            if geometry.kind == ELLIPSE_SHIFT_GEOMETRY
            else "hex_normal_fan_v1"
        ),
        "cavity_layers": spec.cavity_layers,
        "cladding_layers": spec.cladding_layers,
        "layer_semantics": {
            "parent_layer": "one_based_from_center",
            "cladding_shift_layer": "one_based_from_cavity_boundary",
        },
        "bulk_points": cavity_points,
        "cladding_points": cladding_points,
        "bulk_cells": len(cavity_points),
        "cladding_cells": len(cladding_points),
        "bulk_holes": len(cavity_points) * len(cavity_holes),
        "cladding_holes": len(cladding_points) * len(cladding_holes),
        "total_holes": len(records),
        "inner_boundary": boundary_from_lattice_cells(cavity_points, spec.period),
        "finite_boundary": footprint.vertices,
        "finite_hex_a": boundary_a,
        "finite_Lx": footprint.width,
        "finite_Ly": footprint.height,
    }
    return FiniteGeometryPlan(
        spec=spec,
        placements=tuple(placements),
        records=tuple(records),
        cavity_analysis_cells=tuple(cavity_points),
        cladding_parent_cells=tuple(cladding_points),
        footprint=footprint,
        validity_boundary=validity_boundary,
        metadata=metadata,
    )


def _validate_square_hole_groups(holes: Sequence[np.ndarray], label: str) -> None:
    if len(holes) != 6:
        raise ValueError(f"{label} square template must contain exactly six holes")
    centers_x = np.array([np.mean(np.asarray(hole, dtype=float)[:, 0]) for hole in holes])
    if not np.all(centers_x[list(SQUARE_GROUP_HOLE_INDICES[1])] < -_EPSILON):
        raise ValueError(f"{label} square Group 1 holes must lie left of the cell center")
    if not np.all(centers_x[list(SQUARE_GROUP_HOLE_INDICES[2])] > _EPSILON):
        raise ValueError(f"{label} square Group 2 holes must lie right of the cell center")


def _support_inside(
    support: tuple[float, float, float, float],
    half_width: float,
    half_height: float,
) -> bool:
    x_min, x_max, y_min, y_max = support
    return (
        x_min >= -half_width - _EPSILON
        and x_max <= half_width + _EPSILON
        and y_min >= -half_height - _EPSILON
        and y_max <= half_height + _EPSILON
    )


def _square_layer(point: LatticePoint, spec: FiniteGeometrySpec) -> int:
    row_height = math.sqrt(3.0) * spec.period / 2.0
    ux = abs(point.x) / spec.period
    uy = abs(point.y) / row_height
    layer_x = math.ceil(max(0.0, ux - spec.cavity_layers) - 1e-12)
    layer_y = math.ceil(max(0.0, uy - spec.cavity_layers) - 1e-12)
    layer = max(1, layer_x, layer_y)
    if not 1 <= layer <= spec.cladding_layers:
        raise ValueError(
            f"square cladding parent ({point.i}, {point.j}) resolved to invalid layer {layer}"
        )
    return layer


def _axis_sign(primary: float, fallback: float) -> float:
    anchor = primary if abs(primary) > _EPSILON else fallback
    return 0.0 if abs(anchor) <= _EPSILON else math.copysign(1.0, anchor)


def _square_shift(
    point: LatticePoint,
    support: tuple[float, float, float, float],
    layer: int,
    outside_x: bool,
    outside_y: bool,
    spec: FiniteGeometrySpec,
    geometry: CladdingShiftGeometry,
    profile: CladdingShiftProfile,
    x_factor: float,
    y_factor: float,
) -> np.ndarray:
    envelope = profile.envelope(layer)
    if geometry.kind == ELLIPSE_SHIFT_GEOMETRY:
        return np.asarray(
            ellipse_inward_shift(
                point.x,
                point.y,
                spec.period,
                x_factor,
                y_factor,
                envelope=envelope,
            ),
            dtype=float,
        )
    x_min, x_max, y_min, y_max = support
    support_x = (x_min + x_max) / 2.0
    support_y = (y_min + y_max) / 2.0
    return np.array(
        [
            -_axis_sign(point.x, support_x) * spec.period * x_factor * envelope
            if outside_x
            else 0.0,
            -_axis_sign(point.y, support_y) * spec.period * y_factor * envelope
            if outside_y
            else 0.0,
        ],
        dtype=float,
    )


def _build_square_plan(
    spec: FiniteGeometrySpec,
    cavity_holes: Sequence[np.ndarray],
    cladding_holes: Sequence[np.ndarray],
    geometry: CladdingShiftGeometry,
    profile: CladdingShiftProfile,
    x_factor: float,
    y_factor: float,
) -> FiniteGeometryPlan:
    _validate_square_hole_groups(cavity_holes, "cavity")
    _validate_square_hole_groups(cladding_holes, "cladding")
    period = spec.period
    row_height = math.sqrt(3.0) * period / 2.0
    cavity_half_width = spec.cavity_layers * period
    cavity_half_height = (spec.cavity_layers + 0.5) * row_height
    total_layers = spec.cavity_layers + spec.cladding_layers
    structured_half_width = total_layers * period
    structured_half_height = (total_layers + 0.5) * row_height

    placements: list[FragmentPlacement] = []
    records: list[dict[str, object]] = []
    parent_regions: dict[tuple[int, int], list[str]] = {}
    parent_points: dict[tuple[int, int], LatticePoint] = {}
    for n in range(-total_layers, total_layers + 1):
        for m in range(-2 * total_layers, 2 * total_layers + 1):
            x = period * (m + 0.5 * n)
            y = row_height * n
            point = LatticePoint(
                i=m,
                j=n,
                shell=int(index_shell(m, n)),
                x=float(x),
                y=float(y),
            )
            for group_id, (x_lo_offset, x_hi_offset) in (
                (1, (-period / 2.0, 0.0)),
                (2, (0.0, period / 2.0)),
            ):
                support = (
                    x + x_lo_offset,
                    x + x_hi_offset,
                    y - row_height / 2.0,
                    y + row_height / 2.0,
                )
                if not _support_inside(
                    support, structured_half_width, structured_half_height
                ):
                    continue
                in_cavity = _support_inside(
                    support, cavity_half_width, cavity_half_height
                )
                region = "bulk" if in_cavity else "cladding"
                outside_x = (
                    support[0] < -cavity_half_width - _EPSILON
                    or support[1] > cavity_half_width + _EPSILON
                )
                outside_y = (
                    support[2] < -cavity_half_height - _EPSILON
                    or support[3] > cavity_half_height + _EPSILON
                )
                layer = 0 if in_cavity else _square_layer(point, spec)
                shift = (
                    np.zeros(2, dtype=float)
                    if in_cavity
                    else _square_shift(
                        point,
                        support,
                        layer,
                        outside_x,
                        outside_y,
                        spec,
                        geometry,
                        profile,
                        x_factor,
                        y_factor,
                    )
                )
                templates = cavity_holes if in_cavity else cladding_holes
                polygons = tuple(
                    _translate_polygon(
                        templates[hole_idx],
                        (point.x + shift[0], point.y + shift[1]),
                    )
                    for hole_idx in SQUARE_GROUP_HOLE_INDICES[group_id]
                )
                placement = FragmentPlacement(
                    m=m,
                    n=n,
                    point=point,
                    group_id=group_id,
                    region=region,
                    support=support,
                    layer=layer,
                    outside_x=outside_x,
                    outside_y=outside_y,
                    shift=shift,
                    polygons=polygons,
                )
                placements.append(placement)
                parent_regions.setdefault((m, n), []).append(region)
                parent_points[(m, n)] = point
                records.extend(
                    {
                        "region": region,
                        "point": point,
                        "group_id": group_id,
                        "layer": layer,
                        "support": support,
                        "outside_x": outside_x,
                        "outside_y": outside_y,
                        "polygon": polygon,
                        "modulation_shift": shift,
                    }
                    for polygon in polygons
                )

    placements_by_parent: dict[tuple[int, int], list[FragmentPlacement]] = {}
    for placement in placements:
        placements_by_parent.setdefault((placement.m, placement.n), []).append(
            placement
        )
    normalized_placements: list[FragmentPlacement] = []
    normalized_records: list[dict[str, object]] = []
    for placement in placements:
        siblings = placements_by_parent[(placement.m, placement.n)]
        is_full_cladding_parent = (
            len(siblings) == 2
            and all(sibling.region == "cladding" for sibling in siblings)
        )
        outside_x = (
            any(sibling.outside_x for sibling in siblings)
            if is_full_cladding_parent
            else placement.outside_x
        )
        outside_y = (
            any(sibling.outside_y for sibling in siblings)
            if is_full_cladding_parent
            else placement.outside_y
        )
        shift = placement.shift
        if placement.region == "cladding":
            shift = _square_shift(
                placement.point,
                placement.support,
                placement.layer,
                outside_x,
                outside_y,
                spec,
                geometry,
                profile,
                x_factor,
                y_factor,
            )
        templates = cavity_holes if placement.region == "bulk" else cladding_holes
        polygons = tuple(
            _translate_polygon(
                templates[hole_idx],
                (
                    placement.point.x + shift[0],
                    placement.point.y + shift[1],
                ),
            )
            for hole_idx in SQUARE_GROUP_HOLE_INDICES[placement.group_id]
        )
        normalized = FragmentPlacement(
            m=placement.m,
            n=placement.n,
            point=placement.point,
            group_id=placement.group_id,
            region=placement.region,
            support=placement.support,
            layer=placement.layer,
            outside_x=outside_x,
            outside_y=outside_y,
            shift=shift,
            polygons=polygons,
        )
        normalized_placements.append(normalized)
        normalized_records.extend(
            {
                "region": normalized.region,
                "point": normalized.point,
                "group_id": normalized.group_id,
                "layer": normalized.layer,
                "support": normalized.support,
                "outside_x": normalized.outside_x,
                "outside_y": normalized.outside_y,
                "polygon": polygon,
                "modulation_shift": normalized.shift,
            }
            for polygon in normalized.polygons
        )
    placements = normalized_placements
    records = normalized_records

    cavity_analysis_cells = tuple(
        parent_points[key]
        for key, regions in parent_regions.items()
        if regions == ["bulk", "bulk"]
    )
    cladding_parent_cells = tuple(
        parent_points[key]
        for key, regions in parent_regions.items()
        if "cladding" in regions
    )
    buffer = spec.empty_buffer_periods * period
    width = 2.0 * (structured_half_width + buffer)
    height = 2.0 * (structured_half_height + buffer)
    vertices = np.array(
        [
            [width / 2.0, height / 2.0],
            [-width / 2.0, height / 2.0],
            [-width / 2.0, -height / 2.0],
            [width / 2.0, -height / 2.0],
        ],
        dtype=float,
    )
    footprint = FiniteFootprint(
        shape="square",
        width=width,
        height=height,
        vertices=vertices,
        cavity_half_width=cavity_half_width,
        cavity_half_height=cavity_half_height,
        structured_half_width=structured_half_width,
        structured_half_height=structured_half_height,
    )
    validity_margin = spec.cladding_layers // 2
    validity_half_width = (spec.cavity_layers + validity_margin) * period
    validity_half_height = (
        spec.cavity_layers + validity_margin + 0.5
    ) * row_height
    validity_boundary = np.array(
        [
            [validity_half_width, validity_half_height],
            [-validity_half_width, validity_half_height],
            [-validity_half_width, -validity_half_height],
            [validity_half_width, -validity_half_height],
        ],
        dtype=float,
    )
    cavity_fragments = sum(p.region == "bulk" for p in placements)
    mixed_parents = sum(
        len(regions) == 2 and set(regions) == {"bulk", "cladding"}
        for regions in parent_regions.values()
    )
    outer_half_parents = sum(len(regions) == 1 for regions in parent_regions.values())
    metadata = {
        "finite_geometry": "square",
        "finite_geometry_schema": "square_fragments_v1",
        "finite_geometry_shift_formula": (
            "square_radial_ellipse_v1"
            if geometry.kind == ELLIPSE_SHIFT_GEOMETRY
            else "square_axis_components_v1"
        ),
        "lattice_vectors_um": [
            [period, 0.0],
            [period / 2.0, row_height],
        ],
        "organizational_cell_size_um": [period, row_height],
        "group_hole_indices": {
            str(key): list(value) for key, value in SQUARE_GROUP_HOLE_INDICES.items()
        },
        "cavity_fragment_bounds_um": [
            -cavity_half_width,
            cavity_half_width,
            -cavity_half_height,
            cavity_half_height,
        ],
        "structured_fragment_bounds_um": [
            -structured_half_width,
            structured_half_width,
            -structured_half_height,
            structured_half_height,
        ],
        "physical_footprint_um": [width, height],
        "outer_empty_pad_periods": spec.empty_buffer_periods,
        "fragment_counts": {
            "total": len(placements),
            "cavity": cavity_fragments,
            "cladding": len(placements) - cavity_fragments,
        },
        "parent_cell_counts": {
            "total": len(parent_regions),
            "full_cavity": len(cavity_analysis_cells),
            "mixed_cavity_cladding": mixed_parents,
            "outer_half_cell": outer_half_parents,
        },
        "square_layer_formula": "square_parent_chebyshev_v1",
        "analysis_cell_policy": "both_fragments_are_cavity",
        "cavity_layers": spec.cavity_layers,
        "cladding_layers": spec.cladding_layers,
        "bulk_points": list(cavity_analysis_cells),
        "cladding_points": list(cladding_parent_cells),
        "bulk_cells": len(cavity_analysis_cells),
        "cladding_cells": len(cladding_parent_cells),
        "bulk_holes": cavity_fragments * 3,
        "cladding_holes": (len(placements) - cavity_fragments) * 3,
        "total_holes": len(records),
        "inner_boundary": np.array(
            [
                [cavity_half_width, cavity_half_height],
                [-cavity_half_width, cavity_half_height],
                [-cavity_half_width, -cavity_half_height],
                [cavity_half_width, -cavity_half_height],
            ],
            dtype=float,
        ),
        "finite_boundary": vertices,
        "finite_Lx": width,
        "finite_Ly": height,
    }
    return FiniteGeometryPlan(
        spec=spec,
        placements=tuple(placements),
        records=tuple(records),
        cavity_analysis_cells=cavity_analysis_cells,
        cladding_parent_cells=cladding_parent_cells,
        footprint=footprint,
        validity_boundary=validity_boundary,
        metadata=metadata,
    )


def build_finite_geometry_plan(
    spec: FiniteGeometrySpec,
    cavity_holes: Sequence[np.ndarray],
    cladding_holes: Sequence[np.ndarray],
    *,
    shift_geometry: CladdingShiftGeometry,
    shift_profile: CladdingShiftProfile,
    x_shift_factor: float,
    y_shift_factor: float,
) -> FiniteGeometryPlan:
    if not isinstance(spec, FiniteGeometrySpec):
        raise TypeError("spec must be a FiniteGeometrySpec")
    x_factor = float(x_shift_factor)
    y_factor = float(y_shift_factor)
    if not math.isfinite(x_factor) or x_factor < 0.0:
        raise ValueError("x_shift_factor must be finite and non-negative")
    if not math.isfinite(y_factor) or y_factor < 0.0:
        raise ValueError("y_shift_factor must be finite and non-negative")
    if spec.shape == "hex":
        return _build_hex_plan(
            spec,
            cavity_holes,
            cladding_holes,
            shift_geometry,
            shift_profile,
            x_factor,
            y_factor,
        )
    return _build_square_plan(
        spec,
        cavity_holes,
        cladding_holes,
        shift_geometry,
        shift_profile,
        x_factor,
        y_factor,
    )
