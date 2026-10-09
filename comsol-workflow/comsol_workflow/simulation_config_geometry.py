"""Read explicit geometry from a serialized COMSOL simulation config.

The source format contains both geometry and solver settings.  This module
intentionally exposes only the geometry portion: the footprint and the hole
polygons from one selected layer.  Materials, physics, boundary conditions,
frequency shifts, and mode counts remain the responsibility of the active
workflow configuration.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


class SimulationConfigGeometryError(ValueError):
    """Raised when a simulation config cannot provide valid geometry."""


@dataclass(frozen=True, slots=True)
class SimulationConfigGeometry:
    """Geometry extracted from one ``*_simulation_config.json`` file.

    ``footprint`` and every item in ``holes`` are counter-clockwise convex
    polygons in micrometres.  Arrays are marked read-only by the loader so a
    caller cannot accidentally mutate the source geometry in place.
    """

    source_path: Path
    source_sha256: str
    layer_index: int
    footprint: np.ndarray
    holes: tuple[np.ndarray, ...]

    @property
    def hole_count(self) -> int:
        return len(self.holes)

    @property
    def footprint_bounds(self) -> tuple[float, float, float, float]:
        """Return ``(xmin, xmax, ymin, ymax)`` in micrometres."""

        return (
            float(np.min(self.footprint[:, 0])),
            float(np.max(self.footprint[:, 0])),
            float(np.min(self.footprint[:, 1])),
            float(np.max(self.footprint[:, 1])),
        )

    @property
    def footprint_width(self) -> float:
        xmin, xmax, _, _ = self.footprint_bounds
        return xmax - xmin

    @property
    def footprint_height(self) -> float:
        _, _, ymin, ymax = self.footprint_bounds
        return ymax - ymin

    @property
    def is_axis_aligned_rectangle(self) -> bool:
        """Whether the footprint is the four-corner rectangle used by COMSOL."""

        if self.footprint.shape != (4, 2):
            return False
        xmin, xmax, ymin, ymax = self.footprint_bounds
        expected = {
            (xmin, ymin),
            (xmax, ymin),
            (xmax, ymax),
            (xmin, ymax),
        }
        actual = {tuple(point) for point in self.footprint}
        return actual == expected

    @property
    def is_first_quadrant(self) -> bool:
        """Whether the footprint is anchored at the first-quadrant origin."""

        xmin, _, ymin, _ = self.footprint_bounds
        tolerance = 1e-12
        return abs(xmin) <= tolerance and abs(ymin) <= tolerance

    def summary(self) -> dict[str, Any]:
        """Return JSON-friendly provenance and geometry statistics."""

        xmin, xmax, ymin, ymax = self.footprint_bounds
        return {
            "source_path": str(self.source_path),
            "source_sha256": self.source_sha256,
            "layer_index": self.layer_index,
            "length_unit": "um",
            "footprint_vertex_count": int(len(self.footprint)),
            "footprint_bounds_um": {
                "xmin": xmin,
                "xmax": xmax,
                "ymin": ymin,
                "ymax": ymax,
            },
            "footprint_width_um": self.footprint_width,
            "footprint_height_um": self.footprint_height,
            "footprint_kind": (
                "axis_aligned_rectangle"
                if self.is_axis_aligned_rectangle
                else "convex_polygon"
            ),
            "first_quadrant": self.is_first_quadrant,
            "hole_count": self.hole_count,
            "hole_vertex_counts": sorted({int(len(hole)) for hole in self.holes}),
        }


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _coordinate_array(value: object, label: str, *, min_vertices: int) -> np.ndarray:
    if not isinstance(value, list) or len(value) < min_vertices:
        raise SimulationConfigGeometryError(
            f"{label} must be a list with at least {min_vertices} vertices"
        )

    coordinates: list[list[float]] = []
    for vertex_index, vertex in enumerate(value):
        if not isinstance(vertex, list) or len(vertex) != 2:
            raise SimulationConfigGeometryError(
                f"{label}[{vertex_index}] must contain exactly two coordinates"
            )
        if not all(_is_number(component) for component in vertex):
            raise SimulationConfigGeometryError(
                f"{label}[{vertex_index}] coordinates must be finite numbers"
            )
        coordinates.append([float(vertex[0]), float(vertex[1])])

    array = np.asarray(coordinates, dtype=float)
    if not np.all(np.isfinite(array)):
        raise SimulationConfigGeometryError(f"{label} must contain finite coordinates")
    if len(np.unique(array, axis=0)) != len(array):
        raise SimulationConfigGeometryError(f"{label} contains duplicate vertices")
    return array


def _signed_area(polygon: np.ndarray) -> float:
    return 0.5 * float(
        np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
        - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
    )


def _cross_2d(left: np.ndarray, right: np.ndarray) -> float:
    return float(left[0] * right[1] - left[1] * right[0])


def _validate_convex_ccw(polygon: np.ndarray, label: str) -> None:
    tolerance = 1e-12
    area = _signed_area(polygon)
    if area <= tolerance:
        raise SimulationConfigGeometryError(
            f"{label} must be counter-clockwise and non-degenerate; signed area={area:g}"
        )

    edge_crosses = []
    for index in range(len(polygon)):
        a = polygon[index]
        b = polygon[(index + 1) % len(polygon)]
        c = polygon[(index + 2) % len(polygon)]
        edge_crosses.append(_cross_2d(b - a, c - b))
    if min(edge_crosses) <= tolerance:
        raise SimulationConfigGeometryError(
            f"{label} must be strictly convex with no collinear consecutive vertices"
        )


def _inside_convex_polygon(point: np.ndarray, polygon: np.ndarray) -> bool:
    crosses = []
    for index in range(len(polygon)):
        edge = polygon[(index + 1) % len(polygon)] - polygon[index]
        offset = point - polygon[index]
        crosses.append(_cross_2d(edge, offset))
    return min(crosses) >= -1e-10


def _require_mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise SimulationConfigGeometryError(f"{label} must be an object")
    return value


def _read_json(path: Path) -> tuple[Mapping[str, object], str]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise SimulationConfigGeometryError(
            f"cannot read simulation config {path}: {exc}"
        ) from exc
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SimulationConfigGeometryError(
            f"simulation config must be UTF-8 JSON: {path}"
        ) from exc
    return _require_mapping(payload, "simulation config root"), hashlib.sha256(raw).hexdigest()


def load_simulation_config_geometry(
    path: str | Path,
    *,
    layer_index: int = 0,
) -> SimulationConfigGeometry:
    """Load explicit footprint and hole polygons from a simulation config.

    Only ``length_unit``, ``footprint`` and ``layers[layer_index].holes`` are
    consumed.  The source file's material and solver settings are deliberately
    ignored.  The returned coordinates remain in the source's declared
    micrometre system and are suitable for ``LayerSpec.holes``.
    """

    source_path = Path(path).expanduser()
    try:
        source_path = source_path.resolve(strict=True)
    except OSError as exc:
        raise SimulationConfigGeometryError(
            f"simulation config does not exist: {path}"
        ) from exc
    if not source_path.is_file():
        raise SimulationConfigGeometryError(f"simulation config is not a file: {path}")

    payload, source_sha256 = _read_json(source_path)
    if payload.get("length_unit") != "um":
        raise SimulationConfigGeometryError(
            "simulation config length_unit must be 'um'"
        )

    footprint = _coordinate_array(
        payload.get("footprint"),
        "footprint",
        min_vertices=3,
    )
    _validate_convex_ccw(footprint, "footprint")

    if type(layer_index) is not int or layer_index < 0:
        raise SimulationConfigGeometryError("layer_index must be a non-negative integer")
    raw_layers = payload.get("layers")
    if not isinstance(raw_layers, list):
        raise SimulationConfigGeometryError("layers must be a list")
    if layer_index >= len(raw_layers):
        raise SimulationConfigGeometryError(
            f"layer_index {layer_index} is out of range for {len(raw_layers)} layers"
        )
    layer = _require_mapping(raw_layers[layer_index], f"layers[{layer_index}]")
    raw_holes = layer.get("holes")
    if not isinstance(raw_holes, list):
        raise SimulationConfigGeometryError(
            f"layers[{layer_index}].holes must be a list"
        )

    holes: list[np.ndarray] = []
    for hole_index, raw_hole in enumerate(raw_holes):
        hole = _coordinate_array(
            raw_hole,
            f"layers[{layer_index}].holes[{hole_index}]",
            min_vertices=3,
        )
        _validate_convex_ccw(hole, f"layers[{layer_index}].holes[{hole_index}]")
        if not all(_inside_convex_polygon(vertex, footprint) for vertex in hole):
            raise SimulationConfigGeometryError(
                f"layers[{layer_index}].holes[{hole_index}] lies outside footprint"
            )
        hole.setflags(write=False)
        holes.append(hole)

    footprint.setflags(write=False)
    return SimulationConfigGeometry(
        source_path=source_path,
        source_sha256=source_sha256,
        layer_index=layer_index,
        footprint=footprint,
        holes=tuple(holes),
    )


__all__ = [
    "SimulationConfigGeometry",
    "SimulationConfigGeometryError",
    "load_simulation_config_geometry",
]
