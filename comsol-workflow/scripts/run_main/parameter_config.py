"""Load and safely update the shared parameters for primary workflows."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
from typing import Any
from uuid import uuid4

from comsol_workflow.cladding_shift_profile import (
    CladdingShiftGeometry,
    CladdingShiftProfile,
)


SCRIPTS_DIR = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = SCRIPTS_DIR / ".out"
UNIT_CELL_OUTPUT_ROOT = OUTPUT_ROOT / "unit_cell_band"
UNIT_CELL_2D_OUTPUT_ROOT = OUTPUT_ROOT / "unit_cell_2D"
STRIP_1D_OUTPUT_ROOT = OUTPUT_ROOT / "strip_1d"
FINITE_CAVITY_OUTPUT_ROOT = OUTPUT_ROOT / "finite_cavity"
WORKFLOW_OUTPUT_ROOTS = (
    UNIT_CELL_OUTPUT_ROOT,
    UNIT_CELL_2D_OUTPUT_ROOT,
    STRIP_1D_OUTPUT_ROOT,
    FINITE_CAVITY_OUTPUT_ROOT,
)
DEFAULT_PARAMETER_PATH = SCRIPTS_DIR / "parameter.json"
PARAMETER_PATH_ENV = "COMSOL_WORKFLOW_PARAMETER_PATH"


def parameter_path_from_environment() -> Path:
    """Return an optional task-local parameter file without changing defaults."""
    configured = os.environ.get(PARAMETER_PATH_ENV)
    if configured is None or not configured.strip():
        return DEFAULT_PARAMETER_PATH
    return Path(configured).expanduser().resolve()


PARAMETER_PATH = parameter_path_from_environment()
PARAMETER_SCHEMA_VERSION = 2
# Map JSON sections to the existing runtime contract; cache metadata stays stable.
PARAMETER_SECTIONS = {
    "cells": {"cavity": "cavity", "cladding": "cladding"},
    "structure": {name: name for name in (
        "cavity_layers", "cladding_layers", "center_frequency_thz",
        "center_frequency_source", "cladding_shift_factors", "cladding_shift_profile",
    )},
    "unit_cell_band": {"eigenmode_count": "unit_cell_eigenmode_count"},
    "strip_1d": {
        "eigenmode_count": "strip_eigenmode_count",
        "cut_angle_deg": "strip_cut_angle_deg",
        "kx": "strip_kx",
    },
    "finite_cavity": {
        "geometry": "finite_geometry",
        "eigenmode_count": "finite_eigenmode_count",
        "quarter_symmetry_ids": "finite_quarter_symmetry_ids",
        **{name: name for name in (
            "cladding_y_over_x_shift_ratio", "cladding_x_shift_factor",
            "cladding_y_shift_factor", "cladding_shift_geometry",
        )},
    },
}
FOURIER_B0_NM = 230.0
DEFAULT_UNIT_CELL_2D_Q_MAX_OVER_G = 0.1
DEFAULT_UNIT_CELL_2D_Q_POINTS_PER_AXIS = 19
DEFAULT_VALLEY_NORMAL_HALF_WIDTH_OVER_G = 0.003125
DEFAULT_VALLEY_NORMAL_POINTS_PER_ROUND = 3
DEFAULT_VALLEY_ANCHOR_STRIDE = 3
DEFAULT_VALLEY_TANGENT_HALF_WIDTH_OVER_G = 0.0125
DEFAULT_VALLEY_TANGENT_POINTS_PER_ROUND = 5
DEFAULT_VALLEY_LONGITUDINAL_MAX_GAP_OVER_G = 0.00625
DEFAULT_VALLEY_REFINEMENT_ROUNDS = 2
DEFAULT_VALLEY_REFINEMENT_RATIO = 4.0
DEFAULT_VALLEY_JUNCTION_HALF_WIDTH_OVER_G = 0.0125
DEFAULT_VALLEY_JUNCTION_POINTS_PER_AXIS = 7
DEFAULT_VALLEY_MAX_EDGE_RECENTERS = 2
DEFAULT_VALLEY_MAXIMUM_NEW_POINTS = 1000
DEFAULT_VALLEY_DISPLAY_POINTS_PER_AXIS = 577
UNIT_CELL_2D_ANGLE_METHODS = {
    "real", "imag", "complex-axis", "stokes-axis"
}
UNIT_CELL_2D_BAND_LABELS = {"p1", "p2", "d1", "d2"}
UNIT_CELL_2D_K_STRING_PRECISION = 6


def _require_number(
    value: object,
    path: str,
    *,
    positive: bool = False,
    non_negative: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{path} must be a number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{path} must be finite")
    if positive and numeric <= 0.0:
        raise ValueError(f"{path} must be positive")
    if non_negative and numeric < 0.0:
        raise ValueError(f"{path} must be non-negative")
    return numeric


def _format_decimal(value: float, decimal_places: int = 3) -> str:
    return f"{float(value):.{decimal_places}f}".rstrip("0").rstrip(".")


@dataclass(frozen=True)
class CellParameters:
    b0_nm: float
    eta: float
    zeta: float

    @classmethod
    def from_mapping(cls, value: object, path: str) -> "CellParameters":
        if not isinstance(value, dict):
            raise TypeError(f"{path} must be an object")
        required = {"b0_nm", "eta", "zeta"}
        missing = sorted(required.difference(value))
        if missing:
            raise ValueError(f"{path} is missing: {', '.join(missing)}")
        return cls(
            b0_nm=_require_number(value["b0_nm"], f"{path}.b0_nm", positive=True),
            eta=_require_number(value["eta"], f"{path}.eta", positive=True),
            zeta=_require_number(value["zeta"], f"{path}.zeta", positive=True),
        )

    def to_fourier_params(self, name: str) -> dict[str, float | str]:
        return {
            "name": str(name),
            "r_f0": self.eta,
            "b_square_f0": (self.b0_nm / FOURIER_B0_NM) ** 2,
            "b_square_f3": (self.zeta**2 - 1.0) / 2.0,
        }

    def to_dict(self) -> dict[str, float]:
        return {
            "b0_nm": self.b0_nm,
            "eta": self.eta,
            "zeta": self.zeta,
        }

    def label(self, prefix: str) -> str:
        return (
            f"{prefix}({_format_decimal(self.b0_nm)}"
            f"-{_format_decimal(self.eta)}"
            f"-{_format_decimal(self.zeta)})"
        )


@dataclass(frozen=True)
class UnitCell2DValleyRefinementParameters:
    normal_half_width_over_g: float
    normal_points_per_round: int
    anchor_stride: int
    tangent_half_width_over_g: float
    tangent_points_per_round: int
    longitudinal_max_gap_over_g: float
    refinement_rounds: int
    refinement_ratio: float
    junction_half_width_over_g: float
    junction_points_per_axis: int
    max_edge_recenters: int
    maximum_new_points: int
    display_points_per_axis: int

    def to_dict(self) -> dict[str, object]:
        return {
            "normal_half_width_over_G": self.normal_half_width_over_g,
            "normal_points_per_round": self.normal_points_per_round,
            "anchor_stride": self.anchor_stride,
            "tangent_half_width_over_G": self.tangent_half_width_over_g,
            "tangent_points_per_round": self.tangent_points_per_round,
            "longitudinal_max_gap_over_G": self.longitudinal_max_gap_over_g,
            "refinement_rounds": self.refinement_rounds,
            "refinement_ratio": self.refinement_ratio,
            "junction_half_width_over_G": self.junction_half_width_over_g,
            "junction_points_per_axis": self.junction_points_per_axis,
            "max_edge_recenters": self.max_edge_recenters,
            "maximum_new_points": self.maximum_new_points,
            "display_points_per_axis": self.display_points_per_axis,
        }


@dataclass(frozen=True)
class UnitCell2DParameters:
    cell: CellParameters
    eigenmode_pair_count: int
    q_max_over_g: float
    q_points_per_axis: int
    is_quarter: int
    target_bands: tuple[str, ...]
    frequency_tolerance_thz: float
    minimum_overlap: float
    minimum_ambiguity_gap: float
    normalization_kind: str
    field_plane: str
    mask_relative_threshold: float
    angle_methods: tuple[str, ...]
    analysis_enabled: bool
    valley_refinement: UnitCell2DValleyRefinementParameters

    @property
    def q_axis_over_g(self) -> tuple[float, ...]:
        if self.is_quarter:
            step = self.q_max_over_g / (self.q_points_per_axis - 1)
            return tuple(
                0.0 if index == 0 else index * step
                for index in range(self.q_points_per_axis)
            )
        step = 2.0 * self.q_max_over_g / (self.q_points_per_axis - 1)
        return tuple(
            0.0 if index == self.q_points_per_axis // 2
            else -self.q_max_over_g + index * step
            for index in range(self.q_points_per_axis)
        )

    @property
    def plot_q_axis_over_g(self) -> tuple[float, ...]:
        solved = self.q_axis_over_g
        if not self.is_quarter:
            return solved
        return tuple(-value for value in reversed(solved[1:])) + solved

    @property
    def sampling_domain(self) -> str:
        return "quarter" if self.is_quarter else "full"


@dataclass(frozen=True)
class SharedParameters:
    cavity: CellParameters
    cladding: CellParameters
    cavity_layers: int
    cladding_layers: int
    mesh_size: int
    unit_cell_eigenmode_count: int
    finite_eigenmode_count: int
    finite_geometry: str
    center_frequency_thz: float
    cladding_shift_factors: tuple[float, ...]
    cladding_shift_input_kind: str
    cladding_y_over_x_shift_ratio: float | None
    cladding_x_shift_factor: float | None
    cladding_y_shift_factor: float | None
    cladding_shift_geometry: CladdingShiftGeometry
    cladding_shift_profile: CladdingShiftProfile
    unit_cell_2d: UnitCell2DParameters
    source_path: Path
    strip_eigenmode_count: int = 2
    strip_cut_angle_deg: float = 60.0
    strip_kx: float = 0.0
    finite_quarter_symmetry_ids: tuple[int, ...] = (1,)

    @property
    def case_parameter_label(self) -> str:
        return (
            f"{self.cavity.label('cav')}_"
            f"{self.cladding.label('clad')}_mesh{self.mesh_size}"
        )

    @property
    def unit_cell_case_parameter_label(self) -> str:
        """Return the stable cell/mesh label; solver count lives in metadata."""
        return self.case_parameter_label

    @property
    def eigenfrequency_shift(self) -> str:
        return f"{_format_decimal(self.center_frequency_thz, 12)} [THz]"

    @property
    def finite_cladding_shift_factor_pairs(
        self,
    ) -> tuple[tuple[float, float], ...]:
        """Resolve finite x/y targets from ratio scan or one explicit pair."""
        if self.cladding_shift_input_kind == "xy":
            return (
                (
                    float(self.cladding_x_shift_factor),
                    float(self.cladding_y_shift_factor),
                ),
            )
        return tuple(
            (factor, factor * float(self.cladding_y_over_x_shift_ratio))
            for factor in self.cladding_shift_factors
        )

    @property
    def finite_cladding_shift_input(self) -> dict[str, float | str]:
        if self.cladding_shift_input_kind == "xy":
            return {
                "kind": "xy",
                "x_factor": float(self.cladding_x_shift_factor),
                "y_factor": float(self.cladding_y_shift_factor),
            }
        return {
            "kind": "scan_ratio",
            "y_over_x_ratio": float(self.cladding_y_over_x_shift_ratio),
        }

    def structure_series_label(self, prefix: str) -> str:
        """Name a finite-structure series without affecting unit-cell labels."""
        geometry_suffix = (
            self.cladding_shift_geometry.series_suffix
            if prefix in {"finite", "finite_quarter"}
            else ""
        )
        finite_geometry_suffix = (
            "_square"
            if prefix in {"finite", "finite_quarter"}
            and self.finite_geometry == "square"
            else ""
        )
        return (
            f"{prefix}_{self.cavity_layers}-{self.cladding_layers}_"
            f"{self.case_parameter_label}"
            f"{finite_geometry_suffix}"
            f"{self.cladding_shift_profile.series_suffix}"
            f"{geometry_suffix}"
        )

    @property
    def identity(self) -> dict[str, object]:
        return {
            "cavity": self.cavity.to_dict(),
            "cladding": self.cladding.to_dict(),
            "mesh_size": self.mesh_size,
            "unit_cell_eigenmode_count": self.unit_cell_eigenmode_count,
        }

    def to_metadata(self) -> dict[str, object]:
        return {
            **self.identity,
            "finite_eigenmode_count": self.finite_eigenmode_count,
            "finite_geometry": self.finite_geometry,
            "cavity_layers": self.cavity_layers,
            "cladding_layers": self.cladding_layers,
            "center_frequency_thz": self.center_frequency_thz,
            "cladding_shift_factors": list(self.cladding_shift_factors),
            "finite_cladding_shift_input": self.finite_cladding_shift_input,
            "finite_cladding_shift_factor_pairs": [
                {
                    "base_factor": x_factor,
                    "x_factor": x_factor,
                    "y_factor": y_factor,
                }
                for x_factor, y_factor in self.finite_cladding_shift_factor_pairs
            ],
            "cladding_shift_geometry": self.cladding_shift_geometry.to_dict(),
            "cladding_shift_profile": self.cladding_shift_profile.to_dict(),
            "cladding_shift_profile_label": self.cladding_shift_profile.label,
            "case_parameter_label": self.case_parameter_label,
            "unit_cell_case_parameter_label": self.unit_cell_case_parameter_label,
            "unit_cell_2d": unit_cell_2d_metadata(self.unit_cell_2d),
            "parameter_path": str(self.source_path),
        }


def _read_parameter_data(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Shared parameter file does not exist: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Shared parameter file is not valid JSON: {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise TypeError("parameter.json root must be an object")
    return data


def _check_keys(value: object, allowed: set[str], path: str) -> None:
    if not isinstance(value, dict):
        raise TypeError(f"{path} must be an object")
    unknown = sorted(set(value).difference(allowed))
    if unknown:
        raise ValueError(f"{path} has unknown fields: {', '.join(unknown)}")


def parameter_data_to_flat(data: dict[str, Any]) -> dict[str, Any]:
    """Resolve v1 snapshots or v2 sections without mutating the input mapping."""
    version = data.get("schema_version")
    if type(version) is not int or version not in {1, PARAMETER_SCHEMA_VERSION}:
        raise ValueError(
            f"schema_version must be 1 or {PARAMETER_SCHEMA_VERSION}; got {version!r}"
        )
    if version == 1:
        mixed = sorted(set(data).intersection(PARAMETER_SECTIONS))
        if mixed:
            raise ValueError(f"schema_version=1 cannot contain v2 sections: {mixed}")
        return dict(data)

    _check_keys(
        data, {"schema_version", "mesh_size", "unit_cell_2d", *PARAMETER_SECTIONS},
        "parameter.json",
    )
    flat = {"schema_version": version, "mesh_size": data.get("mesh_size")}
    for section, fields in PARAMETER_SECTIONS.items():
        value = data.get(section)
        _check_keys(value, set(fields), section)
        flat.update({target: value[key] for key, target in fields.items() if key in value})
    for cell in ("cavity", "cladding"):
        _check_keys(flat.get(cell), {"b0_nm", "eta", "zeta"}, f"cells.{cell}")
    settings = data.get("unit_cell_2d")
    _check_keys(settings, {
        "b0_nm", "eta", "zeta", "eigenmode_pair_count", "q_max_over_G",
        "q_points_per_axis", "isQuarter", "q_axis_over_G", "target_bands",
        "frequency_tolerance_thz", "minimum_overlap", "minimum_ambiguity_gap",
        "normalization_kind", "field_plane", "mask_relative_threshold",
        "angle_methods", "analysis_enabled", "valley_refinement",
    }, "unit_cell_2d")
    if "valley_refinement" in settings:
        _check_keys(settings["valley_refinement"], {
            "normal_half_width_over_G", "normal_points_per_round", "anchor_stride",
            "tangent_half_width_over_G", "tangent_points_per_round",
            "longitudinal_max_gap_over_G", "refinement_rounds", "refinement_ratio",
            "junction_half_width_over_G", "junction_points_per_axis",
            "max_edge_recenters", "maximum_new_points", "display_points_per_axis",
        }, "unit_cell_2d.valley_refinement")
    flat["unit_cell_2d"] = settings
    for key, target in PARAMETER_SECTIONS["strip_1d"].items():
        if target not in flat:
            raise ValueError(f"strip_1d is missing: {key}")
    if "finite_quarter_symmetry_ids" not in flat:
        raise ValueError("finite_cavity is missing: quarter_symmetry_ids")
    return flat


def _load_unit_cell_2d(
    value: object,
    *,
    legacy_fallback_cell: CellParameters,
) -> UnitCell2DParameters:
    if value is None:
        value = legacy_fallback_cell.to_dict()
    if not isinstance(value, dict):
        raise TypeError("unit_cell_2d must be an object")
    eigenmode_pair_count = value.get("eigenmode_pair_count", 2)
    if type(eigenmode_pair_count) is not int or eigenmode_pair_count <= 0:
        raise ValueError(
            "unit_cell_2d.eigenmode_pair_count must be a positive integer"
        )
    if "q_axis_over_G" in value:
        raise ValueError(
            "unit_cell_2d.q_axis_over_G is no longer accepted; use "
            "q_max_over_G and q_points_per_axis"
        )
    q_max_over_g = _require_number(
        value.get("q_max_over_G", DEFAULT_UNIT_CELL_2D_Q_MAX_OVER_G),
        "unit_cell_2d.q_max_over_G",
        positive=True,
    )
    q_points_per_axis = value.get(
        "q_points_per_axis", DEFAULT_UNIT_CELL_2D_Q_POINTS_PER_AXIS
    )
    is_quarter = value.get("isQuarter", 0)
    if type(is_quarter) is not int or is_quarter not in {0, 1}:
        raise ValueError("unit_cell_2d.isQuarter must be the integer 0 or 1")
    invalid_points = type(q_points_per_axis) is not int
    if is_quarter:
        invalid_points = invalid_points or q_points_per_axis < 2
        points_message = "an integer >= 2 when isQuarter=1"
    else:
        invalid_points = (
            invalid_points
            or q_points_per_axis < 3
            or q_points_per_axis % 2 == 0
        )
        points_message = "an odd integer >= 3 when isQuarter=0"
    if invalid_points:
        raise ValueError(
            f"unit_cell_2d.q_points_per_axis must be {points_message}"
        )
    if is_quarter:
        step = q_max_over_g / (q_points_per_axis - 1)
        axis = tuple(
            0.0 if index == 0 else index * step
            for index in range(q_points_per_axis)
        )
    else:
        step = 2.0 * q_max_over_g / (q_points_per_axis - 1)
        axis = tuple(
            0.0 if index == q_points_per_axis // 2
            else -q_max_over_g + index * step
            for index in range(q_points_per_axis)
        )
    formatted_axis = {
        f"{item:.{UNIT_CELL_2D_K_STRING_PRECISION}f}" for item in axis
    }
    if len(formatted_axis) != len(axis):
        raise ValueError(
            "unit_cell_2d uniform q-axis values must remain unique at "
            f"{UNIT_CELL_2D_K_STRING_PRECISION} decimal places"
        )
    raw_bands = value.get("target_bands", ["p2"])
    if not isinstance(raw_bands, list) or not raw_bands:
        raise ValueError("unit_cell_2d.target_bands must be a non-empty list")
    bands = tuple(str(item).strip() for item in raw_bands)
    if any(not item for item in bands) or len(set(bands)) != len(bands):
        raise ValueError("unit_cell_2d.target_bands must contain unique names")
    unknown_bands = sorted(set(bands).difference(UNIT_CELL_2D_BAND_LABELS))
    if unknown_bands:
        raise ValueError(f"Invalid unit_cell_2d.target_bands: {unknown_bands}")
    raw_methods = value.get(
        "angle_methods", ["real", "imag", "complex-axis", "stokes-axis"]
    )
    if not isinstance(raw_methods, list) or not raw_methods:
        raise ValueError("unit_cell_2d.angle_methods must be a non-empty list")
    methods = tuple(str(item) for item in raw_methods)
    unknown = sorted(set(methods).difference(UNIT_CELL_2D_ANGLE_METHODS))
    if unknown or len(set(methods)) != len(methods):
        raise ValueError(f"Invalid unit_cell_2d.angle_methods: {unknown or methods}")
    normalization = value.get("normalization_kind", "planar_l2")
    field_plane = value.get("field_plane", "air")
    enabled = value.get("analysis_enabled", True)
    if normalization != "planar_l2":
        raise ValueError("unit_cell_2d.normalization_kind must be 'planar_l2'")
    if field_plane != "air":
        raise ValueError("unit_cell_2d.field_plane must be 'air'")
    if type(enabled) is not bool:
        raise TypeError("unit_cell_2d.analysis_enabled must be a boolean")
    if not enabled:
        raise ValueError(
            "unit_cell_2d.analysis_enabled must remain true; "
            "use run_unit_cell_2d.py --skip-analysis for diagnostics"
        )
    overlap = _require_number(
        value.get("minimum_overlap", 0.2),
        "unit_cell_2d.minimum_overlap", non_negative=True,
    )
    gap = _require_number(
        value.get("minimum_ambiguity_gap", 0.05),
        "unit_cell_2d.minimum_ambiguity_gap", non_negative=True,
    )
    mask = _require_number(
        value.get("mask_relative_threshold", 1e-6),
        "unit_cell_2d.mask_relative_threshold", non_negative=True,
    )
    if overlap > 1 or gap > 1 or mask >= 1:
        raise ValueError("unit_cell_2d overlap/gap/mask thresholds are out of range")
    raw_refinement = value.get("valley_refinement", {})
    if not isinstance(raw_refinement, dict):
        raise TypeError("unit_cell_2d.valley_refinement must be an object")

    def refinement_odd_integer(name: str, default: int, minimum: int) -> int:
        item = raw_refinement.get(name, default)
        if type(item) is not int or item < minimum or item % 2 == 0:
            raise ValueError(
                f"unit_cell_2d.valley_refinement.{name} must be an odd "
                f"integer >= {minimum}"
            )
        return item

    def refinement_integer(name: str, default: int, minimum: int) -> int:
        item = raw_refinement.get(name, default)
        if type(item) is not int or item < minimum:
            raise ValueError(
                f"unit_cell_2d.valley_refinement.{name} must be an integer "
                f">= {minimum}"
            )
        return item

    refinement = UnitCell2DValleyRefinementParameters(
        normal_half_width_over_g=_require_number(
            raw_refinement.get(
                "normal_half_width_over_G",
                DEFAULT_VALLEY_NORMAL_HALF_WIDTH_OVER_G,
            ),
            "unit_cell_2d.valley_refinement.normal_half_width_over_G",
            positive=True,
        ),
        normal_points_per_round=refinement_odd_integer(
            "normal_points_per_round",
            DEFAULT_VALLEY_NORMAL_POINTS_PER_ROUND,
            3,
        ),
        anchor_stride=refinement_integer(
            "anchor_stride", DEFAULT_VALLEY_ANCHOR_STRIDE, 1
        ),
        tangent_half_width_over_g=_require_number(
            raw_refinement.get(
                "tangent_half_width_over_G",
                DEFAULT_VALLEY_TANGENT_HALF_WIDTH_OVER_G,
            ),
            "unit_cell_2d.valley_refinement.tangent_half_width_over_G",
            positive=True,
        ),
        tangent_points_per_round=refinement_odd_integer(
            "tangent_points_per_round",
            DEFAULT_VALLEY_TANGENT_POINTS_PER_ROUND,
            3,
        ),
        longitudinal_max_gap_over_g=_require_number(
            raw_refinement.get(
                "longitudinal_max_gap_over_G",
                DEFAULT_VALLEY_LONGITUDINAL_MAX_GAP_OVER_G,
            ),
            "unit_cell_2d.valley_refinement.longitudinal_max_gap_over_G",
            positive=True,
        ),
        refinement_rounds=refinement_integer(
            "refinement_rounds", DEFAULT_VALLEY_REFINEMENT_ROUNDS, 1
        ),
        refinement_ratio=_require_number(
            raw_refinement.get(
                "refinement_ratio", DEFAULT_VALLEY_REFINEMENT_RATIO
            ),
            "unit_cell_2d.valley_refinement.refinement_ratio",
            positive=True,
        ),
        junction_half_width_over_g=_require_number(
            raw_refinement.get(
                "junction_half_width_over_G",
                DEFAULT_VALLEY_JUNCTION_HALF_WIDTH_OVER_G,
            ),
            "unit_cell_2d.valley_refinement.junction_half_width_over_G",
            positive=True,
        ),
        junction_points_per_axis=refinement_odd_integer(
            "junction_points_per_axis",
            DEFAULT_VALLEY_JUNCTION_POINTS_PER_AXIS,
            3,
        ),
        max_edge_recenters=refinement_integer(
            "max_edge_recenters", DEFAULT_VALLEY_MAX_EDGE_RECENTERS, 0
        ),
        maximum_new_points=refinement_integer(
            "maximum_new_points", DEFAULT_VALLEY_MAXIMUM_NEW_POINTS, 1
        ),
        display_points_per_axis=refinement_odd_integer(
            "display_points_per_axis",
            DEFAULT_VALLEY_DISPLAY_POINTS_PER_AXIS,
            3,
        ),
    )
    if refinement.refinement_ratio <= 1.0:
        raise ValueError(
            "unit_cell_2d.valley_refinement.refinement_ratio must be > 1"
        )
    return UnitCell2DParameters(
        cell=CellParameters.from_mapping(value, "unit_cell_2d"),
        eigenmode_pair_count=eigenmode_pair_count,
        q_max_over_g=q_max_over_g,
        q_points_per_axis=q_points_per_axis,
        is_quarter=is_quarter,
        target_bands=bands,
        frequency_tolerance_thz=_require_number(
            value.get("frequency_tolerance_thz", 7.0),
            "unit_cell_2d.frequency_tolerance_thz", positive=True,
        ),
        minimum_overlap=overlap,
        minimum_ambiguity_gap=gap,
        normalization_kind=normalization,
        field_plane=field_plane,
        mask_relative_threshold=mask,
        angle_methods=methods,
        analysis_enabled=enabled,
        valley_refinement=refinement,
    )


def unit_cell_2d_metadata(value: UnitCell2DParameters) -> dict[str, object]:
    return {
        **value.cell.to_dict(),
        "eigenmode_pair_count": value.eigenmode_pair_count,
        "isQuarter": value.is_quarter,
        "sampling_domain": value.sampling_domain,
        "q_grid_kind": (
            "quarter_uniform" if value.is_quarter else "symmetric_uniform"
        ),
        "q_max_over_G": value.q_max_over_g,
        "q_points_per_axis": value.q_points_per_axis,
        "solved_q_axis_over_G": list(value.q_axis_over_g),
        "plot_q_axis_over_G": list(value.plot_q_axis_over_g),
        "solved_point_count": value.q_points_per_axis**2,
        "plot_coordinate_count": len(value.plot_q_axis_over_g) ** 2,
        "target_bands": list(value.target_bands),
        "frequency_tolerance_thz": value.frequency_tolerance_thz,
        "minimum_overlap": value.minimum_overlap,
        "minimum_ambiguity_gap": value.minimum_ambiguity_gap,
        "normalization_kind": value.normalization_kind,
        "field_plane": value.field_plane,
        "mask_relative_threshold": value.mask_relative_threshold,
        "angle_methods": list(value.angle_methods),
        "analysis_enabled": value.analysis_enabled,
        "valley_refinement": value.valley_refinement.to_dict(),
    }


def load_shared_parameters(path: str | os.PathLike[str] = PARAMETER_PATH) -> SharedParameters:
    parameter_path = Path(path).resolve()
    data = parameter_data_to_flat(_read_parameter_data(parameter_path))

    strip_eigenmode_count = data.get("strip_eigenmode_count", 2)
    if type(strip_eigenmode_count) is not int or strip_eigenmode_count <= 0:
        raise ValueError("strip_1d.eigenmode_count must be a positive integer")
    strip_cut_angle_deg = _require_number(
        data.get("strip_cut_angle_deg", 60.0), "strip_1d.cut_angle_deg"
    )
    if strip_cut_angle_deg not in {0.0, 60.0}:
        raise ValueError("strip_1d.cut_angle_deg must be 0 or 60")
    strip_kx = _require_number(data.get("strip_kx", 0.0), "strip_1d.kx")
    symmetry_ids = data.get("finite_quarter_symmetry_ids", [1])
    if (
        not isinstance(symmetry_ids, list) or not symmetry_ids
        or any(type(item) is not int or item not in {1, 2, 3, 4} for item in symmetry_ids)
        or len(set(symmetry_ids)) != len(symmetry_ids)
    ):
        raise ValueError(
            "finite_cavity.quarter_symmetry_ids must be unique integers from 1 to 4"
        )

    mesh_size = data.get("mesh_size")
    if type(mesh_size) is not int or not 1 <= mesh_size <= 9:
        raise ValueError("mesh_size must be an integer from 1 to 9")

    unit_cell_eigenmode_count = data.get("unit_cell_eigenmode_count")
    if type(unit_cell_eigenmode_count) is not int or unit_cell_eigenmode_count <= 0:
        raise ValueError("unit_cell_eigenmode_count must be a positive integer")

    finite_eigenmode_count = data.get("finite_eigenmode_count")
    if type(finite_eigenmode_count) is not int or finite_eigenmode_count <= 0:
        raise ValueError("finite_eigenmode_count must be a positive integer")

    finite_geometry = data.get("finite_geometry", "hex")
    if not isinstance(finite_geometry, str):
        raise TypeError("finite_geometry must be a string")
    if finite_geometry not in {"hex", "square"}:
        raise ValueError("finite_geometry must be 'hex' or 'square'")

    cavity_layers = data.get("cavity_layers")
    if type(cavity_layers) is not int or cavity_layers <= 0:
        raise ValueError("cavity_layers must be a positive integer")

    cladding_layers = data.get("cladding_layers")
    if type(cladding_layers) is not int or cladding_layers <= 0:
        raise ValueError("cladding_layers must be a positive integer")

    raw_shifts = data.get("cladding_shift_factors")
    if not isinstance(raw_shifts, list) or not raw_shifts:
        raise ValueError("cladding_shift_factors must be a non-empty list")
    shifts = tuple(
        _require_number(
            value,
            f"cladding_shift_factors[{idx}]",
            non_negative=True,
        )
        for idx, value in enumerate(raw_shifts)
    )
    has_ratio = "cladding_y_over_x_shift_ratio" in data
    has_x_factor = "cladding_x_shift_factor" in data
    has_y_factor = "cladding_y_shift_factor" in data
    if has_ratio and (has_x_factor or has_y_factor):
        raise ValueError(
            "cladding_y_over_x_shift_ratio conflicts with explicit "
            "cladding_x_shift_factor/cladding_y_shift_factor"
        )
    if has_x_factor != has_y_factor:
        raise ValueError(
            "cladding_x_shift_factor and cladding_y_shift_factor "
            "must be provided together"
        )

    if has_x_factor:
        cladding_shift_input_kind = "xy"
        cladding_y_over_x_shift_ratio = None
        cladding_x_shift_factor = _require_number(
            data["cladding_x_shift_factor"],
            "cladding_x_shift_factor",
            non_negative=True,
        )
        cladding_y_shift_factor = _require_number(
            data["cladding_y_shift_factor"],
            "cladding_y_shift_factor",
            non_negative=True,
        )
    else:
        cladding_shift_input_kind = "scan_ratio"
        cladding_x_shift_factor = None
        cladding_y_shift_factor = None
        raw_ratio = data.get("cladding_y_over_x_shift_ratio", 1.0)
        cladding_y_over_x_shift_ratio = _require_number(
            raw_ratio,
            "cladding_y_over_x_shift_ratio",
            non_negative=True,
        )

    cavity = CellParameters.from_mapping(data.get("cavity"), "cavity")
    return SharedParameters(
        cavity=cavity,
        cladding=CellParameters.from_mapping(data.get("cladding"), "cladding"),
        cavity_layers=cavity_layers,
        cladding_layers=cladding_layers,
        mesh_size=mesh_size,
        unit_cell_eigenmode_count=unit_cell_eigenmode_count,
        finite_eigenmode_count=finite_eigenmode_count,
        finite_geometry=finite_geometry,
        center_frequency_thz=_require_number(
            data.get("center_frequency_thz"),
            "center_frequency_thz",
            positive=True,
        ),
        cladding_shift_factors=shifts,
        cladding_shift_input_kind=cladding_shift_input_kind,
        cladding_y_over_x_shift_ratio=cladding_y_over_x_shift_ratio,
        cladding_x_shift_factor=cladding_x_shift_factor,
        cladding_y_shift_factor=cladding_y_shift_factor,
        cladding_shift_geometry=CladdingShiftGeometry.from_mapping(
            data.get("cladding_shift_geometry")
        ),
        cladding_shift_profile=CladdingShiftProfile.from_mapping(
            data.get("cladding_shift_profile")
        ),
        unit_cell_2d=_load_unit_cell_2d(
            data.get("unit_cell_2d"),
            legacy_fallback_cell=cavity,
        ),
        source_path=parameter_path,
        strip_eigenmode_count=strip_eigenmode_count,
        strip_cut_angle_deg=strip_cut_angle_deg,
        strip_kx=strip_kx,
        finite_quarter_symmetry_ids=tuple(symmetry_ids),
    )


def _identities_match(
    left: dict[str, object],
    right: dict[str, object],
    *,
    include_cladding: bool = True,
) -> bool:
    cell_names = ("cavity", "cladding") if include_cladding else ("cavity",)
    for cell_name in cell_names:
        left_cell = left[cell_name]
        right_cell = right[cell_name]
        if not isinstance(left_cell, dict) or not isinstance(right_cell, dict):
            return False
        for field in ("b0_nm", "eta", "zeta"):
            if not math.isclose(
                float(left_cell[field]),
                float(right_cell[field]),
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                return False
    return (
        int(left["mesh_size"]) == int(right["mesh_size"])
        and int(left["unit_cell_eigenmode_count"])
        == int(right["unit_cell_eigenmode_count"])
    )


def update_center_frequency_from_cavity_p2(
    frequency_thz: float,
    *,
    expected_parameters: SharedParameters,
    path: str | os.PathLike[str] = PARAMETER_PATH,
    workflow: str = "run_band_pair",
    include_cladding: bool = True,
) -> SharedParameters:
    """Update only center_frequency_thz after a complete cavity p2@Gamma run."""
    numeric_frequency = _require_number(
        frequency_thz,
        "cavity p2@Gamma frequency",
        positive=True,
    )
    parameter_path = Path(path).resolve()
    current = load_shared_parameters(parameter_path)
    if workflow not in {"run_band_pair", "run_band_solo"}:
        raise ValueError(f"Unsupported unit-cell workflow: {workflow}")
    if not _identities_match(
        current.identity,
        expected_parameters.identity,
        include_cladding=include_cladding,
    ):
        raise RuntimeError(
            "parameter.json geometry, mesh, or unit-cell mode count changed "
            "during the unit-cell run; "
            "refusing to overwrite center_frequency_thz"
        )

    data = _read_parameter_data(parameter_path)
    structure = data if data["schema_version"] == 1 else data["structure"]
    structure["center_frequency_thz"] = numeric_frequency
    structure["center_frequency_source"] = {
        "workflow": workflow,
        "cell": "cavity",
        "band": "p2",
        "k_point": "Gamma",
    }

    output_root = parameter_path.parent / ".out"
    unit_cell_output_root = output_root / UNIT_CELL_OUTPUT_ROOT.name
    temporary_dir = unit_cell_output_root / ".tmp"
    temporary_dir.mkdir(parents=True, exist_ok=True)
    temporary_path = temporary_dir / (
        f"parameter-{os.getpid()}-{uuid4().hex}.json"
    )
    try:
        temporary_path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary_path, parameter_path)
    finally:
        temporary_path.unlink(missing_ok=True)
        for directory in (temporary_dir, unit_cell_output_root, output_root):
            try:
                directory.rmdir()
            except OSError:
                pass
    return load_shared_parameters(parameter_path)


ACTIVE_PARAMETERS = load_shared_parameters()
