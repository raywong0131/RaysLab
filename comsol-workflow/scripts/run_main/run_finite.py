#!/usr/bin/env python3
"""Build finite-size cavity geometry from the shared cavity and cladding cells."""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis  # noqa: E402
from comsol_workflow.cladding_shift_profile import (  # noqa: E402
    ELLIPSE_SHIFT_GEOMETRY,
    GROUPS_SHIFT_GEOMETRY,
    CladdingShiftGeometry,
    CladdingShiftProfile,
    cladding_shift_group,
    ellipse_inward_shift,
    ellipse_shift_weight,
    finite_shift_profile_metadata,
    relative_cladding_layer,
    save_finite_shift_profile_diagnostics,
    shift_profile_preflight_messages,
)
from comsol_workflow.field_plotting import save_field_plot  # noqa: E402
from comsol_workflow.finite_geometry import (  # noqa: E402
    FiniteGeometryPlan,
    FiniteGeometrySpec,
    build_finite_geometry_plan,
    square_cell_plot_records,
)
from comsol_workflow.finite_geometry_comsol import (  # noqa: E402
    compile_finite_geometry_input,
)
from comsol_workflow.geometry_utils import create_hexagon_design  # noqa: E402
from comsol_workflow.hex_lattice_utils import (  # noqa: E402
    LatticePoint,
    boundary_from_lattice_cells,
    bulk_points as _bulk_points,
    cell_polygon as _cell_polygon,
    lattice_points,
    retained_cladding_points as _retained_cladding_points,
    unit_cell_corners as _unit_cell_corners,
)
from comsol_workflow.lattice_fourier_postprocess import (  # noqa: E402
    run_finite_lattice_fourier_postprocess as run_exported_finite_lattice_fourier_postprocess,
    run_finite_unit_cell_intensity_maps as run_exported_finite_unit_cell_intensity_maps,
    score_finite_modes as score_exported_finite_modes,
    select_finite_fundamental_mode,
    select_finite_mode_closest_to_frequency,
)
from comsol_workflow.farfield_fft import (  # noqa: E402
    run_farfield_fft_case,
    validate_farfield_config,
)
from comsol_workflow.polygon_utils import (  # noqa: E402
    offset_polygon,
    point_inside_or_on_polygon,
)
from scripts.run_main.parameter_config import (  # noqa: E402
    ACTIVE_PARAMETERS,
    FINITE_CAVITY_OUTPUT_ROOT,
)

EXPORT_DIRNAME = "simulation_exports"
OUTPUT_LAYOUT_VERSION = "mode_centric_v6"
MODEL_DIRNAME = "00_model"
RESULTS_DIRNAME = "01_results"
OVERVIEW_DIRNAME = "10_overview"
LOGS_DIRNAME = "80_logs"
CONFIG_DIRNAME = "99_config"
STAGING_DIRNAME = ".staging"
OUTPUT_ROOT = FINITE_CAVITY_OUTPUT_ROOT
DPI = 220

A = 0.82
R_0 = A / 3.0
B_0 = 0.23
D = 0.02
CAVITY_LAYERS = ACTIVE_PARAMETERS.cavity_layers
BULK_RADIUS = CAVITY_LAYERS - 1

FINITE_CAVITY_PRESET = "shared_parameter_json"
FINITE_CAVITY_PRESETS = {
    FINITE_CAVITY_PRESET: {
        "cavity": ACTIVE_PARAMETERS.cavity.to_dict(),
        "cladding": ACTIVE_PARAMETERS.cladding.to_dict(),
        "px_gamma_frequency_thz": ACTIVE_PARAMETERS.center_frequency_thz,
        "cladding_shift_factors": list(
            ACTIVE_PARAMETERS.cladding_shift_factors
        ),
        "finite_cladding_shift_input": (
            ACTIVE_PARAMETERS.finite_cladding_shift_input
        ),
        "cladding_shift_geometry": (
            ACTIVE_PARAMETERS.cladding_shift_geometry.to_dict()
        ),
        "cladding_shift_profile": ACTIVE_PARAMETERS.cladding_shift_profile.to_dict(),
        "finite_geometry": ACTIVE_PARAMETERS.finite_geometry,
    }
}


def _format_decimal(value: float, decimal_places: int = 12) -> str:
    return f"{float(value):.{decimal_places}f}".rstrip("0").rstrip(".")


def _resolve_cell_params(
    cell_label: str,
    cell_name: str,
    values: dict[str, object],
) -> tuple[dict[str, float | str], dict[str, float]]:
    required = {"b0_nm", "eta", "zeta"}
    missing = sorted(required.difference(values))
    if missing:
        raise ValueError(
            f"{cell_label} is missing required fields: {', '.join(missing)}"
        )

    resolved_values = {key: float(values[key]) for key in required}
    for key in ("b0_nm", "eta"):
        value = resolved_values[key]
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{cell_label} {key} must be positive and finite")
    if not np.isfinite(resolved_values["zeta"]):
        raise ValueError(f"{cell_label} zeta must be finite")

    return (
        {
            "name": cell_name,
            "r_f0": resolved_values["eta"],
            "b_square_f0": (resolved_values["b0_nm"] / (B_0 * 1000.0)) ** 2,
            "b_square_f3": (resolved_values["zeta"] ** 2 - 1.0) / 2.0,
        },
        resolved_values,
    )


def resolve_finite_cavity_preset(
    preset_name: str,
    presets: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    registry = FINITE_CAVITY_PRESETS if presets is None else presets
    if preset_name not in registry:
        raise ValueError(f"Unknown finite-cavity preset: {preset_name}")
    preset = registry[preset_name]
    required = {
        "cavity",
        "cladding",
        "px_gamma_frequency_thz",
    }
    missing = sorted(required.difference(preset))
    if missing:
        raise ValueError(f"Preset is missing required fields: {', '.join(missing)}")

    bulk_params, cavity_values = _resolve_cell_params(
        "cavity", "cavity_p_bic", preset["cavity"]
    )
    cladding_params, cladding_values = _resolve_cell_params(
        "cladding", "cladding_bandgap_alignment", preset["cladding"]
    )
    frequency_thz = float(preset["px_gamma_frequency_thz"])
    if not np.isfinite(frequency_thz) or frequency_thz <= 0.0:
        raise ValueError("px_gamma_frequency_thz must be positive and finite")

    raw_shift_factors = preset.get(
        "cladding_shift_factors",
        preset.get("cladding_inward_shift_factors"),
    )
    if raw_shift_factors is None:
        shift_factors = []
    else:
        if not isinstance(raw_shift_factors, (list, tuple)):
            raise TypeError("cladding_shift_factors must be a list")
        shift_factors = [float(value) for value in raw_shift_factors]
        if not all(np.isfinite(value) and value >= 0.0 for value in shift_factors):
            raise ValueError(
                "cladding_shift_factors must be finite and non-negative"
            )

    structured_input = preset.get("finite_cladding_shift_input")
    if structured_input is not None:
        if not isinstance(structured_input, dict):
            raise TypeError("finite_cladding_shift_input must be an object")
        input_kind = structured_input.get("kind")
        if input_kind == "scan_ratio":
            allowed_input_fields = {"kind", "y_over_x_ratio"}
            if "y_over_x_ratio" not in structured_input:
                raise ValueError(
                    "finite_cladding_shift_input is missing y_over_x_ratio"
                )
            y_over_x_ratio = float(structured_input["y_over_x_ratio"])
            explicit_x_factor = None
            explicit_y_factor = None
        elif input_kind == "xy":
            allowed_input_fields = {"kind", "x_factor", "y_factor"}
            missing_input_fields = allowed_input_fields.difference(
                structured_input
            )
            if missing_input_fields:
                raise ValueError(
                    "finite_cladding_shift_input is missing: "
                    + ", ".join(sorted(missing_input_fields))
                )
            y_over_x_ratio = None
            explicit_x_factor = float(structured_input["x_factor"])
            explicit_y_factor = float(structured_input["y_factor"])
        else:
            raise ValueError(
                "finite_cladding_shift_input.kind must be scan_ratio or xy"
            )
        unknown_input_fields = sorted(
            set(structured_input).difference(allowed_input_fields)
        )
        if unknown_input_fields:
            raise ValueError(
                "finite_cladding_shift_input contains unknown fields: "
                + ", ".join(unknown_input_fields)
            )
    else:
        has_ratio = "cladding_y_over_x_shift_ratio" in preset
        has_x_factor = "cladding_x_shift_factor" in preset
        has_y_factor = "cladding_y_shift_factor" in preset
        if has_ratio and (has_x_factor or has_y_factor):
            raise ValueError(
                "cladding_y_over_x_shift_ratio conflicts with explicit x/y"
            )
        if has_x_factor != has_y_factor:
            raise ValueError(
                "cladding_x_shift_factor and cladding_y_shift_factor "
                "must be provided together"
            )
        if has_x_factor:
            input_kind = "xy"
            y_over_x_ratio = None
            explicit_x_factor = float(preset["cladding_x_shift_factor"])
            explicit_y_factor = float(preset["cladding_y_shift_factor"])
        else:
            input_kind = "scan_ratio"
            y_over_x_ratio = float(
                preset.get("cladding_y_over_x_shift_ratio", 1.0)
            )
            explicit_x_factor = None
            explicit_y_factor = None

    if input_kind == "scan_ratio":
        if not shift_factors:
            raise ValueError("cladding_shift_factors must be non-empty")
        if not np.isfinite(y_over_x_ratio) or y_over_x_ratio < 0.0:
            raise ValueError(
                "cladding_y_over_x_shift_ratio must be finite and non-negative"
            )
        finite_shift_input = {
            "kind": "scan_ratio",
            "y_over_x_ratio": y_over_x_ratio,
        }
        shift_factor_pairs = [
            (factor, factor * y_over_x_ratio)
            for factor in shift_factors
        ]
    else:
        if not np.isfinite(explicit_x_factor) or explicit_x_factor < 0.0:
            raise ValueError(
                "cladding_x_shift_factor must be finite and non-negative"
            )
        if not np.isfinite(explicit_y_factor) or explicit_y_factor < 0.0:
            raise ValueError(
                "cladding_y_shift_factor must be finite and non-negative"
            )
        finite_shift_input = {
            "kind": "xy",
            "x_factor": explicit_x_factor,
            "y_factor": explicit_y_factor,
        }
        shift_factor_pairs = [(explicit_x_factor, explicit_y_factor)]
    shift_geometry = CladdingShiftGeometry.from_mapping(
        preset.get("cladding_shift_geometry")
    )
    finite_geometry = preset.get("finite_geometry", "hex")
    if finite_geometry not in {"hex", "square"}:
        raise ValueError("finite_geometry must be 'hex' or 'square'")
    shift_profile = CladdingShiftProfile.from_mapping(
        preset.get("cladding_shift_profile")
    )

    cavity_label = (
        f"cav({_format_decimal(cavity_values['b0_nm'], 3)}"
        f"-{_format_decimal(cavity_values['eta'], 3)}"
        f"-{_format_decimal(cavity_values['zeta'], 3)})"
    )
    cladding_label = (
        f"clad({_format_decimal(cladding_values['b0_nm'], 3)}"
        f"-{_format_decimal(cladding_values['eta'], 3)}"
        f"-{_format_decimal(cladding_values['zeta'], 3)})"
    )
    return {
        "preset_name": preset_name,
        "cavity": cavity_values,
        "cladding": cladding_values,
        "bulk_params": bulk_params,
        "cladding_params": cladding_params,
        "px_gamma_frequency_thz": frequency_thz,
        "eigenfrequency_shift": f"{_format_decimal(frequency_thz)} [THz]",
        "cladding_shift_factors": shift_factors,
        "finite_cladding_shift_input": finite_shift_input,
        "finite_cladding_shift_factor_pairs": shift_factor_pairs,
        "cladding_shift_geometry": shift_geometry.to_dict(),
        "finite_geometry": finite_geometry,
        "cladding_shift_profile": shift_profile.to_dict(),
        "cladding_shift_profile_label": shift_profile.label,
        "case_parameter_label": f"{cavity_label}_{cladding_label}",
    }


FINITE_CAVITY_CONFIG = {
    "preset_name": FINITE_CAVITY_PRESET,
    "cavity": ACTIVE_PARAMETERS.cavity.to_dict(),
    "cladding": ACTIVE_PARAMETERS.cladding.to_dict(),
    "bulk_params": ACTIVE_PARAMETERS.cavity.to_fourier_params("cavity_p_bic"),
    "cladding_params": ACTIVE_PARAMETERS.cladding.to_fourier_params(
        "cladding_bandgap_alignment"
    ),
    "px_gamma_frequency_thz": ACTIVE_PARAMETERS.center_frequency_thz,
    "eigenfrequency_shift": ACTIVE_PARAMETERS.eigenfrequency_shift,
    "cladding_shift_factors": list(ACTIVE_PARAMETERS.cladding_shift_factors),
    "finite_cladding_shift_input": (
        ACTIVE_PARAMETERS.finite_cladding_shift_input
    ),
    "finite_cladding_shift_factor_pairs": [
        list(pair)
        for pair in ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
    ],
    "cladding_shift_geometry": ACTIVE_PARAMETERS.cladding_shift_geometry.to_dict(),
    "cladding_shift_profile": ACTIVE_PARAMETERS.cladding_shift_profile.to_dict(),
    "cladding_shift_profile_label": ACTIVE_PARAMETERS.cladding_shift_profile.label,
    "finite_geometry": ACTIVE_PARAMETERS.finite_geometry,
    "case_parameter_label": (
        f"{ACTIVE_PARAMETERS.cavity.label('cav')}_"
        f"{ACTIVE_PARAMETERS.cladding.label('clad')}"
    ),
}
CLADDING_SHIFT_FACTORS = list(FINITE_CAVITY_CONFIG["cladding_shift_factors"])
CLADDING_SHIFT_INPUT_KIND = ACTIVE_PARAMETERS.cladding_shift_input_kind
CLADDING_Y_OVER_X_SHIFT_RATIO = (
    ACTIVE_PARAMETERS.cladding_y_over_x_shift_ratio
)
CLADDING_X_SHIFT_FACTOR, CLADDING_Y_SHIFT_FACTOR = (
    ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs[0]
)
CLADDING_SHIFT_GEOMETRY = ACTIVE_PARAMETERS.cladding_shift_geometry
CLADDING_SHIFT_PROFILE = ACTIVE_PARAMETERS.cladding_shift_profile
FINITE_GEOMETRY = ACTIVE_PARAMETERS.finite_geometry
CLADDING_LAYERS = ACTIVE_PARAMETERS.cladding_layers
OUT_DIR = OUTPUT_ROOT / ACTIVE_PARAMETERS.structure_series_label("finite")
OUTER_EMPTY_PAD_PERIODS = 0.5
CLADDING_X_INWARD_SHIFT = CLADDING_X_SHIFT_FACTOR * A
CLADDING_Y_INWARD_SHIFT = CLADDING_Y_SHIFT_FACTOR * A
SHIFT_ARROW_PLOT_SCALE = 10.0
STRIP_X_MIN = -A / 2.0
STRIP_X_MAX = A / 2.0

SLAB_HEIGHT = "200 [nm]"
REFRACTIVE_INDEX = "3.3"
EIGENMODE_COUNT = ACTIVE_PARAMETERS.finite_eigenmode_count
MESH_AUTO_SIZE = ACTIVE_PARAMETERS.mesh_size
USE_MANUAL_MESH_CONSTRUCTION = True
MODE_TYPE = "TE"
WAVELENGTH = "1550 [nm]"
EIGENFREQUENCY_SHIFT = str(FINITE_CAVITY_CONFIG["eigenfrequency_shift"])

RESUME_FROM_EXISTING_MPH = False
RUN_FINITE_LATTICE_FOURIER_POSTPROCESS = True
FINITE_LATTICE_FOURIER_MODE_INDICES = None
FINITE_LATTICE_FOURIER_RHO_GRID_SIZE = 41
FINITE_LATTICE_FOURIER_TOP_K = 6
FINITE_SINGLE_MODE_ANALYSIS_ENABLED = True
FINITE_SINGLE_MODE_TARGET_FREQUENCY_THZ = float(
    FINITE_CAVITY_CONFIG["px_gamma_frequency_thz"]
)
RUN_FARFIELD_FFT = True
FARFIELD_GRID_SIZE = 801
FARFIELD_FFT_SIZE = 2001
FARFIELD_NA = 0.9
FARFIELD_H0_UM = None

BULK_PARAMS = dict(FINITE_CAVITY_CONFIG["bulk_params"])
CLADDING_PARAMS = dict(FINITE_CAVITY_CONFIG["cladding_params"])

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


def bulk_points() -> list[LatticePoint]:
    return _bulk_points(BULK_RADIUS, A)


def bulk_boundary(points: list[LatticePoint]) -> np.ndarray:
    return boundary_from_lattice_cells(points, A)


def retained_cladding_points(boundary: np.ndarray) -> list[LatticePoint]:
    _, retained = _retained_cladding_points(A, boundary, CLADDING_LAYERS)
    return retained


def translate_polygon(polygon: np.ndarray, x: float, y: float) -> np.ndarray:
    return np.asarray(polygon, dtype=float) + np.array([x, y], dtype=float)


def hexagon_boundary_vertices(a: float) -> np.ndarray:
    root3 = np.sqrt(3.0)
    return np.array([
        [a / 2.0, 0.0],
        [a / 4.0, root3 * a / 4.0],
        [-a / 4.0, root3 * a / 4.0],
        [-a / 2.0, 0.0],
        [-a / 4.0, -root3 * a / 4.0],
        [a / 4.0, -root3 * a / 4.0],
    ], dtype=float)


def hexagon_boundary_size_for_points(points: np.ndarray, pad: float) -> float:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    if points.shape[0] == 0:
        raise ValueError("points must not be empty")

    root3 = np.sqrt(3.0)
    x = points[:, 0]
    y = points[:, 1]
    min_a = max(
        float(4.0 * np.max(np.abs(y)) / root3),
        float(2.0 * np.max(np.abs(root3 * x + y)) / root3),
        float(2.0 * np.max(np.abs(root3 * x - y)) / root3),
    )
    return min_a + 2.0 * pad


def cladding_inward_shift(point: LatticePoint) -> np.ndarray:
    layer = relative_cladding_layer(point, BULK_RADIUS)
    envelope = CLADDING_SHIFT_PROFILE.envelope(layer)
    if CLADDING_SHIFT_GEOMETRY.kind == ELLIPSE_SHIFT_GEOMETRY:
        return np.asarray(
            ellipse_inward_shift(
                point.x,
                point.y,
                A,
                CLADDING_X_SHIFT_FACTOR,
                CLADDING_Y_SHIFT_FACTOR,
                envelope=envelope,
            ),
            dtype=float,
        )
    kind, idx = normal_fan_feature(point)
    return normal_fan_shift(kind, idx, envelope=envelope)


def current_cladding_y_over_x_shift_ratio() -> float | None:
    """Return configured ratio provenance; explicit x/y mode has no ratio."""
    if CLADDING_SHIFT_INPUT_KIND == "scan_ratio":
        return float(CLADDING_Y_OVER_X_SHIFT_RATIO)
    return None


def current_cladding_shift_base_factor() -> float | None:
    """Return the shared scan value, or None for an explicit x/y pair."""
    if CLADDING_SHIFT_INPUT_KIND == "scan_ratio":
        return CLADDING_X_SHIFT_FACTOR
    return None


def resolved_ellipse_cell_shift(
    point: LatticePoint,
    shift: np.ndarray | None = None,
) -> dict[str, float | int]:
    """Describe the ideal-center radial shift applied to one ellipse-mode cell."""
    layer = relative_cladding_layer(point, BULK_RADIUS)
    envelope = CLADDING_SHIFT_PROFILE.envelope(layer)
    resolved_shift = (
        cladding_inward_shift(point)
        if shift is None
        else np.asarray(shift, dtype=float)
    )
    radius = float(np.hypot(point.x, point.y))
    return {
        "i": int(point.i),
        "j": int(point.j),
        "shell": int(point.shell),
        "layer": layer,
        "ideal_x_um": float(point.x),
        "ideal_y_um": float(point.y),
        "radius_um": radius,
        "theta_degrees": float(np.degrees(np.arctan2(point.y, point.x))),
        "envelope": float(envelope),
        "angular_factor": ellipse_shift_weight(
            point.x,
            point.y,
            CLADDING_X_SHIFT_FACTOR,
            CLADDING_Y_SHIFT_FACTOR,
        ),
        "shift_x_um": float(resolved_shift[0]),
        "shift_y_um": float(resolved_shift[1]),
        "shift_magnitude_um": float(np.linalg.norm(resolved_shift)),
    }


def build_full_lattice_holes() -> tuple[list[dict[str, object]], dict[str, object]]:
    if FINITE_GEOMETRY == "square":
        plan = current_finite_geometry_plan()
        fragment_support_records = [
            {
                "region": placement.region,
                "point": placement.point,
                "group_id": placement.group_id,
                "layer": placement.layer,
                "outside_x": placement.outside_x,
                "outside_y": placement.outside_y,
                "modulation_shift": placement.shift,
                "polygon": np.array(
                    [
                        [placement.support[0], placement.support[2]],
                        [placement.support[1], placement.support[2]],
                        [placement.support[1], placement.support[3]],
                        [placement.support[0], placement.support[3]],
                    ],
                    dtype=float,
                ),
            }
            for placement in plan.placements
        ]
        resolved_cell_shifts = []
        seen_parents: set[tuple[int, int]] = set()
        for placement in plan.placements:
            if placement.region != "cladding":
                continue
            key = (placement.m, placement.n)
            if key in seen_parents:
                continue
            seen_parents.add(key)
            point = placement.point
            resolved_cell_shifts.append(
                {
                    "i": int(point.i),
                    "j": int(point.j),
                    "shell": int(point.shell),
                    "layer": int(placement.layer),
                    "ideal_x_um": float(point.x),
                    "ideal_y_um": float(point.y),
                    "shift_x_um": float(placement.shift[0]),
                    "shift_y_um": float(placement.shift[1]),
                    "shift_magnitude_um": float(np.linalg.norm(placement.shift)),
                }
            )
        metadata = {
            "a": A,
            "cavity_layers": CAVITY_LAYERS,
            "bulk_radius": BULK_RADIUS,
            "cladding_layers": CLADDING_LAYERS,
            "outer_empty_pad_periods": OUTER_EMPTY_PAD_PERIODS,
            "cladding_x_shift_factor": CLADDING_X_SHIFT_FACTOR,
            "cladding_y_shift_factor": CLADDING_Y_SHIFT_FACTOR,
            "cladding_shift_input_kind": CLADDING_SHIFT_INPUT_KIND,
            "cladding_shift_base_factor": current_cladding_shift_base_factor(),
            "cladding_y_over_x_shift_ratio": current_cladding_y_over_x_shift_ratio(),
            "cladding_x_inward_shift": CLADDING_X_INWARD_SHIFT,
            "cladding_y_inward_shift": CLADDING_Y_INWARD_SHIFT,
            **finite_geometry_shift_metadata(),
            **plan.metadata,
            "resolved_cladding_shift_by_cell": resolved_cell_shifts,
            "shift_arrow_plot_scale": SHIFT_ARROW_PLOT_SCALE,
            "bulk_params": BULK_PARAMS,
            "cladding_params": CLADDING_PARAMS,
            "fragment_support_records": fragment_support_records,
            "geometry_plan": plan,
        }
        return list(plan.records), metadata

    bulk = bulk_points()
    inner_boundary = bulk_boundary(bulk)
    cladding = lattice_points(
        BULK_RADIUS + CLADDING_LAYERS,
        A,
        min_shell=BULK_RADIUS + 1,
    )

    bulk_unit_holes = create_unit_cell_holes(BULK_PARAMS)
    cladding_unit_holes = create_unit_cell_holes(CLADDING_PARAMS)

    records: list[dict[str, object]] = []
    resolved_cell_shifts: list[dict[str, float | int]] = []
    for point in bulk:
        for hole in bulk_unit_holes:
            records.append({
                "region": "bulk",
                "point": point,
                "parent_layer": point.shell + 1,
                "cladding_shift_layer": None,
                "polygon": translate_polygon(hole, point.x, point.y),
                "modulation_shift": np.zeros(2, dtype=float),
            })
    for point in cladding:
        shift = cladding_inward_shift(point)
        if CLADDING_SHIFT_GEOMETRY.kind == ELLIPSE_SHIFT_GEOMETRY:
            resolved_cell_shifts.append(resolved_ellipse_cell_shift(point, shift))
        for hole in cladding_unit_holes:
            records.append({
                "region": "cladding",
                "point": point,
                "parent_layer": point.shell + 1,
                "cladding_shift_layer": point.shell - BULK_RADIUS,
                "polygon": translate_polygon(hole, point.x + shift[0], point.y + shift[1]),
                "modulation_shift": shift,
            })

    profile_metadata = finite_shift_profile_metadata(
        CLADDING_SHIFT_GEOMETRY,
        CLADDING_LAYERS,
        CLADDING_X_SHIFT_FACTOR,
        CLADDING_Y_SHIFT_FACTOR,
        A,
        CLADDING_SHIFT_PROFILE,
    )
    geometry_target_metadata: dict[str, object]
    if CLADDING_SHIFT_GEOMETRY.kind == GROUPS_SHIFT_GEOMETRY:
        geometry_target_metadata = {
            "cladding_x_side_inward_shift_target": normal_fan_shift_magnitude(
                "side", target_factor=CLADDING_X_SHIFT_FACTOR
            ),
            "cladding_y_side_inward_shift_target": normal_fan_shift_magnitude(
                "side", target_factor=CLADDING_Y_SHIFT_FACTOR
            ),
            "cladding_x_corner_inward_shift_target": normal_fan_shift_magnitude(
                "corner", target_factor=CLADDING_X_SHIFT_FACTOR
            ),
            "cladding_y_corner_inward_shift_target": normal_fan_shift_magnitude(
                "corner", target_factor=CLADDING_Y_SHIFT_FACTOR
            ),
        }
    else:
        geometry_target_metadata = {
            "cladding_x_axis_inward_shift_target": CLADDING_X_INWARD_SHIFT,
            "cladding_y_axis_inward_shift_target": CLADDING_Y_INWARD_SHIFT,
        }
    metadata = {
        "a": A,
        "finite_geometry": "hex",
        "finite_geometry_schema": "hex_full_cells_v1",
        "finite_geometry_shift_formula": (
            "hex_radial_ellipse_v1"
            if CLADDING_SHIFT_GEOMETRY.kind == ELLIPSE_SHIFT_GEOMETRY
            else "hex_normal_fan_v1"
        ),
        "cavity_layers": CAVITY_LAYERS,
        "bulk_radius": BULK_RADIUS,
        "cladding_layers": CLADDING_LAYERS,
        "outer_empty_pad_periods": OUTER_EMPTY_PAD_PERIODS,
        "cladding_x_shift_factor": CLADDING_X_SHIFT_FACTOR,
        "cladding_y_shift_factor": CLADDING_Y_SHIFT_FACTOR,
        "cladding_shift_input_kind": CLADDING_SHIFT_INPUT_KIND,
        "cladding_shift_base_factor": current_cladding_shift_base_factor(),
        "cladding_y_over_x_shift_ratio": (
            current_cladding_y_over_x_shift_ratio()
        ),
        "cladding_x_inward_shift": CLADDING_X_INWARD_SHIFT,
        "cladding_y_inward_shift": CLADDING_Y_INWARD_SHIFT,
        **geometry_target_metadata,
        **profile_metadata,
        "resolved_cladding_shift_by_cell": resolved_cell_shifts,
        "shift_arrow_plot_scale": SHIFT_ARROW_PLOT_SCALE,
        "bulk_params": BULK_PARAMS,
        "cladding_params": CLADDING_PARAMS,
        "bulk_cells": len(bulk),
        "cladding_cells": len(cladding),
        "bulk_holes": len(bulk) * len(bulk_unit_holes),
        "cladding_holes": len(cladding) * len(cladding_unit_holes),
        "total_holes": len(records),
        "inner_boundary": inner_boundary,
        "bulk_points": bulk,
        "cladding_points": cladding,
    }
    return records, metadata


def current_finite_geometry_plan() -> FiniteGeometryPlan:
    return build_finite_geometry_plan(
        FiniteGeometrySpec(
            shape=FINITE_GEOMETRY,
            period=A,
            cavity_layers=CAVITY_LAYERS,
            cladding_layers=CLADDING_LAYERS,
            empty_buffer_periods=OUTER_EMPTY_PAD_PERIODS,
        ),
        create_unit_cell_holes(BULK_PARAMS),
        create_unit_cell_holes(CLADDING_PARAMS),
        shift_geometry=CLADDING_SHIFT_GEOMETRY,
        shift_profile=CLADDING_SHIFT_PROFILE,
        x_shift_factor=CLADDING_X_SHIFT_FACTOR,
        y_shift_factor=CLADDING_Y_SHIFT_FACTOR,
    )


def finite_geometry_shift_metadata() -> dict[str, object]:
    metadata = finite_shift_profile_metadata(
        CLADDING_SHIFT_GEOMETRY,
        CLADDING_LAYERS,
        CLADDING_X_SHIFT_FACTOR,
        CLADDING_Y_SHIFT_FACTOR,
        A,
        CLADDING_SHIFT_PROFILE,
    )
    if FINITE_GEOMETRY != "square":
        return metadata
    if CLADDING_SHIFT_GEOMETRY.kind == GROUPS_SHIFT_GEOMETRY:
        formula = {
            "kind": "square_axis_components_v1",
            "x_shift": "-sign(x)*A*fx*envelope(layer) when outside_x",
            "y_shift": "-sign(y)*A*fy*envelope(layer) when outside_y",
            "corner_scale": 1.0,
        }
        metadata["cladding_shift_region_groups"] = {
            "x": {"axes": ["left", "right"]},
            "y": {"axes": ["bottom", "top"]},
        }
    else:
        formula = {
            **dict(metadata["cladding_shift_formula"]),
            "kind": "square_radial_ellipse_v1",
        }
    metadata["cladding_shift_formula"] = formula
    metadata["resolved_cladding_shift"] = {
        **dict(metadata["resolved_cladding_shift"]),
        "formula_version": formula["kind"],
    }
    return metadata


def cell_records(points, region: str) -> list[dict[str, object]]:
    return [
        {
            "region": region,
            "point": point,
            "polygon": cell_polygon(point.x, point.y),
        }
        for point in points
    ]


def full_cell_records(metadata: dict[str, object]) -> list[dict[str, object]]:
    if metadata.get("finite_geometry") == "square":
        return list(metadata["fragment_support_records"])
    return [
        *cell_records(metadata["bulk_points"], "bulk"),
        *cell_records(metadata["cladding_points"], "cladding"),
    ]


def strip_validity_layer_count() -> int:
    return CLADDING_LAYERS // 2


def strip_validity_boundary(
    metadata: dict[str, object] | None = None,
) -> np.ndarray:
    if FINITE_GEOMETRY == "square":
        if metadata is not None and isinstance(
            metadata.get("geometry_plan"), FiniteGeometryPlan
        ):
            return metadata["geometry_plan"].validity_boundary
        return current_finite_geometry_plan().validity_boundary
    radius = BULK_RADIUS + strip_validity_layer_count()
    return boundary_from_lattice_cells(lattice_points(radius, A), A)


def case_export_dir(out_dir: Path) -> Path:
    return out_dir / EXPORT_DIRNAME


def finite_case_output_paths(
    case_dir: Path,
    *,
    model_filename: str,
) -> dict[str, Path]:
    case_dir = Path(case_dir)
    model_dir = case_dir / MODEL_DIRNAME
    staging_dir = case_dir / STAGING_DIRNAME
    return {
        "case_dir": case_dir,
        "model_dir": model_dir,
        "results_dir": case_dir / RESULTS_DIRNAME,
        "overview_dir": case_dir / OVERVIEW_DIRNAME,
        "logs_dir": case_dir / LOGS_DIRNAME,
        "config_dir": case_dir / CONFIG_DIRNAME,
        "staging_dir": staging_dir,
        "export_dir": staging_dir / EXPORT_DIRNAME,
        "mph": model_dir / model_filename,
        "progress_log": model_dir / "comsol_progress.log",
        "config": case_dir / CONFIG_DIRNAME / "config.json",
    }


def finite_case_output_metadata(paths: dict[str, Path]) -> dict[str, object]:
    return {
        "case_dir": str(paths["case_dir"]),
        "mph": str(paths["mph"]),
        "mode_output_pattern": str(paths["results_dir"] / "mode{mode_idx}"),
        "model_dir": str(paths["model_dir"]),
        "overview_dir": str(paths["overview_dir"]),
        "logs_dir": str(paths["logs_dir"]),
        "config_dir": str(paths["config_dir"]),
    }


def prepare_finite_case_output(
    case_dir: Path,
    *,
    model_filename: str,
    resume: bool,
    create_staging: bool = True,
) -> dict[str, Path]:
    paths = finite_case_output_paths(case_dir, model_filename=model_filename)
    paths["case_dir"].mkdir(parents=True, exist_ok=True)

    if paths["results_dir"].exists() and any(paths["results_dir"].iterdir()):
        raise FileExistsError(
            "Refusing to overwrite existing finite mode results: "
            f"{paths['results_dir']}"
        )
    if not resume:
        occupied = [
            path
            for path in (paths["mph"], paths["config"])
            if path.exists()
        ]
        if occupied:
            raise FileExistsError(
                "Refusing to overwrite an existing finite case; use the resume "
                f"workflow or a new output directory: {occupied[0]}"
            )

    for key in (
        "model_dir",
        "results_dir",
        "overview_dir",
        "logs_dir",
        "config_dir",
    ):
        paths[key].mkdir(parents=True, exist_ok=True)

    staging_dir = paths["staging_dir"]
    if staging_dir.exists():
        if any(staging_dir.iterdir()):
            raise FileExistsError(
                "Refusing to overwrite non-empty finite staging directory: "
                f"{staging_dir}"
            )
        staging_dir.rmdir()
    if create_staging:
        staging_dir.mkdir(parents=True)
    return paths


def _mode_result_dir(paths: dict[str, Path], mode_idx: int) -> Path:
    return paths["results_dir"] / f"mode{int(mode_idx)}"


def _add_output_move(
    moves: list[tuple[Path, Path]],
    source: Path,
    destination: Path,
) -> None:
    if not source.is_file():
        raise FileNotFoundError(f"Missing staged finite output: {source}")
    moves.append((source, destination))


def _validate_result_mode_rows(
    frame: pd.DataFrame,
    source_path: Path,
    known_modes: set[int],
) -> None:
    if frame.empty and "mode_idx" not in frame.columns:
        return
    if "mode_idx" not in frame.columns:
        raise ValueError(f"Missing mode_idx column in {source_path}")
    mode_indices = frame["mode_idx"].astype(int).tolist()
    unknown_modes = sorted(set(mode_indices).difference(known_modes))
    if unknown_modes:
        raise ValueError(
            f"{source_path} references unknown mode {unknown_modes[0]}"
        )
    if len(mode_indices) != len(set(mode_indices)):
        raise ValueError(f"Duplicate mode_idx values in {source_path}")


def _rewrite_farfield_outputs(
    paths: dict[str, Path],
    mode_indices: list[int],
) -> None:
    for mode_idx in mode_indices:
        mode_dir = _mode_result_dir(paths, mode_idx)
        summary_path = mode_dir / "12_farfield_FFT" / "farfield_summary.json"
        if not summary_path.is_file():
            continue
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        source_path = mode_dir / "11_simulation_exports" / "E_air.parquet"
        if source_path.is_file():
            payload["source_file"] = str(source_path)
        payload["plots"] = [
            str(path)
            for path in sorted(
                (mode_dir / "12_farfield_FFT").glob("*.png"),
                key=lambda item: item.name,
            )
        ]
        write_config(summary_path, payload)

    aggregate_path = paths["overview_dir"] / "farfield_summary.csv"
    if aggregate_path.is_file():
        try:
            summary_df = pd.read_csv(aggregate_path)
        except pd.errors.EmptyDataError:
            return
        for row_idx, row in summary_df.iterrows():
            mode_idx = int(row["mode_idx"])
            mode_dir = _mode_result_dir(paths, mode_idx)
            source_path = mode_dir / "11_simulation_exports" / "E_air.parquet"
            if source_path.is_file():
                summary_df.at[row_idx, "source_file"] = str(source_path)
            plots = [
                str(path)
                for path in sorted(
                    (mode_dir / "12_farfield_FFT").glob("*.png"),
                    key=lambda item: item.name,
                )
            ]
            summary_df.at[row_idx, "plots"] = json.dumps(plots)
        summary_df.to_csv(aggregate_path, index=False)


def _rewrite_lattice_summary(
    paths: dict[str, Path],
) -> None:
    summary_path = paths["overview_dir"] / "finite_lattice_fourier_summary.csv"
    if not summary_path.is_file():
        return
    try:
        summary_df = pd.read_csv(summary_path)
    except pd.errors.EmptyDataError:
        return
    for row_idx, row in summary_df.iterrows():
        mode_idx = int(row["mode_idx"])
        summary_df.at[row_idx, "output_dir"] = str(
            _mode_result_dir(paths, mode_idx) / "13_lattice_fourier_Hz"
        )
    summary_df.to_csv(summary_path, index=False)


def finalize_finite_case_output(
    case_dir: Path,
    *,
    model_filename: str,
) -> dict[str, Path]:
    paths = finite_case_output_paths(case_dir, model_filename=model_filename)
    staging_dir = paths["staging_dir"]
    export_dir = paths["export_dir"]
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.is_file():
        raise FileNotFoundError(
            f"Cannot finalize finite output without {eigen_path}"
        )

    eigen_df = pd.read_csv(eigen_path)
    if "mode_idx" not in eigen_df.columns:
        raise ValueError(f"Missing mode_idx column in {eigen_path}")
    mode_indices = [int(value) for value in eigen_df["mode_idx"].tolist()]
    if len(mode_indices) != len(set(mode_indices)):
        raise ValueError(f"Duplicate mode_idx values in {eigen_path}")
    known_modes = set(mode_indices)

    moves: list[tuple[Path, Path]] = []
    _add_output_move(
        moves,
        eigen_path,
        paths["overview_dir"] / "eigenfrequencies.csv",
    )

    air_metadata = export_dir / "air_field_metadata.json"
    if air_metadata.is_file():
        _add_output_move(
            moves,
            air_metadata,
            paths["config_dir"] / "air_field_metadata.json",
        )

    for source in sorted(export_dir.glob("*"), key=lambda item: item.name):
        if not source.is_file() or source in {eigen_path, air_metadata}:
            continue
        match = re.fullmatch(r"(\d+)_(.+)", source.name)
        if match is None:
            raise ValueError(f"Unknown staged simulation export: {source}")
        mode_idx = int(match.group(1))
        if mode_idx not in known_modes:
            raise ValueError(
                f"Simulation export references unknown mode {mode_idx}: {source}"
            )
        _add_output_move(
            moves,
            source,
            _mode_result_dir(paths, mode_idx)
            / "11_simulation_exports"
            / match.group(2),
        )

    farfield_dir = staging_dir / "farfield_fft"
    if farfield_dir.is_dir():
        farfield_config = farfield_dir / "config.json"
        if farfield_config.is_file():
            _add_output_move(
                moves,
                farfield_config,
                paths["config_dir"] / "farfield_config.json",
            )
        farfield_aggregate = farfield_dir / "farfield_summary.csv"
        if farfield_aggregate.is_file():
            try:
                farfield_frame = pd.read_csv(farfield_aggregate)
            except pd.errors.EmptyDataError:
                farfield_frame = pd.DataFrame()
            _validate_result_mode_rows(
                farfield_frame,
                farfield_aggregate,
                known_modes,
            )
            _add_output_move(
                moves,
                farfield_aggregate,
                paths["overview_dir"] / "farfield_summary.csv",
            )
        for mode_source_dir in sorted(farfield_dir.glob("mode_*")):
            match = re.fullmatch(r"mode_(\d+)", mode_source_dir.name)
            if not mode_source_dir.is_dir() or match is None:
                raise ValueError(f"Unknown staged far-field entry: {mode_source_dir}")
            mode_idx = int(match.group(1))
            if mode_idx not in known_modes:
                raise ValueError(
                    f"Far-field output references unknown mode {mode_idx}: "
                    f"{mode_source_dir}"
                )
            for source in sorted(mode_source_dir.iterdir()):
                if not source.is_file():
                    raise ValueError(f"Unknown staged far-field entry: {source}")
                destination_name = (
                    "farfield_summary.json"
                    if source.name == "summary.json"
                    else source.name
                )
                _add_output_move(
                    moves,
                    source,
                    _mode_result_dir(paths, mode_idx)
                    / "12_farfield_FFT"
                    / destination_name,
                )

    lattice_dir = staging_dir / "finite_lattice_fourier_hz"
    if lattice_dir.is_dir():
        lattice_root_files = {
            "index_r_cavity.png": (
                paths["overview_dir"] / "index_r_cavity.png"
            ),
            "index_k_1stBZ.png": (
                paths["overview_dir"] / "index_k_1stBZ.png"
            ),
            "finite_lattice_fourier_summary.csv": (
                paths["overview_dir"] / "finite_lattice_fourier_summary.csv"
            ),
        }
        for name, destination in lattice_root_files.items():
            source = lattice_dir / name
            if source.is_file():
                if name == "finite_lattice_fourier_summary.csv":
                    try:
                        lattice_frame = pd.read_csv(source)
                    except pd.errors.EmptyDataError:
                        lattice_frame = pd.DataFrame()
                    _validate_result_mode_rows(
                        lattice_frame,
                        source,
                        known_modes,
                    )
                _add_output_move(moves, source, destination)
        for mode_source_dir in sorted(lattice_dir.glob("mode_*")):
            match = re.fullmatch(r"mode_(\d+)", mode_source_dir.name)
            if not mode_source_dir.is_dir() or match is None:
                raise ValueError(f"Unknown staged lattice entry: {mode_source_dir}")
            mode_idx = int(match.group(1))
            if mode_idx not in known_modes:
                raise ValueError(
                    f"Lattice output references unknown mode {mode_idx}: "
                    f"{mode_source_dir}"
                )
            for source in sorted(mode_source_dir.iterdir()):
                if not source.is_file():
                    raise ValueError(f"Unknown staged lattice entry: {source}")
                _add_output_move(
                    moves,
                    source,
                    _mode_result_dir(paths, mode_idx)
                    / "13_lattice_fourier_Hz"
                    / source.name,
                )

    score_path = staging_dir / "mode_scores.csv"
    score_df = None
    if score_path.is_file():
        try:
            score_df = pd.read_csv(score_path)
        except pd.errors.EmptyDataError:
            score_df = pd.DataFrame()
        _validate_result_mode_rows(score_df, score_path, known_modes)
        _add_output_move(
            moves,
            score_path,
            paths["overview_dir"] / "mode_scores.csv",
        )
    mode_selection_path = staging_dir / "analysis_mode_selection.csv"
    if mode_selection_path.is_file():
        mode_selection_df = pd.read_csv(mode_selection_path)
        _validate_result_mode_rows(
            mode_selection_df,
            mode_selection_path,
            known_modes,
        )
        _add_output_move(
            moves,
            mode_selection_path,
            paths["overview_dir"] / "analysis_mode_selection.csv",
        )
    objective_path = staging_dir / "objective.json"
    if objective_path.is_file():
        _add_output_move(
            moves,
            objective_path,
            paths["config_dir"] / "objective.json",
        )

    staged_files = {path for path in staging_dir.rglob("*") if path.is_file()}
    planned_sources = {source for source, _destination in moves}
    unknown_files = sorted(staged_files.difference(planned_sources))
    if unknown_files:
        raise ValueError(f"Unknown staged finite output: {unknown_files[0]}")

    generated_targets = []
    for mode_idx in mode_indices:
        generated_targets.append(
            _mode_result_dir(paths, mode_idx)
            / "10_overview"
            / "eigenfrequency.csv"
        )
    if score_df is not None and "mode_idx" in score_df.columns:
        for mode_idx in score_df["mode_idx"].astype(int).tolist():
            generated_targets.append(
                _mode_result_dir(paths, mode_idx)
                / "10_overview"
                / "mode_score.csv"
            )

    destinations = [destination for _source, destination in moves]
    if len(destinations) != len(set(destinations)):
        raise ValueError("Finite output plan contains duplicate destinations")
    for destination in [*destinations, *generated_targets]:
        if destination.exists():
            raise FileExistsError(
                f"Refusing to overwrite finalized finite output: {destination}"
            )

    for source, destination in moves:
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.replace(destination)

    for mode_idx in mode_indices:
        mode_overview = _mode_result_dir(paths, mode_idx) / "10_overview"
        mode_overview.mkdir(parents=True, exist_ok=True)
        eigen_df.loc[eigen_df["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            mode_overview / "eigenfrequency.csv",
            index=False,
        )
    if score_df is not None and "mode_idx" in score_df.columns:
        for mode_idx in score_df["mode_idx"].astype(int).tolist():
            mode_overview = _mode_result_dir(paths, mode_idx) / "10_overview"
            mode_overview.mkdir(parents=True, exist_ok=True)
            score_df.loc[score_df["mode_idx"].astype(int).eq(mode_idx)].to_csv(
                mode_overview / "mode_score.csv",
                index=False,
            )

    _rewrite_farfield_outputs(paths, mode_indices)
    _rewrite_lattice_summary(paths)

    remaining_files = [path for path in staging_dir.rglob("*") if path.is_file()]
    if remaining_files:
        raise RuntimeError(
            f"Finite staging cleanup found an unhandled file: {remaining_files[0]}"
        )
    for directory in sorted(
        (path for path in staging_dir.rglob("*") if path.is_dir()),
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        directory.rmdir()
    staging_dir.rmdir()
    return paths


def polygons_from_records(records: list[dict[str, object]]) -> list[np.ndarray]:
    return [np.asarray(record["polygon"], dtype=float) for record in records]


def plot_limits(
    polygons: list[np.ndarray],
    extra: list[np.ndarray] | None = None,
    pad_scale: float = 0.04,
) -> tuple[float, float, float, float]:
    arrays = [*polygons, *(extra or [])]
    points = np.vstack(arrays)
    span = max(float(np.ptp(points[:, 0])), float(np.ptp(points[:, 1])), A)
    pad = pad_scale * span
    return (
        float(points[:, 0].min() - pad),
        float(points[:, 0].max() + pad),
        float(points[:, 1].min() - pad),
        float(points[:, 1].max() + pad),
    )


def draw_cells(ax, points, *, facecolor: str, edgecolor: str, alpha: float, zorder: int, label: str) -> None:
    for idx, point in enumerate(points):
        ax.add_patch(Polygon(
            cell_polygon(point.x, point.y),
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=0.32,
            alpha=alpha,
            zorder=zorder,
            label=label if idx == 0 else "_nolegend_",
        ))


def draw_record_cells(
    ax,
    records: list[dict[str, object]],
    *,
    bulk_facecolor: str,
    cladding_facecolor: str,
    bulk_edgecolor: str,
    cladding_edgecolor: str,
    alpha: float,
    zorder: int,
    label_prefix: str,
) -> None:
    seen = set()
    for record in records:
        region = str(record["region"])
        facecolor = bulk_facecolor if region == "bulk" else cladding_facecolor
        edgecolor = bulk_edgecolor if region == "bulk" else cladding_edgecolor
        label = f"{label_prefix} {region} cells"
        if region in seen:
            label = "_nolegend_"
        seen.add(region)
        ax.add_patch(Polygon(
            record["polygon"],
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=0.32,
            alpha=alpha,
            zorder=zorder,
            label=label,
        ))


def draw_record_holes(
    ax,
    records: list[dict[str, object]],
    *,
    bulk_color: str,
    cladding_color: str,
    edgecolor: str,
    alpha: float,
    zorder: int,
    label_prefix: str,
) -> None:
    seen = set()
    for record in records:
        region = str(record["region"])
        color = bulk_color if region == "bulk" else cladding_color
        label = f"{label_prefix} {region} holes"
        if region in seen:
            label = "_nolegend_"
        seen.add(region)
        ax.add_patch(Polygon(
            record["polygon"],
            closed=True,
            facecolor=color,
            edgecolor=edgecolor,
            linewidth=0.35,
            alpha=alpha,
            zorder=zorder,
            label=label,
        ))


def draw_lattice_centers(ax, points: list[LatticePoint], *, color: str, label: str, zorder: int) -> None:
    centers = np.array([[point.x, point.y] for point in points], dtype=float)
    if len(centers) == 0:
        return
    ax.scatter(
        centers[:, 0],
        centers[:, 1],
        s=18,
        facecolors="none",
        edgecolors=color,
        linewidths=0.9,
        alpha=0.74,
        label=label,
        zorder=zorder,
    )


def draw_cladding_shift_arrows(ax, points: list[LatticePoint], *, label: str, zorder: int) -> None:
    centers = []
    shifts = []
    for point in points:
        shift = cladding_inward_shift(point)
        if np.linalg.norm(shift) <= 1e-12:
            continue
        centers.append([point.x, point.y])
        shifts.append(shift)
    if not shifts:
        return
    centers = np.asarray(centers, dtype=float)
    shifts = SHIFT_ARROW_PLOT_SCALE * np.asarray(shifts, dtype=float)
    ax.quiver(
        centers[:, 0],
        centers[:, 1],
        shifts[:, 0],
        shifts[:, 1],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        color="#b45309",
        width=0.0018,
        headwidth=2.4,
        headlength=3.0,
        headaxislength=2.7,
        alpha=0.94,
        label=label,
        zorder=zorder,
    )


def ideal_bulk_hex_vertices() -> np.ndarray:
    axial_vertices = np.array([
        [BULK_RADIUS, 0],
        [0, BULK_RADIUS],
        [-BULK_RADIUS, BULK_RADIUS],
        [-BULK_RADIUS, 0],
        [0, -BULK_RADIUS],
        [BULK_RADIUS, -BULK_RADIUS],
    ], dtype=float)
    x = A * (axial_vertices[:, 0] + 0.5 * axial_vertices[:, 1])
    y = A * np.sqrt(3.0) * axial_vertices[:, 1] / 2.0
    return np.column_stack([x, y])


def normal_fan_feature(point: LatticePoint) -> tuple[str, int]:
    p = np.array([point.x, point.y], dtype=float)
    vertices = ideal_bulk_hex_vertices()
    best_distance = np.inf
    best_kind = "side"
    best_idx = 0
    for idx, start in enumerate(vertices):
        end = vertices[(idx + 1) % len(vertices)]
        edge = end - start
        t = float(np.dot(p - start, edge) / np.dot(edge, edge))
        t_clamped = float(np.clip(t, 0.0, 1.0))
        projection = start + t_clamped * edge
        distance = float(np.sum((p - projection) ** 2))
        if distance < best_distance:
            best_distance = distance
            if 1e-10 < t_clamped < 1.0 - 1e-10:
                best_kind = "side"
                best_idx = idx
            else:
                best_kind = "corner"
                best_idx = idx if t_clamped <= 1e-10 else (idx + 1) % len(vertices)
    return best_kind, best_idx


def normal_fan_side_inward_normal(side_idx: int) -> np.ndarray:
    vertices = ideal_bulk_hex_vertices()
    start = vertices[side_idx]
    end = vertices[(side_idx + 1) % len(vertices)]
    edge = end - start
    normal = np.array([-edge[1], edge[0]], dtype=float)
    return normal / np.linalg.norm(normal)


def normal_fan_direction(kind: str, idx: int) -> np.ndarray:
    if kind == "side":
        return normal_fan_side_inward_normal(idx)
    if kind == "corner":
        previous_normal = normal_fan_side_inward_normal((idx - 1) % 6)
        next_normal = normal_fan_side_inward_normal(idx)
        direction = previous_normal + next_normal
        return direction / np.linalg.norm(direction)
    raise ValueError(f"Unknown normal-fan feature kind: {kind}")


def cladding_shift_factor(kind: str, idx: int) -> float:
    group = cladding_shift_group(kind, idx)
    if group == "x":
        return CLADDING_X_SHIFT_FACTOR
    return CLADDING_Y_SHIFT_FACTOR


def normal_fan_shift_magnitude(
    kind: str,
    *,
    target_factor: float,
    envelope: float = 1.0,
) -> float:
    base_shift = float(target_factor) * A
    if kind == "side":
        return envelope * base_shift
    if kind == "corner":
        return envelope * 2.0 * base_shift / np.sqrt(3.0)
    raise ValueError(f"Unknown normal-fan feature kind: {kind}")


def normal_fan_shift(kind: str, idx: int, *, envelope: float = 1.0) -> np.ndarray:
    target_factor = cladding_shift_factor(kind, idx)
    if target_factor == 0.0 or envelope == 0.0:
        return np.zeros(2, dtype=float)
    return normal_fan_shift_magnitude(
        kind,
        target_factor=target_factor,
        envelope=envelope,
    ) * normal_fan_direction(kind, idx)


def closed_line(vertices: np.ndarray) -> np.ndarray:
    return np.vstack([vertices, vertices[0]])


def limits_from_polygons(polygons: list[np.ndarray], pad_scale: float = 0.04) -> tuple[float, float, float, float]:
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


def draw_simulation_boundary(ax, boundary: np.ndarray, *, color: str, linewidth: float, label: str, zorder: int) -> None:
    ax.add_patch(Polygon(
        boundary,
        closed=True,
        fill=False,
        edgecolor=color,
        linewidth=linewidth,
        label=label,
        zorder=zorder,
    ))


def draw_unit_cells(ax, points: list[LatticePoint], *, facecolor: str, edgecolor: str, alpha: float, zorder: int, label: str) -> None:
    if not points:
        return
    for point in points:
        ax.add_patch(Polygon(cell_polygon(point.x, point.y), closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=0.35, alpha=alpha, zorder=zorder))
    coords = np.array([[point.x, point.y] for point in points], dtype=float)
    ax.scatter(coords[:, 0], coords[:, 1], s=5, c=edgecolor, linewidths=0, alpha=0.7, label=label, zorder=zorder + 1)


def save_full_lattice_normal_fan_plot(path: Path, metadata: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        *[cell_polygon(point.x, point.y) for point in metadata["bulk_points"]],
        *[cell_polygon(point.x, point.y) for point in metadata["cladding_points"]],
    ]
    vertices = ideal_bulk_hex_vertices()
    limits = plot_limits(cells, [vertices])
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)

    draw_cells(ax, metadata["bulk_points"], facecolor="#e0f2fe", edgecolor="#38bdf8", alpha=0.42, zorder=1, label="bulk cells")
    seen = set()
    for point in metadata["cladding_points"]:
        kind, idx = normal_fan_feature(point)
        group = cladding_shift_group(kind, idx)
        if group == "x":
            facecolor = "#dbeafe"
            edgecolor = "#2563eb"
            label = "x-shift group cells"
        else:
            facecolor = "#fed7aa"
            edgecolor = "#ea580c"
            label = "y-shift group cells"
        if group in seen:
            label = "_nolegend_"
        seen.add(group)
        ax.add_patch(Polygon(
            cell_polygon(point.x, point.y),
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=0.34,
            alpha=0.52,
            zorder=2,
            label=label,
        ))

    draw_lattice_centers(ax, metadata["bulk_points"], color="#0369a1", label="bulk lattice centers", zorder=5)
    draw_lattice_centers(ax, metadata["cladding_points"], color="#111827", label="cladding lattice centers", zorder=5)
    ideal_closed = closed_line(vertices)
    ax.plot(ideal_closed[:, 0], ideal_closed[:, 1], color="#111827", linewidth=1.05, label="ideal bulk hexagon", zorder=6)
    inner_closed = closed_line(metadata["inner_boundary"])
    ax.plot(inner_closed[:, 0], inner_closed[:, 1], color="#64748b", linewidth=0.75, linestyle=":", label="bulk cell-union edge", zorder=6)

    centroid = np.mean(vertices, axis=0)
    for idx, start in enumerate(vertices):
        end = vertices[(idx + 1) % len(vertices)]
        midpoint = 0.5 * (start + end)
        side_label_pos = midpoint + 0.42 * (midpoint - centroid) / np.linalg.norm(midpoint - centroid)
        corner_label_pos = start + 0.50 * (start - centroid) / np.linalg.norm(start - centroid)
        side_color = (
            "#2563eb"
            if cladding_shift_group("side", idx) == "x"
            else "#ea580c"
        )
        corner_color = (
            "#2563eb"
            if cladding_shift_group("corner", idx) == "x"
            else "#ea580c"
        )
        ax.text(side_label_pos[0], side_label_pos[1], f"S{idx}", ha="center", va="center", fontsize=9.0, color=side_color, weight="bold", zorder=7)
        ax.text(corner_label_pos[0], corner_label_pos[1], f"C{idx}", ha="center", va="center", fontsize=9.0, color=corner_color, weight="bold", zorder=7)

    apply_limits(ax, limits)
    ax.set_title("Full finite lattice x/y cladding-shift regions")
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def save_full_lattice_ellipse_shift_plot(
    path: Path,
    metadata: dict[str, object],
) -> None:
    """Visualize the continuous angular factor used by radial ellipse shifts."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        *[cell_polygon(point.x, point.y) for point in metadata["bulk_points"]],
        *[
            cell_polygon(point.x, point.y)
            for point in metadata["cladding_points"]
        ],
    ]
    vertices = ideal_bulk_hex_vertices()
    limits = plot_limits(cells, [vertices])
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)
    draw_cells(
        ax,
        metadata["bulk_points"],
        facecolor="#e0f2fe",
        edgecolor="#38bdf8",
        alpha=0.42,
        zorder=1,
        label="bulk cells",
    )

    x_factor = float(metadata["cladding_x_shift_factor"])
    y_factor = float(metadata["cladding_y_shift_factor"])
    scale = max(x_factor, y_factor)
    for point in metadata["cladding_points"]:
        weight = ellipse_shift_weight(
            point.x,
            point.y,
            x_factor,
            y_factor,
        )
        normalized = 0.0 if scale == 0.0 else weight / scale
        color = plt.cm.viridis(float(np.clip(normalized, 0.0, 1.0)))
        ax.add_patch(
            Polygon(
                cell_polygon(point.x, point.y),
                closed=True,
                facecolor=color,
                edgecolor="#475569",
                linewidth=0.30,
                alpha=0.68,
                zorder=2,
            )
        )

    draw_lattice_centers(
        ax,
        metadata["bulk_points"],
        color="#0369a1",
        label="bulk lattice centers",
        zorder=5,
    )
    draw_lattice_centers(
        ax,
        metadata["cladding_points"],
        color="#111827",
        label="ideal cladding centers",
        zorder=5,
    )
    ideal_closed = closed_line(vertices)
    ax.plot(
        ideal_closed[:, 0],
        ideal_closed[:, 1],
        color="#111827",
        linewidth=1.05,
        label="ideal bulk hexagon",
        zorder=6,
    )
    inner_closed = closed_line(metadata["inner_boundary"])
    ax.plot(
        inner_closed[:, 0],
        inner_closed[:, 1],
        color="#64748b",
        linewidth=0.75,
        linestyle=":",
        label="bulk cell-union edge",
        zorder=6,
    )
    ax.text(
        0.012,
        0.018,
        (
            r"$f(\theta)=f_x\cos^2\theta+f_y\sin^2\theta$; "
            f"x={x_factor:.6g}, y={y_factor:.6g}"
        ),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=8,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "#9ca3af",
            "alpha": 0.92,
        },
        zorder=10,
    )
    apply_limits(ax, limits)
    ax.set_title("Full finite lattice radial ellipse-shift weights")
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def save_full_lattice_plot(
    path: Path,
    records: list[dict[str, object]],
    metadata: dict[str, object],
    *,
    annotated: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    polygons = polygons_from_records(records)
    limits = plot_limits(polygons, [metadata["inner_boundary"]])
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)
    if annotated:
        draw_cells(ax, metadata["bulk_points"], facecolor="#dbeafe", edgecolor="#3b82f6", alpha=0.50, zorder=1, label="bulk cells")
        draw_cells(ax, metadata["cladding_points"], facecolor="#dcfce7", edgecolor="#22c55e", alpha=0.32, zorder=0, label="cladding cells")
    draw_record_holes(
        ax,
        records,
        bulk_color="#1d4ed8",
        cladding_color="#15803d",
        edgecolor="#111827",
        alpha=0.78,
        zorder=4,
        label_prefix="uncut" if annotated else "_nolegend_",
    )
    if annotated:
        draw_cladding_shift_arrows(ax, metadata["cladding_points"], label=f"cladding shift arrow x{SHIFT_ARROW_PLOT_SCALE:g}", zorder=6)
    inner_closed = closed_line(metadata["inner_boundary"])
    ax.plot(
        inner_closed[:, 0],
        inner_closed[:, 1],
        color="#111827",
        linewidth=0.85,
        label="bulk cell-union edge" if annotated else "_nolegend_",
        zorder=5,
    )
    ax.axvline(STRIP_X_MIN, color="#dc2626", linewidth=1.1, linestyle="--", label="strip cut" if annotated else "_nolegend_", zorder=7)
    ax.axvline(STRIP_X_MAX, color="#dc2626", linewidth=1.1, linestyle="--", zorder=7)
    apply_limits(ax, limits)
    ax.set_title("Full modulated finite lattice holes" + (" (annotated)" if annotated else ""))
    if annotated:
        ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    else:
        ax.grid(False)
        ax.tick_params(labelbottom=False, labelleft=False)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def save_full_lattice_shift_plot(path: Path, metadata: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        *[cell_polygon(point.x, point.y) for point in metadata["bulk_points"]],
        *[cell_polygon(point.x, point.y) for point in metadata["cladding_points"]],
    ]
    limits = plot_limits(cells, [metadata["inner_boundary"]])
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)
    draw_cells(ax, metadata["bulk_points"], facecolor="#dbeafe", edgecolor="#3b82f6", alpha=0.46, zorder=1, label="bulk cells")
    draw_cells(ax, metadata["cladding_points"], facecolor="#dcfce7", edgecolor="#22c55e", alpha=0.30, zorder=0, label="cladding cells")
    draw_lattice_centers(ax, metadata["bulk_points"], color="#1d4ed8", label="bulk lattice centers", zorder=4)
    draw_lattice_centers(ax, metadata["cladding_points"], color="#15803d", label="cladding lattice centers", zorder=4)
    draw_cladding_shift_arrows(ax, metadata["cladding_points"], label=f"cladding shift arrow x{SHIFT_ARROW_PLOT_SCALE:g}", zorder=5)
    inner_closed = closed_line(metadata["inner_boundary"])
    ax.plot(inner_closed[:, 0], inner_closed[:, 1], color="#111827", linewidth=0.85, label="bulk cell-union edge", zorder=6)
    geometry = CladdingShiftGeometry.from_mapping(
        metadata.get("cladding_shift_geometry")
    )
    target_description = (
        "principal-axis target |shift|"
        if geometry.kind == ELLIPSE_SHIFT_GEOMETRY
        else "target side |shift|"
    )
    ax.text(
        0.012,
        0.018,
        (
            f"{target_description}: x={CLADDING_X_INWARD_SHIFT:.6g} um, "
            f"y={CLADDING_Y_INWARD_SHIFT:.6g} um; "
            f"profile = {metadata['cladding_shift_profile_label']}; "
            f"arrows shown x{SHIFT_ARROW_PLOT_SCALE:g}"
        ),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#9ca3af", "alpha": 0.92},
        zorder=10,
    )
    apply_limits(ax, limits)
    ax.set_title("Full finite lattice modulation vectors")
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def simulation_boundary_for_full_lattice_hexagon(
    cell_records_for_full: list[dict[str, object]],
) -> tuple[np.ndarray, float, float, float]:
    polygons = polygons_from_records(cell_records_for_full)
    points = np.vstack(polygons)
    finite_hex_a = hexagon_boundary_size_for_points(points, pad=A / 2.0)
    boundary = hexagon_boundary_vertices(finite_hex_a)
    finite_lx = finite_hex_a
    finite_ly = np.sqrt(3.0) * finite_hex_a / 2.0
    return boundary, finite_hex_a, finite_lx, finite_ly


def save_finite_simulation_plot(
    path: Path,
    records: list[dict[str, object]],
    cell_records_for_full: list[dict[str, object]],
    boundary: np.ndarray,
    validity_boundary: np.ndarray,
    metadata: dict[str, object],
    *,
    annotated: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    polygons = polygons_from_records(records)
    limits = plot_limits(polygons, [boundary, validity_boundary], pad_scale=0.06)
    fig, ax = plt.subplots(figsize=(8.0, 7.2), constrained_layout=True)
    if annotated:
        draw_record_cells(
            ax,
            cell_records_for_full,
            bulk_facecolor="#dbeafe",
            cladding_facecolor="#dcfce7",
            bulk_edgecolor="#3b82f6",
            cladding_edgecolor="#22c55e",
            alpha=0.32,
            zorder=1,
            label_prefix="full",
        )
    draw_record_holes(
        ax,
        records,
        bulk_color="#2563eb",
        cladding_color="#16a34a",
        edgecolor="#111827",
        alpha=0.88 if annotated else 0.76,
        zorder=4,
        label_prefix="full finite" if annotated else "_nolegend_",
    )

    ax.add_patch(Polygon(
        boundary,
        closed=True,
        fill=False,
        edgecolor="#dc2626",
        linewidth=1.15,
        label="finite simulation boundary" if annotated else "_nolegend_",
        zorder=7,
    ))
    if annotated:
        validity_closed = closed_line(validity_boundary)
        ax.plot(
            validity_closed[:, 0],
            validity_closed[:, 1],
            color="#f59e0b",
            linewidth=1.05,
            linestyle="--",
            label="validity boundary",
            zorder=8,
        )
    apply_limits(ax, limits)
    if annotated:
        ax.set_title("Full finite lattice simulation")
        info = (
            f"Lx = {metadata['finite_Lx']:.6g} um\n"
            f"Ly = {metadata['finite_Ly']:.6g} um\n"
            f"holes = {metadata['simulation_holes']}"
        )
        ax.text(
            0.012,
            0.98,
            info,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#9ca3af", "alpha": 0.92},
            zorder=10,
        )
        ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    else:
        ax.set_title("Full finite lattice simulation")
        ax.grid(False)
        ax.tick_params(labelbottom=False, labelleft=False)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def save_finite_model_geometry_figures(
    model_dir: Path,
    records: list[dict[str, object]],
    cell_records_for_full: list[dict[str, object]],
    boundary: np.ndarray,
    validity_boundary: np.ndarray,
    metadata: dict[str, object],
) -> None:
    """Write the shared, compact model-geometry figure set for finite runs."""

    model_dir.mkdir(parents=True, exist_ok=True)
    save_finite_simulation_plot(
        model_dir / "full_finite_simulation.png",
        records,
        cell_records_for_full,
        boundary,
        validity_boundary,
        metadata,
        annotated=False,
    )
    if metadata.get("finite_geometry") == "square":
        save_square_fragment_shift_plot(
            model_dir / "full_lattice_modulation_vectors.png",
            metadata,
            diagnostic=False,
        )
        geometry = CladdingShiftGeometry.from_mapping(
            metadata.get("cladding_shift_geometry")
        )
        diagnostic_name = (
            "full_lattice_ellipse_shift.png"
            if geometry.kind == ELLIPSE_SHIFT_GEOMETRY
            else "full_lattice_square_axis_regions.png"
        )
        save_square_fragment_shift_plot(
            model_dir / diagnostic_name,
            metadata,
            diagnostic=True,
        )
        return
    save_full_lattice_shift_plot(
        model_dir / "full_lattice_modulation_vectors.png",
        metadata,
    )
    geometry = CladdingShiftGeometry.from_mapping(
        metadata.get("cladding_shift_geometry")
    )
    if geometry.kind == ELLIPSE_SHIFT_GEOMETRY:
        save_full_lattice_ellipse_shift_plot(
            model_dir / "full_lattice_ellipse_shift.png",
            metadata,
        )
    else:
        save_full_lattice_normal_fan_plot(
            model_dir / "full_lattice_normal_fan_regions.png",
            metadata,
        )


def save_square_fragment_shift_plot(
    path: Path,
    metadata: dict[str, object],
    *,
    diagnostic: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plan = metadata.get("geometry_plan")
    if not isinstance(plan, FiniteGeometryPlan) or plan.spec.shape != "square":
        raise ValueError("square cell plots require a square FiniteGeometryPlan")
    cells = list(square_cell_plot_records(plan))
    boundary = np.asarray(metadata["finite_boundary"], dtype=float)
    polygons = [np.asarray(record["polygon"], dtype=float) for record in cells]
    limits = plot_limits(polygons, [boundary], pad_scale=0.04)
    geometry = CladdingShiftGeometry.from_mapping(
        metadata.get("cladding_shift_geometry")
    )
    max_shift = max(
        (
            float(np.linalg.norm(record["modulation_shift"]))
            for record in cells
            if record["region"] == "cladding"
        ),
        default=0.0,
    )
    fig, ax = plt.subplots(figsize=(9.0, 8.2), constrained_layout=True)
    labels_seen: set[str] = set()
    for record in cells:
        region = str(record["region"])
        if region == "bulk":
            facecolor, edgecolor, label = "#dbeafe", "#3b82f6", "cavity cell"
        elif diagnostic and geometry.kind == GROUPS_SHIFT_GEOMETRY:
            outside_x = bool(record["outside_x"])
            outside_y = bool(record["outside_y"])
            if outside_x and outside_y:
                facecolor, edgecolor, label = "#fde68a", "#d97706", "x+y corner cell"
            elif outside_x:
                facecolor, edgecolor, label = "#dcfce7", "#16a34a", "x-shift cell"
            else:
                facecolor, edgecolor, label = "#f3e8ff", "#9333ea", "y-shift cell"
        elif diagnostic:
            magnitude = float(np.linalg.norm(record["modulation_shift"]))
            normalized = 0.0 if max_shift == 0.0 else magnitude / max_shift
            facecolor = plt.cm.viridis(float(np.clip(normalized, 0.0, 1.0)))
            edgecolor, label = "#475569", "ellipse-shift cell"
        else:
            facecolor, edgecolor, label = "#dcfce7", "#16a34a", "cladding cell"
        display_label = label if label not in labels_seen else "_nolegend_"
        labels_seen.add(label)
        polygon = np.asarray(record["polygon"], dtype=float)
        ax.add_patch(
            Polygon(
                polygon,
                closed=True,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.25,
                alpha=0.55,
                label=display_label,
                zorder=1,
            )
        )
        if not diagnostic and region == "cladding":
            center = np.mean(polygon, axis=0)
            shift = np.asarray(record["modulation_shift"], dtype=float)
            ax.arrow(
                center[0],
                center[1],
                SHIFT_ARROW_PLOT_SCALE * shift[0],
                SHIFT_ARROW_PLOT_SCALE * shift[1],
                width=0.003,
                head_width=0.035,
                head_length=0.05,
                length_includes_head=True,
                color="#111827",
                alpha=0.72,
                zorder=4,
            )
    physical = closed_line(boundary)
    ax.plot(
        physical[:, 0],
        physical[:, 1],
        color="#dc2626",
        linewidth=1.1,
        label="physical footprint",
        zorder=6,
    )
    apply_limits(ax, limits)
    ax.set_title(
        "Square fragment shift regions"
        if diagnostic
        else f"Square fragment modulation vectors (x{SHIFT_ARROW_PLOT_SCALE:g})"
    )
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, LatticePoint):
        return value.__dict__
    if isinstance(value, FiniteGeometryPlan):
        return {
            "shape": value.spec.shape,
            "record_count": len(value.records),
            "placement_count": len(value.placements),
        }
    if isinstance(value, dict):
        return {key: to_jsonable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [to_jsonable(item) for item in value]
    return value


def write_config(path: Path, config: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(to_jsonable(config), f, indent=2)


def shift_factor_label(shift_factor: float) -> str:
    return f"{float(shift_factor):.3f}"


def case_stem(
    x_shift_factor: float,
    y_shift_factor: float | None = None,
) -> str:
    if y_shift_factor is None:
        return f"shift{shift_factor_label(x_shift_factor)}"
    return (
        f"shiftx{shift_factor_label(x_shift_factor)}"
        f"_shifty{shift_factor_label(y_shift_factor)}"
    )


def simulation_layers(holes: list[np.ndarray]):
    from comsol_workflow.simulation_utils import LayerSpec

    return [
        LayerSpec(
            height=SLAB_HEIGHT,
            holes=holes,
            refractive_index=REFRACTIVE_INDEX,
            label="patterned_slab",
        )
    ]


def mode_validity(sim_run, mode_count: int, metadata: dict[str, object] | None = None) -> list[dict[str, bool]]:
    validity_rows = []
    z_min = 0.0
    z_max = float(sim_run.model.java.param().evaluate("z_layer_0_hi", "um"))
    tolerance = float(sim_run.model.java.param().evaluate("selection_tol", "um"))
    if metadata is not None and "inner_boundary" in metadata:
        inner_boundary = np.asarray(metadata["inner_boundary"], dtype=float)
    else:
        inner_boundary = bulk_boundary(bulk_points())
    validity_boundary = offset_polygon(inner_boundary, tolerance)

    for mode_idx in range(mode_count):
        norm_e_peak = sim_run.get_3d_field_max_value(mode_idx, "ewfd.normE")
        norm_h_peak = sim_run.get_3d_field_max_value(mode_idx, "ewfd.normH")
        z_valid = (
            z_min - tolerance <= norm_e_peak["z"] <= z_max + tolerance
            and z_min - tolerance <= norm_h_peak["z"] <= z_max + tolerance
        )
        xy_valid = (
            point_inside_or_on_polygon(np.array([norm_e_peak["x"], norm_e_peak["y"]], dtype=float), validity_boundary)
            and point_inside_or_on_polygon(np.array([norm_h_peak["x"], norm_h_peak["y"]], dtype=float), validity_boundary)
        )
        validity_rows.append({
            "is_valid": bool(z_valid and xy_valid),
        })
    return validity_rows


def export_finite_results(sim_run, out_dir: Path, metadata: dict[str, object] | None = None) -> pd.DataFrame:
    from comsol_workflow.energy_recovery import interpolate_field

    eigenfrequencies = sim_run.get_eigenfrequencies()
    df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])
    df = pd.concat([df, pd.DataFrame(mode_validity(sim_run, len(df), metadata))], axis=1)
    df["mode_idx"] = np.arange(len(df), dtype=int)
    df.to_csv(out_dir / "eigenfrequencies.csv", index=False)
    valid_mode_indices = df.loc[df["is_valid"], "mode_idx"].to_numpy(dtype=int).tolist()
    print(f"Valid finite-cavity modes: {valid_mode_indices}")

    for mode_idx in valid_mode_indices:
        coords_re, re_hz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz", "center")
        coords_im, im_hz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz*(-i)", "center")
        im_hz_aligned = interpolate_field(coords_im, im_hz, coords_re)

        field_df = pd.DataFrame({
            "x": coords_re[:, 0],
            "y": coords_re[:, 1],
            "re": np.real(re_hz),
            "im": np.real(im_hz_aligned),
        })
        field_df.to_parquet(out_dir / f"{mode_idx:02d}_Hz_center.parquet")
        sim_run.export_2d_fields(mode_idx, "ewfd.Hz", "ReHz", out_dir, plane="center", export_data=False, export_image=True)
        sim_run.export_2d_fields(mode_idx, "ewfd.Hz*(-i)", "ImHz", out_dir, plane="center", export_data=False, export_image=True)
        sim_run.export_2d_fields(mode_idx, "ewfd.normE", "normE", out_dir, plane="center", export_data=False, export_image=True)
        sim_run.export_3d_fields(mode_idx, "ewfd.normE", "normE", out_dir)
        sim_run.export_3d_fields(mode_idx, "ewfd.normH", "normH", out_dir)

    return df


def finite_simulation_config():
    from comsol_workflow.simulation_utils import SimulationConfig

    return SimulationConfig(
        wavelength=WAVELENGTH,
        eigenfrequency_shift=EIGENFREQUENCY_SHIFT,
        eigenmode_count=EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=MODE_TYPE,
        air_cutplane_z="H_air",
    )


def finite_mesh_builder(layers):
    if not USE_MANUAL_MESH_CONSTRUCTION:
        return None

    from comsol_workflow.mesh_constrcution import construct_manual_mesh

    def build_mesh(jmodel, _config):
        construct_manual_mesh(jmodel, layers, mesh_auto_size=MESH_AUTO_SIZE)

    return build_mesh


def shared_case_metadata(metadata: dict[str, object], validity_boundary: np.ndarray) -> dict[str, object]:
    validity_layers = strip_validity_layer_count()
    validity_region = {
        "type": (
            "rectangular cavity-parent core plus floor(half cladding layers)"
            if FINITE_GEOMETRY == "square"
            else "bulk plus floor(half cladding layers) boundary"
        ),
        "half_cladding_layers": validity_layers,
        "boundary": validity_boundary,
        "finite_geometry": FINITE_GEOMETRY,
    }
    if FINITE_GEOMETRY == "hex":
        validity_region["shell_radius"] = BULK_RADIUS + validity_layers
    return {
        **metadata,
        **finite_geometry_shift_metadata(),
        "finite_cavity_preset": FINITE_CAVITY_PRESET,
        "finite_cavity_preset_config": FINITE_CAVITY_CONFIG,
        "shared_parameters": ACTIVE_PARAMETERS.to_metadata(),
        "strip_x_min": STRIP_X_MIN,
        "strip_x_max": STRIP_X_MAX,
        "strip_rotation_degrees": 0.0,
        "validity_region": validity_region,
        "export_dir": EXPORT_DIRNAME,
        "simulation_common": {
            "slab_height": SLAB_HEIGHT,
            "refractive_index": REFRACTIVE_INDEX,
            "wavelength": WAVELENGTH,
            "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
            "eigenmode_count": EIGENMODE_COUNT,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "use_manual_mesh_construction": USE_MANUAL_MESH_CONSTRUCTION,
            "mode_type": MODE_TYPE,
        },
    }


RESUME_METADATA_PATHS = (
    "a",
    "bulk_radius",
    "cladding_layers",
    "finite_geometry",
    "finite_geometry_schema",
    "finite_geometry_shift_formula",
    "cladding_x_shift_factor",
    "cladding_y_shift_factor",
    "cladding_shift_input_kind",
    "cladding_shift_base_factor",
    "cladding_y_over_x_shift_ratio",
    "cladding_x_inward_shift",
    "cladding_y_inward_shift",
    "cladding_shift_profile",
    "cladding_shift_profile_label",
    "cladding_shift_geometry",
    "outer_envelope",
    "bulk_start_envelope",
    "max_first_difference",
    "max_second_difference",
    "finite_cavity_preset",
    "finite_cavity_preset_config",
    "shared_parameters",
    "bulk_params",
    "cladding_params",
    "simulation_holes",
    "simulation_cells",
    "finite_Lx",
    "finite_Ly",
    "simulation_common",
    "simulation.mode",
    "simulation.boundary",
)
GROUPS_RESUME_METADATA_PATHS = (
    "cladding_shift_region_groups",
    "resolved_cladding_shift_by_group",
)
ELLIPSE_RESUME_METADATA_PATHS = (
    "cladding_shift_formula",
    "resolved_cladding_shift_by_axis",
    "resolved_cladding_shift_by_cell",
)


def _metadata_value(metadata: dict[str, object], path: str) -> object:
    value: object = metadata
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise ValueError(f"Cannot resume finite cavity: missing metadata field {path}")
        value = value[part]
    return value


def _metadata_values_match(saved: object, current: object) -> bool:
    if isinstance(saved, bool) or isinstance(current, bool):
        return type(saved) is type(current) and saved == current
    numeric_types = (int, float, np.integer, np.floating)
    if isinstance(saved, numeric_types) and isinstance(current, numeric_types):
        return bool(np.isclose(float(saved), float(current), rtol=1e-12, atol=1e-12))
    if isinstance(saved, dict) and isinstance(current, dict):
        if saved.keys() != current.keys():
            return False
        return all(_metadata_values_match(saved[key], current[key]) for key in saved)
    if isinstance(saved, list) and isinstance(current, list):
        return len(saved) == len(current) and all(
            _metadata_values_match(left, right)
            for left, right in zip(saved, current)
        )
    return saved == current


def validate_resume_metadata(
    saved: dict[str, object],
    current: dict[str, object],
) -> None:
    def with_compatibility_defaults(
        metadata: dict[str, object],
    ) -> dict[str, object]:
        normalized = dict(metadata)
        normalized.setdefault("finite_geometry", "hex")
        normalized.setdefault("finite_geometry_schema", "hex_full_cells_v1")
        normalized.setdefault(
            "cladding_shift_geometry",
            {"kind": GROUPS_SHIFT_GEOMETRY},
        )
        normalized_shift_geometry = CladdingShiftGeometry.from_mapping(
            normalized["cladding_shift_geometry"]
        )
        normalized.setdefault(
            "finite_geometry_shift_formula",
            (
                "hex_radial_ellipse_v1"
                if normalized_shift_geometry.kind == ELLIPSE_SHIFT_GEOMETRY
                else "hex_normal_fan_v1"
            ),
        )
        normalized.setdefault("cladding_shift_input_kind", "xy")
        normalized.setdefault("cladding_shift_base_factor", None)
        normalized.setdefault("cladding_y_over_x_shift_ratio", None)
        for field in ("finite_cavity_preset_config", "shared_parameters"):
            value = normalized.get(field)
            if isinstance(value, dict):
                nested = dict(value)
                nested.setdefault(
                    "cladding_shift_geometry",
                    {"kind": GROUPS_SHIFT_GEOMETRY},
                )
                nested.setdefault("finite_geometry", "hex")
                normalized[field] = nested
        return normalized

    normalized_saved = with_compatibility_defaults(saved)
    normalized_current = with_compatibility_defaults(current)
    geometry = CladdingShiftGeometry.from_mapping(
        normalized_current["cladding_shift_geometry"]
    )
    geometry_paths = (
        ELLIPSE_RESUME_METADATA_PATHS
        if geometry.kind == ELLIPSE_SHIFT_GEOMETRY
        else GROUPS_RESUME_METADATA_PATHS
    )
    shape_paths = (
        ("finite_hex_a",)
        if normalized_current["finite_geometry"] == "hex"
        else (
            "fragment_counts",
            "parent_cell_counts",
            "group_hole_indices",
            "square_layer_formula",
        )
    )
    for path in (*RESUME_METADATA_PATHS, *shape_paths, *geometry_paths):
        saved_value = _metadata_value(normalized_saved, path)
        current_value = _metadata_value(normalized_current, path)
        if not _metadata_values_match(saved_value, current_value):
            raise ValueError(
                "Cannot resume finite cavity: metadata mismatch at "
                f"{path}; saved={saved_value!r}, current={current_value!r}"
            )


def load_and_validate_resume_metadata(
    config_path: Path,
    current: dict[str, object],
) -> None:
    if not config_path.is_file():
        raise FileNotFoundError(
            f"Cannot resume finite cavity: missing checkpoint config {config_path}"
        )
    saved = json.loads(config_path.read_text(encoding="utf-8"))
    validate_resume_metadata(saved, current)


def mode_validity_for_boundary(
    sim_run,
    mode_count: int,
    validity_boundary: np.ndarray,
) -> list[dict[str, object]]:
    validity_boundary = np.asarray(validity_boundary, dtype=float)
    if len(validity_boundary) < 3:
        raise ValueError("mode validity requires a polygonal boundary")

    validity_rows = []
    z_min = 0.0
    z_max = float(sim_run.model.java.param().evaluate("z_layer_0_hi", "um"))
    tolerance = float(sim_run.model.java.param().evaluate("selection_tol", "um"))

    for mode_idx in range(mode_count):
        norm_e_peak = sim_run.get_3d_field_max_value(mode_idx, "ewfd.normE")
        norm_h_peak = sim_run.get_3d_field_max_value(mode_idx, "ewfd.normH")
        z_valid = (
            z_min - tolerance <= norm_e_peak["z"] <= z_max + tolerance
            and z_min - tolerance <= norm_h_peak["z"] <= z_max + tolerance
        )
        xy_valid = (
            point_inside_or_on_polygon(np.array([norm_e_peak["x"], norm_e_peak["y"]], dtype=float), validity_boundary)
            and point_inside_or_on_polygon(np.array([norm_h_peak["x"], norm_h_peak["y"]], dtype=float), validity_boundary)
        )
        validity_rows.append({
            "is_valid": bool(z_valid and xy_valid),
        })
    return validity_rows


def farfield_square_grid(
    finite_hex_a: float,
    grid_size: int,
    z_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not np.isfinite(finite_hex_a) or float(finite_hex_a) <= 0.0:
        raise ValueError("finite_hex_a must be positive and finite")
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("far-field grid_size must be an odd integer >= 3")
    if not np.isfinite(z_um):
        raise ValueError("far-field plane z must be finite")

    half_width = float(finite_hex_a) / 2.0
    x_axis = np.linspace(-half_width, half_width, grid_size)
    y_axis = np.linspace(-half_width, half_width, grid_size)
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")
    coordinates = np.column_stack(
        [
            xx.ravel(),
            yy.ravel(),
            np.full(xx.size, float(z_um)),
        ]
    )
    return x_axis, y_axis, coordinates


def _csv_bool(value: object) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def export_valid_air_fields(
    sim_run,
    export_dir: Path,
    eigen_df: pd.DataFrame,
    *,
    finite_hex_a: float,
    grid_size: int,
    h0_um: float | None,
    finite_geometry: str = "hex",
    footprint_width: float | None = None,
    footprint_height: float | None = None,
) -> dict[str, object]:
    export_dir.mkdir(parents=True, exist_ok=True)
    air_plane_z_um = float(
        sim_run.model.java.param().evaluate("z_slab_top + H_air", "um")
    )
    window_size = (
        float(finite_hex_a)
        if finite_geometry == "hex"
        else max(float(footprint_width), float(footprint_height))
    )
    x_axis, y_axis, coordinates = farfield_square_grid(
        window_size,
        grid_size,
        air_plane_z_um,
    )
    valid_rows = eigen_df.loc[eigen_df["is_valid"].map(_csv_bool)]
    valid_mode_indices = valid_rows["mode_idx"].to_numpy(dtype=int).tolist()
    tolerance = max(window_size, 1.0) * 1e-12
    if finite_geometry == "square":
        geometric_domain = (
            np.abs(coordinates[:, 0]) <= float(footprint_width) / 2.0 + tolerance
        ) & (
            np.abs(coordinates[:, 1]) <= float(footprint_height) / 2.0 + tolerance
        )
    else:
        root3 = float(np.sqrt(3.0))
        geometric_domain = (
            np.abs(coordinates[:, 1])
            <= root3 * float(finite_hex_a) / 4.0 + tolerance
        ) & (
            root3 * np.abs(coordinates[:, 0]) + np.abs(coordinates[:, 1])
            <= root3 * float(finite_hex_a) / 2.0 + tolerance
        )

    for mode_idx in valid_mode_indices:
        values = np.asarray(
            sim_run.get_fields_at_coordinates(
                mode_idx,
                ["ewfd.Ex", "ewfd.Ey", "ewfd.Ez"],
                coordinates,
                dataset="dset1",
            ),
            dtype=complex,
        )
        expected_shape = (3, len(coordinates))
        if values.shape != expected_shape:
            raise ValueError(
                f"Expected air electric-field shape {expected_shape}, got {values.shape}"
            )
        is_in_domain = geometric_domain & np.all(
            np.isfinite(np.real(values)) & np.isfinite(np.imag(values)),
            axis=0,
        )
        values[:, ~is_in_domain] = 0.0
        field_df = pd.DataFrame(
            {
                "x": coordinates[:, 0],
                "y": coordinates[:, 1],
                "Ex_re": np.real(values[0]),
                "Ex_im": np.imag(values[0]),
                "Ey_re": np.real(values[1]),
                "Ey_im": np.imag(values[1]),
                "Ez_re": np.real(values[2]),
                "Ez_im": np.imag(values[2]),
                "is_in_domain": is_in_domain,
            }
        )
        field_df.to_parquet(export_dir / f"{mode_idx:02d}_E_air.parquet", index=False)

    metadata = {
        "air_plane_expression": "z_slab_top + H_air",
        "air_plane_z_um": air_plane_z_um,
        "configured_h0_um": h0_um,
        "finite_geometry": finite_geometry,
        "footprint_shape": "rectangle" if finite_geometry == "square" else "hexagon",
        "footprint_Lx_um": (
            float(footprint_width) if footprint_width is not None else float(finite_hex_a)
        ),
        "footprint_Ly_um": (
            float(footprint_height)
            if footprint_height is not None
            else float(np.sqrt(3.0) * finite_hex_a / 2.0)
        ),
        "sampling_window_um": window_size,
        "domain_mask_kind": (
            "rectangle_footprint_v1"
            if finite_geometry == "square"
            else "hexagon_footprint_v1"
        ),
        "window_half_width_um": window_size / 2.0,
        "grid_size": int(grid_size),
        "dx_um": float(x_axis[1] - x_axis[0]),
        "dy_um": float(y_axis[1] - y_axis[0]),
        "valid_mode_indices": valid_mode_indices,
    }
    if finite_geometry == "hex":
        metadata["finite_hex_a_um"] = float(finite_hex_a)
    write_config(export_dir / "air_field_metadata.json", metadata)
    return metadata


def export_full_finite_results(
    sim_run,
    export_dir: Path,
    validity_boundary: np.ndarray,
    *,
    label: str,
    finite_hex_a: float | None = None,
    finite_geometry: str = "hex",
    footprint_width: float | None = None,
    footprint_height: float | None = None,
) -> pd.DataFrame:
    from comsol_workflow.energy_recovery import interpolate_field

    export_dir.mkdir(parents=True, exist_ok=True)
    eigenfrequencies = sim_run.get_eigenfrequencies()
    df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])
    validity = mode_validity_for_boundary(sim_run, len(df), validity_boundary)
    df = pd.concat([df, pd.DataFrame(validity)], axis=1)
    df["mode_idx"] = np.arange(len(df), dtype=int)
    df.to_csv(export_dir / "eigenfrequencies.csv", index=False)

    valid_mode_indices = df.loc[df["is_valid"], "mode_idx"].to_numpy(dtype=int).tolist()
    print(f"Valid {label} modes: {valid_mode_indices}")

    for mode_idx in df["mode_idx"].to_numpy(dtype=int).tolist():
        mode_row = df.loc[df["mode_idx"] == mode_idx].iloc[0]
        frequency_thz = float(mode_row["re"])
        quality_factor = float(mode_row["q"])
        coords_re, re_hz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz", "center")
        coords_im, im_hz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz*(-i)", "center")
        im_hz_aligned = interpolate_field(coords_im, im_hz, coords_re)

        field_df = pd.DataFrame({
            "x": coords_re[:, 0],
            "y": coords_re[:, 1],
            "re": np.real(re_hz),
            "im": np.real(im_hz_aligned),
        })
        field_df.to_parquet(export_dir / f"{mode_idx:02d}_Hz_center.parquet")
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Hz_Re_2d.png",
            coords_re,
            np.real(re_hz),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Re(Hz)",
            unit_label="A/m",
            color_scale_mode="linearsymmetric",
            dpi=DPI,
        )
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Hz_Im_2d.png",
            coords_re,
            np.real(im_hz_aligned),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Im(Hz)",
            unit_label="A/m",
            color_scale_mode="linearsymmetric",
            dpi=DPI,
        )
        wem_coordinates, wem = sim_run.get_2d_fields(
            mode_idx,
            "ewfd.Wav",
            "center",
        )
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Wem_2d.png",
            wem_coordinates,
            np.real(wem),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Wem",
            unit_label="J/m³",
            color_scale_mode="linear",
            dpi=DPI,
        )

    has_footprint = finite_hex_a is not None or (
        footprint_width is not None and footprint_height is not None
    )
    if RUN_FARFIELD_FFT and has_footprint:
        export_valid_air_fields(
            sim_run,
            export_dir,
            df,
            finite_hex_a=(
                float(finite_hex_a)
                if finite_hex_a is not None
                else max(float(footprint_width), float(footprint_height))
            ),
            grid_size=FARFIELD_GRID_SIZE,
            h0_um=FARFIELD_H0_UM,
            finite_geometry=finite_geometry,
            footprint_width=footprint_width,
            footprint_height=footprint_height,
        )

    return df


def attach_saved_model(sim_run, model_path: Path) -> None:
    sim_run.client.remove(sim_run.model)
    sim_run.model = sim_run.client.load(str(model_path))
    sim_run._plotting_initialized = False
    sim_run.plane_datasets = {
        "center": "cpl1",
        "air": "cpl2",
        "yz": "cpl3",
        "xz": "cpl4",
    }


def attach_saved_solution(sim_run, solved_path: Path) -> None:
    attach_saved_model(sim_run, solved_path)


def _java_tags(container) -> list[str]:
    return [str(tag) for tag in container.tags()]


def validate_geometry_checkpoint_model(sim_run) -> None:
    jmodel = sim_run.model.java
    try:
        component = jmodel.component("comp1")
        geometry_tags = _java_tags(component.geom("geom1").feature())
        mesh_tags = _java_tags(component.mesh("mesh1").feature())
        study_tags = _java_tags(jmodel.study())
        solution_tags = _java_tags(jmodel.sol())
        dataset_tags = _java_tags(jmodel.result().dataset())
    except Exception as exc:
        raise ValueError(
            "Cannot resume finite cavity: checkpoint is missing required "
            "COMSOL component/geometry/mesh/study structure"
        ) from exc

    if not geometry_tags:
        raise ValueError(
            "Cannot resume finite cavity: geometry feature sequence is empty"
        )
    if study_tags != ["std1"]:
        raise ValueError(
            "Cannot resume finite cavity: expected study tags ['std1'], "
            f"got {study_tags!r}"
        )
    if mesh_tags != ["size"]:
        raise ValueError(
            "Cannot resume finite cavity: expected geometry-only mesh feature tags "
            f"['size'], got {mesh_tags!r}"
        )
    if solution_tags:
        raise ValueError(
            "Cannot resume finite cavity: expected no solution tags, "
            f"got {solution_tags!r}"
        )
    if dataset_tags:
        raise ValueError(
            "Cannot resume finite cavity: expected no dataset tags, "
            f"got {dataset_tags!r}"
        )


def set_cladding_inward_shift_factors(
    x_shift_factor: float,
    y_shift_factor: float,
) -> None:
    global CLADDING_X_INWARD_SHIFT
    global CLADDING_X_SHIFT_FACTOR
    global CLADDING_Y_INWARD_SHIFT
    global CLADDING_Y_SHIFT_FACTOR

    resolved_x_factor = float(x_shift_factor)
    resolved_y_factor = float(y_shift_factor)
    if not np.isfinite(resolved_x_factor) or resolved_x_factor < 0.0:
        raise ValueError("cladding_x_shift_factor must be finite and non-negative")
    if not np.isfinite(resolved_y_factor) or resolved_y_factor < 0.0:
        raise ValueError("cladding_y_shift_factor must be finite and non-negative")
    CLADDING_X_SHIFT_FACTOR = resolved_x_factor
    CLADDING_Y_SHIFT_FACTOR = resolved_y_factor
    CLADDING_X_INWARD_SHIFT = CLADDING_X_SHIFT_FACTOR * A
    CLADDING_Y_INWARD_SHIFT = CLADDING_Y_SHIFT_FACTOR * A


def set_cladding_inward_shift_factor(shift_factor: float) -> None:
    """Compatibility helper that applies one legacy factor to both groups."""
    set_cladding_inward_shift_factors(shift_factor, shift_factor)


def resolve_finite_analysis_mode_indices(
    staging_dir: Path,
    export_dir: Path,
    *,
    target_frequency_thz: float | None = None,
) -> list[int] | None:
    """Return the mode closest to target frequency, or None when disabled."""
    if not FINITE_SINGLE_MODE_ANALYSIS_ENABLED:
        return None
    resolved_target_frequency_thz = (
        FINITE_SINGLE_MODE_TARGET_FREQUENCY_THZ
        if target_frequency_thz is None
        else float(target_frequency_thz)
    )
    selected_mode, _selection = select_finite_mode_closest_to_frequency(
        export_dir,
        staging_dir / "analysis_mode_selection.csv",
        target_frequency_thz=resolved_target_frequency_thz,
        candidate_mode_indices=FINITE_LATTICE_FOURIER_MODE_INDICES,
    )
    return [selected_mode]


def run_finite_full_case(
    out_dir: Path,
    records: list[dict[str, object]],
    metadata: dict[str, object],
    validity_boundary: np.ndarray,
) -> dict[str, Path]:
    from comsol_workflow.simulation_utils import BoundarySpec, SimulationRun

    if RUN_FARFIELD_FFT:
        validate_farfield_config(
            FARFIELD_GRID_SIZE,
            FARFIELD_FFT_SIZE,
            FARFIELD_NA,
            FARFIELD_H0_UM,
        )

    case_dir = out_dir
    output_paths = prepare_finite_case_output(
        case_dir,
        model_filename="finite_cavity.mph",
        resume=RESUME_FROM_EXISTING_MPH,
    )
    export_dir = output_paths["export_dir"]

    cell_records_for_full = full_cell_records(metadata)
    if metadata.get("finite_geometry") == "square":
        boundary = np.asarray(metadata["finite_boundary"], dtype=float)
        finite_hex_a = None
        finite_lx = float(metadata["finite_Lx"])
        finite_ly = float(metadata["finite_Ly"])
    else:
        boundary, finite_hex_a, finite_lx, finite_ly = simulation_boundary_for_full_lattice_hexagon(cell_records_for_full)
    case_metadata = {
        **shared_case_metadata(metadata, validity_boundary),
        "case": "finite_cavity",
        "source": "full finite lattice before strip clipping",
        "simulation_holes": len(records),
        "simulation_cells": len(cell_records_for_full),
        "finite_Lx": finite_lx,
        "finite_Ly": finite_ly,
        "finite_boundary": boundary,
        "simulation": {
            "mode": "finite_size",
            "boundary": (
                "rectangle"
                if metadata.get("finite_geometry") == "square"
                else "hexagon_boundary"
            ),
            "periodic_axis": None,
            "kx_only": False,
        },
        "output_layout": OUTPUT_LAYOUT_VERSION,
        "single_mode_analysis": {
            "enabled": FINITE_SINGLE_MODE_ANALYSIS_ENABLED,
            "selection": "closest_to_target_frequency",
            "target_frequency_thz": FINITE_SINGLE_MODE_TARGET_FREQUENCY_THZ,
            "candidate_mode_indices": FINITE_LATTICE_FOURIER_MODE_INDICES,
        },
        "paths": finite_case_output_metadata(output_paths),
    }
    if finite_hex_a is not None:
        case_metadata["finite_hex_a"] = finite_hex_a
    config_path = output_paths["config"]
    mph_path = output_paths["mph"]
    if RESUME_FROM_EXISTING_MPH:
        if not mph_path.is_file():
            raise FileNotFoundError(
                "Cannot resume finite cavity: missing geometry checkpoint "
                f"{mph_path}"
            )
        load_and_validate_resume_metadata(config_path, case_metadata)

    overview_dir = output_paths["overview_dir"]
    save_finite_shift_profile_diagnostics(
        overview_dir,
        case_metadata,
        dpi=DPI,
    )
    save_finite_model_geometry_figures(
        output_paths["model_dir"],
        records,
        cell_records_for_full,
        boundary,
        validity_boundary,
        case_metadata,
    )
    if not RESUME_FROM_EXISTING_MPH:
        write_config(config_path, case_metadata)

    print(
        f"Prepared finite-full geometry at {case_dir}; holes={len(records)}, "
        f"cells={len(cell_records_for_full)}, Lx={finite_lx:.6g}, Ly={finite_ly:.6g}"
    )

    progress_log_path = output_paths["progress_log"]
    print(f"COMSOL progress log: {progress_log_path}")

    with SimulationRun(config=finite_simulation_config(), progress_log_path=progress_log_path) as sim:
        layers = simulation_layers(polygons_from_records(records))
        if metadata.get("finite_geometry") == "square":
            plan = metadata["geometry_plan"]
            compiled = compile_finite_geometry_input(
                plan,
                layer_height=SLAB_HEIGHT,
                refractive_index=REFRACTIVE_INDEX,
            )
            sim_boundary = compiled.boundary
            layers = list(compiled.layers)
        else:
            sim_boundary = BoundarySpec("hexagon_boundary", a=finite_hex_a)
        if RESUME_FROM_EXISTING_MPH:
            print(f"Resuming geometry-only finite-cavity checkpoint: {mph_path}")
            attach_saved_model(sim, mph_path)
            validate_geometry_checkpoint_model(sim)
            with progress_log_path.open("a", encoding="utf-8") as handle:
                handle.write(f"Resuming geometry-only checkpoint: {mph_path}\n")
            write_config(config_path, case_metadata)
        else:
            sim.build_geometry(sim_boundary, layers, simulation_mode="finite_size", k=None)
            sim.model.save(str(mph_path))
        sim.run_simulation(mesh_save_path=mph_path, mesh_builder=finite_mesh_builder(layers))
        sim.model.save(str(mph_path))
        export_full_finite_results(
            sim,
            export_dir,
            validity_boundary,
            label="finite_cavity",
            finite_hex_a=finite_hex_a,
            finite_geometry=str(metadata.get("finite_geometry", "hex")),
            footprint_width=finite_lx,
            footprint_height=finite_ly,
        )
    analysis_mode_indices = resolve_finite_analysis_mode_indices(
        output_paths["staging_dir"],
        export_dir,
    )
    lattice_mode_indices = (
        analysis_mode_indices
        if FINITE_SINGLE_MODE_ANALYSIS_ENABLED
        else FINITE_LATTICE_FOURIER_MODE_INDICES
    )
    if RUN_FARFIELD_FFT:
        run_farfield_fft_case(
            output_paths["staging_dir"],
            export_dir,
            grid_size=FARFIELD_GRID_SIZE,
            fft_size=FARFIELD_FFT_SIZE,
            na=FARFIELD_NA,
            h0_um=FARFIELD_H0_UM,
            period_um=A,
            mode_indices=analysis_mode_indices,
        )
    if RUN_FINITE_LATTICE_FOURIER_POSTPROCESS:
        run_exported_finite_lattice_fourier_postprocess(
            output_paths["staging_dir"],
            export_dir,
            cells=list(metadata["bulk_points"]),
            cladding_cells=list(metadata["cladding_points"]),
            bulk_radius=BULK_RADIUS,
            period=A,
            rho_grid_size=FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            top_k=FINITE_LATTICE_FOURIER_TOP_K,
            mode_indices=lattice_mode_indices,
            include_unselected_valid_mode_maps=(
                not FINITE_SINGLE_MODE_ANALYSIS_ENABLED
            ),
            dpi=DPI,
            transform_kind=(
                "arbitrary_square_cavity_v1"
                if metadata.get("finite_geometry") == "square"
                else "hex_cyclic_quotient_v1"
            ),
        )
        score_exported_finite_modes(
            output_paths["staging_dir"],
            export_dirname=EXPORT_DIRNAME,
        )
    else:
        run_exported_finite_unit_cell_intensity_maps(
            output_paths["staging_dir"],
            export_dir,
            cells=list(metadata["bulk_points"]),
            cladding_cells=list(metadata["cladding_points"]),
            period=A,
            rho_grid_size=FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            mode_indices=analysis_mode_indices,
            dpi=DPI,
        )
    return finalize_finite_case_output(
        case_dir,
        model_filename="finite_cavity.mph",
    )


def run_finite_cavity_shift(
    x_shift_factor: float,
    y_shift_factor: float | None = None,
) -> dict[str, object]:
    directional_case = y_shift_factor is not None
    resolved_y_factor = (
        float(x_shift_factor)
        if y_shift_factor is None
        else float(y_shift_factor)
    )
    set_cladding_inward_shift_factors(x_shift_factor, resolved_y_factor)
    out_dir = OUT_DIR / case_stem(
        x_shift_factor,
        resolved_y_factor if directional_case else None,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    print(
        "Running finite cavity cladding shift: "
        f"geometry={CLADDING_SHIFT_GEOMETRY.kind}, "
        f"x={shift_factor_label(x_shift_factor)}*A "
        f"({CLADDING_X_INWARD_SHIFT:.6g} um), "
        f"y={shift_factor_label(resolved_y_factor)}*A "
        f"({CLADDING_Y_INWARD_SHIFT:.6g} um); "
        f"y/x={current_cladding_y_over_x_shift_ratio()}; "
        f"profile={CLADDING_SHIFT_PROFILE.label} -> {out_dir}"
    )

    records, metadata = build_full_lattice_holes()
    for message in shift_profile_preflight_messages(metadata):
        print(message)
    validity_boundary = strip_validity_boundary(metadata)
    output_paths = run_finite_full_case(
        out_dir,
        records,
        metadata,
        validity_boundary,
    )

    return {
        "finite_geometry": FINITE_GEOMETRY,
        "cladding_x_shift_factor": CLADDING_X_SHIFT_FACTOR,
        "cladding_y_shift_factor": CLADDING_Y_SHIFT_FACTOR,
        "cladding_shift_input_kind": CLADDING_SHIFT_INPUT_KIND,
        "cladding_shift_base_factor": current_cladding_shift_base_factor(),
        "cladding_y_over_x_shift_ratio": (
            current_cladding_y_over_x_shift_ratio()
        ),
        "cladding_x_inward_shift": CLADDING_X_INWARD_SHIFT,
        "cladding_y_inward_shift": CLADDING_Y_INWARD_SHIFT,
        **finite_geometry_shift_metadata(),
        "out_dir": str(out_dir),
        "mph": str(output_paths["mph"]),
        "mode_output_pattern": str(
            output_paths["results_dir"] / "mode{mode_idx}"
        ),
        "output_layout": OUTPUT_LAYOUT_VERSION,
        "model_dir": str(output_paths["model_dir"]),
        "overview_dir": str(output_paths["overview_dir"]),
        "logs_dir": str(output_paths["logs_dir"]),
        "config_dir": str(output_paths["config_dir"]),
        "finite_cavity_preset": FINITE_CAVITY_PRESET,
        "finite_cavity_preset_config": FINITE_CAVITY_CONFIG,
        "bulk_radius": BULK_RADIUS,
        "cladding_layers": CLADDING_LAYERS,
        "mesh_auto_size": MESH_AUTO_SIZE,
        "eigenmode_count": EIGENMODE_COUNT,
        "case": "finite_cavity",
        "use_manual_mesh_construction": USE_MANUAL_MESH_CONSTRUCTION,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = [
        run_finite_cavity_shift(x_shift_factor, y_shift_factor)
        for x_shift_factor, y_shift_factor in (
            ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
        )
    ]

    write_config(
        OUT_DIR / "run_summary.json",
        {
            "finite_cavity_preset": FINITE_CAVITY_PRESET,
            "finite_cavity_preset_config": FINITE_CAVITY_CONFIG,
            "runs": summaries,
            "output_layout": OUTPUT_LAYOUT_VERSION,
        },
    )
    print(f"Wrote {len(summaries)} finite cavity geometry previews to {OUT_DIR}")
    for item in summaries:
        print(f"  {item['out_dir']}")


if __name__ == "__main__":
    main()
