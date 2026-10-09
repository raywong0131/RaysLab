#!/usr/bin/env python3
"""Build 0- or 60-degree finite-lattice cuts for 1D strip simulations."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis  # noqa: E402
from comsol_workflow.cladding_shift_profile import (  # noqa: E402
    relative_cladding_layer,
    save_shift_profile_diagnostics,
    shift_profile_metadata,
    shift_profile_preflight_messages,
)
from comsol_workflow.geometry_utils import create_hexagon_design  # noqa: E402
from comsol_workflow.field_plotting import save_field_plot  # noqa: E402
from comsol_workflow.hex_lattice_utils import (  # noqa: E402
    LatticePoint,
    boundary_from_lattice_cells,
    bulk_points as _bulk_points,
    cell_polygon as _cell_polygon,
    lattice_points,
    unit_cell_corners,
)
from comsol_workflow.lattice_fourier_postprocess import (  # noqa: E402
    run_finite_lattice_fourier_postprocess,
    run_strip_bulk_fourier_postprocess,
    score_finite_modes,
    score_strip_modes,
)
from comsol_workflow.polygon_utils import point_inside_or_on_polygon  # noqa: E402
from scripts.run_main.parameter_config import (  # noqa: E402
    ACTIVE_PARAMETERS,
    STRIP_1D_OUTPUT_ROOT,
)


EXPORT_DIRNAME = "simulation_exports"
OUTPUT_LAYOUT_VERSION = 3
MODEL_DIRNAME = "00_model"
RESULTS_DIRNAME = "01_results"
OVERVIEW_DIRNAME = "10_overview"
CONFIG_DIRNAME = "99_config"
STAGING_DIRNAME = ".staging"
MODE_OVERVIEW_DIRNAME = "10_overview"
MODE_EXPORT_DIRNAME = "11_simulation_exports"
MODE_FOURIER_DIRNAME = "13_strip_bulk_fourier_Hz"
STRIP_FOURIER_STAGING_DIRNAME = "strip_bulk_fourier_hz"
GAMMA_SUBSPACE_P_VALID_THRESHOLD = 0.9
WEM_CUTLINE_SAMPLE_COUNT = 1001
OUTPUT_ROOT = STRIP_1D_OUTPUT_ROOT
RUN_STRIP_CUT = True
RUN_FINITE_FULL = False
# None selects the default policy: automatically analyze scans with >3 shifts.
# Set explicitly to True/False only when a run needs to override that policy.
RUN_GAMMA_TREND_AFTER_SCAN: bool | None = None
RUN_DIAGNOSTIC_PLOTS = False
REUSE_COMPLETED_CASES = True
DPI = 220
A = 0.82
R_0 = A / 3.0
B_0 = 0.23
CAVITY_LAYERS = ACTIVE_PARAMETERS.cavity_layers
BULK_RADIUS = CAVITY_LAYERS - 1
CLADDING_LAYERS = ACTIVE_PARAMETERS.cladding_layers
OUT_DIR = OUTPUT_ROOT / ACTIVE_PARAMETERS.structure_series_label("strip1d")
OUTER_EMPTY_PAD_PERIODS = 0.5
CLADDING_INWARD_SHIFT_FACTORS = list(ACTIVE_PARAMETERS.cladding_shift_factors)
CLADDING_SHIFT_PROFILE = ACTIVE_PARAMETERS.cladding_shift_profile
CLADDING_INWARD_SHIFT_FACTOR = CLADDING_INWARD_SHIFT_FACTORS[0]
CLADDING_INWARD_SHIFT = CLADDING_INWARD_SHIFT_FACTOR * A
SHIFT_ARROW_PLOT_SCALE = 4.0
SLAB_HEIGHT = "200 [nm]"
REFRACTIVE_INDEX = "3.3"
EIGENMODE_COUNT = ACTIVE_PARAMETERS.strip_eigenmode_count
MESH_AUTO_SIZE = ACTIVE_PARAMETERS.mesh_size
MODE_TYPE = "TE"
WAVELENGTH = "1550 [nm]"
EIGENFREQUENCY_SHIFT = ACTIVE_PARAMETERS.eigenfrequency_shift
STRIP_KX = ACTIVE_PARAMETERS.strip_kx
ALLOWED_STRIP_CUT_ANGLES = (0.0, 60.0)
STRIP_CUT_ANGLE = ACTIVE_PARAMETERS.strip_cut_angle_deg
STRIP_FOURIER_RHO_GRID_SIZE = 41
STRIP_FOURIER_TOP_K = 6
STRIP_FOURIER_MODE_INDICES: list[int] | None = None
FINITE_FOURIER_RHO_GRID_SIZE = STRIP_FOURIER_RHO_GRID_SIZE
FINITE_FOURIER_TOP_K = STRIP_FOURIER_TOP_K
FINITE_FOURIER_MODE_INDICES: list[int] | None = None


def integer_parameter_label(value: int | float) -> str:
    numeric = float(value)
    if not numeric.is_integer():
        raise ValueError(f"Expected an integer-valued parameter; got {value}")
    return str(int(numeric))


def shift_factor_label(value: float) -> str:
    """Use one stable three-decimal shift label across finite workflows."""
    return f"{float(value):.3f}"


def strip_case_stem(shift_factor: float) -> str:
    return f"shift{shift_factor_label(shift_factor)}"


def strip_cut_dirname(angle_degrees: float | None = None) -> str:
    angle = STRIP_CUT_ANGLE if angle_degrees is None else angle_degrees
    angle_int = int(integer_parameter_label(validate_strip_cut_angle(angle)))
    return f"cut_{angle_int:02d}"


def strip_cut_output_dir() -> Path:
    return OUT_DIR / strip_cut_dirname()


def output_dir_for_cladding_shift_factor(shift_factor: float) -> Path:
    return strip_cut_output_dir() / strip_case_stem(shift_factor)


BULK_PARAMS = ACTIVE_PARAMETERS.cavity.to_fourier_params("cavity_p_bic")
CLADDING_PARAMS = ACTIVE_PARAMETERS.cladding.to_fourier_params("bulk_gap_centered")

STRIP_X_MIN = -A / 2.0
STRIP_X_MAX = A / 2.0


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


def cell_polygon(center_x: float, center_y: float, period: float = A) -> np.ndarray:
    return _cell_polygon(center_x, center_y, period)


def bulk_points() -> list[LatticePoint]:
    return _bulk_points(BULK_RADIUS, A)


def bulk_boundary(points: list[LatticePoint]) -> np.ndarray:
    return boundary_from_lattice_cells(points, A)


def translate_polygon(polygon: np.ndarray, x: float, y: float) -> np.ndarray:
    return np.asarray(polygon, dtype=float) + np.array([x, y], dtype=float)


def validate_strip_cut_angle(angle_degrees: float) -> float:
    """Return the canonical supported cut angle or reject the configuration."""
    angle = float(angle_degrees)
    for allowed in ALLOWED_STRIP_CUT_ANGLES:
        if np.isclose(angle, allowed, rtol=0.0, atol=1e-12):
            return allowed
    allowed_text = ", ".join(str(angle) for angle in ALLOWED_STRIP_CUT_ANGLES)
    raise ValueError(f"STRIP_CUT_ANGLE must be one of {allowed_text}; got {angle_degrees}")


def strip_cut_basis(angle_degrees: float) -> tuple[np.ndarray, np.ndarray]:
    """Return global tangent and normal unit vectors for a strip cut."""
    angle = np.radians(validate_strip_cut_angle(angle_degrees))
    tangent = np.array([np.cos(angle), np.sin(angle)], dtype=float)
    normal = np.array([-np.sin(angle), np.cos(angle)], dtype=float)
    return tangent, normal


def strip_fourier_row_index_coefficients(angle_degrees: float) -> tuple[int, int]:
    """Map axial lattice indices to rows along the cut-local finite axis."""
    angle = validate_strip_cut_angle(angle_degrees)
    if angle == 0.0:
        return (0, 1)
    return (-1, 0)


def points_to_strip_frame(points: np.ndarray, angle_degrees: float) -> np.ndarray:
    """Express global xy coordinates in the cut-local tangent/normal frame."""
    tangent, normal = strip_cut_basis(angle_degrees)
    basis = np.column_stack([tangent, normal])
    return np.asarray(points, dtype=float) @ basis


def points_from_strip_frame(points: np.ndarray, angle_degrees: float) -> np.ndarray:
    """Express cut-local tangent/normal coordinates in the global xy frame."""
    tangent, normal = strip_cut_basis(angle_degrees)
    basis = np.column_stack([tangent, normal])
    return np.asarray(points, dtype=float) @ basis.T


def vector_to_strip_frame(vector: np.ndarray, angle_degrees: float) -> np.ndarray:
    return points_to_strip_frame(np.asarray(vector, dtype=float), angle_degrees)


def lattice_point_to_strip_frame(point: LatticePoint, angle_degrees: float) -> LatticePoint:
    local_xy = points_to_strip_frame(np.array([point.x, point.y], dtype=float), angle_degrees)
    return LatticePoint(
        i=point.i,
        j=point.j,
        shell=point.shell,
        x=float(local_xy[0]),
        y=float(local_xy[1]),
    )


def records_to_strip_frame(
    records: list[dict[str, object]],
    angle_degrees: float,
) -> list[dict[str, object]]:
    """Rigidly transform every lattice record into the strip solver frame."""
    validate_strip_cut_angle(angle_degrees)
    transformed = []
    for record in records:
        transformed_record = dict(record)
        transformed_record["polygon"] = points_to_strip_frame(record["polygon"], angle_degrees)
        if "point" in record:
            transformed_record["point"] = lattice_point_to_strip_frame(record["point"], angle_degrees)
        if "modulation_shift" in record:
            transformed_record["modulation_shift"] = vector_to_strip_frame(
                record["modulation_shift"], angle_degrees
            )
        transformed.append(transformed_record)
    return transformed


def metadata_to_strip_frame(metadata: dict[str, object], angle_degrees: float) -> dict[str, object]:
    """Transform geometry-bearing metadata into the strip solver frame."""
    transformed = dict(metadata)
    transformed["inner_boundary"] = points_to_strip_frame(metadata["inner_boundary"], angle_degrees)
    transformed["bulk_points"] = [
        lattice_point_to_strip_frame(point, angle_degrees)
        for point in metadata["bulk_points"]
    ]
    transformed["cladding_points"] = [
        lattice_point_to_strip_frame(point, angle_degrees)
        for point in metadata["cladding_points"]
    ]
    return transformed


def cladding_inward_shift(point: LatticePoint) -> np.ndarray:
    kind, idx = normal_fan_feature(point)
    layer = relative_cladding_layer(point, BULK_RADIUS)
    envelope = CLADDING_SHIFT_PROFILE.envelope(layer)
    return normal_fan_shift(kind, idx, envelope=envelope)


def closed_line(vertices: np.ndarray) -> np.ndarray:
    return np.vstack([vertices, vertices[0]])


def to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, LatticePoint):
        return value.__dict__
    if isinstance(value, dict):
        return {key: to_jsonable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [to_jsonable(item) for item in value]
    return value


def polygon_area(vertices: np.ndarray) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def dedupe_vertices(vertices: np.ndarray, tolerance: float = 1e-12) -> np.ndarray:
    if len(vertices) == 0:
        return vertices
    kept = [vertices[0]]
    for vertex in vertices[1:]:
        if np.linalg.norm(vertex - kept[-1]) > tolerance:
            kept.append(vertex)
    if len(kept) > 1 and np.linalg.norm(kept[0] - kept[-1]) <= tolerance:
        kept.pop()
    return np.asarray(kept, dtype=float)


def clip_polygon_x_leq(polygon: np.ndarray, x_max: float) -> np.ndarray:
    clipped = []
    for idx, current in enumerate(polygon):
        previous = polygon[idx - 1]
        current_inside = current[0] <= x_max + 1e-12
        previous_inside = previous[0] <= x_max + 1e-12
        dx = current[0] - previous[0]
        if current_inside != previous_inside and abs(dx) > 1e-14:
            fraction = (x_max - previous[0]) / dx
            clipped.append(previous + fraction * (current - previous))
        if current_inside:
            clipped.append(current)
    return dedupe_vertices(np.asarray(clipped, dtype=float))


def clip_polygon_x_geq(polygon: np.ndarray, x_min: float) -> np.ndarray:
    clipped = []
    for idx, current in enumerate(polygon):
        previous = polygon[idx - 1]
        current_inside = current[0] >= x_min - 1e-12
        previous_inside = previous[0] >= x_min - 1e-12
        dx = current[0] - previous[0]
        if current_inside != previous_inside and abs(dx) > 1e-14:
            fraction = (x_min - previous[0]) / dx
            clipped.append(previous + fraction * (current - previous))
        if current_inside:
            clipped.append(current)
    return dedupe_vertices(np.asarray(clipped, dtype=float))


def clip_polygon_to_vertical_strip(
    polygon: np.ndarray,
    x_min: float = STRIP_X_MIN,
    x_max: float = STRIP_X_MAX,
) -> np.ndarray | None:
    clipped = clip_polygon_x_geq(np.asarray(polygon, dtype=float), x_min)
    if len(clipped) == 0:
        return None
    clipped = clip_polygon_x_leq(clipped, x_max)
    if len(clipped) < 3:
        return None
    if abs(polygon_area(clipped)) <= 1e-14:
        return None
    return clipped


def build_full_lattice_holes() -> tuple[list[dict[str, object]], dict[str, object]]:
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
        for hole in cladding_unit_holes:
            records.append({
                "region": "cladding",
                "point": point,
                "parent_layer": point.shell + 1,
                "cladding_shift_layer": point.shell - BULK_RADIUS,
                "polygon": translate_polygon(hole, point.x + shift[0], point.y + shift[1]),
                "modulation_shift": shift,
            })

    profile_metadata = shift_profile_metadata(
        CLADDING_LAYERS,
        CLADDING_INWARD_SHIFT_FACTOR,
        A,
        CLADDING_SHIFT_PROFILE,
    )
    metadata = {
        "a": A,
        "bulk_radius": BULK_RADIUS,
        "cladding_layers": CLADDING_LAYERS,
        "outer_empty_pad_periods": OUTER_EMPTY_PAD_PERIODS,
        "cladding_inward_shift_factor": CLADDING_INWARD_SHIFT_FACTOR,
        "cladding_inward_shift": CLADDING_INWARD_SHIFT,
        "cladding_inward_shift_target": CLADDING_INWARD_SHIFT,
        "cladding_side_inward_shift": normal_fan_shift_magnitude("side"),
        "cladding_side_inward_shift_target": normal_fan_shift_magnitude("side"),
        "cladding_corner_inward_shift": normal_fan_shift_magnitude("corner"),
        "cladding_corner_inward_shift_target": normal_fan_shift_magnitude("corner"),
        **profile_metadata,
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


def strip_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    clipped_records = []
    for record in records:
        clipped = clip_polygon_to_vertical_strip(record["polygon"])
        if clipped is None:
            continue
        clipped_records.append({**record, "polygon": clipped})
    return clipped_records


def strip_records_for_cut(
    records: list[dict[str, object]],
    angle_degrees: float,
) -> list[dict[str, object]]:
    return strip_records(records_to_strip_frame(records, angle_degrees))


def cell_records(points, region: str) -> list[dict[str, object]]:
    return [
        {
            "region": region,
            "point": point,
            "polygon": cell_polygon(point.x, point.y),
        }
        for point in points
    ]


def strip_cell_records(metadata: dict[str, object]) -> list[dict[str, object]]:
    return strip_records([
        *cell_records(metadata["bulk_points"], "bulk"),
        *cell_records(metadata["cladding_points"], "cladding"),
    ])


def full_cell_records(metadata: dict[str, object]) -> list[dict[str, object]]:
    return [
        *cell_records(metadata["bulk_points"], "bulk"),
        *cell_records(metadata["cladding_points"], "cladding"),
    ]


def strip_bulk_cell_records(angle_degrees: float = STRIP_CUT_ANGLE) -> list[dict[str, object]]:
    metadata = {
        "bulk_points": [
            lattice_point_to_strip_frame(point, angle_degrees)
            for point in bulk_points()
        ],
        "cladding_points": [],
    }
    return [
        record
        for record in strip_cell_records(metadata)
        if record["region"] == "bulk"
    ]


def strip_validity_layer_count() -> int:
    return CLADDING_LAYERS // 2


def strip_validity_boundary() -> np.ndarray:
    radius = BULK_RADIUS + strip_validity_layer_count()
    return boundary_from_lattice_cells(lattice_points(radius, A), A)


def case_export_dir(out_dir: Path) -> Path:
    return out_dir / EXPORT_DIRNAME


def strip_case_dir(out_dir: Path) -> Path:
    return out_dir


def strip_case_output_paths(case_dir: Path) -> dict[str, Path]:
    case_dir = Path(case_dir)
    model_dir = case_dir / MODEL_DIRNAME
    staging_dir = case_dir / STAGING_DIRNAME
    return {
        "case_dir": case_dir,
        "model_dir": model_dir,
        "results_dir": case_dir / RESULTS_DIRNAME,
        "overview_dir": case_dir / OVERVIEW_DIRNAME,
        "config_dir": case_dir / CONFIG_DIRNAME,
        "staging_dir": staging_dir,
        "export_dir": staging_dir / EXPORT_DIRNAME,
        "mph": model_dir / "strip_1d.mph",
        "progress_log": model_dir / "comsol_progress.log",
        "config": case_dir / CONFIG_DIRNAME / "config.json",
    }


def prepare_strip_case_output(case_dir: Path) -> dict[str, Path]:
    paths = strip_case_output_paths(case_dir)
    paths["case_dir"].mkdir(parents=True, exist_ok=True)

    occupied = []
    for key in ("model_dir", "results_dir", "overview_dir", "config_dir"):
        directory = paths[key]
        if directory.exists() and any(directory.iterdir()):
            occupied.append(directory)
    if occupied:
        raise FileExistsError(
            "Refusing to overwrite existing strip output: "
            f"{occupied[0]}"
        )

    staging_dir = paths["staging_dir"]
    if staging_dir.exists():
        if any(staging_dir.iterdir()):
            raise FileExistsError(
                "Refusing to overwrite non-empty strip staging directory: "
                f"{staging_dir}"
            )
        staging_dir.rmdir()

    for key in ("model_dir", "results_dir", "overview_dir", "config_dir"):
        paths[key].mkdir(parents=True, exist_ok=True)
    staging_dir.mkdir(parents=True)
    return paths


def _strip_mode_dir(paths: dict[str, Path], mode_idx: int) -> Path:
    return paths["results_dir"] / f"mode{int(mode_idx)}"


def _validate_strip_mode_rows(frame, source: Path, known_modes: set[int]) -> None:
    if frame.empty and "mode_idx" not in frame.columns:
        return
    if "mode_idx" not in frame.columns:
        raise ValueError(f"Missing mode_idx column in {source}")
    mode_indices = frame["mode_idx"].astype(int).tolist()
    unknown = sorted(set(mode_indices).difference(known_modes))
    if unknown:
        raise ValueError(f"{source} references unknown mode {unknown[0]}")
    if len(mode_indices) != len(set(mode_indices)):
        raise ValueError(f"Duplicate mode_idx values in {source}")


def finalize_strip_case_output(case_dir: Path) -> dict[str, Path]:
    import pandas as pd

    paths = strip_case_output_paths(case_dir)
    staging_dir = paths["staging_dir"]
    export_dir = paths["export_dir"]
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.is_file():
        raise FileNotFoundError(
            f"Cannot finalize strip output without {eigen_path}"
        )

    eigen_df = pd.read_csv(eigen_path)
    if "mode_idx" not in eigen_df.columns:
        raise ValueError(f"Missing mode_idx column in {eigen_path}")
    mode_indices = eigen_df["mode_idx"].astype(int).tolist()
    if len(mode_indices) != len(set(mode_indices)):
        raise ValueError(f"Duplicate mode_idx values in {eigen_path}")
    known_modes = set(mode_indices)

    score_path = staging_dir / "mode_scores.csv"
    objective_path = staging_dir / "objective.json"
    if not score_path.is_file():
        raise FileNotFoundError(
            f"Cannot finalize strip output without {score_path}"
        )
    if not objective_path.is_file():
        raise FileNotFoundError(
            f"Cannot finalize strip output without {objective_path}"
        )
    score_df = pd.read_csv(score_path)
    _validate_strip_mode_rows(score_df, score_path, known_modes)
    if "gamma_subspace_p" not in score_df.columns:
        raise ValueError(f"Missing gamma_subspace_p column in {score_path}")
    score_df["is_valid"] = score_df["gamma_subspace_p"].gt(
        GAMMA_SUBSPACE_P_VALID_THRESHOLD
    )
    valid_modes = set(
        score_df.loc[score_df["is_valid"], "mode_idx"].astype(int).tolist()
    )

    moves: list[tuple[Path, Path]] = [
        (eigen_path, paths["overview_dir"] / "eigenfrequencies.csv"),
        (score_path, paths["overview_dir"] / "mode_scores.csv"),
        (objective_path, paths["config_dir"] / "objective.json"),
    ]
    discarded_sources: set[Path] = set()
    for source in sorted(export_dir.iterdir(), key=lambda item: item.name):
        if source == eigen_path:
            continue
        if not source.is_file():
            raise ValueError(f"Unknown staged simulation export: {source}")
        match = re.fullmatch(r"(\d+)_(.+)", source.name)
        if match is None:
            raise ValueError(f"Unknown staged simulation export: {source}")
        mode_idx = int(match.group(1))
        if mode_idx not in known_modes:
            raise ValueError(
                f"Simulation export references unknown mode {mode_idx}: {source}"
            )
        if mode_idx in valid_modes:
            moves.append((
                source,
                _strip_mode_dir(paths, mode_idx)
                / MODE_EXPORT_DIRNAME
                / match.group(2),
            ))
        else:
            discarded_sources.add(source)

    fourier_dir = staging_dir / STRIP_FOURIER_STAGING_DIRNAME
    fourier_summary = fourier_dir / "strip_bulk_fourier_summary.csv"
    if not fourier_summary.is_file():
        raise FileNotFoundError(
            f"Cannot finalize strip output without {fourier_summary}"
        )
    fourier_df = pd.read_csv(fourier_summary)
    _validate_strip_mode_rows(fourier_df, fourier_summary, known_modes)
    moves.append((
        fourier_summary,
        paths["overview_dir"] / "strip_bulk_fourier_summary.csv",
    ))
    shared_fourier_files = {
        "bulk_row_indices.png",
        "bulk_k_grid_1d.png",
    }
    for source in sorted(fourier_dir.iterdir(), key=lambda item: item.name):
        if source == fourier_summary:
            continue
        if source.is_file():
            if source.name not in shared_fourier_files:
                raise ValueError(f"Unknown staged strip Fourier output: {source}")
            moves.append((source, paths["overview_dir"] / source.name))
            continue
        match = re.fullmatch(r"mode_(\d+)", source.name)
        if not source.is_dir() or match is None:
            raise ValueError(f"Unknown staged strip Fourier entry: {source}")
        mode_idx = int(match.group(1))
        if mode_idx not in known_modes:
            raise ValueError(
                f"Strip Fourier output references unknown mode {mode_idx}: "
                f"{source}"
            )
        for mode_source in sorted(source.iterdir(), key=lambda item: item.name):
            if not mode_source.is_file():
                raise ValueError(
                    f"Unknown staged strip Fourier entry: {mode_source}"
                )
            if mode_idx in valid_modes:
                moves.append((
                    mode_source,
                    _strip_mode_dir(paths, mode_idx)
                    / MODE_FOURIER_DIRNAME
                    / mode_source.name,
                ))
            else:
                discarded_sources.add(mode_source)

    if not valid_modes.issubset(set(fourier_df["mode_idx"].astype(int))):
        missing = sorted(valid_modes.difference(
            set(fourier_df["mode_idx"].astype(int))
        ))
        raise ValueError(f"Valid strip mode missing Fourier result: {missing[0]}")
    required_exports = {
        "Hz_center.parquet",
        "Hz_Re_2d.png",
        "Hz_Im_2d.png",
        "Wem_2d.png",
        "cutline_Wem.png",
    }
    for mode_idx in valid_modes:
        staged_names = {
            source.name.split("_", 1)[1]
            for source in export_dir.glob(f"{mode_idx:02d}_*")
        }
        missing = sorted(required_exports.difference(staged_names))
        if missing:
            raise FileNotFoundError(
                f"Valid strip mode {mode_idx} missing staged export {missing[0]}"
            )

    staged_files = {path for path in staging_dir.rglob("*") if path.is_file()}
    planned_sources = {source for source, _destination in moves}
    handled_sources = planned_sources.union(discarded_sources)
    unknown_files = sorted(staged_files.difference(handled_sources))
    if unknown_files:
        raise ValueError(f"Unknown staged strip output: {unknown_files[0]}")

    generated_targets = [
        _strip_mode_dir(paths, mode_idx)
        / MODE_OVERVIEW_DIRNAME
        / "eigenfrequency.csv"
        for mode_idx in sorted(valid_modes)
    ]
    generated_targets.extend(
        _strip_mode_dir(paths, mode_idx)
        / MODE_OVERVIEW_DIRNAME
        / "mode_score.csv"
        for mode_idx in sorted(valid_modes)
    )
    destinations = [destination for _source, destination in moves]
    if len(destinations) != len(set(destinations)):
        raise ValueError("Strip output plan contains duplicate destinations")
    for destination in [*destinations, *generated_targets]:
        if destination.exists():
            raise FileExistsError(
                f"Refusing to overwrite finalized strip output: {destination}"
            )

    for source, destination in moves:
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.replace(destination)
    for source in discarded_sources:
        source.unlink()

    final_eigen_path = paths["overview_dir"] / "eigenfrequencies.csv"
    final_eigen_df = pd.read_csv(final_eigen_path)
    final_eigen_df["geometry_is_valid"] = final_eigen_df["is_valid"]
    gamma_by_mode = score_df.set_index("mode_idx")["gamma_subspace_p"]
    final_eigen_df["gamma_subspace_p"] = final_eigen_df["mode_idx"].map(
        gamma_by_mode
    )
    final_eigen_df["is_valid"] = final_eigen_df["mode_idx"].isin(valid_modes)
    final_eigen_df.to_csv(final_eigen_path, index=False)
    score_df.to_csv(paths["overview_dir"] / "mode_scores.csv", index=False)

    for mode_idx in sorted(valid_modes):
        mode_overview = (
            _strip_mode_dir(paths, mode_idx) / MODE_OVERVIEW_DIRNAME
        )
        mode_overview.mkdir(parents=True, exist_ok=True)
        final_eigen_df.loc[
            final_eigen_df["mode_idx"].astype(int).eq(mode_idx)
        ].to_csv(
            mode_overview / "eigenfrequency.csv",
            index=False,
        )
        score_df.loc[score_df["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            mode_overview / "mode_score.csv",
            index=False,
        )

    final_fourier_summary = (
        paths["overview_dir"] / "strip_bulk_fourier_summary.csv"
    )
    final_fourier_df = pd.read_csv(final_fourier_summary)
    final_fourier_df["geometry_is_valid"] = final_fourier_df["is_valid"]
    final_fourier_df["gamma_subspace_p"] = final_fourier_df["mode_idx"].map(
        gamma_by_mode
    )
    final_fourier_df["is_valid"] = final_fourier_df["mode_idx"].isin(
        valid_modes
    )
    for row_idx, row in final_fourier_df.iterrows():
        mode_idx = int(row["mode_idx"])
        final_fourier_df.at[row_idx, "output_dir"] = (
            str(_strip_mode_dir(paths, mode_idx) / MODE_FOURIER_DIRNAME)
            if mode_idx in valid_modes
            else ""
        )
    final_fourier_df.to_csv(final_fourier_summary, index=False)

    remaining_files = [path for path in staging_dir.rglob("*") if path.is_file()]
    if remaining_files:
        raise RuntimeError(
            f"Strip staging cleanup found an unhandled file: {remaining_files[0]}"
        )
    for directory in sorted(
        (path for path in staging_dir.rglob("*") if path.is_dir()),
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        directory.rmdir()
    staging_dir.rmdir()
    return paths


def finite_full_case_dir(out_dir: Path) -> Path:
    return out_dir / "finite_cavity"


def polygons_from_records(records: list[dict[str, object]]) -> list[np.ndarray]:
    return [np.asarray(record["polygon"], dtype=float) for record in records]


def plot_limits(polygons: list[np.ndarray], extra: list[np.ndarray] | None = None, pad_scale: float = 0.04) -> tuple[float, float, float, float]:
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


def apply_limits(ax, limits: tuple[float, float, float, float]) -> None:
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.16)


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


def unique_lattice_points(records: list[dict[str, object]], region: str | None = None) -> list[LatticePoint]:
    points = []
    seen = set()
    for record in records:
        if region is not None and record["region"] != region:
            continue
        point = record["point"]
        key = (point.i, point.j)
        if key in seen:
            continue
        seen.add(key)
        points.append(point)
    return points


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


def normal_fan_shift_magnitude(kind: str, *, envelope: float = 1.0) -> float:
    if kind == "side":
        return envelope * CLADDING_INWARD_SHIFT
    if kind == "corner":
        return envelope * 2.0 * CLADDING_INWARD_SHIFT / np.sqrt(3.0)
    raise ValueError(f"Unknown normal-fan feature kind: {kind}")


def normal_fan_shift(kind: str, idx: int, *, envelope: float = 1.0) -> np.ndarray:
    if CLADDING_INWARD_SHIFT == 0.0 or envelope == 0.0:
        return np.zeros(2, dtype=float)
    return normal_fan_shift_magnitude(
        kind,
        envelope=envelope,
    ) * normal_fan_direction(kind, idx)


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
        kind, _ = normal_fan_feature(point)
        if kind == "side":
            facecolor = "#dbeafe"
            edgecolor = "#2563eb"
            label = "side sector cells"
        else:
            facecolor = "#fed7aa"
            edgecolor = "#ea580c"
            label = "corner sector cells"
        if kind in seen:
            label = "_nolegend_"
        seen.add(kind)
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
        ax.text(side_label_pos[0], side_label_pos[1], f"S{idx}", ha="center", va="center", fontsize=9.0, color="#1d4ed8", weight="bold", zorder=7)
        ax.text(corner_label_pos[0], corner_label_pos[1], f"C{idx}", ha="center", va="center", fontsize=9.0, color="#c2410c", weight="bold", zorder=7)

    ax.axvline(STRIP_X_MIN, color="#dc2626", linewidth=1.1, linestyle="--", label="strip cut", zorder=7)
    ax.axvline(STRIP_X_MAX, color="#dc2626", linewidth=1.1, linestyle="--", zorder=7)
    apply_limits(ax, limits)
    ax.set_title("Full finite lattice side/corner normal-fan regions")
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
    ax.axvline(STRIP_X_MIN, color="#dc2626", linewidth=1.1, linestyle="--", label="strip cut", zorder=7)
    ax.axvline(STRIP_X_MAX, color="#dc2626", linewidth=1.1, linestyle="--", zorder=7)
    ax.text(
        0.012,
        0.018,
        (
            f"target side |shift| = {normal_fan_shift_magnitude('side'):.6g} um; "
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


def save_strip_plot(
    path: Path,
    records: list[dict[str, object]],
    cell_records_for_strip: list[dict[str, object]],
    metadata: dict[str, object],
    *,
    annotated: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    polygons = polygons_from_records(records)
    strip_boundary = np.array([
        [STRIP_X_MIN, min(float(poly[:, 1].min()) for poly in polygons)],
        [STRIP_X_MAX, min(float(poly[:, 1].min()) for poly in polygons)],
        [STRIP_X_MAX, max(float(poly[:, 1].max()) for poly in polygons)],
        [STRIP_X_MIN, max(float(poly[:, 1].max()) for poly in polygons)],
    ], dtype=float)
    limits = plot_limits(polygons, [strip_boundary], pad_scale=0.06)
    fig, ax = plt.subplots(figsize=(4.6, 10.0), constrained_layout=True)
    if annotated:
        draw_record_cells(
            ax,
            cell_records_for_strip,
            bulk_facecolor="#dbeafe",
            cladding_facecolor="#dcfce7",
            bulk_edgecolor="#3b82f6",
            cladding_edgecolor="#22c55e",
            alpha=0.42,
            zorder=1,
            label_prefix="clipped",
        )
    draw_record_holes(
        ax,
        records,
        bulk_color="#2563eb",
        cladding_color="#16a34a",
        edgecolor="#111827",
        alpha=0.86,
        zorder=4,
        label_prefix="clipped" if annotated else "_nolegend_",
    )
    ax.axvline(STRIP_X_MIN, color="#dc2626", linewidth=1.15, linestyle="--", label="x = +/- a/2" if annotated else "_nolegend_", zorder=7)
    ax.axvline(STRIP_X_MAX, color="#dc2626", linewidth=1.15, linestyle="--", zorder=7)
    apply_limits(ax, limits)
    ax.set_title("Cut-local strip holes" + (" (annotated)" if annotated else ""))
    if annotated:
        ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    else:
        ax.grid(False)
        ax.tick_params(labelbottom=False, labelleft=False)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def save_global_cut_preview(
    path: Path,
    local_records: list[dict[str, object]],
    local_boundary: np.ndarray,
    angle_degrees: float,
) -> None:
    """Show the local rectangular simulation domain in the original lattice frame."""
    path.parent.mkdir(parents=True, exist_ok=True)
    global_boundary = points_from_strip_frame(local_boundary, angle_degrees)
    global_records = [
        {
            **record,
            "polygon": points_from_strip_frame(record["polygon"], angle_degrees),
        }
        for record in local_records
    ]
    polygons = polygons_from_records(global_records)
    limits = plot_limits(polygons, [global_boundary], pad_scale=0.04)
    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)
    draw_record_holes(
        ax,
        global_records,
        bulk_color="#2563eb",
        cladding_color="#16a34a",
        edgecolor="#111827",
        alpha=0.58,
        zorder=3,
        label_prefix="global",
    )
    ax.add_patch(Polygon(
        global_boundary,
        closed=True,
        fill=False,
        edgecolor="#dc2626",
        linewidth=1.4,
        label=f"{integer_parameter_label(angle_degrees)} degree cut domain",
        zorder=7,
    ))
    tangent, normal = strip_cut_basis(angle_degrees)
    arrow_length = 2.2 * A
    ax.quiver(
        [0.0, 0.0],
        [0.0, 0.0],
        [arrow_length * tangent[0], arrow_length * normal[0]],
        [arrow_length * tangent[1], arrow_length * normal[1]],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        color=["#7c3aed", "#ea580c"],
        width=0.004,
        zorder=8,
    )
    ax.text(*(arrow_length * tangent), "periodic tangent", color="#7c3aed", fontsize=8)
    ax.text(*(arrow_length * normal), "finite normal", color="#ea580c", fontsize=8)
    apply_limits(ax, limits)
    ax.set_title("Strip cut in global lattice coordinates")
    ax.legend(loc="upper right", framealpha=0.94, fontsize=8)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def simulation_boundary_for_strip(cell_records_for_strip: list[dict[str, object]]) -> tuple[np.ndarray, float, float]:
    polygons = polygons_from_records(cell_records_for_strip)
    points = np.vstack(polygons)
    half_x = A / 2.0
    half_y = float(np.max(np.abs(points[:, 1]))) + OUTER_EMPTY_PAD_PERIODS * A
    boundary = np.array([
        [-half_x, -half_y],
        [half_x, -half_y],
        [half_x, half_y],
        [-half_x, half_y],
    ], dtype=float)
    return boundary, 2.0 * half_x, 2.0 * half_y


def simulation_boundary_for_full_lattice(cell_records_for_full: list[dict[str, object]]) -> tuple[np.ndarray, float, float]:
    polygons = polygons_from_records(cell_records_for_full)
    points = np.vstack(polygons)
    half_x = float(np.max(np.abs(points[:, 0]))) + A / 2.0
    half_y = float(np.max(np.abs(points[:, 1]))) + A / 2.0
    boundary = np.array([
        [-half_x, -half_y],
        [half_x, -half_y],
        [half_x, half_y],
        [-half_x, half_y],
    ], dtype=float)
    return boundary, 2.0 * half_x, 2.0 * half_y


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


def save_simulation_plot(
    path: Path,
    records: list[dict[str, object]],
    boundary: np.ndarray,
    metadata: dict[str, object],
    *,
    annotated: bool,
    cell_records_for_strip: list[dict[str, object]] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    polygons = polygons_from_records(records)
    limits = plot_limits(polygons, [boundary], pad_scale=0.06)
    fig, ax = plt.subplots(figsize=(4.8, 10.4), constrained_layout=True)
    if annotated and cell_records_for_strip is not None:
        draw_record_cells(
            ax,
            cell_records_for_strip,
            bulk_facecolor="#dbeafe",
            cladding_facecolor="#dcfce7",
            bulk_edgecolor="#3b82f6",
            cladding_edgecolor="#22c55e",
            alpha=0.34,
            zorder=1,
            label_prefix="simulation clipped",
        )
    draw_record_holes(
        ax,
        records,
        bulk_color="#2563eb",
        cladding_color="#16a34a",
        edgecolor="#111827",
        alpha=0.88 if annotated else 0.76,
        zorder=4,
        label_prefix="simulation" if annotated else "_nolegend_",
    )
    xmin = float(boundary[:, 0].min())
    xmax = float(boundary[:, 0].max())
    ymin = float(boundary[:, 1].min())
    ymax = float(boundary[:, 1].max())
    ax.add_patch(Rectangle(
        (xmin, ymin),
        xmax - xmin,
        ymax - ymin,
        fill=False,
        edgecolor="#dc2626",
        linewidth=1.15,
        label="strip_1d simulation boundary" if annotated else "_nolegend_",
        zorder=7,
    ))
    apply_limits(ax, limits)
    if annotated:
        ax.set_title("Cut-local strip for strip_1d simulation")
        info = (
            f"Lx = {metadata['strip_Lx']:.6g} um\n"
            f"Ly = {metadata['strip_Ly']:.6g} um\n"
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
        ax.set_title("Cut-local strip for strip_1d")
        ax.grid(False)
        ax.tick_params(labelbottom=False, labelleft=False)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


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


def write_config(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(to_jsonable(data), handle, indent=2)


def wem_cutline_at_x_zero(
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    sample_count: int = WEM_CUTLINE_SAMPLE_COUNT,
) -> tuple[np.ndarray, np.ndarray]:
    from comsol_workflow.energy_recovery import interpolate_field

    points = np.asarray(coordinates, dtype=float)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError("Wem cutline coordinates must have shape (N, >=2)")
    points = points[:, :2]
    wem = np.real(np.asarray(values)).reshape(-1)
    if len(points) != len(wem):
        raise ValueError("Wem cutline coordinates and values must have equal length")
    if sample_count < 2:
        raise ValueError("Wem cutline sample_count must be at least 2")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(wem)):
        raise ValueError("Wem cutline input must contain only finite values")
    if float(np.min(points[:, 0])) > 0.0 or float(np.max(points[:, 0])) < 0.0:
        raise ValueError("Wem field does not span x=0")

    y = np.linspace(
        float(np.min(points[:, 1])),
        float(np.max(points[:, 1])),
        int(sample_count),
    )
    query = np.column_stack([np.zeros_like(y), y])
    cutline = np.real(interpolate_field(points, wem, query))
    scale = float(np.max(np.abs(cutline)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("Wem cutline cannot be normalized from a zero field")
    return y, cutline / scale


def save_wem_cutline_plot(
    path: Path,
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    dpi: int = DPI,
) -> tuple[np.ndarray, np.ndarray]:
    y, normalized_wem = wem_cutline_at_x_zero(coordinates, values)
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    ax.plot(y, normalized_wem, color="#1f2937", linewidth=1.4)
    ax.set_title("Wem cutline at x = 0")
    ax.set_xlabel("y (µm)")
    ax.set_ylabel("Normalized Wem")
    ax.set_xlim(float(y[0]), float(y[-1]))
    ymin = min(0.0, float(np.min(normalized_wem)))
    ymax = max(1.0, float(np.max(normalized_wem)))
    padding = max(0.02, 0.04 * (ymax - ymin))
    ax.set_ylim(ymin - padding, ymax + padding)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return y, normalized_wem


def mode_validity(
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


def export_simulation_results(
    sim_run,
    export_dir: Path,
    validity_boundary: np.ndarray,
    *,
    label: str,
):
    import pandas as pd

    from comsol_workflow.energy_recovery import interpolate_field

    export_dir.mkdir(parents=True, exist_ok=True)
    eigenfrequencies = sim_run.get_eigenfrequencies()
    df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])
    validity = mode_validity(sim_run, len(df), validity_boundary)
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
        coords_wem, wem = sim_run.get_2d_fields(
            mode_idx,
            "ewfd.Wav",
            "center",
        )
        save_wem_cutline_plot(
            export_dir / f"{mode_idx:02d}_cutline_Wem.png",
            coords_wem,
            wem,
            dpi=DPI,
        )
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Wem_2d.png",
            coords_wem,
            np.real(wem),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Wem",
            unit_label="J/m³",
            color_scale_mode="linear",
            dpi=DPI,
        )

    return df


def simulation_layers(records: list[dict[str, object]]):
    from comsol_workflow.simulation_utils import LayerSpec

    return [
        LayerSpec(
            height=SLAB_HEIGHT,
            holes=polygons_from_records(records),
            refractive_index=REFRACTIVE_INDEX,
            label="patterned_slab",
        )
    ]


def simulation_config():
    from comsol_workflow.simulation_utils import SimulationConfig

    return SimulationConfig(
        wavelength=WAVELENGTH,
        eigenfrequency_shift=EIGENFREQUENCY_SHIFT,
        eigenmode_count=EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=MODE_TYPE,
    )


def shared_case_metadata(metadata: dict[str, object], validity_boundary: np.ndarray) -> dict[str, object]:
    validity_layers = strip_validity_layer_count()
    return {
        **metadata,
        **shift_profile_metadata(
            CLADDING_LAYERS,
            CLADDING_INWARD_SHIFT_FACTOR,
            A,
            CLADDING_SHIFT_PROFILE,
        ),
        "strip_x_min": STRIP_X_MIN,
        "strip_x_max": STRIP_X_MAX,
        "validity_region": {
            "type": "bulk plus floor(half cladding layers) boundary",
            "half_cladding_layers": validity_layers,
            "shell_radius": BULK_RADIUS + validity_layers,
            "boundary": validity_boundary,
        },
        "output_layout_version": OUTPUT_LAYOUT_VERSION,
        "strip_cut_directory": strip_cut_dirname(),
        "mode_validity": {
            "metric": "gamma_subspace_p",
            "operator": ">",
            "threshold": GAMMA_SUBSPACE_P_VALID_THRESHOLD,
        },
        "export_dir": (
            f"{RESULTS_DIRNAME}/mode{{mode_idx}}/{MODE_EXPORT_DIRNAME}"
        ),
        "shared_parameters": ACTIVE_PARAMETERS.to_metadata(),
        "simulation_common": {
            "slab_height": SLAB_HEIGHT,
            "refractive_index": REFRACTIVE_INDEX,
            "wavelength": WAVELENGTH,
            "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
            "eigenmode_count": EIGENMODE_COUNT,
            "mesh_auto_size": MESH_AUTO_SIZE,
            "mode_type": MODE_TYPE,
        },
    }


def run_strip_cut_case(
    out_dir: Path,
    records: list[dict[str, object]],
    metadata: dict[str, object],
    validity_boundary: np.ndarray,
) -> None:
    from comsol_workflow.simulation_utils import BoundarySpec, SimulationRun

    case_dir = strip_case_dir(out_dir)
    paths = prepare_strip_case_output(case_dir)
    export_dir = paths["export_dir"]

    angle = validate_strip_cut_angle(STRIP_CUT_ANGLE)
    tangent, normal = strip_cut_basis(angle)
    local_metadata = metadata_to_strip_frame(metadata, angle)
    local_validity_boundary = points_to_strip_frame(validity_boundary, angle)
    clipped = strip_records_for_cut(records, angle)
    clipped_cells = strip_cell_records(local_metadata)
    boundary, strip_lx, strip_ly = simulation_boundary_for_strip(clipped_cells)
    case_metadata = {
        **shared_case_metadata(local_metadata, local_validity_boundary),
        "case": "strip_1d",
        "source": "full finite lattice rigidly transformed to the cut-local frame, then clipped",
        "strip_cut_angle_degrees": angle,
        "solver_frame_rotation_degrees": -angle,
        "global_periodic_tangent": tangent,
        "global_finite_normal": normal,
        "global_inner_boundary": metadata["inner_boundary"],
        "strip_holes": len(clipped),
        "simulation_holes": len(clipped),
        "strip_cells": len(clipped_cells),
        "simulation_cells": len(clipped_cells),
        "strip_Lx": strip_lx,
        "strip_Ly": strip_ly,
        "strip_boundary": boundary,
        "simulation": {
            "mode": "strip_1d",
            "coordinate_frame": "cut-local",
            "finite_axis": "local_y_global_normal",
            "periodic_axis": "local_x_global_tangent",
            "kx_only": True,
            "kx": STRIP_KX,
            "kx_direction_global": tangent,
        },
    }

    if RUN_DIAGNOSTIC_PLOTS:
        save_shift_profile_diagnostics(
            paths["overview_dir"],
            case_metadata,
            dpi=DPI,
        )
        save_global_cut_preview(
            paths["model_dir"] / "strip_cut_global_preview.png",
            clipped,
            boundary,
            angle,
        )
        save_strip_plot(
            paths["model_dir"] / "strip_cut_local_holes.png",
            clipped,
            clipped_cells,
            case_metadata,
            annotated=False,
        )
        save_strip_plot(
            paths["model_dir"] / "strip_cut_local_holes_annotated.png",
            clipped,
            clipped_cells,
            case_metadata,
            annotated=True,
        )
        save_simulation_plot(
            paths["model_dir"] / "strip_cut_local_simulation.png",
            clipped,
            boundary,
            case_metadata,
            annotated=False,
        )
        save_simulation_plot(
            paths["model_dir"] / "strip_cut_local_simulation_annotated.png",
            clipped,
            boundary,
            case_metadata,
            annotated=True,
            cell_records_for_strip=clipped_cells,
        )
    write_config(paths["config"], case_metadata)

    print(
        f"Prepared strip-cut geometry in {case_dir}; holes={len(clipped)}, "
        f"cells={len(clipped_cells)}, Lx={strip_lx:.6g}, Ly={strip_ly:.6g}"
    )

    progress_log_path = paths["progress_log"]
    print(f"COMSOL progress log: {progress_log_path}")

    with SimulationRun(config=simulation_config(), progress_log_path=progress_log_path) as sim:
        mph_path = paths["mph"]
        boundary = BoundarySpec("rectangle", a=strip_lx, b=strip_ly)
        layers = simulation_layers(clipped)
        sim.build_geometry(boundary, layers, simulation_mode="strip_1d", k={"kx": STRIP_KX})
        sim.model.save(str(mph_path))
        sim.run_simulation(mesh_save_path=mph_path)
        sim.model.save(str(mph_path))
        export_simulation_results(sim, export_dir, local_validity_boundary, label="strip_1d")
    run_strip_bulk_fourier_postprocess(
        paths["staging_dir"],
        export_dir,
        bulk_records=strip_bulk_cell_records(angle),
        bulk_radius=BULK_RADIUS,
        period=A,
        kx_fraction=STRIP_KX,
        rho_grid_size=STRIP_FOURIER_RHO_GRID_SIZE,
        top_k=STRIP_FOURIER_TOP_K,
        mode_indices=STRIP_FOURIER_MODE_INDICES,
        dpi=DPI,
        generate_plots=RUN_DIAGNOSTIC_PLOTS,
        row_index_coefficients=strip_fourier_row_index_coefficients(angle),
    )
    score_strip_modes(
        paths["staging_dir"],
        export_dirname=EXPORT_DIRNAME,
        gamma_subspace_p_valid_threshold=(
            GAMMA_SUBSPACE_P_VALID_THRESHOLD
        ),
    )
    finalize_strip_case_output(case_dir)


def run_finite_full_case(
    out_dir: Path,
    records: list[dict[str, object]],
    metadata: dict[str, object],
    validity_boundary: np.ndarray,
) -> None:
    from comsol_workflow.simulation_utils import BoundarySpec, SimulationRun

    case_dir = finite_full_case_dir(out_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    export_dir = case_export_dir(case_dir)

    cell_records_for_full = full_cell_records(metadata)
    boundary, finite_hex_a, finite_lx, finite_ly = simulation_boundary_for_full_lattice_hexagon(cell_records_for_full)
    case_metadata = {
        **shared_case_metadata(metadata, validity_boundary),
        "case": "finite_cavity",
        "source": "full finite lattice before strip clipping",
        "simulation_holes": len(records),
        "simulation_cells": len(cell_records_for_full),
        "finite_hex_a": finite_hex_a,
        "finite_Lx": finite_lx,
        "finite_Ly": finite_ly,
        "finite_boundary": boundary,
        "simulation": {
            "mode": "finite_size",
            "boundary": "hexagon_boundary",
            "periodic_axis": None,
            "kx_only": False,
        },
    }

    save_shift_profile_diagnostics(case_dir, case_metadata, dpi=DPI)
    save_full_lattice_plot(case_dir / "full_lattice_modulated_holes.png", records, case_metadata, annotated=False)
    save_full_lattice_plot(case_dir / "full_lattice_modulated_holes_annotated.png", records, case_metadata, annotated=True)
    save_full_lattice_shift_plot(case_dir / "full_lattice_modulation_vectors.png", case_metadata)
    save_full_lattice_normal_fan_plot(case_dir / "full_lattice_normal_fan_regions.png", case_metadata)
    save_finite_simulation_plot(
        case_dir / "full_finite_simulation.png",
        records,
        cell_records_for_full,
        boundary,
        validity_boundary,
        case_metadata,
        annotated=False,
    )
    save_finite_simulation_plot(
        case_dir / "full_finite_simulation_annotated.png",
        records,
        cell_records_for_full,
        boundary,
        validity_boundary,
        case_metadata,
        annotated=True,
    )
    write_config(case_dir / "config.json", case_metadata)

    print(
        f"Saved finite-full geometry to {case_dir}; holes={len(records)}, "
        f"cells={len(cell_records_for_full)}, Lx={finite_lx:.6g}, Ly={finite_ly:.6g}"
    )

    progress_log_path = case_dir / "comsol_progress.log"
    print(f"COMSOL progress log: {progress_log_path}")

    with SimulationRun(config=simulation_config(), progress_log_path=progress_log_path) as sim:
        boundary = BoundarySpec("hexagon_boundary", a=finite_hex_a)
        layers = simulation_layers(records)
        mph_path = case_dir / "finite_cavity.mph"
        sim.build_geometry(boundary, layers, simulation_mode="finite_size", k=None)
        sim.model.save(str(mph_path))
        sim.run_simulation(mesh_save_path=mph_path)
        sim.model.save(str(mph_path))
        export_simulation_results(sim, export_dir, validity_boundary, label="finite_cavity")
    run_finite_lattice_fourier_postprocess(
        case_dir,
        export_dir,
        cells=bulk_points(),
        cladding_cells=list(metadata["cladding_points"]),
        bulk_radius=BULK_RADIUS,
        period=A,
        rho_grid_size=FINITE_FOURIER_RHO_GRID_SIZE,
        top_k=FINITE_FOURIER_TOP_K,
        mode_indices=FINITE_FOURIER_MODE_INDICES,
        dpi=DPI,
    )
    score_finite_modes(case_dir, export_dirname=EXPORT_DIRNAME)


def set_cladding_inward_shift_factor(shift_factor: float) -> Path:
    global CLADDING_INWARD_SHIFT
    global CLADDING_INWARD_SHIFT_FACTOR

    CLADDING_INWARD_SHIFT_FACTOR = float(shift_factor)
    CLADDING_INWARD_SHIFT = CLADDING_INWARD_SHIFT_FACTOR * A
    return output_dir_for_cladding_shift_factor(shift_factor)


def strip_case_is_complete(out_dir: Path, shift_factor: float) -> bool:
    if not RUN_STRIP_CUT or RUN_FINITE_FULL:
        return False
    paths = strip_case_output_paths(out_dir)
    config_path = paths["config"]
    layout_version = OUTPUT_LAYOUT_VERSION
    if not config_path.is_file():
        layout_version = 1
        config_path = out_dir / "config.json"
    if layout_version == OUTPUT_LAYOUT_VERSION:
        required_paths = [
            config_path,
            paths["mph"],
            paths["overview_dir"] / "eigenfrequencies.csv",
            paths["overview_dir"] / "strip_bulk_fourier_summary.csv",
            paths["overview_dir"] / "mode_scores.csv",
            paths["config_dir"] / "objective.json",
        ]
    else:
        required_paths = [
            config_path,
            out_dir / "strip_1d.mph",
            out_dir / EXPORT_DIRNAME / "eigenfrequencies.csv",
            out_dir / "strip_bulk_fourier_hz"
            / "strip_bulk_fourier_summary.csv",
            out_dir / "mode_scores.csv",
            out_dir / "objective.json",
        ]
    if not all(path.is_file() for path in required_paths):
        return False
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        simulation_common = config["simulation_common"]
        expected_profile = shift_profile_metadata(
            CLADDING_LAYERS,
            shift_factor,
            A,
            CLADDING_SHIFT_PROFILE,
        )
        metadata_matches = (
            int(config["bulk_radius"]) == BULK_RADIUS
            and int(config["cladding_layers"]) == CLADDING_LAYERS
            and np.isclose(float(config["strip_cut_angle_degrees"]), STRIP_CUT_ANGLE)
            and np.isclose(float(config["cladding_inward_shift"]), shift_factor * A)
            and config["cladding_shift_profile"]
            == expected_profile["cladding_shift_profile"]
            and config["resolved_cladding_shift_by_layer"]
            == expected_profile["resolved_cladding_shift_by_layer"]
            and float(simulation_common["mesh_auto_size"]) == float(MESH_AUTO_SIZE)
            and simulation_common["eigenfrequency_shift"] == EIGENFREQUENCY_SHIFT
        )
        if layout_version == OUTPUT_LAYOUT_VERSION:
            metadata_matches = metadata_matches and (
                int(config.get("output_layout_version", 0))
                == OUTPUT_LAYOUT_VERSION
                and config.get("strip_cut_directory") == strip_cut_dirname()
                and config.get("mode_validity")
                == {
                    "metric": "gamma_subspace_p",
                    "operator": ">",
                    "threshold": GAMMA_SUBSPACE_P_VALID_THRESHOLD,
                }
            )
        if not metadata_matches or layout_version == 1:
            return metadata_matches

        import pandas as pd

        eigen_df = pd.read_csv(
            paths["overview_dir"] / "eigenfrequencies.csv"
        )
        score_df = pd.read_csv(paths["overview_dir"] / "mode_scores.csv")
        fourier_df = pd.read_csv(
            paths["overview_dir"] / "strip_bulk_fourier_summary.csv"
        )
        mode_indices = set(eigen_df["mode_idx"].astype(int).tolist())
        score_modes = set(score_df["mode_idx"].astype(int).tolist())
        fourier_modes = set(fourier_df["mode_idx"].astype(int).tolist())
        if score_modes != fourier_modes or score_modes != mode_indices:
            return False
        expected_valid = score_df["gamma_subspace_p"].gt(
            GAMMA_SUBSPACE_P_VALID_THRESHOLD
        )
        if "is_valid" not in score_df.columns:
            return False
        actual_valid = score_df["is_valid"].map(
            lambda value: str(value).strip().lower() in {"true", "1"}
        )
        if not actual_valid.equals(expected_valid):
            return False
        valid_modes = set(
            score_df.loc[expected_valid, "mode_idx"].astype(int).tolist()
        )
        existing_mode_dirs = {
            int(path.name.removeprefix("mode"))
            for path in paths["results_dir"].glob("mode*")
            if path.is_dir() and path.name.removeprefix("mode").isdigit()
        }
        if existing_mode_dirs != valid_modes:
            return False
        for mode_idx in sorted(valid_modes):
            mode_dir = _strip_mode_dir(paths, mode_idx)
            mode_required = [
                mode_dir / MODE_OVERVIEW_DIRNAME / "eigenfrequency.csv",
                mode_dir / MODE_EXPORT_DIRNAME / "Hz_center.parquet",
                mode_dir / MODE_EXPORT_DIRNAME / "Hz_Re_2d.png",
                mode_dir / MODE_EXPORT_DIRNAME / "Hz_Im_2d.png",
                mode_dir / MODE_EXPORT_DIRNAME / "Wem_2d.png",
                mode_dir / MODE_EXPORT_DIRNAME / "cutline_Wem.png",
                mode_dir / MODE_OVERVIEW_DIRNAME / "mode_score.csv",
            ]
            if not all(path.is_file() for path in mode_required):
                return False
            fourier_mode_dir = mode_dir / MODE_FOURIER_DIRNAME
            if not fourier_mode_dir.is_dir() or not any(
                fourier_mode_dir.iterdir()
            ):
                return False
        return True
    except (
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        OSError,
    ):
        return False


def strip_case_summary(shift_factor: float, out_dir: Path) -> dict[str, object]:
    return {
        "output_layout_version": OUTPUT_LAYOUT_VERSION,
        "strip_cut_directory": strip_cut_dirname(),
        "mode_validity_metric": "gamma_subspace_p",
        "mode_validity_operator": ">",
        "mode_validity_threshold": GAMMA_SUBSPACE_P_VALID_THRESHOLD,
        "cladding_inward_shift_factor": shift_factor,
        "cladding_inward_shift": shift_factor * A,
        **shift_profile_metadata(
            CLADDING_LAYERS,
            shift_factor,
            A,
            CLADDING_SHIFT_PROFILE,
        ),
        "strip_cut_angle_degrees": STRIP_CUT_ANGLE,
        "out_dir": str(out_dir),
        "bulk_radius": BULK_RADIUS,
        "cladding_layers": CLADDING_LAYERS,
        "mesh_auto_size": MESH_AUTO_SIZE,
        "eigenmode_count": EIGENMODE_COUNT,
        "case": "strip_1d",
    }


def run_case_for_cladding_shift_factor(shift_factor: float) -> dict[str, object]:
    out_dir = set_cladding_inward_shift_factor(shift_factor)
    if REUSE_COMPLETED_CASES and strip_case_is_complete(out_dir, shift_factor):
        print(f"Reusing completed strip case: {out_dir}")
        return strip_case_summary(shift_factor, out_dir)
    print(
        f"Running target CLADDING_INWARD_SHIFT={shift_factor:.3f}*A "
        f"({CLADDING_INWARD_SHIFT:.6g}); "
        f"profile={CLADDING_SHIFT_PROFILE.label} -> {out_dir}"
    )

    records, metadata = build_full_lattice_holes()
    for message in shift_profile_preflight_messages(metadata):
        print(message)
    validity_boundary = strip_validity_boundary()

    if RUN_STRIP_CUT:
        run_strip_cut_case(out_dir, records, metadata, validity_boundary)

    return strip_case_summary(shift_factor, out_dir)


def completed_strip_case_summaries(cut_dir: Path) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for case_dir in sorted(Path(cut_dir).glob("shift*")):
        if not case_dir.is_dir():
            continue
        match = re.fullmatch(r"shift(-?\d+(?:\.\d+)?)", case_dir.name)
        if match is None:
            continue
        shift_factor = float(match.group(1))
        if strip_case_is_complete(case_dir, shift_factor):
            summaries.append(strip_case_summary(shift_factor, case_dir))
    return sorted(
        summaries,
        key=lambda row: float(row["cladding_inward_shift_factor"]),
    )


def merge_strip_case_summaries(
    *summary_groups: list[dict[str, object]],
) -> list[dict[str, object]]:
    by_shift: dict[float, dict[str, object]] = {}
    for summaries in summary_groups:
        for summary in summaries:
            shift_factor = float(summary["cladding_inward_shift_factor"])
            by_shift[round(shift_factor, 12)] = summary
    return [by_shift[key] for key in sorted(by_shift)]


def should_run_gamma_trend_after_scan(*, shift_count: int | None = None) -> bool:
    if shift_count is None:
        shift_count = len(CLADDING_INWARD_SHIFT_FACTORS)
    trend_enabled = (
        shift_count > 3
        if RUN_GAMMA_TREND_AFTER_SCAN is None
        else RUN_GAMMA_TREND_AFTER_SCAN
    )
    return (
        RUN_STRIP_CUT
        and trend_enabled
        and shift_count > 1
    )


def run_gamma_trend_after_scan(summary_path: Path) -> tuple[Path, list[str], Path]:
    from scripts.analysis.srip1d_trend import (
        TREND_OUTPUT_DIR_NAME,
        generate_trend_from_run_summary,
    )

    target_frequency = float(EIGENFREQUENCY_SHIFT.split()[0])
    return generate_trend_from_run_summary(
        run_summary_path=summary_path,
        output_dir=summary_path.parent / TREND_OUTPUT_DIR_NAME,
        target_frequency=target_frequency,
    )


def run_case() -> list[dict[str, object]]:
    if RUN_FINITE_FULL:
        raise ValueError(
            "RUN_FINITE_FULL is no longer supported by run_strip_1d.py; "
            "use scripts/run_main/run_finite.py so finite results remain "
            "under scripts/.out/finite_cavity/."
        )
    cut_dir = strip_cut_output_dir()
    cut_dir.mkdir(parents=True, exist_ok=True)
    summaries = [
        run_case_for_cladding_shift_factor(shift_factor)
        for shift_factor in CLADDING_INWARD_SHIFT_FACTORS
    ]
    summaries = merge_strip_case_summaries(
        completed_strip_case_summaries(cut_dir),
        summaries,
    )
    summary_path = cut_dir / "run_summary.json"
    write_config(
        summary_path,
        {
            "output_layout_version": OUTPUT_LAYOUT_VERSION,
            "runs": summaries,
        },
    )
    if should_run_gamma_trend_after_scan(shift_count=len(summaries)):
        run_gamma_trend_after_scan(summary_path)
    return summaries


def main() -> None:
    run_case()


if __name__ == "__main__":
    main()
