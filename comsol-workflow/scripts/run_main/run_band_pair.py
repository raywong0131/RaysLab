#!/usr/bin/env python3
"""Calculate and compare cavity/cladding bands with one Model per geometry."""

from __future__ import annotations

import math
import json
import os
import sys
from dataclasses import replace
from pathlib import Path
from typing import Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.band_connector import (  # noqa: E402
    assign_mode_candidates,
    field_overlap,
    load_hz_center,
    load_mode_composition,
)
from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis  # noqa: E402
from comsol_workflow.energy_recovery import interpolate_field  # noqa: E402
from comsol_workflow.geometry_utils import create_hexagon_design, visualize_hexagon_design  # noqa: E402
from comsol_workflow.simulation_spatial.hexagon_unit_cell import (  # noqa: E402
    SimulationConfig,
    SimulationRun as BaseSimulationRun,
)
from scripts.run_main.parameter_config import (  # noqa: E402
    ACTIVE_PARAMETERS,
    UNIT_CELL_OUTPUT_ROOT,
    update_center_frequency_from_cavity_p2,
)

OUT_DIR = UNIT_CELL_OUTPUT_ROOT / (
    f"unitcell_band_{ACTIVE_PARAMETERS.case_parameter_label}"
)
RESULTS_DIRNAME = "01_results"
OVERVIEW_DIRNAME = "10_overview"
CONFIG_DIRNAME = "99_config"
LEGACY_CASE_OVERVIEW_FILES = (
    "geometry.png",
    "connected_bands_gamma_m.csv",
    "connected_bands_gamma_k.csv",
    "selected_bands.csv",
    "mode_composition.csv",
    "q_power_law_fits.csv",
)

A = 0.82
R_0 = A / 3.0
B_0 = 0.23
D = 0.02

SLAB_HEIGHT = "200 [nm]"
REFRACTIVE_INDEX = "3.3"
EIGENMODE_COUNT = ACTIVE_PARAMETERS.unit_cell_eigenmode_count
MESH_AUTO_SIZE = ACTIVE_PARAMETERS.mesh_size
MODE_TYPE = "TE"
WAVELENGTH = "1550 [nm]"
EIGENFREQUENCY_SHIFT = "c_const/1.55[um]"
UPDATE_CENTER_FREQUENCY_AFTER_RUN = True

FREQUENCY_MIN_THZ = 175.0
FREQUENCY_MAX_THZ = 220.0
BAND_POINTS_PER_ARM = 31
Q_FIT_K_MAGNITUDES = (0.002, 0.004, 0.006, 0.010, 0.015, 0.020, 0.030, 0.050)
K_STRING_PRECISION = 6
BAND_FREQUENCY_TOLERANCE_THZ = 7.0
BAND_MINIMUM_OVERLAP = 0.2

BULK_PARAMS = ACTIVE_PARAMETERS.cavity.to_fourier_params("cavity_p_bic")
CLADDING_PARAMS = ACTIVE_PARAMETERS.cladding.to_fourier_params("bulk_gap_centered")

REQUIRED_POINT_COLUMNS = {"re", "is_valid", "p_weight"}
POINT_SOLVER_CONFIG_FILENAME = "solver_config.json"
# Stored px/py weights use the project-internal convention.  User-facing labels
# reverse that mapping: internal px is user py (red), internal py is user px (blue).
P_X_COLOR = "#dc2626"
P_Y_COLOR = "#2563eb"

_MESH_AUTO_SIZE_UNSET = object()


class ReusableSimulationRun(BaseSimulationRun):
    """Own one Model and reuse it for every k point of one fixed geometry."""

    def __init__(self, config=None):
        super().__init__(config)
        self._geometry_signature = None
        self._applied_mesh_auto_size = None
        self._current_k = None
        self.model_build_count = 0
        self.mesh_build_count = 0
        self.solve_count = 0

    @staticmethod
    def geometry_signature(a, holes):
        return (
            float(a),
            tuple(
                tuple(tuple(float(value) for value in vertex) for vertex in hole)
                for hole in holes
            ),
        )

    @staticmethod
    def _normalize_k(k):
        if k is None:
            return {"kx": 0.0, "ky": 0.0}
        return {"kx": float(k["kx"]), "ky": float(k["ky"])}

    @staticmethod
    def _validate_mesh_auto_size(mesh_auto_size):
        if mesh_auto_size is None:
            return None
        if type(mesh_auto_size) is not int or not 1 <= mesh_auto_size <= 9:
            raise ValueError("mesh_auto_size must be None or an integer from 1 to 9")
        return mesh_auto_size

    def _set_k(self, k):
        normalized = self._normalize_k(k)
        jmodel = self.model.java
        jmodel.param().set("kx", f"{normalized['kx']}*G")
        jmodel.param().set("ky", f"{normalized['ky']}*G")
        self._current_k = normalized

    def _rebuild_mesh(self, mesh_auto_size):
        mesh = self.model.java.component("comp1").mesh("mesh1")
        if mesh_auto_size is not None:
            mesh.autoMeshSize(int(mesh_auto_size))
        mesh.run()
        self._applied_mesh_auto_size = mesh_auto_size
        self.mesh_build_count += 1

    def build_and_run(
        self,
        a,
        holes,
        k=None,
        mesh_auto_size=_MESH_AUTO_SIZE_UNSET,
    ):
        signature = self.geometry_signature(a, holes)
        desired_mesh = self._validate_mesh_auto_size(
            self.config.mesh_auto_size
            if mesh_auto_size is _MESH_AUTO_SIZE_UNSET
            else mesh_auto_size
        )
        if self._geometry_signature is None:
            if desired_mesh != self.config.mesh_auto_size:
                self.config = replace(self.config, mesh_auto_size=desired_mesh)
            super().build_and_run(a, holes, k)
            self._geometry_signature = signature
            self._applied_mesh_auto_size = desired_mesh
            self._current_k = self._normalize_k(k)
            self.model_build_count += 1
            self.mesh_build_count += 1
            self.solve_count += 1
            return

        if signature != self._geometry_signature:
            raise ValueError(
                "ReusableSimulationRun already owns a different geometry; "
                "create another runner for the changed geometry"
            )

        self._set_k(k)
        if desired_mesh != self._applied_mesh_auto_size:
            if desired_mesh is None:
                raise ValueError(
                    "Cannot restore the implicit default after an explicit COMSOL "
                    "auto mesh size; create another runner instead"
                )
            self._rebuild_mesh(desired_mesh)

        solution = self.model.java.sol("sol1")
        solution.clearSolutionData()
        solution.runAll()
        self.solve_count += 1

    @property
    def execution_stats(self):
        return {
            "execution_engine": "reused_model",
            "model_reuse_enabled": True,
            "model_build_count": int(self.model_build_count),
            "mesh_build_count": int(self.mesh_build_count),
            "solve_count": int(self.solve_count),
            "mesh_auto_size": self._applied_mesh_auto_size,
            "current_k": None if self._current_k is None else dict(self._current_k),
        }


# Preserve the main module's public runner name for callers that construct or
# replace the per-point runner directly. The canonical runner now reuses Models.
SimulationRun = ReusableSimulationRun


def unit_cell_case_paths(root: Path, case_name: str) -> dict[str, Path]:
    """Return the fixed result, overview, and config locations for one cell."""
    root = Path(root)
    return {
        "result": root / RESULTS_DIRNAME / str(case_name),
        "overview": root / OVERVIEW_DIRNAME / str(case_name),
        "config": root / CONFIG_DIRNAME / str(case_name),
    }


def _tree_file_stats(root: Path) -> tuple[int, int]:
    files = [path for path in Path(root).rglob("*") if path.is_file()]
    return len(files), sum(path.stat().st_size for path in files)


def migrate_legacy_unit_cell_layout(root: Path) -> dict[str, int]:
    """Move one complete legacy result into the durable unit-cell layout."""
    root = Path(root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Unit-cell result does not exist: {root}")
    formal_names = {RESULTS_DIRNAME, OVERVIEW_DIRNAME, CONFIG_DIRNAME}
    existing_formal = sorted(name for name in formal_names if (root / name).exists())
    if existing_formal:
        raise FileExistsError(
            f"Formal unit-cell directories already exist: {existing_formal}"
        )

    required_root_files = {
        "config.json",
        "unit_cell_band_comparison.png",
        "p_bands_q_vs_k_comparison.png",
    }
    allowed_root_names = required_root_files | {"cavity", "cladding", "80_logs"}
    root_names = {path.name for path in root.iterdir()}
    unknown_root = sorted(root_names.difference(allowed_root_names))
    missing_root = sorted(required_root_files.difference(root_names))
    if unknown_root:
        raise RuntimeError(f"Legacy unit-cell root has unknown entries: {unknown_root}")
    if missing_root:
        raise RuntimeError(f"Legacy unit-cell root is incomplete: {missing_root}")

    expected_case_names = {"k_points", "config.json", *LEGACY_CASE_OVERVIEW_FILES}
    for case_name in ("cavity", "cladding"):
        legacy_case = root / case_name
        if not legacy_case.is_dir():
            raise RuntimeError(f"Legacy unit-cell case is missing: {legacy_case}")
        case_names = {path.name for path in legacy_case.iterdir()}
        unknown = sorted(case_names.difference(expected_case_names))
        missing = sorted(expected_case_names.difference(case_names))
        if unknown:
            raise RuntimeError(f"{case_name} has unknown legacy entries: {unknown}")
        if missing:
            raise RuntimeError(f"{case_name} is incomplete: {missing}")

    files_before, bytes_before = _tree_file_stats(root)
    for case_name in ("cavity", "cladding"):
        legacy_case = root / case_name
        paths = unit_cell_case_paths(root, case_name)
        paths["result"].mkdir(parents=True)
        paths["overview"].mkdir(parents=True)
        paths["config"].mkdir(parents=True)
        (legacy_case / "k_points").replace(paths["result"] / "k_points")
        for filename in LEGACY_CASE_OVERVIEW_FILES:
            (legacy_case / filename).replace(paths["overview"] / filename)
        (legacy_case / "config.json").replace(paths["config"] / "config.json")
        legacy_case.rmdir()

    overview_root = root / OVERVIEW_DIRNAME
    config_root = root / CONFIG_DIRNAME
    for filename in (
        "unit_cell_band_comparison.png",
        "p_bands_q_vs_k_comparison.png",
    ):
        (root / filename).replace(overview_root / filename)
    (root / "config.json").replace(config_root / "config.json")

    files_after, bytes_after = _tree_file_stats(root)
    if (files_after, bytes_after) != (files_before, bytes_before):
        raise RuntimeError(
            "Unit-cell layout migration changed file totals: "
            f"before=({files_before}, {bytes_before}), "
            f"after=({files_after}, {bytes_after})"
        )
    return {
        "files_before": files_before,
        "files_after": files_after,
        "bytes_before": bytes_before,
        "bytes_after": bytes_after,
    }


def make_branch_k_points(
    direction: str,
    point_count: int,
    dense_magnitudes: Sequence[float] = (),
) -> list[dict[str, float | str]]:
    """Return an outward Gamma-to-edge branch in dimensionless units of G."""
    endpoints = {"gamma_m": 0.5, "gamma_k": 1.0 / math.sqrt(3.0)}
    if direction not in endpoints:
        raise ValueError(f"Unknown direction: {direction}")
    if point_count < 2:
        raise ValueError("point_count must be at least 2")

    endpoint = endpoints[direction]
    dense = [float(value) for value in dense_magnitudes]
    if any(value < 0.0 or value > endpoint for value in dense):
        raise ValueError(f"Dense k magnitude is outside {direction} endpoint {endpoint:g}")

    magnitudes = np.concatenate([np.linspace(0.0, endpoint, point_count), dense])
    rounded = sorted({round(float(value), K_STRING_PRECISION) for value in magnitudes})
    points = []
    for magnitude in rounded:
        kx = magnitude if direction == "gamma_k" else 0.0
        ky = magnitude if direction == "gamma_m" else 0.0
        points.append(
            {
                "direction": direction,
                "kx": kx,
                "ky": ky,
                "kx_str": f"{kx:.{K_STRING_PRECISION}f}",
                "ky_str": f"{ky:.{K_STRING_PRECISION}f}",
                "k_norm": magnitude,
                "path_fraction": magnitude / endpoint,
            }
        )
    return points


def case_k_points(case_name: str) -> dict[str, list[dict[str, float | str]]]:
    """Return identical directional and Q-fit k points for both cells."""
    if case_name not in {"cavity", "cladding"}:
        raise ValueError(f"Unknown unit-cell case: {case_name}")
    return {
        "gamma_m": make_branch_k_points(
            "gamma_m", BAND_POINTS_PER_ARM, Q_FIT_K_MAGNITUDES
        ),
        "gamma_k": make_branch_k_points(
            "gamma_k", BAND_POINTS_PER_ARM, Q_FIT_K_MAGNITUDES
        ),
    }


def k_point_dir(case_dir: Path, point: dict[str, object]) -> Path:
    return case_dir / "k_points" / f"kx={point['kx_str']}_ky={point['ky_str']}"


def point_solver_identity(config: SimulationConfig) -> dict[str, object]:
    """Return the solver settings that determine k-point cache compatibility."""
    return {
        "eigenmode_count": int(config.eigenmode_count),
        "eigenfrequency_shift": str(config.eigenfrequency_shift),
        "mesh_auto_size": config.mesh_auto_size,
        "mode_type": str(config.mode_type),
        "wavelength": str(config.wavelength),
    }


def k_point_is_complete(
    k_dir: Path,
    simulation_config: SimulationConfig | None = None,
) -> bool:
    """Return whether a saved k point has its table and every valid Hz field."""
    if simulation_config is not None:
        identity_path = Path(k_dir) / POINT_SOLVER_CONFIG_FILENAME
        try:
            identity = json.loads(identity_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            return False
        if identity != point_solver_identity(simulation_config):
            return False
    csv_path = Path(k_dir) / "eigenfrequencies.csv"
    try:
        frame = pd.read_csv(csv_path)
    except (OSError, ValueError, pd.errors.ParserError):
        return False
    if not REQUIRED_POINT_COLUMNS.issubset(frame.columns):
        return False

    field_columns = {"x", "y", "re", "im"}
    for mode_idx in frame.index[frame["is_valid"].astype(bool)]:
        field_path = Path(k_dir) / "eigenmodes" / f"{int(mode_idx):02d}_Hz_center.parquet"
        try:
            field = pd.read_parquet(field_path)
        except (OSError, ValueError, ImportError):
            return False
        if not field_columns.issubset(field.columns):
            return False
    return True


def get_hole_params(params: dict[str, float]) -> np.ndarray:
    """Convert the named sixfold Fourier parameters into six hole rows."""
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

    r_f *= R_0 * float(params["r_f0"])
    b_square_f *= B_0 * B_0 * float(params["b_square_f0"])
    minimum_values = get_fourier_basis_min_nonzero_values(6)
    r = fourier_to_standard_basis(r_f / minimum_values, 6)
    theta = fourier_to_standard_basis(theta_f / minimum_values, 6)
    b_square = fourier_to_standard_basis(b_square_f / minimum_values, 6)
    b = np.sqrt(np.maximum(b_square, 0.0))
    phi = fourier_to_standard_basis(phi_f / minimum_values, 6)
    return np.stack([r, theta, b, phi], axis=1)


def write_csv_atomic(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(item) for item in value]
    return value


def write_json_atomic(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(_to_jsonable(data), handle, indent=2)
    os.replace(temporary, path)


def _field_peak_z(sim_run, mode_idx: int, expression: str, plane: str) -> float:
    coordinates, values = sim_run.get_2d_fields(mode_idx, expression, plane)
    return float(coordinates[np.argmax(np.real(values)), 1])


def _mode_is_valid(sim_run, mode_idx: int) -> tuple[bool, dict[str, float]]:
    peaks = {
        "z_peak_normH_xz": _field_peak_z(sim_run, mode_idx, "ewfd.normH", "xz"),
        "z_peak_normE_xz": _field_peak_z(sim_run, mode_idx, "ewfd.normE", "xz"),
        "z_peak_normH_yz": _field_peak_z(sim_run, mode_idx, "ewfd.normH", "yz"),
        "z_peak_normE_yz": _field_peak_z(sim_run, mode_idx, "ewfd.normE", "yz"),
    }
    return all(value <= 1.0 for value in peaks.values()), peaks


def run_k_point(
    case_dir: Path,
    holes: list[list[list[float]]],
    point: dict[str, object],
    simulation_config: SimulationConfig,
    sim_run: ReusableSimulationRun | None = None,
) -> pd.DataFrame:
    """Run or load one Floquet point and persist all valid center Hz fields."""
    output_dir = k_point_dir(case_dir, point)
    csv_path = output_dir / "eigenfrequencies.csv"
    if k_point_is_complete(output_dir, simulation_config):
        return pd.read_csv(csv_path)

    k = {"kx": float(point["kx"]), "ky": float(point["ky"])}
    if sim_run is not None:
        sim_run.build_and_run(A, holes, k, mesh_auto_size=MESH_AUTO_SIZE)
        frame = persist_k_point_solution(sim_run, output_dir)
        write_json_atomic(
            output_dir / POINT_SOLVER_CONFIG_FILENAME,
            point_solver_identity(simulation_config),
        )
        return frame

    with SimulationRun(simulation_config) as owned_runner:
        owned_runner.build_and_run(A, holes, k)
        frame = persist_k_point_solution(owned_runner, output_dir)
        write_json_atomic(
            output_dir / POINT_SOLVER_CONFIG_FILENAME,
            point_solver_identity(simulation_config),
        )
        return frame


def persist_k_point_solution(
    sim_run: ReusableSimulationRun,
    output_dir: Path,
) -> pd.DataFrame:
    """Persist one already-solved point using the stable unit-cell result layout."""
    output_dir = Path(output_dir)
    csv_path = output_dir / "eigenfrequencies.csv"
    eigenmodes_dir = output_dir / "eigenmodes"
    eigenmodes_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        sim_run.get_eigenfrequencies(),
        columns=["re", "im", "q"],
    )
    rows = []
    for mode_idx in frame.index:
        valid, peaks = _mode_is_valid(sim_run, int(mode_idx))
        row = {"is_valid": bool(valid), **peaks}
        if valid:
            coordinates_re, hz_re = sim_run.get_2d_fields(
                int(mode_idx), "ewfd.Hz", "center"
            )
            coordinates_im, hz_im = sim_run.get_2d_fields(
                int(mode_idx), "ewfd.Hz*(-i)", "center"
            )
            hz_im_aligned = interpolate_field(coordinates_im, hz_im, coordinates_re)
            field_path = eigenmodes_dir / f"{int(mode_idx):02d}_Hz_center.parquet"
            pd.DataFrame(
                {
                    "x": coordinates_re[:, 0],
                    "y": coordinates_re[:, 1],
                    "re": np.real(hz_re),
                    "im": np.real(hz_im_aligned),
                }
            ).to_parquet(field_path)
            row.update(load_mode_composition(output_dir, int(mode_idx)))
        rows.append(row)

    analysis = pd.DataFrame(rows, index=frame.index)
    frame = pd.concat([frame, analysis], axis=1)
    for column in (
        "s_weight",
        "p_weight",
        "d_weight",
        "f_weight",
        "px_weight",
        "py_weight",
        "dx_weight",
        "dy_weight",
        "dominant_subspace",
        "dominant_mode",
    ):
        if column not in frame:
            frame[column] = np.nan
    write_csv_atomic(csv_path, frame)
    return frame


def select_gamma_bands(gamma_modes: pd.DataFrame) -> dict[str, int]:
    """Select two p-dominant and two d-dominant valid modes at Gamma."""
    required = {"re", "is_valid", "dominant_subspace", "p_weight", "d_weight"}
    missing = required.difference(gamma_modes.columns)
    if missing:
        raise ValueError(f"Gamma mode table is missing columns: {sorted(missing)}")

    eligible = gamma_modes[
        gamma_modes["is_valid"].astype(bool)
        & gamma_modes["re"].between(FREQUENCY_MIN_THZ, FREQUENCY_MAX_THZ)
    ]
    selected = {}
    for subspace in ("p", "d"):
        matches = eligible[eligible["dominant_subspace"] == subspace].sort_values("re")
        if len(matches) < 2:
            diagnostic_columns = [
                "re",
                "is_valid",
                "dominant_subspace",
                "p_weight",
                "d_weight",
            ]
            diagnostic = gamma_modes[diagnostic_columns].to_string()
            raise ValueError(
                f"Gamma selection requires two {subspace}-dominant modes in "
                f"[{FREQUENCY_MIN_THZ:g}, {FREQUENCY_MAX_THZ:g}] THz.\n{diagnostic}"
            )
        for rank, mode_idx in enumerate(matches.index[:2], start=1):
            selected[f"{subspace}{rank}"] = int(mode_idx)
    return selected


def _tracked_mode_row(
    point: dict[str, object],
    point_index: int,
    band_label: str,
    mode_idx: int,
    mode: pd.Series,
    match_status: str,
    frequency_difference: float,
    overlap: float,
) -> dict[str, object]:
    sign = -1.0 if point["direction"] == "gamma_m" else 1.0
    return {
        "direction": point["direction"],
        "point_index": point_index,
        "path_coordinate": sign * float(point["path_fraction"]),
        "path_fraction": float(point["path_fraction"]),
        "kx": float(point["kx"]),
        "ky": float(point["ky"]),
        "k_norm": float(point["k_norm"]),
        "kx_str": point["kx_str"],
        "ky_str": point["ky_str"],
        "band_label": band_label,
        "mode_idx": int(mode_idx),
        "match_status": match_status,
        "frequency_difference": float(frequency_difference),
        "field_overlap": float(overlap),
        **mode.to_dict(),
    }


def _unmatched_mode_row(
    point: dict[str, object], point_index: int, band_label: str
) -> dict[str, object]:
    sign = -1.0 if point["direction"] == "gamma_m" else 1.0
    return {
        "direction": point["direction"],
        "point_index": point_index,
        "path_coordinate": sign * float(point["path_fraction"]),
        "path_fraction": float(point["path_fraction"]),
        "kx": float(point["kx"]),
        "ky": float(point["ky"]),
        "k_norm": float(point["k_norm"]),
        "kx_str": point["kx_str"],
        "ky_str": point["ky_str"],
        "band_label": band_label,
        "mode_idx": np.nan,
        "match_status": "unmatched",
        "frequency_difference": np.nan,
        "field_overlap": np.nan,
    }


def track_selected_branch(
    case_dir: Path,
    points: Sequence[dict[str, object]],
    selected_modes: dict[str, int],
) -> pd.DataFrame:
    """Track the four Gamma-selected modes outward along one branch."""
    if not points:
        raise ValueError("Cannot track an empty k-point branch")
    gamma_dir = k_point_dir(case_dir, points[0])
    gamma_modes = pd.read_csv(gamma_dir / "eigenfrequencies.csv")
    rows = []
    states = {}
    labels = list(selected_modes)
    for label in labels:
        mode_idx = int(selected_modes[label])
        mode = gamma_modes.loc[mode_idx]
        field = load_hz_center(gamma_dir, mode_idx)
        states[label] = {
            "active": True,
            "frequency": float(mode["re"]),
            "field": field,
        }
        rows.append(
            _tracked_mode_row(points[0], 0, label, mode_idx, mode, "gamma", 0.0, 1.0)
        )

    for point_index, point in enumerate(points[1:], start=1):
        point_dir = k_point_dir(case_dir, point)
        modes = pd.read_csv(point_dir / "eigenfrequencies.csv")
        valid_modes = modes[modes["is_valid"].astype(bool)]
        candidate_indices = [int(index) for index in valid_modes.index]
        candidate_fields = {
            mode_idx: load_hz_center(point_dir, mode_idx) for mode_idx in candidate_indices
        }
        active_labels = [label for label in labels if states[label]["active"]]
        if not active_labels:
            rows.extend(
                _unmatched_mode_row(point, point_index, label) for label in labels
            )
            continue
        previous_frequencies = [states[label]["frequency"] for label in active_labels]
        candidate_frequencies = [float(modes.loc[index, "re"]) for index in candidate_indices]
        overlaps = np.array(
            [
                [
                    field_overlap(states[label]["field"], candidate_fields[index])
                    for index in candidate_indices
                ]
                for label in active_labels
            ],
            dtype=float,
        )
        assignments = assign_mode_candidates(
            previous_frequencies,
            candidate_frequencies,
            overlaps,
            BAND_FREQUENCY_TOLERANCE_THZ,
            BAND_MINIMUM_OVERLAP,
        )
        assignment_by_label = dict(zip(active_labels, assignments))

        for label in labels:
            if not states[label]["active"]:
                rows.append(_unmatched_mode_row(point, point_index, label))
                continue
            assignment = assignment_by_label[label]
            if not assignment["matched"]:
                states[label]["active"] = False
                rows.append(_unmatched_mode_row(point, point_index, label))
                continue
            mode_idx = candidate_indices[int(assignment["candidate_index"])]
            mode = modes.loc[mode_idx]
            states[label].update(
                {
                    "frequency": float(mode["re"]),
                    "field": candidate_fields[mode_idx],
                }
            )
            rows.append(
                _tracked_mode_row(
                    point,
                    point_index,
                    label,
                    mode_idx,
                    mode,
                    "matched",
                    float(assignment["frequency_difference"]),
                    float(assignment["overlap"]),
                )
            )
    return pd.DataFrame(rows)


def fit_q_power_law(rows: pd.DataFrame) -> dict[str, float | int | str]:
    """Fit Q proportional to |k|**(-alpha) on the configured dense points."""
    k_values = pd.to_numeric(rows["k_norm"], errors="coerce").to_numpy(dtype=float)
    q_values = pd.to_numeric(rows["q"], errors="coerce").to_numpy(dtype=float)
    configured = np.asarray(Q_FIT_K_MAGNITUDES, dtype=float)
    is_configured = (
        np.any(np.isclose(k_values[:, None], configured[None, :], atol=1e-12), axis=1)
        if configured.size
        else np.zeros(len(k_values), dtype=bool)
    )
    usable = is_configured & np.isfinite(k_values) & (k_values > 0.0)
    usable &= np.isfinite(q_values) & (q_values > 0.0)
    k_fit = k_values[usable]
    q_fit = q_values[usable]
    result = {
        "status": "insufficient_points",
        "point_count": int(len(k_fit)),
        "alpha": np.nan,
        "intercept": np.nan,
        "r_squared": np.nan,
        "k_min": float(np.min(k_fit)) if len(k_fit) else np.nan,
        "k_max": float(np.max(k_fit)) if len(k_fit) else np.nan,
    }
    if len(k_fit) < 3:
        return result

    log_k = np.log(k_fit)
    log_q = np.log(q_fit)
    slope, intercept = np.polyfit(log_k, log_q, 1)
    predicted = slope * log_k + intercept
    residual_sum = float(np.sum((log_q - predicted) ** 2))
    total_sum = float(np.sum((log_q - np.mean(log_q)) ** 2))
    result.update(
        {
            "status": "ok",
            "alpha": float(-slope),
            "intercept": float(intercept),
            "r_squared": float(1.0 - residual_sum / total_sum) if total_sum > 0.0 else 1.0,
        }
    )
    return result


def p_x_fraction(px_weight, py_weight) -> np.ndarray:
    """Return the project-defined p_x share within the p_x/p_y pair."""
    px, py = np.broadcast_arrays(
        np.asarray(px_weight, dtype=float),
        np.asarray(py_weight, dtype=float),
    )
    total = px + py
    fraction = np.full(px.shape, np.nan, dtype=float)
    finite = np.isfinite(px) & np.isfinite(py)
    zero = finite & np.isclose(total, 0.0)
    fraction[zero] = 0.5
    positive = finite & (total > 0.0)
    np.divide(px, total, out=fraction, where=positive)
    return np.clip(fraction, 0.0, 1.0)


def p_band_color(frame: pd.DataFrame, band_label: str) -> str:
    """Select red user-p_y or blue user-p_x from the Gamma composition."""
    gamma = frame[
        (frame["band_label"] == band_label)
        & (frame["match_status"] != "unmatched")
        & np.isclose(pd.to_numeric(frame["path_coordinate"]), 0.0)
    ]
    if gamma.empty:
        raise ValueError(f"Missing Gamma composition for {band_label}")
    fractions = p_x_fraction(
        pd.to_numeric(gamma["px_weight"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(gamma["py_weight"], errors="coerce").to_numpy(dtype=float),
    )
    if not np.all(np.isfinite(fractions)):
        raise ValueError(f"Invalid Gamma composition for {band_label}")
    if not np.allclose(fractions, fractions[0], atol=1e-6):
        raise ValueError(f"Inconsistent Gamma composition for {band_label}")
    return P_X_COLOR if float(fractions[0]) >= 0.5 else P_Y_COLOR


def _reference_p2_gamma_frequency(frame: pd.DataFrame) -> float:
    rows = frame[
        (frame["band_label"] == "p2")
        & (frame["match_status"] != "unmatched")
        & np.isclose(pd.to_numeric(frame["path_coordinate"]), 0.0)
    ]
    frequencies = pd.to_numeric(rows["re"], errors="coerce").dropna().to_numpy(
        dtype=float
    )
    if frequencies.size == 0:
        raise ValueError("Missing cavity p2 Gamma frequency")
    if not np.allclose(frequencies, frequencies[0]):
        raise ValueError("Inconsistent cavity p2 Gamma frequencies")
    return float(frequencies[0])


def plot_band_cases(
    cases: Sequence[tuple[str, pd.DataFrame]],
    output_path: Path,
) -> None:
    """Plot one or more unit-cell band cases using a shared cavity reference."""
    cases = tuple(cases)
    if not cases:
        raise ValueError("At least one unit-cell band case is required")
    frequency_half_span = 0.5 * (FREQUENCY_MAX_THZ - FREQUENCY_MIN_THZ)
    center_frequency = _reference_p2_gamma_frequency(cases[0][1])
    frequency_limits = (
        center_frequency - frequency_half_span,
        center_frequency + frequency_half_span,
    )

    p_color_map = LinearSegmentedColormap.from_list(
        "wave", [P_Y_COLOR, "white", P_X_COLOR]
    )
    d_color_map = LinearSegmentedColormap.from_list(
        "d_character", ["#b0b0b0", "black"]
    )
    normalization = Normalize(0.0, 1.0)
    figure, axes_grid = plt.subplots(
        1,
        len(cases),
        figsize=(5.6 * len(cases), 4.8),
        sharey=True,
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes_grid[0]
    for axis in axes[1:]:
        axis.tick_params(labelleft=False)
    for axis, (title, frame) in zip(axes, cases):
        for label in ("p1", "p2", "d1", "d2"):
            band_rows = frame[frame["band_label"] == label]
            band_rows = band_rows.sort_values("path_coordinate").drop_duplicates(
                "path_coordinate", keep="first"
            )
            x = band_rows["path_coordinate"].to_numpy(dtype=float)
            y = band_rows["re"].to_numpy(dtype=float)
            matched = band_rows["match_status"].to_numpy() != "unmatched"
            if label.startswith("p"):
                character = p_x_fraction(
                    band_rows["px_weight"].to_numpy(dtype=float),
                    band_rows["py_weight"].to_numpy(dtype=float),
                )
                color_map = p_color_map
            else:
                first_weight = band_rows["dy_weight"].to_numpy(dtype=float)
                second_weight = band_rows["dx_weight"].to_numpy(dtype=float)
                color_map = d_color_map
                total_weight = first_weight + second_weight
                character = np.full(total_weight.shape, 0.5, dtype=float)
                np.divide(
                    second_weight,
                    total_weight,
                    out=character,
                    where=total_weight > 0.0,
                )
                character = np.clip(character, 0.0, 1.0)
            if len(x) >= 2:
                points = np.column_stack([x, y])
                segments = np.stack([points[:-1], points[1:]], axis=1)
                segment_weights = 0.5 * (character[:-1] + character[1:])
                usable = matched[:-1] & matched[1:]
                usable &= np.all(np.isfinite(segments), axis=(1, 2))
                usable &= np.isfinite(segment_weights)
                segments = segments[usable]
                segment_weights = segment_weights[usable]
                if not len(segments):
                    continue
                collection = LineCollection(
                    segments,
                    cmap=color_map,
                    norm=normalization,
                    linewidth=2.0,
                )
                collection.set_array(segment_weights)
                axis.add_collection(collection)
        axis.set_xlim(-1.0, 1.0)
        axis.set_ylim(*frequency_limits)
        axis.set_xticks([-1.0, 0.0, 1.0], labels=["M", "\u0393", "K"])
        axis.set_xlabel("Wave-vector path")
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
        axis.axvline(0.0, color="0.7", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("Frequency (THz)")
    legend_handles = [
        Line2D([], [], color=P_X_COLOR, linewidth=2.0, label=r"$p_y$"),
        Line2D([], [], color=P_Y_COLOR, linewidth=2.0, label=r"$p_x$"),
        Line2D([], [], color="#b0b0b0", linewidth=2.0, label=r"$d_{xy}$"),
        Line2D([], [], color="black", linewidth=2.0, label=r"$d_{x^2+y^2}$"),
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.96),
        fontsize=8,
        framealpha=0.9,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def plot_band_comparison(
    cavity: pd.DataFrame,
    cladding: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the cavity/cladding pair using the shared multi-case renderer."""
    plot_band_cases(
        (("Cavity cell", cavity), ("Cladding cell", cladding)),
        output_path,
    )


def plot_p_band_q_cases(
    cases: Sequence[tuple[str, pd.DataFrame, pd.DataFrame]],
    output_path: Path,
) -> None:
    """Plot one or more unit-cell p-band Q cases on signed path axes."""
    cases = tuple(cases)
    if not cases:
        raise ValueError("At least one unit-cell Q case is required")
    figure, axes_grid = plt.subplots(
        1,
        len(cases),
        figsize=(5.6 * len(cases), 4.8),
        sharex=True,
        sharey=True,
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes_grid[0]
    direction_styles = {"gamma_m": "-", "gamma_k": "--"}
    direction_names = {"gamma_m": "\u0393\u2013M", "gamma_k": "\u0393\u2013K"}
    for axis, (title, frame, fits) in zip(axes, cases):
        annotations = {"gamma_m": [], "gamma_k": []}
        band_colors = {
            label: p_band_color(frame, label) for label in ("p1", "p2")
        }
        for direction in ("gamma_m", "gamma_k"):
            sign = -1.0 if direction == "gamma_m" else 1.0
            for label in ("p1", "p2"):
                rows = frame[
                    (frame["direction"] == direction)
                    & (frame["band_label"] == label)
                    & (frame["match_status"] != "unmatched")
                    & (frame["k_norm"] <= 0.2)
                ].sort_values("k_norm")
                k_values = rows["k_norm"].to_numpy(dtype=float)
                q_values = rows["q"].to_numpy(dtype=float)
                usable = np.isfinite(q_values) & (q_values > 0.0)
                axis.plot(
                    sign * k_values[usable],
                    np.log10(q_values[usable]),
                    marker="o",
                    markersize=3,
                    linewidth=1.2,
                    linestyle=direction_styles[direction],
                    color=band_colors[label],
                    label=f"{label} {direction_names[direction]}",
                )
                fit_rows = fits[
                    (fits["direction"] == direction)
                    & (fits["band_label"] == label)
                ]
                if fit_rows.empty or fit_rows.iloc[0]["status"] != "ok":
                    annotations[direction].append(
                        f"{direction_names[direction]} {label}: fit unavailable"
                    )
                    continue
                fit = fit_rows.iloc[0]
                k_fit_max = min(float(fit["k_max"]), 0.2)
                k_fit = np.geomspace(float(fit["k_min"]), k_fit_max, 100)
                q_fit = np.exp(float(fit["intercept"])) * k_fit ** (
                    -float(fit["alpha"])
                )
                axis.plot(
                    sign * k_fit,
                    np.log10(q_fit),
                    color=band_colors[label],
                    linestyle=":",
                    linewidth=1.0,
                )
                annotations[direction].append(
                    f"{direction_names[direction]} {label}: \u03b1={fit['alpha']:.3f}, "
                    f"R\u00b2={fit['r_squared']:.3f}\n"
                    f"k/G={fit['k_min']:.3g}\u2013{fit['k_max']:.3g}"
                )
        axis.set_xlim(-0.2, 0.2)
        axis.axvline(0.0, color="0.65", linewidth=0.9, linestyle="--")
        axis.set_title(title)
        axis.set_xlabel("Signed |k|/G  (\u0393\u2013M < 0, \u0393\u2013K > 0)")
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=7, ncol=2, loc="lower center")
        for direction, x_position, alignment in (
            ("gamma_m", 0.02, "left"),
            ("gamma_k", 0.98, "right"),
        ):
            axis.text(
                x_position,
                0.98,
                "\n".join(annotations[direction]),
                transform=axis.transAxes,
                ha=alignment,
                va="top",
                fontsize=7,
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
            )
    axes[0].set_ylabel("log10(Q)")
    for axis in axes[1:]:
        axis.tick_params(labelleft=False)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def plot_p_band_q_comparison(
    cavity: pd.DataFrame,
    cladding: pd.DataFrame,
    cavity_fits: pd.DataFrame,
    cladding_fits: pd.DataFrame,
    output_path: Path,
) -> None:
    """Compare cavity and cladding p-band Q on shared axes."""
    plot_p_band_q_cases(
        (
            ("Cavity cell", cavity, cavity_fits),
            ("Cladding cell", cladding, cladding_fits),
        ),
        output_path,
    )
def _unique_case_points(
    branches: dict[str, list[dict[str, float | str]]],
) -> list[dict[str, float | str]]:
    unique = {}
    for points in branches.values():
        for point in points:
            unique.setdefault((point["kx_str"], point["ky_str"]), point)
    return list(unique.values())


def _mode_composition_table(
    case_name: str,
    case_dir: Path,
    points: Sequence[dict[str, object]],
) -> pd.DataFrame:
    frames = []
    for point in points:
        modes = pd.read_csv(k_point_dir(case_dir, point) / "eigenfrequencies.csv")
        modes.insert(0, "mode_idx", modes.index.astype(int))
        modes.insert(0, "ky", float(point["ky"]))
        modes.insert(0, "kx", float(point["kx"]))
        modes.insert(0, "direction", point["direction"])
        modes.insert(0, "cell", case_name)
        frames.append(modes)
    return pd.concat(frames, ignore_index=True)


def summarize_band_tracking(selected: pd.DataFrame) -> dict[str, object]:
    """Summarize incomplete tracking without treating missing solver modes as fatal."""
    unmatched = selected[selected["match_status"] == "unmatched"]
    detail_columns = [
        "direction",
        "point_index",
        "band_label",
        "kx",
        "ky",
        "k_norm",
    ]
    details = unmatched.reindex(columns=detail_columns).to_dict(orient="records")
    return {
        "status": "complete" if unmatched.empty else "partial",
        "selected_row_count": int(len(selected)),
        "matched_row_count": int(len(selected) - len(unmatched)),
        "unmatched_row_count": int(len(unmatched)),
        "unmatched": details,
    }


def run_case(
    case_name: str,
    params: dict[str, float],
    output_root: Path | None = None,
    branches: dict[str, list[dict[str, float | str]]] | None = None,
) -> pd.DataFrame:
    """Run, track, and save one cavity or cladding unit-cell band case."""
    root = OUT_DIR if output_root is None else Path(output_root)
    paths = unit_cell_case_paths(root, case_name)
    case_dir = paths["result"]
    overview_dir = paths["overview"]
    config_dir = paths["config"]
    case_dir.mkdir(parents=True, exist_ok=True)
    overview_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)
    hole_params = get_hole_params(params)
    hexagon, holes, info = create_hexagon_design(A, hole_params)
    min_distance = float(min(item["min_dist"] for item in info))
    if min_distance < D:
        raise ValueError(
            f"{case_name} unit cell violates clearance: min_dist={min_distance:g}, d={D:g}"
        )
    visualize_hexagon_design(
        hexagon,
        holes,
        info,
        filename=overview_dir / "geometry.png",
        annotate=False,
    )

    branches = case_k_points(case_name) if branches is None else branches
    if not branches:
        raise ValueError("At least one k-space branch is required")
    unique_points = _unique_case_points(branches)
    simulation_config = SimulationConfig(
        slab_height=SLAB_HEIGHT,
        slab_refractive_index=REFRACTIVE_INDEX,
        eigenfrequency_shift=EIGENFREQUENCY_SHIFT,
        eigenmode_count=EIGENMODE_COUNT,
        mesh_auto_size=MESH_AUTO_SIZE,
        mode_type=MODE_TYPE,
        wavelength=WAVELENGTH,
    )
    hole_vertices = [hole.tolist() for hole in holes]
    sim_run = None
    cache_hits = 0
    cache_misses = 0
    try:
        for point_index, point in enumerate(unique_points, start=1):
            output_dir = k_point_dir(case_dir, point)
            is_complete = k_point_is_complete(output_dir, simulation_config)
            state = "cached" if is_complete else "run"
            print(
                f"[{case_name}] {point_index}/{len(unique_points)} "
                f"kx={point['kx_str']} ky={point['ky_str']} ({state})"
            )
            if is_complete:
                cache_hits += 1
            else:
                cache_misses += 1
                if sim_run is None:
                    sim_run = ReusableSimulationRun(simulation_config)
            run_k_point(
                case_dir,
                hole_vertices,
                point,
                simulation_config,
                sim_run=sim_run,
            )
        runner_stats = (
            sim_run.execution_stats
            if sim_run is not None
            else {
                "execution_engine": "reused_model",
                "model_reuse_enabled": True,
                "model_build_count": 0,
                "mesh_build_count": 0,
                "solve_count": 0,
                "mesh_auto_size": MESH_AUTO_SIZE,
                "current_k": None,
            }
        )
    finally:
        if sim_run is not None:
            sim_run.clear()

    execution = {
        **runner_stats,
        "reuse_scope": "one COMSOL Model per fixed unit-cell geometry",
        "cache_hits": cache_hits,
        "cache_misses": cache_misses,
    }

    first_branch = next(iter(branches.values()))
    if not first_branch or float(first_branch[0]["k_norm"]) != 0.0:
        raise ValueError("Every unit-cell dataset must begin with a Gamma point")
    gamma_dir = k_point_dir(case_dir, first_branch[0])
    selected_modes = select_gamma_bands(pd.read_csv(gamma_dir / "eigenfrequencies.csv"))
    tracked = []
    for direction in branches:
        branch = track_selected_branch(case_dir, branches[direction], selected_modes)
        branch.insert(0, "cell", case_name)
        write_csv_atomic(overview_dir / f"connected_bands_{direction}.csv", branch)
        tracked.append(branch)
    selected = pd.concat(tracked, ignore_index=True)
    write_csv_atomic(overview_dir / "selected_bands.csv", selected)
    tracking = summarize_band_tracking(selected)

    compositions = _mode_composition_table(case_name, case_dir, unique_points)
    write_csv_atomic(overview_dir / "mode_composition.csv", compositions)
    write_json_atomic(
        config_dir / "config.json",
        {
            "case": case_name,
            "params": params,
            "hole_params": hole_params,
            "triangles": holes,
            "minimum_clearance": min_distance,
            "selected_gamma_modes": selected_modes,
            "branches": branches,
            "execution": execution,
            "tracking": tracking,
        },
    )
    if tracking["status"] == "partial":
        unmatched = selected[selected["match_status"] == "unmatched"]
        details = unmatched[["direction", "point_index", "band_label"]].to_string(
            index=False
        )
        print(
            f"[{case_name}] WARNING: {tracking['unmatched_row_count']} selected-band "
            "rows are unmatched because the requested eigenspectrum did not contain "
            f"a valid continuation. Continuing with partial tracking:\n{details}"
        )
    return selected


def _build_q_fits(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for direction in ("gamma_m", "gamma_k"):
        for label in ("p1", "p2"):
            data = frame[
                (frame["direction"] == direction) & (frame["band_label"] == label)
            ]
            rows.append(
                {
                    "direction": direction,
                    "band_label": label,
                    **fit_q_power_law(data),
                }
            )
    return pd.DataFrame(rows)


def run_band_workflow(
    output_dir: Path,
    cavity_params: dict[str, float],
    cladding_params: dict[str, float],
) -> dict[str, object]:
    """Run both unit cells and write the complete comparison output set."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    branches_by_case = {name: case_k_points(name) for name in ("cavity", "cladding")}
    write_json_atomic(
        output_dir / CONFIG_DIRNAME / "config.json",
        {
            "execution_engine": "reused_model",
            "model_reuse": {
                "scope": "one COMSOL Model per fixed unit-cell geometry",
                "same_mesh": "reuse without mesh rebuild",
                "changed_mesh": "rebuild mesh1 only",
            },
            "geometry": {"a": A, "r_0": R_0, "b_0": B_0, "d": D},
            "simulation": {
                "slab_height": SLAB_HEIGHT,
                "refractive_index": REFRACTIVE_INDEX,
                "eigenmode_count": EIGENMODE_COUNT,
                "mesh_auto_size": MESH_AUTO_SIZE,
                "mode_type": MODE_TYPE,
                "wavelength": WAVELENGTH,
                "eigenfrequency_shift": EIGENFREQUENCY_SHIFT,
            },
            "frequency_window_thz": [FREQUENCY_MIN_THZ, FREQUENCY_MAX_THZ],
            "band_points_per_arm": BAND_POINTS_PER_ARM,
            "q_fit_k_magnitudes": Q_FIT_K_MAGNITUDES,
            "band_frequency_tolerance_thz": BAND_FREQUENCY_TOLERANCE_THZ,
            "band_minimum_overlap": BAND_MINIMUM_OVERLAP,
            "cases": {"cavity": cavity_params, "cladding": cladding_params},
            "shared_parameters": ACTIVE_PARAMETERS.to_metadata(),
            "sampled_points": branches_by_case,
        },
    )
    cavity = run_case("cavity", cavity_params, output_dir)
    cladding = run_case("cladding", cladding_params, output_dir)
    cavity_fits = _build_q_fits(cavity)
    cladding_fits = _build_q_fits(cladding)
    cavity_overview = unit_cell_case_paths(output_dir, "cavity")["overview"]
    cladding_overview = unit_cell_case_paths(output_dir, "cladding")["overview"]
    overview_dir = output_dir / OVERVIEW_DIRNAME
    write_csv_atomic(cavity_overview / "q_power_law_fits.csv", cavity_fits)
    write_csv_atomic(cladding_overview / "q_power_law_fits.csv", cladding_fits)
    plot_band_comparison(
        cavity,
        cladding,
        overview_dir / "unit_cell_band_comparison.png",
    )
    plot_p_band_q_comparison(
        cavity,
        cladding,
        cavity_fits,
        cladding_fits,
        overview_dir / "p_bands_q_vs_k_comparison.png",
    )
    return {
        "output_dir": output_dir,
        "cavity": cavity,
        "cladding": cladding,
        "cavity_fits": cavity_fits,
        "cladding_fits": cladding_fits,
    }


def cavity_p2_gamma_frequency(cavity: pd.DataFrame) -> float:
    """Return the unique cavity p2 frequency at Gamma from tracked bands."""
    rows = cavity[
        (cavity["band_label"] == "p2")
        & np.isclose(pd.to_numeric(cavity["k_norm"], errors="coerce"), 0.0)
        & (cavity["match_status"] == "gamma")
    ]
    frequencies = pd.to_numeric(rows["re"], errors="coerce").dropna().to_numpy(
        dtype=float
    )
    if frequencies.size == 0 or not np.all(np.isfinite(frequencies)):
        raise RuntimeError("Complete unit-cell result has no valid cavity p2@Gamma frequency")
    if not np.allclose(frequencies, frequencies[0], rtol=0.0, atol=1e-9):
        raise RuntimeError(
            f"Cavity p2@Gamma frequency is inconsistent across branches: {frequencies}"
        )
    return float(frequencies[0])


def main() -> None:
    result = run_band_workflow(OUT_DIR, BULK_PARAMS, CLADDING_PARAMS)
    if UPDATE_CENTER_FREQUENCY_AFTER_RUN:
        frequency_thz = cavity_p2_gamma_frequency(result["cavity"])
        updated = update_center_frequency_from_cavity_p2(
            frequency_thz,
            expected_parameters=ACTIVE_PARAMETERS,
        )
        print(
            "Updated shared center frequency from cavity p2@Gamma: "
            f"{updated.center_frequency_thz:.12g} THz -> {updated.source_path}"
        )


if __name__ == "__main__":
    main()
