#!/usr/bin/env python3
"""Solve all four finite-quarter symmetry sectors in one COMSOL model."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_finite as finite  # noqa: E402
from scripts.run_main import run_finite_quarter as quarter  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

from comsol_workflow.simulation_utils import BoundarySpec, SimulationRun  # noqa: E402


SYMMETRY_IDS = (1, 2, 3, 4)
MODEL_FILENAME = "finite_quarter_all.mph"
SUMMARY_FILENAME = "run_summary.json"
AGGREGATE_TABLES = (
    "eigenfrequencies.csv",
    "farfield_summary.csv",
    "finite_lattice_fourier_summary.csv",
    "mode_scores.csv",
)
SHARED_OVERVIEW_FILES = ("index_r_cavity.png", "index_k_1stBZ.png")
RUN_COMSOL = True
REUSE_MESHED_MODEL: Path | None = None
RESUME_INCOMPLETE = False


def all_symmetry_series_dir(active_parameters=finite.ACTIVE_PARAMETERS) -> Path:
    quarter_label = active_parameters.structure_series_label("finite_quarter")
    if not quarter_label.startswith("finite_quarter_"):
        raise ValueError(
            "Unexpected finite-quarter series label: "
            f"{quarter_label!r}"
        )
    return finite.OUTPUT_ROOT / quarter_label.replace(
        "finite_quarter_",
        "finite_quarter_all_",
        1,
    )


OUT_DIR = all_symmetry_series_dir()


def symmetry_suffix(symmetry_id: int) -> str:
    case = quarter.symmetry_case(symmetry_id)
    return (
        f"_{symmetry_id}"
        f"_x{case['x_boundary']}"
        f"_y{case['y_boundary']}"
    )


def mode_uid(mode_idx: int, symmetry_id: int) -> str:
    return f"mode{int(mode_idx)}{symmetry_suffix(symmetry_id)}"


def solution_identity(symmetry_id: int) -> dict[str, str]:
    return {
        "solution_tag": f"solsym{int(symmetry_id)}",
        "dataset_tag": f"dsetsym{int(symmetry_id)}",
    }


def _partition_dir(root_staging_dir: Path, symmetry_id: int) -> Path:
    return root_staging_dir / f"s{int(symmetry_id)}"


def validate_runtime_config() -> None:
    if type(finite.EIGENMODE_COUNT) is not int or finite.EIGENMODE_COUNT < 1:
        raise ValueError("finite_eigenmode_count must be a positive integer")
    if finite.RESUME_FROM_EXISTING_MPH or quarter.RESUME_FROM_EXISTING_MPH:
        raise ValueError(
            "finite-quarter all-symmetry currently requires fresh output; "
            "disable RESUME_FROM_EXISTING_MPH"
        )
    if RESUME_INCOMPLETE and REUSE_MESHED_MODEL is not None:
        raise ValueError(
            "--resume-incomplete loads the target MPH directly; "
            "do not also pass --reuse-meshed-model"
        )
    if REUSE_MESHED_MODEL is not None:
        pairs = finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
        if len(pairs) != 1:
            raise ValueError(
                "--reuse-meshed-model requires exactly one finite shift pair"
            )
        source_info = reuse_model_source_info(REUSE_MESHED_MODEL, *pairs[0])
        source_path = Path(source_info["mph_path"])
        for case_dir in expected_case_directories():
            target = finite.finite_case_output_paths(
                case_dir,
                model_filename=MODEL_FILENAME,
            )["mph"].resolve()
            if target == source_path:
                raise ValueError("Reused source MPH and target MPH must differ")


def _reuse_config_path(model_path: Path) -> Path:
    if model_path.parent.name != finite.MODEL_DIRNAME:
        raise ValueError(
            "Reused MPH must be stored in a standard 00_model directory"
        )
    return model_path.parent.parent / finite.CONFIG_DIRNAME / "config.json"


def reuse_model_source_info(
    model_path: Path,
    x_shift_factor: float,
    y_shift_factor: float,
) -> dict[str, object]:
    source = Path(model_path).resolve()
    if not source.is_file() or source.suffix.lower() != ".mph":
        raise FileNotFoundError(f"Reused finite-quarter MPH not found: {source}")
    config_path = _reuse_config_path(source)
    if not config_path.is_file():
        raise FileNotFoundError(
            f"Reused finite-quarter config not found: {config_path}"
        )
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    shared = payload.get("shared_parameters", {})
    active = finite.ACTIVE_PARAMETERS.to_metadata()
    checks = {
        "finite_geometry": (payload.get("finite_geometry"), finite.FINITE_GEOMETRY),
        "cavity_layers": (payload.get("cavity_layers"), finite.CAVITY_LAYERS),
        "cladding_layers": (payload.get("cladding_layers"), finite.CLADDING_LAYERS),
        "cavity": (shared.get("cavity"), active["cavity"]),
        "cladding": (shared.get("cladding"), active["cladding"]),
        "mesh_size": (shared.get("mesh_size"), finite.MESH_AUTO_SIZE),
        "cladding_shift_geometry": (
            payload.get("cladding_shift_geometry"),
            active["cladding_shift_geometry"],
        ),
        "cladding_shift_profile": (
            payload.get("cladding_shift_profile"),
            active["cladding_shift_profile"],
        ),
        "cladding_x_shift_factor": (
            payload.get("cladding_x_shift_factor"),
            float(x_shift_factor),
        ),
        "cladding_y_shift_factor": (
            payload.get("cladding_y_shift_factor"),
            float(y_shift_factor),
        ),
    }
    mismatches = [
        f"{name}: source={actual!r}, requested={expected!r}"
        for name, (actual, expected) in checks.items()
        if actual != expected
    ]
    if mismatches:
        raise ValueError(
            "Reused finite-quarter model identity mismatch: "
            + "; ".join(mismatches)
        )
    source_symmetry = payload.get("symmetry")
    if not isinstance(source_symmetry, dict):
        raise ValueError("Reused finite-quarter config is missing symmetry identity")
    stat = source.stat()
    return {
        "mph_path": str(source),
        "mph_size_bytes": int(stat.st_size),
        "mph_last_write_ns": int(stat.st_mtime_ns),
        "config_path": str(config_path),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "source_symmetry": source_symmetry,
        "source_eigenfrequency_shift": payload.get("simulation_common", {}).get(
            "eigenfrequency_shift"
        ),
        "source_eigenmode_count": payload.get("simulation_common", {}).get(
            "eigenmode_count"
        ),
        "quarter_hole_count": payload.get("quarter_hole_count"),
    }


def case_dir_for_shift(
    x_shift_factor: float,
    y_shift_factor: float,
) -> Path:
    return OUT_DIR / finite.case_stem(x_shift_factor, y_shift_factor)


def expected_case_directories() -> list[Path]:
    return [
        case_dir_for_shift(x_factor, y_factor)
        for x_factor, y_factor in (
            finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
        )
    ]


def validate_output_targets() -> None:
    for (x_factor, y_factor), case_dir in zip(
        finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs,
        expected_case_directories(),
    ):
        paths = finite.finite_case_output_paths(
            case_dir,
            model_filename=MODEL_FILENAME,
        )
        if RESUME_INCOMPLETE:
            _load_resume_state(paths, x_factor, y_factor)
            continue
        occupied = [
            path
            for path in (
                paths["mph"],
                paths["config"],
                paths["staging_dir"],
            )
            if path.exists()
        ]
        if paths["results_dir"].is_dir() and any(
            paths["results_dir"].iterdir()
        ):
            occupied.append(paths["results_dir"])
        if occupied:
            raise FileExistsError(
                "Refusing to overwrite an existing all-symmetry case: "
                f"{occupied[0]}"
            )


def _validate_partition_exports(
    partition: dict[str, Path],
    symmetry_id: int,
) -> pd.DataFrame:
    eigen_path = partition["export_dir"] / "eigenfrequencies.csv"
    if not eigen_path.is_file() or not partition["config"].is_file():
        raise FileNotFoundError(
            f"Incomplete saved symmetry partition {symmetry_id}: {partition['case_dir']}"
        )
    frame = pd.read_csv(eigen_path)
    required = {
        "re", "im", "q", "is_valid", "mode_idx",
        "symmetry_id", "x_boundary", "y_boundary",
    }
    if missing := required.difference(frame.columns):
        raise ValueError(
            f"Saved symmetry {symmetry_id} table is missing columns {sorted(missing)}"
        )
    mode_indices = pd.to_numeric(frame["mode_idx"], errors="raise").to_numpy(float)
    if (
        len(frame) < finite.EIGENMODE_COUNT
        or not np.array_equal(mode_indices, np.arange(len(frame), dtype=float))
    ):
        raise ValueError(
            f"Saved symmetry {symmetry_id} must contain at least "
            f"{finite.EIGENMODE_COUNT} modes numbered continuously from zero"
        )
    case = quarter.symmetry_case(symmetry_id)
    if (
        set(frame["symmetry_id"].astype(int)) != {symmetry_id}
        or set(frame["x_boundary"].astype(str)) != {case["x_boundary"]}
        or set(frame["y_boundary"].astype(str)) != {case["y_boundary"]}
    ):
        raise ValueError(f"Saved symmetry {symmetry_id} identity does not match")
    expected_exports = (
        "Hz_center.parquet",
        "Hz_Im_2d.png",
        "Hz_Re_2d.png",
        "Wem_2d.png",
    )
    if finite.RUN_FARFIELD_FFT:
        expected_exports = ("E_air.parquet", *expected_exports)
    valid_indices = frame.loc[
        frame["is_valid"].map(finite._csv_bool), "mode_idx"
    ].astype(int)
    for mode_idx in valid_indices:
        for suffix in expected_exports:
            path = partition["export_dir"] / f"{mode_idx:02d}_{suffix}"
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"Incomplete saved symmetry {symmetry_id} export: {path}"
                )
    air_metadata = partition["export_dir"] / "air_field_metadata.json"
    if finite.RUN_FARFIELD_FFT and not air_metadata.is_file():
        raise FileNotFoundError(
            f"Missing saved symmetry {symmetry_id} air metadata: {air_metadata}"
        )
    return frame


def _validate_resume_identity(
    metadata: dict[str, object],
    x_shift_factor: float,
    y_shift_factor: float,
) -> None:
    active = finite.ACTIVE_PARAMETERS.to_metadata()
    identity_checks = {
        "case": (metadata.get("case"), "finite_quarter_all"),
        "finite_geometry": (metadata.get("finite_geometry"), finite.FINITE_GEOMETRY),
        "cavity_layers": (metadata.get("cavity_layers"), finite.CAVITY_LAYERS),
        "cladding_layers": (metadata.get("cladding_layers"), finite.CLADDING_LAYERS),
        "cavity": (metadata.get("shared_parameters", {}).get("cavity"), active["cavity"]),
        "cladding": (
            metadata.get("shared_parameters", {}).get("cladding"),
            active["cladding"],
        ),
        "mesh_size": (
            metadata.get("shared_parameters", {}).get("mesh_size"),
            finite.MESH_AUTO_SIZE,
        ),
        "center_frequency_thz": (
            metadata.get("shared_parameters", {}).get("center_frequency_thz"),
            finite.ACTIVE_PARAMETERS.center_frequency_thz,
        ),
        "eigenmode_count": (
            metadata.get("simulation", {}).get("eigenmode_count"),
            finite.EIGENMODE_COUNT,
        ),
        "cladding_x_shift_factor": (
            metadata.get("cladding_x_shift_factor"),
            float(x_shift_factor),
        ),
        "cladding_y_shift_factor": (
            metadata.get("cladding_y_shift_factor"),
            float(y_shift_factor),
        ),
        "output_layout": (metadata.get("output_layout"), finite.OUTPUT_LAYOUT_VERSION),
    }
    mismatches = [
        f"{name}: saved={saved!r}, active={active_value!r}"
        for name, (saved, active_value) in identity_checks.items()
        if saved != active_value
    ]
    if mismatches:
        raise ValueError("Resume config identity mismatch: " + "; ".join(mismatches))


def _load_resume_state(
    paths: dict[str, Path],
    x_shift_factor: float,
    y_shift_factor: float,
) -> tuple[dict[str, object], dict[str, object], list[dict[str, Path]]]:
    required = (paths["mph"], paths["config"], paths["staging_dir"])
    missing = [path for path in required if not path.exists()]
    summary_path = paths["config_dir"] / SUMMARY_FILENAME
    if not summary_path.is_file():
        missing.append(summary_path)
    if missing:
        raise FileNotFoundError(f"Cannot resume incomplete case; missing {missing[0]}")
    if paths["results_dir"].is_dir() and any(paths["results_dir"].iterdir()):
        raise ValueError("Cannot resume after final 01_results has been populated")

    metadata = json.loads(paths["config"].read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    _validate_resume_identity(metadata, x_shift_factor, y_shift_factor)
    if summary.get("status") not in {"running", "failed"}:
        raise ValueError(
            f"Cannot resume case with status {summary.get('status')!r}"
        )
    runs = summary.get("symmetry_runs_internal")
    if not isinstance(runs, list):
        raise ValueError("Resume summary is missing symmetry_runs_internal")
    completed_ids = [int(run.get("id", -1)) for run in runs]
    if completed_ids != list(SYMMETRY_IDS[: len(completed_ids)]):
        raise ValueError("Completed symmetry IDs must be a continuous prefix")
    if (
        int(summary.get("solve_count", -1)) != len(completed_ids)
        or int(summary.get("solution_copy_count", -1)) != len(completed_ids)
    ):
        raise ValueError("Resume solve/copy counts do not match completed symmetries")
    partition_dirs = {
        path.name for path in paths["staging_dir"].iterdir() if path.is_dir()
    }
    expected_dirs = {f"s{value}" for value in completed_ids}
    if partition_dirs != expected_dirs:
        raise ValueError(
            "Resume staging partitions do not match completed symmetries: "
            f"saved={sorted(partition_dirs)}, expected={sorted(expected_dirs)}"
        )
    partitions: list[dict[str, Path]] = []
    reference_indices: list[int] | None = None
    for run, symmetry_id in zip(runs, completed_ids):
        case = quarter.symmetry_case(symmetry_id)
        archive = solution_identity(symmetry_id)
        if (
            run.get("status") != "solution_copied_exported_verified"
            or run.get("x_boundary") != case["x_boundary"]
            or run.get("y_boundary") != case["y_boundary"]
            or run.get("solution_tag") != archive["solution_tag"]
            or run.get("dataset_tag") != archive["dataset_tag"]
        ):
            raise ValueError(
                f"Resume summary identity mismatch for symmetry {symmetry_id}"
            )
        partition = finite.finite_case_output_paths(
            _partition_dir(paths["staging_dir"], symmetry_id),
            model_filename="internal_partition.mph",
        )
        frame = _validate_partition_exports(partition, symmetry_id)
        indices = frame["mode_idx"].astype(int).tolist()
        if reference_indices is None:
            reference_indices = indices
        elif indices != reference_indices:
            raise ValueError("Saved symmetry partitions have different mode indices")
        if int(run.get("mode_count", -1)) != len(frame):
            raise ValueError(
                f"Resume summary mode count mismatch for symmetry {symmetry_id}"
            )
        partitions.append(partition)
    return metadata, summary, partitions


def validate_cavity_edge_cell_count(metadata: dict[str, object]) -> None:
    cavity_layers = finite.CAVITY_LAYERS
    points = list(metadata["bulk_points"])
    shells = {int(point.shell) for point in points}
    expected_shells = set(range(cavity_layers))
    if shells != expected_shells:
        raise ValueError(
            "Finite-quarter cavity shell mismatch: "
            f"expected {sorted(expected_shells)}, got {sorted(shells)}"
        )
    expected_radius = cavity_layers - 1
    edge_count = sum(
        1
        for point in points
        if point.i >= 0
        and point.j >= 0
        and point.i + point.j == expected_radius
    )
    if edge_count != cavity_layers:
        raise ValueError(
            "Finite-quarter cavity edge cell mismatch: "
            f"expected {cavity_layers}, got {edge_count}"
        )


def _detect_comsol_processes() -> list[str]:
    if os.name != "nt":
        return []
    try:
        result = subprocess.run(
            ["tasklist", "/FO", "CSV", "/NH"],
            capture_output=True,
            text=True,
            check=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    return [
        line
        for line in result.stdout.splitlines()
        if any(token in line.lower() for token in ("comsol", "mphserver"))
    ]


def print_preflight() -> None:
    active = finite.ACTIVE_PARAMETERS
    reuse_info = None
    if REUSE_MESHED_MODEL is not None:
        reuse_info = reuse_model_source_info(
            REUSE_MESHED_MODEL,
            *active.finite_cladding_shift_factor_pairs[0],
        )
    print("Finite-quarter all-symmetry preflight")
    print(f"  cavity: {active.cavity.to_dict()}")
    print(f"  cladding: {active.cladding.to_dict()}")
    print(
        f"  geometry={finite.FINITE_GEOMETRY}, "
        f"layers={finite.CAVITY_LAYERS}-{finite.CLADDING_LAYERS}, "
        f"mesh={finite.MESH_AUTO_SIZE}"
    )
    print(
        f"  center_frequency_thz={active.center_frequency_thz}, "
        f"requested_eigenmodes_per_symmetry={finite.EIGENMODE_COUNT}, "
        f"minimum_requested_total_modes="
        f"{len(SYMMETRY_IDS) * finite.EIGENMODE_COUNT}"
    )
    print(
        "  symmetry order: "
        + ", ".join(
            (
                f"{case['id']}:x{case['x_boundary']}/"
                f"y{case['y_boundary']}"
            )
            for case in map(quarter.symmetry_case, SYMMETRY_IDS)
        )
    )
    if RESUME_INCOMPLETE:
        paths = finite.finite_case_output_paths(
            expected_case_directories()[0],
            model_filename=MODEL_FILENAME,
        )
        _metadata, resume_summary, _partitions = _load_resume_state(
            paths,
            *active.finite_cladding_shift_factor_pairs[0],
        )
        completed = [
            int(run["id"])
            for run in resume_summary["symmetry_runs_internal"]
        ]
        pending = list(SYMMETRY_IDS[len(completed):])
        print(
            f"  execution: resume target MPH, completed={completed}, "
            f"pending={pending}"
        )
    elif reuse_info is None:
        print("  execution: model=1, mesh=1, solve=4, solution copies=4")
    else:
        print(
            "  execution: model=0, mesh=0, reused model=1, reused mesh=1, "
            "solve=4, solution copies=4"
        )
        print(f"  reused MPH: {reuse_info['mph_path']}")
        print(f"    bytes={reuse_info['mph_size_bytes']}")
        print(f"    config={reuse_info['config_path']}")
        print(f"    config_sha256={reuse_info['config_sha256']}")
        print(
            "    old solver: "
            f"shift={reuse_info['source_eigenfrequency_shift']}, "
            f"eigenmodes={reuse_info['source_eigenmode_count']}"
        )
    for (x_factor, y_factor), case_dir in zip(
        active.finite_cladding_shift_factor_pairs,
        expected_case_directories(),
    ):
        print(
            f"  shift x={x_factor:.6g}, y={y_factor:.6g}: "
            f"{case_dir}"
        )
        print(f"    MPH: {case_dir / '00_model' / MODEL_FILENAME}")
    print(f"  far-field={finite.RUN_FARFIELD_FFT}")
    print(
        "  finite-lattice Fourier="
        f"{finite.RUN_FINITE_LATTICE_FOURIER_POSTPROCESS}"
    )
    comsol_path = shutil.which("comsol")
    print(f"  COMSOL executable on PATH: {comsol_path or 'not found'}")
    processes = _detect_comsol_processes()
    print(f"  active COMSOL/mphserver processes: {len(processes)}")
    for process in processes:
        print(f"    {process}")
    print("  license: deferred to COMSOL startup; not consumed by preflight")


def preflight_geometry() -> None:
    for x_factor, y_factor in (
        finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
    ):
        finite.set_cladding_inward_shift_factors(x_factor, y_factor)
        _, quarter_records, metadata, _, _ = (
            quarter.build_quarter_case_geometry()
        )
        validate_cavity_edge_cell_count(metadata)
        print(
            "  geometry validated: "
            f"x={x_factor:.6g}, y={y_factor:.6g}, "
            f"quarter_holes={len(quarter_records)}"
        )


def _simulation_geometry(
    quarter_records: list[dict[str, object]],
    metadata: dict[str, object],
    finite_extent: float,
):
    layers = finite.simulation_layers(
        finite.polygons_from_records(quarter_records)
    )
    if metadata.get("finite_geometry") == "square":
        compiled = finite.compile_finite_geometry_input(
            metadata["geometry_plan"],
            layer_height=finite.SLAB_HEIGHT,
            refractive_index=finite.REFRACTIVE_INDEX,
            quarter=True,
            holes=finite.polygons_from_records(quarter_records),
        )
        return compiled.boundary, list(compiled.layers)
    return BoundarySpec("quarter_hexagon_boundary", a=finite_extent), layers


def _case_metadata(
    paths: dict[str, Path],
    geometry_metadata: dict[str, object],
    validity_boundary: np.ndarray,
    *,
    x_shift_factor: float,
    y_shift_factor: float,
    reuse_source: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        **finite.shared_case_metadata(geometry_metadata, validity_boundary),
        "case": "finite_quarter_all",
        "source": "four finite-quarter symmetry sectors in one COMSOL model",
        "symmetry_ids": list(SYMMETRY_IDS),
        "symmetry_records_internal": [
            {
                **quarter.symmetry_case(symmetry_id),
                **solution_identity(symmetry_id),
            }
            for symmetry_id in SYMMETRY_IDS
        ],
        "quarter_domain": {"x_min": 0.0, "y_min": 0.0},
        "field_output": "full field reconstructed by symmetry",
        "raw_quarter_exports_saved": False,
        "quarter_boundary": geometry_metadata["quarter_boundary"],
        "simulation_holes": geometry_metadata["quarter_hole_count"],
        "cladding_x_shift_factor": float(x_shift_factor),
        "cladding_y_shift_factor": float(y_shift_factor),
        "simulation": {
            "mode": "finite_quarter",
            "boundary": (
                "quarter_rectangle"
                if geometry_metadata.get("finite_geometry") == "square"
                else "quarter_hexagon_boundary"
            ),
            "eigenfrequency_shift": finite.EIGENFREQUENCY_SHIFT,
            "eigenmode_count": finite.EIGENMODE_COUNT,
            "mesh_auto_size": finite.MESH_AUTO_SIZE,
        },
        "execution": {
            "model_reuse": True,
            "model_build_count": 0 if reuse_source else 1,
            "mesh_build_count": 0 if reuse_source else 1,
            "reused_model_count": 1 if reuse_source else 0,
            "reused_mesh_count": 1 if reuse_source else 0,
            "planned_solve_count": 4,
            "planned_solution_copy_count": 4,
            "reuse_source": reuse_source,
        },
        "analysis": {
            "scope": "all_valid_modes",
            "symmetry_visible_in_figures": False,
        },
        "output_layout": finite.OUTPUT_LAYOUT_VERSION,
        "paths": {
            **finite.finite_case_output_metadata(paths),
            "mode_output_pattern": str(
                paths["results_dir"]
                / "mode{local_mode_idx}_{symmetry_id}_x{x}_y{y}"
            ),
        },
    }


def _prepare_partition(
    root_staging_dir: Path,
    symmetry_id: int,
) -> dict[str, Path]:
    partition_dir = _partition_dir(root_staging_dir, symmetry_id)
    return finite.prepare_finite_case_output(
        partition_dir,
        model_filename="internal_partition.mph",
        resume=False,
    )


def _assert_same_eigenfrequencies(
    working: np.ndarray,
    archived: np.ndarray,
    *,
    symmetry_id: int,
) -> None:
    if working.shape != archived.shape or not np.allclose(
        working,
        archived,
        rtol=1e-11,
        atol=1e-12,
        equal_nan=True,
    ):
        raise ValueError(
            "Copied COMSOL solution eigenfrequencies do not match the "
            f"working solution for symmetry ID {symmetry_id}"
        )


def copy_and_validate_solution(
    sim: SimulationRun,
    symmetry_id: int,
) -> dict[str, object]:
    identity = solution_identity(symmetry_id)
    working = np.asarray(sim.get_eigenfrequencies(), dtype=float)
    sim.copy_solution(
        identity["solution_tag"],
        identity["dataset_tag"],
        label=f"Internal symmetry solution {symmetry_id}",
    )
    archived = np.asarray(
        sim.get_eigenfrequencies(identity["dataset_tag"]),
        dtype=float,
    )
    _assert_same_eigenfrequencies(
        working,
        archived,
        symmetry_id=symmetry_id,
    )
    return {
        **identity,
        "eigenfrequency_row_count": int(len(working)),
        "eigenfrequency_verified": True,
        "field_probe_verified": False,
    }


def validate_solution_field_probe(
    sim: SimulationRun,
    eigen_df: pd.DataFrame,
    archive: dict[str, object],
) -> None:
    valid = eigen_df.loc[eigen_df["is_valid"].map(finite._csv_bool)]
    if eigen_df.empty:
        raise ValueError("Cannot validate solution copy: no mode exists")
    mode_idx = int((valid if not valid.empty else eigen_df).iloc[0]["mode_idx"])
    coordinates, _values = sim.get_2d_fields(
        mode_idx,
        "ewfd.Hz",
        "center",
    )
    point = np.asarray(coordinates, dtype=float)[0, :2]
    probe = np.array([[point[0], point[1], 0.0]], dtype=float)
    working = sim.get_fields_at_coordinates(
        mode_idx,
        ["ewfd.Hz"],
        probe,
        dataset="dset1",
    )
    archived = sim.get_fields_at_coordinates(
        mode_idx,
        ["ewfd.Hz"],
        probe,
        dataset=str(archive["dataset_tag"]),
    )
    if not np.allclose(
        working,
        archived,
        rtol=1e-9,
        atol=1e-12,
        equal_nan=True,
    ):
        raise ValueError(
            "Copied COMSOL solution field probe does not match the working "
            f"solution for mode {mode_idx}"
        )
    archive["field_probe_mode_idx"] = mode_idx
    archive["field_probe_mode_valid"] = bool(not valid.empty)
    archive["field_probe_verified"] = True


def _annotate_eigenfrequencies(
    eigen_df: pd.DataFrame,
    case: dict[str, str | int],
    archive: dict[str, object],
) -> pd.DataFrame:
    frame = eigen_df.copy()
    frame["local_mode_idx"] = frame["mode_idx"].astype(int)
    frame["mode_uid"] = [
        mode_uid(mode_idx, int(case["id"]))
        for mode_idx in frame["mode_idx"].astype(int)
    ]
    frame["solution_tag"] = str(archive["solution_tag"])
    frame["dataset_tag"] = str(archive["dataset_tag"])
    return frame


def _run_partition_postprocess(
    paths: dict[str, Path],
    geometry_metadata: dict[str, object],
) -> None:
    eigen_df = pd.read_csv(paths["export_dir"] / "eigenfrequencies.csv")
    has_valid_modes = bool(eigen_df["is_valid"].map(finite._csv_bool).any())
    if finite.RUN_FARFIELD_FFT:
        finite.run_farfield_fft_case(
            paths["staging_dir"],
            paths["export_dir"],
            grid_size=finite.FARFIELD_GRID_SIZE,
            fft_size=finite.FARFIELD_FFT_SIZE,
            na=finite.FARFIELD_NA,
            h0_um=finite.FARFIELD_H0_UM,
            period_um=finite.A,
            mode_indices=None,
        )
    if finite.RUN_FINITE_LATTICE_FOURIER_POSTPROCESS and has_valid_modes:
        finite.run_exported_finite_lattice_fourier_postprocess(
            paths["staging_dir"],
            paths["export_dir"],
            cells=list(geometry_metadata["bulk_points"]),
            cladding_cells=list(geometry_metadata["cladding_points"]),
            bulk_radius=finite.BULK_RADIUS,
            period=finite.A,
            rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            top_k=finite.FINITE_LATTICE_FOURIER_TOP_K,
            mode_indices=None,
            include_unselected_valid_mode_maps=True,
            dpi=finite.DPI,
            transform_kind=(
                "arbitrary_square_cavity_v1"
                if geometry_metadata.get("finite_geometry") == "square"
                else "hex_cyclic_quotient_v1"
            ),
        )
        finite.score_exported_finite_modes(
            paths["staging_dir"],
            export_dirname=finite.EXPORT_DIRNAME,
        )
    elif has_valid_modes:
        finite.run_exported_finite_unit_cell_intensity_maps(
            paths["staging_dir"],
            paths["export_dir"],
            cells=list(geometry_metadata["bulk_points"]),
            cladding_cells=list(geometry_metadata["cladding_points"]),
            period=finite.A,
            rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            mode_indices=None,
            dpi=finite.DPI,
        )
    finite.finalize_finite_case_output(
        paths["case_dir"],
        model_filename="internal_partition.mph",
    )


def _replace_paths(value, replacements: dict[str, str]):
    if isinstance(value, str):
        for source, destination in replacements.items():
            value = value.replace(source, destination)
        return value
    if isinstance(value, list):
        return [_replace_paths(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            key: _replace_paths(item, replacements)
            for key, item in value.items()
        }
    return value


def _rewrite_mode_tree_paths(
    mode_dir: Path,
    replacements: dict[str, str],
) -> None:
    for path in mode_dir.rglob("*.json"):
        payload = json.loads(path.read_text(encoding="utf-8"))
        finite.write_config(path, _replace_paths(payload, replacements))
    for path in mode_dir.rglob("*.csv"):
        text = path.read_text(encoding="utf-8")
        rewritten = _replace_paths(text, replacements)
        if rewritten != text:
            path.write_text(rewritten, encoding="utf-8")


def _annotate_aggregate_frame(
    frame: pd.DataFrame,
    *,
    symmetry_id: int,
    replacements: dict[str, str],
) -> pd.DataFrame:
    frame = frame.copy()
    case = quarter.symmetry_case(symmetry_id)
    if "mode_idx" in frame.columns:
        frame["local_mode_idx"] = frame["mode_idx"].astype(int)
        frame["mode_uid"] = [
            mode_uid(mode_idx, symmetry_id)
            for mode_idx in frame["mode_idx"].astype(int)
        ]
    frame["symmetry_id"] = symmetry_id
    frame["x_boundary"] = str(case["x_boundary"])
    frame["y_boundary"] = str(case["y_boundary"])
    frame["solution_tag"] = solution_identity(symmetry_id)["solution_tag"]
    frame["dataset_tag"] = solution_identity(symmetry_id)["dataset_tag"]
    for column in frame.select_dtypes(include=["object", "string"]).columns:
        frame[column] = frame[column].map(
            lambda value: _replace_paths(value, replacements)
        )
    return frame


def _remove_empty_tree(root: Path) -> None:
    for directory in sorted(
        (path for path in root.rglob("*") if path.is_dir()),
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        directory.rmdir()
    root.rmdir()


def _merge_partition_outputs(
    paths: dict[str, Path],
) -> dict[str, object]:
    """Merge four finalized internal cases into one mode-centric dataset."""
    aggregate_frames: dict[str, list[pd.DataFrame]] = {
        name: [] for name in AGGREGATE_TABLES
    }
    partition_plans = []
    all_mode_uids: list[str] = []
    reference_indices: set[int] | None = None
    per_symmetry_mode_counts: dict[str, int] = {}
    shared_sources: dict[str, list[Path]] = {
        name: [] for name in SHARED_OVERVIEW_FILES
    }

    for symmetry_id in SYMMETRY_IDS:
        partition_dir = _partition_dir(paths["staging_dir"], symmetry_id)
        partition = finite.finite_case_output_paths(
            partition_dir,
            model_filename="internal_partition.mph",
        )
        eigen_path = partition["overview_dir"] / "eigenfrequencies.csv"
        if not eigen_path.is_file():
            raise FileNotFoundError(
                f"Missing finalized symmetry eigenfrequency table: {eigen_path}"
            )
        eigen_frame = pd.read_csv(eigen_path)
        if "mode_idx" not in eigen_frame.columns:
            raise ValueError(f"Missing mode_idx column in {eigen_path}")
        numeric_indices = pd.to_numeric(
            eigen_frame["mode_idx"], errors="raise"
        ).to_numpy(float)
        actual_indices = set(numeric_indices.astype(int).tolist())
        expected_indices = set(range(len(eigen_frame)))
        if (
            len(eigen_frame) < finite.EIGENMODE_COUNT
            or not np.array_equal(numeric_indices, numeric_indices.astype(int))
            or actual_indices != expected_indices
        ):
            raise ValueError(
                f"Symmetry ID {symmetry_id} must contain at least "
                f"{finite.EIGENMODE_COUNT} unique local modes numbered from zero"
            )
        if reference_indices is None:
            reference_indices = actual_indices
        elif actual_indices != reference_indices:
            raise ValueError("All symmetry sectors must return the same mode indices")
        per_symmetry_mode_counts[str(symmetry_id)] = len(eigen_frame)

        mode_moves = []
        for mode_idx in sorted(expected_indices):
            source = partition["results_dir"] / f"mode{mode_idx}"
            destination = paths["results_dir"] / mode_uid(
                mode_idx,
                symmetry_id,
            )
            if not source.is_dir():
                raise FileNotFoundError(f"Missing finalized mode directory: {source}")
            if destination.exists():
                raise FileExistsError(
                    f"Refusing to overwrite merged mode directory: {destination}"
                )
            mode_moves.append((source, destination))
            all_mode_uids.append(destination.name)

        replacements = {
            str(source): str(destination)
            for source, destination in mode_moves
        }
        replacements[str(partition_dir)] = str(paths["case_dir"])
        for name in AGGREGATE_TABLES:
            source = partition["overview_dir"] / name
            if source.is_file():
                try:
                    frame = pd.read_csv(source)
                except pd.errors.EmptyDataError:
                    frame = pd.DataFrame()
                aggregate_frames[name].append(
                    _annotate_aggregate_frame(
                        frame,
                        symmetry_id=symmetry_id,
                        replacements=replacements,
                    )
                )
        for name in SHARED_OVERVIEW_FILES:
            source = partition["overview_dir"] / name
            if source.is_file():
                shared_sources[name].append(source)

        config_destination = (
            paths["config_dir"] / "internal_symmetry" / str(symmetry_id)
        )
        if config_destination.exists():
            raise FileExistsError(
                f"Refusing to overwrite internal symmetry config: {config_destination}"
            )

        recognized_roots = [
            partition["results_dir"],
            partition["config_dir"],
        ]
        recognized_files = {
            partition["overview_dir"] / name for name in AGGREGATE_TABLES
        }
        recognized_files.update(
            partition["overview_dir"] / name
            for name in SHARED_OVERVIEW_FILES
        )
        unknown = [
            source
            for source in partition_dir.rglob("*")
            if source.is_file()
            and source not in recognized_files
            and not any(source.is_relative_to(root) for root in recognized_roots)
        ]
        if unknown:
            raise ValueError(f"Unknown finalized symmetry output: {unknown[0]}")
        partition_plans.append(
            (
                partition,
                mode_moves,
                replacements,
                config_destination,
            )
        )

    expected_total = sum(per_symmetry_mode_counts.values())
    if (
        len(all_mode_uids) != expected_total
        or len(set(all_mode_uids)) != expected_total
    ):
        raise ValueError(
            f"Merged mode identity check failed: expected {expected_total} unique modes"
        )
    for name, sources in shared_sources.items():
        if len(sources) > 1:
            reference = sources[0].read_bytes()
            if any(source.read_bytes() != reference for source in sources[1:]):
                raise ValueError(
                    f"Internal symmetry copies of {name} are not identical"
                )

    for partition, mode_moves, replacements, config_destination in partition_plans:
        for source, destination in mode_moves:
            _rewrite_mode_tree_paths(source, replacements)
            source.replace(destination)
        if partition["config_dir"].is_dir():
            _rewrite_mode_tree_paths(partition["config_dir"], replacements)
            config_destination.parent.mkdir(parents=True, exist_ok=True)
            partition["config_dir"].replace(config_destination)
        for name in AGGREGATE_TABLES:
            source = partition["overview_dir"] / name
            if source.is_file():
                source.unlink()
        for name in SHARED_OVERVIEW_FILES:
            source = partition["overview_dir"] / name
            if not source.is_file():
                continue
            destination = paths["overview_dir"] / name
            if destination.exists():
                source.unlink()
            else:
                source.replace(destination)
        remaining_files = [
            source for source in partition["case_dir"].rglob("*") if source.is_file()
        ]
        if remaining_files:
            raise RuntimeError(
                f"Symmetry partition cleanup found an unhandled file: {remaining_files[0]}"
            )
        _remove_empty_tree(partition["case_dir"])

    merged_counts = {}
    for name, frames in aggregate_frames.items():
        if not frames:
            continue
        merged = pd.concat(frames, ignore_index=True, sort=False)
        if name == "eigenfrequencies.csv":
            if len(merged) != expected_total:
                raise ValueError(
                    f"Merged eigenfrequency table has {len(merged)} rows; "
                    f"expected {expected_total}"
                )
            if merged["mode_uid"].nunique() != expected_total:
                raise ValueError("Merged eigenfrequency mode_uid values are not unique")
        merged.to_csv(paths["overview_dir"] / name, index=False)
        merged_counts[name] = int(len(merged))

    remaining_files = [
        source for source in paths["staging_dir"].rglob("*") if source.is_file()
    ]
    if remaining_files:
        raise RuntimeError(
            f"All-symmetry staging cleanup found an unhandled file: {remaining_files[0]}"
        )
    _remove_empty_tree(paths["staging_dir"])
    return {
        "mode_count": expected_total,
        "requested_eigenmode_count_per_symmetry": finite.EIGENMODE_COUNT,
        "actual_mode_count_per_symmetry": per_symmetry_mode_counts,
        "table_row_counts": merged_counts,
    }


def save_frequency_q_plot(
    eigenfrequency_csv: Path,
    output_path: Path,
) -> dict[str, object]:
    frame = pd.read_csv(eigenfrequency_csv)
    frequency = pd.to_numeric(frame["re"], errors="coerce").to_numpy(float)
    quality_factor = pd.to_numeric(frame["q"], errors="coerce").to_numpy(float)
    finite_mask = np.isfinite(frequency) & np.isfinite(quality_factor)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(5.2, 4.0), constrained_layout=True)
    axis.scatter(
        frequency[finite_mask],
        quality_factor[finite_mask],
        s=22,
        marker="o",
        color="#1f77b4",
        edgecolors="none",
    )
    axis.set_xscale("linear")
    axis.set_yscale("linear")
    axis.set_xlabel("Frequency (THz)")
    axis.set_ylabel("Q")
    axis.tick_params(direction="in", top=True, right=True)
    figure.savefig(output_path, dpi=finite.DPI)
    plt.close(figure)
    return {
        "source_row_count": int(len(frame)),
        "plotted_row_count": int(finite_mask.sum()),
        "nonfinite_row_count": int((~finite_mask).sum()),
        "x_scale": "linear",
        "y_scale": "linear",
        "symmetry_visible": False,
        "output_path": str(output_path),
    }


def run_all_symmetry_shift(
    x_shift_factor: float,
    y_shift_factor: float,
) -> dict[str, object]:
    finite.set_cladding_inward_shift_factors(
        x_shift_factor,
        y_shift_factor,
    )
    (
        full_records,
        quarter_records,
        geometry_metadata,
        _quarter_boundary,
        finite_extent,
    ) = quarter.build_quarter_case_geometry()
    validate_cavity_edge_cell_count(geometry_metadata)
    reuse_source = (
        reuse_model_source_info(
            REUSE_MESHED_MODEL,
            x_shift_factor,
            y_shift_factor,
        )
        if REUSE_MESHED_MODEL is not None
        else None
    )
    if reuse_source is not None and reuse_source["quarter_hole_count"] != int(
        geometry_metadata["quarter_hole_count"]
    ):
        raise ValueError(
            "Reused finite-quarter hole count mismatch: "
            f"source={reuse_source['quarter_hole_count']}, "
            f"requested={geometry_metadata['quarter_hole_count']}"
        )
    validity_boundary = finite.strip_validity_boundary(geometry_metadata)
    case_dir = case_dir_for_shift(x_shift_factor, y_shift_factor)
    if RESUME_INCOMPLETE:
        paths = finite.finite_case_output_paths(
            case_dir,
            model_filename=MODEL_FILENAME,
        )
        metadata, summary, partitions = _load_resume_state(
            paths,
            x_shift_factor,
            y_shift_factor,
        )
        completed_ids = [
            int(run["id"]) for run in summary["symmetry_runs_internal"]
        ]
        summary.pop("failure", None)
        summary["status"] = "running"
        summary["resume_count"] = int(summary.get("resume_count", 0)) + 1
        summary["resumed_completed_symmetry_ids"] = completed_ids
    else:
        paths = finite.prepare_finite_case_output(
            case_dir,
            model_filename=MODEL_FILENAME,
            resume=False,
        )
        metadata = _case_metadata(
            paths,
            geometry_metadata,
            validity_boundary,
            x_shift_factor=x_shift_factor,
            y_shift_factor=y_shift_factor,
            reuse_source=reuse_source,
        )
        finite.write_config(paths["config"], metadata)
        for message in finite.shift_profile_preflight_messages(metadata):
            print(message)
        finite.save_finite_shift_profile_diagnostics(
            paths["overview_dir"], metadata, dpi=finite.DPI
        )
        finite.save_finite_model_geometry_figures(
            paths["model_dir"],
            full_records,
            finite.full_cell_records(geometry_metadata),
            np.asarray(geometry_metadata["finite_boundary"], dtype=float),
            validity_boundary,
            {**metadata, "simulation_holes": len(full_records)},
        )
        summary = {
            "status": "running",
            "case": "finite_quarter_all",
            "finite_geometry": finite.FINITE_GEOMETRY,
            "cladding_x_shift_factor": float(x_shift_factor),
            "cladding_y_shift_factor": float(y_shift_factor),
            "symmetry_ids": list(SYMMETRY_IDS),
            "model_build_count": 0,
            "mesh_build_count": 0,
            "reused_model_count": 0,
            "reused_mesh_count": 0,
            "solve_count": 0,
            "solution_copy_count": 0,
            "symmetry_runs_internal": [],
            "output_layout": finite.OUTPUT_LAYOUT_VERSION,
            **finite.finite_case_output_metadata(paths),
        }
        partitions = []
        completed_ids = []
    summary_path = paths["config_dir"] / SUMMARY_FILENAME
    finite.write_config(summary_path, summary)

    if not RUN_COMSOL:
        summary["status"] = "geometry_only"
        finite.write_config(summary_path, summary)
        return summary

    if finite.RUN_FARFIELD_FFT:
        finite.validate_farfield_config(
            finite.FARFIELD_GRID_SIZE,
            finite.FARFIELD_FFT_SIZE,
            finite.FARFIELD_NA,
            finite.FARFIELD_H0_UM,
        )

    sim_boundary, layers = _simulation_geometry(
        quarter_records,
        geometry_metadata,
        finite_extent,
    )
    try:
        first_case = quarter.symmetry_case(SYMMETRY_IDS[0])
        with SimulationRun(
            config=quarter.quarter_simulation_config(first_case),
            progress_log_path=paths["progress_log"],
        ) as sim:
            if RESUME_INCOMPLETE:
                finite.attach_saved_model(sim, paths["mph"])
                summary["resume_prepare"] = (
                    sim.prepare_reused_finite_quarter_model()
                )
                jmodel = sim.model.java
                for run, partition in zip(
                    summary["symmetry_runs_internal"], partitions
                ):
                    symmetry_id = int(run["id"])
                    identity = solution_identity(symmetry_id)
                    if (
                        not jmodel.sol().hasTag(identity["solution_tag"])
                        or not jmodel.result().dataset().hasTag(
                            identity["dataset_tag"]
                        )
                    ):
                        raise ValueError(
                            f"Saved MPH is missing symmetry {symmetry_id} solution copy"
                        )
                    staged = pd.read_csv(
                        partition["export_dir"] / "eigenfrequencies.csv"
                    )[["re", "im", "q"]].to_numpy(float)
                    archived = np.asarray(
                        sim.get_eigenfrequencies(identity["dataset_tag"]),
                        dtype=float,
                    )
                    _assert_same_eigenfrequencies(
                        staged,
                        archived,
                        symmetry_id=symmetry_id,
                    )
                summary["resume_solution_validation"] = {
                    "completed_symmetry_ids": completed_ids,
                    "eigenfrequencies_verified": True,
                }
                working_solution_id = None
            elif reuse_source is None:
                sim.build_geometry(
                    sim_boundary,
                    layers,
                    simulation_mode="finite_quarter",
                    k=None,
                )
                summary["model_build_count"] = 1
            else:
                finite.attach_saved_model(sim, Path(reuse_source["mph_path"]))
                summary["reuse_prepare"] = (
                    sim.prepare_reused_finite_quarter_model()
                )
                summary["reused_model_count"] = 1
                summary["reused_mesh_count"] = 1
            if not RESUME_INCOMPLETE:
                sim.set_finite_quarter_symmetry(
                    str(first_case["x_boundary"]),
                    str(first_case["y_boundary"]),
                )
                if reuse_source is None:
                    sim.model.save(str(paths["mph"]))
                    sim.run_simulation(
                        mesh_save_path=paths["mph"],
                        mesh_builder=finite.finite_mesh_builder(layers),
                    )
                    summary["mesh_build_count"] = 1
                else:
                    sim.rerun_eigenfrequency_solution()
                summary["solve_count"] = 1
                working_solution_id = SYMMETRY_IDS[0]
            finite.write_config(summary_path, summary)

            pending_ids = SYMMETRY_IDS[len(completed_ids):]
            for symmetry_id in pending_ids:
                case = quarter.symmetry_case(symmetry_id)
                if working_solution_id != symmetry_id:
                    sim.set_finite_quarter_symmetry(
                        str(case["x_boundary"]),
                        str(case["y_boundary"]),
                    )
                    sim.rerun_eigenfrequency_solution()
                    summary["solve_count"] += 1
                    finite.write_config(summary_path, summary)

                archive = copy_and_validate_solution(sim, symmetry_id)
                sim.model.save(str(paths["mph"]))
                summary["solution_copy_count"] += 1
                partition = _prepare_partition(paths["staging_dir"], symmetry_id)
                partitions.append(partition)
                eigen_df = quarter.export_reconstructed_results(
                    sim,
                    partition["export_dir"],
                    validity_boundary,
                    finite_hex_a=finite_extent,
                    case=case,
                    finite_geometry=str(
                        geometry_metadata.get("finite_geometry", "hex")
                    ),
                    footprint_width=float(geometry_metadata["finite_Lx"]),
                    footprint_height=float(geometry_metadata["finite_Ly"]),
                )
                annotated = _annotate_eigenfrequencies(
                    eigen_df,
                    case,
                    archive,
                )
                annotated.to_csv(
                    partition["export_dir"] / "eigenfrequencies.csv",
                    index=False,
                )
                validate_solution_field_probe(sim, annotated, archive)
                finite.write_config(
                    partition["config"],
                    {
                        "case": "finite_quarter_all_internal_symmetry",
                        "symmetry_internal": case,
                        "solution_archive": archive,
                        "eigenmode_count": finite.EIGENMODE_COUNT,
                        "paths": finite.finite_case_output_metadata(partition),
                    },
                )
                sim.model.save(str(paths["mph"]))
                summary["symmetry_runs_internal"].append(
                    {
                        **case,
                        **archive,
                        "mode_count": int(len(annotated)),
                        "status": "solution_copied_exported_verified",
                    }
                )
                finite.write_config(summary_path, summary)

        for partition in partitions:
            _run_partition_postprocess(partition, geometry_metadata)
        merge_summary = _merge_partition_outputs(paths)
        plot_summary = save_frequency_q_plot(
            paths["overview_dir"] / "eigenfrequencies.csv",
            paths["overview_dir"] / "f_Q.png",
        )
        summary.update(
            {
                "status": "complete",
                "merge": merge_summary,
                "frequency_q_plot": plot_summary,
            }
        )
        actual_counts = (
            summary["model_build_count"],
            summary["mesh_build_count"],
            summary["reused_model_count"],
            summary["reused_mesh_count"],
            summary["solve_count"],
            summary["solution_copy_count"],
        )
        expected_counts = (
            (1, 1, 0, 0, 4, 4)
            if reuse_source is None and not RESUME_INCOMPLETE
            else (0, 0, 1, 1, 4, 4)
        )
        if actual_counts != expected_counts:
            raise RuntimeError("All-symmetry execution counts are incomplete")
        finite.write_config(summary_path, summary)
        return summary
    except Exception as exc:
        summary["status"] = "failed"
        summary["failure"] = {
            "type": type(exc).__name__,
            "message": str(exc),
        }
        finite.write_config(summary_path, summary)
        raise


def run_all_symmetries() -> list[dict[str, object]]:
    validate_runtime_config()
    validate_output_targets()
    return [
        run_all_symmetry_shift(x_shift_factor, y_shift_factor)
        for x_shift_factor, y_shift_factor in (
            finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
        )
    ]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Solve all four finite-quarter symmetry sectors in one COMSOL model"
        )
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="validate parameters, geometry, and output targets without solving",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="override the finite_quarter_all series directory",
    )
    parser.add_argument(
        "--reuse-meshed-model",
        type=Path,
        help="load an existing solved finite-quarter MPH and reuse its mesh",
    )
    parser.add_argument(
        "--resume-incomplete",
        action="store_true",
        help="resume verified symmetry-prefix results from the target output",
    )
    parser.add_argument(
        "--allow-concurrent-comsol",
        action="store_true",
        help="allow startup while another COMSOL process is already active",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    global OUT_DIR, REUSE_MESHED_MODEL, RESUME_INCOMPLETE
    args = parse_args(argv)
    if args.output_dir is not None:
        resolved_output = args.output_dir.resolve()
        finite_root = finite.FINITE_CAVITY_OUTPUT_ROOT.resolve()
        if (
            resolved_output != finite_root
            and finite_root not in resolved_output.parents
        ):
            raise ValueError(
                "--output-dir must be located under "
                f"{finite.FINITE_CAVITY_OUTPUT_ROOT}"
            )
        OUT_DIR = args.output_dir
    REUSE_MESHED_MODEL = args.reuse_meshed_model
    RESUME_INCOMPLETE = args.resume_incomplete
    validate_runtime_config()
    validate_output_targets()
    print_preflight()
    if args.preflight_only:
        preflight_geometry()
        print("Preflight complete; COMSOL was not started.")
        return
    active_processes = _detect_comsol_processes()
    if active_processes and not args.allow_concurrent_comsol:
        raise RuntimeError(
            "Refusing to start while another COMSOL/mphserver process is active; "
            "close it or explicitly pass --allow-concurrent-comsol"
        )
    summaries = run_all_symmetries()
    print(f"Completed {len(summaries)} all-symmetry finite-quarter case(s).")
    for summary in summaries:
        print(f"  {summary['case_dir']}")


if __name__ == "__main__":
    main()
