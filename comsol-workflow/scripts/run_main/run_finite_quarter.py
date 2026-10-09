#!/usr/bin/env python3
"""Solve a finite cavity on its first quadrant and restore the full fields."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_finite as finite  # noqa: E402

from comsol_workflow.field_plotting import save_field_plot  # noqa: E402
from comsol_workflow.polygon_utils import (  # noqa: E402
    clip_polygon_to_convex_region,
    polygon_area,
)
from comsol_workflow.simulation_utils import (  # noqa: E402
    BoundarySpec,
    SimulationConfig,
    SimulationRun,
    hz_reflection_sign,
    reconstruct_quarter_electric_grid,
    reconstruct_quarter_scalar_field,
)


OUTPUT_ROOT = finite.OUTPUT_ROOT
OUT_DIR = OUTPUT_ROOT / finite.ACTIVE_PARAMETERS.structure_series_label(
    "finite_quarter"
)
SYMMETRY_IDS = list(finite.ACTIVE_PARAMETERS.finite_quarter_symmetry_ids)
SYMMETRY_CASES: dict[int, dict[str, str]] = {
    1: {"x_boundary": "PEC", "y_boundary": "PMC"},
    2: {"x_boundary": "PMC", "y_boundary": "PEC"},
    3: {"x_boundary": "PEC", "y_boundary": "PEC"},
    4: {"x_boundary": "PMC", "y_boundary": "PMC"},
}

RUN_COMSOL = True
RESUME_FROM_EXISTING_MPH = False
DPI = finite.DPI


def validate_symmetry_ids(symmetry_ids: list[int]) -> list[int]:
    if not isinstance(symmetry_ids, list) or not symmetry_ids:
        raise ValueError("SYMMETRY_IDS must be a non-empty list")
    normalized = []
    for symmetry_id in symmetry_ids:
        if type(symmetry_id) is not int or symmetry_id not in SYMMETRY_CASES:
            raise ValueError("SYMMETRY_IDS entries must be integers from 1 to 4")
        if symmetry_id in normalized:
            raise ValueError(f"SYMMETRY_IDS contains duplicate ID {symmetry_id}")
        normalized.append(symmetry_id)
    return normalized


def symmetry_case(symmetry_id: int) -> dict[str, str | int]:
    if symmetry_id not in SYMMETRY_CASES:
        raise ValueError(f"Unknown symmetry ID: {symmetry_id}")
    return {"id": symmetry_id, **SYMMETRY_CASES[symmetry_id]}


def symmetry_dirname(symmetry_id: int) -> str:
    case = symmetry_case(symmetry_id)
    return (
        f"symmetry_{symmetry_id}"
        f"_x{case['x_boundary']}"
        f"_y{case['y_boundary']}"
    )


def resolved_target_frequency_thz(target_frequency_thz: float | None) -> float:
    value = (
        finite.FINITE_SINGLE_MODE_TARGET_FREQUENCY_THZ
        if target_frequency_thz is None
        else float(target_frequency_thz)
    )
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("target_frequency_thz must be positive and finite")
    return value


def target_eigenfrequency_shift(target_frequency_thz: float | None) -> str:
    if target_frequency_thz is None:
        return finite.EIGENFREQUENCY_SHIFT
    value = resolved_target_frequency_thz(target_frequency_thz)
    label = f"{value:.12f}".rstrip("0").rstrip(".")
    return f"{label} [THz]"


def quarter_case_stem(
    x_shift_factor: float,
    y_shift_factor: float | None = None,
) -> str:
    return finite.case_stem(x_shift_factor, y_shift_factor)


def quarter_output_case_dir(
    x_shift_factor: float,
    symmetry_id: int,
    symmetry_ids: list[int] | None = None,
    *,
    y_shift_factor: float | None = None,
) -> Path:
    normalized = validate_symmetry_ids(
        list(SYMMETRY_IDS) if symmetry_ids is None else list(symmetry_ids)
    )
    shift_dir = OUT_DIR / quarter_case_stem(
        x_shift_factor,
        y_shift_factor,
    )
    if normalized == [1] and symmetry_id == 1:
        return shift_dir
    return shift_dir / symmetry_dirname(symmetry_id)


def quarter_boundary_vertices(finite_hex_a: float) -> np.ndarray:
    root3 = float(np.sqrt(3.0))
    return np.array(
        [
            [0.0, 0.0],
            [finite_hex_a / 2.0, 0.0],
            [finite_hex_a / 4.0, root3 * finite_hex_a / 4.0],
            [0.0, root3 * finite_hex_a / 4.0],
        ],
        dtype=float,
    )


def quarter_rectangle_vertices(width: float, height: float) -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0],
            [width / 2.0, 0.0],
            [width / 2.0, height / 2.0],
            [0.0, height / 2.0],
        ],
        dtype=float,
    )


def _canonical_polygon(polygon: np.ndarray, decimals: int = 9) -> tuple[tuple[float, float], ...]:
    rounded = np.round(np.asarray(polygon, dtype=float), decimals=decimals)
    return tuple(sorted((float(x), float(y)) for x, y in rounded))


def validate_full_geometry_mirror_symmetry(
    records: list[dict[str, object]],
    *,
    decimals: int = 9,
) -> None:
    polygon_keys = {
        _canonical_polygon(np.asarray(record["polygon"], dtype=float), decimals)
        for record in records
    }
    reflections = (
        np.array([1.0, -1.0]),
        np.array([-1.0, 1.0]),
    )
    for record_idx, record in enumerate(records):
        polygon = np.asarray(record["polygon"], dtype=float)
        for reflection in reflections:
            reflected_key = _canonical_polygon(polygon * reflection, decimals)
            if reflected_key not in polygon_keys:
                axis = "x" if reflection[1] < 0.0 else "y"
                raise ValueError(
                    "Quarter simulation requires mirror-symmetric full geometry; "
                    f"record {record_idx} has no reflection across the {axis} axis"
                )


def clip_records_to_first_quadrant(
    records: list[dict[str, object]],
    boundary: np.ndarray,
) -> list[dict[str, object]]:
    clipped_records = []
    for record in records:
        polygon = clip_polygon_to_convex_region(
            np.asarray(record["polygon"], dtype=float),
            boundary,
        )
        if len(polygon) < 3 or polygon_area(polygon) <= 1e-12:
            continue
        clipped_records.append({**record, "polygon": polygon})
    if not clipped_records:
        raise ValueError("Quarter clipping removed every finite-cavity hole")
    return clipped_records


def quarter_simulation_config(
    case: dict[str, str | int],
    *,
    target_frequency_thz: float | None = None,
) -> SimulationConfig:
    return SimulationConfig(
        wavelength=finite.WAVELENGTH,
        eigenfrequency_shift=target_eigenfrequency_shift(
            target_frequency_thz
        ),
        eigenmode_count=finite.EIGENMODE_COUNT,
        mesh_auto_size=finite.MESH_AUTO_SIZE,
        mode_type=finite.MODE_TYPE,
        air_cutplane_z="H_air",
        symmetry_x_boundary=str(case["x_boundary"]),
        symmetry_y_boundary=str(case["y_boundary"]),
    )


def _reconstruct_hz_field(
    sim_run,
    mode_idx: int,
    case: dict[str, str | int],
) -> tuple[np.ndarray, np.ndarray]:
    from comsol_workflow.energy_recovery import interpolate_field

    coordinates_re, hz_re = sim_run.get_2d_fields(mode_idx, "ewfd.Hz", "center")
    coordinates_im, hz_im = sim_run.get_2d_fields(
        mode_idx,
        "ewfd.Hz*(-i)",
        "center",
    )
    coordinates_re = np.asarray(coordinates_re, dtype=float)[:, :2]
    coordinates_im = np.asarray(coordinates_im, dtype=float)[:, :2]
    hz_im_aligned = interpolate_field(coordinates_im, hz_im, coordinates_re)
    quarter_hz = np.real(hz_re) + 1j * np.real(hz_im_aligned)
    return reconstruct_quarter_scalar_field(
        coordinates_re,
        quarter_hz,
        x_reflection_sign=hz_reflection_sign(str(case["x_boundary"])),
        y_reflection_sign=hz_reflection_sign(str(case["y_boundary"])),
    )


def _reconstruct_wem_field(sim_run, mode_idx: int) -> tuple[np.ndarray, np.ndarray]:
    coordinates, values = sim_run.get_2d_fields(mode_idx, "ewfd.Wav", "center")
    return reconstruct_quarter_scalar_field(
        np.asarray(coordinates, dtype=float)[:, :2],
        np.real(values),
        x_reflection_sign=1.0,
        y_reflection_sign=1.0,
    )


def export_reconstructed_air_fields(
    sim_run,
    export_dir: Path,
    eigen_df: pd.DataFrame,
    *,
    finite_hex_a: float,
    grid_size: int,
    h0_um: float | None,
    case: dict[str, str | int],
    finite_geometry: str = "hex",
    footprint_width: float | None = None,
    footprint_height: float | None = None,
    mode_indices: list[int] | None = None,
) -> dict[str, object]:
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("far-field grid_size must be an odd integer >= 3")
    export_dir.mkdir(parents=True, exist_ok=True)
    air_plane_z_um = float(
        sim_run.model.java.param().evaluate("z_slab_top + H_air", "um")
    )
    window_size = (
        float(finite_hex_a)
        if finite_geometry == "hex"
        else max(float(footprint_width), float(footprint_height))
    )
    half_width = window_size / 2.0
    full_axis = np.linspace(-half_width, half_width, grid_size)
    midpoint = grid_size // 2
    quarter_axis = full_axis[midpoint:]
    xx_quarter, yy_quarter = np.meshgrid(quarter_axis, quarter_axis, indexing="xy")
    quarter_coordinates = np.column_stack(
        [
            xx_quarter.ravel(),
            yy_quarter.ravel(),
            np.full(xx_quarter.size, air_plane_z_um),
        ]
    )
    tolerance = max(window_size, 1.0) * 1e-12
    if finite_geometry == "square":
        geometric_quarter = (
            xx_quarter <= float(footprint_width) / 2.0 + tolerance
        ) & (
            yy_quarter <= float(footprint_height) / 2.0 + tolerance
        )
    else:
        root3 = float(np.sqrt(3.0))
        geometric_quarter = (
            yy_quarter <= root3 * float(finite_hex_a) / 4.0 + tolerance
        ) & (
            root3 * xx_quarter + yy_quarter
            <= root3 * float(finite_hex_a) / 2.0 + tolerance
        )

    valid_rows = eigen_df.loc[eigen_df["is_valid"].map(finite._csv_bool)]
    valid_mode_indices = valid_rows["mode_idx"].to_numpy(dtype=int).tolist()
    if mode_indices is None:
        exported_mode_indices = valid_mode_indices
    else:
        exported_mode_indices = [int(mode_idx) for mode_idx in mode_indices]
        if len(exported_mode_indices) != len(set(exported_mode_indices)):
            raise ValueError("Air-field mode indices must be unique")
        invalid = sorted(
            set(exported_mode_indices).difference(valid_mode_indices)
        )
        if invalid:
            raise ValueError(
                "Air-field modes must be valid exported modes; "
                f"invalid modes: {invalid}"
            )
    for mode_idx in exported_mode_indices:
        values = np.asarray(
            sim_run.get_fields_at_coordinates(
                mode_idx,
                ["ewfd.Ex", "ewfd.Ey", "ewfd.Ez"],
                quarter_coordinates,
                dataset="dset1",
            ),
            dtype=complex,
        )
        expected_shape = (3, len(quarter_coordinates))
        if values.shape != expected_shape:
            raise ValueError(
                f"Expected quarter air electric-field shape {expected_shape}, "
                f"got {values.shape}"
            )
        finite_quarter = np.all(
            np.isfinite(np.real(values)) & np.isfinite(np.imag(values)),
            axis=0,
        ).reshape(geometric_quarter.shape)
        quarter_mask = geometric_quarter & finite_quarter
        quarter_values = values.reshape(3, *geometric_quarter.shape)
        quarter_values[:, ~quarter_mask] = 0.0
        full_values, full_mask = reconstruct_quarter_electric_grid(
            quarter_values,
            quarter_mask,
            x_boundary=str(case["x_boundary"]),
            y_boundary=str(case["y_boundary"]),
        )
        xx_full, yy_full = np.meshgrid(full_axis, full_axis, indexing="xy")
        flattened = full_values.reshape(3, -1)
        field_df = pd.DataFrame(
            {
                "x": xx_full.ravel(),
                "y": yy_full.ravel(),
                "Ex_re": np.real(flattened[0]),
                "Ex_im": np.imag(flattened[0]),
                "Ey_re": np.real(flattened[1]),
                "Ey_im": np.imag(flattened[1]),
                "Ez_re": np.real(flattened[2]),
                "Ez_im": np.imag(flattened[2]),
                "is_in_domain": full_mask.ravel(),
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
        "window_half_width_um": half_width,
        "grid_size": int(grid_size),
        "dx_um": float(full_axis[1] - full_axis[0]),
        "dy_um": float(full_axis[1] - full_axis[0]),
        "valid_mode_indices": valid_mode_indices,
        "exported_mode_indices": exported_mode_indices,
        "field_domain": "full field reconstructed from first-quadrant solution",
        "symmetry": case,
    }
    if finite_geometry == "hex":
        metadata["finite_hex_a_um"] = float(finite_hex_a)
    finite.write_config(export_dir / "air_field_metadata.json", metadata)
    return metadata


def export_reconstructed_candidate_hz_results(
    sim_run,
    export_dir: Path,
    validity_boundary: np.ndarray,
    *,
    case: dict[str, str | int],
) -> pd.DataFrame:
    """Export eigenvalues and the minimum Hz data needed for mode selection."""
    export_dir.mkdir(parents=True, exist_ok=True)
    eigenfrequencies = sim_run.get_eigenfrequencies()
    eigen_df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])
    validity = finite.mode_validity_for_boundary(
        sim_run,
        len(eigen_df),
        validity_boundary,
    )
    eigen_df = pd.concat([eigen_df, pd.DataFrame(validity)], axis=1)
    eigen_df["mode_idx"] = np.arange(len(eigen_df), dtype=int)
    eigen_df["symmetry_id"] = int(case["id"])
    eigen_df["x_boundary"] = str(case["x_boundary"])
    eigen_df["y_boundary"] = str(case["y_boundary"])
    eigen_df.to_csv(export_dir / "eigenfrequencies.csv", index=False)

    valid_mode_indices = eigen_df.loc[
        eigen_df["is_valid"].map(finite._csv_bool),
        "mode_idx",
    ].to_numpy(dtype=int).tolist()
    print(
        f"Valid finite-quarter symmetry {case['id']} modes: "
        f"{valid_mode_indices}"
    )

    for mode_idx in valid_mode_indices:
        hz_coordinates, hz = _reconstruct_hz_field(sim_run, mode_idx, case)
        pd.DataFrame(
            {
                "x": hz_coordinates[:, 0],
                "y": hz_coordinates[:, 1],
                "re": np.real(hz),
                "im": np.imag(hz),
            }
        ).to_parquet(export_dir / f"{mode_idx:02d}_Hz_center.parquet", index=False)
    return eigen_df


def export_reconstructed_selected_results(
    sim_run,
    export_dir: Path,
    eigen_df: pd.DataFrame,
    mode_indices: list[int],
    *,
    finite_hex_a: float,
    case: dict[str, str | int],
    finite_geometry: str = "hex",
    footprint_width: float | None = None,
    footprint_height: float | None = None,
) -> None:
    """Generate the complete field and air-plane exports for selected modes."""
    valid_mode_indices = eigen_df.loc[
        eigen_df["is_valid"].map(finite._csv_bool),
        "mode_idx",
    ].to_numpy(dtype=int).tolist()
    selected_mode_indices = [int(mode_idx) for mode_idx in mode_indices]
    if len(selected_mode_indices) != len(set(selected_mode_indices)):
        raise ValueError("Selected finite-quarter mode indices must be unique")
    invalid = sorted(set(selected_mode_indices).difference(valid_mode_indices))
    if invalid:
        raise ValueError(
            "Selected finite-quarter modes must be valid; "
            f"invalid modes: {invalid}"
        )
    for mode_idx in selected_mode_indices:
        mode_row = eigen_df.loc[eigen_df["mode_idx"] == mode_idx].iloc[0]
        frequency_thz = float(mode_row["re"])
        quality_factor = float(mode_row["q"])
        hz_frame = pd.read_parquet(
            export_dir / f"{mode_idx:02d}_Hz_center.parquet"
        )
        hz_coordinates = hz_frame[["x", "y"]].to_numpy(dtype=float)
        hz = hz_frame["re"].to_numpy(dtype=float) + 1j * hz_frame[
            "im"
        ].to_numpy(dtype=float)
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Hz_Re_2d.png",
            hz_coordinates,
            np.real(hz),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Re(Hz)",
            unit_label="A/m",
            color_scale_mode="linearsymmetric",
            dpi=DPI,
        )
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Hz_Im_2d.png",
            hz_coordinates,
            np.imag(hz),
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Im(Hz)",
            unit_label="A/m",
            color_scale_mode="linearsymmetric",
            dpi=DPI,
        )
        wem_coordinates, wem = _reconstruct_wem_field(sim_run, mode_idx)
        save_field_plot(
            export_dir / f"{mode_idx:02d}_Wem_2d.png",
            wem_coordinates,
            wem,
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label="Wem",
            unit_label="J/m³",
            color_scale_mode="linear",
            dpi=DPI,
        )

    if finite.RUN_FARFIELD_FFT:
        export_reconstructed_air_fields(
            sim_run,
            export_dir,
            eigen_df,
            finite_hex_a=finite_hex_a,
            grid_size=finite.FARFIELD_GRID_SIZE,
            h0_um=finite.FARFIELD_H0_UM,
            case=case,
            finite_geometry=finite_geometry,
            footprint_width=footprint_width,
            footprint_height=footprint_height,
            mode_indices=selected_mode_indices,
        )


def export_reconstructed_results(
    sim_run,
    export_dir: Path,
    validity_boundary: np.ndarray,
    *,
    finite_hex_a: float,
    case: dict[str, str | int],
    finite_geometry: str = "hex",
    footprint_width: float | None = None,
    footprint_height: float | None = None,
) -> pd.DataFrame:
    """Compatibility wrapper that fully exports every valid quarter mode."""
    eigen_df = export_reconstructed_candidate_hz_results(
        sim_run,
        export_dir,
        validity_boundary,
        case=case,
    )
    valid_mode_indices = eigen_df.loc[
        eigen_df["is_valid"].map(finite._csv_bool),
        "mode_idx",
    ].to_numpy(dtype=int).tolist()
    export_reconstructed_selected_results(
        sim_run,
        export_dir,
        eigen_df,
        valid_mode_indices,
        finite_hex_a=finite_hex_a,
        case=case,
        finite_geometry=finite_geometry,
        footprint_width=footprint_width,
        footprint_height=footprint_height,
    )
    return eigen_df


def build_quarter_case_geometry() -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, object],
    np.ndarray,
    float,
]:
    full_records, metadata = finite.build_full_lattice_holes()
    validate_full_geometry_mirror_symmetry(full_records)
    cell_records = finite.full_cell_records(metadata)
    if metadata.get("finite_geometry") == "square":
        full_boundary = np.asarray(metadata["finite_boundary"], dtype=float)
        finite_lx = float(metadata["finite_Lx"])
        finite_ly = float(metadata["finite_Ly"])
        finite_extent = max(finite_lx, finite_ly)
        quarter_boundary = quarter_rectangle_vertices(finite_lx, finite_ly)
    else:
        full_boundary, finite_extent, finite_lx, finite_ly = (
            finite.simulation_boundary_for_full_lattice_hexagon(cell_records)
        )
        quarter_boundary = quarter_boundary_vertices(finite_extent)
    quarter_records = clip_records_to_first_quadrant(
        full_records,
        quarter_boundary,
    )
    case_metadata = {
        **metadata,
        "finite_Lx": finite_lx,
        "finite_Ly": finite_ly,
        "finite_boundary": full_boundary,
        "simulation_cells": len(cell_records),
        "quarter_boundary": quarter_boundary,
        "full_hole_count": len(full_records),
        "quarter_hole_count": len(quarter_records),
    }
    if metadata.get("finite_geometry") != "square":
        case_metadata["finite_hex_a"] = finite_extent
    return (
        full_records,
        quarter_records,
        case_metadata,
        quarter_boundary,
        finite_extent,
    )


def run_symmetry_case(
    x_shift_factor: float,
    symmetry_id: int,
    full_records: list[dict[str, object]],
    quarter_records: list[dict[str, object]],
    geometry_metadata: dict[str, object],
    finite_hex_a: float,
    y_shift_factor: float | None = None,
    *,
    target_frequency_thz: float | None = None,
    target_band_label: str | None = None,
    target_internal_mode: str | None = None,
    unit_cell_source_path: str | Path | None = None,
    symmetry_ids: list[int] | None = None,
) -> dict[str, object]:
    resolved_y_factor = (
        float(x_shift_factor)
        if y_shift_factor is None
        else float(y_shift_factor)
    )
    case = symmetry_case(symmetry_id)
    active_symmetry_ids = validate_symmetry_ids(
        list(SYMMETRY_IDS) if symmetry_ids is None else list(symmetry_ids)
    )
    if symmetry_id not in active_symmetry_ids:
        raise ValueError(
            f"symmetry ID {symmetry_id} is not present in the active batch"
        )
    resolved_target_frequency = resolved_target_frequency_thz(
        target_frequency_thz
    )
    case_dir = quarter_output_case_dir(
        x_shift_factor,
        symmetry_id,
        active_symmetry_ids,
        y_shift_factor=(
            resolved_y_factor if y_shift_factor is not None else None
        ),
    )
    output_paths = finite.prepare_finite_case_output(
        case_dir,
        model_filename="finite_quarter.mph",
        resume=RESUME_FROM_EXISTING_MPH,
        create_staging=RUN_COMSOL,
    )
    export_dir = output_paths["export_dir"]
    mph_path = output_paths["mph"]
    config_path = output_paths["config"]
    progress_log_path = output_paths["progress_log"]
    validity_boundary = finite.strip_validity_boundary(geometry_metadata)

    metadata = {
        **finite.shared_case_metadata(geometry_metadata, validity_boundary),
        "case": "finite_quarter",
        "source": "first-quadrant symmetry reduction of finite cavity",
        "symmetry": case,
        "symmetry_ids": active_symmetry_ids,
        "target_mode": {
            "band_label": target_band_label,
            "internal_mode": target_internal_mode,
            "unit_cell_gamma_frequency_thz": resolved_target_frequency,
            "unit_cell_source_path": (
                None
                if unit_cell_source_path is None
                else str(Path(unit_cell_source_path))
            ),
        },
        "quarter_domain": {"x_min": 0.0, "y_min": 0.0},
        "field_output": "full field reconstructed by symmetry",
        "raw_quarter_exports_saved": False,
        "quarter_boundary": geometry_metadata["quarter_boundary"],
        "simulation_holes": len(quarter_records),
        "simulation": {
            "mode": "finite_quarter",
            "boundary": (
                "quarter_rectangle"
                if geometry_metadata.get("finite_geometry") == "square"
                else "quarter_hexagon_boundary"
            ),
            "x_axis_boundary": case["x_boundary"],
            "y_axis_boundary": case["y_boundary"],
        },
        "output_layout": finite.OUTPUT_LAYOUT_VERSION,
        "single_mode_analysis": {
            "enabled": finite.FINITE_SINGLE_MODE_ANALYSIS_ENABLED,
            "selection": (
                "target_component_centered_nodeless_envelope"
                if target_internal_mode is not None
                else "closest_to_target_frequency"
            ),
            "target_internal_mode": target_internal_mode,
            "target_frequency_thz": resolved_target_frequency,
            "candidate_mode_indices": (
                finite.FINITE_LATTICE_FOURIER_MODE_INDICES
            ),
            "candidate_export": "Hz_center_only",
            "full_processing_scope": "selected_mode_only",
        },
        "paths": finite.finite_case_output_metadata(output_paths),
    }
    for message in finite.shift_profile_preflight_messages(metadata):
        print(message)
    finite.save_finite_shift_profile_diagnostics(
        output_paths["overview_dir"],
        metadata,
        dpi=DPI,
    )
    full_cell_records = finite.full_cell_records(geometry_metadata)
    finite.save_finite_model_geometry_figures(
        output_paths["model_dir"],
        full_records,
        full_cell_records,
        np.asarray(geometry_metadata["finite_boundary"], dtype=float),
        validity_boundary,
        {
            **metadata,
            "simulation_holes": len(full_records),
        },
    )
    finite.write_config(config_path, metadata)
    print(
        f"Prepared finite-quarter symmetry ID {symmetry_id} at {case_dir}; "
        f"holes={len(quarter_records)}, x-axis={case['x_boundary']}, "
        f"y-axis={case['y_boundary']}"
    )
    print(f"COMSOL progress log: {progress_log_path}")

    if not RUN_COMSOL:
        return {
            "finite_geometry": finite.FINITE_GEOMETRY,
            "cladding_x_shift_factor": float(x_shift_factor),
            "cladding_y_shift_factor": resolved_y_factor,
            "cladding_shift_input_kind": finite.CLADDING_SHIFT_INPUT_KIND,
            "cladding_shift_base_factor": (
                finite.current_cladding_shift_base_factor()
            ),
            "cladding_y_over_x_shift_ratio": (
                finite.current_cladding_y_over_x_shift_ratio()
            ),
            "cladding_x_inward_shift": float(x_shift_factor) * finite.A,
            "cladding_y_inward_shift": resolved_y_factor * finite.A,
            **finite.finite_geometry_shift_metadata(),
            "symmetry_id": symmetry_id,
            "target_band_label": target_band_label,
            "target_internal_mode": target_internal_mode,
            "target_frequency_thz": resolved_target_frequency,
            "out_dir": str(case_dir),
            "status": "geometry_only",
            "output_layout": finite.OUTPUT_LAYOUT_VERSION,
            **finite.finite_case_output_metadata(output_paths),
        }

    if finite.RUN_FARFIELD_FFT:
        finite.validate_farfield_config(
            finite.FARFIELD_GRID_SIZE,
            finite.FARFIELD_FFT_SIZE,
            finite.FARFIELD_NA,
            finite.FARFIELD_H0_UM,
        )

    layers = finite.simulation_layers(
        finite.polygons_from_records(quarter_records)
    )
    with SimulationRun(
        config=quarter_simulation_config(
            case,
            target_frequency_thz=resolved_target_frequency,
        ),
        progress_log_path=progress_log_path,
    ) as sim:
        if RESUME_FROM_EXISTING_MPH:
            if not mph_path.is_file():
                raise FileNotFoundError(
                    f"Cannot resume finite quarter: missing checkpoint {mph_path}"
                )
            finite.attach_saved_model(sim, mph_path)
            finite.validate_geometry_checkpoint_model(sim)
        else:
            if geometry_metadata.get("finite_geometry") == "square":
                compiled = finite.compile_finite_geometry_input(
                    geometry_metadata["geometry_plan"],
                    layer_height=finite.SLAB_HEIGHT,
                    refractive_index=finite.REFRACTIVE_INDEX,
                    quarter=True,
                    holes=finite.polygons_from_records(quarter_records),
                )
                sim_boundary = compiled.boundary
                layers = list(compiled.layers)
            else:
                sim_boundary = BoundarySpec(
                    "quarter_hexagon_boundary", a=finite_hex_a
                )
            sim.build_geometry(
                sim_boundary,
                layers,
                simulation_mode="finite_quarter",
                k=None,
            )
            sim.model.save(str(mph_path))
        sim.run_simulation(
            mesh_save_path=mph_path,
            mesh_builder=finite.finite_mesh_builder(layers),
        )
        sim.model.save(str(mph_path))
        eigen_df = export_reconstructed_candidate_hz_results(
            sim,
            export_dir,
            validity_boundary,
            case=case,
        )
        if finite.FINITE_SINGLE_MODE_ANALYSIS_ENABLED:
            if target_internal_mode is None:
                analysis_mode_indices = finite.resolve_finite_analysis_mode_indices(
                    output_paths["staging_dir"],
                    export_dir,
                    target_frequency_thz=resolved_target_frequency,
                )
            else:
                selected_mode, _selection = (
                    finite.select_finite_fundamental_mode(
                        export_dir,
                        output_paths["staging_dir"]
                        / "analysis_mode_selection.csv",
                        cells=list(geometry_metadata["bulk_points"]),
                        period=finite.A,
                        rho_grid_size=(
                            finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE
                        ),
                        target_internal_mode=target_internal_mode,
                        target_frequency_thz=resolved_target_frequency,
                        candidate_mode_indices=(
                            finite.FINITE_LATTICE_FOURIER_MODE_INDICES
                        ),
                    )
                )
                analysis_mode_indices = [selected_mode]
        else:
            analysis_mode_indices = None
        full_export_mode_indices = (
            analysis_mode_indices
            if analysis_mode_indices is not None
            else eigen_df.loc[
                eigen_df["is_valid"].map(finite._csv_bool),
                "mode_idx",
            ].to_numpy(dtype=int).tolist()
        )
        export_reconstructed_selected_results(
            sim,
            export_dir,
            eigen_df,
            full_export_mode_indices,
            finite_hex_a=finite_hex_a,
            case=case,
            finite_geometry=str(
                geometry_metadata.get("finite_geometry", "hex")
            ),
            footprint_width=float(geometry_metadata["finite_Lx"]),
            footprint_height=float(geometry_metadata["finite_Ly"]),
        )
    lattice_mode_indices = (
        analysis_mode_indices
        if finite.FINITE_SINGLE_MODE_ANALYSIS_ENABLED
        else finite.FINITE_LATTICE_FOURIER_MODE_INDICES
    )
    if finite.RUN_FARFIELD_FFT:
        finite.run_farfield_fft_case(
            output_paths["staging_dir"],
            export_dir,
            grid_size=finite.FARFIELD_GRID_SIZE,
            fft_size=finite.FARFIELD_FFT_SIZE,
            na=finite.FARFIELD_NA,
            h0_um=finite.FARFIELD_H0_UM,
            period_um=finite.A,
            mode_indices=analysis_mode_indices,
        )
    if finite.RUN_FINITE_LATTICE_FOURIER_POSTPROCESS:
        finite.run_exported_finite_lattice_fourier_postprocess(
            output_paths["staging_dir"],
            export_dir,
            cells=list(geometry_metadata["bulk_points"]),
            cladding_cells=list(geometry_metadata["cladding_points"]),
            bulk_radius=finite.BULK_RADIUS,
            period=finite.A,
            rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            top_k=finite.FINITE_LATTICE_FOURIER_TOP_K,
            mode_indices=lattice_mode_indices,
            include_unselected_valid_mode_maps=(
                not finite.FINITE_SINGLE_MODE_ANALYSIS_ENABLED
            ),
            dpi=DPI,
            transform_kind=(
                "arbitrary_square_cavity_v1"
                if geometry_metadata.get("finite_geometry") == "square"
                else "hex_cyclic_quotient_v1"
            ),
        )
        finite.score_exported_finite_modes(
            output_paths["staging_dir"],
            export_dirname=finite.EXPORT_DIRNAME,
        )
    else:
        finite.run_exported_finite_unit_cell_intensity_maps(
            output_paths["staging_dir"],
            export_dir,
            cells=list(geometry_metadata["bulk_points"]),
            cladding_cells=list(geometry_metadata["cladding_points"]),
            period=finite.A,
            rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
            mode_indices=analysis_mode_indices,
            dpi=DPI,
        )
    output_paths = finite.finalize_finite_case_output(
        case_dir,
        model_filename="finite_quarter.mph",
    )
    return {
        "finite_geometry": finite.FINITE_GEOMETRY,
        "cladding_x_shift_factor": float(x_shift_factor),
        "cladding_y_shift_factor": resolved_y_factor,
        "cladding_shift_input_kind": finite.CLADDING_SHIFT_INPUT_KIND,
        "cladding_shift_base_factor": (
            finite.current_cladding_shift_base_factor()
        ),
        "cladding_y_over_x_shift_ratio": (
            finite.current_cladding_y_over_x_shift_ratio()
        ),
        "cladding_x_inward_shift": float(x_shift_factor) * finite.A,
        "cladding_y_inward_shift": resolved_y_factor * finite.A,
        **finite.finite_geometry_shift_metadata(),
        "symmetry_id": symmetry_id,
        "target_band_label": target_band_label,
        "target_internal_mode": target_internal_mode,
        "target_frequency_thz": resolved_target_frequency,
        "analysis_mode_indices": analysis_mode_indices,
        "out_dir": str(case_dir),
        "status": "complete",
        "mph": str(mph_path),
        "mode_output_pattern": str(
            output_paths["results_dir"] / "mode{mode_idx}"
        ),
        "output_layout": finite.OUTPUT_LAYOUT_VERSION,
        "model_dir": str(output_paths["model_dir"]),
        "overview_dir": str(output_paths["overview_dir"]),
        "logs_dir": str(output_paths["logs_dir"]),
        "config_dir": str(output_paths["config_dir"]),
    }


def run_finite_quarter_shift(
    x_shift_factor: float,
    y_shift_factor: float | None = None,
) -> list[dict[str, object]]:
    resolved_y_factor = (
        float(x_shift_factor)
        if y_shift_factor is None
        else float(y_shift_factor)
    )
    finite.set_cladding_inward_shift_factors(
        x_shift_factor,
        resolved_y_factor,
    )
    symmetry_ids = validate_symmetry_ids(list(SYMMETRY_IDS))
    (
        full_records,
        quarter_records,
        metadata,
        _boundary,
        finite_hex_a,
    ) = build_quarter_case_geometry()
    return [
        run_symmetry_case(
            x_shift_factor,
            symmetry_id,
            full_records,
            quarter_records,
            metadata,
            finite_hex_a,
            resolved_y_factor if y_shift_factor is not None else None,
        )
        for symmetry_id in symmetry_ids
    ]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []
    for x_shift_factor, y_shift_factor in (
        finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
    ):
        summaries.extend(
            run_finite_quarter_shift(x_shift_factor, y_shift_factor)
        )
    finite.write_config(
        OUT_DIR / "run_summary.json",
        {
            "finite_cavity_preset": finite.FINITE_CAVITY_PRESET,
            "finite_cavity_preset_config": finite.FINITE_CAVITY_CONFIG,
            "symmetry_ids": list(SYMMETRY_IDS),
            "runs": summaries,
            "output_layout": finite.OUTPUT_LAYOUT_VERSION,
        },
    )
    print(f"Wrote {len(summaries)} finite-quarter cases to {OUT_DIR}")
    for item in summaries:
        print(f"  {item['out_dir']}")


if __name__ == "__main__":
    main()
