#!/usr/bin/env python3
"""Recover full post-processing for one mode from a solved quarter MPH."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.hex_lattice_utils import LatticePoint  # noqa: E402
from comsol_workflow.simulation_utils import SimulationRun  # noqa: E402
from scripts.run_main import run_finite as finite  # noqa: E402
from scripts.run_main import run_finite_quarter as quarter  # noqa: E402


STAGING_DIRNAME = ".mode_recovery"


def _csv_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def _read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def _lattice_points(records: object, label: str) -> list[LatticePoint]:
    if not isinstance(records, list) or not records:
        raise ValueError(f"Case config {label} must be a non-empty list")
    points = []
    for idx, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError(f"Case config {label}[{idx}] must be an object")
        points.append(
            LatticePoint(
                i=int(record["i"]),
                j=int(record["j"]),
                shell=int(record["shell"]),
                x=float(record["x"]),
                y=float(record["y"]),
            )
        )
    return points


def _hardlink_or_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def corrected_selection_tables(
    eigen_df: pd.DataFrame,
    previous_selection: pd.DataFrame,
    audit_df: pd.DataFrame,
    *,
    selected_mode_idx: int,
    target_internal_mode: str,
    max_peak_shell: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid_modes = set(
        eigen_df.loc[eigen_df["is_valid"].map(_csv_bool), "mode_idx"]
        .astype(int)
        .tolist()
    )
    if selected_mode_idx not in valid_modes:
        raise ValueError(
            f"Selected mode {selected_mode_idx} is not a valid exported mode"
        )
    if "mode_idx" not in previous_selection or "mode_idx" not in audit_df:
        raise ValueError("Selection and audit tables must contain mode_idx")
    previous = previous_selection.copy()
    previous["mode_idx"] = previous["mode_idx"].astype(int)
    audit = audit_df.copy()
    audit["mode_idx"] = audit["mode_idx"].astype(int)
    missing_audit = sorted(set(previous["mode_idx"]).difference(audit["mode_idx"]))
    if missing_audit:
        raise ValueError(
            f"Fundamental audit is missing candidate modes: {missing_audit}"
        )
    audit_columns = [
        "mode_idx",
        "center_over_peak_cell",
        "center_over_peak_shell",
        "peak_shell",
        "central_energy_shell_le_2",
        "rms_radius_over_a",
        "outward_resurgence",
        "radial_sign_changes",
        "max_significant_shell_phase_degrees",
        "sampling_grid_side",
        "sampling_points_per_cell",
        "max_nearest_sample_distance_um",
    ]
    missing_columns = sorted(set(audit_columns).difference(audit.columns))
    # Native finite-quarter staging tables contain the envelope audit fields,
    # but older/interrupted runs may not include sampling diagnostics.  Those
    # diagnostics are informational and can be filled with NaN during recovery.
    required_audit_columns = {
        "mode_idx",
        "center_over_peak_cell",
        "center_over_peak_shell",
        "peak_shell",
        "outward_resurgence",
        "radial_sign_changes",
    }
    missing_required = sorted(required_audit_columns.difference(audit.columns))
    if missing_required:
        raise ValueError(
            f"Fundamental audit is missing columns: {missing_required}"
        )
    for column in missing_columns:
        audit[column] = np.nan
    selection = previous.drop(
        columns=[column for column in audit_columns if column != "mode_idx"],
        errors="ignore",
    ).merge(audit[audit_columns], on="mode_idx", how="left", validate="one_to_one")
    selection["target_internal_mode"] = str(target_internal_mode)
    selection["target_component_matches"] = (
        selection["dominant_mode"].astype(str).eq(str(target_internal_mode))
        & selection["target_component_weight"].astype(float).ge(0.5)
    )
    selection["fundamental_eligible"] = (
        selection["target_component_matches"]
        & selection["peak_shell"].astype(int).le(int(max_peak_shell))
        & selection["center_over_peak_shell"].astype(float).ge(0.9)
        & selection["radial_sign_changes"].astype(int).eq(0)
        & selection["outward_resurgence"].astype(float).le(0.05)
    )
    eligible_modes = selection.loc[
        selection["fundamental_eligible"], "mode_idx"
    ].astype(int).tolist()
    if eligible_modes != [selected_mode_idx]:
        raise ValueError(
            "Approved mode must be the unique fundamental candidate; "
            f"eligible modes: {eligible_modes}"
        )
    selection["analysis_selected"] = selection["mode_idx"].eq(
        selected_mode_idx
    )
    selection["selection_reason"] = "rejected_higher_order_envelope"
    selection.loc[
        selection["analysis_selected"], "selection_reason"
    ] = (
        "selected_centered_nodeless_boundary_peak_fundamental_envelope"
        if int(selection.loc[selection["analysis_selected"], "peak_shell"].iloc[0]) > 0
        else "selected_centered_nodeless_fundamental_envelope"
    )
    selection = selection.sort_values("mode_idx")

    corrected_eigen = eigen_df.copy()
    corrected_eigen["analysis_selected"] = corrected_eigen["mode_idx"].astype(
        int
    ).eq(selected_mode_idx)
    corrected_eigen["analysis_selection_method"] = (
        "target_component_centered_nodeless_envelope"
    )
    corrected_eigen["analysis_target_internal_mode"] = target_internal_mode
    indexed = selection.set_index("mode_idx")
    for column in (
        "target_component_weight",
        "fundamental_eligible",
        "peak_shell",
        "center_over_peak_shell",
        "outward_resurgence",
        "radial_sign_changes",
        "mode_px",
        "mode_py",
        "mode_dx",
        "mode_dy",
    ):
        corrected_eigen[column] = corrected_eigen["mode_idx"].astype(int).map(
            indexed[column]
        )
    return corrected_eigen, selection


def validate_recovery_inputs(
    case_dir: Path,
    selected_mode_idx: int,
    audit_path: Path,
    backup_dir: Path,
    *,
    resume_staging: bool = False,
    allow_existing_staging: bool = False,
    initial_staging: bool = False,
) -> dict[str, object]:
    case_dir = case_dir.resolve()
    if not case_dir.is_dir():
        raise FileNotFoundError(case_dir)
    if backup_dir.exists():
        raise FileExistsError(f"Backup target already exists: {backup_dir}")
    staging_dir = case_dir / STAGING_DIRNAME
    if initial_staging and not staging_dir.exists():
        native_staging = case_dir / ".staging"
        if native_staging.is_dir():
            staging_dir = native_staging
    if staging_dir.exists() and not (resume_staging or allow_existing_staging):
        raise FileExistsError(f"Recovery staging already exists: {staging_dir}")
    if resume_staging and not staging_dir.is_dir():
        raise FileNotFoundError(
            f"Recovery staging does not exist for --resume-staging: {staging_dir}"
        )
    mph_path = case_dir / "00_model" / "finite_quarter.mph"
    if not mph_path.is_file():
        raise FileNotFoundError(mph_path)
    config_path = case_dir / "99_config" / "config.json"
    config = _read_json(config_path)
    if str(config.get("case")) != "finite_quarter":
        raise ValueError("Recovery requires a finite_quarter case")
    target_mode = config.get("target_mode")
    if not isinstance(target_mode, dict):
        raise ValueError("Case config is missing target_mode")
    eigen_path = case_dir / "10_overview" / "eigenfrequencies.csv"
    selection_path = case_dir / "10_overview" / "analysis_mode_selection.csv"
    if initial_staging:
        eigen_path = staging_dir / finite.EXPORT_DIRNAME / "eigenfrequencies.csv"
        selection_path = staging_dir / "analysis_mode_selection.csv"
    if not eigen_path.is_file() or not selection_path.is_file():
        raise FileNotFoundError("Case is missing eigenfrequency or selection tables")
    mode_dir = case_dir / "01_results" / f"mode{selected_mode_idx}"
    hz_path = mode_dir / "11_simulation_exports" / "Hz_center.parquet"
    if initial_staging:
        hz_path = staging_dir / finite.EXPORT_DIRNAME / f"{selected_mode_idx:02d}_Hz_center.parquet"
    if not hz_path.is_file():
        raise FileNotFoundError(hz_path)
    if not audit_path.is_file():
        raise FileNotFoundError(audit_path)
    for forbidden in (
        mode_dir / "12_farfield_FFT",
        mode_dir / "13_lattice_fourier_Hz",
        mode_dir / "11_simulation_exports" / "E_air.parquet",
        mode_dir / "11_simulation_exports" / "Wem_2d.png",
    ):
        if forbidden.exists():
            raise FileExistsError(
                f"Selected mode already has full processing output: {forbidden}"
            )
    return {
        "case_dir": case_dir,
        "backup_dir": backup_dir.resolve(),
        "staging_dir": staging_dir,
        "mph_path": mph_path,
        "config_path": config_path,
        "config": config,
        "eigen_path": eigen_path,
        "selection_path": selection_path,
        "audit_path": audit_path.resolve(),
        "mode_dir": mode_dir,
        "hz_path": hz_path,
        "initial_staging": bool(initial_staging),
    }


def _move_to_backup(source: Path, backup_root: Path, case_dir: Path) -> None:
    if not source.exists():
        return
    relative = source.relative_to(case_dir)
    destination = backup_root / relative
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))


def recover_selected_mode(
    case_dir: Path,
    selected_mode_idx: int,
    audit_path: Path,
    backup_dir: Path,
    *,
    resume_staging: bool = False,
    allow_existing_staging: bool = False,
    max_peak_shell: int = 0,
    initial_staging: bool = False,
) -> None:
    context = validate_recovery_inputs(
        case_dir,
        selected_mode_idx,
        audit_path,
        backup_dir,
        resume_staging=resume_staging,
        allow_existing_staging=allow_existing_staging,
        initial_staging=initial_staging,
    )
    case_dir = context["case_dir"]
    staging_dir = context["staging_dir"]
    export_dir = staging_dir / finite.EXPORT_DIRNAME
    export_dir.mkdir(
        parents=True,
        exist_ok=resume_staging or allow_existing_staging,
    )

    eigen_df = pd.read_csv(context["eigen_path"])
    previous_selection = pd.read_csv(context["selection_path"])
    audit_df = pd.read_csv(context["audit_path"])
    config = dict(context["config"])
    target_mode = dict(config["target_mode"])
    target_internal_mode = str(target_mode["internal_mode"])
    target_frequency_thz = float(
        target_mode["unit_cell_gamma_frequency_thz"]
    )
    corrected_eigen, corrected_selection = corrected_selection_tables(
        eigen_df,
        previous_selection,
        audit_df,
        selected_mode_idx=selected_mode_idx,
        target_internal_mode=target_internal_mode,
        max_peak_shell=max_peak_shell,
    )
    corrected_eigen.to_csv(export_dir / "eigenfrequencies.csv", index=False)
    corrected_selection.to_csv(
        staging_dir / "analysis_mode_selection.csv",
        index=False,
    )
    staged_hz_path = export_dir / f"{selected_mode_idx:02d}_Hz_center.parquet"
    if not resume_staging and not initial_staging:
        _hardlink_or_copy(context["hz_path"], staged_hz_path)

    symmetry = dict(config["symmetry"])
    finite_geometry = str(config.get("finite_geometry", "hex"))
    finite_hex_a = float(config["finite_hex_a"])
    finite_lx = float(config["finite_Lx"])
    finite_ly = float(config["finite_Ly"])
    if not resume_staging:
        with SimulationRun(
            config=quarter.quarter_simulation_config(
                symmetry,
                target_frequency_thz=target_frequency_thz,
            )
        ) as simulation:
            finite.attach_saved_solution(simulation, context["mph_path"])
            solved = np.asarray(simulation.get_eigenfrequencies(), dtype=float)
            if solved.ndim != 2 or solved.shape[0] != len(eigen_df):
                raise ValueError(
                    "Solved MPH eigenmode count does not match finalized table"
                )
            saved_frequency = float(
                eigen_df.loc[
                    eigen_df["mode_idx"].astype(int).eq(selected_mode_idx), "re"
                ].iloc[0]
            )
            if not np.isclose(
                solved[selected_mode_idx, 0],
                saved_frequency,
                rtol=0.0,
                atol=1e-8,
            ):
                raise ValueError(
                    "Solved MPH mode frequency does not match finalized table"
                )
            quarter.export_reconstructed_selected_results(
                simulation,
                export_dir,
                corrected_eigen,
                [selected_mode_idx],
                finite_hex_a=finite_hex_a,
                case=symmetry,
                finite_geometry=finite_geometry,
                footprint_width=finite_lx,
                footprint_height=finite_ly,
            )

        finite.run_farfield_fft_case(
            staging_dir,
            export_dir,
            grid_size=finite.FARFIELD_GRID_SIZE,
            fft_size=finite.FARFIELD_FFT_SIZE,
            na=finite.FARFIELD_NA,
            h0_um=finite.FARFIELD_H0_UM,
            period_um=finite.A,
            mode_indices=[selected_mode_idx],
        )
    else:
        required_staged = [
            staged_hz_path,
            export_dir / f"{selected_mode_idx:02d}_Hz_Re_2d.png",
            export_dir / f"{selected_mode_idx:02d}_Hz_Im_2d.png",
            export_dir / f"{selected_mode_idx:02d}_Wem_2d.png",
            export_dir / f"{selected_mode_idx:02d}_E_air.parquet",
            export_dir / "air_field_metadata.json",
            staging_dir / "farfield_fft" / "farfield_summary.csv",
            staging_dir / "farfield_fft" / f"mode_{selected_mode_idx:02d}" / "summary.json",
        ]
        missing_staged = [path for path in required_staged if not path.is_file()]
        if missing_staged:
            raise FileNotFoundError(
                f"Incomplete recovery staging: {missing_staged[0]}"
            )
        print("Resuming from staged mode fields; COMSOL will not be started.")
    bulk_points = _lattice_points(config["bulk_points"], "bulk_points")
    cladding_points = _lattice_points(
        config["cladding_points"],
        "cladding_points",
    )
    finite.run_exported_finite_lattice_fourier_postprocess(
        staging_dir,
        export_dir,
        cells=bulk_points,
        cladding_cells=cladding_points,
        bulk_radius=int(config["bulk_radius"]),
        period=float(config["a"]),
        rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
        top_k=finite.FINITE_LATTICE_FOURIER_TOP_K,
        mode_indices=[selected_mode_idx],
        include_unselected_valid_mode_maps=False,
        dpi=quarter.DPI,
        transform_kind=(
            "arbitrary_square_cavity_v1"
            if finite_geometry == "square"
            else "hex_cyclic_quotient_v1"
        ),
    )
    finite.score_exported_finite_modes(
        staging_dir,
        export_dirname=finite.EXPORT_DIRNAME,
    )

    config["single_mode_analysis"] = {
        **dict(config["single_mode_analysis"]),
        "selection": "target_component_centered_nodeless_envelope",
        "approved_selected_mode_idx": int(selected_mode_idx),
        "full_processing_scope": "selected_mode_only",
    }
    finite.write_config(staging_dir / "config.json", config)

    old_selected = eigen_df.loc[
        eigen_df.get("analysis_selected", pd.Series(False, index=eigen_df.index)).map(_csv_bool),
        "mode_idx",
    ].astype(int).tolist()
    if initial_staging:
        old_mode_idx = None
    else:
        if len(old_selected) != 1:
            raise ValueError(
                f"Expected one previous selected mode, got {old_selected}"
            )
        old_mode_idx = old_selected[0]
    backup_root = context["backup_dir"]
    backup_root.mkdir(parents=True)
    if old_mode_idx is not None:
        old_mode_dir = case_dir / "01_results" / f"mode{old_mode_idx}"
        for source in (
            old_mode_dir / "12_farfield_FFT",
            old_mode_dir / "13_lattice_fourier_Hz",
            old_mode_dir / "10_overview" / "mode_score.csv",
        ):
            _move_to_backup(source, backup_root, case_dir)
        for source in (
            old_mode_dir / "11_simulation_exports" / "Hz_Re_2d.png",
            old_mode_dir / "11_simulation_exports" / "Hz_Im_2d.png",
            old_mode_dir / "11_simulation_exports" / "Wem_2d.png",
            old_mode_dir / "11_simulation_exports" / "E_air.parquet",
        ):
            _move_to_backup(source, backup_root, case_dir)
    for source in (
        case_dir / "10_overview" / "eigenfrequencies.csv",
        case_dir / "10_overview" / "analysis_mode_selection.csv",
        case_dir / "10_overview" / "farfield_summary.csv",
        case_dir / "10_overview" / "finite_lattice_fourier_summary.csv",
        case_dir / "10_overview" / "mode_scores.csv",
        case_dir / "99_config" / "config.json",
        case_dir / "99_config" / "air_field_metadata.json",
        case_dir / "99_config" / "farfield_config.json",
        case_dir / "99_config" / "objective.json",
        context["mode_dir"] / "10_overview" / "eigenfrequency.csv",
    ):
        _move_to_backup(source, backup_root, case_dir)

    selected_mode_dir = context["mode_dir"]
    selected_exports = selected_mode_dir / "11_simulation_exports"
    (selected_mode_dir / "10_overview").mkdir(parents=True, exist_ok=True)
    selected_exports.mkdir(parents=True, exist_ok=True)
    for staged_name, final_name in (
        (f"{selected_mode_idx:02d}_Hz_Re_2d.png", "Hz_Re_2d.png"),
        (f"{selected_mode_idx:02d}_Hz_Im_2d.png", "Hz_Im_2d.png"),
        (f"{selected_mode_idx:02d}_Wem_2d.png", "Wem_2d.png"),
        (f"{selected_mode_idx:02d}_E_air.parquet", "E_air.parquet"),
    ):
        shutil.move(
            str(export_dir / staged_name),
            str(selected_exports / final_name),
        )
    shutil.move(
        str(export_dir / "air_field_metadata.json"),
        str(case_dir / "99_config" / "air_field_metadata.json"),
    )
    shutil.move(
        str(export_dir / "eigenfrequencies.csv"),
        str(case_dir / "10_overview" / "eigenfrequencies.csv"),
    )
    shutil.move(
        str(staging_dir / "analysis_mode_selection.csv"),
        str(case_dir / "10_overview" / "analysis_mode_selection.csv"),
    )
    shutil.move(
        str(staging_dir / "config.json"),
        str(case_dir / "99_config" / "config.json"),
    )

    farfield_dir = staging_dir / "farfield_fft"
    shutil.move(
        str(farfield_dir / "config.json"),
        str(case_dir / "99_config" / "farfield_config.json"),
    )
    shutil.move(
        str(farfield_dir / "farfield_summary.csv"),
        str(case_dir / "10_overview" / "farfield_summary.csv"),
    )
    final_farfield = selected_mode_dir / "12_farfield_FFT"
    staged_farfield = farfield_dir / f"mode_{selected_mode_idx:02d}"
    staged_farfield.rename(final_farfield)
    staged_summary = final_farfield / "summary.json"
    if not staged_summary.is_file():
        raise FileNotFoundError(staged_summary)
    staged_summary.rename(final_farfield / "farfield_summary.json")

    lattice_dir = staging_dir / "finite_lattice_fourier_hz"
    shutil.move(
        str(lattice_dir / "finite_lattice_fourier_summary.csv"),
        str(case_dir / "10_overview" / "finite_lattice_fourier_summary.csv"),
    )
    final_lattice = selected_mode_dir / "13_lattice_fourier_Hz"
    (lattice_dir / f"mode_{selected_mode_idx:02d}").rename(final_lattice)
    shutil.move(
        str(staging_dir / "mode_scores.csv"),
        str(case_dir / "10_overview" / "mode_scores.csv"),
    )
    shutil.move(
        str(staging_dir / "objective.json"),
        str(case_dir / "99_config" / "objective.json"),
    )

    corrected_eigen.loc[
        corrected_eigen["mode_idx"].astype(int).eq(selected_mode_idx)
    ].to_csv(
        selected_mode_dir / "10_overview" / "eigenfrequency.csv",
        index=False,
    )
    score_df = pd.read_csv(case_dir / "10_overview" / "mode_scores.csv")
    score_df.loc[
        score_df["mode_idx"].astype(int).eq(selected_mode_idx)
    ].to_csv(
        selected_mode_dir / "10_overview" / "mode_score.csv",
        index=False,
    )
    finite._rewrite_farfield_outputs(
        finite.finite_case_output_paths(
            case_dir,
            model_filename="finite_quarter.mph",
        ),
        [selected_mode_idx],
    )
    finite._rewrite_lattice_summary(
        finite.finite_case_output_paths(
            case_dir,
            model_filename="finite_quarter.mph",
        )
    )

    for root_file in (
        lattice_dir / "index_r_cavity.png",
        lattice_dir / "index_k_1stBZ.png",
    ):
        root_file.unlink(missing_ok=True)
    linked_hz = export_dir / f"{selected_mode_idx:02d}_Hz_center.parquet"
    linked_hz.unlink(missing_ok=True)
    for directory in sorted(
        (path for path in staging_dir.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    ):
        directory.rmdir()
    staging_dir.rmdir()
    print(
        f"Recovered finite-quarter mode {selected_mode_idx} at {case_dir}; "
        f"previous mode {old_mode_idx} outputs backed up to {backup_root}"
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--mode-index", type=int, required=True)
    parser.add_argument("--audit-csv", type=Path, required=True)
    parser.add_argument("--backup-dir", type=Path, required=True)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Load COMSOL and install the recovered selected-mode outputs.",
    )
    parser.add_argument(
        "--resume-staging",
        action="store_true",
        help="Reuse already exported mode fields and resume after COMSOL export.",
    )
    parser.add_argument(
        "--allow-existing-staging",
        action="store_true",
        help="Reuse an existing partial staging while still running COMSOL export.",
    )
    parser.add_argument(
        "--max-peak-shell",
        type=int,
        default=0,
        help="Allow a selected fundamental envelope peak up to this shell.",
    )
    parser.add_argument(
        "--initial-staging",
        action="store_true",
        help="Finalize an interrupted case whose eigen/selection tables exist only in .staging.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    context = validate_recovery_inputs(
        args.case_dir,
        args.mode_index,
        args.audit_csv,
        args.backup_dir,
        resume_staging=args.resume_staging,
        allow_existing_staging=args.allow_existing_staging,
        initial_staging=args.initial_staging,
    )
    print(
        "Validated selected-mode recovery: "
        f"case={context['case_dir']}, mode={args.mode_index}, "
        f"mph={context['mph_path']}, backup={context['backup_dir']}"
    )
    config = dict(context["config"])
    target_mode = dict(config["target_mode"])
    corrected_selection_tables(
        pd.read_csv(context["eigen_path"]),
        pd.read_csv(context["selection_path"]),
        pd.read_csv(context["audit_path"]),
        selected_mode_idx=args.mode_index,
        target_internal_mode=str(target_mode["internal_mode"]),
        max_peak_shell=args.max_peak_shell,
    )
    print(
        "Validated approved fundamental selection: "
        f"mode={args.mode_index}, internal_mode={target_mode['internal_mode']}"
    )
    if not args.apply:
        print("Dry run completed; COMSOL was not started.")
        return
    recover_selected_mode(
        args.case_dir,
        args.mode_index,
        args.audit_csv,
        args.backup_dir,
        resume_staging=args.resume_staging,
        allow_existing_staging=args.allow_existing_staging,
        max_peak_shell=args.max_peak_shell,
        initial_staging=args.initial_staging,
    )


if __name__ == "__main__":
    main()
