#!/usr/bin/env python3
"""Build complete finite-cavity analyses for every mode selected in overview."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from comsol_workflow.simulation_utils import SimulationRun
from scripts.analysis.recover_finite_quarter_selected_mode import _lattice_points
from scripts.run_main import run_finite as finite
from scripts.run_main import run_finite_quarter as quarter


TAG_PATTERN = re.compile(r"_(f\d+\.\d+_e\d+_\d{8}(?:_retry\d+)?)$")
REQUIRED_FILES = (
    "10_overview/eigenfrequency.csv",
    "10_overview/mode_score.csv",
    "11_simulation_exports/Hz_center.parquet",
    "11_simulation_exports/Hz_Re_2d.png",
    "11_simulation_exports/Hz_Im_2d.png",
    "11_simulation_exports/Wem_2d.png",
    "11_simulation_exports/E_air.parquet",
    "12_farfield_FFT/A_overview.png",
    "12_farfield_FFT/B_polarization.png",
    "12_farfield_FFT/C_cutlines.png",
    "12_farfield_FFT/D_gauss_fit.png",
    "12_farfield_FFT/farfield_summary.json",
    "13_lattice_fourier_Hz/finite_lattice_fourier_hz.npz",
    "13_lattice_fourier_Hz/intensity_maps_normH.png",
    "13_lattice_fourier_Hz/k_weight_Hz.png",
    "13_lattice_fourier_Hz/k_weight_tops_norm.png",
    "13_lattice_fourier_Hz/k_weight_tops_phase.png",
    "13_lattice_fourier_Hz/p_subspace_mode_decomposition.csv",
    "13_lattice_fourier_Hz/top_k_peaks.csv",
)


def _tag(source_run: str) -> str:
    match = TAG_PATTERN.search(source_run)
    if not match:
        raise ValueError(f"Cannot derive source tag from {source_run}")
    return match.group(1)


def _read_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def _copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def _mode_record(target: Path, row: pd.Series) -> dict[str, object]:
    source_run = str(row["source_run"])
    tag = _tag(source_run)
    symmetry_name = (
        f"symmetry_{int(row['symmetry_id'])}_"
        f"x{row['x_boundary']}_y{row['y_boundary']}"
    )
    legacy = source_run.startswith("fq4m_retry_")
    if legacy:
        mode_dir = target / "01_results" / tag / symmetry_name / f"mode{int(row['mode_idx'])}"
        config_dir = target / "99_config" / "source_runs" / tag / symmetry_name
        model_dir = target / "00_model" / tag / symmetry_name
        overview_dir = target / "10_overview" / "source_runs" / tag / symmetry_name
        model_name = "finite_quarter.mph"
    else:
        mode_dir = target / "01_results" / tag / str(row["source_mode_uid"])
        config_dir = target / "99_config" / "source_runs" / tag
        model_dir = target / "00_model" / tag
        overview_dir = target / "10_overview" / "source_runs" / tag
        model_name = "finite_quarter_all.mph"
    record = {
        "name": str(row["overview_prefix"]),
        "mode_idx": int(row["mode_idx"]),
        "frequency_thz": float(row["re"]),
        "q": float(row["q"]),
        "tag": tag,
        "legacy": legacy,
        "symmetry_name": symmetry_name,
        "source_mode_dir": mode_dir,
        "config_path": config_dir / "config.json",
        "model_path": model_dir / model_name,
        "eigen_path": overview_dir / "eigenfrequencies.csv",
        "output_dir": target / "10_overview" / str(row["overview_prefix"]),
    }
    for key in ("source_mode_dir", "config_path", "model_path", "eigen_path"):
        if not Path(record[key]).exists():
            raise FileNotFoundError(record[key])
    record["complete"] = all((mode_dir / relative).is_file() for relative in REQUIRED_FILES)
    if not record["complete"] and not legacy:
        raise ValueError(f"Incomplete all-mode archive is not supported: {mode_dir}")
    return record


def load_records(target: Path) -> list[dict[str, object]]:
    selection_path = target / "10_overview" / "q_gt_1000_modes.csv"
    frame = pd.read_csv(selection_path)
    required = {
        "overview_prefix", "source_run", "source_mode_uid", "mode_idx",
        "symmetry_id", "x_boundary", "y_boundary", "re", "q",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Selection table is missing columns: {missing}")
    if frame.empty or not np.all(np.isfinite(frame["q"]) & (frame["q"] > 1000.0)):
        raise ValueError("Overview selection must contain only finite Q > 1000 modes")
    if frame["overview_prefix"].duplicated().any():
        raise ValueError("overview_prefix values must be unique")
    records = [_mode_record(target, row) for _, row in frame.iterrows()]
    conflicts = [record["output_dir"] for record in records if Path(record["output_dir"]).exists()]
    if conflicts:
        raise FileExistsError(conflicts[0])
    return records


def _validate_mode_dir(path: Path, frequency_thz: float, q: float) -> dict[str, object]:
    missing = [relative for relative in REQUIRED_FILES if not (path / relative).is_file()]
    if missing:
        raise FileNotFoundError(path / missing[0])
    if any((path / relative).stat().st_size == 0 for relative in REQUIRED_FILES):
        raise IOError(f"Empty analysis file in {path}")
    eigen = pd.read_csv(path / "10_overview" / "eigenfrequency.csv")
    if len(eigen) != 1 or not np.isclose(float(eigen.iloc[0]["re"]), frequency_thz, atol=1e-8, rtol=0.0):
        raise ValueError(f"Eigenfrequency mismatch in {path}")
    if not np.isclose(float(eigen.iloc[0]["q"]), q, atol=1e-8, rtol=1e-10):
        raise ValueError(f"Q mismatch in {path}")
    farfield = _read_json(path / "12_farfield_FFT" / "farfield_summary.json")
    if str(farfield.get("status")) != "ok":
        raise ValueError(f"Far-field analysis is not complete in {path}")
    pd.read_parquet(path / "11_simulation_exports" / "Hz_center.parquet")
    pd.read_parquet(path / "11_simulation_exports" / "E_air.parquet")
    with np.load(path / "13_lattice_fourier_Hz" / "finite_lattice_fourier_hz.npz") as archive:
        if not archive.files:
            raise ValueError(f"Empty Fourier archive in {path}")
    pngs = list(path.rglob("*.png"))
    for png in pngs:
        with Image.open(png) as image:
            image.verify()
    files = [item for item in path.rglob("*") if item.is_file()]
    return {"file_count": len(files), "bytes": sum(item.stat().st_size for item in files), "png_count": len(pngs)}


def _copy_complete(record: dict[str, object], outputs: Path) -> None:
    shutil.copytree(
        record["source_mode_dir"],
        outputs / str(record["name"]),
        copy_function=shutil.copyfile,
    )


def _farfield_contract(target: Path, tag: str) -> dict[str, object]:
    candidates = sorted((target / "99_config" / "source_runs" / tag).rglob("farfield_config.json"))
    if not candidates:
        raise FileNotFoundError(f"No completed far-field config for {tag}")
    contract = _read_json(candidates[0])
    expected = {
        "grid_size": finite.FARFIELD_GRID_SIZE,
        "fft_size": finite.FARFIELD_FFT_SIZE,
        "na": finite.FARFIELD_NA,
        "period_um": finite.A,
    }
    for key, value in expected.items():
        if not np.isclose(float(contract[key]), float(value), atol=0.0, rtol=0.0):
            raise ValueError(f"Current {key} differs from saved far-field contract")
    return contract


def _analyze_group(target: Path, records: list[dict[str, object]], staging: Path) -> None:
    first = records[0]
    config = _read_json(Path(first["config_path"]))
    case = dict(config["symmetry"])
    mode_indices = [int(record["mode_idx"]) for record in records]
    group_dir = staging / "groups" / str(first["symmetry_name"])
    export_dir = group_dir / finite.EXPORT_DIRNAME
    export_dir.mkdir(parents=True)
    eigen_df = pd.read_csv(first["eigen_path"])
    eigen_df.to_csv(export_dir / "eigenfrequencies.csv", index=False)
    valid_modes = set(
        eigen_df.loc[eigen_df["is_valid"].map(finite._csv_bool), "mode_idx"].astype(int)
    )
    if not set(mode_indices).issubset(valid_modes):
        raise ValueError(f"Selected invalid mode in {first['symmetry_name']}")
    for record in records:
        _copy_file(
            Path(record["source_mode_dir"]) / "11_simulation_exports" / "Hz_center.parquet",
            export_dir / f"{int(record['mode_idx']):02d}_Hz_center.parquet",
        )

    target_frequency = float(dict(config["single_mode_analysis"])["target_frequency_thz"])
    with SimulationRun(config=quarter.quarter_simulation_config(case, target_frequency_thz=target_frequency)) as simulation:
        finite.attach_saved_solution(simulation, Path(first["model_path"]))
        solved = np.asarray(simulation.get_eigenfrequencies(), dtype=float)
        if solved.ndim != 2 or solved.shape[0] != len(eigen_df):
            raise ValueError(f"Saved solution count mismatch for {first['symmetry_name']}")
        for record in records:
            mode_idx = int(record["mode_idx"])
            if not np.isclose(solved[mode_idx, 0], float(record["frequency_thz"]), atol=1e-8, rtol=0.0):
                raise ValueError(f"Saved solution frequency mismatch for {record['name']}")
        quarter.export_reconstructed_selected_results(
            simulation,
            export_dir,
            eigen_df,
            mode_indices,
            finite_hex_a=float(config["finite_hex_a"]),
            case=case,
            finite_geometry=str(config.get("finite_geometry", "hex")),
            footprint_width=float(config["finite_Lx"]),
            footprint_height=float(config["finite_Ly"]),
        )

    farfield = _farfield_contract(target, str(first["tag"]))
    farfield_summary = finite.run_farfield_fft_case(
        group_dir,
        export_dir,
        grid_size=int(farfield["grid_size"]),
        fft_size=int(farfield["fft_size"]),
        na=float(farfield["na"]),
        h0_um=None,
        period_um=float(farfield["period_um"]),
        mode_indices=mode_indices,
    )
    statuses = dict(zip(farfield_summary["mode_idx"].astype(int), farfield_summary["status"].astype(str)))
    if any(statuses.get(mode_idx) != "ok" for mode_idx in mode_indices):
        raise RuntimeError(f"Far-field failure in {first['symmetry_name']}: {statuses}")
    finite.run_exported_finite_lattice_fourier_postprocess(
        group_dir,
        export_dir,
        cells=_lattice_points(config["bulk_points"], "bulk_points"),
        cladding_cells=_lattice_points(config["cladding_points"], "cladding_points"),
        bulk_radius=int(config["bulk_radius"]),
        period=float(config["a"]),
        rho_grid_size=finite.FINITE_LATTICE_FOURIER_RHO_GRID_SIZE,
        top_k=finite.FINITE_LATTICE_FOURIER_TOP_K,
        mode_indices=mode_indices,
        include_unselected_valid_mode_maps=False,
        dpi=quarter.DPI,
        transform_kind="hex_cyclic_quotient_v1",
    )
    finite.score_exported_finite_modes(group_dir, export_dirname=finite.EXPORT_DIRNAME)
    score_df = pd.read_csv(group_dir / "mode_scores.csv")

    for record in records:
        mode_idx = int(record["mode_idx"])
        output = staging / "outputs" / str(record["name"])
        for section in ("10_overview", "11_simulation_exports"):
            (output / section).mkdir(parents=True, exist_ok=True)
        eigen_df.loc[eigen_df["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            output / "10_overview" / "eigenfrequency.csv", index=False
        )
        score_df.loc[score_df["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            output / "10_overview" / "mode_score.csv", index=False
        )
        for staged_name, final_name in (
            (f"{mode_idx:02d}_Hz_center.parquet", "Hz_center.parquet"),
            (f"{mode_idx:02d}_Hz_Re_2d.png", "Hz_Re_2d.png"),
            (f"{mode_idx:02d}_Hz_Im_2d.png", "Hz_Im_2d.png"),
            (f"{mode_idx:02d}_Wem_2d.png", "Wem_2d.png"),
            (f"{mode_idx:02d}_E_air.parquet", "E_air.parquet"),
        ):
            _copy_file(export_dir / staged_name, output / "11_simulation_exports" / final_name)
        shutil.move(str(group_dir / "farfield_fft" / f"mode_{mode_idx:02d}"), output / "12_farfield_FFT")
        summary_path = output / "12_farfield_FFT" / "summary.json"
        summary_path.rename(output / "12_farfield_FFT" / "farfield_summary.json")
        shutil.move(
            str(group_dir / "finite_lattice_fourier_hz" / f"mode_{mode_idx:02d}"),
            output / "13_lattice_fourier_Hz",
        )


def run(target: Path, apply: bool) -> dict[str, object]:
    target = target.resolve(strict=True)
    records = load_records(target)
    staging = target.parent / ".fq_mode_analysis_staging"
    if staging.exists():
        raise FileExistsError(staging)
    report: dict[str, object] = {
        "workflow": "finite_quarter_overview_selected_mode_analysis",
        "status": "dry_run" if not apply else "complete",
        "target": str(target),
        "selection": str(target / "10_overview" / "q_gt_1000_modes.csv"),
        "mode_count": len(records),
        "reused_modes": [record["name"] for record in records if record["complete"]],
        "generated_modes": [record["name"] for record in records if not record["complete"]],
        "groups_to_load": sorted({record["symmetry_name"] for record in records if not record["complete"]}),
    }
    if not apply:
        return report

    (staging / "outputs").mkdir(parents=True)
    for record in records:
        if record["complete"]:
            print(f"Reusing complete analysis: {record['name']}", flush=True)
            _copy_complete(record, staging / "outputs")
    incomplete = [record for record in records if not record["complete"]]
    groups: dict[str, list[dict[str, object]]] = {}
    for record in incomplete:
        groups.setdefault(str(record["symmetry_name"]), []).append(record)
    for index, (name, group_records) in enumerate(groups.items(), start=1):
        print(f"[{index}/{len(groups)}] Loading saved solution for {name}", flush=True)
        _analyze_group(target, group_records, staging)

    validation = {}
    for record in records:
        validation[str(record["name"])] = _validate_mode_dir(
            staging / "outputs" / str(record["name"]),
            float(record["frequency_thz"]),
            float(record["q"]),
        )
    report["validation"] = validation
    report["created_at"] = datetime.now(timezone.utc).isoformat()
    manifest = staging / "finite_mode_analysis_manifest.json"
    manifest.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    overview = target / "10_overview"
    manifest_target = overview / manifest.name
    if manifest_target.exists():
        raise FileExistsError(manifest_target)
    for record in records:
        os.replace(staging / "outputs" / str(record["name"]), record["output_dir"])
    os.replace(manifest, manifest_target)
    shutil.rmtree(staging)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.target, args.apply), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
