#!/usr/bin/env python3
"""Migrate completed strip results from layout v2 to layout v3.

The migration loads saved MPH models only to evaluate ``ewfd.Wav`` for the
new x=0 cutline. It does not mesh or solve the model again.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from dataclasses import dataclass, replace
from pathlib import Path

import pandas as pd

from comsol_workflow.simulation_utils import SimulationRun
from comsol_workflow.hex_lattice_utils import LatticePoint
from comsol_workflow.lattice_fourier_postprocess import (
    run_strip_bulk_fourier_postprocess,
    score_strip_modes,
)
from scripts.run_main import run_strip_1d


LAYOUT_VERSION = 3
VALIDITY_THRESHOLD = 0.9


@dataclass(frozen=True)
class CasePlan:
    source: Path
    destination: Path
    shift_name: str
    valid_modes: tuple[int, ...]
    all_modes: tuple[int, ...]
    unscored_modes: tuple[int, ...]


class _SimulationView:
    """Borrow field extraction without owning or clearing a loaded model."""

    get_2d_fields = SimulationRun.get_2d_fields

    def __init__(self, model):
        self.model = model
        self.plane_datasets = {"center": "cpl1"}


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def discover_cases(series_dir: Path, cut_angle: int) -> list[CasePlan]:
    cut_dir = series_dir / f"cut_{cut_angle:02d}"
    if cut_dir.exists():
        raise FileExistsError(f"Migration destination already exists: {cut_dir}")
    pattern = re.compile(rf"cut{cut_angle}_shift(\d+\.\d{{3}})$")
    plans: list[CasePlan] = []
    for source in sorted(series_dir.iterdir(), key=lambda path: path.name):
        if not source.is_dir():
            continue
        match = pattern.fullmatch(source.name)
        if match is None:
            continue
        shift_name = f"shift{match.group(1)}"
        score_path = source / "10_overview" / "mode_scores.csv"
        eigen_path = source / "10_overview" / "eigenfrequencies.csv"
        fourier_path = source / "10_overview" / "strip_bulk_fourier_summary.csv"
        model_path = source / "00_model" / "strip_1d.mph"
        config_path = source / "99_config" / "config.json"
        objective_path = source / "99_config" / "objective.json"
        required = (
            score_path,
            eigen_path,
            fourier_path,
            model_path,
            config_path,
            objective_path,
            source / "01_results",
        )
        missing = [path for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Incomplete legacy strip case: {missing[0]}")

        scores = pd.read_csv(score_path)
        if not {"mode_idx", "gamma_subspace_p"}.issubset(scores.columns):
            raise ValueError(f"Missing strip validity columns in {score_path}")
        eigen = pd.read_csv(eigen_path)
        all_modes = tuple(eigen["mode_idx"].astype(int).tolist())
        scored_modes = set(scores["mode_idx"].astype(int).tolist())
        if not scored_modes.issubset(set(all_modes)):
            raise ValueError(f"mode_scores.csv contains an unknown mode in {source}")
        unscored_modes = tuple(sorted(set(all_modes).difference(scored_modes)))
        valid_modes = tuple(
            scores.loc[
                scores["gamma_subspace_p"].gt(VALIDITY_THRESHOLD),
                "mode_idx",
            ].astype(int).tolist()
        )
        if not valid_modes:
            raise ValueError(f"No gamma-valid mode in {source}")
        existing_modes = {
            int(path.name.removeprefix("mode"))
            for path in (source / "01_results").iterdir()
            if path.is_dir() and re.fullmatch(r"mode\d+", path.name)
        }
        if existing_modes != set(all_modes):
            raise ValueError(
                f"Legacy mode directories differ from mode_scores.csv in {source}"
            )
        for mode_idx in valid_modes:
            exports = source / "01_results" / f"mode{mode_idx}" / "11_simulation_exports"
            required_exports = {
                "Hz_center.parquet",
                "Hz_Re_2d.png",
                "Hz_Im_2d.png",
                "Wem_2d.png",
            }
            missing_exports = sorted(
                required_exports.difference(path.name for path in exports.iterdir())
            )
            if missing_exports:
                raise FileNotFoundError(
                    f"Valid mode {mode_idx} missing {missing_exports[0]} in {source}"
                )
        plans.append(
            CasePlan(
                source=source,
                destination=cut_dir / shift_name,
                shift_name=shift_name,
                valid_modes=valid_modes,
                all_modes=all_modes,
                unscored_modes=unscored_modes,
            )
        )
    if not plans:
        raise FileNotFoundError(
            f"No cut{cut_angle}_shiftX.XXX cases found in {series_dir}"
        )
    return plans


def stage_missing_fourier(plan: CasePlan, staging_dir: Path) -> CasePlan:
    if not plan.unscored_modes:
        return plan
    config = _read_json(plan.source / "99_config" / "config.json")
    backfill_dir = staging_dir / "backfill" / plan.source.name
    export_dir = backfill_dir / "simulation_exports"
    export_dir.mkdir(parents=True)
    shutil.copy2(
        plan.source / "10_overview" / "eigenfrequencies.csv",
        export_dir / "eigenfrequencies.csv",
    )
    for mode_idx in plan.unscored_modes:
        shutil.copy2(
            plan.source
            / "01_results"
            / f"mode{mode_idx}"
            / "11_simulation_exports"
            / "Hz_center.parquet",
            export_dir / f"{mode_idx:02d}_Hz_center.parquet",
        )

    points = [LatticePoint(**row) for row in config["bulk_points"]]
    bulk_records = run_strip_1d.strip_records(
        run_strip_1d.cell_records(points, "bulk")
    )
    angle = float(config["strip_cut_angle_degrees"])
    run_strip_bulk_fourier_postprocess(
        backfill_dir,
        export_dir,
        bulk_records=bulk_records,
        bulk_radius=int(config["bulk_radius"]),
        period=float(config["a"]),
        kx_fraction=float(config["simulation"]["kx"]),
        rho_grid_size=run_strip_1d.STRIP_FOURIER_RHO_GRID_SIZE,
        top_k=run_strip_1d.STRIP_FOURIER_TOP_K,
        mode_indices=list(plan.unscored_modes),
        dpi=run_strip_1d.DPI,
        generate_plots=False,
        row_index_coefficients=run_strip_1d.strip_fourier_row_index_coefficients(angle),
    )
    _, backfill_scores = score_strip_modes(
        backfill_dir,
        export_dirname="simulation_exports",
        gamma_subspace_p_valid_threshold=VALIDITY_THRESHOLD,
    )
    added_valid = set(
        backfill_scores.loc[backfill_scores["is_valid"], "mode_idx"]
        .astype(int)
        .tolist()
    )
    print(
        f"Backfilled {plan.source.name}: modes={list(plan.unscored_modes)}, "
        f"gamma-valid={sorted(added_valid)}",
        flush=True,
    )
    return replace(
        plan,
        valid_modes=tuple(sorted(set(plan.valid_modes).union(added_valid))),
    )


def install_staged_backfill(plan: CasePlan, staging_dir: Path) -> None:
    if not plan.unscored_modes:
        return
    backfill_dir = staging_dir / "backfill" / plan.source.name
    overview = plan.destination / "10_overview"
    scores = pd.concat(
        [
            pd.read_csv(overview / "mode_scores.csv"),
            pd.read_csv(backfill_dir / "mode_scores.csv"),
        ],
        ignore_index=True,
    ).sort_values("mode_idx")
    scores.to_csv(overview / "mode_scores.csv", index=False)
    summary = pd.concat(
        [
            pd.read_csv(overview / "strip_bulk_fourier_summary.csv"),
            pd.read_csv(
                backfill_dir
                / "strip_bulk_fourier_hz"
                / "strip_bulk_fourier_summary.csv"
            ),
        ],
        ignore_index=True,
    ).sort_values("mode_idx")
    summary.to_csv(overview / "strip_bulk_fourier_summary.csv", index=False)
    for mode_idx in plan.unscored_modes:
        if mode_idx not in set(plan.valid_modes):
            continue
        source = backfill_dir / "strip_bulk_fourier_hz" / f"mode_{mode_idx:02d}"
        destination = (
            plan.destination
            / "01_results"
            / f"mode{mode_idx}"
            / "13_strip_bulk_fourier_Hz"
        )
        if destination.exists():
            raise FileExistsError(f"Refusing to overwrite {destination}")
        shutil.copytree(source, destination)


def stage_cutlines(client, plans: list[CasePlan], staging_dir: Path) -> None:
    for case_number, plan in enumerate(plans, start=1):
        model_path = plan.source / "00_model" / "strip_1d.mph"
        print(
            f"Loading model {case_number}/{len(plans)}: {model_path}",
            flush=True,
        )
        model = client.load(str(model_path))
        try:
            sim_run = _SimulationView(model)
            for mode_idx in plan.valid_modes:
                coordinates, values = sim_run.get_2d_fields(
                    mode_idx,
                    "ewfd.Wav",
                    "center",
                )
                output = staging_dir / plan.source.name / f"mode{mode_idx}" / "cutline_Wem.png"
                run_strip_1d.save_wem_cutline_plot(output, coordinates, values)
                print(f"  staged mode{mode_idx}: {output}", flush=True)
        finally:
            client.remove(model)


def _rewrite_case_tables(case_dir: Path, valid_modes: set[int]) -> None:
    overview = case_dir / "10_overview"
    score_path = overview / "mode_scores.csv"
    scores = pd.read_csv(score_path)
    scores["is_valid"] = scores["gamma_subspace_p"].gt(VALIDITY_THRESHOLD)
    scores.to_csv(score_path, index=False)
    gamma_by_mode = scores.set_index("mode_idx")["gamma_subspace_p"]

    eigen_path = overview / "eigenfrequencies.csv"
    eigen = pd.read_csv(eigen_path)
    if "geometry_is_valid" not in eigen.columns:
        eigen["geometry_is_valid"] = eigen["is_valid"]
    eigen["gamma_subspace_p"] = eigen["mode_idx"].map(gamma_by_mode)
    eigen["is_valid"] = eigen["mode_idx"].astype(int).isin(valid_modes)
    eigen.to_csv(eigen_path, index=False)

    fourier_path = overview / "strip_bulk_fourier_summary.csv"
    fourier = pd.read_csv(fourier_path)
    if "geometry_is_valid" not in fourier.columns:
        fourier["geometry_is_valid"] = fourier["is_valid"]
    fourier["gamma_subspace_p"] = fourier["mode_idx"].map(gamma_by_mode)
    fourier["is_valid"] = fourier["mode_idx"].astype(int).isin(valid_modes)
    for row_idx, row in fourier.iterrows():
        mode_idx = int(row["mode_idx"])
        fourier.at[row_idx, "output_dir"] = (
            str(case_dir / "01_results" / f"mode{mode_idx}" / "13_strip_bulk_fourier_Hz")
            if mode_idx in valid_modes
            else ""
        )
    fourier.to_csv(fourier_path, index=False)

    for mode_idx in sorted(valid_modes):
        mode_overview = case_dir / "01_results" / f"mode{mode_idx}" / "10_overview"
        eigen.loc[eigen["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            mode_overview / "eigenfrequency.csv",
            index=False,
        )
        scores.loc[scores["mode_idx"].astype(int).eq(mode_idx)].to_csv(
            mode_overview / "mode_score.csv",
            index=False,
        )


def _rewrite_case_json(case_dir: Path, valid_modes: set[int]) -> None:
    config_path = case_dir / "99_config" / "config.json"
    config = _read_json(config_path)
    config.update(
        {
            "output_layout_version": LAYOUT_VERSION,
            "strip_cut_directory": case_dir.parent.name,
            "mode_validity": {
                "metric": "gamma_subspace_p",
                "operator": ">",
                "threshold": VALIDITY_THRESHOLD,
            },
            "export_dir": "01_results/mode{mode_idx}/11_simulation_exports",
        }
    )
    _write_json(config_path, config)

    score_df = pd.read_csv(case_dir / "10_overview" / "mode_scores.csv")
    eligible = score_df.loc[score_df["mode_idx"].astype(int).isin(valid_modes)]
    best = eligible.loc[eligible["mode_score"].idxmax()].to_dict()
    best = {
        key: value.item() if hasattr(value, "item") else value
        for key, value in best.items()
    }
    score = float(best["mode_score"])
    attrs = {
        "valid_mode_count": len(valid_modes),
        "scored_mode_count": len(score_df),
        "best_mode_idx": int(best["mode_idx"]),
        "best_frequency": float(best["frequency"]),
        "best_q": float(best["q"]),
        "gamma_k_weight_fraction": float(best["gamma_k_weight_fraction"]),
        "gamma_mode_px": float(best["gamma_mode_px"]),
        "score_px_coefficient": float(best["score_px_coefficient"]),
        "raw_gamma_p_px_weight_fraction": float(best["raw_gamma_p_px_weight_fraction"]),
        "gamma_p_px_weight_fraction": score,
        "mode_score": score,
        "validity_metric": "gamma_subspace_p",
        "validity_operator": ">",
        "validity_threshold": VALIDITY_THRESHOLD,
    }
    _write_json(
        case_dir / "99_config" / "objective.json",
        {"score": score, "best_mode": best, "attrs": attrs},
    )


def _rewrite_run_summary(series_dir: Path, cut_dir: Path, plans: list[CasePlan]) -> None:
    source = series_dir / "run_summary.json"
    if source.exists():
        summary = _read_json(source)
    else:
        summary = {"runs": []}
    runs_by_factor = {
        f"{float(row['cladding_inward_shift_factor']):.3f}": row
        for row in summary.get("runs", [])
    }
    updated_runs = []
    for plan in plans:
        factor_key = plan.shift_name.removeprefix("shift")
        row = dict(runs_by_factor.get(factor_key, {}))
        if not row:
            config = _read_json(plan.destination / "99_config" / "config.json")
            row = {
                "cladding_inward_shift_factor": float(factor_key),
                "cladding_inward_shift": float(config["cladding_inward_shift"]),
                "strip_cut_angle_degrees": float(config["strip_cut_angle_degrees"]),
                "bulk_radius": int(config["bulk_radius"]),
                "cladding_layers": int(config["cladding_layers"]),
                "mesh_auto_size": int(config["simulation_common"]["mesh_auto_size"]),
                "eigenmode_count": int(config["simulation_common"]["eigenmode_count"]),
                "case": "strip_1d",
            }
        row.update(
            {
                "output_layout_version": LAYOUT_VERSION,
                "strip_cut_directory": cut_dir.name,
                "mode_validity_metric": "gamma_subspace_p",
                "mode_validity_operator": ">",
                "mode_validity_threshold": VALIDITY_THRESHOLD,
                "out_dir": str(plan.destination),
            }
        )
        updated_runs.append(row)
    _write_json(
        cut_dir / "run_summary.json",
        {"output_layout_version": LAYOUT_VERSION, "runs": updated_runs},
    )
    if source.exists():
        source.unlink()


def _move_trend(series_dir: Path, cut_dir: Path) -> None:
    source = series_dir / "strip1d_trend"
    if not source.exists():
        return
    destination = cut_dir / "strip1d_trend"
    source.rename(destination)
    csv_path = destination / "strip1d_trend.csv"
    if csv_path.exists():
        trend = pd.read_csv(csv_path)
        if "folder" in trend.columns:
            trend["folder"] = trend["folder"].astype(str).str.replace(
                r"^cut\d+_", "", regex=True
            )
        if "mode_scores_path" in trend.columns:
            trend["mode_scores_path"] = trend["mode_scores_path"].astype(str)
            for child in cut_dir.iterdir():
                if child.is_dir() and child.name.startswith("shift"):
                    old_name = f"cut{int(cut_dir.name[-2:])}_{child.name}"
                    trend["mode_scores_path"] = trend["mode_scores_path"].str.replace(
                        old_name,
                        str(Path(cut_dir.name) / child.name),
                        regex=False,
                    )
        trend.to_csv(csv_path, index=False)


def commit_series(
    series_dir: Path,
    plans: list[CasePlan],
    staging_dir: Path,
    cut_angle: int,
) -> None:
    cut_dir = series_dir / f"cut_{cut_angle:02d}"
    cut_dir.mkdir()
    for plan in plans:
        plan.source.rename(plan.destination)
        install_staged_backfill(plan, staging_dir)
        valid_modes = set(plan.valid_modes)
        for mode_idx in plan.valid_modes:
            staged = staging_dir / plan.source.name / f"mode{mode_idx}" / "cutline_Wem.png"
            destination = (
                plan.destination
                / "01_results"
                / f"mode{mode_idx}"
                / "11_simulation_exports"
                / "cutline_Wem.png"
            )
            if destination.exists():
                raise FileExistsError(f"Refusing to overwrite {destination}")
            shutil.move(str(staged), str(destination))
        for mode_idx in set(plan.all_modes).difference(valid_modes):
            mode_dir = plan.destination / "01_results" / f"mode{mode_idx}"
            backup = staging_dir / "removed_modes" / plan.shift_name / mode_dir.name
            backup.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(mode_dir), str(backup))
        _rewrite_case_tables(plan.destination, valid_modes)
        _rewrite_case_json(plan.destination, valid_modes)
        print(f"Migrated {plan.destination}", flush=True)

    _rewrite_run_summary(series_dir, cut_dir, plans)
    _move_trend(series_dir, cut_dir)


def verify_series(series_dir: Path, plans: list[CasePlan], cut_angle: int) -> None:
    cut_dir = series_dir / f"cut_{cut_angle:02d}"
    summary = _read_json(cut_dir / "run_summary.json")
    if summary.get("output_layout_version") != LAYOUT_VERSION:
        raise ValueError(f"Wrong run summary layout in {cut_dir}")
    for plan in plans:
        case_dir = cut_dir / plan.shift_name
        actual_modes = {
            int(path.name.removeprefix("mode"))
            for path in (case_dir / "01_results").iterdir()
            if path.is_dir()
        }
        if actual_modes != set(plan.valid_modes):
            raise ValueError(f"Wrong retained modes in {case_dir}: {actual_modes}")
        for mode_idx in plan.valid_modes:
            cutline = (
                case_dir
                / "01_results"
                / f"mode{mode_idx}"
                / "11_simulation_exports"
                / "cutline_Wem.png"
            )
            if not cutline.is_file() or cutline.stat().st_size == 0:
                raise FileNotFoundError(f"Missing cutline: {cutline}")
        config = _read_json(case_dir / "99_config" / "config.json")
        if config.get("output_layout_version") != LAYOUT_VERSION:
            raise ValueError(f"Wrong config layout in {case_dir}")
    print(f"Verified {len(plans)} cases in {cut_dir}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series", nargs="+", type=Path)
    parser.add_argument("--cut-angle", type=int, default=60)
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    all_plans = {
        series.resolve(): discover_cases(series.resolve(), args.cut_angle)
        for series in args.series
    }
    for series, plans in all_plans.items():
        print(f"{series}: {len(plans)} cases")
        for plan in plans:
            print(
                f"  {plan.source.name} -> {plan.destination.parent.name}/{plan.shift_name}; "
                f"valid={list(plan.valid_modes)}; backfill={list(plan.unscored_modes)}"
            )
    if not args.execute:
        print("Dry run only; pass --execute to migrate.")
        return

    import mph

    client = mph.start(cores=1)
    try:
        for series, plans in all_plans.items():
            staging_dir = series / ".layout_v3_migration"
            if staging_dir.exists():
                raise FileExistsError(f"Migration staging already exists: {staging_dir}")
            staging_dir.mkdir()
            plans = [stage_missing_fourier(plan, staging_dir) for plan in plans]
            stage_cutlines(client, plans, staging_dir)
            commit_series(series, plans, staging_dir, args.cut_angle)
            verify_series(series, plans, args.cut_angle)
            shutil.rmtree(staging_dir)
    finally:
        client.disconnect()


if __name__ == "__main__":
    main()
