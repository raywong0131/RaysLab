#!/usr/bin/env python3
"""Merge legacy and all-symmetry finite-quarter result tables without COMSOL."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.run_main.run_finite_quarter_all_symmetries import save_frequency_q_plot


CASE_NAME = "shiftx0.000_shifty0.000"


def _mode_uid(row: pd.Series) -> str:
    return (
        f"mode{int(row['mode_idx'])}_{int(row['symmetry_id'])}"
        f"_x{row['x_boundary']}_y{row['y_boundary']}"
    )


def _annotate(frame: pd.DataFrame, run: Path, mode_dirs: list[Path]) -> pd.DataFrame:
    result = frame.copy()
    source_uids = result["mode_uid"].astype(str) if "mode_uid" in result else result.apply(_mode_uid, axis=1)
    result.insert(0, "source_mode_dir", [str(path) for path in mode_dirs])
    result.insert(0, "source_mode_uid", source_uids)
    result.insert(0, "source_run", run.name)
    result.insert(0, "record_uid", [f"{run.name}::{uid}" for uid in source_uids])
    return result


def load_legacy(run: Path) -> pd.DataFrame:
    case = run / CASE_NAME
    frames: list[pd.DataFrame] = []
    for symmetry_dir in sorted(case.glob("symmetry_*")):
        csv_path = symmetry_dir / "10_overview" / "eigenfrequencies.csv"
        if not csv_path.is_file():
            raise FileNotFoundError(csv_path)
        frame = pd.read_csv(csv_path)
        mode_dirs = [symmetry_dir / "01_results" / f"mode{int(idx)}" for idx in frame["mode_idx"]]
        frames.append(_annotate(frame, run, mode_dirs))
    if len(frames) != 4:
        raise ValueError(f"Expected four legacy symmetry cases in {case}, found {len(frames)}")
    return pd.concat(frames, ignore_index=True, sort=False)


def load_all_symmetry(run: Path) -> pd.DataFrame:
    case = run / CASE_NAME
    csv_path = case / "10_overview" / "eigenfrequencies.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    frame = pd.read_csv(csv_path)
    if "mode_uid" not in frame:
        raise ValueError(f"Missing mode_uid in {csv_path}")
    mode_dirs = [case / "01_results" / str(uid) for uid in frame["mode_uid"]]
    return _annotate(frame, run, mode_dirs)


def load_target_dataset(run: Path) -> pd.DataFrame:
    csv_path = run / "overview" / "eigenfrequencies.csv"
    if not csv_path.is_file():
        return load_legacy(run)
    frame = pd.read_csv(csv_path)
    required = {"record_uid", "source_run", "source_mode_uid", "source_mode_dir"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Existing merged table is missing columns: {sorted(missing)}")
    return frame


def merge(target_run: Path, additional_run: Path, q_threshold: float = 1000.0) -> dict[str, object]:
    target_run = target_run.resolve(strict=True)
    additional_run = additional_run.resolve(strict=True)
    overview = target_run / "overview"
    if not overview.is_dir():
        raise FileNotFoundError(overview)

    frames = [load_target_dataset(target_run), load_all_symmetry(additional_run)]
    merged = pd.concat(frames, ignore_index=True, sort=False)
    if merged["record_uid"].duplicated().any():
        raise ValueError("Merged record_uid values are not unique")
    source_dirs = [Path(path) for path in merged["source_mode_dir"]]
    missing = [path for path in source_dirs if not path.is_dir()]
    if missing:
        raise FileNotFoundError(f"Missing source mode directory: {missing[0]}")
    merged = merged.sort_values(["re", "record_uid"], kind="stable").reset_index(drop=True)

    q_values = pd.to_numeric(merged["q"], errors="coerce")
    high_q = merged.loc[np.isfinite(q_values) & (q_values > q_threshold)].copy()

    staging_parent = Path(__file__).resolve().parents[2] / 'tmp'
    staging_parent.mkdir(exist_ok=True)
    if staging_parent.drive.lower() != target_run.drive.lower():
        raise ValueError('Merge staging and target run must be on the same drive')
    staging = Path(tempfile.mkdtemp(prefix=f'fqmerge_{uuid.uuid4().hex[:8]}_', dir=staging_parent))
    try:
        merged_csv = staging / "eigenfrequencies.csv"
        merged.to_csv(merged_csv, index=False)
        plot_summary = save_frequency_q_plot(merged_csv, staging / "f_Q.png")
        plot_summary['output_path'] = str(overview / 'f_Q.png')

        staged_images: list[tuple[Path, Path]] = []
        overview_image_count = 0
        new_image_count = 0
        existing_image_count = 0
        used_prefixes: set[str] = set()
        for index, row in high_q.iterrows():
            prefix = (
                f"symmetry_{int(row['symmetry_id'])}"
                f"_x{row['x_boundary']}_y{row['y_boundary']}"
                f"_mode{int(row['mode_idx'])}"
            )
            if prefix in used_prefixes:
                raise ValueError(f"Duplicate high-Q overview prefix: {prefix}")
            used_prefixes.add(prefix)
            high_q.loc[index, "overview_prefix"] = prefix
            source_dir = Path(row["source_mode_dir"]) / "11_simulation_exports"
            for field, source_name in (("hz", "Hz_Im_2d.png"), ("wem", "Wem_2d.png")):
                destination = overview / f"{prefix}_{source_name}"
                source = source_dir / source_name
                if destination.is_file():
                    if destination.stat().st_size == 0:
                        raise IOError(f"Empty existing overview image: {destination}")
                    existing_image_count += 1
                else:
                    if not source.is_file():
                        raise FileNotFoundError(source)
                    staged = staging / destination.name
                    shutil.copy2(source, staged)
                    if staged.stat().st_size != source.stat().st_size:
                        raise IOError(f"Copied image verification failed: {source}")
                    staged_images.append((staged, destination))
                    new_image_count += 1
                overview_image_count += 1
                high_q.loc[index, f"overview_{field}_image"] = str(destination)

        high_q.to_csv(staging / "q_gt_1000_modes.csv", index=False)
        summary = {
            "status": "complete",
            "target_run": str(target_run),
            "additional_run": str(additional_run),
            "source_row_counts": [int(len(frame)) for frame in frames],
            "merged_row_count": int(len(merged)),
            "unique_record_count": int(merged["record_uid"].nunique()),
            "q_threshold_strictly_greater_than": q_threshold,
            "high_q_mode_count": int(len(high_q)),
            "high_q_overview_image_count": overview_image_count,
            "new_overview_image_count": new_image_count,
            "existing_overview_image_count": existing_image_count,
            "frequency_q_plot": plot_summary,
        }
        (staging / "merge_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        for source, destination in staged_images:
            os.replace(source, destination)
        for name in ("eigenfrequencies.csv", "f_Q.png", "q_gt_1000_modes.csv", "merge_summary.json"):
            os.replace(staging / name, overview / name)
        return summary
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_run", type=Path)
    parser.add_argument("additional_run", type=Path)
    parser.add_argument("--q-threshold", type=float, default=1000.0)
    args = parser.parse_args()
    print(json.dumps(merge(args.target_run, args.additional_run, args.q_threshold), indent=2))


if __name__ == "__main__":
    main()
