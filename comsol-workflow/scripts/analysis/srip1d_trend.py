#!/usr/bin/env python3
"""Plot strip-1D gamma trends across cladding-shift scan results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.run_main.parameter_config import STRIP_1D_OUTPUT_ROOT  # noqa: E402
from scripts.analysis._gamma_trend_common import (  # noqa: E402
    load_trend_csv,
    parse_frequency_thz,
    save_weight_frequency_plot,
    select_closest_mode_by_weight,
)

DEFAULT_OUT_ROOT = STRIP_1D_OUTPUT_ROOT
TREND_OUTPUT_DIR_NAME = "strip1d_trend"
TREND_OUTPUT_STEM = "strip1d_trend"
SHIFT_FACTORS = [idx / 100.0 for idx in range(10)]
GAMMA_SUBSPACE_P_THRESHOLD = 0.9
REQUIRED_COLUMNS = {"frequency", "gamma_k_weight_fraction", "gamma_subspace_p"}
PLOT_COLUMNS = {"shift_factor", "frequency", "gamma_k_weight_fraction"}


def mesh_auto_size_label(mesh_auto_size: int | float) -> str:
    mesh_value = int(mesh_auto_size) if float(mesh_auto_size).is_integer() else mesh_auto_size
    return f"mesh_auto_size_{mesh_value}"


def analysis_dir_for_mesh_auto_size(mesh_auto_size: int | float, *, out_root: Path = DEFAULT_OUT_ROOT) -> Path:
    del mesh_auto_size
    return out_root / TREND_OUTPUT_DIR_NAME


def output_stem_for_mesh_auto_size(mesh_auto_size: int | float) -> str:
    del mesh_auto_size
    return TREND_OUTPUT_STEM


def numeric_label(value: int | float) -> str:
    numeric = float(value)
    return str(int(numeric)) if numeric.is_integer() else f"{numeric:g}"


def scan_output_stem(mesh_auto_size: int | float, strip_cut_angle: int | float) -> str:
    del mesh_auto_size, strip_cut_angle
    return TREND_OUTPUT_STEM


def legacy_scan_output_stem(mesh_auto_size: int | float, bulk_rotation: int | float) -> str:
    del mesh_auto_size, bulk_rotation
    return TREND_OUTPUT_STEM


def resolve_out_root(
    configured_out_root: Path | None,
    *,
    current_series_dir: Path,
) -> Path:
    return Path(configured_out_root) if configured_out_root is not None else Path(current_series_dir)


def selected_modes_csv_path(output_dir: Path, mesh_auto_size: int | float) -> Path:
    del mesh_auto_size
    return output_dir / f"{TREND_OUTPUT_STEM}.csv"


def select_closest_mode(
    mode_scores: pd.DataFrame,
    *,
    target_frequency: float,
    source: str,
) -> pd.Series:
    return select_closest_mode_by_weight(
        mode_scores,
        target_frequency=target_frequency,
        source=source,
        weight_column="gamma_subspace_p",
        weight_threshold=GAMMA_SUBSPACE_P_THRESHOLD,
        additional_numeric_columns=("gamma_k_weight_fraction",),
    )


def mode_scores_path_for_case(case_dir: Path) -> Path:
    """Prefer layout v2 while retaining read-only legacy compatibility."""
    case_dir = Path(case_dir)
    current = case_dir / "10_overview" / "mode_scores.csv"
    if current.is_file():
        return current
    return case_dir / "mode_scores.csv"


def _trend_row(
    selected: pd.Series,
    *,
    shift_factor: float,
    folder: str,
    target_frequency: float,
    mode_scores_path: Path,
) -> dict[str, object]:
    """Keep the exported row schema identical for current and legacy layouts."""
    return {
        "shift_factor": shift_factor,
        "folder": folder,
        "mode_idx": int(selected["mode_idx"]) if "mode_idx" in selected and pd.notna(selected["mode_idx"]) else None,
        "frequency": float(selected["frequency"]),
        "target_frequency": target_frequency,
        "frequency_distance_to_target": float(selected["frequency_distance_to_target"]),
        "gamma_k_weight_fraction": float(selected["gamma_k_weight_fraction"]),
        "gamma_subspace_p": float(selected["gamma_subspace_p"]),
        "mode_scores_path": str(mode_scores_path),
    }


def build_trend_rows(
    *,
    out_root: Path,
    shift_factors: list[float],
    target_frequency: float,
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    issues: list[str] = []

    for shift_factor in shift_factors:
        folder = f"hexagon_cavity_{shift_factor:.2f}"
        mode_scores_path = mode_scores_path_for_case(
            out_root / folder / "strip_1d"
        )
        if not mode_scores_path.exists():
            issues.append(f"{folder}: missing {mode_scores_path.relative_to(out_root)}")
            continue

        try:
            mode_scores = pd.read_csv(mode_scores_path)
            selected = select_closest_mode(
                mode_scores,
                target_frequency=target_frequency,
                source=folder,
            )
        except Exception as exc:
            issues.append(f"{folder}: {exc}")
            continue

        rows.append(_trend_row(
            selected,
            shift_factor=shift_factor,
            folder=folder,
            target_frequency=target_frequency,
            mode_scores_path=mode_scores_path,
        ))

    return rows, issues


def load_run_summary(path: Path) -> list[dict[str, object]]:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    runs = data.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ValueError(f"{path} must contain a non-empty runs list")

    required = {
        "cladding_inward_shift_factor",
        "mesh_auto_size",
        "out_dir",
    }
    for idx, run in enumerate(runs):
        if not isinstance(run, dict):
            raise ValueError(f"{path} run {idx} is not an object")
        missing = required - set(run)
        if missing:
            raise ValueError(f"{path} run {idx} missing required fields: {sorted(missing)}")
        if "strip_cut_angle_degrees" not in run and "bulk_rotation_degrees" not in run:
            raise ValueError(
                f"{path} run {idx} missing strip_cut_angle_degrees "
                "(or legacy bulk_rotation_degrees)"
            )
    return runs


def build_trend_rows_from_run_summary(
    run_summary_path: Path,
    *,
    target_frequency: float,
) -> tuple[list[dict[str, object]], list[str], dict[str, int | float]]:
    run_summary_path = Path(run_summary_path)
    runs = load_run_summary(run_summary_path)
    mesh_values = {float(run["mesh_auto_size"]) for run in runs}
    current_cut_metadata = ["strip_cut_angle_degrees" in run for run in runs]
    if any(current_cut_metadata) and not all(current_cut_metadata):
        raise ValueError(f"{run_summary_path} mixes current cut and legacy rotation metadata")
    orientation_field = (
        "strip_cut_angle_degrees" if all(current_cut_metadata) else "bulk_rotation_degrees"
    )
    orientation_values = {float(run[orientation_field]) for run in runs}
    if len(mesh_values) != 1:
        raise ValueError(f"{run_summary_path} contains inconsistent mesh_auto_size values")
    if len(orientation_values) != 1:
        raise ValueError(f"{run_summary_path} contains inconsistent {orientation_field} values")

    rows: list[dict[str, object]] = []
    issues: list[str] = []
    for run in runs:
        shift_factor = float(run["cladding_inward_shift_factor"])
        case_dir = Path(str(run["out_dir"]))
        if not case_dir.is_absolute():
            case_dir = run_summary_path.parent / case_dir
        elif not case_dir.exists():
            relocated_case_dir = run_summary_path.parent / case_dir.name
            if relocated_case_dir.exists():
                case_dir = relocated_case_dir
        folder = case_dir.name
        mode_scores_path = mode_scores_path_for_case(case_dir)
        if not mode_scores_path.exists():
            issues.append(f"{folder}: missing {mode_scores_path}")
            continue

        try:
            mode_scores = pd.read_csv(mode_scores_path)
            selected = select_closest_mode(
                mode_scores,
                target_frequency=target_frequency,
                source=folder,
            )
        except Exception as exc:
            issues.append(f"{folder}: {exc}")
            continue

        rows.append(_trend_row(
            selected,
            shift_factor=shift_factor,
            folder=folder,
            target_frequency=target_frequency,
            mode_scores_path=mode_scores_path,
        ))

    mesh_auto_size = next(iter(mesh_values))
    orientation_value = next(iter(orientation_values))
    scan_metadata = {
        "mesh_auto_size": int(mesh_auto_size) if mesh_auto_size.is_integer() else mesh_auto_size,
        orientation_field: (
            int(orientation_value) if orientation_value.is_integer() else orientation_value
        ),
    }
    return rows, issues, scan_metadata


def save_trend_plot(selected_modes: pd.DataFrame, path: Path, *, target_frequency: float) -> None:
    save_weight_frequency_plot(
        selected_modes,
        path,
        target_frequency=target_frequency,
        weight_column="gamma_k_weight_fraction",
        title="Strip Gamma Weight And Frequency Closest To Target",
        target_label=f"EIGENFREQUENCY_SHIFT = {target_frequency:g} THz",
    )


def load_selected_modes_csv(path: Path, *, target_frequency: float | None = None) -> tuple[pd.DataFrame, float]:
    return load_trend_csv(
        path,
        required_plot_columns=PLOT_COLUMNS,
        target_frequency=target_frequency,
        allowed_shift_factors=SHIFT_FACTORS,
    )


def plot_from_selected_modes_csv(
    selected_modes_csv: Path,
    *,
    output_dir: Path | None = None,
    target_frequency: float | None = None,
) -> Path:
    selected_modes_csv = Path(selected_modes_csv)
    selected_modes, target_frequency = load_selected_modes_csv(
        selected_modes_csv,
        target_frequency=target_frequency,
    )
    output_dir = output_dir or selected_modes_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_png = output_dir / f"{TREND_OUTPUT_STEM}.png"
    save_trend_plot(selected_modes, plot_png, target_frequency=target_frequency)
    return plot_png


def write_selected_modes_csv(
    *,
    rows: list[dict[str, object]],
    output_dir: Path,
    mesh_auto_size: int | float,
) -> Path:
    if not rows:
        raise ValueError("No selected mode rows to write.")

    output_dir.mkdir(parents=True, exist_ok=True)
    selected_modes = pd.DataFrame(rows).sort_values("shift_factor")
    selected_csv = selected_modes_csv_path(output_dir, mesh_auto_size)

    selected_modes.to_csv(selected_csv, index=False)
    return selected_csv


def write_outputs(
    *,
    rows: list[dict[str, object]],
    issues: list[str],
    output_dir: Path,
    mesh_auto_size: int | float,
) -> tuple[Path, list[str], Path]:
    selected_csv = write_selected_modes_csv(
        rows=rows,
        output_dir=output_dir,
        mesh_auto_size=mesh_auto_size,
    )
    plot_png = plot_from_selected_modes_csv(selected_csv, output_dir=output_dir)
    return selected_csv, issues, plot_png


def generate_trend_plot(
    *,
    out_root: Path,
    output_dir: Path,
    mesh_auto_size: int | float,
    build_target_frequency: float,
    plot_target_frequency: float | None = None,
    shift_factors: list[float] = SHIFT_FACTORS,
) -> tuple[Path, list[str], Path, bool]:
    selected_csv = selected_modes_csv_path(output_dir, mesh_auto_size)
    used_existing_csv = selected_csv.exists()
    issues: list[str] = []

    if not used_existing_csv:
        rows, issues = build_trend_rows(
            out_root=out_root,
            shift_factors=shift_factors,
            target_frequency=build_target_frequency,
        )
        if not rows:
            issue_text = "; ".join(issues)
            raise ValueError(f"No usable mode_scores.csv rows were found. {issue_text}".strip())
        selected_csv = write_selected_modes_csv(
            rows=rows,
            output_dir=output_dir,
            mesh_auto_size=mesh_auto_size,
        )

    plot_png = plot_from_selected_modes_csv(
        selected_csv,
        output_dir=output_dir,
        target_frequency=plot_target_frequency,
    )
    return selected_csv, issues, plot_png, used_existing_csv


def generate_trend_from_run_summary(
    *,
    run_summary_path: Path,
    output_dir: Path,
    target_frequency: float,
) -> tuple[Path, list[str], Path]:
    rows, issues, scan_metadata = build_trend_rows_from_run_summary(
        run_summary_path,
        target_frequency=target_frequency,
    )
    if not rows:
        issue_text = "; ".join(issues)
        raise ValueError(f"No usable mode_scores.csv rows were found. {issue_text}".strip())

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if "strip_cut_angle_degrees" in scan_metadata:
        output_stem = scan_output_stem(
            scan_metadata["mesh_auto_size"],
            scan_metadata["strip_cut_angle_degrees"],
        )
    else:
        output_stem = legacy_scan_output_stem(
            scan_metadata["mesh_auto_size"],
            scan_metadata["bulk_rotation_degrees"],
        )
    selected_csv = output_dir / f"{output_stem}.csv"
    selected_modes = pd.DataFrame(rows).sort_values("shift_factor")
    selected_modes.to_csv(selected_csv, index=False)
    plot_png = output_dir / f"{output_stem}.png"
    save_trend_plot(selected_modes, plot_png, target_frequency=target_frequency)
    return selected_csv, issues, plot_png


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--target-frequency", type=float, default=None)
    parser.add_argument(
        "--selected-modes-csv",
        type=Path,
        default=None,
        help="Plot directly from a previously generated strip1d trend CSV file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.selected_modes_csv is not None:
        try:
            plot_png = plot_from_selected_modes_csv(
                args.selected_modes_csv,
                output_dir=args.output_dir,
                target_frequency=args.target_frequency,
            )
        except Exception as exc:
            print(f"Could not plot selected modes CSV: {exc}", file=sys.stderr)
            return 1
        print(f"Selected modes CSV: {args.selected_modes_csv}")
        print(f"Trend plot: {plot_png}")
        return 0

    from scripts.run_main import run_strip_1d

    target_frequency = (
        args.target_frequency
        if args.target_frequency is not None
        else parse_frequency_thz(run_strip_1d.EIGENFREQUENCY_SHIFT)
    )
    out_root = resolve_out_root(
        args.out_root,
        current_series_dir=run_strip_1d.strip_cut_output_dir(),
    )
    output_dir = args.output_dir or out_root / TREND_OUTPUT_DIR_NAME
    try:
        selected_csv, issues, plot_png = generate_trend_from_run_summary(
            run_summary_path=out_root / "run_summary.json",
            output_dir=output_dir,
            target_frequency=target_frequency,
        )
    except Exception as exc:
        print(f"Could not generate trend plot: {exc}", file=sys.stderr)
        return 1

    print(f"Target frequency: {target_frequency:g} THz")
    print(f"Selected modes: {selected_csv}")
    print(f"Trend plot: {plot_png}")
    print("Selected modes CSV source: generated from run_summary.json")
    for issue in issues:
        print(f"Warning: {issue}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
