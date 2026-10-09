#!/usr/bin/env python3
"""Plot finite-cavity py-mode, frequency, and Q trends across shift results."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.analysis._gamma_trend_common import (  # noqa: E402
    load_trend_csv,
    parse_frequency_thz,
    save_weight_frequency_plot,
    select_closest_mode_by_weight,
)

TREND_OUTPUT_DIR_NAME = "finite_trend"
TREND_OUTPUT_STEM = "finite_trend"
Q_TREND_OUTPUT_STEM = "finite_q_trend"
PY_MODE_WEIGHT_COLUMN = "gamma_p_px_weight_fraction"
PY_MODE_WEIGHT_THRESHOLD = 0.80
REQUIRED_SCORE_COLUMNS = {
    "mode_idx",
    "frequency",
    PY_MODE_WEIGHT_COLUMN,
    "q",
}
PLOT_COLUMNS = {
    "shift_factor",
    "frequency",
    PY_MODE_WEIGHT_COLUMN,
    "q",
}
NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
SHIFT_PATH_PATTERN = re.compile(
    rf"^shift(?P<value>{NUMBER_PATTERN})$"
)
DIRECTIONAL_SHIFT_PATH_PATTERN = re.compile(
    rf"^shiftx(?P<x>{NUMBER_PATTERN})_shifty(?P<y>{NUMBER_PATTERN})$"
)


def resolve_out_root(
    configured_out_root: Path | None,
    *,
    current_series_dir: Path,
) -> Path:
    return (
        Path(configured_out_root)
        if configured_out_root is not None
        else Path(current_series_dir)
    )


def load_run_summary(path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    runs = payload.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ValueError(f"{path} must contain a non-empty runs list")
    for idx, run in enumerate(runs):
        if not isinstance(run, dict):
            raise ValueError(f"{path} run {idx} is not an object")
        if "out_dir" not in run:
            raise ValueError(f"{path} run {idx} missing required field: out_dir")
    return payload, runs


def detect_source_kind(
    payload: dict[str, object],
    runs: list[dict[str, object]],
) -> str:
    if "symmetry_ids" in payload or any("symmetry_id" in run for run in runs):
        return "finite_quarter"
    if any(run.get("case") == "finite_cavity" for run in runs):
        return "finite"
    raise ValueError(
        "Could not identify the series as finite or finite_quarter from run_summary.json"
    )


def _shift_from_path(path: Path) -> float | None:
    for part in reversed(Path(path).parts):
        match = SHIFT_PATH_PATTERN.fullmatch(part)
        if match is not None:
            return float(match.group("value"))
    return None


def _directional_shift_from_path(path: Path) -> tuple[float, float] | None:
    for part in reversed(Path(path).parts):
        match = DIRECTIONAL_SHIFT_PATH_PATTERN.fullmatch(part)
        if match is not None:
            return float(match.group("x")), float(match.group("y"))
    return None


def shift_pair_for_run(run: dict[str, object]) -> tuple[float, float] | None:
    x_value = run.get("cladding_x_shift_factor")
    y_value = run.get("cladding_y_shift_factor")
    if (x_value is None) != (y_value is None):
        raise ValueError(
            "cladding_x_shift_factor and cladding_y_shift_factor must be present together"
        )
    explicit_pair = (
        None
        if x_value is None
        else (float(x_value), float(y_value))
    )
    path_pair = _directional_shift_from_path(Path(str(run["out_dir"])))
    if explicit_pair is not None and path_pair is not None:
        if not all(
            math.isclose(explicit, named, rel_tol=0.0, abs_tol=1e-12)
            for explicit, named in zip(explicit_pair, path_pair)
        ):
            raise ValueError(
                "cladding x/y shift factors disagree with out_dir: "
                f"({explicit_pair[0]:g}, {explicit_pair[1]:g}) != "
                f"({path_pair[0]:g}, {path_pair[1]:g})"
            )
    pair = explicit_pair if explicit_pair is not None else path_pair
    if pair is not None and not all(math.isfinite(value) for value in pair):
        raise ValueError("cladding x/y shift factors must be finite")
    return pair


def shift_factor_for_run(run: dict[str, object]) -> float:
    pair = shift_pair_for_run(run)
    base_value = run.get("cladding_shift_base_factor")
    base = None if base_value is None else float(base_value)
    if base is not None:
        if not math.isfinite(base):
            raise ValueError("cladding_shift_base_factor must be finite")
        if pair is not None and not math.isclose(
            base, pair[0], rel_tol=0.0, abs_tol=1e-12
        ):
            raise ValueError(
                "cladding_shift_base_factor disagrees with x shift: "
                f"{base:g} != {pair[0]:g}"
            )
        return base
    if pair is not None:
        return pair[0]

    explicit_value = run.get("cladding_inward_shift_factor")
    explicit = None if explicit_value is None else float(explicit_value)
    path_value = _shift_from_path(Path(str(run["out_dir"])))
    if explicit is not None and path_value is not None:
        if not math.isclose(explicit, path_value, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                "cladding_inward_shift_factor disagrees with out_dir: "
                f"{explicit:g} != {path_value:g}"
            )
    if explicit is not None:
        return explicit
    if path_value is not None:
        return path_value
    raise ValueError(
        "missing scalar or directional cladding shift metadata/path"
    )


def _shift_identity_for_run(
    run: dict[str, object],
    symmetry_id: int | None,
) -> tuple[float, float | None, int | None]:
    pair = shift_pair_for_run(run)
    if pair is not None:
        return round(pair[0], 12), round(pair[1], 12), symmetry_id
    return round(shift_factor_for_run(run), 12), None, symmetry_id


def _shift_pair_matches_ratio(
    shift_pair: tuple[float, float] | None,
    y_over_x_ratio: float | None,
) -> bool:
    if y_over_x_ratio is None:
        return True
    if shift_pair is None:
        return False
    x_factor, y_factor = shift_pair
    if math.isclose(x_factor, 0.0, abs_tol=1e-12):
        return math.isclose(y_factor, 0.0, abs_tol=1e-12)
    return math.isclose(
        y_factor / x_factor,
        y_over_x_ratio,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _series_relative_tail(path: Path) -> Path | None:
    parts = Path(path).parts
    for idx, part in enumerate(parts):
        if (
            SHIFT_PATH_PATTERN.fullmatch(part)
            or DIRECTIONAL_SHIFT_PATH_PATTERN.fullmatch(part)
        ):
            return Path(*parts[idx:])
    return None


def resolve_case_dir(run: dict[str, object], run_summary_path: Path) -> Path:
    run_summary_path = Path(run_summary_path)
    declared = Path(str(run["out_dir"]))
    if not declared.is_absolute():
        relative_candidate = run_summary_path.parent / declared
        if relative_candidate.exists():
            return relative_candidate
    elif declared.exists():
        return declared

    tail = _series_relative_tail(declared)
    if tail is not None:
        relocated = run_summary_path.parent / tail
        if relocated.exists():
            return relocated

    sibling = run_summary_path.parent / declared.name
    if sibling.exists():
        return sibling
    return run_summary_path.parent / declared if not declared.is_absolute() else declared


def resolve_mode_scores_path(
    run: dict[str, object],
    *,
    case_dir: Path,
    run_summary_path: Path,
) -> Path:
    candidates: list[Path] = []
    overview_value = run.get("overview_dir")
    if overview_value is not None:
        overview_dir = Path(str(overview_value))
        if not overview_dir.is_absolute():
            overview_dir = Path(run_summary_path).parent / overview_dir
        candidates.append(overview_dir / "mode_scores.csv")
    candidates.extend(
        [
            Path(case_dir) / "10_overview" / "mode_scores.csv",
            Path(case_dir) / "mode_scores.csv",
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[-2]


def resolve_target_frequency(
    payload: dict[str, object],
    *,
    target_frequency: float | None,
    fallback_target_frequency: float | str | None,
) -> float:
    if target_frequency is not None:
        resolved = float(target_frequency)
    else:
        preset = payload.get("finite_cavity_preset_config")
        summary_value = (
            preset.get("eigenfrequency_shift")
            if isinstance(preset, dict)
            else None
        )
        if summary_value is not None:
            resolved = parse_frequency_thz(str(summary_value))
        elif fallback_target_frequency is not None:
            resolved = (
                float(fallback_target_frequency)
                if isinstance(fallback_target_frequency, (int, float))
                else parse_frequency_thz(str(fallback_target_frequency))
            )
        else:
            raise ValueError(
                "No target frequency in run_summary.json; pass --target-frequency"
            )
    if not math.isfinite(resolved):
        raise ValueError("target frequency must be finite")
    return resolved


def select_closest_mode(
    mode_scores: pd.DataFrame,
    *,
    target_frequency: float,
    source: str,
) -> pd.Series:
    missing = REQUIRED_SCORE_COLUMNS - set(mode_scores.columns)
    if missing:
        raise ValueError(f"{source} missing required columns: {sorted(missing)}")
    return select_closest_mode_by_weight(
        mode_scores,
        target_frequency=target_frequency,
        source=source,
        weight_column=PY_MODE_WEIGHT_COLUMN,
        weight_threshold=PY_MODE_WEIGHT_THRESHOLD,
    )


def _available_symmetry_ids(runs: list[dict[str, object]]) -> set[int]:
    return {
        int(run["symmetry_id"])
        for run in runs
        if run.get("symmetry_id") is not None
    }


def _validate_symmetry_selection(
    *,
    source_kind: str,
    runs: list[dict[str, object]],
    symmetry_id: int | None,
) -> int | None:
    available = _available_symmetry_ids(runs)
    if source_kind == "finite":
        if symmetry_id is not None:
            raise ValueError("--symmetry-id is only valid for finite_quarter series")
        return None
    if symmetry_id is not None:
        if symmetry_id not in available:
            raise ValueError(
                f"symmetry_id {symmetry_id} is not present; available IDs: "
                f"{sorted(available)}"
            )
        return int(symmetry_id)
    if len(available) > 1:
        raise ValueError(
            "run_summary.json contains multiple symmetry IDs; pass --symmetry-id"
        )
    return next(iter(available)) if available else None


def discover_series_runs(
    run_summary_path: Path,
    *,
    source_kind: str,
    runs: list[dict[str, object]],
    symmetry_id: int | None,
) -> list[dict[str, object]]:
    """Merge summary runs with every shift case physically present in the series."""
    run_summary_path = Path(run_summary_path)
    merged_runs = list(runs)
    known_keys: set[tuple[float, float | None, int | None]] = set()
    for run in runs:
        try:
            key = _shift_identity_for_run(
                run,
                int(run["symmetry_id"])
                if run.get("symmetry_id") is not None
                else None,
            )
        except (TypeError, ValueError):
            continue
        known_keys.add(key)

    for shift_dir in sorted(run_summary_path.parent.iterdir(), key=lambda path: path.name):
        if not shift_dir.is_dir():
            continue
        scalar_match = SHIFT_PATH_PATTERN.fullmatch(shift_dir.name)
        directional_match = DIRECTIONAL_SHIFT_PATH_PATTERN.fullmatch(shift_dir.name)
        if scalar_match is None and directional_match is None:
            continue
        shift_pair = (
            None
            if directional_match is None
            else (
                float(directional_match.group("x")),
                float(directional_match.group("y")),
            )
        )
        shift_factor = (
            float(scalar_match.group("value"))
            if scalar_match is not None
            else shift_pair[0]
        )

        discovered_symmetry = None
        case_dir = shift_dir
        if source_kind == "finite_quarter":
            discovered_symmetry = symmetry_id
            if symmetry_id is not None:
                nested = sorted(
                    path
                    for path in shift_dir.glob(f"symmetry_{symmetry_id}_*")
                    if path.is_dir()
                )
                if len(nested) > 1:
                    raise ValueError(
                        f"{shift_dir} has multiple directories for symmetry_id "
                        f"{symmetry_id}"
                    )
                if nested:
                    case_dir = nested[0]

        key = (
            round(shift_factor, 12),
            None if shift_pair is None else round(shift_pair[1], 12),
            discovered_symmetry,
        )
        if key in known_keys:
            continue
        discovered_run: dict[str, object] = {
            "cladding_inward_shift_factor": shift_factor,
            "symmetry_id": discovered_symmetry,
            "out_dir": str(case_dir),
            "overview_dir": str(case_dir / "10_overview"),
            "case": source_kind,
            "discovered_from_series": True,
        }
        if shift_pair is not None:
            discovered_run.pop("cladding_inward_shift_factor")
            discovered_run.update({
                "cladding_x_shift_factor": shift_pair[0],
                "cladding_y_shift_factor": shift_pair[1],
            })
        merged_runs.append(discovered_run)
        known_keys.add(key)
    return merged_runs


def build_trend_rows_from_run_summary(
    run_summary_path: Path,
    *,
    target_frequency: float,
    symmetry_id: int | None = None,
    y_over_x_ratio: float | None = None,
) -> tuple[list[dict[str, object]], list[str], dict[str, object]]:
    run_summary_path = Path(run_summary_path)
    payload, runs = load_run_summary(run_summary_path)
    source_kind = detect_source_kind(payload, runs)
    selected_symmetry = _validate_symmetry_selection(
        source_kind=source_kind,
        runs=runs,
        symmetry_id=symmetry_id,
    )
    analysis_runs = discover_series_runs(
        run_summary_path,
        source_kind=source_kind,
        runs=runs,
        symmetry_id=selected_symmetry,
    )
    if y_over_x_ratio is not None:
        y_over_x_ratio = float(y_over_x_ratio)
        if not math.isfinite(y_over_x_ratio) or y_over_x_ratio < 0.0:
            raise ValueError("y_over_x_ratio must be finite and non-negative")

    rows: list[dict[str, object]] = []
    issues: list[str] = []
    seen_keys: set[tuple[float, float | None, int | None]] = set()
    for run in analysis_runs:
        run_symmetry = (
            int(run["symmetry_id"])
            if run.get("symmetry_id") is not None
            else None
        )
        if source_kind == "finite_quarter" and run_symmetry != selected_symmetry:
            continue
        try:
            shift_factor = shift_factor_for_run(run)
            shift_pair = shift_pair_for_run(run)
        except Exception as exc:
            issues.append(f"{run.get('out_dir')}: {exc}")
            continue
        if not _shift_pair_matches_ratio(shift_pair, y_over_x_ratio):
            continue

        key = _shift_identity_for_run(run, run_symmetry)
        if key in seen_keys:
            raise ValueError(
                "run_summary.json has duplicate runs for "
                f"shift={key[:2]}, symmetry_id={run_symmetry}"
            )
        seen_keys.add(key)

        case_dir = resolve_case_dir(run, run_summary_path)
        mode_scores_path = resolve_mode_scores_path(
            run,
            case_dir=case_dir,
            run_summary_path=run_summary_path,
        )
        if not mode_scores_path.is_file():
            issues.append(f"{case_dir.name}: missing {mode_scores_path}")
            continue

        try:
            mode_scores = pd.read_csv(mode_scores_path)
            selected = select_closest_mode(
                mode_scores,
                target_frequency=target_frequency,
                source=case_dir.name,
            )
            q_value = float(selected["q"])
            if not math.isfinite(q_value):
                raise ValueError("selected mode has non-finite q")
        except Exception as exc:
            issues.append(f"{case_dir.name}: {exc}")
            continue

        row: dict[str, object] = {
            "shift_factor": shift_factor,
            "shift_kind": "directional" if shift_pair is not None else "scalar",
            "cladding_x_shift_factor": (
                shift_pair[0] if shift_pair is not None else None
            ),
            "cladding_y_shift_factor": (
                shift_pair[1] if shift_pair is not None else None
            ),
            "cladding_y_over_x_shift_ratio": (
                round(shift_pair[1] / shift_pair[0], 12)
                if shift_pair is not None and shift_pair[0] != 0.0
                else None
            ),
            "folder": case_dir.name,
            "source_kind": source_kind,
            "symmetry_id": run_symmetry,
            "mode_idx": int(selected["mode_idx"]),
            "frequency": float(selected["frequency"]),
            "target_frequency": float(target_frequency),
            "frequency_distance_to_target": float(
                selected["frequency_distance_to_target"]
            ),
            PY_MODE_WEIGHT_COLUMN: float(selected[PY_MODE_WEIGHT_COLUMN]),
            "q": q_value,
            "mode_scores_path": str(mode_scores_path),
        }
        for audit_column in ("gamma_k_weight_fraction", "gamma_subspace_p"):
            if audit_column in selected and pd.notna(selected[audit_column]):
                row[audit_column] = float(selected[audit_column])
        rows.append(row)

    metadata = {
        "source_kind": source_kind,
        "symmetry_id": selected_symmetry,
    }
    if y_over_x_ratio is not None:
        metadata["cladding_y_over_x_shift_ratio"] = y_over_x_ratio
    return rows, issues, metadata


def save_trend_plot(
    selected_modes: pd.DataFrame,
    path: Path,
    *,
    target_frequency: float,
    source_kind: str = "finite",
) -> None:
    source_label = "Finite Quarter" if source_kind == "finite_quarter" else "Finite"
    directional_groups = _directional_plot_groups(selected_modes)
    if directional_groups is not None:
        _save_directional_weight_frequency_plot(
            directional_groups,
            path,
            target_frequency=target_frequency,
            source_label=source_label,
        )
        return
    save_weight_frequency_plot(
        selected_modes,
        path,
        target_frequency=target_frequency,
        weight_column=PY_MODE_WEIGHT_COLUMN,
        title=f"{source_label} Py-Mode Weight And Frequency Closest To Target",
        target_label=f"TARGET FREQUENCY = {target_frequency:g} THz",
    )


def _directional_plot_groups(
    selected_modes: pd.DataFrame,
) -> list[tuple[str, pd.DataFrame]] | None:
    directional_columns = {
        "cladding_x_shift_factor",
        "cladding_y_shift_factor",
    }
    if not directional_columns.issubset(selected_modes.columns):
        return None

    directional = selected_modes.copy()
    for column in sorted(directional_columns):
        directional[column] = pd.to_numeric(directional[column], errors="coerce")
    if directional[list(directional_columns)].isna().any(axis=None):
        return None

    x_values = directional["cladding_x_shift_factor"].to_numpy(dtype=float)
    y_values = directional["cladding_y_shift_factor"].to_numpy(dtype=float)
    origin_mask = np.isclose(x_values, 0.0) & np.isclose(y_values, 0.0)
    non_origin = directional.loc[~origin_mask].copy()
    if non_origin.empty:
        return [("x=0, y=0", directional)]

    group_keys: list[tuple[str, float | None]] = []
    for x_factor, y_factor in non_origin[
        ["cladding_x_shift_factor", "cladding_y_shift_factor"]
    ].itertuples(index=False, name=None):
        key = (
            ("ratio", round(float(y_factor) / float(x_factor), 12))
            if not math.isclose(float(x_factor), 0.0, abs_tol=1e-12)
            else ("y_only", None)
        )
        if key not in group_keys:
            group_keys.append(key)

    origin = directional.loc[origin_mask].head(1)
    groups: list[tuple[str, pd.DataFrame]] = []
    for kind, ratio in sorted(
        group_keys,
        key=lambda item: (item[0] == "y_only", float("inf") if item[1] is None else item[1]),
    ):
        if kind == "ratio":
            ratios = (
                non_origin["cladding_y_shift_factor"]
                / non_origin["cladding_x_shift_factor"]
            )
            group = non_origin.loc[np.isclose(ratios, ratio)].copy()
            label = f"y/x={ratio:g}"
        else:
            group = non_origin.loc[
                np.isclose(non_origin["cladding_x_shift_factor"], 0.0)
            ].copy()
            label = "x=0"
        group = pd.concat([origin, group], ignore_index=True)
        group = group.sort_values(
            ["cladding_x_shift_factor", "cladding_y_shift_factor"]
        )
        groups.append((label, group))
    return groups


def _set_directional_x_axis(ax: plt.Axes, groups: list[tuple[str, pd.DataFrame]]) -> None:
    x_values = sorted({
        float(value)
        for _label, group in groups
        for value in group["cladding_x_shift_factor"]
    })
    ax.set_xlabel("CLADDING X SHIFT / A")
    ax.set_xticks(x_values)
    if len(x_values) == 1:
        ax.set_xlim(x_values[0] - 0.005, x_values[0] + 0.005)
    else:
        padding = max(0.005, 0.05 * (x_values[-1] - x_values[0]))
        ax.set_xlim(x_values[0] - padding, x_values[-1] + padding)


def _save_directional_weight_frequency_plot(
    groups: list[tuple[str, pd.DataFrame]],
    path: Path,
    *,
    target_frequency: float,
    source_label: str,
) -> None:
    line_styles = ["-", "--", "-.", ":"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    lines = []
    for idx, (group_label, group) in enumerate(groups):
        style = line_styles[idx % len(line_styles)]
        weight_line, = ax.plot(
            group["cladding_x_shift_factor"],
            group[PY_MODE_WEIGHT_COLUMN],
            marker="o",
            linestyle=style,
            linewidth=1.6,
            markersize=5.5,
            color="#2563eb",
            label=f"{PY_MODE_WEIGHT_COLUMN} ({group_label})",
        )
        lines.append(weight_line)
    ax.set_ylabel(PY_MODE_WEIGHT_COLUMN)
    ax.tick_params(axis="y", labelcolor="#2563eb")
    ax.set_title(
        f"{source_label} Py-Mode Weight And Frequency Closest To Target"
    )
    _set_directional_x_axis(ax, groups)
    ax.grid(True, alpha=0.28, linestyle="--", linewidth=0.7)

    freq_ax = ax.twinx()
    for idx, (group_label, group) in enumerate(groups):
        style = line_styles[idx % len(line_styles)]
        frequency_line, = freq_ax.plot(
            group["cladding_x_shift_factor"],
            group["frequency"],
            marker="s",
            linestyle=style,
            linewidth=1.5,
            markersize=5.0,
            color="#dc2626",
            label=f"frequency ({group_label})",
        )
        lines.append(frequency_line)
    target_line = freq_ax.axhline(
        target_frequency,
        color="#dc2626",
        linestyle=(0, (4, 3)),
        linewidth=1.1,
        alpha=0.8,
        label=f"TARGET FREQUENCY = {target_frequency:g} THz",
    )
    lines.append(target_line)
    freq_ax.set_ylabel("frequency (THz)")
    freq_ax.tick_params(axis="y", labelcolor="#dc2626")
    ax.legend(
        lines,
        [line.get_label() for line in lines],
        loc="best",
        fontsize=7.5,
    )
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_q_trend_plot(
    selected_modes: pd.DataFrame,
    path: Path,
    *,
    source_kind: str = "finite",
) -> None:
    source_label = "Finite Quarter" if source_kind == "finite_quarter" else "Finite"
    directional_groups = _directional_plot_groups(selected_modes)
    if directional_groups is not None:
        line_styles = ["-", "--", "-.", ":"]
        fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
        for idx, (group_label, group) in enumerate(directional_groups):
            ax.plot(
                group["cladding_x_shift_factor"],
                group["q"],
                marker="o",
                linestyle=line_styles[idx % len(line_styles)],
                linewidth=1.6,
                markersize=5.5,
                color="#7c3aed",
                label=f"Q ({group_label})",
            )
        ax.set_ylabel("Q")
        ax.set_title(f"{source_label} Target-Mode Q Versus Cladding Shift")
        _set_directional_x_axis(ax, directional_groups)
        ax.grid(True, alpha=0.28, linestyle="--", linewidth=0.7)
        ax.legend(loc="best", fontsize=8)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return

    shift_factors = sorted(selected_modes["shift_factor"].unique())
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.plot(
        selected_modes["shift_factor"],
        selected_modes["q"],
        marker="o",
        linewidth=1.6,
        markersize=5.5,
        color="#7c3aed",
        label="Q",
    )
    ax.set_xlabel("CLADDING_INWARD_SHIFT / A")
    ax.set_ylabel("Q")
    ax.set_title(f"{source_label} Target-Mode Q Versus Cladding Shift")
    ax.set_xticks(shift_factors)
    ax.set_xlim(min(shift_factors) - 0.005, max(shift_factors) + 0.005)
    ax.grid(True, alpha=0.28, linestyle="--", linewidth=0.7)
    ax.legend(loc="best", fontsize=8)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def load_selected_modes_csv(
    path: Path,
    *,
    target_frequency: float | None = None,
) -> tuple[pd.DataFrame, float, str]:
    selected_modes, resolved_target = load_trend_csv(
        path,
        required_plot_columns=PLOT_COLUMNS,
        target_frequency=target_frequency,
        require_unique_target=True,
    )
    numeric = selected_modes[list(PLOT_COLUMNS)].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError(f"{path} has non-finite values in required plot columns")
    if not (selected_modes[PY_MODE_WEIGHT_COLUMN] > PY_MODE_WEIGHT_THRESHOLD).all():
        raise ValueError(
            f"{path} contains {PY_MODE_WEIGHT_COLUMN} values that are not "
            f"> {PY_MODE_WEIGHT_THRESHOLD:g}"
        )
    source_kind = "finite"
    if "source_kind" in selected_modes.columns:
        source_values = selected_modes["source_kind"].dropna().astype(str).unique()
        if len(source_values) != 1:
            raise ValueError(f"{path} has inconsistent source_kind values")
        source_kind = source_values[0]
    return selected_modes, resolved_target, source_kind


def plot_from_selected_modes_csv(
    selected_modes_csv: Path,
    *,
    output_dir: Path | None = None,
    target_frequency: float | None = None,
) -> tuple[Path, Path]:
    selected_modes_csv = Path(selected_modes_csv)
    selected_modes, resolved_target, source_kind = load_selected_modes_csv(
        selected_modes_csv,
        target_frequency=target_frequency,
    )
    output_dir = Path(output_dir) if output_dir is not None else selected_modes_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    trend_png = output_dir / f"{TREND_OUTPUT_STEM}.png"
    q_trend_png = output_dir / f"{Q_TREND_OUTPUT_STEM}.png"
    save_trend_plot(
        selected_modes,
        trend_png,
        target_frequency=resolved_target,
        source_kind=source_kind,
    )
    save_q_trend_plot(selected_modes, q_trend_png, source_kind=source_kind)
    return trend_png, q_trend_png


def write_selected_modes_csv(
    *,
    rows: list[dict[str, object]],
    output_dir: Path,
) -> Path:
    if not rows:
        raise ValueError("No selected mode rows to write.")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_csv = output_dir / f"{TREND_OUTPUT_STEM}.csv"
    selected = pd.DataFrame(rows)
    sort_columns = ["shift_factor"]
    if "cladding_y_shift_factor" in selected.columns:
        sort_columns.append("cladding_y_shift_factor")
    selected.sort_values(sort_columns).to_csv(selected_csv, index=False)
    return selected_csv


def generate_trend_from_run_summary(
    *,
    run_summary_path: Path,
    output_dir: Path,
    target_frequency: float | None = None,
    fallback_target_frequency: float | str | None = None,
    symmetry_id: int | None = None,
    y_over_x_ratio: float | None = None,
) -> tuple[Path, list[str], Path, Path, float]:
    payload, _runs = load_run_summary(run_summary_path)
    resolved_target = resolve_target_frequency(
        payload,
        target_frequency=target_frequency,
        fallback_target_frequency=fallback_target_frequency,
    )
    rows, issues, _metadata = build_trend_rows_from_run_summary(
        run_summary_path,
        target_frequency=resolved_target,
        symmetry_id=symmetry_id,
        y_over_x_ratio=y_over_x_ratio,
    )
    if not rows:
        issue_text = "; ".join(issues)
        raise ValueError(
            f"No usable mode_scores.csv rows were found. {issue_text}".strip()
        )
    selected_csv = write_selected_modes_csv(rows=rows, output_dir=output_dir)
    trend_png, q_trend_png = plot_from_selected_modes_csv(
        selected_csv,
        output_dir=output_dir,
    )
    return selected_csv, issues, trend_png, q_trend_png, resolved_target


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--target-frequency", type=float, default=None)
    parser.add_argument("--symmetry-id", type=int, default=None)
    parser.add_argument(
        "--y-over-x-ratio",
        type=float,
        default=None,
        help=(
            "For directional cases, keep only this y/x ratio; the (0, 0) "
            "baseline is retained."
        ),
    )
    parser.add_argument(
        "--selected-modes-csv",
        type=Path,
        default=None,
        help="Plot directly from a previously generated finite trend CSV file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.selected_modes_csv is not None:
        try:
            trend_png, q_trend_png = plot_from_selected_modes_csv(
                args.selected_modes_csv,
                output_dir=args.output_dir,
                target_frequency=args.target_frequency,
            )
        except Exception as exc:
            print(f"Could not plot selected modes CSV: {exc}", file=sys.stderr)
            return 1
        print(f"Selected modes CSV: {args.selected_modes_csv}")
        print(f"Finite trend plot: {trend_png}")
        print(f"Finite Q trend plot: {q_trend_png}")
        return 0

    from scripts.run_main import run_finite

    out_root = resolve_out_root(
        args.out_root,
        current_series_dir=run_finite.OUT_DIR,
    )
    output_dir = args.output_dir or out_root / TREND_OUTPUT_DIR_NAME
    try:
        selected_csv, issues, trend_png, q_trend_png, target_frequency = (
            generate_trend_from_run_summary(
                run_summary_path=out_root / "run_summary.json",
                output_dir=output_dir,
                target_frequency=args.target_frequency,
                fallback_target_frequency=run_finite.EIGENFREQUENCY_SHIFT,
                symmetry_id=args.symmetry_id,
                y_over_x_ratio=args.y_over_x_ratio,
            )
        )
    except Exception as exc:
        print(f"Could not generate finite trend plots: {exc}", file=sys.stderr)
        return 1

    print(f"Target frequency: {target_frequency:g} THz")
    print(f"Selected modes: {selected_csv}")
    print(f"Finite trend plot: {trend_png}")
    print(f"Finite Q trend plot: {q_trend_png}")
    print(
        "Selected modes CSV source: generated from run_summary.json "
        "and series shift directories"
    )
    for issue in issues:
        print(f"Warning: {issue}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
