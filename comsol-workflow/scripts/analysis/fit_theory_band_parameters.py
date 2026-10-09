#!/usr/bin/env python3
"""Fit S17 parameters from a completed theory-fit dataset without COMSOL."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.c2v_band_model import (  # noqa: E402
    INTERNAL_TO_THEORY_MODE,
    STANDING_WAVE_LABELS,
    reciprocal_scale_inv_um,
)
from comsol_workflow.theory_band_fitting import (  # noqa: E402
    PARAMETER_NAMES,
    calibrate_c6v,
    calibration_covariance,
    fit_c2v_parameters,
    propagate_calibration_uncertainty,
    tracked_s17_bands,
)


OUTPUT_ROOT = SCRIPT_DIR / ".out" / "unit_cell_band"
DEFAULT_DATASET = OUTPUT_ROOT / "theoryfit_dataset_b245_eta0.94-1.04_target_eta0.96_zeta1.156_mesh5"
DEFAULT_RESULT = OUTPUT_ROOT / "theoryfit_result_b245_eta0.96_zeta1.156_S17"
BAND_COMPARISON_FILENAME = "band_fit_comparison.png"
FULL_BAND_COMPARISON_FILENAME = "full_band_fit_comparison.png"
FULL_BAND_SECOND_ORDER_FILENAME = "full_band_second_order_comparison.png"
BAND_PANEL_TITLES = (r"$\Gamma$–K ($+k_x$)", r"$\Gamma$–M ($+k_y$)")
BAND_MODE_COLORS = {
    "d_xy": "#b8b8b8",
    "d_x2_minus_y2": "#111111",
    "p_y": "#d62728",
    "p_x": "#1f77b4",
}
BAND_MODE_DISPLAY = {
    "d_xy": r"$d_{xy}$",
    "d_x2_minus_y2": r"$d_{x^2-y^2}$",
    "p_y": r"$p_y$",
    "p_x": r"$p_x$",
}


def load_selected(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "direction", "point_index", "band_label", "kx", "ky", "re",
        "is_valid", "match_status", "dominant_mode",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    valid = frame[frame["is_valid"].astype(bool)].copy()
    return valid[valid["match_status"] != "unmatched"]


def gamma_runtime_mapping(frame: pd.DataFrame) -> dict[str, str]:
    gamma = frame[np.isclose(frame["kx"], 0.0) & np.isclose(frame["ky"], 0.0)]
    gamma = gamma.drop_duplicates("band_label")
    if len(gamma) != 4:
        raise ValueError(f"Expected four unique Gamma bands, found {len(gamma)}")
    mapping = {
        str(row.band_label): INTERNAL_TO_THEORY_MODE[str(row.dominant_mode)]
        for row in gamma.itertuples()
    }
    if set(mapping.values()) != set(STANDING_WAVE_LABELS):
        raise ValueError(f"Gamma basis mapping is incomplete or ambiguous: {mapping}")
    return mapping


def labeled_ray(frame: pd.DataFrame, direction: str) -> tuple[np.ndarray, np.ndarray]:
    """Return one Gamma-first ray and modes in canonical standing-wave order."""

    mapping = gamma_runtime_mapping(frame)
    ray = frame[frame["direction"] == direction].copy()
    ray["theory_mode"] = ray["band_label"].map(mapping)
    frequencies = ray.pivot(index="point_index", columns="theory_mode", values="re")
    coordinates = ray.drop_duplicates("point_index").set_index("point_index")[["kx", "ky"]]
    indices = sorted(set(frequencies.index).intersection(coordinates.index))
    frequencies = frequencies.loc[indices, list(STANDING_WAVE_LABELS)]
    coordinates = coordinates.loc[indices]
    if frequencies.isna().any().any():
        raise ValueError(f"Ray {direction} contains incomplete four-band rows")
    return coordinates.to_numpy(float), frequencies.to_numpy(float)


def evaluate_fit_directions(
    frame: pd.DataFrame, *, directions: tuple[str, ...], lattice_constant_um: float,
    calibration: dict[str, float], parameters: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate labeled rays and return point residuals plus split metrics."""

    records = []
    scale = reciprocal_scale_inv_um(lattice_constant_um)
    kwargs = {
        "t0_thz": float(calibration["t0_thz"]),
        "t1_thz": float(calibration["t1_thz"]),
        "v_thz_um": float(calibration["velocity_thz_um"]),
        "alpha": float(parameters["alpha"]),
        "beta": float(parameters["beta"]),
        "gamma": float(parameters["gamma"]),
        "mu_thz": float(parameters["mu_thz"]),
    }
    for direction in directions:
        points, observed = labeled_ray(frame, direction)
        center = float(np.mean(observed[0]))
        observed = observed - center
        predicted, _ = tracked_s17_bands(points * scale, **kwargs)
        split = (
            "fit" if direction.startswith("fit_")
            else "validation" if direction.startswith("validation_") else "qa"
        )
        for point_index, (point, observed_row, predicted_row) in enumerate(
            zip(points, observed, predicted)
        ):
            for mode_index, mode in enumerate(STANDING_WAVE_LABELS):
                records.append({
                    "split": split, "direction": direction,
                    "point_index": point_index, "is_gamma": point_index == 0,
                    "kx_over_G": float(point[0]), "ky_over_G": float(point[1]),
                    "theory_mode": mode,
                    "observed_centered_thz": float(observed_row[mode_index]),
                    "predicted_centered_thz": float(predicted_row[mode_index]),
                    "residual_thz": float(predicted_row[mode_index] - observed_row[mode_index]),
                })
    residuals = pd.DataFrame(records)
    metric_rows = []
    for direction, group in residuals.groupby("direction", sort=False):
        for scope, scoped in (
            ("all", group), ("non_gamma", group[~group["is_gamma"]])
        ):
            values = scoped["residual_thz"].to_numpy(float)
            metric_rows.append({
                "split": str(group["split"].iloc[0]), "direction": direction,
                "scope": scope, "sample_count": int(values.size),
                "rmse_thz": float(np.sqrt(np.mean(values**2))),
                "max_abs_thz": float(np.max(np.abs(values))),
            })
    return residuals, pd.DataFrame(metric_rows)


def calibrate_structure(selected_path: Path, lattice_constant_um: float):
    frame = load_selected(selected_path)
    gamma = frame[np.isclose(frame["kx"], 0.0) & np.isclose(frame["ky"], 0.0)]
    gamma = gamma.drop_duplicates("band_label")
    p_gamma = gamma[gamma["band_label"].astype(str).str.startswith("p")]["re"]
    d_gamma = gamma[gamma["band_label"].astype(str).str.startswith("d")]["re"]
    k_values, lower, upper = [], [], []
    scale = reciprocal_scale_inv_um(lattice_constant_um)
    for direction in ("fit_x", "fit_y"):
        coordinates, frequencies = labeled_ray(frame, direction)
        for point, row in zip(coordinates, frequencies):
            k_values.append(float(np.linalg.norm(point) * scale))
            ordered = np.sort(row)
            lower.append(ordered[:2])
            upper.append(ordered[2:])
    return calibrate_c6v(
        gamma_p_frequencies_thz=p_gamma,
        gamma_d_frequencies_thz=d_gamma,
        k_abs_inv_um=k_values,
        lower_doublet_frequencies_thz=lower,
        upper_doublet_frequencies_thz=upper,
        lattice_constant_um=lattice_constant_um,
    )


def fit_dataset(dataset_dir: Path, result_dir: Path) -> dict[str, object]:
    manifest = json.loads(
        (dataset_dir / "99_config" / "dry_run_manifest.json").read_text(encoding="utf-8")
    )
    if manifest.get("status") != "complete":
        raise RuntimeError("Dataset manifest is not complete; refusing a partial fit")
    config = manifest["fit_config_snapshot"]
    lattice_constant = float(config["lattice_constant_um"])
    calibration_rows = []
    calibration_estimates = {}
    for structure in manifest["structures"]:
        if structure["role"] != "calibration":
            continue
        selected = (
            dataset_dir / structure["relative_output"] /
            "10_overview" / "cavity" / "selected_bands.csv"
        )
        estimate = calibrate_structure(selected, lattice_constant)
        calibration_estimates[float(structure["cell"]["eta"])] = estimate
        calibration_rows.append({"eta": structure["cell"]["eta"], **asdict(estimate)})
    calibration = pd.DataFrame(calibration_rows).sort_values("eta")
    target_eta = float(config["target"]["eta"])
    target_calibration = calibration[np.isclose(calibration["eta"], target_eta)]
    if len(target_calibration) != 1:
        raise ValueError(f"Expected exactly one calibration row at eta={target_eta}")
    fixed = target_calibration.iloc[0]
    target_structure = next(item for item in manifest["structures"] if item["role"] == "target")
    target_path = (
        dataset_dir / target_structure["relative_output"] /
        "10_overview" / "cavity" / "selected_bands.csv"
    )
    target = load_selected(target_path)
    points, modes, slices = [], [], []
    cursor = 0
    for direction in ("fit_x", "fit_y"):
        ray_points, ray_modes = labeled_ray(target, direction)
        points.append(ray_points * reciprocal_scale_inv_um(lattice_constant))
        center = float(np.mean(ray_modes[0]))
        modes.append(ray_modes - center)
        slices.append(slice(cursor, cursor + len(ray_points)))
        cursor += len(ray_points)
    fit_options = config["fit"]
    result = fit_c2v_parameters(
        k_points_inv_um=np.vstack(points), observed_modes_thz=np.vstack(modes),
        t0_thz=float(fixed.t0_thz), t1_thz=float(fixed.t1_thz),
        v_thz_um=float(fixed.velocity_thz_um), ray_slices=tuple(slices),
        lower_bounds=fit_options["lower_bounds"], upper_bounds=fit_options["upper_bounds"],
        multistarts=int(fit_options["multistarts"]),
        random_seed=int(fit_options["random_seed"]), loss=str(fit_options["loss"]),
        frequency_scale_thz=float(fit_options["frequency_scale_thz"]),
        condition_warning=float(fit_options["condition_warning"]),
    )
    calibration_estimate = calibration_estimates[target_eta]
    propagated = propagate_calibration_uncertainty(
        calibration_mean=[fixed.t0_thz, fixed.t1_thz, fixed.velocity_thz_um],
        calibration_covariance_matrix=calibration_covariance(
            calibration_estimate, lattice_constant
        ),
        samples=int(fit_options["calibration_samples"]),
        random_seed=int(fit_options["random_seed"]) + 1,
        fit_arguments={
            "k_points_inv_um": np.vstack(points),
            "observed_modes_thz": np.vstack(modes),
            "ray_slices": tuple(slices),
            "lower_bounds": fit_options["lower_bounds"],
            "upper_bounds": fit_options["upper_bounds"],
            "multistarts": 3,
            "random_seed": int(fit_options["random_seed"]),
            "loss": str(fit_options["loss"]),
            "frequency_scale_thz": float(fit_options["frequency_scale_thz"]),
            "condition_warning": float(fit_options["condition_warning"]),
        },
    )
    result_dir.mkdir(parents=True, exist_ok=False)
    calibration.to_csv(result_dir / "eta_calibration.csv", index=False)
    pd.DataFrame(result.diagnostics.covariance, index=result.parameters, columns=result.parameters).to_csv(
        result_dir / "covariance.csv"
    )
    pd.DataFrame(result.diagnostics.correlation, index=result.parameters, columns=result.parameters).to_csv(
        result_dir / "correlation.csv"
    )
    pd.DataFrame(propagated, columns=PARAMETER_NAMES).to_csv(
        result_dir / "calibration_propagation_samples.csv", index=False
    )
    quantiles = np.quantile(propagated, [0.025, 0.5, 0.975], axis=0)
    summary = {
        "model": config["model"], "parameters": result.parameters,
        "fixed_calibration": fixed.to_dict(),
        "cost": result.cost, "success": result.success, "message": result.message,
        "diagnostics": {
            "singular_values": result.diagnostics.singular_values.tolist(),
            "condition_number": result.diagnostics.condition_number,
            "rank": result.diagnostics.rank,
            "degrees_of_freedom": result.diagnostics.degrees_of_freedom,
            "residual_variance": result.diagnostics.residual_variance,
            "ill_conditioned": result.diagnostics.ill_conditioned,
        },
        "calibration_propagation": {
            name: {
                "q025": float(quantiles[0, index]),
                "median": float(quantiles[1, index]),
                "q975": float(quantiles[2, index]),
            }
            for index, name in enumerate(PARAMETER_NAMES)
        },
    }
    (result_dir / "fit_result.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_additional_diagnostics(dataset_dir, result_dir, summary=summary)
    plot_band_fit_comparison(dataset_dir, result_dir)
    return summary


def write_additional_diagnostics(
    dataset_dir: Path, result_dir: Path, *, summary: dict[str, object] | None = None,
) -> tuple[Path, ...]:
    """Add residual and parameter audit tables without overwriting existing files."""

    targets = tuple(
        result_dir / name for name in (
            "parameters.csv", "residuals.csv", "validation_metrics.csv",
            "jacobian_singular_values.csv",
        )
    )
    existing = [path for path in targets if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite diagnostic files: {existing}")
    if summary is None:
        summary = json.loads((result_dir / "fit_result.json").read_text(encoding="utf-8"))
    manifest = json.loads(
        (dataset_dir / "99_config" / "dry_run_manifest.json").read_text(encoding="utf-8")
    )
    target_structure = next(item for item in manifest["structures"] if item["role"] == "target")
    target = load_selected(
        dataset_dir / target_structure["relative_output"] /
        "10_overview" / "cavity" / "selected_bands.csv"
    )
    residuals, metrics = evaluate_fit_directions(
        target,
        directions=("fit_x", "fit_y", "validation_1", "validation_2", "qa_1", "qa_2"),
        lattice_constant_um=float(manifest["fit_config_snapshot"]["lattice_constant_um"]),
        calibration=summary["fixed_calibration"],
        parameters=summary["parameters"],
    )
    residuals.to_csv(result_dir / "residuals.csv", index=False)
    metrics.to_csv(result_dir / "validation_metrics.csv", index=False)
    singular_values = summary["diagnostics"]["singular_values"]
    pd.DataFrame({
        "index": np.arange(1, len(singular_values) + 1),
        "singular_value": singular_values,
    }).to_csv(result_dir / "jacobian_singular_values.csv", index=False)
    covariance = pd.read_csv(result_dir / "covariance.csv", index_col=0)
    propagation = pd.read_csv(result_dir / "calibration_propagation_samples.csv")
    parameter_rows = []
    for name in PARAMETER_NAMES:
        parameter_rows.append({
            "parameter": name, "estimate": float(summary["parameters"][name]),
            "conditional_std": float(np.sqrt(max(covariance.loc[name, name], 0.0))),
            "propagated_q025": float(propagation[name].quantile(0.025)),
            "propagated_median": float(propagation[name].quantile(0.5)),
            "propagated_q975": float(propagation[name].quantile(0.975)),
        })
    pd.DataFrame(parameter_rows).to_csv(result_dir / "parameters.csv", index=False)
    return targets


def plot_band_fit_comparison(
    dataset_dir: Path, result_dir: Path, *, output_path: Path | None = None,
    overwrite: bool = False,
) -> dict[str, object]:
    """Plot COMSOL and fitted S17 bands along the two fitted positive axes."""

    output = result_dir / BAND_COMPARISON_FILENAME if output_path is None else output_path
    if output.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing figure: {output}")
    manifest = json.loads(
        (dataset_dir / "99_config" / "dry_run_manifest.json").read_text(encoding="utf-8")
    )
    summary = json.loads((result_dir / "fit_result.json").read_text(encoding="utf-8"))
    target_structure = next(item for item in manifest["structures"] if item["role"] == "target")
    target = load_selected(
        dataset_dir / target_structure["relative_output"] /
        "10_overview" / "cavity" / "selected_bands.csv"
    )
    lattice_constant = float(manifest["fit_config_snapshot"]["lattice_constant_um"])
    scale = reciprocal_scale_inv_um(lattice_constant)
    calibration = summary["fixed_calibration"]
    parameters = summary["parameters"]
    model_kwargs = {
        "t0_thz": float(calibration["t0_thz"]),
        "t1_thz": float(calibration["t1_thz"]),
        "v_thz_um": float(calibration["velocity_thz_um"]),
        "alpha": float(parameters["alpha"]),
        "beta": float(parameters["beta"]),
        "gamma": float(parameters["gamma"]),
        "mu_thz": float(parameters["mu_thz"]),
    }
    panel_data = []
    all_frequencies = []
    for direction in ("fit_x", "fit_y"):
        points, observed = labeled_ray(target, direction)
        center = float(np.mean(observed[0]))
        magnitudes = np.linalg.norm(points, axis=1)
        smooth_magnitudes = np.linspace(0.0, float(np.max(magnitudes)), 301)
        smooth_points = (
            np.column_stack([smooth_magnitudes, np.zeros_like(smooth_magnitudes)])
            if direction == "fit_x"
            else np.column_stack([np.zeros_like(smooth_magnitudes), smooth_magnitudes])
        )
        predicted_at_data, _ = tracked_s17_bands(points * scale, **model_kwargs)
        smooth_predicted, _ = tracked_s17_bands(smooth_points * scale, **model_kwargs)
        predicted_absolute = predicted_at_data + center
        smooth_absolute = smooth_predicted + center
        residual = predicted_absolute[1:] - observed[1:]
        rmse = float(np.sqrt(np.mean(residual**2)))
        panel_data.append((magnitudes, observed, smooth_magnitudes, smooth_absolute, rmse))
        all_frequencies.extend([observed.ravel(), smooth_absolute.ravel()])

    frequency_values = np.concatenate(all_frequencies)
    padding = max(0.04 * float(np.ptp(frequency_values)), 0.1)
    y_limits = (
        float(np.min(frequency_values) - padding),
        float(np.max(frequency_values) + padding),
    )
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), sharey=True)
    for axis, title, data in zip(axes, BAND_PANEL_TITLES, panel_data):
        magnitudes, observed, smooth_magnitudes, smooth_absolute, rmse = data
        for mode_index, mode in enumerate(STANDING_WAVE_LABELS):
            color = BAND_MODE_COLORS[mode]
            axis.plot(magnitudes, observed[:, mode_index], color=color, lw=2.5, alpha=0.42)
            axis.plot(
                smooth_magnitudes, smooth_absolute[:, mode_index],
                color=color, lw=1.35, linestyle=(0, (3.0, 2.0)),
            )
        axis.set_title(title, fontsize=10.5, pad=6)
        axis.set_xlabel(r"$k/G$", fontsize=9.5)
        axis.set_xlim(0.0, max(float(np.max(magnitudes)), 1e-12))
        ticks = np.linspace(0.0, float(np.max(magnitudes)), 5)
        axis.set_xticks(ticks)
        axis.set_xticklabels([f"{value:.3f}" for value in ticks])
        axis.set_ylim(*y_limits)
        axis.text(
            0.04, 0.05, f"RMSE = {rmse:.3f} THz", transform=axis.transAxes,
            fontsize=7.8, ha="left", va="bottom",
        )
        axis.tick_params(direction="out", length=3, width=0.8, labelsize=8.5)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[0].set_ylabel("Frequency (THz)", fontsize=9.5)
    source_handles = [
        Line2D([0], [0], color="#333333", lw=2.5, alpha=0.42, label="COMSOL"),
        Line2D([0], [0], color="#333333", lw=1.35, linestyle=(0, (3.0, 2.0)), label="S17 fit"),
    ]
    mode_handles = [
        Line2D([0], [0], color=BAND_MODE_COLORS[mode], lw=2.0, label=BAND_MODE_DISPLAY[mode])
        for mode in STANDING_WAVE_LABELS
    ]
    figure.legend(
        handles=source_handles + mode_handles, loc="upper center", ncol=6,
        frameon=False, fontsize=7.6, handlelength=2.2, columnspacing=1.25,
        bbox_to_anchor=(0.5, 1.02),
    )
    figure.subplots_adjust(left=0.095, right=0.985, bottom=0.17, top=0.80, wspace=0.10)
    figure.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return {
        "output": str(output), "panel_titles": list(BAND_PANEL_TITLES),
        "x_limits_over_G": [0.0, float(max(data[0].max() for data in panel_data))],
        "y_limits_thz": list(y_limits),
        "rmse_thz": [float(data[4]) for data in panel_data],
        "line_count": 16,
    }


def plot_full_band_fit_comparison(
    full_band_selected_path: Path, result_dir: Path, *,
    output_path: Path | None = None, overwrite: bool = False,
    include_second_order: bool = False,
) -> dict[str, object]:
    """Compare COMSOL with first-order S17 and optionally Appendix C5."""

    default_name = (
        FULL_BAND_SECOND_ORDER_FILENAME
        if include_second_order else FULL_BAND_COMPARISON_FILENAME
    )
    output = result_dir / default_name if output_path is None else output_path
    if output.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing figure: {output}")
    frame = load_selected(full_band_selected_path)
    summary = json.loads((result_dir / "fit_result.json").read_text(encoding="utf-8"))
    source_config_path = full_band_selected_path.parents[2] / "99_config" / "config.json"
    source_config = json.loads(source_config_path.read_text(encoding="utf-8"))
    lattice_constant = float(source_config["geometry"]["a"])
    frequency_window = source_config["frequency_window_thz"]
    half_span = 0.5 * (float(frequency_window[1]) - float(frequency_window[0]))
    gamma_p2 = frame[
        (frame["band_label"] == "p2")
        & np.isclose(frame["kx"], 0.0)
        & np.isclose(frame["ky"], 0.0)
    ]["re"].drop_duplicates()
    if len(gamma_p2) != 1:
        raise ValueError("Full-band source has inconsistent p2@Gamma frequencies")
    frequency_limits = (
        float(gamma_p2.iloc[0]) - half_span,
        float(gamma_p2.iloc[0]) + half_span,
    )
    gamma_center = float(
        frame[np.isclose(frame["kx"], 0.0) & np.isclose(frame["ky"], 0.0)]
        .drop_duplicates("band_label")["re"].mean()
    )
    calibration = summary["fixed_calibration"]
    parameters = summary["parameters"]
    model_kwargs = {
        "t0_thz": float(calibration["t0_thz"]),
        "t1_thz": float(calibration["t1_thz"]),
        "v_thz_um": float(calibration["velocity_thz_um"]),
        "alpha": float(parameters["alpha"]),
        "beta": float(parameters["beta"]),
        "gamma": float(parameters["gamma"]),
        "mu_thz": float(parameters["mu_thz"]),
    }
    branch_contract = (
        ("gamma_m", -1.0, 0.5),
        ("gamma_k", 1.0, 1.0 / np.sqrt(3.0)),
    )
    scale = reciprocal_scale_inv_um(lattice_constant)
    figure_width = 7.2 if include_second_order else 6.5
    figure, axis = plt.subplots(figsize=(figure_width, 4.45))
    first_order_residuals = {}
    second_order_residuals = {}
    for direction, sign, endpoint in branch_contract:
        points, observed = labeled_ray(frame, direction)
        magnitudes = np.linalg.norm(points, axis=1)
        path_coordinates = sign * magnitudes / endpoint
        smooth_magnitudes = np.linspace(0.0, endpoint, 601)
        smooth_points = (
            np.column_stack([np.zeros_like(smooth_magnitudes), smooth_magnitudes])
            if direction == "gamma_m"
            else np.column_stack([smooth_magnitudes, np.zeros_like(smooth_magnitudes)])
        )
        predicted_at_data, _ = tracked_s17_bands(points * scale, **model_kwargs)
        predicted_at_data += gamma_center
        smooth_predicted, _ = tracked_s17_bands(smooth_points * scale, **model_kwargs)
        smooth_predicted += gamma_center
        residual = predicted_at_data - observed
        first_order_residuals[direction] = {
            "rmse_thz": float(np.sqrt(np.mean(residual**2))),
            "max_abs_thz": float(np.max(np.abs(residual))),
        }
        if include_second_order:
            quadratic_at_data, _ = tracked_s17_bands(
                points * scale, **model_kwargs, include_quadratic=True,
                lattice_constant_um=lattice_constant,
            )
            quadratic_at_data += gamma_center
            quadratic_smooth, _ = tracked_s17_bands(
                smooth_points * scale, **model_kwargs, include_quadratic=True,
                lattice_constant_um=lattice_constant,
            )
            quadratic_smooth += gamma_center
            quadratic_residual = quadratic_at_data - observed
            second_order_residuals[direction] = {
                "rmse_thz": float(np.sqrt(np.mean(quadratic_residual**2))),
                "max_abs_thz": float(np.max(np.abs(quadratic_residual))),
            }
        for mode_index, mode in enumerate(STANDING_WAVE_LABELS):
            color = BAND_MODE_COLORS[mode]
            axis.plot(
                path_coordinates, observed[:, mode_index],
                color=color, lw=2.4, alpha=0.42,
            )
            axis.plot(
                sign * smooth_magnitudes / endpoint,
                smooth_predicted[:, mode_index],
                color=color, lw=1.3,
                linestyle=(0, (1.0, 1.7)) if include_second_order else (0, (3.0, 2.0)),
            )
            if include_second_order:
                axis.plot(
                    sign * smooth_magnitudes / endpoint,
                    quadratic_smooth[:, mode_index],
                    color=color, lw=1.35, linestyle=(0, (4.0, 2.0)),
                )
    fit_left = -0.02 / 0.5
    fit_right = float(0.02 / (1.0 / np.sqrt(3.0)))
    axis.axvspan(fit_left, fit_right, color="#eeeeee", zorder=-5)
    axis.axvline(0.0, color="#b8b8b8", lw=0.8)
    axis.set_xlim(-1.0, 1.0)
    axis.set_ylim(*frequency_limits)
    axis.set_xticks([-1.0, 0.0, 1.0], labels=["M", r"$\Gamma$", "K"])
    axis.set_xlabel("Wave-vector path", fontsize=9.5)
    axis.set_ylabel("Frequency (THz)", fontsize=9.5)
    axis.set_title("(245, 0.96, 1.156) full band", fontsize=10.5, pad=8)
    axis.text(
        0.5 * (fit_left + fit_right), frequency_limits[0] + 0.035 * (frequency_limits[1] - frequency_limits[0]),
        r"fit window: $|k|/G\leq0.02$", ha="center", va="bottom", fontsize=7.5,
    )
    axis.grid(axis="y", color="#d9d9d9", lw=0.6, alpha=0.55)
    axis.tick_params(direction="out", length=3, width=0.8, labelsize=8.5)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    source_handles = [
        Line2D([0], [0], color="#333333", lw=2.4, alpha=0.42, label="COMSOL"),
        Line2D(
            [0], [0], color="#333333", lw=1.3,
            linestyle=(0, (1.0, 1.7)) if include_second_order else (0, (3.0, 2.0)),
            label="S17 first order" if include_second_order else "S17 fit",
        ),
    ]
    if include_second_order:
        source_handles.append(
            Line2D(
                [0], [0], color="#333333", lw=1.35,
                linestyle=(0, (4.0, 2.0)), label="C5 second order",
            )
        )
    mode_handles = [
        Line2D([0], [0], color=BAND_MODE_COLORS[mode], lw=2.0, label=BAND_MODE_DISPLAY[mode])
        for mode in STANDING_WAVE_LABELS
    ]
    figure.legend(
        handles=source_handles + mode_handles, loc="upper center",
        ncol=7 if include_second_order else 6,
        frameon=False, fontsize=7.6, handlelength=2.2,
        columnspacing=1.0 if include_second_order else 1.25,
        bbox_to_anchor=(0.5, 1.01),
    )
    figure.subplots_adjust(left=0.11, right=0.985, bottom=0.13, top=0.82)
    figure.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    full_residuals = (
        {
            "first_order": first_order_residuals,
            "second_order": second_order_residuals,
        }
        if include_second_order else first_order_residuals
    )
    return {
        "output": str(output), "source": str(full_band_selected_path),
        "x_limits": [-1.0, 1.0], "frequency_limits_thz": list(frequency_limits),
        "fit_window_path_coordinates": [fit_left, fit_right],
        "full_path_residuals": full_residuals,
        "line_count": 24 if include_second_order else 16,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = fit_dataset(args.dataset, args.output_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Result: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
