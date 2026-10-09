from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit

from .finite_lattice_fourier import first_bz_polygon, reciprocal_lattice_basis


AIR_FIELD_COLUMNS = [
    "x",
    "y",
    "Ex_re",
    "Ex_im",
    "Ey_re",
    "Ey_im",
    "Ez_re",
    "Ez_im",
    "is_in_domain",
]

_FARFIELD_INTENSITY_CMAP = "magma"
_NORMALIZED_K_PLOT_LIMIT = 0.2
_FWHM_ARROW_LENGTH_NORMALIZED = 0.035
_COLORBAR_LABEL = "Normalized Intensity"
_COLORBAR_LABEL_ROTATION = 270
_COLORBAR_LABEL_PAD = 18
_LOG_INTENSITY_MIN = -5.0
_NA_CIRCLE_COLOR = "#9ca3af"


def brillouin_zone_polygons(period_um: float) -> list[np.ndarray]:
    if not np.isfinite(period_um) or float(period_um) <= 0.0:
        raise ValueError("period_um must be positive and finite")
    period_um = float(period_um)
    central = first_bz_polygon(period_um)
    reciprocal = reciprocal_lattice_basis(period_um)
    candidates = []
    for h0 in range(-2, 3):
        for h1 in range(-2, 3):
            if h0 == 0 and h1 == 0:
                continue
            vector = reciprocal @ np.array([h0, h1], dtype=float)
            candidates.append(vector)
    nearest_norm = min(float(np.linalg.norm(vector)) for vector in candidates)
    neighbours = [
        vector
        for vector in candidates
        if np.isclose(np.linalg.norm(vector), nearest_norm, rtol=1e-10, atol=1e-12)
    ]
    neighbours.sort(key=lambda vector: float(np.arctan2(vector[1], vector[0])))
    if len(neighbours) != 6:
        raise RuntimeError(f"Expected six nearest reciprocal-lattice neighbours, got {len(neighbours)}")
    return [central, *[central + vector for vector in neighbours]]


def polygon_axis_interval(polygon: np.ndarray, axis: str) -> tuple[float, float]:
    polygon = np.asarray(polygon, dtype=float)
    if polygon.ndim != 2 or polygon.shape[1] != 2 or len(polygon) < 3:
        raise ValueError("polygon must have shape (n, 2) with n >= 3")
    if axis not in {"kx", "ky"}:
        raise ValueError("axis must be 'kx' or 'ky'")
    coordinate_idx = 0 if axis == "kx" else 1
    transverse_idx = 1 - coordinate_idx
    intersections: list[float] = []
    tolerance = 1e-12
    for current, following in zip(polygon, np.roll(polygon, -1, axis=0), strict=True):
        current_t = float(current[transverse_idx])
        following_t = float(following[transverse_idx])
        if abs(current_t) <= tolerance:
            intersections.append(float(current[coordinate_idx]))
        if current_t * following_t < 0.0:
            fraction = -current_t / (following_t - current_t)
            intersections.append(
                float(current[coordinate_idx] + fraction * (following[coordinate_idx] - current[coordinate_idx]))
            )
    if len(intersections) < 2:
        raise ValueError("polygon does not cross the requested central axis")
    return min(intersections), max(intersections)


def validate_farfield_config(
    grid_size: int,
    fft_size: int,
    na: float,
    h0_um: float | None,
) -> None:
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("FARFIELD_GRID_SIZE must be an odd integer >= 3")
    if fft_size % 2 == 0:
        raise ValueError("FARFIELD_FFT_SIZE must be odd")
    if fft_size < grid_size:
        raise ValueError("FARFIELD_FFT_SIZE must be >= FARFIELD_GRID_SIZE")
    if (fft_size - grid_size) % 2:
        raise ValueError("Far-field zero padding must be symmetric")
    if not np.isfinite(na) or not 0.0 < float(na) <= 1.0:
        raise ValueError("FARFIELD_NA must satisfy 0 < NA <= 1")
    if h0_um is not None and (not np.isfinite(h0_um) or float(h0_um) < 0.0):
        raise ValueError("FARFIELD_H0_UM must be finite and nonnegative")


def symmetric_pad(field: np.ndarray, fft_size: int) -> np.ndarray:
    field = np.asarray(field, dtype=complex)
    if field.ndim != 2 or field.shape[0] != field.shape[1]:
        raise ValueError("Far-field input must be a square 2D array")
    if fft_size < field.shape[0] or (fft_size - field.shape[0]) % 2:
        raise ValueError("FFT size must permit symmetric padding")
    pad = (fft_size - field.shape[0]) // 2
    return np.pad(field, ((pad, pad), (pad, pad)), mode="constant")


def centered_fft2(field: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field)))


def centered_ifft2(field: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(field)))


def _uniform_spacing(axis: np.ndarray, name: str) -> float:
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or len(axis) < 2:
        raise ValueError(f"{name} axis must contain at least two points")
    diffs = np.diff(axis)
    spacing = float(diffs[0])
    if spacing <= 0.0 or not np.allclose(diffs, spacing, rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} axis must be strictly increasing and uniform")
    return spacing


def load_air_field(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    df = pd.read_parquet(path)
    missing = [column for column in AIR_FIELD_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing air-field columns in {path}: {missing}")

    x_axis = np.sort(df["x"].unique().astype(float))
    y_axis = np.sort(df["y"].unique().astype(float))
    if len(x_axis) != len(y_axis) or len(df) != len(x_axis) * len(y_axis):
        raise ValueError(f"Air field must be a complete square grid: {path}")
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")
    if not np.allclose(df["x"].to_numpy(dtype=float), xx.ravel()):
        raise ValueError(f"Air-field x coordinates are not in C-order grid layout: {path}")
    if not np.allclose(df["y"].to_numpy(dtype=float), yy.ravel()):
        raise ValueError(f"Air-field y coordinates are not in C-order grid layout: {path}")

    shape = (len(y_axis), len(x_axis))
    result: dict[str, Any] = {
        "path": path,
        "x_axis": x_axis,
        "y_axis": y_axis,
        "is_in_domain": df["is_in_domain"].astype(bool).to_numpy().reshape(shape),
    }
    for component in ("Ex", "Ey", "Ez"):
        result[component] = (
            df[f"{component}_re"].to_numpy(dtype=float)
            + 1j * df[f"{component}_im"].to_numpy(dtype=float)
        ).reshape(shape)
    return result


def compute_farfield_fft(
    ex: np.ndarray,
    ey: np.ndarray,
    ez: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    *,
    frequency_thz: float,
    h0_um: float,
    na: float,
    fft_size: int,
) -> dict[str, Any]:
    ex = np.asarray(ex, dtype=complex)
    ey = np.asarray(ey, dtype=complex)
    ez = np.asarray(ez, dtype=complex)
    if ex.shape != ey.shape or ex.shape != ez.shape or ex.ndim != 2 or ex.shape[0] != ex.shape[1]:
        raise ValueError("Ex, Ey, and Ez must share one square 2D grid")
    validate_farfield_config(ex.shape[0], fft_size, na, h0_um)
    dx = _uniform_spacing(np.asarray(x_axis), "x")
    dy = _uniform_spacing(np.asarray(y_axis), "y")
    if len(x_axis) != ex.shape[1] or len(y_axis) != ex.shape[0]:
        raise ValueError("Field shape does not match coordinate axes")
    if not np.isfinite(frequency_thz) or float(frequency_thz) <= 0.0:
        raise ValueError("frequency_thz must be positive and finite")

    lambda_um = 299.792458 / float(frequency_thz)
    k0 = 2.0 * np.pi / lambda_um
    norm_e0 = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2

    fex = centered_fft2(symmetric_pad(ex, fft_size))
    fey = centered_fft2(symmetric_pad(ey, fft_size))
    fez = centered_fft2(symmetric_pad(ez, fft_size))
    norm_fex = np.abs(fex) ** 2
    norm_fe0 = norm_fex + np.abs(fey) ** 2 + np.abs(fez) ** 2

    kx = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(fft_size, d=dx))
    ky = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(fft_size, d=dy))
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing="xy")
    kr2 = kx_grid**2 + ky_grid**2
    na_mask = kr2 <= (k0 * float(na)) ** 2
    kz = np.sqrt(np.maximum(k0**2 - kr2, 0.0))
    propagation = np.exp(1j * kz * float(h0_um)) * na_mask
    fex *= propagation
    fey *= propagation
    fez *= propagation
    del propagation

    norm_ife0 = np.abs(centered_ifft2(fex)) ** 2
    norm_ife0 += np.abs(centered_ifft2(fey)) ** 2
    norm_ife0 += np.abs(centered_ifft2(fez)) ** 2

    return {
        "x_axis": np.asarray(x_axis, dtype=float),
        "y_axis": np.asarray(y_axis, dtype=float),
        "kx": kx,
        "ky": ky,
        "normE0": norm_e0,
        "normfE0": norm_fe0,
        "normfEx": norm_fex,
        "normifE0": norm_ife0,
        "ffEx": fex,
        "ffEy": fey,
        "na_mask": na_mask,
        "lambda_um": lambda_um,
        "k0": k0,
        "dx_um": dx,
        "dy_um": dy,
        "grid_size": int(ex.shape[0]),
        "fft_size": int(fft_size),
    }


def gaussian_fwhm_metrics(fit_b: float, fit_c: float, k0: float) -> dict[str, float]:
    unavailable = {
        "fwhm_k_um_inv": np.nan,
        "lower_k_um_inv": np.nan,
        "upper_k_um_inv": np.nan,
        "divergence_full_angle_deg": np.nan,
    }
    if not np.all(np.isfinite([fit_b, fit_c, k0])) or float(k0) <= 0.0:
        return unavailable

    half_width = abs(float(fit_c)) * np.sqrt(np.log(2.0))
    lower_k = float(fit_b) - half_width
    upper_k = float(fit_b) + half_width
    lower_sine = lower_k / float(k0)
    upper_sine = upper_k / float(k0)
    divergence_deg = np.nan
    if -1.0 <= lower_sine <= 1.0 and -1.0 <= upper_sine <= 1.0:
        divergence_deg = float(
            np.degrees(np.arcsin(upper_sine) - np.arcsin(lower_sine))
        )
    return {
        "fwhm_k_um_inv": float(2.0 * half_width),
        "lower_k_um_inv": float(lower_k),
        "upper_k_um_inv": float(upper_k),
        "divergence_full_angle_deg": divergence_deg,
    }


def summarize_farfield(
    result: dict[str, Any],
    *,
    source_file: str,
    mode_idx: int,
    frequency_thz: float,
    q: float,
    h0_um: float,
    na: float,
    period_um: float,
) -> dict[str, Any]:
    angles = np.linspace(0.0, 2.0 * np.pi, 73)
    totals = np.array(
        [
            np.sum(np.abs(result["ffEy"] * np.cos(angle) - result["ffEx"] * np.sin(angle)) ** 2)
            for angle in angles
        ],
        dtype=float,
    )
    maximum = float(np.max(totals)) if totals.size else 0.0
    minimum = float(np.min(totals)) if totals.size else 0.0
    denominator = abs(maximum + minimum)
    polarization_degree = (maximum - minimum) / denominator if denominator else 0.0
    dominant_theta_deg = float(np.degrees(angles[int(np.argmax(totals))])) if totals.size else 0.0

    center = int(result["fft_size"]) // 2
    cut_ky = np.asarray(result["normfE0"], dtype=float)[:, center]
    ky = np.asarray(result["ky"], dtype=float)
    ky_min, ky_max = polygon_axis_interval(brillouin_zone_polygons(period_um)[0], "ky")
    fit_mask = (ky >= ky_min) & (ky <= ky_max)
    fit_ky = ky[fit_mask]
    fit_cut = cut_ky[fit_mask]
    fit_a = fit_b = fit_c = np.nan
    if len(fit_ky) >= 4 and float(np.max(fit_cut)) > 0.0:
        try:
            normalized = fit_cut / np.max(fit_cut)

            def gaussian(x, a, b, c):
                return a * np.exp(-((x - b) / c) ** 2)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                fit, _ = curve_fit(
                    gaussian,
                    fit_ky,
                    normalized,
                    p0=[1.0, 0.0, max(0.1, 0.1 * (ky_max - ky_min))],
                    maxfev=10000,
                )
            fit_a, fit_b, fit_c = map(float, fit)
        except (RuntimeError, ValueError, FloatingPointError):
            pass

    width_metrics = gaussian_fwhm_metrics(fit_b, fit_c, float(result["k0"]))
    k_scale = 2.0 * np.pi / float(period_um)

    return {
        "status": "ok",
        "mode_idx": int(mode_idx),
        "source_file": source_file,
        "frequency_thz": float(frequency_thz),
        "q": float(q),
        "lambda_um": float(result["lambda_um"]),
        "h0_um": float(h0_um),
        "grid_size": int(result["grid_size"]),
        "fft_size": int(result["fft_size"]),
        "dx_um": float(result["dx_um"]),
        "dy_um": float(result["dy_um"]),
        "na": float(na),
        "period_um": float(period_um),
        "polarization_degree": float(polarization_degree),
        "dominant_theta_deg": dominant_theta_deg,
        "gauss_a": fit_a,
        "gauss_b": fit_b,
        "gauss_c": fit_c,
        "gauss_fwhm_k_um_inv": width_metrics["fwhm_k_um_inv"],
        "gauss_fwhm_k_normalized": width_metrics["fwhm_k_um_inv"] / k_scale,
        "gauss_halfmax_lower_k_um_inv": width_metrics["lower_k_um_inv"],
        "gauss_halfmax_upper_k_um_inv": width_metrics["upper_k_um_inv"],
        "gauss_divergence_full_angle_deg": width_metrics["divergence_full_angle_deg"],
    }


def _normalized(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=float)
    maximum = float(np.nanmax(image)) if image.size else 0.0
    return image / maximum if maximum > 0.0 else image


def plot_farfield_figures(
    result: dict[str, Any],
    summary: dict[str, Any],
    output_dir: str | Path,
    *,
    period_um: float,
    dpi: int = 150,
) -> list[str]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []

    def save(fig, name: str) -> None:
        path = output_dir / name
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths.append(str(path))

    norm_fe = np.asarray(result["normfE0"], dtype=float)
    norm_if = np.asarray(result["normifE0"], dtype=float)
    x_axis = np.asarray(result["x_axis"], dtype=float)
    y_axis = np.asarray(result["y_axis"], dtype=float)
    kx = np.asarray(result["kx"], dtype=float)
    ky = np.asarray(result["ky"], dtype=float)
    fft_size = int(result["fft_size"])
    center = fft_size // 2
    propagated_x = (np.arange(fft_size) - center) * float(result["dx_um"])
    propagated_y = (np.arange(fft_size) - center) * float(result["dy_um"])
    k_scale = 2.0 * np.pi / float(period_um)
    bz_polygons = brillouin_zone_polygons(period_um)
    central_bz = bz_polygons[0]

    def crop_image(
        image: np.ndarray,
        x_values: np.ndarray,
        y_values: np.ndarray,
        bounds: tuple[float, float, float, float],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x_min, x_max, y_min, y_max = bounds
        x_keep = (x_values >= x_min) & (x_values <= x_max)
        y_keep = (y_values >= y_min) & (y_values <= y_max)
        if np.count_nonzero(x_keep) < 2 or np.count_nonzero(y_keep) < 2:
            raise ValueError("FFT grid does not cover the requested Brillouin-zone plot range")
        return image[np.ix_(y_keep, x_keep)], x_values[x_keep], y_values[y_keep]

    def centered_square_bounds(
        bounds: tuple[float, float, float, float],
    ) -> tuple[float, float, float, float]:
        half_width = 1.05 * max(abs(float(value)) for value in bounds)
        return -half_width, half_width, -half_width, half_width

    central_bounds = (
        float(np.min(central_bz[:, 0])),
        float(np.max(central_bz[:, 0])),
        float(np.min(central_bz[:, 1])),
        float(np.max(central_bz[:, 1])),
    )
    linear_bounds = centered_square_bounds(central_bounds)
    normalized_farfield = _normalized(norm_fe)
    log_farfield = np.log10(normalized_farfield + 1e-30)
    zoom_limit = _NORMALIZED_K_PLOT_LIMIT * k_scale
    zoom_bounds = (-zoom_limit, zoom_limit, -zoom_limit, zoom_limit)
    zoom_linear_image, _zoom_kx, _zoom_ky = crop_image(
        normalized_farfield, kx, ky, zoom_bounds
    )
    first_bz_log_image, _first_bz_kx, _first_bz_ky = crop_image(
        log_farfield, kx, ky, linear_bounds
    )

    all_vertices = np.vstack(bz_polygons)
    neighbouring_bounds = (
        float(np.min(all_vertices[:, 0])),
        float(np.max(all_vertices[:, 0])),
        float(np.min(all_vertices[:, 1])),
        float(np.max(all_vertices[:, 1])),
    )
    log_bounds = centered_square_bounds(neighbouring_bounds)
    log_image, _log_kx, _log_ky = crop_image(
        log_farfield, kx, ky, log_bounds
    )
    zoom_extent = [-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT] * 2
    linear_extent = [value / k_scale for value in linear_bounds]
    log_extent = [value / k_scale for value in log_bounds]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout="constrained")
    propagated_image = axes[0, 0].imshow(
        _normalized(norm_if),
        origin="lower",
        aspect="equal",
        extent=[propagated_x[0], propagated_x[-1], propagated_y[0], propagated_y[-1]],
        cmap=_FARFIELD_INTENSITY_CMAP,
        vmin=0.0,
        vmax=1.0,
    )
    axes[0, 0].set_title(f"Near field intensity (NA={float(summary['na']):.2f})")
    axes[0, 0].set_xlabel("x (μm)")
    axes[0, 0].set_ylabel("y (μm)")
    axes[0, 0].set_xlim(float(x_axis[0]), float(x_axis[-1]))
    axes[0, 0].set_ylim(float(y_axis[0]), float(y_axis[-1]))
    propagated_colorbar = fig.colorbar(propagated_image, ax=axes[0, 0])
    propagated_colorbar.set_label(
        _COLORBAR_LABEL,
        rotation=_COLORBAR_LABEL_ROTATION,
        labelpad=_COLORBAR_LABEL_PAD,
    )

    zoom_linear_plot = axes[0, 1].imshow(
        zoom_linear_image,
        origin="lower",
        aspect="equal",
        extent=zoom_extent,
        cmap=_FARFIELD_INTENSITY_CMAP,
        vmin=0.0,
        vmax=1.0,
    )
    axes[0, 1].set_title("Far field Intensity (k-space, linear)")
    axes[0, 1].set_xlabel("kx/(2π/a)")
    axes[0, 1].set_ylabel("ky/(2π/a)")
    axes[0, 1].set_xlim(-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT)
    axes[0, 1].set_ylim(-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT)
    zoom_linear_colorbar = fig.colorbar(zoom_linear_plot, ax=axes[0, 1])
    zoom_linear_colorbar.set_label(
        _COLORBAR_LABEL,
        rotation=_COLORBAR_LABEL_ROTATION,
        labelpad=_COLORBAR_LABEL_PAD,
    )

    first_bz_log_plot = axes[1, 0].imshow(
        first_bz_log_image,
        origin="lower",
        aspect="equal",
        extent=[
            linear_extent[0],
            linear_extent[1],
            linear_extent[2],
            linear_extent[3],
        ],
        cmap=_FARFIELD_INTENSITY_CMAP,
        vmin=_LOG_INTENSITY_MIN,
        vmax=0.0,
    )
    closed_central = np.vstack([central_bz, central_bz[0]]) / k_scale
    axes[1, 0].plot(closed_central[:, 0], closed_central[:, 1], color="white", linewidth=1.2)
    na_radius = float(result["k0"]) * float(summary["na"]) / k_scale
    axes[1, 0].add_patch(
        Circle(
            (0.0, 0.0),
            na_radius,
            fill=False,
            edgecolor=_NA_CIRCLE_COLOR,
            linewidth=1.0,
            linestyle="--",
        )
    )
    axes[1, 0].set_title("Far field intensity (1st BZ, log)")
    axes[1, 0].set_xlabel("kx/(2π/a)")
    axes[1, 0].set_ylabel("ky/(2π/a)")
    axes[1, 0].set_xlim(linear_extent[:2])
    axes[1, 0].set_ylim(linear_extent[2:])
    first_bz_log_colorbar = fig.colorbar(
        first_bz_log_plot,
        ax=axes[1, 0],
        ticks=np.arange(int(_LOG_INTENSITY_MIN), 1),
    )
    first_bz_log_colorbar.ax.set_yticklabels(
        [rf"$10^{{{power}}}$" for power in range(int(_LOG_INTENSITY_MIN), 1)]
    )
    first_bz_log_colorbar.set_label(
        _COLORBAR_LABEL,
        rotation=_COLORBAR_LABEL_ROTATION,
        labelpad=_COLORBAR_LABEL_PAD,
    )

    log_plot = axes[1, 1].imshow(
        log_image,
        origin="lower",
        aspect="equal",
        extent=[
            log_extent[0],
            log_extent[1],
            log_extent[2],
            log_extent[3],
        ],
        cmap=_FARFIELD_INTENSITY_CMAP,
        vmin=_LOG_INTENSITY_MIN,
        vmax=0.0,
    )
    for polygon in bz_polygons:
        closed = np.vstack([polygon, polygon[0]]) / k_scale
        axes[1, 1].plot(closed[:, 0], closed[:, 1], color="white", linewidth=0.9)
    axes[1, 1].add_patch(
        Circle(
            (0.0, 0.0),
            na_radius,
            fill=False,
            edgecolor=_NA_CIRCLE_COLOR,
            linewidth=1.0,
            linestyle="--",
        )
    )
    axes[1, 1].set_title("Far field Intensity (k-space, log)")
    axes[1, 1].set_xlabel("kx/(2π/a)")
    axes[1, 1].set_ylabel("ky/(2π/a)")
    axes[1, 1].set_xlim(log_extent[:2])
    axes[1, 1].set_ylim(log_extent[2:])
    log_colorbar = fig.colorbar(
        log_plot,
        ax=axes[1, 1],
        ticks=np.arange(int(_LOG_INTENSITY_MIN), 1),
    )
    log_colorbar.ax.set_yticklabels(
        [rf"$10^{{{power}}}$" for power in range(int(_LOG_INTENSITY_MIN), 1)]
    )
    log_colorbar.set_label(
        _COLORBAR_LABEL,
        rotation=_COLORBAR_LABEL_ROTATION,
        labelpad=_COLORBAR_LABEL_PAD,
    )

    for axis in axes.ravel():
        axis.set_box_aspect(1.0)
    fig.suptitle(
        f"f={float(summary['frequency_thz']):.3f} THz, "
        f"Q={float(summary['q']):.2f}"
    )
    save(fig, "A_overview.png")

    angles = np.linspace(0.0, 2.0 * np.pi, 73)
    totals = np.array(
        [
            np.sum(np.abs(result["ffEy"] * np.cos(angle) - result["ffEx"] * np.sin(angle)) ** 2)
            for angle in angles
        ]
    )
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "polar"})
    ax.plot(angles, _normalized(totals))
    ax.set_theta_zero_location("N")
    ax.set_title(f"Polarization degree={summary['polarization_degree']:.3f}")
    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("Normalized transmitted power")
    ax.set_ylim(0.0, 1.05)
    save(fig, "B_polarization.png")

    kx_min, kx_max = polygon_axis_interval(central_bz, "kx")
    ky_min, ky_max = polygon_axis_interval(central_bz, "ky")
    kx_keep = (kx >= kx_min) & (kx <= kx_max)
    ky_keep = (ky >= ky_min) & (ky <= ky_max)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(ky[ky_keep] / k_scale, _normalized(norm_fe[ky_keep, center]))
    axes[0].set_title("ky cut at kx=0 (central BZ)")
    axes[0].set_xlabel("ky/(2π/a)")
    axes[0].set_ylabel("Normalized |E(k)|²")
    axes[0].set_xlim(-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT)
    axes[1].plot(kx[kx_keep] / k_scale, _normalized(norm_fe[center, kx_keep]))
    axes[1].set_title("kx cut at ky=0 (central BZ)")
    axes[1].set_xlabel("kx/(2π/a)")
    axes[1].set_ylabel("Normalized |E(k)|²")
    axes[1].set_xlim(-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT)
    fig.tight_layout()
    save(fig, "C_cutlines.png")

    fig, ax = plt.subplots(figsize=(6, 4))
    fit_ky = ky[ky_keep]
    fit_data = _normalized(norm_fe[ky_keep, center])
    ax.plot(fit_ky / k_scale, fit_data, label="FFT data")
    fit_a = float(summary.get("gauss_a", np.nan))
    fit_b = float(summary.get("gauss_b", np.nan))
    fit_c = float(summary.get("gauss_c", np.nan))
    if np.all(np.isfinite([fit_a, fit_b, fit_c])) and fit_c != 0.0:
        fitted = fit_a * np.exp(-((fit_ky - fit_b) / fit_c) ** 2)
        ax.plot(fit_ky / k_scale, fitted, "--", label="Gaussian fit")
    fwhm_k = float(summary.get("gauss_fwhm_k_um_inv", np.nan))
    fwhm_normalized = float(summary.get("gauss_fwhm_k_normalized", np.nan))
    lower_k = float(summary.get("gauss_halfmax_lower_k_um_inv", np.nan))
    upper_k = float(summary.get("gauss_halfmax_upper_k_um_inv", np.nan))
    divergence_deg = float(summary.get("gauss_divergence_full_angle_deg", np.nan))
    if np.all(np.isfinite([fit_a, fwhm_k, fwhm_normalized, lower_k, upper_k, divergence_deg])):
        half_maximum = 0.5 * fit_a
        lower_normalized = lower_k / k_scale
        upper_normalized = upper_k / k_scale
        for endpoint, tail in (
            (lower_normalized, lower_normalized - _FWHM_ARROW_LENGTH_NORMALIZED),
            (upper_normalized, upper_normalized + _FWHM_ARROW_LENGTH_NORMALIZED),
        ):
            ax.annotate(
                "",
                xy=(endpoint, half_maximum),
                xytext=(tail, half_maximum),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "tab:red",
                    "linewidth": 1.3,
                    "shrinkA": 0.0,
                    "shrinkB": 0.0,
                },
            )
        ax.text(
            0.02,
            0.98,
            (
                f"Fit FWHM_k = {fwhm_k:.4g} μm⁻¹\n"
                f"= {fwhm_normalized:.4g} × (2π/a)\n"
                f" FWHM divergence = {divergence_deg:.3f}°"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "alpha": 0.72,
            },
        )
    else:
        ax.text(
            0.03,
            0.97,
            "Gaussian FWHM unavailable",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "alpha": 0.72,
            },
        )
    ax.set_title("Gaussian fit of ky cut (central BZ)")
    ax.set_xlabel("ky/(2π/a)")
    ax.set_ylabel("Normalized |E(kx=0, ky)|²")
    ax.set_xlim(-_NORMALIZED_K_PLOT_LIMIT, _NORMALIZED_K_PLOT_LIMIT)
    ax.legend(fontsize=8.0)
    fig.tight_layout()
    save(fig, "D_gauss_fit.png")
    return paths


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")


def run_farfield_fft_file(
    path: str | Path,
    output_dir: str | Path,
    *,
    mode_idx: int,
    frequency_thz: float,
    q: float,
    h0_um: float,
    na: float,
    fft_size: int,
    period_um: float,
    plot: bool = True,
) -> dict[str, Any]:
    field = load_air_field(path)
    result = compute_farfield_fft(
        field["Ex"],
        field["Ey"],
        field["Ez"],
        field["x_axis"],
        field["y_axis"],
        frequency_thz=frequency_thz,
        h0_um=h0_um,
        na=na,
        fft_size=fft_size,
    )
    summary = summarize_farfield(
        result,
        source_file=str(path),
        mode_idx=mode_idx,
        frequency_thz=frequency_thz,
        q=q,
        h0_um=h0_um,
        na=na,
        period_um=period_um,
    )
    output_dir = Path(output_dir)
    summary["plots"] = (
        plot_farfield_figures(result, summary, output_dir, period_um=period_um) if plot else []
    )
    _write_json(output_dir / "summary.json", summary)
    return summary


def _csv_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def run_farfield_fft_case(
    case_dir: str | Path,
    export_dir: str | Path,
    *,
    grid_size: int,
    fft_size: int,
    na: float,
    h0_um: float | None,
    period_um: float,
    mode_indices: list[int] | None = None,
    plot: bool = True,
) -> pd.DataFrame:
    validate_farfield_config(grid_size, fft_size, na, h0_um)
    if not np.isfinite(period_um) or float(period_um) <= 0.0:
        raise ValueError("period_um must be positive and finite")
    case_dir = Path(case_dir)
    export_dir = Path(export_dir)
    output_dir = case_dir / "farfield_fft"
    output_dir.mkdir(parents=True, exist_ok=True)

    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(f"Missing finite eigenfrequency table: {eigen_path}")
    metadata_path = export_dir / "air_field_metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing air-field metadata: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    automatic_h0 = float(metadata["air_plane_z_um"])
    effective_h0 = automatic_h0 if h0_um is None else float(h0_um)

    config = {
        "grid_size": int(grid_size),
        "fft_size": int(fft_size),
        "na": float(na),
        "configured_h0_um": h0_um,
        "automatic_h0_um": automatic_h0,
        "effective_h0_um": effective_h0,
        "period_um": float(period_um),
        "mode_indices": (
            None if mode_indices is None else [int(value) for value in mode_indices]
        ),
    }
    _write_json(output_dir / "config.json", config)

    eigen_df = pd.read_csv(eigen_path)
    valid = eigen_df.loc[eigen_df["is_valid"].map(_csv_bool)].copy()
    if mode_indices is not None:
        requested = [int(mode_idx) for mode_idx in mode_indices]
        if len(requested) != len(set(requested)):
            raise ValueError("mode_indices must not contain duplicates")
        valid_modes = set(valid["mode_idx"].astype(int).tolist())
        unknown = sorted(set(requested).difference(valid_modes))
        if unknown:
            raise ValueError(
                "Far-field mode_indices must reference valid modes; "
                f"invalid modes: {unknown}"
            )
        order = {mode_idx: idx for idx, mode_idx in enumerate(requested)}
        valid = valid.loc[valid["mode_idx"].astype(int).isin(requested)].copy()
        valid["_requested_order"] = valid["mode_idx"].astype(int).map(order)
        valid = valid.sort_values("_requested_order").drop(
            columns="_requested_order"
        )
    rows: list[dict[str, Any]] = []
    for _, eigen_row in valid.iterrows():
        mode_idx = int(eigen_row["mode_idx"])
        field_path = export_dir / f"{mode_idx:02d}_E_air.parquet"
        try:
            rows.append(
                run_farfield_fft_file(
                    field_path,
                    output_dir / f"mode_{mode_idx:02d}",
                    mode_idx=mode_idx,
                    frequency_thz=float(eigen_row["re"]),
                    q=float(eigen_row["q"]),
                    h0_um=effective_h0,
                    na=na,
                    fft_size=fft_size,
                    period_um=period_um,
                    plot=plot,
                )
            )
        except Exception as exc:
            rows.append(
                {
                    "status": "error",
                    "mode_idx": mode_idx,
                    "source_file": str(field_path),
                    "frequency_thz": float(eigen_row["re"]),
                    "q": float(eigen_row["q"]),
                    "h0_um": effective_h0,
                    "grid_size": int(grid_size),
                    "fft_size": int(fft_size),
                    "na": float(na),
                    "period_um": float(period_um),
                    "error": str(exc),
                }
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "farfield_summary.csv", index=False)
    if len(valid) > 0 and not any(row.get("status") == "ok" for row in rows):
        raise RuntimeError(f"Far-field FFT failed for all {len(valid)} valid modes")
    return summary
