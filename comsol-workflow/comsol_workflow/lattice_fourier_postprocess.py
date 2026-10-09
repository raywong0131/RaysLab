from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

from .basis_utils import standard_to_fourier_basis
from .energy_recovery import field_to_vector, fourier_subspace_energies_field
from .finite_lattice_fourier import (
    arbitrary_finite_lattice_fourier_components,
    finite_lattice_fourier_components,
    first_bz_polygon,
    fold_xis_to_first_bz,
    hex_cyclic_indices,
    hex_quotient_cell_count,
    hex_spiral_indices,
    hex_xi_values,
    unit_cell_rho_grid,
)
from .hex_lattice_utils import LatticePoint, cell_polygon, unit_cell_corners
from .polygon_utils import point_inside_or_on_polygon, polygon_area


INDEX_LABEL_SHELL_COUNT = 4
INTENSITY_MAPS_NORMH_FILENAME = "intensity_maps_normH.png"
CENTER_NORMALIZED_INTENSITY_CMAP = "RdYlBu_r"
CENTER_NORMALIZED_INTENSITY_VMIN = 0.0
CENTER_NORMALIZED_INTENSITY_VCENTER = 1.0
CENTER_NORMALIZED_INTENSITY_VMAX = 2.0
CENTER_NORMALIZED_INTENSITY_TICKS = (0.0, 0.5, 1.0, 1.5, 2.0)
CENTER_NORMALIZED_INTENSITY_TICK_LABELS = (
    "0",
    "0.5",
    "1.0",
    "1.5",
    "2+",
)
CENTER_NORMALIZED_INTENSITY_TITLE = (
    "Unit-Cell Mean Hz Intensity\n"
    r"Linear scale; colorbar midpoint marks $I_R/I_0=1$"
)
CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL = (
    r"$I_R / I_0$ (linear scale)"
)


def spiral_index_annotation_mask(
    indices: np.ndarray,
    shell_count: int = INDEX_LABEL_SHELL_COUNT,
) -> np.ndarray:
    if shell_count < 1:
        raise ValueError("shell_count must be at least 1")
    outer_shell = shell_count - 1
    max_labeled_index = 3 * outer_shell * (outer_shell + 1)
    return np.asarray(indices, dtype=int) <= max_labeled_index


def normalize_profile_magnitudes(profiles: np.ndarray) -> np.ndarray:
    values = np.abs(np.asarray(profiles))
    if values.size == 0:
        return values.astype(float)
    maximum = float(np.max(values))
    if not np.isfinite(maximum) or maximum <= 0.0:
        return np.zeros_like(values, dtype=float)
    return values / maximum


def closed_line(vertices: np.ndarray) -> np.ndarray:
    return np.vstack([vertices, vertices[0]])


def to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: to_jsonable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [to_jsonable(item) for item in value]
    return value


def write_json(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(to_jsonable(data), handle, indent=2)


def csv_bool(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def read_hz_center_field(field_path: Path) -> tuple[np.ndarray, np.ndarray]:
    field_df = pd.read_parquet(field_path)
    points = field_df[["x", "y"]].to_numpy(dtype=float)
    field = field_df["re"].to_numpy(dtype=float) + 1j * field_df["im"].to_numpy(dtype=float)
    return points, field


def interpolate_complex_field(points: np.ndarray, field: np.ndarray, query_points: np.ndarray) -> np.ndarray:
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    re_linear = LinearNDInterpolator(points, np.real(field))
    im_linear = LinearNDInterpolator(points, np.imag(field))
    re = re_linear(query_points)
    im = im_linear(query_points)

    missing = ~np.isfinite(re) | ~np.isfinite(im)
    if np.any(missing):
        re_nearest = NearestNDInterpolator(points, np.real(field))
        im_nearest = NearestNDInterpolator(points, np.imag(field))
        re[missing] = re_nearest(query_points[missing])
        im[missing] = im_nearest(query_points[missing])
    return re + 1j * im


def decompose_profiles(
    rho_points: np.ndarray,
    profiles: np.ndarray,
    weights: np.ndarray,
    weight_fraction: np.ndarray,
    index_name: str,
) -> pd.DataFrame:
    subspace_idx_to_name = {0: "s", 1: "p", 2: "d", 3: "f"}
    mode_idx_to_name = {0: "s", 1: "px", 2: "py", 3: "dx", 4: "dy", 5: "f"}
    subspace_to_mode_names = {"s": ["s"], "p": ["px", "py"], "d": ["dx", "dy"], "f": ["f"]}
    rows = []

    for idx, profile in enumerate(profiles):
        df_re = pd.DataFrame({0: rho_points[:, 0], 1: rho_points[:, 1], 2: np.real(profile)})
        df_im = pd.DataFrame({0: rho_points[:, 0], 1: rho_points[:, 1], 2: np.imag(profile)})

        subspace_res = fourier_subspace_energies_field(df_re, df_im, N=6)
        subspace_values = np.concatenate([
            [subspace_res["E_dc"]],
            subspace_res["E_pairs"],
            [subspace_res["E_nyquist"]] if subspace_res["E_nyquist"] is not None else [],
        ])
        subspace_total = float(subspace_res["total_energy"])
        if subspace_total > 0.0 and np.isfinite(subspace_total):
            subspace_fraction = subspace_values / subspace_total
        else:
            subspace_fraction = np.full_like(subspace_values, np.nan, dtype=float)
        dominant_subspace = int(np.nanargmax(subspace_fraction))
        dominant_subspace_name = subspace_idx_to_name[dominant_subspace]

        standard_coeffs = field_to_vector(df_re, df_im, N=6)
        fourier_coeffs = standard_to_fourier_basis(standard_coeffs, N=6)
        mode_values = np.abs(fourier_coeffs) ** 2
        mode_total = float(np.sum(mode_values))
        if mode_total > 0.0 and np.isfinite(mode_total):
            mode_fraction = mode_values / mode_total
        else:
            mode_fraction = np.full_like(mode_values, np.nan, dtype=float)
        dominant_mode = int(np.nanargmax(mode_fraction))
        dominant_mode_name = mode_idx_to_name[dominant_mode]

        rows.append({
            index_name: int(idx),
            "k_weight": float(weights[idx]),
            "k_weight_fraction": float(weight_fraction[idx]),
            "subspace_s": float(subspace_fraction[0]),
            "subspace_p": float(subspace_fraction[1]),
            "subspace_d": float(subspace_fraction[2]),
            "subspace_f": float(subspace_fraction[3]),
            "dominant_subspace": dominant_subspace_name,
            "dominant_subspace_fraction": float(subspace_fraction[dominant_subspace]),
            "mode_s": float(mode_fraction[0]),
            "mode_px": float(mode_fraction[1]),
            "mode_py": float(mode_fraction[2]),
            "mode_dx": float(mode_fraction[3]),
            "mode_dy": float(mode_fraction[4]),
            "mode_f": float(mode_fraction[5]),
            "dominant_mode": dominant_mode_name,
            "dominant_mode_fraction": float(mode_fraction[dominant_mode]),
            "subspace_matches_mode": bool(dominant_mode_name in subspace_to_mode_names[dominant_subspace_name]),
        })

    return pd.DataFrame(rows)


def strip_bulk_fourier_dir(out_dir: Path) -> Path:
    return out_dir / "strip_bulk_fourier_hz"


def finite_lattice_fourier_dir(out_dir: Path) -> Path:
    return out_dir / "finite_lattice_fourier_hz"


def simulation_export_dir(out_dir: Path, export_dirname: str) -> Path:
    return out_dir / export_dirname


def strip_bulk_rows(
    records: list[dict[str, Any]],
    bulk_radius: int,
    row_index_coefficients: tuple[int, int] = (0, 1),
) -> list[int]:
    i_coefficient, j_coefficient = row_index_coefficients
    rows = sorted({
        int(i_coefficient * record["point"].i + j_coefficient * record["point"].j)
        for record in records
    })
    expected = list(range(-bulk_radius, bulk_radius + 1))
    if rows != expected:
        raise ValueError(f"Expected bulk strip rows {expected}, got {rows}")
    return rows


def strip_bulk_row_centers(rows: list[int], period: float) -> np.ndarray:
    centers = []
    for row in rows:
        centers.append(np.array([0.5 * period * row, np.sqrt(3.0) * period * row / 2.0], dtype=float))
    return np.asarray(centers, dtype=float)


def strip_bulk_rho_grid(period: float, grid_size: int) -> dict[str, np.ndarray]:
    base_grid = unit_cell_rho_grid(period, grid_size)
    return {
        "points": np.asarray(base_grid["points"], dtype=float),
        "cell": unit_cell_corners(period),
    }


def wrap_periodic_x(points: np.ndarray, period: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(points, dtype=float)
    wrapped = points.copy()
    half_period = period / 2.0
    wrapped_x = ((wrapped[:, 0] + half_period) % period) - half_period
    shifts = np.rint((wrapped[:, 0] - wrapped_x) / period).astype(int)
    wrapped[:, 0] = wrapped_x
    return wrapped, shifts


def extract_strip_bulk_row_fields(
    field_points: np.ndarray,
    field: np.ndarray,
    rho_points: np.ndarray,
    row_centers: np.ndarray,
    kx_fraction: float,
    period: float,
) -> np.ndarray:
    row_centers = np.asarray(row_centers, dtype=float)
    rho_points = np.asarray(rho_points, dtype=float)
    query_unwrapped = row_centers[:, None, :] + rho_points[None, :, :]
    query_shape = query_unwrapped.shape[:2]
    query_wrapped, shifts = wrap_periodic_x(query_unwrapped.reshape(-1, 2), period)
    sampled = interpolate_complex_field(field_points, field, query_wrapped)
    bloch_phase = np.exp(-2j * np.pi * kx_fraction * shifts)
    return np.asarray(bloch_phase * sampled, dtype=complex).reshape(query_shape)


def strip_k_values(n_rows: int, kx_fraction: float, period: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    modes = np.arange(n_rows, dtype=int)
    xi = modes.astype(float) / float(n_rows)
    xi_folded = ((xi + 0.5) % 1.0) - 0.5
    b_perp = np.array([0.0, 4.0 * np.pi / (np.sqrt(3.0) * period)], dtype=float)
    k_parallel = np.array([kx_fraction * 2.0 * np.pi / period, 0.0], dtype=float)
    raw_k = k_parallel[None, :] + xi[:, None] * b_perp[None, :]
    folded_k = k_parallel[None, :] + xi_folded[:, None] * b_perp[None, :]
    return xi, raw_k, folded_k


def strip_fourier_components(
    row_fields_by_residue: np.ndarray,
    rho_points: np.ndarray,
    raw_k: np.ndarray,
    period: float,
) -> dict[str, np.ndarray | float]:
    f_components = np.fft.ifft(row_fields_by_residue, axis=0, norm="ortho")
    reconstructed = np.fft.fft(f_components, axis=0, norm="ortho")
    field_norm = float(np.linalg.norm(row_fields_by_residue))
    reconstruction_error = float(np.linalg.norm(reconstructed - row_fields_by_residue) / field_norm)

    phase = np.exp(1j * (raw_k @ rho_points.T))
    g_components = f_components * phase

    sample_weight = polygon_area(unit_cell_corners(period)) / len(rho_points)
    original_norm = float(sample_weight * np.sum(np.abs(row_fields_by_residue) ** 2))
    fourier_norm = float(sample_weight * np.sum(np.abs(f_components) ** 2))
    parseval_error = abs(original_norm - fourier_norm) / original_norm
    weights = sample_weight * np.sum(np.abs(g_components) ** 2, axis=1)
    weight_fraction = weights / float(np.sum(weights))

    p_components = np.zeros_like(g_components)
    positive_weight = weights > 0.0
    p_components[positive_weight] = g_components[positive_weight] / np.sqrt(weights[positive_weight])[:, None]
    profile_norms = sample_weight * np.sum(np.abs(p_components) ** 2, axis=1)
    profile_norm_error = float(np.max(np.abs(profile_norms[positive_weight] - 1.0))) if np.any(positive_weight) else np.nan

    return {
        "F": f_components,
        "G": g_components,
        "P": p_components,
        "weights": weights,
        "weight_fraction": weight_fraction,
        "profile_norms": profile_norms,
        "reconstruction_error": reconstruction_error,
        "parseval_error": float(parseval_error),
        "profile_norm_error": profile_norm_error,
    }


def plot_strip_bulk_rows(output_path: Path, records: list[dict[str, Any]], rows: list[int], *, period: float, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap("coolwarm")
    norm = plt.Normalize(min(rows), max(rows))

    fig, ax = plt.subplots(figsize=(4.8, 10.0), constrained_layout=True)
    for record in records:
        point = record["point"]
        color = cmap(norm(point.j))
        ax.add_patch(Polygon(
            record["polygon"],
            closed=True,
            facecolor=color,
            edgecolor="#111827",
            linewidth=0.45,
            alpha=0.72,
        ))
        center = np.mean(record["polygon"], axis=0)
        ax.text(center[0], center[1], str(point.j), ha="center", va="center", fontsize=6.5)

    polygons = [np.asarray(record["polygon"], dtype=float) for record in records]
    points = np.vstack(polygons)
    span = max(float(np.ptp(points[:, 0])), float(np.ptp(points[:, 1])), period)
    pad = 0.08 * span
    ax.set_xlim(float(points[:, 0].min() - pad), float(points[:, 0].max() + pad))
    ax.set_ylim(float(points[:, 1].min() - pad), float(points[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.grid(True, alpha=0.16)
    ax.set_title("strip bulk row representatives")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="bulk row j")
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_strip_k_grid(output_path: Path, folded_k: np.ndarray, *, period: float, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scale = 2.0 * np.pi / period
    y = folded_k[:, 1] / scale

    fig, ax = plt.subplots(figsize=(6.8, 2.8), constrained_layout=True)
    ax.axhline(0.0, color="#9ca3af", linewidth=0.7)
    ax.scatter(y, np.zeros(len(y)), s=48, c="#2563eb", edgecolor="white", linewidth=0.5, zorder=3)
    for mode_idx, y_value in enumerate(y):
        ax.text(y_value, 0.0, str(mode_idx), fontsize=7.0, ha="center", va="bottom")
    ax.set_ylim(-0.08, 0.12)
    ax.set_yticks([])
    ax.set_xlabel(r"$k_\perp / (2\pi/a)$")
    ax.set_title(f"strip bulk Z_{len(folded_k)} Fourier k points")
    ax.grid(True, axis="x", alpha=0.18)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_strip_mode_k_weights(
    output_path: Path,
    mode_idx: int,
    eigen_row,
    folded_k: np.ndarray,
    weight_fraction: np.ndarray,
    *,
    period: float,
    dpi: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scale = 2.0 * np.pi / period
    x = folded_k[:, 1] / scale
    order = np.argsort(x)

    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    ax.vlines(x[order], 0.0, weight_fraction[order], color="#64748b", linewidth=1.0, alpha=0.8)
    ax.scatter(x, weight_fraction, s=70, c=np.log10(np.maximum(weight_fraction, np.finfo(float).tiny)), cmap="magma", edgecolor="white", linewidth=0.5, zorder=3)
    for m_idx, x_value in enumerate(x):
        ax.text(x_value, weight_fraction[m_idx], f"m={m_idx}", fontsize=7.0, ha="center", va="bottom")

    title = f"mode {mode_idx:02d}: strip bulk Hz weight"
    if "re" in eigen_row and np.isfinite(float(eigen_row["re"])):
        title += f", f={float(eigen_row['re']):.3f} THz"
    if "q" in eigen_row and np.isfinite(float(eigen_row["q"])):
        title += f", Q={float(eigen_row['q']):.2f}"
    ax.set_title(title)
    ax.set_xlabel(r"$k_\perp / (2\pi/a)$")
    ax.set_ylabel("weight fraction")
    ax.set_ylim(0.0, max(1e-12, float(weight_fraction.max())) * 1.18)
    ax.grid(True, alpha=0.18)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_top_profiles(
    output_path: Path,
    mode_idx: int,
    rho_grid: dict[str, np.ndarray],
    profiles: np.ndarray,
    weight_fraction: np.ndarray,
    top_order: np.ndarray,
    *,
    period: float,
    top_k: int,
    phase: bool,
    title_prefix: str,
    dpi: int,
    display_indices: np.ndarray | None = None,
) -> None:
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    output_path.parent.mkdir(parents=True, exist_ok=True)
    top = top_order[: min(top_k, len(top_order))]
    if len(top) == 0:
        raise ValueError("top_order and top_k must select at least one profile")
    if display_indices is None:
        display_indices = np.arange(len(profiles), dtype=int)
    else:
        display_indices = np.asarray(display_indices, dtype=int)
        if display_indices.shape != (len(profiles),):
            raise ValueError("display_indices must align with profiles")
    rho_points = np.asarray(rho_grid["points"], dtype=float)
    cell = np.asarray(rho_grid.get("cell", unit_cell_corners(period)), dtype=float)
    values = np.angle(profiles[top]) if phase else normalize_profile_magnitudes(profiles[top])
    cmap = "twilight_shifted" if phase else "plasma"
    vmin = -np.pi if phase else 0.0
    vmax = np.pi if phase else 1.0

    x_axis = np.linspace(float(cell[:, 0].min()), float(cell[:, 0].max()), 3 * int(np.sqrt(len(rho_points))))
    y_axis = np.linspace(float(cell[:, 1].min()), float(cell[:, 1].max()), 3 * int(np.sqrt(len(rho_points))))
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")
    dense_points = np.column_stack([xx.ravel(), yy.ravel()])
    mask = np.array([point_inside_or_on_polygon(point, cell) for point in dense_points], dtype=bool).reshape(xx.shape)

    n_cols = min(3, len(top))
    n_rows = int(np.ceil(len(top) / n_cols))
    fig = plt.figure(
        figsize=(3.6 * n_cols + 0.8, 3.3 * n_rows),
        constrained_layout=True,
    )
    grid = fig.add_gridspec(
        n_rows,
        n_cols + 1,
        width_ratios=[1.0] * n_cols + [0.08],
    )
    axes_flat = np.array(
        [fig.add_subplot(grid[row, col]) for row in range(n_rows) for col in range(n_cols)],
        dtype=object,
    )
    colorbar_ax = fig.add_subplot(grid[:, -1])
    cell_closed = closed_line(cell)
    last_mesh = None
    for ax, k_idx, component in zip(axes_flat, top, values, strict=False):
        linear = LinearNDInterpolator(rho_points, component)
        image_flat = linear(dense_points)
        missing = ~np.isfinite(image_flat)
        if np.any(missing):
            nearest = NearestNDInterpolator(rho_points, component)
            image_flat[missing] = nearest(dense_points[missing])
        image = image_flat.reshape(xx.shape)
        image = np.where(mask, image, np.nan)
        last_mesh = ax.pcolormesh(x_axis, y_axis, image, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.plot(cell_closed[:, 0], cell_closed[:, 1], color="#111827", linewidth=0.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(float(cell[:, 0].min()), float(cell[:, 0].max()))
        ax.set_ylim(float(cell[:, 1].min()), float(cell[:, 1].max()))
        ax.set_title(
            f"k-index={int(display_indices[k_idx])}, "
            f"P(k)={weight_fraction[k_idx]:.3g}",
            fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes_flat[len(top):]:
        ax.axis("off")

    if last_mesh is not None:
        colorbar_ticks = [-np.pi, 0.0, np.pi] if phase else None
        colorbar = fig.colorbar(last_mesh, cax=colorbar_ax, ticks=colorbar_ticks)
        if phase:
            colorbar.ax.set_yticklabels([r"$-\pi$", r"$0$", r"$+\pi$"])
            colorbar_label = r"$\arg\,u_k(r)$"
        else:
            colorbar_label = r"Normalized $|u_k(r)|$"
        colorbar.set_label(colorbar_label, rotation=270, labelpad=20)
    figure_title = (
        "Top weight finite-cavity profiles: Phase"
        if phase
        else "Top weight finite-cavity profiles: Intensity"
    )
    fig.suptitle(figure_title, fontsize=12)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def strip_bulk_fourier_mode(
    export_dir: Path,
    output_root: Path,
    mode_idx: int,
    eigen_row,
    rho_grid: dict[str, np.ndarray],
    rows: list[int],
    row_centers: np.ndarray,
    xi_values: np.ndarray,
    raw_k: np.ndarray,
    folded_k: np.ndarray,
    kx_fraction: float,
    *,
    period: float,
    top_k: int,
    dpi: int,
    generate_plots: bool = True,
) -> dict[str, float | int | str]:
    field_path = export_dir / f"{mode_idx:02d}_Hz_center.parquet"
    if not field_path.exists():
        raise FileNotFoundError(f"Missing strip Hz field parquet: {field_path}")

    field_points, field = read_hz_center_field(field_path)
    rho_points = np.asarray(rho_grid["points"], dtype=float)
    row_fields = extract_strip_bulk_row_fields(field_points, field, rho_points, row_centers, kx_fraction, period)

    n_rows = len(rows)
    residues = np.mod(np.asarray(rows, dtype=int), n_rows)
    if set(residues.tolist()) != set(range(n_rows)):
        raise ValueError(f"Bulk strip rows are not Z_{n_rows} representatives: {rows}")
    row_fields_by_residue = np.empty_like(row_fields)
    row_fields_by_residue[residues] = row_fields

    fourier = strip_fourier_components(row_fields_by_residue, rho_points, raw_k, period)
    f_components = fourier["F"]
    g_components = fourier["G"]
    p_components = fourier["P"]
    profile_norms = fourier["profile_norms"]
    weights = fourier["weights"]
    weight_fraction = fourier["weight_fraction"]
    reconstruction_error = float(fourier["reconstruction_error"])
    parseval_error = float(fourier["parseval_error"])
    profile_norm_error = float(fourier["profile_norm_error"])
    top_order = np.argsort(weight_fraction)[::-1]

    output_dir = strip_bulk_fourier_dir(output_root) / f"mode_{mode_idx:02d}"
    output_dir.mkdir(parents=True, exist_ok=True)

    decomposition_df = decompose_profiles(rho_points, p_components, weights, weight_fraction, "m")
    decomposition_df["xi"] = xi_values
    decomposition_df["kx"] = raw_k[:, 0]
    decomposition_df["ky"] = raw_k[:, 1]
    decomposition_df["kx_folded"] = folded_k[:, 0]
    decomposition_df["ky_folded"] = folded_k[:, 1]
    decomposition_df["kx_folded_over_2pi_a"] = folded_k[:, 0] / (2.0 * np.pi / period)
    decomposition_df["ky_folded_over_2pi_a"] = folded_k[:, 1] / (2.0 * np.pi / period)
    decomposition_df.to_csv(output_dir / "p_subspace_mode_decomposition.csv", index=False)

    peaks = []
    for rank, k_idx in enumerate(top_order, start=1):
        decomposition_row = decomposition_df.iloc[k_idx]
        peaks.append({
            "rank": rank,
            "m": int(k_idx),
            "xi": float(xi_values[k_idx]),
            "kx": float(raw_k[k_idx, 0]),
            "ky": float(raw_k[k_idx, 1]),
            "kx_folded": float(folded_k[k_idx, 0]),
            "ky_folded": float(folded_k[k_idx, 1]),
            "kx_folded_over_2pi_a": float(folded_k[k_idx, 0] / (2.0 * np.pi / period)),
            "ky_folded_over_2pi_a": float(folded_k[k_idx, 1] / (2.0 * np.pi / period)),
            "weight": float(weights[k_idx]),
            "weight_fraction": float(weight_fraction[k_idx]),
            "dominant_subspace": decomposition_row["dominant_subspace"],
            "dominant_subspace_fraction": float(decomposition_row["dominant_subspace_fraction"]),
            "dominant_mode": decomposition_row["dominant_mode"],
            "dominant_mode_fraction": float(decomposition_row["dominant_mode_fraction"]),
            "subspace_matches_mode": bool(decomposition_row["subspace_matches_mode"]),
        })
    pd.DataFrame(peaks).to_csv(output_dir / "top_k_peaks.csv", index=False)

    np.savez_compressed(
        output_dir / "strip_bulk_fourier_hz.npz",
        mode_idx=int(mode_idx),
        kx_fraction=float(kx_fraction),
        rows=np.asarray(rows, dtype=int),
        residues=residues,
        row_centers=row_centers,
        rho_points=rho_points,
        xi_values=xi_values,
        k_raw=raw_k,
        k_folded=folded_k,
        row_fields=row_fields,
        row_fields_by_residue=row_fields_by_residue,
        F=f_components,
        G=g_components,
        P=p_components,
        profile_norms=profile_norms,
        weights=weights,
        weight_fraction=weight_fraction,
        reconstruction_error=reconstruction_error,
        parseval_error=parseval_error,
        profile_norm_error=profile_norm_error,
    )

    if generate_plots:
        plot_strip_mode_k_weights(output_dir / "k_weight_1d.png", mode_idx, eigen_row, folded_k, weight_fraction, period=period, dpi=dpi)
        plot_top_profiles(output_dir / "top_p_abs.png", mode_idx, rho_grid, p_components, weight_fraction, top_order, period=period, top_k=top_k, phase=False, title_prefix="strip bulk", dpi=dpi)
        plot_top_profiles(output_dir / "top_p_phase.png", mode_idx, rho_grid, p_components, weight_fraction, top_order, period=period, top_k=top_k, phase=True, title_prefix="strip bulk", dpi=dpi)

    top_idx = int(top_order[0])
    top_row = decomposition_df.iloc[top_idx]
    return {
        "mode_idx": int(mode_idx),
        "frequency": float(eigen_row["re"]) if "re" in eigen_row and np.isfinite(float(eigen_row["re"])) else np.nan,
        "q": float(eigen_row["q"]) if "q" in eigen_row and np.isfinite(float(eigen_row["q"])) else np.nan,
        "is_valid": bool(eigen_row["is_valid"]) if "is_valid" in eigen_row else True,
        "n_rows": int(n_rows),
        "n_rho": int(len(rho_points)),
        "kx_fraction": float(kx_fraction),
        "reconstruction_error": reconstruction_error,
        "parseval_error": parseval_error,
        "profile_norm_error": profile_norm_error,
        "top_m": top_idx,
        "top_kx_folded_over_2pi_a": float(folded_k[top_idx, 0] / (2.0 * np.pi / period)),
        "top_ky_folded_over_2pi_a": float(folded_k[top_idx, 1] / (2.0 * np.pi / period)),
        "top_weight_fraction": float(weight_fraction[top_idx]),
        "top_dominant_subspace": top_row["dominant_subspace"],
        "top_dominant_subspace_fraction": float(top_row["dominant_subspace_fraction"]),
        "top_dominant_mode": top_row["dominant_mode"],
        "top_dominant_mode_fraction": float(top_row["dominant_mode_fraction"]),
        "output_dir": str(output_dir),
    }


def selected_mode_indices(
    eigen_df: pd.DataFrame,
    mode_indices: list[int] | None,
    eigen_path: Path,
    *,
    respect_validity: bool = True,
) -> list[int]:
    if mode_indices is not None:
        return list(mode_indices)
    if not respect_validity:
        if "mode_idx" in eigen_df.columns:
            return eigen_df["mode_idx"].to_numpy(dtype=int).tolist()
        return eigen_df.index.to_numpy(dtype=int).tolist()
    if "is_valid" not in eigen_df.columns:
        raise ValueError(f"Missing is_valid column in {eigen_path}; rerun field export first.")
    valid_mask = eigen_df["is_valid"].map(csv_bool)
    if "mode_idx" in eigen_df.columns:
        selected = eigen_df.loc[valid_mask, "mode_idx"].to_numpy(dtype=int).tolist()
    else:
        selected = eigen_df.index[valid_mask].to_numpy(dtype=int).tolist()
    if not selected:
        raise ValueError(f"No valid modes found in {eigen_path}")
    return selected


def eigen_row_for_mode(eigen_df: pd.DataFrame, mode_idx: int, eigen_path: Path):
    if "mode_idx" in eigen_df.columns:
        matches = eigen_df.loc[eigen_df["mode_idx"] == mode_idx]
        if matches.empty:
            raise ValueError(f"Mode index {mode_idx} not found in {eigen_path}")
        return matches.iloc[0]
    return eigen_df.iloc[mode_idx]


def run_strip_bulk_fourier_postprocess(
    out_dir: Path,
    export_dir: Path,
    *,
    bulk_records: list[dict[str, Any]],
    bulk_radius: int,
    period: float,
    kx_fraction: float,
    rho_grid_size: int,
    top_k: int,
    mode_indices: list[int] | None = None,
    dpi: int = 220,
    generate_plots: bool = True,
    row_index_coefficients: tuple[int, int] = (0, 1),
) -> list[dict[str, float | int | str]]:
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(f"Missing strip eigenfrequency table: {eigen_path}")

    rows = strip_bulk_rows(
        bulk_records,
        bulk_radius,
        row_index_coefficients=row_index_coefficients,
    )
    row_centers = strip_bulk_row_centers(rows, period)
    rho_grid = strip_bulk_rho_grid(period, rho_grid_size)
    xi_values, raw_k, folded_k = strip_k_values(len(rows), kx_fraction, period)

    output_dir = strip_bulk_fourier_dir(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if generate_plots:
        plot_strip_bulk_rows(output_dir / "bulk_row_indices.png", bulk_records, rows, period=period, dpi=dpi)
        plot_strip_k_grid(output_dir / "bulk_k_grid_1d.png", folded_k, period=period, dpi=dpi)

    eigen_df = pd.read_csv(eigen_path)
    selected = selected_mode_indices(
        eigen_df,
        mode_indices,
        eigen_path,
        respect_validity=False,
    )
    print(f"Strip bulk Hz Fourier analyzed modes: {selected}")
    summary_rows = []
    for mode_idx in selected:
        summary_rows.append(strip_bulk_fourier_mode(
            export_dir,
            out_dir,
            mode_idx,
            eigen_row_for_mode(eigen_df, mode_idx, eigen_path),
            rho_grid,
            rows,
            row_centers,
            xi_values,
            raw_k,
            folded_k,
            kx_fraction,
            period=period,
            top_k=top_k,
            dpi=dpi,
            generate_plots=generate_plots,
        ))

    summary_path = output_dir / "strip_bulk_fourier_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"Wrote strip bulk Hz Fourier outputs to {output_dir}")
    return summary_rows


def plot_bulk_cyclic_index(output_path: Path, cells: list[LatticePoint], cyclic_indices: np.ndarray, *, bulk_radius: int, period: float, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    max_index = hex_quotient_cell_count(bulk_radius) - 1
    norm = plt.Normalize(0, max_index)
    cmap = plt.get_cmap("viridis")
    if len(cells) != len(cyclic_indices):
        raise ValueError("cells and cyclic_indices must have the same length")
    if not np.array_equal(
        np.sort(np.asarray(cyclic_indices, dtype=int)),
        np.arange(len(cells), dtype=int),
    ):
        raise ValueError("cyclic_indices must remain a complete FFT-index permutation")
    coords = np.array([[point.x, point.y] for point in cells], dtype=float)
    lattice_indices = hex_spiral_indices(coords)
    annotate = spiral_index_annotation_mask(lattice_indices)

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    for point, lattice_index, show_label in zip(
        cells,
        lattice_indices,
        annotate,
        strict=True,
    ):
        polygon = cell_polygon(point.x, point.y, period)
        ax.add_patch(Polygon(
            polygon,
            closed=True,
            facecolor=cmap(norm(lattice_index)),
            edgecolor="#111827",
            linewidth=0.45,
            alpha=0.86,
        ))
        if show_label:
            ax.text(
                point.x,
                point.y,
                str(int(lattice_index)),
                ha="center",
                va="center",
                fontsize=5.0,
                color="white",
            )

    pad = period
    ax.set_xlim(float(coords[:, 0].min() - pad), float(coords[:, 0].max() + pad))
    ax.set_ylim(float(coords[:, 1].min() - pad), float(coords[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"bulk lattice index, L={bulk_radius}")
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array(lattice_indices)
    colorbar_ax = make_axes_locatable(ax).append_axes("right", size="4.2%", pad=0.12)
    colorbar = fig.colorbar(
        scalar_mappable,
        cax=colorbar_ax,
        ticks=[0, max_index],
    )
    colorbar.set_label("lattice index", rotation=270, labelpad=18)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_k_grid(output_path: Path, folded_k: np.ndarray, bz_polygon: np.ndarray, *, period: float, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scale = 2.0 * np.pi / period
    bz = closed_line(bz_polygon) / scale
    points = folded_k / scale
    k_grid_indices = hex_spiral_indices(points)
    max_index = len(k_grid_indices) - 1
    norm = plt.Normalize(0, max_index)
    cmap = plt.get_cmap("viridis")
    annotate = spiral_index_annotation_mask(k_grid_indices)

    fig, ax = plt.subplots(figsize=(6.6, 5.9))
    ax.plot(bz[:, 0], bz[:, 1], color="#111827", linewidth=1.3, label="first BZ")
    scatter = ax.scatter(
        points[:, 0],
        points[:, 1],
        s=32,
        c=k_grid_indices,
        cmap=cmap,
        norm=norm,
        edgecolor="#111827",
        linewidth=0.25,
        zorder=3,
    )
    for k_grid_index, point, show_label in zip(
        k_grid_indices,
        points,
        annotate,
        strict=True,
    ):
        if show_label:
            ax.text(
                point[0],
                point[1],
                str(int(k_grid_index)),
                fontsize=5.0,
                ha="center",
                va="center",
                color="white",
                zorder=4,
            )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$k_x / (2\pi/a)$")
    ax.set_ylabel(r"$k_y / (2\pi/a)$")
    ax.set_title("finite k-grid index")
    ax.grid(True, alpha=0.18)
    ax.legend(loc="best", fontsize=8)
    colorbar_ax = make_axes_locatable(ax).append_axes("right", size="4.2%", pad=0.12)
    colorbar = fig.colorbar(
        scatter,
        cax=colorbar_ax,
        ticks=[0, max_index],
    )
    colorbar.set_label("k-grid index", rotation=270, labelpad=18)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_arbitrary_index_grid(
    output_path: Path,
    points: np.ndarray,
    boundary: np.ndarray | None,
    *,
    title: str,
    period: float,
    dpi: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scale = 2.0 * np.pi / period if boundary is not None else 1.0
    plotted_points = np.asarray(points, dtype=float) / scale
    indices = np.arange(len(plotted_points), dtype=int)
    fig, ax = plt.subplots(figsize=(6.6, 5.9))
    if boundary is not None:
        closed = closed_line(np.asarray(boundary, dtype=float)) / scale
        ax.plot(closed[:, 0], closed[:, 1], color="#111827", linewidth=1.3)
    scatter = ax.scatter(
        plotted_points[:, 0],
        plotted_points[:, 1],
        c=indices,
        s=28,
        cmap="viridis",
        edgecolor="#111827",
        linewidth=0.2,
    )
    ax.scatter(
        plotted_points[0, 0],
        plotted_points[0, 1],
        marker="x",
        s=44,
        color="#dc2626",
        linewidth=1.2,
        label="index 0",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.set_xlabel(r"$k_x/(2\pi/a)$" if boundary is not None else "x (um)")
    ax.set_ylabel(r"$k_y/(2\pi/a)$" if boundary is not None else "y (um)")
    ax.legend(loc="best", fontsize=8)
    fig.colorbar(scatter, ax=ax, label="arbitrary index")
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_finite_mode_k_weights(
    output_path: Path,
    eigen_row,
    folded_k: np.ndarray,
    bz_polygon: np.ndarray,
    weight_fraction: np.ndarray,
    rho_points: np.ndarray,
    profiles: np.ndarray,
    *,
    period: float,
    dpi: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scale = 2.0 * np.pi / period
    bz = closed_line(bz_polygon) / scale
    points = np.asarray(folded_k, dtype=float) / scale
    fractions = np.asarray(weight_fraction, dtype=float)
    rho = np.asarray(rho_points, dtype=float) / period
    mode_profiles = np.asarray(profiles)
    if len(fractions) == 0:
        raise ValueError("weight_fraction must be non-empty")
    if mode_profiles.ndim != 2 or mode_profiles.shape != (len(fractions), len(rho)):
        raise ValueError(
            "profiles must have shape (len(weight_fraction), len(rho_points))"
        )
    log_weight = np.log10(np.maximum(fractions, np.finfo(float).tiny))
    gamma_idx = int(np.argmax(fractions))

    fig, ax = plt.subplots(figsize=(10.5, 8.0))
    ax.set_facecolor("black")
    ax.plot(bz[:, 0], bz[:, 1], color="white", linewidth=1.35, zorder=2)
    scatter = ax.scatter(
        points[:, 0],
        points[:, 1],
        s=21.0,
        c=log_weight,
        cmap="magma",
        vmin=-5.0,
        vmax=0.0,
        linewidths=0.0,
        zorder=3,
    )

    profile = np.abs(mode_profiles[gamma_idx])
    profile /= max(float(np.max(profile)), np.finfo(float).tiny)
    inset = ax.inset_axes([0.81, 0.80, 0.18, 0.18])
    inset.set_facecolor("none")
    inset.patch.set_alpha(0.0)
    inset.tricontourf(
        rho[:, 0],
        rho[:, 1],
        profile,
        levels=np.linspace(0.0, 1.0, 101),
        cmap="plasma",
        vmin=0.0,
        vmax=1.0,
    )
    inset.set_aspect("equal", adjustable="box")
    inset.set_xticks([])
    inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_visible(False)
    ax.text(
        0.90,
        0.785,
        rf"$P(\Gamma)={100.0 * fractions[gamma_idx]:.1f}\%$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10.5,
        color="white",
    )

    title = r"k-weight ($H_z$)"
    if "re" in eigen_row and pd.notna(eigen_row["re"]):
        title += f", f={float(eigen_row['re']):.3f} THz"
    if "q" in eigen_row and pd.notna(eigen_row["q"]):
        title += f", Q={float(eigen_row['q']):.2f}"
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.73, 0.73)
    ax.set_ylim(-0.635, 0.635)
    ax.set_xlabel(r"$k_x/(2\pi/a)$")
    ax.set_ylabel(r"$k_y/(2\pi/a)$")
    ax.grid(False)

    divider = make_axes_locatable(ax)
    colorbar_ax = divider.append_axes("right", size="4.2%", pad=0.12)
    colorbar = fig.colorbar(
        scatter,
        cax=colorbar_ax,
        ticks=[-5, -4, -3, -2, -1, 0],
    )
    colorbar.ax.set_yticklabels([
        r"$10^{-5}$",
        r"$10^{-4}$",
        r"$10^{-3}$",
        r"$10^{-2}$",
        r"$10^{-1}$",
        r"$10^{0}$",
    ])
    colorbar.set_label(
        r"$P(k)=W(k) / \sum_k W(k)$",
        rotation=270,
        labelpad=20,
    )
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def center_normalized_cell_intensity(
    cells: list[LatticePoint],
    cell_fields: np.ndarray,
) -> np.ndarray:
    """Return each complete cell's mean |Hz|^2 normalized by the center cell."""
    fields = np.asarray(cell_fields, dtype=complex)
    if fields.ndim != 2 or fields.shape[0] != len(cells):
        raise ValueError(
            "cell_fields must have shape (len(cells), n_rho)"
        )
    if fields.shape[1] == 0:
        raise ValueError("cell_fields must contain at least one rho sample")
    center_positions = [
        idx
        for idx, point in enumerate(cells)
        if point.i == 0 and point.j == 0
    ]
    if len(center_positions) != 1:
        raise ValueError("cells must contain exactly one center cell (i, j)=(0, 0)")
    intensities = np.mean(np.abs(fields) ** 2, axis=1)
    if not np.isfinite(intensities).all():
        raise ValueError("cell intensities must be finite")
    center_intensity = float(intensities[center_positions[0]])
    if center_intensity <= 0.0:
        raise ValueError("center-cell intensity must be positive")
    return intensities / center_intensity


def plot_center_normalized_cell_intensity(
    output_path: Path,
    cells: list[LatticePoint],
    cell_fields: np.ndarray,
    *,
    period: float,
    dpi: int,
) -> None:
    """Plot the standard full finite unit-cell Hz intensity map."""
    if not cells:
        raise ValueError("cells must not be empty")
    lattice_keys = [(point.i, point.j) for point in cells]
    if len(lattice_keys) != len(set(lattice_keys)):
        raise ValueError("cells contain duplicate lattice indices")
    normalized = center_normalized_cell_intensity(cells, cell_fields)
    norm = TwoSlopeNorm(
        vmin=CENTER_NORMALIZED_INTENSITY_VMIN,
        vcenter=CENTER_NORMALIZED_INTENSITY_VCENTER,
        vmax=CENTER_NORMALIZED_INTENSITY_VMAX,
    )
    cmap = plt.get_cmap(CENTER_NORMALIZED_INTENSITY_CMAP)
    patches = [
        Polygon(cell_polygon(point.x, point.y, period), closed=True)
        for point in cells
    ]

    fig, ax = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    collection = PatchCollection(
        patches,
        cmap=cmap,
        norm=norm,
        edgecolor="none",
    )
    collection.set_array(
        np.clip(
            normalized,
            CENTER_NORMALIZED_INTENSITY_VMIN,
            CENTER_NORMALIZED_INTENSITY_VMAX,
        )
    )
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(CENTER_NORMALIZED_INTENSITY_TITLE)

    colorbar_ax = make_axes_locatable(ax).append_axes(
        "right",
        size="4.2%",
        pad=0.12,
    )
    colorbar = fig.colorbar(collection, cax=colorbar_ax)
    colorbar.set_ticks(CENTER_NORMALIZED_INTENSITY_TICKS)
    colorbar.set_ticklabels(CENTER_NORMALIZED_INTENSITY_TICK_LABELS)
    colorbar.set_label(
        CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL,
        rotation=270,
        labelpad=20,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def sample_finite_cell_fields(
    export_dir: Path,
    mode_idx: int,
    rho_grid: dict[str, np.ndarray],
    cells: list[LatticePoint],
) -> np.ndarray:
    """Sample one exported Hz field over the requested complete cells."""
    field_path = export_dir / f"{mode_idx:02d}_Hz_center.parquet"
    if not field_path.exists():
        raise FileNotFoundError(f"Missing finite Hz field parquet: {field_path}")
    if not cells:
        raise ValueError("cells must not be empty")
    field_points, field = read_hz_center_field(field_path)
    rho_points = np.asarray(rho_grid["points"], dtype=float)
    cell_centers = np.asarray(
        [[point.x, point.y] for point in cells],
        dtype=float,
    )
    query_points = (
        cell_centers[:, None, :] + rho_points[None, :, :]
    ).reshape(-1, 2)
    return interpolate_complex_field(
        field_points,
        field,
        query_points,
    ).reshape(len(cells), len(rho_points))


def finite_mode_envelope_metrics(
    cell_fields: np.ndarray,
    cells: list[LatticePoint],
    *,
    period: float,
    significant_envelope_fraction: float = 0.05,
) -> dict[str, object]:
    """Measure the radial finite-cavity envelope of one internal mode."""
    fields = np.asarray(cell_fields, dtype=complex)
    if fields.ndim != 2 or fields.shape[0] != len(cells):
        raise ValueError(
            "cell_fields must have shape (len(cells), n_rho)"
        )
    if fields.shape[1] == 0:
        raise ValueError("cell_fields must contain rho samples")
    if not np.all(np.isfinite(fields)):
        raise ValueError("cell_fields must contain only finite values")
    significant_envelope_fraction = float(significant_envelope_fraction)
    if not 0.0 < significant_envelope_fraction < 1.0:
        raise ValueError(
            "significant_envelope_fraction must be between 0 and 1"
        )

    shells = np.asarray([int(point.shell) for point in cells], dtype=int)
    if np.any(shells < 0):
        raise ValueError("cell shells must be non-negative")
    shell_count = int(shells.max()) + 1
    if set(shells.tolist()) != set(range(shell_count)):
        raise ValueError("cell shells must form a contiguous range from zero")
    center_positions = np.flatnonzero(shells == 0)
    if len(center_positions) != 1:
        raise ValueError("cells must contain exactly one center cell")
    center_idx = int(center_positions[0])

    cell_intensity = np.mean(np.abs(fields) ** 2, axis=1)
    total_intensity = float(np.sum(cell_intensity))
    if not np.isfinite(total_intensity) or total_intensity <= 0.0:
        raise ValueError("finite mode has zero or non-finite cell intensity")
    shell_intensity = np.asarray(
        [
            float(np.mean(cell_intensity[shells == shell]))
            for shell in range(shell_count)
        ],
        dtype=float,
    )
    peak_shell = int(np.argmax(shell_intensity))
    peak_cell_intensity = float(np.max(cell_intensity))
    peak_shell_intensity = float(np.max(shell_intensity))
    normalized_shell_intensity = shell_intensity / peak_shell_intensity
    outward_resurgence = float(
        np.sum(
            np.maximum(
                np.diff(normalized_shell_intensity[peak_shell:]),
                0.0,
            )
        )
    )

    reference = fields[center_idx]
    reference_norm = float(np.vdot(reference, reference).real)
    if not np.isfinite(reference_norm) or reference_norm <= 0.0:
        raise ValueError("finite mode has a zero center-cell profile")
    cell_envelope = fields @ np.conj(reference) / reference_norm
    shell_envelope = np.asarray(
        [
            np.mean(cell_envelope[shells == shell])
            for shell in range(shell_count)
        ],
        dtype=complex,
    )
    if abs(shell_envelope[0]) <= np.finfo(float).tiny:
        raise ValueError("finite mode has a zero center-shell envelope")
    shell_envelope = shell_envelope / shell_envelope[0]
    significant = np.abs(shell_envelope) >= significant_envelope_fraction
    significant_real = np.real(shell_envelope[significant])
    nonzero_signs = np.sign(significant_real)
    nonzero_signs = nonzero_signs[nonzero_signs != 0.0]
    radial_sign_changes = int(
        np.sum(nonzero_signs[1:] * nonzero_signs[:-1] < 0.0)
    )
    max_phase_degrees = float(
        np.max(np.abs(np.angle(shell_envelope[significant], deg=True)))
    )

    centers = np.asarray(
        [[float(point.x), float(point.y)] for point in cells],
        dtype=float,
    )
    radius2_over_a2 = np.sum(centers**2, axis=1) / float(period) ** 2
    rms_radius_over_a = float(
        np.sqrt(
            np.sum(cell_intensity * radius2_over_a2) / total_intensity
        )
    )
    return {
        "center_over_peak_cell": float(
            cell_intensity[center_idx] / peak_cell_intensity
        ),
        "center_over_peak_shell": float(
            shell_intensity[0] / peak_shell_intensity
        ),
        "peak_shell": peak_shell,
        "central_energy_shell_le_2": float(
            np.sum(cell_intensity[shells <= 2]) / total_intensity
        ),
        "rms_radius_over_a": rms_radius_over_a,
        "outward_resurgence": outward_resurgence,
        "radial_sign_changes": radial_sign_changes,
        "max_significant_shell_phase_degrees": max_phase_degrees,
        "significant_envelope_fraction": significant_envelope_fraction,
        "radial_shell_intensity_normalized": json.dumps(
            normalized_shell_intensity.tolist()
        ),
        "radial_shell_envelope_real": json.dumps(
            np.real(shell_envelope).tolist()
        ),
        "radial_shell_envelope_imag": json.dumps(
            np.imag(shell_envelope).tolist()
        ),
    }


def select_finite_fundamental_mode(
    export_dir: Path,
    output_path: Path,
    *,
    cells: list[LatticePoint],
    period: float,
    rho_grid_size: int,
    target_internal_mode: str,
    target_frequency_thz: float,
    candidate_mode_indices: list[int] | None = None,
    target_component_minimum: float = 0.5,
    center_over_peak_shell_minimum: float = 0.9,
    outward_resurgence_maximum: float = 0.05,
    significant_envelope_fraction: float = 0.05,
) -> tuple[int, pd.DataFrame]:
    """Select the centered, nodeless envelope in one symmetry sector.

    The Gamma component verifies the internal ``px/py/dx/dy`` identity.  The
    finite-cavity fundamental is then determined from the complex radial
    envelope, with target frequency used only as a late tie-break.
    """
    component_columns = {
        "px": "mode_px",
        "py": "mode_py",
        "dx": "mode_dx",
        "dy": "mode_dy",
    }
    target_internal_mode = str(target_internal_mode).strip().lower()
    if target_internal_mode not in component_columns:
        raise ValueError(
            "target_internal_mode must be one of px, py, dx, or dy"
        )
    target_frequency_thz = float(target_frequency_thz)
    if not np.isfinite(target_frequency_thz):
        raise ValueError("target_frequency_thz must be finite")
    if not cells:
        raise ValueError("cells must not be empty")
    target_component_minimum = float(target_component_minimum)
    center_over_peak_shell_minimum = float(
        center_over_peak_shell_minimum
    )
    outward_resurgence_maximum = float(outward_resurgence_maximum)
    if not 0.0 <= target_component_minimum <= 1.0:
        raise ValueError("target_component_minimum must be between 0 and 1")
    if not 0.0 <= center_over_peak_shell_minimum <= 1.0:
        raise ValueError(
            "center_over_peak_shell_minimum must be between 0 and 1"
        )
    if outward_resurgence_maximum < 0.0:
        raise ValueError("outward_resurgence_maximum must be non-negative")

    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(
            f"Missing finite eigenfrequency table: {eigen_path}"
        )
    eigen_df = pd.read_csv(eigen_path).drop(
        columns=["fundamental_envelope_score"],
        errors="ignore",
    )
    valid_modes = selected_mode_indices(eigen_df, None, eigen_path)
    if candidate_mode_indices is None:
        candidates = valid_modes
    else:
        candidates = [int(mode_idx) for mode_idx in candidate_mode_indices]
        invalid = sorted(set(candidates).difference(valid_modes))
        if invalid:
            raise ValueError(
                "Gamma-component candidates must be valid exported modes; "
                f"invalid candidates: {invalid}"
            )
        if len(candidates) != len(set(candidates)):
            raise ValueError("Gamma-component candidates must be unique")
    if not candidates:
        raise ValueError("Gamma-component selection has no candidate modes")

    rho_grid = unit_cell_rho_grid(period, rho_grid_size)
    rho_points = np.asarray(rho_grid["points"], dtype=float)
    rows: list[dict[str, object]] = []
    for mode_idx in candidates:
        eigen_row = eigen_row_for_mode(eigen_df, mode_idx, eigen_path)
        frequency = float(eigen_row["re"])
        if not np.isfinite(frequency):
            raise ValueError(f"Mode {mode_idx} has non-finite frequency")
        cell_fields = sample_finite_cell_fields(
            export_dir,
            mode_idx,
            rho_grid,
            cells,
        )
        envelope_metrics = finite_mode_envelope_metrics(
            cell_fields,
            cells,
            period=period,
            significant_envelope_fraction=(
                significant_envelope_fraction
            ),
        )
        gamma_profile = np.sum(cell_fields, axis=0) / np.sqrt(
            float(len(cells))
        )
        gamma_weight = float(np.sum(np.abs(gamma_profile) ** 2))
        if not np.isfinite(gamma_weight) or gamma_weight <= 0.0:
            raise ValueError(
                f"Mode {mode_idx} has a zero or non-finite Gamma profile"
            )
        decomposition = decompose_profiles(
            rho_points,
            gamma_profile[None, :],
            np.array([gamma_weight], dtype=float),
            np.array([1.0], dtype=float),
            "t",
        ).iloc[0]
        row = {
            "mode_idx": int(mode_idx),
            "frequency": frequency,
            "q": (
                float(eigen_row["q"])
                if "q" in eigen_row and pd.notna(eigen_row["q"])
                else np.nan
            ),
            "target_frequency_thz": target_frequency_thz,
            "frequency_distance_to_target_thz": abs(
                frequency - target_frequency_thz
            ),
            "target_internal_mode": target_internal_mode,
            "target_component_column": component_columns[
                target_internal_mode
            ],
            "gamma_profile_weight": gamma_weight,
        }
        for column in (
            "subspace_s",
            "subspace_p",
            "subspace_d",
            "subspace_f",
            "mode_s",
            "mode_px",
            "mode_py",
            "mode_dx",
            "mode_dy",
            "mode_f",
            "dominant_subspace",
            "dominant_subspace_fraction",
            "dominant_mode",
            "dominant_mode_fraction",
            "subspace_matches_mode",
        ):
            row[column] = decomposition[column]
        row["target_component_weight"] = float(
            decomposition[component_columns[target_internal_mode]]
        )
        row.update(envelope_metrics)
        row["target_component_matches"] = bool(
            row["dominant_mode"] == target_internal_mode
            and row["target_component_weight"]
            >= target_component_minimum
        )
        row["fundamental_eligible"] = bool(
            row["target_component_matches"]
            and int(row["peak_shell"]) == 0
            and float(row["center_over_peak_shell"])
            >= center_over_peak_shell_minimum
            and int(row["radial_sign_changes"]) == 0
            and float(row["outward_resurgence"])
            <= outward_resurgence_maximum
        )
        rows.append(row)

    selection_df = pd.DataFrame(rows)
    eligible = selection_df.loc[selection_df["fundamental_eligible"]].copy()
    selection_df["analysis_selected"] = False
    selection_df["selection_reason"] = "rejected_higher_order_envelope"
    if eligible.empty:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        selection_df.sort_values("mode_idx").to_csv(
            output_path,
            index=False,
        )
        raise ValueError(
            "No valid centered nodeless fundamental mode found; "
            f"inspect {output_path}"
        )
    eligible = eligible.sort_values(
        [
            "radial_sign_changes",
            "peak_shell",
            "outward_resurgence",
            "center_over_peak_shell",
            "target_component_weight",
            "frequency_distance_to_target_thz",
            "mode_idx",
        ],
        ascending=[True, True, True, False, False, True, True],
    )
    selected_mode = int(eligible.iloc[0]["mode_idx"])
    selection_df["analysis_selected"] = selection_df["mode_idx"].eq(
        selected_mode
    )
    selection_df.loc[
        selection_df["analysis_selected"], "selection_reason"
    ] = "selected_centered_nodeless_fundamental_envelope"
    selection_df = selection_df.sort_values("mode_idx")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    selection_df.to_csv(output_path, index=False)

    eigen_df["analysis_selected"] = eigen_df["mode_idx"].astype(int).eq(
        selected_mode
    )
    eigen_df["analysis_selection_method"] = (
        "target_component_centered_nodeless_envelope"
    )
    eigen_df["analysis_target_internal_mode"] = target_internal_mode
    indexed_selection = selection_df.set_index("mode_idx")
    for column in (
        "frequency_distance_to_target_thz",
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
        eigen_df[column] = eigen_df["mode_idx"].astype(int).map(
            indexed_selection[column]
        )
    eigen_df.to_csv(eigen_path, index=False)
    selected_row = selection_df.loc[
        selection_df["analysis_selected"]
    ].iloc[0]
    print(
        "Selected finite fundamental mode: "
        f"mode {selected_mode}, {target_internal_mode}="
        f"{float(selected_row['target_component_weight']):.9g}, "
        f"radial nodes={int(selected_row['radial_sign_changes'])}, "
        f"f={float(selected_row['frequency']):.9g} THz "
        f"(target {target_frequency_thz:.9g} THz)",
        flush=True,
    )
    return selected_mode, selection_df


def select_finite_mode_closest_to_frequency(
    export_dir: Path,
    output_path: Path,
    *,
    target_frequency_thz: float,
    candidate_mode_indices: list[int] | None = None,
) -> tuple[int, pd.DataFrame]:
    """Select the valid finite mode closest to the configured target frequency."""
    target_frequency_thz = float(target_frequency_thz)
    if not np.isfinite(target_frequency_thz):
        raise ValueError("target_frequency_thz must be finite")
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(f"Missing finite eigenfrequency table: {eigen_path}")
    eigen_df = pd.read_csv(eigen_path)
    eigen_df = eigen_df.drop(
        columns=["fundamental_envelope_score"],
        errors="ignore",
    )
    valid_modes = selected_mode_indices(eigen_df, None, eigen_path)
    if candidate_mode_indices is None:
        candidates = valid_modes
    else:
        candidates = [int(mode_idx) for mode_idx in candidate_mode_indices]
        invalid = sorted(set(candidates).difference(valid_modes))
        if invalid:
            raise ValueError(
                "Single-mode candidates must be valid exported modes; "
                f"invalid candidates: {invalid}"
            )
        if len(candidates) != len(set(candidates)):
            raise ValueError("Single-mode candidates must be unique")
    if not candidates:
        raise ValueError("Single-mode selection has no candidate modes")

    rows: list[dict[str, object]] = []
    for mode_idx in candidates:
        eigen_row = eigen_row_for_mode(eigen_df, mode_idx, eigen_path)
        frequency = float(eigen_row["re"])
        if not np.isfinite(frequency):
            raise ValueError(f"Mode {mode_idx} has non-finite frequency")
        rows.append(
            {
                "mode_idx": mode_idx,
                "frequency": frequency,
                "q": (
                    float(eigen_row["q"])
                    if "q" in eigen_row and pd.notna(eigen_row["q"])
                    else np.nan
                ),
                "target_frequency_thz": target_frequency_thz,
                "frequency_distance_to_target_thz": abs(
                    frequency - target_frequency_thz
                ),
            }
        )

    selection_df = pd.DataFrame(rows).sort_values(
        ["frequency_distance_to_target_thz", "mode_idx"],
        ascending=[True, True],
    )
    selected_mode = int(selection_df.iloc[0]["mode_idx"])
    selection_df["analysis_selected"] = selection_df["mode_idx"].eq(
        selected_mode
    )
    selection_df["selection_reason"] = "not_closest_to_target_frequency"
    selection_df.loc[
        selection_df["analysis_selected"], "selection_reason"
    ] = "selected_closest_to_target_frequency"
    selection_df = selection_df.sort_values("mode_idx")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    selection_df.to_csv(output_path, index=False)

    eigen_df["analysis_selected"] = eigen_df["mode_idx"].astype(int).eq(
        selected_mode
    )
    distance_by_mode = selection_df.set_index("mode_idx")[
        "frequency_distance_to_target_thz"
    ]
    eigen_df["frequency_distance_to_target_thz"] = (
        eigen_df["mode_idx"].astype(int).map(distance_by_mode)
    )
    eigen_df.to_csv(eigen_path, index=False)
    print(
        "Selected finite mode closest to target frequency: "
        f"mode {selected_mode} at "
        f"{selection_df.loc[selection_df['analysis_selected'], 'frequency'].iloc[0]:.9g} "
        f"THz (target {target_frequency_thz:.9g} THz)",
        flush=True,
    )
    return selected_mode, selection_df


def sample_finite_unit_cell_fields(
    export_dir: Path,
    mode_idx: int,
    rho_grid: dict[str, np.ndarray],
    cells: list[LatticePoint],
    cladding_cells: list[LatticePoint],
) -> tuple[list[LatticePoint], np.ndarray]:
    """Sample one exported Hz field over every cavity and cladding cell."""
    if not cladding_cells:
        raise ValueError("cladding_cells must not be empty")
    all_cells = [*cells, *cladding_cells]
    return all_cells, sample_finite_cell_fields(
        export_dir,
        mode_idx,
        rho_grid,
        all_cells,
    )


def write_finite_unit_cell_intensity_map(
    export_dir: Path,
    output_root: Path,
    mode_idx: int,
    rho_grid: dict[str, np.ndarray],
    cells: list[LatticePoint],
    cladding_cells: list[LatticePoint],
    *,
    period: float,
    dpi: int,
) -> Path:
    """Write the standard center-normalized map for one finite mode."""
    all_cells, all_cell_fields = sample_finite_unit_cell_fields(
        export_dir,
        mode_idx,
        rho_grid,
        cells,
        cladding_cells,
    )
    output_path = (
        finite_lattice_fourier_dir(output_root)
        / f"mode_{mode_idx:02d}"
        / INTENSITY_MAPS_NORMH_FILENAME
    )
    plot_center_normalized_cell_intensity(
        output_path,
        all_cells,
        all_cell_fields,
        period=period,
        dpi=dpi,
    )
    return output_path


def run_finite_unit_cell_intensity_maps(
    out_dir: Path,
    export_dir: Path,
    *,
    cells: list[LatticePoint],
    cladding_cells: list[LatticePoint],
    period: float,
    rho_grid_size: int,
    mode_indices: list[int] | None = None,
    dpi: int = 220,
) -> list[Path]:
    """Write the standard map for all selected valid finite modes."""
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(f"Missing finite eigenfrequency table: {eigen_path}")
    eigen_df = pd.read_csv(eigen_path)
    selected = selected_mode_indices(eigen_df, mode_indices, eigen_path)
    rho_grid = unit_cell_rho_grid(period, rho_grid_size)
    outputs = [
        write_finite_unit_cell_intensity_map(
            export_dir,
            out_dir,
            mode_idx,
            rho_grid,
            cells,
            cladding_cells,
            period=period,
            dpi=dpi,
        )
        for mode_idx in selected
    ]
    print(f"Finite unit-cell Hz intensity maps: {selected}")
    return outputs


def finite_lattice_fourier_mode(
    export_dir: Path,
    output_root: Path,
    mode_idx: int,
    eigen_row,
    rho_grid: dict[str, np.ndarray],
    cells: list[LatticePoint],
    cladding_cells: list[LatticePoint],
    cyclic_indices: np.ndarray,
    xi_values: np.ndarray,
    folded_k: np.ndarray,
    bz_polygon: np.ndarray,
    *,
    bulk_radius: int,
    period: float,
    top_k: int,
    dpi: int,
    transform_kind: str = "hex_cyclic_quotient_v1",
) -> dict[str, float | int | str]:
    rho_points = np.asarray(rho_grid["points"], dtype=float)
    all_cells, all_cell_fields = sample_finite_unit_cell_fields(
        export_dir,
        mode_idx,
        rho_grid,
        cells,
        cladding_cells,
    )
    all_cell_centers = np.array(
        [[point.x, point.y] for point in all_cells],
        dtype=float,
    )
    cell_centers = all_cell_centers[:len(cells)]
    if transform_kind == "arbitrary_square_cavity_v1":
        lattice_indices = np.arange(len(cell_centers), dtype=int)
        k_grid_indices = np.arange(len(folded_k), dtype=int)
    else:
        lattice_indices = hex_spiral_indices(cell_centers)
        k_grid_indices = hex_spiral_indices(folded_k)
    cell_fields = all_cell_fields[:len(cells)]

    n_cells = len(cells)
    if transform_kind == "arbitrary_square_cavity_v1":
        fourier = arbitrary_finite_lattice_fourier_components(
            cell_fields,
            cell_centers,
            rho_points,
            folded_k,
            period,
        )
    else:
        expected_cells = hex_quotient_cell_count(bulk_radius)
        if n_cells != expected_cells:
            raise ValueError(f"Expected {expected_cells} bulk cells, got {n_cells}")
        if set(cyclic_indices.tolist()) != set(range(n_cells)):
            raise ValueError("Bulk cells are not a complete hex quotient representative set.")
        fields_by_x = np.empty_like(cell_fields)
        fields_by_x[cyclic_indices] = cell_fields
        fourier = finite_lattice_fourier_components(
            fields_by_x, rho_points, xi_values, period
        )
    f_components = fourier["F"]
    g_components = fourier["G"]
    p_components = fourier["P"]
    profile_norms = fourier["profile_norms"]
    weights = fourier["weights"]
    weight_fraction = fourier["weight_fraction"]
    reconstruction_error = float(fourier["reconstruction_error"])
    parseval_error = float(fourier["parseval_error"])
    profile_norm_error = float(fourier["profile_norm_error"])
    top_order = np.argsort(weight_fraction)[::-1]

    output_dir = finite_lattice_fourier_dir(output_root) / f"mode_{mode_idx:02d}"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_center_normalized_cell_intensity(
        output_dir / INTENSITY_MAPS_NORMH_FILENAME,
        all_cells,
        all_cell_fields,
        period=period,
        dpi=dpi,
    )

    decomposition_df = decompose_profiles(rho_points, p_components, weights, weight_fraction, "t")
    decomposition_df["k_index"] = k_grid_indices
    decomposition_df["xi_1"] = xi_values[:, 0]
    decomposition_df["xi_2"] = xi_values[:, 1]
    decomposition_df["kx_folded"] = folded_k[:, 0]
    decomposition_df["ky_folded"] = folded_k[:, 1]
    decomposition_df["kx_folded_over_2pi_a"] = folded_k[:, 0] / (2.0 * np.pi / period)
    decomposition_df["ky_folded_over_2pi_a"] = folded_k[:, 1] / (2.0 * np.pi / period)
    decomposition_df.to_csv(output_dir / "p_subspace_mode_decomposition.csv", index=False)

    peaks = []
    for rank, k_idx in enumerate(top_order[:12], start=1):
        decomposition_row = decomposition_df.iloc[k_idx]
        peaks.append({
            "rank": rank,
            "t": int(k_idx),
            "k_index": int(k_grid_indices[k_idx]),
            "xi_1": float(xi_values[k_idx, 0]),
            "xi_2": float(xi_values[k_idx, 1]),
            "kx_folded": float(folded_k[k_idx, 0]),
            "ky_folded": float(folded_k[k_idx, 1]),
            "kx_folded_over_2pi_a": float(folded_k[k_idx, 0] / (2.0 * np.pi / period)),
            "ky_folded_over_2pi_a": float(folded_k[k_idx, 1] / (2.0 * np.pi / period)),
            "weight": float(weights[k_idx]),
            "weight_fraction": float(weight_fraction[k_idx]),
            "dominant_subspace": decomposition_row["dominant_subspace"],
            "dominant_subspace_fraction": float(decomposition_row["dominant_subspace_fraction"]),
            "dominant_mode": decomposition_row["dominant_mode"],
            "dominant_mode_fraction": float(decomposition_row["dominant_mode_fraction"]),
            "subspace_matches_mode": bool(decomposition_row["subspace_matches_mode"]),
        })
    pd.DataFrame(peaks).to_csv(output_dir / "top_k_peaks.csv", index=False)

    np.savez_compressed(
        output_dir / "finite_lattice_fourier_hz.npz",
        mode_idx=int(mode_idx),
        cell_indices=np.array([[point.i, point.j] for point in cells], dtype=int),
        cell_centers=cell_centers,
        lattice_indices=lattice_indices,
        cyclic_indices=cyclic_indices,
        rho_points=rho_points,
        xi_values=xi_values,
        k_folded=folded_k,
        k_grid_indices=k_grid_indices,
        F=f_components,
        G=g_components,
        P=p_components,
        profile_norms=profile_norms,
        weights=weights,
        weight_fraction=weight_fraction,
        reconstruction_error=reconstruction_error,
        parseval_error=parseval_error,
        profile_norm_error=profile_norm_error,
        transform_kind=transform_kind,
        index_kind=transform_kind,
    )

    (output_dir / "k_weight_first_bz.png").unlink(missing_ok=True)
    plot_finite_mode_k_weights(
        output_dir / "k_weight_Hz.png",
        eigen_row,
        folded_k,
        bz_polygon,
        weight_fraction,
        rho_points,
        p_components,
        period=period,
        dpi=dpi,
    )
    for stale_name in ("top_p_abs.png", "top_p_phase.png"):
        (output_dir / stale_name).unlink(missing_ok=True)
    plot_top_profiles(
        output_dir / "k_weight_tops_norm.png",
        mode_idx,
        rho_grid,
        p_components,
        weight_fraction,
        top_order,
        period=period,
        top_k=top_k,
        phase=False,
        title_prefix="finite-lattice",
        dpi=dpi,
        display_indices=k_grid_indices,
    )
    plot_top_profiles(
        output_dir / "k_weight_tops_phase.png",
        mode_idx,
        rho_grid,
        p_components,
        weight_fraction,
        top_order,
        period=period,
        top_k=top_k,
        phase=True,
        title_prefix="finite-lattice",
        dpi=dpi,
        display_indices=k_grid_indices,
    )

    top_idx = int(top_order[0])
    top_row = decomposition_df.iloc[top_idx]
    return {
        "mode_idx": int(mode_idx),
        "frequency": float(eigen_row["re"]) if "re" in eigen_row and pd.notna(eigen_row["re"]) else np.nan,
        "q": float(eigen_row["q"]) if "q" in eigen_row and pd.notna(eigen_row["q"]) else np.nan,
        "is_valid": bool(eigen_row["is_valid"]) if "is_valid" in eigen_row else True,
        "transform_kind": transform_kind,
        "index_kind": transform_kind,
        "n_cells": int(n_cells),
        "n_rho": int(len(rho_points)),
        "reconstruction_error": reconstruction_error,
        "parseval_error": float(parseval_error),
        "profile_norm_error": profile_norm_error,
        "top_t": top_idx,
        "top_k_index": int(k_grid_indices[top_idx]),
        "top_kx_folded_over_2pi_a": float(folded_k[top_idx, 0] / (2.0 * np.pi / period)),
        "top_ky_folded_over_2pi_a": float(folded_k[top_idx, 1] / (2.0 * np.pi / period)),
        "top_weight_fraction": float(weight_fraction[top_idx]),
        "top_dominant_subspace": top_row["dominant_subspace"],
        "top_dominant_subspace_fraction": float(top_row["dominant_subspace_fraction"]),
        "top_dominant_mode": top_row["dominant_mode"],
        "top_dominant_mode_fraction": float(top_row["dominant_mode_fraction"]),
        "output_dir": str(output_dir),
    }


def run_finite_lattice_fourier_postprocess(
    out_dir: Path,
    export_dir: Path,
    *,
    cells: list[LatticePoint],
    cladding_cells: list[LatticePoint],
    bulk_radius: int,
    period: float,
    rho_grid_size: int,
    top_k: int,
    mode_indices: list[int] | None = None,
    include_unselected_valid_mode_maps: bool = True,
    dpi: int = 220,
    transform_kind: str = "hex_cyclic_quotient_v1",
) -> list[dict[str, float | int | str]]:
    eigen_path = export_dir / "eigenfrequencies.csv"
    if not eigen_path.exists():
        raise FileNotFoundError(f"Missing finite eigenfrequency table: {eigen_path}")

    output_dir = finite_lattice_fourier_dir(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if transform_kind == "arbitrary_square_cavity_v1":
        side = 2 * bulk_radius + 1
        axis = np.arange(-bulk_radius, bulk_radius + 1, dtype=float) / float(side)
        xi_1, xi_2 = np.meshgrid(axis, axis, indexing="xy")
        xi_values = np.column_stack([xi_1.ravel(), xi_2.ravel()])
        gamma_idx = int(np.flatnonzero(np.all(np.isclose(xi_values, 0.0), axis=1))[0])
        xi_values = np.concatenate(
            [xi_values[gamma_idx : gamma_idx + 1], np.delete(xi_values, gamma_idx, axis=0)],
            axis=0,
        )
        cyclic_indices = np.array([], dtype=int)
    elif transform_kind == "hex_cyclic_quotient_v1":
        indices = np.array([[point.i, point.j] for point in cells], dtype=int)
        cyclic_indices = hex_cyclic_indices(indices, bulk_radius)
        xi_values = hex_xi_values(bulk_radius)
    else:
        raise ValueError(
            "transform_kind must be 'hex_cyclic_quotient_v1' or "
            "'arbitrary_square_cavity_v1'"
        )
    folded_k = fold_xis_to_first_bz(xi_values, period)
    bz_polygon = first_bz_polygon(period)
    rho_grid = unit_cell_rho_grid(period, rho_grid_size)

    if transform_kind == "hex_cyclic_quotient_v1":
        plot_bulk_cyclic_index(output_dir / "index_r_cavity.png", cells, cyclic_indices, bulk_radius=bulk_radius, period=period, dpi=dpi)
        plot_k_grid(output_dir / "index_k_1stBZ.png", folded_k, bz_polygon, period=period, dpi=dpi)
    else:
        cell_centers = np.array([[point.x, point.y] for point in cells], dtype=float)
        plot_arbitrary_index_grid(
            output_dir / "index_r_cavity.png",
            cell_centers,
            None,
            title="square cavity analysis-cell index",
            period=period,
            dpi=dpi,
        )
        plot_arbitrary_index_grid(
            output_dir / "index_k_1stBZ.png",
            folded_k,
            bz_polygon,
            title="arbitrary finite-lattice k grid",
            period=period,
            dpi=dpi,
        )

    eigen_df = pd.read_csv(eigen_path)
    selected = selected_mode_indices(eigen_df, mode_indices, eigen_path)
    all_valid = selected_mode_indices(eigen_df, None, eigen_path)
    map_only_modes = (
        [mode_idx for mode_idx in all_valid if mode_idx not in selected]
        if include_unselected_valid_mode_maps
        else []
    )
    for mode_idx in map_only_modes:
        write_finite_unit_cell_intensity_map(
            export_dir,
            out_dir,
            mode_idx,
            rho_grid,
            cells,
            cladding_cells,
            period=period,
            dpi=dpi,
        )
    print(f"Finite-lattice Hz Fourier valid modes: {selected}")
    summary_rows = []
    for mode_idx in selected:
        summary_rows.append(finite_lattice_fourier_mode(
            export_dir,
            out_dir,
            mode_idx,
            eigen_row_for_mode(eigen_df, mode_idx, eigen_path),
            rho_grid,
            cells,
            cladding_cells,
            cyclic_indices,
            xi_values,
            folded_k,
            bz_polygon,
            bulk_radius=bulk_radius,
            period=period,
            top_k=top_k,
            dpi=dpi,
            transform_kind=transform_kind,
        ))

    summary_path = output_dir / "finite_lattice_fourier_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"Wrote finite-lattice Hz Fourier outputs to {output_dir}")
    return summary_rows


def score_gamma_modes(
    out_dir: Path,
    *,
    fourier_dir: Path,
    summary_filename: str,
    component_column: str,
    decomposition_dir_prefix: str = "mode_",
    export_dirname: str = "simulation_exports",
    respect_summary_validity: bool = True,
    gamma_subspace_p_valid_threshold: float | None = None,
) -> tuple[float, pd.DataFrame]:
    summary_path = fourier_dir / summary_filename
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing Fourier summary: {summary_path}")

    eigen_path = simulation_export_dir(out_dir, export_dirname) / "eigenfrequencies.csv"
    valid_mode_count = None
    if eigen_path.exists():
        eigen_df = pd.read_csv(eigen_path)
        if "is_valid" in eigen_df.columns:
            valid_mode_count = int(eigen_df["is_valid"].map(csv_bool).sum())

    summary_df = pd.read_csv(summary_path)
    rows = []
    for _, summary_row in summary_df.iterrows():
        if (
            respect_summary_validity
            and "is_valid" in summary_row
            and not csv_bool(summary_row["is_valid"])
        ):
            continue
        mode_idx = int(summary_row["mode_idx"])
        decomposition_path = fourier_dir / f"{decomposition_dir_prefix}{mode_idx:02d}" / "p_subspace_mode_decomposition.csv"
        if not decomposition_path.exists():
            raise FileNotFoundError(f"Missing mode decomposition: {decomposition_path}")

        decomposition_df = pd.read_csv(decomposition_path)
        gamma_rows = decomposition_df.loc[decomposition_df[component_column].eq(0)]
        if gamma_rows.empty:
            raise ValueError(f"Missing Gamma component in {decomposition_path}")
        gamma_row = gamma_rows.iloc[0]

        gamma_k_weight_fraction = float(gamma_row["k_weight_fraction"])
        gamma_mode_px = float(gamma_row["mode_px"])
        score_px_coefficient = min(1.0, gamma_mode_px + 0.02)
        subspace_matches_mode = csv_bool(gamma_row["subspace_matches_mode"])
        raw_gamma_p_px_weight_fraction = gamma_k_weight_fraction * gamma_mode_px
        mode_score = gamma_k_weight_fraction * score_px_coefficient if subspace_matches_mode else 0.0
        rows.append({
            "mode_idx": mode_idx,
            "frequency": float(summary_row["frequency"]),
            "q": float(summary_row["q"]),
            "gamma_k_weight_fraction": gamma_k_weight_fraction,
            "gamma_subspace_p": float(gamma_row["subspace_p"]),
            "gamma_mode_px": gamma_mode_px,
            "score_px_coefficient": score_px_coefficient,
            "raw_gamma_p_px_weight_fraction": raw_gamma_p_px_weight_fraction,
            "gamma_p_px_weight_fraction": mode_score,
            "mode_score": mode_score,
            "dominant_subspace": gamma_row["dominant_subspace"],
            "dominant_mode": gamma_row["dominant_mode"],
            "subspace_matches_mode": subspace_matches_mode,
        })

    score_df = pd.DataFrame(rows)
    if gamma_subspace_p_valid_threshold is not None:
        if score_df.empty:
            score_df["is_valid"] = pd.Series(dtype=bool)
        else:
            score_df["is_valid"] = score_df["gamma_subspace_p"].gt(
                float(gamma_subspace_p_valid_threshold)
            )
        eligible_df = score_df.loc[score_df["is_valid"]]
        valid_mode_count = int(score_df["is_valid"].sum())
    else:
        eligible_df = score_df
    score_df.to_csv(out_dir / "mode_scores.csv", index=False)
    if eligible_df.empty:
        score = 0.0
        best = None
    else:
        score = float(eligible_df["mode_score"].max())
        best = eligible_df.loc[eligible_df["mode_score"].idxmax()].to_dict()

    attrs = {
        "valid_mode_count": valid_mode_count,
        "scored_mode_count": int(len(score_df)),
        "best_mode_idx": None if best is None else int(best["mode_idx"]),
        "best_frequency": None if best is None else float(best["frequency"]),
        "best_q": None if best is None else float(best["q"]),
        "gamma_k_weight_fraction": None if best is None else float(best["gamma_k_weight_fraction"]),
        "gamma_mode_px": None if best is None else float(best["gamma_mode_px"]),
        "score_px_coefficient": None if best is None else float(best["score_px_coefficient"]),
        "raw_gamma_p_px_weight_fraction": None if best is None else float(best["raw_gamma_p_px_weight_fraction"]),
        "gamma_p_px_weight_fraction": score,
        "mode_score": score,
    }
    if gamma_subspace_p_valid_threshold is not None:
        attrs.update({
            "validity_metric": "gamma_subspace_p",
            "validity_operator": ">",
            "validity_threshold": float(gamma_subspace_p_valid_threshold),
        })
    write_json(out_dir / "objective.json", {
        "score": score,
        "best_mode": best,
        "attrs": attrs,
    })
    print(f"Gamma score: {score:.6g}")
    return score, score_df


def score_strip_modes(
    out_dir: Path,
    *,
    export_dirname: str = "simulation_exports",
    gamma_subspace_p_valid_threshold: float = 0.9,
) -> tuple[float, pd.DataFrame]:
    return score_gamma_modes(
        out_dir,
        fourier_dir=strip_bulk_fourier_dir(out_dir),
        summary_filename="strip_bulk_fourier_summary.csv",
        component_column="m",
        export_dirname=export_dirname,
        respect_summary_validity=False,
        gamma_subspace_p_valid_threshold=gamma_subspace_p_valid_threshold,
    )


def score_finite_modes(out_dir: Path, *, export_dirname: str = "simulation_exports") -> tuple[float, pd.DataFrame]:
    return score_gamma_modes(
        out_dir,
        fourier_dir=finite_lattice_fourier_dir(out_dir),
        summary_filename="finite_lattice_fourier_summary.csv",
        component_column="t",
        export_dirname=export_dirname,
    )
