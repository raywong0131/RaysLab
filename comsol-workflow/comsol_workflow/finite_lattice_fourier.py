from __future__ import annotations

import numpy as np

from .hex_lattice_utils import unit_cell_corners
from .polygon_utils import point_inside_or_on_polygon, polygon_area


def direct_lattice_basis(period: float) -> np.ndarray:
    return np.array([
        [period, 0.5 * period],
        [0.0, np.sqrt(3.0) * period / 2.0],
    ], dtype=float)


def reciprocal_lattice_basis(period: float) -> np.ndarray:
    return 2.0 * np.pi * np.linalg.inv(direct_lattice_basis(period)).T


def hex_quotient_cell_count(shells: int) -> int:
    return 3 * shells * (shells + 1) + 1


def hex_cyclic_indices(indices: np.ndarray, shells: int) -> np.ndarray:
    n_cells = hex_quotient_cell_count(shells)
    indices = np.asarray(indices, dtype=int)
    return ((2 * shells + 1) * indices[:, 0] + shells * indices[:, 1]) % n_cells


def hex_xi_values(shells: int) -> np.ndarray:
    n_cells = hex_quotient_cell_count(shells)
    direction = np.array([2 * shells + 1, shells], dtype=float)
    return np.mod(np.arange(n_cells, dtype=float)[:, None] * direction[None, :] / n_cells, 1.0)


def hex_spiral_indices(points: np.ndarray) -> np.ndarray:
    """Return zero-based concentric-hexagon labels aligned with ``points``.

    The center is 0.  Each shell starts at the upper point in the first
    quadrant and proceeds counterclockwise.  The input order is irrelevant;
    the returned array contains the display label for each input row.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) == 0:
        raise ValueError("points must be a non-empty array with shape (n, 2)")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must contain only finite coordinates")

    n_points = len(points)
    shells = int(round((np.sqrt(12.0 * n_points - 3.0) - 3.0) / 6.0))
    if hex_quotient_cell_count(shells) != n_points:
        raise ValueError(
            f"Expected a complete hexagonal point set, got {n_points} points"
        )
    if n_points == 1:
        return np.zeros(1, dtype=int)

    radii_from_origin = np.linalg.norm(points, axis=1)
    center_idx = int(np.argmin(radii_from_origin))
    relative = points - points[center_idx]
    radii = np.linalg.norm(relative, axis=1)
    nonzero_radii = radii[radii > np.finfo(float).eps]
    if len(nonzero_radii) == 0:
        raise ValueError("points do not contain a nonzero hexagonal shell")
    spacing = float(np.min(nonzero_radii))
    tolerance = max(spacing * 1e-7, np.finfo(float).eps * 100.0)
    if radii_from_origin[center_idx] > tolerance:
        raise ValueError("hexagonal point set must contain the origin")

    nearest = np.flatnonzero(np.isclose(radii, spacing, rtol=1e-7, atol=tolerance))
    if len(nearest) != 6:
        raise ValueError(f"Expected six nearest neighbours, got {len(nearest)}")

    angles = np.mod(np.arctan2(relative[:, 1], relative[:, 0]), 2.0 * np.pi)
    first_quadrant = nearest[
        (relative[nearest, 0] > tolerance) & (relative[nearest, 1] > tolerance)
    ]
    if len(first_quadrant) != 1:
        raise ValueError(
            "Expected one nearest neighbour above the origin in the first quadrant"
        )
    start_idx = int(first_quadrant[0])
    start_angle = float(angles[start_idx])

    nearest_delta = np.mod(angles[nearest] - start_angle, 2.0 * np.pi)
    nearest_order = nearest[np.argsort(nearest_delta)]
    next_idx = int(nearest_order[1])
    basis = np.column_stack([relative[start_idx], relative[next_idx]])
    if abs(float(np.linalg.det(basis))) <= tolerance**2:
        raise ValueError("Failed to recover a non-degenerate triangular-lattice basis")

    axial_float = np.linalg.solve(basis, relative.T).T
    axial = np.rint(axial_float).astype(int)
    reconstructed = axial @ basis.T
    reconstruction_error = float(
        np.max(np.linalg.norm(relative - reconstructed, axis=1))
    )
    if reconstruction_error > tolerance:
        raise ValueError(
            "points are not a complete regular triangular-lattice hexagon "
            f"(reconstruction error {reconstruction_error:g})"
        )
    if len(np.unique(axial, axis=0)) != n_points:
        raise ValueError("points map to duplicate triangular-lattice coordinates")

    shell_values = np.maximum.reduce(
        [np.abs(axial[:, 0]), np.abs(axial[:, 1]), np.abs(axial[:, 0] + axial[:, 1])]
    )
    labels = np.empty(n_points, dtype=int)
    next_label = 0
    angular_offset = np.mod(angles - start_angle, 2.0 * np.pi)
    for shell in range(shells + 1):
        members = np.flatnonzero(shell_values == shell)
        expected_count = 1 if shell == 0 else 6 * shell
        if len(members) != expected_count:
            raise ValueError(
                f"Expected {expected_count} points on shell {shell}, got {len(members)}"
            )
        ordered = members[np.argsort(angular_offset[members])]
        labels[ordered] = np.arange(next_label, next_label + expected_count, dtype=int)
        next_label += expected_count

    if not np.array_equal(np.sort(labels), np.arange(n_points, dtype=int)):
        raise RuntimeError("hexagonal spiral labels are not a complete permutation")
    return labels


def clip_polygon_halfspace(polygon: np.ndarray, normal: np.ndarray, offset: float) -> np.ndarray:
    clipped = []
    for idx, current in enumerate(polygon):
        previous = polygon[idx - 1]
        current_inside = float(np.dot(current, normal)) <= offset + 1e-12
        previous_inside = float(np.dot(previous, normal)) <= offset + 1e-12
        edge = current - previous
        denom = float(np.dot(edge, normal))
        if current_inside != previous_inside and abs(denom) > 1e-14:
            t = (offset - float(np.dot(previous, normal))) / denom
            clipped.append(previous + t * edge)
        if current_inside:
            clipped.append(current)
    return np.asarray(clipped, dtype=float)


def first_bz_polygon(period: float) -> np.ndarray:
    basis = reciprocal_lattice_basis(period)
    vectors = []
    for h0 in range(-3, 4):
        for h1 in range(-3, 4):
            if h0 == 0 and h1 == 0:
                continue
            vectors.append(basis @ np.array([h0, h1], dtype=float))
    vectors.sort(key=np.linalg.norm)
    bound = 2.0 * max(float(np.linalg.norm(vector)) for vector in vectors)
    polygon = np.array([
        [-bound, -bound],
        [bound, -bound],
        [bound, bound],
        [-bound, bound],
    ], dtype=float)
    for vector in vectors:
        polygon = clip_polygon_halfspace(polygon, vector, 0.5 * float(np.dot(vector, vector)))
        if len(polygon) == 0:
            raise RuntimeError("Failed to construct the first Brillouin zone polygon.")
    center = polygon.mean(axis=0)
    angles = np.arctan2(polygon[:, 1] - center[1], polygon[:, 0] - center[0])
    return polygon[np.argsort(angles)]


def fold_xis_to_first_bz(xis: np.ndarray, period: float) -> np.ndarray:
    basis = reciprocal_lattice_basis(period)
    folded = []
    for xi in np.asarray(xis, dtype=float):
        candidates = []
        for h0 in range(-3, 4):
            for h1 in range(-3, 4):
                shifted = xi + np.array([h0, h1], dtype=float)
                physical = basis @ shifted
                candidates.append((float(np.dot(physical, physical)), physical))
        folded.append(min(candidates, key=lambda item: item[0])[1])
    return np.asarray(folded, dtype=float)


def unit_cell_rho_grid(period: float, grid_size: int) -> dict[str, np.ndarray]:
    cell = unit_cell_corners(period)
    x_axis = np.linspace(float(cell[:, 0].min()), float(cell[:, 0].max()), grid_size)
    y_axis = np.linspace(float(cell[:, 1].min()), float(cell[:, 1].max()), grid_size)
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")
    points = np.column_stack([xx.ravel(), yy.ravel()])
    keep = np.array([point_inside_or_on_polygon(point, cell) for point in points], dtype=bool)
    return {
        "x_axis": x_axis,
        "y_axis": y_axis,
        "xx": xx,
        "yy": yy,
        "mask": keep.reshape(xx.shape),
        "points": points[keep],
    }


def finite_lattice_fourier_components(
    cell_fields: np.ndarray,
    rho_points: np.ndarray,
    xi_values: np.ndarray,
    period: float,
) -> dict[str, np.ndarray | float]:
    cell_fields = np.asarray(cell_fields, dtype=complex)
    rho_points = np.asarray(rho_points, dtype=float)
    xi_values = np.asarray(xi_values, dtype=float)

    f_components = np.fft.ifft(cell_fields, axis=0, norm="ortho")
    reconstructed = np.fft.fft(f_components, axis=0, norm="ortho")
    reconstruction_error = float(np.linalg.norm(reconstructed - cell_fields) / np.linalg.norm(cell_fields))

    reciprocal = reciprocal_lattice_basis(period)
    raw_k = xi_values @ reciprocal.T
    phase = np.exp(1j * (raw_k @ rho_points.T))
    g_components = f_components * phase

    sample_weight = polygon_area(unit_cell_corners(period)) / len(rho_points)
    original_norm = float(sample_weight * np.sum(np.abs(cell_fields) ** 2))
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


def arbitrary_finite_lattice_fourier_components(
    cell_fields: np.ndarray,
    cell_centers: np.ndarray,
    rho_points: np.ndarray,
    k_values: np.ndarray,
    period: float,
) -> dict[str, np.ndarray | float]:
    """Evaluate a direct finite-lattice DFT without a cyclic quotient."""
    fields = np.asarray(cell_fields, dtype=complex)
    centers = np.asarray(cell_centers, dtype=float)
    rho = np.asarray(rho_points, dtype=float)
    wavevectors = np.asarray(k_values, dtype=float)
    if fields.ndim != 2:
        raise ValueError("cell_fields must have shape (n_cells, n_rho)")
    if centers.shape != (fields.shape[0], 2):
        raise ValueError("cell_centers must have shape (n_cells, 2)")
    if rho.shape != (fields.shape[1], 2):
        raise ValueError("rho_points must have shape (n_rho, 2)")
    if wavevectors.ndim != 2 or wavevectors.shape[1] != 2 or len(wavevectors) == 0:
        raise ValueError("k_values must be a non-empty array with shape (n_k, 2)")

    phase_centers = np.exp(-1j * (wavevectors @ centers.T))
    f_components = phase_centers @ fields / np.sqrt(float(len(centers)))
    phase_rho = np.exp(1j * (wavevectors @ rho.T))
    g_components = f_components * phase_rho

    sample_weight = polygon_area(unit_cell_corners(period)) / len(rho)
    weights = sample_weight * np.sum(np.abs(g_components) ** 2, axis=1)
    total_weight = float(np.sum(weights))
    weight_fraction = (
        weights / total_weight if total_weight > 0.0 else np.zeros_like(weights)
    )
    p_components = np.zeros_like(g_components)
    positive_weight = weights > 0.0
    p_components[positive_weight] = (
        g_components[positive_weight] / np.sqrt(weights[positive_weight])[:, None]
    )
    profile_norms = sample_weight * np.sum(np.abs(p_components) ** 2, axis=1)
    profile_norm_error = (
        float(np.max(np.abs(profile_norms[positive_weight] - 1.0)))
        if np.any(positive_weight)
        else np.nan
    )
    return {
        "F": f_components,
        "G": g_components,
        "P": p_components,
        "weights": weights,
        "weight_fraction": weight_fraction,
        "profile_norms": profile_norms,
        "reconstruction_error": np.nan,
        "parseval_error": np.nan,
        "profile_norm_error": profile_norm_error,
    }
