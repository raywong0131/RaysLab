#!/usr/bin/env python3
"""Focused checks for shared finite hex-lattice, polygon, and Fourier helpers."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]

from comsol_workflow.finite_lattice_fourier import (
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
from comsol_workflow.hex_lattice_utils import boundary_from_lattice_cells, bulk_points, unit_cell_corners
from comsol_workflow.lattice_fourier_postprocess import (
    normalize_profile_magnitudes,
    plot_finite_mode_k_weights,
    plot_top_profiles,
    spiral_index_annotation_mask,
)
from comsol_workflow.polygon_utils import (
    clip_polygon_to_convex_region,
    clip_polygon_to_outside_region,
    offset_polygon,
    polygon_area,
    validate_clipped_hole,
)


def test_hex_lattice_boundary_and_quotient() -> None:
    period = 0.82
    shells = 3
    cells = bulk_points(shells, period)
    assert len(cells) == hex_quotient_cell_count(shells)

    boundary = boundary_from_lattice_cells(cells, period)
    assert len(boundary) > 6
    assert polygon_area(boundary) > 0.0

    indices = np.array([[point.i, point.j] for point in cells], dtype=int)
    cyclic = hex_cyclic_indices(indices, shells)
    assert set(cyclic.tolist()) == set(range(hex_quotient_cell_count(shells)))


def test_hex_spiral_indices_start_at_center_and_proceed_counterclockwise() -> None:
    cells = bulk_points(2, 1.0)
    coordinates = np.array([[point.x, point.y] for point in cells], dtype=float)
    labels = hex_spiral_indices(coordinates)
    labels_by_axial_index = {
        (point.i, point.j): int(label)
        for point, label in zip(cells, labels, strict=True)
    }

    expected_order = [
        (0, 0),
        (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1), (1, 0),
        (0, 2), (-1, 2), (-2, 2), (-2, 1), (-2, 0), (-1, -1),
        (0, -2), (1, -2), (2, -2), (2, -1), (2, 0), (1, 1),
    ]
    assert [labels_by_axial_index[index] for index in expected_order] == list(range(19))


def test_hex_spiral_indices_are_input_order_independent() -> None:
    cells = bulk_points(3, 0.82)
    coordinates = np.array([[point.x, point.y] for point in cells], dtype=float)
    labels = hex_spiral_indices(coordinates)
    permutation = np.random.default_rng(42).permutation(len(coordinates))
    shuffled_labels = hex_spiral_indices(coordinates[permutation])

    restored = np.empty_like(shuffled_labels)
    restored[permutation] = shuffled_labels
    assert np.array_equal(restored, labels)


def test_finite_k_grid_uses_the_same_concentric_counterclockwise_indices() -> None:
    period = 0.82
    shells = 2
    folded = fold_xis_to_first_bz(hex_xi_values(shells), period)
    folded_before = folded.copy()
    labels = hex_spiral_indices(folded)
    ordered = folded[np.argsort(labels)]

    assert np.array_equal(np.sort(labels), np.arange(hex_quotient_cell_count(shells)))
    assert np.allclose(ordered[0], 0.0)
    for start, stop in ((1, 7), (7, 19)):
        shell = ordered[start:stop]
        assert shell[0, 0] > 0.0
        assert shell[0, 1] > 0.0
        angles = np.mod(
            np.arctan2(shell[:, 1], shell[:, 0])
            - np.arctan2(shell[0, 1], shell[0, 0]),
            2.0 * np.pi,
        )
        assert np.all(np.diff(angles) > 0.0)
    assert np.array_equal(folded, folded_before)


def test_spiral_index_annotations_are_limited_to_first_four_shells() -> None:
    labels = np.arange(61, dtype=int)
    mask = spiral_index_annotation_mask(labels)

    assert np.array_equal(labels[mask], np.arange(37, dtype=int))
    assert not np.any(mask[37:])


def test_polygon_offset_and_outside_clipping() -> None:
    boundary = np.array([
        [-1.0, -1.0],
        [1.0, -1.0],
        [1.0, 1.0],
        [-1.0, 1.0],
    ], dtype=float)
    expanded = offset_polygon(boundary, 0.2)
    assert polygon_area(expanded) > polygon_area(boundary)

    triangle = np.array([
        [0.5, 0.0],
        [1.5, 0.25],
        [1.5, -0.25],
    ], dtype=float)
    pieces = clip_polygon_to_outside_region(triangle, boundary)
    validate_clipped_hole(triangle, pieces, boundary)
    assert len(pieces) == 1
    assert np.isclose(pieces[0][:, 0].max(), 1.5)
    assert np.isclose(pieces[0][:, 0].min(), 1.0)


def test_polygon_clipping_to_first_quadrant_convex_region() -> None:
    subject = np.array(
        [
            [-0.5, -0.5],
            [1.5, -0.5],
            [1.5, 1.5],
            [-0.5, 1.5],
        ],
        dtype=float,
    )
    first_quadrant_window = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        dtype=float,
    )

    clipped = clip_polygon_to_convex_region(subject, first_quadrant_window)

    assert polygon_area(clipped) == 1.0
    assert np.min(clipped[:, 0]) == 0.0
    assert np.min(clipped[:, 1]) == 0.0
    assert np.max(clipped[:, 0]) == 1.0
    assert np.max(clipped[:, 1]) == 1.0


def test_finite_lattice_fourier_components() -> None:
    period = 0.82
    shells = 2
    n_cells = hex_quotient_cell_count(shells)
    rho_grid = unit_cell_rho_grid(period, 9)
    rho_points = rho_grid["points"]
    xis = hex_xi_values(shells)

    cell_idx = np.arange(n_cells, dtype=float)[:, None]
    rho_idx = np.arange(len(rho_points), dtype=float)[None, :]
    field = np.exp(0.17j * cell_idx) * (1.0 + 0.03 * rho_idx)

    result = finite_lattice_fourier_components(field, rho_points, xis, period)
    assert result["reconstruction_error"] < 1e-12
    assert result["parseval_error"] < 1e-12
    assert result["profile_norm_error"] < 1e-12
    assert np.isclose(np.sum(result["weight_fraction"]), 1.0)

    folded = fold_xis_to_first_bz(xis, period)
    bz = first_bz_polygon(period)
    assert folded.shape == (n_cells, 2)
    assert bz.shape[1] == 2
    assert len(bz) >= 4
    assert polygon_area(unit_cell_corners(period)) > 0.0


def test_arbitrary_finite_lattice_dft_accepts_non_hex_cell_set() -> None:
    period = 0.82
    rho_points = unit_cell_rho_grid(period, 7)["points"]
    centers = np.array(
        [[-period, 0.0], [0.0, 0.0], [period, 0.0], [period / 2.0, 0.7]],
        dtype=float,
    )
    fields = np.ones((len(centers), len(rho_points)), dtype=complex)
    k_values = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    result = arbitrary_finite_lattice_fourier_components(
        fields,
        centers,
        rho_points,
        k_values,
        period,
    )
    assert result["P"].shape == (len(k_values), len(rho_points))
    assert np.isclose(np.sum(result["weight_fraction"]), 1.0)
    assert np.isnan(result["reconstruction_error"])


def test_finite_k_weight_plot_uses_approved_standard(tmp_path: Path) -> None:
    period = 0.82
    shells = 1
    rho_grid = unit_cell_rho_grid(period, 9)
    rho_points = rho_grid["points"]
    folded = fold_xis_to_first_bz(hex_xi_values(shells), period)
    bz = first_bz_polygon(period)
    n_modes = len(folded)

    fractions = np.full(n_modes, 0.15 / (n_modes - 1), dtype=float)
    fractions[0] = 0.85
    base_profile = 1.0 + 0.4 * np.cos(2.0 * np.pi * rho_points[:, 0] / period)
    profiles = np.array([
        base_profile * np.exp(1j * mode_idx / n_modes)
        for mode_idx in range(n_modes)
    ])
    output_path = tmp_path / "k_weight_Hz.png"

    plot_finite_mode_k_weights(
        output_path,
        {"re": 198.106790548, "q": 5070.723304},
        folded,
        bz,
        fractions,
        rho_points,
        profiles,
        period=period,
        dpi=60,
    )

    assert output_path.is_file()
    assert output_path.stat().st_size > 0


def test_top_profile_magnitudes_use_one_shared_normalization() -> None:
    profiles = np.array(
        [
            [1.0 + 0.0j, 2.0j],
            [0.5 + 0.0j, 1.0j],
        ]
    )

    normalized = normalize_profile_magnitudes(profiles)

    assert np.isclose(np.max(normalized), 1.0)
    assert np.allclose(normalized, [[0.5, 1.0], [0.25, 0.5]])
    assert np.allclose(normalize_profile_magnitudes(np.zeros((2, 3))), 0.0)


def test_top_profile_plots_use_approved_titles_labels_and_palettes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    import matplotlib.pyplot as plt

    period = 0.82
    rho_grid = unit_cell_rho_grid(period, 9)
    rho_points = rho_grid["points"]
    profile_count = 6
    base = 1.0 + 0.3 * np.cos(2.0 * np.pi * rho_points[:, 0] / period)
    profiles = np.array([
        base * np.exp(1j * (index + 1) * rho_points[:, 1] / period)
        for index in range(profile_count)
    ])
    fractions = np.linspace(0.6, 0.1, profile_count)
    order = np.arange(profile_count)
    display_indices = np.array([0, 6, 3, 1, 4, 2])
    captured = []
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda figure: captured.append(figure))

    plot_top_profiles(
        tmp_path / "norm.png",
        1,
        rho_grid,
        profiles,
        fractions,
        order,
        period=period,
        top_k=profile_count,
        phase=False,
        title_prefix="finite-lattice",
        dpi=30,
        display_indices=display_indices,
    )
    norm_figure = captured.pop()
    assert norm_figure._suptitle.get_text() == (
        "Top weight finite-cavity profiles: Intensity"
    )
    assert norm_figure.axes[0].get_title() == "k-index=0, P(k)=0.6"
    assert norm_figure.axes[0].collections[0].cmap.name == "plasma"

    plot_top_profiles(
        tmp_path / "phase.png",
        1,
        rho_grid,
        profiles,
        fractions,
        order,
        period=period,
        top_k=profile_count,
        phase=True,
        title_prefix="finite-lattice",
        dpi=30,
        display_indices=display_indices,
    )
    phase_figure = captured.pop()
    assert phase_figure._suptitle.get_text() == (
        "Top weight finite-cavity profiles: Phase"
    )
    assert phase_figure.axes[0].get_title() == "k-index=0, P(k)=0.6"
    assert phase_figure.axes[0].collections[0].cmap.name == "twilight_shifted"

    original_close(norm_figure)
    original_close(phase_figure)


def main() -> None:
    test_hex_lattice_boundary_and_quotient()
    test_hex_spiral_indices_start_at_center_and_proceed_counterclockwise()
    test_hex_spiral_indices_are_input_order_independent()
    test_finite_k_grid_uses_the_same_concentric_counterclockwise_indices()
    test_spiral_index_annotations_are_limited_to_first_four_shells()
    test_polygon_offset_and_outside_clipping()
    test_polygon_clipping_to_first_quadrant_convex_region()
    test_finite_lattice_fourier_components()
    test_top_profile_magnitudes_use_one_shared_normalization()
    print("shared finite lattice/polygon/Fourier utility checks passed")


if __name__ == "__main__":
    main()
