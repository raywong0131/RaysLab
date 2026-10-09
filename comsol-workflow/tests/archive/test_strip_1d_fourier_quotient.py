#!/usr/bin/env python3
"""Verify the exact 1D finite-Fourier quotient for the strip geometry."""

from __future__ import annotations

from tests._paths import output_root
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow-tests")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]

from comsol_workflow.hex_lattice_utils import (
    bulk_points as _bulk_points,
    cell_polygon,
    lattice_points,
)


OUTPUT_DIR = output_root() / "strip_1d_fourier_quotient"
DPI = 220
TOL = 1e-10

A = 0.82
BULK_RADIUS = 3
CLADDING_LAYERS = 4
STRIP_X_MIN = -A / 2.0
STRIP_X_MAX = A / 2.0


def polygon_area(vertices: np.ndarray) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def dedupe_vertices(vertices: np.ndarray, tolerance: float = 1e-12) -> np.ndarray:
    if len(vertices) == 0:
        return vertices
    kept = [vertices[0]]
    for vertex in vertices[1:]:
        if np.linalg.norm(vertex - kept[-1]) > tolerance:
            kept.append(vertex)
    if len(kept) > 1 and np.linalg.norm(kept[0] - kept[-1]) <= tolerance:
        kept.pop()
    return np.asarray(kept, dtype=float)


def clip_polygon_x_leq(polygon: np.ndarray, x_max: float) -> np.ndarray:
    clipped = []
    for idx, current in enumerate(polygon):
        previous = polygon[idx - 1]
        current_inside = current[0] <= x_max + 1e-12
        previous_inside = previous[0] <= x_max + 1e-12
        dx = current[0] - previous[0]
        if current_inside != previous_inside and abs(dx) > 1e-14:
            fraction = (x_max - previous[0]) / dx
            clipped.append(previous + fraction * (current - previous))
        if current_inside:
            clipped.append(current)
    return dedupe_vertices(np.asarray(clipped, dtype=float))


def clip_polygon_x_geq(polygon: np.ndarray, x_min: float) -> np.ndarray:
    clipped = []
    for idx, current in enumerate(polygon):
        previous = polygon[idx - 1]
        current_inside = current[0] >= x_min - 1e-12
        previous_inside = previous[0] >= x_min - 1e-12
        dx = current[0] - previous[0]
        if current_inside != previous_inside and abs(dx) > 1e-14:
            fraction = (x_min - previous[0]) / dx
            clipped.append(previous + fraction * (current - previous))
        if current_inside:
            clipped.append(current)
    return dedupe_vertices(np.asarray(clipped, dtype=float))


def clip_polygon_to_strip(polygon: np.ndarray) -> np.ndarray | None:
    clipped = clip_polygon_x_geq(np.asarray(polygon, dtype=float), STRIP_X_MIN)
    if len(clipped) == 0:
        return None
    clipped = clip_polygon_x_leq(clipped, STRIP_X_MAX)
    if len(clipped) < 3:
        return None
    if abs(polygon_area(clipped)) <= 1e-14:
        return None
    return clipped


def strip_cell_records() -> list[dict[str, object]]:
    bulk = _bulk_points(BULK_RADIUS, A)
    cladding = lattice_points(
        BULK_RADIUS + CLADDING_LAYERS,
        A,
        min_shell=BULK_RADIUS + 1,
    )

    records = []
    for region, points in (("bulk", bulk), ("cladding", cladding)):
        for point in points:
            clipped = clip_polygon_to_strip(cell_polygon(point.x, point.y, A))
            if clipped is None:
                continue
            records.append({
                "region": region,
                "point": point,
                "polygon": clipped,
            })
    return records


def row_groups(records: list[dict[str, object]], region: str | None = None) -> dict[int, list[dict[str, object]]]:
    groups: dict[int, list[dict[str, object]]] = {}
    for record in records:
        if region is not None and record["region"] != region:
            continue
        point = record["point"]
        groups.setdefault(point.j, []).append(record)
    return groups


def one_dimensional_unitary(rows: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows, dtype=int)
    n_rows = len(rows)
    modes = np.arange(n_rows, dtype=int)
    return np.exp(-2j * np.pi * rows[:, None] * modes[None, :] / n_rows) / np.sqrt(n_rows)


def assert_unitary(U: np.ndarray, label: str, tolerance: float = TOL) -> float:
    gram_error = float(np.linalg.norm(U.conj().T @ U - np.eye(U.shape[1])))
    assert gram_error < tolerance, f"{label}: U^dagger U error {gram_error:.3e}"
    row_error = float(np.linalg.norm(U @ U.conj().T - np.eye(U.shape[0])))
    assert row_error < tolerance, f"{label}: U U^dagger error {row_error:.3e}"
    return gram_error


def assert_one_dimensional_quotient(rows: list[int], label: str) -> tuple[float, float, float]:
    rows_array = np.asarray(rows, dtype=int)
    n_rows = len(rows_array)
    residues = np.mod(rows_array, n_rows)
    assert set(residues.tolist()) == set(range(n_rows)), f"{label}: rows are not Z_{n_rows} representatives"

    U = one_dimensional_unitary(rows_array)
    gram_error = assert_unitary(U, label)

    rng = np.random.default_rng(20260702 + n_rows)
    field = rng.normal(size=(n_rows, 9)) + 1j * rng.normal(size=(n_rows, 9))
    coeffs_dense = U.conj().T @ field
    reconstructed = U @ coeffs_dense
    recon_error = float(np.linalg.norm(reconstructed - field) / np.linalg.norm(field))
    assert recon_error < TOL, f"{label}: dense reconstruction error {recon_error:.3e}"

    parseval_error = abs(float(np.vdot(field, field).real - np.vdot(coeffs_dense, coeffs_dense).real))
    parseval_error /= float(np.vdot(field, field).real)
    assert parseval_error < TOL, f"{label}: Parseval error {parseval_error:.3e}"

    field_by_residue = np.empty_like(field)
    field_by_residue[residues] = field
    coeffs_fft = np.fft.ifft(field_by_residue, axis=0, norm="ortho")
    fft_error = float(np.linalg.norm(coeffs_fft - coeffs_dense) / np.linalg.norm(coeffs_dense))
    assert fft_error < TOL, f"{label}: dense/FFT coefficient error {fft_error:.3e}"

    return gram_error, recon_error, fft_error


def assert_bloch_phase_gluing(groups: dict[int, list[dict[str, object]]]) -> None:
    rng = np.random.default_rng(410)
    profiles = {
        row: rng.normal(size=5) + 1j * rng.normal(size=5)
        for row in groups
    }

    for kx_fraction in (0.0, 0.17, 0.5):
        for row, records in groups.items():
            canonical_i = min(record["point"].i for record in records)
            for record in records:
                point = record["point"]
                periodic_shift = point.i - canonical_i
                physical = np.exp(-2j * np.pi * kx_fraction * periodic_shift) * profiles[row]
                periodic_gauge = np.exp(+2j * np.pi * kx_fraction * periodic_shift) * physical
                error = float(np.linalg.norm(periodic_gauge - profiles[row]))
                assert error < TOL, (
                    f"kx={kx_fraction:g}, row={row}, i={point.i}: "
                    f"Bloch-phase gluing error {error:.3e}"
                )


def plot_current_strip_rows(records: list[dict[str, object]]) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "current_strip_row_quotient.png"

    rows = sorted({record["point"].j for record in records})
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
            linewidth=0.42,
            alpha=0.72 if record["region"] == "bulk" else 0.38,
        ))
        center = np.mean(record["polygon"], axis=0)
        ax.text(center[0], center[1], str(point.j), ha="center", va="center", fontsize=6.5)

    ax.axvline(STRIP_X_MIN, color="#dc2626", linestyle="--", linewidth=1.0)
    ax.axvline(STRIP_X_MAX, color="#dc2626", linestyle="--", linewidth=1.0)
    all_points = np.vstack([record["polygon"] for record in records])
    pad = 0.25 * A
    ax.set_xlim(float(all_points[:, 0].min() - pad), float(all_points[:, 0].max() + pad))
    ax.set_ylim(float(all_points[:, 1].min() - pad), float(all_points[:, 1].max() + pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("simulation x (um)")
    ax.set_ylabel("simulation y (um)")
    ax.set_title("vertical strip_1d quotient rows")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="row j")
    fig.savefig(output_path, dpi=DPI)
    plt.close(fig)
    return output_path


def test_current_strip_cell_structure() -> list[dict[str, object]]:
    records = strip_cell_records()
    bulk_groups = row_groups(records, region="bulk")
    all_groups = row_groups(records)

    assert len(records) == 23, f"expected 23 clipped cell pieces, got {len(records)}"
    assert sum(record["region"] == "bulk" for record in records) == 11
    assert sum(record["region"] == "cladding" for record in records) == 12
    assert sorted(bulk_groups) == list(range(-BULK_RADIUS, BULK_RADIUS + 1))
    assert sorted(all_groups) == list(range(-(BULK_RADIUS + CLADDING_LAYERS), BULK_RADIUS + CLADDING_LAYERS + 1))
    assert any(len(group) > 1 for group in bulk_groups.values()), "expected some bulk rows to be split at the seam"
    assert len(bulk_groups) < sum(record["region"] == "bulk" for record in records)

    for groups in (bulk_groups, all_groups):
        for row, group in groups.items():
            indices = sorted(record["point"].i for record in group)
            assert indices == list(range(indices[0], indices[-1] + 1)), f"row {row}: noncontiguous seam pieces"

    print(
        "[Pass] current strip cells: "
        f"pieces={len(records)}, bulk_rows={sorted(bulk_groups)}, all_rows={sorted(all_groups)}"
    )
    return records


def test_current_bulk_strip_quotient() -> tuple[float, float, float]:
    records = strip_cell_records()
    rows = sorted(row_groups(records, region="bulk"))
    result = assert_one_dimensional_quotient(rows, "current strip bulk rows")
    print(
        "[Pass] current strip bulk Z_7 quotient: "
        f"unitary={result[0]:.2e}, recon={result[1]:.2e}, dense_fft={result[2]:.2e}"
    )
    return result


def test_current_all_strip_quotient() -> tuple[float, float, float]:
    records = strip_cell_records()
    rows = sorted(row_groups(records))
    result = assert_one_dimensional_quotient(rows, "current strip all rows")
    print(
        "[Pass] current strip all-row Z_15 quotient: "
        f"unitary={result[0]:.2e}, recon={result[1]:.2e}, dense_fft={result[2]:.2e}"
    )
    return result


def test_bloch_phase_gluing_for_current_strip() -> None:
    records = strip_cell_records()
    assert_bloch_phase_gluing(row_groups(records, region="bulk"))
    assert_bloch_phase_gluing(row_groups(records))
    print("[Pass] current strip seam pieces glue with arbitrary Bloch phase")


def main() -> None:
    records = test_current_strip_cell_structure()
    test_current_bulk_strip_quotient()
    test_current_all_strip_quotient()
    test_bloch_phase_gluing_for_current_strip()
    figure_path = plot_current_strip_rows(records)
    print(f"\nGenerated figure:\n  {figure_path}")
    print("\nAll strip 1D finite-Fourier quotient checks passed.")


if __name__ == "__main__":
    main()
