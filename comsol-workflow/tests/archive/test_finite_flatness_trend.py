from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from comsol_workflow.hex_lattice_utils import LatticePoint
from comsol_workflow.lattice_fourier_postprocess import (
    CENTER_NORMALIZED_INTENSITY_CMAP,
    CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL,
    CENTER_NORMALIZED_INTENSITY_TICKS,
    CENTER_NORMALIZED_INTENSITY_TICK_LABELS,
    CENTER_NORMALIZED_INTENSITY_TITLE,
    CENTER_NORMALIZED_INTENSITY_VCENTER,
    CENTER_NORMALIZED_INTENSITY_VMAX,
    CENTER_NORMALIZED_INTENSITY_VMIN,
    INTENSITY_MAPS_NORMH_FILENAME,
    center_normalized_cell_intensity,
    plot_center_normalized_cell_intensity,
    run_finite_unit_cell_intensity_maps,
)
from scripts.analysis.finite_flatness_trend import (
    CANONICAL_FULL_CELL_MAP_FILENAME,
    COMPATIBILITY_FULL_CELL_MAP_FILENAME,
    LINEAR_CENTER_MAP_CMAP,
    LINEAR_CENTER_MAP_COLORBAR_LABEL,
    LINEAR_CENTER_MAP_TICKS,
    LINEAR_CENTER_MAP_TICK_LABELS,
    LINEAR_CENTER_MAP_TITLE,
    LINEAR_CENTER_MAP_VCENTER,
    LINEAR_CENTER_MAP_VMAX,
    LINEAR_CENTER_MAP_VMIN,
    _make_shift_panel_axes,
    _shift_axis_label,
    cladding_guardrails,
    lattice_shell,
    normalize_cells_by_cavity_center,
    normalized_intensity_metrics,
    save_full_cell_maps_linear_center_normalized,
)


def test_standard_finite_intensity_map_contract(tmp_path, monkeypatch) -> None:
    cells = [
        LatticePoint(i=0, j=0, shell=0, x=0.0, y=0.0),
        LatticePoint(i=1, j=0, shell=1, x=1.0, y=0.0),
        LatticePoint(i=2, j=0, shell=2, x=2.0, y=0.0),
    ]
    fields = np.array([
        [1.0, 1.0],
        [np.sqrt(0.5), np.sqrt(0.5)],
        [2.0, 2.0],
    ])
    captured = {}

    def capture_savefig(fig, path, *args, **kwargs):
        captured.update(fig=fig, path=path, kwargs=kwargs)

    monkeypatch.setattr(Figure, "savefig", capture_savefig)
    output = tmp_path / INTENSITY_MAPS_NORMH_FILENAME

    plot_center_normalized_cell_intensity(
        output,
        cells,
        fields,
        period=1.0,
        dpi=220,
    )

    assert center_normalized_cell_intensity(cells, fields).tolist() == [
        1.0,
        0.5000000000000001,
        4.0,
    ]
    fig = captured["fig"]
    assert captured["path"] == output
    assert captured["kwargs"] == {"dpi": 220, "bbox_inches": "tight"}
    assert len(fig.axes) == 2
    ax, colorbar_ax = fig.axes
    assert ax.get_title() == CENTER_NORMALIZED_INTENSITY_TITLE
    assert ax.get_xlabel() == "x (um)"
    assert ax.get_ylabel() == "y (um)"
    assert ax.get_aspect() == 1.0
    collection = ax.collections[0]
    assert collection.cmap.name == CENTER_NORMALIZED_INTENSITY_CMAP
    assert collection.norm.vmin == CENTER_NORMALIZED_INTENSITY_VMIN
    assert collection.norm.vcenter == CENTER_NORMALIZED_INTENSITY_VCENTER
    assert collection.norm.vmax == CENTER_NORMALIZED_INTENSITY_VMAX
    assert collection.get_array().tolist() == [1.0, 0.5000000000000001, 2.0]
    assert len(collection.get_edgecolors()) == 0
    assert colorbar_ax.get_ylabel() == CENTER_NORMALIZED_INTENSITY_COLORBAR_LABEL
    assert tuple(colorbar_ax.get_yticks()) == CENTER_NORMALIZED_INTENSITY_TICKS
    assert tuple(
        label.get_text() for label in colorbar_ax.get_yticklabels()
    ) == CENTER_NORMALIZED_INTENSITY_TICK_LABELS


def test_standard_intensity_map_runner_covers_every_valid_mode(
    tmp_path,
    monkeypatch,
) -> None:
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame({
        "mode_idx": [1, 2, 3],
        "is_valid": [True, False, True],
    }).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    cells = [LatticePoint(i=0, j=0, shell=0, x=0.0, y=0.0)]
    cladding_cells = [LatticePoint(i=1, j=0, shell=1, x=1.0, y=0.0)]
    calls = []
    monkeypatch.setattr(
        "comsol_workflow.lattice_fourier_postprocess.unit_cell_rho_grid",
        lambda period, size: {"period": period, "size": size},
    )

    def capture_map(
        _export_dir,
        _output_root,
        mode_idx,
        rho_grid,
        _cells,
        _cladding_cells,
        **kwargs,
    ):
        calls.append((mode_idx, rho_grid, kwargs))
        return tmp_path / f"mode{mode_idx}.png"

    monkeypatch.setattr(
        "comsol_workflow.lattice_fourier_postprocess."
        "write_finite_unit_cell_intensity_map",
        capture_map,
    )

    outputs = run_finite_unit_cell_intensity_maps(
        tmp_path / "staging",
        export_dir,
        cells=cells,
        cladding_cells=cladding_cells,
        period=0.82,
        rho_grid_size=31,
        dpi=180,
    )

    assert outputs == [tmp_path / "mode1.png", tmp_path / "mode3.png"]
    assert [call[0] for call in calls] == [1, 3]
    assert calls[0][1] == {"period": 0.82, "size": 31}
    assert calls[0][2] == {"period": 0.82, "dpi": 180}


def test_directional_flatness_uses_x_shift_axis_label() -> None:
    directional = pd.DataFrame({"shift_kind": ["directional", "directional"]})
    scalar = pd.DataFrame({"shift_kind": ["scalar"]})

    assert _shift_axis_label(directional) == "CLADDING X SHIFT / A"
    assert _shift_axis_label(scalar) == "CLADDING_INWARD_SHIFT / A"


def test_canonical_linear_center_map_contract(tmp_path, monkeypatch) -> None:
    cells = pd.DataFrame({
        "shift_factor": [0.0, 0.0, 0.1, 0.1],
        "region": ["cavity", "cladding", "cavity", "cladding"],
        "i": [0, 1, 0, 1],
        "j": [0, 0, 0, 0],
        "x": [0.0, 1.0, 0.0, 1.0],
        "y": [0.0, 0.0, 0.0, 0.0],
        "intensity": [2.0, 5.0, 4.0, 2.0],
    })
    summary = pd.DataFrame({
        "shift_factor": [0.0, 0.1],
        "period": [1.0, 1.0],
    })
    captured = {}

    def capture_savefig(fig, path, *args, **kwargs):
        captured.update(fig=fig, path=path, kwargs=kwargs)

    monkeypatch.setattr(Figure, "savefig", capture_savefig)
    output = tmp_path / CANONICAL_FULL_CELL_MAP_FILENAME

    save_full_cell_maps_linear_center_normalized(cells, summary, output)

    fig = captured["fig"]
    assert captured["path"] == output
    assert captured["kwargs"]["dpi"] == 220
    assert fig._suptitle.get_text() == LINEAR_CENTER_MAP_TITLE
    assert len(fig.axes) == 3
    for ax in fig.axes[:2]:
        collection = ax.collections[0]
        assert collection.cmap.name == LINEAR_CENTER_MAP_CMAP
        assert collection.norm.vmin == LINEAR_CENTER_MAP_VMIN
        assert collection.norm.vcenter == LINEAR_CENTER_MAP_VCENTER
        assert collection.norm.vmax == LINEAR_CENTER_MAP_VMAX
        assert ax.get_aspect() == 1.0
    colorbar_ax = fig.axes[-1]
    assert colorbar_ax.get_ylabel() == LINEAR_CENTER_MAP_COLORBAR_LABEL
    assert tuple(colorbar_ax.get_yticks()) == LINEAR_CENTER_MAP_TICKS
    assert tuple(label.get_text() for label in colorbar_ax.get_yticklabels()) == (
        LINEAR_CENTER_MAP_TICK_LABELS
    )
    assert CANONICAL_FULL_CELL_MAP_FILENAME == "cavity_cladding_intensity_maps.png"
    assert COMPATIBILITY_FULL_CELL_MAP_FILENAME == "intensity_maps_normH.png"


def test_shift_panel_layout_supports_non_ten_point_series() -> None:
    fig, axes = _make_shift_panel_axes(7)
    try:
        assert len(axes) == 7
        assert len(fig.axes) == 8
        assert sum(ax.get_visible() for ax in fig.axes) == 7
    finally:
        fig.clear()


def test_lattice_shell_uses_hexagonal_distance() -> None:
    indices = np.array([[0, 0], [2, -1], [-2, 1], [1, 1], [-1, -1]])
    assert lattice_shell(indices).tolist() == [0, 2, 2, 2, 2]


def test_uniform_intensity_is_perfectly_flat() -> None:
    metrics = normalized_intensity_metrics(np.ones(7))
    assert metrics["cv"] == 0.0
    assert metrics["flatness"] == 1.0
    assert metrics["dmax"] == 0.0


def test_flatness_score_matches_cv_identity() -> None:
    metrics = normalized_intensity_metrics(np.array([0.5, 1.5]))
    assert np.isclose(metrics["cv"], 0.5)
    assert np.isclose(metrics["flatness"], 0.8)
    assert np.isclose(metrics["dmax"], 0.5)


def test_cladding_guardrails_detect_outward_increase() -> None:
    intensity = np.array([0.8] * 6 + [0.4] * 12 + [0.6] * 18)
    shells = np.array([2] * 6 + [3] * 12 + [4] * 18)
    metrics = cladding_guardrails(
        cladding_intensity=intensity,
        cladding_shells=shells,
        cavity_mean=1.0,
        cavity_outer_mean_ratio=0.9,
        bulk_radius=1,
    )
    assert np.isclose(metrics["cladding_upward_variation"], 0.2)
    assert np.isclose(metrics["cladding_hotspot_max"], 0.8)
    assert np.isclose(metrics["first_cladding_over_cavity_outer"], 8.0 / 9.0)


def test_center_cell_normalization_is_per_shift() -> None:
    cells = pd.DataFrame(
        {
            "shift_factor": [0.0, 0.0, 0.1, 0.1],
            "region": ["cavity", "cladding", "cavity", "cladding"],
            "i": [0, 1, 0, 1],
            "j": [0, 0, 0, 0],
            "intensity": [2.0, 5.0, 4.0, 2.0],
        }
    )
    normalized = normalize_cells_by_cavity_center(cells)
    assert normalized["intensity_over_center"].tolist() == [1.0, 2.5, 1.0, 0.5]
