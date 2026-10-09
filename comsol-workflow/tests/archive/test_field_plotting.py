from pathlib import Path
import re

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pytest

from comsol_workflow.field_plotting import (
    HEAT_CAMERA_COLORMAP,
    WAVE_COLORMAP,
    _create_field_figure,
    save_field_plot,
    save_field_plot_hybrid_svg,
    save_field_plot_svg,
)


FIELD_COORDINATES = np.array(
    [
        [-1.0, 0.0],
        [-0.5, -0.866],
        [0.5, -0.866],
        [1.0, 0.0],
        [0.5, 0.866],
        [-0.5, 0.866],
        [0.0, 0.0],
    ]
)
FIELD_VALUES = np.array([-1.0, -0.5, 0.0, 1.0, 0.5, -0.25, 0.75])


def _svg_viewbox_size(svg: str) -> tuple[float, float]:
    match = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg)
    assert match is not None
    return float(match.group(1)), float(match.group(2))


def test_comsol_reference_colormaps_match_supplied_colorbars():
    assert WAVE_COLORMAP(0.0) == pytest.approx(to_rgba("#3D1F82"), abs=0.01)
    assert WAVE_COLORMAP(0.5) == pytest.approx(to_rgba("#F2F1F2"), abs=0.01)
    assert WAVE_COLORMAP(1.0) == pytest.approx(to_rgba("#680A22"), abs=0.01)
    assert HEAT_CAMERA_COLORMAP(0.0) == pytest.approx(
        to_rgba("#1B0959"),
        abs=0.01,
    )
    assert HEAT_CAMERA_COLORMAP(1.0) == pytest.approx(
        to_rgba("#FEF9D4"),
        abs=0.01,
    )


def test_field_plot_layout_has_metadata_unit_and_exact_colorbar_height():
    fig, axis, colorbar_axis = _create_field_figure(
        FIELD_COORDINATES,
        FIELD_VALUES,
        frequency_thz=194.912345,
        quality_factor=1234.5,
        quantity_label="Re(Hz)",
        unit_label="A/m",
        color_scale_mode="linearsymmetric",
    )
    try:
        plot_position = axis.get_position()
        colorbar_position = colorbar_axis.get_position()
        assert colorbar_position.y0 == pytest.approx(plot_position.y0, abs=1e-12)
        assert colorbar_position.height == pytest.approx(
            plot_position.height,
            abs=1e-12,
        )
        assert {text.get_text() for text in fig.texts} == {
            "Q = 1234.5",
            "frequency = 194.91 THz",
            "Re(Hz)",
            "1",
        }
        assert len({text.get_position()[1] for text in fig.texts}) == 1
        q_text = next(text for text in fig.texts if text.get_text().startswith("Q ="))
        assert q_text.get_position()[0] == pytest.approx(
            plot_position.x0 + 0.60 * plot_position.width,
            abs=1e-12,
        )
        assert colorbar_position.x0 - plot_position.x1 == pytest.approx(
            0.035,
            abs=1e-12,
        )
        assert axis.get_xlabel() == r"$x$ ($\mu$m)"
        assert axis.get_ylabel() == r"$y$ ($\mu$m)"
        assert axis.get_xticks().size > 0
        assert axis.get_yticks().size > 0
        assert set(axis.xaxis.get_ticks_direction()) == {"in"}
        assert set(axis.yaxis.get_ticks_direction()) == {"in"}
        assert all(spine.get_edgecolor() == (0.0, 0.0, 0.0, 1.0) for spine in axis.spines.values())
    finally:
        plt.close(fig)


def test_saved_field_plot_has_transparent_canvas(tmp_path: Path):
    path = tmp_path / "field.png"
    save_field_plot(
        path,
        FIELD_COORDINATES,
        np.abs(FIELD_VALUES),
        frequency_thz=194.9,
        quality_factor=1000.0,
        quantity_label="Wem",
        unit_label="J/m³",
        color_scale_mode="linear",
        dpi=80,
    )

    image = mpimg.imread(path)
    assert image.shape[2] == 4
    assert image.shape[0] < 6 * 80
    assert image.shape[1] < 7 * 80
    assert image[0, 0, 3] == pytest.approx(0.0)
    assert np.any(image[:, :, 3] == 0.0)
    assert np.any(image[:, :, 3] == 1.0)


def test_saved_svg_field_plot_is_fully_vector(tmp_path: Path):
    path = tmp_path / "field.svg"
    save_field_plot_svg(
        path,
        FIELD_COORDINATES,
        FIELD_VALUES,
        frequency_thz=194.9,
        quality_factor=1000.0,
        quantity_label="Re(Hz)",
        unit_label="A/m",
        color_scale_mode="linearsymmetric",
    )

    svg = path.read_text(encoding="utf-8")
    assert "<svg" in svg
    assert "<image" not in svg
    width, height = _svg_viewbox_size(svg)
    assert width < 7 * 72
    assert height < 6 * 72
    assert "frequency = 194.90 THz" in svg
    assert "Q = 1000" in svg


def test_saved_hybrid_svg_rasterizes_only_the_field(tmp_path: Path):
    path = tmp_path / "field_hybrid.svg"
    save_field_plot_hybrid_svg(
        path,
        FIELD_COORDINATES,
        FIELD_VALUES,
        frequency_thz=194.9,
        quality_factor=1000.0,
        quantity_label="Re(Hz)",
        unit_label="A/m",
        color_scale_mode="linearsymmetric",
        dpi=300,
    )

    svg = path.read_text(encoding="utf-8")
    assert svg.count("<image") == 1
    width, height = _svg_viewbox_size(svg)
    assert width < 7 * 72
    assert height < 6 * 72
    assert "frequency = 194.90 THz" in svg
    assert "Q = 1000" in svg


@pytest.mark.parametrize("label", ["Re(Hz)", "Im(Hz)"])
@pytest.mark.parametrize("values", [FIELD_VALUES * 37, np.zeros(7), np.ones(7) * 4])
def test_hz_reference_colorbar_normalizes_without_mutating_data(label, values):
    original = values.copy()
    fig, axis, bar = _create_field_figure(
        FIELD_COORDINATES, values, frequency_thz=194.9, quality_factor=1000,
        quantity_label=label, unit_label="A/m", color_scale_mode="linearsymmetric",
    )
    try:
        field = axis.collections[0]
        assert field.cmap.name == "RdBu_r"
        assert field.cmap(0.0) == pytest.approx(to_rgba("#053061"))
        assert field.cmap(1.0) == pytest.approx(to_rgba("#67001f"))
        assert field.get_clim() == (-1.0, 1.0)
        np.testing.assert_allclose(bar.get_yticks(), np.linspace(-1, 1, 9))
        assert [tick.get_text().replace("−", "-") for tick in bar.get_yticklabels()] == [
            f"{tick:.2f}" for tick in np.linspace(-1, 1, 9)
        ]
        peak = np.max(np.abs(original))
        np.testing.assert_allclose(field.get_array(), original / peak if peak else original)
        np.testing.assert_array_equal(values, original)
        assert label in {text.get_text() for text in fig.texts}
        assert bar.get_ylabel() == "normalized"
        assert "A/m" not in {text.get_text() for text in fig.texts}
    finally:
        plt.close(fig)


def test_hz_nonfinite_samples_are_excluded_before_normalization():
    coordinates = np.vstack([FIELD_COORDINATES, [np.nan, 0], [2, 2]])
    values = np.r_[FIELD_VALUES * 5, 1e8, np.nan]
    fig, axis, _bar = _create_field_figure(
        coordinates, values, frequency_thz=194.9, quality_factor=1000,
        quantity_label="Im(Hz)", unit_label="A/m", color_scale_mode="linearsymmetric",
    )
    try:
        np.testing.assert_allclose(axis.collections[0].get_array(), FIELD_VALUES)
    finally:
        plt.close(fig)
