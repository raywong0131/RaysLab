import json
import inspect
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Bbox

from comsol_workflow.unit_cell_2d_analysis import (
    C2vQuarterFieldView,
    _cell_edges,
    _minor_axis_arrow_geometry,
    _polarization_ellipse_curve,
    compute_winding,
    derive_polarization_quantities,
    field_angles,
    ordered_square_loop,
    plot_component_magnitude_maps,
    plot_cy_magnitude_log_map,
    plot_intensity_maps,
    plot_polarization_angle_map,
    plot_vector_map,
    plot_winding_detail,
    redraw_cy_magnitude_log,
    run_analysis,
    scan_shrinking_loops,
    write_json_atomic,
)


AXIS = np.array([-0.1, -0.02, -0.001, 0.0, 0.001, 0.02, 0.1])


def quarter_radial_grid() -> pd.DataFrame:
    rows = []
    for qy in (0.0, 0.05, 0.1):
        for qx in (0.0, 0.05, 0.1):
            rows.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "band_label": "p2",
                    "match_status": "gamma" if qx == qy == 0.0 else "matched",
                    "cx_raw_re": qx,
                    "cx_raw_im": 0.25 * qx,
                    "cy_raw_re": qy,
                    "cy_raw_im": 0.25 * qy,
                    "cx_norm_re": qx,
                    "cx_norm_im": 0.25 * qx,
                    "cy_norm_re": qy,
                    "cy_norm_im": 0.25 * qy,
                }
            )
    return derive_polarization_quantities(pd.DataFrame(rows))


def test_quarter_view_applies_c2v_without_materializing_rows():
    source = quarter_radial_grid()
    original = source.copy(deep=True)
    view = C2vQuarterFieldView(source)
    assert len(view.source) == 9
    assert view.qx_axis == pytest.approx(np.linspace(-0.1, 0.1, 5))
    cx, cy = view.complex_grids()
    assert cx[0, 0] == pytest.approx(-0.1 - 0.025j)
    assert cy[0, 0] == pytest.approx(-0.1 - 0.025j)
    assert cx[-1, 0] == pytest.approx(-0.1 - 0.025j)
    assert cy[-1, 0] == pytest.approx(0.1 + 0.025j)
    assert view.grid_matrix("S0")[0, 0] == pytest.approx(
        view.grid_matrix("S0")[-1, -1]
    )
    assert view.grid_matrix("S2")[0, -1] == pytest.approx(
        -view.grid_matrix("S2")[-1, -1]
    )
    pd.testing.assert_frame_equal(source, original)


def test_quarter_stokes_angle_uses_axis_branch_not_interior_mirror_sign():
    view = C2vQuarterFieldView(quarter_radial_grid())
    angles = view.grid_matrix("psi_stokes")
    origin = len(view.qx_axis) // 2

    # Interior points retain the odd director-angle parity of one mirror.
    assert angles[-1, 0] == pytest.approx(-angles[-1, -1])
    assert angles[0, -1] == pytest.approx(-angles[-1, -1])

    # On either invariant axis, the negative half-axis uses the same director
    # branch as the independently solved positive half-axis.  This prevents
    # cyclic interpolation from turning an axis-only branch choice into a
    # visible phase seam.
    assert angles[:origin, origin] == pytest.approx(angles[:origin:-1, origin])
    assert angles[origin, :origin] == pytest.approx(angles[origin, :origin:-1])


def test_quarter_view_builds_closed_winding_without_expanded_table():
    view = C2vQuarterFieldView(quarter_radial_grid())
    loop = ordered_square_loop(view, 0.1)
    result = compute_winding(loop, 0.1, "complex-axis")
    assert result.winding == pytest.approx(1.0)
    assert len(view.source) == 9
    assert len(loop) == 16


def test_quarter_run_analysis_keeps_only_real_source_rows(tmp_path):
    series = tmp_path / "quarter-series"
    (series / "80_logs").mkdir(parents=True)
    (series / "99_config").mkdir(parents=True)
    source = quarter_radial_grid()
    source.to_csv(series / "80_logs" / "polarization_grid.csv", index=False)
    (series / "99_config" / "config.json").write_text(
        json.dumps(
            {
                "cache_identity": {
                    "isQuarter": 1,
                    "unit_cell_2d": {"isQuarter": 1},
                }
            }
        ),
        encoding="utf-8",
    )
    summary = run_analysis(series)
    saved = pd.read_csv(series / "80_logs" / "polarization_grid.csv")
    render = json.loads(
        (series / "99_config" / "symmetry_render_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert len(saved) == 9
    assert render["solved_coordinate_count"] == 9
    assert render["plot_coordinate_count"] == 25
    assert render["materialized_symmetry_row_count"] == 0
    assert summary["status"] == "ok"


def synthetic_grid(*, reverse_y=False, constant=False):
    rows = []
    for qy in AXIS:
        for qx in AXIS:
            if constant:
                cx, cy = 1.0 + 1.0j, 0.25 + 0.25j
            else:
                cx = (1.0 + 1.0j) * qx
                cy = (1.0 + 1.0j) * (-qy if reverse_y else qy)
            rows.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "band_label": "p2",
                    "match_status": "gamma" if qx == qy == 0 else "matched",
                    "cx_raw_re": cx.real,
                    "cx_raw_im": cx.imag,
                    "cy_raw_re": cy.real,
                    "cy_raw_im": cy.imag,
                    "cx_norm_re": cx.real,
                    "cx_norm_im": cx.imag,
                    "cy_norm_re": cy.real,
                    "cy_norm_im": cy.imag,
                }
            )
    return derive_polarization_quantities(pd.DataFrame(rows))


def uniform_grid(point_count=5):
    axis = np.linspace(-0.1, 0.1, point_count)
    rows = []
    for qy in axis:
        for qx in axis:
            cx = 1.0 + 0.1j
            cy = (0.12 + 0.6 * (qx - 0.25 * qy) ** 2) + 0.02j
            rows.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "band_label": "p2",
                    "match_status": "gamma" if qx == qy == 0 else "matched",
                    "cx_raw_re": cx.real,
                    "cx_raw_im": cx.imag,
                    "cy_raw_re": cy.real,
                    "cy_raw_im": cy.imag,
                    "cx_norm_re": cx.real,
                    "cx_norm_im": cx.imag,
                    "cy_norm_re": cy.real,
                    "cy_norm_im": cy.imag,
                }
            )
    return derive_polarization_quantities(pd.DataFrame(rows))


@pytest.mark.parametrize("method", ["real", "imag", "complex-axis", "stokes-axis"])
def test_radial_field_has_positive_unit_winding(method):
    frame = synthetic_grid()
    loop = ordered_square_loop(frame, 0.1)
    result = compute_winding(loop, 0.1, method)
    assert result.nearest_integer == 1
    assert result.winding == pytest.approx(1.0, abs=1e-12)


def test_reverse_and_constant_fields_have_expected_winding():
    reverse = compute_winding(
        ordered_square_loop(synthetic_grid(reverse_y=True), 0.1),
        0.1,
        "stokes-axis",
    )
    constant = compute_winding(
        ordered_square_loop(synthetic_grid(constant=True), 0.1),
        0.1,
        "stokes-axis",
    )
    assert reverse.nearest_integer == -1
    assert constant.nearest_integer == 0


def test_stokes_angle_is_invariant_under_global_complex_phase():
    frame = synthetic_grid()
    cx = frame["cx_norm_re"].to_numpy() + 1j * frame["cx_norm_im"].to_numpy()
    cy = frame["cy_norm_re"].to_numpy() + 1j * frame["cy_norm_im"].to_numpy()
    original, *_ = field_angles(cx, cy, "stokes-axis")
    phase = np.exp(0.713j)
    shifted, *_ = field_angles(cx * phase, cy * phase, "stokes-axis")
    assert np.allclose(original, shifted, atol=1e-12)


def test_primary_intensity_uses_raw_projection_while_stokes_uses_normalized():
    frame = synthetic_grid()
    frame[["cx_raw_re", "cx_raw_im", "cy_raw_re", "cy_raw_im"]] *= 3.0
    derived = derive_polarization_quantities(frame)
    assert np.allclose(derived["C"], 9.0 * derived["C_norm"])


def test_square_loop_is_counter_clockwise_and_shrinks_stably():
    frame = synthetic_grid()
    loop = ordered_square_loop(frame, 0.1)
    x = loop["qx_over_G"].to_numpy()
    y = loop["qy_over_G"].to_numpy()
    signed_area = 0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))
    assert signed_area > 0
    scan, stable = scan_shrinking_loops(frame, "stokes-axis")
    assert stable is not None
    assert stable.half_width == pytest.approx(0.1)
    assert stable.nearest_integer == 1
    assert set(scan["status"]) == {"ok"}


def test_isolated_outer_winding_is_not_reported_as_multiscale_stable():
    frame = synthetic_grid()
    inner = np.maximum(np.abs(frame["qx_over_G"]), np.abs(frame["qy_over_G"])) <= 0.02
    for prefix, value in (("cx", 1.0 + 1.0j), ("cy", 0.25 + 0.25j)):
        for kind in ("raw", "norm"):
            frame.loc[inner, f"{prefix}_{kind}_re"] = value.real
            frame.loc[inner, f"{prefix}_{kind}_im"] = value.imag
    frame = derive_polarization_quantities(frame)
    scan, stable = scan_shrinking_loops(frame, "stokes-axis")
    assert scan.iloc[0]["nearest_integer"] == 1
    assert stable is None


def test_nonuniform_cell_edges_use_coordinate_midpoints():
    assert np.allclose(
        _cell_edges(np.array([-0.1, -0.01, 0.0, 0.03])),
        [-0.145, -0.055, -0.005, 0.015, 0.045],
    )


def test_log_map_uses_equal_aspect_magma_and_power_of_ten_ticks(
    tmp_path, monkeypatch
):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    plot_intensity_maps(synthetic_grid(), tmp_path, "p2")
    figure = plt.figure(max(plt.get_fignums()))
    figure.canvas.draw()
    image = figure.axes[0].collections[0]
    labels = [label.get_text() for label in figure.axes[1].get_yticklabels()]
    map_box = figure.axes[0].get_position()
    colorbar_box = figure.axes[1].get_position()
    assert figure.axes[0].get_aspect() == 1.0
    assert colorbar_box.y0 == pytest.approx(map_box.y0, abs=1e-12)
    assert colorbar_box.y1 == pytest.approx(map_box.y1, abs=1e-12)
    expected_smooth_axis_count = (len(AXIS) - 1) * 24 + 1
    assert image.get_array().size == expected_smooth_axis_count**2
    label_box = figure.axes[1].yaxis.label.get_window_extent(
        renderer=figure.canvas.get_renderer()
    )
    assert label_box.x0 >= figure.bbox.x0
    assert label_box.x1 <= figure.bbox.x1
    assert "py" in figure.axes[0].get_title()
    assert "p2" not in figure.axes[0].get_title()
    assert figure.axes[0].get_title() == (
        r"Normalized $|c_x|^2+|c_y|^2$ (py)"
    )
    assert figure.axes[1].yaxis.label.get_text() == "Normalized intensity"
    assert image.cmap.name == "magma"
    assert labels
    assert all("10^{" in label for label in labels if label)
    for number in list(plt.get_fignums()):
        original_close(plt.figure(number))


def test_quarter_a_maps_plot_normalized_c_magnitude(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    outputs = plot_intensity_maps(view, tmp_path, "p2")
    figures = [plt.figure(number) for number in sorted(plt.get_fignums())]
    linear_values = np.ma.asarray(
        figures[0].axes[0].collections[0].get_array()
    ).reshape(97, 97)
    assert [output.name for output in outputs] == [
        "A1_ck_intensity_linear.png",
        "A2_ck_intensity_log.png",
    ]
    assert linear_values[48, 72] == pytest.approx(1.0 / math.sqrt(8.0))
    for figure in figures:
        assert figure.axes[0].get_title() == r"Normalized $|c(k)|$ (py)"
        assert figure.axes[1].yaxis.label.get_text() == r"Normalized $|c(k)|$"
        original_close(figure)


def test_component_maps_share_max_c_denominator_without_changing_color_range(
    tmp_path, monkeypatch
):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    plot_component_magnitude_maps(synthetic_grid(), tmp_path, "p2")
    figures = [plt.figure(number) for number in sorted(plt.get_fignums())]
    cx_linear = figures[0].axes[0].collections[0]
    cy_linear = figures[2].axes[0].collections[0]
    expected_component_max = 1.0 / math.sqrt(2.0)
    assert cx_linear.norm.vmax == pytest.approx(expected_component_max)
    assert cy_linear.norm.vmax == pytest.approx(expected_component_max)
    assert np.max(cx_linear.norm(cx_linear.get_array())) == pytest.approx(1.0)
    assert np.max(cy_linear.norm(cy_linear.get_array())) == pytest.approx(1.0)
    assert "Normalized $|c_x(k)|$" in figures[0].axes[0].get_title()
    assert "Normalized $|c_y(k)|$" in figures[2].axes[0].get_title()
    for figure in figures:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        image = figure.axes[0].collections[0]
        colorbar_axis = figure.axes[1]
        assert max(colorbar_axis.get_yticks()) == pytest.approx(image.norm.vmax)
        tick_boxes = [
            label.get_window_extent(renderer=renderer)
            for label in colorbar_axis.get_yticklabels()
            if label.get_visible() and label.get_text()
        ]
        tick_box = Bbox.union(tick_boxes)
        label_box = figure.axes[1].yaxis.label.get_window_extent(
            renderer=renderer
        )
        assert label_box.x0 - tick_box.x1 >= 8.0
        assert label_box.x1 <= figure.bbox.x1
        decorated = Bbox.union(
            [
                figure.axes[0].get_tightbbox(renderer),
                colorbar_axis.get_tightbbox(renderer),
            ]
        )
        assert 0.5 * (decorated.x0 + decorated.x1) == pytest.approx(
            0.5 * figure.bbox.width, abs=1.0
        )
        original_close(figure)


def test_b4_uses_same_raw_component_contract_as_other_b_panels(
    tmp_path, monkeypatch
):
    axis = np.linspace(-0.15, 0.15, 25)
    rows = []
    for qy in axis:
        for qx in axis:
            rows.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "band_label": "p2",
                    "match_status": "gamma" if qx == qy == 0 else "matched",
                    "cx_raw_re": 1.0,
                    "cx_raw_im": 0.0,
                    "cy_raw_re": abs(qx + qy) + 1e-6,
                    "cy_raw_im": 0.0,
                    "cx_norm_re": 1.0,
                    "cx_norm_im": 0.0,
                    "cy_norm_re": abs(qx + qy) + 1e-6,
                    "cy_norm_im": 0.0,
                }
            )
    frame = derive_polarization_quantities(pd.DataFrame(rows))
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    output = plot_cy_magnitude_log_map(frame, tmp_path, "p2")
    figure = plt.figure(max(plt.get_fignums()))
    image = figure.axes[0].collections[0]
    assert output.name == "B4_cy_magnitude_log.png"
    assert image.get_array().size == 577**2
    assert image.cmap.name == "magma"
    common_maximum = math.hypot(1.0, 0.300001)
    component_maximum = 0.300001
    expected_maximum = component_maximum / common_maximum
    expected_floor = 1e-6 * expected_maximum
    assert image.norm.vmax == pytest.approx(expected_maximum)
    assert image.norm.vmin == pytest.approx(expected_floor)
    assert max(figure.axes[1].get_yticks()) == pytest.approx(image.norm.vmax)
    assert figure.axes[0].get_title() == "Normalized $|c_y(k)|$ (py)"
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tick_boxes = [
        label.get_window_extent(renderer=renderer)
        for label in figure.axes[1].get_yticklabels()
        if label.get_visible() and label.get_text()
    ]
    assert all(
        lower.y1 + 2.0 <= upper.y0
        for lower, upper in zip(tick_boxes, tick_boxes[1:])
    )
    original_close(figure)


def test_quarter_b4_uses_fixed_log_colorbar_scale(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    output = plot_cy_magnitude_log_map(view, tmp_path, "p2")
    figure = plt.figure(max(plt.get_fignums()))
    image = figure.axes[0].collections[0]
    colorbar_axis = figure.axes[1]
    assert output.name == "B4_cy_magnitude_log.png"
    assert image.norm.vmin == pytest.approx(1e-4)
    assert image.norm.vmax == pytest.approx(1e-1)
    assert colorbar_axis.get_yticks() == pytest.approx(
        [1e-4, 1e-3, 1e-2, 1e-1]
    )
    assert [label.get_text() for label in colorbar_axis.get_yticklabels()] == [
        r"$10^{-4}$",
        r"$10^{-3}$",
        r"$10^{-2}$",
        r"$10^{-1}$",
    ]
    original_close(figure)


def test_quarter_b2_uses_fixed_log_colorbar_scale(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    outputs = plot_component_magnitude_maps(
        view,
        tmp_path,
        "p2",
        panels={"B2"},
    )
    figure = plt.figure(max(plt.get_fignums()))
    image = figure.axes[0].collections[0]
    colorbar_axis = figure.axes[1]
    assert [output.name for output in outputs] == [
        "B2_cx_magnitude_log.png"
    ]
    assert image.norm.vmin == pytest.approx(1e-4)
    assert image.norm.vmax == pytest.approx(1.0)
    assert colorbar_axis.get_yticks() == pytest.approx(
        [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    )
    assert [label.get_text() for label in colorbar_axis.get_yticklabels()] == [
        r"$10^{-4}$",
        r"$10^{-3}$",
        r"$10^{-2}$",
        r"$10^{-1}$",
        r"$10^{0}$",
    ]
    original_close(figure)


def test_quarter_b1_uses_fixed_linear_colorbar_scale(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    outputs = plot_component_magnitude_maps(
        view,
        tmp_path,
        "p2",
        panels={"B1"},
    )
    figure = plt.figure(max(plt.get_fignums()))
    image = figure.axes[0].collections[0]
    colorbar_axis = figure.axes[1]
    assert [output.name for output in outputs] == [
        "B1_cx_magnitude_linear.png"
    ]
    assert image.norm.vmin == pytest.approx(0.0)
    assert image.norm.vmax == pytest.approx(1.0)
    assert colorbar_axis.get_yticks() == pytest.approx(
        [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    )
    original_close(figure)


def test_b4_ignores_norm_columns_like_the_other_b_panels(tmp_path, monkeypatch):
    frame = synthetic_grid()
    changed = frame.copy()
    changed["cx_norm_re"] = 1000.0
    changed["cx_norm_im"] = -500.0
    changed["cy_norm_re"] = -750.0
    changed["cy_norm_im"] = 250.0
    changed = derive_polarization_quantities(changed)
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    plot_cy_magnitude_log_map(frame, tmp_path / "one", "p2")
    first = plt.figure(max(plt.get_fignums()))
    first_values = np.asarray(first.axes[0].collections[0].get_array()).copy()
    first_norm = first.axes[0].collections[0].norm
    plot_cy_magnitude_log_map(changed, tmp_path / "two", "p2")
    second = plt.figure(max(plt.get_fignums()))
    second_values = np.asarray(second.axes[0].collections[0].get_array()).copy()
    second_norm = second.axes[0].collections[0].norm
    assert second_values == pytest.approx(first_values)
    assert second_norm.vmin == pytest.approx(first_norm.vmin)
    assert second_norm.vmax == pytest.approx(first_norm.vmax)
    for figure in (first, second):
        original_close(figure)


def test_redraw_b4_preserves_named_reference_file(tmp_path):
    series = tmp_path / "series"
    logs = series / "80_logs"
    overview = series / "10_overview"
    logs.mkdir(parents=True)
    overview.mkdir(parents=True)
    frame = synthetic_grid().drop(
        columns=[
            "C", "C_norm", "S0", "S1", "S2", "S3",
            "psi_stokes", "analysis_valid",
        ]
    )
    frame.to_csv(logs / "polarization_grid.csv", index=False)
    reference = overview / "B4_cy_magnitude_log.pre_regenerate_20260824.png"
    reference.write_bytes(b"preserved reference")
    outputs = redraw_cy_magnitude_log(series)
    assert [path.name for path in outputs] == ["B4_cy_magnitude_log.png"]
    assert reference.read_bytes() == b"preserved reference"
    assert not (overview / "B4_cy_magnitude_log_smoothing_only.png").exists()


def test_b4_implementation_contains_no_nodal_curve_prior():
    source = inspect.getsource(plot_cy_magnitude_log_map)
    source += inspect.getsource(plot_component_magnitude_maps)
    forbidden = (
        "branch_extraction",
        "signed_distance",
        "nodal_distance",
        "gamma_floor",
        "center_restore",
    )
    assert all(token not in source for token in forbidden)


def test_phase_map_uses_smooth_cyclic_rendering(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    plot_polarization_angle_map(synthetic_grid(), tmp_path, "p2", 1e-6)
    figure = plt.figure(max(plt.get_fignums()))
    image = figure.axes[0].collections[0]
    expected_smooth_axis_count = (len(AXIS) - 1) * 24 + 1
    assert image.get_array().size == expected_smooth_axis_count**2
    assert image.cmap.name == "twilight_shifted"
    assert image.cmap(0.0) == pytest.approx(image.cmap(1.0), abs=0.02)
    assert "Polarization orientation angle" in figure.axes[0].get_title()
    assert "py mode" in figure.axes[0].get_title()
    assert "stokes-axis" not in figure.axes[0].get_title()
    assert figure.axes[1].get_yticks() == pytest.approx(
        [-math.pi / 2, -math.pi / 4, 0, math.pi / 4, math.pi / 2]
    )
    original_close(figure)


def test_c1_and_d_use_arrowed_polarization_ellipses(tmp_path, monkeypatch):
    ellipse_x, ellipse_y = _polarization_ellipse_curve(1.0 + 0.0j, 1.0j, 9)
    assert ellipse_x[0] > 0.0
    assert ellipse_y[1] > ellipse_y[0]
    line_endpoint, line_tangent = _minor_axis_arrow_geometry(1.0, 0.0)
    assert line_endpoint == pytest.approx([0.0, 0.0], abs=1e-12)
    assert np.linalg.norm(line_tangent) == pytest.approx(1.0)
    minor_endpoint, tangent = _minor_axis_arrow_geometry(1.0, 0.5j)
    assert np.linalg.norm(minor_endpoint) == pytest.approx(
        0.5 / math.sqrt(1.0 + 0.5**2)
    )
    signed_rotation = (
        minor_endpoint[0] * tangent[1]
        - minor_endpoint[1] * tangent[0]
    )
    assert signed_rotation > 0.0

    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    frame = synthetic_grid()
    for kind in ("raw", "norm"):
        frame[f"cx_{kind}_re"] = 1.0
        frame[f"cx_{kind}_im"] = 0.0
        frame[f"cy_{kind}_re"] = 0.0
        frame[f"cy_{kind}_im"] = 0.5
    frame = derive_polarization_quantities(frame)
    plot_vector_map(frame, tmp_path, "p2", "complex-axis", 1e-6, None)
    c1 = plt.figure(max(plt.get_fignums()))
    assert c1.axes[0].get_xlim() == pytest.approx((-0.18, 0.18))
    assert c1.axes[0].get_ylim() == pytest.approx((-0.18, 0.18))
    assert len(c1.axes[0].lines) == len(frame)
    assert sum(
        isinstance(patch, FancyArrowPatch) for patch in c1.axes[0].patches
    ) == len(frame)
    assert "Polarization ellipses" in c1.axes[0].get_title()

    _scan, stable = scan_shrinking_loops(frame, "complex-axis")
    assert stable is not None
    plot_winding_detail(
        frame, tmp_path, "p2", "complex-axis", stable, 1e-6
    )
    winding = plt.figure(max(plt.get_fignums()))
    assert len(winding.axes) == 4
    assert "polarization ellipses" in winding.axes[0].get_title()
    assert winding.axes[1].get_title() == "Winding number = 0"
    assert winding.axes[2].get_title() == (
        r"Winding loop and polarization ellipses ($c_y \times 10$)"
    )
    assert winding.axes[3].get_title() == "Winding number = 0"
    assert all(
        "residual" not in axis.get_title().lower()
        for axis in winding.axes
    )
    assert sum(
        isinstance(patch, FancyArrowPatch)
        for patch in winding.axes[0].patches
    ) >= len(frame)
    assert sum(
        isinstance(patch, FancyArrowPatch)
        for patch in winding.axes[2].patches
    ) >= len(frame)
    assert winding.axes[0].get_xlim() == pytest.approx((-0.18, 0.18))
    assert winding.axes[0].get_ylim() == pytest.approx((-0.18, 0.18))
    assert winding.axes[2].get_xlim() == pytest.approx((-0.18, 0.18))
    assert winding.axes[2].get_ylim() == pytest.approx((-0.18, 0.18))
    for number in list(plt.get_fignums()):
        original_close(plt.figure(number))


def test_quarter_d_omits_cy_times_ten_panels(tmp_path, monkeypatch):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    stable = compute_winding(
        ordered_square_loop(view, 0.1),
        0.1,
        "complex-axis",
    )
    output = plot_winding_detail(
        view,
        tmp_path,
        "p2",
        "complex-axis",
        stable,
        1e-6,
    )
    figure = plt.figure(max(plt.get_fignums()))
    assert output.name == "D_BIC_winding.png"
    assert len(figure.axes) == 2
    assert figure.axes[0].get_title() == (
        "Winding loop and polarization ellipses"
    )
    assert figure.axes[1].get_title() == "Winding number = 1"
    assert all("times 10" not in axis.get_title() for axis in figure.axes)
    assert all(r"\times 10" not in axis.get_title() for axis in figure.axes)
    original_close(figure)


def test_quarter_c1_uses_source_axis_count_and_handedness_colors(
    tmp_path, monkeypatch
):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    source = quarter_radial_grid()
    for kind in ("raw", "norm"):
        source[f"cx_{kind}_re"] = 1.0
        source[f"cx_{kind}_im"] = 0.0
        source[f"cy_{kind}_re"] = 0.0
        source[f"cy_{kind}_im"] = 0.5
    view = C2vQuarterFieldView(derive_polarization_quantities(source))

    plot_vector_map(view, tmp_path, "p2", "complex-axis", 1e-6, None)
    figure = plt.figure(max(plt.get_fignums()))
    axis = figure.axes[0]
    colors = {line.get_color().lower() for line in axis.lines}
    assert len(axis.lines) == 3 * 3
    assert sum(isinstance(patch, FancyArrowPatch) for patch in axis.patches) == 9
    assert "#d62728" in colors
    assert "#1f77b4" in colors
    assert [text.get_text() for text in axis.get_legend().get_texts()] == [
        r"Positive ($S_3>0$, CCW)",
        r"Negative ($S_3<0$, CW)",
    ]
    original_close(figure)


def test_all_quarter_k_space_panels_use_005_ticks_with_endpoints(
    tmp_path, monkeypatch
):
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda _figure: None)
    view = C2vQuarterFieldView(quarter_radial_grid())
    plot_intensity_maps(view, tmp_path, "p2")
    plot_component_magnitude_maps(view, tmp_path, "p2")
    plot_polarization_angle_map(view, tmp_path, "p2", 1e-6)
    plot_vector_map(view, tmp_path, "p2", "complex-axis", 1e-6, None)
    stable = compute_winding(
        ordered_square_loop(view, 0.1), 0.1, "complex-axis"
    )
    plot_winding_detail(
        view, tmp_path, "p2", "complex-axis", stable, 1e-6
    )
    figures = [plt.figure(number) for number in sorted(plt.get_fignums())]
    expected = np.asarray([-0.1, -0.05, 0.0, 0.05, 0.1])
    k_axes = [figure.axes[0] for figure in figures]
    assert len(k_axes) == 9
    for axis in k_axes:
        assert axis.get_xticks() == pytest.approx(expected)
        assert axis.get_yticks() == pytest.approx(expected)
        axis.figure.canvas.draw()
        assert [label.get_text() for label in axis.get_xticklabels()] == [
            "-0.1", "-0.05", "0", "0.05", "0.1"
        ]
        assert [label.get_text() for label in axis.get_yticklabels()] == [
            "-0.1", "-0.05", "0", "0.05", "0.1"
        ]
    for figure in figures:
        original_close(figure)


def test_json_writer_replaces_non_finite_values(tmp_path):
    target = tmp_path / "summary.json"
    write_json_atomic(target, {"nan": float("nan"), "inf": np.float64(np.inf)})
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload == {"nan": None, "inf": None}


def test_run_analysis_generates_complete_figure_contract(tmp_path):
    series = tmp_path / "series"
    results = series / "01_results"
    results.mkdir(parents=True)
    synthetic_grid().drop(columns=["C", "C_norm", "S0", "S1", "S2", "S3", "psi_stokes", "analysis_valid"]).to_csv(
        results / "polarization_grid.csv", index=False
    )
    summary = run_analysis(series)
    overview = series / "10_overview"
    expected = {
        "A1_ck_intensity_linear.png",
        "A2_ck_intensity_log.png",
        "B1_cx_magnitude_linear.png",
        "B2_cx_magnitude_log.png",
        "B3_cy_magnitude_linear.png",
        "B4_cy_magnitude_log.png",
        "C1_ck_vector.png",
        "C2_ck_phase.png",
        "D_BIC_winding.png",
    }
    assert expected == {path.name for path in overview.iterdir()}
    assert summary["status"] == "ok"
    logs = series / "80_logs"
    winding = pd.read_csv(logs / "winding_summary.csv")
    assert set(winding["nearest_integer"]) == {1}
    assert all((overview / name).stat().st_size > 0 for name in expected)
    assert {
        "winding_summary.csv",
        "winding_loop_scan_py_real.csv",
        "winding_loop_scan_py_imag.csv",
        "winding_loop_scan_py_complex_axis.csv",
        "winding_loop_scan_py_stokes_axis.csv",
    } == {path.name for path in logs.iterdir()}
    assert (series / "99_config" / "winding_summary.json").is_file()


def test_run_analysis_plots_outer_loop_without_multiscale_stability(tmp_path):
    frame = synthetic_grid()
    inner = np.maximum(
        np.abs(frame["qx_over_G"]), np.abs(frame["qy_over_G"])
    ) <= 0.02
    for prefix, value in (
        ("cx", 1.0 + 1.0j),
        ("cy", 0.25 + 0.25j),
    ):
        for kind in ("raw", "norm"):
            frame.loc[inner, f"{prefix}_{kind}_re"] = value.real
            frame.loc[inner, f"{prefix}_{kind}_im"] = value.imag
    frame = derive_polarization_quantities(frame)
    _scan, stable = scan_shrinking_loops(frame, "complex-axis")
    assert stable is None

    series = tmp_path / "series"
    logs = series / "80_logs"
    logs.mkdir(parents=True)
    frame.drop(
        columns=[
            "C", "C_norm", "S0", "S1", "S2", "S3",
            "psi_stokes", "analysis_valid",
        ]
    ).to_csv(logs / "polarization_grid.csv", index=False)
    summary = run_analysis(series, angle_methods=("complex-axis",))
    assert (series / "10_overview" / "D_BIC_winding.png").is_file()
    assert summary["winding"][0]["status"] == "no_stable_loop"


def test_offline_entry_does_not_reference_mph_or_simulation_runner():
    source = Path("scripts/analysis/unit_cell_2D.py").read_text(encoding="utf-8")
    assert "import mph" not in source
    assert "SimulationRun" not in source
