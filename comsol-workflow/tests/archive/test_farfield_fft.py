import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import comsol_workflow.farfield_fft as farfield_fft
from comsol_workflow.farfield_fft import (
    brillouin_zone_polygons,
    compute_farfield_fft,
    gaussian_fwhm_metrics,
    load_air_field,
    plot_farfield_figures,
    polygon_axis_interval,
    run_farfield_fft_case,
    summarize_farfield,
    symmetric_pad,
    validate_farfield_config,
)


def test_default_odd_sizes_keep_exact_centers():
    validate_farfield_config(801, 2001, 0.9, None)
    padded = symmetric_pad(np.ones((801, 801), dtype=complex), 2001)

    assert padded.shape == (2001, 2001)
    assert padded[1000, 1000] == 1
    assert np.count_nonzero(padded[:600]) == 0
    assert np.count_nonzero(padded[:, :600]) == 0


@pytest.mark.parametrize(
    ("grid_size", "fft_size", "na", "h0_um"),
    [
        (800, 2001, 0.9, None),
        (801, 2000, 0.9, None),
        (803, 801, 0.9, None),
        (801, 2001, 0.0, None),
        (801, 2001, 1.1, None),
        (801, 2001, 0.9, -1.0),
    ],
)
def test_invalid_farfield_config_is_rejected(grid_size, fft_size, na, h0_um):
    with pytest.raises(ValueError):
        validate_farfield_config(grid_size, fft_size, na, h0_um)


def test_centered_plane_wave_hits_exact_fft_bin():
    size = 9
    axis = np.linspace(-2.0, 2.0, size)
    dx = axis[1] - axis[0]
    target_bin = 2
    k = 2 * np.pi * target_bin / (size * dx)
    xx, _yy = np.meshgrid(axis, axis, indexing="xy")
    ex = np.exp(1j * k * xx)
    zero = np.zeros_like(ex)

    result = compute_farfield_fft(
        ex,
        zero,
        zero,
        axis,
        axis,
        frequency_thz=200.0,
        h0_um=0.0,
        na=1.0,
        fft_size=size,
    )

    peak = np.unravel_index(np.argmax(result["normfE0"]), result["normfE0"].shape)
    assert peak == (4, 6)
    assert result["kx"][4] == pytest.approx(0.0)
    assert result["ky"][4] == pytest.approx(0.0)
    assert result["na_mask"][4, 4]


def test_brillouin_zone_polygons_include_central_and_six_neighbours():
    period_um = 0.82

    polygons = brillouin_zone_polygons(period_um)

    assert len(polygons) == 7
    assert np.mean(polygons[0], axis=0) == pytest.approx([0.0, 0.0], abs=1e-12)
    neighbour_centres = np.array([np.mean(polygon, axis=0) for polygon in polygons[1:]])
    expected_radius = 4.0 * np.pi / (np.sqrt(3.0) * period_um)
    assert np.linalg.norm(neighbour_centres, axis=1) == pytest.approx(
        np.full(6, expected_radius)
    )


def test_central_bz_axis_intervals_match_triangular_lattice():
    period_um = 0.82
    central_bz = brillouin_zone_polygons(period_um)[0]

    kx_min, kx_max = polygon_axis_interval(central_bz, "kx")
    ky_min, ky_max = polygon_axis_interval(central_bz, "ky")

    assert (kx_min, kx_max) == pytest.approx(
        (-4.0 * np.pi / (3.0 * period_um), 4.0 * np.pi / (3.0 * period_um))
    )
    assert (ky_min, ky_max) == pytest.approx(
        (-2.0 * np.pi / (np.sqrt(3.0) * period_um), 2.0 * np.pi / (np.sqrt(3.0) * period_um))
    )


def test_gaussian_fwhm_metrics_convert_k_width_to_full_angle():
    fit_b = 0.2
    fit_c = 0.4
    k0 = 4.0
    half_width = abs(fit_c) * np.sqrt(np.log(2.0))

    metrics = gaussian_fwhm_metrics(fit_b, fit_c, k0)

    assert metrics["fwhm_k_um_inv"] == pytest.approx(2.0 * half_width)
    assert metrics["lower_k_um_inv"] == pytest.approx(fit_b - half_width)
    assert metrics["upper_k_um_inv"] == pytest.approx(fit_b + half_width)
    expected_angle = np.degrees(
        np.arcsin((fit_b + half_width) / k0) - np.arcsin((fit_b - half_width) / k0)
    )
    assert metrics["divergence_full_angle_deg"] == pytest.approx(expected_angle)


def test_farfield_plots_label_physical_and_normalized_axes(tmp_path, monkeypatch):
    from matplotlib.colors import to_rgba

    size = 33
    axis = np.linspace(-2.0, 2.0, size)
    xx, yy = np.meshgrid(axis, axis, indexing="xy")
    ex = np.exp(-0.5 * (xx**2 + yy**2)).astype(complex)
    ey = 0.7 * np.exp(-0.5 * ((xx - 0.5) ** 2 + (yy + 0.25) ** 2)).astype(complex)
    zero = np.zeros_like(ex)
    result = compute_farfield_fft(
        ex,
        ey,
        zero,
        axis,
        axis,
        frequency_thz=200.0,
        h0_um=0.0,
        na=0.9,
        fft_size=size,
    )
    summary = summarize_farfield(
        result,
        source_file="synthetic",
        mode_idx=0,
        frequency_thz=200.0,
        q=1000.0,
        h0_um=0.0,
        na=0.9,
        period_um=0.82,
    )
    original_close = farfield_fft.plt.close
    monkeypatch.setattr(farfield_fft.plt, "close", lambda _fig: None)

    plot_farfield_figures(result, summary, tmp_path, period_um=0.82, dpi=40)
    figures = [farfield_fft.plt.figure(number) for number in farfield_fft.plt.get_fignums()]
    xlabels = {axis.get_xlabel() for figure in figures for axis in figure.axes}
    ylabels = {axis.get_ylabel() for figure in figures for axis in figure.axes}
    overview_axes = figures[0].axes[:4]
    overview_colorbar_axes = figures[0].axes[4:]
    near_axis = overview_axes[0]
    zoom_linear_axis = overview_axes[1]
    first_bz_log_axis = overview_axes[2]
    neighbouring_bz_log_axis = overview_axes[3]

    assert {"x (μm)", "kx/(2π/a)", "ky/(2π/a)", "θ (deg)"} <= xlabels
    assert {
        "y (μm)",
        "ky/(2π/a)",
        "Normalized transmitted power",
        "Normalized |E(k)|²",
        "Normalized |E(kx=0, ky)|²",
    } <= ylabels
    assert [axis.get_title() for axis in overview_axes] == [
        "Near field intensity (NA=0.90)",
        "Far field Intensity (k-space, linear)",
        "Far field intensity (1st BZ, log)",
        "Far field Intensity (k-space, log)",
    ]
    assert figures[0]._suptitle.get_text() == "f=200.000 THz, Q=1000.00"
    assert [axis.images[0].get_cmap().name for axis in overview_axes] == ["magma"] * 4
    assert [axis.get_ylabel() for axis in overview_colorbar_axes] == [
        "Normalized Intensity"
    ] * 4
    colorbar_label_rotations = [
        axis.yaxis.label.get_rotation() for axis in overview_colorbar_axes
    ]
    assert colorbar_label_rotations == pytest.approx([270.0] * 4)
    assert near_axis.get_xlim() == pytest.approx(
        (result["x_axis"][0], result["x_axis"][-1])
    )
    assert near_axis.get_ylim() == pytest.approx(
        (result["y_axis"][0], result["y_axis"][-1])
    )
    assert near_axis.images[0].get_clim() == pytest.approx((0.0, 1.0))
    assert zoom_linear_axis.get_xlim() == pytest.approx((-0.2, 0.2))
    assert zoom_linear_axis.get_ylim() == pytest.approx((-0.2, 0.2))
    assert zoom_linear_axis.images[0].get_clim() == pytest.approx((0.0, 1.0))
    assert not np.ma.getmaskarray(zoom_linear_axis.images[0].get_array()).any()
    k_scale = 2.0 * np.pi / 0.82
    image_extent = zoom_linear_axis.images[0].get_extent()
    kx_keep = (result["kx"] >= image_extent[0] * k_scale) & (
        result["kx"] <= image_extent[1] * k_scale
    )
    ky_keep = (result["ky"] >= image_extent[2] * k_scale) & (
        result["ky"] <= image_extent[3] * k_scale
    )
    expected_linear = farfield_fft._normalized(result["normfE0"])[np.ix_(ky_keep, kx_keep)]
    assert np.asarray(zoom_linear_axis.images[0].get_array()) == pytest.approx(expected_linear)
    assert not zoom_linear_axis.lines
    assert not zoom_linear_axis.patches
    assert len(first_bz_log_axis.lines) == 1
    assert len(neighbouring_bz_log_axis.lines) == 7
    for axis in (first_bz_log_axis, neighbouring_bz_log_axis):
        assert len(axis.patches) == 1
        assert axis.patches[0].get_linestyle() == "--"
        assert axis.patches[0].get_edgecolor() == pytest.approx(to_rgba("#9ca3af"))
        assert axis.images[0].get_clim() == pytest.approx((-5.0, 0.0))
    expected_log_ticklabels = [
        r"$10^{-5}$",
        r"$10^{-4}$",
        r"$10^{-3}$",
        r"$10^{-2}$",
        r"$10^{-1}$",
        r"$10^{0}$",
    ]
    for colorbar_axis in overview_colorbar_axes[2:]:
        assert [tick.get_text() for tick in colorbar_axis.get_yticklabels()] == (
            expected_log_ticklabels
        )
    assert [axis.get_box_aspect() for axis in overview_axes] == pytest.approx([1.0] * 4)
    for axis in (first_bz_log_axis, neighbouring_bz_log_axis):
        image_extent = axis.images[0].get_extent()
        assert axis.get_xlim() == pytest.approx(image_extent[:2])
        assert axis.get_ylim() == pytest.approx(image_extent[2:])
        assert axis.get_aspect() == pytest.approx(1.0)
        assert np.ptp(axis.get_xlim()) == pytest.approx(np.ptp(axis.get_ylim()))
        margin = 0.01 * np.ptp(axis.get_xlim())
        for boundary in axis.lines:
            assert min(boundary.get_xdata()) > axis.get_xlim()[0] + margin
            assert max(boundary.get_xdata()) < axis.get_xlim()[1] - margin
            assert min(boundary.get_ydata()) > axis.get_ylim()[0] + margin
            assert max(boundary.get_ydata()) < axis.get_ylim()[1] - margin
    polar_axis = figures[1].axes[0]
    assert polar_axis.get_theta_offset() == pytest.approx(np.pi / 2.0)
    for axis in figures[2].axes[:2]:
        assert axis.get_xlim() == pytest.approx((-0.2, 0.2))
    gauss_annotation = "\n".join(text.get_text() for text in figures[3].axes[0].texts)
    assert "FWHM_k" in gauss_annotation
    assert "FWHM divergence" in gauss_annotation
    gauss_axis = figures[3].axes[0]
    assert gauss_axis.get_xlim() == pytest.approx((-0.2, 0.2))
    metrics_text = next(text for text in gauss_axis.texts if "FWHM_k" in text.get_text())
    assert metrics_text.get_fontsize() <= 8.0
    fwhm_arrows = [
        annotation
        for annotation in gauss_axis.texts
        if getattr(annotation, "arrow_patch", None) is not None
    ]
    assert len(fwhm_arrows) == 2
    fwhm_arrows.sort(key=lambda annotation: annotation.xyann[0])
    left_arrow, right_arrow = fwhm_arrows
    assert left_arrow.xyann[0] < left_arrow.xy[0] < 0.0
    assert right_arrow.xyann[0] > right_arrow.xy[0] > 0.0
    assert left_arrow.xy[0] - left_arrow.xyann[0] == pytest.approx(0.035)
    assert right_arrow.xyann[0] - right_arrow.xy[0] == pytest.approx(0.035)
    assert left_arrow.xyann[1] == pytest.approx(left_arrow.xy[1])
    assert right_arrow.xyann[1] == pytest.approx(right_arrow.xy[1])
    original_close("all")


def _write_air_field(path: Path, size: int = 3) -> None:
    axis = np.linspace(-1.0, 1.0, size)
    xx, yy = np.meshgrid(axis, axis, indexing="xy")
    ones = np.ones(xx.size)
    pd.DataFrame(
        {
            "x": xx.ravel(),
            "y": yy.ravel(),
            "Ex_re": ones,
            "Ex_im": np.zeros(xx.size),
            "Ey_re": np.zeros(xx.size),
            "Ey_im": np.zeros(xx.size),
            "Ez_re": np.zeros(xx.size),
            "Ez_im": np.zeros(xx.size),
            "is_in_domain": np.ones(xx.size, dtype=bool),
        }
    ).to_parquet(path, index=False)


def test_load_air_field_restores_c_order_square_grid(tmp_path):
    path = tmp_path / "00_E_air.parquet"
    _write_air_field(path)

    field = load_air_field(path)

    assert field["Ex"].shape == (3, 3)
    assert np.all(field["Ex"] == 1.0)
    assert field["x_axis"].tolist() == [-1.0, 0.0, 1.0]
    assert field["y_axis"].tolist() == [-1.0, 0.0, 1.0]


def test_case_with_no_valid_modes_writes_empty_summary(tmp_path):
    case_dir = tmp_path / "case"
    export_dir = case_dir / "simulation_exports"
    export_dir.mkdir(parents=True)
    pd.DataFrame(
        {"re": [194.0], "im": [0.1], "q": [970.0], "is_valid": [False], "mode_idx": [0]}
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    (export_dir / "air_field_metadata.json").write_text(
        json.dumps({"air_plane_z_um": 1.75, "grid_size": 3}),
        encoding="utf-8",
    )

    summary = run_farfield_fft_case(
        case_dir,
        export_dir,
        grid_size=3,
        fft_size=3,
        na=0.9,
        h0_um=None,
        period_um=0.82,
        plot=False,
    )

    assert summary.empty
    assert (case_dir / "farfield_fft" / "farfield_summary.csv").exists()


def test_case_mode_subset_only_runs_selected_valid_mode(tmp_path, monkeypatch):
    case_dir = tmp_path / "case"
    export_dir = case_dir / "simulation_exports"
    export_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "re": [194.0, 195.0],
            "im": [0.1, 0.1],
            "q": [970.0, 980.0],
            "is_valid": [True, True],
            "mode_idx": [0, 1],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    (export_dir / "air_field_metadata.json").write_text(
        json.dumps({"air_plane_z_um": 1.75, "grid_size": 3}),
        encoding="utf-8",
    )
    calls = []

    def capture_file(_path, _output_dir, **kwargs):
        calls.append(kwargs["mode_idx"])
        return {"status": "ok", "mode_idx": kwargs["mode_idx"]}

    monkeypatch.setattr(
        "comsol_workflow.farfield_fft.run_farfield_fft_file",
        capture_file,
    )

    summary = run_farfield_fft_case(
        case_dir,
        export_dir,
        grid_size=3,
        fft_size=3,
        na=0.9,
        h0_um=None,
        period_um=0.82,
        mode_indices=[1],
        plot=False,
    )

    assert calls == [1]
    assert summary["mode_idx"].tolist() == [1]


def test_case_records_error_before_raising_when_all_valid_modes_fail(tmp_path):
    case_dir = tmp_path / "case"
    export_dir = case_dir / "simulation_exports"
    export_dir.mkdir(parents=True)
    pd.DataFrame(
        {"re": [194.0], "im": [0.1], "q": [970.0], "is_valid": [True], "mode_idx": [0]}
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    (export_dir / "air_field_metadata.json").write_text(
        json.dumps({"air_plane_z_um": 1.75, "grid_size": 3}),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="failed for all"):
        run_farfield_fft_case(
            case_dir,
            export_dir,
            grid_size=3,
            fft_size=3,
            na=0.9,
            h0_um=None,
            period_um=0.82,
            plot=False,
        )

    summary = pd.read_csv(case_dir / "farfield_fft" / "farfield_summary.csv")
    assert summary.loc[0, "status"] == "error"
    assert "00_E_air.parquet" in summary.loc[0, "error"]
