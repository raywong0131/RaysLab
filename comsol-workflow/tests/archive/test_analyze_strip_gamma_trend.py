import os
import sys
import tempfile
from pathlib import Path

import matplotlib.axes
import pandas as pd

from scripts.analysis import srip1d_trend as trend


def test_shift_factors_cover_zero_to_point_zero_nine():
    expected = [idx / 100.0 for idx in range(10)]

    assert trend.SHIFT_FACTORS == expected
    assert trend.SHIFT_FACTORS[0] == 0.0
    assert trend.SHIFT_FACTORS[-1] == 0.09
    assert 0.10 not in trend.SHIFT_FACTORS


def test_select_closest_mode_can_select_frequency_above_target():
    df = pd.DataFrame({
        "mode_idx": [0, 1, 2],
        "frequency": [194.1, 194.7, 195.0],
        "gamma_k_weight_fraction": [0.2, 0.7, 0.9],
        "gamma_subspace_p": [0.995, 0.996, 0.997],
    })

    row = trend.select_closest_mode(df, target_frequency=194.93, source="case")

    assert int(row["mode_idx"]) == 2
    assert float(row["frequency"]) == 195.0
    assert abs(float(row["frequency_distance_to_target"]) - 0.07) < 1e-12
    assert float(row["gamma_k_weight_fraction"]) == 0.9


def test_select_closest_mode_filters_gamma_subspace_before_frequency_distance():
    df = pd.DataFrame({
        "mode_idx": [1, 2, 3],
        "frequency": [195.3927163474688, 195.9291336575822, 196.05676913004],
        "gamma_k_weight_fraction": [0.2711032133103852, 0.1199954982106016, 0.9765314757633616],
        "gamma_subspace_p": [0.9951073664285848, 0.0002762894715036, 0.9955824744603105],
    })

    row = trend.select_closest_mode(df, target_frequency=195.902486611, source="shift_0.09")

    assert int(row["mode_idx"]) == 3
    assert float(row["gamma_subspace_p"]) > 0.99


def test_select_closest_mode_requires_gamma_subspace_strictly_above_threshold():
    df = pd.DataFrame({
        "mode_idx": [0, 1],
        "frequency": [194.93, 194.8],
        "gamma_k_weight_fraction": [0.9, 0.8],
        "gamma_subspace_p": [0.9, 0.9001],
    })

    row = trend.select_closest_mode(df, target_frequency=194.93, source="case")

    assert int(row["mode_idx"]) == 1


def test_select_closest_mode_rejects_case_without_eligible_p_mode():
    df = pd.DataFrame({
        "mode_idx": [0, 1],
        "frequency": [194.9, 195.0],
        "gamma_k_weight_fraction": [0.9, 0.8],
        "gamma_subspace_p": [0.9, 0.2],
    })

    try:
        trend.select_closest_mode(df, target_frequency=194.93, source="case")
    except ValueError as exc:
        assert "gamma_subspace_p > 0.9" in str(exc)
    else:
        raise AssertionError("Expected no eligible p-mode error")


def test_select_closest_mode_breaks_equal_distance_tie_by_mode_idx():
    df = pd.DataFrame({
        "mode_idx": [4, 2],
        "frequency": [195.03, 194.83],
        "gamma_k_weight_fraction": [0.8, 0.7],
        "gamma_subspace_p": [0.995, 0.995],
    })

    row = trend.select_closest_mode(df, target_frequency=194.93, source="case")

    assert int(row["mode_idx"]) == 2


def test_build_trend_rows_reads_shift_directories():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        case_dir = out_root / "hexagon_cavity_0.01" / "strip_1d"
        case_dir.mkdir(parents=True)
        pd.DataFrame({
            "mode_idx": [0, 1],
            "frequency": [193.0, 194.5],
            "gamma_k_weight_fraction": [0.1, 0.6],
            "gamma_subspace_p": [0.995, 0.996],
        }).to_csv(case_dir / "mode_scores.csv", index=False)

        rows, issues = trend.build_trend_rows(
            out_root=out_root,
            shift_factors=[0.01, 0.02],
            target_frequency=194.93,
        )

    assert len(rows) == 1
    assert rows[0]["shift_factor"] == 0.01
    assert rows[0]["folder"] == "hexagon_cavity_0.01"
    assert rows[0]["mode_idx"] == 1
    assert rows[0]["frequency"] == 194.5
    assert abs(rows[0]["frequency_distance_to_target"] - 0.43) < 1e-12
    assert "frequency_delta_below_target" not in rows[0]
    assert rows[0]["gamma_k_weight_fraction"] == 0.6
    assert rows[0]["gamma_subspace_p"] == 0.996
    assert len(issues) == 1
    assert "hexagon_cavity_0.02" in issues[0]


def test_save_trend_plot_accepts_target_frequency_for_reference_line():
    with tempfile.TemporaryDirectory() as tmpdir:
        plot_path = Path(tmpdir) / "trend.png"
        selected_modes = pd.DataFrame({
            "shift_factor": [0.01, 0.02],
            "frequency": [194.4, 194.6],
            "gamma_k_weight_fraction": [0.5, 0.7],
        })

        trend.save_trend_plot(selected_modes, plot_path, target_frequency=194.93)

        assert plot_path.exists()
        assert plot_path.stat().st_size > 0


def test_save_trend_plot_anchors_legend_below_top_left():
    calls = []
    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        calls.append(kwargs.copy())
        return original_legend(self, *args, **kwargs)

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_path = Path(tmpdir) / "trend.png"
        selected_modes = pd.DataFrame({
            "shift_factor": [0.01, 0.02],
            "frequency": [194.4, 194.6],
            "gamma_k_weight_fraction": [0.5, 0.7],
        })

        matplotlib.axes.Axes.legend = capture_legend
        try:
            trend.save_trend_plot(selected_modes, plot_path, target_frequency=194.93)
        finally:
            matplotlib.axes.Axes.legend = original_legend

    assert calls
    assert calls[-1]["loc"] == "upper left"
    assert calls[-1]["bbox_to_anchor"] == (0.02, 0.82)


def test_output_paths_include_mesh_auto_size():
    assert trend.mesh_auto_size_label(9) == "mesh_auto_size_9"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = trend.analysis_dir_for_mesh_auto_size(9, out_root=Path(tmpdir))
        rows = [{
            "shift_factor": 0.01,
            "folder": "hexagon_cavity_0.01",
            "mode_idx": 1,
            "frequency": 194.5,
            "target_frequency": 194.93,
            "frequency_distance_to_target": 0.43,
            "gamma_k_weight_fraction": 0.6,
            "mode_scores_path": "mode_scores.csv",
        }]

        selected_csv, issues, plot_png = trend.write_outputs(
            rows=rows,
            issues=[],
            output_dir=output_dir,
            mesh_auto_size=9,
        )

        assert output_dir.name == "strip1d_trend"
        assert selected_csv.name == "strip1d_trend.csv"
        assert issues == []
        assert plot_png.name == "strip1d_trend.png"
        assert not list(output_dir.glob("*_issues.txt"))


def test_plot_from_selected_modes_csv_uses_csv_data_and_stem():
    calls = []

    def capture_plot(selected_modes, path, *, target_frequency):
        calls.append((selected_modes.copy(), path, target_frequency))
        path.write_text("plot", encoding="utf-8")

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "strip1d_trend"
        output_dir.mkdir()
        selected_csv = output_dir / "strip1d_trend.csv"
        pd.DataFrame({
            "shift_factor": [0.02, 0.01],
            "frequency": [194.6, 194.4],
            "target_frequency": [194.93, 194.93],
            "gamma_k_weight_fraction": [0.7, 0.5],
        }).to_csv(selected_csv, index=False)

        original_plot = trend.save_trend_plot
        trend.save_trend_plot = capture_plot
        try:
            plot_png = trend.plot_from_selected_modes_csv(selected_csv)
        finally:
            trend.save_trend_plot = original_plot

        assert plot_png.name == "strip1d_trend.png"
        assert plot_png.read_text(encoding="utf-8") == "plot"
        assert len(calls) == 1
        plotted, path, target_frequency = calls[0]
        assert path == plot_png
        assert target_frequency == 194.93
        assert list(plotted["shift_factor"]) == [0.01, 0.02]


def test_load_selected_modes_csv_filters_rows_to_shift_factors():
    with tempfile.TemporaryDirectory() as tmpdir:
        selected_csv = Path(tmpdir) / "selected_modes.csv"
        pd.DataFrame({
            "shift_factor": [0.00, 0.09, 0.10],
            "frequency": [194.2, 194.7, 194.8],
            "target_frequency": [194.93, 194.93, 194.93],
            "gamma_k_weight_fraction": [0.4, 0.9, 0.95],
        }).to_csv(selected_csv, index=False)

        selected_modes, target_frequency = trend.load_selected_modes_csv(selected_csv)

    assert target_frequency == 194.93
    assert list(selected_modes["shift_factor"]) == [0.00, 0.09]


def test_generate_trend_plot_reuses_existing_selected_modes_csv():
    calls = []

    def fail_build(*args, **kwargs):
        raise AssertionError("raw mode_scores should not be rebuilt when selected CSV exists")

    def capture_plot(selected_modes, path, *, target_frequency):
        calls.append((selected_modes.copy(), path, target_frequency))
        path.write_text("plot", encoding="utf-8")

    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        output_dir = trend.analysis_dir_for_mesh_auto_size(9, out_root=out_root)
        output_dir.mkdir(parents=True)
        selected_csv = output_dir / "strip1d_trend.csv"
        pd.DataFrame({
            "shift_factor": [0.02, 0.01],
            "frequency": [194.6, 194.4],
            "target_frequency": [194.93, 194.93],
            "gamma_k_weight_fraction": [0.7, 0.5],
        }).to_csv(selected_csv, index=False)

        original_build = trend.build_trend_rows
        original_plot = trend.save_trend_plot
        trend.build_trend_rows = fail_build
        trend.save_trend_plot = capture_plot
        try:
            selected_csv_out, issues, plot_png, used_existing_csv = trend.generate_trend_plot(
                out_root=out_root,
                output_dir=output_dir,
                mesh_auto_size=9,
                build_target_frequency=999.0,
            )
        finally:
            trend.build_trend_rows = original_build
            trend.save_trend_plot = original_plot

        assert used_existing_csv is True
        assert selected_csv_out == selected_csv
        assert issues == []
        assert plot_png.name == "strip1d_trend.png"
        assert plot_png.read_text(encoding="utf-8") == "plot"
        assert len(calls) == 1
        plotted, path, target_frequency = calls[0]
        assert path == plot_png
        assert target_frequency == 194.93
        assert list(plotted["shift_factor"]) == [0.01, 0.02]


def test_generate_trend_plot_creates_selected_csv_before_plotting_when_missing():
    calls = []

    def fake_build(*, out_root, shift_factors, target_frequency):
        assert target_frequency == 194.93
        return [
            {
                "shift_factor": 0.02,
                "folder": "hexagon_cavity_0.02",
                "mode_idx": 2,
                "frequency": 194.6,
                "target_frequency": target_frequency,
                "frequency_distance_to_target": 0.33,
                "gamma_k_weight_fraction": 0.7,
                "mode_scores_path": "mode_scores.csv",
            },
            {
                "shift_factor": 0.01,
                "folder": "hexagon_cavity_0.01",
                "mode_idx": 1,
                "frequency": 194.4,
                "target_frequency": target_frequency,
                "frequency_distance_to_target": 0.53,
                "gamma_k_weight_fraction": 0.5,
                "mode_scores_path": "mode_scores.csv",
            },
        ], ["hexagon_cavity_0.03: missing mode_scores.csv"]

    def capture_plot(selected_modes, path, *, target_frequency):
        calls.append((selected_modes.copy(), path, target_frequency))
        path.write_text("plot", encoding="utf-8")

    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        output_dir = trend.analysis_dir_for_mesh_auto_size(9, out_root=out_root)

        original_build = trend.build_trend_rows
        original_plot = trend.save_trend_plot
        trend.build_trend_rows = fake_build
        trend.save_trend_plot = capture_plot
        try:
            selected_csv, issues, plot_png, used_existing_csv = trend.generate_trend_plot(
                out_root=out_root,
                output_dir=output_dir,
                mesh_auto_size=9,
                build_target_frequency=194.93,
            )
        finally:
            trend.build_trend_rows = original_build
            trend.save_trend_plot = original_plot

        assert used_existing_csv is False
        assert selected_csv.exists()
        assert selected_csv.name == "strip1d_trend.csv"
        assert any("hexagon_cavity_0.03" in issue for issue in issues)
        assert not list(output_dir.glob("*_issues.txt"))
        assert plot_png.read_text(encoding="utf-8") == "plot"
        assert len(calls) == 1
        plotted, path, target_frequency = calls[0]
        assert path == plot_png
        assert target_frequency == 194.93
        assert list(plotted["shift_factor"]) == [0.01, 0.02]


if __name__ == "__main__":
    test_shift_factors_cover_zero_to_point_zero_nine()
    test_select_closest_mode_can_select_frequency_above_target()
    test_select_closest_mode_filters_gamma_subspace_before_frequency_distance()
    test_select_closest_mode_requires_gamma_subspace_strictly_above_threshold()
    test_select_closest_mode_rejects_case_without_eligible_p_mode()
    test_select_closest_mode_breaks_equal_distance_tie_by_mode_idx()
    test_build_trend_rows_reads_shift_directories()
    test_save_trend_plot_accepts_target_frequency_for_reference_line()
    test_save_trend_plot_anchors_legend_below_top_left()
    test_output_paths_include_mesh_auto_size()
    test_plot_from_selected_modes_csv_uses_csv_data_and_stem()
    test_load_selected_modes_csv_filters_rows_to_shift_factors()
    test_generate_trend_plot_reuses_existing_selected_modes_csv()
    test_generate_trend_plot_creates_selected_csv_before_plotting_when_missing()
    print("All strip gamma trend analysis tests passed.")
