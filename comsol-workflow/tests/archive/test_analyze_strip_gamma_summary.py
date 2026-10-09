import json
import os
import sys
import tempfile
import types
from pathlib import Path

import pandas as pd

from scripts.analysis import srip1d_trend as trend


def _summary_row(case_dir, shift_factor, *, cut_angle=60.0, mesh=5):
    return {
        "cladding_inward_shift_factor": shift_factor,
        "cladding_inward_shift": shift_factor * 0.82,
        "strip_cut_angle_degrees": cut_angle,
        "out_dir": str(case_dir),
        "bulk_radius": 10,
        "cladding_layers": 10,
        "mesh_auto_size": mesh,
        "eigenmode_count": 2,
        "case": "strip_1d",
    }


def test_build_trend_rows_reads_case_paths_from_run_summary():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        valid_case = out_root / "strip_1d_cav10_clad10_mesh5_cut60_0.01"
        missing_case = out_root / "strip_1d_cav10_clad10_mesh5_cut60_0.02"
        valid_case.mkdir()
        pd.DataFrame({
            "mode_idx": [0, 1],
            "frequency": [193.0, 194.5],
            "gamma_k_weight_fraction": [0.1, 0.6],
            "gamma_subspace_p": [0.995, 0.996],
        }).to_csv(valid_case / "mode_scores.csv", index=False)
        summary_path = out_root / "run_summary.json"
        summary_path.write_text(json.dumps({"runs": [
            _summary_row(valid_case.name, 0.01),
            _summary_row(missing_case.name, 0.02),
        ]}), encoding="utf-8")

        rows, issues, scan_metadata = trend.build_trend_rows_from_run_summary(
            summary_path,
            target_frequency=194.93,
        )

    assert len(rows) == 1
    assert rows[0]["shift_factor"] == 0.01
    assert rows[0]["folder"] == valid_case.name
    assert rows[0]["mode_idx"] == 1
    assert rows[0]["frequency"] == 194.5
    assert abs(rows[0]["frequency_distance_to_target"] - 0.43) < 1e-12
    assert rows[0]["gamma_subspace_p"] == 0.996
    assert "frequency_delta_below_target" not in rows[0]
    assert len(issues) == 1
    assert "missing" in issues[0]
    assert scan_metadata == {"mesh_auto_size": 5, "strip_cut_angle_degrees": 60.0}


def test_build_trend_rows_recovers_case_moved_beside_run_summary():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        case_name = "strip_1d_cav10_clad10_mesh5_cut60_0.09"
        actual_case = out_root / case_name
        actual_case.mkdir()
        pd.DataFrame({
            "mode_idx": [2, 3],
            "frequency": [195.9291336575822, 196.05676913004],
            "gamma_k_weight_fraction": [0.1199954982106016, 0.9765314757633616],
            "gamma_subspace_p": [0.0002762894715036, 0.9955824744603105],
        }).to_csv(actual_case / "mode_scores.csv", index=False)
        stale_case = out_root / "old-layout" / case_name
        summary_path = out_root / "run_summary.json"
        summary_path.write_text(
            json.dumps({"runs": [_summary_row(stale_case, 0.09)]}),
            encoding="utf-8",
        )

        rows, issues, _ = trend.build_trend_rows_from_run_summary(
            summary_path,
            target_frequency=195.902486611,
        )

    assert not issues
    assert len(rows) == 1
    assert rows[0]["mode_idx"] == 3
    assert Path(rows[0]["mode_scores_path"]).parent == actual_case


def test_mode_scores_path_prefers_layout_v2_and_falls_back_to_legacy(tmp_path):
    case_dir = tmp_path / "cut60_shift0.030"
    case_dir.mkdir()
    legacy = case_dir / "mode_scores.csv"
    legacy.write_text("mode_idx,frequency\n0,198\n", encoding="utf-8")

    assert trend.mode_scores_path_for_case(case_dir) == legacy

    current = case_dir / "10_overview" / "mode_scores.csv"
    current.parent.mkdir()
    current.write_text("mode_idx,frequency\n0,198\n", encoding="utf-8")

    assert trend.mode_scores_path_for_case(case_dir) == current


def test_scan_output_stem_is_stable():
    assert trend.scan_output_stem(5, 0) == "strip1d_trend"
    assert trend.scan_output_stem(5, 60) == "strip1d_trend"


def test_resolve_out_root_uses_current_series_by_default_and_preserves_override():
    current = Path("current-series")
    explicit = Path("explicit-series")

    assert trend.resolve_out_root(None, current_series_dir=current) == current
    assert trend.resolve_out_root(explicit, current_series_dir=current) == explicit


def test_main_uses_current_strip_output_directory_by_default():
    from scripts.run_main import run_strip_1d

    calls = []
    original_parse_args = trend.parse_args
    original_generate = trend.generate_trend_from_run_summary
    original_out_dir = run_strip_1d.OUT_DIR

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir) / "strip_1d"

        def fake_generate(*, run_summary_path, output_dir, target_frequency):
            calls.append((run_summary_path, output_dir, target_frequency))
            return output_dir / "trend.csv", [], output_dir / "trend.png"

        trend.parse_args = lambda: types.SimpleNamespace(
            out_root=None,
            output_dir=None,
            target_frequency=194.93,
            selected_modes_csv=None,
        )
        trend.generate_trend_from_run_summary = fake_generate
        run_strip_1d.OUT_DIR = out_dir
        try:
            result = trend.main()
        finally:
            trend.parse_args = original_parse_args
            trend.generate_trend_from_run_summary = original_generate
            run_strip_1d.OUT_DIR = original_out_dir

    assert result == 0
    assert calls == [(
        out_dir / "cut_60" / "run_summary.json",
        out_dir / "cut_60" / "strip1d_trend",
        194.93,
    )]


def test_generate_trend_from_summary_writes_strip1d_trend_artifacts():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        case_dir = out_root / "strip_1d_cav10_clad10_mesh5_cut60_0.01"
        case_dir.mkdir()
        pd.DataFrame({
            "mode_idx": [0, 1],
            "frequency": [193.0, 194.5],
            "gamma_k_weight_fraction": [0.1, 0.6],
            "gamma_subspace_p": [0.995, 0.996],
        }).to_csv(case_dir / "mode_scores.csv", index=False)
        summary_path = out_root / "run_summary.json"
        summary_path.write_text(
            json.dumps({"runs": [_summary_row(case_dir, 0.01)]}),
            encoding="utf-8",
        )

        output_dir = out_root / "strip1d_trend"
        selected_csv, issues, plot_png = trend.generate_trend_from_run_summary(
            run_summary_path=summary_path,
            output_dir=output_dir,
            target_frequency=194.93,
        )

        assert selected_csv == output_dir / "strip1d_trend.csv"
        assert plot_png == output_dir / "strip1d_trend.png"
        assert selected_csv.exists()
        assert issues == []
        assert plot_png.exists()
        assert not list(output_dir.glob("*_issues.txt"))


def test_generate_trend_from_summary_rejects_no_usable_cases():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        missing_case = out_root / "strip_1d_cav10_clad10_mesh5_cut60_0.01"
        summary_path = out_root / "run_summary.json"
        summary_path.write_text(
            json.dumps({"runs": [_summary_row(missing_case, 0.01)]}),
            encoding="utf-8",
        )

        try:
            trend.generate_trend_from_run_summary(
                run_summary_path=summary_path,
                output_dir=out_root,
                target_frequency=194.93,
            )
        except ValueError as exc:
            assert "No usable mode_scores.csv rows were found" in str(exc)
        else:
            raise AssertionError("Expected an error when no summary cases are usable")


def test_generate_trend_from_summary_supports_configured_shift_outside_legacy_grid():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_root = Path(tmpdir)
        case_dir = out_root / "strip_1d_cav10_clad10_mesh5_cut60_0.15"
        case_dir.mkdir()
        pd.DataFrame({
            "mode_idx": [0],
            "frequency": [194.5],
            "gamma_k_weight_fraction": [0.6],
            "gamma_subspace_p": [0.995],
        }).to_csv(case_dir / "mode_scores.csv", index=False)
        summary_path = out_root / "run_summary.json"
        summary_path.write_text(
            json.dumps({"runs": [_summary_row(case_dir, 0.15)]}),
            encoding="utf-8",
        )

        selected_csv, _, plot_png = trend.generate_trend_from_run_summary(
            run_summary_path=summary_path,
            output_dir=out_root,
            target_frequency=194.93,
        )

        selected = pd.read_csv(selected_csv)
        assert selected["shift_factor"].tolist() == [0.15]
        assert plot_png.exists()


if __name__ == "__main__":
    test_build_trend_rows_reads_case_paths_from_run_summary()
    test_build_trend_rows_recovers_case_moved_beside_run_summary()
    test_scan_output_stem_is_stable()
    test_resolve_out_root_uses_current_series_by_default_and_preserves_override()
    test_main_uses_current_strip_output_directory_by_default()
    test_generate_trend_from_summary_writes_strip1d_trend_artifacts()
    test_generate_trend_from_summary_rejects_no_usable_cases()
    test_generate_trend_from_summary_supports_configured_shift_outside_legacy_grid()
    print("All strip gamma trend summary tests passed.")
