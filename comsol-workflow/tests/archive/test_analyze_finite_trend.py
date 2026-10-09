import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.analysis import finite_trend as trend


def _scores(rows):
    return pd.DataFrame(rows, columns=[
        "mode_idx",
        "frequency",
        "q",
        "gamma_p_px_weight_fraction",
        "gamma_k_weight_fraction",
        "gamma_subspace_p",
    ])


def _write_case(case_dir: Path, rows) -> Path:
    overview = case_dir / "10_overview"
    overview.mkdir(parents=True)
    path = overview / "mode_scores.csv"
    _scores(rows).to_csv(path, index=False)
    return path


def _write_summary(series_dir: Path, payload: dict) -> Path:
    path = series_dir / "run_summary.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_select_closest_mode_uses_strict_py_weight_threshold():
    scores = _scores([
        [0, 198.4, 1000.0, 0.80, 0.99, 0.999],
        [1, 198.6, 2000.0, 0.8001, 0.10, 0.10],
    ])

    selected = trend.select_closest_mode(
        scores,
        target_frequency=198.4,
        source="shift0.010",
    )

    assert int(selected["mode_idx"]) == 1
    assert float(selected["gamma_p_px_weight_fraction"]) > 0.80


def test_select_closest_mode_ignores_other_gamma_metrics_for_threshold():
    scores = _scores([
        [0, 198.4, 1000.0, 0.79, 0.999, 0.999],
        [1, 198.7, 1200.0, 0.90, 0.01, 0.01],
    ])

    selected = trend.select_closest_mode(
        scores,
        target_frequency=198.4,
        source="shift0.010",
    )

    assert int(selected["mode_idx"]) == 1


def test_build_full_finite_rows_reads_every_shift_present_in_series(tmp_path):
    series = tmp_path / "finite_series"
    listed = series / "shift0.010"
    unlisted = series / "shift0.030"
    score_path = _write_case(listed, [
        [0, 198.41, 100.0, 0.90, 0.2, 0.5],
        [1, 198.60, 9999.0, 0.95, 0.9, 0.9],
    ])
    _write_case(unlisted, [[2, 198.40, 5000.0, 0.99, 0.9, 0.9]])
    summary = _write_summary(series, {
        "finite_cavity_preset_config": {
            "eigenfrequency_shift": "198.424 [THz]",
        },
        "runs": [{
            "case": "finite_cavity",
            "cladding_inward_shift_factor": 0.01,
            "out_dir": str(listed),
            "overview_dir": str(listed / "10_overview"),
        }],
    })

    rows, issues, metadata = trend.build_trend_rows_from_run_summary(
        summary,
        target_frequency=198.424,
    )

    assert issues == []
    assert metadata == {"source_kind": "finite", "symmetry_id": None}
    assert len(rows) == 2
    assert [row["shift_factor"] for row in rows] == [0.01, 0.03]
    assert rows[0]["mode_idx"] == 0
    assert rows[0]["q"] == 100.0
    assert Path(rows[0]["mode_scores_path"]) == score_path
    assert rows[1]["mode_idx"] == 2


def test_build_quarter_rows_recovers_shift_and_stale_absolute_case_path(tmp_path):
    series = tmp_path / "quarter_series"
    actual_case = series / "shift0.030"
    _write_case(actual_case, [[1, 198.3, 777.0, 0.91, 0.2, 0.3]])
    stale_case = tmp_path / "old_machine" / "shift0.030"
    summary = _write_summary(series, {
        "symmetry_ids": [1],
        "runs": [{
            "symmetry_id": 1,
            "out_dir": str(stale_case),
            "overview_dir": str(stale_case / "10_overview"),
        }],
    })

    rows, issues, metadata = trend.build_trend_rows_from_run_summary(
        summary,
        target_frequency=198.4,
    )

    assert issues == []
    assert metadata == {"source_kind": "finite_quarter", "symmetry_id": 1}
    assert rows[0]["shift_factor"] == 0.03
    assert Path(rows[0]["mode_scores_path"]).parent == actual_case / "10_overview"


def test_build_directional_rows_distinguishes_same_x_with_different_y(tmp_path):
    series = tmp_path / "quarter_series"
    listed = series / "shiftx0.050_shifty0.100"
    supplemental = series / "shiftx0.050_shifty0.080"
    _write_case(listed, [[1, 198.3, 777.0, 0.91, 0.2, 0.3]])
    _write_case(supplemental, [[2, 198.4, 888.0, 0.92, 0.2, 0.3]])
    summary = _write_summary(series, {
        "symmetry_ids": [1],
        "runs": [{
            "symmetry_id": 1,
            "cladding_shift_base_factor": 0.05,
            "cladding_x_shift_factor": 0.05,
            "cladding_y_shift_factor": 0.10,
            "out_dir": str(listed),
        }],
    })

    rows, issues, metadata = trend.build_trend_rows_from_run_summary(
        summary,
        target_frequency=198.4,
    )

    assert issues == []
    assert metadata == {"source_kind": "finite_quarter", "symmetry_id": 1}
    assert len(rows) == 2
    assert {
        (row["cladding_x_shift_factor"], row["cladding_y_shift_factor"])
        for row in rows
    } == {(0.05, 0.08), (0.05, 0.10)}
    assert {row["shift_factor"] for row in rows} == {0.05}
    assert {row["cladding_y_over_x_shift_ratio"] for row in rows} == {1.6, 2.0}

    ratio_rows, ratio_issues, ratio_metadata = (
        trend.build_trend_rows_from_run_summary(
            summary,
            target_frequency=198.4,
            y_over_x_ratio=2.0,
        )
    )
    assert ratio_issues == []
    assert len(ratio_rows) == 1
    assert ratio_rows[0]["cladding_y_shift_factor"] == 0.10
    assert ratio_metadata["cladding_y_over_x_shift_ratio"] == 2.0


def test_directional_ratio_filter_keeps_zero_shift_baseline():
    assert trend._shift_pair_matches_ratio((0.0, 0.0), 2.0)
    assert trend._shift_pair_matches_ratio((0.05, 0.10), 2.0)
    assert not trend._shift_pair_matches_ratio((0.05, 0.08), 2.0)
    assert not trend._shift_pair_matches_ratio(None, 2.0)


def test_directional_plot_groups_share_origin_without_joining_ratios():
    selected = pd.DataFrame({
        "cladding_x_shift_factor": [0.00, 0.01, 0.05, 0.05],
        "cladding_y_shift_factor": [0.00, 0.02, 0.08, 0.10],
        "frequency": [198.0, 198.1, 198.2, 198.3],
        "q": [100.0, 200.0, 300.0, 400.0],
        "gamma_p_px_weight_fraction": [0.9, 0.9, 0.9, 0.9],
    })

    groups = trend._directional_plot_groups(selected)

    assert groups is not None
    assert [label for label, _group in groups] == ["y/x=1.6", "y/x=2"]
    assert [len(group) for _label, group in groups] == [2, 3]
    assert all(group.iloc[0]["cladding_x_shift_factor"] == 0.0 for _, group in groups)


def test_selected_mode_with_invalid_q_is_not_replaced_by_farther_mode(tmp_path):
    series = tmp_path / "finite_series"
    case = series / "shift0.010"
    _write_case(case, [
        [0, 198.40, float("nan"), 0.95, 0.9, 0.9],
        [1, 198.80, 5000.0, 0.99, 0.9, 0.9],
    ])
    summary = _write_summary(series, {
        "runs": [{
            "case": "finite_cavity",
            "cladding_inward_shift_factor": 0.01,
            "out_dir": str(case),
        }],
    })

    rows, issues, _ = trend.build_trend_rows_from_run_summary(
        summary,
        target_frequency=198.40,
    )

    assert rows == []
    assert len(issues) == 1
    assert "non-finite q" in issues[0]


def test_multiple_quarter_symmetries_require_explicit_selection(tmp_path):
    series = tmp_path / "quarter_series"
    case1 = series / "shift0.010" / "symmetry_1"
    case2 = series / "shift0.010" / "symmetry_2"
    _write_case(case1, [[0, 198.4, 100.0, 0.9, 0.2, 0.3]])
    _write_case(case2, [[0, 198.5, 200.0, 0.9, 0.2, 0.3]])
    summary = _write_summary(series, {
        "symmetry_ids": [1, 2],
        "runs": [
            {"symmetry_id": 1, "out_dir": str(case1)},
            {"symmetry_id": 2, "out_dir": str(case2)},
        ],
    })

    with pytest.raises(ValueError, match="multiple symmetry IDs"):
        trend.build_trend_rows_from_run_summary(
            summary,
            target_frequency=198.4,
        )

    rows, issues, metadata = trend.build_trend_rows_from_run_summary(
        summary,
        target_frequency=198.4,
        symmetry_id=2,
    )
    assert issues == []
    assert metadata["symmetry_id"] == 2
    assert rows[0]["symmetry_id"] == 2
    assert rows[0]["frequency"] == 198.5


def test_resolve_target_frequency_prefers_explicit_then_summary_then_fallback():
    payload = {
        "finite_cavity_preset_config": {
            "eigenfrequency_shift": "198.424 [THz]",
        }
    }

    assert trend.resolve_target_frequency(
        payload,
        target_frequency=199.0,
        fallback_target_frequency="197 [THz]",
    ) == 199.0
    assert trend.resolve_target_frequency(
        payload,
        target_frequency=None,
        fallback_target_frequency="197 [THz]",
    ) == 198.424
    assert trend.resolve_target_frequency(
        {},
        target_frequency=None,
        fallback_target_frequency="197 [THz]",
    ) == 197.0


def test_generate_trend_writes_csv_and_two_plots_for_same_modes(tmp_path):
    series = tmp_path / "quarter_series"
    for shift, mode_idx, frequency, q_value, weight in (
        (0.01, 1, 198.2, 100.0, 0.81),
        (0.03, 2, 198.3, 300.0, 0.95),
    ):
        _write_case(
            series / f"shift{shift:.3f}",
            [[mode_idx, frequency, q_value, weight, 0.2, 0.3]],
        )
    summary = _write_summary(series, {
        "finite_cavity_preset_config": {
            "eigenfrequency_shift": "198.25 [THz]",
        },
        "symmetry_ids": [1],
        "runs": [
            {"symmetry_id": 1, "out_dir": str(series / "shift0.010")},
            {"symmetry_id": 1, "out_dir": str(series / "shift0.030")},
        ],
    })
    output_dir = series / "finite_trend"

    selected_csv, issues, finite_png, q_png, target = (
        trend.generate_trend_from_run_summary(
            run_summary_path=summary,
            output_dir=output_dir,
        )
    )

    assert issues == []
    assert target == 198.25
    assert selected_csv == output_dir / "finite_trend.csv"
    assert finite_png == output_dir / "finite_trend.png"
    assert q_png == output_dir / "finite_q_trend.png"
    assert selected_csv.is_file()
    assert finite_png.stat().st_size > 0
    assert q_png.stat().st_size > 0
    selected = pd.read_csv(selected_csv)
    assert selected["mode_idx"].tolist() == [1, 2]
    assert selected["q"].tolist() == [100.0, 300.0]


def test_plot_from_existing_csv_rejects_weight_at_threshold(tmp_path):
    selected_csv = tmp_path / "finite_trend.csv"
    pd.DataFrame({
        "shift_factor": [0.01],
        "frequency": [198.2],
        "target_frequency": [198.25],
        "gamma_p_px_weight_fraction": [0.80],
        "q": [100.0],
        "source_kind": ["finite"],
    }).to_csv(selected_csv, index=False)

    with pytest.raises(ValueError, match="not > 0.8"):
        trend.plot_from_selected_modes_csv(selected_csv)
