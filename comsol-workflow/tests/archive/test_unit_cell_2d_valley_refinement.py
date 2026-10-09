from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from comsol_workflow.unit_cell_2d_valley_refinement import (
    SourceScan,
    build_initial_sampling_plan,
    build_valley_segments,
    detect_junction_centers,
    extract_valley_candidates,
    linear_interpolate_display,
    load_source_scan,
    merge_polarization_points,
)
from scripts.run_main.parameter_config import load_shared_parameters
from scripts.run_sweep.unit_cell_2d_valley_solver import (
    _combine_approved_repair_plan,
    build_convergence_repair_plan,
    build_convergence_repair_targets,
    audit_refinement_convergence,
    _edge_recenter_points,
    _next_round_points,
)


def _source_from_values(values: np.ndarray) -> SourceScan:
    axis = np.linspace(-0.1, 0.1, values.shape[0])
    records = []
    for iy, qy in enumerate(axis):
        for ix, qx in enumerate(axis):
            records.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "abs_cy_norm": values[iy, ix],
                    "cx_norm_re": 1.0,
                    "cx_norm_im": 0.0,
                    "cy_norm_re": values[iy, ix],
                    "cy_norm_im": 0.0,
                    "cx_raw_re": 1.0,
                    "cx_raw_im": 0.0,
                    "cy_raw_re": values[iy, ix],
                    "cy_raw_im": 0.0,
                }
            )
    return SourceScan(
        series_dir=Path("source"),
        frame=pd.DataFrame(records),
        q_axis=axis,
        coarse_step=float(axis[1] - axis[0]),
        q_max=0.1,
        source_table=Path("source.csv"),
        source_hashes={},
        config={},
    )


def test_parameter_snapshot_without_refinement_uses_reviewed_defaults(tmp_path):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    data["unit_cell_2d"].pop("valley_refinement")
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    settings = load_shared_parameters(path).unit_cell_2d.valley_refinement
    assert settings.normal_half_width_over_g == pytest.approx(0.003125)
    assert settings.normal_points_per_round == 3
    assert settings.anchor_stride == 3
    assert settings.tangent_half_width_over_g == pytest.approx(0.0125)
    assert settings.tangent_points_per_round == 5
    assert settings.longitudinal_max_gap_over_g == pytest.approx(0.00625)
    assert settings.refinement_rounds == 2
    assert settings.junction_points_per_axis == 7
    assert settings.maximum_new_points == 1000


def test_load_source_scan_requires_complete_planar_l2_grid(tmp_path):
    series = tmp_path / "series"
    (series / "99_config").mkdir(parents=True)
    (series / "80_logs").mkdir(parents=True)
    config = {"cache_identity": {"expected_eigenvalue_count": 4}}
    (series / "99_config" / "config.json").write_text(json.dumps(config), encoding="utf-8")
    rows = []
    for qy in (-0.1, 0.0, 0.1):
        for qx in (-0.1, 0.0, 0.1):
            rows.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "band_label": "p2",
                    "match_status": "matched",
                    "mode_idx": 1,
                    "cx_norm_re": 1.0,
                    "cx_norm_im": 0.0,
                    "cy_norm_re": qx * qx + qy * qy,
                    "cy_norm_im": 0.0,
                    "cx_raw_re": 1.0,
                    "cx_raw_im": 0.0,
                    "cy_raw_re": qx * qx + qy * qy,
                    "cy_raw_im": 0.0,
                    "normalization_kind": "planar_l2",
                    "field_plane": "air",
                    "polarization_definition_version": "unit-cell-2d-v1",
                }
            )
            point_dir = series / "01_results" / "k_points" / f"kx={qx:.6f}_ky={qy:.6f}"
            (point_dir / "eigenmodes").mkdir(parents=True, exist_ok=True)
            (point_dir / "point_metadata.json").write_text("{}", encoding="utf-8")
            (point_dir / "eigenmodes" / "01_Hz_center.parquet").touch()
    pd.DataFrame(rows).to_csv(series / "80_logs" / "polarization_grid.csv", index=False)
    source = load_source_scan(series)
    assert source.frame.shape[0] == 9
    assert source.coarse_step == pytest.approx(0.1)
    assert source.source_hashes


def test_candidates_are_strict_minima_and_position_fit_does_not_change_value():
    axis = np.linspace(-0.1, 0.1, 5)
    qx, qy = np.meshgrid(axis, axis)
    values = (qx - 0.2 * qy) ** 2 + 0.1
    source = _source_from_values(values)
    candidates = extract_valley_candidates(source)
    row = candidates.loc[candidates.axis == "row"]
    assert len(row) == 5
    assert np.allclose(row.abs_cy_raw_sample, values[:, 2])
    # The coordinate can move below the original grid spacing, but no fitted
    # amplitude column exists in the output contract.
    assert np.any(np.abs(row.qx_over_G - row.sample_qx_over_G) > 1e-6)
    assert "fitted_value" not in candidates.columns


def test_quarter_valley_detection_and_sampling_remain_canonical():
    axis = np.asarray([0.0, 0.05, 0.1])
    records = []
    for qy in axis:
        for qx in axis:
            value = qx**2 + 0.1
            records.append(
                {
                    "qx_over_G": qx,
                    "qy_over_G": qy,
                    "cx_norm_re": 1.0,
                    "cx_norm_im": 0.0,
                    "cy_norm_re": value,
                    "cy_norm_im": 0.0,
                    "cx_raw_re": 1.0,
                    "cx_raw_im": 0.0,
                    "cy_raw_re": value,
                    "cy_raw_im": 0.0,
                    "abs_cy_raw": value,
                }
            )
    source = SourceScan(
        series_dir=Path("source"),
        frame=pd.DataFrame(records),
        q_axis=axis,
        coarse_step=0.05,
        q_max=0.1,
        source_table=Path("source.csv"),
        source_hashes={},
        config={},
        is_quarter=1,
    )
    candidates = extract_valley_candidates(source)
    assert len(candidates) == 3
    assert np.allclose(candidates.qx_over_G, 0.0)
    assert (candidates[["qx_over_G", "qy_over_G"]] >= 0.0).all().all()
    segments = build_valley_segments(candidates, source.coarse_step)
    settings = SimpleNamespace(
        normal_half_width_over_g=0.01,
        normal_points_per_round=3,
        anchor_stride=1,
        tangent_half_width_over_g=0.02,
        tangent_points_per_round=3,
        longitudinal_max_gap_over_g=0.02,
        junction_half_width_over_g=0.02,
        junction_points_per_axis=3,
        maximum_new_points=100,
    )
    points = build_initial_sampling_plan(
        source, candidates, segments, pd.DataFrame(), settings
    )
    assert not points.empty
    assert (points[["qx_over_G", "qy_over_G"]] >= 0.0).all().all()
    assert not points[["kx_str", "ky_str"]].duplicated().any()


def test_slice_count_change_creates_data_driven_junction():
    source = _source_from_values(np.ones((5, 5)))
    candidates = pd.DataFrame(
        [
            {"candidate_id": 0, "axis": "row", "slice_index": 1, "qx_over_G": -0.03, "qy_over_G": -0.05},
            {"candidate_id": 1, "axis": "row", "slice_index": 1, "qx_over_G": 0.03, "qy_over_G": -0.05},
            {"candidate_id": 2, "axis": "row", "slice_index": 2, "qx_over_G": 0.0, "qy_over_G": 0.0},
            {"candidate_id": 3, "axis": "row", "slice_index": 3, "qx_over_G": -0.03, "qy_over_G": 0.05},
            {"candidate_id": 4, "axis": "row", "slice_index": 3, "qx_over_G": 0.03, "qy_over_G": 0.05},
        ]
    )
    segments = build_valley_segments(candidates, source.coarse_step)
    junctions = detect_junction_centers(candidates, segments, source.coarse_step)
    assert len(junctions) == 1
    assert abs(junctions.iloc[0].qx_over_G) < source.coarse_step
    assert abs(junctions.iloc[0].qy_over_G) < source.coarse_step


def test_initial_plan_deduplicates_coarse_and_six_decimal_points():
    axis = np.linspace(-0.1, 0.1, 5)
    qx, _qy = np.meshgrid(axis, axis)
    source = _source_from_values(qx**2)
    candidates = extract_valley_candidates(source)
    segments = build_valley_segments(candidates, source.coarse_step)
    settings = SimpleNamespace(
        normal_half_width_over_g=0.01,
        normal_points_per_round=3,
        anchor_stride=3,
        tangent_half_width_over_g=0.0375,
        tangent_points_per_round=5,
        longitudinal_max_gap_over_g=0.02,
        junction_half_width_over_g=0.02,
        junction_points_per_axis=5,
        maximum_new_points=100,
    )
    points = build_initial_sampling_plan(
        source, candidates, segments, pd.DataFrame(), settings
    )
    keys = set(zip(points.kx_str, points.ky_str))
    assert len(keys) == len(points)
    coarse = {
        (f"{qx:.6f}", f"{qy:.6f}")
        for qx, qy in zip(source.frame.qx_over_G, source.frame.qy_over_G)
    }
    assert keys.isdisjoint(coarse)
    assert {"normal", "tangent"}.issubset(set(points.stencil_axis))


def test_merge_rejects_conflict_and_uses_one_common_normalization():
    coarse = pd.DataFrame(
        [{"qx_over_G": 0.0, "qy_over_G": 0.0, "cx_raw_re": 3.0, "cx_raw_im": 0.0, "cy_raw_re": 4.0, "cy_raw_im": 0.0, "cx_norm_re": 30.0, "cx_norm_im": 0.0, "cy_norm_re": 40.0, "cy_norm_im": 0.0}]
    )
    refined = pd.DataFrame(
        [{"qx_over_G": 0.1, "qy_over_G": 0.0, "cx_raw_re": 0.0, "cx_raw_im": 0.0, "cy_raw_re": 1.0, "cy_raw_im": 0.0, "cx_norm_re": 0.0, "cx_norm_im": 0.0, "cy_norm_re": 100.0, "cy_norm_im": 0.0}]
    )
    merged = merge_polarization_points(coarse, refined)
    assert merged.common_c_magnitude_max.unique().tolist() == [5.0]
    assert merged.normalized_cy_magnitude.tolist() == pytest.approx([0.8, 0.2])
    conflict = coarse.copy()
    conflict["cy_raw_re"] = 9.0
    with pytest.raises(ValueError, match="Conflicting COMSOL values"):
        merge_polarization_points(coarse, conflict)


def test_merge_normalizes_mixed_k_labels_and_writes_parquet(tmp_path):
    coarse = pd.DataFrame(
        [
            {
                "qx_over_G": 0.000001,
                "qy_over_G": 0.0,
                "kx_str": 0.000001,
                "ky_str": 0.0,
                "cx_norm_re": 1.0,
                "cx_norm_im": 0.0,
                "cy_norm_re": 0.2,
                "cy_norm_im": 0.0,
                "cx_raw_re": 1.0,
                "cx_raw_im": 0.0,
                "cy_raw_re": 0.2,
                "cy_raw_im": 0.0,
            }
        ]
    )
    refined = pd.DataFrame(
        [
            {
                "qx_over_G": 0.000002,
                "qy_over_G": 0.0,
                "kx_str": "0.000002",
                "ky_str": "0.000000",
                "cx_norm_re": 1.0,
                "cx_norm_im": 0.0,
                "cy_norm_re": 0.1,
                "cy_norm_im": 0.0,
                "cx_raw_re": 1.0,
                "cx_raw_im": 0.0,
                "cy_raw_re": 0.1,
                "cy_raw_im": 0.0,
            }
        ]
    )
    merged = merge_polarization_points(coarse, refined)
    assert merged.kx_str.tolist() == ["0.000001", "0.000002"]
    assert merged.ky_str.tolist() == ["0.000000", "0.000000"]
    merged.to_parquet(tmp_path / "merged.parquet", index=False)


def test_plot_only_preserves_existing_dynamic_plan(tmp_path, monkeypatch):
    from scripts.analysis import unit_cell_2D_valley_refine as cli

    source_series = tmp_path / "source"
    source_series.mkdir()
    parameters = SimpleNamespace()
    source = SimpleNamespace(series_dir=source_series)
    plotted = []

    monkeypatch.setattr(cli, "load_shared_parameters", lambda _path: parameters)
    monkeypatch.setattr(cli, "load_source_scan", lambda _path: source)
    monkeypatch.setattr(cli, "refined_series_name", lambda _params, _source: "refined")
    monkeypatch.setattr(cli, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(
        cli,
        "prepare_valley_refinement",
        lambda *_args, **_kwargs: pytest.fail(
            "plot-only must not rebuild the sampling plan"
        ),
    )
    fake_solver = SimpleNamespace(
        plot_completed_refinement=lambda series, _params: (
            plotted.append(Path(series)),
            Path(series) / "10_overview" / "figure.png",
        )[1]
    )
    monkeypatch.setitem(
        sys.modules, "scripts.run_sweep.unit_cell_2d_valley_solver", fake_solver
    )

    assert cli.main([str(source_series), "--plot-only"]) == 0
    assert plotted == [(tmp_path / "refined").resolve()]


def test_run_approved_repair_preserves_plan_and_uses_explicit_solver(
    tmp_path, monkeypatch
):
    from scripts.analysis import unit_cell_2D_valley_refine as cli

    source_series = tmp_path / "source"
    source_series.mkdir()
    parameters = SimpleNamespace()
    source = SimpleNamespace(series_dir=source_series)
    calls = []
    monkeypatch.setattr(cli, "load_shared_parameters", lambda _path: parameters)
    monkeypatch.setattr(cli, "load_source_scan", lambda _path: source)
    monkeypatch.setattr(cli, "refined_series_name", lambda _params, _source: "refined")
    monkeypatch.setattr(cli, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(
        cli,
        "prepare_valley_refinement",
        lambda *_args, **_kwargs: pytest.fail(
            "approved repair must not rebuild the sampling plan"
        ),
    )
    fake_solver = SimpleNamespace(
        solve_approved_convergence_repair=lambda series, _params, **kwargs: (
            calls.append((Path(series), kwargs)),
            {"status": "complete", "solve_count": 2},
        )[1],
        plot_completed_refinement=lambda series, _params: Path(series)
        / "10_overview"
        / "figure.png",
    )
    monkeypatch.setitem(
        sys.modules, "scripts.run_sweep.unit_cell_2d_valley_solver", fake_solver
    )

    assert cli.main([str(source_series), "--run-approved-repair"]) == 0
    assert calls == [
        (
            (tmp_path / "refined").resolve(),
            {"resume": True, "allow_concurrent_comsol": False},
        )
    ]


def test_piecewise_linear_interpolation_is_exact_at_nodes_and_no_overshoot():
    points = pd.DataFrame(
        {
            "qx_over_G": [0.0, 1.0, 0.0, 1.0],
            "qy_over_G": [0.0, 0.0, 1.0, 1.0],
            "normalized_cy_magnitude": [0.0, 1.0, 1.0, 0.0],
        }
    )
    qx, qy, values = linear_interpolate_display(points, np.asarray([0.0, 0.5, 1.0]))
    assert values.min() >= 0.0
    assert values.max() <= 1.0
    for x, y, expected in zip(points.qx_over_G, points.qy_over_G, points.normalized_cy_magnitude):
        ix = int(np.where(np.isclose(qx[0], x))[0][0])
        iy = int(np.where(np.isclose(qy[:, 0], y))[0][0])
        assert float(values[iy, ix]) == pytest.approx(expected)


def _adaptive_settings():
    return SimpleNamespace(
        refinement_ratio=4.0,
        normal_half_width_over_g=0.01,
        normal_points_per_round=3,
        tangent_half_width_over_g=0.02,
        tangent_points_per_round=5,
        junction_half_width_over_g=0.02,
        junction_points_per_axis=5,
    )


def test_next_round_center_uses_normal_samples_not_lower_tangent_endpoint():
    completed = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "stencil_axis": "normal", "qx_over_G": -0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "stencil_axis": "normal", "qx_over_G": 0.0, "qy_over_G": 0.0, "cy_raw_re": 0.1, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "stencil_axis": "normal", "qx_over_G": 0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "stencil_axis": "tangent", "qx_over_G": 0.0, "qy_over_G": 0.02, "cy_raw_re": 0.001, "cy_raw_im": 0.0, "match_status": "matched"},
        ]
    )
    previous = pd.DataFrame(
        [{"plan_source": "valley_cross_stencil", "source_id": 1, "qx_over_G": -0.01, "qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0}]
    )
    source = SimpleNamespace(
        frame=pd.DataFrame({"qx_over_G": [], "qy_over_G": []}), q_max=0.1
    )
    points = _next_round_points(
        completed, previous, source, _adaptive_settings(), round_number=2
    )
    assert not points.empty
    assert np.allclose(points.center_qx_over_G, 0.0)
    assert np.allclose(points.center_qy_over_G, 0.0)


def test_edge_recenter_reuses_deduplicated_center_value():
    completed = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 1, "qx_over_G": 0.0, "qy_over_G": 0.0, "cy_raw_re": 0.1, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 2, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 0, "qx_over_G": -0.0025, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 2, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 2, "qx_over_G": 0.0025, "qy_over_G": 0.0, "cy_raw_re": 0.3, "cy_raw_im": 0.0, "match_status": "matched"},
        ]
    )
    previous = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 2, "recenter_index": 0, "qx_over_G": -0.0025, "qy_over_G": 0.0, "center_qx_over_G": 0.0, "center_qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0}
        ]
    )
    source = SimpleNamespace(
        frame=pd.DataFrame(
            columns=["qx_over_G", "qy_over_G", "cy_raw_re", "cy_raw_im"]
        ),
        q_max=0.1,
    )
    points = _edge_recenter_points(
        completed,
        previous,
        source,
        _adaptive_settings(),
        round_number=2,
        recenter_index=1,
    )
    assert points.empty


def test_convergence_repair_is_unique_and_respects_total_point_limit():
    completed = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 0, "qx_over_G": -0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 1, "qx_over_G": 0.0, "qy_over_G": 0.0, "cy_raw_re": 0.1, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 2, "qx_over_G": 0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "tangent", "axis_index": 0, "qx_over_G": 0.0, "qy_over_G": -0.02, "cy_raw_re": 0.001, "cy_raw_im": 0.0, "match_status": "matched"},
        ]
    )
    previous = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "qx_over_G": -0.01, "qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "qx_over_G": 0.0, "qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "qx_over_G": 0.01, "qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0},
        ]
    )
    source = SimpleNamespace(
        frame=pd.DataFrame({"qx_over_G": [0.0], "qy_over_G": [0.0]}),
        q_max=0.1,
    )
    settings = _adaptive_settings()
    settings.maximum_new_points = 100
    repair = build_convergence_repair_plan(completed, previous, source, settings)
    targets = build_convergence_repair_targets(completed, previous, settings)
    assert len(targets) == 1
    repair_keys = set(zip(repair.kx_str, repair.ky_str))
    previous_keys = {
        (f"{qx:.6f}", f"{qy:.6f}")
        for qx, qy in zip(previous.qx_over_G, previous.qy_over_G)
    }
    assert len(repair_keys) == len(repair)
    assert repair_keys.isdisjoint(previous_keys)
    assert len(previous) + len(repair) <= settings.maximum_new_points

    settings.maximum_new_points = len(previous)
    with pytest.raises(RuntimeError, match="maximum_new_points"):
        build_convergence_repair_plan(completed, previous, source, settings)


def test_repair_target_survives_when_all_requested_coordinates_already_exist():
    completed = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 0, "qx_over_G": -0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 1, "qx_over_G": 0.0, "qy_over_G": 0.0, "cy_raw_re": 0.1, "cy_raw_im": 0.0, "match_status": "matched"},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_axis": "normal", "axis_index": 2, "qx_over_G": 0.01, "qy_over_G": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "match_status": "matched"},
        ]
    )
    existing_coordinates = {
        (-0.01, 0.0),
        (-0.0025, 0.0),
        (0.0025, 0.0),
        (0.0, -0.005),
        (0.0, -0.0025),
        (0.0, 0.0025),
        (0.0, 0.005),
    }
    previous = pd.DataFrame(
        [
            {
                "plan_source": "valley_cross_stencil",
                "source_id": 1,
                "round": 1,
                "recenter_index": 0,
                "qx_over_G": qx,
                "qy_over_G": qy,
                "normal_x": 1.0,
                "normal_y": 0.0,
                "tangent_x": 0.0,
                "tangent_y": 1.0,
            }
            for qx, qy in existing_coordinates
        ]
    )
    source = SimpleNamespace(
        frame=pd.DataFrame({"qx_over_G": [0.0], "qy_over_G": [0.0]}),
        q_max=0.1,
    )
    settings = _adaptive_settings()
    settings.maximum_new_points = 100
    repair = build_convergence_repair_plan(completed, previous, source, settings)
    targets = build_convergence_repair_targets(completed, previous, settings)
    assert repair.empty
    assert len(targets) == 1


def test_approved_repair_append_is_unique_and_resume_safe():
    previous = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 1, "recenter_index": 0, "stencil_index": 0, "stencil_axis": "normal", "axis_index": 0, "qx_over_G": -0.01, "qy_over_G": 0.0, "center_qx_over_G": 0.0, "center_qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0, "point_id": 0, "solve_order": 0},
        ]
    )
    repair = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "round": 2, "recenter_index": 1, "stencil_index": 0, "stencil_axis": "normal", "axis_index": 0, "qx_over_G": -0.0025, "qy_over_G": 0.0, "center_qx_over_G": 0.0, "center_qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0, "tangent_x": 0.0, "tangent_y": 1.0, "point_id": 1, "solve_order": 1},
        ]
    )
    source = SimpleNamespace(
        frame=pd.DataFrame({"qx_over_G": [0.0], "qy_over_G": [0.0]}),
        q_max=0.1,
    )
    settings = _adaptive_settings()
    settings.maximum_new_points = 10
    combined, new = _combine_approved_repair_plan(
        previous, repair, source, settings
    )
    assert len(combined) == 2
    assert len(new) == 1
    resumed, resumed_new = _combine_approved_repair_plan(
        combined, repair, source, settings
    )
    assert len(resumed) == 2
    assert resumed_new.empty


def test_convergence_audit_requires_an_interior_normal_minimum():
    axis = np.linspace(-0.1, 0.1, 5)
    qx, qy = np.meshgrid(axis, axis)
    source = _source_from_values(qx**2 + qy**2 + 0.1)
    repair = pd.DataFrame(
        [
            {"plan_source": "valley_cross_stencil", "source_id": 1, "repair_kind": "round2_normal_refine", "qx_over_G": -0.0025, "qy_over_G": 0.0, "center_qx_over_G": 0.0, "center_qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0},
            {"plan_source": "valley_cross_stencil", "source_id": 1, "repair_kind": "round2_normal_refine", "qx_over_G": 0.0025, "qy_over_G": 0.0, "center_qx_over_G": 0.0, "center_qy_over_G": 0.0, "normal_x": 1.0, "normal_y": 0.0},
        ]
    )

    def completed(left: float) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"qx_over_G": -0.0025, "qy_over_G": 0.0, "match_status": "matched", "cx_raw_re": 1.0, "cx_raw_im": 0.0, "cy_raw_re": left, "cy_raw_im": 0.0, "cx_norm_re": 1.0, "cx_norm_im": 0.0, "cy_norm_re": left, "cy_norm_im": 0.0},
                {"qx_over_G": 0.0025, "qy_over_G": 0.0, "match_status": "matched", "cx_raw_re": 1.0, "cx_raw_im": 0.0, "cy_raw_re": 0.2, "cy_raw_im": 0.0, "cx_norm_re": 1.0, "cx_norm_im": 0.0, "cy_norm_re": 0.2, "cy_norm_im": 0.0},
            ]
        )

    settings = _adaptive_settings()
    settings.maximum_new_points = 100
    converged = audit_refinement_convergence(
        source, completed(0.2), repair, settings
    )
    assert converged["converged"] is True
    assert converged["edge_minimum_count"] == 0
    edge = audit_refinement_convergence(source, completed(0.01), repair, settings)
    assert edge["converged"] is False
    assert edge["edge_minimum_count"] == 1
