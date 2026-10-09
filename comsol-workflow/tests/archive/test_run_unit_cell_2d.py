import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from comsol_workflow.simulation_spatial.hexagon_unit_cell import (
    SimulationRun,
    _split_numerical_data,
)
from scripts.run_main.parameter_config import load_shared_parameters
from scripts.run_main.run_unit_cell_2d import (
    _cache_identity,
    _match_passes_ambiguity,
    _nearest_predecessors,
    _save_runner_model,
    build_parser,
    make_grid_points,
    resolved_series_dir,
    select_point_modes,
)


def test_default_axis_builds_625_points_with_gamma_first():
    parameters = load_shared_parameters()
    settings = parameters.unit_cell_2d
    points = make_grid_points(
        settings.q_axis_over_g, sampling_domain=settings.sampling_domain
    )
    assert settings.q_max_over_g == pytest.approx(0.15)
    assert settings.q_points_per_axis == 25
    assert settings.is_quarter == 1
    assert settings.q_axis_over_g == pytest.approx(np.linspace(0.0, 0.15, 25))
    assert settings.plot_q_axis_over_g == pytest.approx(
        np.linspace(-0.15, 0.15, 49)
    )
    assert np.diff(settings.q_axis_over_g) == pytest.approx(
        np.full(24, 0.15 / 24)
    )
    assert len(points) == 625
    assert points[0]["qx_over_G"] == 0.0
    assert points[0]["qy_over_G"] == 0.0
    assert sum(point["k_norm"] == 0.0 for point in points) == 1
    assert all(point["sampling_domain"] == "quarter" for point in points)
    assert all(point["qx_over_G"] >= 0.0 for point in points)
    assert all(point["qy_over_G"] >= 0.0 for point in points)
    assert [point["solve_order"] for point in points] == list(range(625))


def test_series_and_cache_identity_use_unit_cell_2d_geometry():
    parameters = load_shared_parameters()
    geometry = parameters.unit_cell_2d.cell
    name = resolved_series_dir(parameters).name
    assert f"({geometry.b0_nm:g}-" in name
    assert "clad" not in name.lower()
    assert "_uniform25_" in name
    assert "_quarter_band-py" in name
    assert "multiscale" not in name
    identity = _cache_identity(parameters)
    assert parameters.unit_cell_eigenmode_count == 10
    assert parameters.unit_cell_2d.eigenmode_pair_count == 2
    assert identity["eigenmode_pair_count"] == 2
    assert identity["expected_eigenvalue_count"] == 4
    assert "unit_cell_eigenmode_count" not in identity
    assert identity["unit_cell_2d"]["b0_nm"] == geometry.b0_nm
    assert identity["isQuarter"] == 1
    assert identity["solved_point_count"] == 625
    assert identity["plot_coordinate_count"] == 2401
    assert "cavity" not in identity
    assert "cladding" not in identity


def test_nearest_predecessors_are_deterministic_and_can_be_outward():
    point = {"kx": 0.01, "ky": 0.0}
    states = {
        (0.02, 0.0): {"name": "outward"},
        (0.0, 0.0): {"name": "gamma"},
        (0.01, 0.02): {"name": "vertical"},
    }
    selected = _nearest_predecessors(point, states, limit=3)
    assert [item["name"] for item in selected] == ["gamma", "outward", "vertical"]


def test_ambiguity_gate_is_inclusive():
    settings = SimpleNamespace(minimum_ambiguity_gap=0.05)
    assert _match_passes_ambiguity({"ambiguity_gap": 0.05}, settings)
    assert not _match_passes_ambiguity({"ambiguity_gap": 0.049}, settings)


def test_tracking_averages_multiple_predecessors_and_ignores_out_of_gate_runner_up(
    monkeypatch, tmp_path
):
    import scripts.run_main.run_unit_cell_2d as workflow

    modes = pd.DataFrame(
        {"re": [200.0, 250.0], "q": [100.0, 200.0], "is_valid": [True, True]}
    )
    states = {
        "p2": {
            (0.0, 0.0): {"frequency_thz": 199.8, "field": "pred-a", "qx_over_G": 0.0, "qy_over_G": 0.0},
            (0.01, 0.0): {"frequency_thz": 200.2, "field": "pred-b", "qx_over_G": 0.01, "qy_over_G": 0.0},
        }
    }
    monkeypatch.setattr(workflow, "load_hz_center", lambda _path, idx: f"candidate-{idx}")
    monkeypatch.setattr(
        workflow,
        "field_overlap",
        lambda _previous, candidate: 0.8 if candidate == "candidate-0" else 0.99,
    )
    settings = SimpleNamespace(
        frequency_tolerance_thz=7.0,
        minimum_overlap=0.2,
        minimum_ambiguity_gap=0.05,
    )
    selected = select_point_modes(
        {"kx": 0.02, "ky": 0.0, "k_norm": 0.02},
        modes,
        tmp_path,
        ["p2"],
        states,
        settings,
    )["p2"]
    assert selected["mode_idx"] == 0
    assert selected["predecessors"] == [[0.01, 0.0], [0.0, 0.0]]
    assert selected["selected_score"] == pytest.approx(0.8)
    assert selected["runner_up_score"] == 0.0


def test_analysis_is_on_by_default_and_only_cli_flag_skips_it():
    assert build_parser().parse_args([]).skip_analysis is False
    assert build_parser().parse_args(["--skip-analysis"]).skip_analysis is True
    source = Path("scripts/run_main/run_unit_cell_2d.py").read_text(encoding="utf-8")
    assert "if not skip_analysis:" in source
    assert "parameters.unit_cell_eigenmode_count" not in source
    assert "parameters.cavity" not in source
    assert "parameters.cladding" not in source


def test_runner_model_is_saved_to_00_model(tmp_path):
    saved = []

    class FakeModel:
        def save(self, path):
            saved.append(path)

    target = tmp_path / "00_model" / "unit_cell_2D.mph"
    _save_runner_model(SimpleNamespace(model=FakeModel()), target)
    assert target.parent.is_dir()
    assert saved == [str(target)]


@pytest.mark.parametrize(
    ("raw", "count", "expected"),
    [
        (np.arange(12).reshape(3, 4), 3, [np.arange(4), np.arange(4, 8), np.arange(8, 12)]),
        (np.arange(12).reshape(4, 3), 3, [np.array([0, 3, 6, 9]), np.array([1, 4, 7, 10]), np.array([2, 5, 8, 11])]),
    ],
)
def test_split_numerical_data_supports_expression_axes(raw, count, expected):
    actual = _split_numerical_data(raw, count)
    assert all(np.array_equal(left, right) for left, right in zip(actual, expected))


def test_compute_polarization_data_reconstructs_complex_fields_and_normalizes():
    class FakeRun:
        def get_2d_fields_batch(self, eigenmode_idx, expressions, plane):
            assert eigenmode_idx == 2
            assert plane == "air"
            points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
            values = {
                expressions[0]: np.full(4, 2.0),
                expressions[1]: np.full(4, 3.0),
                expressions[2]: np.full(4, 4.0),
                expressions[3]: np.full(4, 5.0),
            }
            return points, values

    result = SimulationRun.compute_polarization_data(FakeRun(), 2, "air")
    assert result["cx_raw"] == pytest.approx(2.0 + 3.0j)
    assert result["cy_raw"] == pytest.approx(4.0 + 5.0j)
    assert result["normalization_value"] == pytest.approx(np.sqrt(54.0))
    assert result["cx_norm"] == pytest.approx((2.0 + 3.0j) / np.sqrt(54.0))
    assert result["cy_norm"] == pytest.approx((4.0 + 5.0j) / np.sqrt(54.0))


def test_parameter_file_defines_unit_cell_2d_geometry_and_enables_analysis():
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    assert data["unit_cell_band"]["eigenmode_count"] == load_shared_parameters().unit_cell_eigenmode_count
    assert data["unit_cell_2d"]["eigenmode_pair_count"] == 2
    assert data["unit_cell_2d"]["analysis_enabled"] is True
    assert {"b0_nm", "eta", "zeta"} <= data["unit_cell_2d"].keys()
    assert data["unit_cell_2d"]["q_max_over_G"] == 0.15
    assert data["unit_cell_2d"]["q_points_per_axis"] == 25
    assert data["unit_cell_2d"]["isQuarter"] == 1
    assert "q_axis_over_G" not in data["unit_cell_2d"]


def test_unit_cell_2d_geometry_is_independent_from_top_level_cavity(tmp_path):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    data["cells"]["cavity"] = {"b0_nm": 999.0, "eta": 0.5, "zeta": 0.75}
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    parameters = load_shared_parameters(path)
    assert parameters.unit_cell_2d.cell.to_dict() == {
        key: data["unit_cell_2d"][key] for key in ("b0_nm", "eta", "zeta")
    }


def test_missing_is_quarter_preserves_legacy_full_scan(tmp_path):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    data["unit_cell_2d"].pop("isQuarter")
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    settings = load_shared_parameters(path).unit_cell_2d
    assert settings.is_quarter == 0
    assert settings.sampling_domain == "full"
    assert settings.q_axis_over_g == pytest.approx(np.linspace(-0.15, 0.15, 25))
    assert settings.plot_q_axis_over_g == pytest.approx(settings.q_axis_over_g)


def test_quarter_allows_even_axis_count(tmp_path):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    data["unit_cell_2d"]["q_points_per_axis"] = 18
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    settings = load_shared_parameters(path).unit_cell_2d
    assert len(settings.q_axis_over_g) == 18
    assert len(settings.plot_q_axis_over_g) == 35


def test_present_unit_cell_2d_requires_complete_geometry(tmp_path):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    del data["unit_cell_2d"]["zeta"]
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match=r"unit_cell_2d is missing: zeta"):
        load_shared_parameters(path)


@pytest.mark.parametrize(
    ("update", "message"),
    [
        ({"q_axis_over_G": [-0.1, 0.0, 0.1]}, "no longer accepted"),
        ({"q_max_over_G": 0.0}, "positive"),
        ({"q_points_per_axis": 1}, "integer >= 2"),
        ({"q_points_per_axis": 3.5}, "integer >= 2"),
        ({"isQuarter": 2}, "integer 0 or 1"),
        ({"isQuarter": True}, "integer 0 or 1"),
        ({"isQuarter": "1"}, "integer 0 or 1"),
        ({"isQuarter": 0, "q_points_per_axis": 18}, "odd integer"),
        ({"isQuarter": 0, "q_points_per_axis": 2}, "odd integer"),
        ({"angle_methods": ["stokes-axis", "unknown"]}, "Invalid"),
        ({"target_bands": ["p3"]}, "Invalid"),
        (
            {"q_max_over_G": 0.0000001, "q_points_per_axis": 3},
            "decimal places",
        ),
        ({"analysis_enabled": False}, "must remain true"),
        ({"b0_nm": 0.0}, "positive"),
        ({"eigenmode_pair_count": 0}, "positive integer"),
    ],
)
def test_invalid_unit_cell_2d_configuration_is_rejected(tmp_path, update, message):
    data = json.loads(Path("scripts/parameter.json").read_text(encoding="utf-8"))
    data["unit_cell_2d"].update(update)
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises((TypeError, ValueError), match=message):
        load_shared_parameters(path)
