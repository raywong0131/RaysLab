from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = (
    REPO_ROOT / "scripts" / "run_sweep" / "optimize_cladding_p2_d2_alignment.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "optimize_cladding_p2_d2_alignment",
        SCRIPT_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_solver_settings_are_fixed_to_requested_mesh_and_safe_mode_count():
    optimize = load_module()
    config = optimize.simulation_config()

    assert optimize.FIXED_CAVITY.to_dict() == {
        "b0_nm": 245.0,
        "eta": 0.960,
        "zeta": 1.156,
    }
    assert optimize.START_CLADDING.to_dict() == {
        "b0_nm": 242.0,
        "eta": 0.980,
        "zeta": 0.930,
    }
    assert config.mesh_auto_size == 5
    assert config.eigenmode_count == optimize.ACTIVE_PARAMETERS.unit_cell_eigenmode_count
    assert config.eigenmode_count > 0
    assert config.eigenfrequency_shift == "c_const/1.55[um]"
    assert optimize.OUT_DIR.name.endswith("_mesh5")


def test_summary_uses_cladding_p2_and_d2_against_requested_cavity_modes():
    optimize = load_module()
    frame = pd.DataFrame(
        {
            "re": [197.0, 198.41, 199.5, 200.48],
            "is_valid": [True, True, True, True],
        }
    )
    selected = {"p1": 0, "p2": 1, "d1": 2, "d2": 3}
    reference = {
        "cavity_p2_frequency_thz": 198.40,
        "cavity_d1_frequency_thz": 200.50,
    }

    result = optimize.summarize_candidate(
        optimize.START_CLADDING,
        "test",
        frame,
        selected,
        reference,
    )

    assert result["p2_error_thz"] == pytest.approx(0.01)
    assert result["d1_d2_soft_error_thz"] == pytest.approx(0.02)
    assert result["p2_aligned"] is True
    assert result["p2_mode_idx"] == 1
    assert result["d2_mode_idx"] == 3


def test_gamma_runner_skips_field_analysis_outside_selection_window(
    monkeypatch,
    tmp_path,
):
    optimize = load_module()
    checked = []
    field_calls = []

    class FakeSimulationRun:
        def __init__(self, _config):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def build_and_run(self, *_args, **_kwargs):
            return None

        def get_eigenfrequencies(self):
            return [
                (180.0, 0.0, 1.0),
                (197.0, 0.0, 1.0),
                (198.0, 0.0, 1.0),
                (200.0, 0.0, 1.0),
                (201.0, 0.0, 1.0),
                (230.0, 0.0, 1.0),
            ]

        def get_2d_fields(self, mode_idx, expression, plane):
            field_calls.append((mode_idx, expression, plane))
            coordinates = np.array([[0.0, 0.0], [1.0, 1.0]])
            return coordinates, np.array([1.0, 0.0])

    def fake_valid(_simulation, mode_idx):
        checked.append(mode_idx)
        return True, {"z_peak_normH_xz": 0.0}

    monkeypatch.setattr(optimize.band, "SimulationRun", FakeSimulationRun)
    monkeypatch.setattr(optimize.band, "_mode_is_valid", fake_valid)
    monkeypatch.setattr(
        optimize.band,
        "interpolate_field",
        lambda _source_points, values, _target_points: values,
    )
    monkeypatch.setattr(
        optimize.band,
        "load_mode_composition",
        lambda _root, mode_idx: {
            "p_weight": float(mode_idx < 3),
            "d_weight": float(mode_idx >= 3),
            "dominant_subspace": "p" if mode_idx < 3 else "d",
        },
    )

    frame = optimize.run_gamma_point(tmp_path, [])

    frequencies = [180.0, 197.0, 198.0, 200.0, 201.0, 230.0]
    expected_checked = [
        mode_idx
        for mode_idx, frequency in enumerate(frequencies)
        if optimize.band.FREQUENCY_MIN_THZ
        <= frequency
        <= optimize.band.FREQUENCY_MAX_THZ
    ]
    assert checked == expected_checked
    for mode_idx in range(len(frequencies)):
        assert bool(frame.loc[mode_idx, "is_valid"]) == (mode_idx in expected_checked)
    assert {
        mode_idx for mode_idx, _expression, _plane in field_calls
    } == set(expected_checked)


def test_ranking_applies_d_soft_constraint_only_inside_p2_tolerance():
    optimize = load_module()
    rows = pd.DataFrame(
        [
            {
                "name": "best_p_but_larger_soft_error",
                "feasible": True,
                "p2_aligned": True,
                "p2_error_thz": 0.001,
                "d1_d2_soft_error_thz": 0.5,
                "evaluation_order": 1,
            },
            {
                "name": "soft_constraint_winner",
                "feasible": True,
                "p2_aligned": True,
                "p2_error_thz": 0.019,
                "d1_d2_soft_error_thz": 0.1,
                "evaluation_order": 2,
            },
            {
                "name": "not_p_aligned",
                "feasible": True,
                "p2_aligned": False,
                "p2_error_thz": 0.021,
                "d1_d2_soft_error_thz": 0.0,
                "evaluation_order": 3,
            },
        ]
    )

    ranked = optimize.rank_candidates(rows)

    assert ranked["name"].tolist() == [
        "soft_constraint_winner",
        "best_p_but_larger_soft_error",
        "not_p_aligned",
    ]


def test_search_changes_coordinates_in_zeta_eta_b0_priority_order(
    monkeypatch,
    tmp_path,
):
    optimize = load_module()
    monkeypatch.setattr(optimize, "OUT_DIR", tmp_path / "scan")
    calls = []

    def fake_evaluator(cell, stage, _reference):
        calls.append((stage, cell.b0_nm, cell.eta, cell.zeta))
        d_error = (
            abs(cell.zeta - 0.938)
            + abs(cell.eta - 0.984)
            + abs(cell.b0_nm - 242.5) / 100.0
        )
        return {
            **cell.to_dict(),
            "stage": stage,
            "status": "p2_aligned",
            "feasible": True,
            "p2_aligned": True,
            "p2_error_thz": 0.01,
            "d1_d2_soft_error_thz": d_error,
        }

    _ranked, best = optimize.run_search(
        {
            "cavity_p2_frequency_thz": 198.4,
            "cavity_d1_frequency_thz": 200.5,
        },
        evaluator=fake_evaluator,
    )

    stages = [stage for stage, *_values in calls]
    first_eta = min(index for index, stage in enumerate(stages) if stage.startswith("eta"))
    first_b0 = min(index for index, stage in enumerate(stages) if stage.startswith("b0"))
    assert all(stage in {"start", "zeta", "zeta_fine"} for stage in stages[:first_eta])
    assert all(stage.startswith("eta") for stage in stages[first_eta:first_b0])
    assert best["zeta"] == pytest.approx(0.938)
    assert best["eta"] == pytest.approx(0.984)
    assert best["b0_nm"] == pytest.approx(242.5)


def test_dry_run_does_not_load_cavity_reference(monkeypatch, tmp_path):
    optimize = load_module()
    monkeypatch.setattr(optimize, "OUT_DIR", tmp_path / "scan")
    monkeypatch.setattr(
        optimize,
        "load_cavity_reference",
        lambda: pytest.fail("dry-run must not start the cavity solve"),
    )

    assert optimize.main(["--dry-run"]) is None
    assert (optimize.OUT_DIR / "scan_config.json").is_file()
