from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.run_main.parameter_config import (
    DEFAULT_PARAMETER_PATH,
    PARAMETER_PATH_ENV,
    parameter_path_from_environment,
)
from scripts.run_sweep.optimize_finite_quarter_q import (
    Candidate,
    build_candidate_parameter_data,
    candidate_key,
    directional_axis_search,
    plot_q_history,
    result_row,
    select_high_q_py_mode,
    write_history,
)


def _modes(q0: float = 400.0, q1: float = 7000.0):
    return [
        {"mode_idx": 0, "frequency_thz": 197.5, "q": q0, "is_valid": True},
        {"mode_idx": 1, "frequency_thz": 198.1, "q": q1, "is_valid": True},
    ]


def test_parameter_path_environment_override_is_optional(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv(PARAMETER_PATH_ENV, raising=False)
    assert parameter_path_from_environment() == DEFAULT_PARAMETER_PATH
    override = tmp_path / "candidate.json"
    monkeypatch.setenv(PARAMETER_PATH_ENV, str(override))
    assert parameter_path_from_environment() == override.resolve()


def test_high_q_mode_is_selected_without_hardcoding_mode1() -> None:
    selected = select_high_q_py_mode(_modes(q0=7200.0, q1=5000.0))
    assert selected["mode_idx"] == 0
    assert selected["q"] == 7200.0


def test_mode_selection_requires_exactly_two_solutions() -> None:
    with pytest.raises(ValueError, match="exactly two modes"):
        select_high_q_py_mode(_modes()[:1])


def test_candidate_parameter_data_changes_only_approved_cladding_fields() -> None:
    baseline = {
        "schema_version": 1,
        "cavity": {"b0_nm": 245.0, "eta": 0.96, "zeta": 1.155},
        "cladding": {"b0_nm": 243.7, "eta": 0.98, "zeta": 0.928},
        "cavity_layers": 20,
        "cladding_layers": 20,
        "mesh_size": 9,
        "unit_cell_eigenmode_count": 6,
        "finite_eigenmode_count": 2,
        "center_frequency_thz": 198.424,
        "cladding_shift_factors": [0.0],
        "cladding_shift_profile": {"kind": "uniform"},
    }
    candidate = Candidate(244.2, 0.98, 0.932)
    updated = build_candidate_parameter_data(baseline, candidate)
    assert updated["cavity"] == baseline["cavity"]
    assert updated["mesh_size"] == 9
    assert updated["finite_eigenmode_count"] == 2
    assert updated["cladding"] == {
        "b0_nm": 244.2,
        "eta": 0.98,
        "zeta": 0.932,
    }
    assert baseline["cladding"]["b0_nm"] == 243.7


def test_directional_search_finds_local_zeta_peak_with_fine_refinement() -> None:
    seen: list[float] = []

    def evaluate(value: float, stage: str):
        seen.append(value)
        q = 7100.0 - 2.0e7 * (value - 0.935) ** 2
        return {
            "stage": stage,
            "b0_nm": 243.7,
            "eta": 0.98,
            "zeta": value,
            "py_q": q,
        }

    best = directional_axis_search(
        start_value=0.928,
        coarse_step=0.004,
        fine_step=0.001,
        bounds=(0.912, 0.944),
        evaluate_value=evaluate,
        stage_prefix="zeta",
    )
    assert best["zeta"] == pytest.approx(0.935)
    assert best["py_q"] == pytest.approx(7100.0)
    assert 0.924 in seen
    assert 0.932 in seen


def test_result_csv_and_plot_keep_both_modes(tmp_path: Path) -> None:
    row = result_row(
        iteration=1,
        stage="zeta_coarse",
        source="worker",
        candidate=Candidate(243.7, 0.98, 0.932),
        modes=_modes(),
        previous_best_q=5070.0,
        elapsed_seconds=12.0,
    )
    csv_path = tmp_path / "q_optimization.csv"
    plot_path = tmp_path / "q_vs_iteration.png"
    write_history(csv_path, [row])
    plot_q_history([row], plot_path)
    with csv_path.open(newline="", encoding="utf-8") as stream:
        saved = next(csv.DictReader(stream))
    assert saved["mode0_q"] == "400.0"
    assert saved["mode1_q"] == "7000.0"
    assert saved["selected_mode_idx"] == "1"
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0


def test_candidate_key_has_stable_three_decimal_identity() -> None:
    assert candidate_key(243.7, 0.98, 0.928) == (
        "b0=243.700|eta=0.980|zeta=0.928"
    )
