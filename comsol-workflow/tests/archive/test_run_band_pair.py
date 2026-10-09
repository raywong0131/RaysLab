from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_main" / "run_band_pair.py"
SCAN_SCRIPT_PATH = REPO_ROOT / "scripts" / "run_sweep" / "scan_unit_cell_p2_gamma.py"
ALIGN_SCRIPT_PATH = REPO_ROOT / "scripts" / "run_sweep" / "cladding_bandgap_alignment.py"
PAIR_OPTIMIZATION_SCRIPT_PATH = (
    REPO_ROOT / "scripts" / "run_sweep" / "optimize_cavity_cladding_cells.py"
)
DOUBLET_ALIGNMENT_SCRIPT_PATH = (
    REPO_ROOT
    / "scripts"
    / "run_sweep"
    / "optimize_cladding_doublet_alignment.py"
)
FINITE_SCRIPT_PATH = REPO_ROOT / "scripts" / "run_main" / "run_finite.py"
STRIP_SCRIPT_PATH = REPO_ROOT / "scripts" / "run_main" / "run_strip_1d.py"


def load_script_module():
    spec = importlib.util.spec_from_file_location("run_band_pair", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_scan_script_module():
    spec = importlib.util.spec_from_file_location("scan_unit_cell_p2_gamma", SCAN_SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_unique_case_points_accepts_custom_branch_names():
    module = load_script_module()
    gamma = {"kx_str": "0.000000", "ky_str": "0.000000"}
    x_point = {"kx_str": "0.010000", "ky_str": "0.000000"}
    y_point = {"kx_str": "0.000000", "ky_str": "0.010000"}
    unique = module._unique_case_points({
        "fit_x": [gamma, x_point],
        "fit_y": [gamma.copy(), y_point],
    })
    assert unique == [gamma, x_point, y_point]


def load_alignment_script_module():
    spec = importlib.util.spec_from_file_location(
        "cladding_bandgap_alignment", ALIGN_SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_pair_optimization_script_module():
    spec = importlib.util.spec_from_file_location(
        "optimize_cavity_cladding_cells", PAIR_OPTIMIZATION_SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_doublet_alignment_script_module():
    spec = importlib.util.spec_from_file_location(
        "optimize_cladding_doublet_alignment", DOUBLET_ALIGNMENT_SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_finite_script_module():
    spec = importlib.util.spec_from_file_location("run_finite", FINITE_SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_strip_script_module():
    spec = importlib.util.spec_from_file_location("run_strip_1d", STRIP_SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


band = load_script_module()

import comsol_workflow.band_connector as band_connector_module
from comsol_workflow.simulation_spatial import hexagon_unit_cell as hexagon_unit_cell_module
from scripts.run_main import parameter_config


def test_unit_cell_p2_gamma_scan_script_exists():
    assert SCAN_SCRIPT_PATH.exists()


def test_cladding_bandgap_alignment_script_exists():
    assert ALIGN_SCRIPT_PATH.exists()


def test_cladding_doublet_alignment_script_exists():
    assert DOUBLET_ALIGNMENT_SCRIPT_PATH.exists()


def test_cladding_doublet_alignment_does_not_auto_run_full_band():
    optimize = load_doublet_alignment_script_module()

    assert optimize.RUN_FULL_BAND_VALIDATION is False
    assert optimize.FIXED_CAVITY.to_dict() == {
        "b0_nm": 245.0,
        "eta": 0.960,
        "zeta": 1.156,
    }
    assert optimize.START_CLADDING.to_dict() == {
        "b0_nm": 240.2,
        "eta": 0.967,
        "zeta": 1.0,
    }
    assert optimize.CLADDING_ZETA == 1.0
    assert optimize.simulation_config().mesh_auto_size == 5
    assert optimize.OUT_DIR.name.endswith("cav(245-0.96-1.156)_zeta1_mesh5")


def test_cavity_cladding_pair_optimization_script_exists():
    assert PAIR_OPTIMIZATION_SCRIPT_PATH.exists()


def test_all_workflow_output_roots_use_four_established_directories():
    expected_names = [
        "unit_cell_band",
        "unit_cell_2D",
        "strip_1d",
        "finite_cavity",
    ]
    assert [path.name for path in parameter_config.WORKFLOW_OUTPUT_ROOTS] == (
        expected_names
    )
    assert all(
        path.parent == parameter_config.OUTPUT_ROOT
        for path in parameter_config.WORKFLOW_OUTPUT_ROOTS
    )
    assert band.OUT_DIR.parent == parameter_config.UNIT_CELL_OUTPUT_ROOT

    scan = load_scan_script_module()
    align = load_alignment_script_module()
    optimize = load_pair_optimization_script_module()
    assert scan.OUT_DIR.parent == parameter_config.UNIT_CELL_OUTPUT_ROOT
    assert align.OUT_DIR.parent == parameter_config.UNIT_CELL_OUTPUT_ROOT
    assert optimize.OUT_DIR.parent == parameter_config.UNIT_CELL_OUTPUT_ROOT
    assert scan.OUT_DIR.name == "unit_cell_p2_gamma_scan"
    assert align.OUT_DIR.name == "cladding_bandgap_alignment"
    assert optimize.OUT_DIR.name == "cell_pair_opt_mesh9"


def test_pair_optimization_uses_user_mode_mapping_and_mesh9():
    optimize = load_pair_optimization_script_module()
    gamma = pd.DataFrame(
        {
            "re": [198.0, 199.0],
            "q": [1.0e6, 2.0e6],
            "px_weight": [0.95, 0.05],
            "py_weight": [0.05, 0.95],
        }
    )
    selected = {"p1": 0, "p2": 1}

    py_label, py_idx, _, py_share = optimize.identify_user_p_mode(
        gamma, selected, "py"
    )
    px_label, px_idx, _, px_share = optimize.identify_user_p_mode(
        gamma, selected, "px"
    )

    assert (py_label, py_idx) == ("p1", 0)
    assert py_share == pytest.approx(0.95)
    assert (px_label, px_idx) == ("p2", 1)
    assert px_share == pytest.approx(0.95)
    assert optimize.simulation_config().mesh_auto_size == 9


def test_pair_optimization_output_path_leaves_windows_headroom():
    optimize = load_pair_optimization_script_module()
    deepest = (
        optimize.cavity_candidate_dir(0.950, 1.146)
        / "k_points"
        / "kx=0.000000_ky=0.000000"
        / "eigenmodes"
        / "00_Hz_center.parquet"
    )

    assert len(str(deepest.resolve())) < 240


def test_pair_optimization_cavity_frequency_constraint_is_hard():
    optimize = load_pair_optimization_script_module()
    target = pd.DataFrame(
        {
            "direction": ["gamma_m", "gamma_m", "gamma_k"],
            "k_norm": [0.0, 0.01, 0.01],
            "q": [1.0e7, 1.0e6, 1.0e6],
            "p_weight": [0.99, 0.99, 0.99],
            "is_valid": [True, True, True],
            "match_status": ["gamma", "matched", "matched"],
        }
    )
    metadata = {
        "target_gamma_frequency_thz": (
            optimize.CAVITY_REFERENCE_FREQUENCY_THZ
            + optimize.CAVITY_FREQUENCY_TOLERANCE_THZ
            + 0.001
        ),
        "target_gamma_q": 1.0e7,
        "target_user_mode_share": 0.99,
        "target_band_label": "p2",
    }

    summary = optimize.summarize_cavity_candidate(
        0.96, 1.156, target, metadata, "test"
    )

    assert summary["feasible"] is False
    assert summary["status"] == "frequency_drift_outside_tolerance"


def test_pair_validation_is_written_to_99_config(tmp_path, monkeypatch):
    optimize = load_pair_optimization_script_module()
    monkeypatch.setattr(optimize, "OUT_DIR", tmp_path)

    def gamma_frame(frequencies, px_weights, py_weights):
        return pd.DataFrame(
            {
                "direction": ["gamma_m", "gamma_m"],
                "band_label": ["p1", "p2"],
                "mode_idx": [0, 1],
                "k_norm": [0.0, 0.0],
                "re": frequencies,
                "q": [1.0e6, 2.0e6],
                "px_weight": px_weights,
                "py_weight": py_weights,
            }
        )

    monkeypatch.setattr(
        optimize.band,
        "run_band_workflow",
        lambda *_args, **_kwargs: {
            "cavity": gamma_frame([198.0, 199.0], [0.95, 0.05], [0.05, 0.95]),
            "cladding": gamma_frame([197.9, 198.1], [0.05, 0.95], [0.95, 0.05]),
        },
    )

    optimize.validate_best_pair(
        {"eta": 0.96, "zeta": 1.155},
        {"b0_nm": 243.7, "eta": 0.98, "zeta": 0.928},
    )

    root = tmp_path / "best_pair_full_band"
    assert (root / "99_config" / "pair_validation.json").is_file()
    assert not (root / "pair_validation.json").exists()


def test_pair_optimization_loads_user_accepted_cavity_without_rerun(
    tmp_path, monkeypatch
):
    optimize = load_pair_optimization_script_module()
    monkeypatch.setattr(optimize, "OUT_DIR", tmp_path)
    candidate = optimize.cavity_candidate_dir(
        optimize.ACCEPTED_CAVITY_ETA, optimize.ACCEPTED_CAVITY_ZETA
    )
    candidate.mkdir(parents=True)
    summary = {
        "eta": optimize.ACCEPTED_CAVITY_ETA,
        "zeta": optimize.ACCEPTED_CAVITY_ZETA,
        "feasible": True,
        "frequency_drift_thz": -0.025,
        "gamma_frequency_thz": 198.424,
        "gamma_q": 1.17e7,
        "peak_at_gamma": False,
    }
    (candidate / "config.json").write_text(
        json.dumps(
            {
                "summary": summary,
                "eigenmode_count": optimize.band.EIGENMODE_COUNT,
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame({"k_norm": [0.0], "q": [1.17e7]}).to_csv(
        candidate / "selected_user_py.csv", index=False
    )

    selected, rows = optimize.load_accepted_cavity_candidate()

    assert selected["selection_status"] == "user_accepted_closest_to_gamma"
    assert selected["gamma_frequency_thz"] == pytest.approx(198.424)
    assert rows["q"].iloc[0] == pytest.approx(1.17e7)
    assert (tmp_path / "cavity" / "best" / "best_config.json").is_file()


def make_alignment_gamma_modes(
    p_frequencies,
    p_py_shares,
    d_frequencies,
):
    rows = []
    for frequency, py_share in zip(p_frequencies, p_py_shares):
        rows.append(
            {
                "re": frequency,
                "is_valid": True,
                "dominant_subspace": "p",
                "p_weight": 1.0,
                "d_weight": 0.0,
                "px_weight": 1.0 - py_share,
                "py_weight": py_share,
                "dx_weight": 0.0,
                "dy_weight": 0.0,
            }
        )
    for mode_number, frequency in enumerate(d_frequencies):
        rows.append(
            {
                "re": frequency,
                "is_valid": True,
                "dominant_subspace": "d",
                "p_weight": 0.0,
                "d_weight": 1.0,
                "px_weight": 0.0,
                "py_weight": 0.0,
                "dx_weight": float(mode_number == 0),
                "dy_weight": float(mode_number == 1),
            }
        )
    return pd.DataFrame(rows)


def test_cladding_alignment_parameter_translation_and_grids():
    align = load_alignment_script_module()

    params = align.make_cladding_params(234.2, 0.958, 0.837)

    assert params["name"] == "cladding_bandgap_alignment"
    assert params["r_f0"] == pytest.approx(0.958)
    assert params["b_square_f0"] == pytest.approx((234.2 / 230.0) ** 2)
    assert params["b_square_f3"] == pytest.approx((0.837**2 - 1.0) / 2.0)
    assert align.coarse_b0_values() == pytest.approx(tuple(range(230, 241)))
    assert align.coarse_zeta_values() == pytest.approx(
        (0.830, 0.835, 0.840, 0.845, 0.850)
    )
    assert align.fine_zeta_values(0.840) == pytest.approx(
        tuple(value / 1000 for value in range(835, 846))
    )
    assert align.coarse_eta_values() == pytest.approx(
        (0.950, 0.955, 0.960, 0.965, 0.970)
    )
    assert align.fine_eta_values(0.950) == pytest.approx(
        tuple(value / 1000 for value in range(950, 956))
    )


def test_cladding_alignment_uses_py_composition_and_frequency_gap_edges():
    align = load_alignment_script_module()
    modes = make_alignment_gamma_modes(
        p_frequencies=(194.8, 194.5),
        p_py_shares=(0.1, 0.99),
        d_frequencies=(198.2, 197.4),
    )
    selected = {"p1": 0, "p2": 1, "d1": 2, "d2": 3}

    summary = align.summarize_candidate(
        235.0,
        0.960,
        0.840,
        "test",
        modes,
        selected,
        cavity_py_frequency_thz=194.45,
    )

    assert summary["cladding_py_mode_idx"] == 1
    assert summary["cladding_py_frequency_thz"] == pytest.approx(194.5)
    assert summary["gap_lower_edge_thz"] == pytest.approx(194.8)
    assert summary["gap_upper_edge_thz"] == pytest.approx(197.4)
    assert summary["gap_width_thz"] == pytest.approx(2.6)
    assert summary["midgap_thz"] == pytest.approx(196.1)
    assert summary["alignment_error_thz"] == pytest.approx(0.05)
    assert summary["feasible"]
    assert summary["within_alignment_tolerance"]


def test_cladding_alignment_rejects_gap_below_two_thz():
    align = load_alignment_script_module()
    modes = make_alignment_gamma_modes(
        p_frequencies=(195.0, 194.5),
        p_py_shares=(0.1, 0.99),
        d_frequencies=(196.8, 196.5),
    )

    summary = align.summarize_candidate(
        235.0,
        0.960,
        0.840,
        "test",
        modes,
        {"p1": 0, "p2": 1, "d1": 2, "d2": 3},
        cavity_py_frequency_thz=194.5,
    )

    assert summary["gap_width_thz"] == pytest.approx(1.5)
    assert not summary["feasible"]
    assert summary["status"] == "gap_below_minimum"


def test_cladding_alignment_ranking_uses_midgap_only_after_alignment():
    align = load_alignment_script_module()
    selected = {"p1": 0, "p2": 1, "d1": 2, "d2": 3}

    def candidate(b0, py_frequency, other_p, lower_d, stage="b0"):
        modes = make_alignment_gamma_modes(
            p_frequencies=(other_p, py_frequency),
            p_py_shares=(0.05, 0.99),
            d_frequencies=(lower_d + 0.5, lower_d),
        )
        return align.summarize_candidate(
            b0,
            0.960,
            0.840,
            stage,
            modes,
            selected,
            cavity_py_frequency_thz=194.5,
        )

    aligned_worse_midgap = candidate(235.0, 194.52, 194.52, 197.0)
    aligned_better_midgap = candidate(235.1, 194.54, 194.60, 196.60)
    unaligned = candidate(235.2, 194.60, 194.60, 196.80)
    invalid_gap = candidate(235.3, 194.50, 195.00, 196.50)

    ranked = align.rank_candidates(
        pd.DataFrame(
            [aligned_worse_midgap, aligned_better_midgap, unaligned, invalid_gap]
        )
    )

    assert ranked["b0_nm"].tolist() == pytest.approx([235.1, 235.0, 235.2, 235.3])
    assert ranked["rank"].tolist() == [1, 2, 3, 4]


def test_cladding_alignment_candidate_reuses_gamma_cache(tmp_path, monkeypatch):
    align = load_alignment_script_module()
    assert hasattr(align, "evaluate_candidate")
    monkeypatch.setattr(align, "OUT_DIR", tmp_path / "alignment")
    monkeypatch.setattr(align.band, "get_hole_params", lambda _params: np.ones((6, 4)))
    holes = [np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]]) for _ in range(6)]
    monkeypatch.setattr(
        align.band,
        "create_hexagon_design",
        lambda *_args: (np.zeros((6, 2)), holes, [{"min_dist": 0.03}] * 6),
    )
    monkeypatch.setattr(
        align.band,
        "visualize_hexagon_design",
        lambda *_args, filename, **_kwargs: Path(filename).write_bytes(b"geometry"),
    )
    solve_calls = []

    def fake_run_k_point(case_dir, _holes, point, _config):
        point_dir = align.band.k_point_dir(case_dir, point)
        csv_path = point_dir / "eigenfrequencies.csv"
        if csv_path.exists():
            return pd.read_csv(csv_path)
        solve_calls.append((point["kx_str"], point["ky_str"]))
        point_dir.mkdir(parents=True, exist_ok=True)
        frame = make_alignment_gamma_modes(
            p_frequencies=(194.0, 194.5),
            p_py_shares=(0.05, 0.99),
            d_frequencies=(197.0, 198.0),
        )
        frame.to_csv(csv_path, index=False)
        return frame

    monkeypatch.setattr(align.band, "run_k_point", fake_run_k_point)
    cavity_reference = {"py_frequency_thz": 194.48}

    first = align.evaluate_candidate(235.0, 0.960, 0.840, "b0", cavity_reference)
    second = align.evaluate_candidate(235.0, 0.960, 0.840, "b0", cavity_reference)

    expected_dir = (
        align.scan_root()
        / "candidates"
        / "clad_b0_235.0_eta_0.960_zeta_0.840"
    )
    assert solve_calls == [("0.000000", "0.000000")]
    assert first["status"] == second["status"] == "ok_aligned"
    assert expected_dir.joinpath("config.json").exists()
    config = pd.read_json(expected_dir / "config.json", typ="series")
    assert config["params"]["r_f0"] == pytest.approx(0.960)
    assert config["selected_gamma_modes"] == {"p1": 0, "p2": 1, "d1": 2, "d2": 3}


def fake_alignment_summary(b0_nm, eta, zeta, stage, signed_error, midgap_error):
    alignment_error = abs(float(signed_error))
    return {
        "b0_nm": round(float(b0_nm), 1),
        "eta": round(float(eta), 3),
        "zeta": round(float(zeta), 3),
        "stage": stage,
        "status": "ok_aligned" if alignment_error <= 0.05 else "ok_unaligned",
        "feasible": True,
        "within_alignment_tolerance": alignment_error <= 0.05,
        "cavity_py_frequency_thz": 194.5,
        "cladding_py_mode_idx": 1,
        "cladding_py_frequency_thz": 194.5 + signed_error,
        "cladding_py_share": 0.99,
        "signed_alignment_error_thz": signed_error,
        "alignment_error_thz": alignment_error,
        "gap_lower_edge_thz": 195.0,
        "gap_upper_edge_thz": 197.5,
        "gap_width_thz": 2.5,
        "midgap_thz": 196.25,
        "midgap_error_thz": float(midgap_error),
    }


def test_cladding_alignment_nested_search_refines_b0_and_zeta_without_eta(
    tmp_path, monkeypatch
):
    align = load_alignment_script_module()
    assert hasattr(align, "run_search")
    monkeypatch.setattr(align, "OUT_DIR", tmp_path / "alignment")
    calls = []

    def fake_evaluate(b0_nm, eta, zeta, stage, _cavity_reference):
        key = (round(b0_nm, 1), round(eta, 3), round(zeta, 3))
        calls.append(key)
        target_b0 = 234.0 + 10.0 * (zeta - 0.842)
        signed_error = 0.2 * (b0_nm - target_b0)
        return fake_alignment_summary(
            b0_nm,
            eta,
            zeta,
            stage,
            signed_error,
            abs(zeta - 0.842),
        )

    summary, best = align.run_search(
        evaluate=fake_evaluate,
        cavity_reference={"py_frequency_thz": 194.5},
    )

    assert best["b0_nm"] == pytest.approx(234.0, abs=0.1)
    assert best["zeta"] == pytest.approx(0.842)
    assert best["eta"] == pytest.approx(0.960)
    assert len(calls) == len(set(calls))
    assert len(calls) < 100
    assert summary["eta"].unique().tolist() == pytest.approx([0.960])
    assert not summary["stage"].str.startswith("eta_").any()
    assert align.scan_root().joinpath("scan_summary.csv").exists()
    assert align.scan_root().joinpath("optimization_trace.png").exists()


def test_cladding_alignment_uses_eta_only_when_b0_zeta_cannot_align(
    tmp_path, monkeypatch
):
    align = load_alignment_script_module()
    assert hasattr(align, "run_search")
    monkeypatch.setattr(align, "OUT_DIR", tmp_path / "alignment")

    def fake_evaluate(b0_nm, eta, zeta, stage, _cavity_reference):
        signed_error = 100.0 * (eta - 0.958)
        return fake_alignment_summary(
            b0_nm,
            eta,
            zeta,
            stage,
            signed_error,
            abs(zeta - 0.840),
        )

    summary, best = align.run_search(
        evaluate=fake_evaluate,
        cavity_reference={"py_frequency_thz": 194.5},
    )

    assert best["eta"] == pytest.approx(0.958)
    assert best["within_alignment_tolerance"]
    assert summary["stage"].str.startswith("eta_").any()


def test_cladding_alignment_trace_has_parameter_axes_and_both_constraints(
    tmp_path, monkeypatch
):
    align = load_alignment_script_module()
    rows = pd.DataFrame(
        [
            {
                "evaluation_order": 1,
                "stage": "b0",
                "b0_nm": 235.0,
                "zeta": 0.840,
                "eta": 0.960,
                "alignment_error_thz": 0.04,
                "midgap_error_thz": 1.2,
                "gap_width_thz": 2.5,
            }
        ]
    )
    original_close = align.plt.close
    monkeypatch.setattr(align.plt, "close", lambda *_args, **_kwargs: None)

    align.plot_optimization_trace(rows, tmp_path / "trace.png")

    figure = align.plt.gcf()
    assert [axis.get_xlabel() for axis in figure.axes] == [
        r"Cladding $b_0$ (nm)",
        r"Cladding $\zeta$",
        r"Cladding $\eta$",
    ]
    legend_labels = [text.get_text() for text in figure.legends[0].get_texts()]
    assert "alignment tolerance" in legend_labels
    assert "minimum gap" in legend_labels
    original_close(figure)


def make_alignment_full_band_frame(cell, close_gap=False):
    rows = []
    gamma_frequencies = {
        "p1": 194.5 if cell == "cavity" else 194.0,
        "p2": 195.5 if cell == "cavity" else 194.52,
        "d1": 198.0 if cell == "cavity" else 197.0,
        "d2": 200.0 if cell == "cavity" else 198.0,
    }
    for direction in ("gamma_m", "gamma_k"):
        for path_coordinate in (-1.0, 0.0, 1.0):
            if direction == "gamma_m" and path_coordinate > 0.0:
                continue
            if direction == "gamma_k" and path_coordinate < 0.0:
                continue
            for label, gamma_frequency in gamma_frequencies.items():
                frequency = gamma_frequency
                if not np.isclose(path_coordinate, 0.0):
                    frequency += -0.2 if label.startswith("p") else 0.2
                if close_gap and label == "d1" and path_coordinate > 0.0:
                    frequency = 194.4
                is_py = label == "p1" if cell == "cavity" else label == "p2"
                rows.append(
                    {
                        "cell": cell,
                        "direction": direction,
                        "path_coordinate": path_coordinate,
                        "band_label": label,
                        "match_status": "gamma" if np.isclose(path_coordinate, 0.0) else "matched",
                        "re": frequency,
                        "is_valid": True,
                        "dominant_subspace": label[0],
                        "p_weight": float(label.startswith("p")),
                        "d_weight": float(label.startswith("d")),
                        "px_weight": float(label.startswith("p") and not is_py),
                        "py_weight": float(label.startswith("p") and is_py),
                        "dx_weight": float(label == "d1"),
                        "dy_weight": float(label == "d2"),
                    }
                )
    return pd.DataFrame(rows)


def test_cladding_alignment_full_band_handoff_and_path_gap(tmp_path, monkeypatch):
    align = load_alignment_script_module()
    assert hasattr(align, "validate_full_band")
    monkeypatch.setattr(align, "OUT_DIR", tmp_path / "alignment")
    cavity = make_alignment_full_band_frame("cavity")
    cladding = make_alignment_full_band_frame("cladding")
    calls = []

    def fake_run_band_workflow(output_dir, cavity_params, cladding_params):
        calls.append((Path(output_dir), cavity_params, cladding_params))
        return {"cavity": cavity, "cladding": cladding}

    monkeypatch.setattr(align.band, "run_band_workflow", fake_run_band_workflow)
    best = {"b0_nm": 234.2, "eta": 0.960, "zeta": 0.842, "stage": "zeta_fine"}

    validation = align.validate_full_band(
        best,
        cavity_reference={"py_frequency_thz": 194.5},
    )

    expected_dir = align.scan_root() / "best" / "full_band"
    assert calls[0][0] == expected_dir
    assert calls[0][1] == align.band.BULK_PARAMS
    assert calls[0][2] == align.make_cladding_params(234.2, 0.960, 0.842)
    expected_gap = (
        cladding.loc[cladding["band_label"].str.startswith("d"), "re"].min()
        - cladding.loc[cladding["band_label"].str.startswith("p"), "re"].max()
    )
    assert validation["full_path_gap_thz"] == pytest.approx(expected_gap)
    assert validation["gamma_summary"]["alignment_error_thz"] == pytest.approx(0.02)
    assert align.scan_root().joinpath("best", "validation.json").exists()


def test_cladding_alignment_full_band_rejects_closed_path_gap(tmp_path, monkeypatch):
    align = load_alignment_script_module()
    assert hasattr(align, "validate_full_band")
    monkeypatch.setattr(align, "OUT_DIR", tmp_path / "alignment")
    monkeypatch.setattr(
        align.band,
        "run_band_workflow",
        lambda *_args, **_kwargs: {
            "cavity": make_alignment_full_band_frame("cavity"),
            "cladding": make_alignment_full_band_frame("cladding", close_gap=True),
        },
    )

    with pytest.raises(RuntimeError, match="Full-path cladding p-d gap closed"):
        align.validate_full_band(
            {"b0_nm": 234.2, "eta": 0.960, "zeta": 0.842, "stage": "zeta_fine"},
            cavity_reference={"py_frequency_thz": 194.5},
        )


def test_zeta_scan_parameter_translation_uses_existing_fourier_definition():
    scan = load_scan_script_module()

    params = scan.make_cavity_params(1.1)

    assert params["r_f0"] == pytest.approx(0.96)
    assert params["b_square_f0"] == pytest.approx((235.0 / 230.0) ** 2)
    assert params["b_square_f3"] == pytest.approx((1.1**2 - 1.0) / 2.0)


def test_zeta_scan_builds_inclusive_coarse_and_fine_grids():
    scan = load_scan_script_module()

    coarse = scan.coarse_zeta_values()
    fine = scan.fine_zeta_values(1.15)
    edge = scan.fine_zeta_values(1.10)

    assert coarse == pytest.approx(tuple(value / 1000 for value in range(1100, 1201, 10)))
    assert fine == pytest.approx(tuple(value / 1000 for value in range(1140, 1161)))
    assert edge == pytest.approx(tuple(value / 1000 for value in range(1100, 1111)))
    assert all(round(value, 3) == value for value in coarse + fine + edge)


def test_zeta_scan_direction_points_use_signed_plot_directions():
    scan = load_scan_script_module()
    magnitudes = (0.0, 0.01, 0.2)

    gamma_m = scan.make_direction_points("gamma_m", magnitudes)
    gamma_k = scan.make_direction_points("gamma_k", magnitudes)

    assert [point["kx"] for point in gamma_m] == pytest.approx([0.0, 0.0, 0.0])
    assert [point["ky"] for point in gamma_m] == pytest.approx(magnitudes)
    assert [point["kx"] for point in gamma_k] == pytest.approx(magnitudes)
    assert [point["ky"] for point in gamma_k] == pytest.approx([0.0, 0.0, 0.0])
    assert [point["path_coordinate"] for point in gamma_m] == pytest.approx(
        [0.0, -0.01, -0.2]
    )
    assert [point["path_coordinate"] for point in gamma_k] == pytest.approx(magnitudes)


def make_p2_scan_rows(gamma_q=1.0e7, gamma_frequency=194.0, gamma_p_weight=0.95):
    rows = []
    for direction, q_values in (
        ("gamma_m", (gamma_q, 0.8 * gamma_q, 0.5 * gamma_q)),
        ("gamma_k", (gamma_q, 0.9 * gamma_q, 0.6 * gamma_q)),
    ):
        for k_norm, q_value in zip((0.0, 0.01, 0.02), q_values):
            rows.append(
                {
                    "direction": direction,
                    "band_label": "p2",
                    "k_norm": k_norm,
                    "q": q_value,
                    "re": gamma_frequency if k_norm == 0.0 else gamma_frequency - k_norm,
                    "p_weight": gamma_p_weight,
                    "is_valid": True,
                    "match_status": "gamma" if k_norm == 0.0 else "matched",
                }
            )
    return pd.DataFrame(rows)


def test_zeta_scan_summary_accepts_constrained_gamma_peak():
    scan = load_scan_script_module()

    summary = scan.summarize_p2_candidate(1.15, make_p2_scan_rows(), "coarse")

    assert summary["status"] == "ok"
    assert summary["feasible"]
    assert summary["peak_at_gamma"]
    assert summary["peak_excess_dex"] < 0.0
    assert summary["peak_k_signed"] == 0.0
    assert summary["gamma_frequency_thz"] == pytest.approx(194.0)
    assert summary["gamma_p_weight"] == pytest.approx(0.95)


def test_zeta_scan_summary_records_off_gamma_peak_and_direction():
    scan = load_scan_script_module()
    rows = make_p2_scan_rows()
    rows.loc[(rows["direction"] == "gamma_m") & (rows["k_norm"] == 0.01), "q"] = 2.0e7

    summary = scan.summarize_p2_candidate(1.15, rows, "fine")

    assert summary["status"] == "ok"
    assert not summary["peak_at_gamma"]
    assert summary["peak_excess_dex"] == pytest.approx(np.log10(2.0))
    assert summary["peak_k_signed"] == pytest.approx(-0.01)
    assert summary["peak_direction"] == "gamma_m"


@pytest.mark.parametrize(
    ("mutator", "expected_status"),
    [
        (lambda rows: rows.assign(re=190.0), "gamma_frequency_outside_gap"),
        (lambda rows: rows.assign(p_weight=0.85), "gamma_p_weight_below_min"),
        (
            lambda rows: rows.assign(
                match_status=np.where(rows["k_norm"] == 0.02, "unmatched", rows["match_status"])
            ),
            "unmatched_p2",
        ),
        (lambda rows: rows.assign(q=np.nan), "invalid_q"),
    ],
)
def test_zeta_scan_summary_rejects_invalid_candidates(mutator, expected_status):
    scan = load_scan_script_module()

    summary = scan.summarize_p2_candidate(1.15, mutator(make_p2_scan_rows()), "coarse")

    assert summary["status"] == expected_status
    assert not summary["feasible"]


def test_zeta_scan_ranking_prefers_gamma_peaks_then_higher_gamma_q():
    scan = load_scan_script_module()
    summaries = pd.DataFrame(
        [
            scan.summarize_p2_candidate(1.10, make_p2_scan_rows(gamma_q=1.0e7), "coarse"),
            scan.summarize_p2_candidate(1.11, make_p2_scan_rows(gamma_q=2.0e7), "coarse"),
            scan.summarize_p2_candidate(
                1.12,
                make_p2_scan_rows().assign(
                    q=lambda frame: np.where(frame["k_norm"] == 0.01, 1.1e7, frame["q"])
                ),
                "coarse",
            ),
        ]
    )

    ranked = scan.rank_candidates(summaries)

    assert ranked["zeta"].tolist() == pytest.approx([1.11, 1.10, 1.12])
    assert ranked["rank"].tolist() == [1, 2, 3]


def test_zeta_scan_run_candidate_reuses_complete_k_points(tmp_path, monkeypatch):
    scan = load_scan_script_module()
    monkeypatch.setattr(scan, "OUT_DIR", tmp_path / "unit_cell_p2_gamma_scan")
    monkeypatch.setattr(
        scan.band,
        "visualize_hexagon_design",
        lambda *_args, filename, **_kwargs: Path(filename).write_bytes(b"geometry"),
    )
    monkeypatch.setattr(
        scan.band,
        "load_hz_center",
        lambda k_dir, mode_idx: (str(k_dir), int(mode_idx)),
    )
    monkeypatch.setattr(
        scan.band,
        "field_overlap",
        lambda previous, candidate: 0.99 if previous[1] == candidate[1] else 0.01,
    )
    solve_calls = []

    def fake_run_k_point(case_dir, _holes, point, _simulation_config):
        output_dir = scan.band.k_point_dir(case_dir, point)
        if scan.band.k_point_is_complete(output_dir):
            return pd.read_csv(output_dir / "eigenfrequencies.csv")
        solve_calls.append((point["kx_str"], point["ky_str"]))
        output_dir.joinpath("eigenmodes").mkdir(parents=True, exist_ok=True)
        frame = make_tracked_modes([193.0, 194.0, 196.0, 197.0])
        frame.loc[1, "p_weight"] = 0.95
        frame.loc[1, "q"] = 2.0e7 if point["k_norm"] == 0.0 else 1.5e7
        frame.to_csv(output_dir / "eigenfrequencies.csv", index=False)
        for mode_idx in frame.index:
            pd.DataFrame(
                {
                    "x": [0.0, 0.1, 0.0],
                    "y": [0.0, 0.0, 0.1],
                    "re": [1.0, 0.5, 0.25],
                    "im": [0.0, 0.1, 0.2],
                }
            ).to_parquet(output_dir / "eigenmodes" / f"{mode_idx:02d}_Hz_center.parquet")
        return frame

    monkeypatch.setattr(scan.band, "run_k_point", fake_run_k_point)

    first_rows, first_summary = scan.run_candidate(1.15, "coarse", (0.0, 0.01))
    first_solve_count = len(solve_calls)
    second_rows, second_summary = scan.run_candidate(1.15, "coarse", (0.0, 0.01))

    case_dir = scan.scan_root() / "zeta_1.150"
    assert first_solve_count == 3
    assert len(solve_calls) == first_solve_count
    assert set(first_rows["band_label"]) == {"p2"}
    assert first_summary["status"] == "ok"
    assert first_summary["peak_at_gamma"]
    assert second_summary == first_summary
    pd.testing.assert_frame_equal(first_rows, second_rows)
    for relative_path in ("geometry.png", "config.json", "selected_p2.csv"):
        assert (case_dir / relative_path).exists()


def test_zeta_scan_main_runs_coarse_fine_and_dense_validation(tmp_path, monkeypatch):
    scan = load_scan_script_module()
    monkeypatch.setattr(scan, "OUT_DIR", tmp_path / "unit_cell_p2_gamma_scan")
    monkeypatch.setattr(scan, "coarse_zeta_values", lambda: (1.100, 1.110, 1.120))
    monkeypatch.setattr(
        scan,
        "fine_zeta_values",
        lambda _center: (1.100, 1.110, 1.111, 1.120),
    )
    monkeypatch.setattr(scan, "FINAL_K_MAGNITUDES", (0.0, 0.01, 0.02))
    calls = []

    def fake_run_candidate(zeta, stage, magnitudes):
        calls.append((round(float(zeta), 3), stage, tuple(magnitudes)))
        gamma_q = {1.100: 1.0e7, 1.110: 2.0e7, 1.111: 3.0e7, 1.120: 0.8e7}[
            round(float(zeta), 3)
        ]
        rows = make_p2_scan_rows(gamma_q=gamma_q)
        rows.insert(0, "zeta", round(float(zeta), 3))
        rows.insert(1, "stage", stage)
        return rows, scan.summarize_p2_candidate(zeta, rows, stage)

    monkeypatch.setattr(scan, "run_candidate", fake_run_candidate)

    best = scan.main()

    root = scan.scan_root()
    assert best["zeta"] == pytest.approx(1.111)
    assert [zeta for zeta, stage, _ in calls if stage == "coarse"] == pytest.approx(
        [1.100, 1.110, 1.120]
    )
    assert [zeta for zeta, stage, _ in calls if stage == "fine"] == pytest.approx([1.111])
    assert len([1 for _, stage, _ in calls if stage == "final"]) == 3
    for relative_path in (
        "scan_config.json",
        "scan_summary.csv",
        "zeta_scan.png",
        "best/best_config.json",
        "best/selected_p2.csv",
        "best/cavity_p2_q_vs_k.png",
    ):
        path = root / relative_path
        assert path.exists(), relative_path
        assert path.stat().st_size > 0, relative_path


def test_make_branch_k_points_uses_full_hexagonal_endpoints():
    gamma_m = band.make_branch_k_points("gamma_m", 31)
    gamma_k = band.make_branch_k_points("gamma_k", 31)

    assert gamma_m[0]["kx"] == 0.0
    assert gamma_m[0]["ky"] == 0.0
    assert gamma_m[-1]["ky"] == pytest.approx(0.5)
    assert gamma_k[-1]["kx"] == pytest.approx(1.0 / np.sqrt(3.0))
    assert gamma_k[-1]["ky"] == 0.0


def test_cavity_sampling_adds_dense_points_without_duplicates():
    branches = band.case_k_points("cavity")

    for direction, axis in (("gamma_m", "ky"), ("gamma_k", "kx")):
        values = [point[axis] for point in branches[direction]]
        for expected in band.Q_FIT_K_MAGNITUDES:
            assert any(value == pytest.approx(expected) for value in values)
        keys = [(point["kx_str"], point["ky_str"]) for point in branches[direction]]
        assert len(keys) == len(set(keys))


def test_default_cavity_and_cladding_c2_parameters_and_output_directory():
    assert band.OUT_DIR.name == (
        f"unitcell_band_{band.ACTIVE_PARAMETERS.unit_cell_case_parameter_label}"
    )
    assert band.MESH_AUTO_SIZE == band.ACTIVE_PARAMETERS.mesh_size
    assert band.EIGENMODE_COUNT == band.ACTIVE_PARAMETERS.unit_cell_eigenmode_count
    assert band.EIGENMODE_COUNT > 0
    assert band.BULK_PARAMS == band.ACTIVE_PARAMETERS.cavity.to_fourier_params(
        "cavity_p_bic"
    )
    assert band.CLADDING_PARAMS == band.ACTIVE_PARAMETERS.cladding.to_fourier_params(
        "bulk_gap_centered"
    )


def _shared_parameter_payload(center_frequency_thz=198.0):
    return {
        "schema_version": 1,
        "cavity": {"b0_nm": 245.0, "eta": 0.96, "zeta": 1.156},
        "cladding": {"b0_nm": 242.0, "eta": 0.98, "zeta": 0.93},
        "cavity_layers": 10,
        "cladding_layers": 10,
        "mesh_size": 9,
        "unit_cell_eigenmode_count": 6,
        "finite_eigenmode_count": 1,
        "center_frequency_thz": center_frequency_thz,
        "cladding_shift_factors": [0.03],
    }


def test_shared_parameter_json_resolves_geometry_mesh_frequency_and_label(tmp_path):
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(_shared_parameter_payload()), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.case_parameter_label == (
        "cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9"
    )
    assert shared.unit_cell_eigenmode_count == 6
    assert shared.finite_eigenmode_count == 1
    assert shared.finite_geometry == "hex"
    assert shared.unit_cell_case_parameter_label == (
        "cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9"
    )
    assert shared.eigenfrequency_shift == "198 [THz]"
    assert shared.cladding_shift_factors == (0.03,)
    assert shared.cladding_y_over_x_shift_ratio == 1.0
    assert shared.finite_cladding_shift_factor_pairs == ((0.03, 0.03),)
    assert shared.cavity_layers == 10
    assert shared.cladding_layers == 10
    assert shared.structure_series_label("finite") == (
        "finite_10-10_"
        "cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9"
    )
    assert shared.cavity.to_fourier_params("cavity")["b_square_f0"] == pytest.approx(
        (245.0 / 230.0) ** 2
    )


def test_shared_parameter_json_selects_square_finite_geometry(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["finite_geometry"] = "square"
    path.write_text(json.dumps(payload), encoding="utf-8")
    shared = parameter_config.load_shared_parameters(path)
    assert shared.finite_geometry == "square"
    assert shared.structure_series_label("finite") == (
        "finite_10-10_"
        "cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9_square"
    )
    assert shared.structure_series_label("strip1d") == (
        "strip1d_10-10_"
        "cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9"
    )


@pytest.mark.parametrize("value", ["Square", "triangle", 1])
def test_shared_parameter_json_rejects_invalid_finite_geometry(tmp_path, value):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["finite_geometry"] = value
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises((TypeError, ValueError), match="finite_geometry"):
        parameter_config.load_shared_parameters(path)


def test_cavity_p2_update_changes_only_frequency_and_cleans_temporary_file(tmp_path):
    path = tmp_path / "parameter.json"
    original = _shared_parameter_payload()
    original["custom_future_field"] = {"keep": True}
    path.write_text(json.dumps(original), encoding="utf-8")
    expected = parameter_config.load_shared_parameters(path)

    updated = parameter_config.update_center_frequency_from_cavity_p2(
        198.4493726243485,
        expected_parameters=expected,
        path=path,
    )

    saved = json.loads(path.read_text(encoding="utf-8"))
    assert updated.center_frequency_thz == pytest.approx(198.4493726243485)
    assert saved["center_frequency_thz"] == pytest.approx(198.4493726243485)
    assert saved["custom_future_field"] == {"keep": True}
    assert saved["center_frequency_source"] == {
        "workflow": "run_band_pair",
        "cell": "cavity",
        "band": "p2",
        "k_point": "Gamma",
    }
    assert not (tmp_path / ".out").exists()


def test_cavity_p2_update_rejects_geometry_changed_during_run(tmp_path):
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(_shared_parameter_payload()), encoding="utf-8")
    expected = parameter_config.load_shared_parameters(path)
    changed = _shared_parameter_payload()
    changed["cavity"]["b0_nm"] = 246.0
    path.write_text(json.dumps(changed), encoding="utf-8")

    with pytest.raises(RuntimeError, match="geometry, mesh, or unit-cell mode count changed"):
        parameter_config.update_center_frequency_from_cavity_p2(
            198.5,
            expected_parameters=expected,
            path=path,
        )

    assert json.loads(path.read_text(encoding="utf-8"))["center_frequency_thz"] == 198.0


def test_cavity_p2_update_rejects_mode_count_changed_during_run(tmp_path):
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(_shared_parameter_payload()), encoding="utf-8")
    expected = parameter_config.load_shared_parameters(path)
    changed = _shared_parameter_payload()
    changed["unit_cell_eigenmode_count"] = 8
    path.write_text(json.dumps(changed), encoding="utf-8")

    with pytest.raises(RuntimeError, match="geometry, mesh, or unit-cell mode count changed"):
        parameter_config.update_center_frequency_from_cavity_p2(
            198.5,
            expected_parameters=expected,
            path=path,
        )

    assert json.loads(path.read_text(encoding="utf-8"))["center_frequency_thz"] == 198.0


def test_cavity_p2_update_preserves_layer_changes_not_used_by_unit_cell(tmp_path):
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(_shared_parameter_payload()), encoding="utf-8")
    expected = parameter_config.load_shared_parameters(path)
    changed = _shared_parameter_payload()
    changed["cavity_layers"] = 20
    changed["cladding_layers"] = 30
    path.write_text(json.dumps(changed), encoding="utf-8")

    parameter_config.update_center_frequency_from_cavity_p2(
        198.5,
        expected_parameters=expected,
        path=path,
    )

    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["cavity_layers"] == 20
    assert saved["cladding_layers"] == 30
    assert saved["center_frequency_thz"] == pytest.approx(198.5)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("cavity_layers", 0, "cavity_layers must be a positive integer"),
        ("cladding_layers", 1.5, "cladding_layers must be a positive integer"),
    ],
)
def test_shared_parameter_json_rejects_invalid_structure_layers(
    tmp_path, field, value, message
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload[field] = value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        parameter_config.load_shared_parameters(path)


@pytest.mark.parametrize("value", [None, 0, -1, 1.5, True])
def test_shared_parameter_json_rejects_invalid_unit_cell_mode_count(
    tmp_path,
    value,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["unit_cell_eigenmode_count"] = value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        ValueError,
        match="unit_cell_eigenmode_count must be a positive integer",
    ):
        parameter_config.load_shared_parameters(path)


@pytest.mark.parametrize("value", [None, 0, -1, 1.5, True])
def test_shared_parameter_json_rejects_invalid_finite_mode_count(
    tmp_path,
    value,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["finite_eigenmode_count"] = value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        ValueError,
        match="finite_eigenmode_count must be a positive integer",
    ):
        parameter_config.load_shared_parameters(path)


def test_unit_cell_solver_shift_remains_fixed_and_cavity_p2_gamma_is_strict():
    assert band.EIGENFREQUENCY_SHIFT == "c_const/1.55[um]"
    cavity = pd.DataFrame(
        {
            "band_label": ["p1", "p2", "p2"],
            "k_norm": [0.0, 0.0, 0.0],
            "match_status": ["gamma", "gamma", "gamma"],
            "re": [197.0, 198.5, 198.5],
        }
    )
    assert band.cavity_p2_gamma_frequency(cavity) == pytest.approx(198.5)


def test_unit_cell_main_updates_shared_frequency_from_cavity_p2(monkeypatch):
    cavity = pd.DataFrame(
        {
            "band_label": ["p2", "p2"],
            "k_norm": [0.0, 0.0],
            "match_status": ["gamma", "gamma"],
            "re": [198.75, 198.75],
        }
    )
    captured = {}
    monkeypatch.setattr(
        band,
        "run_band_workflow",
        lambda *_args: {"cavity": cavity},
    )

    def fake_update(frequency_thz, *, expected_parameters):
        captured["frequency_thz"] = frequency_thz
        captured["expected_parameters"] = expected_parameters
        return band.ACTIVE_PARAMETERS

    monkeypatch.setattr(band, "update_center_frequency_from_cavity_p2", fake_update)
    monkeypatch.setattr(band, "UPDATE_CENTER_FREQUENCY_AFTER_RUN", True)

    band.main()

    assert captured["frequency_thz"] == pytest.approx(198.75)
    assert captured["expected_parameters"] is band.ACTIVE_PARAMETERS


def test_spatial_hexagon_applies_optional_auto_mesh_size():
    calls = []

    class FakeMesh:
        def autoMeshSize(self, value):
            calls.append(value)

    class FakeComponent:
        def mesh(self, tag):
            assert tag == "mesh1"
            return FakeMesh()

    class FakeModel:
        def component(self, tag):
            assert tag == "comp1"
            return FakeComponent()

    config = hexagon_unit_cell_module.SimulationConfig(mesh_auto_size=9)
    hexagon_unit_cell_module._configure_auto_mesh(FakeModel(), config)

    assert calls == [9]


def test_spatial_hexagon_keeps_comsol_default_mesh_when_size_is_none():
    class FakeComponent:
        def mesh(self, _tag):
            raise AssertionError("COMSOL mesh should not be configured")

    class FakeModel:
        def component(self, _tag):
            return FakeComponent()

    config = hexagon_unit_cell_module.SimulationConfig()
    hexagon_unit_cell_module._configure_auto_mesh(FakeModel(), config)


def test_cladding_sampling_matches_cavity_sampling():
    cavity = band.case_k_points("cavity")
    cladding = band.case_k_points("cladding")

    assert cladding == cavity
    for direction, axis in (("gamma_m", "ky"), ("gamma_k", "kx")):
        values = [point[axis] for point in cladding[direction]]
        for expected in band.Q_FIT_K_MAGNITUDES:
            assert any(value == pytest.approx(expected) for value in values)


def test_make_branch_k_points_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="direction"):
        band.make_branch_k_points("gamma_x", 31)
    with pytest.raises(ValueError, match="point_count"):
        band.make_branch_k_points("gamma_m", 1)
    with pytest.raises(ValueError, match="outside"):
        band.make_branch_k_points("gamma_m", 31, dense_magnitudes=(0.6,))


def test_assign_mode_candidates_is_global_and_one_to_one():
    result = band_connector_module.assign_mode_candidates(
        previous_frequencies=[193.0, 193.1],
        candidate_frequencies=[193.08, 193.02],
        overlap_matrix=np.array([[0.80, 0.79], [0.99, 0.10]]),
        frequency_tolerance=1.0,
        minimum_overlap=0.2,
    )

    assert [row["candidate_index"] for row in result] == [1, 0]
    assert len({row["candidate_index"] for row in result}) == 2


def test_assign_mode_candidates_marks_disallowed_match():
    result = band_connector_module.assign_mode_candidates(
        previous_frequencies=[190.0],
        candidate_frequencies=[200.0],
        overlap_matrix=np.array([[0.99]]),
        frequency_tolerance=1.0,
        minimum_overlap=0.2,
    )

    assert result[0]["matched"] is False
    assert result[0]["candidate_index"] is None


def test_band_tracking_tolerance_accepts_high_overlap_six_thz_step():
    result = band_connector_module.assign_mode_candidates(
        previous_frequencies=[240.6],
        candidate_frequencies=[234.6, 246.1],
        overlap_matrix=np.array([[0.95, 0.98]]),
        frequency_tolerance=band.BAND_FREQUENCY_TOLERANCE_THZ,
        minimum_overlap=band.BAND_MINIMUM_OVERLAP,
    )

    assert band.BAND_FREQUENCY_TOLERANCE_THZ == pytest.approx(7.0)
    assert result[0]["matched"] is True
    assert result[0]["candidate_index"] == 1


def test_mode_composition_from_frames_returns_normalized_weights(monkeypatch):
    monkeypatch.setattr(
        band_connector_module,
        "fourier_subspace_energies_field",
        lambda *_args, **_kwargs: {
            "total_energy": 4.0,
            "E_dc": 0.4,
            "E_pairs": np.array([3.2, 0.4]),
            "E_nyquist": 0.0,
        },
    )
    monkeypatch.setattr(
        band_connector_module,
        "field_to_vector",
        lambda *_args, **_kwargs: np.ones(6),
    )
    monkeypatch.setattr(
        band_connector_module,
        "standard_to_fourier_basis",
        lambda *_args, **_kwargs: np.sqrt([0.10, 0.50, 0.30, 0.05, 0.03, 0.02]),
    )
    frame_re = pd.DataFrame({0: [0.0], 1: [0.0], 2: [1.0]})
    frame_im = pd.DataFrame({0: [0.0], 1: [0.0], 2: [0.0]})

    result = band_connector_module.mode_composition_from_frames(frame_re, frame_im)

    assert result["dominant_subspace"] == "p"
    assert result["dominant_mode"] == "px"
    assert result["p_weight"] == pytest.approx(0.8)
    assert result["px_weight"] == pytest.approx(0.5)
    assert sum(result[f"{name}_weight"] for name in ("s", "p", "d", "f")) == pytest.approx(1.0)


def test_mode_composition_keeps_s_and_f_mode_weights_distinct(monkeypatch):
    monkeypatch.setattr(
        band_connector_module,
        "fourier_subspace_energies_field",
        lambda *_args, **_kwargs: {
            "total_energy": 1.0,
            "E_dc": 0.1,
            "E_pairs": np.array([0.8, 0.1]),
            "E_nyquist": 0.0,
        },
    )
    monkeypatch.setattr(band_connector_module, "field_to_vector", lambda *_args, **_kwargs: np.ones(6))
    monkeypatch.setattr(
        band_connector_module,
        "standard_to_fourier_basis",
        lambda *_args, **_kwargs: np.sqrt([0.6, 0.1, 0.1, 0.05, 0.05, 0.1]),
    )
    frame = pd.DataFrame({0: [0.0], 1: [0.0], 2: [1.0]})

    result = band_connector_module.mode_composition_from_frames(frame, frame)

    assert result["dominant_subspace"] == "p"
    assert result["dominant_mode"] == "s"
    assert result["s_weight"] == pytest.approx(0.1)
    assert result["s_mode_weight"] == pytest.approx(0.6)


def test_k_point_is_complete_requires_all_valid_mode_fields(tmp_path):
    k_dir = tmp_path / "kx=0.000000_ky=0.000000"
    eigenmodes_dir = k_dir / "eigenmodes"
    eigenmodes_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "re": [194.0, 195.0],
            "is_valid": [True, False],
            "p_weight": [0.9, np.nan],
        }
    ).to_csv(k_dir / "eigenfrequencies.csv", index=False)

    assert not band.k_point_is_complete(k_dir)

    pd.DataFrame(
        {"x": [0.0], "y": [0.0], "re": [1.0], "im": [0.0]}
    ).to_parquet(eigenmodes_dir / "00_Hz_center.parquet")

    assert band.k_point_is_complete(k_dir)

    config = band.SimulationConfig(eigenmode_count=6, mesh_auto_size=5)
    assert not band.k_point_is_complete(k_dir, config)
    band.write_json_atomic(
        k_dir / band.POINT_SOLVER_CONFIG_FILENAME,
        band.point_solver_identity(config),
    )
    assert band.k_point_is_complete(k_dir, config)
    changed = band.SimulationConfig(eigenmode_count=8, mesh_auto_size=5)
    assert not band.k_point_is_complete(k_dir, changed)


@pytest.mark.parametrize(
    "frame",
    [
        pd.DataFrame({"re": [194.0]}),
        pd.DataFrame({"re": [194.0], "is_valid": [True]}),
    ],
)
def test_k_point_is_complete_rejects_incomplete_schema(tmp_path, frame):
    k_dir = tmp_path / "point"
    (k_dir / "eigenmodes").mkdir(parents=True)
    frame.to_csv(k_dir / "eigenfrequencies.csv", index=False)
    if "is_valid" in frame:
        pd.DataFrame(
            {"x": [0.0], "y": [0.0], "re": [1.0], "im": [0.0]}
        ).to_parquet(k_dir / "eigenmodes" / "00_Hz_center.parquet")

    assert not band.k_point_is_complete(k_dir)


def test_get_hole_params_matches_cladding_c2_reference():
    result = band.get_hole_params(band.CLADDING_PARAMS)
    cladding = band.ACTIVE_PARAMETERS.cladding
    b0 = cladding.b0_nm / 1000.0
    transverse_side = b0 * np.sqrt((3.0 - cladding.zeta**2) / 2.0)

    assert result.shape == (6, 4)
    assert result[:, 0] == pytest.approx(np.full(6, band.R_0 * cladding.eta))
    assert result[:, 1] == pytest.approx(np.zeros(6))
    assert result[:, 2] == pytest.approx(
        np.array(
            [
                b0 * cladding.zeta,
                transverse_side,
                transverse_side,
                b0 * cladding.zeta,
                transverse_side,
                transverse_side,
            ]
        )
    )
    assert result[:, 3] == pytest.approx(np.zeros(6))


def test_all_active_hole_converters_preserve_exact_zeta_one_geometry():
    params = parameter_config.CellParameters(
        b0_nm=243.7,
        eta=0.98,
        zeta=1.0,
    ).to_fourier_params("exact_c6")

    for module in (band, load_strip_script_module(), load_finite_script_module()):
        result = module.get_hole_params(params)
        assert result.shape == (6, 4)
        assert result[:, 0] == pytest.approx(np.full(6, module.R_0 * 0.98))
        assert result[:, 1] == pytest.approx(np.zeros(6), abs=1.0e-15)
        assert result[:, 2] == pytest.approx(np.full(6, 0.2437), abs=1.0e-12)
        assert result[:, 3] == pytest.approx(np.zeros(6), abs=1.0e-15)


def test_unit_cell_hole_converter_preserves_small_fourier_coefficients(monkeypatch):
    captured: list[np.ndarray] = []

    def capture(alpha, _basis_size):
        values = np.asarray(alpha, dtype=float)
        captured.append(values.copy())
        return values

    monkeypatch.setattr(band, "fourier_to_standard_basis", capture)
    params = {
        "name": "small_coefficients",
        "r_f0": 1.0,
        "b_square_f0": 1.0,
        "r_f1": 0.001,
        "theta_f1": -0.002,
        "b_square_f3": 0.003,
        "phi_f1": -0.004,
    }
    minimum_values = band.get_fourier_basis_min_nonzero_values(6)

    band.get_hole_params(params)

    assert captured[0][1] == pytest.approx(band.R_0 * 0.001 / minimum_values[1])
    assert captured[1][1] == pytest.approx(-0.002 / minimum_values[1])
    assert captured[2][3] == pytest.approx(
        band.B_0**2 * 0.003 / minimum_values[3]
    )
    assert captured[3][1] == pytest.approx(-0.004 / minimum_values[1])


def test_doublet_alignment_uses_exact_c6_cladding_mapping():
    optimize = load_doublet_alignment_script_module()

    params = optimize.make_cladding_params(247.2, 0.951)
    holes = band.get_hole_params(params)

    assert params["b_square_f3"] == 0.0
    assert holes[:, 2] == pytest.approx(np.full(6, 0.2472), abs=1.0e-12)


def test_doublet_alignment_summary_uses_pair_means_and_requested_cavity_targets():
    optimize = load_doublet_alignment_script_module()
    targets = {
        "cavity_p1_frequency_thz": 196.645,
        "cavity_d1_frequency_thz": 200.502,
    }
    gamma_modes = pd.DataFrame(
        {
            "re": [196.640, 196.648, 200.496, 200.504],
            "is_valid": [True, True, True, True],
            "dominant_subspace": ["p", "p", "d", "d"],
            "p_weight": [0.99, 0.99, 0.01, 0.01],
            "d_weight": [0.01, 0.01, 0.99, 0.99],
        }
    )

    summary = optimize.summarize_candidate(
        247.2, 0.951, gamma_modes, targets, "test"
    )

    assert summary["cladding_p_doublet_frequency_thz"] == pytest.approx(196.644)
    assert summary["cladding_d_doublet_frequency_thz"] == pytest.approx(200.500)
    assert summary["p_splitting_thz"] == pytest.approx(0.008)
    assert summary["d_splitting_thz"] == pytest.approx(0.008)
    assert summary["p_signed_error_thz"] == pytest.approx(-0.001)
    assert summary["d_signed_error_thz"] == pytest.approx(-0.002)
    assert summary["accepted"] is True


def test_doublet_alignment_adaptive_search_reaches_requested_precision(
    tmp_path, monkeypatch
):
    optimize = load_doublet_alignment_script_module()
    monkeypatch.setattr(optimize, "OUT_DIR", tmp_path / "optimization")
    target_b0 = 247.2
    target_eta = 0.951
    targets = {
        "cavity_p1_frequency_thz": 196.645,
        "cavity_d1_frequency_thz": 200.502,
    }

    def fake_evaluator(b0_nm, eta, _targets, stage):
        p_error = 0.10 * (b0_nm - target_b0) + 10.0 * (eta - target_eta)
        d_error = 0.05 * (b0_nm - target_b0) - 15.0 * (eta - target_eta)
        maximum = max(abs(p_error), abs(d_error))
        return {
            "stage": stage,
            "b0_nm": b0_nm,
            "eta": eta,
            "zeta": 1.0,
            "status": "accepted" if maximum <= 0.02 else "ok",
            "feasible": True,
            "accepted": maximum <= 0.02,
            "p_signed_error_thz": p_error,
            "d_signed_error_thz": d_error,
            "p_error_thz": abs(p_error),
            "d_error_thz": abs(d_error),
            "maximum_alignment_error_thz": maximum,
            "rms_alignment_error_thz": np.sqrt(0.5 * (p_error**2 + d_error**2)),
            "p_splitting_thz": 0.001,
            "d_splitting_thz": 0.001,
        }

    best, ranked = optimize.adaptive_search(targets, evaluator=fake_evaluator)

    assert best["b0_nm"] == pytest.approx(target_b0, abs=0.1)
    assert best["eta"] == pytest.approx(target_eta, abs=0.001)
    assert best["accepted"] is True
    assert ranked.iloc[0]["maximum_alignment_error_thz"] <= 0.02


def test_run_k_point_writes_fields_and_reuses_complete_result(tmp_path, monkeypatch):
    calls = []

    class FakeSimulationRun:
        def __init__(self, config):
            calls.append(("init", config))

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def build_and_run(self, a, holes, k):
            calls.append(("build", a, holes, k))

        def get_eigenfrequencies(self):
            return [(194.5, -0.001, 1.0e5)]

        def get_2d_fields(self, _mode_idx, expression, _plane):
            coords = np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.5]])
            if expression in {"ewfd.normH", "ewfd.normE"}:
                return coords, np.array([1.0, 1.0, 2.0])
            if expression == "ewfd.Hz":
                return coords[:, :2], np.array([1.0, 0.5, 0.25])
            if expression == "ewfd.Hz*(-i)":
                return coords[:, :2], np.array([0.0, 0.1, 0.2])
            raise AssertionError(expression)

    composition = {
        "s_weight": 0.0,
        "p_weight": 1.0,
        "d_weight": 0.0,
        "f_weight": 0.0,
        "px_weight": 0.7,
        "py_weight": 0.3,
        "dx_weight": 0.0,
        "dy_weight": 0.0,
        "dominant_subspace": "p",
        "dominant_mode": "px",
    }
    monkeypatch.setattr(band, "SimulationRun", FakeSimulationRun)
    monkeypatch.setattr(band, "load_mode_composition", lambda *_args: composition)
    point = band.make_branch_k_points("gamma_k", 2)[0]
    case_dir = tmp_path / "cavity"
    simulation_config = band.SimulationConfig(eigenmode_count=6, mesh_auto_size=5)

    first = band.run_k_point(
        case_dir,
        [[[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]]],
        point,
        simulation_config,
    )
    second = band.run_k_point(case_dir, [], point, simulation_config)

    assert len(first) == 1
    assert first.loc[0, "is_valid"]
    assert first.loc[0, "dominant_subspace"] == "p"
    assert second.loc[0, "re"] == pytest.approx(194.5)
    assert sum(item[0] == "build" for item in calls) == 1
    point_dir = band.k_point_dir(case_dir, point)
    assert band.k_point_is_complete(point_dir, simulation_config)
    assert json.loads(
        (point_dir / band.POINT_SOLVER_CONFIG_FILENAME).read_text(encoding="utf-8")
    ) == band.point_solver_identity(simulation_config)


def make_gamma_modes():
    return pd.DataFrame(
        {
            "re": [190.0, 191.0, 192.0, 194.0, 196.0, 193.0, 206.0, 198.0],
            "is_valid": [True, False, True, True, True, True, True, True],
            "dominant_subspace": ["s", "p", "p", "s", "p", "d", "d", "d"],
            "p_weight": [0.0, 0.9, 0.8, 0.1, 0.7, 0.1, 0.0, 0.2],
            "d_weight": [0.0, 0.1, 0.2, 0.1, 0.3, 0.9, 1.0, 0.8],
        }
    )


def test_select_gamma_bands_chooses_two_p_and_two_d_by_frequency():
    selected = band.select_gamma_bands(make_gamma_modes())

    assert selected == {"p1": 2, "p2": 4, "d1": 5, "d2": 7}


@pytest.mark.parametrize("missing_subspace", ["p", "d"])
def test_select_gamma_bands_reports_diagnostics_when_modes_missing(missing_subspace):
    modes = make_gamma_modes()
    modes.loc[modes["dominant_subspace"] == missing_subspace, "dominant_subspace"] = "s"

    with pytest.raises(ValueError) as exc_info:
        band.select_gamma_bands(modes)

    message = str(exc_info.value)
    assert "Gamma" in message
    assert "dominant_subspace" in message
    assert "p_weight" in message
    assert "d_weight" in message


def make_tracked_modes(frequencies):
    count = len(frequencies)
    labels = ["p", "p", "d", "d"][:count]
    return pd.DataFrame(
        {
            "re": frequencies,
            "im": np.full(count, -0.001),
            "q": np.linspace(1.0e5, 4.0e5, count),
            "is_valid": np.ones(count, dtype=bool),
            "dominant_subspace": labels,
            "dominant_mode": ["px", "py", "dx", "dy"][:count],
            "s_weight": np.zeros(count),
            "p_weight": [0.9, 0.8, 0.1, 0.2][:count],
            "d_weight": [0.1, 0.2, 0.9, 0.8][:count],
            "f_weight": np.zeros(count),
            "px_weight": [0.7, 0.2, 0.0, 0.0][:count],
            "py_weight": [0.3, 0.8, 0.0, 0.0][:count],
            "dx_weight": [0.0, 0.0, 0.7, 0.2][:count],
            "dy_weight": [0.0, 0.0, 0.3, 0.8][:count],
        }
    )


def test_track_selected_branch_preserves_labels_and_marks_unmatched(tmp_path, monkeypatch):
    points = band.make_branch_k_points("gamma_k", 2)
    case_dir = tmp_path / "cavity"
    gamma_dir = band.k_point_dir(case_dir, points[0])
    next_dir = band.k_point_dir(case_dir, points[1])
    gamma_dir.mkdir(parents=True)
    next_dir.mkdir(parents=True)
    make_tracked_modes([192.0, 193.0, 197.0, 198.0]).to_csv(
        gamma_dir / "eigenfrequencies.csv", index=False
    )
    make_tracked_modes([192.1, 193.1, 197.1, 210.0]).to_csv(
        next_dir / "eigenfrequencies.csv", index=False
    )

    monkeypatch.setattr(band, "load_hz_center", lambda k_dir, mode_idx: (str(k_dir), mode_idx))
    monkeypatch.setattr(
        band,
        "field_overlap",
        lambda previous, candidate: 0.95 if previous[1] == candidate[1] else 0.05,
    )

    result = band.track_selected_branch(
        case_dir,
        points,
        {"p1": 0, "p2": 1, "d1": 2, "d2": 3},
    )

    assert len(result) == 8
    assert set(result["band_label"]) == {"p1", "p2", "d1", "d2"}
    second_point = result[result["point_index"] == 1].set_index("band_label")
    assert second_point.loc["p1", "mode_idx"] == 0
    assert second_point.loc["p2", "mode_idx"] == 1
    assert second_point.loc["d1", "mode_idx"] == 2
    assert second_point.loc["d2", "match_status"] == "unmatched"
    assert pd.isna(second_point.loc["d2", "mode_idx"])


def test_track_selected_branch_keeps_reporting_after_all_bands_are_lost(tmp_path, monkeypatch):
    points = band.make_branch_k_points("gamma_k", 3)
    case_dir = tmp_path / "cavity"
    for point_index, point in enumerate(points):
        point_dir = band.k_point_dir(case_dir, point)
        point_dir.mkdir(parents=True)
        frequencies = [192.0, 193.0, 197.0, 198.0] if point_index == 0 else [210.0] * 4
        make_tracked_modes(frequencies).to_csv(point_dir / "eigenfrequencies.csv", index=False)
    monkeypatch.setattr(band, "load_hz_center", lambda k_dir, mode_idx: (str(k_dir), mode_idx))
    monkeypatch.setattr(band, "field_overlap", lambda *_args: 0.99)

    result = band.track_selected_branch(
        case_dir,
        points,
        {"p1": 0, "p2": 1, "d1": 2, "d2": 3},
    )

    assert len(result) == 12
    assert set(result[result["point_index"] == 2]["match_status"]) == {"unmatched"}


def test_fit_q_power_law_recovers_exponent_and_excludes_invalid_rows():
    k_values = np.asarray(band.Q_FIT_K_MAGNITUDES)
    q_values = 7.5e3 * k_values ** -2.25
    rows = pd.DataFrame(
        {
            "k_norm": np.concatenate([[0.0], k_values, [0.01, 0.02]]),
            "q": np.concatenate([[1.0e12], q_values, [np.nan, -1.0]]),
        }
    )

    fit = band.fit_q_power_law(rows)

    assert fit["status"] == "ok"
    assert fit["alpha"] == pytest.approx(2.25, rel=1e-6)
    assert fit["point_count"] == len(band.Q_FIT_K_MAGNITUDES)
    assert fit["r_squared"] == pytest.approx(1.0)


def test_fit_q_power_law_reports_insufficient_points():
    rows = pd.DataFrame({"k_norm": [0.0, 0.002, 0.004], "q": [1e12, 1e8, np.nan]})

    fit = band.fit_q_power_law(rows)

    assert fit["status"] == "insufficient_points"
    assert fit["point_count"] == 1


def test_p_x_fraction_uses_project_fourier_definition():
    result = band.p_x_fraction(
        np.array([1.0, 0.0, 0.0, np.nan]),
        np.array([0.0, 1.0, 0.0, 1.0]),
    )

    np.testing.assert_allclose(result[:3], np.array([1.0, 0.0, 0.5]))
    assert np.isnan(result[3])


def make_selected_band_rows(cell):
    rows = []
    base_frequency = {"p1": 192.0, "p2": 193.0, "d1": 197.0, "d2": 198.0}
    p_weight = {"p1": 0.9, "p2": 0.8, "d1": 0.1, "d2": 0.2}
    for direction, sign in (("gamma_m", -1.0), ("gamma_k", 1.0)):
        for path_fraction, k_norm in ((0.0, 0.0), (0.5, 0.02), (1.0, 0.05)):
            for label in ("p1", "p2", "d1", "d2"):
                rows.append(
                    {
                        "cell": cell,
                        "direction": direction,
                        "path_fraction": path_fraction,
                        "path_coordinate": sign * path_fraction,
                        "k_norm": k_norm,
                        "band_label": label,
                        "re": base_frequency[label] + sign * 0.3 * path_fraction,
                        "q": 1e7 * max(k_norm, 0.002) ** (-2.0),
                        "p_weight": p_weight[label],
                        "px_weight": 0.7 if label == "p1" else 0.2,
                        "py_weight": 0.3 if label == "p1" else 0.8,
                        "dx_weight": 0.3 if label == "d1" else 0.8,
                        "dy_weight": 0.7 if label == "d1" else 0.2,
                        "match_status": "matched" if path_fraction else "gamma",
                    }
                )
    return pd.DataFrame(rows)


def test_p_band_color_rejects_missing_gamma_composition():
    frame = make_selected_band_rows("cavity")
    frame = frame[
        ~(
            (frame["band_label"] == "p1")
            & np.isclose(frame["path_coordinate"], 0.0)
        )
    ]

    with pytest.raises(ValueError, match="Missing Gamma composition for p1"):
        band.p_band_color(frame, "p1")


def test_p_band_color_rejects_inconsistent_directional_gamma_composition():
    frame = make_selected_band_rows("cavity")
    gamma_k_p1 = (
        (frame["direction"] == "gamma_k")
        & (frame["band_label"] == "p1")
        & np.isclose(frame["path_coordinate"], 0.0)
    )
    frame.loc[gamma_k_p1, ["px_weight", "py_weight"]] = [0.0, 1.0]

    with pytest.raises(ValueError, match="Inconsistent Gamma composition for p1"):
        band.p_band_color(frame, "p1")


def test_plot_band_comparison_centers_shared_limits_on_cavity_p2_gamma(
    tmp_path, monkeypatch
):
    original_close = band.plt.close
    monkeypatch.setattr(band.plt, "close", lambda *_args, **_kwargs: None)
    output = tmp_path / "unit_cell_band_comparison.png"
    cavity = make_selected_band_rows("cavity")

    band.plot_band_comparison(
        cavity,
        make_selected_band_rows("cladding"),
        output,
    )

    assert output.exists()
    assert output.stat().st_size > 0
    figure = band.plt.figure(band.plt.get_fignums()[-1])
    assert len(figure.axes) == 2
    plot_axes = figure.axes
    cavity_p2_gamma = cavity[
        (cavity["band_label"] == "p2")
        & np.isclose(cavity["path_coordinate"], 0.0)
    ]
    center_frequency = float(cavity_p2_gamma.iloc[0]["re"])
    half_span = 0.5 * (band.FREQUENCY_MAX_THZ - band.FREQUENCY_MIN_THZ)
    expected_limits = (center_frequency - half_span, center_frequency + half_span)
    assert all(axis.get_ylim() == pytest.approx(expected_limits) for axis in plot_axes)
    assert [tick.get_text() for tick in plot_axes[0].get_xticklabels()] == ["M", "Γ", "K"]
    assert all(not axis.texts for axis in plot_axes)
    assert len(figure.legends) == 1
    assert [text.get_text() for text in figure.legends[0].get_texts()] == [
        r"$p_y$",
        r"$p_x$",
        r"$d_{xy}$",
        r"$d_{x^2+y^2}$",
    ]
    legend_colors = [
        to_rgba(handle.get_color()) for handle in figure.legends[0].legend_handles[:2]
    ]
    assert legend_colors == pytest.approx(
        [to_rgba("#dc2626"), to_rgba("#2563eb")]
    )
    line_collections = [
        collection
        for collection in plot_axes[0].collections
        if isinstance(collection, LineCollection)
    ]
    assert len(line_collections) == 4
    assert len(plot_axes[0].collections) == 4
    p_collection = line_collections[0]
    d_collection = line_collections[2]
    assert p_collection.cmap(0.0) == pytest.approx(to_rgba("#2563eb"), abs=0.01)
    assert p_collection.cmap(0.5) == pytest.approx(to_rgba("white"), abs=0.01)
    assert p_collection.cmap(1.0) == pytest.approx(to_rgba("#dc2626"), abs=0.01)
    np.testing.assert_allclose(p_collection.get_array(), np.full(4, 0.7))
    assert d_collection.cmap(0.0) == pytest.approx(to_rgba("#b0b0b0"), abs=0.01)
    assert d_collection.cmap(1.0) == pytest.approx(to_rgba("black"), abs=0.01)
    np.testing.assert_allclose(d_collection.get_array(), np.full(4, 0.3))
    original_close(figure)


def test_plot_band_comparison_rejects_inconsistent_cavity_p2_gamma(tmp_path):
    cavity = make_selected_band_rows("cavity")
    cavity.loc[cavity["band_label"] == "p1", ["px_weight", "py_weight"]] = [
        0.0,
        1.0,
    ]
    cavity.loc[cavity["band_label"] == "p2", ["px_weight", "py_weight"]] = [
        1.0,
        0.0,
    ]
    gamma_p2 = (
        (cavity["band_label"] == "p2")
        & np.isclose(cavity["path_coordinate"], 0.0)
    )
    cavity.loc[cavity.index[gamma_p2][-1], "re"] += 0.1

    with pytest.raises(ValueError, match="cavity p2 Gamma"):
        band.plot_band_comparison(
            cavity,
            make_selected_band_rows("cladding"),
            tmp_path / "unit_cell_band_comparison.png",
        )


def test_plot_p_band_q_comparison_combines_both_cells_on_shared_axes(
    tmp_path, monkeypatch
):
    original_close = band.plt.close
    monkeypatch.setattr(band.plt, "close", lambda *_args, **_kwargs: None)
    output = tmp_path / "p_bands_q_vs_k_comparison.png"
    fits = pd.DataFrame(
        [
            {
                "direction": direction,
                "band_label": label,
                "status": "ok",
                "alpha": 2.0,
                "intercept": np.log(1e7),
                "r_squared": 0.99,
                "k_min": 0.002,
                "k_max": 0.05,
                "point_count": 8,
            }
            for direction in ("gamma_m", "gamma_k")
            for label in ("p1", "p2")
        ]
    )

    cavity = make_selected_band_rows("cavity")
    cavity.loc[cavity["band_label"] == "p1", ["px_weight", "py_weight"]] = [
        0.0,
        1.0,
    ]
    cavity.loc[cavity["band_label"] == "p2", ["px_weight", "py_weight"]] = [
        1.0,
        0.0,
    ]
    outside = cavity[cavity["direction"] == "gamma_k"].iloc[[0]].copy()
    outside["k_norm"] = 0.3
    outside["path_fraction"] = 0.6
    outside["path_coordinate"] = 0.6
    cavity = pd.concat([cavity, outside], ignore_index=True)

    cladding = make_selected_band_rows("cladding")
    cladding.loc[cladding["band_label"] == "p1", ["px_weight", "py_weight"]] = [
        1.0,
        0.0,
    ]
    cladding.loc[cladding["band_label"] == "p2", ["px_weight", "py_weight"]] = [
        0.0,
        1.0,
    ]
    cladding["q"] *= 0.1
    band.plot_p_band_q_comparison(cavity, cladding, fits, fits, output)

    assert output.exists()
    assert output.stat().st_size > 0
    figure = band.plt.figure(band.plt.get_fignums()[-1])
    assert len(figure.axes) == 2
    assert [axis.get_title() for axis in figure.axes] == ["Cavity cell", "Cladding cell"]
    assert all(axis.get_xlim() == pytest.approx((-0.2, 0.2)) for axis in figure.axes)
    assert figure.axes[0].get_ylim() == pytest.approx(figure.axes[1].get_ylim())
    assert all(
        axis.get_xlabel() == "Signed |k|/G  (Γ–M < 0, Γ–K > 0)"
        for axis in figure.axes
    )
    assert figure.axes[0].get_ylabel() == "log10(Q)"
    expected_colors = (
        ["#2563eb", "#dc2626", "#2563eb", "#dc2626"],
        ["#dc2626", "#2563eb", "#dc2626", "#2563eb"],
    )
    for axis, expected in zip(figure.axes, expected_colors):
        plotted_x = np.concatenate([line.get_xdata() for line in axis.lines])
        assert np.all(np.abs(plotted_x) <= 0.2 + 1e-12)
        assert np.any(plotted_x < 0.0)
        assert np.any(plotted_x > 0.0)
        raw_lines = [line for line in axis.lines if str(line.get_label()).startswith("p")]
        fit_lines = [line for line in axis.lines if line.get_linestyle() == ":"]
        assert [to_rgba(line.get_color()) for line in raw_lines] == pytest.approx(
            [to_rgba(color) for color in expected]
        )
        assert [to_rgba(line.get_color()) for line in fit_lines] == pytest.approx(
            [to_rgba(color) for color in expected]
        )
        annotation = " ".join(text.get_text() for text in axis.texts)
        assert "Γ–M p1: α" in annotation
        assert "Γ–K p1: α" in annotation
    original_close(figure)


def test_main_orchestrates_outputs_and_reuses_complete_points(tmp_path, monkeypatch):
    monkeypatch.setattr(band, "OUT_DIR", tmp_path / "unit_cell_band")
    monkeypatch.setattr(band, "BAND_POINTS_PER_ARM", 2)
    monkeypatch.setattr(band, "Q_FIT_K_MAGNITUDES", (0.002, 0.004, 0.006))
    monkeypatch.setattr(band, "UPDATE_CENTER_FREQUENCY_AFTER_RUN", False)
    monkeypatch.setattr(
        band,
        "visualize_hexagon_design",
        lambda *_args, filename, **_kwargs: Path(filename).write_bytes(b"geometry"),
    )
    monkeypatch.setattr(
        band,
        "load_hz_center",
        lambda k_dir, mode_idx: (str(k_dir), int(mode_idx)),
    )
    monkeypatch.setattr(
        band,
        "field_overlap",
        lambda previous, candidate: 0.99 if previous[1] == candidate[1] else 0.01,
    )
    solve_calls = []

    class FakeReusableSimulationRun:
        def __init__(self, config):
            self.execution_stats = {
                "execution_engine": "fake_reused_model",
                "model_reuse_enabled": True,
                "model_build_count": 1,
                "mesh_build_count": 1,
                "solve_count": 0,
                "mesh_auto_size": config.mesh_auto_size,
                "current_k": None,
            }

        def clear(self):
            return None

    monkeypatch.setattr(band, "ReusableSimulationRun", FakeReusableSimulationRun)

    def fake_run_k_point(
        case_dir,
        _holes,
        point,
        _simulation_config,
        sim_run=None,
    ):
        output_dir = band.k_point_dir(case_dir, point)
        if band.k_point_is_complete(output_dir):
            return pd.read_csv(output_dir / "eigenfrequencies.csv")
        assert sim_run is not None
        solve_calls.append((case_dir.name, point["kx_str"], point["ky_str"]))
        output_dir.joinpath("eigenmodes").mkdir(parents=True, exist_ok=True)
        k_norm = max(float(point["k_norm"]), 0.002)
        frame = make_tracked_modes([192.0, 193.0, 197.0, 198.0])
        if case_dir.name == "cavity" and point["ky_str"] == "0.500000":
            frame.loc[0, "re"] = 220.0
        frame["q"] = [1e4 * k_norm**-2.0, 2e4 * k_norm**-2.2, 2e4, 1e4]
        frame.to_csv(output_dir / "eigenfrequencies.csv", index=False)
        for mode_idx in frame.index:
            pd.DataFrame(
                {
                    "x": [0.0, 0.1, 0.0],
                    "y": [0.0, 0.0, 0.1],
                    "re": [1.0, 0.5, 0.25],
                    "im": [0.0, 0.1, 0.2],
                }
            ).to_parquet(output_dir / "eigenmodes" / f"{mode_idx:02d}_Hz_center.parquet")
        return frame

    monkeypatch.setattr(band, "run_k_point", fake_run_k_point)

    band.main()
    first_solve_count = len(solve_calls)
    band.main()

    expected = [
        "99_config/config.json",
        "99_config/cavity/config.json",
        "99_config/cladding/config.json",
        "10_overview/unit_cell_band_comparison.png",
        "10_overview/p_bands_q_vs_k_comparison.png",
        "10_overview/cavity/selected_bands.csv",
        "10_overview/cavity/mode_composition.csv",
        "10_overview/cavity/q_power_law_fits.csv",
        "10_overview/cladding/selected_bands.csv",
        "10_overview/cladding/mode_composition.csv",
        "10_overview/cladding/q_power_law_fits.csv",
    ]
    for relative_path in expected:
        path = band.OUT_DIR / relative_path
        assert path.exists(), relative_path
        assert path.stat().st_size > 0, relative_path
    config = json.loads(
        (band.OUT_DIR / "99_config" / "config.json").read_text(encoding="utf-8")
    )
    assert config["simulation"]["mesh_auto_size"] == band.MESH_AUTO_SIZE
    assert config["simulation"]["eigenmode_count"] == band.EIGENMODE_COUNT
    cavity_config = json.loads(
        (band.OUT_DIR / "99_config" / "cavity" / "config.json").read_text(
            encoding="utf-8"
        )
    )
    assert cavity_config["tracking"]["status"] == "partial"
    assert cavity_config["tracking"]["unmatched_row_count"] == 1
    assert cavity_config["tracking"]["unmatched"][0]["band_label"] == "p1"
    assert first_solve_count > 0
    assert len(solve_calls) == first_solve_count
    assert (band.OUT_DIR / "01_results" / "cavity" / "k_points").is_dir()
    assert (band.OUT_DIR / "01_results" / "cladding" / "k_points").is_dir()
    assert not (band.OUT_DIR / "cavity").exists()
    assert not (band.OUT_DIR / "cladding").exists()
    assert not (band.OUT_DIR / "config.json").exists()


def test_migrate_legacy_unit_cell_layout_preserves_files_and_bytes(tmp_path):
    root = tmp_path / "legacy"
    root.mkdir()
    root_files = {
        "config.json": b"config",
        "unit_cell_band_comparison.png": b"bands",
        "p_bands_q_vs_k_comparison.png": b"q",
    }
    for name, payload in root_files.items():
        (root / name).write_bytes(payload)
    (root / "80_logs").mkdir()
    (root / "80_logs" / "run_stdout.log").write_bytes(b"log")
    for case_name in ("cavity", "cladding"):
        case_dir = root / case_name
        (case_dir / "k_points" / "gamma").mkdir(parents=True)
        (case_dir / "k_points" / "gamma" / "eigenfrequencies.csv").write_bytes(
            b"re\n198\n"
        )
        (case_dir / "config.json").write_bytes(b"case-config")
        for filename in band.LEGACY_CASE_OVERVIEW_FILES:
            (case_dir / filename).write_bytes(filename.encode("utf-8"))

    before = [path for path in root.rglob("*") if path.is_file()]
    before_bytes = sum(path.stat().st_size for path in before)
    result = band.migrate_legacy_unit_cell_layout(root)
    after = [path for path in root.rglob("*") if path.is_file()]

    assert result == {
        "files_before": len(before),
        "files_after": len(before),
        "bytes_before": before_bytes,
        "bytes_after": before_bytes,
    }
    assert (root / "01_results" / "cavity" / "k_points").is_dir()
    assert (root / "01_results" / "cladding" / "k_points").is_dir()
    assert (root / "10_overview" / "cavity" / "selected_bands.csv").is_file()
    assert (root / "10_overview" / "cladding" / "selected_bands.csv").is_file()
    assert (root / "99_config" / "cavity" / "config.json").is_file()
    assert (root / "99_config" / "cladding" / "config.json").is_file()
    assert (root / "99_config" / "config.json").is_file()
    assert (root / "80_logs" / "run_stdout.log").is_file()
    assert not (root / "cavity").exists()
    assert not (root / "cladding").exists()


def test_run_band_workflow_accepts_explicit_output_and_params(tmp_path, monkeypatch):
    cavity_params = {"name": "cavity", "r_f0": 0.96, "b_square_f0": 1.0}
    cladding_params = {"name": "cladding", "r_f0": 0.95, "b_square_f0": 1.0}
    calls = []

    def fake_run_case(case_name, params, output_root=None):
        calls.append((case_name, params, Path(output_root)))
        return pd.DataFrame({"cell": [case_name]})

    monkeypatch.setattr(band, "run_case", fake_run_case)
    monkeypatch.setattr(band, "_build_q_fits", lambda _frame: pd.DataFrame())
    monkeypatch.setattr(band, "plot_band_comparison", lambda *_args: None)
    monkeypatch.setattr(band, "plot_p_band_q_comparison", lambda *_args: None)

    output_dir = tmp_path / "full"
    result = band.run_band_workflow(output_dir, cavity_params, cladding_params)

    assert calls == [
        ("cavity", cavity_params, output_dir),
        ("cladding", cladding_params, output_dir),
    ]
    assert result["output_dir"] == output_dir
    assert result["cavity"].iloc[0]["cell"] == "cavity"
    assert result["cladding"].iloc[0]["cell"] == "cladding"


def test_finite_cavity_preset_resolves_current_geometry_and_frequency():
    finite = load_finite_script_module()

    assert finite.BULK_PARAMS == finite.ACTIVE_PARAMETERS.cavity.to_fourier_params(
        "cavity_p_bic"
    )
    assert finite.CLADDING_PARAMS == (
        finite.ACTIVE_PARAMETERS.cladding.to_fourier_params(
            "cladding_bandgap_alignment"
        )
    )
    assert finite.EIGENFREQUENCY_SHIFT == finite.ACTIVE_PARAMETERS.eigenfrequency_shift
    assert finite.CLADDING_SHIFT_FACTORS == list(
        finite.ACTIVE_PARAMETERS.cladding_shift_factors
    )
    assert finite.CLADDING_Y_OVER_X_SHIFT_RATIO == (
        finite.ACTIVE_PARAMETERS.cladding_y_over_x_shift_ratio
    )
    assert finite.CLADDING_SHIFT_INPUT_KIND == (
        finite.ACTIVE_PARAMETERS.cladding_shift_input_kind
    )
    assert finite.CLADDING_X_SHIFT_FACTOR == finite.CLADDING_SHIFT_FACTORS[0]
    assert finite.CLADDING_Y_SHIFT_FACTOR == pytest.approx(
        finite.CLADDING_SHIFT_FACTORS[0]
        * finite.CLADDING_Y_OVER_X_SHIFT_RATIO
    )


@pytest.mark.parametrize(
    ("preset_name", "presets", "message"),
    [
        ("missing", {}, "Unknown finite-cavity preset"),
        (
            "bad",
            {
                "bad": {
                    "cavity": {"b0_nm": 235.0, "eta": 0.96},
                    "cladding": {"b0_nm": 231.8, "eta": 0.96, "zeta": 0.83},
                    "px_gamma_frequency_thz": 195.9,
                    "cladding_inward_shift_factors": [0.0],
                }
            },
            "cavity is missing required fields",
        ),
        (
            "bad",
            {
                "bad": {
                    "cavity": {"b0_nm": 0.0, "eta": 0.96, "zeta": 1.156},
                    "cladding": {"b0_nm": 231.8, "eta": 0.96, "zeta": 0.83},
                    "px_gamma_frequency_thz": 195.9,
                    "cladding_inward_shift_factors": [0.0],
                }
            },
            "cavity b0_nm must be positive and finite",
        ),
        (
            "bad",
            {
                "bad": {
                    "cavity": {"b0_nm": 235.0, "eta": 0.0, "zeta": 1.156},
                    "cladding": {"b0_nm": 231.8, "eta": 0.96, "zeta": 0.83},
                    "px_gamma_frequency_thz": 195.9,
                    "cladding_inward_shift_factors": [0.0],
                }
            },
            "cavity eta must be positive and finite",
        ),
        (
            "bad",
            {
                "bad": {
                    "cavity": {"b0_nm": 235.0, "eta": 0.96, "zeta": 1.156},
                    "cladding": {"b0_nm": 231.8, "eta": 0.96, "zeta": 0.83},
                    "px_gamma_frequency_thz": 195.9,
                    "cladding_inward_shift_factors": [],
                }
            },
            "cladding_shift_factors must be non-empty",
        ),
    ],
)
def test_finite_cavity_preset_rejects_invalid_records(preset_name, presets, message):
    finite = load_finite_script_module()

    with pytest.raises(ValueError, match=message):
        finite.resolve_finite_cavity_preset(preset_name, presets)


def test_finite_cavity_preset_controls_case_name_and_metadata():
    finite = load_finite_script_module()

    assert finite.OUTPUT_ROOT == parameter_config.FINITE_CAVITY_OUTPUT_ROOT
    assert finite.OUT_DIR.name == finite.ACTIVE_PARAMETERS.structure_series_label(
        "finite"
    )
    assert finite.case_stem(0.0) == "shift0.000"
    assert finite.case_stem(0.05, 0.12) == "shiftx0.050_shifty0.120"
    assert finite.BULK_RADIUS == finite.ACTIVE_PARAMETERS.cavity_layers - 1
    assert finite.CLADDING_LAYERS == finite.ACTIVE_PARAMETERS.cladding_layers
    assert finite.MESH_AUTO_SIZE == finite.ACTIVE_PARAMETERS.mesh_size
    assert finite.USE_MANUAL_MESH_CONSTRUCTION is True
    assert finite.RESUME_FROM_EXISTING_MPH is False
    metadata = finite.shared_case_metadata(
        {}, np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    )
    assert metadata["finite_cavity_preset"] == finite.FINITE_CAVITY_PRESET
    config = metadata["finite_cavity_preset_config"]
    assert config["cavity"] == finite.ACTIVE_PARAMETERS.cavity.to_dict()
    assert config["cladding"] == finite.ACTIVE_PARAMETERS.cladding.to_dict()
    assert config["bulk_params"] == finite.BULK_PARAMS
    assert config["cladding_params"] == finite.CLADDING_PARAMS
    assert config["eigenfrequency_shift"] == finite.ACTIVE_PARAMETERS.eigenfrequency_shift
    assert config["cladding_shift_factors"] == list(
        finite.ACTIVE_PARAMETERS.cladding_shift_factors
    )
    assert config["finite_cladding_shift_input"] == (
        finite.ACTIVE_PARAMETERS.finite_cladding_shift_input
    )
    assert config["finite_cladding_shift_factor_pairs"] == [
        list(pair)
        for pair in finite.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs
    ]


def test_finite_cavity_preset_supports_explicit_single_axis_zero_pair():
    finite = load_finite_script_module()
    resolved = finite.resolve_finite_cavity_preset(
        "explicit",
        {
            "explicit": {
                "cavity": {"b0_nm": 245.0, "eta": 0.96, "zeta": 1.156},
                "cladding": {"b0_nm": 240.1, "eta": 0.967, "zeta": 1.0},
                "px_gamma_frequency_thz": 198.4,
                "cladding_x_shift_factor": 0.0,
                "cladding_y_shift_factor": 0.08,
                "cladding_shift_geometry": {"kind": "ellipse"},
                "cladding_shift_profile": {"kind": "uniform"},
            }
        },
    )

    assert resolved["finite_cladding_shift_input"] == {
        "kind": "xy",
        "x_factor": 0.0,
        "y_factor": 0.08,
    }
    assert resolved["finite_cladding_shift_factor_pairs"] == [(0.0, 0.08)]
