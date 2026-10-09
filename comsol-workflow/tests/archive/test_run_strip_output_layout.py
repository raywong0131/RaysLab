import json
import sys
import types
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

from scripts.run_main import run_strip_1d
from scripts.run_main.parameter_config import STRIP_1D_OUTPUT_ROOT


def test_strip_output_root_is_stable():
    assert run_strip_1d.OUTPUT_ROOT == STRIP_1D_OUTPUT_ROOT
    assert run_strip_1d.OUT_DIR.name == (
        run_strip_1d.ACTIVE_PARAMETERS.structure_series_label("strip1d")
    )


def test_strip_cut_directory_and_shift_case_are_not_redundant(monkeypatch):
    monkeypatch.setattr(run_strip_1d, "STRIP_CUT_ANGLE", 0.0)
    assert run_strip_1d.strip_cut_dirname() == "cut_00"
    assert run_strip_1d.strip_case_stem(0.03) == "shift0.030"

    monkeypatch.setattr(run_strip_1d, "STRIP_CUT_ANGLE", 60.0)
    assert run_strip_1d.strip_cut_dirname() == "cut_60"
    assert run_strip_1d.strip_case_stem(0.03) == "shift0.030"


def test_strip_artifacts_are_written_directly_in_parameterized_case_dir():
    case_dir = Path(
        "strip1d_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9"
    ) / "cut_60" / "shift0.030"
    assert run_strip_1d.strip_case_dir(case_dir) == case_dir


def _write_staged_strip_outputs(case_dir: Path) -> dict[str, Path]:
    paths = run_strip_1d.prepare_strip_case_output(case_dir)
    export_dir = paths["export_dir"]
    export_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {"re": 198.0, "im": 0.1, "q": 1000.0, "is_valid": True, "mode_idx": 0},
            {"re": 199.0, "im": 0.2, "q": 500.0, "is_valid": False, "mode_idx": 1},
        ]
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    for mode_idx in (0, 1):
        for name in (
            "Hz_center.parquet",
            "Hz_Re_2d.png",
            "Hz_Im_2d.png",
            "Wem_2d.png",
            "cutline_Wem.png",
        ):
            (export_dir / f"{mode_idx:02d}_{name}").write_bytes(
                f"mode={mode_idx};name={name}".encode("utf-8")
            )

    fourier_dir = paths["staging_dir"] / "strip_bulk_fourier_hz"
    for mode_idx in (0, 1):
        mode_fourier = fourier_dir / f"mode_{mode_idx:02d}"
        mode_fourier.mkdir(parents=True)
        for name in (
            "p_subspace_mode_decomposition.csv",
            "strip_bulk_fourier_hz.npz",
            "top_k_peaks.csv",
        ):
            (mode_fourier / name).write_text(name, encoding="utf-8")
    pd.DataFrame(
        [
            {
                "mode_idx": mode_idx,
                "frequency": 198.0 + mode_idx,
                "is_valid": mode_idx == 0,
                "output_dir": str(fourier_dir / f"mode_{mode_idx:02d}"),
            }
            for mode_idx in (0, 1)
        ]
    ).to_csv(fourier_dir / "strip_bulk_fourier_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "mode_idx": 0,
                "frequency": 198.0,
                "gamma_subspace_p": 0.9001,
                "mode_score": 0.9,
            },
            {
                "mode_idx": 1,
                "frequency": 199.0,
                "gamma_subspace_p": 0.9,
                "mode_score": 0.8,
            },
        ]
    ).to_csv(paths["staging_dir"] / "mode_scores.csv", index=False)
    (paths["staging_dir"] / "objective.json").write_text(
        json.dumps({"score": 0.9, "best_mode": {"mode_idx": 0}}),
        encoding="utf-8",
    )
    return paths


def test_strip_finalization_keeps_only_strict_gamma_valid_modes(tmp_path):
    case_dir = tmp_path / "cut_60" / "shift0.030"
    paths = _write_staged_strip_outputs(case_dir)

    run_strip_1d.finalize_strip_case_output(case_dir)

    assert not paths["staging_dir"].exists()
    assert not list(case_dir.glob("*.*"))
    assert (case_dir / "10_overview/eigenfrequencies.csv").is_file()
    assert (case_dir / "10_overview/mode_scores.csv").is_file()
    assert (case_dir / "99_config/objective.json").is_file()
    valid_mode = case_dir / "01_results/mode0"
    invalid_mode = case_dir / "01_results/mode1"
    assert not list(valid_mode.glob("*.*"))
    assert (valid_mode / "10_overview/eigenfrequency.csv").is_file()
    assert (valid_mode / "10_overview/mode_score.csv").is_file()
    assert (valid_mode / "11_simulation_exports/Hz_center.parquet").is_file()
    assert (valid_mode / "11_simulation_exports/Hz_Re_2d.png").is_file()
    assert (valid_mode / "11_simulation_exports/Hz_Im_2d.png").is_file()
    assert (valid_mode / "11_simulation_exports/Wem_2d.png").is_file()
    assert (valid_mode / "11_simulation_exports/cutline_Wem.png").is_file()
    assert (valid_mode / "13_strip_bulk_fourier_Hz/top_k_peaks.csv").is_file()
    assert not invalid_mode.exists()

    summary = pd.read_csv(
        case_dir / "10_overview/strip_bulk_fourier_summary.csv"
    )
    assert Path(summary.loc[0, "output_dir"]) == (
        valid_mode / "13_strip_bulk_fourier_Hz"
    )
    assert pd.isna(summary.loc[1, "output_dir"])
    assert summary["is_valid"].tolist() == [True, False]
    eigen = pd.read_csv(case_dir / "10_overview/eigenfrequencies.csv")
    assert eigen["geometry_is_valid"].tolist() == [True, False]
    assert eigen["is_valid"].tolist() == [True, False]


def test_strip_finalization_rejects_unknown_staging_before_moving(tmp_path):
    case_dir = tmp_path / "cut_60" / "shift0.030"
    paths = _write_staged_strip_outputs(case_dir)
    unknown = paths["staging_dir"] / "unexpected.bin"
    unknown.write_bytes(b"unexpected")

    with pytest.raises(ValueError, match="Unknown staged strip output"):
        run_strip_1d.finalize_strip_case_output(case_dir)

    assert unknown.is_file()
    assert (paths["export_dir"] / "eigenfrequencies.csv").is_file()
    assert not (case_dir / "10_overview/eigenfrequencies.csv").exists()


def test_strip_calculation_defaults_to_noninteractive_no_diagnostic_plots():
    assert matplotlib.get_backend().lower() == "agg"
    assert run_strip_1d.RUN_DIAGNOSTIC_PLOTS is False


def test_wem_cutline_plot_contract(tmp_path, monkeypatch):
    x = np.array([-1.0, 0.0, 1.0])
    y = np.array([-2.0, 0.0, 2.0])
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack([xx.ravel(), yy.ravel()])
    wem = (1.0 + yy.ravel() ** 2) * (1.0 - 0.1 * xx.ravel() ** 2)
    captured = {}
    original_subplots = run_strip_1d.plt.subplots

    def capture_subplots(*args, **kwargs):
        fig, ax = original_subplots(*args, **kwargs)
        captured["ax"] = ax
        return fig, ax

    monkeypatch.setattr(run_strip_1d.plt, "subplots", capture_subplots)
    output = tmp_path / "cutline_Wem.png"
    cut_y, normalized = run_strip_1d.save_wem_cutline_plot(
        output,
        points,
        wem,
    )

    assert output.is_file()
    assert len(cut_y) == run_strip_1d.WEM_CUTLINE_SAMPLE_COUNT
    assert np.isclose(np.max(np.abs(normalized)), 1.0)
    ax = captured["ax"]
    assert ax.get_title() == "Wem cutline at x = 0"
    assert ax.get_xlabel() == "y (µm)"
    assert ax.get_ylabel() == "Normalized Wem"
    assert len(ax.lines) == 1


def test_run_case_writes_current_run_summary(tmp_path, monkeypatch):
    out_root = tmp_path / "strip_1d"
    factors = [0.01, 0.02]

    def fake_run(shift_factor):
        return {
            "cladding_inward_shift_factor": shift_factor,
            "cladding_inward_shift": shift_factor * run_strip_1d.A,
            "strip_cut_angle_degrees": run_strip_1d.STRIP_CUT_ANGLE,
            "out_dir": str(out_root / run_strip_1d.strip_case_stem(shift_factor)),
            "bulk_radius": run_strip_1d.BULK_RADIUS,
            "cladding_layers": run_strip_1d.CLADDING_LAYERS,
            "mesh_auto_size": run_strip_1d.MESH_AUTO_SIZE,
            "eigenmode_count": run_strip_1d.EIGENMODE_COUNT,
            "case": "strip_1d",
        }

    monkeypatch.setattr(run_strip_1d, "OUT_DIR", out_root)
    monkeypatch.setattr(run_strip_1d, "CLADDING_INWARD_SHIFT_FACTORS", factors)
    monkeypatch.setattr(run_strip_1d, "run_case_for_cladding_shift_factor", fake_run)
    monkeypatch.setattr(run_strip_1d, "RUN_STRIP_CUT", False)

    summaries = run_strip_1d.run_case()
    cut_root = out_root / "cut_60"
    summary = json.loads((cut_root / "run_summary.json").read_text(encoding="utf-8"))

    assert summaries == summary["runs"]
    assert summary["output_layout_version"] == 3
    assert [row["cladding_inward_shift_factor"] for row in summary["runs"]] == factors
    assert all(
        row["strip_cut_angle_degrees"] == run_strip_1d.STRIP_CUT_ANGLE
        for row in summary["runs"]
    )


def test_gamma_trend_auto_run_condition_uses_more_than_three_shifts_by_default(monkeypatch):
    assert run_strip_1d.RUN_GAMMA_TREND_AFTER_SCAN is None
    monkeypatch.setattr(run_strip_1d, "RUN_STRIP_CUT", True)
    monkeypatch.setattr(
        run_strip_1d,
        "CLADDING_INWARD_SHIFT_FACTORS",
        [0.0, 0.01, 0.02, 0.03],
    )
    assert run_strip_1d.should_run_gamma_trend_after_scan() is True

    monkeypatch.setattr(
        run_strip_1d,
        "CLADDING_INWARD_SHIFT_FACTORS",
        [0.0, 0.01, 0.02],
    )
    assert run_strip_1d.should_run_gamma_trend_after_scan() is False

    assert run_strip_1d.should_run_gamma_trend_after_scan(shift_count=4) is True

    monkeypatch.setattr(run_strip_1d, "RUN_GAMMA_TREND_AFTER_SCAN", True)
    assert run_strip_1d.should_run_gamma_trend_after_scan() is True

    monkeypatch.setattr(run_strip_1d, "RUN_GAMMA_TREND_AFTER_SCAN", False)
    assert run_strip_1d.should_run_gamma_trend_after_scan() is False

    monkeypatch.setattr(run_strip_1d, "RUN_STRIP_CUT", False)
    monkeypatch.setattr(run_strip_1d, "RUN_GAMMA_TREND_AFTER_SCAN", True)
    assert run_strip_1d.should_run_gamma_trend_after_scan() is False


def test_merge_strip_case_summaries_keeps_completed_supplemental_points():
    earlier = [
        {"cladding_inward_shift_factor": 0.0, "out_dir": "shift0.000"},
        {"cladding_inward_shift_factor": 0.01, "out_dir": "shift0.010"},
    ]
    supplemental = [
        {"cladding_inward_shift_factor": 0.01, "out_dir": "new-shift0.010"},
        {"cladding_inward_shift_factor": 0.02, "out_dir": "shift0.020"},
    ]

    merged = run_strip_1d.merge_strip_case_summaries(earlier, supplemental)

    assert [row["cladding_inward_shift_factor"] for row in merged] == [
        0.0,
        0.01,
        0.02,
    ]
    assert merged[1]["out_dir"] == "new-shift0.010"


def test_run_case_invokes_gamma_trend_after_summary_is_written(tmp_path, monkeypatch):
    out_root = tmp_path / "strip_1d"
    calls = []

    monkeypatch.setattr(run_strip_1d, "OUT_DIR", out_root)
    monkeypatch.setattr(run_strip_1d, "CLADDING_INWARD_SHIFT_FACTORS", [0.01, 0.02])
    monkeypatch.setattr(
        run_strip_1d,
        "run_case_for_cladding_shift_factor",
        lambda shift: {"cladding_inward_shift_factor": shift},
    )
    monkeypatch.setattr(run_strip_1d, "RUN_STRIP_CUT", True)
    monkeypatch.setattr(run_strip_1d, "RUN_GAMMA_TREND_AFTER_SCAN", True)

    def fake_analyze(summary_path):
        assert len(json.loads(summary_path.read_text(encoding="utf-8"))["runs"]) == 2
        calls.append(summary_path)

    monkeypatch.setattr(run_strip_1d, "run_gamma_trend_after_scan", fake_analyze)
    run_strip_1d.run_case()

    assert calls == [out_root / "cut_60" / "run_summary.json"]


def test_run_gamma_trend_after_scan_writes_outputs_in_trend_subdirectory(monkeypatch):
    calls = []
    module_name = "scripts.analysis.srip1d_trend"

    def fake_generate(*, run_summary_path, output_dir, target_frequency):
        calls.append((run_summary_path, output_dir, target_frequency))
        return output_dir / "trend.csv", [], output_dir / "trend.png"

    monkeypatch.setitem(
        sys.modules,
        module_name,
        types.SimpleNamespace(
            TREND_OUTPUT_DIR_NAME="strip1d_trend",
            generate_trend_from_run_summary=fake_generate,
        ),
    )
    monkeypatch.setattr(run_strip_1d, "OUT_DIR", Path("strip-output"))
    summary_path = Path("strip-output/cut_60") / "run_summary.json"

    result = run_strip_1d.run_gamma_trend_after_scan(summary_path)

    target = float(run_strip_1d.EIGENFREQUENCY_SHIFT.split()[0])
    trend_dir = Path("strip-output/cut_60/strip1d_trend")
    assert calls == [(summary_path, trend_dir, target)]
    assert result == (
        trend_dir / "trend.csv",
        [],
        trend_dir / "trend.png",
    )


def test_strip_main_rejects_finite_full_output(monkeypatch, tmp_path):
    monkeypatch.setattr(run_strip_1d, "OUT_DIR", tmp_path / "strip")
    monkeypatch.setattr(run_strip_1d, "RUN_FINITE_FULL", True)

    with pytest.raises(ValueError, match="run_finite.py"):
        run_strip_1d.run_case()
