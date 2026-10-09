import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts.run_main import run_finite
from scripts.run_main import run_finite_quarter
from scripts.run_main import run_finite_quarter_four_modes as four_modes


def test_all_symmetry_names_preserve_series_suffix(monkeypatch, tmp_path):
    active = SimpleNamespace(
        structure_series_label=lambda prefix: (
            f"{prefix}_20-20_cell_mesh5_shiftprof-power2"
        )
    )
    monkeypatch.setattr(run_finite, "OUTPUT_ROOT", tmp_path)

    assert four_modes.all_symmetry_series_dir(active) == (
        tmp_path
        / "finite_quarter_all_20-20_cell_mesh5_shiftprof-power2"
    )
    assert four_modes.mode_uid(7, 1) == "mode7_1_xPEC_yPMC"
    assert four_modes.mode_uid(7, 2) == "mode7_2_xPMC_yPEC"
    assert four_modes.mode_uid(7, 3) == "mode7_3_xPEC_yPEC"
    assert four_modes.mode_uid(7, 4) == "mode7_4_xPMC_yPMC"


def test_copy_solution_is_verified_before_return():
    events = []
    working = np.array([[198.0, -0.01, 9900.0]])

    class FakeSimulation:
        def get_eigenfrequencies(self, dataset="dset1"):
            events.append(("read", dataset))
            return working.copy()

        def copy_solution(self, solution_tag, dataset_tag, *, label):
            events.append(("copy", solution_tag, dataset_tag, label))

    result = four_modes.copy_and_validate_solution(FakeSimulation(), 3)

    assert events == [
        ("read", "dset1"),
        ("copy", "solsym3", "dsetsym3", "Internal symmetry solution 3"),
        ("read", "dsetsym3"),
    ]
    assert result["eigenfrequency_verified"] is True


def test_solution_field_probe_falls_back_to_first_invalid_mode():
    class FakeSimulation:
        def get_2d_fields(self, *_args):
            return np.array([[0.0, 0.0]]), np.array([1.0])

        def get_fields_at_coordinates(self, *_args, **_kwargs):
            return np.array([[1.0 + 0.0j]])

    archive = {"dataset_tag": "dsetsym4"}
    four_modes.validate_solution_field_probe(
        FakeSimulation(),
        pd.DataFrame({"mode_idx": [0], "is_valid": [False]}),
        archive,
    )

    assert archive["field_probe_verified"] is True
    assert archive["field_probe_mode_valid"] is False


def test_load_resume_state_accepts_verified_continuous_prefix(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setattr(run_finite, "EIGENMODE_COUNT", 2)
    monkeypatch.setattr(four_modes, "_validate_resume_identity", lambda *args: None)
    paths = run_finite.prepare_finite_case_output(
        tmp_path / "shiftx0.000_shifty0.000",
        model_filename=four_modes.MODEL_FILENAME,
        resume=False,
    )
    paths["mph"].write_bytes(b"saved-model")
    paths["config"].write_text("{}", encoding="utf-8")
    runs = []
    for symmetry_id in (1, 2, 3):
        case = run_finite_quarter.symmetry_case(symmetry_id)
        partition = four_modes._prepare_partition(
            paths["staging_dir"], symmetry_id
        )
        partition["export_dir"].mkdir(parents=True)
        pd.DataFrame(
            {
                "re": [203.0, 204.0],
                "im": [-0.1, -0.1],
                "q": [1000.0, 1100.0],
                "is_valid": [False, False],
                "mode_idx": [0, 1],
                "symmetry_id": [symmetry_id, symmetry_id],
                "x_boundary": [case["x_boundary"]] * 2,
                "y_boundary": [case["y_boundary"]] * 2,
            }
        ).to_csv(
            partition["export_dir"] / "eigenfrequencies.csv",
            index=False,
        )
        (partition["export_dir"] / "air_field_metadata.json").write_text(
            "{}", encoding="utf-8"
        )
        partition["config"].write_text("{}", encoding="utf-8")
        runs.append(
            {
                **case,
                **four_modes.solution_identity(symmetry_id),
                "mode_count": 2,
                "status": "solution_copied_exported_verified",
            }
        )
    summary = {
        "status": "running",
        "solve_count": 3,
        "solution_copy_count": 3,
        "symmetry_runs_internal": runs,
    }
    (paths["config_dir"] / four_modes.SUMMARY_FILENAME).write_text(
        json.dumps(summary), encoding="utf-8"
    )

    _metadata, loaded, partitions = four_modes._load_resume_state(
        paths, 0.0, 0.0
    )

    assert [run["id"] for run in loaded["symmetry_runs_internal"]] == [1, 2, 3]
    assert len(partitions) == 3


def _make_finalized_partition(root_paths, symmetry_id, mode_count):
    partition_dir = four_modes._partition_dir(
        root_paths["staging_dir"], symmetry_id
    )
    paths = run_finite.finite_case_output_paths(
        partition_dir,
        model_filename="internal_partition.mph",
    )
    for key in (
        "results_dir",
        "overview_dir",
        "model_dir",
        "logs_dir",
        "config_dir",
    ):
        paths[key].mkdir(parents=True, exist_ok=True)
    rows = []
    for mode_idx in range(mode_count):
        mode_dir = paths["results_dir"] / f"mode{mode_idx}"
        mode_dir.mkdir()
        (mode_dir / "paths.json").write_text(
            json.dumps({"output": str(mode_dir)}),
            encoding="utf-8",
        )
        rows.append(
            {
                "mode_idx": mode_idx,
                "re": 197.0 + symmetry_id + mode_idx / 10,
                "im": -0.01,
                "q": 1000.0 + mode_idx,
                "is_valid": True,
            }
        )
    pd.DataFrame(rows).to_csv(
        paths["overview_dir"] / "eigenfrequencies.csv",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "mode_idx": 0,
                "output_dir": str(paths["results_dir"] / "mode0"),
            }
        ]
    ).to_csv(
        paths["overview_dir"] / "finite_lattice_fourier_summary.csv",
        index=False,
    )
    (paths["overview_dir"] / "index_r_cavity.png").write_bytes(b"same-r")
    (paths["overview_dir"] / "index_k_1stBZ.png").write_bytes(b"same-k")
    (paths["config_dir"] / "config.json").write_text(
        json.dumps({"case_dir": str(partition_dir)}),
        encoding="utf-8",
    )


def test_merge_partitions_builds_one_dataset(monkeypatch, tmp_path):
    monkeypatch.setattr(run_finite, "EIGENMODE_COUNT", 2)
    root_paths = run_finite.prepare_finite_case_output(
        tmp_path / "shiftx0.000_shifty0.000",
        model_filename=four_modes.MODEL_FILENAME,
        resume=False,
    )
    for symmetry_id in four_modes.SYMMETRY_IDS:
        _make_finalized_partition(root_paths, symmetry_id, 2)

    result = four_modes._merge_partition_outputs(root_paths)

    expected = {
        four_modes.mode_uid(mode_idx, symmetry_id)
        for symmetry_id in four_modes.SYMMETRY_IDS
        for mode_idx in range(2)
    }
    assert {path.name for path in root_paths["results_dir"].iterdir()} == expected
    merged = pd.read_csv(
        root_paths["overview_dir"] / "eigenfrequencies.csv"
    )
    assert len(merged) == 8
    assert merged["mode_uid"].nunique() == 8
    assert result["mode_count"] == 8
    assert not root_paths["staging_dir"].exists()
    assert (
        root_paths["config_dir"]
        / "internal_symmetry"
        / "1"
        / "config.json"
    ).is_file()
    assert (
        root_paths["overview_dir"] / "index_r_cavity.png"
    ).read_bytes() == b"same-r"
    mode_payload = json.loads(
        (
            root_paths["results_dir"]
            / "mode0_1_xPEC_yPMC"
            / "paths.json"
        ).read_text(encoding="utf-8")
    )
    assert str(root_paths["staging_dir"]) not in mode_payload["output"]
    assert mode_payload["output"].endswith("mode0_1_xPEC_yPMC")


def test_merge_partitions_uses_actual_returned_mode_count(monkeypatch, tmp_path):
    monkeypatch.setattr(run_finite, "EIGENMODE_COUNT", 2)
    root_paths = run_finite.prepare_finite_case_output(
        tmp_path / "shiftx0.000_shifty0.000",
        model_filename=four_modes.MODEL_FILENAME,
        resume=False,
    )
    for symmetry_id in four_modes.SYMMETRY_IDS:
        _make_finalized_partition(root_paths, symmetry_id, 3)

    result = four_modes._merge_partition_outputs(root_paths)

    assert result["mode_count"] == 12
    assert result["requested_eigenmode_count_per_symmetry"] == 2
    assert result["actual_mode_count_per_symmetry"] == {
        str(symmetry_id): 3 for symmetry_id in four_modes.SYMMETRY_IDS
    }
    assert len(pd.read_csv(root_paths["overview_dir"] / "eigenfrequencies.csv")) == 12


def test_frequency_q_plot_uses_one_linear_scatter(monkeypatch, tmp_path):
    source = tmp_path / "eigenfrequencies.csv"
    pd.DataFrame(
        {
            "re": [196.0, 197.0, np.nan],
            "q": [1000.0, 1200.0, 900.0],
            "symmetry_id": [1, 2, 3],
        }
    ).to_csv(source, index=False)
    calls = []

    class FakeAxis:
        def scatter(self, x, y, **kwargs):
            calls.append(("scatter", list(x), list(y), kwargs))

        def set_xscale(self, value):
            calls.append(("xscale", value))

        def set_yscale(self, value):
            calls.append(("yscale", value))

        def set_xlabel(self, value):
            calls.append(("xlabel", value))

        def set_ylabel(self, value):
            calls.append(("ylabel", value))

        def tick_params(self, **kwargs):
            calls.append(("ticks", kwargs))

    class FakeFigure:
        def savefig(self, path, **kwargs):
            calls.append(("savefig", path, kwargs))

    monkeypatch.setattr(
        four_modes.plt,
        "subplots",
        lambda **kwargs: (FakeFigure(), FakeAxis()),
    )
    monkeypatch.setattr(four_modes.plt, "close", lambda figure: None)

    result = four_modes.save_frequency_q_plot(source, tmp_path / "f_Q.png")

    scatter = next(call for call in calls if call[0] == "scatter")
    assert scatter[1:3] == ([196.0, 197.0], [1000.0, 1200.0])
    assert ("xscale", "linear") in calls
    assert ("yscale", "linear") in calls
    assert ("xlabel", "Frequency (THz)") in calls
    assert ("ylabel", "Q") in calls
    assert result["plotted_row_count"] == 2
    assert result["nonfinite_row_count"] == 1
    assert not any(call[0] in {"legend", "title"} for call in calls)


@pytest.mark.parametrize("reuse_existing_mesh", [False, True])
def test_run_shift_reuses_one_model_and_copies_before_next_boundary(
    monkeypatch,
    tmp_path,
    reuse_existing_mesh,
):
    events = []
    monkeypatch.setattr(four_modes, "OUT_DIR", tmp_path / "series")
    monkeypatch.setattr(run_finite, "EIGENMODE_COUNT", 2)
    monkeypatch.setattr(run_finite, "RUN_FARFIELD_FFT", False)
    monkeypatch.setattr(four_modes, "RUN_COMSOL", True)
    reuse_path = tmp_path / "source" / "00_model" / "finite_quarter.mph"
    monkeypatch.setattr(
        four_modes,
        "REUSE_MESHED_MODEL",
        reuse_path if reuse_existing_mesh else None,
    )
    monkeypatch.setattr(
        four_modes,
        "reuse_model_source_info",
        lambda *args: {
            "mph_path": str(reuse_path),
            "quarter_hole_count": 1,
        },
    )
    monkeypatch.setattr(
        run_finite,
        "attach_saved_model",
        lambda sim, path: events.append(("attach", path)),
    )
    monkeypatch.setattr(
        run_finite,
        "set_cladding_inward_shift_factors",
        lambda x, y: events.append(("shift", x, y)),
    )
    metadata = {
        "finite_geometry": "hex",
        "finite_boundary": np.zeros((3, 2)),
        "finite_Lx": 2.0,
        "finite_Ly": 2.0,
        "quarter_boundary": np.zeros((3, 2)),
        "quarter_hole_count": 1,
        "bulk_points": [],
        "cladding_points": [],
    }
    monkeypatch.setattr(
        run_finite_quarter,
        "build_quarter_case_geometry",
        lambda: ([{"polygon": np.zeros((3, 2))}], [], metadata, None, 2.0),
    )
    monkeypatch.setattr(four_modes, "validate_cavity_edge_cell_count", lambda value: None)
    monkeypatch.setattr(
        four_modes,
        "_case_metadata",
        lambda *args, **kwargs: {"case": "finite_quarter_all"},
    )
    monkeypatch.setattr(run_finite, "strip_validity_boundary", lambda value: np.zeros((3, 2)))
    monkeypatch.setattr(run_finite, "shift_profile_preflight_messages", lambda value: [])
    monkeypatch.setattr(run_finite, "save_finite_shift_profile_diagnostics", lambda *a, **k: None)
    monkeypatch.setattr(run_finite, "save_finite_model_geometry_figures", lambda *a, **k: None)
    monkeypatch.setattr(run_finite, "full_cell_records", lambda value: [])
    monkeypatch.setattr(four_modes, "_simulation_geometry", lambda *a: ("boundary", ["layer"]))
    monkeypatch.setattr(run_finite, "finite_mesh_builder", lambda layers: "mesh")
    monkeypatch.setattr(
        run_finite_quarter,
        "quarter_simulation_config",
        lambda case: "config",
    )

    class FakeModel:
        def save(self, path):
            events.append(("save", path))

    class FakeSimulation:
        def __init__(self, **kwargs):
            self.model = FakeModel()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def build_geometry(self, *args, **kwargs):
            events.append(("build",))

        def prepare_reused_finite_quarter_model(self):
            events.append(("prepare_reuse",))
            return {"eigenfrequency_shift": "203.6 [THz]", "eigenmode_count": 2}

        def set_finite_quarter_symmetry(self, x, y):
            events.append(("boundary", x, y))

        def run_simulation(self, **kwargs):
            events.append(("initial_solve",))

        def rerun_eigenfrequency_solution(self):
            events.append(("rerun",))

    monkeypatch.setattr(four_modes, "SimulationRun", FakeSimulation)

    def fake_copy(sim, symmetry_id):
        events.append(("copy", symmetry_id))
        return {
            **four_modes.solution_identity(symmetry_id),
            "eigenfrequency_verified": True,
            "field_probe_verified": False,
        }

    monkeypatch.setattr(four_modes, "copy_and_validate_solution", fake_copy)
    monkeypatch.setattr(four_modes, "validate_solution_field_probe", lambda *a: None)

    def fake_export(sim, export_dir, *args, **kwargs):
        export_dir.mkdir(parents=True, exist_ok=True)
        return pd.DataFrame(
            {
                "re": [197.0, 198.0],
                "im": [-0.01, -0.01],
                "q": [1000.0, 1100.0],
                "is_valid": [True, True],
                "mode_idx": [0, 1],
                "symmetry_id": [kwargs["case"]["id"]] * 2,
                "x_boundary": [kwargs["case"]["x_boundary"]] * 2,
                "y_boundary": [kwargs["case"]["y_boundary"]] * 2,
            }
        )

    monkeypatch.setattr(run_finite_quarter, "export_reconstructed_results", fake_export)
    monkeypatch.setattr(four_modes, "_run_partition_postprocess", lambda *a: None)
    monkeypatch.setattr(
        four_modes,
        "_merge_partition_outputs",
        lambda paths: {"mode_count": 8},
    )
    monkeypatch.setattr(
        four_modes,
        "save_frequency_q_plot",
        lambda *a: {"x_scale": "linear", "y_scale": "linear"},
    )

    summary = four_modes.run_all_symmetry_shift(0.0, 0.0)

    assert events.count(("build",)) == (0 if reuse_existing_mesh else 1)
    assert events.count(("initial_solve",)) == (
        0 if reuse_existing_mesh else 1
    )
    assert events.count(("rerun",)) == (4 if reuse_existing_mesh else 3)
    assert events.count(("prepare_reuse",)) == (
        1 if reuse_existing_mesh else 0
    )
    assert [event for event in events if event[0] == "boundary"] == [
        ("boundary", "PEC", "PMC"),
        ("boundary", "PMC", "PEC"),
        ("boundary", "PEC", "PEC"),
        ("boundary", "PMC", "PMC"),
    ]
    assert [event for event in events if event[0] == "copy"] == [
        ("copy", 1),
        ("copy", 2),
        ("copy", 3),
        ("copy", 4),
    ]
    for symmetry_id in (1, 2, 3):
        copy_position = events.index(("copy", symmetry_id))
        next_boundary = events.index(
            [
                ("boundary", "PMC", "PEC"),
                ("boundary", "PEC", "PEC"),
                ("boundary", "PMC", "PMC"),
            ][symmetry_id - 1]
        )
        assert any(
            event[0] == "save"
            for event in events[copy_position + 1 : next_boundary]
        )
    assert summary["status"] == "complete"
    counts = (
        summary["model_build_count"],
        summary["mesh_build_count"],
        summary["reused_model_count"],
        summary["reused_mesh_count"],
        summary["solve_count"],
        summary["solution_copy_count"],
    )
    assert counts == (
        (0, 0, 1, 1, 4, 4)
        if reuse_existing_mesh
        else (1, 1, 0, 0, 4, 4)
    )


def test_single_quarter_finalizer_keeps_plain_mode_name(tmp_path):
    paths = run_finite.prepare_finite_case_output(
        tmp_path / "single",
        model_filename="finite_quarter.mph",
        resume=False,
    )
    paths["export_dir"].mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"re": 198.0, "im": -0.01, "q": 9900.0, "mode_idx": 0}]
    ).to_csv(paths["export_dir"] / "eigenfrequencies.csv", index=False)

    run_finite.finalize_finite_case_output(
        paths["case_dir"],
        model_filename="finite_quarter.mph",
    )

    assert (paths["results_dir"] / "mode0").is_dir()
    assert not any("_xPEC_" in path.name for path in paths["results_dir"].iterdir())


def test_partition_postprocess_allows_zero_valid_modes(monkeypatch, tmp_path):
    paths = run_finite.prepare_finite_case_output(
        tmp_path / "partition",
        model_filename="internal_partition.mph",
        resume=False,
    )
    paths["export_dir"].mkdir(parents=True)
    pd.DataFrame(
        [{"mode_idx": 0, "re": 203.6, "q": 1000.0, "is_valid": False}]
    ).to_csv(paths["export_dir"] / "eigenfrequencies.csv", index=False)
    calls = []
    monkeypatch.setattr(
        run_finite,
        "run_farfield_fft_case",
        lambda *args, **kwargs: calls.append("farfield"),
    )
    monkeypatch.setattr(
        run_finite,
        "run_exported_finite_lattice_fourier_postprocess",
        lambda *args, **kwargs: calls.append("lattice"),
    )
    monkeypatch.setattr(
        run_finite,
        "score_exported_finite_modes",
        lambda *args, **kwargs: calls.append("score"),
    )
    monkeypatch.setattr(
        run_finite,
        "finalize_finite_case_output",
        lambda *args, **kwargs: calls.append("finalize"),
    )

    four_modes._run_partition_postprocess(paths, {})

    assert calls == ["farfield", "finalize"]


def test_output_dir_argument_is_available():
    args = four_modes.parse_args(
        [
            "--output-dir",
            "scripts/.out/finite_cavity/custom",
            "--reuse-meshed-model",
            "source/finite_quarter.mph",
            "--resume-incomplete",
        ]
    )
    assert args.output_dir == Path("scripts/.out/finite_cavity/custom")
    assert args.reuse_meshed_model == Path("source/finite_quarter.mph")
    assert args.resume_incomplete is True
