from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_main" / "run_band_pair.py"


def load_band():
    spec = importlib.util.spec_from_file_location("run_band_pair_model_reuse", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


band = load_band()


class FakeParam:
    def __init__(self):
        self.values = {}

    def set(self, name, value):
        self.values[name] = value


class FakeMesh:
    def __init__(self):
        self.auto_sizes = []
        self.run_count = 0

    def autoMeshSize(self, value):
        self.auto_sizes.append(value)

    def run(self):
        self.run_count += 1


class FakeSolution:
    def __init__(self):
        self.clear_count = 0
        self.run_count = 0

    def clearSolutionData(self):
        self.clear_count += 1

    def runAll(self):
        self.run_count += 1


class FakeComponent:
    def __init__(self, mesh):
        self._mesh = mesh

    def mesh(self, name):
        assert name == "mesh1"
        return self._mesh


class FakeJava:
    def __init__(self):
        self._param = FakeParam()
        self._mesh = FakeMesh()
        self._solution = FakeSolution()

    def param(self):
        return self._param

    def component(self, name):
        assert name == "comp1"
        return FakeComponent(self._mesh)

    def sol(self, name):
        assert name == "sol1"
        return self._solution


class FakeModel:
    def __init__(self):
        self.java = FakeJava()


class FakeClient:
    def remove(self, _model):
        return None


def make_runner(mesh_auto_size=5):
    runner = band.ReusableSimulationRun.__new__(band.ReusableSimulationRun)
    runner.config = band.SimulationConfig(mesh_auto_size=mesh_auto_size)
    runner.client = FakeClient()
    runner.model = FakeModel()
    runner._plotting_initialized = False
    runner._geometry_signature = None
    runner._applied_mesh_auto_size = None
    runner._current_k = None
    runner.model_build_count = 0
    runner.mesh_build_count = 0
    runner.solve_count = 0
    return runner


def test_canonical_entry_uses_reused_model_and_standard_output():
    source = SCRIPT_PATH.read_text(encoding="utf-8")

    assert SCRIPT_PATH.exists()
    assert sorted(path.name for path in SCRIPT_PATH.parent.glob("run_band_pair*.py")) == [
        "run_band_pair.py"
    ]
    assert "class ReusableSimulationRun" in source
    assert band.OUT_DIR.parent == band.UNIT_CELL_OUTPUT_ROOT
    assert band.OUT_DIR.name == (
        f"unitcell_band_{band.ACTIVE_PARAMETERS.unit_cell_case_parameter_label}"
    )
    assert band.EIGENMODE_COUNT == band.ACTIVE_PARAMETERS.unit_cell_eigenmode_count
    assert band.EIGENMODE_COUNT > 0
    assert band.OUT_DIR.name.endswith(f"_mesh{band.MESH_AUTO_SIZE}")


def test_reusable_runner_reuses_model_and_rebuilds_only_changed_mesh(monkeypatch):
    runner = make_runner()
    holes = [[[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]]]
    first_builds = []

    def fake_first_build(self, a, received_holes, k):
        first_builds.append((a, received_holes, dict(k)))

    monkeypatch.setattr(band.BaseSimulationRun, "build_and_run", fake_first_build)

    runner.build_and_run(0.82, holes, {"kx": 0.0, "ky": 0.0})
    runner.build_and_run(0.82, holes, {"kx": 0.1, "ky": 0.0})
    runner.build_and_run(
        0.82,
        holes,
        {"kx": 0.2, "ky": 0.0},
        mesh_auto_size=6,
    )

    assert len(first_builds) == 1
    assert runner.model_build_count == 1
    assert runner.mesh_build_count == 2
    assert runner.solve_count == 3
    assert runner.model.java._mesh.auto_sizes == [6]
    assert runner.model.java._mesh.run_count == 1
    assert runner.model.java._solution.clear_count == 2
    assert runner.model.java._solution.run_count == 2
    assert runner.model.java._param.values == {"kx": "0.2*G", "ky": "0.0*G"}
    assert runner.execution_stats["current_k"] == {"kx": 0.2, "ky": 0.0}


def test_reusable_runner_rejects_changed_geometry(monkeypatch):
    runner = make_runner()
    holes = [[[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]]]
    monkeypatch.setattr(band.BaseSimulationRun, "build_and_run", lambda *_args: None)
    runner.build_and_run(0.82, holes, {"kx": 0.0, "ky": 0.0})

    changed = [[[0.0, 0.0], [0.2, 0.0], [0.0, 0.1]]]
    with pytest.raises(ValueError, match="different geometry"):
        runner.build_and_run(0.82, changed, {"kx": 0.1, "ky": 0.0})


def test_run_k_point_uses_caller_owned_runner(monkeypatch, tmp_path):
    point = {
        "direction": "gamma_m",
        "point_index": 1,
        "kx": 0.1,
        "ky": 0.0,
        "kx_str": "0.100000",
        "ky_str": "0.000000",
        "k_norm": 0.1,
    }
    calls = []

    class Runner:
        def build_and_run(self, a, holes, k, mesh_auto_size):
            calls.append((a, holes, k, mesh_auto_size))

    expected = pd.DataFrame({"re": [198.0]})
    monkeypatch.setattr(
        band,
        "k_point_is_complete",
        lambda _path, _config=None: False,
    )
    monkeypatch.setattr(
        band,
        "persist_k_point_solution",
        lambda runner, output_dir: expected,
    )

    result = band.run_k_point(
        tmp_path,
        [[[0.0, 0.0]]],
        point,
        band.SimulationConfig(mesh_auto_size=5),
        sim_run=Runner(),
    )

    assert result is expected
    assert calls == [
        (
            band.A,
            [[[0.0, 0.0]]],
            {"kx": 0.1, "ky": 0.0},
            band.MESH_AUTO_SIZE,
        )
    ]


def test_mesh_size_validation():
    assert band.ReusableSimulationRun._validate_mesh_auto_size(None) is None
    assert band.ReusableSimulationRun._validate_mesh_auto_size(5) == 5
    for invalid in (0, 10, 5.0, True):
        with pytest.raises(ValueError, match="mesh_auto_size"):
            band.ReusableSimulationRun._validate_mesh_auto_size(invalid)


@pytest.mark.parametrize(
    ("cache_complete", "expected_runner_count"),
    ((False, 1), (True, 0)),
)
def test_run_case_owns_at_most_one_runner(
    monkeypatch,
    tmp_path,
    cache_complete,
    expected_runner_count,
):
    points = band.make_branch_k_points("gamma_m", 2)
    runners = []
    received_runners = []

    class Runner:
        def __init__(self, _config):
            runners.append(self)
            self.cleared = False

        @property
        def execution_stats(self):
            return {
                "execution_engine": "reused_model",
                "model_reuse_enabled": True,
                "model_build_count": 1,
                "mesh_build_count": 1,
                "solve_count": len(received_runners),
                "mesh_auto_size": band.MESH_AUTO_SIZE,
                "current_k": {"kx": 0.0, "ky": 0.1},
            }

        def clear(self):
            self.cleared = True

    def fake_run_k_point(case_dir, _holes, point, _config, sim_run=None):
        received_runners.append(sim_run)
        output_dir = band.k_point_dir(case_dir, point)
        output_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"re": [198.0]}).to_csv(
            output_dir / "eigenfrequencies.csv",
            index=False,
        )
        return pd.DataFrame({"re": [198.0]})

    tracked = pd.DataFrame(
        {
            "direction": ["gamma_m"],
            "point_index": [0],
            "band_label": ["p1"],
            "match_status": ["gamma"],
        }
    )
    monkeypatch.setattr(band, "ReusableSimulationRun", Runner)
    monkeypatch.setattr(
        band,
        "create_hexagon_design",
        lambda *_args: (
            object(),
            [np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]])],
            [{"min_dist": band.D + 0.01}],
        ),
    )
    monkeypatch.setattr(band, "visualize_hexagon_design", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        band,
        "case_k_points",
        lambda _name: {"gamma_m": points, "gamma_k": points},
    )
    monkeypatch.setattr(
        band,
        "k_point_is_complete",
        lambda _path, _config=None: cache_complete,
    )
    monkeypatch.setattr(band, "run_k_point", fake_run_k_point)
    monkeypatch.setattr(band, "select_gamma_bands", lambda _frame: {"p1": 0})
    monkeypatch.setattr(
        band,
        "track_selected_branch",
        lambda *_args, **_kwargs: tracked.copy(),
    )
    monkeypatch.setattr(
        band,
        "_mode_composition_table",
        lambda *_args, **_kwargs: pd.DataFrame({"mode": []}),
    )

    band.run_case("cavity", band.BULK_PARAMS, tmp_path)

    assert len(runners) == expected_runner_count
    if runners:
        assert runners[0].cleared is True
        assert all(runner is runners[0] for runner in received_runners)
    else:
        assert all(runner is None for runner in received_runners)
    config = json.loads(
        (tmp_path / "99_config" / "cavity" / "config.json").read_text(
            encoding="utf-8"
        )
    )
    assert config["execution"]["cache_hits"] == (2 if cache_complete else 0)
    assert config["execution"]["cache_misses"] == (0 if cache_complete else 2)
    assert config["execution"]["model_build_count"] == expected_runner_count
