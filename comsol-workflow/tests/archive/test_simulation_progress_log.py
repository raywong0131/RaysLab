import os
import sys
import tempfile
import types
from pathlib import Path

import pytest

import comsol_workflow.mesh_constrcution as mesh_constrcution
import comsol_workflow.simulation_utils as simulation_utils


class _NoopJavaObject:
    def __init__(self, events=None):
        self.events = events if events is not None else []
        self.auto_mesh_sizes = []

    def __getattr__(self, _name):
        return self

    def __call__(self, *_args, **_kwargs):
        return self

    def __getitem__(self, _key):
        return self

    def set(self, *_args, **_kwargs):
        return self

    def create(self, *_args, **_kwargs):
        return self

    def run(self, *_args, **_kwargs):
        self.events.append(("run", None))
        return self

    def runAll(self, *_args, **_kwargs):
        self.events.append(("runAll", None))
        return self

    def autoMeshSize(self, value):
        self.auto_mesh_sizes.append(value)
        self.events.append(("autoMeshSize", value))
        return self

    def setResult(self, *_args, **_kwargs):
        return self

    def label(self, *_args, **_kwargs):
        return self


class _FakeModel:
    def __init__(self):
        self.events = []
        self.java = _NoopJavaObject(self.events)

    def save(self, path):
        self.events.append(("save", Path(path).name))


class _FakeClient:
    def create(self, _name):
        return _FakeModel()

    def remove(self, _model):
        return None


class _FakeMeshSelection:
    def __init__(self, events, tag):
        self.events = events
        self.tag = tag

    def named(self, name):
        self.events.append(("selection_named", self.tag, name))
        return self


class _FakeMeshFeature:
    def __init__(self, events, tag):
        self.events = events
        self.tag = tag
        self.children = {}
        self._selection = _FakeMeshSelection(events, tag)

    def selection(self):
        return self._selection

    def set(self, key, value):
        self.events.append(("feature_set", self.tag, key, value))
        return self

    def create(self, tag, kind):
        child_tag = f"{self.tag}/{tag}"
        self.children[tag] = _FakeMeshFeature(self.events, child_tag)
        self.events.append(("child_create", self.tag, tag, kind))
        return self.children[tag]

    def feature(self, tag):
        return self.children[tag]


class _FakeMesh:
    def __init__(self):
        self.events = []
        self.features = {"size": _FakeMeshFeature(self.events, "size")}

    def create(self, tag, kind):
        self.features[tag] = _FakeMeshFeature(self.events, tag)
        self.events.append(("feature_create", tag, kind))
        return self.features[tag]

    def feature(self, tag):
        return self.features[tag]


class _FakeMeshJModel:
    def __init__(self):
        self.mesh1 = _FakeMesh()

    def component(self, name):
        assert name == "comp1"
        return self

    def mesh(self, tag):
        assert tag == "mesh1"
        return self.mesh1


@pytest.mark.parametrize(
    ("mesh_auto_size", "expected"),
    [
        (1, (1.3, 0.2, 1.0)),
        (2, (1.35, 0.3, 0.85)),
        (3, (1.4, 0.4, 0.7)),
        (4, (1.45, 0.5, 0.6)),
        (5, (1.5, 0.6, 0.5)),
        (6, (1.6, 0.7, 0.4)),
        (7, (1.7, 0.8, 0.3)),
        (8, (1.85, 0.9, 0.2)),
        (9, (2.0, 1.0, 0.1)),
    ],
)
def test_manual_mesh_size_level_mapping(mesh_auto_size, expected):
    assert (
        mesh_constrcution.resolve_manual_mesh_size_parameters(mesh_auto_size)
        == expected
    )


@pytest.mark.parametrize("mesh_auto_size", [0, 10, 1.5, True, "5", None])
def test_manual_mesh_size_level_rejects_invalid_values(mesh_auto_size):
    with pytest.raises(ValueError, match="integer from 1 to 9"):
        mesh_constrcution.resolve_manual_mesh_size_parameters(mesh_auto_size)


def test_manual_mesh_construction_uses_material_named_selections_and_wavelength_sizes():
    jmodel = _FakeMeshJModel()
    layers = [types.SimpleNamespace(refractive_index="3.3")]
    original_jint = getattr(mesh_constrcution, "JInt", None)
    mesh_constrcution.JInt = lambda value: ("JInt", value)
    try:
        mesh_constrcution.construct_manual_mesh(jmodel, layers, mesh_auto_size=9)
    finally:
        if original_jint is None:
            delattr(mesh_constrcution, "JInt")
        else:
            mesh_constrcution.JInt = original_jint
    events = jmodel.mesh1.events

    assert ("feature_set", "size", "hmax", "lda0/5") in events, events
    assert ("feature_set", "size", "hmin", "(lda0/5)/33.3") in events, events
    for tag in ("size", "size_layer0"):
        assert ("feature_set", tag, "hgrad", 2.0) in events, events
        assert ("feature_set", tag, "hcurve", 1.0) in events, events
        assert ("feature_set", tag, "hnarrow", 0.1) in events, events
    for property_name in ("hgradactive", "hcurveactive", "hnarrowactive"):
        assert ("feature_set", "size", property_name, True) not in events, events
        assert ("feature_set", "size_layer0", property_name, True) in events, events
    assert ("feature_create", "size_layer0", "Size") in events, events
    assert ("selection_named", "size_layer0", "sel_layer_0_mat_dom") in events, events
    assert ("feature_set", "size_layer0", "hmax", "lda0/(3.3)/5") in events, events
    assert ("feature_set", "size_layer0", "hmin", "(lda0/(3.3)/5)/33.3") in events, events
    assert ("feature_create", "ftet1", "FreeTet") in events, events
    assert ("selection_named", "ftet1", "sel_physical_dom") in events, events
    assert ("feature_create", "swe1", "Sweep") in events, events
    assert ("selection_named", "swe1", "sel_pml_dom") in events, events
    assert ("child_create", "swe1", "dis1", "Distribution") in events, events
    assert (
        "feature_set",
        "swe1/dis1",
        "numelem",
        ("JInt", 8),
    ) in jmodel.mesh1.events


def test_manual_mesh_construction_defaults_to_level_five():
    jmodel = _FakeMeshJModel()
    original_jint = getattr(mesh_constrcution, "JInt", None)
    mesh_constrcution.JInt = lambda value: ("JInt", value)
    try:
        mesh_constrcution.construct_manual_mesh(jmodel, [])
    finally:
        if original_jint is None:
            delattr(mesh_constrcution, "JInt")
        else:
            mesh_constrcution.JInt = original_jint

    events = jmodel.mesh1.events
    assert ("feature_set", "size", "hgrad", 1.5) in events, events
    assert ("feature_set", "size", "hcurve", 0.6) in events, events
    assert ("feature_set", "size", "hnarrow", 0.5) in events, events


def test_run_simulation_writes_comsol_progress_to_requested_file():
    calls = []
    original_start = simulation_utils.mph.start
    original_jclass = getattr(simulation_utils, "JClass", None)

    class _FakeModelUtil:
        @staticmethod
        def showProgress(value):
            calls.append(value)

    try:
        simulation_utils.mph.start = lambda: _FakeClient()
        simulation_utils.JClass = lambda name: _FakeModelUtil if name == "com.comsol.model.util.ModelUtil" else None

        with tempfile.TemporaryDirectory() as tmpdir:
            progress_log = Path(tmpdir) / "comsol_progress.log"
            sim = simulation_utils.SimulationRun(progress_log_path=progress_log)
            try:
                sim.run_simulation()
            finally:
                sim.clear()

            assert calls == [str(progress_log), False], calls
    finally:
        simulation_utils.mph.start = original_start
        if original_jclass is None:
            delattr(simulation_utils, "JClass")
        else:
            simulation_utils.JClass = original_jclass


def test_run_simulation_sets_auto_mesh_size_to_9():
    original_start = simulation_utils.mph.start

    try:
        simulation_utils.mph.start = lambda: _FakeClient()

        sim = simulation_utils.SimulationRun()
        try:
            sim.run_simulation()
            assert sim.model.java.auto_mesh_sizes == [9], sim.model.java.auto_mesh_sizes
        finally:
            sim.clear()
    finally:
        simulation_utils.mph.start = original_start


def test_run_simulation_uses_manual_mesh_builder_without_auto_mesh_size():
    original_start = simulation_utils.mph.start

    try:
        simulation_utils.mph.start = lambda: _FakeClient()

        sim = simulation_utils.SimulationRun()
        try:
            def build_manual_mesh(jmodel, config):
                jmodel.events.append(("manual_mesh_builder", config.mesh_auto_size))

            sim.run_simulation(mesh_builder=build_manual_mesh)
            events = sim.model.events
            assert sim.model.java.auto_mesh_sizes == [], sim.model.java.auto_mesh_sizes
            assert events.index(("manual_mesh_builder", 9)) < events.index(("run", None))
            assert events.index(("run", None)) < events.index(("runAll", None))
        finally:
            sim.clear()
    finally:
        simulation_utils.mph.start = original_start


def test_run_simulation_saves_requested_model_after_mesh_before_solving():
    original_start = simulation_utils.mph.start

    try:
        simulation_utils.mph.start = lambda: _FakeClient()

        with tempfile.TemporaryDirectory() as tmpdir:
            mesh_save_path = Path(tmpdir) / "finite_cavity.mph"
            sim = simulation_utils.SimulationRun()
            try:
                sim.run_simulation(mesh_save_path=mesh_save_path)
                events = sim.model.events
                save_event = ("save", "finite_cavity.mph")
                assert save_event in events, events
                assert events.index(("run", None)) < events.index(save_event)
                assert events.index(save_event) < events.index(("runAll", None))
            finally:
                sim.clear()
    finally:
        simulation_utils.mph.start = original_start


if __name__ == "__main__":
    test_manual_mesh_construction_uses_material_named_selections_and_wavelength_sizes()
    test_run_simulation_writes_comsol_progress_to_requested_file()
    test_run_simulation_sets_auto_mesh_size_to_9()
    test_run_simulation_uses_manual_mesh_builder_without_auto_mesh_size()
    test_run_simulation_saves_requested_model_after_mesh_before_solving()
    print("All tests passed.")
