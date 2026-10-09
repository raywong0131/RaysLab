import sys
import importlib

import numpy as np
import pytest

import comsol_workflow.field_plotting as field_plotting
from pathlib import Path


import comsol_workflow.simulation_utils as simulation_utils
from comsol_workflow.simulation_utils import SimulationRun


class RecordingNode:
    def __init__(self):
        self.settings = []
        self.children = {}
        self.run_count = 0

    def feature(self, tag):
        return self.children[tag]

    def export(self, tag=None):
        export_root = self.children["export"]
        return export_root if tag is None else export_root.children[tag]

    def set(self, name, value):
        self.settings.append((name, value))

    def run(self):
        self.run_count += 1


def make_run():
    result = RecordingNode()
    pg1 = RecordingNode()
    surf1 = RecordingNode()
    img1 = RecordingNode()
    plot1 = RecordingNode()
    export_root = RecordingNode()
    pg1.children["surf1"] = surf1
    export_root.children.update({"img1": img1, "plot1": plot1})
    result.children.update({"pg1": pg1, "export": export_root})

    class JavaModel:
        def result(self, tag=None):
            return result if tag is None else result.children[tag]

    run = SimulationRun.__new__(SimulationRun)
    run.client = type("Client", (), {"remove": lambda self, model: None})()
    run.model = type("Model", (), {"java": JavaModel()})()
    run.plane_datasets = {"center": "cpl1"}
    run.setup_plotting = lambda: None
    return run, surf1, img1


def disable_jpype_array_conversion(monkeypatch):
    monkeypatch.setattr(simulation_utils, "JInt", int)
    monkeypatch.setattr(
        simulation_utils,
        "JArray",
        lambda _value_type, _dimensions: lambda values: values,
    )


def test_export_2d_fields_applies_custom_filename_and_color_table(monkeypatch, tmp_path):
    disable_jpype_array_conversion(monkeypatch)
    run, surf1, img1 = make_run()

    run.export_2d_fields(
        0,
        "ewfd.normE",
        "unused",
        tmp_path,
        plane="center",
        export_data=False,
        image_filename="00_Hz_Re_2d.png",
        color_table="Wave",
        color_scale_mode="linearsymmetric",
        show_colorbar=True,
    )

    assert ("colortable", "Wave") in surf1.settings
    assert ("colorscalemode", "linearsymmetric") in surf1.settings
    assert ("colorlegend", "on") in surf1.settings
    assert ("options2d", "on") in img1.settings
    assert ("legend2d", "on") in img1.settings
    assert ("pngfilename", str(tmp_path / "00_Hz_Re_2d.png")) in img1.settings


@pytest.mark.parametrize("module_name", [
    "comsol_workflow.simulation_utils",
    "comsol_workflow.simulation_spatial.hexagon_unit_cell",
    "comsol_workflow.simulation_spatial.rectangle_finite_size",
    "comsol_workflow.simulation_spatial.square_unit_cell",
])
@pytest.mark.parametrize("expr,label", [("ewfd.Hz", "Re(Hz)"), ("ewfd.Hz*(-i)", "Im(Hz)")])
def test_hz_exports_use_shared_renderer_and_preserve_raw_data(monkeypatch, tmp_path, module_name, expr, label):
    module = importlib.import_module(module_name)
    monkeypatch.setattr(module, "JInt", int)
    monkeypatch.setattr(module, "JArray", lambda *_args: lambda values: values)
    run, surface, image = make_run()
    coordinates = np.array([[0, 0], [1, 0], [0, 1]])
    values = np.array([-3., 1., 2.])
    run.get_2d_fields = lambda *args: (coordinates, values)
    run.get_eigenfrequencies = lambda: [[194.9, 0.1, 1000.]] * 4
    captured = []
    monkeypatch.setattr(field_plotting, "save_field_plot", lambda *args, **kwargs: captured.append((args, kwargs)))
    module.SimulationRun.export_2d_fields(run, 3, expr, "component", tmp_path)
    assert len(captured) == 1
    args, kwargs = captured[0]
    assert args[0] == tmp_path / "03_component_center_2d.png"
    np.testing.assert_array_equal(args[2], values)
    assert kwargs["quantity_label"] == label
    assert ("expr", expr) in surface.settings
    assert run.model.java.result().export("plot1").run_count == 1
    assert image.run_count == 0


def test_hz_export_preserves_custom_filename(monkeypatch, tmp_path):
    disable_jpype_array_conversion(monkeypatch)
    run, _surface, _image = make_run()
    paths = []
    monkeypatch.setattr(field_plotting, "save_hz_simulation_plot", lambda _run, _idx, _expr, path, _plane: paths.append(path) or True)
    run.export_2d_fields(0, "ewfd.Hz", "Hz", tmp_path, export_data=False, image_filename="Hz_Re_2d.png")
    assert paths == [str(tmp_path / "Hz_Re_2d.png")]
