import sys
from pathlib import Path

import numpy as np
import pytest


import comsol_workflow.simulation_utils as simulation_utils
from comsol_workflow.simulation_utils import SimulationRun


class FakeInterp:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag
        self.settings = {}
        self.coordinates = None
        self.run_count = 0

    def set(self, name, value):
        self.settings[name] = value

    def setInterpolationCoordinates(self, coordinates):
        self.coordinates = np.asarray(coordinates, dtype=float)

    def run(self):
        self.run_count += 1

    def getData(self):
        return self.real

    def getImagData(self):
        return self.imag


class FakeNumericalRoot:
    def __init__(self, interp):
        self.interp = interp

    def __call__(self, tag):
        assert tag == "int_grid"
        return self.interp


class FakeResults:
    def __init__(self, interp):
        self._numerical = FakeNumericalRoot(interp)

    def numerical(self, tag=None):
        return self._numerical if tag is None else self._numerical(tag)


class FakeJavaModel:
    def __init__(self, interp):
        self._result = FakeResults(interp)

    def result(self):
        return self._result


def make_fake_run(real, imag):
    interp = FakeInterp(real, imag)
    run = SimulationRun.__new__(SimulationRun)
    run.client = type("Client", (), {"remove": lambda self, model: None})()
    run.model = type("Model", (), {"java": FakeJavaModel(interp)})()
    return run, interp


def disable_jpype_conversion(monkeypatch):
    monkeypatch.setattr(simulation_utils, "JString", str)
    monkeypatch.setattr(simulation_utils, "JInt", int)
    monkeypatch.setattr(simulation_utils, "JDouble", float)
    monkeypatch.setattr(
        simulation_utils,
        "JArray",
        lambda _value_type, _dimensions: lambda values: values,
    )


def test_get_fields_at_coordinates_sets_grid_and_returns_complex_values(monkeypatch):
    disable_jpype_conversion(monkeypatch)
    run, interp = make_fake_run(
        real=[[[1.0, 2.0]], [[3.0, 4.0]]],
        imag=[[[0.5, 0.0]], [[0.0, -0.5]]],
    )
    coordinates = np.array([[0.0, 0.0, 1.2], [1.0, 0.0, 1.2]])

    values = run.get_fields_at_coordinates(
        0,
        ["ewfd.Ex", "ewfd.Ey"],
        coordinates,
    )

    assert interp.settings["data"] == "dset1"
    assert interp.settings["expr"] == ["ewfd.Ex", "ewfd.Ey"]
    assert "innerinput" not in interp.settings
    assert interp.settings["solnum"] == [1]
    assert interp.coordinates == pytest.approx(coordinates.T)
    assert interp.run_count == 1
    assert values.shape == (2, 2)
    assert values[0, 0] == pytest.approx(1.0 + 0.5j)
    assert values[1, 1] == pytest.approx(4.0 - 0.5j)


def test_get_fields_at_coordinates_rejects_wrong_coordinate_shape(monkeypatch):
    disable_jpype_conversion(monkeypatch)
    run, _interp = make_fake_run(real=[], imag=[])

    with pytest.raises(ValueError, match="shape"):
        run.get_fields_at_coordinates(0, ["ewfd.Ex"], np.zeros((3, 2)))
