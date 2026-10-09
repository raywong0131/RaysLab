import json

import numpy as np
import pandas as pd
import pytest

from scripts.run_main import run_finite


class FakeParameters:
    def __init__(self, air_plane_z_um):
        self.air_plane_z_um = air_plane_z_um
        self.calls = []

    def evaluate(self, expression, unit):
        self.calls.append((expression, unit))
        return self.air_plane_z_um


class FakeJavaModel:
    def __init__(self, parameters):
        self._parameters = parameters

    def param(self):
        return self._parameters


class FakeSimulation:
    def __init__(self, air_plane_z_um=1.75):
        self.parameters = FakeParameters(air_plane_z_um)
        self.model = type("Model", (), {"java": FakeJavaModel(self.parameters)})()
        self.calls = []

    def get_fields_at_coordinates(self, mode_idx, expressions, coordinates, dataset="dset1"):
        self.calls.append((mode_idx, list(expressions), coordinates.copy(), dataset))
        values = np.vstack(
            [
                np.ones(len(coordinates), dtype=complex),
                np.full(len(coordinates), 2.0 + 3.0j, dtype=complex),
                np.full(len(coordinates), -1.0j, dtype=complex),
            ]
        )
        values[:, 0] = np.nan + 1j * np.nan
        return values


def test_farfield_defaults_and_air_plane():
    assert run_finite.RUN_FARFIELD_FFT is True
    assert run_finite.FARFIELD_GRID_SIZE == 801
    assert run_finite.FARFIELD_FFT_SIZE == 2001
    assert run_finite.FARFIELD_NA == pytest.approx(0.9)
    assert run_finite.FARFIELD_H0_UM is None
    assert run_finite.finite_simulation_config().air_cutplane_z == "H_air"


def test_farfield_square_grid_is_centered_minimum_hex_square():
    x_axis, y_axis, coordinates = run_finite.farfield_square_grid(34.0, 801, 1.75)

    assert x_axis[[0, 400, 800]].tolist() == [-17.0, 0.0, 17.0]
    assert y_axis[[0, 400, 800]].tolist() == [-17.0, 0.0, 17.0]
    assert coordinates.shape == (801 * 801, 3)
    assert coordinates[400 * 801 + 400].tolist() == [0.0, 0.0, 1.75]


def test_export_valid_air_fields_exports_only_valid_modes_and_zeros_domain_exterior(tmp_path):
    sim = FakeSimulation()
    eigen_df = pd.DataFrame(
        {
            "re": [194.8, 195.1],
            "im": [0.1, 0.2],
            "q": [974.0, 488.0],
            "is_valid": [True, False],
            "mode_idx": [0, 1],
        }
    )

    metadata = run_finite.export_valid_air_fields(
        sim,
        tmp_path,
        eigen_df,
        finite_hex_a=2.0,
        grid_size=3,
        h0_um=None,
    )

    assert len(sim.calls) == 1
    assert sim.calls[0][0] == 0
    assert sim.calls[0][1] == ["ewfd.Ex", "ewfd.Ey", "ewfd.Ez"]
    assert not (tmp_path / "01_E_air.parquet").exists()
    field = pd.read_parquet(tmp_path / "00_E_air.parquet")
    assert len(field) == 9
    assert field.loc[0, "is_in_domain"] == False  # noqa: E712
    assert field.loc[0, ["Ex_re", "Ex_im", "Ey_re", "Ey_im", "Ez_re", "Ez_im"]].eq(0.0).all()
    assert field.loc[1, "is_in_domain"] == False  # noqa: E712
    assert field.loc[1, ["Ex_re", "Ex_im", "Ey_re", "Ey_im", "Ez_re", "Ez_im"]].eq(0.0).all()
    assert int(field["is_in_domain"].sum()) == 3
    assert field.loc[4, ["x", "y"]].tolist() == [0.0, 0.0]
    assert metadata["air_plane_z_um"] == pytest.approx(1.75)
    saved = json.loads((tmp_path / "air_field_metadata.json").read_text(encoding="utf-8"))
    assert saved["grid_size"] == 3
    assert saved["valid_mode_indices"] == [0]
    assert sim.parameters.calls == [("z_slab_top + H_air", "um")]


def test_export_valid_air_fields_uses_square_rectangle_mask(tmp_path):
    sim = FakeSimulation()
    eigen_df = pd.DataFrame(
        {"re": [195.0], "im": [0.1], "q": [900.0], "is_valid": [True], "mode_idx": [0]}
    )
    metadata = run_finite.export_valid_air_fields(
        sim,
        tmp_path,
        eigen_df,
        finite_hex_a=4.0,
        grid_size=5,
        h0_um=None,
        finite_geometry="square",
        footprint_width=4.0,
        footprint_height=2.0,
    )
    field = pd.read_parquet(tmp_path / "00_E_air.parquet")
    assert int(field["is_in_domain"].sum()) == 15
    assert metadata["domain_mask_kind"] == "rectangle_footprint_v1"
    assert metadata["sampling_window_um"] == pytest.approx(4.0)
    assert "finite_hex_a_um" not in metadata
