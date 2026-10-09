"""Configuration-only checks: no COMSOL or formal result writes."""

from copy import deepcopy
from dataclasses import replace
import importlib.util
import json
from pathlib import Path

import pytest

from scripts.run_main import parameter_config as config
from scripts.run_sweep.optimize_finite_quarter_q import (
    Candidate, build_candidate_parameter_data,
)


def current_data():
    return json.loads(config.DEFAULT_PARAMETER_PATH.read_text(encoding="utf-8"))


def save(tmp_path, data):
    path = tmp_path / "parameter.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_v2_and_legacy_snapshot_resolve_same_runtime_contract(tmp_path):
    data = current_data()
    path = save(tmp_path, data)
    grouped = config.load_shared_parameters(path)
    legacy = {
        "schema_version": 1,
        "mesh_size": data["mesh_size"],
        **data["cells"], **data["structure"],
        "unit_cell_eigenmode_count": data["unit_cell_band"]["eigenmode_count"],
        "finite_geometry": data["finite_cavity"]["geometry"],
        "finite_eigenmode_count": data["finite_cavity"]["eigenmode_count"],
        "unit_cell_2d": data["unit_cell_2d"],
        **{key: value for key, value in data["finite_cavity"].items()
           if key.startswith("cladding_")},
        "strip_eigenmode_count": data["strip_1d"]["eigenmode_count"],
        "strip_cut_angle_deg": data["strip_1d"]["cut_angle_deg"],
        "strip_kx": data["strip_1d"]["kx"],
        "finite_quarter_symmetry_ids": data["finite_cavity"]["quarter_symmetry_ids"],
    }
    assert config.load_shared_parameters(save(tmp_path, legacy)) == grouped


@pytest.mark.parametrize(("section", "key", "value", "message"), [
    (None, "cavity", {}, "unknown fields"),
    (None, "schema_version", True, "schema_version"),
    (None, "schema_version", 99, "schema_version"),
    (None, "schema_version", 1, "v2 sections"),
    (None, "unit_cell_2d", None, "unit_cell_2d must be an object"),
    ("unit_cell_band", "eigenmode_count", 0, "positive integer"),
    ("unit_cell_band", "eigenmode_counts", 8, "unknown fields"),
    ("finite_cavity", "geometry", "triangle", "finite_geometry"),
    ("finite_cavity", "eigenmode_count", True, "positive integer"),
    ("finite_cavity", "quarter_symmetry_ids", [], "quarter_symmetry_ids"),
    ("finite_cavity", "quarter_symmetry_ids", [1, 1], "quarter_symmetry_ids"),
    ("finite_cavity", "quarter_symmetry_ids", [True], "quarter_symmetry_ids"),
    ("finite_cavity", "quarter_symmetry_ids", [5], "quarter_symmetry_ids"),
    ("strip_1d", "eigenmode_count", 0, "strip_1d.eigenmode_count"),
    ("strip_1d", "eigenmode_count", 2.5, "strip_1d.eigenmode_count"),
    ("strip_1d", "cut_angle_deg", 30, "cut_angle_deg"),
    ("strip_1d", "kx", float("nan"), "must be finite"),
    ("unit_cell_2d", "valley_refinement", {"maximum_new_point": 1}, "unknown fields"),
])
def test_invalid_v2_inputs_fail_before_any_solver(tmp_path, section, key, value, message):
    data = current_data()
    (data if section is None else data[section])[key] = value
    with pytest.raises((TypeError, ValueError), match=message):
        config.load_shared_parameters(save(tmp_path, data))


def test_v2_explicit_xy_stays_finite_only_and_rejects_conflicting_ratio(tmp_path):
    data = current_data()
    finite = data["finite_cavity"]
    finite.pop("cladding_y_over_x_shift_ratio", None)
    finite.update(cladding_x_shift_factor=0.02, cladding_y_shift_factor=0.04)
    parsed = config.load_shared_parameters(save(tmp_path, data))
    assert parsed.finite_cladding_shift_factor_pairs == ((0.02, 0.04),)
    assert list(parsed.cladding_shift_factors) == data["structure"]["cladding_shift_factors"]
    finite["cladding_y_over_x_shift_ratio"] = 1.0
    with pytest.raises(ValueError, match="conflicts"):
        config.load_shared_parameters(save(tmp_path, data))


def test_v2_frequency_update_preserves_other_sections_and_concurrent_layers(tmp_path):
    data = current_data()
    path = save(tmp_path, data)
    expected = config.load_shared_parameters(path)
    data["structure"]["cavity_layers"] += 1
    save(tmp_path, data)
    frequency = expected.center_frequency_thz + 0.1
    config.update_center_frequency_from_cavity_p2(
        frequency, expected_parameters=expected, path=path,
    )
    data["structure"].update(center_frequency_thz=frequency, center_frequency_source={
        "workflow": "run_band_pair", "cell": "cavity", "band": "p2", "k_point": "Gamma",
    })
    assert json.loads(path.read_text(encoding="utf-8")) == data
    assert not (tmp_path / ".out").exists()


@pytest.mark.parametrize("section,key", [("cells", "cavity"), ("unit_cell_band", "eigenmode_count")])
def test_v2_frequency_update_rejects_changed_band_identity(tmp_path, section, key):
    data = current_data()
    path = save(tmp_path, data)
    expected = config.load_shared_parameters(path)
    if section == "cells":
        data[section][key]["b0_nm"] += 1
    else:
        data[section][key] += 1
    save(tmp_path, data)
    before = path.read_bytes()
    with pytest.raises(RuntimeError, match="refusing to overwrite"):
        config.update_center_frequency_from_cavity_p2(
            expected.center_frequency_thz + 0.1, expected_parameters=expected, path=path,
        )
    assert path.read_bytes() == before


def test_v2_optimizer_only_updates_candidate_snapshot(tmp_path):
    baseline = current_data()
    original = deepcopy(baseline)
    candidate = Candidate(244.2, 0.98, 0.932)
    updated = build_candidate_parameter_data(baseline, candidate)
    assert baseline == original
    expected = deepcopy(original)
    expected["cells"]["cladding"] = {"b0_nm": 244.2, "eta": 0.98, "zeta": 0.932}
    expected["structure"].update(cladding_shift_factors=[0.0], cladding_shift_profile={"kind": "uniform"})
    assert updated == expected
    assert config.load_shared_parameters(save(tmp_path, updated)).cladding.b0_nm == 244.2


def test_primary_runners_bind_strip_and_quarter_configuration(monkeypatch):
    from scripts.run_main import run_finite

    parameters = replace(config.ACTIVE_PARAMETERS, strip_eigenmode_count=7,
                         strip_cut_angle_deg=0.0, strip_kx=0.125,
                         finite_quarter_symmetry_ids=(2, 4))
    monkeypatch.setattr(config, "ACTIVE_PARAMETERS", parameters)
    monkeypatch.setattr(run_finite, "ACTIVE_PARAMETERS", parameters)
    modules = []
    for name in ("run_strip_1d", "run_finite_quarter"):
        path = Path(config.__file__).with_name(name + ".py")
        spec = importlib.util.spec_from_file_location("parameter_test_" + name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        modules.append(module)
    strip, quarter = modules
    assert (strip.EIGENMODE_COUNT, strip.STRIP_CUT_ANGLE, strip.STRIP_KX) == (7, 0.0, 0.125)
    assert quarter.SYMMETRY_IDS == [2, 4]
    assert run_finite.EIGENMODE_COUNT == parameters.finite_eigenmode_count
