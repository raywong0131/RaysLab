import copy
import os
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import pytest

from comsol_workflow.cladding_shift_profile import (
    CladdingShiftProfile,
    directional_shift_profile_metadata,
)
from scripts.run_main import run_finite
from scripts.run_main import run_strip_1d


def _resume_metadata():
    return {
        **directional_shift_profile_metadata(
            20,
            0.06,
            0.06,
            0.82,
            CladdingShiftProfile(),
        ),
        "a": 0.82,
        "bulk_radius": 20,
        "cladding_layers": 20,
        "cladding_x_shift_factor": 0.06,
        "cladding_y_shift_factor": 0.06,
        "cladding_x_inward_shift": 0.0492,
        "cladding_y_inward_shift": 0.0492,
        "finite_cavity_preset": "preset",
        "finite_cavity_preset_config": {
            "cavity": {"b0_nm": 235.0, "eta": 0.96, "zeta": 1.156},
            "cladding": {"b0_nm": 231.8, "eta": 0.96, "zeta": 0.83},
        },
        "shared_parameters": {"source": "test"},
        "bulk_params": {"r_f0": 0.96, "b_square_f0": 1.0},
        "cladding_params": {"r_f0": 0.96, "b_square_f0": 1.0},
        "simulation_holes": 29526,
        "simulation_cells": 4921,
        "finite_hex_a": 67.51333333333334,
        "finite_Lx": 67.51333333333334,
        "finite_Ly": 58.4682617608334,
        "simulation_common": {
            "slab_height": "200 [nm]",
            "refractive_index": "3.3",
            "wavelength": "1550 [nm]",
            "eigenfrequency_shift": "195.950356405206 [THz]",
            "eigenmode_count": 2,
            "mesh_auto_size": 9,
            "use_manual_mesh_construction": True,
            "mode_type": "TE",
        },
        "simulation": {
            "mode": "finite_size",
            "boundary": "hexagon_boundary",
        },
    }


def test_resume_metadata_accepts_matching_checkpoint():
    metadata = _resume_metadata()

    run_finite.validate_resume_metadata(metadata, metadata)


def test_resume_metadata_reports_nested_mismatch():
    saved = _resume_metadata()
    current = copy.deepcopy(saved)
    current["simulation_common"]["mesh_auto_size"] = 8

    with pytest.raises(
        ValueError,
        match=r"simulation_common.*saved=.*mesh_auto_size.*9.*current=.*mesh_auto_size.*8",
    ):
        run_finite.validate_resume_metadata(saved, current)


def test_resume_metadata_reports_missing_field():
    saved = _resume_metadata()
    del saved["finite_cavity_preset_config"]

    with pytest.raises(ValueError, match="finite_cavity_preset_config"):
        run_finite.validate_resume_metadata(saved, _resume_metadata())


def test_resume_metadata_reports_cladding_shift_profile_mismatch():
    saved = _resume_metadata()
    current = copy.deepcopy(saved)
    current["cladding_shift_profile"] = {
        "kind": "tanh_power",
        "scale_layers": 4.0,
        "power": 2.0,
    }

    with pytest.raises(ValueError, match="cladding_shift_profile"):
        run_finite.validate_resume_metadata(saved, current)


class _FakeTags:
    def __init__(self, tags):
        self._tags = tags

    def tags(self):
        return self._tags


class _FakeCheckpointGeometry:
    def __init__(self, feature_tags):
        self._features = _FakeTags(feature_tags)

    def feature(self):
        return self._features


class _FakeCheckpointMesh(_FakeCheckpointGeometry):
    pass


class _FakeCheckpointComponent:
    def __init__(self, geometry_tags, mesh_tags):
        self._geometry = _FakeCheckpointGeometry(geometry_tags)
        self._mesh = _FakeCheckpointMesh(mesh_tags)

    def geom(self, tag):
        if tag != "geom1":
            raise KeyError(tag)
        return self._geometry

    def mesh(self, tag):
        if tag != "mesh1":
            raise KeyError(tag)
        return self._mesh


class _FakeCheckpointResult:
    def __init__(self, dataset_tags):
        self._datasets = _FakeTags(dataset_tags)

    def dataset(self):
        return self._datasets


class _FakeCheckpointJava:
    def __init__(
        self,
        *,
        geometry_tags,
        mesh_tags,
        study_tags,
        solution_tags,
        dataset_tags,
    ):
        self._component = _FakeCheckpointComponent(geometry_tags, mesh_tags)
        self._studies = _FakeTags(study_tags)
        self._solutions = _FakeTags(solution_tags)
        self._result = _FakeCheckpointResult(dataset_tags)

    def component(self, tag):
        if tag != "comp1":
            raise KeyError(tag)
        return self._component

    def study(self):
        return self._studies

    def sol(self):
        return self._solutions

    def result(self):
        return self._result


def _checkpoint_simulation_run(
    *,
    mesh_tags,
    solution_tags,
    dataset_tags,
):
    java = _FakeCheckpointJava(
        geometry_tags=["blk1"],
        mesh_tags=mesh_tags,
        study_tags=["std1"],
        solution_tags=solution_tags,
        dataset_tags=dataset_tags,
    )
    return types.SimpleNamespace(model=types.SimpleNamespace(java=java))


@pytest.mark.parametrize(
    ("mesh_tags", "solution_tags", "dataset_tags", "message"),
    [
        (["size", "ftet1"], [], [], "mesh feature tags"),
        (["size"], ["sol1"], [], "solution tags"),
        (["size"], [], ["dset1"], "dataset tags"),
    ],
)
def test_geometry_checkpoint_rejects_non_geometry_only_state(
    mesh_tags,
    solution_tags,
    dataset_tags,
    message,
):
    sim = _checkpoint_simulation_run(
        mesh_tags=mesh_tags,
        solution_tags=solution_tags,
        dataset_tags=dataset_tags,
    )

    with pytest.raises(ValueError, match=message):
        run_finite.validate_geometry_checkpoint_model(sim)


def test_geometry_checkpoint_accepts_expected_state():
    sim = _checkpoint_simulation_run(
        mesh_tags=["size"],
        solution_tags=[],
        dataset_tags=[],
    )

    run_finite.validate_geometry_checkpoint_model(sim)


class _FakeBoundarySpec:
    def __init__(self, shape, a, b=None):
        self.shape = shape
        self.a = a
        self.b = b


class _FakeLayerSpec:
    def __init__(self, height, holes, refractive_index, label=None):
        self.height = height
        self.holes = holes
        self.refractive_index = refractive_index
        self.label = label


class _FakeSimulationConfig:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class _FakeModel:
    def __init__(self, events):
        self._events = events
        self.java = object()

    def save(self, path):
        self._events.append(("save", Path(path).name))


class _FakeSimulationRun:
    events = []

    def __init__(self, config=None, progress_log_path=None):
        self.config = config
        self.progress_log_path = progress_log_path
        self.model = _FakeModel(self.events)

    def __enter__(self):
        self.events.append(("enter", Path(self.progress_log_path).name))
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.events.append(("exit", None))
        return False

    def build_geometry(self, boundary, layers, simulation_mode="unit_cell", k=None):
        self.events.append(("build_geometry", simulation_mode))

    def run_simulation(self, mesh_save_path=None, mesh_builder=None):
        self.events.append(("run_simulation", None if mesh_save_path is None else Path(mesh_save_path).name))
        if mesh_builder is not None:
            mesh_builder(self.model.java, self.config)
        if mesh_save_path is not None:
            self.model.save(mesh_save_path)

    def build_and_run(self, boundary, layers, simulation_mode="unit_cell", k=None):
        self.events.append(("build_and_run", simulation_mode))


def _fake_simulation_utils_module():
    return types.SimpleNamespace(
        BoundarySpec=_FakeBoundarySpec,
        LayerSpec=_FakeLayerSpec,
        SimulationConfig=_FakeSimulationConfig,
        SimulationRun=_FakeSimulationRun,
    )


def _patch_attr(module, originals, name, value):
    originals[name] = getattr(module, name)
    setattr(module, name, value)


def _run_finite_full_case_with_fakes(
    module,
    *,
    manual_mesh=False,
    resume=False,
    create_resume_mph=True,
    resume_metadata_matches=True,
    resume_model_valid=True,
):
    originals = {}
    simulation_module_name = "comsol_workflow.simulation_utils"
    mesh_module_name = "comsol_workflow.mesh_constrcution"
    original_simulation_utils = sys.modules.get(simulation_module_name)
    original_mesh_constrcution = sys.modules.get(mesh_module_name)
    _FakeSimulationRun.events = []

    _patch_attr(module, originals, "full_cell_records", lambda _metadata: [{"polygon": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])}])
    _patch_attr(
        module,
        originals,
        "simulation_boundary_for_full_lattice_hexagon",
        lambda _cells: (np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]), 2.0, 2.0, 1.0),
    )
    _patch_attr(module, originals, "save_full_lattice_plot", lambda *args, **kwargs: None)
    _patch_attr(module, originals, "save_full_lattice_shift_plot", lambda *args, **kwargs: None)
    _patch_attr(module, originals, "save_full_lattice_normal_fan_plot", lambda *args, **kwargs: None)
    _patch_attr(module, originals, "save_finite_simulation_plot", lambda *args, **kwargs: None)
    if hasattr(module, "save_finite_model_geometry_figures"):
        _patch_attr(
            module,
            originals,
            "save_finite_model_geometry_figures",
            lambda model_dir, *_args, **_kwargs: _FakeSimulationRun.events.append(
                ("save_finite_model_geometry_figures", Path(model_dir).name)
            ),
        )
    for diagnostic_name in (
        "save_shift_profile_diagnostics",
        "save_directional_shift_profile_diagnostics",
    ):
        if hasattr(module, diagnostic_name):
            _patch_attr(
                module,
                originals,
                diagnostic_name,
                lambda *_args, **_kwargs: None,
            )
    _patch_attr(
        module,
        originals,
        "write_config",
        lambda path, *_args, **_kwargs: _FakeSimulationRun.events.append(
            ("write_config", Path(path).name)
        ),
    )
    _patch_attr(module, originals, "bulk_points", lambda: [])

    if hasattr(module, "export_full_finite_results"):
        _patch_attr(module, originals, "export_full_finite_results", lambda *args, **kwargs: None)
    if hasattr(module, "export_simulation_results"):
        _patch_attr(module, originals, "export_simulation_results", lambda *args, **kwargs: None)
    if hasattr(module, "resolve_finite_analysis_mode_indices"):
        _patch_attr(
            module,
            originals,
            "resolve_finite_analysis_mode_indices",
            lambda *_args, **_kwargs: (
                _FakeSimulationRun.events.append(
                    ("resolve_finite_analysis_mode_indices", [1])
                )
                or [1]
            ),
        )
    if hasattr(module, "run_exported_finite_lattice_fourier_postprocess"):
        _patch_attr(
            module,
            originals,
            "run_exported_finite_lattice_fourier_postprocess",
            lambda *args, **kwargs: _FakeSimulationRun.events.append(
                ("run_exported_finite_lattice_fourier_postprocess", kwargs)
            ),
        )
    if hasattr(module, "run_finite_lattice_fourier_postprocess"):
        _patch_attr(module, originals, "run_finite_lattice_fourier_postprocess", lambda *args, **kwargs: None)
    if hasattr(module, "run_farfield_fft_case"):
        _patch_attr(
            module,
            originals,
            "run_farfield_fft_case",
            lambda *args, **kwargs: _FakeSimulationRun.events.append(
                (
                    "run_farfield_fft_case",
                    (
                        kwargs["grid_size"],
                        kwargs["fft_size"],
                        kwargs["na"],
                        kwargs["h0_um"],
                        kwargs["period_um"],
                        kwargs.get("mode_indices"),
                    ),
                )
            ),
        )
    if hasattr(module, "score_exported_finite_modes"):
        _patch_attr(module, originals, "score_exported_finite_modes", lambda *args, **kwargs: None)
    if hasattr(module, "score_finite_modes"):
        _patch_attr(module, originals, "score_finite_modes", lambda *args, **kwargs: None)
    if hasattr(module, "finalize_finite_case_output"):
        _patch_attr(
            module,
            originals,
            "finalize_finite_case_output",
            lambda case_dir, *, model_filename: module.finite_case_output_paths(
                Path(case_dir),
                model_filename=model_filename,
            ),
        )
    if hasattr(module, "USE_MANUAL_MESH_CONSTRUCTION"):
        _patch_attr(module, originals, "USE_MANUAL_MESH_CONSTRUCTION", manual_mesh)
    if hasattr(module, "RESUME_FROM_EXISTING_MPH"):
        _patch_attr(module, originals, "RESUME_FROM_EXISTING_MPH", resume)
        _patch_attr(
            module,
            originals,
            "load_and_validate_resume_metadata",
            lambda *_args, **_kwargs: (
                _FakeSimulationRun.events.append(("validate_resume_metadata", None))
                if resume_metadata_matches
                else (_ for _ in ()).throw(ValueError("metadata mismatch"))
            ),
        )
        _patch_attr(
            module,
            originals,
            "attach_saved_model",
            lambda _sim, path: _FakeSimulationRun.events.append(
                ("attach_saved_model", Path(path).name)
            ),
        )
        _patch_attr(
            module,
            originals,
            "validate_geometry_checkpoint_model",
            lambda _sim: (
                _FakeSimulationRun.events.append(
                    ("validate_geometry_checkpoint_model", None)
                )
                if resume_model_valid
                else (_ for _ in ()).throw(ValueError("invalid checkpoint model"))
            ),
        )

    try:
        sys.modules[simulation_module_name] = _fake_simulation_utils_module()
        if manual_mesh:
            sys.modules[mesh_module_name] = types.SimpleNamespace(
                construct_manual_mesh=lambda *_args, **kwargs: _FakeSimulationRun.events.append(
                    ("construct_manual_mesh", kwargs.get("mesh_auto_size"))
                )
            )
        with tempfile.TemporaryDirectory() as tmpdir:
            if resume and create_resume_mph:
                model_dir = Path(tmpdir) / run_finite.MODEL_DIRNAME
                model_dir.mkdir(parents=True)
                (model_dir / "finite_cavity.mph").touch()
            records = [{"polygon": np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]])}]
            validity_boundary = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
            metadata = {
                "bulk_points": [],
                "cladding_points": [],
                "inner_boundary": validity_boundary,
                "a": 0.82,
                "bulk_radius": 1,
                "cladding_layers": 1,
            }
            module.run_finite_full_case(Path(tmpdir), records, metadata, validity_boundary)
        return list(_FakeSimulationRun.events)
    finally:
        for name, value in originals.items():
            setattr(module, name, value)
        if original_simulation_utils is None:
            sys.modules.pop(simulation_module_name, None)
        else:
            sys.modules[simulation_module_name] = original_simulation_utils
        if original_mesh_constrcution is None:
            sys.modules.pop(mesh_module_name, None)
        else:
            sys.modules[mesh_module_name] = original_mesh_constrcution


def _assert_single_mph_saved_at_three_stages(events):
    assert ("build_geometry", "finite_size") in events, events
    assert ("run_simulation", "finite_cavity.mph") in events, events

    save_events = [value for event, value in events if event == "save"]
    assert save_events == ["finite_cavity.mph", "finite_cavity.mph", "finite_cavity.mph"], events
    assert "finite_cavity_geometry.mph" not in save_events
    assert "finite_cavity_solved.mph" not in save_events

    geometry_save = events.index(("save", "finite_cavity.mph"))
    run = events.index(("run_simulation", "finite_cavity.mph"))
    mesh_save = events.index(("save", "finite_cavity.mph"), run)
    final_save = events.index(("save", "finite_cavity.mph"), mesh_save + 1)

    assert events.index(("build_geometry", "finite_size")) < geometry_save
    assert geometry_save < run
    assert run < mesh_save
    assert mesh_save < final_save


def _run_strip_cut_case_with_fakes():
    originals = {}
    simulation_module_name = "comsol_workflow.simulation_utils"
    original_simulation_utils = sys.modules.get(simulation_module_name)
    _FakeSimulationRun.events = []

    _patch_attr(run_strip_1d, originals, "strip_records_for_cut", lambda records, _angle: records)
    _patch_attr(run_strip_1d, originals, "metadata_to_strip_frame", lambda metadata, _angle: metadata)
    _patch_attr(run_strip_1d, originals, "points_to_strip_frame", lambda points, _angle: points)
    _patch_attr(run_strip_1d, originals, "strip_cell_records", lambda _metadata: [{"polygon": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])}])
    _patch_attr(
        run_strip_1d,
        originals,
        "simulation_boundary_for_strip",
        lambda _cells: (np.array([[-0.5, -1.0], [0.5, -1.0], [0.5, 1.0], [-0.5, 1.0]]), 1.0, 2.0),
    )
    _patch_attr(run_strip_1d, originals, "save_global_cut_preview", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "save_strip_plot", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "save_simulation_plot", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "write_config", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "export_simulation_results", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "run_strip_bulk_fourier_postprocess", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "score_strip_modes", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "finalize_strip_case_output", lambda *args, **kwargs: None)
    _patch_attr(run_strip_1d, originals, "strip_bulk_cell_records", lambda *_args: [])

    try:
        sys.modules[simulation_module_name] = _fake_simulation_utils_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            records = [{"polygon": np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]])}]
            validity_boundary = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
            metadata = {
                "bulk_points": [],
                "cladding_points": [],
                "inner_boundary": validity_boundary,
                "a": 0.82,
                "bulk_radius": 1,
                "cladding_layers": 1,
            }
            run_strip_1d.run_strip_cut_case(Path(tmpdir), records, metadata, validity_boundary)
        return list(_FakeSimulationRun.events)
    finally:
        for name, value in originals.items():
            setattr(run_strip_1d, name, value)
        if original_simulation_utils is None:
            sys.modules.pop(simulation_module_name, None)
        else:
            sys.modules[simulation_module_name] = original_simulation_utils


def _assert_strip_mph_saved_at_three_stages(events):
    assert ("build_geometry", "strip_1d") in events, events
    assert ("run_simulation", "strip_1d.mph") in events, events

    save_events = [value for event, value in events if event == "save"]
    assert save_events == ["strip_1d.mph", "strip_1d.mph", "strip_1d.mph"], events
    assert "strip_1d_solved.mph" not in save_events

    geometry_save = events.index(("save", "strip_1d.mph"))
    run = events.index(("run_simulation", "strip_1d.mph"))
    mesh_save = events.index(("save", "strip_1d.mph"), run)
    final_save = events.index(("save", "strip_1d.mph"), mesh_save + 1)

    assert events.index(("build_geometry", "strip_1d")) < geometry_save
    assert geometry_save < run
    assert run < mesh_save
    assert mesh_save < final_save


def test_cladding_inward_shift_factors_come_from_shared_parameters():
    assert run_strip_1d.CLADDING_INWARD_SHIFT_FACTORS == list(
        run_strip_1d.ACTIVE_PARAMETERS.cladding_shift_factors
    )


def test_strip_structure_layers_come_from_shared_parameters():
    assert (
        run_strip_1d.BULK_RADIUS
        == run_strip_1d.ACTIVE_PARAMETERS.cavity_layers - 1
    )
    assert (
        run_strip_1d.CLADDING_LAYERS
        == run_strip_1d.ACTIVE_PARAMETERS.cladding_layers
    )


def test_shift_factor_label_always_uses_three_decimals():
    assert run_strip_1d.shift_factor_label(0.00) == "0.000"
    assert run_strip_1d.shift_factor_label(0.07) == "0.070"
    assert run_strip_1d.shift_factor_label(0.075) == "0.075"
    assert run_strip_1d.shift_factor_label(0.09) == "0.090"


def test_output_dir_for_supplemental_shift_is_distinct_from_point_zero_seven():
    supplemental = run_strip_1d.output_dir_for_cladding_shift_factor(0.075)
    existing = run_strip_1d.output_dir_for_cladding_shift_factor(0.07)

    assert supplemental.parent.name == run_strip_1d.strip_cut_dirname()
    assert supplemental.name == "shift0.075"
    assert existing.name == "shift0.070"
    assert supplemental != existing


def test_run_finite_shift_label_always_uses_three_decimals():
    assert run_finite.shift_factor_label(0.01) == "0.010"
    assert run_finite.shift_factor_label(0.02) == "0.020"
    assert run_finite.shift_factor_label(0.05) == "0.050"
    assert run_finite.shift_factor_label(0.07) == "0.070"
    assert run_finite.shift_factor_label(0.075) == "0.075"
    case_stem = run_finite.case_stem(0.075)
    assert case_stem == "shift0.075"
    assert run_finite.case_stem(0.075) != run_finite.case_stem(0.07)


def test_run_strip_finite_full_case_overwrites_one_mph_at_three_stages():
    events = _run_finite_full_case_with_fakes(run_strip_1d)

    _assert_single_mph_saved_at_three_stages(events)


def test_run_finite_full_case_overwrites_one_mph_at_three_stages():
    events = _run_finite_full_case_with_fakes(run_finite)

    _assert_single_mph_saved_at_three_stages(events)
    assert ("save_finite_model_geometry_figures", "00_model") in events
    assert (
        "run_farfield_fft_case",
        (801, 2001, 0.9, None, 0.82, [1]),
    ) in events
    lattice_events = [
        value
        for event, value in events
        if event == "run_exported_finite_lattice_fourier_postprocess"
    ]
    assert lattice_events[0]["mode_indices"] == [1]
    assert lattice_events[0]["include_unselected_valid_mode_maps"] is False


def test_finite_single_mode_analysis_can_be_disabled(monkeypatch, tmp_path):
    monkeypatch.setattr(
        run_finite,
        "FINITE_SINGLE_MODE_ANALYSIS_ENABLED",
        False,
    )
    monkeypatch.setattr(
        run_finite,
        "select_finite_mode_closest_to_frequency",
        lambda *_args, **_kwargs: pytest.fail("selector should not run"),
    )

    assert run_finite.resolve_finite_analysis_mode_indices(
        tmp_path,
        tmp_path,
    ) is None


def test_run_finite_full_case_can_use_manual_mesh_construction():
    events = _run_finite_full_case_with_fakes(run_finite, manual_mesh=True)

    mesh_event = ("construct_manual_mesh", run_finite.MESH_AUTO_SIZE)
    assert mesh_event in events, events
    assert events.index(("run_simulation", "finite_cavity.mph")) < events.index(mesh_event)


def test_run_finite_resume_validates_before_write_and_skips_geometry():
    events = _run_finite_full_case_with_fakes(
        run_finite,
        manual_mesh=True,
        resume=True,
    )

    assert ("validate_resume_metadata", None) in events
    assert ("attach_saved_model", "finite_cavity.mph") in events
    assert ("validate_geometry_checkpoint_model", None) in events
    assert ("build_geometry", "finite_size") not in events
    assert ("construct_manual_mesh", run_finite.MESH_AUTO_SIZE) in events
    assert ("run_simulation", "finite_cavity.mph") in events

    metadata_validation = events.index(("validate_resume_metadata", None))
    model_validation = events.index(("validate_geometry_checkpoint_model", None))
    config_write = events.index(("write_config", "config.json"))
    run = events.index(("run_simulation", "finite_cavity.mph"))
    assert metadata_validation < model_validation < config_write < run


def test_run_finite_resume_requires_existing_mph():
    with pytest.raises(FileNotFoundError, match="missing geometry checkpoint"):
        _run_finite_full_case_with_fakes(
            run_finite,
            resume=True,
            create_resume_mph=False,
        )


def test_run_finite_resume_stops_on_metadata_mismatch():
    with pytest.raises(ValueError, match="metadata mismatch"):
        _run_finite_full_case_with_fakes(
            run_finite,
            resume=True,
            resume_metadata_matches=False,
        )

    assert not any(
        event in {"attach_saved_model", "run_simulation"}
        for event, _value in _FakeSimulationRun.events
    )


def test_run_finite_resume_stops_on_invalid_model_state():
    with pytest.raises(ValueError, match="invalid checkpoint model"):
        _run_finite_full_case_with_fakes(
            run_finite,
            resume=True,
            resume_model_valid=False,
        )

    assert not any(
        event == "run_simulation"
        for event, _value in _FakeSimulationRun.events
    )


def test_run_strip_cut_case_overwrites_one_mph_at_three_stages():
    events = _run_strip_cut_case_with_fakes()

    _assert_strip_mph_saved_at_three_stages(events)


if __name__ == "__main__":
    test_cladding_inward_shift_factors_come_from_shared_parameters()
    test_shift_factor_label_always_uses_three_decimals()
    test_output_dir_for_supplemental_shift_is_distinct_from_point_zero_seven()
    test_run_finite_shift_label_always_uses_three_decimals()
    test_run_strip_finite_full_case_overwrites_one_mph_at_three_stages()
    test_run_finite_full_case_overwrites_one_mph_at_three_stages()
    test_run_finite_full_case_can_use_manual_mesh_construction()
    test_run_strip_cut_case_overwrites_one_mph_at_three_stages()
    print("All run_strip_1d shift scan tests passed.")
