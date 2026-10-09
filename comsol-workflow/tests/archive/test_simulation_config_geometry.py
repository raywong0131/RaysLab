from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from comsol_workflow.simulation_config_geometry import (
    SimulationConfigGeometryError,
    load_simulation_config_geometry,
)


REFERENCE_QUARTER_CONFIG = (
    Path(__file__).resolve().parents[2]
    / "results"
    / "finite_cavity_reference_results"
    / "finite_cavity_xy_fourier26_bulk5_cladding6_wide_bounds_mma_move0p01_eval0031"
    / "quarter_simulation_config.json"
)


def _write_config(path: Path, *, footprint=None, holes=None, length_unit="um") -> None:
    payload = {
        "length_unit": length_unit,
        "frequency_unit": "THz",
        "footprint": (
            [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]
            if footprint is None
            else footprint
        ),
        "materials": {"air": {"refractive_index": 1.0}},
        "layers": [
            {
                "thickness": 0.1,
                "material": "silicon",
                "holes": (
                    [[[0.2, 0.2], [0.5, 0.2], [0.2, 0.5]]]
                    if holes is None
                    else holes
                ),
                "hole_material": "air",
            },
            {"thickness": 1.55, "material": "air"},
        ],
        "physics": [{"type": "maxwell", "field": "E", "order": 2}],
        "boundary_conditions": [{"type": "pec", "boundaries": ["side0"]}],
        "shift_frequency": 194.9,
        "mode_count": 6,
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_loads_real_quarter_config_geometry_only() -> None:
    geometry = load_simulation_config_geometry(REFERENCE_QUARTER_CONFIG)

    assert geometry.hole_count == 607
    assert geometry.footprint.shape == (4, 2)
    assert geometry.footprint_bounds == pytest.approx((0.0, 9.84, 0.0, 8.69497636287113))
    assert geometry.footprint_width == pytest.approx(9.84)
    assert geometry.footprint_height == pytest.approx(8.69497636287113)
    assert geometry.is_axis_aligned_rectangle
    assert geometry.is_first_quadrant
    assert {hole.shape for hole in geometry.holes} == {(3, 2)}
    assert all(not hole.flags.writeable for hole in geometry.holes)
    assert not geometry.footprint.flags.writeable

    summary = geometry.summary()
    assert summary["length_unit"] == "um"
    assert summary["hole_count"] == 607
    assert summary["footprint_kind"] == "axis_aligned_rectangle"


def test_solver_fields_are_not_part_of_geometry_result(tmp_path: Path) -> None:
    source = tmp_path / "simulation_config.json"
    _write_config(source)

    geometry = load_simulation_config_geometry(source)

    assert not hasattr(geometry, "mode_count")
    assert not hasattr(geometry, "shift_frequency")
    assert geometry.hole_count == 1


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (lambda data: data.update(length_unit="nm"), "length_unit"),
        (lambda data: data.update(footprint=[[0, 0], [1, 0]]), "footprint"),
        (lambda data: data.update(layers=[]), "layer_index"),
        (
            lambda data: data["layers"][0].update(holes=[[[0.2, 0.2], [0.2, 0.5], [0.5, 0.2]]]),
            "counter-clockwise",
        ),
        (
            lambda data: data["layers"][0].update(
                holes=[[[0.2, 0.2], [2.5, 0.2], [0.2, 0.5]]]
            ),
            "outside footprint",
        ),
    ],
)
def test_rejects_invalid_geometry_contract(
    tmp_path: Path,
    mutator,
    message: str,
) -> None:
    source = tmp_path / "simulation_config.json"
    _write_config(source)
    data = json.loads(source.read_text(encoding="utf-8"))
    mutator(data)
    source.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(SimulationConfigGeometryError, match=message):
        load_simulation_config_geometry(source)


def test_rejects_non_numeric_and_non_convex_holes(tmp_path: Path) -> None:
    source = tmp_path / "simulation_config.json"
    _write_config(source)
    data = json.loads(source.read_text(encoding="utf-8"))
    data["layers"][0]["holes"] = [[[0.2, 0.2], [0.5, 0.2], ["bad", 0.5]]]
    source.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(SimulationConfigGeometryError, match="finite numbers"):
        load_simulation_config_geometry(source)

    _write_config(source)
    data = json.loads(source.read_text(encoding="utf-8"))
    data["layers"][0]["holes"] = [
        [[0.2, 0.2], [1.0, 0.2], [1.0, 1.0], [0.6, 0.4], [0.2, 1.0]]
    ]
    source.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(SimulationConfigGeometryError, match="strictly convex"):
        load_simulation_config_geometry(source)


def test_summary_hash_is_source_byte_hash(tmp_path: Path) -> None:
    source = tmp_path / "simulation_config.json"
    _write_config(source)
    geometry = load_simulation_config_geometry(source)

    import hashlib

    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    assert geometry.source_sha256 == expected
    assert np.allclose(geometry.holes[0], [[0.2, 0.2], [0.5, 0.2], [0.2, 0.5]])
