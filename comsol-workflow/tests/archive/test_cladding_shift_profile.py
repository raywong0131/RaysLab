import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from comsol_workflow.cladding_shift_profile import (
    ELLIPSE_SHIFT_GEOMETRY,
    GROUPS_SHIFT_GEOMETRY,
    CladdingShiftGeometry,
    CladdingShiftProfile,
    cladding_shift_group,
    directional_shift_profile_metadata,
    ellipse_inward_shift,
    ellipse_shift_profile_metadata,
    ellipse_shift_weight,
    finite_shift_profile_metadata,
    relative_cladding_layer,
    resolved_shift_profile,
    save_directional_shift_profile_diagnostics,
    save_ellipse_shift_profile_diagnostics,
    save_shift_profile_diagnostics,
    shift_profile_metadata,
    shift_profile_preflight_messages,
)
from comsol_workflow.hex_lattice_utils import lattice_points
from scripts.run_main import parameter_config
from scripts.run_main import run_finite
from scripts.run_main import run_finite_quarter
from scripts.run_main import run_strip_1d


def _shared_parameter_payload() -> dict[str, object]:
    return {
        "schema_version": 1,
        "cavity": {"b0_nm": 245.0, "eta": 0.96, "zeta": 1.155},
        "cladding": {"b0_nm": 240.2, "eta": 0.967, "zeta": 1.0},
        "cavity_layers": 10,
        "cladding_layers": 10,
        "mesh_size": 9,
        "unit_cell_eigenmode_count": 6,
        "finite_eigenmode_count": 1,
        "center_frequency_thz": 198.4,
        "cladding_shift_factors": [0.04],
    }


def test_active_profile_matches_parameter_file():
    profile = parameter_config.ACTIVE_PARAMETERS.cladding_shift_profile
    payload = parameter_config.parameter_data_to_flat(
        json.loads(parameter_config.PARAMETER_PATH.read_text(encoding="utf-8"))
    )
    expected = CladdingShiftProfile.from_mapping(
        payload.get("cladding_shift_profile")
    )

    assert profile == expected
    assert profile.label == expected.label
    series_label = parameter_config.ACTIVE_PARAMETERS.structure_series_label(
        "finite"
    )
    if parameter_config.ACTIVE_PARAMETERS.cladding_shift_geometry.kind == "ellipse":
        assert series_label.endswith("_ellipse-shift")
    else:
        assert not series_label.endswith("_ellipse-shift")
    expected_ratio = float(payload["cladding_y_over_x_shift_ratio"])
    assert (
        parameter_config.ACTIVE_PARAMETERS.cladding_y_over_x_shift_ratio
        == expected_ratio
    )
    assert parameter_config.ACTIVE_PARAMETERS.cladding_shift_input_kind == (
        "scan_ratio"
    )
    expected_pairs = tuple(
        (factor, factor * expected_ratio)
        for factor in parameter_config.ACTIVE_PARAMETERS.cladding_shift_factors
    )
    np.testing.assert_allclose(
        parameter_config.ACTIVE_PARAMETERS.finite_cladding_shift_factor_pairs,
        expected_pairs,
    )
    expected_last_factor = float(payload["cladding_shift_factors"][-1])
    assert expected_pairs[-1] == pytest.approx(
        (expected_last_factor, expected_last_factor * expected_ratio)
    )
    assert parameter_config.ACTIVE_PARAMETERS.cladding_shift_geometry == (
        CladdingShiftGeometry.from_mapping(
            payload.get("cladding_shift_geometry")
        )
    )


def test_missing_profile_reproduces_legacy_uniform_series_name(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    path.write_text(json.dumps(payload), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.cladding_shift_profile == CladdingShiftProfile()
    assert shared.cladding_shift_profile.to_dict() == {"kind": "uniform"}
    assert shared.cladding_y_over_x_shift_ratio == 1.0
    assert shared.cladding_shift_input_kind == "scan_ratio"
    assert shared.finite_cladding_shift_factor_pairs == ((0.04, 0.04),)
    assert shared.cladding_shift_geometry == CladdingShiftGeometry()
    assert not shared.structure_series_label("finite").endswith("shiftprof-uniform")
    assert not shared.structure_series_label("finite").endswith("_ellipse-shift")


def test_explicit_directional_factors_select_one_xy_case_without_affecting_identity(
    tmp_path,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["cladding_x_shift_factor"] = 0.05
    payload["cladding_y_shift_factor"] = 0.12
    path.write_text(json.dumps(payload), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.cladding_shift_input_kind == "xy"
    assert shared.cladding_y_over_x_shift_ratio is None
    assert shared.finite_cladding_shift_factor_pairs == ((0.05, 0.12),)
    assert shared.finite_cladding_shift_input == {
        "kind": "xy",
        "x_factor": 0.05,
        "y_factor": 0.12,
    }
    assert "finite_cladding_shift_input" not in shared.identity


def test_ellipse_series_suffix_applies_only_to_finite_workflows(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["cladding_shift_geometry"] = {"kind": "ellipse"}
    path.write_text(json.dumps(payload), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.structure_series_label("finite").endswith("_ellipse-shift")
    assert shared.structure_series_label("finite_quarter").endswith(
        "_ellipse-shift"
    )
    assert not shared.structure_series_label("strip1d").endswith(
        "_ellipse-shift"
    )


def test_ellipse_parameter_configuration_accepts_zero_y_over_x_ratio(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["cladding_y_over_x_shift_ratio"] = 0.0
    payload["cladding_shift_geometry"] = {"kind": "ellipse"}
    path.write_text(json.dumps(payload), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.finite_cladding_shift_factor_pairs == ((0.04, 0.0),)
    assert shared.cladding_shift_geometry.kind == ELLIPSE_SHIFT_GEOMETRY


def test_explicit_xy_configuration_accepts_x_zero_y_positive(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["cladding_x_shift_factor"] = 0.0
    payload["cladding_y_shift_factor"] = 0.08
    payload["cladding_shift_geometry"] = {"kind": "ellipse"}
    path.write_text(json.dumps(payload), encoding="utf-8")

    shared = parameter_config.load_shared_parameters(path)

    assert shared.cladding_shift_input_kind == "xy"
    assert shared.finite_cladding_shift_factor_pairs == ((0.0, 0.08),)
    assert shared.cladding_y_over_x_shift_ratio is None


@pytest.mark.parametrize(
    ("mapping", "error_type", "message"),
    [
        ([], TypeError, "must be an object"),
        ({}, TypeError, "kind must be a string"),
        ({"kind": "radial"}, ValueError, "must be one of"),
        (
            {"kind": "ellipse", "axis": "x"},
            ValueError,
            "unknown fields",
        ),
    ],
)
def test_shift_geometry_configuration_rejects_invalid_records(
    mapping,
    error_type,
    message,
):
    with pytest.raises(error_type, match=message):
        CladdingShiftGeometry.from_mapping(mapping)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("cladding_y_over_x_shift_ratio", float("nan")),
        ("cladding_y_over_x_shift_ratio", float("inf")),
        ("cladding_y_over_x_shift_ratio", -0.01),
        ("cladding_shift_factors", [-0.01]),
    ],
)
def test_shift_scan_ratio_configuration_rejects_invalid_values(
    tmp_path,
    field,
    value,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload[field] = value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=field):
        parameter_config.load_shared_parameters(path)


@pytest.mark.parametrize(
    "payload_updates",
    [
        {
            "cladding_y_over_x_shift_ratio": 1.6,
            "cladding_x_shift_factor": 0.05,
            "cladding_y_shift_factor": 0.08,
        },
        {"cladding_x_shift_factor": 0.05},
        {"cladding_y_shift_factor": 0.08},
    ],
)
def test_shift_input_modes_reject_conflicts_and_incomplete_pairs(
    tmp_path,
    payload_updates,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload.update(payload_updates)
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="conflicts|provided together"):
        parameter_config.load_shared_parameters(path)


@pytest.mark.parametrize(
    ("x_factor", "y_factor", "field"),
    [
        (-0.01, 0.08, "cladding_x_shift_factor"),
        (0.05, float("inf"), "cladding_y_shift_factor"),
    ],
)
def test_explicit_xy_configuration_rejects_invalid_factors(
    tmp_path,
    x_factor,
    y_factor,
    field,
):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    payload["cladding_x_shift_factor"] = x_factor
    payload["cladding_y_shift_factor"] = y_factor
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=field):
        parameter_config.load_shared_parameters(path)


def test_unit_cell_frequency_update_ignores_profile_change(tmp_path):
    path = tmp_path / "parameter.json"
    payload = _shared_parameter_payload()
    path.write_text(json.dumps(payload), encoding="utf-8")
    expected = parameter_config.load_shared_parameters(path)
    payload["cladding_shift_profile"] = {
        "kind": "tanh_power",
        "scale_layers": 4.0,
        "power": 2.0,
    }
    path.write_text(json.dumps(payload), encoding="utf-8")

    updated = parameter_config.update_center_frequency_from_cavity_p2(
        199.0,
        expected_parameters=expected,
        path=path,
    )

    assert updated.center_frequency_thz == 199.0
    assert updated.cladding_shift_profile.label == "tanhpow-l4-p2"


@pytest.mark.parametrize(
    ("mapping", "error_type", "message"),
    [
        ([], TypeError, "must be an object"),
        ({}, TypeError, "kind must be a string"),
        ({"kind": "unknown"}, ValueError, "must be one of"),
        (
            {"kind": "uniform", "power": 2.0},
            ValueError,
            "unknown fields",
        ),
        (
            {"kind": "tanh_power", "scale_layers": 4.0},
            ValueError,
            "is missing: power",
        ),
        (
            {"kind": "tanh_power", "scale_layers": 0.0, "power": 2.0},
            ValueError,
            "scale_layers must be positive",
        ),
        (
            {"kind": "tanh_power", "scale_layers": 4.0, "power": 1.0},
            ValueError,
            "power must be greater than 1",
        ),
    ],
)
def test_profile_configuration_rejects_invalid_records(mapping, error_type, message):
    with pytest.raises(error_type, match=message):
        CladdingShiftProfile.from_mapping(mapping)


def test_tanh_power_profile_has_expected_layer_values_and_saturation():
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    values = [profile.envelope(layer) for layer in range(11)]

    assert values[0] == 0.0
    assert values[1] == pytest.approx(math.tanh(1.0 / 16.0))
    assert values[4] == pytest.approx(math.tanh(1.0))
    assert values[10] == pytest.approx(0.9999925467214317)
    assert all(left < right for left, right in zip(values, values[1:]))


def test_profile_preflight_reports_default_saturation_without_warning():
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )

    messages = shift_profile_preflight_messages(
        shift_profile_metadata(10, 0.05, 0.82, profile)
    )

    assert len(messages) == 1
    assert "outer envelope=0.999993" in messages[0]
    assert "layer 8 envelope=0.999329" in messages[0]


def test_profile_preflight_warns_without_normalizing_short_cladding():
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    metadata = shift_profile_metadata(4, 0.05, 0.82, profile)

    messages = shift_profile_preflight_messages(metadata)

    assert len(messages) == 2
    assert messages[1].startswith("WARNING:")
    assert metadata["outer_envelope"] == pytest.approx(math.tanh(1.0))


def test_relative_layer_uses_hex_shell_offset():
    assert relative_cladding_layer(SimpleNamespace(shell=10), 10) == 0
    assert relative_cladding_layer(SimpleNamespace(shell=11), 10) == 1
    assert relative_cladding_layer(SimpleNamespace(shell=20), 10) == 10

    with pytest.raises(ValueError, match="lies inside"):
        relative_cladding_layer(SimpleNamespace(shell=9), 10)


def test_resolved_profile_uses_absolute_not_cumulative_shift():
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    rows = resolved_shift_profile(10, 0.09, 0.82, profile)

    assert rows[1]["side_shift_um"] == pytest.approx(
        0.09 * 0.82 * profile.envelope(1)
    )
    assert rows[-1]["side_shift_um"] == pytest.approx(
        0.09 * 0.82 * profile.envelope(10)
    )
    assert rows[-1]["corner_shift_um"] / rows[-1][
        "side_shift_um"
    ] == pytest.approx(2.0 / math.sqrt(3.0))
    assert rows[-1]["side_shift_um"] < 0.09 * 0.82


def test_zero_target_produces_zero_shift_for_every_layer():
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    rows = resolved_shift_profile(10, 0.0, 0.82, profile)

    assert all(row["effective_factor"] == 0.0 for row in rows)
    assert all(row["side_shift_um"] == 0.0 for row in rows)
    assert all(row["corner_shift_um"] == 0.0 for row in rows)


def test_directional_region_partition_covers_exactly_six_x_and_six_y_regions():
    x_regions = {
        (kind, idx)
        for kind in ("side", "corner")
        for idx in range(6)
        if cladding_shift_group(kind, idx) == "x"
    }
    y_regions = {
        (kind, idx)
        for kind in ("side", "corner")
        for idx in range(6)
        if cladding_shift_group(kind, idx) == "y"
    }

    assert x_regions == {
        ("side", 0),
        ("side", 2),
        ("side", 3),
        ("side", 5),
        ("corner", 0),
        ("corner", 3),
    }
    assert y_regions == {
        ("side", 1),
        ("side", 4),
        ("corner", 1),
        ("corner", 2),
        ("corner", 4),
        ("corner", 5),
    }
    assert x_regions.isdisjoint(y_regions)


def test_finite_and_strip_apply_the_same_layer_resolved_vectors(monkeypatch):
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    target = 0.09
    monkeypatch.setattr(
        run_finite,
        "CLADDING_SHIFT_GEOMETRY",
        CladdingShiftGeometry(kind=GROUPS_SHIFT_GEOMETRY),
    )
    monkeypatch.setattr(run_finite, "CLADDING_SHIFT_PROFILE", profile)
    monkeypatch.setattr(run_finite, "CLADDING_X_SHIFT_FACTOR", target)
    monkeypatch.setattr(run_finite, "CLADDING_Y_SHIFT_FACTOR", target)
    monkeypatch.setattr(run_finite, "CLADDING_X_INWARD_SHIFT", target * run_finite.A)
    monkeypatch.setattr(run_finite, "CLADDING_Y_INWARD_SHIFT", target * run_finite.A)
    monkeypatch.setattr(run_strip_1d, "CLADDING_SHIFT_PROFILE", profile)
    monkeypatch.setattr(run_strip_1d, "CLADDING_INWARD_SHIFT_FACTOR", target)
    monkeypatch.setattr(run_strip_1d, "CLADDING_INWARD_SHIFT", target * run_strip_1d.A)

    max_shell = run_finite.BULK_RADIUS + min(run_finite.CLADDING_LAYERS, 3)
    points = lattice_points(
        max_shell,
        run_finite.A,
        min_shell=run_finite.BULK_RADIUS + 1,
    )
    for point in points:
        finite_shift = run_finite.cladding_inward_shift(point)
        strip_shift = run_strip_1d.cladding_inward_shift(point)
        layer = point.shell - run_finite.BULK_RADIUS
        kind, _ = run_finite.normal_fan_feature(point)
        expected_magnitude = target * run_finite.A * profile.envelope(layer)
        if kind == "corner":
            expected_magnitude *= 2.0 / math.sqrt(3.0)

        assert strip_shift == pytest.approx(finite_shift)
        assert np.linalg.norm(finite_shift) == pytest.approx(expected_magnitude)


def test_uniform_profile_keeps_every_actual_cladding_layer_at_target(monkeypatch):
    target = 0.04
    monkeypatch.setattr(
        run_finite,
        "CLADDING_SHIFT_GEOMETRY",
        CladdingShiftGeometry(kind=GROUPS_SHIFT_GEOMETRY),
    )
    monkeypatch.setattr(run_finite, "CLADDING_SHIFT_PROFILE", CladdingShiftProfile())
    monkeypatch.setattr(run_finite, "CLADDING_X_SHIFT_FACTOR", target)
    monkeypatch.setattr(run_finite, "CLADDING_Y_SHIFT_FACTOR", target)
    monkeypatch.setattr(run_finite, "CLADDING_X_INWARD_SHIFT", target * run_finite.A)
    monkeypatch.setattr(run_finite, "CLADDING_Y_INWARD_SHIFT", target * run_finite.A)

    max_shell = run_finite.BULK_RADIUS + min(run_finite.CLADDING_LAYERS, 3)
    points = lattice_points(
        max_shell,
        run_finite.A,
        min_shell=run_finite.BULK_RADIUS + 1,
    )
    for point in points:
        shift = run_finite.cladding_inward_shift(point)
        kind, _ = run_finite.normal_fan_feature(point)
        expected = target * run_finite.A
        if kind == "corner":
            expected *= 2.0 / math.sqrt(3.0)
        assert np.linalg.norm(shift) == pytest.approx(expected)


def test_finite_uses_distinct_x_and_y_group_target_factors(monkeypatch):
    x_target = 0.05
    y_target = 0.12
    monkeypatch.setattr(
        run_finite,
        "CLADDING_SHIFT_GEOMETRY",
        CladdingShiftGeometry(kind=GROUPS_SHIFT_GEOMETRY),
    )
    monkeypatch.setattr(run_finite, "CLADDING_SHIFT_PROFILE", CladdingShiftProfile())
    monkeypatch.setattr(run_finite, "CLADDING_X_SHIFT_FACTOR", x_target)
    monkeypatch.setattr(run_finite, "CLADDING_Y_SHIFT_FACTOR", y_target)
    monkeypatch.setattr(run_finite, "CLADDING_X_INWARD_SHIFT", x_target * run_finite.A)
    monkeypatch.setattr(run_finite, "CLADDING_Y_INWARD_SHIFT", y_target * run_finite.A)

    points = lattice_points(
        run_finite.BULK_RADIUS + run_finite.CLADDING_LAYERS,
        run_finite.A,
        min_shell=run_finite.BULK_RADIUS + 1,
    )
    representatives = {}
    for point in points:
        representatives.setdefault(run_finite.normal_fan_feature(point), point)
    assert len(representatives) == 12

    for (kind, idx), point in representatives.items():
        group = cladding_shift_group(kind, idx)
        target = x_target if group == "x" else y_target
        expected = target * run_finite.A
        if kind == "corner":
            expected *= 2.0 / math.sqrt(3.0)
        assert np.linalg.norm(run_finite.cladding_inward_shift(point)) == pytest.approx(
            expected
        )


def test_ellipse_weight_anchors_axes_and_interpolates_smoothly():
    x_factor = 0.05
    y_factor = 0.08

    assert ellipse_shift_weight(1.0, 0.0, x_factor, y_factor) == pytest.approx(
        x_factor
    )
    assert ellipse_shift_weight(0.0, 1.0, x_factor, y_factor) == pytest.approx(
        y_factor
    )
    assert ellipse_shift_weight(1.0, 1.0, x_factor, y_factor) == pytest.approx(
        0.065
    )


def test_ellipse_shift_uses_radial_direction_and_supports_one_zero_axis():
    period = 0.82

    assert ellipse_inward_shift(2.0, 0.0, period, 0.05, 0.0) == pytest.approx(
        (-period * 0.05, 0.0)
    )
    assert ellipse_inward_shift(0.0, 2.0, period, 0.05, 0.0) == pytest.approx(
        (0.0, 0.0)
    )
    diagonal = np.asarray(
        ellipse_inward_shift(1.0, 1.0, period, 0.05, 0.0)
    )
    assert np.linalg.norm(diagonal) == pytest.approx(period * 0.025)
    assert diagonal[0] == pytest.approx(diagonal[1])
    assert diagonal[0] < 0.0


def test_finite_ellipse_shift_matches_principal_axes_and_c2v(monkeypatch):
    monkeypatch.setattr(
        run_finite,
        "CLADDING_SHIFT_GEOMETRY",
        CladdingShiftGeometry(kind=ELLIPSE_SHIFT_GEOMETRY),
    )
    monkeypatch.setattr(run_finite, "CLADDING_SHIFT_PROFILE", CladdingShiftProfile())
    monkeypatch.setattr(run_finite, "CLADDING_X_SHIFT_FACTOR", 0.05)
    monkeypatch.setattr(run_finite, "CLADDING_Y_SHIFT_FACTOR", 0.08)

    def point(x, y):
        return SimpleNamespace(
            x=float(x),
            y=float(y),
            shell=run_finite.BULK_RADIUS + 1,
        )

    assert run_finite.cladding_inward_shift(point(2.0, 0.0)) == pytest.approx(
        (-run_finite.A * 0.05, 0.0)
    )
    assert run_finite.cladding_inward_shift(point(0.0, 2.0)) == pytest.approx(
        (0.0, -run_finite.A * 0.08)
    )

    reference = run_finite.cladding_inward_shift(point(2.0, 1.0))
    mirror_x = run_finite.cladding_inward_shift(point(-2.0, 1.0))
    mirror_y = run_finite.cladding_inward_shift(point(2.0, -1.0))
    assert mirror_x == pytest.approx((-reference[0], reference[1]))
    assert mirror_y == pytest.approx((reference[0], -reference[1]))
    assert 2.0 * reference[1] - reference[0] == pytest.approx(0.0)
    assert np.dot(np.array([2.0, 1.0]), reference) < 0.0


def test_ellipse_metadata_uses_axes_without_corner_compensation():
    geometry = CladdingShiftGeometry(kind=ELLIPSE_SHIFT_GEOMETRY)
    metadata = finite_shift_profile_metadata(
        geometry,
        10,
        0.05,
        0.08,
        0.82,
        CladdingShiftProfile(),
    )

    assert metadata["cladding_shift_geometry"] == {"kind": "ellipse"}
    assert metadata["cladding_shift_formula"]["kind"] == (
        "radial_cos2_sin2_v1"
    )
    assert set(metadata["resolved_cladding_shift_by_axis"]) == {"x", "y"}
    assert "resolved_cladding_shift_by_group" not in metadata
    assert "cladding_shift_region_groups" not in metadata
    assert metadata["resolved_cladding_shift_by_axis"]["x"][-1][
        "axis_shift_um"
    ] == pytest.approx(0.82 * 0.05)
    assert metadata["resolved_cladding_shift_by_axis"]["y"][-1][
        "axis_shift_um"
    ] == pytest.approx(0.82 * 0.08)


def test_full_geometry_preserves_hole_rigidity_and_quarter_mirror_symmetry(monkeypatch):
    x_target = 0.05
    y_target = 0.12
    monkeypatch.setattr(
        run_finite,
        "CLADDING_SHIFT_GEOMETRY",
        CladdingShiftGeometry(kind=ELLIPSE_SHIFT_GEOMETRY),
    )
    monkeypatch.setattr(run_finite, "CLADDING_X_SHIFT_FACTOR", x_target)
    monkeypatch.setattr(run_finite, "CLADDING_Y_SHIFT_FACTOR", y_target)
    monkeypatch.setattr(run_finite, "CLADDING_X_INWARD_SHIFT", x_target * run_finite.A)
    monkeypatch.setattr(run_finite, "CLADDING_Y_INWARD_SHIFT", y_target * run_finite.A)

    records, metadata = run_finite.build_full_lattice_holes()
    first_point = metadata["cladding_points"][0]
    point_records = [
        record
        for record in records
        if record["region"] == "cladding" and record["point"] == first_point
    ]

    assert len(point_records) > 1
    assert all(
        record["modulation_shift"] == pytest.approx(
            point_records[0]["modulation_shift"]
        )
        for record in point_records
    )
    assert len(metadata["resolved_cladding_shift_by_cell"]) == len(
        metadata["cladding_points"]
    )
    run_finite_quarter.validate_full_geometry_mirror_symmetry(records)


def test_profile_diagnostics_write_resolved_csv_and_png(tmp_path):
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    metadata = shift_profile_metadata(10, 0.05, 0.82, profile)

    csv_path, png_path = save_shift_profile_diagnostics(tmp_path, metadata, dpi=80)

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [int(row["layer"]) for row in rows] == list(range(11))
    assert float(rows[-1]["effective_factor"]) == pytest.approx(
        0.05 * profile.envelope(10)
    )
    assert png_path.is_file()
    assert png_path.stat().st_size > 0


def test_directional_profile_diagnostics_write_both_groups(tmp_path):
    metadata = directional_shift_profile_metadata(
        10,
        0.05,
        0.12,
        0.82,
        CladdingShiftProfile(),
    )

    csv_path, png_path = save_directional_shift_profile_diagnostics(
        tmp_path,
        metadata,
        dpi=80,
    )

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["group"] for row in rows} == {"x", "y"}
    assert len(rows) == 22
    x_outer = next(
        row for row in rows if row["group"] == "x" and row["layer"] == "10"
    )
    y_outer = next(
        row for row in rows if row["group"] == "y" and row["layer"] == "10"
    )
    assert float(x_outer["target_factor"]) == 0.05
    assert float(y_outer["target_factor"]) == 0.12
    assert png_path.stat().st_size > 0


def test_ellipse_profile_diagnostics_write_both_principal_axes(tmp_path):
    metadata = ellipse_shift_profile_metadata(
        10,
        0.05,
        0.08,
        0.82,
        CladdingShiftProfile(),
    )

    csv_path, png_path = save_ellipse_shift_profile_diagnostics(
        tmp_path,
        metadata,
        dpi=80,
    )

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["axis"] for row in rows} == {"x", "y"}
    assert len(rows) == 22
    assert "corner_shift_um" not in rows[0]
    assert png_path.stat().st_size > 0


def test_strip_completion_rejects_profile_mismatch(tmp_path, monkeypatch):
    target = 0.04
    profile = CladdingShiftProfile(
        kind="tanh_power",
        scale_layers=4.0,
        power=2.0,
    )
    monkeypatch.setattr(run_strip_1d, "CLADDING_SHIFT_PROFILE", profile)
    monkeypatch.setattr(run_strip_1d, "CLADDING_INWARD_SHIFT_FACTOR", target)
    monkeypatch.setattr(
        run_strip_1d,
        "CLADDING_INWARD_SHIFT",
        target * run_strip_1d.A,
    )
    required = [
        tmp_path / "strip_1d.mph",
        tmp_path / run_strip_1d.EXPORT_DIRNAME / "eigenfrequencies.csv",
        tmp_path
        / "strip_bulk_fourier_hz"
        / "strip_bulk_fourier_summary.csv",
        tmp_path / "mode_scores.csv",
        tmp_path / "objective.json",
    ]
    for path in required:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    metadata = run_strip_1d.shared_case_metadata(
        {
            "bulk_radius": run_strip_1d.BULK_RADIUS,
            "cladding_layers": run_strip_1d.CLADDING_LAYERS,
            "cladding_inward_shift": target * run_strip_1d.A,
        },
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
    )
    metadata["strip_cut_angle_degrees"] = run_strip_1d.STRIP_CUT_ANGLE
    (tmp_path / "config.json").write_text(
        json.dumps(run_strip_1d.to_jsonable(metadata)),
        encoding="utf-8",
    )

    assert run_strip_1d.strip_case_is_complete(tmp_path, target)

    config = json.loads((tmp_path / "config.json").read_text(encoding="utf-8"))
    config["cladding_shift_profile"] = {"kind": "uniform"}
    (tmp_path / "config.json").write_text(json.dumps(config), encoding="utf-8")
    assert not run_strip_1d.strip_case_is_complete(tmp_path, target)


def test_resume_metadata_contract_tracks_layer_profile_identity():
    saved = finite_shift_profile_metadata(
        CladdingShiftGeometry(kind=GROUPS_SHIFT_GEOMETRY),
        10,
        0.05,
        0.12,
        0.82,
        CladdingShiftProfile(),
    )
    current = finite_shift_profile_metadata(
        CladdingShiftGeometry(kind=GROUPS_SHIFT_GEOMETRY),
        10,
        0.05,
        0.12,
        0.82,
        CladdingShiftProfile(
            kind="tanh_power",
            scale_layers=4.0,
            power=2.0,
        ),
    )
    paths = (
        *run_finite.RESUME_METADATA_PATHS,
        *run_finite.GROUPS_RESUME_METADATA_PATHS,
    )
    monkey_paths = (
        "cladding_x_shift_factor",
        "cladding_y_shift_factor",
        "cladding_shift_profile",
        "cladding_shift_profile_label",
        "cladding_shift_geometry",
        "cladding_shift_region_groups",
        "resolved_cladding_shift_by_group",
    )
    assert all(path in paths for path in monkey_paths)
    assert saved["resolved_cladding_shift_by_group"] != current[
        "resolved_cladding_shift_by_group"
    ]
    assert set(run_finite.ELLIPSE_RESUME_METADATA_PATHS) == {
        "cladding_shift_formula",
        "resolved_cladding_shift_by_axis",
        "resolved_cladding_shift_by_cell",
    }
