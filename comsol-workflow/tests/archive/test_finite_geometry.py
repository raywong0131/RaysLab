from __future__ import annotations

import numpy as np
import pytest

from comsol_workflow.cladding_shift_profile import (
    CladdingShiftGeometry,
    CladdingShiftProfile,
)
from comsol_workflow.finite_geometry import (
    FiniteGeometrySpec,
    build_finite_geometry_plan,
    square_cell_plot_records,
    square_fragment_cell_polygon,
)
from comsol_workflow.finite_geometry_comsol import compile_finite_geometry_input
from comsol_workflow.geometry_utils import create_hexagon_design
from scripts.run_main import run_finite as finite


def _templates(params: dict[str, float]) -> list[np.ndarray]:
    _, holes, _ = create_hexagon_design(finite.A, finite.get_hole_params(params))
    return list(holes)


def _square_plan(
    *,
    cavity_layers: int = 5,
    cladding_layers: int = 5,
    geometry: str = "groups",
    x_factor: float = 0.0,
    y_factor: float = 0.0,
    profile: CladdingShiftProfile | None = None,
):
    return build_finite_geometry_plan(
        FiniteGeometrySpec(
            "square",
            finite.A,
            cavity_layers,
            cladding_layers,
        ),
        _templates(finite.BULK_PARAMS),
        _templates(finite.CLADDING_PARAMS),
        shift_geometry=CladdingShiftGeometry(geometry),
        shift_profile=profile or CladdingShiftProfile("uniform"),
        x_shift_factor=x_factor,
        y_shift_factor=y_factor,
    )


def test_square_cav5_clad5_approved_counts() -> None:
    plan = _square_plan()
    assert plan.metadata["fragment_counts"] == {
        "total": 840,
        "cavity": 220,
        "cladding": 620,
    }
    assert plan.metadata["parent_cell_counts"] == {
        "total": 431,
        "full_cavity": 105,
        "mixed_cavity_cladding": 10,
        "outer_half_cell": 22,
    }
    assert len(plan.records) == 2520
    assert len({(p.m, p.n, p.group_id) for p in plan.placements}) == 840


def test_square_footprint_adds_half_period_buffer() -> None:
    plan = _square_plan()
    footprint = plan.footprint
    assert footprint.width / 2.0 - footprint.structured_half_width == pytest.approx(
        finite.A / 2.0
    )
    assert footprint.height / 2.0 - footprint.structured_half_height == pytest.approx(
        finite.A / 2.0
    )


def test_square_overview_uses_true_hex_cells_and_hex_half_cells() -> None:
    plan = _square_plan(geometry="ellipse", x_factor=0.05, y_factor=0.08)
    records = square_cell_plot_records(plan)

    assert len(records) == plan.metadata["parent_cell_counts"]["total"] + plan.metadata[
        "parent_cell_counts"
    ]["mixed_cavity_cladding"]
    assert {len(np.asarray(record["polygon"])) for record in records} == {4, 6}

    full = next(record for record in records if record["group_id"] == 0)
    expected_full = finite.cell_polygon(full["point"].x, full["point"].y)
    np.testing.assert_allclose(full["polygon"], expected_full)

    half = next(record for record in records if record["group_id"] in {1, 2})
    expected_half = square_fragment_cell_polygon(
        half["point"], half["group_id"], finite.A
    )
    np.testing.assert_allclose(half["polygon"], expected_half)
    assert np.ptp(np.asarray(half["polygon"])[:, 1]) == pytest.approx(
        2.0 * finite.A / np.sqrt(3.0)
    )


def test_square_overview_merges_rigid_cladding_parent_to_one_arrow_record() -> None:
    plan = _square_plan(geometry="ellipse", x_factor=0.05, y_factor=0.08)
    records = square_cell_plot_records(plan)
    placement_counts: dict[tuple[int, int], int] = {}
    record_counts: dict[tuple[int, int], int] = {}
    for placement in plan.placements:
        if placement.region == "cladding":
            key = (placement.m, placement.n)
            placement_counts[key] = placement_counts.get(key, 0) + 1
    for record in records:
        if record["region"] == "cladding":
            point = record["point"]
            key = (point.i, point.j)
            record_counts[key] = record_counts.get(key, 0) + 1

    rigid_parent = next(key for key, count in placement_counts.items() if count == 2)
    assert record_counts[rigid_parent] == 1


@pytest.mark.parametrize("geometry", ["groups", "ellipse"])
def test_square_shift_keeps_full_cladding_parents_rigid(geometry: str) -> None:
    plan = _square_plan(geometry=geometry, x_factor=0.05, y_factor=0.08)
    by_parent: dict[tuple[int, int], list] = {}
    for placement in plan.placements:
        by_parent.setdefault((placement.m, placement.n), []).append(placement)
    full_cladding = [
        placements
        for placements in by_parent.values()
        if len(placements) == 2
        and all(placement.region == "cladding" for placement in placements)
    ]
    assert full_cladding
    for placements in full_cladding:
        np.testing.assert_allclose(placements[0].shift, placements[1].shift)
        assert placements[0].layer == placements[1].layer


def test_square_groups_corner_uses_independent_axis_components() -> None:
    plan = _square_plan(x_factor=0.05, y_factor=0.08)
    corner = next(
        placement
        for placement in plan.placements
        if placement.region == "cladding"
        and placement.outside_x
        and placement.outside_y
        and placement.point.x > 0.0
        and placement.point.y > 0.0
        and placement.layer == 1
    )
    np.testing.assert_allclose(corner.shift, [-0.05 * finite.A, -0.08 * finite.A])
    assert plan.metadata["finite_geometry_shift_formula"] == "square_axis_components_v1"


def test_square_mixed_parent_moves_only_cladding_fragment() -> None:
    plan = _square_plan(x_factor=0.05, y_factor=0.08)
    by_parent: dict[tuple[int, int], list] = {}
    for placement in plan.placements:
        by_parent.setdefault((placement.m, placement.n), []).append(placement)
    mixed = next(
        placements
        for placements in by_parent.values()
        if len(placements) == 2
        and {placement.region for placement in placements} == {"bulk", "cladding"}
    )
    cavity = next(p for p in mixed if p.region == "bulk")
    cladding = next(p for p in mixed if p.region == "cladding")
    np.testing.assert_allclose(cavity.shift, [0.0, 0.0])
    assert np.linalg.norm(cladding.shift) > 0.0


def test_hex_dispatch_preserves_existing_cell_and_hole_counts() -> None:
    cavity_layers = 2
    plan = build_finite_geometry_plan(
        FiniteGeometrySpec("hex", finite.A, cavity_layers, 2),
        _templates(finite.BULK_PARAMS),
        _templates(finite.CLADDING_PARAMS),
        shift_geometry=CladdingShiftGeometry("groups"),
        shift_profile=CladdingShiftProfile("uniform"),
        x_shift_factor=0.05,
        y_shift_factor=0.08,
    )
    assert len(plan.cavity_analysis_cells) == 7
    assert len(plan.cladding_parent_cells) == 30
    assert len(plan.records) == 37 * 6
    assert {point.shell for point in plan.cavity_analysis_cells} == {0, 1}
    assert {point.shell for point in plan.cladding_parent_cells} == {2, 3}
    assert plan.metadata["cavity_layers"] == cavity_layers
    edge_cells = [
        point
        for point in plan.cavity_analysis_cells
        if point.i >= 0
        and point.j >= 0
        and point.i + point.j == cavity_layers - 1
    ]
    assert len(edge_cells) == cavity_layers
    assert plan.footprint.boundary_shape == "hexagon_boundary"


def test_square_comsol_compiler_uses_rectangle_boundaries() -> None:
    plan = _square_plan(cavity_layers=1, cladding_layers=1)
    full = compile_finite_geometry_input(
        plan,
        layer_height="200 [nm]",
        refractive_index="3.3",
    )
    quarter = compile_finite_geometry_input(
        plan,
        layer_height="200 [nm]",
        refractive_index="3.3",
        quarter=True,
    )
    assert full.boundary.shape == "rectangle"
    assert quarter.boundary.shape == "quarter_rectangle"
    assert full.boundary.a == pytest.approx(plan.footprint.width)
    assert full.boundary.b == pytest.approx(plan.footprint.height)


def test_square_nonuniform_profile_uses_absolute_parent_layer_envelope() -> None:
    profile = CladdingShiftProfile("tanh_power", scale_layers=3.0, power=2.0)
    plan = _square_plan(
        x_factor=0.05,
        y_factor=0.08,
        profile=profile,
    )
    placement = next(
        item
        for item in plan.placements
        if item.region == "cladding"
        and item.outside_x
        and not item.outside_y
        and item.point.x > 0.0
        and item.layer == 2
    )
    assert placement.shift[0] == pytest.approx(
        -finite.A * 0.05 * profile.envelope(2)
    )
    assert placement.shift[1] == pytest.approx(0.0)
