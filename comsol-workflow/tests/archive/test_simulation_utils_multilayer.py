import os
import sys

import numpy as np

from comsol_workflow.simulation_utils import (
    BoundarySpec,
    LayerSpec,
    SimulationRun,
    _boundary_definition,
    _make_layer_stack,
    _validate_build_inputs,
    _wall_seed_points,
    electric_reflection_matrix,
    hz_reflection_sign,
    reconstruct_quarter_electric_grid,
    reconstruct_quarter_scalar_field,
)


def triangle():
    return np.array([
        [0.0, 0.0],
        [0.1, 0.0],
        [0.0, 0.1],
    ])


def test_boundary_rectangle_validation():
    boundary = BoundarySpec("rectangle", a=2.0, b=3.0)
    assert boundary.shape == "rectangle"
    assert boundary.a == 2.0
    assert boundary.b == 3.0

    try:
        BoundarySpec("rectangle", a=2.0)
        assert False, "Expected rectangle without b to fail"
    except ValueError as exc:
        assert "b is required" in str(exc)

    print("[Pass] rectangle boundary validation")


def test_hexagon_and_square_reject_b():
    for shape in ["hexagon", "hexagon_boundary", "quarter_hexagon_boundary", "square"]:
        try:
            BoundarySpec(shape, a=1.0, b=2.0)
            assert False, f"Expected {shape} with b to fail"
        except ValueError as exc:
            assert "b must be None" in str(exc)

    print("[Pass] hexagon/hexagon_boundary/square reject b")


def test_hexagon_boundary_footprint_uses_left_right_vertices():
    definition = _boundary_definition(BoundarySpec("hexagon_boundary", a=2.0))

    assert definition.footprint_label == "HexagonBoundary"
    assert definition.footprint_table == [
        ["a/2", "0"],
        ["a/4", "sqrt(3)*a/4"],
        ["-a/4", "sqrt(3)*a/4"],
        ["-a/2", "0"],
        ["-a/4", "-sqrt(3)*a/4"],
        ["a/4", "-sqrt(3)*a/4"],
    ]
    assert definition.gb_expr is None

    print("[Pass] hexagon_boundary footprint uses left/right vertices")


def test_quarter_hexagon_boundary_footprint_and_wall_seeds():
    boundary = BoundarySpec("quarter_hexagon_boundary", a=4.0)
    definition = _boundary_definition(boundary)
    root3 = float(np.sqrt(3.0))

    assert definition.footprint_table == [
        ["0", "0"],
        ["a/2", "0"],
        ["a/4", "sqrt(3)*a/4"],
        ["0", "sqrt(3)*a/4"],
    ]
    assert _wall_seed_points(boundary, z=1.5) == {
        "symmetry_x_axis": [(1.0, 0.0, 1.5)],
        "symmetry_y_axis": [(0.0, root3 / 2.0, 1.5)],
        "outer_hex_ne": [(1.5, root3 / 2.0, 1.5)],
        "outer_hex_n": [(0.5, root3, 1.5)],
    }


def test_quarter_rectangle_footprint_and_wall_seeds():
    boundary = BoundarySpec("quarter_rectangle", a=4.0, b=6.0)
    definition = _boundary_definition(boundary)
    assert definition.footprint_table == [
        ["0", "0"],
        ["a/2", "0"],
        ["a/2", "b/2"],
        ["0", "b/2"],
    ]
    assert _wall_seed_points(boundary, z=1.5) == {
        "symmetry_x_axis": [(1.0, 0.0, 1.5)],
        "symmetry_y_axis": [(0.0, 1.5, 1.5)],
        "outer_rect_x": [(2.0, 1.5, 1.5)],
        "outer_rect_y": [(1.0, 3.0, 1.5)],
    }


def test_rectangle_wall_seeds_and_k_scale():
    boundary = BoundarySpec("rectangle", a=2.0, b=4.0)
    seeds = _wall_seed_points(boundary, z=1.5)
    assert seeds["side_x"] == [(1.0, 0.0, 1.5), (-1.0, 0.0, 1.5)]
    assert seeds["side_y"] == [(0.0, 2.0, 1.5), (0.0, -2.0, 1.5)]

    definition = _boundary_definition(boundary)
    assert definition.ga_expr == "2*pi/a"
    assert definition.gb_expr == "2*pi/b"

    print("[Pass] rectangle wall seeds and k scale")


def test_hexagon_boundary_wall_seeds_group_opposite_edges():
    boundary = BoundarySpec("hexagon_boundary", a=4.0)
    seeds = _wall_seed_points(boundary, z=1.5)
    root3 = float(np.sqrt(3.0))

    assert set(seeds) == {"side_hex_ne_sw", "side_hex_n_s", "side_hex_nw_se"}
    assert seeds["side_hex_ne_sw"] == [(1.5, root3 / 2.0, 1.5), (-1.5, -root3 / 2.0, 1.5)]
    assert seeds["side_hex_n_s"] == [(0.0, root3, 1.5), (0.0, -root3, 1.5)]
    assert seeds["side_hex_nw_se"] == [(-1.5, root3 / 2.0, 1.5), (1.5, -root3 / 2.0, 1.5)]

    print("[Pass] hexagon_boundary wall seeds group opposite edges")


def test_single_parameter_boundaries_have_only_ga():
    for shape, expected_ga in [
        ("hexagon", "4*pi/sqrt(3)/a"),
        ("hexagon_boundary", "2*pi/a"),
        ("quarter_hexagon_boundary", "2*pi/a"),
        ("square", "2*pi/a"),
    ]:
        definition = _boundary_definition(BoundarySpec(shape, a=2.0))
        assert definition.ga_expr == expected_ga
        assert definition.gb_expr is None

    print("[Pass] hexagon/square use Ga only")


def test_layer_stack_first_layer_half_height():
    layers = [
        LayerSpec("200 [nm]", [triangle()], "3.3"),
        LayerSpec("50 [nm]", [], "2.0"),
        LayerSpec("25 [nm]", [], "1.5"),
    ]
    stack = _make_layer_stack(layers)

    assert stack.entries[0].height_param == "H_layer_0"
    assert stack.entries[0].model_height_expr == "H_layer_0/2"
    assert stack.entries[0].z_lo_expr == "0 [um]"
    assert stack.entries[0].z_hi_expr == "H_layer_0/2"
    assert stack.entries[1].model_height_expr == "H_layer_1"
    assert stack.entries[1].z_lo_expr == "H_layer_0/2"
    assert stack.entries[1].z_hi_expr == "H_layer_0/2 + H_layer_1"
    assert stack.z_slab_top_expr == "H_layer_0/2 + H_layer_1 + H_layer_2"

    print("[Pass] layer stack uses half first layer and full later layers")


def test_finite_size_rejects_k():
    boundary = BoundarySpec("square", a=1.0)
    layers = [LayerSpec("200 [nm]", [], "3.3")]
    try:
        _validate_build_inputs(boundary, layers, "finite_size", {"kx": 0.0, "ky": 0.0})
        assert False, "Expected finite_size with k to fail"
    except ValueError as exc:
        assert "finite_size" in str(exc)

    print("[Pass] finite_size rejects k")


def test_finite_quarter_requires_quarter_boundary():
    layers = [LayerSpec("200 [nm]", [], "3.3")]
    boundary, _, mode, k, use_floquet = _validate_build_inputs(
        BoundarySpec("quarter_hexagon_boundary", a=1.0),
        layers,
        "finite_quarter",
        None,
    )
    assert boundary.shape == "quarter_hexagon_boundary"
    assert mode == "finite_quarter"
    assert k == {"kx": 0.0, "ky": 0.0}
    assert use_floquet is False
    rectangle, _, mode, k, use_floquet = _validate_build_inputs(
        BoundarySpec("quarter_rectangle", a=2.0, b=3.0),
        layers,
        "finite_quarter",
        None,
    )
    assert rectangle.shape == "quarter_rectangle"
    assert mode == "finite_quarter"
    assert k == {"kx": 0.0, "ky": 0.0}
    assert use_floquet is False

    try:
        _validate_build_inputs(
            BoundarySpec("hexagon_boundary", a=1.0),
            layers,
            "finite_quarter",
            None,
        )
        assert False, "Expected finite_quarter with full boundary to fail"
    except ValueError as exc:
        assert "quarter_hexagon_boundary" in str(exc)


def test_quarter_field_reconstruction_uses_pec_pmc_component_parity():
    coordinates = np.array([[0.0, 0.0], [1.0, 1.0]])
    values = np.array([3.0, 5.0])
    full_coordinates, full_values = reconstruct_quarter_scalar_field(
        coordinates,
        values,
        x_reflection_sign=hz_reflection_sign("PEC"),
        y_reflection_sign=hz_reflection_sign("PMC"),
    )
    field = {tuple(point): value for point, value in zip(full_coordinates, full_values)}
    assert field[(1.0, 1.0)] == 5.0
    assert field[(1.0, -1.0)] == 5.0
    assert field[(-1.0, 1.0)] == -5.0
    assert field[(-1.0, -1.0)] == -5.0
    assert field[(0.0, 0.0)] == 0.0

    quarter_e = np.ones((3, 2, 2), dtype=complex)
    quarter_mask = np.ones((2, 2), dtype=bool)
    full_e, full_mask = reconstruct_quarter_electric_grid(
        quarter_e,
        quarter_mask,
        x_boundary="PEC",
        y_boundary="PMC",
    )
    assert full_e.shape == (3, 3, 3)
    assert full_mask.all()
    assert np.allclose(
        full_e[:, 0, 2],
        electric_reflection_matrix("x", "PEC") @ np.ones(3),
    )
    assert np.allclose(
        full_e[:, 2, 0],
        electric_reflection_matrix("y", "PMC") @ np.ones(3),
    )


def test_old_api_shape_rejected_before_model_access():
    runner = SimulationRun.__new__(SimulationRun)
    try:
        runner.build_geometry(1.0, [triangle()], "unit_cell")
        assert False, "Expected old build_geometry(a, holes, ...) shape to fail"
    except TypeError as exc:
        assert "old build_and_run" in str(exc)

    try:
        _validate_build_inputs(BoundarySpec("hexagon", a=1.0), [triangle()], "unit_cell", None)
        assert False, "Expected old holes list as layers to fail"
    except TypeError as exc:
        assert "LayerSpec" in str(exc)

    print("[Pass] old API shape rejected before model access")


def test_finite_quarter_boundary_features_toggle_without_rebuild():
    states = {}

    class Feature:
        def __init__(self, tag):
            self.tag = tag

        def active(self, value):
            states[self.tag] = value

    class Physics:
        def feature(self, tag):
            return Feature(tag)

    class Component:
        def physics(self, tag):
            assert tag == "ewfd"
            return Physics()

    class Java:
        def component(self, tag):
            assert tag == "comp1"
            return Component()

    runner = SimulationRun.__new__(SimulationRun)
    runner.model = type("Model", (), {"java": Java()})()
    runner.set_finite_quarter_symmetry("PMC", "PEC")

    assert states == {
        "bc_symmetry_x_axis_pec": False,
        "bc_symmetry_x_axis_pmc": True,
        "bc_symmetry_y_axis_pec": True,
        "bc_symmetry_y_axis_pmc": False,
    }
    runner.model = None


def test_reused_finite_quarter_model_upgrades_boundaries_and_solver(monkeypatch):
    events = []
    monkeypatch.setattr("comsol_workflow.simulation_utils.JInt", int)

    class Selection:
        def __init__(self, tag):
            self.tag = tag

        def named(self, value):
            events.append(("selection", self.tag, value))

    class Feature:
        def __init__(self, tag):
            self.tag = tag

        def active(self, value):
            events.append(("active", self.tag, value))

        def selection(self):
            return Selection(self.tag)

        def set(self, key, value):
            events.append(("set", self.tag, key, value))

    class Features:
        def __init__(self, tags):
            self.items = {tag: Feature(tag) for tag in tags}

        def __call__(self, tag=None):
            return self if tag is None else self.items[tag]

        def feature(self, tag=None):
            return self(tag)

        def hasTag(self, tag):
            return tag in self.items

        def create(self, tag, kind, dimension):
            events.append(("create", tag, kind, dimension))
            self.items[tag] = Feature(tag)

    physics = Features(("bc_symmetry_x_axis", "bc_symmetry_y_axis"))
    solution_features = Features(("e1",))

    class Component:
        def physics(self, tag):
            assert tag == "ewfd"
            return physics

    class Study:
        def feature(self, tag):
            assert tag == "eig"
            return Feature(tag)

    class Solution:
        def feature(self, tag=None):
            return solution_features(tag)

    class Tags:
        def __init__(self, tags):
            self.tags = set(tags)

        def hasTag(self, tag):
            return tag in self.tags

    class Result:
        def dataset(self):
            return Tags(("dset1", "cpl1", "cpl2", "cpl3", "cpl4"))

        def numerical(self):
            return Tags(("gev1", "int1", "int_grid", "max3d"))

    class Java:
        def component(self, tag):
            assert tag == "comp1"
            return Component()

        def study(self, tag=None):
            return Tags(("std1",)) if tag is None else Study()

        def sol(self, tag=None):
            return Tags(("sol1",)) if tag is None else Solution()

        def result(self):
            return Result()

    runner = SimulationRun.__new__(SimulationRun)
    runner.config = type(
        "Config",
        (),
        {"eigenfrequency_shift": "203.6 [THz]", "eigenmode_count": 6},
    )()
    runner.model = type("Model", (), {"java": Java()})()

    result = runner.prepare_reused_finite_quarter_model()

    created_tags = {event[1] for event in events if event[0] == "create"}
    assert created_tags == {
        "bc_symmetry_x_axis_pec",
        "bc_symmetry_x_axis_pmc",
        "bc_symmetry_y_axis_pec",
        "bc_symmetry_y_axis_pmc",
    }
    assert ("active", "bc_symmetry_x_axis", False) in events
    assert ("active", "bc_symmetry_y_axis", False) in events
    assert ("set", "eig", "shift", "203.6 [THz]") in events
    assert ("set", "e1", "neigs", 6) in events
    assert result["eigenmode_count"] == 6
    runner.model = None


def test_solution_rerun_and_copy_use_working_solution():
    events = []

    class WorkingSolution:
        def clearSolutionData(self):
            events.append("clear")

        def runAll(self):
            events.append("run")

        def copySolution(self, destination):
            events.append(("copy_solution", "sol1", destination))

        def label(self, value):
            events.append(("solution_label", value))

    class SolutionList:
        def __init__(self):
            self.working = WorkingSolution()

        def __call__(self, tag=None):
            return self if tag is None else self.working

        def hasTag(self, tag):
            return False

    class Dataset:
        def set(self, key, value):
            events.append(("dataset_set", key, value))

        def label(self, value):
            events.append(("dataset_label", value))

    class DatasetList:
        def __init__(self):
            self.dataset = Dataset()

        def __call__(self, tag=None):
            return self if tag is None else self.dataset

        def hasTag(self, tag):
            return False

        def create(self, tag, kind):
            events.append(("dataset_create", tag, kind))

    class Result:
        def __init__(self):
            self.datasets = DatasetList()

        def __call__(self):
            return self

        def dataset(self, tag=None):
            return self.datasets(tag)

    class Java:
        def __init__(self):
            self.solutions = SolutionList()
            self.results = Result()

        def sol(self, tag=None):
            return self.solutions(tag)

        def result(self):
            return self.results

    runner = SimulationRun.__new__(SimulationRun)
    runner.model = type("Model", (), {"java": Java()})()
    runner.rerun_eigenfrequency_solution()
    copied = runner.copy_solution(
        "solsym2",
        "dsetsym2",
        label="Internal symmetry solution 2",
    )

    assert events[:2] == ["clear", "run"]
    assert ("copy_solution", "sol1", "solsym2") in events
    assert ("dataset_create", "dsetsym2", "Solution") in events
    assert ("dataset_set", "solution", "solsym2") in events
    assert copied == {
        "solution_tag": "solsym2",
        "dataset_tag": "dsetsym2",
    }
    runner.model = None


if __name__ == "__main__":
    test_boundary_rectangle_validation()
    test_hexagon_and_square_reject_b()
    test_hexagon_boundary_footprint_uses_left_right_vertices()
    test_rectangle_wall_seeds_and_k_scale()
    test_hexagon_boundary_wall_seeds_group_opposite_edges()
    test_single_parameter_boundaries_have_only_ga()
    test_layer_stack_first_layer_half_height()
    test_finite_size_rejects_k()
    test_old_api_shape_rejected_before_model_access()
    print("\nAll simulation_utils multilayer tests passed.")
