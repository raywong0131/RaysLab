import sys
import numpy as np
import pytest

from comsol_workflow.s4_boundary import (
    EPS0, area_rule, complex_columns, edge_integrals, group_totals, hole_edges,
    maxwell_terms, parity_errors, phase_and_norm, polygon_area, polygon_edges,
    length_only_prediction, segment_minimum,
)


def geometry():
    theta = np.pi / 6 + np.arange(6) * np.pi / 3
    outer = 0.82e-6 / np.sqrt(3) * np.column_stack([np.cos(theta), np.sin(theta)])
    holes = []
    for k in range(6):
        t = np.arange(3) * 2 * np.pi / 3 + np.pi
        p = 0.07e-6 * np.column_stack([np.cos(t), np.sin(t)]) + [0.25e-6, 0]
        a = k * np.pi / 3
        rotation = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        holes.append(p @ rotation.T)
    return outer, holes


def test_material_area_partition():
    outer, holes = geometry()
    xy, w, air = area_rule(outer, holes, 4, 1)
    assert np.all(w > 0)
    assert np.isclose(sum(w), polygon_area(outer), rtol=1e-12, atol=0)
    assert np.isclose(sum(w[air]), sum(polygon_area(h) for h in holes), rtol=1e-12, atol=0)
    assert np.isclose(w @ xy[:, 0], 0, atol=1e-32)


def test_normals_numbering_and_constant_field_closure():
    _, holes = geometry()
    for h in (holes, [p[::-1] for p in holes]):
        edges = hole_edges(h)
        assert len(edges) == 18
        rows = edge_integrals(edges, lambda xy: np.full(len(xy), 1 + 2j), 8)
        for k in range(1, 7):
            assert abs(sum(r["qx"] for r in rows if r["hole"] == k)) < 1e-20
            assert abs(sum(r["qy"] for r in rows if r["hole"] == k)) < 1e-20
        assert np.isclose(sum(r["length_m"] for r in rows if r["group"] == "B"),
                          2 * sum(r["length_m"] for r in rows if r["group"] == "A"))
    with pytest.raises(ValueError, match="numbering"):
        hole_edges(holes[1:] + holes[:1])


def test_full_area_boundary_identity_complex_field_and_3d_terms():
    outer, holes = geometry()
    xy, w, air = area_rule(outer, holes, 5)
    length = 0.82e-6
    hz = lambda p: (1 + 0.7j) * (p[:, 0] / length + 0.3 * p[:, 1] / length) + 0.4j
    eps = np.where(air, 1.0, 3.3**2)
    omega = 2 * np.pi * (198e12 - 1e8j)
    dx = np.full(len(xy), (0.2 + 0.1j) / length)
    dy = np.full(len(xy), (-0.4 + 0.2j) / length)
    ex = 1j / (omega * EPS0 * eps) * ((1 + 0.7j) * 0.3 / length - dy)
    ey = 1j / (omega * EPS0 * eps) * (dx - (1 + 0.7j) / length)
    q = group_totals(edge_integrals(hole_edges(holes), hz, 8))
    boundary = edge_integrals(polygon_edges(outer), hz, 8)
    ox = sum(r["qx"] for r in boundary) / 3.3**2
    oy = sum(r["qy"] for r in boundary) / 3.3**2
    terms = maxwell_terms(q["Qx"], q["Qy"], ox, oy, w @ (dx / eps), w @ (dy / eps),
                          omega, sum(w), 1, 3.3**2)
    assert np.isclose(w @ ex / sum(w), terms["predicted_Ex"], rtol=1e-12)
    assert np.isclose(w @ ey / sum(w), terms["predicted_Ey"], rtol=1e-12)
    assert abs(terms["Ox"]) > 0 and abs(terms["Zx"]) > 0


def test_periodic_outer_cancellation_independent_of_holes():
    size = 1e-6
    square = size * np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
    h = lambda p: np.sin(2 * np.pi * p[:, 0] / size) + 1j * np.cos(2 * np.pi * p[:, 1] / size)
    b = edge_integrals(polygon_edges(square), h, 32)
    assert abs(sum(r["qx"] for r in b)) < 1e-20
    assert abs(sum(r["qy"] for r in b)) < 1e-20


def test_phase_and_same_normalization_for_both_fields():
    ref = np.array([1, -2, 3], complex)
    phase = np.exp(1.2j)
    h = 7 * phase * ref
    w = np.array([1, 2, 1], float)
    factor = phase_and_norm(h, w, ref)
    normalized = h * factor
    assert np.isclose(w @ abs(normalized)**2 / sum(w), 1)
    assert np.max(abs(normalized.imag)) < 1e-14
    # A radiation zero is never the normalization anchor.
    assert 0j * factor == 0j


def test_parity_is_a_measurement_not_projection():
    h = np.array([1 + .1j, 2 - .3j])
    assert max(parity_errors(h, -h, h).values()) == 0
    assert parity_errors(h, h, h)["odd_x_error"] == 2
    assert np.array_equal(h, [1 + .1j, 2 - .3j])


def test_group_sum_uses_six_and_twelve_edges():
    _, holes = geometry()
    rows = [{"group": e.group, "qx": (-2 + 1j if e.group == "A" else 1 - .5j), "qy": 0j}
            for e in hole_edges(holes)]
    g = group_totals(rows)
    assert g["Qx"] == 0 and g["cancellation_ratio"] == 0
    assert np.isclose(abs(g["AB_phase_rad"]), np.pi)
    with pytest.raises(ValueError):
        group_totals(rows[:-1])


def test_serialized_complex_fields_keep_both_components():
    assert complex_columns({"q": 2 - 3j}) == {"q_re": 2, "q_im": -3}


def test_pure_module_does_not_load_mph():
    assert "mph" not in sys.modules


def test_length_only_uses_actual_reference_not_zeta_equal_one():
    _, holes = geometry()
    edges = hole_edges(holes)
    reference = edge_integrals(edges, lambda xy: 1 + xy[:, 0] / 1e-6, 8)
    enlarged = [h.mean(axis=0) + 1.16 * (h - h.mean(axis=0)) for h in holes]
    predicted = length_only_prediction(reference, hole_edges(enlarged))
    for r, p in zip(reference, predicted):
        assert np.isclose(p["qx"], 1.16 * r["qx"], atol=1e-22)


def test_real_sign_change_is_not_a_complex_zero():
    estimate = segment_minimum(1, 2, -1 + .1j, 1 + .1j)
    assert estimate["zeta_estimate"] == 1.5
    assert estimate["real_sign_change"]
    assert estimate["complex_residual"] == .1
    exact = segment_minimum(1, 2, -1 - 2j, 1 + 2j)
    assert exact["complex_residual"] < 1e-14


def test_preparation_preserves_shared_parameter_and_imports_no_solver(tmp_path, monkeypatch):
    from scripts.run_main import run_s4_boundary as entry
    monkeypatch.setattr(entry, "TASK_ROOT", tmp_path / "tmp")
    monkeypatch.setattr(entry, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path / "formal")
    parameter = tmp_path / "parameter.json"
    data = entry.json_read(entry.parameter_path_from_environment())
    data["unit_cell_2d"]["zeta"] = 1.15
    entry.write_json(parameter, data)
    before = entry.digest(parameter)
    entry.main(["prepare", "--run-id", "test-s4", "--parameter-file", str(parameter)])
    directory = tmp_path / "tmp/s4_test-s4"
    manifest = entry.prepared_manifest(directory)
    active = entry.load_shared_parameters(parameter)
    assert manifest["cases"][0]["zeta"] == 1.156
    assert manifest["source_zeta"] == active.unit_cell_2d.cell.zeta
    assert manifest["cases"][0]["mesh"] == active.mesh_size
    assert manifest["execution"]["solve_count"] == 1
    assert "mph" not in sys.modules
    assert entry.digest(parameter) == before
    entry.report(directory)
    assert (directory / "80_logs/acceptance.md").is_file()
    assert not (directory / "acceptance.md").exists()
    assert not (tmp_path / "formal").exists()
    assert "mph" not in sys.modules
    with pytest.raises(FileExistsError):
        entry.main(["prepare", "--run-id", "test-s4"])


@pytest.mark.parametrize("case_data_dir", ["01_results", "80_logs"])
def test_run_report_keeps_tables_in_run_logs(tmp_path, case_data_dir):
    from scripts.run_main import run_s4_boundary as entry
    import pandas as pd

    cases = [{"case_id": f"case{i}", "eta": .96, "mesh": 5, "zeta": z}
             for i, z in enumerate((1.15, 1.16))]
    entry.write_json(tmp_path / "99_config/s4_manifest.json", {"run_id": "layout", "cases": cases})
    protected = {}
    for i, case in enumerate(cases):
        case_dir = tmp_path / case["case_id"]
        entry.write_json(case_dir / "99_config/s4_summary.json", case)
        path = case_dir / case_data_dir / "group_totals.csv"
        entry.write_csv(path, [{"refinement": 2, "z_over_H": 0., "Qx": complex(2 * i - 1, .1)}])
        protected[path] = entry.digest(path)
    entry.report(tmp_path)
    logs = tmp_path / "80_logs"
    assert {p.name for p in logs.iterdir()} == {"case_table.csv", "zero_estimates.csv", "acceptance.md"}
    assert len(pd.read_csv(logs / "case_table.csv")) == 2
    estimate = pd.read_csv(logs / "zero_estimates.csv").iloc[0]
    assert np.isclose(estimate.zeta_estimate, 1.155)
    assert np.isclose(estimate.complex_residual, .1)
    assert not list(tmp_path.glob("*.csv")) and not list(tmp_path.glob("*.md"))
    assert all(entry.digest(path) == sha for path, sha in protected.items())
    assert "mph" not in sys.modules


def test_prepare_budget_and_section_guard(tmp_path, monkeypatch):
    from scripts.run_main import run_s4_boundary as entry
    monkeypatch.setattr(entry, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(entry, "TASK_ROOT", tmp_path / "tmp")
    with pytest.raises(ValueError, match="budget"):
        entry.main(["prepare", "--zeta", "1.154", "1.158"])
    with pytest.raises(SystemExit):
        entry.main(["prepare", "--z-nm", "50"])


def test_s4_extraction_preserves_nonzero_vector_field_thickness_term(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import pandas as pd
    from scripts.run_main import run_s4_boundary as entry

    outer, holes = geometry()
    omega = 2*np.pi*(198e12-1e8j)
    derivative = 1e6*(1+.2j)
    ey = 1j*derivative/(omega*EPS0)

    class AnalyticSampler:
        def plane(self, mode, xy, z, expressions):
            values = {"x": xy[:, 0], "y": xy[:, 1], "z": z, "dom": 1,
                      "ewfd.Hz": 1, "ewfd.Ex": 0j, "ewfd.Ey": ey,
                      "d(ewfd.Hx,z)": 0, "d(ewfd.Hy,z)": 0,
                      "d(laginterp(2,ewfd.Hx),z)/(1[A/m^2])": derivative,
                      "d(laginterp(2,ewfd.Hy),z)/(1[A/m^2])": 0}
            return np.column_stack([np.broadcast_to(values[e], len(xy)) for e in expressions])

    params = {"a": .82e-6, "H": 200e-9, "H_air": 1.55e-6}
    runner = SimpleNamespace(model=SimpleNamespace(java=SimpleNamespace(param=lambda: SimpleNamespace(evaluate=params.get))),
                             config=SimpleNamespace(air_refractive_index=1, slab_refractive_index=1))
    monkeypatch.setattr(entry, "export_area_field", lambda *args: None)
    settings = {"z_nm": [0], "area_order": 2, "area_subdivisions": 1, "edge_order": 4,
                "trace_offset_over_a": 1e-5, "identity_tolerance": 1e-3}
    result = entry.evaluate_case(runner, AnalyticSampler(), 1, omega, outer, holes, settings, 1+0j, tmp_path)
    rows = pd.read_csv(tmp_path / "area_comparison.csv")
    assert np.allclose(rows.Zy_re+1j*rows.Zy_im, ey, rtol=1e-12)
    assert np.allclose(rows.Ey_area_re+1j*rows.Ey_area_im, ey, rtol=1e-12)
    assert result["derivative_evaluation"] == "elementwise_laginterp2"
    assert result["identity_gate"]


def test_native_s31_integral_preserves_units_mode_and_complex_phase(monkeypatch):
    from types import SimpleNamespace
    from scripts.analysis.audit_s31_boundary import native_integral

    class Node:
        def __init__(self):
            self.settings = {}
            self.real = [[2.0], [3.0]]

        def set(self, key, value=None):
            if value is None:
                self.entities = key
            else:
                self.settings[key] = value

        def selection(self):
            return self

        def getReal(self):
            return self.real

        def getImag(self):
            return [[0.0], [4.0]]

    node = Node()
    numerical = SimpleNamespace(tags=lambda: [], create=lambda tag, kind: node)
    model = SimpleNamespace(java=SimpleNamespace(result=lambda: SimpleNamespace(numerical=lambda: numerical)))
    monkeypatch.setitem(sys.modules, "jpype.types", SimpleNamespace(
        JArray=lambda kind: lambda values: values, JInt=int, JString=str))
    values = native_integral(model, "audit", "IntSurface", [3, 16], ["1", "ewfd.Ey"], ["m^2", "V*m"], 1, 8)
    assert np.array_equal(values, [2, 3-4j])
    assert node.entities == [3, 16]
    assert node.settings["unit"] == ["m^2", "V*m"]
    assert node.settings["solnum"] == [2]
    assert node.settings["intorder"] == 8
    node.real = [[float("nan")], [3.0]]
    with pytest.raises(ValueError, match="invalid shape or values"):
        native_integral(model, "audit", "IntSurface", [3], ["1", "ewfd.Ey"], ["m^2", "V*m"], 1, 8)


def test_sampler_preserves_imaginary_parts_and_explicit_si_units():
    from scripts.run_main.run_s4_boundary import Sampler

    class Node:
        def set(self, key, value):
            setattr(self, key, value)

        def setInterpolationCoordinates(self, value):
            self.coordinates = np.asarray(value)

        def getData(self):
            return np.array([[[1.0, 2.0]], [[3.0, 4.0]]])

        def getImagData(self):
            return np.array([[[.1, .2]], [[.3, .4]]])

        def isComplex(self):
            return True

    sample = Sampler.__new__(Sampler)
    sample.node = Node()
    sample.JArray = lambda kind, dimensions: lambda values: values
    sample.JString, sample.JInt, sample.JDouble = str, int, float
    xyz = np.array([[1e-7, 2e-7, 3e-8], [-1e-7, 2e-7, 3e-8]])
    values = sample(0, xyz, ["ewfd.Hz", "d(ewfd.Hx,z)"])
    assert np.allclose(values, [[1 - .1j, 3 - .3j], [2 - .2j, 4 - .4j]])
    assert sample.node.expr == ["(ewfd.Hz)/(1[A/m])", "(d(ewfd.Hx,z))/(1[A/m^2])"]
    assert np.allclose(sample.node.coordinates, xyz.T / 1e-6)


def test_run_rejects_tampering_conflicts_and_failed_preflight(tmp_path, monkeypatch):
    from scripts.run_main import run_s4_boundary as entry
    monkeypatch.setattr(entry, "TASK_ROOT", tmp_path / "tmp")
    monkeypatch.setattr(entry, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path / "formal")
    entry.main(["prepare", "--run-id", "guards"])
    directory = entry.TASK_ROOT / "s4_guards"
    manifest = entry.prepared_manifest(directory)
    output = entry.Path(manifest["output_directory"])
    calls = []

    def unavailable():
        calls.append(True)
        raise RuntimeError("test unavailable backend")

    monkeypatch.setattr(entry, "runtime_preflight", unavailable)
    with pytest.raises(RuntimeError, match="unavailable backend"):
        entry.main(["run", str(directory)])
    assert calls == [True]
    assert not output.parent.exists()
    assert entry.json_read(directory / "99_config/preflight.json")["status"] == "failed"
    snapshot = directory / manifest["cases"][0]["snapshot"]
    original = snapshot.read_bytes()
    snapshot.write_bytes(original + b" ")
    with pytest.raises(RuntimeError, match="snapshot changed"):
        entry.main(["run", str(directory)])
    snapshot.write_bytes(original)
    manifest_path = directory / "99_config/s4_manifest.json"
    original_manifest = manifest_path.read_bytes()
    manifest_path.write_bytes(original_manifest + b" ")
    with pytest.raises(RuntimeError, match="Manifest changed"):
        entry.main(["run", str(directory)])
    manifest_path.write_bytes(original_manifest)
    with pytest.raises(ValueError, match="task directory"):
        entry.main(["run", str(output)])
    output.mkdir(parents=True)
    marker = output / "keep.txt"
    marker.write_text("existing result")
    with pytest.raises(FileExistsError):
        entry.main(["run", str(directory)])
    assert marker.read_text() == "existing result"
    assert list(output.iterdir()) == [marker]
    assert calls == [True]


def test_check_model_keeps_audit_in_tmp_and_never_saves_or_solves(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from scripts.run_main import run_s4_boundary as entry
    monkeypatch.setattr(entry, "TASK_ROOT", tmp_path / "tmp")
    monkeypatch.setattr(entry, "UNIT_CELL_2D_OUTPUT_ROOT", tmp_path / "formal")
    source = tmp_path / "old.mph"
    source.write_bytes(b"protected checkpoint")
    before = entry.digest(source)
    model = SimpleNamespace(java=SimpleNamespace(param=lambda: SimpleNamespace(
        evaluate=lambda name: {"a": .82e-6, "H": 200e-9, "kx": .15, "ky": .15}[name])))
    removed = []
    client = SimpleNamespace(load=lambda path: model, remove=removed.append)
    monkeypatch.setitem(sys.modules, "mph", SimpleNamespace(start=lambda **kw: client))
    monkeypatch.setattr(entry, "runtime_preflight", lambda: {"status": "test"})
    monkeypatch.setattr(entry, "Sampler", lambda runner: lambda mode, xyz, expressions:
                        np.column_stack([xyz, np.ones((len(xyz), 6), complex) * (1 + 2j)]))
    entry.main(["check-model", str(source)])
    audit, = entry.TASK_ROOT.glob("s4_backend_check_*/backend_check.json")
    assert entry.json_read(audit)["status"] == "passed"
    assert entry.json_read(audit)["source_unchanged"]
    assert removed == [model]
    assert entry.digest(source) == before
    assert not (tmp_path / "formal").exists()


def test_s4_plot_contract_complex_sums_and_missing_edges(tmp_path):
    from types import SimpleNamespace
    from matplotlib.image import imread
    from scripts.run_main import run_s4_boundary as entry
    from scripts.analysis.plot_s4_boundary import main as plot_main
    from comsol_workflow.s4_plotting import build_s4_figures

    outer, holes = geometry()
    edges = hole_edges(holes)
    output = tmp_path / "01_results"
    geometry_rows = [{"hole": e.hole, "edge": e.edge, "group": e.group,
                      "x0_m": e.start[0], "y0_m": e.start[1], "x1_m": e.end[0], "y1_m": e.end[1],
                      "length_m": e.length, "nx_dielectric_to_air": e.normal[0],
                      "ny_dielectric_to_air": e.normal[1]} for e in edges]
    edge_rows, group_rows = [], []
    for section in (0., .25, .5):
        for refinement in (1, 2):
            rows = edge_integrals(edges, lambda xy: (1 + section) * (
                xy[:, 0] / 1e-6 + 1j * .03 * xy[:, 1] / 1e-6), 8 * refinement)
            edge_rows.extend({"z_over_H": section, "refinement": refinement, **r} for r in rows)
            group_rows.append({"z_over_H": section, "refinement": refinement, **group_totals(rows)})
    entry.write_csv(output / "edge_geometry.csv", geometry_rows[::-1])
    entry.write_csv(output / "edge_integrals.csv", edge_rows[::-1])
    entry.write_csv(output / "group_totals.csv", group_rows[::-1])
    area_m2 = abs(polygon_area(outer))
    prefactor, jump = (2 + 3j) * 1e7, .9
    area_rows = [{"z_over_H": z, "refinement": r, "Ex_area": 1 + 2j, "Ey_area": 3 - 4j,
                  "area_samples": 256 * r,
                  "area_m2": area_m2, "Ex_integral": (1 + 2j) * area_m2,
                  "Ey_integral": (3 - 4j) * area_m2, "normalization_factor": 1 + 0j,
                  "prefactor": prefactor, "jump": jump,
                  "Lx": .3 + .2j, "Ly": .4j, "Ox": 1j, "Oy": 2j,
                  "Zx": 0j, "Zy": 0j, "predicted_Ex": .3 + 1.2j, "predicted_Ey": 2.4j}
                 for z in (0., .25, .5) for r in (1, 2)]
    entry.write_csv(output / "area_comparison.csv", area_rows)
    entry.write_json(tmp_path / "99_config/s4_summary.json", {
        "frequency_thz": 198.12345, "mode_idx": 1, "b0_nm": 245., "eta": .96, "zeta": 1.156, "mesh": 5})
    with pytest.raises(ValueError, match="Missing area_field_z0"):
        build_s4_figures(tmp_path)
    xy, weights, _ = area_rule(outer, holes, 4, 2)
    def electric_samples(mode, points, z, expressions):
        assert mode == 1 and z == 0 and expressions == ['ewfd.Ex', 'ewfd.Ey']
        return np.tile([1 + 2j, 3 - 4j], (len(points), 1))
    field_path = entry.export_area_field(output, SimpleNamespace(plane=electric_samples), 1, outer,
                                         xy, weights, electric_samples(1, xy, 0, ['ewfd.Ex', 'ewfd.Ey']), 1, 2)
    before = {p: entry.digest(p) for p in output.iterdir() if p.suffix in ('.csv', '.npz')}
    figures = build_s4_figures(tmp_path)
    assert set(figures) == {"s4_edges_py", "s4_edges_py_Ex_reference", "s4_group_contributions", "s4_area_comparison"}
    scale = max(abs(r["integral_hz"]) for r in edge_rows if r["refinement"] == 2 and r["z_over_H"] == 0)
    for axis in ("x", "y"):
        fig = figures['s4_edges_py' if axis == 'x' else 's4_edges_py_Ex_reference']
        channel_label = r"$p_y$: $E_y$ integral" if axis == "x" else r"$p_y$: $E_x$ reference"
        assert channel_label in fig._suptitle.get_text()
        field_ax, spatial, contribution_ax = fig.axes
        fine = [r for r in edge_rows if r["refinement"] == 2 and r["z_over_H"] == 0]
        assert spatial.get_aspect() == 1
        assert spatial.collections[0].get_cmap().name == "RdBu_r"
        conversion = (1 if axis == 'x' else -1) * prefactor * jump / (3 - 4j)
        expected_edges = np.array([conversion * r['q'+axis] for r in fine])
        assert spatial.collections[0].get_clim() == contribution_ax.get_ylim()
        edge_limit = contribution_ax.get_ylim()[1]
        all_edge_max = max(abs((sign * prefactor * jump / (3-4j) * r['q'+key]).real)
                           for key, sign in [('x', 1), ('y', -1)] for r in fine)
        assert all_edge_max < edge_limit < 2*all_edge_max
        assert spatial.child_axes[0].get_ylim() == contribution_ax.get_ylim()
        assert fig.texts == [fig._suptitle]
        assert np.allclose(spatial.collections[0].get_array(), expected_edges.real)
        heights = np.array([p.get_height() for p in contribution_ax.patches])
        assert np.allclose(heights, expected_edges.real)
        assert heights.min() < 0 < heights.max()
        assert np.isclose(heights.sum(), expected_edges.sum().real)
        assert contribution_ax.get_ylim()[0] == -contribution_ax.get_ylim()[1]
        assert len(contribution_ax.get_yticks()) <= 5
        assert contribution_ax.get_ylabel() == r"Re($c_e$) (1)"
        assert len(fig.axes) == 3
        fig.canvas.draw()
        left_box, right_box, bar_box = [ax.get_window_extent() for ax in fig.axes]
        assert np.isclose(bar_box.width / bar_box.height, 4 / 3)
        assert np.isclose(right_box.x0 - left_box.x1, bar_box.x0 - right_box.x1)
        assert np.allclose([left_box.y0, left_box.y1], [bar_box.y0, bar_box.y1])
        assert np.allclose([left_box.y0, left_box.y1], [right_box.y0, right_box.y1])
        electric_axis, expected_field = ('y', 1) if axis == 'x' else ('x', -.2)
        assert field_ax.get_aspect() == 1
        assert field_ax.get_title() == rf"(a) Area field $E_{electric_axis}$"
        assert field_ax.collections[0].get_cmap().name == 'RdBu_r'
        assert np.allclose(field_ax.collections[0].get_array().compressed(), expected_field)
        expected_integral = ((3 - 4j) if axis == 'x' else (1 + 2j)) / (3 - 4j)
        value_label = '1' if axis == 'x' else f'{expected_integral.real:.4g} {expected_integral.imag:+.4g}i'
        assert value_label in field_ax.texts[0].get_text()
        assert r'\int_\Omega' in field_ax.texts[0].get_text()
        total = expected_edges.sum()
        assert f'{total.real:.4g} {total.imag:+.4g}i' in contribution_ax.texts[-1].get_text()
        assert r'\sum_{e=1}^{18}c_e' in contribution_ax.texts[-1].get_text()
        samples = electric_samples(1, xy, 0, ['ewfd.Ex', 'ewfd.Ey'])[:, 1 if axis == 'x' else 0]
        assert np.isclose((weights / area_m2) @ (samples / (3 - 4j)), expected_integral)
        assert field_ax.child_axes[0].get_ylabel() == rf"Re($F_{electric_axis}=S E_{electric_axis}/I_y$) (1)"
        assert spatial.get_title() == "(b) Boundary geometry"
        assert contribution_ax.get_title() == "(c) Edge contributions"
    assert figures['s4_edges_py'].axes[2].get_ylim() == figures['s4_edges_py_Ex_reference'].axes[2].get_ylim()
    from matplotlib.text import Text
    for fig in figures.values():
        labels = ' '.join(text.get_text() for text in fig.findobj(Text))
        assert not any(old in labels for old in ('qx', 'qy', 'Qx', 'Qy', 'q_{e,', 'Q_x', 'Q_y'))
    fine_group = next(r for r in group_rows if r["refinement"] == 2 and r["z_over_H"] == 0)
    for component, offset in (("x", 0), ("y", 2)):
        for part, col in ((np.real, 0), (np.imag, 1)):
            ax = figures['s4_group_contributions'].axes[offset + col]
            assert np.allclose([p.get_height() for p in ax.patches],
                               [part(fine_group['Q'+g+component]) / scale for g in ('A', 'B', '')])
    assert np.allclose([p.get_height() for p in figures['s4_area_comparison'].axes[0].patches], [1, .3, 0, 0, .3])
    assert np.allclose([p.get_height() for p in figures['s4_area_comparison'].axes[1].patches], [2, .2, 1, 0, 1.2])
    for fig in figures.values():
        fig.clear()
    plot_main([str(tmp_path)])
    assert len(list(output.glob("*.png"))) == len(list(output.glob("*.svg"))) == 5
    image = imread(output / "s4_edges_py.png")
    assert image.shape[1] > 1000 and np.std(image[:, :, :3]) > .01
    untouched = {p: (entry.digest(p), p.stat().st_mtime_ns) for p in output.glob('s4_*')
                 if not p.name.startswith('s4_edges_')}
    plot_main([str(tmp_path), '--edges-only'])
    assert all((entry.digest(p), p.stat().st_mtime_ns) == saved for p, saved in untouched.items())
    assert all(entry.digest(p) == sha for p, sha in before.items())
    assert "mph" not in sys.modules
    original_field = field_path.read_bytes()
    with np.load(field_path) as saved:
        altered_field = dict(saved)
    altered_field['electric_vm'] *= 2
    np.savez_compressed(field_path, **altered_field)
    with pytest.raises(ValueError, match="quadrature disagrees"):
        build_s4_figures(tmp_path)
    field_path.write_bytes(original_field)
    for row in area_rows:
        row['Ey_area'] = row['Ey_integral'] = 0j
    entry.write_csv(output / 'area_comparison.csv', area_rows)
    altered_field['electric_vm'] = np.tile([1 + 2j, 0j], (len(xy), 1))
    altered_field['grid_electric_vm'][:, :, 1] = 0j
    np.savez_compressed(field_path, **altered_field)
    with pytest.raises(ValueError, match='zero or numerically unresolved'):
        build_s4_figures(tmp_path)
    field_path.write_bytes(original_field)
    for row in area_rows:
        row['Ey_area'], row['Ey_integral'] = 3 - 4j, (3 - 4j) * area_m2
    entry.write_csv(output / 'area_comparison.csv', area_rows)
    entry.write_csv(output / "edge_integrals.csv", edge_rows[:-1])
    with pytest.raises(ValueError, match="Missing edge"):
        build_s4_figures(tmp_path)
    entry.write_csv(output / "edge_integrals.csv", edge_rows)
    group_rows[-1]["Qx"] += scale
    entry.write_csv(output / "group_totals.csv", group_rows)
    with pytest.raises(ValueError, match="disagrees"):
        build_s4_figures(tmp_path)

    # The compact layout must remain reproducible without restoring auxiliary plots.
    group_rows[-1]["Qx"] -= scale
    entry.write_csv(output / "group_totals.csv", group_rows)
    logs = tmp_path / "80_logs"
    output.rename(logs)
    saved = {p: entry.digest(p) for p in logs.iterdir() if p.is_file()}
    overview = tmp_path / "10_overview"
    overview.mkdir()
    diagnostic = overview / "s31_diagnostic_py.png"
    diagnostic.write_bytes(b"preserve independent audit figure")
    assert entry.data_directory(tmp_path) == logs
    plot_main([str(tmp_path)])
    assert {p.name for p in overview.iterdir()} == {
        's4_edges_py.png', 's4_edges_py.svg', 's31_diagnostic_py.png'}
    assert diagnostic.read_bytes() == b"preserve independent audit figure"
    assert not output.exists()
    assert all(entry.digest(p) == sha for p, sha in saved.items())
    plot_main([str(tmp_path), '--validation-only'])
    assert (overview / 's4_validation_py.png').is_file()
    assert 'mph' not in sys.modules


def test_s31_diagnostic_preserves_common_reference_and_complex_residual(tmp_path):
    from scripts.run_main import run_s4_boundary as entry
    from scripts.analysis.plot_s4_boundary import main as plot_main
    from comsol_workflow.s4_plotting import build_s4_validation_figure
    reference = 1 + 2j
    rows = []
    for refinement, a, b in ((1, -reference, 3 * reference), (2, reference, (-2 + .3j) * reference)):
        outer = .4 * reference
        rows.append({'z_over_H': 0, 'refinement': refinement, 'area_samples': 100 * refinement,
                     'Ey_area': a, 'Ly': b, 'Oy': outer, 'Zy': 0j, 'predicted_Ey': b + outer})
    path = tmp_path / '01_results/area_comparison.csv'
    entry.write_csv(path, rows)
    entry.write_json(tmp_path / '99_config/s4_summary.json', {'frequency_thz': 198., 'mesh': 5})
    before = entry.digest(path)
    figure = build_s4_validation_figure(tmp_path)
    comparison, precision, balance = figure.axes
    assert np.allclose([p.get_height() for p in comparison.patches], [1, -2])
    lines = {line.get_label(): line.get_ydata() for line in precision.lines}
    assert np.allclose(lines['Area A'], [-1, 1])  # coarse must not be independently normalized to 1
    assert np.allclose(lines['Boundary B'], [3, -2])
    assert np.allclose([p.get_height() for p in balance.patches], [1, -2, .4, 0, -1.6])
    assert f'{abs(3 - .3j):.4g}' in comparison.texts[-1].get_text()
    assert 'changes sign' in precision.texts[-1].get_text()
    assert 'derivative unvalidated' in balance.texts[-1].get_text()
    assert f'{abs(2.6 - .3j):.4g}' in balance.texts[-1].get_text()
    figure.canvas.draw()
    boxes = [ax.get_window_extent() for ax in figure.axes]
    assert all(np.isclose(box.width / box.height, 4 / 3) for box in boxes)
    assert np.isclose(boxes[1].x0 - boxes[0].x1, boxes[2].x0 - boxes[1].x1)
    figure.clear()
    plot_main([str(tmp_path), '--validation-only'])
    assert sorted(p.name for p in path.parent.glob('*.png')) == ['s4_validation_py.png']
    assert entry.digest(path) == before and 'mph' not in sys.modules
    rows[-1]['Ey_area'] = 0j
    entry.write_csv(path, rows)
    figure = build_s4_validation_figure(tmp_path)
    assert figure.axes[0].get_ylabel() == 'Re(value) (V/m)'
    assert figure.axes[0].patches[0].get_height() == 0
    figure.clear()
    rows[-1]['predicted_Ey'] += 1
    entry.write_csv(path, rows)
    with pytest.raises(ValueError, match='reconstruction disagrees'):
        build_s4_validation_figure(tmp_path)


def test_s4_closed_half_slab_sections_and_material_side():
    from types import SimpleNamespace
    from scripts.run_main.run_s4_boundary import Sampler, section_coordinates
    sections = section_coordinates({"z_nm": [0, 50, 100]}, 200e-9)
    assert np.allclose([s[1] for s in sections], [0, 50e-9, 100e-9], atol=1e-20)
    for invalid in ([1, 50], [0, 101], [0, float('nan')], [0, 0]):
        with pytest.raises(ValueError):
            section_coordinates({"z_nm": invalid}, 200e-9)
    selected = []
    selection = SimpleNamespace(set=lambda ids: selected.append(list(ids)), all=lambda: selected.append('all'))
    class PlaneProbe(Sampler):
        def __call__(self, mode, xyz, expressions):
            return xyz
    sampler = PlaneProbe.__new__(PlaneProbe)
    sampler.node = SimpleNamespace(selection=lambda: selection)
    sampler.slab_top, sampler.slab_domains = 100e-9, [1, 4, 5]
    for z in [0., 50e-9, 100e-9, 900e-9]:
        xyz = sampler.plane(1, np.array([[1e-7, 2e-7]]), z, ['z'])
        assert xyz[0, 2] == z
    assert selected == [[1, 4, 5], [1, 4, 5], [1, 4, 5], 'all']


def test_s4_reprocess_preflight_failure_preserves_existing_case(tmp_path, monkeypatch):
    from scripts.analysis.reprocess_s4_boundary import reprocess
    from scripts.run_main import run_s4_boundary as entry
    monkeypatch.setattr(entry, 'TASK_ROOT', tmp_path / 'tasks')
    case = tmp_path / 'run/case0'
    snapshot = case / '99_config/parameter.json'
    entry.write_json(snapshot, {})
    entry.write_json(case / '99_config/realized_geometry.json', {})
    entry.write_json(case / '99_config/s4_summary.json', {})
    entry.write_json(case.parent / '99_config/s4_manifest.json', {
        'cases': [{'case_id': 'case0', 'snapshot': 'case0/99_config/parameter.json',
                   'snapshot_sha256': entry.digest(snapshot)}], 's4': {}})
    (case / '00_model').mkdir()
    (case / '00_model/s4_gamma.mph').write_bytes(b'protected saved model')
    before = {p: entry.digest(p) for p in case.parent.rglob('*') if p.is_file()}
    def unavailable():
        raise RuntimeError('backend unavailable')
    monkeypatch.setattr(entry, 'runtime_preflight', unavailable)
    with pytest.raises(RuntimeError, match='backend unavailable'):
        reprocess(case, replace_results=True)
    assert before == {p: entry.digest(p) for p in case.parent.rglob('*') if p.is_file()}
    audit_path, = entry.TASK_ROOT.glob('s4_postprocess_*/audit.json')
    audit = entry.json_read(audit_path)
    assert audit['status'] == 'failed' and audit['source_mph_unchanged']
    assert audit['eigensolves'] == 0 and 'mph' not in sys.modules
