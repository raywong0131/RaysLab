import sys
import numpy as np
import pytest

from scripts.analysis.decompose_s4_c6v import symmetry_coefficients, weighted_parity, scaled_terms
from comsol_workflow.s4_boundary import maxwell_terms


def basis(xy):
    x, y = xy.T
    return np.column_stack([x*np.exp(-x*x-y*y), y*np.exp(-x*x-y*y)])


def test_constant_complex_subspace_recovers_target_on_new_points():
    rng = np.random.default_rng(1917)
    xy = rng.uniform(-1, 1, (251, 2))
    mix = np.array([[1+.3j, -.8+.1j], [.13-.2j, 1+.5j]]) @ np.diag([1e-4, 1e3])
    h, hx, hy = [basis(xy*r) @ mix for r in ([1, 1], [-1, 1], [1, -1])]
    coeff = symmetry_coefficients(h, hx, hy, np.ones(len(xy)))
    new = rng.uniform(-1, 1, (617, 2))
    fields = [basis(new*r) @ mix @ coeff for r in ([1, 1], [-1, 1], [1, -1])]
    assert max(weighted_parity(*fields, np.ones(len(new))).values()) < 1e-12
    assert abs((mix @ coeff)[1]) < 1e-12
    assert "mph" not in sys.modules


def test_rank_deficient_subspace_rejected():
    x = np.arange(1, 11)
    h = np.column_stack([x, 2j*x])
    with pytest.raises(ValueError, match="Rank-deficient"):
        symmetry_coefficients(h, -h, h, np.ones(len(x)))


def test_validation_does_not_force_parity_to_zero():
    rng = np.random.default_rng(52)
    xy = rng.uniform(-1, 1, (401, 2))
    h = basis(xy)
    coeff = symmetry_coefficients(h, basis(xy*[-1, 1]), basis(xy*[1, -1]), np.ones(len(xy)))
    combined = h @ coeff
    errors = weighted_parity(combined, -combined+0.1, combined, np.ones(len(xy)))
    assert errors["odd_x_error"] > 0.05


def test_maxwell_sum_keeps_individual_complex_frequencies():
    coefficients = np.array([.8+.2j, -.17+.07j])
    frequencies = np.array([1e15-1e12j, 1.001e15-1.1e12j])
    data = [np.array([1+.1j, -.2+.3j])*scale for scale in (1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6)]
    actual = scaled_terms(*data, frequencies, 6e-13, 1, 10.89, coefficients)
    parts = [maxwell_terms(*(d[i] for d in data), frequencies[i], 6e-13, 1, 10.89) for i in range(2)]
    for key in actual:
        assert np.allclose(actual[key], sum(coefficients[i]*parts[i][key] for i in range(2)))
    common = maxwell_terms(*(np.dot(coefficients, d) for d in data), frequencies.mean(), 6e-13, 1, 10.89)
    assert abs(common["Ly"]-actual["Ly"]) > 1e-6


def test_gradient_average_preserves_sign_units_and_each_frequency():
    from scripts.analysis.plot_s4_control import gradient_average
    from comsol_workflow.s4_boundary import EPS0

    omega = np.array([1e15-1e12j, 1.2e15-2e12j])
    coefficients = np.array([.8+.2j, -.17+.07j])
    electric = np.array([3-4j, -2+1j])
    area = 6e-13
    integrals = 1j*omega*EPS0*area*electric
    actual = gradient_average(integrals, coefficients, omega, area)
    assert np.isclose(actual, coefficients @ electric, rtol=1e-12)
    assert not np.isclose(actual, gradient_average(integrals, coefficients, np.full(2, omega.mean()), area))
    with pytest.raises(ValueError, match="Invalid"):
        gradient_average(integrals, coefficients, [0, omega[1]], area)
    assert "mph" not in sys.modules


def test_dimensional_panels_do_not_normalize_or_rotate_fields():
    import pandas as pd
    from comsol_workflow.s4_plotting import build_edge_figures
    from comsol_workflow.s4_boundary import hole_edges, polygon_area

    outer = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * 1e-6
    triangle = np.array([[-.05, -.05], [.05, -.05], [0, .05]])
    triangle -= triangle.mean(axis=0)
    holes = [(triangle + .6*np.array([np.cos(angle), np.sin(angle)]))*1e-6
             for angle in np.arange(6)*np.pi/3]
    rows = [{"hole": e.hole, "edge": e.edge, "group": e.group,
             "x0_m": e.start[0], "y0_m": e.start[1], "x1_m": e.end[0], "y1_m": e.end[1],
             "nx_dielectric_to_air": e.normal[0], "ny_dielectric_to_air": e.normal[1]}
            for e in hole_edges(holes)]
    frame = pd.DataFrame(rows)
    field = {"grid_electric_vm": np.full((3, 3, 2), 3+4j), "inside": np.ones((3, 3), bool),
             "grid_x_m": np.linspace(-1e-6, 1e-6, 3), "grid_y_m": np.linspace(-1e-6, 1e-6, 3), "outer_m": outer}
    boundary = {"x": np.linspace(-2, 3, 18)+2j, "y": np.ones(18)*(1-1j)}
    figures = build_edge_figures(frame, frame.assign(z_over_H=0.),
        pd.DataFrame({"area_m2": [abs(polygon_area(outer))]}), {"frequency_thz": 197.},
        field, boundary, 2-3j, dimensional=True)
    try:
        fig = figures["s4_edges_py"]
        area_ax, _, bars = fig.axes
        assert np.allclose(area_ax.collections[0].get_array(), 3.)
        assert np.allclose([p.get_height() for p in bars.patches], boundary["x"].real)
        assert 'V/m' in area_ax.child_axes[0].get_ylabel()
        assert 'V/m' in bars.get_ylabel()
        assert fig.texts == [fig._suptitle]
    finally:
        for fig in figures.values():
            fig.clear()
