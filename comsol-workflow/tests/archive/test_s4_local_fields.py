"""Pure numerical tests: no COMSOL model, solver or official results access."""
import sys

import numpy as np

from scripts.analysis.validate_s4_local_fields import (
    central_difference, coherent_pair, contains, convert_derivatives,
    distance, fd_mask, material_ids, metrics,
)
from comsol_workflow.s4_boundary import EPS0


def geometry():
    outer = np.array([[-5., -5.], [5., -5.], [5., 5.], [-5., 5.]])*1e-9
    hole = np.array([[-1., -1.], [1., -1.], [0., 1.]])*1e-9
    return outer, [hole]


def test_convex_membership_and_distance():
    outer, holes = geometry()
    xy = np.array([[0., 0.], [2., 0.], [6., 0.]])*1e-9
    np.testing.assert_array_equal(contains(xy, outer), [True, True, False])
    np.testing.assert_array_equal(contains(xy, outer[::-1]), contains(xy, outer))
    np.testing.assert_array_equal(material_ids(xy, holes), [1, 0, 0])
    np.testing.assert_allclose(distance(xy[1:2], [outer]), [3e-9])


def test_fd_mask_keeps_material_and_cell():
    outer, holes = geometry()
    xy = np.array([[0., 0.], [2., 0.], [4.95, 0.], [.49, 0.]])*1e-9
    np.testing.assert_array_equal(fd_mask(xy, outer, holes, .1e-9), [True, True, False, False])


def test_central_difference():
    x, h = np.array([1., 2., 3.]), .01
    np.testing.assert_allclose(central_difference(2*(x+h)+3, 2*(x-h)+3, h), 2)
    np.testing.assert_allclose(central_difference((x+h)**3, (x-h)**3, h), 3*x*x+h*h)


def test_each_mode_keeps_own_frequency():
    gradients = np.array([[[2+1j, 3-2j]], [[4-1j, 5+3j]]])
    omegas, eps = np.array([2-0.1j, 3-0.2j]), np.array([1., 10.89])
    coefficients = np.array([1+.2j, -.3+.1j])
    converted = convert_derivatives(gradients, omegas, eps)
    expected = sum(-1j*coefficients[m]*gradients[m]/(omegas[m]*EPS0*eps) for m in range(2))
    np.testing.assert_allclose(coherent_pair(converted, coefficients), expected)


def test_region_means_add_without_renormalization():
    ref = np.array([2+1j, 3-2j, -4+1j])
    pred, weights = 2*ref, np.array([1., 2., 3.])
    masks = [np.array([True, False, False]), np.array([False, True, True])]
    whole = metrics(ref, pred, weights, np.ones(3, bool), 6.)
    parts = [metrics(ref, pred, weights, m, 6.) for m in masks]
    for field in ("Ey_mean_fullcell", "Hz_Ey_mean_fullcell", "residual_mean_fullcell", "area_fraction"):
        np.testing.assert_allclose(sum(v[field] for v in parts), whole[field])
    np.testing.assert_allclose(whole["relative_L2"], 1.)


def test_import_does_not_load_mph():
    assert "mph" not in sys.modules
