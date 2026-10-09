import numpy as np

from scripts.analysis.plot_s5_gamma_mechanism import quadrant_integrals


def test_quadrants_keep_complex_phase_and_share_axis_samples():
    axis = np.array([-1., 0., 1.])
    x, y = np.meshgrid(axis, axis)
    odd, even = 2+3j, 1-2j
    air = {'x_axis': axis, 'y_axis': axis, 'Ex': odd*x*y,
           'Ey': np.full((3, 3), even)}
    parts, total = quadrant_integrals(air)
    np.testing.assert_allclose(parts[:, 0], odd*np.array([1, -1, 1, -1])*1e-12)
    np.testing.assert_allclose(parts[:, 1], np.full(4, 2.25*even)*1e-12)
    np.testing.assert_allclose(total, [0, 9*even*1e-12])
