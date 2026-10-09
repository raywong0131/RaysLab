import numpy as np
import pytest
from scipy.constants import epsilon_0, mu_0

from scripts.analysis.plot_dipolar_angle_filter import cone_fields, energy, Z0


def test_polar_cone_preserves_complex_plane_waves_and_maxwell_relation():
    q = np.array([-1., 0., 1.])
    k0 = 4.
    electric = np.zeros((3, 3, 3), complex)
    electric[1, 1, 1] = 1+2j
    electric[1, 2, 1] = -3+.5j  # Off-axis, transverse y polarization.
    e, h, mask = cone_fields(q, electric, k0, 5.)
    assert mask.sum() == 1
    np.testing.assert_array_equal(e[1, 1], electric[1, 1])
    np.testing.assert_allclose(h[1, 1], [-(1+2j)/Z0, 0, 0])
    assert not e[1, 2].any()
    e, h, mask = cone_fields(q, electric, k0, 15.)
    assert mask[1, 2] and not mask[2, 2]  # Circular cone, not a square crop.
    np.testing.assert_array_equal(e, electric)
    np.testing.assert_allclose(epsilon_0*np.sum(abs(e)**2, -1), mu_0*np.sum(abs(h)**2, -1), atol=1e-25)
    np.testing.assert_allclose(energy(e, h), epsilon_0/2*np.sum(abs(e)**2, -1), atol=1e-25)
    for theta in (0., -1., 91., np.nan):
        with pytest.raises(ValueError):
            cone_fields(q, electric, k0, theta)
