import numpy as np
import pytest
from comsol_workflow.s4_boundary import signed_indicator


def test_signed_projection_invariants_and_zero():
    expected = signed_indicator(2 + 3j, -4 - 5j)
    assert expected['s'] < 0 and expected['inv_s'] < 0
    for c in (1e-100, 1e100, np.exp(1.7j)):
        result = signed_indicator(c * (2 + 3j), c * (-4 - 5j))
        np.testing.assert_allclose(list(result.values()), list(expected.values()))
        assert np.isclose(result['s']**2 + result['q_perp']**2, result['r']**2)
    cancelled = signed_indicator(1, -1)
    assert cancelled['s'] == cancelled['r'] == 0 and np.isnan(cancelled['inv_s'])
    perpendicular = signed_indicator(1, -1 + 1j)
    assert perpendicular['s'] == 0 and perpendicular['r'] > 0
    for a, b in ((0, 1), (np.nan, 1), (1, np.inf)):
        with pytest.raises(ValueError):
            signed_indicator(a, b)
