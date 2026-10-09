import sys
import numpy as np

from scripts.run_main.run_s4_mesh_convergence import subspace_scores, complex_change, snapshot_for_mesh
from scripts.run_main.run_s4_mesh_convergence import explicit_frequencies
from types import SimpleNamespace


def test_projection_is_invariant_to_complex_pair_rotation():
    rng = np.random.default_rng(42)
    basis = rng.normal(size=(100, 2))+1j*rng.normal(size=(100, 2))
    weights = rng.uniform(.1, 1., 100)
    rotation = np.array([[1, 1j], [1j, 2]])
    np.testing.assert_allclose(subspace_scores(basis, basis@rotation, weights), [1., 1.], atol=1e-12)
    test = rng.normal(size=(100, 4))+1j*rng.normal(size=(100, 4))
    np.testing.assert_allclose(subspace_scores(basis, test, weights),
                               subspace_scores(basis@rotation, test, weights), atol=1e-12)


def test_orthogonal_mode_has_zero_projection():
    basis = np.eye(4, dtype=complex)[:, :2]
    candidates = np.eye(4, dtype=complex)
    np.testing.assert_allclose(subspace_scores(basis, candidates, np.ones(4)), [1., 1., 0., 0.])


def test_change_uses_full_complex_difference():
    absolute, relative = complex_change(1., 1j)
    np.testing.assert_allclose([absolute, relative], [np.sqrt(2), np.sqrt(2)])
    assert complex_change(0j, 0j) == (0., 0.)


def test_snapshot_changes_only_mesh_and_preserves_input():
    source = {"mesh_size": 5, "unit_cell_2d": {"b0_nm": 245, "eta": .95, "zeta": 1.}}
    target = snapshot_for_mesh(source, 3)
    assert target["mesh_size"] == 3 and source["mesh_size"] == 5
    assert target["unit_cell_2d"] == source["unit_cell_2d"]
    assert target["unit_cell_2d"] is not source["unit_cell_2d"]


def test_import_has_no_comsol_side_effect():
    assert "mph" not in sys.modules


def test_frequency_reader_uses_explicit_dimensionless_units():
    class Node:
        def __init__(self):
            self.settings = {}
        def set(self, key, val):
            self.settings[key] = val
        def getReal(self):
            return [[200., 201.], [.1, .2], [1000., 502.5]]
    node = Node()
    numerical = SimpleNamespace(create=lambda *args: node)
    model = SimpleNamespace(java=SimpleNamespace(result=lambda: SimpleNamespace(numerical=lambda: numerical)))
    f = explicit_frequencies(model)
    np.testing.assert_allclose(f[:, 1], [.1, .2])
    assert node.settings["unit"] == ["1", "1", "1"]
    assert all("1[THz]" in e for e in node.settings["expr"][:2])


def test_frequency_readers_keep_negative_damping_and_check_q_magnitude():
    import pytest
    from scripts.run_main.run_s4_surface_zeta import frequencies
    values=np.array([[197.,200.],[.1,-2.7e-11],[985.,200./(2*2.7e-11)]])
    class Node:
        def set(self,*args):pass
        def getReal(self):return values
    node=Node()
    model=SimpleNamespace(java=SimpleNamespace(result=lambda:SimpleNamespace(numerical=lambda:SimpleNamespace(create=lambda *args:node))))
    for read in (explicit_frequencies,lambda m:frequencies(m).to_numpy()):
        np.testing.assert_array_equal(read(model),values.T)
    values[2,1]*=2
    with pytest.raises(ValueError,match="inconsistent"):explicit_frequencies(model)
    with pytest.raises(AssertionError):frequencies(model)
