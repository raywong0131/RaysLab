import os
import numpy as np
import sys

from comsol_workflow.energy_recovery import (
    compute_autocorrelation,
    fourier_subspace_energies,
    fourier_subspace_energies_direct
)


def test_autocorrelation_properties():
    """Test autocorrelation properties for both real and complex vectors"""
    N = 8
    
    # Real case
    x_real = np.random.randn(N)
    r_real = compute_autocorrelation(x_real)
    
    # Check r_0 = ||x||^2
    assert np.isclose(r_real[0], np.dot(x_real, x_real)), \
        f"Real: r_0={r_real[0]:.4f} != ||x||^2={np.dot(x_real, x_real):.4f}"
    
    # Check symmetry: r_m = r_{N-m}
    for m in range(1, N // 2 + 1):
        assert np.isclose(r_real[m], r_real[N - m]), \
            f"Real: Symmetry broken at m={m}"
    
    # Complex case
    x_complex = np.random.randn(N) + 1j * np.random.randn(N)
    r_complex = compute_autocorrelation(x_complex)
    
    # Check r_0 = ||x||^2 and is real
    x_norm_sq = np.vdot(x_complex, x_complex).real
    assert np.isclose(r_complex[0].real, x_norm_sq) and np.abs(r_complex[0].imag) < 1e-10, \
        f"Complex: r_0={r_complex[0]} != ||x||^2={x_norm_sq:.4f}"
    
    # Check Hermitian symmetry: r_{N-m} = conj(r_m)
    for m in range(1, N // 2 + 1):
        assert np.isclose(r_complex[N - m], np.conj(r_complex[m])), \
            f"Complex: Hermitian symmetry broken at m={m}"


def test_dft_properties():
    """Test DFT of autocorrelation is real and non-negative"""
    N = 8
    
    # Real case
    x_real = np.random.randn(N)
    energies_real = fourier_subspace_energies(x_real)
    R_real = energies_real['R']
    
    assert np.all(np.abs(np.imag(R_real)) < 1e-10), \
        f"Real: DFT has imaginary components: max={np.max(np.abs(np.imag(R_real))):.2e}"
    assert np.all(R_real >= -1e-10), \
        f"Real: DFT has negative values: min={np.min(R_real):.4f}"
    
    # Complex case
    x_complex = np.random.randn(N) + 1j * np.random.randn(N)
    energies_complex = fourier_subspace_energies(x_complex)
    R_complex = energies_complex['R']
    
    assert np.all(np.abs(np.imag(R_complex)) < 1e-10), \
        f"Complex: DFT has imaginary components: max={np.max(np.abs(np.imag(R_complex))):.2e}"
    assert np.all(R_complex >= -1e-10), \
        f"Complex: DFT has negative values: min={np.min(R_complex):.4f}"


def test_documentation_example_real():
    """Test against documentation example (N=6, real)"""
    x_doc = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    r_doc = compute_autocorrelation(x_doc)
    
    r_expected = np.array([91, 76, 67, 64, 67, 76], dtype=float)
    assert np.allclose(r_doc, r_expected), \
        f"Autocorrelation mismatch: expected {r_expected}, got {r_doc}"
    
    energies_doc = fourier_subspace_energies(x_doc)
    R_doc = energies_doc['R']
    R_expected = np.array([441, 36, 12, 9, 12, 36], dtype=float)
    assert np.allclose(R_doc, R_expected, atol=1e-10), \
        f"DFT mismatch: expected {R_expected}, got {R_doc}"
    
    assert np.isclose(energies_doc['E_dc'], 73.5), "E_dc mismatch"
    assert np.isclose(energies_doc['E_pairs'][0], 12), "E_1 mismatch"
    assert np.isclose(energies_doc['E_pairs'][1], 4), "E_2 mismatch"
    assert np.isclose(energies_doc['E_nyquist'], 1.5), "E_nyquist mismatch"


def test_documentation_example_complex():
    """Test against documentation example (N=6, complex)"""
    x_doc = np.array([1, 1j, 0, 0, 0, 0], dtype=complex)
    r_doc = compute_autocorrelation(x_doc)
    
    r_expected = np.array([2, -1j, 0, 0, 0, 1j], dtype=complex)
    assert np.allclose(r_doc, r_expected), \
        f"Autocorrelation mismatch: expected {r_expected}, got {r_doc}"
    
    energies_doc = fourier_subspace_energies(x_doc)
    R_doc = energies_doc['R']
    sqrt3 = np.sqrt(3)
    R_expected = np.array([2, 2-sqrt3, 2-sqrt3, 2, 2+sqrt3, 2+sqrt3])
    assert np.allclose(R_doc, R_expected, atol=1e-10), \
        f"DFT mismatch: expected {R_expected}, got {R_doc}"
    
    assert np.isclose(energies_doc['E_dc'], 1/3), "E_dc mismatch"
    assert np.isclose(energies_doc['E_pairs'][0], 2/3), "E_1 mismatch"
    assert np.isclose(energies_doc['E_pairs'][1], 2/3), "E_2 mismatch"
    assert np.isclose(energies_doc['E_nyquist'], 1/3), "E_nyquist mismatch"


def test_energy_conservation():
    """Test Parseval's theorem (energy conservation)"""
    for N in [5, 6, 8, 12, 15]:
        # Real case
        x_real = np.random.randn(N)
        energies_real = fourier_subspace_energies(x_real)
        true_total_real = np.dot(x_real, x_real)
        assert np.isclose(energies_real['total_energy'], true_total_real, rtol=1e-10), \
            f"Real N={N}: Energy not conserved"
        
        # Complex case
        x_complex = np.random.randn(N) + 1j * np.random.randn(N)
        energies_complex = fourier_subspace_energies(x_complex)
        true_total_complex = np.vdot(x_complex, x_complex).real
        assert np.isclose(energies_complex['total_energy'], true_total_complex, rtol=1e-10), \
            f"Complex N={N}: Energy not conserved"


def test_comparison_with_direct_method():
    """Compare cyclic-shift method with direct Fourier decomposition"""
    for N in [5, 6, 8, 12]:
        # Real case
        x_real = np.random.randn(N)
        energies_shift = fourier_subspace_energies(x_real)
        energies_direct = fourier_subspace_energies_direct(x_real)
        
        assert np.isclose(energies_shift['E_dc'], energies_direct['E_dc']), \
            f"Real N={N}: DC mismatch"
        assert np.allclose(energies_shift['E_pairs'], energies_direct['E_pairs']), \
            f"Real N={N}: Pairs mismatch"
        if N % 2 == 0:
            assert np.isclose(energies_shift['E_nyquist'], energies_direct['E_nyquist']), \
                f"Real N={N}: Nyquist mismatch"
        
        # Complex case
        x_complex = np.random.randn(N) + 1j * np.random.randn(N)
        energies_shift_c = fourier_subspace_energies(x_complex)
        energies_direct_c = fourier_subspace_energies_direct(x_complex)
        
        assert np.isclose(energies_shift_c['E_dc'], energies_direct_c['E_dc']), \
            f"Complex N={N}: DC mismatch"
        assert np.allclose(energies_shift_c['E_pairs'], energies_direct_c['E_pairs']), \
            f"Complex N={N}: Pairs mismatch"
        if N % 2 == 0:
            assert np.isclose(energies_shift_c['E_nyquist'], energies_direct_c['E_nyquist']), \
                f"Complex N={N}: Nyquist mismatch"


def test_invariance_global_phase():
    """Test invariance under global sign (real) or phase (complex) transformations"""
    N = 10
    
    # Real: sign flip
    x_real = np.random.randn(N)
    x_neg = -x_real
    energies_pos = fourier_subspace_energies(x_real)
    energies_neg = fourier_subspace_energies(x_neg)
    
    assert np.isclose(energies_pos['E_dc'], energies_neg['E_dc']), \
        "Real: DC not invariant under sign flip"
    assert np.allclose(energies_pos['E_pairs'], energies_neg['E_pairs']), \
        "Real: Pairs not invariant under sign flip"
    if N % 2 == 0:
        assert np.isclose(energies_pos['E_nyquist'], energies_neg['E_nyquist']), \
            "Real: Nyquist not invariant under sign flip"
    
    # Complex: global phase
    x_complex = np.random.randn(N) + 1j * np.random.randn(N)
    phi = np.random.rand() * 2 * np.pi
    x_phase = np.exp(1j * phi) * x_complex
    
    energies_orig = fourier_subspace_energies(x_complex)
    energies_phase = fourier_subspace_energies(x_phase)
    
    assert np.isclose(energies_orig['E_dc'], energies_phase['E_dc']), \
        "Complex: DC not invariant under phase"
    assert np.allclose(energies_orig['E_pairs'], energies_phase['E_pairs']), \
        "Complex: Pairs not invariant under phase"
    if N % 2 == 0:
        assert np.isclose(energies_orig['E_nyquist'], energies_phase['E_nyquist']), \
            "Complex: Nyquist not invariant under phase"


def test_invariance_cyclic_shift():
    """Test invariance under cyclic shifts"""
    N = 10
    
    # Real case
    x_real = np.random.randn(N)
    x_shifted = np.roll(x_real, 3)
    energies_orig = fourier_subspace_energies(x_real)
    energies_shift = fourier_subspace_energies(x_shifted)
    
    assert np.isclose(energies_orig['E_dc'], energies_shift['E_dc']), \
        "Real: DC changes under shift"
    assert np.allclose(energies_orig['E_pairs'], energies_shift['E_pairs']), \
        "Real: Pairs change under shift"
    if N % 2 == 0:
        assert np.isclose(energies_orig['E_nyquist'], energies_shift['E_nyquist']), \
            "Real: Nyquist changes under shift"
    
    # Complex case
    x_complex = np.random.randn(N) + 1j * np.random.randn(N)
    x_shifted_c = np.roll(x_complex, 3)
    energies_orig_c = fourier_subspace_energies(x_complex)
    energies_shift_c = fourier_subspace_energies(x_shifted_c)
    
    assert np.isclose(energies_orig_c['E_dc'], energies_shift_c['E_dc']), \
        "Complex: DC changes under shift"
    assert np.allclose(energies_orig_c['E_pairs'], energies_shift_c['E_pairs']), \
        "Complex: Pairs change under shift"
    if N % 2 == 0:
        assert np.isclose(energies_orig_c['E_nyquist'], energies_shift_c['E_nyquist']), \
            "Complex: Nyquist changes under shift"


def test_edge_cases():
    """Test edge cases: constant, alternating, and pure imaginary vectors"""
    # Constant vector (all energy in DC)
    x_const = np.ones(8)
    energies_const = fourier_subspace_energies(x_const)
    
    assert np.isclose(energies_const['E_dc'], 8), "Constant: DC energy incorrect"
    assert np.allclose(energies_const['E_pairs'], 0, atol=1e-10), \
        "Constant: Non-DC energy should be zero"
    assert np.isclose(energies_const['E_nyquist'], 0, atol=1e-10), \
        "Constant: Nyquist energy should be zero"
    
    # Alternating vector (all energy in Nyquist)
    x_alt = np.array([1, -1, 1, -1, 1, -1, 1, -1], dtype=float)
    energies_alt = fourier_subspace_energies(x_alt)
    
    assert np.isclose(energies_alt['E_dc'], 0, atol=1e-10), \
        "Alternating: DC energy should be zero"
    assert np.allclose(energies_alt['E_pairs'], 0, atol=1e-10), \
        "Alternating: Pair energy should be zero"
    assert np.isclose(energies_alt['E_nyquist'], 8), \
        "Alternating: Nyquist energy incorrect"
    
    # Complex: pure imaginary
    x_imag = 1j * np.random.randn(8)
    energies_imag = fourier_subspace_energies(x_imag)
    true_energy = np.vdot(x_imag, x_imag).real
    
    assert np.isclose(energies_imag['total_energy'], true_energy), \
        "Pure imaginary: energy not conserved"


if __name__ == "__main__":
    print("Running tests...")
    
    test_autocorrelation_properties()
    print("✓ test_autocorrelation_properties")
    
    test_dft_properties()
    print("✓ test_dft_properties")
    
    test_documentation_example_real()
    print("✓ test_documentation_example_real")
    
    test_documentation_example_complex()
    print("✓ test_documentation_example_complex")
    
    test_energy_conservation()
    print("✓ test_energy_conservation")
    
    test_comparison_with_direct_method()
    print("✓ test_comparison_with_direct_method")
    
    test_invariance_global_phase()
    print("✓ test_invariance_global_phase")
    
    test_invariance_cyclic_shift()
    print("✓ test_invariance_cyclic_shift")
    
    test_edge_cases()
    print("✓ test_edge_cases")
    
    print("\n=== All Tests Passed ===")
