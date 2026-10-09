import os
import numpy as np

from comsol_workflow.basis_utils import (
    get_fourier_basis_matrix, 
    get_fourier_basis_min_nonzero_values,
    standard_to_fourier_basis, 
    fourier_to_standard_basis
)

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)

    print("=== TEST 1: N=6 Specific Vector Checks (from Doc Section 1.4) ===")
    N_test = 6
    U_6 = get_fourier_basis_matrix(N_test)

    print("Fourier Basis Matrix U (N=6):")
    print(U_6)
    
    # Check dimensions
    assert U_6.shape == (6, 6), "Shape mismatch"
    
    # 1. Check DC (Column 0)
    # Doc: 1/sqrt(6) * [1, 1, 1, 1, 1, 1]
    v0_expected = np.ones(6) / np.sqrt(6)
    if np.allclose(U_6[:, 0], v0_expected):
        print("[Pass] DC vector matches")
    else:
        print("[Fail] DC vector mismatch")
        
    # 2. Check c_1 (Column 1)
    # Doc: 1/(2*sqrt(3)) * [2, 1, -1, -2, -1, 1]
    factor_c1 = 1.0 / (2.0 * np.sqrt(3))
    c1_expected = factor_c1 * np.array([2, 1, -1, -2, -1, 1])
    if np.allclose(U_6[:, 1], c1_expected):
        print("[Pass] c_1 vector matches")
    else:
        print("[Fail] c_1 vector mismatch")
        print("Calculated:", U_6[:, 1])
        print("Expected:  ", c1_expected)

    # 3. Check s_1 (Column 2)
    # Doc: 1/2 * [0, 1, 1, 0, -1, -1]
    s1_expected = 0.5 * np.array([0, 1, 1, 0, -1, -1])
    if np.allclose(U_6[:, 2], s1_expected):
        print("[Pass] s_1 vector matches")
    else:
        print("[Fail] s_1 vector mismatch")
        
    # 4. Check Nyquist (Column 5, last)
    # Doc: 1/sqrt(6) * [1, -1, 1, -1, 1, -1]
    nyq_expected = (1.0 / np.sqrt(6)) * np.array([1, -1, 1, -1, 1, -1])
    if np.allclose(U_6[:, -1], nyq_expected):
        print("[Pass] Nyquist vector matches")
    else:
        print("[Fail] Nyquist vector mismatch")

    print("\n=== TEST 2: Orthogonality (U^T U = I) ===")
    # Test for odd N as well
    for n_dim in [5, 6, 12]:
        U = get_fourier_basis_matrix(n_dim)
        gram = U.T @ U
        identity = np.eye(n_dim)
        err = np.linalg.norm(gram - identity)
        if err < 1e-14:
            print(f"[Pass] N={n_dim} is orthonormal (Error: {err:.2e})")
        else:
            print(f"[Fail] N={n_dim} matrix is not orthonormal. Error: {err}")

    print("\n=== TEST 3: Reconstruction (Forward -> Inverse) ===")
    x_original = np.random.randn(N_test)
    alpha = standard_to_fourier_basis(x_original, N_test)
    x_recon = fourier_to_standard_basis(alpha, N_test)
    
    recon_err = np.linalg.norm(x_original - x_recon)
    print(f"Original vector: {x_original[:3]}...")
    print(f"Reconstructed:   {x_recon[:3]}...")
    if recon_err < 1e-14:
        print(f"[Pass] Perfect reconstruction (Error: {recon_err:.2e})")
    else:
        print(f"[Fail] Reconstruction error high: {recon_err}")

    print("\n=== TEST 4: Parseval's Theorem (Norm Preservation) ===")
    norm_x = np.linalg.norm(x_original)
    norm_alpha = np.linalg.norm(alpha)
    if np.isclose(norm_x, norm_alpha):
        print(f"[Pass] Norm preserved: {norm_x:.4f} == {norm_alpha:.4f}")
    else:
        print(f"[Fail] Norms differ: {norm_x} vs {norm_alpha}")

    print("\n=== TEST 5: Analytical Min Non-Zero Values ===")
    test_5_passed = True
    # Test a mix of dimensions: Even, Odd, Multiples of 4
    n_values = [4, 5, 6, 7, 8, 12, 15, 16]
    
    for n in n_values:
        U = get_fourier_basis_matrix(n)
        mins_analytic = get_fourier_basis_min_nonzero_values(n)
        mins_numeric = np.zeros(n)
        
        for i in range(n):
            # Extract column vector
            vec = U[:, i]
            abs_vec = np.abs(vec)
            # Find min non-zero value (filtering out floating point near-zeros)
            non_zeros = abs_vec[abs_vec > 1e-10]
            
            if non_zeros.size > 0:
                mins_numeric[i] = np.min(non_zeros)
            else:
                mins_numeric[i] = 0.0 # Should not occur for valid basis vectors
        
        if not np.allclose(mins_analytic, mins_numeric, atol=1e-10):
            print(f"[Fail] N={n} mismatch")
            print(f"  Analytic: {mins_analytic}")
            print(f"  Numeric:  {mins_numeric}")
            test_5_passed = False
            
    if test_5_passed:
        print(f"[Pass] Analytical values match numerical search for N={n_values}")
