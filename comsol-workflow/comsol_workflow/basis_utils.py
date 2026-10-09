import math
import numpy as np

def get_fourier_basis_matrix(N):
    """
    Constructs the symmetry-adapted orthonormal basis matrix U for dimension N.
    
    The columns are ordered according to the convention:
      - Column 0: DC vector
      - Columns 1..(N_odd)/2 pairs: c_k, s_k
      - Final Column (only if N even): Nyquist vector v_{N/2}
      
    Returns:
        U (np.ndarray): N x N orthogonal matrix.
    """
    U = np.zeros((N, N))
    
    # 1. Constant (DC) vector
    # v0 = 1/sqrt(N) * [1, 1, ..., 1]
    U[:, 0] = 1.0 / np.sqrt(N)
    
    # 2. Cosine/Sine pairs
    # k ranges from 1 to floor((N-1)/2)
    max_k = (N - 1) // 2
    col_idx = 1
    
    # Use standard 0-indexed grid: theta_n = 2*pi*n/N
    n_grid = np.arange(N)
    
    for k in range(1, max_k + 1):
        theta = (2.0 * np.pi * k * n_grid) / N
        
        # c_k
        c_k = np.sqrt(2.0 / N) * np.cos(theta)
        U[:, col_idx] = c_k
        col_idx += 1
        
        # s_k
        s_k = np.sqrt(2.0 / N) * np.sin(theta)
        U[:, col_idx] = s_k
        col_idx += 1
        
    # 3. Alternating (Nyquist) vector: only if N is even
    if N % 2 == 0:
        # v_{N/2} = 1/sqrt(N) * [1, -1, 1, -1, ...]
        v_nyq = (1.0 / np.sqrt(N)) * ((-1.0)**n_grid)
        U[:, col_idx] = v_nyq
        
    return U

def get_fourier_basis_min_nonzero_values(N):
    """
    Computes the analytical minimum non-zero absolute value for each column 
    of the orthonormal Fourier basis matrix U.
    
    Returns:
        mins (np.ndarray): 1D array of length N.
    """
    mins = np.zeros(N)
    
    # 1. DC Vector [1, 1, ..., 1] * 1/sqrt(N)
    mins[0] = 1.0 / np.sqrt(N)
    
    max_k = (N - 1) // 2
    idx = 1
    
    for k in range(1, max_k + 1):
        # Calculate properties of the grid for frequency k
        g = math.gcd(k, N)
        N_prime = N // g
        
        # Analytic min for Cosine vector
        if N_prime % 4 == 0:
            val_c = np.sin(2.0 * np.pi / N_prime)
        elif N_prime % 2 == 0:
            val_c = np.sin(np.pi / N_prime)
        else:
            val_c = np.sin(np.pi / (2.0 * N_prime))
            
        mins[idx] = val_c * np.sqrt(2.0 / N)
        idx += 1
        
        # Analytic min for Sine vector
        if N_prime % 2 == 0:
            val_s = np.sin(2.0 * np.pi / N_prime)
        else:
            val_s = np.sin(np.pi / N_prime)
            
        mins[idx] = val_s * np.sqrt(2.0 / N)
        idx += 1
        
    # 3. Nyquist Vector (if N is even) [1, -1, ...] * 1/sqrt(N)
    if N % 2 == 0:
        mins[idx] = 1.0 / np.sqrt(N)
        
    return mins

def standard_to_fourier_basis(x, N):
    """
    Transforms vector x (standard basis) to alpha (Fourier basis).
    alpha = U.T * x
    """
    U = get_fourier_basis_matrix(N)
    return U.T @ x

def fourier_to_standard_basis(alpha, N):
    """
    Transforms coefficients alpha (Fourier basis) to x (standard basis).
    x = U * alpha
    """
    U = get_fourier_basis_matrix(N)
    return U @ alpha
