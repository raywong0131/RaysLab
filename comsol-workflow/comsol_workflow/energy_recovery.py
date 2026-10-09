import numpy as np
import pandas as pd

from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay

from .basis_utils import standard_to_fourier_basis
from .geometry_utils import rotate_points

def interpolate_field(points, f, query_points):
    lin_interp = LinearNDInterpolator(points, f)
    interpolated_values = lin_interp(query_points)
    
    nan_mask = np.isnan(interpolated_values)
    if np.any(nan_mask):
        nn_interp = NearestNDInterpolator(points, f)
        interpolated_values[nan_mask] = nn_interp(query_points[nan_mask])
    
    return interpolated_values

def integrate_field(points, f):
    points = np.asarray(points)      # (N,2)
    f = np.asarray(f)          # (N,)
    tri = Delaunay(points)
    simplices = tri.simplices         # (M,3) indices

    p0 = points[simplices[:, 0]]
    p1 = points[simplices[:, 1]]
    p2 = points[simplices[:, 2]]

    # Triangle areas in 2D: 0.5 * |(p1-p0) x (p2-p0)|
    v1 = p1 - p0
    v2 = p2 - p0
    areas = 0.5 * np.abs(v1[:, 0]*v2[:, 1] - v1[:, 1]*v2[:, 0]) # (M,)

    f_tri_mean = f[simplices].mean(axis=1)   # (f0+f1+f2)/3 (M,)
    return np.sum(areas * f_tri_mean)

def get_complex_field(df_re, df_im):
    pts_re = df_re.iloc[:, :2].values
    re = df_re.iloc[:, 2].values
    pts_im = df_im.iloc[:, :2].values
    im = df_im.iloc[:, 2].values
    
    x = re + 1j * interpolate_field(pts_im, im, pts_re)
    return pts_re, x

def inner_product_rotated_field(points, x, angle_deg):
    # Evaluate rotated field at the same physical points:
    # x_rot(p) = x(R^{-1} p) = x( rotate(p, -theta) )
    q = rotate_points(points, -angle_deg)
    x_rot = interpolate_field(points, np.real(x), q) + 1j * interpolate_field(points, np.imag(x), q)

    inner_product = integrate_field(points, np.conj(x) * x_rot)
    return inner_product

def compute_autocorrelation_field(df_re, df_im, N):
    pts, x = get_complex_field(df_re, df_im)
    angle_deg = 360.0 / N
    r = [inner_product_rotated_field(pts, x, m * angle_deg) for m in range(N)]
    return np.array(r)

def fourier_subspace_energies_field(df_re, df_im, N):
    """
    Compute Fourier subspace energies from field data using cyclic rotations.
    
    Args:
        df_re: DataFrame with real part of field (x, y, Re[field])
        df_im: DataFrame with imaginary part of field (x, y, Im[field])
        N: Rotational symmetry order
        
    Returns:
        dict: Dictionary with keys:
            - 'E_dc': DC energy (scalar)
            - 'E_pairs': Array of energies for (c_k, s_k) pairs, k=1,...,floor((N-1)/2)
            - 'E_nyquist': Nyquist energy (scalar, only if N is even, else None)
            - 'total_energy': Total energy (should equal ||x||^2)
            - 'R': Raw DFT coefficients (for debugging)
            - 'autocorrelation': Autocorrelation sequence
    """
    # Compute autocorrelation from field rotations
    r = compute_autocorrelation_field(df_re, df_im, N)
    
    # Compute DFT of autocorrelation
    R = np.fft.fft(r).real  # Should be real and non-negative
    
    # Map to real Fourier-subspace energies
    E_dc = R[0] / N
    
    # Cosine/Sine pair energies
    max_k = (N - 1) // 2
    E_pairs = np.zeros(max_k)
    
    for k in range(1, max_k + 1):
        E_pairs[k - 1] = (R[k] + R[N - k]) / N
    
    # Nyquist energy (only if N is even)
    E_nyquist = None
    if N % 2 == 0:
        E_nyquist = R[N // 2] / N
    
    # Compute total energy for verification
    total_energy = E_dc + np.sum(E_pairs)
    if E_nyquist is not None:
        total_energy += E_nyquist
    
    return {
        'E_dc': E_dc,
        'E_pairs': E_pairs,
        'E_nyquist': E_nyquist,
        'total_energy': total_energy,
        'R': R,
    }

def get_standard_basis_field(points, N):
    # Const function -> const in all points
    x = points[:, 0]
    y = points[:, 1]
    theta = np.arctan2(y, x)  # (-pi, pi]
    width = 2 * np.pi / N
    idx = np.floor(((theta + width/2) % (2*np.pi)) / width).astype(int)

    basis = [np.zeros(len(points)) for _ in range(N)]
    for m in range(N):
        basis[m][idx == m] = 1.0
        basis[m][(x==0) & (y==0)] = 1.0 / N
    return basis

def field_to_vector(df_re, df_im, N):
    pts, x = get_complex_field(df_re, df_im)
    standard_basis = get_standard_basis_field(pts, N)
    
    coeffs = []
    for m in range(N):
        # Compute <standard_basis_m, x> to get repr in standard basis
        coeff = integrate_field(pts, np.conj(standard_basis[m]) * x)
        coeffs.append(coeff)
    
    return np.array(coeffs)

def compute_autocorrelation(x):
    """
    Compute cyclic autocorrelation sequence r_m = <x, R^m x> for m=0,...,N-1.
    Uses Hermitian inner product for complex (reduces to dot product for real).
    
    Args:
        x (np.ndarray): Input vector of length N (real or complex)
        
    Returns:
        np.ndarray: Autocorrelation sequence [r_0, r_1, ..., r_{N-1}]
                   (real for real x, complex for complex x)
    """
    N = len(x)
    r = [0] * N
    
    for m in range(N):
        shifted_x = np.roll(x, m)  # R^m x
        r[m] = np.vdot(x, shifted_x)  # Hermitian inner product
    
    return np.array(r)

def fourier_subspace_energies(x):
    """
    Recover energies in each real Fourier subspace from cyclic shift inner products.
    
    This implements the algorithm from inner_product_decomposition*.md:
    1. Measure autocorrelation r_m = <x, R^m x>
    2. Compute DFT to get R_k = |X_k|^2
    3. Map to real Fourier-subspace energies E_k
    
    Works for both real and complex vectors. For:
    - Real x: E_k = 2*R_k / N for pair subspaces
    - Complex x: E_k = (R_k + R_{N-k}) / N for pair subspaces
    
    Args:
        x (np.ndarray): Input vector of length N (real or complex)
        
    Returns:
        dict: Dictionary with keys:
            - 'E_dc': DC energy (scalar)
            - 'E_pairs': Array of energies for (c_k, s_k) pairs, k=1,...,floor((N-1)/2)
            - 'E_nyquist': Nyquist energy (scalar, only if N is even, else None)
            - 'total_energy': Total energy (should equal ||x||^2)
            - 'R': Raw DFT coefficients (for debugging)
    """
    N = len(x)
    
    # Step 1: Compute autocorrelation
    r = compute_autocorrelation(x)
    
    # Step 2: Compute DFT of autocorrelation
    R = np.fft.fft(r).real  # Should be real and non-negative due to properties of autocorrelation
    
    # Step 3: Map to real Fourier-subspace energies
    
    # DC energy: E_0 = R_0 / N
    E_dc = R[0] / N
    
    # Cosine/Sine pair energies
    max_k = (N - 1) // 2
    E_pairs = np.zeros(max_k)
    
    for k in range(1, max_k + 1):
        E_pairs[k - 1] = (R[k] + R[N - k]) / N
    
    # Nyquist energy (only if N is even): E_{N/2} = R_{N/2} / N
    E_nyquist = None
    if N % 2 == 0:
        E_nyquist = R[N // 2] / N
    
    # Compute total energy for verification
    total_energy = E_dc + np.sum(E_pairs)
    if E_nyquist is not None:
        total_energy += E_nyquist
    
    return {
        'E_dc': E_dc,
        'E_pairs': E_pairs,
        'E_nyquist': E_nyquist,
        'total_energy': total_energy,
        'R': R  # Include raw DFT for debugging
    }

def fourier_subspace_energies_direct(x):
    """
    Directly compute Fourier-subspace energies from the Fourier basis decomposition.
    This serves as ground truth for testing the cyclic-shift method.
    
    Works for both real and complex vectors.
    
    Args:
        x (np.ndarray): Input vector of length N (real or complex)
        
    Returns:
        dict: Same format as fourier_subspace_energies()
    """
    
    N = len(x)
    alpha = standard_to_fourier_basis(x, N)
    
    # DC energy
    E_dc = np.abs(alpha[0])**2
    
    # Pair energies
    max_k = (N - 1) // 2
    E_pairs = np.zeros(max_k)
    idx = 1
    for k in range(max_k):
        # Each pair: |alpha_c|^2 + |alpha_s|^2
        E_pairs[k] = np.abs(alpha[idx])**2 + np.abs(alpha[idx + 1])**2
        idx += 2
    
    # Nyquist energy
    E_nyquist = None
    if N % 2 == 0:
        E_nyquist = np.abs(alpha[-1])**2
    
    # Total energy
    total_energy = E_dc + np.sum(E_pairs)
    if E_nyquist is not None:
        total_energy += E_nyquist
    
    return {
        'E_dc': E_dc,
        'E_pairs': E_pairs,
        'E_nyquist': E_nyquist,
        'total_energy': total_energy
    }
