"""Four-band C2v Hamiltonian from S10--S17 and Appendix C of the SI.

The canonical representation is S17 in the standing-wave basis
``(d_xy, d_x2_minus_y2, p_y, p_x)``. Frequencies and ``t0/t1/mu`` are
ordinary-frequency THz values, ``v`` is THz um, and physical wave vectors
are in 1/um. COMSOL coordinates are converted from k/G with
``G = 4*pi/(sqrt(3)*a)`` before evaluating the Hamiltonian.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


STANDING_WAVE_LABELS = ("d_xy", "d_x2_minus_y2", "p_y", "p_x")
SPIN_LABELS = ("d_plus", "p_plus", "d_minus", "p_minus")
INTERNAL_TO_THEORY_MODE = {
    "px": "p_y",
    "py": "p_x",
    "dx": "d_x2_minus_y2",
    "dy": "d_xy",
}

SW_FROM_SPIN = np.asarray(
    [
        [1j, 0.0, -1j, 0.0],
        [1.0, 0.0, 1.0, 0.0],
        [0.0, 1j, 0.0, -1j],
        [0.0, -1.0, 0.0, -1.0],
    ],
    dtype=np.complex128,
) / np.sqrt(2.0)


@dataclass(frozen=True)
class EffectiveParameters:
    """S13--S14 parameter combinations in ordinary-frequency units."""

    mx_thz: float
    my_thz: float
    chi: float
    mean_mass_thz: float
    delta_mass_thz: float
    v_plus_thz_um: float
    v_minus_thz_um: float


@dataclass(frozen=True)
class QuadraticCorrections:
    """Appendix C2 direct-projection coefficients in THz."""

    lambda_x_thz: float
    lambda_y_thz: float
    lambda_xy_thz: float
    lambda_0_thz: float
    lambda_2_thz: complex


@dataclass(frozen=True)
class BandSolution:
    """Eigenpairs evaluated at a sequence of COMSOL k/G points."""

    k_points_over_g: NDArray[np.float64]
    k_points_inv_um: NDArray[np.float64]
    eigenvalues_thz: NDArray[np.float64]
    eigenvectors: NDArray[np.complex128]


def reciprocal_scale_inv_um(lattice_constant_um: float) -> float:
    """Return ``G = 4*pi/(sqrt(3)*a)`` in 1/um."""

    if not np.isfinite(lattice_constant_um) or lattice_constant_um <= 0:
        raise ValueError("lattice_constant_um must be finite and positive")
    return float(4.0 * np.pi / (np.sqrt(3.0) * lattice_constant_um))


def k_over_g_to_inv_um(
    k_points_over_g: ArrayLike, lattice_constant_um: float
) -> NDArray[np.float64]:
    """Convert an ``(..., 2)`` array of dimensionless k/G coordinates."""

    points = np.asarray(k_points_over_g, dtype=float)
    if points.shape == (2,):
        points = points[None, :]
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("k_points_over_g must have shape (n, 2) or (2,)")
    if not np.all(np.isfinite(points)):
        raise ValueError("k_points_over_g must contain only finite values")
    return points * reciprocal_scale_inv_um(lattice_constant_um)


def effective_parameters(
    *, t0_thz: float, t1_thz: float, v_thz_um: float,
    alpha: float, beta: float, gamma: float,
) -> EffectiveParameters:
    """Evaluate the combinations defined by S13 and S14."""

    mx = gamma * (t0_thz - t1_thz)
    my = ((4.0 * beta - gamma) * t0_thz - (2.0 * alpha + gamma) * t1_thz) / 3.0
    chi = (4.0 * alpha - gamma) / 3.0
    return EffectiveParameters(
        mx_thz=float(mx), my_thz=float(my), chi=float(chi),
        mean_mass_thz=float((mx + my) / 2.0),
        delta_mass_thz=float((my - mx) / 2.0),
        v_plus_thz_um=float((chi + gamma) * v_thz_um / 2.0),
        v_minus_thz_um=float((chi - gamma) * v_thz_um / 2.0),
    )


def quadratic_corrections(
    kx_inv_um: float, ky_inv_um: float, *, lattice_constant_um: float,
    t1_thz: float, alpha: float, gamma: float,
) -> QuadraticCorrections:
    """Evaluate the direct quadratic Bloch-phase terms in Appendix C2."""

    a2t1 = lattice_constant_um**2 * t1_thz
    lambda_x = a2t1 * gamma * (kx_inv_um**2 + 3.0 * ky_inv_um**2) / 8.0
    lambda_y = a2t1 * (
        (8.0 * alpha + gamma) * kx_inv_um**2
        + 3.0 * gamma * ky_inv_um**2
    ) / 24.0
    lambda_xy = a2t1 * gamma * kx_inv_um * ky_inv_um / 4.0
    lambda_0 = (lambda_x + lambda_y) / 2.0
    lambda_2 = (lambda_y - lambda_x) / 2.0 - 1j * lambda_xy
    return QuadraticCorrections(
        lambda_x_thz=float(lambda_x), lambda_y_thz=float(lambda_y),
        lambda_xy_thz=float(lambda_xy), lambda_0_thz=float(lambda_0),
        lambda_2_thz=complex(lambda_2),
    )


def hamiltonian_s17(
    kx_inv_um: float, ky_inv_um: float, *, t0_thz: float, t1_thz: float,
    v_thz_um: float, alpha: float, beta: float, gamma: float, mu_thz: float,
) -> NDArray[np.complex128]:
    """Return canonical S17 in the standing-wave basis."""

    p = effective_parameters(
        t0_thz=t0_thz, t1_thz=t1_thz, v_thz_um=v_thz_um,
        alpha=alpha, beta=beta, gamma=gamma,
    )
    gx = gamma * v_thz_um * kx_inv_um
    gy = gamma * v_thz_um * ky_inv_um
    cx = p.chi * v_thz_um * kx_inv_um
    return np.asarray([
        [p.mx_thz - mu_thz, 0.0, -1j * gy, 1j * gx],
        [0.0, p.my_thz + mu_thz, 1j * cx, 1j * gy],
        [1j * gy, -1j * cx, -p.my_thz + mu_thz, 0.0],
        [-1j * gx, -1j * gy, 0.0, -p.mx_thz - mu_thz],
    ], dtype=np.complex128)


def hamiltonian_s16(
    kx_inv_um: float, ky_inv_um: float, *, t0_thz: float, t1_thz: float,
    v_thz_um: float, alpha: float, beta: float, gamma: float, mu_thz: float,
) -> NDArray[np.complex128]:
    """Return S16 in ``(d+, p+, d-, p-)`` order for verification."""

    p = effective_parameters(
        t0_thz=t0_thz, t1_thz=t1_thz, v_thz_um=v_thz_um,
        alpha=alpha, beta=beta, gamma=gamma,
    )
    vp_x = p.v_plus_thz_um * kx_inv_um
    vm_x = p.v_minus_thz_um * kx_inv_um
    gy = gamma * v_thz_um * ky_inv_um
    mass = p.mean_mass_thz
    delta = p.delta_mass_thz
    return np.asarray([
        [mass, -vp_x - 1j * gy, delta + mu_thz, vm_x],
        [-vp_x + 1j * gy, -mass, -vm_x, delta - mu_thz],
        [delta + mu_thz, -vm_x, mass, vp_x - 1j * gy],
        [vm_x, delta - mu_thz, vp_x + 1j * gy, -mass],
    ], dtype=np.complex128)


def hamiltonian_c5(
    kx_inv_um: float, ky_inv_um: float, *, lattice_constant_um: float,
    t0_thz: float, t1_thz: float, v_thz_um: float,
    alpha: float, beta: float, gamma: float, mu_thz: float,
) -> NDArray[np.complex128]:
    """Return Appendix C5 through direct order k^2 in standing-wave order."""

    matrix = hamiltonian_s17(
        kx_inv_um, ky_inv_um, t0_thz=t0_thz, t1_thz=t1_thz,
        v_thz_um=v_thz_um, alpha=alpha, beta=beta, gamma=gamma,
        mu_thz=mu_thz,
    )
    q = quadratic_corrections(
        kx_inv_um, ky_inv_um, lattice_constant_um=lattice_constant_um,
        t1_thz=t1_thz, alpha=alpha, gamma=gamma,
    )
    matrix[0, 0] += q.lambda_x_thz
    matrix[1, 1] += q.lambda_y_thz
    matrix[2, 2] -= q.lambda_y_thz
    matrix[3, 3] -= q.lambda_x_thz
    matrix[0, 1] = matrix[1, 0] = q.lambda_xy_thz
    matrix[2, 3] = matrix[3, 2] = q.lambda_xy_thz
    return matrix


def hamiltonian_c4(
    kx_inv_um: float, ky_inv_um: float, *, lattice_constant_um: float,
    t0_thz: float, t1_thz: float, v_thz_um: float,
    alpha: float, beta: float, gamma: float, mu_thz: float,
) -> NDArray[np.complex128]:
    """Return Appendix C4 in ``(d+, p+, d-, p-)`` order for verification."""

    matrix = hamiltonian_s16(
        kx_inv_um, ky_inv_um, t0_thz=t0_thz, t1_thz=t1_thz,
        v_thz_um=v_thz_um, alpha=alpha, beta=beta, gamma=gamma,
        mu_thz=mu_thz,
    )
    q = quadratic_corrections(
        kx_inv_um, ky_inv_um, lattice_constant_um=lattice_constant_um,
        t1_thz=t1_thz, alpha=alpha, gamma=gamma,
    )
    matrix[0, 0] += q.lambda_0_thz
    matrix[1, 1] -= q.lambda_0_thz
    matrix[2, 2] += q.lambda_0_thz
    matrix[3, 3] -= q.lambda_0_thz
    matrix[0, 2] += q.lambda_2_thz
    matrix[2, 0] += q.lambda_2_thz.conjugate()
    matrix[1, 3] += q.lambda_2_thz.conjugate()
    matrix[3, 1] += q.lambda_2_thz
    return matrix


def hamiltonian_c6v_s12(
    kx_inv_um: float, ky_inv_um: float, *, mass_thz: float, v_thz_um: float,
) -> NDArray[np.complex128]:
    """Return the isotropic standing-wave Hamiltonian S12."""

    vx = v_thz_um * kx_inv_um
    vy = v_thz_um * ky_inv_um
    return np.asarray([
        [mass_thz, 0.0, -1j * vy, 1j * vx],
        [0.0, mass_thz, 1j * vx, 1j * vy],
        [1j * vy, -1j * vx, -mass_thz, 0.0],
        [-1j * vx, -1j * vy, 0.0, -mass_thz],
    ], dtype=np.complex128)


def calculate_bands(
    k_points_over_g: ArrayLike, *, lattice_constant_um: float,
    t0_thz: float, t1_thz: float, v_thz_um: float,
    alpha: float, beta: float, gamma: float, mu_thz: float,
    include_quadratic: bool = False,
) -> BandSolution:
    """Evaluate eigenpairs using S17 or, explicitly, Appendix C5."""

    points_over_g = np.asarray(k_points_over_g, dtype=float)
    if points_over_g.shape == (2,):
        points_over_g = points_over_g[None, :]
    points = k_over_g_to_inv_um(points_over_g, lattice_constant_um)
    values = np.empty((len(points), 4), dtype=float)
    vectors = np.empty((len(points), 4, 4), dtype=np.complex128)
    kwargs = dict(
        t0_thz=t0_thz, t1_thz=t1_thz, v_thz_um=v_thz_um,
        alpha=alpha, beta=beta, gamma=gamma, mu_thz=mu_thz,
    )
    for index, (kx, ky) in enumerate(points):
        matrix = (
            hamiltonian_c5(
                kx, ky, lattice_constant_um=lattice_constant_um, **kwargs
            )
            if include_quadratic else hamiltonian_s17(kx, ky, **kwargs)
        )
        if not np.allclose(matrix, matrix.conj().T, rtol=0.0, atol=1e-12):
            raise ValueError("Hamiltonian is not Hermitian")
        values[index], vectors[index] = np.linalg.eigh(matrix)
    return BandSolution(points_over_g.copy(), points, values, vectors)
