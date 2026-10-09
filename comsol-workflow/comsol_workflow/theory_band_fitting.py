"""Pure numerical calibration and fitting for the S17 four-band model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares
from scipy.optimize import linear_sum_assignment

from .c2v_band_model import hamiltonian_c5, hamiltonian_s17


PARAMETER_NAMES = ("alpha", "beta", "gamma", "mu_thz")


@dataclass(frozen=True)
class C6vCalibration:
    mass_thz: float
    mass_std_thz: float
    velocity_thz_um: float
    t0_thz: float
    t1_thz: float
    frequency_center_thz: float
    rmse_thz: float
    velocity_std_thz_um: float


@dataclass(frozen=True)
class FitDiagnostics:
    covariance: NDArray[np.float64]
    correlation: NDArray[np.float64]
    singular_values: NDArray[np.float64]
    condition_number: float
    rank: int
    degrees_of_freedom: int
    residual_variance: float
    ill_conditioned: bool


@dataclass(frozen=True)
class C2vFitResult:
    parameters: dict[str, float]
    residuals_thz: NDArray[np.float64]
    predicted_thz: NDArray[np.float64]
    cost: float
    success: bool
    message: str
    multistart_parameters: NDArray[np.float64]
    diagnostics: FitDiagnostics


def center_gamma_frequencies(frequencies_thz: ArrayLike) -> tuple[NDArray[np.float64], float]:
    """Center a four-mode dataset on its Gamma four-mode mean."""

    values = np.asarray(frequencies_thz, dtype=float)
    if values.shape[-1] != 4:
        raise ValueError("frequencies_thz must have four modes on its last axis")
    center = float(np.mean(values[0])) if values.ndim == 2 else float(np.mean(values))
    return values - center, center


def calibrate_c6v(
    *, gamma_p_frequencies_thz: ArrayLike, gamma_d_frequencies_thz: ArrayLike,
    k_abs_inv_um: ArrayLike, lower_doublet_frequencies_thz: ArrayLike,
    upper_doublet_frequencies_thz: ArrayLike, lattice_constant_um: float,
) -> C6vCalibration:
    """Extract signed mass and velocity from one C6v eta point using S12."""

    p_gamma = np.asarray(gamma_p_frequencies_thz, dtype=float).reshape(-1)
    d_gamma = np.asarray(gamma_d_frequencies_thz, dtype=float).reshape(-1)
    if p_gamma.size != 2 or d_gamma.size != 2:
        raise ValueError("Gamma p and d inputs must each contain exactly two modes")
    center = float(np.mean(np.concatenate([p_gamma, d_gamma])))
    mass = float((np.mean(d_gamma) - np.mean(p_gamma)) / 2.0)
    p_mean_variance = float(np.var(p_gamma, ddof=1) / p_gamma.size)
    d_mean_variance = float(np.var(d_gamma, ddof=1) / d_gamma.size)
    mass_std = float(np.sqrt((p_mean_variance + d_mean_variance) / 4.0))
    k_abs = np.asarray(k_abs_inv_um, dtype=float).reshape(-1)
    lower = np.asarray(lower_doublet_frequencies_thz, dtype=float)
    upper = np.asarray(upper_doublet_frequencies_thz, dtype=float)
    if lower.shape[0] != k_abs.size or upper.shape[0] != k_abs.size:
        raise ValueError("Each dispersion array must have one row per k value")
    lower_mean = np.mean(lower.reshape(k_abs.size, -1), axis=1) - center
    upper_mean = np.mean(upper.reshape(k_abs.size, -1), axis=1) - center

    def residual(x: NDArray[np.float64]) -> NDArray[np.float64]:
        energy = np.sqrt(mass * mass + (x[0] * k_abs) ** 2)
        return np.concatenate([lower_mean + energy, upper_mean - energy])

    nonzero = k_abs > 0
    if not np.any(nonzero):
        raise ValueError("At least one nonzero k point is required to determine v")
    scale = np.maximum(k_abs[nonzero], np.finfo(float).eps)
    energy = (upper_mean[nonzero] - lower_mean[nonzero]) / 2.0
    velocity_guess = np.median(np.sqrt(np.maximum(energy**2 - mass**2, 0.0)) / scale)
    fit = least_squares(residual, [max(float(velocity_guess), 1e-9)], bounds=(0.0, np.inf))
    res = residual(fit.x)
    dof = max(res.size - 1, 1)
    jtj = float((fit.jac.T @ fit.jac).item())
    variance = float(res @ res / dof)
    velocity_std = float(np.sqrt(variance / jtj)) if jtj > 0 else float("inf")
    velocity = float(fit.x[0])
    t1 = 2.0 * velocity / lattice_constant_um
    return C6vCalibration(
        mass_thz=mass, mass_std_thz=mass_std,
        velocity_thz_um=velocity, t0_thz=t1 + mass,
        t1_thz=t1, frequency_center_thz=center,
        rmse_thz=float(np.sqrt(np.mean(res**2))), velocity_std_thz_um=velocity_std,
    )


def tracked_s17_bands(
    k_points_inv_um: ArrayLike, *, t0_thz: float, t1_thz: float,
    v_thz_um: float, alpha: float, beta: float, gamma: float, mu_thz: float,
    include_quadratic: bool = False, lattice_constant_um: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """Track S17 or Appendix C5 eigenvectors by standing-wave overlap.

    The point sequence must start at Gamma and then follow one continuous ray.
    Returned columns retain ``(d_xy, d_x2_minus_y2, p_y, p_x)`` labels.
    """

    points = np.asarray(k_points_inv_um, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) == 0:
        raise ValueError("k_points_inv_um must have nonempty shape (n, 2)")
    if not np.allclose(points[0], 0.0, atol=1e-14, rtol=0.0):
        raise ValueError("A tracked ray must begin at Gamma")
    if include_quadratic and lattice_constant_um is None:
        raise ValueError("lattice_constant_um is required for quadratic tracking")
    values = np.empty((len(points), 4), dtype=float)
    vectors = np.empty((len(points), 4, 4), dtype=np.complex128)
    previous = np.eye(4, dtype=np.complex128)
    kwargs = dict(
        t0_thz=t0_thz, t1_thz=t1_thz, v_thz_um=v_thz_um,
        alpha=alpha, beta=beta, gamma=gamma, mu_thz=mu_thz,
    )
    for index, (kx, ky) in enumerate(points):
        matrix = (
            hamiltonian_c5(
                kx, ky, lattice_constant_um=float(lattice_constant_um), **kwargs
            )
            if include_quadratic else hamiltonian_s17(kx, ky, **kwargs)
        )
        raw_values, raw_vectors = np.linalg.eigh(matrix)
        overlap = np.abs(previous.conj().T @ raw_vectors) ** 2
        rows, columns = linear_sum_assignment(-overlap)
        permutation = columns[np.argsort(rows)]
        current = raw_vectors[:, permutation]
        phases = np.sum(previous.conj() * current, axis=0)
        nonzero = np.abs(phases) > 0
        current[:, nonzero] *= np.exp(-1j * np.angle(phases[nonzero]))
        values[index] = raw_values[permutation]
        vectors[index] = current
        previous = current
    return values, vectors


def _svd_diagnostics(
    jacobian: NDArray[np.float64], residuals: NDArray[np.float64],
    *, rcond: float, condition_warning: float,
) -> FitDiagnostics:
    u, singular_values, vt = np.linalg.svd(jacobian, full_matrices=False)
    del u
    cutoff = rcond * singular_values[0] if singular_values.size else 0.0
    keep = singular_values > cutoff
    rank = int(np.count_nonzero(keep))
    dof = max(int(residuals.size - rank), 0)
    variance = float(residuals @ residuals / dof) if dof else float("nan")
    inverse_squared = np.zeros_like(singular_values)
    inverse_squared[keep] = 1.0 / singular_values[keep] ** 2
    covariance = (vt.T * inverse_squared) @ vt
    covariance *= variance
    std = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(std, std)
    correlation = np.divide(
        covariance, denominator, out=np.zeros_like(covariance), where=denominator > 0
    )
    condition = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > 0 else float("inf")
    )
    return FitDiagnostics(
        covariance=covariance, correlation=correlation,
        singular_values=singular_values, condition_number=condition, rank=rank,
        degrees_of_freedom=dof, residual_variance=variance,
        ill_conditioned=rank < jacobian.shape[1] or condition > condition_warning,
    )


def gamma_initial_guess(
    gamma_modes_thz: ArrayLike, *, t0_thz: float, t1_thz: float,
) -> NDArray[np.float64]:
    """Build a bounded-fit seed from Gamma modes in standing-wave order."""

    dxy, dx2, py, px = np.asarray(gamma_modes_thz, dtype=float).reshape(4)
    mx = (dxy - px) / 2.0
    mu = max(0.0, (dx2 + py - dxy - px) / 4.0)
    mass = t0_thz - t1_thz
    gamma = mx / mass if abs(mass) > 1e-10 else 1.0
    return np.asarray([1.0, 1.0, max(gamma, 1e-6), mu], dtype=float)


def _predict_rays(
    points: NDArray[np.float64], ray_slices: Sequence[slice], *,
    t0_thz: float, t1_thz: float, v_thz_um: float,
    parameters: NDArray[np.float64],
) -> NDArray[np.float64]:
    predicted = np.empty((len(points), 4), dtype=float)
    alpha, beta, gamma, mu = parameters
    for ray_slice in ray_slices:
        ray = points[ray_slice]
        values, _ = tracked_s17_bands(
            ray, t0_thz=t0_thz, t1_thz=t1_thz, v_thz_um=v_thz_um,
            alpha=alpha, beta=beta, gamma=gamma, mu_thz=mu,
        )
        predicted[ray_slice] = values
    return predicted


def fit_c2v_parameters(
    *, k_points_inv_um: ArrayLike, observed_modes_thz: ArrayLike,
    t0_thz: float, t1_thz: float, v_thz_um: float,
    ray_slices: Sequence[slice] | None = None, weights: ArrayLike | None = None,
    initial: ArrayLike | None = None,
    lower_bounds: ArrayLike = (1e-6, 1e-6, 1e-6, 0.0),
    upper_bounds: ArrayLike = (5.0, 5.0, 5.0, 10.0),
    multistarts: int = 12, random_seed: int = 20260819,
    loss: str = "soft_l1", frequency_scale_thz: float = 1.0,
    svd_rcond: float = 1e-10, condition_warning: float = 1e8,
) -> C2vFitResult:
    """Fit ``alpha,beta,gamma,mu`` with deterministic multistart least squares."""

    points = np.asarray(k_points_inv_um, dtype=float)
    observed = np.asarray(observed_modes_thz, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or observed.shape != (len(points), 4):
        raise ValueError("Expected k_points (n,2) and observed_modes (n,4)")
    if frequency_scale_thz <= 0 or multistarts < 1:
        raise ValueError("frequency_scale_thz and multistarts must be positive")
    slices = tuple(ray_slices or (slice(0, len(points)),))
    for ray_slice in slices:
        if not np.allclose(points[ray_slice][0], 0.0, atol=1e-14, rtol=0.0):
            raise ValueError("Every tracked ray must begin at Gamma")
    point_weights = np.ones_like(observed) if weights is None else np.broadcast_to(
        np.asarray(weights, dtype=float), observed.shape
    ).copy()
    if np.any(point_weights < 0):
        raise ValueError("weights must be nonnegative squared-loss weights")
    sqrt_weight = np.sqrt(point_weights)
    lower = np.asarray(lower_bounds, dtype=float)
    upper = np.asarray(upper_bounds, dtype=float)
    seed = (
        gamma_initial_guess(observed[slices[0]][0], t0_thz=t0_thz, t1_thz=t1_thz)
        if initial is None else np.asarray(initial, dtype=float)
    )
    seed = np.clip(seed, lower + 1e-9, upper - 1e-9)

    def residual(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        predicted = _predict_rays(
            points, slices, t0_thz=t0_thz, t1_thz=t1_thz,
            v_thz_um=v_thz_um, parameters=parameters,
        )
        return ((predicted - observed) * sqrt_weight / frequency_scale_thz).ravel()

    rng = np.random.default_rng(random_seed)
    starts = [seed]
    for _ in range(multistarts - 1):
        starts.append(lower + rng.random(4) * (upper - lower))
    fits = [
        least_squares(residual, start, bounds=(lower, upper), loss=loss)
        for start in starts
    ]
    best = min(fits, key=lambda result: float(result.cost))
    predicted = _predict_rays(
        points, slices, t0_thz=t0_thz, t1_thz=t1_thz,
        v_thz_um=v_thz_um, parameters=best.x,
    )
    unscaled_residuals = predicted - observed
    diagnostics = _svd_diagnostics(
        np.asarray(best.jac, dtype=float), np.asarray(best.fun, dtype=float),
        rcond=svd_rcond, condition_warning=condition_warning,
    )
    return C2vFitResult(
        parameters={name: float(value) for name, value in zip(PARAMETER_NAMES, best.x)},
        residuals_thz=unscaled_residuals, predicted_thz=predicted,
        cost=float(best.cost), success=bool(best.success), message=str(best.message),
        multistart_parameters=np.vstack([fit.x for fit in fits]), diagnostics=diagnostics,
    )


def calibration_covariance(
    calibration: C6vCalibration, lattice_constant_um: float,
) -> NDArray[np.float64]:
    """Return covariance of ``(t0,t1,v)`` from independent mass/v estimates."""

    jacobian = np.asarray([
        [1.0, 2.0 / lattice_constant_um],
        [0.0, 2.0 / lattice_constant_um],
        [0.0, 1.0],
    ])
    source = np.diag([
        calibration.mass_std_thz**2,
        calibration.velocity_std_thz_um**2,
    ])
    return jacobian @ source @ jacobian.T


def propagate_calibration_uncertainty(
    *, calibration_mean: ArrayLike, calibration_covariance_matrix: ArrayLike,
    samples: int, random_seed: int, fit_arguments: dict[str, object],
) -> NDArray[np.float64]:
    """Refit Gaussian calibration draws and return alpha/beta/gamma/mu samples."""

    if samples < 1:
        raise ValueError("samples must be positive")
    mean = np.asarray(calibration_mean, dtype=float).reshape(3)
    covariance = np.asarray(calibration_covariance_matrix, dtype=float).reshape(3, 3)
    rng = np.random.default_rng(random_seed)
    draws = rng.multivariate_normal(mean, covariance, size=samples)
    fitted = []
    for t0, t1, velocity in draws:
        if min(t0, t1, velocity) <= 0:
            continue
        result = fit_c2v_parameters(
            **fit_arguments, t0_thz=float(t0), t1_thz=float(t1),
            v_thz_um=float(velocity),
        )
        if result.success:
            fitted.append([result.parameters[name] for name in PARAMETER_NAMES])
    if not fitted:
        raise RuntimeError("No valid calibration-uncertainty refits converged")
    return np.asarray(fitted, dtype=float)
