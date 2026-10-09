import numpy as np

from comsol_workflow.theory_band_fitting import (
    calibrate_c6v,
    fit_c2v_parameters,
    tracked_s17_bands,
)


def test_c6v_calibration_recovers_signed_mass_velocity_and_hoppings():
    lattice_constant = 0.82
    t0, t1 = 11.0, 10.0
    velocity = lattice_constant * t1 / 2.0
    mass = t0 - t1
    center = 200.0
    k_abs = np.asarray([0.0, 0.04, 0.08, 0.12])
    energy = np.sqrt(mass**2 + (velocity * k_abs) ** 2)
    calibration = calibrate_c6v(
        gamma_p_frequencies_thz=[center - mass, center - mass],
        gamma_d_frequencies_thz=[center + mass, center + mass],
        k_abs_inv_um=k_abs,
        lower_doublet_frequencies_thz=np.column_stack(
            [center - energy, center - energy]
        ),
        upper_doublet_frequencies_thz=np.column_stack(
            [center + energy, center + energy]
        ),
        lattice_constant_um=lattice_constant,
    )
    assert np.isclose(calibration.mass_thz, mass)
    assert np.isclose(calibration.velocity_thz_um, velocity)
    assert np.isclose(calibration.t0_thz, t0)
    assert np.isclose(calibration.t1_thz, t1)
    assert calibration.rmse_thz < 1e-12


def _synthetic_two_rays(parameters):
    ray_x = np.asarray([[0.0, 0.0], [0.04, 0.0], [0.08, 0.0], [0.12, 0.0]])
    ray_y = np.asarray([[0.0, 0.0], [0.0, 0.04], [0.0, 0.08], [0.0, 0.12]])
    common = dict(t0_thz=11.0, t1_thz=10.0, v_thz_um=4.1)
    bands_x, _ = tracked_s17_bands(ray_x, **common, **parameters)
    bands_y, _ = tracked_s17_bands(ray_y, **common, **parameters)
    return np.vstack([ray_x, ray_y]), np.vstack([bands_x, bands_y]), common


def test_c2v_multistart_fit_recovers_synthetic_parameters():
    truth = dict(alpha=1.20, beta=0.91, gamma=0.82, mu_thz=0.27)
    points, observed, common = _synthetic_two_rays(truth)
    result = fit_c2v_parameters(
        k_points_inv_um=points,
        observed_modes_thz=observed,
        ray_slices=(slice(0, 4), slice(4, 8)),
        multistarts=6,
        loss="linear",
        **common,
    )
    assert result.success
    for name, expected in truth.items():
        assert np.isclose(result.parameters[name], expected, atol=1e-7)
    assert np.max(np.abs(result.residuals_thz)) < 1e-8
    assert not result.diagnostics.ill_conditioned


def test_gamma_only_fit_is_flagged_as_ill_conditioned():
    truth = dict(alpha=1.20, beta=0.91, gamma=0.82, mu_thz=0.27)
    points, observed, common = _synthetic_two_rays(truth)
    result = fit_c2v_parameters(
        k_points_inv_um=points[:1], observed_modes_thz=observed[:1],
        multistarts=2, loss="linear", **common,
    )
    assert result.diagnostics.ill_conditioned
    assert result.diagnostics.rank < 4


def test_tracking_keeps_standing_wave_labels_across_frequency_reordering():
    parameters = dict(alpha=1.3, beta=0.8, gamma=0.7, mu_thz=0.5)
    ray = np.column_stack([np.linspace(0.0, 0.5, 40), np.zeros(40)])
    values, vectors = tracked_s17_bands(
        ray, t0_thz=11.0, t1_thz=10.0, v_thz_um=4.1, **parameters
    )
    consecutive = np.abs(
        np.einsum("nij,nij->nj", vectors[:-1].conj(), vectors[1:])
    )
    assert values.shape == (40, 4)
    assert np.min(consecutive) > 0.9
