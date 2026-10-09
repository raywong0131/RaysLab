import numpy as np

from comsol_workflow.c2v_band_model import (
    INTERNAL_TO_THEORY_MODE,
    SW_FROM_SPIN,
    calculate_bands,
    effective_parameters,
    hamiltonian_c4,
    hamiltonian_c5,
    hamiltonian_c6v_s12,
    hamiltonian_s16,
    hamiltonian_s17,
    k_over_g_to_inv_um,
    quadratic_corrections,
    reciprocal_scale_inv_um,
)


PARAMS = dict(
    t0_thz=12.0,
    t1_thz=10.0,
    v_thz_um=4.1,
    alpha=1.15,
    beta=0.91,
    gamma=0.83,
    mu_thz=0.27,
)


def test_s17_is_hermitian_and_s16_transforms_to_it():
    for kx, ky in [(0.0, 0.0), (0.13, 0.0), (0.07, -0.11)]:
        spin = hamiltonian_s16(kx, ky, **PARAMS)
        standing = hamiltonian_s17(kx, ky, **PARAMS)
        assert np.allclose(standing, standing.conj().T)
        assert np.allclose(SW_FROM_SPIN @ spin @ SW_FROM_SPIN.conj().T, standing)


def test_gamma_diagonal_and_c6v_limit():
    p = effective_parameters(**{k: PARAMS[k] for k in PARAMS if k != "mu_thz"})
    expected = [
        p.mx_thz - PARAMS["mu_thz"],
        p.my_thz + PARAMS["mu_thz"],
        -p.my_thz + PARAMS["mu_thz"],
        -p.mx_thz - PARAMS["mu_thz"],
    ]
    assert np.allclose(np.diag(hamiltonian_s17(0.0, 0.0, **PARAMS)), expected)

    c6v = dict(PARAMS, alpha=1.0, beta=1.0, gamma=1.0, mu_thz=0.0)
    mass = c6v["t0_thz"] - c6v["t1_thz"]
    assert np.allclose(
        hamiltonian_s17(0.13, -0.11, **c6v),
        hamiltonian_c6v_s12(
            0.13, -0.11, mass_thz=mass, v_thz_um=c6v["v_thz_um"]
        ),
    )


def test_time_reversal_spectrum_and_k_conversion():
    points = np.asarray([[0.01, -0.02], [-0.01, 0.02]])
    solution = calculate_bands(points, lattice_constant_um=0.82, **PARAMS)
    assert np.allclose(solution.eigenvalues_thz[0], solution.eigenvalues_thz[1])
    assert np.allclose(
        solution.k_points_inv_um,
        points * reciprocal_scale_inv_um(0.82),
    )
    assert k_over_g_to_inv_um([0.0, 0.0], 0.82).shape == (1, 2)


def test_internal_fourier_label_mapping_is_locked():
    assert INTERNAL_TO_THEORY_MODE == {
        "px": "p_y",
        "py": "p_x",
        "dx": "d_x2_minus_y2",
        "dy": "d_xy",
    }


def test_appendix_c2_quadratic_coefficients_and_axis_mixing():
    kx, ky, lattice = 0.13, -0.07, 0.82
    q = quadratic_corrections(
        kx, ky, lattice_constant_um=lattice,
        t1_thz=PARAMS["t1_thz"], alpha=PARAMS["alpha"], gamma=PARAMS["gamma"],
    )
    a2t1 = lattice**2 * PARAMS["t1_thz"]
    assert np.isclose(
        q.lambda_x_thz,
        a2t1 * PARAMS["gamma"] * (kx**2 + 3.0 * ky**2) / 8.0,
    )
    assert np.isclose(
        q.lambda_y_thz,
        a2t1 * ((8.0 * PARAMS["alpha"] + PARAMS["gamma"]) * kx**2
                + 3.0 * PARAMS["gamma"] * ky**2) / 24.0,
    )
    assert np.isclose(q.lambda_xy_thz, a2t1 * PARAMS["gamma"] * kx * ky / 4.0)
    assert quadratic_corrections(
        kx, 0.0, lattice_constant_um=lattice,
        t1_thz=PARAMS["t1_thz"], alpha=PARAMS["alpha"], gamma=PARAMS["gamma"],
    ).lambda_xy_thz == 0.0
    assert q.lambda_xy_thz != 0.0


def test_appendix_c4_transforms_to_c5_and_reduces_to_s17_at_gamma():
    kwargs = dict(PARAMS, lattice_constant_um=0.82)
    for kx, ky in [(0.0, 0.0), (0.13, 0.0), (0.07, -0.11)]:
        spin = hamiltonian_c4(kx, ky, **kwargs)
        standing = hamiltonian_c5(kx, ky, **kwargs)
        assert np.allclose(standing, standing.conj().T)
        assert np.allclose(SW_FROM_SPIN @ spin @ SW_FROM_SPIN.conj().T, standing)
    assert np.allclose(
        hamiltonian_c5(0.0, 0.0, **kwargs), hamiltonian_s17(0.0, 0.0, **PARAMS)
    )


def test_quadratic_spectrum_obeys_time_reversal_and_explicit_switch():
    points = np.asarray([[0.01, -0.02], [-0.01, 0.02]])
    quadratic = calculate_bands(
        points, lattice_constant_um=0.82, include_quadratic=True, **PARAMS
    )
    linear = calculate_bands(points, lattice_constant_um=0.82, **PARAMS)
    assert np.allclose(quadratic.eigenvalues_thz[0], quadratic.eigenvalues_thz[1])
    assert not np.allclose(quadratic.eigenvalues_thz, linear.eigenvalues_thz)
