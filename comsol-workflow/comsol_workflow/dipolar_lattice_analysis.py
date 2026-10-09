"""Saved hexagonal-quotient fields and finite-cell aperture analysis."""
import numpy as np
from .finite_lattice_fourier import direct_lattice_basis, hex_spiral_indices


def load_data(source):
    with np.load(source, allow_pickle=False) as saved:
        d = {key: saved[key] for key in (
            "cell_indices", "cell_centers", "rho_points", "k_folded", "k_grid_indices",
            "F", "weights", "weight_fraction", "mode_idx", "transform_kind",
            "cyclic_indices", "xi_values")}
    if str(d["transform_kind"]) != "hex_cyclic_quotient_v1":
        raise ValueError("This calculation requires the saved complete hexagonal quotient.")
    np.testing.assert_array_equal(hex_spiral_indices(d["k_folded"]), d["k_grid_indices"])
    d["selected"] = np.argsort(d["k_grid_indices"])[:19]
    np.testing.assert_array_equal(d["k_grid_indices"][d["selected"]], np.arange(19))
    basis = np.linalg.lstsq(d["cell_indices"], d["cell_centers"], rcond=None)[0].T
    d["period"] = float(basis[0, 0])
    np.testing.assert_allclose(basis, direct_lattice_basis(d["period"]), atol=1e-12)
    np.testing.assert_allclose(d["cell_indices"] @ basis.T, d["cell_centers"], atol=1e-12)
    d["rows"] = []
    for i in np.unique(d["cell_indices"][:, 0]):
        js = np.sort(d["cell_indices"][d["cell_indices"][:, 0] == i, 1])
        np.testing.assert_array_equal(js, np.arange(js.min(), js.max() + 1))
        d["rows"].append((int(i), int(js.min()), int(js.max())))
    d["basis"] = basis
    d["area"] = float(np.linalg.det(basis))
    d["drho"] = d["area"] / len(d["rho_points"])
    d["N"] = len(d["cell_centers"])
    d["scale"] = 2 * np.pi / d["period"]
    d["F19"] = d["F"][d["selected"]]
    d["k19"] = d["k_folded"][d["selected"]]
    d["P19"] = d["weight_fraction"][d["selected"]]
    d["Cstar"] = float(np.sqrt(d["N"] * d["area"] * d["weights"].sum()))
    d["gram"] = d["drho"] * (d["F19"] @ d["F19"].conj().T) / d["weights"].sum()
    if not all(np.isfinite(d[key]).all() for key in ("F", "k_folded", "weights")):
        raise ValueError("Non-finite source data")
    return d


def window(q, d):
    """Exact S(q)=mean_R exp(i q.R), summing each contiguous lattice row.

    This is the existing finite-cell aperture, not a second applied window.
    Geometric sums avoid allocating a (pixels, 1141) phase matrix.
    """
    uv = np.asarray(q) @ d["basis"]
    u = uv[..., 0]
    v = (uv[..., 1] + np.pi) % (2 * np.pi) - np.pi
    result = np.zeros_like(u, dtype=complex)
    for i, lo, hi in d["rows"]:
        n = hi - lo + 1
        row = n * np.sinc(n * v / (2 * np.pi)) / np.sinc(v / (2 * np.pi))
        result += np.exp(1j * (i * u + (lo + hi) * v / 2)) * row
    return result / d["N"]


def components(q, d):
    return window(np.asarray(q)[:, None, :] - d["k19"][None, :, :], d)


def scalar_spectrum(q, s, d):
    """C_n(K)/Cstar, including the sampled intracell Fourier integral."""
    out = np.empty_like(s)
    factor = np.sqrt(d["N"]) * d["drho"] / d["Cstar"]
    for start in range(0, len(q), 512):
        section = slice(start, start + 512)
        intracell = np.exp(1j * (q[section] @ d["rho_points"].T)) @ d["F19"].T
        out[section] = factor * s[section] * intracell
    return out


def coherent_terms(amplitudes):
    total = np.sum(amplitudes, axis=-1)
    incoherent = np.sum(abs(amplitudes)**2, axis=-1)
    coherent = abs(total)**2
    return total, incoherent, coherent, coherent - incoherent
