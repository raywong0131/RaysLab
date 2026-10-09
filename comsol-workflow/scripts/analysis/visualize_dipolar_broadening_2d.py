"""Offline 19-point, two-dimensional extension of the saved lattice Hz transform.

Run: python -B -m scripts.analysis.visualize_dipolar_broadening_2d [--check]
The +i K.r transform convention preserves the saved IFFT's k labels.
"""
from __future__ import annotations

from comsol_workflow.dipolar_lattice_analysis import load_data, window, components, scalar_spectrum, coherent_terms

from comsol_workflow.output_paths import output_path

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import brentq
from scipy.interpolate import RegularGridInterpolator

from comsol_workflow.finite_lattice_fourier import (
    direct_lattice_basis, hex_spiral_indices, finite_lattice_fourier_components,
)
from comsol_workflow.farfield_fft import load_air_field


ROOT = Path(__file__).resolve().parents[2]
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
from comsol_workflow.figure_output import save_figure_formats

INPUTS = load_dipolar_inputs()
SERIES = INPUTS.series_dir
SOURCE = INPUTS.fourier_source
OUTPUT = INPUTS.output_dir










def checks(d):
    q = np.random.default_rng(42).uniform(-12, 12, (27, 2))
    brute = np.exp(1j * (q @ d["cell_centers"].T)).mean(axis=1)
    kernel_error = float(np.max(abs(window(q, d) - brute)))
    assert kernel_error < 5e-12
    knots = components(d["k19"], d)
    orth_error = float(np.max(abs(knots - np.eye(19))))
    assert orth_error < 5e-12
    knot_weights = np.einsum("qi,ij,qj->q", knots, d["gram"], knots.conj()).real
    np.testing.assert_allclose(knot_weights, d["P19"], atol=1e-12)
    s0 = components(np.zeros((1, 2)), d)
    c0 = scalar_spectrum(np.zeros((1, 2)), s0, d)[0]
    reconstructed = np.fft.fft(d["F"], axis=0, norm="ortho")
    direct_zero = reconstructed.sum() * d["drho"] / d["Cstar"]
    zero_error = float(abs(direct_zero - c0.sum()))
    assert zero_error < 5e-12
    # Independent direct sum of the 19 reconstructed components away from K=0.
    qcheck = np.array([[0.067, -0.093], [-0.211, 0.138]])
    for point in qcheck:
        direct = 0j
        for k, f in zip(d["k19"], d["F19"]):
            direct += (np.exp(1j * ((point - k) @ d["cell_centers"].T)).sum()
                       * np.sum(f * np.exp(1j * (d["rho_points"] @ point))))
        direct *= d["drho"] / (np.sqrt(d["N"]) * d["Cstar"])
        calc = scalar_spectrum(point[None], components(point[None], d), d).sum()
        np.testing.assert_allclose(calc, direct, atol=1e-12)
    return {"kernel_direct_error": kernel_error, "knot_orthogonality_error": orth_error,
            "K0_reconstruction_error_over_Cstar": zero_error}, s0[0], c0


def colorbar(fig, ax, artist, label, ticks=None):
    cax = make_axes_locatable(ax).append_axes("right", size="4.5%", pad=0.09)
    bar = fig.colorbar(artist, cax=cax, ticks=ticks)
    bar.set_label(label, rotation=270, labelpad=17)
    bar.ax.tick_params(direction="in")


def k_axes(ax, limit):
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit),
           xlabel=r"$K_x/(2\pi/a)$", ylabel=r"$K_y/(2\pi/a)$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([-0.08, -0.04, 0, 0.04, 0.08])
    ax.set_yticks([-0.08, -0.04, 0, 0.04, 0.08])
    ax.tick_params(direction="in", which="both")


def save(fig, output, stem):
    save_figure_formats(fig, output, stem, dpi=190, bbox_inches="tight")
    plt.close(fig)


def plot(d, axis, kernels, broad, spectrum, c0, widths, output):
    points = d["k19"] / d["scale"]
    limit = axis[-1]
    extent = [-limit, limit, -limit, limit]
    norm = LogNorm(1e-5, 1)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "xtick.direction": "in",
                         "ytick.direction": "in", "axes.linewidth": 0.8, "pdf.fonttype": 42})
    fig, axs = plt.subplots(1, 3, figsize=(16.6, 5.6))
    for ax in axs:
        k_axes(ax, limit)
        ax.set_facecolor("black")
    dots = axs[0].scatter(*points.T, c=d["P19"], s=70, cmap="magma", norm=norm, lw=0)
    for n, point in enumerate(points):
        axs[0].annotate(str(n), point, xytext=(0, 9), textcoords="offset points",
                        ha="center", color="white", fontsize=9)
    axs[0].set_title("a  Discrete k-weight: 1 + 6 + 12 points", loc="left")
    colorbar(fig, axs[0], dots, r"$P_n=W_n/\sum_{\mathrm{all}\ k}W_k$", [1e-5, 1e-3, 1e-1, 1])
    im = axs[1].imshow(abs(kernels[:, :, 1])**2, origin="lower", extent=extent,
                       cmap="magma", norm=norm, interpolation="nearest")
    axs[1].contour(axis, axis, abs(kernels[:, :, 1])**2, levels=[0.1, 0.5],
                   colors=["#9199a4", "white"], linewidths=[0.65, 1.05])
    axs[1].scatter(*points.T, s=16, facecolors="none", edgecolors="#b7c1cf", lw=0.65)
    axs[1].plot(0, 0, "+", color="#68d7e9", ms=9)
    axs[1].annotate("K = 0", (0, 0), xytext=(6, -15), textcoords="offset points", color="#68d7e9")
    axs[1].annotate("1", points[1], xytext=(6, 6), textcoords="offset points", color="white")
    axs[1].set_title("b  2D envelope of point 1", loc="left")
    colorbar(fig, axs[1], im, r"$|S(\mathbf{K}-\mathbf{k}_1)|^2$", [1e-5, 1e-3, 1e-1, 1])
    im = axs[2].imshow(np.maximum(broad, 1e-20), origin="lower", extent=extent,
                       cmap="magma", norm=norm, interpolation="nearest")
    for n in range(19):
        axs[2].contour(axis, axis, abs(kernels[:, :, n])**2, levels=[0.5],
                       colors="white", linewidths=0.45, alpha=0.40)
    axs[2].scatter(*points.T, s=12, facecolors="none", edgecolors="white", lw=0.6)
    axs[2].set_title("c  Coherent broadening: all 19 points", loc="left")
    colorbar(fig, axs[2], im, r"$\|\sum_n S(\mathbf{K}-\mathbf{k}_n)F_n\|_\rho^2/\sum_k W_k$",
             [1e-5, 1e-3, 1e-1, 1])
    fig.suptitle(r"$H_z$ | First three layers, original positions and indices", y=0.995)
    fig.text(0.5, 0.055, f"Retained weight: {100*d['P19'].sum():.3f}% of all 1141 points.  "
             f"Envelope intensity FWHM: x = {widths[0]:.5f}, y = {widths[1]:.5f} in 2pi/a.",
             ha="center", fontsize=10)
    fig.text(0.5, 0.015, "Center = layer 1 (index 0).  Panel b: 10% / 50% contours; "
             "panel c: 50% contours for each envelope.  1141-cell cavity aperture; no extra window.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.055, right=0.97, top=0.91, bottom=0.19, wspace=0.48)
    save(fig, output, "04_k_3layers_broadening_2d")

    fig, axs = plt.subplots(1, 3, figsize=(16.6, 5.7))
    component_scale = float(max(abs(spectrum.real).max(), abs(spectrum.imag).max()))
    for ax, part, title in zip(axs[:2], (spectrum.real, spectrum.imag),
                               (r"a  Re $C_{19}(\mathbf{K})$", r"b  Im $C_{19}(\mathbf{K})$")):
        k_axes(ax, limit)
        im = ax.imshow(part / component_scale, origin="lower", extent=extent,
                       cmap="RdBu_r", norm=Normalize(-1, 1), interpolation="nearest")
        ax.scatter(*points.T, s=12, facecolors="none", edgecolors="#616773", lw=0.5)
        ax.plot(0, 0, "+", color="black", ms=9)
        ax.set_title(title, loc="left")
        colorbar(fig, ax, im, "normalized (1)", [-1, -.5, 0, .5, 1])
    ax = axs[2]
    spacer = make_axes_locatable(ax).append_axes("right", size="4.5%", pad=0.09)
    spacer.set_axis_off()
    span = 1.45 * abs(c0.sum())
    ax.set(xlim=(-span, span), ylim=(-span, span), xlabel=r"Re $[C_n(0)/C_*]$",
           ylabel=r"Im $[C_n(0)/C_*]$")
    ax.set_aspect("equal", adjustable="box")
    ax.axhline(0, color="#dadde2", lw=0.7)
    ax.axvline(0, color="#dadde2", lw=0.7)
    ax.scatter(c0[1:].real, c0[1:].imag, c="#526a87", s=24)
    ax.annotate("", (c0[0].real, c0[0].imag), (0, 0),
                arrowprops={"arrowstyle": "->", "color": "#bc4141", "lw": 2})
    ax.annotate("0 = total\n(sampled residual)", (c0[0].real, c0[0].imag),
                xytext=(-8, 18), textcoords="offset points", ha="center", color="#9a3030")
    ax.annotate("1-18: zero analytically", (0, 0), xytext=(0, -25),
                textcoords="offset points", ha="center", color="#526a87", fontsize=9)
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0, 0), useMathText=True)
    ax.set_title("c  Complex contributions at K = 0", loc="left")
    fig.suptitle(r"Spatial Fourier amplitude $C_{19}(\mathbf{K})=\int H_z^{(19)}(\mathbf{r})e^{+i\mathbf{K}\cdot\mathbf{r}}\,d^2r$", y=0.995)
    fig.text(0.5, 0.060, f"Panels a-b share one scale: {component_scale:.5e} C*.  "
             f"C19(0)/C* = {c0.sum().real:.6e} {c0.sum().imag:+.2e}i.", ha="center", fontsize=10)
    fig.text(0.5, 0.020, "C* = sqrt(cavity area x integrated |Hz|^2).  "
             "The small origin residual is not a resolved radiation signal; Hz is not the electric far field.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=0.06, right=0.97, top=0.90, bottom=0.19, wspace=0.48)
    save(fig, output, "05_k0_complex_amplitude_2d")
    return component_scale



def plot_all_envelopes(d, axis, kernels, output, *, component="Hz", weights=None, plane_note=""):
    """Compare P_n |S(K-k_n)|^2 using the original all-k weight denominator."""
    stem = {"Hz": "06_k_19_individual_envelopes", "Ex": "07_k_19_Ex_envelopes",
            "Ey": "08_k_19_Ey_envelopes"}[component]
    # An explicit envelope export updates only these two named figures.
    axis = np.asarray(axis)
    if kernels.shape != (len(axis), len(axis), 19) or not np.isfinite(kernels).all():
        raise ValueError("Expected finite cached S19 with shape (ny, nx, 19)")
    np.testing.assert_allclose(np.diff(axis), np.diff(axis)[0], atol=1e-12)
    assert np.all(np.diff(axis) > 0)
    assert np.max(abs(kernels)) <= 1 + 1e-12
    # Compare cached values to the source-defined window before plotting.
    samples = [(0, 0), (len(axis)//2, len(axis)//2), (17, 41), (len(axis)-1, len(axis)-1)]
    q = np.array([[axis[x], axis[y]] for y, x in samples]) * d["scale"]
    np.testing.assert_allclose(np.array([kernels[y, x] for y, x in samples]),
                               components(q, d), atol=1e-12)
    points = d["k19"] / d["scale"]
    weights = np.asarray(d["P19"] if weights is None else weights, dtype=float)
    if weights.shape != (19,) or not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("Expected 19 finite nonnegative original k weights")
    norm = LogNorm(1e-5, 1)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.linewidth": 0.7,
                         "xtick.direction": "in", "ytick.direction": "in", "pdf.fonttype": 42})
    fig, axs = plt.subplots(4, 5, figsize=(15.2, 12.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.885, bottom=0.095, top=0.92, wspace=0.13, hspace=0.25)
    extent = [axis[0], axis[-1], axis[0], axis[-1]]
    images = []
    for slot, ax in enumerate(axs.flat):
        ax.set_facecolor("black")
        ax.set(xlim=extent[:2], ylim=extent[2:])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([-0.08, 0, 0.08])
        ax.set_yticks([-0.08, 0, 0.08])
        ax.tick_params(which="both", direction="in", labelsize=9, pad=3)
        if slot == 0:
            ax.scatter(*points.T, s=18, c=weights, cmap="magma", norm=norm, lw=0)
            for idx, xy in enumerate(points):
                ax.annotate(str(idx), xy, xytext=(0, 5), textcoords="offset points",
                            ha="center", color="white", fontsize=8)
            ax.set_title("Index map", loc="left", pad=7)
            continue
        n = slot - 1
        intensity = weights[n] * abs(kernels[:, :, n])**2
        images.append(ax.imshow(intensity, origin="lower", extent=extent, cmap="magma",
                                norm=norm, interpolation="nearest"))
        ax.scatter(*points.T, s=6, facecolors="none", edgecolors="#b7c1cf", lw=0.4)
        ax.plot(0, 0, "+", color="#68d7e9", ms=6, mew=0.9)
        layer = 1 if n == 0 else (2 if n <= 6 else 3)
        ax.set_title(f"{n} | Layer {layer} | P = {100 * weights[n]:.3g}%",
                     loc="left", pad=7, fontsize=10)
    field_label = {"Hz": r"$H_z$", "Ex": r"$E_x$", "Ey": r"$E_y$"}[component]
    fig.suptitle(field_label + " | Weighted 2D envelopes of all 19 k points", y=0.978, fontsize=17)
    if component == "Hz":
        definition = r"$P_n=W_n/\sum_{\mathrm{all}\ k}W_k$"
    else:
        definition = r"$P_{\alpha n}=W_{\alpha n}/\sum_{\mathrm{all}\ k}(W_{xk}+W_{yk})$"
    fig.text(0.475, 0.947, r"$I_n(\mathbf{K})=P_n|S(\mathbf{K}-\mathbf{k}_n)|^2$"
             "   |   " + definition, ha="center", fontsize=11)
    fig.text(0.475, 0.060, r"$K_x/(2\pi/a)$", ha="center", fontsize=12)
    fig.text(0.012, 0.510, r"$K_y/(2\pi/a)$", va="center", rotation=90, fontsize=12)
    footer = ("Cyan +: K = 0.  " + plane_note + "  Ex and Ey share one normalization."
              if component != "Hz" else
              f"Cyan +: K = 0.  Hz norm weights; denominator includes all {len(d['weights'])} k points.")
    fig.text(0.475, 0.023, footer, ha="center", fontsize=10)
    fig.canvas.draw()
    boxes = [ax.get_position() for ax in axs.flat]
    bottom, top = min(b.y0 for b in boxes), max(b.y1 for b in boxes)
    cax = fig.add_axes([max(b.x1 for b in boxes) + 0.020, bottom, 0.015, top-bottom])
    bar = fig.colorbar(images[0], cax=cax, ticks=10.0**np.arange(-5, 1))
    bar.set_label(r"$P_n|S(\mathbf{K}-\mathbf{k}_n)|^2$", rotation=270, labelpad=21, fontsize=13)
    bar.ax.tick_params(which="both", direction="in", labelsize=11)
    # The shared scale and aligned colorbar are part of this figure's contract.
    assert axs.flat[0].get_title(loc="left") == "Index map"
    assert len(images) == 19 and all(im.norm is norm for im in images)
    assert all(im.get_clim() == (1e-5, 1.0) for im in images)
    # Detect accidental per-panel normalization; grid maxima can miss exact centers.
    sampled_peaks = np.array([np.max(im.get_array()) for im in images])
    np.testing.assert_allclose(sampled_peaks, weights, rtol=0.005, atol=1e-15)
    positive = weights > 0
    print(component, "maximum sampled-peak relative error:",
          float(np.max(abs(sampled_peaks[positive] / weights[positive] - 1), initial=0)))
    np.testing.assert_allclose([cax.get_position().y0, cax.get_position().y1], [bottom, top], atol=1e-12)
    save(fig, output, stem)
    return [str(output_path(output, f"{stem}.{ext}")) for ext in ("png", "pdf")]



def electric_envelopes(d, source, output, axis, kernels):
    """Independently decompose saved air-plane Ex/Ey on the same cavity footprint."""
    air_path = source.parent.parent / "11_simulation_exports/E_air.parquet"
    field_hash = hashlib.sha256(air_path.read_bytes()).hexdigest()
    metadata = json.loads((air_path.parent.parent / "12_farfield_FFT/farfield_summary.json").read_text(encoding="utf-8"))
    assert int(metadata["mode_idx"]) == int(d["mode_idx"])
    np.testing.assert_allclose(metadata["period_um"], d["period"], atol=1e-12)
    z = float(metadata["h0_um"])
    fields = load_air_field(air_path)
    query = (d["cell_centers"][:, None, :] + d["rho_points"][None, :, :]).reshape(-1, 2)
    grid = (fields["y_axis"], fields["x_axis"])
    inside = RegularGridInterpolator(grid, fields["is_in_domain"].astype(float), bounds_error=True)(query[:, ::-1])
    if not np.all(inside >= 1 - 1e-12):
        raise ValueError("Cavity sample points include invalid air-field grid cells")
    cyclic = np.asarray(d["cyclic_indices"], dtype=int)
    np.testing.assert_array_equal(np.sort(cyclic), np.arange(d["N"]))
    # Same phase and indexing as the original Hz decomposition, checked directly.
    phase = np.exp(1j * (d["k19"] @ d["cell_centers"].T)) / np.sqrt(d["N"])
    profiles, all_weights, diagnostics, direct_zero = {}, {}, {}, {}
    for comp in ("Ex", "Ey"):
        values = RegularGridInterpolator(grid, fields[comp], bounds_error=True)(query[:, ::-1])
        cell_fields = values.reshape(d["N"], len(d["rho_points"]))
        if not np.isfinite(cell_fields).all():
            raise ValueError(f"Nonfinite {comp} samples")
        ordered = np.empty_like(cell_fields)
        ordered[cyclic] = cell_fields
        transformed = finite_lattice_fourier_components(ordered, d["rho_points"], d["xi_values"], d["period"])
        profiles[comp] = transformed["F"][d["selected"]].copy()
        all_weights[comp] = transformed["weights"].copy()
        direct = phase @ cell_fields
        dft_error = float(np.linalg.norm(direct - profiles[comp]) / np.linalg.norm(profiles[comp]))
        diagnostics[comp] = {"reconstruction_error": float(transformed["reconstruction_error"]),
                             "parseval_error": float(transformed["parseval_error"]),
                             "selected_direct_DFT_error": dft_error}
        assert max(diagnostics[comp].values()) < 1e-10
        direct_zero[comp] = complex(cell_fields.sum() * d["drho"])
    denominator = float(sum(w.sum() for w in all_weights.values()))
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("Nonpositive joint electric-field norm")
    fractions = {comp: w / denominator for comp, w in all_weights.items()}
    np.testing.assert_allclose(sum(w.sum() for w in fractions.values()), 1, atol=1e-12)
    cstar = float(np.sqrt(d["N"] * d["area"] * denominator))
    s0 = components(np.zeros((1, 2)), d)[0]
    amplitudes, gram, paths = {}, {}, []
    for comp in ("Ex", "Ey"):
        amplitudes[comp] = np.sqrt(d["N"]) * d["drho"] * profiles[comp].sum(axis=1) * s0
        error = float(abs(amplitudes[comp].sum() - direct_zero[comp]) / cstar)
        diagnostics[comp]["K0_vs_direct_integral_over_Cstar"] = error
        assert error < 1e-12
        gram[comp] = d["drho"] * (profiles[comp] @ profiles[comp].conj().T) / denominator
        np.testing.assert_allclose(gram[comp].diagonal().real, fractions[comp][d["selected"]], atol=1e-12)
        paths.extend(plot_all_envelopes(
            d, axis, kernels, output, component=comp, weights=fractions[comp][d["selected"]],
            plane_note=f"Air plane z = {z:g} um; {d['N']}-cell cavity footprint."))
    np.savez_compressed(output_path(output, "electric_19_components.npz"),
                        F_Ex_19=profiles["Ex"], F_Ey_19=profiles["Ey"],
                        W_Ex_all=all_weights["Ex"], W_Ey_all=all_weights["Ey"],
                        P_Ex_all_joint=fractions["Ex"], P_Ey_all_joint=fractions["Ey"],
                        Gram_Ex_19=gram["Ex"], Gram_Ey_19=gram["Ey"],
                        C0_Ex_19_V_m=amplitudes["Ex"]*1e-12, C0_Ey_19_V_m=amplitudes["Ey"]*1e-12,
                        rho_points_um=d["rho_points"], k19_um_inv=d["k19"], selected_rows=d["selected"],
                        joint_norm_denominator=denominator, Cstar_V_m=cstar*1e-12)
    with (output_path(output, "electric_19_weights_and_K0.csv")).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "layer", "kx_over_2pi_a", "ky_over_2pi_a", "Px_joint", "Py_joint",
                         "Cx0_V_m_re", "Cx0_V_m_im", "Cy0_V_m_re", "Cy0_V_m_im"])
        for n, t in enumerate(d["selected"]):
            cx, cy = amplitudes["Ex"][n]*1e-12, amplitudes["Ey"][n]*1e-12
            writer.writerow([n, 1 if n == 0 else (2 if n <= 6 else 3), *(d["k19"][n]/d["scale"]),
                             fractions["Ex"][t], fractions["Ey"][t], cx.real, cx.imag, cy.real, cy.imag])
    assert hashlib.sha256(air_path.read_bytes()).hexdigest() == field_hash
    summary = {
        "source": str(air_path.resolve()), "sha256": field_hash, "plane_z_um": z,
        "mode_idx": int(d["mode_idx"]), "source_grid_size": list(fields["Ex"].shape),
        "source_grid_step_um": [float(np.diff(fields["x_axis"])[0]), float(np.diff(fields["y_axis"])[0])],
        "cavity_cells": d["N"], "rho_samples_per_cell": len(d["rho_points"]),
        "interpolation": "complex bilinear, within valid source grid; no extrapolation",
        "Fourier_convention": "+i K.r; same original cavity quotient and k labels",
        "normalization": "D=sum_all_k(Wx_k+Wy_k); no per-polarization or per-panel normalization",
        "quantity": "electric-component squared field norm, not total electromagnetic energy",
        "scope": "air-plane field independently decomposed over the cavity footprint; not an internal-Hz-to-radiation map",
        "interference": "images are individual-component intensities; complex F and Gram matrices retained for coherent sums",
        "joint_norm_D_units": "(V/m)^2 um^2", "joint_norm_D": denominator,
        "component_fractions_all_k": {c: float(v.sum()) for c, v in fractions.items()},
        "component_fractions_19_k": {c: float(v[d["selected"]].sum()) for c, v in fractions.items()},
        "retained_fraction_within_each_component": {c: float(v[d["selected"]].sum()/v.sum()) for c,v in fractions.items()},
        "Gamma_weights_joint": {c: float(v[d["selected"][0]]) for c, v in fractions.items()},
        "C0_cavity_footprint_V_m": {c: [float(v.real*1e-12), float(v.imag*1e-12)] for c,v in direct_zero.items()},
        "Cstar_V_m": cstar*1e-12, "checks": diagnostics, "source_unchanged": True, "plots": paths}
    (output_path(output, "electric_19_summary.json")).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary





def plot_electric_interference(d, axis, result, scales, component, output, z):
    label = {"Ex": r"E_x", "Ey": r"E_y"}[component]
    c = "x" if component == "Ex" else "y"
    total, inc, coh, cross = [result[key] for key in ("total", "incoherent", "coherent", "cross")]
    limit = float(axis[-1])
    extent = [axis[0], axis[-1], axis[0], axis[-1]]
    amplitude_norm = Normalize(-scales["amplitude"], scales["amplitude"])
    intensity_norm = LogNorm(scales["intensity_min"], scales["intensity_max"])
    contrast = np.divide(cross, coh+inc, out=np.zeros_like(cross), where=(coh+inc)>0)
    contrast = np.ma.masked_where(coh+inc < scales["phase_intensity_threshold"], contrast)
    contrast_norm = Normalize(-1, 1)
    contrast_cmap = plt.get_cmap("RdBu_r").copy()
    contrast_cmap.set_bad("#777777")
    phase = np.ma.masked_where(coh < scales["phase_intensity_threshold"], np.angle(total))
    phase_cmap = plt.get_cmap("twilight_shifted").copy()
    phase_cmap.set_bad("#777777")
    panels = [
        (total.real, "RdBu_r", amplitude_norm, rf"a  Re $C_{c}$", rf"Re $(C_{c}/C_*)$"),
        (total.imag, "RdBu_r", amplitude_norm, rf"b  Im $C_{c}$", rf"Im $(C_{c}/C_*)$"),
        (phase, phase_cmap, Normalize(-np.pi, np.pi), rf"c  Phase of $C_{c}$", "Phase (rad)"),
        (np.maximum(inc, np.finfo(float).tiny), "magma", intensity_norm,
         r"d  Incoherent reference: $\sum_n|C_n|^2$", r"$\sum_n|C_n/C_*|^2$"),
        (np.maximum(coh, np.finfo(float).tiny), "magma", intensity_norm,
         r"e  Coherent sum: $|\sum_n C_n|^2$", r"$|C/C_*|^2$"),
        (contrast, contrast_cmap, contrast_norm, "f  Relative interference",
         r"$(I_{\rm coh}-I_{\rm inc})/(I_{\rm coh}+I_{\rm inc})$"),
    ]
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.linewidth": .8,
                         "xtick.direction": "in", "ytick.direction": "in", "pdf.fonttype": 42})
    fig, axs = plt.subplots(2, 3, figsize=(16, 10.2))
    images = []
    for n, (ax, (data, cmap, norm, title, bar_label)) in enumerate(zip(axs.flat, panels)):
        k_axes(ax, limit)
        ax.set_xticks([-.08, -.04, 0, .04, .08])
        ax.set_yticks([-.08, -.04, 0, .04, .08])
        im = ax.imshow(data, origin="lower", extent=extent, interpolation="nearest", cmap=cmap, norm=norm)
        images.append(im)
        points = d["k19"] / d["scale"]
        marker_color = "#838c99" if n in (3, 4) else "#565f6a"
        ax.scatter(*points.T, s=8, facecolors="none", edgecolors=marker_color, lw=.45, alpha=.65)
        ax.plot(0, 0, "+", color="#5eabbc" if n in (3, 4) else "#252c34", ms=7, mew=.9)
        ax.set_title(title, loc="left", pad=9)
        if n == 2:
            colorbar(fig, ax, im, bar_label, [-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
            fig.axes[-1].set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
        else:
            ticks = (np.linspace(-scales["amplitude"], scales["amplitude"], 5) if n in (0, 1)
                     else [-1, -.5, 0, .5, 1] if n == 5 else None)
            colorbar(fig, ax, im, bar_label, ticks)
    assert images[0].norm is images[1].norm and images[3].norm is images[4].norm
    fig.suptitle(rf"${label}$ | Complex sum of the 19 k components", fontsize=17, y=.990)
    fig.text(.5, .945, rf"$C_{c}(\mathbf{{K}})=C_{{{c},0}}(\mathbf{{K}})+\cdots+C_{{{c},18}}(\mathbf{{K}})$"
             rf"   |   Air plane $z={z:g}\,\mu$m; cavity footprint", ha="center", fontsize=11)
    fig.text(.5, .055, "Common scales across Ex and Ey.  Red in f: constructive; blue: destructive.  "
             "Gray: below the shared intensity threshold.", ha="center", fontsize=10)
    fig.text(.5, .021, r"Spatial Fourier amplitude includes the intracell form factor.  "
             r"$C_*$ is one common electric-field reference; these are not total radiated-power maps.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=.06, right=.967, top=.892, bottom=.14, wspace=.43, hspace=.43)
    stem = f"{'09' if component == 'Ex' else '10'}_{component}_19_coherent_interference"
    save(fig, output, stem)
    return [str(output_path(output, f"{stem}.{ext}")) for ext in ("png", "pdf")]


def electric_interference(d, output, axis, kernels):
    source = output_path(output, "electric_19_components.npz")
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    meta = json.loads((output_path(output, "electric_19_summary.json")).read_text(encoding="utf-8"))
    if hashlib.sha256(Path(meta["source"]).read_bytes()).hexdigest() != meta["sha256"]:
        raise ValueError("Original electric field changed since the saved decomposition")
    with np.load(source, allow_pickle=False) as saved:
        electric = {key: saved[key] for key in saved.files}
    np.testing.assert_array_equal(electric["selected_rows"], d["selected"])
    np.testing.assert_allclose(electric["k19_um_inv"], d["k19"], atol=1e-12)
    np.testing.assert_allclose(electric["rho_points_um"], d["rho_points"], atol=1e-12)
    _, test_inc, test_coh, test_cross = coherent_terms(np.array([[1,1],[1,-1],[1,1j]], dtype=complex))
    np.testing.assert_allclose(test_inc, [2,2,2])
    np.testing.assert_allclose(test_coh, [4,0,2])
    np.testing.assert_allclose(test_cross, [2,-2,0], atol=1e-14)
    xx, yy = np.meshgrid(axis, axis)
    q = np.column_stack((xx.ravel(), yy.ravel()))*d["scale"]
    s = kernels.reshape(-1,19)
    mid = len(axis)//2
    origin = mid*len(axis)+mid
    result, statistics = {}, {}
    cstar = float(electric["Cstar_V_m"])
    np.testing.assert_allclose(cstar, np.sqrt(d["N"]*d["area"]*electric["joint_norm_denominator"])*1e-12)
    for comp in ("Ex", "Ey"):
        print(f"Summing full complex {comp} amplitudes...", flush=True)
        local = dict(d, F19=electric[f"F_{comp}_19"], Cstar=cstar/1e-12)
        amplitudes = scalar_spectrum(q, s, local)
        total, inc, coh, cross = coherent_terms(amplitudes)
        for idx in [origin, 1001, 11117, len(q)-1]:
            a = amplitudes[idx]
            pairwise = 2*np.real(sum(a[i]*a[j].conjugate() for i in range(19) for j in range(i+1,19)))
            np.testing.assert_allclose(cross[idx], pairwise, atol=1e-13)
        np.testing.assert_allclose(amplitudes[origin]*cstar, electric[f"C0_{comp}_19_V_m"], atol=1e-16)
        # Check the factorization against the reconstructed spatial field at off-grid K.
        fields19 = np.exp(-1j*(d["cell_centers"]@d["k19"].T)) @ local["F19"] / np.sqrt(d["N"])
        offgrid = np.array([[.067,-.093],[-.211,.138]])
        factored = scalar_spectrum(offgrid, components(offgrid,d), local).sum(axis=1)
        direct = np.array([np.exp(1j*(d["cell_centers"]@v)) @ fields19 @ np.exp(1j*(d["rho_points"]@v))
                           for v in offgrid])*d["drho"]/local["Cstar"]
        error = float(np.max(abs(factored-direct)))
        assert error < 1e-12
        result[comp] = {key: arr.reshape(xx.shape) for key,arr in
                        zip(("total","incoherent","coherent","cross"),(total,inc,coh,cross))}
        def location(values):
            idx = int(np.argmax(values))
            return [float(v) for v in q[idx]/d["scale"]]
        c0 = complex(total[origin]*cstar)
        statistics[comp] = {
            "C0_V_m": [float(c0.real),float(c0.imag)], "C0_abs_V_m": float(abs(c0)),
            "C0_over_Cstar": [float(total[origin].real),float(total[origin].imag)],
            "I0_over_Cstar_squared": float(coh[origin]),
            "coherent_peak": float(coh.max()), "coherent_peak_K_over_2pi_a": location(coh),
            "incoherent_peak": float(inc.max()), "cross_max": float(cross.max()), "cross_min": float(cross.min()),
            "cross_max_K_over_2pi_a": location(cross), "cross_min_K_over_2pi_a": location(-cross),
            "display_window_integrated_coherent_to_incoherent": float(np.trapezoid(np.trapezoid(coh.reshape(xx.shape),x=axis,axis=1),x=axis)
                / np.trapezoid(np.trapezoid(inc.reshape(xx.shape),x=axis,axis=1),x=axis)),
            "direct_spatial_factorization_error_over_Cstar": error}
    max_intensity = max(float(r[k].max()) for r in result.values() for k in ("coherent","incoherent"))
    top = float(10**np.ceil(np.log10(max_intensity)))
    def nice_bound(value):
        unit = 10**np.floor(np.log10(value))
        return float(np.ceil(value/unit)*unit)
    scales = {"amplitude": nice_bound(max(float(max(abs(r["total"].real).max(),abs(r["total"].imag).max())) for r in result.values())),
              "intensity_min": top*1e-5, "intensity_max": top,
              "cross": nice_bound(max(float(abs(r["cross"]).max()) for r in result.values())),
              "phase_intensity_threshold": max_intensity*1e-10}
    plots = []
    for comp in ("Ex", "Ey"):
        plots += plot_electric_interference(d,axis,result[comp],scales,comp,output,float(meta["plane_z_um"]))
    arrays = {"K_axis_over_2pi_a": axis, "Cstar_V_m": cstar}
    for comp,r in result.items():
        arrays[f"C_{comp}_over_Cstar"] = r["total"]
        for key in ("incoherent","coherent","cross"):
            arrays[f"I_{comp}_{key}_over_Cstar_squared"] = r[key]
        arrays[f"relative_interference_{comp}"] = np.divide(r["cross"], r["coherent"]+r["incoherent"],
            out=np.zeros_like(r["cross"]), where=(r["coherent"]+r["incoherent"])>0)
        assert np.max(abs(arrays[f"relative_interference_{comp}"])) <= 1 + 1e-12
    np.savez_compressed(output_path(output, "electric_19_interference.npz"), **arrays)
    with (output_path(output, "electric_19_coherent_K0.csv")).open("w",newline="",encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["component","C0_V_m_re","C0_V_m_im","abs_C0_V_m","C0_over_Cstar_re","C0_over_Cstar_im","I0_over_Cstar_squared"])
        for comp,v in statistics.items():
            writer.writerow([comp,*v["C0_V_m"],v["C0_abs_V_m"],*v["C0_over_Cstar"],v["I0_over_Cstar_squared"]])
    assert hashlib.sha256(source.read_bytes()).hexdigest() == source_hash
    summary = {"source": str(source.resolve()), "sha256": source_hash, "plane_z_um": meta["plane_z_um"],
               "points": 19, "Cstar_V_m": cstar, "scales": scales, "results": statistics,
               "display_window_K_over_2pi_a": [float(axis[0]),float(axis[-1])],
               "quantity": "C_alpha(K)=sum_n sqrt(N) S(K-kn) integral_cell F_alpha,n(rho) exp(+iK.rho) d2rho",
               "interference": "abs(sum Cn)^2 - sum abs(Cn)^2, evaluated separately for Ex and Ey",
               "relative_interference_display": "eta=(Icoh-Iinc)/(Icoh+Iinc), shared [-1,1]; raw cross term also saved",
               "difference_from_previous_envelopes": "includes the intracell Fourier integral; previous plots used the intracell squared norm",
               "K0_interpretation": "nonGamma terms are zero analytically for this exact aperture; Ex total is a numerical residual, not 19 nonzero origin phasors cancelling",
               "phase_reference": "original stored eigenmode phase; no separate rotations; phase masked below one joint intensity threshold",
               "scope": "19-component air-plane cavity-footprint reconstruction; not a full-mode far-field/power calculation",
               "retained_fraction_within_each_component": meta["retained_fraction_within_each_component"],
               "inputs_unchanged": True, "plots": plots}
    (output_path(output, "electric_19_interference_summary.json")).write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--envelopes-only", action="store_true",
                        help="Update the 19 weighted envelope plots from existing output cache only.")
    parser.add_argument("--electric-envelopes", action="store_true",
                        help="Plot separate Ex/Ey air-plane envelopes with a joint denominator.")
    parser.add_argument("--electric-interference", action="store_true",
                        help="Sum the saved 19 Ex/Ey complex spatial Fourier components.")
    args = parser.parse_args()
    source_hash = hashlib.sha256(args.source.read_bytes()).hexdigest()
    d = load_data(args.source)
    check, s0, c0 = checks(d)
    if args.check:
        print(json.dumps(check, indent=2))
        return
    if args.envelopes_only or args.electric_envelopes or args.electric_interference:
        metadata = json.loads((output_path(args.output, "broadening_2d_summary.json")).read_text(encoding="utf-8"))
        if metadata["sha256"] != source_hash:
            raise ValueError("Cached envelopes do not match the requested source hash.")
        with np.load(output_path(args.output, "broadening_2d_data.npz"), allow_pickle=False) as cached:
            np.testing.assert_array_equal(cached["selected_rows"], d["selected"])
            if args.electric_interference:
                summary = electric_interference(d, args.output, cached["K_axis_over_2pi_a"], cached["S19"])
                print(json.dumps(summary, indent=2))
                return
            if args.electric_envelopes:
                summary = electric_envelopes(d, args.source, args.output,
                                             cached["K_axis_over_2pi_a"], cached["S19"])
                print(json.dumps(summary, indent=2))
                return
            paths = plot_all_envelopes(d, cached["K_axis_over_2pi_a"], cached["S19"], args.output)
        print(json.dumps({"outputs": paths, "checks": check}, indent=2))
        return
    outputs = ["04_k_3layers_broadening_2d.png", "04_k_3layers_broadening_2d.pdf",
               "05_k0_complex_amplitude_2d.png", "05_k0_complex_amplitude_2d.pdf",
               "k0_19_contributions.csv", "broadening_2d_summary.json", "broadening_2d_data.npz"]
    if any((output_path(args.output, name)).exists() for name in outputs):
        raise FileExistsError("Refusing to overwrite existing output; choose a new --output.")
    args.output.mkdir(parents=True, exist_ok=True)
    spacing = float(np.linalg.norm(d["k19"][1]) / d["scale"])
    limit = 2.8 * spacing
    axis = np.linspace(-limit, limit, 241)
    xx, yy = np.meshgrid(axis, axis)
    q = np.column_stack((xx.ravel(), yy.ravel())) * d["scale"]
    print("Evaluating 19 finite-aperture envelopes...", flush=True)
    kernels = components(q, d)
    broad = np.einsum("qi,ij,qj->q", kernels, d["gram"], kernels.conj(), optimize=True).real
    print("Integrating the saved complex intracell profiles...", flush=True)
    spectrum = scalar_spectrum(q, kernels, d).sum(axis=1).reshape(xx.shape)
    mid = len(axis) // 2
    np.testing.assert_allclose(spectrum[mid, mid], c0.sum(), atol=1e-12)
    widths = []
    for direction in np.eye(2):
        root = brentq(lambda t: abs(window(t * direction, d))**2 - 0.5,
                      0, spacing * d["scale"])
        widths.append(float(2 * root / d["scale"]))
    scale = plot(d, axis, kernels.reshape(*xx.shape, 19), broad.reshape(xx.shape),
                 spectrum, c0, widths, args.output)
    with (output_path(args.output, "k0_19_contributions.csv")).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["index", "layer", "array_row", "kx_over_2pi_a", "ky_over_2pi_a", "P_all_k",
                         "S0_re", "S0_im", "Cn0_over_Cstar_re", "Cn0_over_Cstar_im", "Cn0_A_m_re", "Cn0_A_m_im"])
        for n, t in enumerate(d["selected"]):
            writer.writerow([n, 1 if n == 0 else (2 if n <= 6 else 3), int(t),
                             *(d["k19"][n] / d["scale"]), d["P19"][n], s0[n].real, s0[n].imag,
                             c0[n].real, c0[n].imag, *(np.array([c0[n].real, c0[n].imag]) * d["Cstar"] * 1e-12)])
    np.savez_compressed(output_path(args.output, "broadening_2d_data.npz"), K_axis_over_2pi_a=axis,
                        S19=kernels.reshape(*xx.shape, 19), coherent_weight=broad.reshape(xx.shape),
                        C19_over_Cstar=spectrum, Cn0_over_Cstar=c0, selected_rows=d["selected"])
    assert hashlib.sha256(args.source.read_bytes()).hexdigest() == source_hash
    summary = {"source": str(args.source.resolve()), "sha256": source_hash,
               "mode": int(d["mode_idx"]), "period_um": d["period"], "cavity_cells": d["N"],
               "layers": [1, 6, 12], "display_indices": list(range(19)),
               "retained_weight": float(d["P19"].sum()), "spacing_over_2pi_a": spacing,
               "FWHM_intensity_xy_over_2pi_a": widths,
               "Cstar_A_m": d["Cstar"] * 1e-12,
               "C19_zero_over_Cstar": [float(c0.sum().real), float(c0.sum().imag)],
               "max_nonGamma_zero_over_Cstar": float(abs(c0[1:]).max()),
               "max_nonGamma_S0": float(abs(s0[1:]).max()),
               "Re_Im_shared_scale_over_Cstar": scale,
               "definitions": {"Fourier_sign": "+i K.r, matching saved IFFT labels",
                   "window": "S(q) = mean_R exp(i q.R), exact 1141-cell aperture; no extra taper",
                   "continuous_lattice_profile": "T19(K,rho) = sum_n S(K-kn) F_n(rho)",
                   "coherent_weight": "integral_rho abs(T19)^2 / sum_all_k Wk",
                   "spatial_amplitude": "C19(K) = sqrt(N) integral_rho exp(i K.rho) T19(K,rho)",
                   "rho_quadrature": "equal area weights, same 1241 samples as saved source",
                   "Cstar": "sqrt(N * cell_area * sum_all_k Wk), converted to A m",
                   "K0_interpretation": "nonGamma terms vanish analytically; Gamma scalar integral has small sampling residual; global eigenmode phase arbitrary",
                   "scope": "19-component Hz reconstruction on cavity cells only; excludes cladding, not electric far field"},
               "checks": check, "source_unchanged": True}
    (output_path(args.output, "broadening_2d_summary.json")).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
