"""Offline window illustration and saved mode19 air-plane K=0 diagnostics.

Run from the repository root with python -B -m scripts.analysis.visualize_dipolar_broadening.
No solver imports, model loads, or source-output rewrites.
"""
from __future__ import annotations

from comsol_workflow.output_paths import output_path

import argparse
import base64
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from comsol_workflow.farfield_fft import centered_fft2, load_air_field


ROOT = Path(__file__).resolve().parents[2]
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
from comsol_workflow.figure_output import save_figure_formats

INPUTS = load_dipolar_inputs()
SERIES = INPUTS.series_dir
MODE = INPUTS.field_mode_dir
DEFAULT_OUT = INPUTS.output_dir
COLORS = ["#3477b5", "#ce6a36", "#299782"]


def kernel(s, taper):
    """FT / window area for w(x)=1-t/2+(t/2)cos(2*pi*x/L), |x|<=L/2."""
    s = np.asarray(s, dtype=float)
    if not 0 <= taper <= 1:
        raise ValueError("taper must be between 0 and 1")
    return ((1 - taper / 2) * np.sinc(s)
            + taper / 4 * (np.sinc(s - 1) + np.sinc(s + 1))) / (1 - taper / 2)


def toy_components(s, taper=1.0, left_phase=60.0, right_phase=-60.0, side_amplitude=1.0):
    amplitudes = np.array([side_amplitude * np.exp(1j * np.deg2rad(left_phase)),
                           1.0, side_amplitude * np.exp(1j * np.deg2rad(right_phase))])
    return amplitudes[:, None] * kernel(np.asarray(s)[None, :] - np.array([-1, 0, 1])[:, None], taper)


def self_check():
    x = np.linspace(-0.5, 0.5, 8001)
    s = np.array([-2.31, -1, 0, 0.41, 1, 2])
    max_error = 0.0
    for t in [0.0, 0.37, 1.0]:
        w = 1 - t / 2 + t / 2 * np.cos(2 * np.pi * x)
        numerical = np.trapezoid(w * np.exp(-2j * np.pi * s[:, None] * x), x=x, axis=1) / np.trapezoid(w, x=x)
        max_error = max(max_error, float(np.max(abs(numerical - kernel(s, t)))))
    assert max_error < 1e-7
    np.testing.assert_allclose(kernel([-2, -1, 0, 1, 2], 0), [0, 0, 1, 0, 0], atol=1e-14)
    np.testing.assert_allclose(kernel([-2, -1, 0, 1, 2], 1), [0, 0.5, 1, 0.5, 0], atol=1e-14)
    assert abs(toy_components([0], left_phase=180, right_phase=180).sum()) < 1e-14
    np.testing.assert_allclose(toy_components([0], left_phase=0, right_phase=0).sum(), 2)
    return {"analytic_kernel_quadrature_max_error": max_error, "toy_checks": "passed"}


def complex_pair(z):
    return [float(np.real(z)), float(np.imag(z))]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def real_diagnostics():
    source = MODE / "11_simulation_exports/E_air.parquet"
    source_hash = digest(source)
    fields = load_air_field(source)
    x, y = fields["x_axis"], fields["y_axis"]
    dx, dy = float(np.diff(x).mean()), float(np.diff(y).mean())
    np.testing.assert_allclose(np.diff(x), dx, atol=1e-12)
    np.testing.assert_allclose(np.diff(y), dy, atol=1e-12)
    e = np.stack([fields[c] for c in ["Ex", "Ey", "Ez"]])
    if not np.isfinite(e).all():
        raise ValueError("Non-finite source field")
    area = dx * dy * 1e-12  # Source coordinates in um; integrate over m^2.
    total = e.sum(axis=(1, 2)) * area
    if abs(total[1]) == 0:
        raise ValueError("Ey(0) is zero; cannot use it as a display reference")
    scale = total[1]
    quadrant = []
    partition = np.zeros_like(e[0].real)
    for sx, sy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
        # Share axis samples instead of duplicating or dropping them.
        wx = np.where(np.isclose(x, 0, atol=1e-10), 0.5, (sx * x > 0).astype(float))
        wy = np.where(np.isclose(y, 0, atol=1e-10), 0.5, (sy * y > 0).astype(float))
        weight = wy[:, None] * wx[None, :]
        partition += weight
        quadrant.append((e * weight).sum(axis=(1, 2)) * area)
    quadrant = np.asarray(quadrant)
    np.testing.assert_allclose(partition, 1, atol=1e-15)
    q_error = float(np.max(abs(quadrant.sum(axis=0) - total)) / abs(scale))
    mid = len(x) // 2
    fft_zero = np.array([centered_fft2(component)[mid, mid] * area for component in e])
    fft_error = float(np.max(abs(fft_zero - total)) / abs(scale))
    assert max(q_error, fft_error) < 1e-12
    eigen_path = MODE / "10_overview/eigenfrequency.csv"
    eigen = pd.read_csv(eigen_path).iloc[0]
    if int(eigen.mode_idx) != 19 or abs(eigen.re - 197.9530009044521) > 1e-8:
        raise ValueError("Unexpected source mode; refuse a silent substitution")
    summary_path = MODE / "12_farfield_FFT/farfield_summary.json"
    saved_summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    assert digest(source) == source_hash
    return {
        "source": str(source), "source_sha256": source_hash,
        "source_bytes": source.stat().st_size,
        "source_mode_id": (MODE.relative_to(SERIES / "01_results").as_posix()
                           if MODE.is_relative_to(SERIES / "01_results") else str(MODE)),
        "frequency_thz": float(eigen.re), "frequency_imag_thz": float(eigen.im), "Q": float(eigen.q),
        "grid_size": len(x), "dx_um": dx, "dy_um": dy, "plane_z_um": saved_summary["h0_um"],
        "definition": "C_i(0)=dx_m*dy_m*sum E_i(x,y,z0), spatial Fourier amplitude on the saved air plane",
        "units": "V m, assuming COMSOL E in V/m; arbitrary eigenmode normalization",
        "display_gauge": "one common multiplier 1/C_y(0) for Ex, Ey, Ez and every quadrant",
        "phase_convention": "raw saved COMSOL complex fields, no conjugation, no propagation/reference-distance prefactor",
        "total_raw": [complex_pair(z) for z in total],
        "total_normalized": [complex_pair(z / scale) for z in total],
        "quadrants_raw": [[complex_pair(z) for z in row] for row in quadrant],
        "quadrants_normalized": [[complex_pair(z / scale) for z in row] for row in quadrant],
        "abs_Ex_over_Ey": float(abs(total[0] / scale)),
        "quadrant_sum_relative_error": q_error,
        "fft_center_relative_error": fft_error,
        "limitations": [
            "Full plane was reconstructed from one quarter solution: Ex cancellation is not an independent full-device check.",
            "Quadrants are spatial partitions, not internal Bloch-k contributions.",
            "Air-plane angular spectrum only: no calibrated far-field units or independent radiation-surface calculation.",
            "Raw phase is arbitrary; the phase of a near-zero component is not physically resolved.",
            "No basis radiation f_bk, transfer T, or S52 validation is inferred."
        ],
    }


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in", "axes.titlepad": 12,
        "savefig.dpi": 170})


def save_figure(fig, out, stem):
    save_figure_formats(fig, out, stem, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def toy_figure(out):
    s = np.linspace(-3.5, 3.5, 1201)
    parts = toy_components(s)
    values = toy_components([0])[:, 0]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), layout="constrained")
    fig.suptitle("Finite-window broadening and coherent addition | scalar illustration", fontsize=14)
    x = np.linspace(-0.65, 0.65, 700)
    window = np.where(abs(x) <= .5, .5 * (1 + np.cos(2 * np.pi * x)), 0)
    ax = axes[0, 0]
    ax.fill_between(x, window, color="#dfe8ee")
    ax.plot(x, window, color="#344b62")
    ax.set(xlabel="x / L", ylabel="Window amplitude", title="a  A specified finite Hann window", ylim=(-.05, 1.1))
    ax = axes[0, 1]
    for n, c, part in zip([-1, 0, 1], COLORS, parts):
        ax.plot(s, abs(part) ** 2, color=c, label=f"n = {n:+d}")
    ax.axvline(0, color="#777", lw=.8, ls=":")
    ax.set(xlabel=r"$s=KL/(2\pi)$", ylabel="Component intensity (common scale)", title="b  Broadened components overlap", xlim=(-3.5, 3.5))
    ax.legend(frameon=False)
    ax = axes[1, 0]
    ax.plot(s, abs(parts.sum(axis=0)) ** 2, color="#152d43", label="Coherent: |sum| squared", lw=2)
    ax.plot(s, (abs(parts) ** 2).sum(axis=0), color="#969da5", ls="--", label="Incoherent: sum of intensities")
    ax.axvline(0, color="#777", lw=.8, ls=":")
    ax.set(xlabel=r"$s=KL/(2\pi)$", ylabel="Intensity (same common scale)", title="c  Relative phase changes the result", xlim=(-3.5, 3.5))
    ax.legend(frameon=False, fontsize=9)
    ax = axes[1, 1]
    start = 0j
    for n, color, z in zip([-1, 0, 1], COLORS, values):
        end = start + z
        ax.annotate("", (end.real, end.imag), (start.real, start.imag), arrowprops={"arrowstyle": "->", "color": color, "lw": 2.3})
        ax.text((start.real + end.real) / 2, (start.imag + end.imag) / 2 + .07, f"n={n:+d}", color=color, ha="center")
        start = end
    ax.scatter([start.real], [start.imag], color="#152d43", zorder=4)
    ax.text(start.real, start.imag - .14, f"sum = {start.real:.2f} + {start.imag:.2f}i", ha="center")
    ax.axhline(0, color="#bbb", lw=.7); ax.axvline(0, color="#bbb", lw=.7)
    ax.set(xlabel="Real amplitude", ylabel="Imaginary amplitude", title="d  Complex contributions at K = 0", xlim=(-.3, 1.9), ylim=(-.5, .9), aspect="equal")
    fig.supxlabel("Specified carriers n = -1, 0, +1; phases +60, 0, -60 deg. All amplitudes divided by window area.\nThis illustration is not a fitted decomposition of the cavity mode.", fontsize=9)
    save_figure(fig, out, "01_broadening_illustration")


def real_figure(out, data):
    quad = np.array(data["quadrants_normalized"])
    q = quad[:, :, 0] + 1j * quad[:, :, 1]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.1), layout="constrained")
    fig.suptitle(f"Saved mode19 | {data['frequency_thz']:.6f} THz | air-plane K = 0", fontsize=14)
    ax = axes[0]
    for idx, c in enumerate(["#3477b5", "#ce6a36", "#299782", "#9463a8"]):
        z = q[idx, 0]
        ax.annotate("", (z.real, z.imag), (0, 0), arrowprops={"arrowstyle": "->", "color": c, "lw": 2, "alpha": .7})
    ax.text(-.053, .020, "Q1, Q3", ha="center")
    ax.text(.053, -.023, "Q2, Q4", ha="center")
    ax.scatter([0], [0], color="#152d43", s=30, zorder=5)
    ax.text(0, -.007, "sum ~ 0", ha="center")
    ax.set(title="a  Ex: opposite complex contributions", xlim=(-.079, .079), ylim=(-.052, .052))
    ax = axes[1]
    start = 0j
    for idx, c in enumerate(["#3477b5", "#ce6a36", "#299782", "#9463a8"]):
        z = q[idx, 1]; end = start + z
        ax.annotate("", (end.real, end.imag), (start.real, start.imag), arrowprops={"arrowstyle": "->", "color": c, "lw": 2.7})
        ax.text((start.real + end.real) / 2, .065, f"Q{idx+1}", color=c, ha="center")
        start = end
    ax.scatter([1], [0], color="#152d43", zorder=4)
    ax.text(1, -.10, "sum = 1", ha="center")
    ax.set(title="b  Ey: four quarters add in phase", xlim=(-.15, 1.15), ylim=(-.43, .43))
    for ax in axes:
        ax.axhline(0, color="#bbb", lw=.7, zorder=0); ax.axvline(0, color="#bbb", lw=.7, zorder=0)
        ax.set(xlabel="Re [C / Cy(0)]", ylabel="Im [C / Cy(0)]", aspect="equal")
    fig.supxlabel("One common complex normalization Cy(0); axis ranges differ. Q1-Q4 are spatial quadrants, not internal k states.\nThe source is symmetry-reconstructed: cancellation is a bookkeeping check, not independent symmetry validation.", fontsize=9)
    save_figure(fig, out, "02_mode19_normal_phasors")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Run numerical checks only; create no files")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    checks = self_check()
    if args.check:
        print(json.dumps(checks)); return
    out = args.output.resolve()
    allowed = (ROOT / "results").resolve()
    if not out.is_relative_to(allowed) or out == allowed:
        raise ValueError("Output must be a dedicated folder below results")
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite an existing visualization: {out}")
    data = real_diagnostics()
    out.mkdir(parents=True)
    style(); toy_figure(out); real_figure(out, data)
    checks.update({k: data[k] for k in ["quadrant_sum_relative_error", "fft_center_relative_error"]})
    payload = {"real": data, "checks": checks,
        "toy": {"default_taper": 1, "left_phase_deg": 60, "right_phase_deg": -60,
                "side_amplitude": 1, "definition": "Three specified scalar carriers times one finite window; not the saved mode's DFT/Bloch basis."}}
    with (output_path(out, "normal_contributions.csv")).open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["region", "component", "raw_re_Vm", "raw_im_Vm", "relative_re", "relative_im", "relative_abs", "phase_deg_if_resolved"])
        for region, raw, normalized in [(f"Q{i+1}", r, n) for i, (r, n) in enumerate(zip(data["quadrants_raw"], data["quadrants_normalized"]))] + [("total", data["total_raw"], data["total_normalized"])]:
            for component, r, n in zip(["Ex", "Ey", "Ez"], raw, normalized):
                z = complex(*n)
                writer.writerow([region, component, *r, *n, abs(z), np.angle(z, deg=True) if abs(z)>1e-12 else "unresolved"])
    payload["script_sha256"] = digest(Path(__file__))
    (output_path(out, "summary.json")).write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    template = Path(__file__).with_name("visualize_s5_broadening.html").read_text(encoding="utf-8")
    html = template.replace("__S5_DATA__", json.dumps(payload, ensure_ascii=False).replace("<", "\\u003c"))
    for path in out.rglob("*"):
        mime = {".png": "image/png", ".pdf": "application/pdf", ".csv": "text/csv", ".json": "application/json"}.get(path.suffix)
        if mime and (f'href="{path.name}"' in html or f'src="{path.name}"' in html):
            uri = f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")
            html = html.replace(f'href="{path.name}"', f'href="{uri}" download="{path.name}"')
            html = html.replace(f'src="{path.name}"', f'src="{uri}"')
    (output_path(out, "index.html")).write_text(html, encoding="utf-8")
    assert digest(Path(data["source"])) == data["source_sha256"]
    print(json.dumps({"output": str(out), "checks": checks, "raw_K0": data["total_raw"], "Ex_over_Ey": data["abs_Ex_over_Ey"]}))


if __name__ == "__main__":
    main()
