"""Offline S4 boundary sums from the five existing eta=.96 Gamma caches."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import LinearNDInterpolator

from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.boundary_integrals import area_rule, complex_columns, group_totals, hole_edges, parity_errors, phase_and_norm, polygon_area

ROOT = Path(__file__).resolve().parents[2] / "scripts/.out/unit_cell_2D"
ZETAS = [1.154, 1.155, 1.156, 1.157, 1.158]


def sampled_edges(edges, xy_m, hz, method):
    """Integrate only cached points lying on each actual edge, including endpoints."""
    rows = []
    for e in edges:
        delta = e.end - e.start
        t = (xy_m - e.start) @ delta / e.length**2
        distance = np.linalg.norm(xy_m - e.start - t[:, None]*delta, axis=1)
        keep = (distance < 1e-15) & (t >= -1e-9) & (t <= 1+1e-9)
        samples = pd.DataFrame({"t": np.round(t[keep], 12), "hz": hz[keep]})
        samples = samples.groupby("t", as_index=False).mean().sort_values("t")
        if len(samples) < 3 or abs(samples.t.iloc[0]) > 1e-9 or abs(samples.t.iloc[-1]-1) > 1e-9:
            raise ValueError(f"Incomplete cached edge {e.hole}/{e.edge}")
        integral = method(samples.hz.to_numpy(), x=samples.t.to_numpy()*e.length)
        rows.append({"hole": e.hole, "edge": e.edge, "group": e.group,
                     "sample_count": len(samples), "length_m": e.length,
                     "nx": e.normal[0], "ny": e.normal[1], "integral_hz": integral,
                     "qx": e.normal[0]*integral, "qy": e.normal[1]*integral})
    return rows


def self_check():
    outer, holes, _ = create_hexagon_design(.82, [[.2624, 0, .24, 0]]*6)
    edges = hole_edges(np.asarray(holes)*1e-6)
    points = np.concatenate([e.start + np.linspace(0, 1, 13)[:, None]*(e.end-e.start)
                             for e in edges])
    values = (2+3j)*points[:, 0]
    expected = -(2+3j)*sum(abs(polygon_area(h*1e-6)) for h in np.asarray(holes))
    for method in (trapezoid, simpson):
        rows = sampled_edges(edges, points, values, method)
        totals = group_totals(rows)
        assert len(rows) == 18 and all(r["sample_count"] == 13 for r in rows)
        np.testing.assert_allclose(totals["Qx"], expected, rtol=1e-10, atol=1e-25)
        np.testing.assert_allclose(totals["QAx"]+totals["QBx"], expected, rtol=1e-10, atol=1e-25)
    h = np.array([1+2j, 3-1j, -2+1j])
    reference = h*phase_and_norm(h, np.ones(3))
    target = h*4*np.exp(.7j)
    np.testing.assert_allclose(target*phase_and_norm(target, np.ones(3), reference), reference)
    print("Analytic boundary and phase checks passed.")


def run(output):
    output = output.resolve()
    if not output.is_relative_to(ROOT.resolve()) or output == ROOT.resolve():
        raise ValueError("Output must be a new task directory under unit_cell_2D")
    if output.exists():
        raise FileExistsError(output)
    sources, cases, mode_rows = {}, [], []

    def source(path):
        sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    for zeta in ZETAS:
        series = ROOT / f"unit_cell_2D_(245-0.96-{zeta:g})_mesh5_kmax0.15_uniform25_quarter_band-py"
        config = json.loads(source(series/"99_config/config.json").read_text(encoding="utf-8"))
        identity = config["cache_identity"]
        geom = config["unit_cell_2d_geometry"]
        if geom != {"b0_nm": 245., "eta": .96, "zeta": zeta} or identity["mesh_size"] != 5:
            raise ValueError("Unexpected source identity")
        point = series / "01_results/k_points/kx=0.000000_ky=0.000000"
        metadata = json.loads(source(point/"point_metadata.json").read_text(encoding="utf-8"))
        if metadata["cache_identity"] != identity:
            raise ValueError("Point/series identity mismatch")
        a_m = config["lattice_constant_um"]*1e-6
        outer, holes, _ = create_hexagon_design(config["lattice_constant_um"], config["hole_params"])
        if zeta == 1.156:
            direct = ROOT / "S4_b245_eta0.96_zeta1.156_mesh5_Gamma_1case_20260916-gamma-v2/case000_eta0.96_zeta1.156_mesh5/99_config/realized_geometry.json"
            realized = json.loads(source(direct).read_text(encoding="utf-8"))
            np.testing.assert_allclose(holes, realized["holes"], rtol=0, atol=1e-14)
        if not cases:
            probes, weights, _ = area_rule(np.asarray(outer)*1e-6, [], 12, 4)
        else:
            np.testing.assert_allclose(outer, cases[0]["outer"], rtol=0, atol=1e-14)
        frequencies = pd.read_csv(source(point/"eigenfrequencies.csv"))
        candidates = []
        for idx, row in frequencies.iterrows():
            if row.p_weight < .8:
                continue
            data = pd.read_parquet(source(point/f"eigenmodes/{idx:02d}_Hz_center.parquet"))
            if not np.isfinite(data[["x", "y", "re", "im"]].to_numpy()).all():
                raise ValueError("Nonfinite cached field")
            grouped = data.groupby(["x", "y"], as_index=False)[["re", "im"]].mean()
            xy = grouped[["x", "y"]].to_numpy()*1e-6  # cached um -> library m
            hz = grouped.re.to_numpy()+1j*grouped.im.to_numpy()
            interp = LinearNDInterpolator(xy/1e-6, hz)
            h = interp(probes/1e-6)
            errors = parity_errors(h, interp(probes/1e-6*[-1, 1]), interp(probes/1e-6*[1, -1]))
            if not np.isfinite(h).all() or not all(np.isfinite(list(errors.values()))):
                raise ValueError("Missing parity/normalization samples")
            audit = {"zeta": zeta, "mode_idx": idx, "p_weight": row.p_weight,
                     "internal_px_weight": row.px_weight, **errors}
            mode_rows.append(audit)
            candidates.append({"mode_idx": idx, "xy": xy, "hz": hz, "h": h,
                               "audit": audit, "score": sum(errors.values()),
                               "frequency_thz": row.re, "Q": row.q,
                               "cache_rows": len(data), "unique_coordinates": len(grouped)})
        chosen = min(candidates, key=lambda c: c["score"])
        if chosen["score"] > .4 or chosen["audit"]["internal_px_weight"] < .9:
            raise ValueError("User py mode not reliably identified")
        for c in candidates:
            c["audit"]["selected"] = c is chosen
        cases.append({**chosen, "zeta": zeta, "a_m": a_m, "outer": outer,
                      "edges": hole_edges(np.asarray(holes)*1e-6)})

    ref = next(c for c in cases if c["zeta"] == 1.156)["h"]
    ref = ref*phase_and_norm(ref, weights)
    totals, edge_rows = [], []
    for case in cases:
        h = case["h"]
        overlap = abs(np.dot(weights, ref.conj()*h))/np.sqrt(
            np.dot(weights, abs(ref)**2)*np.dot(weights, abs(h)**2))
        if overlap < .9:
            raise ValueError("Insufficient cross-zeta mode overlap")
        factor = phase_and_norm(h, weights, ref)/case["a_m"]
        case["audit"]["overlap_reference"] = overlap
        for name, method in (("trapezoid", trapezoid), ("simpson", simpson)):
            rows = sampled_edges(case["edges"], case["xy"], case["hz"], method)
            total = group_totals(rows)
            for r in rows:
                edge_rows.append(complex_columns({"zeta": case["zeta"], "method": name,
                                                  **r, "Fx": r["qx"]*factor}))
            totals.append(complex_columns({"zeta": case["zeta"], "eta": .96, "method": name,
                "mode_idx": case["mode_idx"], "frequency_thz": case["frequency_thz"], "Q": case["Q"],
                "cache_rows": case["cache_rows"], "unique_coordinates": case["unique_coordinates"],
                "Hz_rms": np.sqrt(np.dot(weights, abs(h)**2)/weights.sum()),
                "phase_rotation_rad": np.angle(factor), "overlap_reference": overlap,
                "f_A_raw": total["QAx"], "f_B_raw": total["QBx"], "f_raw": total["Qx"],
                "F_A": total["QAx"]*factor, "F_B": total["QBx"]*factor, "F": total["Qx"]*factor,
                "cancellation_ratio": total["cancellation_ratio"]}))
    if any(hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha for p, sha in sources.items()):
        raise ValueError("Source changed during analysis")
    output.mkdir(parents=True)
    table = pd.DataFrame(totals)
    table.to_csv(output/"boundary_sums.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(output/"edge_integrals.csv", index=False)
    pd.DataFrame(mode_rows).to_csv(output/"mode_audit.csv", index=False)
    manifest = {"source_sha256": sources, "source_unchanged": True, "comsol_started": False,
                "eta": .96, "zetas": ZETAS, "reference_zeta": 1.156,
                "normalization": "F=f/(a*Hz_RMS), area-weighted RMS on the same full-cell quadrature grid",
                "phase": "Hz overlap with reference; never rotate individual f to real",
                "boundary": "dielectric to air; A=1,4; B=2,3,5,6; all 18 edges unsymmetrized",
                "sampling": "exact duplicate xy complex mean; on-edge distance tolerance 1e-15 m; t rounded to 12 decimals",
                "limitations": ["center cut-plane actual z and boundary trace not independently re-exported",
                                "trapezoid/Simpson difference is quadrature sensitivity, not FEM uncertainty",
                                "raw cached fields have arbitrary eigenmode amplitude",
                                "Q is cached eigenfrequency Q, not a new far-field validation"]}
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    plot(table, output)
    print(table[["zeta", "method", "F_re", "F_im", "cancellation_ratio", "Q", "overlap_reference"]].to_string(index=False))
    print(output)


def plot(table, output, *, methods=("trapezoid", "simpson"),
         method_labels=("Trapezoid", "Simpson"),
         section_label="existing center-$H_z$ cache", q_title=r"d  Cached $Q$ at $\Gamma$"):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "xtick.direction": "in", "ytick.direction": "in",
                         "pdf.fonttype": 42, "axes.linewidth": .8})
    main = table[table.method == methods[0]]
    alt = table[table.method == methods[1]]
    z = main.zeta
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.9), layout="constrained")
    a, b, c, d = axes.flat
    for key, label, color in [("F_A_re", r"Re $F_A$", "#c24a42"),
                              ("F_B_re", r"Re $F_B$", "#387ba8"), ("F_re", r"Re $F$", "#252525")]:
        a.plot(z, main[key], "o-", color=color, label=label, ms=4)
    a.set(ylabel="Normalized boundary sum (1)", title="a  A/B contributions")
    a.legend(frameon=False, ncol=3, fontsize=9, loc="center right", bbox_to_anchor=(1, .7))
    b.plot(z, main.F_re, "o-", color="#252525", label=r"Re $F$", ms=4)
    b.plot(z, main.F_im, "s-", color="#b16d2e", label=r"Im $F$", ms=4)
    b.set(ylabel="Normalized boundary sum (1)", title="b  Complex total")
    b.legend(frameon=False)
    for ax in (a, b):
        ax.axhline(0, color=".8", linewidth=.7, zorder=0)
    c.plot(z, main.cancellation_ratio, "o-", color="#252525", label=method_labels[0], ms=4)
    c.plot(z, alt.cancellation_ratio, "s-", color="#89939c", label=method_labels[1], ms=4, mfc="white")
    c.set(ylabel=r"$|f|/(|f_A|+|f_B|)$ (1)", title="c  Cancellation residual", ylim=(0, None))
    c.legend(frameon=False)
    d.plot(z, main.Q/1e6, "o-", color="#854f99", ms=4)
    d.set(ylabel=r"Eigenfrequency $Q$ ($10^6$)", title=q_title, ylim=(0, None))
    for ax in axes.flat:
        ax.set_xlabel(r"$\zeta$ (1)")
        ax.set_xticks(ZETAS, [f"{v:.3f}" for v in ZETAS])
        ax.tick_params(which="both", direction="in")
        ax.margins(x=.08)
    fig.suptitle(r"$\eta=0.96$  |  $b=245$ nm  |  user $p_y$  |  " + section_label, fontsize=13)
    fig.supxlabel(r"$F=f/(aH_{z,\mathrm{RMS}})$; Hz-overlap phase alignment. Points joined without a fitted zero.", fontsize=9)
    fig.savefig(output/"boundary_sum_vs_zeta.png", dpi=200)
    fig.savefig(output/"boundary_sum_vs_zeta.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    elif not args.self_check:
        parser.error("--output or --self-check required")
