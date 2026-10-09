"""Plot saved S4 edge and group integrals without importing the solver."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .output_paths import data_directory

GROUP_COLORS = {"A": "#d55e00", "B": "#0072b2", "A+B": "#222222"}


def _complex(frame, name):
    values = frame[name + "_re"].to_numpy(float) + 1j * frame[name + "_im"].to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite S4 values: {name}")
    return values


def _load_case(case_dir):
    case_dir = Path(case_dir)
    results = data_directory(case_dir)
    geometry = pd.read_csv(results / "edge_geometry.csv").sort_values(["hole", "edge"])
    edges = pd.read_csv(results / "edge_integrals.csv")
    groups = pd.read_csv(results / "group_totals.csv")
    summary = json.loads((case_dir / "99_config/s4_summary.json").read_text(encoding="utf8"))
    expected = [(hole, edge) for hole in range(1, 7) for edge in range(1, 4)]
    if list(zip(geometry.hole, geometry.edge)) != expected:
        raise ValueError("Expected exactly 18 unique numbered edges")
    if not np.isfinite(geometry.select_dtypes(include="number")).all().all():
        raise ValueError("Nonfinite edge geometry")
    if edges.duplicated(["z_over_H", "refinement", "hole", "edge"]).any() or groups.duplicated(["z_over_H", "refinement"]).any():
        raise ValueError("Duplicate S4 samples")
    refinement = int(edges.refinement.max())
    sections = sorted(edges.z_over_H.unique())
    if not np.isfinite(sections).all() or not sections:
        raise ValueError("Missing/nonfinite S4 sections")
    edges = edges[edges.refinement == refinement].sort_values(["z_over_H", "hole", "edge"])
    groups = groups[groups.refinement == refinement].sort_values("z_over_H")
    if list(groups.z_over_H) != sections:
        raise ValueError("Missing group section at selected refinement")
    expected_groups = np.array(["A" if hole in (1, 4) else "B" for hole, _ in expected])
    if not np.array_equal(geometry.group, expected_groups):
        raise ValueError("Geometry A/B numbering mismatch")
    for section in sections:
        rows = edges[edges.z_over_H == section]
        if list(zip(rows.hole, rows.edge)) != expected or not np.array_equal(rows.group, expected_groups):
            raise ValueError("Missing edge or A/B numbering mismatch at selected refinement")
    raw_scale = float(np.max(abs(_complex(edges, "integral_hz"))))
    scale = raw_scale if raw_scale else 1.0
    for axis in ("x", "y"):
        for section in sections:
            rows = edges[edges.z_over_H == section]
            total = groups[groups.z_over_H == section]
            for group in ("A", "B", ""):
                selected = rows if not group else rows[rows.group == group]
                actual = _complex(selected, "q" + axis).sum()
                if not np.allclose(actual, _complex(total, "Q" + group + axis)[0], rtol=1e-9, atol=scale * 1e-12):
                    raise ValueError("Saved group total disagrees with edge contributions")
    if 0.0 not in sections:
        raise ValueError("Missing z=0 data; re-extract the saved model before plotting")
    edges = edges[edges.z_over_H == 0]
    groups = groups[groups.z_over_H == 0]
    area = pd.read_csv(results / "area_comparison.csv")
    area = area[(area.z_over_H == 0) & (area.refinement == refinement)]
    if len(area) != 1:
        raise ValueError("Expected one z=0 area comparison at selected refinement")
    raw_scale = float(np.max(abs(_complex(edges, "integral_hz"))))
    scale = raw_scale if raw_scale else 1.0
    return geometry, edges, groups, area, summary, [0.0], refinement, scale, raw_scale


def _figure(size):
    fig = Figure(figsize=size)
    FigureCanvasAgg(fig)
    return fig


def _load_area_field(case_dir, area, meta, refinement):
    path = data_directory(case_dir) / "area_field_z0.npz"
    if not path.is_file():
        raise ValueError("Missing area_field_z0.npz; extract the saved model with --area-field-only")
    with np.load(path, allow_pickle=False) as saved:
        field = {key: saved[key] for key in saved.files}
    xy, weights, electric = (field[key] for key in ("xy_m", "weights_m2", "electric_vm"))
    x, y, grid, inside = (field[key] for key in ("grid_x_m", "grid_y_m", "grid_electric_vm", "inside"))
    if (field["z_nm"] != 0 or field["refinement"] != refinement or field["mode_idx"] != meta["mode_idx"]
            or xy.shape != (len(weights), 2) or electric.shape != (len(weights), 2)
            or grid.shape != (len(y), len(x), 2) or inside.shape != grid.shape[:2]
            or inside.dtype != bool or not inside.any() or np.any(weights <= 0)
            or not all(np.isfinite(v).all() for v in (xy, weights, electric, x, y, grid[inside], field["outer_m"]))
            or np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0)):
        raise ValueError("Invalid z=0 area field identity, coordinates or samples")
    if not np.isclose(weights.sum(), area.area_m2.iloc[0], rtol=1e-10, atol=0):
        raise ValueError("Area field weights disagree with cell area")
    factor = _complex(area, "normalization_factor")[0]
    if not np.isclose(field["normalization_factor"], factor, rtol=1e-12, atol=0):
        raise ValueError("Area field normalization disagrees with saved integral")
    for index, axis in enumerate(("x", "y")):
        integral = weights @ electric[:, index]
        tolerance = max(float(weights @ abs(electric[:, index])) * 1e-10, 1e-30)
        expected = _complex(area, "E" + axis + "_integral")[0]
        if (not np.isclose(integral, expected, rtol=1e-9, atol=tolerance)
                or not np.isclose(expected, _complex(area, "E" + axis + "_area")[0] * area.area_m2.iloc[0],
                                  rtol=1e-9, atol=tolerance)):
            raise ValueError("Area field quadrature disagrees with saved integral")
    return field


def build_s4_validation_figure(case_dir):
    """Expose S31 mismatch and integration sensitivity using saved z=0 data only."""
    case_dir = Path(case_dir)
    meta = json.loads((case_dir / "99_config/s4_summary.json").read_text(encoding="utf8"))
    area = pd.read_csv(data_directory(case_dir) / "area_comparison.csv")
    return build_validation_figure(area, meta)


def build_validation_figure(area, meta):
    """Render integral sensitivity without assuming how the fields were obtained."""
    section = area["z_nm"] if "z_nm" in area else area["z_over_H"]
    area = area[section == 0].sort_values("refinement")
    if len(area) < 2 or area.refinement.duplicated().any():
        raise ValueError("S31 diagnostics require at least two unique z=0 integration levels")
    counts = area.area_samples.to_numpy(float)
    if not np.isfinite(counts).all() or np.any(counts <= 0) or np.any(np.diff(counts) <= 0):
        raise ValueError("Invalid area quadrature sample counts")
    values = np.column_stack([_complex(area, key) for key in ("Ey_area", "Ly", "Oy", "Zy")])
    reconstructed = values[:, 1:].sum(axis=1)
    if not np.allclose(reconstructed, _complex(area, "predicted_Ey"), rtol=1e-9, atol=1e-12):
        raise ValueError("Saved reconstruction disagrees with B+O+Z")
    fine = values[-1]
    relative = fine[0] != 0
    reference = fine[0] if relative else 1 + 0j
    normalized = values / reference
    s31_error = abs(fine[0] - fine[1]) / abs(reference)
    full_error = abs(fine[0] - reconstructed[-1]) / abs(reference)
    area_change = abs(values[-1, 0] - values[-2, 0]) / abs(reference)
    ylabel = r"Re(value / $A_f$) (1)" if relative else "Re(value) (V/m)"
    metric_unit = "" if relative else " V/m"
    with mpl.rc_context({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.direction": "in", "ytick.direction": "in", "svg.fonttype": "none"}):
        fig = _figure((14, 5.8))
        axes = fig.subplots(1, 3)
        fig.subplots_adjust(left=.065, right=.975, top=.73, bottom=.29, wspace=.36)
        frequency_label = meta.get("frequency_label") or f"{meta['frequency_thz']:.5f} THz"
        fig.suptitle(rf"S31 check | $p_y$ / $E_y$ | z=0 nm | {frequency_label} | mesh {meta['mesh']}",
                     y=.975, fontsize=12)
        fig.text(.5, .905, meta.get("validation_formula", r"$A=\langle E_y\rangle_{area},\quad B=\mathcal{K}_\Gamma\oint H_z n_x\,dl$")
                 + r"    |    S31: $A=B$    |    3D balance: $A=B+O+Z$", ha="center", fontsize=10)
        ref_note = "One complex reference for every panel and level: fine area = 1. Ratios are sensitive near cancellation."
        if not relative:
            ref_note = "Fine area is zero: absolute values are shown; no division by the area integral."
        fig.text(.5, .842, ref_note, ha="center", color=".3")
        for ax in axes:
            ax.set_box_aspect(3 / 4)
            ax.axhline(0, color=".65", lw=.7, zorder=0)
            ax.set_ylabel(ylabel)
        comparison, precision, balance = axes
        comparison.bar([0, 1], normalized[-1, :2].real, color=["#222222", "#d55e00"], width=.55)
        comparison.set(xticks=[0, 1], xticklabels=["Area A", "Boundary B"], title="(a) Does S31 agree?")
        comparison.axhline(normalized[-1, 0].real, color=".35", lw=.8, ls="--")
        comparison.margins(y=.25)
        for x, value in enumerate(normalized[-1, :2].real):
            comparison.annotate(f"{value:+.3f}", (x, value), xytext=(0, 7 if value >= 0 else -7),
                                textcoords="offset points", ha="center", va="bottom" if value >= 0 else "top")
        ratio = normalized[-1, 1]
        comparison.text(.5, -.26, f"S31 residual = {s31_error:.4g}{metric_unit}\n"
                        + (f"B / A_f = {ratio.real:.4g} {ratio.imag:+.4g}i" if relative else "Fine area = 0; absolute residual shown"),
                        transform=comparison.transAxes, ha="center", va="top", fontsize=9)
        for col, label, color in ((0, "Area A", "#222222"), (1, "Boundary B", "#d55e00")):
            precision.plot(np.arange(len(area)), normalized[:, col].real, "o-", color=color, label=label, lw=1.4, ms=5)
        precision.set(xticks=np.arange(len(area)), xticklabels=[f"{int(n):,}" for n in counts],
                      xlabel="Area quadrature points", title="(b) Sensitivity to integration precision")
        precision.legend(loc="best", frameon=False, fontsize=8)
        precision.margins(x=.2, y=.25)
        sign_change = normalized[-2, 0].real * normalized[-1, 0].real < 0
        precision.text(.5, -.28, f"Area change = {area_change:.4g}{metric_unit}\n"
                       + ("Area projection changes sign" if sign_change else "Two levels do not establish convergence"),
                       transform=precision.transAxes, ha="center", va="top", fontsize=9,
                       color="#a63603" if sign_change else ".3")
        terms = np.r_[normalized[-1], normalized[-1, 1:].sum()]
        balance.bar(np.arange(5), terms.real, color=["#222222", "#d55e00", "#0072b2", "#009e73", "#9467bd"], width=.58)
        balance.set(xticks=np.arange(5), xticklabels=["A", "B", "O", "Z", "B+O+Z"], title="(c) Does the 3D balance close?")
        balance.axhline(normalized[-1, 0].real, color=".35", lw=.8, ls="--")
        balance.margins(y=.25)
        for x, value in enumerate(terms.real):
            balance.annotate(f"{value:+.2f}", (x, value), xytext=(0, 6 if value > 0 else -6),
                             textcoords="offset points", ha="center", va="bottom" if value > 0 else "top", fontsize=8)
        z_note = "Z exported as zero; derivative unvalidated" if fine[3] == 0 else "Z: exported thickness-derivative term"
        balance.text(.5, -.26, f"Full residual = {full_error:.4g}{metric_unit}\n{z_note}",
                     transform=balance.transAxes, ha="center", va="top", fontsize=9)
        fig.text(.5, .085, "B: hole boundary; O: cell outer boundary; Z: thickness term. Residuals use full complex differences.", ha="center", fontsize=9)
        fig.text(.5, .037, "Integration levels use the same FEM solution; this figure is not a mesh-convergence test.", ha="center", color=".35", fontsize=9)
    return fig


def build_edge_figures(geometry, edges, area, meta, field, boundary, reference, *, dimensional=False):
    """Render panels from converted edge terms, normalized or in V/m."""
    cell_area = float(area.area_m2.iloc[0])
    normalized_field = field["grid_electric_vm"] if dimensional else field["grid_electric_vm"] * (cell_area / reference)
    maximum = max(1.0, float(np.max(abs(normalized_field[field["inside"]].real))),
                  *(float(np.max(abs(values.real))) for values in boundary.values()))
    display_ticks = mpl.ticker.MaxNLocator(nbins=4, symmetric=True).tick_values(-maximum, maximum)
    display_limit = float(display_ticks[-1])
    # Keep the same area normalization, but fit boundary panels to boundary
    # values rather than the much larger local electric-field amplitudes.
    edge_maximum = max(float(np.max(abs(values.real))) for values in boundary.values())
    edge_ticks = mpl.ticker.MaxNLocator(nbins=6, symmetric=True).tick_values(
        -1.08 * (edge_maximum or 1.0), 1.08 * (edge_maximum or 1.0))
    edge_limit = float(edge_ticks[-1])
    colors = [GROUP_COLORS[group] for group in geometry.group]
    segments = geometry[["x0_m", "y0_m", "x1_m", "y1_m"]].to_numpy().reshape(18, 2, 2) * 1e6
    midpoints = segments.mean(axis=1)
    normals = geometry[["nx_dielectric_to_air", "ny_dielectric_to_air"]].to_numpy()
    extent = np.max(abs(segments)) * 1.16
    span = 2 * extent
    figures = {}
    with mpl.rc_context({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.direction": "in", "ytick.direction": "in", "svg.fonttype": "none"}):
        for component in ("x", "y"):
            channel_label = (r"$p_y$: $E_y$ integral" if component == "x"
                             else r"$p_y$: $E_x$ reference")
            if dimensional and component == "x":
                channel_label = r"$p_y$: $E_y$ average"
            fig = _figure((13.5, 5.2))
            axes = fig.subplots(1, 3, squeeze=False, width_ratios=[1, 1, 4 / 3])
            fig.subplots_adjust(left=.05, right=.92, top=.80, bottom=.32, wspace=.65)
            frequency_label = meta.get("frequency_label") or f"{meta['frequency_thz']:.5f} THz"
            fig.suptitle(rf"{channel_label} | $\Gamma$ | z=0 nm | {frequency_label}", y=.975, fontsize=11)
            electric_axis, electric_index = ("y", 1) if component == "x" else ("x", 0)
            sign = "" if component == "x" else "-"
            for index, section in enumerate([0.0]):
                field_ax, spatial, contribution_ax = axes[index]
                contribution_ax.set_box_aspect(3 / 4)
                rows = edges[edges.z_over_H == section]
                values = boundary[component]
                collection = LineCollection(segments, cmap="RdBu_r", norm=mpl.colors.Normalize(-edge_limit, edge_limit), linewidths=3)
                collection.set_array(values.real)
                spatial.add_collection(collection)
                spatial.quiver(midpoints[:, 0], midpoints[:, 1], normals[:, 0] * span * .027,
                               normals[:, 1] * span * .027, angles="xy", scale_units="xy", scale=1,
                               width=.004, color="#555555", zorder=3)
                for hole in range(1, 7):
                    rows_geometry = geometry[geometry.hole == hole]
                    center = rows_geometry[["x0_m", "y0_m"]].mean().to_numpy() * 1e6
                    group = rows_geometry.group.iloc[0]
                    spatial.text(*center, f"{group}{hole}", ha="center", va="center",
                                 color=GROUP_COLORS[group], fontsize=8, fontweight="bold")
                for point, normal, edge in zip(midpoints, normals, geometry.edge):
                    spatial.text(*(point - normal * span * .025), str(edge), ha="center", va="center", fontsize=6)
                spatial.set(xlim=(-extent, extent), ylim=(-extent, extent), aspect="equal",
                            xlabel=r"x ($\mu$m)", ylabel=r"y ($\mu$m)",
                            title="(b) Boundary geometry")
                cax = spatial.inset_axes([1.03, 0, .045, 1])
                colorbar = fig.colorbar(collection, cax=cax, ticks=edge_ticks)
                edge_label = rf"Re($E_{{{electric_axis},e}}$) (V/m)" if dimensional else r"Re($c_e$) (1)"
                colorbar.set_label(edge_label, fontsize=8)
                colorbar.ax.tick_params(direction="in", labelsize=7)
                # Complex division fixes the phase by the primary Ey area integral.
                contribution_ax.bar(np.arange(18), values.real, color=colors, width=.9)
                contribution_ax.axhline(0, color=".5", lw=.6)
                contribution_ax.set(xticks=np.arange(18), xticklabels=[f"{h}:{e}" for h, e in zip(rows.hole, rows.edge)],
                                    xlim=(-.6, 17.6), ylim=(-edge_limit, edge_limit), xlabel="Hole:edge",
                                    ylabel=edge_label,
                                    title="(c) Edge contributions")
                contribution_ax.set_yticks(np.linspace(-edge_limit, edge_limit, 5))
                contribution_ax.legend(handles=[Line2D([], [], color=GROUP_COLORS[g], lw=3, label=g)
                                               for g in ('A', 'B')], loc="upper right", ncol=2, frameon=False, fontsize=8)
                contribution_ax.tick_params(axis="x", labelrotation=60, labelsize=7)
                contribution_ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
                displayed = np.ma.array(normalized_field[:, :, electric_index].real,
                                        mask=~field["inside"])
                mesh = field_ax.pcolormesh(field["grid_x_m"] * 1e6, field["grid_y_m"] * 1e6,
                                          displayed, shading="nearest", cmap="RdBu_r",
                                          vmin=-display_limit, vmax=display_limit, rasterized=True)
                field_ax.add_collection(LineCollection(segments, colors=".35", linewidths=.5))
                outer = field["outer_m"] * 1e6
                closed = np.vstack([outer, outer[0]])
                field_ax.plot(closed[:, 0], closed[:, 1], color=".25", lw=.7)
                field_extent = np.max(abs(outer)) * 1.035
                field_ax.set(xlim=(-field_extent, field_extent), ylim=(-field_extent, field_extent),
                             aspect="equal", xlabel=r"x ($\mu$m)", ylabel=r"y ($\mu$m)",
                             title=rf"(a) Area field $E_{electric_axis}$")
                field_cax = field_ax.inset_axes([1.03, 0, .045, 1])
                field_bar = fig.colorbar(mesh, cax=field_cax, ticks=display_ticks)
                field_label = (rf"Re($E_{electric_axis}$) (V/m)" if dimensional
                               else rf"Re($F_{electric_axis}=S E_{electric_axis}/I_y$) (1)")
                field_bar.set_label(field_label, fontsize=8)
                field_bar.ax.tick_params(direction="in", labelsize=7)
                if dimensional:
                    continue  # The caller supplies the three method-specific formulas.
                integral = _complex(area, "E" + electric_axis + "_integral")[0] / reference
                value_label = "1" if component == "x" else f"{integral.real:.4g} {integral.imag:+.4g}i"
                field_ax.text(.5, -.29, rf"$\frac{{1}}{{S}}\int_\Omega F_{electric_axis}\,dS"
                              + rf"=\frac{{\int_\Omega E_{electric_axis}\,dS}}{{I_y}}={value_label}$",
                              transform=field_ax.transAxes, ha="center", va="top", fontsize=10)
                total = values.sum()
                total_label = f"{total.real:.4g} {total.imag:+.4g}i"
                contribution_ax.text(.5, -.29, meta.get("boundary_formula_" + component,
                             rf"$C_\partial={sign}\frac{{i\,\Delta(1/\epsilon)}}{{\omega\epsilon_0 I_y}}"
                             + rf"\oint H_z n_{component}\,dl$")
                             + "\n" + rf"$=\sum_{{e=1}}^{{18}}c_e={total_label}$",
                             transform=contribution_ax.transAxes, ha="center", va="top", fontsize=10)
            name = "s4_edges_py" if component == "x" else "s4_edges_py_Ex_reference"
            figures[name] = fig

    return figures


def build_s4_figures(case_dir):
    """Plot only z=0; other saved sections remain numerical references."""
    geometry, edges, groups, area, meta, sections, refinement, scale, raw_scale = _load_case(case_dir)
    field = _load_area_field(case_dir, area, meta, refinement)
    reference = _complex(area, "Ey_integral")[0]
    absolute_integral = field["weights_m2"] @ abs(field["electric_vm"][:, 1])
    if abs(reference) <= 32 * np.finfo(float).eps * absolute_integral:
        raise ValueError("Ey area integral is zero or numerically unresolved; cannot normalize")
    cell_area = float(area.area_m2.iloc[0])
    jump = float(area.jump.iloc[0])
    if not np.isfinite(jump):
        raise ValueError("Nonfinite permittivity jump")
    # prefactor=i/(omega*eps0*S). Convert Hz line integrals to E area integrals first.
    conversion = cell_area * _complex(area, "prefactor")[0] * jump / reference
    boundary = {"x": conversion * _complex(edges, "qx"), "y": -conversion * _complex(edges, "qy")}
    title = (rf"S4 | $p_y$ | {meta['frequency_thz']:.5f} THz | mode {meta['mode_idx']}"
             + f" | b0={meta['b0_nm']:g} nm, eta={meta['eta']:g}, zeta={meta['zeta']:g}, mesh={meta['mesh']}")
    normalization = (rf"$C_{{ref}}=\max_e|\int_e H_z\,ds|_{{z=0}}={raw_scale:.3e}$ A; "
                     + "one x/y scale; Hz RMS normalization at z=0.")
    if raw_scale == 0:
        normalization = "All edge integrals are zero; display divisor is 1 A."
    status = "Boundary contributions only; numerical acceptance is not established by these plots."
    figures = build_edge_figures(geometry, edges, area, meta, field, boundary, reference)
    with mpl.rc_context({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.direction": "in", "ytick.direction": "in", "svg.fonttype": "none"}):
        fig = _figure((10.5, 7.6))
        axes = fig.subplots(2, 2)
        fig.subplots_adjust(left=.1, right=.975, bottom=.17, top=.82, hspace=.52, wspace=.33)
        fig.suptitle(title, y=.98, fontsize=10)
        fig.text(.5, .92, r"$C_A=\sum_{e\in A}c_e,\quad C_B=\sum_{e\in B}c_e,\quad C=C_A+C_B$ | z=0 nm", ha="center")
        for row, component in enumerate(("x", "y")):
            for col, (part, label) in enumerate(((np.real, "Re"), (np.imag, "Im"))):
                ax = axes[row, col]
                values = [part(_complex(groups, "Q" + key + component)[0]) / scale for key in ("A", "B", "")]
                ax.bar([0, 1, 2], values, color=list(GROUP_COLORS.values()), width=.58)
                ax.axhline(0, color=".65", lw=.6)
                channel = "py: Ey integral" if component == "x" else "py: Ex reference"
                ax.set(xlabel="Group", ylabel=rf"{label}($C$) / $C_{{ref}}$ (1)",
                       title=f"{'Real' if col == 0 else 'Imaginary'} group contributions | {channel}",
                       xticks=[0, 1, 2], xticklabels=["A", "B", "A+B"])
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
        fig.text(.5, .092, normalization, ha="center", fontsize=8)
        fig.text(.5, .058, "Only z=0 nm is shown; other sections are data references.", ha="center", fontsize=9)
        fig.text(.5, .028, status, ha="center", fontsize=9, color=".35")
        figures["s4_group_contributions"] = fig

        fig = _figure((12, 8))
        axes = fig.subplots(2, 2)
        fig.subplots_adjust(left=.085, right=.98, top=.82, bottom=.18, hspace=.62, wspace=.32)
        fig.suptitle(title, y=.98, fontsize=10)
        fig.text(.5, .919, r"z=0 nm | $\langle E\rangle_{area}$ versus $L+O+Z$ | same complex field normalization", ha="center")
        labels = ["Area", "Hole L", "Outer O", "3D Z", "L+O+Z"]
        for row, component in enumerate(("x", "y")):
            keys = ["E" + component + "_area", "L" + component, "O" + component,
                    "Z" + component, "predicted_E" + component]
            values = np.array([_complex(area, key)[0] for key in keys])
            for col, (part, label) in enumerate(((np.real, "Re"), (np.imag, "Im"))):
                ax = axes[row, col]
                ax.bar(np.arange(5), part(values), color=["#222222", "#d55e00", "#0072b2", "#009e73", "#9467bd"], width=.6)
                ax.axhline(0, color=".65", lw=.6)
                ax.set(xticks=np.arange(5), xticklabels=labels,
                       ylabel=rf"{label}($E_{component}$) (V/m)", title=f"{'Real' if col == 0 else 'Imaginary'} area / boundary comparison | E{component}")
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
        fig.text(.5, .10, "Area is the direct electric-field integral divided by cell area; all terms use the z=0 section.", ha="center")
        fig.text(.5, .060, "Raw area integrals and z=50/100 nm reference values remain in area_comparison.csv.", ha="center")
        fig.text(.5, .025, "Zero exported 3D terms require derivative validation; agreement is not assumed.", ha="center", color=".35")
        figures["s4_area_comparison"] = fig
    return figures


def plot_s4_case(case_dir, *, edges_only=False, validation_only=False):
    """Save figures in the case layout; never alter numerical results."""
    if edges_only and validation_only:
        raise ValueError("Choose edges-only or validation-only, not both")
    figures = ({"s4_validation_py": build_s4_validation_figure(case_dir)} if validation_only
               else build_s4_figures(case_dir))
    compact = data_directory(case_dir).name == "80_logs"
    if not compact and not edges_only and not validation_only:
        figures["s4_validation_py"] = build_s4_validation_figure(case_dir)
    output = Path(case_dir) / ("10_overview" if compact else "01_results")
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    try:
        for name, figure in figures.items():
            if compact and not validation_only and name != "s4_edges_py":
                continue
            if edges_only and not name.startswith("s4_edges_"):
                continue
            for extension in ("png", "svg"):
                path = output / f"{name}.{extension}"
                with mpl.rc_context({"svg.fonttype": "none"}):
                    figure.savefig(path, dpi=300, facecolor="white")
                paths.append(path)
    finally:
        for figure in figures.values():
            figure.clear()
    return paths
