"""Consistent transparent field-image rendering for finite simulations."""

from __future__ import annotations

import base64
from html import escape
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.tri as mtri
import numpy as np


# Sampled bottom-to-top from the COMSOL 6.3 reference colorbars supplied for
# this workflow.  Keeping the nodes here makes reconstructed and full-model
# exports visually identical without requiring the reference PNGs at runtime.
WAVE_COLORS = (
    "#3D1F82",
    "#402699",
    "#422DB0",
    "#4539CE",
    "#4852ED",
    "#4D7AFC",
    "#5DA9FE",
    "#96D3F8",
    "#F2F1F2",
    "#F8B8D4",
    "#F67E8C",
    "#E24643",
    "#BE211E",
    "#9D1219",
    "#870E1C",
    "#770C1F",
    "#680A22",
)
HEAT_CAMERA_COLORS = (
    "#1B0959",
    "#2B0A71",
    "#400B84",
    "#570D90",
    "#6F1093",
    "#8A158D",
    "#A51A7E",
    "#C12267",
    "#D82F4D",
    "#E94335",
    "#F55E1D",
    "#F87A0E",
    "#F6990B",
    "#F5B51C",
    "#F6CD45",
    "#FAE484",
    "#FEF9D4",
)

WAVE_COLORMAP = LinearSegmentedColormap.from_list("COMSOL Wave", WAVE_COLORS)
HEAT_CAMERA_COLORMAP = LinearSegmentedColormap.from_list(
    "COMSOL HeatCamera",
    HEAT_CAMERA_COLORS,
)

# Approved Hz reference: negative blue, zero white, positive red.
HZ_COLORMAP = matplotlib.colormaps["RdBu_r"]
HZ_COLORBAR_TICKS = np.linspace(-1.0, 1.0, 9)


def save_hz_simulation_plot(simulation, mode_idx, expression, path, plane) -> bool:
    """Render signed 2D Hz images; leave raw COMSOL exports untouched."""
    expression_key = expression.replace(" ", "")
    if expression_key in {"ewfd.Hz", "real(ewfd.Hz)"}:
        label = "Re(Hz)"
    elif expression_key in {"ewfd.Hz*(-i)", "imag(ewfd.Hz)"}:
        label = "Im(Hz)"
    else:
        return False
    coordinates, values = simulation.get_2d_fields(mode_idx, expression, plane)
    frequency, _imaginary, quality = simulation.get_eigenfrequencies()[mode_idx]
    save_field_plot(
        Path(path), coordinates, values,
        frequency_thz=frequency, quality_factor=quality,
        quantity_label=label, unit_label="A/m",
        color_scale_mode="linearsymmetric", dpi=300,
    )
    return True


# Crop the transparent canvas to the actual annotations while retaining a
# small safety margin for antialiasing and font ascenders/descenders.
_TIGHT_SAVEFIG_OPTIONS = {
    "transparent": True,
    "facecolor": "none",
    "edgecolor": "none",
    "bbox_inches": "tight",
    "pad_inches": 0.03,
}


def _finite_field_data(
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    quantity_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    coordinates = np.asarray(coordinates, dtype=float)
    values = np.real(np.asarray(values)).reshape(-1)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("coordinates must have shape (N, 2) or (N, 3)")
    if len(coordinates) != len(values):
        raise ValueError("coordinates and values must contain the same number of rows")
    coordinates = coordinates[:, :2]
    finite_rows = np.isfinite(values) & np.all(np.isfinite(coordinates), axis=1)
    coordinates = coordinates[finite_rows]
    values = values[finite_rows]
    if len(coordinates) < 3:
        raise ValueError(
            f"Cannot plot {quantity_label}: fewer than three finite field points"
        )
    return coordinates, values


def _format_frequency(frequency_thz: float) -> str:
    value = float(frequency_thz)
    return f"frequency = {value:.2f} THz" if np.isfinite(value) else "frequency = n/a"


def _format_quality_factor(quality_factor: float) -> str:
    value = float(quality_factor)
    return f"Q = {value:.6g}" if np.isfinite(value) else "Q = n/a"


def _create_field_figure(
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    frequency_thz: float,
    quality_factor: float,
    quantity_label: str,
    unit_label: str,
    color_scale_mode: str,
    rasterize_field: bool = False,
):
    coordinates, values = _finite_field_data(
        coordinates,
        values,
        quantity_label=quantity_label,
    )
    is_hz = quantity_label in {"Re(Hz)", "Im(Hz)", "ReHz", "ImHz", "Hz"}
    if is_hz:
        limit = float(np.max(np.abs(values)))
        values = values / limit if limit > 0.0 else values.copy()
        unit_label = "1"
        color_scale_mode = "linearsymmetric"
    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1])

    fig = plt.figure(figsize=(7.0, 6.0))
    fig.patch.set_alpha(0.0)
    axis = fig.add_axes((0.065, 0.075, 0.745, 0.715))
    axis.patch.set_alpha(0.0)

    if color_scale_mode == "linearsymmetric":
        limit = float(np.max(np.abs(values), initial=0.0))
        if limit <= 0.0:
            limit = 1.0
        image = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap=HZ_COLORMAP if is_hz else WAVE_COLORMAP,
            vmin=-1.0 if is_hz else -limit,
            vmax=1.0 if is_hz else limit,
            rasterized=rasterize_field,
        )
    elif color_scale_mode == "linear":
        lower = min(float(np.min(values)), 0.0)
        upper = float(np.max(values))
        if upper <= lower:
            upper = lower + 1.0
        image = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap=HEAT_CAMERA_COLORMAP,
            vmin=lower,
            vmax=upper,
            rasterized=rasterize_field,
        )
    else:
        raise ValueError(f"Unsupported color scale mode: {color_scale_mode}")

    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$x$ ($\mu$m)", fontsize=10.0, labelpad=5.0)
    axis.set_ylabel(r"$y$ ($\mu$m)", fontsize=10.0, labelpad=5.0)
    axis.tick_params(
        axis="both",
        direction="in",
        labelsize=8.5,
        colors="black",
        length=3.5,
    )
    for spine in axis.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.9)

    # Equal-aspect adjustment can shrink the originally requested axes box.
    # Draw first, then use its final position verbatim for the colorbar.
    fig.canvas.draw()
    plot_position = axis.get_position()
    colorbar_left = plot_position.x1 + 0.035
    colorbar_axis = fig.add_axes(
        (colorbar_left, plot_position.y0, 0.022, plot_position.height)
    )
    colorbar = fig.colorbar(
        image, cax=colorbar_axis,
        **({"ticks": HZ_COLORBAR_TICKS, "format": "%.2f"} if is_hz else {}),
    )
    if is_hz:
        colorbar.set_label("normalized", rotation=270, labelpad=12)
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    colorbar.outline.set_edgecolor("black")
    colorbar.outline.set_linewidth(0.9)
    colorbar_axis.tick_params(labelsize=8.0, colors="black", length=3.0)

    metadata_y = plot_position.y1 + 0.012
    fig.text(
        plot_position.x0,
        metadata_y,
        _format_frequency(frequency_thz),
        ha="left",
        va="bottom",
        fontsize=10.5,
        color="black",
    )
    fig.text(
        plot_position.x0 + 0.60 * plot_position.width,
        metadata_y,
        _format_quality_factor(quality_factor),
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="black",
    )
    fig.text(
        plot_position.x1,
        metadata_y,
        quantity_label,
        ha="right",
        va="bottom",
        fontsize=10.5,
        color="black",
    )
    fig.text(
        colorbar_axis.get_position().x0 + colorbar_axis.get_position().width / 2.0,
        metadata_y,
        unit_label,
        ha="center",
        va="bottom",
        fontsize=10.0,
        color="black",
    )
    return fig, axis, colorbar_axis


def save_field_plot(
    path: Path,
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    frequency_thz: float,
    quality_factor: float,
    quantity_label: str,
    unit_label: str,
    color_scale_mode: str,
    dpi: int,
) -> None:
    """Save one field plot with transparent canvas and aligned colorbar."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, _axis, _colorbar_axis = _create_field_figure(
        coordinates,
        values,
        frequency_thz=frequency_thz,
        quality_factor=quality_factor,
        quantity_label=quantity_label,
        unit_label=unit_label,
        color_scale_mode=color_scale_mode,
    )
    try:
        fig.savefig(
            path,
            dpi=dpi,
            **_TIGHT_SAVEFIG_OPTIONS,
        )
    finally:
        plt.close(fig)


def save_png_plot_html(
    html_path: Path,
    png_path: Path,
    *,
    document_title: str,
) -> None:
    """Wrap a transparent PNG plot in a self-contained responsive HTML file."""

    html_path = Path(html_path)
    png_path = Path(png_path)
    if not png_path.is_file():
        raise FileNotFoundError(f"PNG plot does not exist: {png_path}")
    encoded_png = base64.b64encode(png_path.read_bytes()).decode("ascii")
    safe_title = escape(str(document_title))
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{safe_title}</title>
  <style>
    html, body {{
      margin: 0;
      min-height: 100%;
      background: transparent;
    }}
    main {{
      display: flex;
      justify-content: center;
      align-items: flex-start;
      width: 100%;
    }}
    img {{
      display: block;
      max-width: 100%;
      width: auto;
      height: auto;
    }}
  </style>
</head>
<body>
  <main>
    <img src="data:image/png;base64,{encoded_png}" alt="{safe_title}">
  </main>
</body>
</html>
"""
    html_path.write_text(html, encoding="utf-8", newline="\n")


def save_field_plot_svg(
    path: Path,
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    frequency_thz: float,
    quality_factor: float,
    quantity_label: str,
    unit_label: str,
    color_scale_mode: str,
) -> None:
    """Save a fully vector SVG field plot with editable text."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, _axis, _colorbar_axis = _create_field_figure(
            coordinates,
            values,
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label=quantity_label,
            unit_label=unit_label,
            color_scale_mode=color_scale_mode,
            rasterize_field=False,
        )
        try:
            fig.savefig(
                path,
                format="svg",
                **_TIGHT_SAVEFIG_OPTIONS,
            )
        finally:
            plt.close(fig)


def save_field_plot_hybrid_svg(
    path: Path,
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    frequency_thz: float,
    quality_factor: float,
    quantity_label: str,
    unit_label: str,
    color_scale_mode: str,
    dpi: int = 600,
) -> None:
    """Save an SVG with a rasterized field and vector annotations/colorbar."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, _axis, _colorbar_axis = _create_field_figure(
            coordinates,
            values,
            frequency_thz=frequency_thz,
            quality_factor=quality_factor,
            quantity_label=quantity_label,
            unit_label=unit_label,
            color_scale_mode=color_scale_mode,
            rasterize_field=True,
        )
        try:
            fig.savefig(
                path,
                format="svg",
                dpi=dpi,
                **_TIGHT_SAVEFIG_OPTIONS,
            )
        finally:
            plt.close(fig)
