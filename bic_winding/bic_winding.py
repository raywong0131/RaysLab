#!/usr/bin/env python3
"""Compute the winding number of a BIC polarization singularity.

The input table is expected to contain kx, ky, cx, cy. By default this script
uses Excel-style columns B, C, D, E and treats the file as having no header,
which matches the example workbook described by the user.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


# Optional manual parameters.
# Uncomment and edit these lines when you want fixed values without typing CLI flags.
# If a line is commented out, the program uses the default behavior.
#
# half_width = 0.05
# scale_x = 1.0
scale_y = 1.0


@dataclass(frozen=True)
class WindingResult:
    winding: float
    nearest_integer: int
    residual: float
    total_angle: float
    min_norm: float
    max_norm: float
    path_points: int
    field_mode: str
    path_mode: str
    center: tuple[float, float]


@dataclass(frozen=True)
class LoopScanPoint:
    half_width: float
    winding: float
    nearest_integer: int
    residual: float
    path_points: int
    min_norm: float


@dataclass(frozen=True)
class LoopScanResult:
    reference_integer: int
    stable_until: float | None
    changed_at: LoopScanPoint | None
    points: tuple[LoopScanPoint, ...]


def manual_parameter(name: str, default: float | None) -> float | None:
    value = globals().get(name, default)
    if value is None:
        return None
    return float(value)


def excel_col_to_index(col: str) -> int:
    text = col.strip()
    if text.isdigit():
        idx = int(text)
        if idx < 0:
            raise ValueError(f"Column index must be non-negative: {col!r}")
        return idx

    if not re.fullmatch(r"[A-Za-z]+", text):
        raise ValueError(f"Column must be an Excel letter or zero-based index: {col!r}")

    value = 0
    for ch in text.upper():
        value = value * 26 + (ord(ch) - ord("A") + 1)
    return value - 1


def parse_complex(value: object) -> complex:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return complex(np.nan, np.nan)
    if isinstance(value, complex):
        return value
    if isinstance(value, (int, float, np.integer, np.floating)):
        return complex(float(value), 0.0)

    text = str(value).strip()
    text = text.replace(" ", "")
    text = text.replace("I", "j").replace("i", "j")
    text = text.replace("−", "-")

    if text in {"", "nan", "NaN"}:
        return complex(np.nan, np.nan)
    if text in {"j", "+j"}:
        return 1j
    if text == "-j":
        return -1j

    try:
        return complex(text)
    except ValueError as exc:
        raise ValueError(f"Could not parse complex value {value!r}") from exc


def load_table(
    path: Path,
    sheet: str | int,
    header: str,
    kx_col: str,
    ky_col: str,
    cx_col: str,
    cy_col: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    read_header = None if header == "none" else 0
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls", ".xlsm"}:
        df = pd.read_excel(path, sheet_name=sheet, header=read_header)
    elif suffix in {".csv", ".txt"}:
        df = pd.read_csv(path, header=read_header)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t", header=read_header)
    else:
        raise ValueError(f"Unsupported input type: {path.suffix}")

    def get_col(selector: str) -> pd.Series:
        if header == "first" and selector in df.columns:
            return df[selector]
        return df.iloc[:, excel_col_to_index(selector)]

    kx = pd.to_numeric(get_col(kx_col), errors="coerce").to_numpy(float)
    ky = pd.to_numeric(get_col(ky_col), errors="coerce").to_numpy(float)
    cx = get_col(cx_col).map(parse_complex).to_numpy(np.complex128)
    cy = get_col(cy_col).map(parse_complex).to_numpy(np.complex128)

    good = np.isfinite(kx) & np.isfinite(ky) & np.isfinite(cx.real) & np.isfinite(cy.real)
    good &= np.isfinite(cx.imag) & np.isfinite(cy.imag)
    if not np.any(good):
        raise ValueError("No valid rows were found after parsing kx, ky, cx, cy.")

    return kx[good], ky[good], cx[good], cy[good]


def default_center(kx: np.ndarray, ky: np.ndarray) -> tuple[float, float]:
    return float((np.min(kx) + np.max(kx)) / 2), float((np.min(ky) + np.max(ky)) / 2)


def is_rectangular_grid(kx: np.ndarray, ky: np.ndarray) -> bool:
    ux = np.unique(kx)
    uy = np.unique(ky)
    return len(ux) * len(uy) == len(kx) and len(ux) > 1 and len(uy) > 1


def boundary_path_indices(kx: np.ndarray, ky: np.ndarray, tol: float) -> np.ndarray:
    xmin, xmax = float(np.min(kx)), float(np.max(kx))
    ymin, ymax = float(np.min(ky)), float(np.max(ky))

    pieces: list[int] = []

    bottom = np.where(np.abs(ky - ymin) <= tol)[0]
    pieces.extend(sorted(bottom, key=lambda i: kx[i]))

    right = np.where((np.abs(kx - xmax) <= tol) & (ky > ymin + tol))[0]
    pieces.extend(sorted(right, key=lambda i: ky[i]))

    top = np.where((np.abs(ky - ymax) <= tol) & (kx < xmax - tol))[0]
    pieces.extend(sorted(top, key=lambda i: -kx[i]))

    left = np.where(
        (np.abs(kx - xmin) <= tol) & (ky < ymax - tol) & (ky > ymin + tol)
    )[0]
    pieces.extend(sorted(left, key=lambda i: -ky[i]))

    if len(pieces) < 4:
        raise ValueError("Could not extract a rectangular boundary loop.")

    return np.array(pieces, dtype=int)


def angle_sorted_indices(
    kx: np.ndarray, ky: np.ndarray, center: tuple[float, float]
) -> np.ndarray:
    angles = np.arctan2(ky - center[1], kx - center[0])
    radii = np.hypot(kx - center[0], ky - center[1])
    return np.lexsort((radii, angles))


def input_path_indices(kx: np.ndarray, _ky: np.ndarray) -> np.ndarray:
    return np.arange(len(kx), dtype=int)


def select_path(
    kx: np.ndarray,
    ky: np.ndarray,
    path_mode: str,
    center: tuple[float, float],
    tol: float,
) -> tuple[np.ndarray, str]:
    if path_mode == "auto":
        if is_rectangular_grid(kx, ky):
            return boundary_path_indices(kx, ky, tol), "boundary"
        return angle_sorted_indices(kx, ky, center), "angle"
    if path_mode == "boundary":
        return boundary_path_indices(kx, ky, tol), "boundary"
    if path_mode == "angle":
        return angle_sorted_indices(kx, ky, center), "angle"
    if path_mode == "input":
        return input_path_indices(kx, ky), "input"
    raise ValueError(f"Unknown path mode: {path_mode}")


def close_path(
    kx: np.ndarray,
    ky: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    tol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    same_k = abs(kx[0] - kx[-1]) <= tol and abs(ky[0] - ky[-1]) <= tol
    same_c = abs(cx[0] - cx[-1]) <= tol and abs(cy[0] - cy[-1]) <= tol
    if same_k and same_c:
        return kx, ky, cx, cy
    return (
        np.r_[kx, kx[0]],
        np.r_[ky, ky[0]],
        np.r_[cx, cx[0]],
        np.r_[cy, cy[0]],
    )


def has_significant_imag(cx: np.ndarray, cy: np.ndarray, rel_tol: float = 1e-8) -> bool:
    scale = max(float(np.max(np.abs(cx))), float(np.max(np.abs(cy))), np.finfo(float).tiny)
    imag = max(float(np.max(np.abs(cx.imag))), float(np.max(np.abs(cy.imag))))
    return imag > rel_tol * scale


def complex_axis_vectors(cx: np.ndarray, cy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    field = np.column_stack([cx, cy])
    phase_seed = np.sum(field * field, axis=1)
    alpha = 0.5 * np.angle(phase_seed)
    vec = np.real(field * np.exp(-1j * alpha)[:, None])

    for i in range(1, len(vec)):
        if float(np.dot(vec[i - 1], vec[i])) < 0:
            vec[i] *= -1

    return vec[:, 0], vec[:, 1]


def field_angles(
    cx: np.ndarray,
    cy: np.ndarray,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, float, float]:
    resolved = mode
    if mode == "auto":
        resolved = "complex-axis" if has_significant_imag(cx, cy) else "real"

    if resolved == "real":
        vx, vy = cx.real, cy.real
        angles = np.unwrap(np.arctan2(vy, vx))
    elif resolved == "imag":
        vx, vy = cx.imag, cy.imag
        angles = np.unwrap(np.arctan2(vy, vx))
    elif resolved == "complex-axis":
        vx, vy = complex_axis_vectors(cx, cy)
        angles = np.unwrap(np.arctan2(vy, vx))
    elif resolved == "stokes-axis":
        vx = np.abs(cx) ** 2 - np.abs(cy) ** 2
        vy = 2 * np.real(cx * np.conj(cy))
        angles = 0.5 * np.unwrap(np.angle(vx + 1j * vy))
    else:
        raise ValueError(f"Unknown field mode: {mode}")

    norm = np.hypot(vx, vy)
    if np.any(norm == 0):
        raise ValueError("The selected loop contains a zero polarization vector.")

    return angles, vx, vy, resolved, float(np.min(norm)), float(np.max(norm))


def compute_winding(
    kx: np.ndarray,
    ky: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    path_mode: str,
    field_mode: str,
    center: tuple[float, float],
    tol: float,
    reverse: bool = False,
) -> WindingResult:
    indices, resolved_path = select_path(kx, ky, path_mode, center, tol)
    if reverse:
        indices = indices[::-1]

    pkx, pky = kx[indices], ky[indices]
    pcx, pcy = cx[indices], cy[indices]
    pkx, pky, pcx, pcy = close_path(pkx, pky, pcx, pcy, tol)

    angles, _vx, _vy, resolved_field, min_norm, max_norm = field_angles(pcx, pcy, field_mode)
    total_angle = float(angles[-1] - angles[0])
    winding = total_angle / (2 * math.pi)
    nearest = int(round(winding))

    return WindingResult(
        winding=winding,
        nearest_integer=nearest,
        residual=winding - nearest,
        total_angle=total_angle,
        min_norm=min_norm,
        max_norm=max_norm,
        path_points=len(pkx),
        field_mode=resolved_field,
        path_mode=resolved_path,
        center=center,
    )


def _contains_value(values: np.ndarray, target: float, tol: float) -> bool:
    return bool(np.any(np.abs(values - target) <= tol))


def symmetric_half_widths(values: np.ndarray, center_value: float, tol: float) -> list[float]:
    unique = np.unique(values)
    widths = sorted(
        {
            float(abs(value - center_value))
            for value in unique
            if abs(float(value - center_value)) > tol
        },
        reverse=True,
    )
    return [
        width
        for width in widths
        if _contains_value(unique, center_value - width, tol)
        and _contains_value(unique, center_value + width, tol)
    ]


def common_half_widths(kx: np.ndarray, ky: np.ndarray, center: tuple[float, float], tol: float) -> list[float]:
    x_widths = symmetric_half_widths(kx, center[0], tol)
    y_widths = symmetric_half_widths(ky, center[1], tol)

    common: list[float] = []
    for x_width in x_widths:
        match = next((y_width for y_width in y_widths if abs(y_width - x_width) <= tol), None)
        if match is not None:
            common.append(float((x_width + match) / 2))
    return common


def square_loop_mask(
    kx: np.ndarray,
    ky: np.ndarray,
    center: tuple[float, float],
    half_width: float,
    tol: float,
) -> np.ndarray:
    dx = kx - center[0]
    dy = ky - center[1]
    inside = (np.abs(dx) <= half_width + tol) & (np.abs(dy) <= half_width + tol)
    on_boundary = (np.abs(np.abs(dx) - half_width) <= tol) | (
        np.abs(np.abs(dy) - half_width) <= tol
    )
    return inside & on_boundary


def scan_shrinking_square_loops(
    kx: np.ndarray,
    ky: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    field_mode: str,
    center: tuple[float, float],
    tol: float,
    reverse: bool,
) -> LoopScanResult | None:
    if not is_rectangular_grid(kx, ky):
        return None

    points: list[LoopScanPoint] = []
    reference_integer: int | None = None
    stable_until: float | None = None
    changed_at: LoopScanPoint | None = None

    for half_width in common_half_widths(kx, ky, center, tol):
        mask = square_loop_mask(kx, ky, center, half_width, tol)
        if int(np.sum(mask)) < 4:
            continue

        try:
            result = compute_winding(
                kx[mask],
                ky[mask],
                cx[mask],
                cy[mask],
                "boundary",
                field_mode,
                center,
                tol,
                reverse,
            )
        except Exception:
            continue

        point = LoopScanPoint(
            half_width=half_width,
            winding=result.winding,
            nearest_integer=result.nearest_integer,
            residual=result.residual,
            path_points=result.path_points,
            min_norm=result.min_norm,
        )
        points.append(point)

        if reference_integer is None:
            reference_integer = point.nearest_integer
            stable_until = point.half_width
            continue

        if (
            point.nearest_integer != reference_integer
            or abs(point.winding - reference_integer) > LOOP_SCAN_INTEGER_TOL
        ):
            changed_at = point
            break

        stable_until = point.half_width

    if reference_integer is None:
        return None
    return LoopScanResult(reference_integer, stable_until, changed_at, tuple(points))


def threshold_plot_half_width(loop_scan: LoopScanResult | None) -> float | None:
    if loop_scan is None:
        return None
    if loop_scan.changed_at is not None:
        return loop_scan.stable_until
    if loop_scan.points:
        return loop_scan.points[-1].half_width
    return None


def select_calculation_data(
    kx: np.ndarray,
    ky: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    center: tuple[float, float],
    tol: float,
    half_width: float | None,
    path_mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    if half_width is None:
        return kx, ky, cx, cy, path_mode

    mask = square_loop_mask(kx, ky, center, half_width, tol)
    if int(np.sum(mask)) < 4:
        raise ValueError(f"half_width={half_width:.12g} does not select a valid square loop.")
    return kx[mask], ky[mask], cx[mask], cy[mask], "boundary"


def make_plot(
    output: Path,
    kx: np.ndarray,
    ky: np.ndarray,
    cx: np.ndarray,
    cy: np.ndarray,
    path_mode: str,
    field_mode: str,
    center: tuple[float, float],
    tol: float,
    reverse: bool,
    plot_half_width: float | None = None,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib, which is not available in this Python environment."
        ) from exc

    _all_angles, all_vx, all_vy, resolved_field, _all_min_norm, _all_max_norm = field_angles(
        cx,
        cy,
        field_mode,
    )
    all_norm = np.hypot(all_vx, all_vy)
    all_ux, all_uy = all_vx / all_norm, all_vy / all_norm

    loop_kx, loop_ky, loop_cx, loop_cy = kx, ky, cx, cy
    loop_path_mode = path_mode
    plot_title_extra = ""
    if plot_half_width is not None:
        mask = square_loop_mask(kx, ky, center, plot_half_width, tol)
        if int(np.sum(mask)) >= 4:
            loop_kx, loop_ky = kx[mask], ky[mask]
            loop_cx, loop_cy = cx[mask], cy[mask]
            loop_path_mode = "boundary"
            plot_title_extra = f", half_width={plot_half_width:.6g}"

    indices, _resolved_path = select_path(loop_kx, loop_ky, loop_path_mode, center, tol)
    if reverse:
        indices = indices[::-1]
    pkx, pky = loop_kx[indices], loop_ky[indices]
    pcx, pcy = loop_cx[indices], loop_cy[indices]
    pkx, pky, pcx, pcy = close_path(pkx, pky, pcx, pcy, tol)
    angles, vx, vy, _loop_field, _min_norm, _max_norm = field_angles(pcx, pcy, field_mode)

    norm = np.hypot(vx, vy)
    ux, uy = vx / norm, vy / norm

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    quiver_scale = 180
    axes[0].quiver(
        kx,
        ky,
        all_ux,
        all_uy,
        angles="xy",
        scale_units="xy",
        scale=quiver_scale,
        width=0.0018,
        color="0.55",
        alpha=0.6,
        pivot="middle",
    )
    axes[0].plot(pkx, pky, "-", color="tab:red", linewidth=1.0)
    axes[0].quiver(
        pkx,
        pky,
        ux,
        uy,
        angles="xy",
        scale_units="xy",
        scale=quiver_scale,
        width=0.0032,
        color="tab:blue",
        pivot="middle",
    )
    axes[0].scatter([center[0]], [center[1]], c="tab:red", s=24)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("kx")
    axes[0].set_ylabel("ky")
    axes[0].set_title(f"Loop and c-vector ({resolved_field}{plot_title_extra})")

    path_index = np.arange(len(angles))
    axes[1].plot(path_index, angles, color="tab:blue", label="unwrapped angle")

    if len(angles) > 1:
        trend = np.linspace(float(angles[0]), float(angles[-1]), len(angles))
        total_angle = float(angles[-1] - angles[0])
        winding = total_angle / (2 * math.pi)
        direction = "decreasing" if total_angle < 0 else "increasing"
        axes[1].plot(path_index, trend, "--", color="tab:red", linewidth=1.4, label="trend")

        start = max(0, len(angles) // 5)
        end = min(len(angles) - 1, start + max(1, len(angles) // 4))
        axes[1].annotate(
            "",
            xy=(path_index[end], trend[end]),
            xytext=(path_index[start], trend[start]),
            arrowprops={"arrowstyle": "->", "color": "tab:red", "lw": 1.8},
        )
        axes[1].text(
            0.03,
            0.05,
            f"{direction}\nDelta angle = {total_angle:.3g} rad\nwinding = {winding:.3g}",
            transform=axes[1].transAxes,
            va="bottom",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.8", "alpha": 0.9},
        )

    axes[1].set_xlabel("path index")
    axes[1].set_ylabel("unwrapped angle (rad)")
    axes[1].set_title("Accumulated polarization angle")
    axes[1].legend(loc="best")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


SUPPORTED_EXTS = {".xlsx", ".xls", ".xlsm", ".csv", ".tsv", ".txt"}
LOOP_SCAN_INTEGER_TOL = 0.25


def discover_data_files(script_dir: Path) -> list[Path]:
    """Return all supported data files in the same directory as the script."""
    found = sorted(
        p for p in script_dir.iterdir()
        if p.is_file() and p.suffix.lower() in SUPPORTED_EXTS
    )
    return found


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute a BIC winding number from kx, ky, cx, cy table data. "
            "If no input file is given, all supported data files in the script's "
            "directory are processed automatically."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        default=None,
        help="Input .xlsx, .xls, .csv, or .tsv file. Omit to auto-discover.",
    )
    parser.add_argument("--sheet", default=0, help="Excel sheet name or zero-based sheet index.")
    parser.add_argument("--header", choices=["none", "first"], default="none")
    parser.add_argument("--kx-col", default="B", help="kx column, default B.")
    parser.add_argument("--ky-col", default="C", help="ky column, default C.")
    parser.add_argument("--cx-col", default="D", help="cx column, default D.")
    parser.add_argument("--cy-col", default="E", help="cy column, default E.")
    parser.add_argument(
        "--half-width",
        type=float,
        default=manual_parameter("half_width", None),
        help="Use a square loop with this half width around the center.",
    )
    parser.add_argument(
        "--scale-x",
        "--cx-scale",
        dest="scale_x",
        type=float,
        default=manual_parameter("scale_x", 1.0),
        help="Multiply cx by this factor.",
    )
    parser.add_argument(
        "--scale-y",
        "--cy-scale",
        dest="scale_y",
        type=float,
        default=manual_parameter("scale_y", 1.0),
        help="Multiply cy by this factor.",
    )
    parser.add_argument(
        "--path",
        choices=["auto", "boundary", "angle", "input"],
        default="auto",
        help="Loop extraction mode. auto uses rectangular boundary for grid data.",
    )
    parser.add_argument(
        "--field-mode",
        choices=["auto", "real", "imag", "complex-axis", "stokes-axis"],
        default="auto",
        help="How cx,cy are converted to a polarization angle.",
    )
    parser.add_argument(
        "--center",
        nargs=2,
        type=float,
        metavar=("KX0", "KY0"),
        help="Center used for angle sorting and reporting. Defaults to data midpoint.",
    )
    parser.add_argument("--tol", type=float, default=1e-10, help="Coordinate tolerance.")
    parser.add_argument("--reverse", action="store_true", help="Reverse loop direction.")
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="Directory to save PNG diagnostic plots (one per file). Defaults to the script directory.",
    )
    parser.add_argument("--plot", type=Path, help="Output PNG plot (single-file mode only).")
    return parser.parse_args(argv)


def default_plot_path(input_path: Path, args: argparse.Namespace, script_dir: Path) -> Path:
    if args.plot is not None:
        return args.plot
    output_dir = args.plot_dir if args.plot_dir is not None else script_dir
    return output_dir / (input_path.stem + "_winding.png")


def print_loop_scan(scan: LoopScanResult | None) -> None:
    if scan is None:
        print("  shrink scan    : skipped (requires a rectangular grid)")
        return

    print("  shrink scan    :")
    print(f"    stable if |winding - reference| <= {LOOP_SCAN_INTEGER_TOL:g}")
    print("    half_width       winding   nearest   min |c|")
    for point in scan.points:
        print(
            f"    {point.half_width:>10.6g}  "
            f"{point.winding:>12.6g}  "
            f"{point.nearest_integer:>7d}  "
            f"{point.min_norm:.3e}"
        )

    if scan.changed_at is None:
        print(
            f"  threshold      : no change found; winding stays {scan.reference_integer} "
            f"down to half_width {scan.stable_until:.6g}"
        )
    else:
        stable = scan.stable_until
        changed = scan.changed_at
        if stable is None:
            print(
                f"  threshold      : first checked half_width {changed.half_width:.6g} "
                f"already gives {changed.nearest_integer}"
            )
        else:
            print(
                f"  threshold      : winding stays {scan.reference_integer} down to "
                f"half_width {stable:.6g}; changes at {changed.half_width:.6g} "
                f"to {changed.winding:.6g} (nearest {changed.nearest_integer})"
            )
            print(
                f"                   transition interval: "
                f"({changed.half_width:.6g}, {stable:.6g}]"
            )


def print_result(
    path: Path,
    kx: np.ndarray,
    result: WindingResult,
    plot_message: str | None,
    loop_scan: LoopScanResult | None,
    component_scale: tuple[float, float],
    half_width: float | None,
) -> None:
    print(f"\n{'=' * 60}")
    print(f"File             : {path.name}")
    print(f"  input rows     : {len(kx)}")
    print(f"  component scale: scale_x = {component_scale[0]:.12g}, scale_y = {component_scale[1]:.12g}")
    if half_width is None:
        print("  manual loop    : none")
    else:
        print(f"  manual loop    : half_width = {half_width:.12g}")
    print(f"  path mode      : {result.path_mode}")
    print(f"  path points    : {result.path_points}")
    print(f"  field mode     : {result.field_mode}")
    print(f"  center         : ({result.center[0]:.12g}, {result.center[1]:.12g})")
    print(f"  total angle    : {result.total_angle:.12g} rad")
    print(f"  winding number : {result.winding:.12g}")
    print(f"  nearest integer: {result.nearest_integer}")
    print(f"  residual       : {result.residual:.3e}")
    print(f"  min/max |c|    : {result.min_norm:.6e} / {result.max_norm:.6e}")
    print_loop_scan(loop_scan)
    if plot_message:
        print(f"  plot           : {plot_message}")


def process_file(
    path: Path,
    args: argparse.Namespace,
    sheet: str | int,
    plot_path: Path | None,
) -> bool:
    """Process a single file. Returns True on success, False on error."""
    try:
        kx, ky, cx, cy = load_table(
            path, sheet, args.header,
            args.kx_col, args.ky_col, args.cx_col, args.cy_col,
        )
    except Exception as exc:
        print(f"\n[SKIP] {path.name}: {exc}")
        return False

    cx = cx * args.scale_x
    cy = cy * args.scale_y
    center = tuple(args.center) if args.center is not None else default_center(kx, ky)

    try:
        calc_kx, calc_ky, calc_cx, calc_cy, calc_path = select_calculation_data(
            kx,
            ky,
            cx,
            cy,
            center,
            args.tol,
            args.half_width,
            args.path,
        )
        result = compute_winding(
            calc_kx,
            calc_ky,
            calc_cx,
            calc_cy,
            calc_path,
            args.field_mode,
            center,
            args.tol,
            args.reverse,
        )
    except Exception as exc:
        print(f"\n[ERROR] {path.name}: {exc}")
        return False

    loop_scan = scan_shrinking_square_loops(
        kx,
        ky,
        cx,
        cy,
        args.field_mode,
        center,
        args.tol,
        args.reverse,
    )

    plot_message = None
    if plot_path:
        try:
            make_plot(
                plot_path,
                kx,
                ky,
                cx,
                cy,
                args.path,
                args.field_mode,
                center,
                args.tol,
                args.reverse,
                args.half_width if args.half_width is not None else threshold_plot_half_width(loop_scan),
            )
            plot_message = str(plot_path)
        except RuntimeError as exc:
            plot_message = f"skipped ({exc})"

    print_result(path, kx, result, plot_message, loop_scan, (args.scale_x, args.scale_y), args.half_width)
    return True


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    sheet: str | int = args.sheet
    if isinstance(sheet, str) and sheet.isdigit():
        sheet = int(sheet)

    script_dir = Path(__file__).parent

    if args.input is not None:
        # Single-file mode
        plot_path = default_plot_path(args.input, args, script_dir)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        success = process_file(args.input, args, sheet, plot_path)
        return 0 if success else 1

    # Auto-discovery mode
    files = discover_data_files(script_dir)
    if not files:
        print(f"No data files found in {script_dir}")
        print(f"Supported extensions: {', '.join(sorted(SUPPORTED_EXTS))}")
        return 1

    print(f"Found {len(files)} file(s) in {script_dir}:")
    for f in files:
        print(f"  {f.name}")

    errors = 0
    for path in files:
        plot_path = default_plot_path(path, args, script_dir)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        ok = process_file(path, args, sheet, plot_path)
        if not ok:
            errors += 1

    print(f"\n{'=' * 60}")
    print(f"Done: {len(files) - errors}/{len(files)} file(s) processed successfully.")
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
