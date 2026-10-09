"""Offline polarization, color-map, and winding analysis for unit-cell scans."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import FuncFormatter, LogFormatterMathtext, MaxNLocator
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
from scipy.interpolate import (
    CloughTocher2DInterpolator,
    RectBivariateSpline,
)


ANGLE_METHODS = ("real", "imag", "complex-axis", "stokes-axis")
POLARIZATION_DEFINITION_VERSION = "unit-cell-2d-v1"
VALID_MATCH_STATUSES = {"gamma", "matched", "repaired"}
FORMAL_MODE_LABELS = {"p1": "px", "p2": "py"}
SYMMETRY_GAUGE_VERSION = "symmetry_gauge_v1"
SYMMETRY_RENDER_VERSION = "c2v-quarter-render-v1"
QUARTER_B4_COLORBAR_MINIMUM = 1e-4
QUARTER_B4_COLORBAR_MAXIMUM = 1e-1
QUARTER_B2_COLORBAR_MINIMUM = 1e-4
QUARTER_B2_COLORBAR_MAXIMUM = 1.0
QUARTER_K_TICK_INTERVAL = 0.05
POSITIVE_HANDEDNESS_COLOR = "#D62728"
NEGATIVE_HANDEDNESS_COLOR = "#1F77B4"
ZERO_HANDEDNESS_COLOR = "0.5"
REQUIRED_COLUMNS = {
    "qx_over_G",
    "qy_over_G",
    "band_label",
    "match_status",
    "cx_raw_re",
    "cx_raw_im",
    "cy_raw_re",
    "cy_raw_im",
    "cx_norm_re",
    "cx_norm_im",
    "cy_norm_re",
    "cy_norm_im",
}


def _format_decimal(value: float, places: int = 6) -> str:
    return f"{float(value):.{places}f}".rstrip("0").rstrip(".")


def unit_cell_2d_series_name(
    *,
    b0_nm: float,
    eta: float,
    zeta: float,
    mesh_size: int,
    q_max_over_g: float,
    q_points_per_axis: int,
    target_bands: Sequence[str],
    is_quarter: int = 0,
) -> str:
    if type(is_quarter) is not int or is_quarter not in {0, 1}:
        raise ValueError("is_quarter must be the integer 0 or 1")
    domain = "quarter" if is_quarter else "full"
    return (
        f"unit_cell_2D_({_format_decimal(b0_nm)}-"
        f"{_format_decimal(eta)}-{_format_decimal(zeta)})"
        f"_mesh{int(mesh_size)}_kmax{_format_decimal(q_max_over_g)}"
        f"_uniform{int(q_points_per_axis)}"
        f"_{domain}"
        f"_band-{'-'.join(formal_mode_label(item) for item in target_bands)}"
    )


def formal_mode_label(internal_label: str) -> str:
    return FORMAL_MODE_LABELS.get(str(internal_label), str(internal_label))


@dataclass(frozen=True)
class WindingResult:
    half_width: float
    winding: float
    nearest_integer: int
    residual: float
    total_angle: float
    min_norm: float
    path_points: int
    angle_method: str
    status: str = "ok"


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_polarization_grid(source: str | Path) -> pd.DataFrame:
    source = Path(source)
    if source.is_dir():
        logged = source / "80_logs" / "polarization_grid.csv"
        source = logged if logged.is_file() else source / "01_results" / "polarization_grid.csv"
    frame = pd.read_csv(source)
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Polarization grid is missing columns: {missing}")
    duplicated = frame.duplicated(["band_label", "qx_over_G", "qy_over_G"])
    if bool(duplicated.any()):
        raise ValueError("Polarization grid contains duplicate band/k rows")
    return derive_polarization_quantities(frame)


def complex_coefficients(
    frame: pd.DataFrame | PolarizationFieldView | PolarizationLoop,
    *,
    prefix: str = "cx_norm",
) -> tuple[np.ndarray, np.ndarray]:
    if prefix not in {"cx_norm", "cx_raw"}:
        raise ValueError(f"Unknown coefficient prefix: {prefix}")
    if isinstance(frame, PolarizationLoop):
        if prefix != "cx_norm":
            raise ValueError("PolarizationLoop only carries normalized coefficients")
        return frame.cx.copy(), frame.cy.copy()
    if isinstance(frame, PolarizationFieldView):
        cx, cy = frame.complex_grids(prefix=prefix)
        return cx.ravel(), cy.ravel()
    cy_prefix = "cy" + prefix.removeprefix("cx")
    cx = (
        pd.to_numeric(frame[f"{prefix}_re"], errors="coerce").to_numpy(float)
        + 1j * pd.to_numeric(frame[f"{prefix}_im"], errors="coerce").to_numpy(float)
    )
    cy = (
        pd.to_numeric(frame[f"{cy_prefix}_re"], errors="coerce").to_numpy(float)
        + 1j * pd.to_numeric(frame[f"{cy_prefix}_im"], errors="coerce").to_numpy(float)
    )
    return cx, cy


def derive_polarization_quantities(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    cx, cy = complex_coefficients(result)
    cx_raw, cy_raw = complex_coefficients(result, prefix="cx_raw")
    s0 = np.abs(cx) ** 2 + np.abs(cy) ** 2
    raw_s0 = np.abs(cx_raw) ** 2 + np.abs(cy_raw) ** 2
    s1 = np.abs(cx) ** 2 - np.abs(cy) ** 2
    s2 = 2.0 * np.real(cx * np.conj(cy))
    s3 = -2.0 * np.imag(cx * np.conj(cy))
    result["C"] = raw_s0
    result["C_norm"] = s0
    result["S0"] = s0
    result["S1"] = s1
    result["S2"] = s2
    result["S3"] = s3
    result["psi_stokes"] = 0.5 * np.angle(s1 + 1j * s2)
    finite = np.isfinite(cx.real) & np.isfinite(cx.imag)
    finite &= np.isfinite(cy.real) & np.isfinite(cy.imag)
    finite &= np.isfinite(cx_raw.real) & np.isfinite(cx_raw.imag)
    finite &= np.isfinite(cy_raw.real) & np.isfinite(cy_raw.imag)
    result["analysis_valid"] = (
        result["match_status"].isin(VALID_MATCH_STATUSES).to_numpy() & finite
    )
    return result


@dataclass(frozen=True)
class PolarizationLoop:
    """One transient closed-path source without materialized symmetry rows."""

    qx: np.ndarray
    qy: np.ndarray
    cx: np.ndarray
    cy: np.ndarray
    analysis_valid: np.ndarray

    def closed(self) -> "PolarizationLoop":
        return PolarizationLoop(
            qx=np.r_[self.qx, self.qx[0]],
            qy=np.r_[self.qy, self.qy[0]],
            cx=np.r_[self.cx, self.cx[0]],
            cy=np.r_[self.cy, self.cy[0]],
            analysis_valid=np.r_[self.analysis_valid, self.analysis_valid[0]],
        )

    def __len__(self) -> int:
        return len(self.qx)


class PolarizationFieldView:
    """Read one real COMSOL table through a full-frame coordinate interface."""

    is_quarter = False

    def __init__(self, source: pd.DataFrame, *, cy_scale: float = 1.0):
        self.source = source.copy()
        self.cy_scale = float(cy_scale)
        if not math.isfinite(self.cy_scale) or self.cy_scale <= 0.0:
            raise ValueError("cy_scale must be finite and positive")
        if self.source.empty:
            raise ValueError("Polarization field view requires at least one row")
        self._source_qx = np.sort(
            self.source["qx_over_G"].astype(float).unique()
        )
        self._source_qy = np.sort(
            self.source["qy_over_G"].astype(float).unique()
        )
        if (
            self._source_qx.size != self._source_qy.size
            or not np.allclose(self._source_qx, self._source_qy, atol=1e-12)
        ):
            raise ValueError("Polarization source must use one common square axis")
        if self._source_qx.size < 2:
            raise ValueError("Polarization source needs at least two points per axis")
        if not np.all(np.diff(self._source_qx) > 0.0):
            raise ValueError("Polarization source axis must be strictly increasing")
        for band_label, band in self.source.groupby("band_label", sort=False):
            expected = self._source_qx.size * self._source_qy.size
            if len(band) != expected:
                raise ValueError(
                    f"Band {band_label!r} is not a complete Cartesian source grid"
                )
            if band[["qx_over_G", "qy_over_G"]].duplicated().any():
                raise ValueError(f"Band {band_label!r} contains duplicate source points")
        if self.is_quarter:
            tolerance = 0.5 * 10 ** -6
            if self._source_qx[0] < -tolerance or not math.isclose(
                float(self._source_qx[0]), 0.0, abs_tol=tolerance
            ):
                raise ValueError("Quarter source axis must start at zero")

    @property
    def qx_axis(self) -> np.ndarray:
        if not self.is_quarter:
            return self._source_qx.copy()
        return np.r_[-self._source_qx[:0:-1], self._source_qx]

    @property
    def qy_axis(self) -> np.ndarray:
        if not self.is_quarter:
            return self._source_qy.copy()
        return np.r_[-self._source_qy[:0:-1], self._source_qy]

    @property
    def band_labels(self) -> list[str]:
        return [str(value) for value in self.source["band_label"].drop_duplicates()]

    @property
    def sampling_domain(self) -> str:
        return "quarter" if self.is_quarter else "full"

    def for_band(self, band_label: str) -> "PolarizationFieldView":
        band = self.source.loc[
            self.source["band_label"].astype(str) == str(band_label)
        ].copy()
        if band.empty:
            raise ValueError(f"Polarization grid contains no band {band_label!r}")
        return type(self)(band, cy_scale=self.cy_scale)

    def scaled_cy(self, factor: float) -> "PolarizationFieldView":
        return type(self)(self.source, cy_scale=self.cy_scale * float(factor))

    def _source_matrix(self, column: str) -> np.ndarray:
        if self.source["band_label"].nunique() != 1:
            raise ValueError("Select one band before requesting a grid matrix")
        table = self.source.pivot(
            index="qy_over_G", columns="qx_over_G", values=column
        ).sort_index(axis=0).sort_index(axis=1)
        if table.shape != (len(self._source_qy), len(self._source_qx)):
            raise ValueError(f"Column {column!r} does not form a complete source grid")
        return table.to_numpy()

    def _map_source_matrix(self, matrix: np.ndarray) -> np.ndarray:
        if not self.is_quarter:
            return np.asarray(matrix).copy()
        y_indices = np.r_[
            np.arange(len(self._source_qy) - 1, 0, -1),
            np.arange(len(self._source_qy)),
        ]
        x_indices = np.r_[
            np.arange(len(self._source_qx) - 1, 0, -1),
            np.arange(len(self._source_qx)),
        ]
        return np.asarray(matrix)[np.ix_(y_indices, x_indices)].copy()

    def complex_grids(
        self, *, prefix: str = "cx_norm"
    ) -> tuple[np.ndarray, np.ndarray]:
        if prefix not in {"cx_norm", "cx_raw"}:
            raise ValueError(f"Unknown coefficient prefix: {prefix}")
        suffix = prefix.removeprefix("cx")
        cx_source = self._source_matrix(f"cx{suffix}_re").astype(float)
        cx_source = cx_source + 1j * self._source_matrix(
            f"cx{suffix}_im"
        ).astype(float)
        cy_source = self._source_matrix(f"cy{suffix}_re").astype(float)
        cy_source = cy_source + 1j * self._source_matrix(
            f"cy{suffix}_im"
        ).astype(float)
        cx = self._map_source_matrix(cx_source)
        cy = self._map_source_matrix(cy_source)
        if self.is_quarter:
            cx[:, self.qx_axis < 0.0] *= -1.0
            cy[self.qy_axis < 0.0, :] *= -1.0
        if prefix == "cx_norm":
            cy *= self.cy_scale
        return cx, cy

    def grid_matrix(self, column: str) -> np.ndarray:
        coefficient_columns = {
            f"{component}_{kind}_{part}"
            for component in ("cx", "cy")
            for kind in ("raw", "norm")
            for part in ("re", "im")
        }
        if column in coefficient_columns:
            component, kind, part = column.split("_")
            cx, cy = self.complex_grids(prefix=f"cx_{kind}")
            selected = cx if component == "cx" else cy
            return selected.real if part == "re" else selected.imag
        if column in {"C", "C_norm", "S0", "S1", "S2", "S3", "psi_stokes"}:
            cx, cy = self.complex_grids(
                prefix="cx_raw" if column == "C" else "cx_norm"
            )
            s0 = np.abs(cx) ** 2 + np.abs(cy) ** 2
            if column in {"C", "C_norm", "S0"}:
                return s0
            s1 = np.abs(cx) ** 2 - np.abs(cy) ** 2
            if column == "S1":
                return s1
            s2 = 2.0 * np.real(cx * np.conj(cy))
            if column == "S2":
                return s2
            if column == "S3":
                return -2.0 * np.imag(cx * np.conj(cy))
            angles = 0.5 * np.angle(s1 + 1j * s2)
            if self.is_quarter:
                # A single mirror reverses the director angle at interior
                # points.  Points on the invariant coordinate axes require
                # the axis branch itself, however: applying the interior sign
                # rule there creates an artificial phase seam that is absent
                # from an independently solved full-space grid.
                origin_x = int(np.flatnonzero(np.isclose(self.qx_axis, 0.0))[0])
                origin_y = int(np.flatnonzero(np.isclose(self.qy_axis, 0.0))[0])
                angles[:origin_y, origin_x] = angles[:origin_y:-1, origin_x]
                angles[origin_y, :origin_x] = angles[origin_y, :origin_x:-1]
            return angles
        return self._map_source_matrix(self._source_matrix(column))

    def flat_samples(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        mesh_x, mesh_y = np.meshgrid(self.qx_axis, self.qy_axis)
        cx, cy = self.complex_grids()
        valid = self.grid_matrix("analysis_valid").astype(bool)
        return (
            mesh_x.ravel(),
            mesh_y.ravel(),
            cx.ravel(),
            cy.ravel(),
            valid.ravel(),
        )

    def ordered_square_loop(self, half_width: float) -> PolarizationLoop:
        qx_axis = self.qx_axis
        qy_axis = self.qy_axis
        tolerance = max(abs(float(half_width)) * 1e-9, 1e-13)
        x_candidates = np.flatnonzero(np.isclose(qx_axis, half_width, atol=tolerance))
        y_candidates = np.flatnonzero(np.isclose(qy_axis, half_width, atol=tolerance))
        negative_x = np.flatnonzero(np.isclose(qx_axis, -half_width, atol=tolerance))
        negative_y = np.flatnonzero(np.isclose(qy_axis, -half_width, atol=tolerance))
        if not all(len(item) == 1 for item in (x_candidates, y_candidates, negative_x, negative_y)):
            raise ValueError(f"half_width={half_width:g} is not on the plot grid")
        ix0, ix1 = int(negative_x[0]), int(x_candidates[0])
        iy0, iy1 = int(negative_y[0]), int(y_candidates[0])
        coordinates = []
        coordinates.extend((ix, iy0) for ix in range(ix0, ix1 + 1))
        coordinates.extend((ix1, iy) for iy in range(iy0 + 1, iy1 + 1))
        coordinates.extend((ix, iy1) for ix in range(ix1 - 1, ix0 - 1, -1))
        coordinates.extend((ix0, iy) for iy in range(iy1 - 1, iy0, -1))
        cx, cy = self.complex_grids()
        valid = self.grid_matrix("analysis_valid").astype(bool)
        return PolarizationLoop(
            qx=np.asarray([qx_axis[ix] for ix, _iy in coordinates], dtype=float),
            qy=np.asarray([qy_axis[iy] for _ix, iy in coordinates], dtype=float),
            cx=np.asarray([cx[iy, ix] for ix, iy in coordinates], dtype=complex),
            cy=np.asarray([cy[iy, ix] for ix, iy in coordinates], dtype=complex),
            analysis_valid=np.asarray(
                [valid[iy, ix] for ix, iy in coordinates], dtype=bool
            ),
        )


class C2vQuarterFieldView(PolarizationFieldView):
    """C2v view whose source rows remain restricted to the first quadrant."""

    is_quarter = True


class FullFieldView(PolarizationFieldView):
    """Direct view of a legacy or explicit full-space COMSOL table."""

    is_quarter = False


def _series_is_quarter(series_dir: Path) -> int:
    config_path = Path(series_dir) / "99_config" / "config.json"
    if not config_path.is_file():
        return 0
    config = json.loads(config_path.read_text(encoding="utf-8"))
    identity = config.get("cache_identity", {})
    unit_cell = identity.get("unit_cell_2d", {}) if isinstance(identity, dict) else {}
    value = identity.get("isQuarter", unit_cell.get("isQuarter", 0))
    if type(value) is not int or value not in {0, 1}:
        raise ValueError(f"Invalid isQuarter in {config_path}: {value!r}")
    return value


def load_polarization_field_view(source: str | Path) -> PolarizationFieldView:
    source_path = Path(source)
    frame = load_polarization_grid(source_path)
    is_quarter = _series_is_quarter(source_path) if source_path.is_dir() else 0
    view_type = C2vQuarterFieldView if is_quarter else FullFieldView
    return view_type(frame)


def _coerce_field_view(
    source: pd.DataFrame | PolarizationFieldView,
) -> PolarizationFieldView:
    if isinstance(source, PolarizationFieldView):
        return source
    return FullFieldView(derive_polarization_quantities(source))


def complex_axis_vectors(
    cx: np.ndarray,
    cy: np.ndarray,
    *,
    continuous: bool,
) -> tuple[np.ndarray, np.ndarray]:
    field = np.column_stack([cx, cy])
    phase = 0.5 * np.angle(np.sum(field * field, axis=1))
    vectors = np.real(field * np.exp(-1j * phase)[:, None])
    if continuous:
        for index in range(1, len(vectors)):
            if float(np.dot(vectors[index - 1], vectors[index])) < 0.0:
                vectors[index] *= -1.0
    return vectors[:, 0], vectors[:, 1]


def field_angles(
    cx: np.ndarray,
    cy: np.ndarray,
    method: str,
    *,
    continuous: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if method == "real":
        vx, vy = cx.real, cy.real
        angles = np.unwrap(np.arctan2(vy, vx))
    elif method == "imag":
        vx, vy = cx.imag, cy.imag
        angles = np.unwrap(np.arctan2(vy, vx))
    elif method == "complex-axis":
        vx, vy = complex_axis_vectors(cx, cy, continuous=continuous)
        angles = np.unwrap(np.arctan2(vy, vx))
    elif method == "stokes-axis":
        vx = np.abs(cx) ** 2 - np.abs(cy) ** 2
        vy = 2.0 * np.real(cx * np.conj(cy))
        angles = 0.5 * np.unwrap(np.angle(vx + 1j * vy))
    else:
        raise ValueError(f"Unknown angle method: {method}")
    return angles, np.asarray(vx, float), np.asarray(vy, float)


def common_half_widths(
    frame: pd.DataFrame | PolarizationFieldView,
) -> list[float]:
    if isinstance(frame, PolarizationFieldView):
        qx, qy = frame.qx_axis, frame.qy_axis
    else:
        qx = np.sort(frame["qx_over_G"].unique().astype(float))
        qy = np.sort(frame["qy_over_G"].unique().astype(float))
    x_widths = {round(abs(value), 14) for value in qx if value > 0}
    y_widths = {round(abs(value), 14) for value in qy if value > 0}
    return sorted(x_widths.intersection(y_widths), reverse=True)


def ordered_square_loop(
    frame: pd.DataFrame | PolarizationFieldView,
    half_width: float,
) -> pd.DataFrame | PolarizationLoop:
    if isinstance(frame, PolarizationFieldView):
        return frame.ordered_square_loop(half_width)
    qx = pd.to_numeric(frame["qx_over_G"], errors="coerce").to_numpy(float)
    qy = pd.to_numeric(frame["qy_over_G"], errors="coerce").to_numpy(float)
    tolerance = max(abs(float(half_width)) * 1e-9, 1e-13)
    inside = (np.abs(qx) <= half_width + tolerance) & (
        np.abs(qy) <= half_width + tolerance
    )
    boundary = np.isclose(np.abs(qx), half_width, atol=tolerance) | np.isclose(
        np.abs(qy), half_width, atol=tolerance
    )
    selected = frame.loc[inside & boundary].copy()
    if len(selected) < 4:
        raise ValueError(f"half_width={half_width:g} has fewer than four points")
    bottom = selected[np.isclose(selected["qy_over_G"], -half_width, atol=tolerance)]
    right = selected[
        np.isclose(selected["qx_over_G"], half_width, atol=tolerance)
        & (selected["qy_over_G"] > -half_width + tolerance)
    ]
    top = selected[
        np.isclose(selected["qy_over_G"], half_width, atol=tolerance)
        & (selected["qx_over_G"] < half_width - tolerance)
    ]
    left = selected[
        np.isclose(selected["qx_over_G"], -half_width, atol=tolerance)
        & (selected["qy_over_G"] < half_width - tolerance)
        & (selected["qy_over_G"] > -half_width + tolerance)
    ]
    pieces = [
        bottom.sort_values("qx_over_G"),
        right.sort_values("qy_over_G"),
        top.sort_values("qx_over_G", ascending=False),
        left.sort_values("qy_over_G", ascending=False),
    ]
    loop = pd.concat(pieces, ignore_index=True)
    if len(loop) != len(selected):
        raise ValueError(f"half_width={half_width:g} does not form one square loop")
    return loop


def compute_winding(
    loop: pd.DataFrame | PolarizationLoop,
    half_width: float,
    angle_method: str,
) -> WindingResult:
    loop_valid = (
        loop.analysis_valid
        if isinstance(loop, PolarizationLoop)
        else loop["analysis_valid"].to_numpy(bool)
    )
    if not bool(np.all(loop_valid)):
        raise ValueError("Winding loop contains invalid tracked points")
    cx, cy = complex_coefficients(loop)
    cx = np.r_[cx, cx[0]]
    cy = np.r_[cy, cy[0]]
    angles, vx, vy = field_angles(cx, cy, angle_method, continuous=True)
    norm = np.hypot(vx, vy)
    if not bool(np.all(np.isfinite(norm))) or bool(np.any(norm <= 0.0)):
        raise ValueError("Winding loop contains a zero or invalid polarization vector")
    total = float(angles[-1] - angles[0])
    winding = total / (2.0 * math.pi)
    nearest = int(round(winding))
    return WindingResult(
        half_width=float(half_width),
        winding=winding,
        nearest_integer=nearest,
        residual=winding - nearest,
        total_angle=total,
        min_norm=float(np.min(norm)),
        path_points=len(loop) + 1,
        angle_method=angle_method,
    )


def scan_shrinking_loops(
    frame: pd.DataFrame | PolarizationFieldView,
    angle_method: str,
) -> tuple[pd.DataFrame, WindingResult | None]:
    rows: list[dict[str, object]] = []
    stable: WindingResult | None = None
    reference_integer: int | None = None
    stable_candidates: list[tuple[int, WindingResult]] = []
    for half_width in common_half_widths(frame):
        try:
            result = compute_winding(
                ordered_square_loop(frame, half_width),
                half_width,
                angle_method,
            )
            row = asdict(result)
            if reference_integer is None:
                reference_integer = result.nearest_integer
            if (
                result.nearest_integer == reference_integer
                and abs(result.residual) <= 0.25
            ):
                stable_candidates.append((len(rows), result))
            else:
                row["status"] = "unstable"
        except (ValueError, FloatingPointError) as exc:
            row = {
                "half_width": half_width,
                "angle_method": angle_method,
                "status": "invalid",
                "error": str(exc),
            }
        rows.append(row)
    if len(rows) == 1 and len(stable_candidates) == 1:
        stable = stable_candidates[0][1]
    else:
        for (left_index, left), (right_index, _right) in zip(
            stable_candidates, stable_candidates[1:]
        ):
            if right_index == left_index + 1:
                stable = left
                break
    return pd.DataFrame(rows), stable


def _cell_edges(values: np.ndarray) -> np.ndarray:
    values = np.sort(np.asarray(values, float))
    if len(values) < 2:
        raise ValueError("A color map requires at least two coordinates per axis")
    middle = 0.5 * (values[:-1] + values[1:])
    return np.r_[
        values[0] - 0.5 * (values[1] - values[0]),
        middle,
        values[-1] + 0.5 * (values[-1] - values[-2]),
    ]


def _band_frame(frame: pd.DataFrame, band_label: str) -> pd.DataFrame:
    band = frame[frame["band_label"] == band_label].copy()
    if band.empty:
        raise ValueError(f"Polarization grid contains no band {band_label!r}")
    x_count = band["qx_over_G"].nunique()
    y_count = band["qy_over_G"].nunique()
    if x_count * y_count != len(band):
        raise ValueError(f"Band {band_label!r} is not a complete Cartesian grid")
    return band


def _grid_matrix(
    frame: pd.DataFrame | PolarizationFieldView,
    column: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if isinstance(frame, PolarizationFieldView):
        return frame.qx_axis, frame.qy_axis, frame.grid_matrix(column)
    table = frame.pivot(index="qy_over_G", columns="qx_over_G", values=column)
    table = table.sort_index(axis=0).sort_index(axis=1)
    return (
        table.columns.to_numpy(float),
        table.index.to_numpy(float),
        table.to_numpy(float),
    )


def _formal_mode_label(internal_label: str) -> str:
    return formal_mode_label(internal_label)


def _smooth_grid(
    qx: np.ndarray,
    qy: np.ndarray,
    values: np.ma.MaskedArray,
    *,
    subdivisions_per_interval: int = 24,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Densify a complete nonuniform Cartesian grid for smooth rendering."""
    if np.ma.getmaskarray(values).any():
        raise ValueError("Smooth color maps require a complete valid Cartesian grid")
    source = np.asarray(values, dtype=float)
    dense_qx = np.concatenate(
        [
            np.linspace(left, right, subdivisions_per_interval, endpoint=False)
            for left, right in zip(qx[:-1], qx[1:])
        ]
        + [np.asarray([qx[-1]], dtype=float)]
    )
    dense_qy = np.concatenate(
        [
            np.linspace(left, right, subdivisions_per_interval, endpoint=False)
            for left, right in zip(qy[:-1], qy[1:])
        ]
        + [np.asarray([qy[-1]], dtype=float)]
    )
    spline = RectBivariateSpline(
        np.asarray(qy, dtype=float),
        np.asarray(qx, dtype=float),
        source,
        kx=min(3, len(qy) - 1),
        ky=min(3, len(qx) - 1),
        s=0.0,
    )
    dense = spline(dense_qy, dense_qx)
    dense = np.clip(dense, float(np.min(source)), float(np.max(source)))
    return dense_qx, dense_qy, dense


def _valid_intensity_mask(
    frame: pd.DataFrame,
    relative_threshold: float,
) -> np.ndarray:
    intensity = pd.to_numeric(frame["C"], errors="coerce").to_numpy(float)
    valid = frame["analysis_valid"].to_numpy(bool) & np.isfinite(intensity)
    if not bool(np.any(valid)):
        return valid
    threshold = float(np.nanmax(intensity[valid])) * float(relative_threshold)
    return valid & (intensity > threshold)


def _smooth_cyclic_angles(
    frame: pd.DataFrame | PolarizationFieldView,
    relative_threshold: float,
    *,
    subdivisions_per_interval: int = 24,
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    """Smooth a pi-periodic axis through its doubled-angle unit vector."""
    field = _coerce_field_view(frame)
    qx, qy, intensity = _grid_matrix(field, "C")
    _qx, _qy, tracked = _grid_matrix(field, "analysis_valid")
    _qx, _qy, angles = _grid_matrix(field, "psi_stokes")
    finite_intensity = np.isfinite(intensity) & tracked.astype(bool)
    if not np.any(finite_intensity):
        raise ValueError("Cyclic phase smoothing has no valid intensity")
    threshold = float(np.max(intensity[finite_intensity])) * float(
        relative_threshold
    )
    valid = finite_intensity & (intensity > threshold) & np.isfinite(angles)
    if int(np.count_nonzero(valid)) < 4:
        raise ValueError("Cyclic phase smoothing requires at least four valid points")
    mesh_x, mesh_y = np.meshgrid(qx, qy)
    points = np.column_stack([mesh_x[valid], mesh_y[valid]])
    doubled_x = np.cos(2.0 * angles[valid])
    doubled_y = np.sin(2.0 * angles[valid])
    dense_qx = np.concatenate(
        [
            np.linspace(left, right, subdivisions_per_interval, endpoint=False)
            for left, right in zip(qx[:-1], qx[1:])
        ]
        + [np.asarray([qx[-1]], dtype=float)]
    )
    dense_qy = np.concatenate(
        [
            np.linspace(left, right, subdivisions_per_interval, endpoint=False)
            for left, right in zip(qy[:-1], qy[1:])
        ]
        + [np.asarray([qy[-1]], dtype=float)]
    )
    mesh_x, mesh_y = np.meshgrid(dense_qx, dense_qy)
    smooth_x = CloughTocher2DInterpolator(points, doubled_x)(mesh_x, mesh_y)
    smooth_y = CloughTocher2DInterpolator(points, doubled_y)(mesh_x, mesh_y)
    smooth_norm = np.hypot(smooth_x, smooth_y)
    invalid = ~np.isfinite(smooth_norm) | (smooth_norm <= 1e-12)
    smooth_angles = 0.5 * np.arctan2(smooth_y, smooth_x)
    return dense_qx, dense_qy, np.ma.masked_where(invalid, smooth_angles)


def _aligned_colorbar(figure, axis, image, **kwargs):
    """Create a colorbar axis whose top and bottom exactly match the map axis."""
    divider = make_axes_locatable(axis)
    colorbar_axis = divider.append_axes("right", size="4.5%", pad=0.10)
    colorbar = figure.colorbar(image, cax=colorbar_axis, **kwargs)
    colorbar._layout_divider = divider
    return colorbar


def _color_map_figure():
    """Reserve fixed room for an equal-height colorbar and its vertical label."""
    figure, axis = plt.subplots(figsize=(5.2, 4.4))
    figure.subplots_adjust(left=0.15, right=0.78, bottom=0.15, top=0.88)
    return figure, axis


def _center_color_map_group(figure, axis, colorbar) -> None:
    """Center the fully decorated map/colorbar group on the canvas."""
    for _iteration in range(2):
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        decorated = Bbox.union(
            [
                axis.get_tightbbox(renderer),
                colorbar.ax.get_tightbbox(renderer),
            ]
        )
        shift = (
            0.5 * figure.bbox.width
            - 0.5 * (decorated.x0 + decorated.x1)
        ) / figure.bbox.width
        divider = colorbar._layout_divider
        left, bottom, width, height = divider.get_position()
        divider.set_position((left + shift, bottom, width, height))
    figure.canvas.draw()


def _save_color_map(figure, axis, colorbar, path: Path) -> None:
    _center_color_map_group(figure, axis, colorbar)
    figure.savefig(path, dpi=300)


def _set_component_colorbar_ticks(
    colorbar,
    norm,
    *,
    maximum_as_power: bool = False,
) -> None:
    """Include the exact displayed maximum while keeping readable log decades."""
    maximum = float(norm.vmax)
    if isinstance(norm, LogNorm):
        first_power = int(math.ceil(math.log10(float(norm.vmin))))
        last_power = int(math.floor(math.log10(maximum)))
        decades = np.power(
            10.0, np.arange(first_power, last_power + 1, dtype=int)
        )
        decades = decades[decades < maximum * (1.0 - 1e-12)]
        ticks = np.r_[decades, maximum]
        labels = [
            rf"$10^{{{int(round(math.log10(value)))}}}$"
            for value in decades
        ] + [
            rf"$10^{{{int(round(math.log10(maximum)))}}}$"
            if maximum_as_power
            else f"{maximum:.4g}"
        ]
        while len(ticks) > 1:
            colorbar.set_ticks(ticks)
            colorbar.ax.set_yticklabels(labels)
            colorbar.ax.figure.canvas.draw()
            renderer = colorbar.ax.figure.canvas.get_renderer()
            tick_labels = [
                label
                for label in colorbar.ax.get_yticklabels()
                if label.get_visible() and label.get_text()
            ]
            if len(tick_labels) < 2:
                break
            lower_box = tick_labels[-2].get_window_extent(renderer=renderer)
            maximum_box = tick_labels[-1].get_window_extent(renderer=renderer)
            if lower_box.y1 + 2.0 <= maximum_box.y0:
                break
            ticks = np.delete(ticks, -2)
            del labels[-2]
    else:
        ticks = np.linspace(float(norm.vmin), maximum, 6)
        labels = [f"{value:.4g}" for value in ticks]
    colorbar.set_ticks(ticks)
    colorbar.ax.set_yticklabels(labels)


def _set_k_space_ticks(
    axis,
    qx: np.ndarray,
    qy: np.ndarray,
    *,
    is_quarter: bool,
) -> None:
    if not is_quarter:
        axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        return

    def ticks(values: np.ndarray) -> np.ndarray:
        lower = float(np.min(values))
        upper = float(np.max(values))
        first = math.ceil((lower - 1e-12) / QUARTER_K_TICK_INTERVAL)
        last = math.floor((upper + 1e-12) / QUARTER_K_TICK_INTERVAL)
        result = QUARTER_K_TICK_INTERVAL * np.arange(first, last + 1)
        result = np.round(result, 12)
        if result.size == 0 or not math.isclose(
            float(result[0]), lower, abs_tol=1e-12
        ):
            result = np.r_[lower, result]
        if not math.isclose(float(result[-1]), upper, abs_tol=1e-12):
            result = np.r_[result, upper]
        return result

    axis.set_xticks(ticks(np.asarray(qx, dtype=float)))
    axis.set_yticks(ticks(np.asarray(qy, dtype=float)))
    concise = FuncFormatter(
        lambda value, _position: (
            "0"
            if math.isclose(float(value), 0.0, abs_tol=1e-12)
            else f"{float(value):.2f}".rstrip("0").rstrip(".")
        )
    )
    axis.xaxis.set_major_formatter(concise)
    axis.yaxis.set_major_formatter(concise)


def plot_intensity_maps(
    frame: pd.DataFrame | PolarizationFieldView,
    output_dir: Path,
    band_label: str,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = _coerce_field_view(frame)
    qx, qy, intensity = _grid_matrix(frame, "C")
    _qx_valid, _qy_valid, valid = _grid_matrix(frame, "analysis_valid")
    intensity = np.ma.masked_where(~valid.astype(bool), intensity)
    if frame.is_quarter:
        intensity = np.ma.sqrt(np.ma.maximum(intensity, 0.0))
    compressed = intensity.compressed()
    if compressed.size == 0 or not bool(np.all(np.isfinite(compressed))):
        raise ValueError(f"Band {band_label!r} has no finite tracked intensity")
    maximum = float(np.max(compressed))
    positive = compressed[compressed > 0.0]
    if maximum <= 0.0 or positive.size == 0:
        raise ValueError(f"Band {band_label!r} has no positive tracked intensity")
    intensity = intensity / maximum
    plot_qx, plot_qy, plot_intensity = _smooth_grid(qx, qy, intensity)
    positive = intensity.compressed()[intensity.compressed() > 0.0]
    log_minimum = 10.0 ** math.floor(math.log10(float(np.min(positive))))
    log_maximum = 1.0
    if math.isclose(log_minimum, log_maximum, rel_tol=1e-12, abs_tol=0.0):
        log_minimum /= 10.0
    x_edges, y_edges = _cell_edges(qx), _cell_edges(qy)
    outputs: list[Path] = []
    specifications = (
        ("A1_ck_intensity_linear", Normalize(vmin=0.0, vmax=1.0)),
        (
            "A2_ck_intensity_log",
            LogNorm(
                vmin=log_minimum,
                vmax=log_maximum,
            ),
        ),
    )
    for stem, norm in specifications:
        figure, axis = _color_map_figure()
        image = axis.pcolormesh(
            _cell_edges(plot_qx),
            _cell_edges(plot_qy),
            plot_intensity,
            cmap="magma",
            norm=norm,
            shading="flat",
        )
        axis.set_aspect("equal", adjustable="box")
        _set_k_space_ticks(
            axis, qx, qy, is_quarter=frame.is_quarter
        )
        axis.set_xlabel(r"$k_x/G$")
        axis.set_ylabel(r"$k_y/G$")
        if frame.is_quarter:
            title = rf"Normalized $|c(k)|$ ({_formal_mode_label(band_label)})"
            colorbar_label = r"Normalized $|c(k)|$"
        else:
            title = (
                rf"Normalized $|c_x|^2+|c_y|^2$ "
                rf"({_formal_mode_label(band_label)})"
            )
            colorbar_label = "Normalized intensity"
        axis.set_title(title)
        colorbar = _aligned_colorbar(figure, axis, image)
        if isinstance(norm, LogNorm):
            first_power = int(math.ceil(math.log10(float(norm.vmin))))
            last_power = int(math.floor(math.log10(float(norm.vmax))))
            powers = np.arange(first_power, last_power + 1, dtype=int)
            if powers.size == 0:
                powers = np.asarray([int(round(math.log10(float(norm.vmax))))])
            colorbar.set_ticks(np.power(10.0, powers))
            colorbar.formatter = LogFormatterMathtext(
                base=10.0, labelOnlyBase=True
            )
            colorbar.minorticks_off()
            colorbar.update_ticks()
        colorbar.set_label(colorbar_label, rotation=270, labelpad=24)
        path = output_dir / f"{stem}.png"
        _save_color_map(figure, axis, colorbar, path)
        plt.close(figure)
        outputs.append(path)
    return outputs


def plot_component_magnitude_maps(
    frame: pd.DataFrame | PolarizationFieldView,
    output_dir: Path,
    band_label: str,
    *,
    panels: set[str] | None = None,
) -> list[Path]:
    """Plot selected B panels through one shared data and rendering contract."""
    output_dir.mkdir(parents=True, exist_ok=True)
    field = _coerce_field_view(frame)
    qx, qy = field.qx_axis, field.qy_axis
    cx, cy = field.complex_grids(prefix="cx_raw")
    valid = field.grid_matrix("analysis_valid").astype(bool)
    vector_magnitude = np.hypot(np.abs(cx), np.abs(cy))
    common_maximum = float(np.max(vector_magnitude[valid]))
    if common_maximum <= 0.0:
        raise ValueError(f"Band {band_label!r} has no positive |c(k)|")
    outputs: list[Path] = []
    specifications = (
        ("B1", "B2", "cx", np.abs(cx), "c_x"),
        ("B3", "B4", "cy", np.abs(cy), "c_y"),
    )
    for linear_panel, log_panel, component, component_values, symbol in specifications:
        values = component_values
        values = np.ma.masked_where(~valid.astype(bool), values)
        component_maximum = float(np.max(values.compressed()))
        if component_maximum <= 0.0:
            raise ValueError(f"Band {band_label!r} has no positive {component}")
        relative_values = values / component_maximum
        values = values / common_maximum
        colorbar_maximum = component_maximum / common_maximum
        plot_qx, plot_qy, plot_values = _smooth_grid(qx, qy, values)
        relative_positive = relative_values.compressed()[
            relative_values.compressed() > 0.0
        ]
        relative_floor = 10.0 ** math.floor(
            math.log10(float(np.min(relative_positive)))
        )
        log_floor = relative_floor * colorbar_maximum
        log_norm = LogNorm(log_floor, colorbar_maximum)
        if field.is_quarter:
            if log_panel == "B2":
                log_norm = LogNorm(
                    QUARTER_B2_COLORBAR_MINIMUM,
                    QUARTER_B2_COLORBAR_MAXIMUM,
                )
            elif log_panel == "B4":
                log_norm = LogNorm(
                    QUARTER_B4_COLORBAR_MINIMUM,
                    QUARTER_B4_COLORBAR_MAXIMUM,
                )
        linear_norm = Normalize(0.0, colorbar_maximum)
        if field.is_quarter and linear_panel == "B1":
            linear_norm = Normalize(0.0, 1.0)
        for panel, scale, norm in (
            (linear_panel, "linear", linear_norm),
            (log_panel, "log", log_norm),
        ):
            if panels is not None and panel not in panels:
                continue
            figure, axis = _color_map_figure()
            image = axis.pcolormesh(
                _cell_edges(plot_qx), _cell_edges(plot_qy), plot_values,
                cmap="magma", norm=norm, shading="flat",
            )
            axis.set_aspect("equal", adjustable="box")
            _set_k_space_ticks(
                axis, qx, qy, is_quarter=field.is_quarter
            )
            axis.set_xlabel(r"$k_x/G$")
            axis.set_ylabel(r"$k_y/G$")
            axis.set_title(
                rf"Normalized $|{symbol}(k)|$ ({_formal_mode_label(band_label)})"
            )
            colorbar = _aligned_colorbar(figure, axis, image)
            _set_component_colorbar_ticks(
                colorbar,
                norm,
                maximum_as_power=(
                    field.is_quarter and panel in {"B2", "B4"}
                ),
            )
            if isinstance(norm, LogNorm):
                colorbar.minorticks_off()
            colorbar.set_label(
                rf"Normalized $|{symbol}(k)|$", rotation=270, labelpad=24
            )
            path = output_dir / f"{panel}_{component}_magnitude_{scale}.png"
            _save_color_map(figure, axis, colorbar, path)
            plt.close(figure)
            outputs.append(path)
    return outputs


def plot_cy_magnitude_log_map(
    frame: pd.DataFrame,
    output_dir: Path,
    band_label: str,
) -> Path:
    """Plot only B4 through the same shared contract as B1--B3."""
    outputs = plot_component_magnitude_maps(
        frame,
        output_dir,
        band_label,
        panels={"B4"},
    )
    if len(outputs) != 1:
        raise RuntimeError(f"Expected one B4 output, got {len(outputs)}")
    return outputs[0]


def redraw_cy_magnitude_log(series_dir: Path) -> list[Path]:
    """Redraw only standard B4 for every band in one existing series."""
    series_dir = Path(series_dir)
    frame = load_polarization_field_view(series_dir)
    output_dir = series_dir / "10_overview"
    outputs = []
    for band_label in frame.band_labels:
        band = frame.for_band(str(band_label))
        outputs.append(
            plot_cy_magnitude_log_map(band, output_dir, str(band_label))
        )
    return outputs


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def plot_polarization_angle_map(
    frame: pd.DataFrame | PolarizationFieldView,
    output_dir: Path,
    band_label: str,
    relative_threshold: float,
) -> Path:
    frame = _coerce_field_view(frame)
    qx, qy, angles = _smooth_cyclic_angles(frame, relative_threshold)
    figure, axis = _color_map_figure()
    image = axis.pcolormesh(
        _cell_edges(qx),
        _cell_edges(qy),
        angles,
        cmap="twilight_shifted",
        norm=Normalize(-0.5 * math.pi, 0.5 * math.pi),
        shading="flat",
    )
    axis.set_aspect("equal", adjustable="box")
    _set_k_space_ticks(axis, qx, qy, is_quarter=frame.is_quarter)
    axis.set_xlabel(r"$k_x/G$")
    axis.set_ylabel(r"$k_y/G$")
    axis.set_title(
        rf"Polarization orientation angle $\psi(k)$ "
        rf"({_formal_mode_label(band_label)} mode)"
    )
    colorbar = _aligned_colorbar(
        figure,
        axis,
        image,
        ticks=[-math.pi / 2, -math.pi / 4, 0, math.pi / 4, math.pi / 2],
    )
    colorbar.ax.set_yticklabels(
        [r"$-\pi/2$", r"$-\pi/4$", "0", r"$\pi/4$", r"$\pi/2$"]
    )
    colorbar.set_label(
        r"Polarization angle $\psi(k)$", rotation=270, labelpad=24
    )
    path = output_dir / "C2_ck_phase.png"
    _save_color_map(figure, axis, colorbar, path)
    plt.close(figure)
    return path


def _plot_vectors(
    axis,
    frame: pd.DataFrame,
    angle_method: str,
    color: str,
    *,
    width: float,
    alpha: float,
) -> None:
    cx, cy = complex_coefficients(frame)
    angles, _vx, _vy = field_angles(
        cx, cy, angle_method, continuous=False
    )
    qx_values = np.sort(frame["qx_over_G"].unique().astype(float))
    qy_values = np.sort(frame["qy_over_G"].unique().astype(float))
    x_edges = _cell_edges(qx_values)
    y_edges = _cell_edges(qy_values)
    x_widths = dict(zip(qx_values, np.diff(x_edges)))
    y_widths = dict(zip(qy_values, np.diff(y_edges)))
    local_length = np.asarray(
        [
            0.58 * min(x_widths[float(qx)], y_widths[float(qy)])
            for qx, qy in zip(frame["qx_over_G"], frame["qy_over_G"])
        ],
        dtype=float,
    )
    axis.quiver(
        frame["qx_over_G"],
        frame["qy_over_G"],
        local_length * np.cos(angles),
        local_length * np.sin(angles),
        angles="xy",
        scale_units="xy",
        scale=1.0,
        pivot="middle",
        color=color,
        width=width,
        alpha=alpha,
        headwidth=0,
        headlength=0,
        headaxislength=0,
    )


def _polarization_ellipse_curve(
    cx: complex, cy: complex, sample_count: int = 81
) -> tuple[np.ndarray, np.ndarray]:
    """Return Re([cx, cy] exp(-it)); positive S3 traverses counter-clockwise."""
    magnitude = math.sqrt(abs(cx) ** 2 + abs(cy) ** 2)
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("Polarization ellipse requires a non-zero Jones vector")
    jones = np.asarray([cx, cy], dtype=complex) / magnitude
    parameter = np.linspace(0.0, 2.0 * math.pi, sample_count)
    curve = np.real(jones[:, None] * np.exp(-1j * parameter)[None, :])
    return curve[0], curve[1]


def _minor_axis_arrow_geometry(
    cx: complex, cy: complex
) -> tuple[np.ndarray, np.ndarray]:
    """Return the minor-axis endpoint and forward tangent of the ellipse."""
    magnitude = math.sqrt(abs(cx) ** 2 + abs(cy) ** 2)
    if not math.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("Polarization ellipse requires a non-zero Jones vector")
    jones = np.asarray([cx, cy], dtype=complex) / magnitude
    ellipse_matrix = np.column_stack([jones.real, jones.imag])
    _left, _singular_values, right_transpose = np.linalg.svd(
        ellipse_matrix, full_matrices=True
    )
    minor_parameter = right_transpose[-1]
    minor_endpoint = ellipse_matrix @ minor_parameter
    forward_parameter = np.asarray(
        [-minor_parameter[1], minor_parameter[0]], dtype=float
    )
    tangent = ellipse_matrix @ forward_parameter
    tangent_norm = float(np.linalg.norm(tangent))
    if tangent_norm <= 0.0:
        raise ValueError("Polarization ellipse has no defined tangent")
    return minor_endpoint, tangent / tangent_norm


def _plot_polarization_ellipses(
    axis,
    frame: pd.DataFrame | PolarizationFieldView | PolarizationLoop,
    *,
    color: str,
    alpha: float,
    sample_axis_count: int | None = None,
    color_by_handedness: bool = False,
) -> int:
    """Draw one locally scaled, arrowed polarization ellipse per valid k point."""
    if isinstance(frame, PolarizationFieldView):
        qx_values, qy_values = frame.qx_axis, frame.qy_axis
        cx_grid, cy_grid = frame.complex_grids()
        valid_grid = frame.grid_matrix("analysis_valid").astype(bool)
        if sample_axis_count is not None:
            sample_axis_count = int(sample_axis_count)
            if sample_axis_count < 2 or sample_axis_count > len(qx_values):
                raise ValueError(
                    "sample_axis_count must be between 2 and the display-axis size"
                )
            x_indices = np.rint(
                np.linspace(0, len(qx_values) - 1, sample_axis_count)
            ).astype(int)
            y_indices = np.rint(
                np.linspace(0, len(qy_values) - 1, sample_axis_count)
            ).astype(int)
            if len(np.unique(x_indices)) != sample_axis_count or len(
                np.unique(y_indices)
            ) != sample_axis_count:
                raise ValueError("sample_axis_count does not produce unique samples")
            qx_values = qx_values[x_indices]
            qy_values = qy_values[y_indices]
            cx_grid = cx_grid[np.ix_(y_indices, x_indices)]
            cy_grid = cy_grid[np.ix_(y_indices, x_indices)]
            valid_grid = valid_grid[np.ix_(y_indices, x_indices)]
        mesh_x, mesh_y = np.meshgrid(qx_values, qy_values)
        sample_qx = mesh_x.ravel()
        sample_qy = mesh_y.ravel()
        sample_cx = cx_grid.ravel()
        sample_cy = cy_grid.ravel()
        valid = valid_grid.ravel()
    elif isinstance(frame, PolarizationLoop):
        sample_qx, sample_qy = frame.qx, frame.qy
        sample_cx, sample_cy = frame.cx, frame.cy
        valid = frame.analysis_valid
        qx_values = np.sort(np.unique(frame.qx))
        qy_values = np.sort(np.unique(frame.qy))
    else:
        sample_qx = frame["qx_over_G"].to_numpy(float)
        sample_qy = frame["qy_over_G"].to_numpy(float)
        sample_cx, sample_cy = complex_coefficients(frame)
        valid = frame["analysis_valid"].to_numpy(bool)
        qx_values = np.sort(frame["qx_over_G"].unique().astype(float))
        qy_values = np.sort(frame["qy_over_G"].unique().astype(float))
    sample_qx = np.asarray(sample_qx)[valid]
    sample_qy = np.asarray(sample_qy)[valid]
    cx = np.asarray(sample_cx)[valid]
    cy = np.asarray(sample_cy)[valid]
    x_widths = dict(zip(qx_values, np.diff(_cell_edges(qx_values))))
    y_widths = dict(zip(qy_values, np.diff(_cell_edges(qy_values))))
    count = 0
    for qx, qy, coefficient_x, coefficient_y in zip(
        sample_qx, sample_qy, cx, cy
    ):
        try:
            ellipse_x, ellipse_y = _polarization_ellipse_curve(
                coefficient_x, coefficient_y
            )
            minor_endpoint, arrow_tangent = _minor_axis_arrow_geometry(
                coefficient_x, coefficient_y
            )
        except ValueError:
            continue
        local_half_size = 0.34 * min(
            x_widths[float(qx)],
            y_widths[float(qy)],
        )
        radius = float(np.max(np.hypot(ellipse_x, ellipse_y)))
        scale = local_half_size / radius
        plotted_x = float(qx) + scale * ellipse_x
        plotted_y = float(qy) + scale * ellipse_y
        ellipse_color = color
        if color_by_handedness:
            s3 = -2.0 * float(
                np.imag(coefficient_x * np.conj(coefficient_y))
            )
            ellipse_color = (
                POSITIVE_HANDEDNESS_COLOR
                if s3 > 0.0
                else NEGATIVE_HANDEDNESS_COLOR
                if s3 < 0.0
                else ZERO_HANDEDNESS_COLOR
            )
        axis.plot(
            plotted_x,
            plotted_y,
            color=ellipse_color,
            linewidth=0.55,
            alpha=alpha,
            solid_capstyle="round",
        )
        arrow_tip = np.asarray(
            [float(qx), float(qy)]
        ) + scale * minor_endpoint
        arrow_length = 0.24 * local_half_size
        arrow_start = arrow_tip - arrow_length * arrow_tangent
        arrow = FancyArrowPatch(
            tuple(arrow_start),
            tuple(arrow_tip),
            arrowstyle="-|>",
            mutation_scale=4.5,
            linewidth=0.5,
            color=ellipse_color,
            alpha=alpha,
            shrinkA=0.0,
            shrinkB=0.0,
        )
        axis.add_patch(arrow)
        count += 1
    return count


def _axis_limits_with_one_step(values: pd.Series | np.ndarray) -> tuple[float, float]:
    raw = values.to_numpy(dtype=float) if isinstance(values, pd.Series) else values
    coordinates = np.sort(np.unique(np.asarray(raw, dtype=float)))
    if coordinates.size < 2:
        raise ValueError("At least two coordinates are required for axis limits.")
    lower_step = float(coordinates[1] - coordinates[0])
    upper_step = float(coordinates[-1] - coordinates[-2])
    return (
        float(coordinates[0]) - lower_step,
        float(coordinates[-1]) + upper_step,
    )


def plot_vector_map(
    frame: pd.DataFrame | PolarizationFieldView,
    output_dir: Path,
    band_label: str,
    angle_method: str,
    relative_threshold: float,
    stable: WindingResult | None,
) -> Path:
    frame = _coerce_field_view(frame)
    figure, axis = plt.subplots(figsize=(5.2, 4.6), constrained_layout=True)
    axis.set_xlim(*_axis_limits_with_one_step(frame.qx_axis))
    axis.set_ylim(*_axis_limits_with_one_step(frame.qy_axis))
    sample_axis_count = (
        int(frame.source["qx_over_G"].nunique()) if frame.is_quarter else None
    )
    _plot_polarization_ellipses(
        axis,
        frame,
        color="0.2",
        alpha=0.9,
        sample_axis_count=sample_axis_count,
        color_by_handedness=frame.is_quarter,
    )
    if frame.is_quarter:
        axis.legend(
            handles=[
                Line2D(
                    [0], [0], color=POSITIVE_HANDEDNESS_COLOR, linewidth=1.4,
                    label=r"Positive ($S_3>0$, CCW)",
                ),
                Line2D(
                    [0], [0], color=NEGATIVE_HANDEDNESS_COLOR, linewidth=1.4,
                    label=r"Negative ($S_3<0$, CW)",
                ),
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=2,
            frameon=False,
            fontsize=7.5,
            handlelength=1.8,
        )
    axis.set_aspect("equal", adjustable="box")
    _set_k_space_ticks(
        axis,
        frame.qx_axis,
        frame.qy_axis,
        is_quarter=frame.is_quarter,
    )
    axis.set_xlabel(r"$k_x/G$")
    axis.set_ylabel(r"$k_y/G$")
    axis.set_title(
        f"Polarization ellipses ({_formal_mode_label(band_label)} mode)",
        pad=34 if frame.is_quarter else None,
    )
    path = output_dir / "C1_ck_vector.png"
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return path


def _plot_winding_row(
    left_axis,
    right_axis,
    frame: pd.DataFrame | PolarizationFieldView,
    angle_method: str,
    stable: WindingResult,
    *,
    left_title: str,
) -> None:
    loop = ordered_square_loop(frame, stable.half_width)
    closed = (
        loop.closed()
        if isinstance(loop, PolarizationLoop)
        else pd.concat([loop, loop.iloc[[0]]], ignore_index=True)
    )
    cx, cy = complex_coefficients(closed)
    angles, _vx, _vy = field_angles(cx, cy, angle_method, continuous=True)
    _plot_polarization_ellipses(
        left_axis, frame, color="0.55", alpha=0.65
    )
    closed_qx = closed.qx if isinstance(closed, PolarizationLoop) else closed["qx_over_G"]
    closed_qy = closed.qy if isinstance(closed, PolarizationLoop) else closed["qy_over_G"]
    left_axis.plot(
        closed_qx,
        closed_qy,
        color="0.25",
        linewidth=1.0,
    )
    _plot_polarization_ellipses(
        left_axis, loop, color="tab:blue", alpha=1.0
    )
    left_axis.set_aspect("equal", adjustable="box")
    left_axis.set_xlabel(r"$k_x/G$")
    left_axis.set_ylabel(r"$k_y/G$")
    left_axis.set_title(left_title)
    if isinstance(frame, PolarizationFieldView):
        qx_values, qy_values = frame.qx_axis, frame.qy_axis
    else:
        qx_values, qy_values = frame["qx_over_G"], frame["qy_over_G"]
    _set_k_space_ticks(
        left_axis,
        qx_values,
        qy_values,
        is_quarter=(
            isinstance(frame, PolarizationFieldView) and frame.is_quarter
        ),
    )
    left_axis.set_xlim(*_axis_limits_with_one_step(qx_values))
    left_axis.set_ylim(*_axis_limits_with_one_step(qy_values))

    index = np.arange(len(angles))
    right_axis.plot(index, angles, color="tab:blue", linewidth=1.5)
    right_axis.plot(
        index,
        np.linspace(float(angles[0]), float(angles[-1]), len(angles)),
        color="0.55",
        linewidth=1.0,
        linestyle="--",
    )
    right_axis.set_xlabel("Counter-clockwise loop index")
    right_axis.set_ylabel("Unwrapped angle (rad)")
    right_axis.set_title(
        f"Winding number = {stable.nearest_integer}"
    )


def plot_winding_detail(
    frame: pd.DataFrame | PolarizationFieldView,
    output_dir: Path,
    band_label: str,
    angle_method: str,
    stable: WindingResult,
    relative_threshold: float,
) -> Path:
    frame = _coerce_field_view(frame)
    if frame.is_quarter:
        figure, axes = plt.subplots(
            1, 2, figsize=(10.4, 4.4), constrained_layout=True
        )
        _plot_winding_row(
            axes[0],
            axes[1],
            frame,
            angle_method,
            stable,
            left_title="Winding loop and polarization ellipses",
        )
        path = output_dir / "D_BIC_winding.png"
        figure.savefig(path, dpi=300)
        plt.close(figure)
        return path

    scaled = frame.scaled_cy(10.0)
    scaled_winding = compute_winding(
        ordered_square_loop(scaled, stable.half_width),
        stable.half_width,
        angle_method,
    )

    figure, axes = plt.subplots(
        2, 2, figsize=(10.4, 8.8), constrained_layout=True
    )
    _plot_winding_row(
        axes[0, 0],
        axes[0, 1],
        frame,
        angle_method,
        stable,
        left_title="Winding loop and polarization ellipses",
    )
    _plot_winding_row(
        axes[1, 0],
        axes[1, 1],
        scaled,
        angle_method,
        scaled_winding,
        left_title=(
            r"Winding loop and polarization ellipses ($c_y \times 10$)"
        ),
    )
    path = output_dir / "D_BIC_winding.png"
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return path


def plot_winding_scan(
    scan: pd.DataFrame,
    output_dir: Path,
    band_label: str,
    angle_method: str,
) -> Path:
    usable = scan[scan["status"].isin(["ok", "unstable"])].copy()
    figure, axes = plt.subplots(
        2, 1, figsize=(5.4, 5.6), sharex=True, constrained_layout=True
    )
    if not usable.empty:
        axes[0].plot(usable["half_width"], usable["winding"], "o-", color="tab:blue")
        axes[1].plot(
            usable["half_width"],
            np.abs(usable["residual"]),
            "o-",
            color="tab:red",
        )
    axes[0].set_ylabel("Winding")
    axes[0].set_title(f"Shrinking-loop scan ({band_label}, {angle_method})")
    axes[1].set_ylabel("|Residual|")
    axes[1].set_xlabel(r"Square half-width $k/G$")
    suffix = angle_method.replace("-", "_")
    path = output_dir / (
        f"winding_loop_scan_{formal_mode_label(band_label)}_{suffix}.png"
    )
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return path


def run_analysis(
    series_dir: str | Path,
    *,
    angle_methods: Sequence[str] = ANGLE_METHODS,
    mask_relative_threshold: float = 1e-6,
) -> dict[str, object]:
    series_dir = Path(series_dir)
    frame = load_polarization_field_view(series_dir)
    output_dir = series_dir / "10_overview"
    logs_dir = series_dir / "80_logs"
    config_dir = series_dir / "99_config"
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)
    unknown = sorted(set(angle_methods).difference(ANGLE_METHODS))
    if unknown:
        raise ValueError(f"Unknown angle methods: {unknown}")
    outputs: list[str] = []
    winding_rows: list[dict[str, object]] = []
    for band_label in frame.band_labels:
        band = frame.for_band(str(band_label))
        display_band_label = formal_mode_label(str(band_label))
        outputs.extend(
            str(path) for path in plot_intensity_maps(band, output_dir, str(band_label))
        )
        outputs.extend(
            str(path)
            for path in plot_component_magnitude_maps(
                band, output_dir, str(band_label)
            )
        )
        outputs.append(
            str(
                plot_polarization_angle_map(
                    band,
                    output_dir,
                    str(band_label),
                    mask_relative_threshold,
                )
            )
        )
        outputs.append(
            str(
                plot_vector_map(
                    band,
                    output_dir,
                    str(band_label),
                    "complex-axis",
                    mask_relative_threshold,
                    None,
                )
            )
        )
        for angle_method in angle_methods:
            scan, stable = scan_shrinking_loops(band, angle_method)
            suffix = angle_method.replace("-", "_")
            scan_path = logs_dir / (
                f"winding_loop_scan_{display_band_label}_{suffix}.csv"
            )
            scan.to_csv(scan_path, index=False)
            if angle_method == "complex-axis":
                display_half_width = common_half_widths(band)[0]
                display_winding = compute_winding(
                    ordered_square_loop(band, display_half_width),
                    display_half_width,
                    angle_method,
                )
                outputs.append(
                    str(
                        plot_winding_detail(
                            band,
                            output_dir,
                            str(band_label),
                            angle_method,
                            display_winding,
                            mask_relative_threshold,
                        )
                    )
                )
            if stable is not None:
                winding_rows.append(
                    {
                        "band_label": display_band_label,
                        "internal_band_label": band_label,
                        **asdict(stable),
                    }
                )
            else:
                winding_rows.append(
                    {
                        "band_label": display_band_label,
                        "internal_band_label": band_label,
                        "angle_method": angle_method,
                        "status": "no_stable_loop",
                    }
                )
    summary_frame = pd.DataFrame(winding_rows)
    summary_csv = logs_dir / "winding_summary.csv"
    summary_frame.to_csv(summary_csv, index=False)
    primary_missing = summary_frame[
        (summary_frame["angle_method"] == "stokes-axis")
        & (summary_frame["status"] != "ok")
    ] if "stokes-axis" in angle_methods else summary_frame.iloc[0:0]
    summary = {
        "status": "incomplete" if not primary_missing.empty else "ok",
        "series_dir": str(series_dir),
        "polarization_definition_version": POLARIZATION_DEFINITION_VERSION,
        "angle_methods": list(angle_methods),
        "mask_relative_threshold": mask_relative_threshold,
        "outputs": outputs + [str(summary_csv)],
        "winding": summary_frame.to_dict(orient="records"),
    }
    source_csv = series_dir / "80_logs" / "polarization_grid.csv"
    source_parquet = series_dir / "01_results" / "polarization_grid.parquet"
    source_table = source_csv if source_csv.is_file() else source_parquet
    source_coordinate_count = int(
        frame.source[["qx_over_G", "qy_over_G"]].drop_duplicates().shape[0]
    )
    render_summary = {
        "status": "ok",
        "sampling_domain": frame.sampling_domain,
        "isQuarter": int(frame.is_quarter),
        "symmetry_group": "C2v" if frame.is_quarter else None,
        "symmetry_render_version": (
            SYMMETRY_RENDER_VERSION if frame.is_quarter else None
        ),
        "symmetry_gauge_version": (
            SYMMETRY_GAUGE_VERSION if frame.is_quarter else None
        ),
        "solved_q_axis_over_G": frame._source_qx.tolist(),
        "plot_q_axis_over_G": frame.qx_axis.tolist(),
        "solved_coordinate_count": source_coordinate_count,
        "plot_coordinate_count": int(len(frame.qx_axis) * len(frame.qy_axis)),
        "plot_coordinate_count_is_data_point_count": False,
        "materialized_symmetry_row_count": 0,
        "source_table": str(source_table),
        "source_sha256": _sha256_file(source_table) if source_table.is_file() else None,
    }
    write_json_atomic(
        config_dir / "symmetry_render_summary.json", render_summary
    )
    summary["symmetry_render_summary"] = str(
        config_dir / "symmetry_render_summary.json"
    )
    write_json_atomic(config_dir / "winding_summary.json", summary)
    if not primary_missing.empty:
        raise ValueError(
            "No stable stokes-axis winding loop for: "
            + ", ".join(primary_missing["band_label"].astype(str))
        )
    return summary
