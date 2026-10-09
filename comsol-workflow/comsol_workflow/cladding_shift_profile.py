"""Layer-resolved cladding-shift profiles shared by finite workflows."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import math
from numbers import Integral
from pathlib import Path
from typing import Mapping, Protocol


UNIFORM_PROFILE = "uniform"
TANH_POWER_PROFILE = "tanh_power"
SUPPORTED_PROFILE_KINDS = (UNIFORM_PROFILE, TANH_POWER_PROFILE)
GROUPS_SHIFT_GEOMETRY = "groups"
ELLIPSE_SHIFT_GEOMETRY = "ellipse"
SUPPORTED_SHIFT_GEOMETRY_KINDS = (
    GROUPS_SHIFT_GEOMETRY,
    ELLIPSE_SHIFT_GEOMETRY,
)
X_SHIFT_SIDE_INDICES = (0, 2, 3, 5)
X_SHIFT_CORNER_INDICES = (0, 3)
Y_SHIFT_SIDE_INDICES = (1, 4)
Y_SHIFT_CORNER_INDICES = (1, 2, 4, 5)


def cladding_shift_group(kind: str, idx: int) -> str:
    """Map one of the 12 normal-fan regions to the left/right or top/bottom group."""
    if type(idx) is not int or not 0 <= idx < 6:
        raise ValueError("normal-fan region index must be an integer from 0 to 5")
    if kind == "side":
        return "y" if idx in Y_SHIFT_SIDE_INDICES else "x"
    if kind == "corner":
        return "y" if idx in Y_SHIFT_CORNER_INDICES else "x"
    raise ValueError(f"Unknown normal-fan feature kind: {kind}")


def _require_finite_number(
    value: object,
    path: str,
    *,
    positive: bool = False,
    greater_than_one: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{path} must be a number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{path} must be finite")
    if positive and numeric <= 0.0:
        raise ValueError(f"{path} must be positive")
    if greater_than_one and numeric <= 1.0:
        raise ValueError(f"{path} must be greater than 1")
    return numeric


def _format_label_number(value: float) -> str:
    return f"{float(value):.12f}".rstrip("0").rstrip(".")


@dataclass(frozen=True)
class CladdingShiftGeometry:
    """Select how finite cladding-cell target shifts vary with polar angle."""

    kind: str = GROUPS_SHIFT_GEOMETRY

    def __post_init__(self) -> None:
        if self.kind not in SUPPORTED_SHIFT_GEOMETRY_KINDS:
            raise ValueError(
                "cladding_shift_geometry.kind must be one of "
                f"{', '.join(SUPPORTED_SHIFT_GEOMETRY_KINDS)}; got {self.kind!r}"
            )

    @classmethod
    def from_mapping(cls, value: object | None) -> "CladdingShiftGeometry":
        """Load a geometry selector; a missing field preserves groups behavior."""
        if value is None:
            return cls()
        if not isinstance(value, Mapping):
            raise TypeError("cladding_shift_geometry must be an object")
        kind = value.get("kind")
        if not isinstance(kind, str):
            raise TypeError("cladding_shift_geometry.kind must be a string")
        unknown = sorted(set(value).difference({"kind"}))
        if unknown:
            raise ValueError(
                "cladding_shift_geometry contains unknown fields: "
                + ", ".join(unknown)
            )
        return cls(kind=kind)

    @property
    def series_suffix(self) -> str:
        if self.kind == ELLIPSE_SHIFT_GEOMETRY:
            return "_ellipse-shift"
        return ""

    def to_dict(self) -> dict[str, str]:
        return {"kind": self.kind}


def ellipse_shift_weight(
    x: float,
    y: float,
    x_factor: float,
    y_factor: float,
) -> float:
    """Return f_x*cos(theta)^2 + f_y*sin(theta)^2 at an ideal center."""
    x_position = _require_finite_number(x, "ellipse shift x")
    y_position = _require_finite_number(y, "ellipse shift y")
    x_target = _require_finite_number(x_factor, "ellipse shift x_factor")
    y_target = _require_finite_number(y_factor, "ellipse shift y_factor")
    if x_target < 0.0 or y_target < 0.0:
        raise ValueError("ellipse shift factors must be non-negative")
    radius_squared = x_position * x_position + y_position * y_position
    if radius_squared == 0.0:
        return 0.0
    return (
        x_target * x_position * x_position
        + y_target * y_position * y_position
    ) / radius_squared


def ellipse_inward_shift(
    x: float,
    y: float,
    period: float,
    x_factor: float,
    y_factor: float,
    *,
    envelope: float = 1.0,
) -> tuple[float, float]:
    """Return the smooth radial inward shift of one ideal cladding-cell center."""
    x_position = _require_finite_number(x, "ellipse shift x")
    y_position = _require_finite_number(y, "ellipse shift y")
    lattice_period = _require_finite_number(
        period,
        "ellipse shift period",
        positive=True,
    )
    resolved_envelope = _require_finite_number(
        envelope,
        "ellipse shift envelope",
    )
    if resolved_envelope < 0.0:
        raise ValueError("ellipse shift envelope must be non-negative")
    radius = math.hypot(x_position, y_position)
    if radius == 0.0 or resolved_envelope == 0.0:
        return (0.0, 0.0)
    weight = ellipse_shift_weight(
        x_position,
        y_position,
        x_factor,
        y_factor,
    )
    magnitude = lattice_period * resolved_envelope * weight
    return (
        -magnitude * x_position / radius,
        -magnitude * y_position / radius,
    )


@dataclass(frozen=True)
class CladdingShiftProfile:
    """Validated scalar envelope applied to a far-field cladding-shift target."""

    kind: str = UNIFORM_PROFILE
    scale_layers: float | None = None
    power: float | None = None

    def __post_init__(self) -> None:
        if self.kind not in SUPPORTED_PROFILE_KINDS:
            raise ValueError(
                "cladding_shift_profile.kind must be one of "
                f"{', '.join(SUPPORTED_PROFILE_KINDS)}; got {self.kind!r}"
            )
        if self.kind == UNIFORM_PROFILE:
            if self.scale_layers is not None or self.power is not None:
                raise ValueError(
                    "uniform cladding_shift_profile does not accept "
                    "scale_layers or power"
                )
            return

        scale_layers = _require_finite_number(
            self.scale_layers,
            "cladding_shift_profile.scale_layers",
            positive=True,
        )
        power = _require_finite_number(
            self.power,
            "cladding_shift_profile.power",
            greater_than_one=True,
        )
        object.__setattr__(self, "scale_layers", scale_layers)
        object.__setattr__(self, "power", power)

    @classmethod
    def from_mapping(cls, value: object | None) -> "CladdingShiftProfile":
        """Load a profile; a missing field reproduces the legacy uniform shift."""
        if value is None:
            return cls()
        if not isinstance(value, Mapping):
            raise TypeError("cladding_shift_profile must be an object")
        kind = value.get("kind")
        if not isinstance(kind, str):
            raise TypeError("cladding_shift_profile.kind must be a string")

        allowed = {"kind"}
        if kind == TANH_POWER_PROFILE:
            allowed.update({"scale_layers", "power"})
        unknown = sorted(set(value).difference(allowed))
        if unknown:
            raise ValueError(
                "cladding_shift_profile contains unknown fields: "
                + ", ".join(unknown)
            )
        if kind == TANH_POWER_PROFILE:
            missing = sorted({"scale_layers", "power"}.difference(value))
            if missing:
                raise ValueError(
                    "cladding_shift_profile is missing: " + ", ".join(missing)
                )
            return cls(
                kind=kind,
                scale_layers=value["scale_layers"],
                power=value["power"],
            )
        return cls(kind=kind)

    @property
    def label(self) -> str:
        if self.kind == UNIFORM_PROFILE:
            return UNIFORM_PROFILE
        return (
            "tanhpow-l"
            f"{_format_label_number(float(self.scale_layers))}"
            "-p"
            f"{_format_label_number(float(self.power))}"
        )

    @property
    def series_suffix(self) -> str:
        if self.kind == UNIFORM_PROFILE:
            return ""
        return f"_shiftprof-{self.label}"

    def to_dict(self) -> dict[str, float | str]:
        if self.kind == UNIFORM_PROFILE:
            return {"kind": self.kind}
        return {
            "kind": self.kind,
            "scale_layers": float(self.scale_layers),
            "power": float(self.power),
        }

    def envelope(self, layer: int) -> float:
        if not isinstance(layer, Integral):
            raise TypeError("cladding shift layer must be an integer")
        normalized_layer = int(layer)
        if normalized_layer < 0:
            raise ValueError("cladding shift layer must be non-negative")
        if normalized_layer == 0:
            return 0.0
        if self.kind == UNIFORM_PROFILE:
            return 1.0
        argument = (normalized_layer / float(self.scale_layers)) ** float(self.power)
        return math.tanh(argument)


class _ShellPoint(Protocol):
    shell: int


def relative_cladding_layer(point: _ShellPoint, bulk_radius: int) -> int:
    """Return the C6-symmetric shell offset from the outer cavity shell."""
    if type(bulk_radius) is not int or bulk_radius < 0:
        raise ValueError("bulk_radius must be a non-negative integer")
    if not isinstance(point.shell, Integral):
        raise TypeError("lattice point shell must be an integer")
    layer = int(point.shell) - bulk_radius
    if layer < 0:
        raise ValueError(
            f"lattice point shell {point.shell} lies inside bulk radius {bulk_radius}"
        )
    return layer


def resolved_shift_profile(
    cladding_layers: int,
    target_factor: float,
    period: float,
    profile: CladdingShiftProfile,
) -> list[dict[str, float | int]]:
    """Resolve absolute shifts and discrete differences from layer 0 through N."""
    if type(cladding_layers) is not int or cladding_layers <= 0:
        raise ValueError("cladding_layers must be a positive integer")
    target = _require_finite_number(target_factor, "target_factor")
    lattice_period = _require_finite_number(period, "period", positive=True)

    rows: list[dict[str, float | int]] = []
    effective_factors: list[float] = []
    for layer in range(cladding_layers + 1):
        envelope = profile.envelope(layer)
        effective_factor = target * envelope
        effective_factors.append(effective_factor)
        first_difference = (
            0.0 if layer == 0 else effective_factor - effective_factors[layer - 1]
        )
        second_difference = (
            0.0
            if layer < 2
            else effective_factor
            - 2.0 * effective_factors[layer - 1]
            + effective_factors[layer - 2]
        )
        side_shift = lattice_period * effective_factor
        rows.append(
            {
                "layer": layer,
                "envelope": envelope,
                "target_factor": target,
                "effective_factor": effective_factor,
                "side_shift_um": side_shift,
                "corner_shift_um": 2.0 * side_shift / math.sqrt(3.0),
                "first_difference_factor": first_difference,
                "second_difference_factor": second_difference,
                "side_shift_first_difference_um": lattice_period
                * first_difference,
            }
        )
    return rows


def shift_profile_metadata(
    cladding_layers: int,
    target_factor: float,
    period: float,
    profile: CladdingShiftProfile,
) -> dict[str, object]:
    rows = resolved_shift_profile(
        cladding_layers,
        target_factor,
        period,
        profile,
    )
    bulk_start_layer = max(1, cladding_layers - 2)
    return {
        "cladding_shift_profile": profile.to_dict(),
        "cladding_shift_profile_label": profile.label,
        "resolved_cladding_shift_by_layer": rows,
        "outer_envelope": float(rows[-1]["envelope"]),
        "bulk_start_layer": bulk_start_layer,
        "bulk_start_envelope": float(rows[bulk_start_layer]["envelope"]),
        "max_first_difference": max(
            abs(float(row["first_difference_factor"])) for row in rows
        ),
        "max_second_difference": max(
            abs(float(row["second_difference_factor"])) for row in rows
        ),
        "max_side_shift_first_difference_um": max(
            abs(float(row["side_shift_first_difference_um"])) for row in rows
        ),
    }


def directional_shift_profile_metadata(
    cladding_layers: int,
    x_target_factor: float,
    y_target_factor: float,
    period: float,
    profile: CladdingShiftProfile,
) -> dict[str, object]:
    """Resolve the shared layer envelope for left/right and top/bottom groups."""
    rows_by_group = {
        "x": resolved_shift_profile(
            cladding_layers,
            x_target_factor,
            period,
            profile,
        ),
        "y": resolved_shift_profile(
            cladding_layers,
            y_target_factor,
            period,
            profile,
        ),
    }
    bulk_start_layer = max(1, cladding_layers - 2)

    def group_max(field: str) -> dict[str, float]:
        return {
            group: max(abs(float(row[field])) for row in rows)
            for group, rows in rows_by_group.items()
        }

    max_first_by_group = group_max("first_difference_factor")
    max_second_by_group = group_max("second_difference_factor")
    max_side_first_by_group = group_max("side_shift_first_difference_um")
    reference_rows = rows_by_group["x"]
    return {
        "cladding_x_shift_factor": float(x_target_factor),
        "cladding_y_shift_factor": float(y_target_factor),
        "cladding_x_inward_shift": float(period) * float(x_target_factor),
        "cladding_y_inward_shift": float(period) * float(y_target_factor),
        "cladding_shift_profile": profile.to_dict(),
        "cladding_shift_profile_label": profile.label,
        "cladding_shift_region_groups": {
            "x": {
                "sides": list(X_SHIFT_SIDE_INDICES),
                "corners": list(X_SHIFT_CORNER_INDICES),
            },
            "y": {
                "sides": list(Y_SHIFT_SIDE_INDICES),
                "corners": list(Y_SHIFT_CORNER_INDICES),
            },
        },
        "resolved_cladding_shift_by_group": rows_by_group,
        "outer_envelope": float(reference_rows[-1]["envelope"]),
        "bulk_start_layer": bulk_start_layer,
        "bulk_start_envelope": float(
            reference_rows[bulk_start_layer]["envelope"]
        ),
        "max_first_difference_by_group": max_first_by_group,
        "max_second_difference_by_group": max_second_by_group,
        "max_side_shift_first_difference_um_by_group": max_side_first_by_group,
        "max_first_difference": max(max_first_by_group.values()),
        "max_second_difference": max(max_second_by_group.values()),
        "max_side_shift_first_difference_um": max(
            max_side_first_by_group.values()
        ),
    }


def ellipse_shift_profile_metadata(
    cladding_layers: int,
    x_target_factor: float,
    y_target_factor: float,
    period: float,
    profile: CladdingShiftProfile,
) -> dict[str, object]:
    """Resolve layer-dependent shifts on the two principal ellipse axes."""
    x_target = _require_finite_number(
        x_target_factor,
        "cladding_x_shift_factor",
    )
    y_target = _require_finite_number(
        y_target_factor,
        "cladding_y_shift_factor",
    )
    if x_target < 0.0 or y_target < 0.0:
        raise ValueError("cladding ellipse shift factors must be non-negative")

    def axis_rows(target: float) -> list[dict[str, float | int]]:
        return [
            {
                "layer": row["layer"],
                "envelope": row["envelope"],
                "target_factor": row["target_factor"],
                "effective_factor": row["effective_factor"],
                "axis_shift_um": row["side_shift_um"],
                "first_difference_factor": row["first_difference_factor"],
                "second_difference_factor": row["second_difference_factor"],
                "axis_shift_first_difference_um": row[
                    "side_shift_first_difference_um"
                ],
            }
            for row in resolved_shift_profile(
                cladding_layers,
                target,
                period,
                profile,
            )
        ]

    rows_by_axis = {
        "x": axis_rows(x_target),
        "y": axis_rows(y_target),
    }
    bulk_start_layer = max(1, cladding_layers - 2)

    def axis_max(field: str) -> dict[str, float]:
        return {
            axis: max(abs(float(row[field])) for row in rows)
            for axis, rows in rows_by_axis.items()
        }

    max_first_by_axis = axis_max("first_difference_factor")
    max_second_by_axis = axis_max("second_difference_factor")
    max_shift_first_by_axis = axis_max("axis_shift_first_difference_um")
    reference_rows = rows_by_axis["x"]
    return {
        "cladding_x_shift_factor": x_target,
        "cladding_y_shift_factor": y_target,
        "cladding_x_inward_shift": float(period) * x_target,
        "cladding_y_inward_shift": float(period) * y_target,
        "cladding_shift_profile": profile.to_dict(),
        "cladding_shift_profile_label": profile.label,
        "resolved_cladding_shift_by_axis": rows_by_axis,
        "outer_envelope": float(reference_rows[-1]["envelope"]),
        "bulk_start_layer": bulk_start_layer,
        "bulk_start_envelope": float(
            reference_rows[bulk_start_layer]["envelope"]
        ),
        "max_first_difference_by_axis": max_first_by_axis,
        "max_second_difference_by_axis": max_second_by_axis,
        "max_axis_shift_first_difference_um_by_axis": max_shift_first_by_axis,
        "max_first_difference": max(max_first_by_axis.values()),
        "max_second_difference": max(max_second_by_axis.values()),
        "max_axis_shift_first_difference_um": max(
            max_shift_first_by_axis.values()
        ),
    }


def finite_shift_profile_metadata(
    geometry: CladdingShiftGeometry,
    cladding_layers: int,
    x_target_factor: float,
    y_target_factor: float,
    period: float,
    profile: CladdingShiftProfile,
) -> dict[str, object]:
    """Return mode-aware finite cladding-shift metadata."""
    if not isinstance(geometry, CladdingShiftGeometry):
        raise TypeError("geometry must be a CladdingShiftGeometry")
    if geometry.kind == GROUPS_SHIFT_GEOMETRY:
        metadata = directional_shift_profile_metadata(
            cladding_layers,
            x_target_factor,
            y_target_factor,
            period,
            profile,
        )
        resolved = metadata["resolved_cladding_shift_by_group"]
        formula = {
            "kind": "normal_fan_groups_v1",
            "side_shift": "A * envelope(layer) * group_factor",
            "corner_scale": 2.0 / math.sqrt(3.0),
        }
    else:
        metadata = ellipse_shift_profile_metadata(
            cladding_layers,
            x_target_factor,
            y_target_factor,
            period,
            profile,
        )
        resolved = metadata["resolved_cladding_shift_by_axis"]
        formula = {
            "kind": "radial_cos2_sin2_v1",
            "origin": [0.0, 0.0],
            "angular_weight": "fx*cos(theta)^2 + fy*sin(theta)^2",
            "shift": "-A*envelope(layer)*angular_weight*p/|p|",
        }
    return {
        "cladding_shift_geometry": geometry.to_dict(),
        "cladding_shift_formula": formula,
        "resolved_cladding_shift": {
            "kind": geometry.kind,
            "values": resolved,
        },
        **metadata,
    }


def shift_profile_preflight_messages(
    metadata: Mapping[str, object],
    *,
    minimum_outer_envelope: float = 0.999,
    minimum_bulk_start_envelope: float = 0.99,
) -> list[str]:
    """Describe far-field saturation and warn without changing the profile."""
    profile_label = str(metadata.get("cladding_shift_profile_label", "unknown"))
    outer_envelope = _require_finite_number(
        metadata.get("outer_envelope"),
        "outer_envelope",
    )
    bulk_start_layer = metadata.get("bulk_start_layer")
    if not isinstance(bulk_start_layer, Integral):
        raise TypeError("bulk_start_layer must be an integer")
    bulk_start_envelope = _require_finite_number(
        metadata.get("bulk_start_envelope"),
        "bulk_start_envelope",
    )
    max_first_difference = _require_finite_number(
        metadata.get("max_first_difference"),
        "max_first_difference",
    )
    max_second_difference = _require_finite_number(
        metadata.get("max_second_difference"),
        "max_second_difference",
    )
    messages = [
        "Cladding-shift profile preflight: "
        f"profile={profile_label}, outer envelope={outer_envelope:.6f}, "
        f"layer {int(bulk_start_layer)} envelope={bulk_start_envelope:.6f}, "
        f"max |first difference|={max_first_difference:.6g}, "
        f"max |second difference|={max_second_difference:.6g}."
    ]
    if (
        outer_envelope < minimum_outer_envelope
        or bulk_start_envelope < minimum_bulk_start_envelope
    ):
        messages.append(
            "WARNING: cladding-shift profile has not reached the configured "
            "far-field saturation guard; inputs were left unchanged."
        )
    return messages


def save_shift_profile_diagnostics(
    directory: str | Path,
    metadata: Mapping[str, object],
    *,
    dpi: int = 220,
) -> tuple[Path, Path]:
    """Write the resolved layer table and a compact geometry-only profile plot."""
    rows = metadata.get("resolved_cladding_shift_by_layer")
    if not isinstance(rows, list) or not rows or not all(
        isinstance(row, Mapping) for row in rows
    ):
        raise ValueError("resolved_cladding_shift_by_layer must be a non-empty list")
    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "cladding_shift_profile.csv"
    png_path = output_dir / "cladding_shift_profile.png"

    fieldnames = [
        "layer",
        "envelope",
        "target_factor",
        "effective_factor",
        "side_shift_um",
        "corner_shift_um",
        "first_difference_factor",
        "second_difference_factor",
        "side_shift_first_difference_um",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({key: row[key] for key in fieldnames} for row in rows)

    import matplotlib.pyplot as plt

    layers = [int(row["layer"]) for row in rows]
    side_shifts = [float(row["side_shift_um"]) for row in rows]
    corner_shifts = [float(row["corner_shift_um"]) for row in rows]
    first_differences = [
        float(row["side_shift_first_difference_um"]) for row in rows
    ]

    fig, (ax_shift, ax_difference) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        sharex=True,
        constrained_layout=True,
    )
    ax_shift.plot(layers, side_shifts, marker="o", label="side absolute shift")
    ax_shift.plot(
        layers,
        corner_shifts,
        marker="s",
        label="corner absolute shift",
    )
    ax_shift.set_ylabel("absolute shift (um)")
    ax_shift.grid(alpha=0.25)
    ax_shift.legend()

    ax_difference.bar(layers, first_differences, color="#0f766e", alpha=0.82)
    ax_difference.axhline(0.0, color="#111827", linewidth=0.8)
    ax_difference.set_xlabel("relative cladding layer")
    ax_difference.set_ylabel("side first difference (um)")
    ax_difference.set_xticks(layers)
    ax_difference.grid(axis="y", alpha=0.25)

    profile_label = str(metadata.get("cladding_shift_profile_label", "unknown"))
    target_factor = float(rows[-1]["target_factor"])
    fig.suptitle(
        f"Cladding shift profile: {profile_label}; target={target_factor:.6g} A"
    )
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    return csv_path, png_path


def save_directional_shift_profile_diagnostics(
    directory: str | Path,
    metadata: Mapping[str, object],
    *,
    dpi: int = 220,
) -> tuple[Path, Path]:
    """Write x/y group layer tables and a combined geometry-only profile plot."""
    rows_by_group = metadata.get("resolved_cladding_shift_by_group")
    if not isinstance(rows_by_group, Mapping) or set(rows_by_group) != {"x", "y"}:
        raise ValueError("resolved_cladding_shift_by_group must contain x and y")
    if not all(
        isinstance(rows_by_group[group], list)
        and rows_by_group[group]
        and all(isinstance(row, Mapping) for row in rows_by_group[group])
        for group in ("x", "y")
    ):
        raise ValueError("directional shift rows must be non-empty lists")

    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "cladding_shift_profile.csv"
    png_path = output_dir / "cladding_shift_profile.png"
    value_fields = [
        "layer",
        "envelope",
        "target_factor",
        "effective_factor",
        "side_shift_um",
        "corner_shift_um",
        "first_difference_factor",
        "second_difference_factor",
        "side_shift_first_difference_um",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group", *value_fields])
        writer.writeheader()
        for group in ("x", "y"):
            for row in rows_by_group[group]:
                writer.writerow(
                    {"group": group, **{key: row[key] for key in value_fields}}
                )

    import matplotlib.pyplot as plt

    fig, (ax_shift, ax_difference) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        sharex=True,
        constrained_layout=True,
    )
    colors = {"x": "#2563eb", "y": "#dc2626"}
    for group in ("x", "y"):
        rows = rows_by_group[group]
        layers = [int(row["layer"]) for row in rows]
        side_shifts = [float(row["side_shift_um"]) for row in rows]
        corner_shifts = [float(row["corner_shift_um"]) for row in rows]
        first_differences = [
            float(row["side_shift_first_difference_um"]) for row in rows
        ]
        ax_shift.plot(
            layers,
            side_shifts,
            marker="o",
            color=colors[group],
            label=f"{group} group side",
        )
        ax_shift.plot(
            layers,
            corner_shifts,
            linestyle="--",
            color=colors[group],
            label=f"{group} group corner",
        )
        ax_difference.plot(
            layers,
            first_differences,
            marker="o",
            color=colors[group],
            label=f"{group} group",
        )

    layers = [int(row["layer"]) for row in rows_by_group["x"]]
    ax_shift.set_ylabel("absolute shift (um)")
    ax_shift.grid(alpha=0.25)
    ax_shift.legend()
    ax_difference.axhline(0.0, color="#111827", linewidth=0.8)
    ax_difference.set_xlabel("relative cladding layer")
    ax_difference.set_ylabel("side first difference (um)")
    ax_difference.set_xticks(layers)
    ax_difference.grid(alpha=0.25)
    ax_difference.legend()

    profile_label = str(metadata.get("cladding_shift_profile_label", "unknown"))
    x_target = float(rows_by_group["x"][-1]["target_factor"])
    y_target = float(rows_by_group["y"][-1]["target_factor"])
    fig.suptitle(
        "Directional cladding shift: "
        f"{profile_label}; x={x_target:.6g} A; y={y_target:.6g} A"
    )
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    return csv_path, png_path


def save_ellipse_shift_profile_diagnostics(
    directory: str | Path,
    metadata: Mapping[str, object],
    *,
    dpi: int = 220,
) -> tuple[Path, Path]:
    """Write principal-axis layer tables for the smooth radial ellipse mode."""
    rows_by_axis = metadata.get("resolved_cladding_shift_by_axis")
    if not isinstance(rows_by_axis, Mapping) or set(rows_by_axis) != {"x", "y"}:
        raise ValueError("resolved_cladding_shift_by_axis must contain x and y")
    if not all(
        isinstance(rows_by_axis[axis], list)
        and rows_by_axis[axis]
        and all(isinstance(row, Mapping) for row in rows_by_axis[axis])
        for axis in ("x", "y")
    ):
        raise ValueError("ellipse shift rows must be non-empty lists")

    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "cladding_shift_profile.csv"
    png_path = output_dir / "cladding_shift_profile.png"
    value_fields = [
        "layer",
        "envelope",
        "target_factor",
        "effective_factor",
        "axis_shift_um",
        "first_difference_factor",
        "second_difference_factor",
        "axis_shift_first_difference_um",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["axis", *value_fields])
        writer.writeheader()
        for axis in ("x", "y"):
            for row in rows_by_axis[axis]:
                writer.writerow(
                    {"axis": axis, **{key: row[key] for key in value_fields}}
                )

    import matplotlib.pyplot as plt

    fig, (ax_shift, ax_difference) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        sharex=True,
        constrained_layout=True,
    )
    colors = {"x": "#2563eb", "y": "#dc2626"}
    for axis in ("x", "y"):
        rows = rows_by_axis[axis]
        layers = [int(row["layer"]) for row in rows]
        shifts = [float(row["axis_shift_um"]) for row in rows]
        first_differences = [
            float(row["axis_shift_first_difference_um"]) for row in rows
        ]
        ax_shift.plot(
            layers,
            shifts,
            marker="o",
            color=colors[axis],
            label=f"{axis} principal axis",
        )
        ax_difference.plot(
            layers,
            first_differences,
            marker="o",
            color=colors[axis],
            label=f"{axis} principal axis",
        )

    layers = [int(row["layer"]) for row in rows_by_axis["x"]]
    ax_shift.set_ylabel("principal-axis shift (um)")
    ax_shift.grid(alpha=0.25)
    ax_shift.legend()
    ax_difference.axhline(0.0, color="#111827", linewidth=0.8)
    ax_difference.set_xlabel("relative cladding layer")
    ax_difference.set_ylabel("first difference (um)")
    ax_difference.set_xticks(layers)
    ax_difference.grid(alpha=0.25)
    ax_difference.legend()

    profile_label = str(metadata.get("cladding_shift_profile_label", "unknown"))
    x_target = float(rows_by_axis["x"][-1]["target_factor"])
    y_target = float(rows_by_axis["y"][-1]["target_factor"])
    fig.suptitle(
        "Radial ellipse cladding shift: "
        f"{profile_label}; x={x_target:.6g} A; y={y_target:.6g} A"
    )
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    return csv_path, png_path


def save_finite_shift_profile_diagnostics(
    directory: str | Path,
    metadata: Mapping[str, object],
    *,
    dpi: int = 220,
) -> tuple[Path, Path]:
    """Dispatch finite profile diagnostics according to the geometry kind."""
    geometry = CladdingShiftGeometry.from_mapping(
        metadata.get("cladding_shift_geometry")
    )
    if geometry.kind == ELLIPSE_SHIFT_GEOMETRY:
        return save_ellipse_shift_profile_diagnostics(
            directory,
            metadata,
            dpi=dpi,
        )
    return save_directional_shift_profile_diagnostics(
        directory,
        metadata,
        dpi=dpi,
    )
