"""Data-driven local refinement support for the unit-cell 2D ``py`` valley.

The geometry extracted here is navigation data only: it selects new k points.
It never supplies or modifies a polarization magnitude.  Final magnitudes are
accepted only from the coarse COMSOL scan or a refined COMSOL solve.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd


VALID_MATCH_STATUS = frozenset({"gamma", "matched", "repaired"})
K_PRECISION = 6
SOURCE_BAND_LABEL = "p2"
DISPLAY_BAND_LABEL = "py"
REQUIRED_COLUMNS = frozenset(
    {
        "qx_over_G",
        "qy_over_G",
        "band_label",
        "match_status",
        "mode_idx",
        "cx_raw_re",
        "cx_raw_im",
        "cy_raw_re",
        "cy_raw_im",
        "cx_norm_re",
        "cx_norm_im",
        "cy_norm_re",
        "cy_norm_im",
        "normalization_kind",
        "field_plane",
        "polarization_definition_version",
    }
)


@dataclass(frozen=True)
class SourceScan:
    series_dir: Path
    frame: pd.DataFrame
    q_axis: np.ndarray
    coarse_step: float
    q_max: float
    source_table: Path
    source_hashes: dict[str, str]
    config: dict[str, object]
    is_quarter: int = 0

    @property
    def plot_q_axis(self) -> np.ndarray:
        if not self.is_quarter:
            return self.q_axis.copy()
        return np.r_[-self.q_axis[:0:-1], self.q_axis]

    def canonical_coordinate(self, qx: float, qy: float) -> tuple[float, float]:
        if not self.is_quarter:
            return float(qx), float(qy)
        return abs(float(qx)), abs(float(qy))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary.replace(path)


def load_source_scan(
    series_dir: str | Path,
    *,
    band_label: str = SOURCE_BAND_LABEL,
) -> SourceScan:
    """Load and validate one complete Cartesian coarse COMSOL scan."""
    series = Path(series_dir).resolve()
    config_path = series / "99_config" / "config.json"
    csv_path = series / "80_logs" / "polarization_grid.csv"
    parquet_path = series / "01_results" / "polarization_grid.parquet"
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing source configuration: {config_path}")
    if csv_path.is_file():
        table_path = csv_path
        frame = pd.read_csv(csv_path)
    elif parquet_path.is_file():
        table_path = parquet_path
        frame = pd.read_parquet(parquet_path)
    else:
        raise FileNotFoundError(
            "Source series has neither polarization_grid.csv nor parquet"
        )
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Source polarization table is missing columns: {missing}")
    frame = frame.loc[frame["band_label"].astype(str) == band_label].copy()
    if frame.empty:
        raise ValueError(f"Source polarization table has no band {band_label!r}")
    invalid = sorted(set(frame["match_status"].astype(str)) - VALID_MATCH_STATUS)
    if invalid:
        raise ValueError(f"Source scan contains invalid tracked points: {invalid}")
    if set(frame["normalization_kind"].astype(str)) != {"planar_l2"}:
        raise ValueError("Source scan must use planar_l2 normalization")
    if set(frame["field_plane"].astype(str)) != {"air"}:
        raise ValueError("Source scan must use the air field plane")
    if frame[["qx_over_G", "qy_over_G"]].duplicated().any():
        raise ValueError("Source scan contains duplicate Cartesian coordinates")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    identity = config.get("cache_identity", {})
    unit_cell_identity = (
        identity.get("unit_cell_2d", {}) if isinstance(identity, dict) else {}
    )
    is_quarter = identity.get(
        "isQuarter", unit_cell_identity.get("isQuarter", 0)
    )
    if type(is_quarter) is not int or is_quarter not in {0, 1}:
        raise ValueError(f"Source isQuarter must be 0 or 1; got {is_quarter!r}")
    qx = np.sort(frame["qx_over_G"].astype(float).unique())
    qy = np.sort(frame["qy_over_G"].astype(float).unique())
    minimum_axis_points = 2 if is_quarter else 3
    if (
        qx.size < minimum_axis_points
        or qx.size != qy.size
        or not np.allclose(qx, qy, atol=1e-12)
    ):
        raise ValueError("Source scan must use one common square Cartesian axis")
    if is_quarter and (
        qx[0] < -0.5 * 10**-K_PRECISION
        or not math.isclose(float(qx[0]), 0.0, abs_tol=0.5 * 10**-K_PRECISION)
    ):
        raise ValueError("Quarter source scan must start at zero")
    differences = np.diff(qx)
    if not np.allclose(differences, differences[0], atol=1e-12, rtol=1e-10):
        raise ValueError("Source scan axis is not uniformly spaced")
    if len(frame) != qx.size * qy.size:
        raise ValueError("Source scan does not contain a complete Cartesian grid")
    numeric = frame[
        [
            "cx_raw_re",
            "cx_raw_im",
            "cy_raw_re",
            "cy_raw_im",
            "cx_norm_re",
            "cx_norm_im",
            "cy_norm_re",
            "cy_norm_im",
        ]
    ].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("Source normalized polarization contains non-finite values")
    frame.loc[:, numeric.columns] = numeric
    frame["abs_cy_norm"] = np.hypot(frame["cy_norm_re"], frame["cy_norm_im"])
    frame["abs_cy_raw"] = np.hypot(frame["cy_raw_re"], frame["cy_raw_im"])
    expected = int(
        config.get("cache_identity", {}).get("expected_eigenvalue_count", -1)
    )
    if expected != 4:
        raise ValueError(f"Source scan must expect four eigenvalues; got {expected}")
    for row in frame.itertuples(index=False):
        point_dir = (
            series
            / "01_results"
            / "k_points"
            / f"kx={float(row.qx_over_G):.{K_PRECISION}f}_ky={float(row.qy_over_G):.{K_PRECISION}f}"
        )
        metadata = point_dir / "point_metadata.json"
        field = point_dir / "eigenmodes" / f"{int(row.mode_idx):02d}_Hz_center.parquet"
        if not metadata.is_file() or not field.is_file():
            raise FileNotFoundError(
                "Source point cache is incomplete at "
                f"({float(row.qx_over_G):g}, {float(row.qy_over_G):g}): "
                f"metadata={metadata.is_file()}, Hz={field.is_file()}"
            )
    return SourceScan(
        series_dir=series,
        frame=frame.sort_values(["qy_over_G", "qx_over_G"]).reset_index(drop=True),
        q_axis=qx,
        coarse_step=float(differences[0]),
        q_max=float(max(abs(qx[0]), abs(qx[-1]))),
        source_table=table_path,
        source_hashes={
            str(config_path): _sha256(config_path),
            str(table_path): _sha256(table_path),
        },
        config=config,
        is_quarter=is_quarter,
    )


def _parabolic_coordinate(
    coordinates: np.ndarray, values: np.ndarray, index: int
) -> float:
    """Refine a minimum coordinate only; the fitted value is intentionally unused."""
    left, center, right = values[index - 1 : index + 2]
    denominator = left - 2.0 * center + right
    if not np.isfinite(denominator) or denominator <= 0.0:
        return float(coordinates[index])
    offset = 0.5 * (left - right) / denominator
    offset = float(np.clip(offset, -1.0, 1.0))
    return float(coordinates[index] + offset * (coordinates[index + 1] - coordinates[index]))


def extract_valley_candidates(source: SourceScan) -> pd.DataFrame:
    """Find every strict row/column local minimum without an amplitude cutoff."""
    frame = source.frame
    if "abs_cy_raw" not in frame.columns:
        required = {"cy_raw_re", "cy_raw_im"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(
                "Raw cy magnitude cannot be derived; missing columns: "
                + ", ".join(sorted(missing))
            )
        frame = frame.assign(
            abs_cy_raw=np.hypot(frame["cy_raw_re"], frame["cy_raw_im"])
        )
    grid = frame.pivot(
        index="qy_over_G", columns="qx_over_G", values="abs_cy_raw"
    ).loc[source.q_axis, source.q_axis]
    values = grid.to_numpy(dtype=float)
    scan_axis = source.q_axis
    if source.is_quarter:
        mirror_indices = np.r_[
            np.arange(len(source.q_axis) - 1, 0, -1),
            np.arange(len(source.q_axis)),
        ]
        values = values[np.ix_(mirror_indices, mirror_indices)]
        scan_axis = source.plot_q_axis
    rows: list[dict[str, object]] = []
    candidate_id = 0
    for axis_kind in ("row", "column"):
        for slice_index, slice_coordinate in enumerate(scan_axis):
            profile = values[slice_index, :] if axis_kind == "row" else values[:, slice_index]
            for sample_index in range(1, scan_axis.size - 1):
                if not (
                    profile[sample_index] < profile[sample_index - 1]
                    and profile[sample_index] < profile[sample_index + 1]
                ):
                    continue
                refined = _parabolic_coordinate(scan_axis, profile, sample_index)
                sample_qx = (
                    scan_axis[sample_index]
                    if axis_kind == "row"
                    else slice_coordinate
                )
                sample_qy = (
                    slice_coordinate
                    if axis_kind == "row"
                    else scan_axis[sample_index]
                )
                rows.append(
                    {
                        "candidate_id": candidate_id,
                        "axis": axis_kind,
                        "slice_index": slice_index,
                        "slice_coordinate": float(slice_coordinate),
                        "sample_index": sample_index,
                        "sample_qx_over_G": float(sample_qx),
                        "sample_qy_over_G": float(sample_qy),
                        "qx_over_G": float(refined if axis_kind == "row" else slice_coordinate),
                        "qy_over_G": float(slice_coordinate if axis_kind == "row" else refined),
                        "abs_cy_raw_sample": float(profile[sample_index]),
                    }
                )
                candidate_id += 1
    result = pd.DataFrame(rows)
    if source.is_quarter and not result.empty:
        tolerance = 0.5 * 10**-K_PRECISION
        result = result.loc[
            (result["sample_qx_over_G"] >= -tolerance)
            & (result["sample_qy_over_G"] >= -tolerance)
        ].copy()
        for column in (
            "sample_qx_over_G",
            "sample_qy_over_G",
            "qx_over_G",
            "qy_over_G",
        ):
            result[column] = result[column].abs()
        result = result.reset_index(drop=True)
        result["candidate_id"] = np.arange(len(result), dtype=int)
    return result


def build_valley_segments(
    candidates: pd.DataFrame, coarse_step: float
) -> pd.DataFrame:
    """Connect candidates in adjacent slices using deterministic mutual nearests."""
    segments: list[dict[str, object]] = []
    segment_id = 0
    for axis_kind in ("row", "column"):
        subset = candidates.loc[candidates["axis"] == axis_kind]
        by_slice = {
            int(key): group.sort_values("candidate_id")
            for key, group in subset.groupby("slice_index")
        }
        for slice_index in sorted(by_slice):
            if slice_index + 1 not in by_slice:
                continue
            first = by_slice[slice_index]
            second = by_slice[slice_index + 1]
            varying = "qx_over_G" if axis_kind == "row" else "qy_over_G"
            nearest_forward = {
                int(row.candidate_id): int(
                    second.iloc[
                        np.argmin(np.abs(second[varying].to_numpy() - float(getattr(row, varying))))
                    ].candidate_id
                )
                for row in first.itertuples(index=False)
            }
            nearest_reverse = {
                int(row.candidate_id): int(
                    first.iloc[
                        np.argmin(np.abs(first[varying].to_numpy() - float(getattr(row, varying))))
                    ].candidate_id
                )
                for row in second.itertuples(index=False)
            }
            indexed = candidates.set_index("candidate_id")
            for start_id, end_id in nearest_forward.items():
                if nearest_reverse.get(end_id) != start_id:
                    continue
                start = indexed.loc[start_id]
                end = indexed.loc[end_id]
                dx = float(end.qx_over_G - start.qx_over_G)
                dy = float(end.qy_over_G - start.qy_over_G)
                length = math.hypot(dx, dy)
                if length > math.hypot(coarse_step, 2.0 * coarse_step) + 1e-12:
                    continue
                tx, ty = dx / length, dy / length
                segments.append(
                    {
                        "segment_id": segment_id,
                        "axis": axis_kind,
                        "start_candidate_id": start_id,
                        "end_candidate_id": end_id,
                        "qx0_over_G": float(start.qx_over_G),
                        "qy0_over_G": float(start.qy_over_G),
                        "qx1_over_G": float(end.qx_over_G),
                        "qy1_over_G": float(end.qy_over_G),
                        "tangent_x": tx,
                        "tangent_y": ty,
                        "normal_x": -ty,
                        "normal_y": tx,
                        "length_over_G": length,
                    }
                )
                segment_id += 1
    return pd.DataFrame(segments)


def _candidate_tangent_map(
    candidates: pd.DataFrame, segments: pd.DataFrame
) -> dict[int, np.ndarray]:
    vectors: dict[int, list[np.ndarray]] = {
        int(item): [] for item in candidates["candidate_id"]
    }
    for segment in segments.itertuples(index=False):
        vector = np.asarray([segment.tangent_x, segment.tangent_y], dtype=float)
        vectors[int(segment.start_candidate_id)].append(vector)
        vectors[int(segment.end_candidate_id)].append(vector)
    result = {}
    for candidate_id, items in vectors.items():
        if not items:
            continue
        reference = items[0]
        aligned = [item if np.dot(item, reference) >= 0 else -item for item in items]
        mean = np.mean(aligned, axis=0)
        norm = np.linalg.norm(mean)
        if norm > 1e-12:
            result[candidate_id] = mean / norm
    return result


def detect_junction_centers(
    candidates: pd.DataFrame,
    segments: pd.DataFrame,
    coarse_step: float,
) -> pd.DataFrame:
    """Detect sampling ambiguities from row/column evidence, without Γ assumptions."""
    tangents = _candidate_tangent_map(candidates, segments)
    indexed = candidates.set_index("candidate_id")
    row_ids = [int(value) for value in candidates.loc[candidates.axis == "row", "candidate_id"]]
    column_ids = [int(value) for value in candidates.loc[candidates.axis == "column", "candidate_id"]]
    evidence: list[tuple[float, float, int, int]] = []
    for row_id in row_ids:
        if row_id not in tangents:
            continue
        row = indexed.loc[row_id]
        for column_id in column_ids:
            if column_id not in tangents:
                continue
            column = indexed.loc[column_id]
            distance = math.hypot(
                float(row.qx_over_G - column.qx_over_G),
                float(row.qy_over_G - column.qy_over_G),
            )
            if distance > 0.45 * coarse_step:
                continue
            # Parallel row/column detections describe one stable branch.  A
            # directional disagreement marks a region where a 2-D stencil is safer.
            alignment = abs(float(np.dot(tangents[row_id], tangents[column_id])))
            if alignment >= math.cos(math.radians(30.0)):
                continue
            evidence.append(
                (
                    0.5 * float(row.qx_over_G + column.qx_over_G),
                    0.5 * float(row.qy_over_G + column.qy_over_G),
                    row_id,
                    column_id,
                )
            )
    # A change in the number of strict minima between adjacent slices is direct
    # evidence that a one-dimensional normal stencil is not well-defined.  The
    # rule is topology-neutral: it also catches branch birth/death, near
    # splitting, and coarse-grid ambiguity away from Gamma.
    for axis_kind in ("row", "column"):
        subset = candidates.loc[candidates.axis == axis_kind]
        by_slice = {
            int(index): group
            for index, group in subset.groupby("slice_index")
        }
        for slice_index, group in sorted(by_slice.items()):
            adjacent_counts = [
                len(by_slice[other])
                for other in (slice_index - 1, slice_index + 1)
                if other in by_slice
            ]
            if not adjacent_counts or all(count == len(group) for count in adjacent_counts):
                continue
            for row in group.itertuples(index=False):
                candidate_id = int(row.candidate_id)
                evidence.append(
                    (
                        float(row.qx_over_G),
                        float(row.qy_over_G),
                        candidate_id,
                        candidate_id,
                    )
                )
    if not evidence:
        return pd.DataFrame(
            columns=["junction_id", "qx_over_G", "qy_over_G", "evidence_count", "candidate_ids"]
        )
    # Connected components in evidence-coordinate space avoid hard-coding the
    # number or location of ambiguous regions.
    coordinates = np.asarray([[item[0], item[1]] for item in evidence])
    unassigned = set(range(len(evidence)))
    clusters: list[list[int]] = []
    while unassigned:
        seed = min(unassigned)
        component = {seed}
        frontier = [seed]
        unassigned.remove(seed)
        while frontier:
            current = frontier.pop()
            neighbours = [
                other
                for other in sorted(unassigned)
                if np.linalg.norm(coordinates[current] - coordinates[other]) <= 2.1 * coarse_step
            ]
            for other in neighbours:
                unassigned.remove(other)
                component.add(other)
                frontier.append(other)
        clusters.append(sorted(component))
    records = []
    for junction_id, cluster in enumerate(clusters):
        ids = sorted({identifier for index in cluster for identifier in evidence[index][2:]})
        records.append(
            {
                "junction_id": junction_id,
                "qx_over_G": float(np.mean(coordinates[cluster, 0])),
                "qy_over_G": float(np.mean(coordinates[cluster, 1])),
                "evidence_count": len(cluster),
                "candidate_ids": json.dumps(ids),
            }
        )
    return pd.DataFrame(records)


def _coordinate_key(qx: float, qy: float) -> tuple[str, str]:
    return f"{qx:.{K_PRECISION}f}", f"{qy:.{K_PRECISION}f}"


def _select_sampling_anchor_ids(
    candidates: pd.DataFrame,
    segments: pd.DataFrame,
    tangent_map: dict[int, np.ndarray],
    junction_xy: np.ndarray,
    junction_radius: float,
    anchor_stride: int,
    coarse_step: float,
) -> list[int]:
    """Subsample connected valley branches without fitting their geometry."""
    indexed = candidates.set_index("candidate_id")
    eligible = set()
    for row in candidates.itertuples(index=False):
        candidate_id = int(row.candidate_id)
        if candidate_id not in tangent_map:
            continue
        coordinate = np.asarray([row.qx_over_G, row.qy_over_G], dtype=float)
        if junction_xy.size and np.min(np.linalg.norm(junction_xy - coordinate, axis=1)) <= junction_radius:
            continue
        eligible.add(candidate_id)
    adjacency = {candidate_id: set() for candidate_id in eligible}
    for segment in segments.itertuples(index=False):
        start = int(segment.start_candidate_id)
        end = int(segment.end_candidate_id)
        if start in eligible and end in eligible:
            adjacency[start].add(end)
            adjacency[end].add(start)
    components = []
    remaining = set(eligible)
    while remaining:
        seed = min(remaining)
        component = {seed}
        frontier = [seed]
        remaining.remove(seed)
        while frontier:
            current = frontier.pop()
            for neighbour in sorted(adjacency[current] & remaining):
                remaining.remove(neighbour)
                component.add(neighbour)
                frontier.append(neighbour)
        components.append(component)
    selected = []
    for component in components:
        ordered = sorted(
            component,
            key=lambda candidate_id: (
                int(indexed.loc[candidate_id].slice_index), candidate_id
            ),
        )
        chosen = ordered[:: int(anchor_stride)]
        if ordered[-1] not in chosen:
            chosen.append(ordered[-1])
        selected.extend(chosen)
    # Row and column evidence may select the same physical location. Keep one
    # deterministic anchor when they lie within half a coarse interval.
    accepted: list[int] = []
    for candidate_id in sorted(selected):
        row = indexed.loc[candidate_id]
        coordinate = np.asarray([row.qx_over_G, row.qy_over_G], dtype=float)
        if any(
            np.linalg.norm(
                coordinate
                - indexed.loc[other, ["qx_over_G", "qy_over_G"]].to_numpy(dtype=float)
            )
            <= 0.5 * coarse_step
            for other in accepted
        ):
            continue
        accepted.append(candidate_id)
    return accepted


def build_initial_sampling_plan(
    source: SourceScan,
    candidates: pd.DataFrame,
    segments: pd.DataFrame,
    junctions: pd.DataFrame,
    settings,
) -> pd.DataFrame:
    """Build immutable round-1 coordinates; later rounds require real COMSOL data."""
    tangent_map = _candidate_tangent_map(candidates, segments)
    existing = {
        _coordinate_key(float(row.qx_over_G), float(row.qy_over_G))
        for row in source.frame.itertuples(index=False)
    }
    generated: list[dict[str, object]] = []
    junction_radius = max(
        float(settings.junction_half_width_over_g), 1.25 * source.coarse_step
    )
    junction_xy = (
        junctions[["qx_over_G", "qy_over_G"]].to_numpy(dtype=float)
        if not junctions.empty
        else np.empty((0, 2))
    )
    normal_offsets = np.linspace(
        -float(settings.normal_half_width_over_g),
        float(settings.normal_half_width_over_g),
        int(settings.normal_points_per_round),
    )
    tangent_offsets = np.linspace(
        -float(settings.tangent_half_width_over_g),
        float(settings.tangent_half_width_over_g),
        int(settings.tangent_points_per_round),
    )
    anchor_ids = _select_sampling_anchor_ids(
        candidates,
        segments,
        tangent_map,
        junction_xy,
        junction_radius,
        int(settings.anchor_stride),
        source.coarse_step,
    )
    indexed = candidates.set_index("candidate_id")
    for candidate_id in anchor_ids:
        row = indexed.loc[candidate_id]
        tangent = tangent_map[candidate_id]
        center = np.asarray([row.qx_over_G, row.qy_over_G], dtype=float)
        normal = np.asarray([-tangent[1], tangent[0]])
        for stencil_axis, direction, offsets in (
            ("normal", normal, normal_offsets),
            ("tangent", tangent, tangent_offsets),
        ):
            for axis_index, offset in enumerate(offsets):
                qx, qy = center + offset * direction
                if abs(qx) > source.q_max + 1e-12 or abs(qy) > source.q_max + 1e-12:
                    continue
                generated.append(
                    {
                        "plan_source": "valley_cross_stencil",
                        "source_id": candidate_id,
                        "round": 1,
                        "recenter_index": 0,
                        "stencil_axis": stencil_axis,
                        "axis_index": axis_index,
                        "stencil_index": axis_index,
                        "qx_over_G": float(qx),
                        "qy_over_G": float(qy),
                        "center_qx_over_G": float(center[0]),
                        "center_qy_over_G": float(center[1]),
                        "offset_over_G": float(offset),
                        "normal_x": float(normal[0]),
                        "normal_y": float(normal[1]),
                        "tangent_x": float(tangent[0]),
                        "tangent_y": float(tangent[1]),
                        "status": "planned",
                    }
                )
    junction_offsets = np.linspace(
        -float(settings.junction_half_width_over_g),
        float(settings.junction_half_width_over_g),
        int(settings.junction_points_per_axis),
    )
    for junction in junctions.itertuples(index=False):
        for y_index, dy in enumerate(junction_offsets):
            for x_index, dx in enumerate(junction_offsets):
                qx = float(junction.qx_over_G + dx)
                qy = float(junction.qy_over_G + dy)
                if abs(qx) > source.q_max + 1e-12 or abs(qy) > source.q_max + 1e-12:
                    continue
                generated.append(
                    {
                        "plan_source": "junction_grid",
                        "source_id": int(junction.junction_id),
                        "round": 1,
                        "recenter_index": 0,
                        "stencil_axis": "junction",
                        "axis_index": y_index * len(junction_offsets) + x_index,
                        "stencil_index": y_index * len(junction_offsets) + x_index,
                        "qx_over_G": qx,
                        "qy_over_G": qy,
                        "center_qx_over_G": float(junction.qx_over_G),
                        "center_qy_over_G": float(junction.qy_over_G),
                        "offset_over_G": float(math.hypot(dx, dy)),
                        "normal_x": np.nan,
                        "normal_y": np.nan,
                        "tangent_x": np.nan,
                        "tangent_y": np.nan,
                        "status": "planned",
                    }
                )
    # Fill only longitudinal holes left by sparse cross-stencil anchors.  The
    # piecewise segment midpoint is a scan coordinate, never an inferred
    # amplitude.  Existing coarse coordinates are excluded before insertion.
    planned_xy = np.asarray(
        [[item["qx_over_G"], item["qy_over_G"]] for item in generated],
        dtype=float,
    )
    for segment in segments.itertuples(index=False):
        midpoint = np.asarray(
            [
                0.5 * (segment.qx0_over_G + segment.qx1_over_G),
                0.5 * (segment.qy0_over_G + segment.qy1_over_G),
            ],
            dtype=float,
        )
        key = _coordinate_key(float(midpoint[0]), float(midpoint[1]))
        if key in existing:
            continue
        nearest = (
            float(np.min(np.linalg.norm(planned_xy - midpoint, axis=1)))
            if planned_xy.size
            else math.inf
        )
        if nearest <= float(settings.longitudinal_max_gap_over_g) + 1e-12:
            continue
        tangent = np.asarray([segment.tangent_x, segment.tangent_y], dtype=float)
        normal = np.asarray([segment.normal_x, segment.normal_y], dtype=float)
        generated.append(
            {
                "plan_source": "segment_gap_fill",
                "source_id": int(segment.segment_id),
                "round": 1,
                "recenter_index": 0,
                "stencil_axis": "tangent",
                "axis_index": 0,
                "stencil_index": 0,
                "qx_over_G": float(midpoint[0]),
                "qy_over_G": float(midpoint[1]),
                "center_qx_over_G": float(midpoint[0]),
                "center_qy_over_G": float(midpoint[1]),
                "offset_over_G": 0.0,
                "normal_x": float(normal[0]),
                "normal_y": float(normal[1]),
                "tangent_x": float(tangent[0]),
                "tangent_y": float(tangent[1]),
                "status": "planned_gap_fill",
            }
        )
        planned_xy = np.vstack([planned_xy, midpoint])
    # Six-decimal point identity is shared with the unit-cell solver.  Keep the
    # first deterministic provenance record and never average coordinates.
    unique: dict[tuple[str, str], dict[str, object]] = {}
    for record in generated:
        if source.is_quarter:
            qx, qy = source.canonical_coordinate(
                record["qx_over_G"], record["qy_over_G"]
            )
            center_qx, center_qy = source.canonical_coordinate(
                record["center_qx_over_G"], record["center_qy_over_G"]
            )
            record["qx_over_G"], record["qy_over_G"] = qx, qy
            record["center_qx_over_G"] = center_qx
            record["center_qy_over_G"] = center_qy
        key = _coordinate_key(record["qx_over_G"], record["qy_over_G"])
        if key in existing or key in unique:
            continue
        record["kx_str"], record["ky_str"] = key
        unique[key] = record
    records = sorted(
        unique.values(),
        key=lambda item: (
            math.hypot(float(item["qx_over_G"]), float(item["qy_over_G"])),
            float(item["qy_over_G"]),
            float(item["qx_over_G"]),
        ),
    )
    for point_id, record in enumerate(records):
        record["point_id"] = point_id
        record["solve_order"] = point_id
    if len(records) > int(settings.maximum_new_points):
        raise ValueError(
            f"Sampling plan has {len(records)} new points, exceeding maximum_new_points="
            f"{settings.maximum_new_points}"
        )
    return pd.DataFrame(records)


def refined_series_name(parameters, source: SourceScan) -> str:
    settings = parameters.unit_cell_2d
    cell = settings.cell
    source_points = len(source.q_axis)
    return (
        "unit_cell_2D_valley_refined_"
        f"({cell.b0_nm:g}-{cell.eta:g}-{cell.zeta:g})_"
        f"mesh{parameters.mesh_size}_kmax{source.q_max:g}_"
        f"source-uniform{source_points}_"
        f"{'quarter' if source.is_quarter else 'full'}_band-{DISPLAY_BAND_LABEL}"
    )


def _sampling_plan_figure(
    source: SourceScan,
    candidates: pd.DataFrame,
    segments: pd.DataFrame,
    junctions: pd.DataFrame,
    points: pd.DataFrame,
    path: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7.0, 6.2))
    axis.scatter(
        source.frame.qx_over_G,
        source.frame.qy_over_G,
        s=5,
        c="#c8c8c8",
        linewidths=0,
        label="coarse COMSOL",
        zorder=1,
    )
    for segment in segments.itertuples(index=False):
        axis.plot(
            [segment.qx0_over_G, segment.qx1_over_G],
            [segment.qy0_over_G, segment.qy1_over_G],
            color="#2b6cb0" if segment.axis == "row" else "#2f855a",
            linewidth=0.8,
            alpha=0.8,
            zorder=2,
        )
    if not points.empty:
        normal = points.loc[points.plan_source == "valley_cross_stencil"]
        junction = points.loc[points.plan_source == "junction_grid"]
        gap_fill = points.loc[points.plan_source == "segment_gap_fill"]
        axis.scatter(normal.qx_over_G, normal.qy_over_G, s=8, c="#dd6b20", label="normal + tangent stencil", zorder=3)
        axis.scatter(gap_fill.qx_over_G, gap_fill.qy_over_G, s=16, c="#c05621", marker="D", label="longitudinal gap fill", zorder=4)
        axis.scatter(junction.qx_over_G, junction.qy_over_G, s=8, c="#805ad5", label="junction grid", zorder=3)
    if not junctions.empty:
        axis.scatter(junctions.qx_over_G, junctions.qy_over_G, s=42, facecolors="none", edgecolors="#553c9a", linewidths=1.0, zorder=4)
    axis.set(xlabel=r"$k_x/G$", ylabel=r"$k_y/G$", xlim=(-source.q_max, source.q_max), ylim=(-source.q_max, source.q_max))
    axis.set_aspect("equal", adjustable="box")
    axis.legend(frameon=False, fontsize=8, loc="upper right")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_convergence_repair_sampling_plan(
    source: SourceScan,
    existing_refined: pd.DataFrame,
    repair_points: pd.DataFrame,
    path: str | Path,
) -> None:
    """Plot a review-only correction pass without modifying solved data."""
    figure, axis = plt.subplots(figsize=(7.0, 6.2))
    axis.scatter(
        source.frame.qx_over_G,
        source.frame.qy_over_G,
        s=7,
        c="#c7c7c7",
        label="coarse COMSOL",
        zorder=1,
    )
    axis.scatter(
        existing_refined.qx_over_G,
        existing_refined.qy_over_G,
        s=5,
        c="#dd6b20",
        alpha=0.32,
        label="existing refined COMSOL",
        zorder=2,
    )
    if not repair_points.empty:
        axis.scatter(
            repair_points.qx_over_G,
            repair_points.qy_over_G,
            s=14,
            c="#2563eb",
            label="proposed correction points",
            zorder=3,
        )
    axis.set(
        xlabel=r"$k_x/G$",
        ylabel=r"$k_y/G$",
        xlim=(-source.q_max, source.q_max),
        ylim=(-source.q_max, source.q_max),
    )
    axis.set_aspect("equal", adjustable="box")
    axis.legend(frameon=False, fontsize=8, loc="upper right")
    figure.tight_layout()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def prepare_valley_refinement(
    source_series: str | Path,
    parameters,
    output_root: str | Path,
) -> dict[str, object]:
    source = load_source_scan(source_series)
    settings = parameters.unit_cell_2d
    source_cell = source.config.get("unit_cell_2d_geometry", {})
    if source_cell != settings.cell.to_dict():
        raise ValueError(
            "Parameter file unit_cell_2d geometry does not match the source series"
        )
    if parameters.mesh_size != source.config.get("cache_identity", {}).get("mesh_size"):
        raise ValueError("Parameter file mesh_size does not match the source series")
    if settings.eigenmode_pair_count != 2:
        raise ValueError("Valley refinement requires eigenmode_pair_count=2 (four eigenvalues)")
    if not math.isclose(settings.q_max_over_g, source.q_max, abs_tol=1e-12):
        raise ValueError("Parameter q_max_over_G does not match the source series")
    if settings.q_points_per_axis != len(source.q_axis):
        raise ValueError("Parameter q_points_per_axis does not match the source series")
    if settings.is_quarter != source.is_quarter:
        raise ValueError("Parameter isQuarter does not match the source series")
    candidates = extract_valley_candidates(source)
    segments = build_valley_segments(candidates, source.coarse_step)
    junctions = detect_junction_centers(candidates, segments, source.coarse_step)
    points = build_initial_sampling_plan(
        source, candidates, segments, junctions, settings.valley_refinement
    )
    series_dir = Path(output_root).resolve() / refined_series_name(parameters, source)
    model_dir = series_dir / "00_model"
    logs_dir = series_dir / "80_logs"
    config_dir = series_dir / "99_config"
    for directory in (model_dir, logs_dir, config_dir):
        directory.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(logs_dir / "valley_seed_candidates.csv", index=False)
    segments.to_csv(logs_dir / "valley_position_segments.csv", index=False)
    junctions.to_csv(logs_dir / "valley_junctions.csv", index=False)
    points.to_csv(logs_dir / "refinement_points.csv", index=False)
    sampling_figure = model_dir / "valley_refinement_sampling_plan.png"
    _sampling_plan_figure(source, candidates, segments, junctions, points, sampling_figure)
    maximum_round_points = min(
        int(settings.valley_refinement.maximum_new_points),
        len(points)
        * int(settings.valley_refinement.refinement_rounds)
        * (int(settings.valley_refinement.max_edge_recenters) + 1),
    )
    plan_payload = {
        "workflow": "unit_cell_2d_valley_refinement",
        "status": "prepared",
        "source_series": str(source.series_dir),
        "source_table": str(source.source_table),
        "source_hashes": source.source_hashes,
        "band_internal": SOURCE_BAND_LABEL,
        "band_display": DISPLAY_BAND_LABEL,
        "coarse_step_over_G": source.coarse_step,
        "q_max_over_G": source.q_max,
        "coarse_points_per_axis": len(source.q_axis),
        "isQuarter": source.is_quarter,
        "sampling_domain": "quarter" if source.is_quarter else "full",
        "plot_points_per_axis": len(source.plot_q_axis),
        "candidate_count": len(candidates),
        "segment_count": len(segments),
        "junction_count": len(junctions),
        "round_1_unique_new_point_count": len(points),
        "estimated_maximum_points_for_configured_rounds": maximum_round_points,
        "sampling_parameters": settings.valley_refinement.to_dict(),
        "points": points.to_dict(orient="records"),
        "amplitude_contract": "coordinates-only; no fitted amplitude is used",
    }
    _atomic_json(config_dir / "sampling_plan.json", plan_payload)
    _atomic_json(
        config_dir / "config.json",
        {
            "workflow": "unit_cell_2d_valley_refinement",
            "parameter_path": str(parameters.source_path),
            "source_config": str(source.series_dir / "99_config" / "config.json"),
            "source_hashes": source.source_hashes,
            "unit_cell_2d": {
                **settings.cell.to_dict(),
                "eigenmode_pair_count": settings.eigenmode_pair_count,
                "expected_eigenvalue_count": 2 * settings.eigenmode_pair_count,
                "isQuarter": settings.is_quarter,
                "sampling_domain": settings.sampling_domain,
                "normalization_kind": settings.normalization_kind,
                "field_plane": settings.field_plane,
                "valley_refinement": settings.valley_refinement.to_dict(),
            },
        },
    )
    summary = {
        **{key: plan_payload[key] for key in (
            "status", "source_series", "coarse_step_over_G", "candidate_count",
            "segment_count", "junction_count", "round_1_unique_new_point_count",
            "estimated_maximum_points_for_configured_rounds",
        )},
        "output_dir": str(series_dir),
        "sampling_plan_figure": str(sampling_figure),
        "comsol_started": False,
    }
    _atomic_json(config_dir / "run_summary.json", summary)
    return summary


def merge_polarization_points(
    coarse: pd.DataFrame, refined: pd.DataFrame
) -> pd.DataFrame:
    """Merge real COMSOL nodes and reject coordinate/value conflicts."""
    required = {
        "qx_over_G",
        "qy_over_G",
        "cx_raw_re",
        "cx_raw_im",
        "cy_raw_re",
        "cy_raw_im",
        "cx_norm_re",
        "cx_norm_im",
        "cy_norm_re",
        "cy_norm_im",
    }
    for name, frame in (("coarse", coarse), ("refined", refined)):
        missing = sorted(required.difference(frame.columns))
        if missing:
            raise ValueError(f"{name} points are missing columns: {missing}")
    combined = pd.concat(
        [coarse.assign(point_source="coarse"), refined.assign(point_source="refined")],
        ignore_index=True,
    )
    combined["coordinate_key"] = [
        ":".join(_coordinate_key(qx, qy))
        for qx, qy in zip(combined.qx_over_G, combined.qy_over_G)
    ]
    value_columns = [
        "cx_raw_re",
        "cx_raw_im",
        "cy_raw_re",
        "cy_raw_im",
        "cx_norm_re",
        "cx_norm_im",
        "cy_norm_re",
        "cy_norm_im",
    ]
    records = []
    for key, group in combined.groupby("coordinate_key", sort=False):
        values = group[value_columns].to_numpy(dtype=float)
        if not np.allclose(values, values[0], atol=1e-10, rtol=1e-8):
            raise ValueError(f"Conflicting COMSOL values at coordinate {key}")
        records.append(group.iloc[0].to_dict())
    merged = pd.DataFrame(records).drop(columns="coordinate_key")
    # CSV loading may infer coarse k labels such as ``0.000001`` as floats,
    # while the refined parquet cache preserves the same labels as strings.
    # Rebuild both labels from the authoritative numeric coordinates so the
    # merged table has a stable schema and remains parquet-serializable.
    merged["kx_str"] = [
        f"{float(value):.{K_PRECISION}f}" for value in merged["qx_over_G"]
    ]
    merged["ky_str"] = [
        f"{float(value):.{K_PRECISION}f}" for value in merged["qy_over_G"]
    ]
    # The refined figure is a B-panel and therefore uses the same saved raw
    # coefficient contract as standard B1--B4.  One scalar denominator is
    # computed over the complete coarse+refined set; no per-k renormalization
    # enters the displayed magnitude.
    cx = np.hypot(merged.cx_raw_re, merged.cx_raw_im)
    cy = np.hypot(merged.cy_raw_re, merged.cy_raw_im)
    common_max = float(np.max(np.hypot(cx, cy)))
    if not np.isfinite(common_max) or common_max <= 0:
        raise ValueError("Merged normalized polarization has no positive finite scale")
    merged["normalized_cy_magnitude"] = cy / common_max
    merged["common_c_magnitude_max"] = common_max
    return merged


def linear_interpolate_display(
    points: pd.DataFrame,
    q_axis: np.ndarray,
    *,
    is_quarter: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    if type(is_quarter) is not int or is_quarter not in {0, 1}:
        raise ValueError("is_quarter must be 0 or 1")
    sample_x = points.qx_over_G.to_numpy(dtype=float)
    sample_y = points.qy_over_G.to_numpy(dtype=float)
    sample_values = points.normalized_cy_magnitude.to_numpy(dtype=float)
    if is_quarter:
        render_samples: dict[tuple[str, str], tuple[float, float, float]] = {}
        for sign_x, sign_y in ((1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)):
            for qx_value, qy_value, magnitude in zip(
                sign_x * sample_x, sign_y * sample_y, sample_values
            ):
                key = _coordinate_key(qx_value, qy_value)
                render_samples.setdefault(
                    key, (float(qx_value), float(qy_value), float(magnitude))
                )
        sample_x = np.asarray([item[0] for item in render_samples.values()])
        sample_y = np.asarray([item[1] for item in render_samples.values()])
        sample_values = np.asarray([item[2] for item in render_samples.values()])
    triangulation = mtri.Triangulation(
        sample_x,
        sample_y,
    )
    interpolator = mtri.LinearTriInterpolator(
        triangulation, sample_values
    )
    qx, qy = np.meshgrid(q_axis, q_axis)
    return qx, qy, interpolator(qx, qy)


def plot_comsol_refined_b4(
    merged: pd.DataFrame,
    output_path: str | Path,
    *,
    q_max: float,
    display_points_per_axis: int,
    is_quarter: int = 0,
) -> None:
    """Render only a Delaunay piecewise-linear field from real COMSOL nodes."""
    axis_values = np.linspace(-q_max, q_max, int(display_points_per_axis))
    qx, qy, image_values = linear_interpolate_display(
        merged, axis_values, is_quarter=is_quarter
    )
    figure, axis = plt.subplots(figsize=(6.4, 5.4))
    image = axis.pcolormesh(
        qx,
        qy,
        image_values,
        shading="auto",
        cmap="magma",
        norm=LogNorm(vmin=1e-5, vmax=1e-1),
        rasterized=True,
    )
    axis.set_title(r"Normalized $|c_y(k)|$ (py)")
    axis.set_xlabel(r"$k_x/G$")
    axis.set_ylabel(r"$k_y/G$")
    axis.set_xlim(-q_max, q_max)
    axis.set_ylim(-q_max, q_max)
    axis.set_aspect("equal", adjustable="box")
    ticks = np.arange(math.ceil(-q_max / 0.05), math.floor(q_max / 0.05) + 1) * 0.05
    axis.set_xticks(ticks)
    axis.set_yticks(ticks)
    divider = make_axes_locatable(axis)
    colorbar_axis = divider.append_axes("right", size="4.5%", pad=0.16)
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_ticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    colorbar.set_ticklabels([r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$"])
    colorbar.set_label(r"Normalized $|c_y(k)|$", rotation=270, labelpad=24)
    figure.tight_layout()
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)
