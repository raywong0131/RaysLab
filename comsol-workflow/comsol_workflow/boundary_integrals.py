"""Pure numerical operations for S4; coordinates in metres, no COMSOL import."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss

EPS0 = 8.8541878128e-12


def polygon_area(p):
    p = np.asarray(p, float)
    return float(np.sum(p[:, 0] * np.roll(p[:, 1], -1)
                        - p[:, 1] * np.roll(p[:, 0], -1)) / 2)


@dataclass(frozen=True)
class Edge:
    hole: int
    edge: int
    group: str
    start: np.ndarray
    end: np.ndarray
    normal: np.ndarray

    @property
    def length(self):
        return float(np.linalg.norm(self.end - self.start))


def polygon_edges(polygon, *, hole=0, group="outer", inward=False):
    """Retain input vertex IDs. Choose material normal independent of winding."""
    p = np.asarray(polygon, float)
    area = polygon_area(p)
    if not np.isfinite(p).all() or area == 0:
        raise ValueError("Invalid polygon")
    edges = []
    for j, (start, end) in enumerate(zip(p, np.roll(p, -1, axis=0)), 1):
        tangent = end - start
        length = np.linalg.norm(tangent)
        if length <= 0:
            raise ValueError("Zero-length edge")
        normal = np.sign(area) * np.array([tangent[1], -tangent[0]]) / length
        if inward:
            normal = -normal
        midpoint = (start + end) / 2
        if inward and np.dot(normal, p.mean(axis=0) - midpoint) <= 0:
            raise ValueError("Hole normal does not point into air")
        edges.append(Edge(hole, j, group, start, end, normal))
    return edges


def hole_edges(holes):
    if len(holes) != 6 or any(np.shape(h) != (3, 2) for h in holes):
        raise ValueError("S4 requires exactly six triangular holes")
    # The existing geometry generator orders holes counter-clockwise from +x.
    centers = np.asarray([np.mean(h, axis=0) for h in holes])
    expected = np.arange(6) * np.pi / 3
    angles = np.arctan2(centers[:, 1], centers[:, 0])
    if np.max(np.abs(np.angle(np.exp(1j * (angles - expected))))) > 1e-8:
        raise ValueError("Unexpected hole numbering; verify manuscript coordinates")
    return [e for n, h in enumerate(holes, 1)
            for e in polygon_edges(h, hole=n, group="A" if n in (1, 4) else "B",
                                   inward=True)]


def line_rule(edge, order):
    nodes, weights = leggauss(order)
    t = (nodes + 1) / 2
    return (edge.start + t[:, None] * (edge.end - edge.start),
            weights * edge.length / 2)


def _vertical_span(polygon, x):
    ys = []
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        if min(a[0], b[0]) < x < max(a[0], b[0]):
            ys.append(a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]))
    if not ys:
        return None
    if len(ys) != 2:
        raise ValueError("Expected convex, nonoverlapping polygons")
    return min(ys), max(ys)


def area_rule(outer, holes, order, subdivisions=1):
    """Material-conforming quadrature: split at every polygon x vertex.

    Each vertical interval is either dielectric or one air hole. No quadrature
    point lies on an interface; weights have m² units. Uniform subintervals
    control sampling of the FEM field, independently of the geometry split.
    """
    outer = np.asarray(outer, float)
    holes = [np.asarray(h, float) for h in holes]
    xbreaks = np.unique(np.concatenate([p[:, 0] for p in [outer, *holes]]))
    nodes, weights = leggauss(order)
    points, areas, air = [], [], []
    for xa, xb in zip(xbreaks[:-1], xbreaks[1:]):
        if xb - xa < np.ptp(outer[:, 0]) * 1e-12:
            continue
        xsplits = np.linspace(xa, xb, subdivisions + 1)
        for x0, x1 in zip(xsplits[:-1], xsplits[1:]):
            for x, wx in zip((x0 + x1) / 2 + (x1 - x0) / 2 * nodes,
                             (x1 - x0) / 2 * weights):
                lower, upper = _vertical_span(outer, x)
                spans = sorted(s for h in holes if (s := _vertical_span(h, x)))
                intervals = []
                cursor = lower
                for lo, hi in spans:
                    if lo < cursor - 1e-15 or hi > upper + 1e-15:
                        raise ValueError("Hole overlap or hole outside the unit cell")
                    if lo > cursor:
                        intervals.append((cursor, lo, False))
                    intervals.append((lo, hi, True))
                    cursor = hi
                if cursor < upper:
                    intervals.append((cursor, upper, False))
                for lo, hi, is_air in intervals:
                    ysplits = np.linspace(lo, hi, subdivisions + 1)
                    for y0, y1 in zip(ysplits[:-1], ysplits[1:]):
                        ys = (y0 + y1) / 2 + (y1 - y0) / 2 * nodes
                        points.extend(zip(np.full(order, x), ys))
                        areas.extend(wx * (y1 - y0) / 2 * weights)
                        air.extend([is_air] * order)
    return np.asarray(points), np.asarray(areas), np.asarray(air, bool)


def edge_integrals(edges, sample_hz, order):
    rows = []
    for e in edges:
        xy, w = line_rule(e, order)
        hz = np.asarray(sample_hz(xy), complex)
        if hz.shape != w.shape or not np.isfinite(hz).all():
            raise ValueError("Missing/nonfinite Hz samples")
        integral = np.dot(w, hz)
        rows.append({"hole": e.hole, "edge": e.edge, "group": e.group,
                     "length_m": e.length, "nx": float(e.normal[0]),
                     "ny": float(e.normal[1]), "integral_hz": integral,
                     "mean_hz": integral / e.length,
                     "qx": integral * e.normal[0],
                     "qy": integral * e.normal[1]})
    return rows


def group_totals(rows):
    result = {}
    for g in ("A", "B"):
        selected = [r for r in rows if r["group"] == g]
        if len(selected) != (6 if g == "A" else 12):
            raise ValueError("Incomplete A/B edge set")
        for axis in ("x", "y"):
            result[f"Q{g}{axis}"] = sum(r[f"q{axis}"] for r in selected)
    result["Qx"] = result["QAx"] + result["QBx"]
    result["Qy"] = result["QAy"] + result["QBy"]
    result["Dq"] = sum(abs(r["qx"]) + abs(r["qy"]) for r in rows)
    den = abs(result["QAx"]) + abs(result["QBx"])
    result["cancellation_ratio"] = abs(result["Qx"]) / den if den else float("nan")
    result["AB_phase_rad"] = float(np.angle(result["QAx"] * np.conj(result["QBx"])))
    return result


def maxwell_terms(qx, qy, outer_x, outer_y, dz_hx, dz_hy,
                  omega, area, epsilon_air, epsilon_slab):
    """e^(-iωt), n: dielectric→air. Derivative inputs include 1/epsilon.

    The returned quantities all have V/m units with unnormalised H in A/m.
    Outer and dz terms are integrals, not averages. No term is discarded.
    """
    p = 1j / (omega * EPS0 * area)
    jump = 1 / epsilon_air - 1 / epsilon_slab
    lx, ly = -p * jump * qy, p * jump * qx
    ox, oy = p * outer_y, -p * outer_x
    zx, zy = -p * dz_hy, p * dz_hx
    return {"Lx": lx, "Ly": ly, "Ox": ox, "Oy": oy, "Zx": zx, "Zy": zy,
            "predicted_Ex": lx + ox + zx, "predicted_Ey": ly + oy + zy,
            "prefactor": p, "jump": jump}


def phase_and_norm(hz, weights, reference=None):
    """One factor for all E/H. Hz has unit RMS on the reference section."""
    hz = np.asarray(hz, complex)
    rms = np.sqrt(np.dot(weights, abs(hz) ** 2) / np.sum(weights))
    if not np.isfinite(rms) or rms <= 0:
        raise ValueError("Zero or nonfinite cell Hz norm")
    overlap = np.dot(weights, np.conj(reference) * hz) if reference is not None else 0
    if reference is not None and abs(overlap) == 0:
        raise ValueError("Reference field is orthogonal to the selected mode")
    anchor = overlap if reference is not None else hz[np.argmax(abs(hz))]
    return np.exp(-1j * np.angle(anchor)) / rms


def parity_errors(h, mirror_x, mirror_y):
    scale = np.linalg.norm(h)
    if scale == 0:
        return {"odd_x_error": float("inf"), "even_y_error": float("inf")}
    return {"odd_x_error": float(np.linalg.norm(mirror_x + h) / scale),
            "even_y_error": float(np.linalg.norm(mirror_y - h) / scale)}


def complex_columns(record):
    result = {}
    for key, value in record.items():
        if isinstance(value, (complex, np.complexfloating)):
            result[key + "_re"] = float(np.real(value))
            result[key + "_im"] = float(np.imag(value))
        elif isinstance(value, np.generic):
            result[key] = value.item()
        else:
            result[key] = value
    return result


def length_only_prediction(reference_rows, target_edges):
    """Keep each complex edge-averaged Hz; use the actual target lengths/normals."""
    lookup = {(r["hole"], r["edge"]): r for r in reference_rows}
    rows = []
    for e in target_edges:
        ref = lookup[(e.hole, e.edge)]
        h = ref["mean_hz"]
        rows.append({"hole": e.hole, "edge": e.edge, "group": e.group,
                     "qx": h * e.length * e.normal[0],
                     "qy": h * e.length * e.normal[1]})
    return rows


def segment_minimum(z0, z1, value0, value1):
    """Complex linear interpolation diagnostic, never a solved radiation zero."""
    dv = value1 - value0
    if abs(dv) == 0:
        return None
    t = float(np.clip(-np.real(np.conj(dv) * value0) / abs(dv)**2, 0, 1))
    value = value0 + t * dv
    return {"zeta_estimate": z0 + t * (z1 - z0), "complex_residual": abs(value),
            "value_estimate": value, "interior": 0 < t < 1,
            "real_sign_change": value0.real * value1.real <= 0,
            "method": "linear_complex_segment_not_an_eigensolve"}

def signed_indicator(fa, fb):
    """Phase/scale invariant A-directed projection and cancellation diagnostics."""
    if not np.isfinite([fa, fb]).all() or abs(fa) == 0:
        raise ValueError("Finite integrals and a nonzero A reference are required")
    scale = max(abs(fa), abs(fb))
    a, b = fa / scale, fb / scale
    projection = ((a + b) * np.conj(a)) / (abs(a) * (abs(a) + abs(b)))
    s, q = float(projection.real), float(projection.imag)
    return {'s': s, 'inv_s': 1 / s if s != 0 else float('nan'),
            'q_perp': q, 'r': float(abs(a + b) / (abs(a) + abs(b)))}
