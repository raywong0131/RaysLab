import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def rotate_points(points, angle_deg, origin=(0.0, 0.0)):
    """Rotate Nx2 points by angle_deg about origin."""
    pts = np.asarray(points, dtype=float)
    ox, oy = origin
    ang = np.radians(angle_deg)
    # R[2,2] x P[2,N] -> P_rot[2,N]
    R = np.array([[np.cos(ang), -np.sin(ang)],
                  [np.sin(ang),  np.cos(ang)]], dtype=float)
    # P[N,2] x R[2,2].T -> P_rot[N,2]
    return (pts - np.array([ox, oy], dtype=float)) @ R.T + np.array([ox, oy], dtype=float)


def create_outer_triangle(a):
    """
    equilateral triangle with apex at origin.
    Two edges from apex are at angles ±30° to x-axis.
    Given: height (apex -> opposite side) is a/2.

    Returns vertices as (3,2): [apex, outer_down, outer_up].
    """
    height = a / 2.0
    edge = a / np.sqrt(3.0)  # because height = edge*sqrt(3)/2 = a/2
    tri = np.array([
        [0.0, 0.0],
        [height, -edge / 2.0],
        [height,  edge / 2.0],
    ], dtype=float)
    return tri


def create_inner_triangle(b, center=(0.0, 0.0), rotation=0.0):
    """
    equilateral triangle defined by edge length b, centroid at `center`.
    Base orientation has one side vertical (parallel to outer triangle's opposite side),
    and the other two edges parallel to outer triangle's ±30° edges.
    Then rotate by `rotation` (deg) about its centroid.
    """
    h = b * np.sqrt(3.0) / 2.0

    # Start with same orientation as outer triangle (side opposite "apex" vertical), with centroid at origin.
    tri = np.array([
        [-2.0 * h / 3.0, 0.0],
        [ h / 3.0, -b / 2.0],
        [ h / 3.0,  b / 2.0],
    ], dtype=float)

    tri = rotate_points(tri, rotation, origin=(0.0, 0.0))

    return tri + np.array(center, dtype=float)


def point_to_segment_distance(point, seg_start, seg_end):
    """Minimum distance from a point to a line segment, returning distance and both points."""
    p = np.asarray(point, dtype=float)
    a = np.asarray(seg_start, dtype=float)
    b = np.asarray(seg_end, dtype=float)

    ab = b - a
    if np.allclose(ab, 0.0):
        return float(np.linalg.norm(p - a)), p, a

    t = float(np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0.0, 1.0))
    nearest = a + t * ab
    return float(np.linalg.norm(p - nearest)), p, nearest


def point_in_triangle(point, triangle, eps=1e-10):
    """Check if point is inside (or on edge) of triangle using barycentric coordinates."""
    x, y = point
    x1, y1 = triangle[0]
    x2, y2 = triangle[1]
    x3, y3 = triangle[2]

    denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    if abs(denom) < eps:
        return False

    a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denom
    b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denom
    c = 1.0 - a - b

    return (a >= -eps) and (b >= -eps) and (c >= -eps)


def get_min_dist(tri_inner, tri_outer):
    """Minimum distance between two triangles, returning distance and closest points."""
    tri_inner = np.asarray(tri_inner, dtype=float)
    tri_outer = np.asarray(tri_outer, dtype=float)
    min_dist = float("inf")
    closest_point_inner = None
    closest_point_outer = None

    # Check tri_inner vertices to tri_outer edges
    for i in range(3): # tri_inner vertices
        vertex_inner = tri_inner[i]
        is_inside = point_in_triangle(vertex_inner, tri_outer)
        vertex_min_dist_abs = float("inf")
        vertex_closest_point_inner = None
        vertex_closest_point_outer = None
        for j in range(3): # tri_outer edges
            dist_abs, pt_inner, pt_outer = point_to_segment_distance(
                vertex_inner, tri_outer[j], tri_outer[(j + 1) % 3]
            )

            if dist_abs < vertex_min_dist_abs:
                vertex_min_dist_abs = dist_abs
                vertex_closest_point_inner = pt_inner
                vertex_closest_point_outer = pt_outer

        vertex_min_dist = (1 if is_inside else -1) * vertex_min_dist_abs
        if vertex_min_dist < min_dist:
            min_dist = vertex_min_dist
            closest_point_inner = vertex_closest_point_inner
            closest_point_outer = vertex_closest_point_outer

    return min_dist, closest_point_inner, closest_point_outer


def create_single_sector(a, r, theta, b, phi):
    """
    Create one sector's geometry:
      - outer triangle with apex at origin and edges ±30° to x-axis
      - inner triangular hole centered at (r, r * tan(theta)), then rotated by phi about its centroid

    Returns: (hole_vertices, info)
    """
    tri_outer = create_outer_triangle(a)

    cx, cy = r, r * np.tan(np.radians(theta))
    tri_inner = create_inner_triangle(b, center=(cx, cy), rotation=phi)

    # signed distance from inner triangle vertices to outer triangle (negative if outside)
    min_dist, pt1, pt2 = get_min_dist(tri_inner, tri_outer)

    info = {
        "min_dist": min_dist,
        "closest_points": [pt1, pt2],
    }

    return tri_inner, tri_outer, info


def create_hexagon_design(a, hole_params):
    """
    Full design with 6 sectors rotated by 0, 60, ..., 300 degrees.

    Input: a, hole_params = [(r, theta, b, phi)] * 6  (theta is local sector angle in [-30, +30])
    Output: (hexagon_vertices (6,2), hole_vertices (6,3,2), sector_info)
    """
    if len(hole_params) != 6:
        raise ValueError("Need exactly 6 hole parameters")

    hole_vertices = []
    sector_info = []

    for i, (r, theta, b, phi) in enumerate(hole_params):
        hole_coord, _, info = create_single_sector(a, r, theta, b, phi)

        rot = 60.0 * i
        hole_coord_r = rotate_points(hole_coord, rot, origin=(0.0, 0.0))

        # Rotate the closest points as well
        pt1_r = rotate_points([info["closest_points"][0]], rot, origin=(0.0, 0.0))[0]
        pt2_r = rotate_points([info["closest_points"][1]], rot, origin=(0.0, 0.0))[0]
        rotated_closest = [pt1_r.tolist(), pt2_r.tolist()]

        hole_vertices.append(hole_coord_r)
        
        # Validation for this sector
        sector_result = {
            "sector_index": i,
            "rotation": rot,
            "min_dist": info["min_dist"],
            "closest_points": rotated_closest
        }
        sector_info.append(sector_result)

    # The outer boundary is a regular hexagon; corners are the "outer up" vertex of each sector.
    height = a / 2.0
    edge = a / np.sqrt(3.0)
    corner0 = np.array([height, edge / 2.0], dtype=float)  # angle +30°, radius=edge
    hexagon_vertices = np.vstack([rotate_points(corner0[None, :], 60.0 * i, origin=(0.0, 0.0))[0] for i in range(6)])

    return hexagon_vertices, hole_vertices, sector_info


def visualize_hexagon_design(
    hexagon_vertices,
    hole_vertices,
    sector_info,
    filename,
    annotate,
):
    """Visualize the hexagon design (hexagon outline + triangular holes)."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Draw hexagon boundary
    hex_patch = Polygon(hexagon_vertices, closed=True, fill=False, edgecolor="blue", linewidth=2, label="Hexagon")
    ax.add_patch(hex_patch)

    # Draw lines connecting opposite vertices (diagonals)
    for i in range(3):
        p1 = hexagon_vertices[i]
        p2 = hexagon_vertices[i + 3]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="blue", linestyle="--", linewidth=1, alpha=0.5)

    if annotate:
        for i, v in enumerate(hexagon_vertices):
            ax.plot(v[0], v[1], "bo", markersize=4)
            ax.text(v[0], v[1], f"  ({v[0]:.3f}, {v[1]:.3f})", fontsize=8, ha="left", va="center", color="blue")

    # Draw holes and min distance annotations
    for i, hole_coord in enumerate(hole_vertices):
        ax.add_patch(Polygon(hole_coord, closed=True, fill=False, edgecolor="red", linewidth=1.5))

        if annotate:
            min_dist = sector_info[i]["min_dist"]

            # Annotate triangle vertices
            for v in hole_coord:
                ax.plot(v[0], v[1], "ro", markersize=3)
                ax.text(v[0], v[1], f"  ({v[0]:.3f}, {v[1]:.3f})", fontsize=7, ha="left", va="center")

            # Draw line between closest points
            pt1, pt2 = sector_info[i]["closest_points"]
            ax.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'g--', linewidth=1, alpha=0.7)
            ax.plot(pt1[0], pt1[1], 'go', markersize=4)
            ax.plot(pt2[0], pt2[1], 'go', markersize=4)
            
            # Show min distance next to the line
            mid_x = (pt1[0] + pt2[0]) / 2
            mid_y = (pt1[1] + pt2[1]) / 2
            color = "green" if min_dist > 0 else "red"
            ax.text(mid_x, mid_y, f"{min_dist:.4f}", fontsize=8, ha="center", color=color, weight="bold")

    # Explicit limits + padding so geometry fills the frame consistently
    pts = [np.asarray(hexagon_vertices, float)]
    for t in hole_vertices:
        pts.append(np.asarray(t, float))
    P = np.vstack(pts)

    xmin, ymin = P.min(axis=0)
    xmax, ymax = P.max(axis=0)
    dx, dy = (xmax - xmin), (ymax - ymin)
    pad = 0.08 * max(dx, dy)

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)

    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    
    title = "Unit Cell Geometry"
    ax.set_title(title)

    plt.tight_layout()
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
