from __future__ import annotations

import numpy as np


def rounded_vertex(point: np.ndarray, decimals: int = 10) -> tuple[float, float]:
    return tuple(np.round(point, decimals))


def polygon_area(vertices: np.ndarray) -> float:
    vertices = np.asarray(vertices, dtype=float)
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon_min_width(vertices: np.ndarray) -> float:
    vertices = np.asarray(vertices, dtype=float)
    if vertices.shape[0] < 3:
        return 0.0

    min_width = float("inf")
    for idx, start in enumerate(vertices):
        end = vertices[(idx + 1) % vertices.shape[0]]
        edge = end - start
        length = float(np.linalg.norm(edge))
        if length < 1e-12:
            return 0.0
        normal = np.array([-edge[1], edge[0]], dtype=float) / length
        projected = vertices @ normal
        min_width = min(min_width, float(projected.max() - projected.min()))
    return min_width


def exterior_boundary_segments(polygons: list[np.ndarray]) -> list[tuple[np.ndarray, np.ndarray]]:
    edge_counts: dict[tuple[tuple[float, float], tuple[float, float]], int] = {}
    vertex_values: dict[tuple[float, float], np.ndarray] = {}
    for polygon in polygons:
        vertices = [rounded_vertex(vertex) for vertex in polygon]
        for key, vertex in zip(vertices, polygon):
            vertex_values[key] = np.asarray(vertex, dtype=float)
        for idx, start in enumerate(vertices):
            end = vertices[(idx + 1) % len(vertices)]
            edge = tuple(sorted([start, end]))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

    return [(vertex_values[start], vertex_values[end]) for (start, end), count in edge_counts.items() if count == 1]


def ordered_boundary_loop(segments: list[tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    vertices: dict[tuple[float, float], np.ndarray] = {}
    adjacency: dict[tuple[float, float], list[tuple[float, float]]] = {}
    for start, end in segments:
        start_key = rounded_vertex(start)
        end_key = rounded_vertex(end)
        vertices[start_key] = start
        vertices[end_key] = end
        adjacency.setdefault(start_key, []).append(end_key)
        adjacency.setdefault(end_key, []).append(start_key)

    start = min(adjacency, key=lambda key: (vertices[key][0], vertices[key][1]))
    loop = [start]
    previous = None
    current = start
    while True:
        neighbors = adjacency[current]
        candidates = neighbors if previous is None else [key for key in neighbors if key != previous]
        if not candidates:
            break
        next_key = candidates[0]
        if next_key == start:
            break
        loop.append(next_key)
        previous, current = current, next_key

    ordered = np.array([vertices[key] for key in loop], dtype=float)
    if polygon_area(ordered) < 0.0:
        ordered = ordered[::-1]
    return ordered


def line_intersection(point_a: np.ndarray, direction_a: np.ndarray, point_b: np.ndarray, direction_b: np.ndarray) -> np.ndarray:
    matrix = np.column_stack([direction_a, -direction_b])
    rhs = point_b - point_a
    det = float(np.linalg.det(matrix))
    if abs(det) < 1e-12:
        return point_a
    t, _ = np.linalg.solve(matrix, rhs)
    return point_a + t * direction_a


def offset_polygon(vertices: np.ndarray, distance: float) -> np.ndarray:
    vertices = np.asarray(vertices, dtype=float)
    if polygon_area(vertices) < 0.0:
        vertices = vertices[::-1]

    shifted = []
    for idx in range(len(vertices)):
        prev_point = vertices[(idx - 1) % len(vertices)]
        point = vertices[idx]
        next_point = vertices[(idx + 1) % len(vertices)]

        prev_direction = point - prev_point
        next_direction = next_point - point
        prev_direction = prev_direction / np.linalg.norm(prev_direction)
        next_direction = next_direction / np.linalg.norm(next_direction)

        prev_normal = np.array([prev_direction[1], -prev_direction[0]], dtype=float)
        next_normal = np.array([next_direction[1], -next_direction[0]], dtype=float)
        shifted.append(line_intersection(point + distance * prev_normal, prev_direction, point + distance * next_normal, next_direction))
    return np.asarray(shifted, dtype=float)


def cross_2d(vector_a: np.ndarray, vector_b: np.ndarray) -> float:
    return float(vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0])


def point_on_segment(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> bool:
    if abs(cross_2d(point - start, end - start)) > 1e-10:
        return False
    return float(np.dot(point - start, end - start)) >= -1e-10 and float(np.dot(point - end, start - end)) >= -1e-10


def point_in_polygon_winding(point: np.ndarray, polygon: np.ndarray) -> bool:
    x, y = point
    winding = 0
    for idx, start in enumerate(polygon):
        end = polygon[(idx + 1) % len(polygon)]
        if start[1] <= y:
            if end[1] > y and cross_2d(end - start, point - start) > 1e-10:
                winding += 1
        else:
            if end[1] <= y and cross_2d(end - start, point - start) < -1e-10:
                winding -= 1
    return winding != 0


def point_inside_or_on_polygon(point: np.ndarray, polygon: np.ndarray) -> bool:
    if point_in_polygon_winding(point, polygon):
        return True
    return any(point_on_segment(point, polygon[idx], polygon[(idx + 1) % len(polygon)]) for idx in range(len(polygon)))


def point_strictly_inside_polygon(point: np.ndarray, polygon: np.ndarray) -> bool:
    if any(point_on_segment(point, polygon[idx], polygon[(idx + 1) % len(polygon)]) for idx in range(len(polygon))):
        return False
    return point_in_polygon_winding(point, polygon)


def segment_intersection_parameter(
    start_a: np.ndarray,
    end_a: np.ndarray,
    start_b: np.ndarray,
    end_b: np.ndarray,
) -> tuple[float, float, np.ndarray] | None:
    direction_a = end_a - start_a
    direction_b = end_b - start_b
    matrix = np.column_stack([direction_a, -direction_b])
    det = float(np.linalg.det(matrix))
    if abs(det) < 1e-12:
        return None
    t, u = np.linalg.solve(matrix, start_b - start_a)
    if -1e-10 <= t <= 1.0 + 1e-10 and -1e-10 <= u <= 1.0 + 1e-10:
        t = min(1.0, max(0.0, float(t)))
        u = min(1.0, max(0.0, float(u)))
        return t, u, start_a + t * direction_a
    return None


def point_segment_parameter(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float | None:
    direction = end - start
    length_sq = float(np.dot(direction, direction))
    if length_sq < 1e-14:
        return 0.0 if np.linalg.norm(point - start) <= 1e-10 else None
    if abs(cross_2d(point - start, direction)) > 1e-10:
        return None
    t = float(np.dot(point - start, direction) / length_sq)
    if -1e-10 <= t <= 1.0 + 1e-10:
        return min(1.0, max(0.0, t))
    return None


def segment_split_points(edge_start: np.ndarray, edge_end: np.ndarray, cutter_start: np.ndarray, cutter_end: np.ndarray) -> list[tuple[float, np.ndarray]]:
    split_points = []
    result = segment_intersection_parameter(edge_start, edge_end, cutter_start, cutter_end)
    if result is not None:
        t, _, point = result
        split_points.append((t, point))

    for point in [cutter_start, cutter_end]:
        t = point_segment_parameter(point, edge_start, edge_end)
        if t is not None:
            split_points.append((t, edge_start + t * (edge_end - edge_start)))

    return split_points


def segments_intersect(start_a: np.ndarray, end_a: np.ndarray, start_b: np.ndarray, end_b: np.ndarray) -> bool:
    if segment_intersection_parameter(start_a, end_a, start_b, end_b) is not None:
        return True
    return (
        point_on_segment(start_a, start_b, end_b)
        or point_on_segment(end_a, start_b, end_b)
        or point_on_segment(start_b, start_a, end_a)
        or point_on_segment(end_b, start_a, end_a)
    )


def segments_properly_intersect(start_a: np.ndarray, end_a: np.ndarray, start_b: np.ndarray, end_b: np.ndarray) -> bool:
    direction_a = end_a - start_a
    direction_b = end_b - start_b
    orient_c = cross_2d(direction_a, start_b - start_a)
    orient_d = cross_2d(direction_a, end_b - start_a)
    orient_a = cross_2d(direction_b, start_a - start_b)
    orient_b = cross_2d(direction_b, end_a - start_b)
    return orient_c * orient_d < -1e-20 and orient_a * orient_b < -1e-20


def split_edge_by_polygon(edge_start: np.ndarray, edge_end: np.ndarray, polygon: np.ndarray) -> list[np.ndarray]:
    splits = [(0.0, edge_start), (1.0, edge_end)]
    for idx in range(len(polygon)):
        splits.extend(segment_split_points(edge_start, edge_end, polygon[idx], polygon[(idx + 1) % len(polygon)]))

    deduped = []
    for t, point in sorted(splits, key=lambda item: item[0]):
        if not deduped or np.linalg.norm(point - deduped[-1][1]) > 1e-10:
            deduped.append((t, point))
    return [point for _, point in deduped]


def polygons_intersect(polygon_a: np.ndarray, polygon_b: np.ndarray) -> bool:
    if any(point_inside_or_on_polygon(point, polygon_a) for point in polygon_b):
        return True
    if any(point_inside_or_on_polygon(point, polygon_b) for point in polygon_a):
        return True
    for idx_a in range(len(polygon_a)):
        start_a = polygon_a[idx_a]
        end_a = polygon_a[(idx_a + 1) % len(polygon_a)]
        for idx_b in range(len(polygon_b)):
            if segment_intersection_parameter(start_a, end_a, polygon_b[idx_b], polygon_b[(idx_b + 1) % len(polygon_b)]) is not None:
                return True
    return False


def segment_midpoint(start: np.ndarray, end: np.ndarray) -> np.ndarray:
    return 0.5 * (start + end)


def add_segment(segments: list[tuple[np.ndarray, np.ndarray]], start: np.ndarray, end: np.ndarray) -> None:
    if np.linalg.norm(end - start) > 1e-10:
        segments.append((np.asarray(start, dtype=float), np.asarray(end, dtype=float)))


def directed_loops_from_segments(segments: list[tuple[np.ndarray, np.ndarray]]) -> list[np.ndarray]:
    edges = [(rounded_vertex(start), rounded_vertex(end), start, end) for start, end in segments if np.linalg.norm(end - start) > 1e-10]
    starts: dict[tuple[float, float], list[int]] = {}
    for idx, (start_key, _, _, _) in enumerate(edges):
        starts.setdefault(start_key, []).append(idx)

    used = set()
    loops = []
    for edge_idx in range(len(edges)):
        if edge_idx in used:
            continue
        start_key, end_key, start, end = edges[edge_idx]
        used.add(edge_idx)
        loop_points = [start, end]
        current_key = end_key

        while current_key != start_key:
            candidates = [idx for idx in starts.get(current_key, []) if idx not in used]
            if not candidates:
                break
            next_idx = candidates[0]
            _, next_end_key, _, next_end = edges[next_idx]
            used.add(next_idx)
            loop_points.append(next_end)
            current_key = next_end_key

        if current_key == start_key and len(loop_points) >= 4:
            loop = np.asarray(loop_points[:-1], dtype=float)
            area = polygon_area(loop)
            if abs(area) > 1e-12:
                loops.append(loop if area > 0.0 else loop[::-1])
    return loops


def clip_polygon_to_outside_region(subject: np.ndarray, boundary: np.ndarray) -> list[np.ndarray]:
    subject = np.asarray(subject, dtype=float)
    boundary = np.asarray(boundary, dtype=float)
    if polygon_area(subject) < 0.0:
        subject = subject[::-1]
    if polygon_area(boundary) < 0.0:
        boundary = boundary[::-1]

    subject_vertices_inside = [point_inside_or_on_polygon(point, boundary) for point in subject]
    if all(subject_vertices_inside) and not polygons_intersect(subject, boundary):
        return []
    if not any(subject_vertices_inside) and not polygons_intersect(subject, boundary):
        return [subject]

    segments: list[tuple[np.ndarray, np.ndarray]] = []

    for idx in range(len(subject)):
        points = split_edge_by_polygon(subject[idx], subject[(idx + 1) % len(subject)], boundary)
        for start, end in zip(points[:-1], points[1:]):
            if not point_inside_or_on_polygon(segment_midpoint(start, end), boundary):
                add_segment(segments, start, end)

    for idx in range(len(boundary)):
        points = split_edge_by_polygon(boundary[idx], boundary[(idx + 1) % len(boundary)], subject)
        for start, end in zip(points[:-1], points[1:]):
            if point_inside_or_on_polygon(segment_midpoint(start, end), subject):
                add_segment(segments, end, start)

    return directed_loops_from_segments(segments)


def clipped_pieces_contain_point(pieces: list[np.ndarray], point: np.ndarray, tolerance: float = 1e-8) -> bool:
    return any(np.linalg.norm(piece - point, axis=1).min(initial=np.inf) <= tolerance for piece in pieces if len(piece) > 0)


def validate_clipped_hole(original: np.ndarray, clipped: list[np.ndarray], boundary: np.ndarray) -> None:
    for vertex in original:
        is_outside = not point_inside_or_on_polygon(vertex, boundary)
        is_inside = point_strictly_inside_polygon(vertex, boundary)
        is_retained = clipped_pieces_contain_point(clipped, vertex)
        if is_outside and not is_retained:
            raise AssertionError(f"Outside triangle vertex was not retained: {vertex.tolist()}")
        if is_inside and is_retained:
            raise AssertionError(f"Inside triangle vertex was retained: {vertex.tolist()}")

    for piece in clipped:
        for idx, start in enumerate(piece):
            end = piece[(idx + 1) % len(piece)]
            midpoint = segment_midpoint(start, end)
            if point_strictly_inside_polygon(midpoint, boundary):
                raise AssertionError(f"Clipped triangle edge still lies inside boundary: {midpoint.tolist()}")


def polygon_fully_inside_region(polygon: np.ndarray, boundary: np.ndarray) -> bool:
    if not all(point_inside_or_on_polygon(point, boundary) for point in polygon):
        return False
    for idx in range(len(polygon)):
        start = polygon[idx]
        end = polygon[(idx + 1) % len(polygon)]
        for boundary_idx in range(len(boundary)):
            if segments_properly_intersect(start, end, boundary[boundary_idx], boundary[(boundary_idx + 1) % len(boundary)]):
                return False
    return True


def clip_polygon_to_convex_region(
    subject: np.ndarray,
    boundary: np.ndarray,
    *,
    tolerance: float = 1e-10,
) -> np.ndarray:
    """Return the part of a polygon lying inside a convex CCW boundary.

    The implementation is the Sutherland-Hodgman algorithm.  It is used by
    symmetry-reduced finite models to clip holes that cross a mirror plane.
    Empty intersections are returned with shape ``(0, 2)``.
    """
    subject = np.asarray(subject, dtype=float)
    boundary = np.asarray(boundary, dtype=float)
    if subject.ndim != 2 or subject.shape[1] != 2 or len(subject) < 3:
        raise ValueError("subject must have shape (N, 2) with N >= 3")
    if boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 3:
        raise ValueError("boundary must have shape (N, 2) with N >= 3")
    if polygon_area(boundary) <= tolerance:
        raise ValueError("boundary must be a non-degenerate counter-clockwise polygon")

    output = subject.copy()
    if polygon_area(output) < 0.0:
        output = output[::-1]

    def inside(point: np.ndarray, edge_start: np.ndarray, edge_end: np.ndarray) -> bool:
        return cross_2d(edge_end - edge_start, point - edge_start) >= -tolerance

    def intersection(
        segment_start: np.ndarray,
        segment_end: np.ndarray,
        edge_start: np.ndarray,
        edge_end: np.ndarray,
    ) -> np.ndarray:
        segment = segment_end - segment_start
        edge = edge_end - edge_start
        denominator = cross_2d(segment, edge)
        if abs(denominator) <= tolerance:
            return 0.5 * (segment_start + segment_end)
        parameter = cross_2d(edge_start - segment_start, edge) / denominator
        return segment_start + parameter * segment

    for edge_idx, edge_start in enumerate(boundary):
        if len(output) == 0:
            break
        edge_end = boundary[(edge_idx + 1) % len(boundary)]
        input_vertices = output
        clipped: list[np.ndarray] = []
        previous = input_vertices[-1]
        previous_inside = inside(previous, edge_start, edge_end)
        for current in input_vertices:
            current_inside = inside(current, edge_start, edge_end)
            if current_inside:
                if not previous_inside:
                    clipped.append(intersection(previous, current, edge_start, edge_end))
                clipped.append(current)
            elif previous_inside:
                clipped.append(intersection(previous, current, edge_start, edge_end))
            previous = current
            previous_inside = current_inside

        deduplicated: list[np.ndarray] = []
        for point in clipped:
            if not deduplicated or np.linalg.norm(point - deduplicated[-1]) > tolerance:
                deduplicated.append(np.asarray(point, dtype=float))
        if len(deduplicated) > 1 and np.linalg.norm(deduplicated[0] - deduplicated[-1]) <= tolerance:
            deduplicated.pop()
        output = np.asarray(deduplicated, dtype=float).reshape(-1, 2)

    if len(output) < 3 or abs(polygon_area(output)) <= tolerance:
        return np.empty((0, 2), dtype=float)
    if polygon_area(output) < 0.0:
        output = output[::-1]
    return output


def cell_has_outside_region_part(polygon: np.ndarray, boundary: np.ndarray) -> bool:
    return not polygon_fully_inside_region(polygon, boundary)
