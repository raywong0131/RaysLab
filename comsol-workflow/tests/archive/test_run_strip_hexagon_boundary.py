import os
import sys

import numpy as np

from scripts.run_main import run_strip_1d


def _inside_hexagon_boundary(points: np.ndarray, a: float, tolerance: float = 1e-12) -> np.ndarray:
    root3 = np.sqrt(3.0)
    x = points[:, 0]
    y = points[:, 1]
    return (
        (np.abs(y) <= root3 * a / 4.0 + tolerance)
        & (np.abs(root3 * x + y) <= root3 * a / 2.0 + tolerance)
        & (np.abs(root3 * x - y) <= root3 * a / 2.0 + tolerance)
    )


def test_hexagon_boundary_vertices_have_left_right_tips():
    vertices = run_strip_1d.hexagon_boundary_vertices(2.0)
    expected = np.array([
        [1.0, 0.0],
        [0.5, np.sqrt(3.0) / 2.0],
        [-0.5, np.sqrt(3.0) / 2.0],
        [-1.0, 0.0],
        [-0.5, -np.sqrt(3.0) / 2.0],
        [0.5, -np.sqrt(3.0) / 2.0],
    ])

    assert np.allclose(vertices, expected)


def test_hexagon_boundary_size_contains_points_with_padding():
    points = np.array([
        [1.0, 0.0],
        [0.0, np.sqrt(3.0) / 2.0],
        [-1.0, 0.0],
        [0.0, -np.sqrt(3.0) / 2.0],
    ])

    a = run_strip_1d.hexagon_boundary_size_for_points(points, pad=0.25)
    vertices = run_strip_1d.hexagon_boundary_vertices(a)

    assert a > 2.0
    assert np.all(_inside_hexagon_boundary(points, a))
    assert np.all(_inside_hexagon_boundary(vertices, a))


if __name__ == "__main__":
    test_hexagon_boundary_vertices_have_left_right_tips()
    test_hexagon_boundary_size_contains_points_with_padding()
    print("All run_strip_1d hexagon boundary tests passed.")
