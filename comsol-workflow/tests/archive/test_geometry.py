from tests._paths import output_root
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

from comsol_workflow.geometry_utils import (
    rotate_points,
    create_outer_triangle,
    create_inner_triangle,
    point_to_segment_distance,
    point_in_triangle,
    get_min_dist,
    create_single_sector,
    create_hexagon_design,
    visualize_hexagon_design
)

def test_rotation_logic_static():
    print("\n=== TEST: Static Rotation Logic ===")
    p = np.array([[1.0, 0.0]])
    
    # 1. Rotate 90 degrees
    p_90 = rotate_points(p, 90.0)
    expected_90 = np.array([[0.0, 1.0]])
    assert np.allclose(p_90, expected_90, atol=1e-7), f"Rotate 90 failed. Got {p_90}"
    print("[Pass] Rotate (1,0) by 90 deg -> (0,1)")

    # 2. Rotate 180 degrees
    p_180 = rotate_points(p, 180.0)
    expected_180 = np.array([[-1.0, 0.0]])
    assert np.allclose(p_180, expected_180, atol=1e-7), f"Rotate 180 failed. Got {p_180}"
    print("[Pass] Rotate (1,0) by 180 deg -> (-1,0)")

def test_rotation_random_properties():
    print("\n=== TEST: Random Rotation Properties ===")
    # Generate 100 random points
    points = np.random.rand(100, 2) * 10 - 5  # Range [-5, 5]
    angle = np.random.uniform(0, 360)
    
    rotated_points = rotate_points(points, angle)
    
    # Property 1: Conservation of distance from origin (norm)
    norms_orig = np.linalg.norm(points, axis=1)
    norms_rot = np.linalg.norm(rotated_points, axis=1)
    
    if np.allclose(norms_orig, norms_rot, atol=1e-7):
        print(f"[Pass] Norms conserved for random angle {angle:.2f}")
    else:
        print(f"[Fail] Norms changed after rotation!")

    # Property 2: 360 degree rotation returns to start
    full_circle = rotate_points(points, 360.0)
    if np.allclose(points, full_circle, atol=1e-7):
        print("[Pass] 360 degree rotation returns original points")
    else:
        print("[Fail] 360 degree idenity check failed")

def test_triangle_geometry_random():
    print("\n=== TEST: Triangle Geometry Scaling (Random) ===")
    a_val = np.random.uniform(1.0, 10.0)
    tri_outer = create_outer_triangle(a_val)
    
    # Expected height for equilateral triangle side 'side' where 'a' usually relates to height or side
    # Based on existing code logic: create_outer_triangle(a) -> height = a/2
    # So side length = height / (sqrt(3)/2) = (a/2) * (2/sqrt(3)) = a/sqrt(3)
    
    height_x = tri_outer[1, 0] # Base x-coord
    edge_y_top = tri_outer[2, 1]
    edge_y_bot = tri_outer[1, 1]
    
    expected_height = a_val / 2.0
    if np.isclose(height_x, expected_height):
        print(f"[Pass] Height scales correctly with a={a_val:.2f}")
    else:
        print(f"[Fail] Height mismatch. Expected {expected_height}, Got {height_x}")

    # Check symmetry of base
    if np.isclose(abs(edge_y_top), abs(edge_y_bot)):
        print("[Pass] Base is symmetric around x-axis")
    else:
        print("[Fail] Base asymmetry detected")

def test_point_checks_random():
    print("\n=== TEST: Point Containment Random Sampling ===")
    # Triangle vertices: (0,0), (2,0), (0,2)
    tri = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
    
    # Generate points inside valid bounds [0,2]x[0,2]
    pts = np.random.rand(50, 2) * 2.5
    
    passed = True
    for p in pts:
        # Mathematical check: x > 0, y > 0, x + y < 2
        math_check = (p[0] >= 0) and (p[1] >= 0) and (p[0] + p[1] <= 2.0)
        func_check = point_in_triangle(p, tri)
        
        # Allow small epsilon for boundary cases
        if not np.isclose(p[0]+p[1], 2.0, atol=1e-5): 
            if math_check != func_check:
                passed = False
                print(f"[Fail] Point {p} -> Math: {math_check}, Func: {func_check}")
    
    if passed:
        print("[Pass] 50 random points validated against analytical constraints")

def test_full_flow_visualization():
    print("\n=== TEST: Full Flow & Visualization ===")
    
    # 1. Valid Configuration (Safe)
    a = 0.82
    # Params: radius(r), theta, size(b), rotation(phi)
    # Ensure r < height is satisfied easily
    valid_hole_params = [
        [0.150, 5.0, 0.100, 0.0],
        [0.175, 0.0, 0.125, 10.0],
        [0.200, 0.0, 0.150, -30.0],
        [0.225, -5.0, 0.175, 0.0],
        [0.250, 0.0, 0.200, 0.0],
        [0.275, 0.0, 0.225, 0.0],
    ]
    
    print("Generating Valid Design...")
    hex_v, holes_v_valid, info_valid = create_hexagon_design(a, valid_hole_params)
    
    # Check minimum distance
    min_dists_valid = [x['min_dist'] for x in info_valid]
    print(f"Valid Design Min Distances: {min_dists_valid}")
    if all(d > 0 for d in min_dists_valid):
        print("[Pass] Valid design has positive distances (No collision)")
    else:
        print("[Fail] Valid design detected collision incorrectly")

    visualize_hexagon_design(
        hex_v, holes_v_valid, info_valid, 
        filename=str(output_root() / 'test_design_valid.png'),
        annotate=False
    )
    visualize_hexagon_design(
        hex_v, holes_v_valid, info_valid, 
        filename=str(output_root() / 'test_design_valid_annotated.png'),
        annotate=True
    )

    # 2. Invalid Configuration (Collision/Outside)
    print("\nGenerating Invalid Design (Outside/Collision)...")
    invalid_hole_params = [
        [0.050, -5.0, 0.200, 0.0],
        [0.410, 0.0, 0.125, 10.0],
        [0.200, 0.0, 0.250, -30.0],
        [0.225, 25.0, 0.175, 0.0],
        [0.600, 0.0, 0.200, 0.0],
        [0.275, 0.0, 0.800, 0.0],
    ]
    
    hex_v, holes_v_invalid, info_invalid = create_hexagon_design(a, invalid_hole_params)
    
    min_dists_invalid = [x['min_dist'] for x in info_invalid]
    print(f"Invalid Design Min Distances: {min_dists_invalid}")
    
    if any(d < 0 for d in min_dists_invalid):
        print("[Pass] Invalid design correctly reports negative distance")
    else:
        print("[Fail] Invalid design reported positive distance")

    visualize_hexagon_design(
        hex_v, holes_v_invalid, info_invalid, 
        filename=str(output_root() / 'test_design_invalid.png'),
        annotate=False
    )
    visualize_hexagon_design(
        hex_v, holes_v_invalid, info_invalid, 
        filename=str(output_root() / 'test_design_invalid_annotated.png'),
        annotate=True
    )


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    
    test_rotation_logic_static()
    test_rotation_random_properties()
    test_triangle_geometry_random()
    test_point_checks_random()
    test_full_flow_visualization()
    print("\nAll Tests Completed.")
