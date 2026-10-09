# Explicit Parametrization: Equilateral Tri2 inside fixed Equilateral Tri1 (Check-based)

This document describes an **explicit** parametrization for positioning an equilateral triangle **Tri2** inside a fixed equilateral triangle **Tri1**.

Unlike the *feasible-by-construction* method, this approach allows arbitrary input parameters for position and size, generating the geometry first and then **validating** if the placement satisfies the containment and clearance constraints ($d$).

---

## Part 1 — Explanation

### 1) Tri1 Geometry (Fixed Outer)

Tri1 is defined by a scaling parameter $a$.
-   **Height**: $H_{1} = a/2$.
-   **Side Length**: $L_{1} = a / \sqrt{3}$.
-   **Apex**: At the origin $(0,0)$.
-   **Orientation**: Points to the left. The two legs extending from the apex form angles of $\pm 30^\circ$ with the positive $x$-axis. The base is a vertical line on the right.

Vertices (CCW):
1.  **Apex**: $A = (0, 0)$
2.  **Bottom-Right**: $B = ( H_1, -L_1/2 )$
3.  **Top-Right**: $C = ( H_1, +L_1/2 )$

Tri1 Interior inequalities:
-   (upper) $y \le x / \sqrt{3}$
-   (lower) $y \ge -x / \sqrt{3}$
-   (right) $x \le a/2$

---

### 2) Tri2 Geometry (Inner Hole)

Tri2 is an equilateral triangle defined by side length $b$.
-   **Side Length**: $L_2 = b$.
-   **Height**: $h_{2} = b\sqrt{3}/2$.
-   **Radius** (centroid to vertex): $R_{2} = b / \sqrt{3}$.

**Orientation**:
Initially (before rotation), Tri2 is aligned with Tri1 (pointing left).
Its vertices relative to its centroid are:
-   $V_0 = (-R_{2}, 0)$  (Left vertex)
-   $V_1 = (R_{2}/2, -b/2)$ (Bottom-right)
-   $V_2 = (R_{2}/2, +b/2)$ (Top-right)

This shape is then rotated by an angle $\phi$ around its centroid.

---

### 3) Positioning Tri2 via $(r, \theta)$

The centroid $(c_x, c_y)$ of Tri2 is parametrized by two variables: $r$ and $\theta$.
In this specific construction, $r$ represents the **x-coordinate projection**, not the Euclidean radius.

Mapping:
-   **x-coordinate**: $c_x = r$
-   **y-coordinate**: $c_y = r \cdot \tan(\theta)$

This places the centroid at a horizontal distance $r$ from the origin, lying on the ray at angle $\theta$ from the x-axis.

---

### 4) Validation (Clearance $d$)

Since the inputs $(r, \theta, b, \phi)$ are arbitrary, the resulting Tri2 might intersect Tri1 or likely be too close to the boundary.

A valid design requires:
1.  **Strict Containment**: Every vertex of Tri2 must be inside Tri1.
2.  **Clearance**: The minimum Euclidean distance between the sets of line segments (Tri1 edges vs Tri2 edges) must be $\ge d$.

---

## Part 2 — Implementation Recipe

### Inputs
-   $a > 0$ (scales Tri1)
-   $r$ (x-position of centroid)
-   $\theta \in [-30^\circ, 30^\circ]$ (angle of position ray)
-   $b > 0$ (Tri2 side length)
-   $\phi$ (Tri2 rotation in degrees)
-   $d \ge 0$ (required clearance threshold)

### Algorithm

1.  **Compute Tri1 Vertices**
    -   $H_1 = a / 2$
    -   $L_1 = a / \sqrt{3}$
    -   $V_{\text{out}} = \{ (0,0), \; (H_1, -L_1/2), \; (H_1, L_1/2) \}$

2.  **Compute Tri2 Center**
    -   $c_x = r$
    -   $c_y = r \cdot \tan(\theta_{\text{rad}})$

3.  **Compute Tri2 Local Vertices (aligned)**
    -   $R_2 = b / \sqrt{3}$
    -   $v_0 = (-R_2, 0)$
    -   $v_1 = (R_2/2, -b/2)$
    -   $v_2 = (R_2/2, +b/2)$

4.  **Rotate and Translate Tri2**
    For each local vertex $v$ in $\{v_0, v_1, v_2\}$:
    -   Rotate by $\phi$:
        $$ x_{\text{rot}} = v_x \cos\phi - v_y \sin\phi $$
        $$ y_{\text{rot}} = v_x \sin\phi + v_y \cos\phi $$
    -   Translate to center:
        $$ V_{\text{in}} = ( x_{\text{rot}} + c_x, \; y_{\text{rot}} + c_y ) $$

5.  **Calculate Minimum Distance**
    Compute the shortest distance $D_{\min}$ between the boundary of Tri1 and the boundary of Tri2 using segment-to-segment distance checks.
    -   Check distances from all Tri2 vertices to all Tri1 edges.
    -   Check distances from all Tri1 vertices to all Tri2 edges.

6.  **Validate**
    -   **is_inside**: Check if all $V_{\text{in}}$ are strictly inside Tri1 bounds.
    -   **is_valid**: `is_inside` AND (`min_dist` $\ge d$).