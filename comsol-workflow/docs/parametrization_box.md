# Feasible-by-construction parametrization: equilateral Tri2 inside fixed equilateral Tri1 (with clearance $d$)

This document gives a **box-constrained** parametrization for an equilateral triangle **Tri2** that is guaranteed to be:
1.  **Inside** a fixed equilateral triangle **Tri1**.
2.  At **minimum distance $\ge d$** from the boundary of Tri1 (a clearance / safety margin).

It also enforces **counter-clockwise (CCW)** vertex ordering for **both** triangles.

---

## Part 1 — Explanation (clear reasoning)

### 1) Coordinate frame and CCW ordering for Tri1 (fixed)

Let Tri1 have side length $L > 0$. You want:
-   Vertex $A$ at the origin $(0,0)$.
-   The edges $AB$ and $AC$ along the rays at angles $-30^\circ$ and $+30^\circ$ from the $x$-axis.
-   **CCW vertex order**.

Set Tri1 vertices (CCW):
-   $A = (0, 0)$
-   $B = ( \frac{\sqrt{3}}{2}L , -L/2 )$   (on the $-30^\circ$ ray)
-   $C = ( \frac{\sqrt{3}}{2}L , +L/2 )$   (on the $+30^\circ$ ray)

This ordering $A \to B \to C$ is CCW.

The three sides of Tri1 correspond to these lines:
-   Upper slanted side ($AC$): $y = x / \sqrt{3}$
-   Lower slanted side ($AB$): $y = -x / \sqrt{3}$
-   Vertical side ($BC$): $x = \frac{\sqrt{3}}{2}L$

So the *interior* of Tri1 is the set of points $(x,y)$ satisfying:
-   (upper)  $y \le \frac{1}{\sqrt{3}}x$
-   (lower)  $y \ge -\frac{1}{\sqrt{3}}x$
-   (right)  $x \le \frac{\sqrt{3}}{2}L$

---

### 2) Adding a minimum clearance $d$ to Tri1’s boundary

We want Tri2 not only inside Tri1, but also with **minimum distance $\ge d$** from Tri1’s boundary.

Because Tri1’s sides are straight lines, "distance $\ge d$ from a side" means the point must lie inside an **inset** half-plane obtained by shifting each side **inward** by distance $d$.

For these three sides, the inward-shifted inequalities become:

-   (upper, shifted down)
    $$ y \le \frac{1}{\sqrt{3}}x - \frac{2}{\sqrt{3}}d $$

-   (lower, shifted up)
    $$ y \ge -\frac{1}{\sqrt{3}}x + \frac{2}{\sqrt{3}}d $$

-   (right, shifted left)
    $$ x \le \frac{\sqrt{3}}{2}L - d $$

If $d$ is too large, the inset triangle becomes empty. The largest possible is Tri1’s inradius:
$$ d < \frac{\sqrt{3}}{6}L $$
(At equality, the inset collapses to a point).

From now on, whenever we say "Tri2 is feasible", we mean **Tri2 lies inside this inset region**, which guarantees the boundary-to-boundary distance between Tri2 and Tri1 is at least $d$.

---

### 3) Parametrizing Tri2 by centroid + size + rotation

Tri2 is equilateral, so it is completely described by:
-   $(p,q)$ = centroid (center) of Tri2
-   $\ell$ = side length of Tri2
-   $\theta$ = rotation about the $z$-axis through the centroid

Key equilateral geometry:
-   The distance from the centroid to each vertex is $r = \ell / \sqrt{3}$.

#### "Initially aligned with Tri1" at $\theta = 0$
At $\theta = 0$, Tri2 should "point" like Tri1: one vertex to the left, two to the right (one up, one down).

Define base vertex directions (angles) around the centroid in **CCW order**:
-   $\alpha_0 = \pi$      (left)
-   $\alpha_1 = -\pi/3$    (down-right)
-   $\alpha_2 = +\pi/3$    (up-right)

Then apply rotation $\theta$:
-   $\beta_k = \alpha_k + \theta$ for $k = 0,1,2$

Tri2 vertices (CCW) are:
-   $V_k = (p,q) + r \cdot (\cos \beta_k, \sin \beta_k)$, for $k = 0,1,2$

So:
-   $\theta = 0$ gives Tri2 aligned with Tri1.
-   Increasing $\theta$ rotates Tri2 around its centroid.
-   The vertex order stays CCW by construction.

Because equilateral triangles repeat every $120^\circ$, you can restrict:
$$ \theta \in [0, 2\pi/3] $$

---

### 4) From "all vertices satisfy inset inequalities" to constraints on the centroid

A direct containment check would enforce the inset inequalities for **each** of the 3 vertices $\to$ 9 inequalities.

Instead, we do this side-by-side:
> For each side of Tri1, find which vertex of Tri2 is "worst" (closest to violating that side).
> If the centroid satisfies the inequality for that worst vertex, then all vertices satisfy it.

This works because:
-   Each side inequality is linear in $(x,y)$.
-   Tri2 is convex.
-   The maximum (or minimum) of a linear function over a polygon occurs at a **vertex**.

#### Define "worst-vertex" terms $M_1, M_2, M_3$ (depend only on $\theta$)

Let $\beta_k = \alpha_k + \theta$.

For the upper inset inequality $y \le \frac{1}{\sqrt{3}}x - \frac{2}{\sqrt{3}}d$, rearrange to:
$$ y - \frac{1}{\sqrt{3}}x \le -\frac{2}{\sqrt{3}}d $$

When substituting $V_k = (p,q) + r(\cos \beta_k, \sin \beta_k)$, the vertex-dependent part becomes:
$$ \sin \beta_k - \frac{1}{\sqrt{3}}\cos \beta_k $$

So define:
$$ M_1(\theta) = \max_{k} \left[ \sin(\beta_k) - \frac{1}{\sqrt{3}}\cos(\beta_k) \right] $$

For the lower inset inequality $y \ge -\frac{1}{\sqrt{3}}x + \frac{2}{\sqrt{3}}d$, rewrite as:
$$ -y - \frac{1}{\sqrt{3}}x \le -\frac{2}{\sqrt{3}}d $$

Vertex-dependent part becomes:
$$ -\sin \beta_k - \frac{1}{\sqrt{3}}\cos \beta_k $$

Define:
$$ M_2(\theta) = \max_{k} \left[ -\sin(\beta_k) - \frac{1}{\sqrt{3}}\cos(\beta_k) \right] $$

For the right inset inequality $x \le \frac{\sqrt{3}}{2}L - d$, the vertex-dependent part is:
$$ \cos \beta_k $$

Define:
$$ M_3(\theta) = \max_{k} \left[ \cos(\beta_k) \right] $$

#### Centroid inequalities (guarantee all Tri2 vertices satisfy clearance)

With these definitions, it is sufficient to enforce **three** inequalities on the centroid $(p,q)$:

-   (upper)
    $$ q \le \frac{1}{\sqrt{3}}p - \frac{2}{\sqrt{3}}d - r \cdot M_1(\theta) $$

-   (lower)
    $$ q \ge -\frac{1}{\sqrt{3}}p + \frac{2}{\sqrt{3}}d + r \cdot M_2(\theta) $$

-   (right)
    $$ p \le \frac{\sqrt{3}}{2}L - d - r \cdot M_3(\theta) $$

If these hold, then *every* vertex $V_k$ lies inside the inset triangle, so Tri2 is inside Tri1 with clearance $\ge d$.

---

### 5) The feasible centroid region is a triangle (for fixed $r, \theta$)

For fixed $r$ and $\theta$, the three centroid inequalities above define an intersection of three half-planes, which is a (possibly empty) **triangle** in the $(p,q)$ plane.

Define the right-side bound:
$$ p_B = \frac{\sqrt{3}}{2}L - d - r \cdot M_3 $$

Then the three vertices of the centroid-feasible triangle are:

-   $P_U$ (right bound + upper bound):
    $$ P_U = \left( p_B , \frac{1}{\sqrt{3}}(p_B - 2d) - r M_1 \right) $$

-   $P_L$ (right bound + lower bound):
    $$ P_L = \left( p_B , -\frac{1}{\sqrt{3}}(p_B - 2d) + r M_2 \right) $$

-   $P_A$ (upper bound $\cap$ lower bound):
    $$ P_A = \left( 2d + \frac{\sqrt{3}}{2}r(M_1+M_2) , \frac{1}{2}r(M_2 - M_1) \right) $$

Any centroid inside this triangle is feasible.

---

### 6) Maximum feasible size $r_{\max}(\theta, d)$

If $r$ is too large, the centroid-feasible triangle becomes empty.

A convenient (and tight) way to ensure non-emptiness is to require $P_{A,x} \le p_B$, which yields the maximum feasible radius:

$$ r_{\max} = \frac{L - 2\sqrt{3}d}{ (M_1+M_2) + \frac{2}{\sqrt{3}}M_3 } $$

This formula also shows feasibility requires $L - 2\sqrt{3}d > 0$, equivalent to $d < \frac{\sqrt{3}}{6}L$.

Finally, choose a normalized size parameter $s \in [0,1]$ and set:
-   $r = (1 - \epsilon) \cdot s \cdot r_{\max}$
-   $\ell = \sqrt{3} \cdot r$

$\epsilon$ is a tiny safety margin (e.g., $10^{-9}$) to avoid numerical boundary issues.

---

### 7) Picking the centroid using $c_1, c_2, c_3$ in $[0,1]$

You requested the pattern:
-   choose $c_1, c_2, c_3 \in [0,1]$,
-   normalize them to weights,
-   use those weights to combine the triangle vertices.

Let $S = c_1 + c_2 + c_3$.
-   If $S > 0$: $w_1 = c_1/S,\; w_2 = c_2/S,\; w_3 = c_3/S$.
-   Else (rare corner): $w_1 = w_2 = w_3 = 1/3$.

Then pick centroid:
$$ (p,q) = w_1 P_U + w_2 P_L + w_3 P_A $$

Because $w_i \ge 0$ and $\sum w_i = 1$, $(p,q)$ lies inside the centroid-feasible triangle and is therefore feasible.

---

## Part 2 — Implementation (recipe only)

### Inputs
-   $L > 0$ (Tri1 side length, fixed)
-   $d \ge 0$ (minimum clearance to Tri1 boundary, fixed)
-   $\theta \in [0, 2\pi/3]$
-   $s \in [0, 1]$
-   $c_1, c_2, c_3 \in [0, 1]$
-   $\epsilon$ small, e.g. $10^{-9}$
-   $\sqrt{3} \approx 1.73205$

### Output
-   Tri1 vertices in CCW: $A, B, C$
-   Tri2 vertices in CCW: $V_0, V_1, V_2$
-   (optional) centroid $(p,q)$, side length $\ell$

### Algorithm

1.  **Tri1 (CCW)**
    -   $A = (0, 0)$
    -   $B = ( \frac{\sqrt{3}}{2}L , -L/2 )$
    -   $C = ( \frac{\sqrt{3}}{2}L , +L/2 )$

2.  **Base angles for Tri2 (CCW, aligned when $\theta=0$)**
    -   $\alpha_0 = \pi$
    -   $\alpha_1 = -\pi/3$
    -   $\alpha_2 = +\pi/3$

3.  **Rotated angles**
    -   $\beta_k = \alpha_k + \theta$ for $k=0,1,2$

4.  **Compute $M_1, M_2, M_3$**
    For each $\beta_k$ in $\{\beta_0, \beta_1, \beta_2\}$ compute:
    -   $t_{1,k} = \sin(\beta_k) - \frac{1}{\sqrt{3}}\cos(\beta_k)$
    -   $t_{2,k} = -\sin(\beta_k) - \frac{1}{\sqrt{3}}\cos(\beta_k)$
    -   $t_{3,k} = \cos(\beta_k)$

    Then:
    -   $M_1 = \max_k(t_{1,k})$
    -   $M_2 = \max_k(t_{2,k})$
    -   $M_3 = \max_k(t_{3,k})$

5.  **Compute $r_{\max}$**
    -   $N = L - 2\sqrt{3}d$
    -   $D = (M_1 + M_2) + \frac{2}{\sqrt{3}}M_3$
    -   $r_{\max} = \max(0, N / D)$   (if $N \le 0$, no positive-size triangle fits)

6.  **Choose $r$ and $\ell$**
    -   $r = (1 - \epsilon) \cdot s \cdot r_{\max}$
    -   $\ell = \sqrt{3} \cdot r$

7.  **Centroid-feasible triangle vertices**
    -   $p_B = \frac{\sqrt{3}}{2}L - d - r M_3$

    -   $P_U = ( p_B ,\; \frac{1}{\sqrt{3}}p_B - \frac{2}{\sqrt{3}}d - r M_1 )$
    -   $P_L = ( p_B ,\; -\frac{1}{\sqrt{3}}p_B + \frac{2}{\sqrt{3}}d + r M_2 )$
    -   $P_A = ( 2d + \frac{\sqrt{3}}{2}r(M_1+M_2) ,\; \frac{1}{2}r(M_2 - M_1) )$

8.  **Normalize $c_1, c_2, c_3$ to weights**
    -   $S = c_1 + c_2 + c_3$
    -   if $S > 0$:
        -   $w_1 = c_1/S, \quad w_2 = c_2/S, \quad w_3 = c_3/S$
    -   else:
        -   $w_1 = w_2 = w_3 = 1/3$

9.  **Centroid**
    -   $(p,q) = w_1 P_U + w_2 P_L + w_3 P_A$

10. **Tri2 vertices (CCW)**
    -   $V_0 = (p,q) + r(\cos \beta_0, \sin \beta_0)$
    -   $V_1 = (p,q) + r(\cos \beta_1, \sin \beta_1)$
    -   $V_2 = (p,q) + r(\cos \beta_2, \sin \beta_2)$

### Guarantee
If $d < \frac{\sqrt{3}}{6}L$, then for any inputs in the given bounds, the produced Tri2 lies inside Tri1 and the minimum distance from Tri2 to Tri1’s boundary is at least $d$ (up to the numerical margin controlled by $\epsilon$).