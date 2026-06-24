import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Patch

# ── Grid ─────────────────────────────────────────────────────────────────────
N = 100
x = np.linspace(-2, 2, N)
y = np.linspace(-2, 2, N)
X, Y = np.meshgrid(x, y)
D = 2.0

# ── Phase-transition model ────────────────────────────────────────────────────
#
#   δ  = z1 − z2              (control parameter)
#   αx = tanh(s · δ)          (smooth phase switch, s = sharpness)
#
#   Surface 1:  Z = k·αx·x²  + k·y²  + z1
#   Surface 2:  Z = −k·αx·x² − k·y²  + z2
#
#   δ < 0  →  αx < 0  →  two saddle surfaces  (x: dip/rise, y: rise/dip)
#   δ = 0  →  αx = 0  →  two parabolic cylinders along x (transition)
#   δ > 0  →  αx > 0  →  paraboloid (up) + paraboloid (down), NO intersection
#
#   Why no intersection when δ > 0:
#     Setting Z1 = Z2:  αx·x² + y² = −δ/(2k) < 0  →  no real solution
#
#   Intersection (only when δ ≤ 0):
#     αx·x² + y² = −δ/(2k)   at   z = (z1 + z2)/2
#     δ = 0  →  y = 0  (line along x-axis)
#     δ < 0  →  y² − |αx|·x² = C > 0  →  hyperbola opening along y

def compute_surfaces(k, z1, z2, s):
    delta = z1 - z2
    ax = np.tanh(s * delta)           # phase parameter
    Z1 = k * ax * X**2 + k * Y**2 + z1
    Z2 = -k * ax * X**2 - k * Y**2 + z2
    return Z1, Z2, ax, delta


def compute_intersection(k, z1, z2, delta, ax):
    """
    Analytically solve  αx·x² + y² = C  where C = −δ/(2k).
    Returns curve segments and a label string.
    """
    z_mid = (z1 + z2) / 2.0
    C = -delta / (2.0 * k)     # RHS
    curves = []

    # ── Paraboloid regime: C < 0 → no real solution ───────────────────────
    if delta > 1e-6:
        return [], 'None  (paraboloids do not intersect)'

    # ── Transition: δ = 0, C = 0 → y = 0 (line along x-axis) ────────────
    if abs(delta) <= 1e-6:
        t = np.linspace(-D, D, 400)
        curves.append((t, np.zeros_like(t), np.full_like(t, z_mid)))
        return curves, 'Line  y = 0  (surfaces touch at z = {:.3f})'.format(z_mid)

    # ── Saddle regime: δ < 0, C > 0, αx < 0 ──────────────────────────────
    # y² − |αx|·x² = C  →  hyperbola opening along y
    # y = ±√(C + |αx|·x²)
    abs_ax = abs(ax)

    # Clip to domain: |y| ≤ D → C + |αx|·x² ≤ D²
    if C > D**2:
        # Even at x=0, y = √C > D → entire hyperbola outside domain
        return [], 'Hyperbola (outside domain — move z2 closer to z1)'

    xlim = min(D, np.sqrt((D**2 - C) / abs_ax) if abs_ax > 1e-12 else D)
    t = np.linspace(-xlim, xlim, 600)
    yp = np.sqrt(np.maximum(C + abs_ax * t**2, 0.0))
    curves.append((t,  yp, np.full_like(t, z_mid)))
    curves.append((t, -yp, np.full_like(t, z_mid)))

    y_vertex = np.sqrt(C)
    label = 'Hyperbola (along y,  vertex y\u2080 = \u00b1{:.3f})'.format(y_vertex)
    return curves, label


# ── Visual style ──────────────────────────────────────────────────────────────
BG, PANE, TICK, LC = '#1a1a2e', '#333355', '#cccccc', '#aaaacc'
INTER_COLOR, INTER_GLOW = '#ffff00', '#ff8800'

REGIME_COLORS = {
    'saddle':      '#aaaaff',
    'transition':  '#ffffff',
    'paraboloid':  '#aaffaa',
}

def style_ax(ax):
    ax.set_facecolor(BG)
    ax.tick_params(colors=TICK, labelsize=8)
    for p in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        p.fill = False
        p.set_edgecolor(PANE)
    ax.set_xlabel('X', color=LC, labelpad=5)
    ax.set_ylabel('Y', color=LC, labelpad=5)
    ax.set_zlabel('Z', color=LC, labelpad=5)

legend_elements = [
    Patch(facecolor='#66ccff', alpha=0.85, label='Surface 1'),
    Patch(facecolor='#ff8844', alpha=0.85, label='Surface 2'),
    Patch(facecolor=INTER_COLOR, alpha=1.0,  label='Intersection'),
]

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 10))
fig.patch.set_facecolor(BG)

ax3d = fig.add_axes([0.0, 0.32, 1.0, 0.66], projection='3d')

info_ax = fig.add_axes([0.08, 0.285, 0.84, 0.030])
info_ax.axis('off')
info_txt = info_ax.text(0.5, 0.5, '', transform=info_ax.transAxes,
                         ha='center', va='center', fontsize=9,
                         family='monospace')

SBG = '#2a2a4a'
ax_k  = fig.add_axes([0.15, 0.230, 0.70, 0.024], facecolor=SBG)
ax_s  = fig.add_axes([0.15, 0.178, 0.70, 0.024], facecolor=SBG)
ax_z1 = fig.add_axes([0.15, 0.114, 0.70, 0.024], facecolor=SBG)
ax_z2 = fig.add_axes([0.15, 0.052, 0.70, 0.024], facecolor=SBG)

s_k  = Slider(ax_k,  'Curvature k',   0.1,  3.0, valinit=1.0,  color='#7788ff')
s_s  = Slider(ax_s,  'Sharpness',     0.2,  6.0, valinit=2.0,  color='#88aaff')
# Default: z1 < z2 → saddle regime with visible intersection
s_z1 = Slider(ax_z1, 'Surface 1  Z', -3.0,  3.0, valinit=-0.5, color='#66ccff')
s_z2 = Slider(ax_z2, 'Surface 2  Z', -3.0,  3.0, valinit=0.5,  color='#ff8844')

for sl in (s_k, s_s, s_z1, s_z2):
    sl.label.set_color('white')
    sl.valtext.set_color('white')


# ── Main draw function ────────────────────────────────────────────────────────
def draw(k, z1, z2, sharpness):
    ax3d.clear()
    style_ax(ax3d)

    Z1, Z2, ax_val, delta = compute_surfaces(k, z1, z2, sharpness)

    ax3d.plot_surface(X, Y, Z1, cmap='cool',   alpha=0.68, linewidth=0, antialiased=True)
    ax3d.plot_surface(X, Y, Z2, cmap='autumn', alpha=0.68, linewidth=0, antialiased=True)

    curves, curve_label = compute_intersection(k, z1, z2, delta, ax_val)
    for xi, yi, zi in curves:
        ax3d.plot(xi, yi, zi, color=INTER_GLOW,  linewidth=6,   alpha=0.35)
        ax3d.plot(xi, yi, zi, color=INTER_COLOR, linewidth=2.2, alpha=1.0)

    if abs(delta) <= 1e-6:
        regime_key, regime_str = 'transition', '── TRANSITION  (z1 = z2,  surfaces touch) ──'
    elif delta > 0:
        regime_key, regime_str = 'paraboloid', 'PARABOLOID REGIME  (z1 > z2,  no intersection)'
    else:
        regime_key, regime_str = 'saddle', 'SADDLE REGIME  (z1 < z2,  intersection exists)'

    info_txt.set_color(REGIME_COLORS[regime_key])
    info_txt.set_text(
        '{:s}   |   \u03b1x = {:+.3f}   |   {:s}   |   z\u209a = {:.3f}'.format(
            regime_str, ax_val, curve_label, (z1 + z2) / 2)
    )

    ax3d.set_title('Phase Transition: Saddle Surfaces \u2194 Paraboloids',
                   color='white', fontsize=12, pad=8)
    ax3d.legend(handles=legend_elements, loc='upper left',
                facecolor='#2a2a4a', edgecolor='#555577',
                labelcolor='white', fontsize=9)
    fig.canvas.draw_idle()


# ── Init & callbacks ──────────────────────────────────────────────────────────
draw(s_k.valinit, s_z1.valinit, s_z2.valinit, s_s.valinit)

def on_change(_):
    draw(s_k.val, s_z1.val, s_z2.val, s_s.val)

for sl in (s_k, s_s, s_z1, s_z2):
    sl.on_changed(on_change)

plt.show()
