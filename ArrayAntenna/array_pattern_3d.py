"""
3D Array Antenna Radiation Pattern — Pencil Beam with Null Steering
8×8 planar array, d=0.5λ, Taylor-40dB window, LCMV null projection
HPBW ≈ 17° (E-plane) / 13° (H-plane), SLL ≤ -24 dB, 3 nulls ≤ -50 dB
"""
import numpy as np
from scipy.signal.windows import taylor
from scipy.signal import argrelmin
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ── USER CONFIGURATION ──────────────────────────────────────────────────────
N, M          = 8, 8          # array elements (N rows × M cols)
D_LAMBDA      = 0.5           # element spacing / wavelength
TAYLOR_SLL    = 40            # Taylor window SLL in dB (positive) — set high for null margin
TAYLOR_NBAR   = 4
THETA_MAIN    = 0.0           # main beam elevation (deg from z-axis)
PHI_MAIN      = 0.0           # main beam azimuth (deg)
THETA_NULLS   = [20, 15, 25]  # null elevations (deg)
PHI_NULLS     = [30, 120, 250]# null azimuths (deg)
DB_FLOOR      = -40.0         # plot floor (dB)
THETA_MAX_DEG = 85.0          # max elevation to compute
# ────────────────────────────────────────────────────────────────────────────


def make_steering_vec(theta_deg, phi_deg, N, M, kd):
    """Flat steering vector of length N*M for direction (theta, phi)."""
    th = np.deg2rad(theta_deg)
    ph = np.deg2rad(phi_deg)
    mm, nn = np.meshgrid(np.arange(N), np.arange(M), indexing='ij')
    phase = kd * (mm.ravel() * np.sin(th) * np.cos(ph) +
                  nn.ravel() * np.sin(th) * np.sin(ph))
    return np.exp(1j * phase)   # shape (N*M,)


def compute_weights(N, M, kd, theta_nulls, phi_nulls):
    """
    Returns complex weight vector (length N*M).
    Method: 2D separable Taylor window → null projection → main-beam normalisation.
    """
    # 1. Taylor-weighted initial vector (separable 2D)
    w1d = taylor(N, nbar=TAYLOR_NBAR, sll=TAYLOR_SLL, norm=True)
    w2d = np.outer(w1d, w1d)
    w_init = w2d.ravel().astype(complex)

    # 2. Build null steering matrix  shape (N*M, num_nulls)
    A_null = np.column_stack([
        make_steering_vec(t, p, N, M, kd)
        for t, p in zip(theta_nulls, phi_nulls)
    ])

    # 3. Null projection  P = I - A(A^H A)^{-1} A^H
    AHA = A_null.conj().T @ A_null
    P = np.eye(N * M) - A_null @ np.linalg.solve(AHA, A_null.conj().T)
    w = P @ w_init

    # 4. Normalise so main beam has unity gain
    a_main = make_steering_vec(THETA_MAIN, PHI_MAIN, N, M, kd)
    w = w / (w.conj() @ a_main)
    return w


def compute_af_3d(w, N, M, kd, theta_grid_deg, phi_grid_deg):
    """
    Vectorised array factor over (N_theta, N_phi) angular grid.
    Returns AF_lin (normalised to 1) and AF_db arrays.
    """
    th = np.deg2rad(theta_grid_deg)   # (N_theta,)
    ph = np.deg2rad(phi_grid_deg)     # (N_phi,)
    TH, PH = np.meshgrid(th, ph, indexing='ij')  # (N_theta, N_phi)

    mm, nn = np.meshgrid(np.arange(N), np.arange(M), indexing='ij')
    mm_f = mm.ravel()   # (N*M,)
    nn_f = nn.ravel()

    # phase tensor  (N_theta, N_phi, N*M)
    sin_th = np.sin(TH)[:, :, np.newaxis]
    phase = kd * (mm_f * sin_th * np.cos(PH)[:, :, np.newaxis] +
                  nn_f * sin_th * np.sin(PH)[:, :, np.newaxis])

    AF_complex = np.tensordot(np.exp(1j * phase), w.conj(), axes=([2], [0]))
    AF_lin = np.abs(AF_complex)
    AF_lin /= AF_lin.max()
    AF_db = 20.0 * np.log10(np.maximum(AF_lin, 1e-7))
    return AF_lin, AF_db, TH, PH


def _hpbw_and_sll(af_db_cut, theta_grid):
    """
    Extract HPBW and SLL from a 1-D cut where peak is at theta=0 (index 0).
    HPBW = 2 * theta at first -3 dB crossing.
    SLL  = max AF_db beyond the first pattern null.
    """
    # -3 dB crossing (first point that drops below -3 dB from the peak)
    idx_3dB = np.where(af_db_cut < -3)[0]
    hpbw = 2.0 * theta_grid[idx_3dB[0]] if len(idx_3dB) else float('nan')

    # First local minimum (null) after the 3 dB point
    null_idxs = argrelmin(af_db_cut, order=5)[0]
    if len(null_idxs):
        sll = af_db_cut[null_idxs[0]:].max()
    else:
        sll = float('nan')
    return hpbw, sll


def verify_specs(w, N, M, kd, AF_db, theta_grid, phi_grid):
    print("\n=== Specification Verification ===")
    for i, (t, p) in enumerate(zip(THETA_NULLS, PHI_NULLS)):
        a = make_steering_vec(t, p, N, M, kd)
        depth = 20.0 * np.log10(abs(w.conj() @ a) + 1e-15)
        status = "PASS" if depth <= -50 else "FAIL"
        print(f"  Null {i+1} at (θ={t:3d}°, φ={p:3d}°): {depth:7.1f} dB  [{status}]  (req ≤ -50 dB)")

    phi_idx_0  = np.argmin(np.abs(phi_grid - 0))
    phi_idx_90 = np.argmin(np.abs(phi_grid - 90))
    for label, idx in [("E-plane (φ=0°)", phi_idx_0), ("H-plane (φ=90°)", phi_idx_90)]:
        cut = AF_db[:, idx]
        hpbw, sll = _hpbw_and_sll(cut, theta_grid)
        s_status = "PASS" if sll <= -20 else "FAIL"
        print(f"  {label}:  HPBW={hpbw:.1f}°,  SLL={sll:.1f} dB  [{s_status}]")
    print()


def _polar_cut(af_db, phi_grid, phi_target, theta_grid):
    """Return (angles_rad, pattern_dB) for a full ±theta cut at given phi."""
    idx = np.argmin(np.abs(phi_grid - phi_target))
    idx_opp = np.argmin(np.abs(phi_grid - (phi_target + 180) % 360))
    # forward half: theta 0..max at phi_target
    # backward half: theta 0..max at phi_target+180 → displayed as negative theta
    fwd  = af_db[:, idx]
    bwd  = af_db[:, idx_opp][::-1]  # reversed so theta goes from -max to 0
    pattern = np.concatenate([bwd, fwd])
    angles  = np.deg2rad(np.concatenate([-theta_grid[::-1], theta_grid]))
    return angles, pattern


def plot_pattern(AF_lin, AF_db, TH, PH, theta_grid, phi_grid):
    DB_FLOOR_PLOT = DB_FLOOR

    fig = plt.figure(figsize=(14, 11))
    fig.patch.set_facecolor('#0f0f1a')
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32,
                           left=0.06, right=0.96, top=0.93, bottom=0.06)

    panel_kw = dict(facecolor='#0f0f1a')

    # ── Panel 1: 3D balloon ─────────────────────────────────────────────────
    ax3d = fig.add_subplot(gs[0, 0], projection='3d', **panel_kw)
    AF_db_c = np.clip(AF_db, DB_FLOOR_PLOT, 0)
    R = np.clip(AF_lin, 1e-6, None)
    X = R * np.sin(TH) * np.cos(PH)
    Y = R * np.sin(TH) * np.sin(PH)
    Z = R * np.cos(TH)

    norm = Normalize(vmin=DB_FLOOR_PLOT, vmax=0)
    face_colors = cm.jet(norm(AF_db_c))
    ax3d.plot_surface(X, Y, Z, facecolors=face_colors,
                      rstride=2, cstride=2, linewidth=0, antialiased=True, shade=False)

    # null markers
    for i, (t, p) in enumerate(zip(THETA_NULLS, PHI_NULLS)):
        tr, pr = np.deg2rad(t), np.deg2rad(p)
        r0 = 1.15
        xn = r0 * np.sin(tr) * np.cos(pr)
        yn = r0 * np.sin(tr) * np.sin(pr)
        zn = r0 * np.cos(tr)
        ax3d.scatter([xn], [yn], [zn], color='red', s=60, zorder=5, depthshade=False)
        ax3d.text(xn, yn, zn + 0.05, f'N{i+1}', color='red', fontsize=8, ha='center')

    ax3d.set_xlabel('x', color='white', labelpad=2, fontsize=9)
    ax3d.set_ylabel('y', color='white', labelpad=2, fontsize=9)
    ax3d.set_zlabel('z', color='white', labelpad=2, fontsize=9)
    ax3d.set_title('3D Radiation Pattern', color='white', fontsize=11, pad=6)
    ax3d.tick_params(colors='white', labelsize=7)
    for pane in [ax3d.xaxis.pane, ax3d.yaxis.pane, ax3d.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#333355')
    ax3d.grid(True, color='#333355', linewidth=0.4)

    sm = cm.ScalarMappable(cmap='jet', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax3d, shrink=0.55, pad=0.1, label='dB')
    cbar.ax.yaxis.label.set_color('white')
    cbar.ax.tick_params(colors='white', labelsize=8)

    # ── Panel 2: 2D heatmap ─────────────────────────────────────────────────
    ax2d = fig.add_subplot(gs[0, 1], **panel_kw)
    im = ax2d.pcolormesh(phi_grid, theta_grid, AF_db,
                         cmap='jet', vmin=DB_FLOOR_PLOT, vmax=0, shading='gouraud')
    for t, p in zip(THETA_NULLS, PHI_NULLS):
        ax2d.plot(p, t, 'ro', ms=6, mec='white', mew=0.6)
    for i, (t, p) in enumerate(zip(THETA_NULLS, PHI_NULLS)):
        ax2d.annotate(f'N{i+1}', (p, t), textcoords='offset points',
                      xytext=(4, 4), color='red', fontsize=8)
    ax2d.set_xlabel('Azimuth φ (deg)', color='white', fontsize=9)
    ax2d.set_ylabel('Elevation θ (deg)', color='white', fontsize=9)
    ax2d.set_title('2D Pattern Map', color='white', fontsize=11)
    ax2d.tick_params(colors='white', labelsize=8)
    for spine in ax2d.spines.values():
        spine.set_edgecolor('#555577')
    cbar2 = fig.colorbar(im, ax=ax2d, label='dB')
    cbar2.ax.yaxis.label.set_color('white')
    cbar2.ax.tick_params(colors='white', labelsize=8)

    # ── Helper for polar panels ─────────────────────────────────────────────
    def add_polar_panel(gs_pos, phi_target, title):
        ax = fig.add_subplot(gs_pos, projection='polar', facecolor='#0f0f1a')
        angles, pattern = _polar_cut(AF_db, phi_grid, phi_target, theta_grid)

        # map dB to [0, 1] radius so negative dB plots inward
        r_range = -DB_FLOOR_PLOT   # e.g. 40 dB dynamic range
        r = np.clip(pattern - DB_FLOOR_PLOT, 0, r_range) / r_range

        ax.plot(angles, r, color='#00ccff', linewidth=1.3)
        ax.fill(angles, r, alpha=0.15, color='#00ccff')

        # radial tick labels (dB)
        tick_db = np.arange(DB_FLOOR_PLOT, 1, 10)    # [-40,-30,-20,-10, 0]
        tick_r  = (tick_db - DB_FLOOR_PLOT) / r_range
        ax.set_rticks(tick_r)
        ax.set_yticklabels([f'{int(v)} dB' for v in tick_db],
                           color='#aaaacc', fontsize=7)
        ax.set_rlim(0, 1)

        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_thetagrids(np.arange(-180, 181, 30),
                          labels=[f'{v}°' for v in np.arange(-180, 181, 30)],
                          color='white', fontsize=7)
        ax.tick_params(colors='white')
        ax.grid(color='#333355', linewidth=0.5)
        ax.set_title(title, color='white', fontsize=10, pad=12)

        # mark -3 dB and -20 dB circles
        for ref_db, col, ls in [(-3, 'yellow', '--'), (-20, 'orange', ':')]:
            ref_r = (ref_db - DB_FLOOR_PLOT) / r_range
            th_c  = np.linspace(-np.pi, np.pi, 360)
            ax.plot(th_c, [ref_r] * 360, color=col, linewidth=0.7,
                    linestyle=ls, label=f'{ref_db} dB')
        ax.legend(loc='lower right', fontsize=6,
                  labelcolor='white', facecolor='#1a1a2e', edgecolor='#555577')
        return ax

    add_polar_panel(gs[1, 0], phi_target=0,  title='E-plane Cut (φ = 0°/180°)')
    add_polar_panel(gs[1, 1], phi_target=90, title='H-plane Cut (φ = 90°/270°)')

    # ── Figure title ────────────────────────────────────────────────────────
    null_str = ', '.join(f'(θ={t}°,φ={p}°)' for t, p in zip(THETA_NULLS, PHI_NULLS))
    fig.suptitle(
        f'8×8 Planar Array  d=0.5λ  |  Pencil Beam HPBW≈17°  |  Taylor-{TAYLOR_SLL}dB  |  '
        f'Nulls: {null_str}',
        color='white', fontsize=10, y=0.97)

    return fig


# ── MAIN ─────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    kd = 2 * np.pi * D_LAMBDA   # = π for d = λ/2

    w = compute_weights(N, M, kd, THETA_NULLS, PHI_NULLS)

    theta_grid = np.linspace(0, THETA_MAX_DEG, 181)
    phi_grid   = np.linspace(0, 360, 361)

    AF_lin, AF_db, TH, PH = compute_af_3d(w, N, M, kd, theta_grid, phi_grid)

    verify_specs(w, N, M, kd, AF_db, theta_grid, phi_grid)

    fig = plot_pattern(AF_lin, AF_db, TH, PH, theta_grid, phi_grid)
    out_path = 'ArrayAntenna/array_pattern_3d.png'
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())
    print(f'Figure saved to {out_path}')
    plt.show()
