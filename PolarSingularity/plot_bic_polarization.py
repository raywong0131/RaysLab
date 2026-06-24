import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.colors import TwoSlopeNorm
import sys
import os

# ── complex string parser ──────────────────────────────────────────────────────
_cpx_re = re.compile(
    r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)'   # real part
    r'([+-](?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i'   # imaginary part (required)
    r'|([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i' # pure imaginary
    r'|([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)'  # pure real
)

def parse_complex(s):
    s = str(s).strip()
    m = re.fullmatch(
        r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)'
        r'([+-](?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i',
        s)
    if m:
        return complex(float(m.group(1)), float(m.group(2)))
    m = re.fullmatch(r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i', s)
    if m:
        return complex(0, float(m.group(1)))
    try:
        return complex(float(s), 0)
    except ValueError:
        return complex(np.nan, np.nan)

# ── load data ──────────────────────────────────────────────────────────────────
def load_file(path):
    df = pd.read_excel(path, header=None)
    df.columns = ['param', 'kx', 'ky', 'cx_str', 'cy_str']
    df['cx'] = df['cx_str'].apply(parse_complex)
    df['cy'] = df['cy_str'].apply(parse_complex)
    return df

# ── Stokes / ellipse parameters from Jones vector ─────────────────────────────
def jones_to_ellipse(cx, cy):
    """Returns (theta_deg, a, b, s3_sign) for each (cx, cy) pair.
    theta: orientation angle of major axis (degrees)
    a, b : semi-axes (normalised so a²+b²=1)
    s3_sign: +1 left-hand, -1 right-hand, 0 linear
    """
    S0 = np.abs(cx)**2 + np.abs(cy)**2
    S1 = np.abs(cx)**2 - np.abs(cy)**2
    S2 = 2 * np.real(cx * np.conj(cy))
    S3 = 2 * np.imag(cx * np.conj(cy))

    with np.errstate(invalid='ignore', divide='ignore'):
        S_lin = np.sqrt(S1**2 + S2**2)
        a = np.sqrt(np.maximum((S0 + S_lin) / 2, 0))
        b = np.sqrt(np.maximum((S0 - S_lin) / 2, 0))
        # normalise so major axis is 1
        norm = np.where(a > 0, a, 1)
        b_norm = b / norm
        theta = 0.5 * np.degrees(np.arctan2(S2, S1))

    s3_sign = np.sign(S3)
    return theta, np.ones_like(theta), b_norm, s3_sign

# ── rotation arrow helper ─────────────────────────────────────────────────────
def _add_rotation_arrow(ax, x0, y0, width, height, ang_deg, sg, color, lw=0.9):
    """Draw a small tangent arrow at the major-axis tip to indicate rotation sense."""
    if sg == 0:
        return
    ang_rad = np.radians(ang_deg)
    a = width / 2
    b = height / 2
    # At t=π/2 the tangent is exactly parallel to the major axis.
    # As b→0 (near-linear), this point collapses to the centre (midpoint of segment).
    px = x0 - b * np.sin(ang_rad)
    py = y0 + b * np.cos(ang_rad)
    # CCW tangent at t=π/2: (-cos θ, -sin θ) = negative major-axis direction
    tx = -np.cos(ang_rad)
    ty = -np.sin(ang_rad)
    if sg < 0:          # RCP → CW: reverse tangent
        tx, ty = -tx, -ty
    half = 0.4 * a      # half arrow length
    ax.annotate('',
                xy=(px + tx * half, py + ty * half),
                xytext=(px - tx * half, py - ty * half),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                mutation_scale=5))

# ── main plotting routine ──────────────────────────────────────────────────────
def plot_bic(filepath, output_dir=None, subsample=2):
    df = load_file(filepath)
    param_val = df['param'].iloc[0]
    fname_base = os.path.splitext(os.path.basename(filepath))[0]
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(filepath), 'Plot_output')
    os.makedirs(output_dir, exist_ok=True)

    kx_vals = np.sort(df['kx'].unique())
    ky_vals = np.sort(df['ky'].unique())
    Nx, Ny = len(kx_vals), len(ky_vals)

    # pivot to 2D grids
    df_piv = df.set_index(['kx', 'ky'])
    CX = np.array([[df_piv.loc[(kx, ky), 'cx'] for ky in ky_vals]
                   for kx in kx_vals], dtype=complex)  # shape (Nx, Ny)
    CY = np.array([[df_piv.loc[(kx, ky), 'cy'] for ky in ky_vals]
                   for kx in kx_vals], dtype=complex)

    KX, KY = np.meshgrid(kx_vals, ky_vals, indexing='ij')

    cx_int = np.abs(CX)**2
    cy_int = np.abs(CY)**2
    norm_factor = cx_int.max()
    cx_int = cx_int / norm_factor
    cy_int = cy_int / norm_factor

    extent = [kx_vals.min(), kx_vals.max(), ky_vals.min(), ky_vals.max()]

    def _smooth_colormap(data, label, title, outpath, cmap='RdBu_r'):
        fig, ax = plt.subplots(figsize=(5, 4.5))
        # transpose so rows=ky, cols=kx (imshow origin='lower')
        im = ax.imshow(data.T, origin='lower', extent=extent,
                       aspect='equal', cmap=cmap,
                       interpolation='bicubic')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(label, fontsize=12)
        ax.set_xlabel(r'$k_x \ (2\pi/a)$', fontsize=12)
        ax.set_ylabel(r'$k_y \ (2\pi/a)$', fontsize=12)
        ax.set_title(title, fontsize=12)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        print(f'Saved: {outpath}')

    out1 = os.path.join(output_dir, fname_base + '_cx.png')
    _smooth_colormap(cx_int, r'$|c_x|^2$',
                     rf'$|c_x|^2$ — param = {param_val}', out1)

    out2 = os.path.join(output_dir, fname_base + '_cy.png')
    _smooth_colormap(cy_int, r'$|c_y|^2$',
                     rf'$|c_y|^2$ — param = {param_val}', out2)

    # ── Figure 3: polarisation ellipses ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5.5, 5))

    # subsample grid
    idx_x = np.arange(0, Nx, subsample)
    idx_y = np.arange(0, Ny, subsample)

    kx_step = kx_vals[1] - kx_vals[0]
    ky_step = ky_vals[1] - ky_vals[0]
    # semi-major axis = 45% of subsampled grid spacing → ~10% gap between neighbours
    ellipse_scale = 0.45 * subsample * min(kx_step, ky_step)

    cx_sub = CX[np.ix_(idx_x, idx_y)]
    cy_sub = CY[np.ix_(idx_x, idx_y)]
    theta, a_ax, b_ax, s3_sign = jones_to_ellipse(cx_sub, cy_sub)

    for ii, ix in enumerate(idx_x):
        for jj, jy in enumerate(idx_y):
            x0 = kx_vals[ix]
            y0 = ky_vals[jy]
            ang = theta[ii, jj]
            b_n = b_ax[ii, jj]  # normalised minor axis (0 to 1)
            sg  = s3_sign[ii, jj]

            width  = 2 * ellipse_scale          # major axis length
            height = 2 * ellipse_scale * b_n    # minor axis length

            color = '#d62728' if sg > 0 else ('#1f77b4' if sg < 0 else 'k')
            e = Ellipse(xy=(x0, y0),
                        width=width, height=height,
                        angle=ang,
                        edgecolor=color, facecolor='none', linewidth=0.9)
            ax.add_patch(e)
            _add_rotation_arrow(ax, x0, y0, width, height, ang, sg, color)

    ax.set_xlim(kx_vals.min() - kx_step, kx_vals.max() + kx_step)
    ax.set_ylim(ky_vals.min() - ky_step, ky_vals.max() + ky_step)
    ax.set_xlabel(r'$k_x \ (2\pi/a)$', fontsize=12)
    ax.set_ylabel(r'$k_y \ (2\pi/a)$', fontsize=12)
    ax.set_title(rf'Far-field polarisation — param = {param_val}', fontsize=12)
    ax.set_aspect('equal')
    fig.tight_layout()
    out3 = os.path.join(output_dir, fname_base + '_ellipses.png')
    fig.savefig(out3, dpi=150)
    plt.close(fig)
    print(f'Saved: {out3}')

    # ── Figure 4: polarisation ellipses with intensity-scaled size (log) ──────
    fig, ax = plt.subplots(figsize=(5.5, 5))

    # intensity at subsampled points (use unnormalised raw values for log scaling)
    I_sub = np.abs(cx_sub)**2 + np.abs(cy_sub)**2
    log_I = np.log10(np.maximum(I_sub, 1e-30))
    log_min, log_max = log_I.min(), log_I.max()
    if log_max > log_min:
        size_norm = (log_I - log_min) / (log_max - log_min)  # 0 … 1
    else:
        size_norm = np.ones_like(log_I)
    # map to [0.15, 1.0] so even the weakest points remain visible
    size_factor = 0.15 + 0.85 * size_norm

    for ii, ix in enumerate(idx_x):
        for jj, jy in enumerate(idx_y):
            x0 = kx_vals[ix]
            y0 = ky_vals[jy]
            ang = theta[ii, jj]
            b_n = b_ax[ii, jj]
            sg  = s3_sign[ii, jj]
            sf  = size_factor[ii, jj]

            width  = 2 * ellipse_scale * sf
            height = 2 * ellipse_scale * sf * b_n

            color = '#d62728' if sg > 0 else ('#1f77b4' if sg < 0 else 'k')
            e = Ellipse(xy=(x0, y0),
                        width=width, height=height,
                        angle=ang,
                        edgecolor=color, facecolor='none', linewidth=0.9)
            ax.add_patch(e)
            _add_rotation_arrow(ax, x0, y0, width, height, ang, sg, color)

    ax.set_xlim(kx_vals.min() - kx_step, kx_vals.max() + kx_step)
    ax.set_ylim(ky_vals.min() - ky_step, ky_vals.max() + ky_step)
    ax.set_xlabel(r'$k_x \ (2\pi/a)$', fontsize=12)
    ax.set_ylabel(r'$k_y \ (2\pi/a)$', fontsize=12)
    ax.set_title(rf'Far-field polarisation (log-$I$ size) — param = {param_val}', fontsize=12)
    ax.set_aspect('equal')
    fig.tight_layout()
    out4 = os.path.join(output_dir, fname_base + '_ellipses_logI.png')
    fig.savefig(out4, dpi=150)
    plt.close(fig)
    print(f'Saved: {out4}')

# ── entry point ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    folder = r'D:\Claude Code\PolarSingularity'
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        # process first file by default for preview
        import glob
        files = sorted(glob.glob(os.path.join(folder, '*.xlsx')))[:]

    for f in files:
        print(f'Processing: {f}')
        plot_bic(f, subsample=2)
