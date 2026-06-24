"""
2D Electromagnetic Field Far-Field Distribution Calculator via FFT

Input format (COMSOL exported text):
  Header lines start with '%'.  The parser reads:
    % Length unit:   nm
    % X  Y  [Z]  Ex  Ey  Ez      <- last header line = column names
  Data columns: x, y, [z], Ex, Ey, Ez
  Ex / Ey / Ez are complex strings, e.g. '1.2-3.4i'.

Supports rectangular (N×M) grids.  Coordinates are converted to µm
internally; plots display in the original data unit.
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time, os, glob, re

# ╔══════════════════════════════════════════════════════════════╗
# ║                    USER CONFIGURATION                        ║
# ╚══════════════════════════════════════════════════════════════╝

AUTO_PARAM = True    # True → read FREQ_THz and H0_UM from each data file
                     # False → use the values below
FREQ_THz = 193.97   # optical frequency [THz]      (used only when AUTO_PARAM=False)
H0_UM    = 1.55     # propagation distance [µm]     (used only when AUTO_PARAM=False)
A_UM     = 0.82     # photonic crystal lattice constant [µm]
NA       = 0.9      # numerical aperture

EXPAND   = 1000# zero-padding width (pixels added on each side, per axis)
# File selection
#   None  → all *.txt files in the same directory as this script
#   list  → specific filenames relative to the script directory
SELECT_FILES = None  # e.g. ['file_a.txt', 'file_b.txt']

# Output directory for PNG files (relative to script directory)
OUTPUT = 'Results'

# ── Plot toggles ──────────────────────────────────────────────────────────────
# Fig A: 2×2 overview  — r-space, k-space log, k-space |Ex|², iFFT near-field
PLOT_OVERVIEW     = True
# Fig B: polarization polar plot
PLOT_POLARIZATION = True
# Fig C: 3×2 cutlines  — k-cuts, θ-cuts, radiation patterns
PLOT_CUTLINES     = True
# Fig D: Gaussian fit on ky-cutline
PLOT_GAUSS        = True

# ╔══════════════════════════════════════════════════════════════╗
# ║                    UNIT HELPERS                              ║
# ╚══════════════════════════════════════════════════════════════╝

_TO_UM = {
    'pm': 1e-6, 'nm': 1e-3, 'um': 1.0, 'µm': 1.0,
    'mm': 1e3,  'cm': 1e4,  'm':  1e6,
}

def to_um(unit: str) -> float:
    """Return scale factor: 1 [unit] = ? µm."""
    return _TO_UM.get(unit.strip().lower(), 1.0)

# ╔══════════════════════════════════════════════════════════════╗
# ║                    FILE LOADER                               ║
# ╚══════════════════════════════════════════════════════════════╝

def _parse_complex(s: str) -> complex:
    """Parse 'a+bi', 'bi', or 'a' into Python complex."""
    s = s.strip()
    # a±bi
    m = re.fullmatch(
        r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)'
        r'([+-](?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i', s)
    if m:
        return complex(float(m.group(1)), float(m.group(2)))
    # pure imaginary
    m = re.fullmatch(r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)i', s)
    if m:
        return complex(0.0, float(m.group(1)))
    # pure real
    try:
        return complex(float(s), 0.0)
    except ValueError:
        return complex(np.nan, np.nan)


def _read_lines(path: str) -> list[str]:
    """
    Read all lines robustly.  COMSOL exports are ASCII text; UTF-16 variants
    just insert a \\x00 after every byte.  Strip all null bytes first, then
    decode as UTF-8 (which is identical to ASCII for COMSOL content).
    """
    with open(path, 'rb') as f:
        raw = f.read()
    # Remove UTF-16 BOM bytes if present, then strip all null bytes.
    raw = raw.lstrip(b'\xff\xfe\xfe\xff')
    clean = raw.replace(b'\x00', b'')
    return clean.decode('utf-8', errors='replace').splitlines()


def load_comsol(path: str):
    """
    Load a COMSOL-exported text file.

    Returns
    -------
    data : dict  {x, y, z (or None), Ex, Ey, Ez}  – arrays in original units
    meta : dict  {length_unit, col_names, has_z}
    """
    meta = {'length_unit': 'um', 'col_names': None, 'has_z': False,
            'freq_thz': None, 'nodes': None}
    data_lines = []

    for raw in _read_lines(path):
        line = raw.rstrip('\n')
        if line.startswith('%'):
            content = line[1:].strip()
            low = content.lower()
            if low.startswith('length unit'):
                meta['length_unit'] = content.split(':', 1)[1].strip()
            elif low.startswith('nodes'):
                try:
                    meta['nodes'] = int(content.split(':', 1)[1].strip())
                except ValueError:
                    pass
            elif content and content[0].upper() in 'XYZ':
                # Last header line describes columns, e.g. "X  Y  Z  Ex  Ey  Ez"
                meta['col_names'] = content.split()
                # Extract frequency: look for token before 'THz'
                # e.g. "... @ 198.74+0.47958i THz ..."  → real part = 198.74
                tokens = meta['col_names']
                for i, tok in enumerate(tokens):
                    if tok.lower() == 'thz' and i > 0:
                        m = re.match(r'([+-]?[\d.]+(?:[Ee][+-]?\d+)?)', tokens[i-1])
                        if m:
                            meta['freq_thz'] = float(m.group(1))
                            break
        elif line.strip():
            # Skip bare-integer lines (stray metadata without '%' prefix,
            # e.g. a lone node-count such as "-27575").
            if not re.fullmatch(r'[+-]?\d+', line.strip()):
                data_lines.append(line)

    # ── Detect coordinate columns from col_names header ──────────────
    # col_names may contain long COMSOL descriptions; only look at leading
    # tokens that are single letters X / Y / Z.
    if meta['col_names']:
        coord_tokens = [c for c in meta['col_names']
                        if c.upper() in ('X', 'Y', 'Z') and len(c) == 1]
        meta['has_z'] = 'Z' in [c.upper() for c in coord_tokens]
    # Fall back to column count if no header found
    if not meta['col_names'] and data_lines:
        ncols = len(data_lines[0].split())
        meta['has_z'] = ncols in (7, 9)   # 2+3E or 3+3E (real+imag pairs)

    nc = 3 if meta['has_z'] else 2   # number of coordinate columns

    # ── Collect valid rows (must have at least nc + 6 or nc + 3 cols) ─
    all_rows = [ln.split() for ln in data_lines]
    # Accept rows with at least nc+3 cols (complex strings) or nc+6 (real/imag pairs)
    rows = [r for r in all_rows if len(r) >= nc + 3]
    bad  = len(all_rows) - len(rows)
    if bad:
        print(f'  WARNING: skipped {bad} short lines '
              f'(first: {next(r for r in all_rows if len(r) < nc + 3)})')
    if not rows:
        print('  First 5 data lines after decode:')
        for l in data_lines[:5]:
            print('   ', repr(l))
        raise ValueError(f'No valid data rows found in {path}')

    # ── Auto-detect real/imag split vs complex-string format ─────────
    # Sample the first valid row: if col[nc] looks like a plain float
    # (no 'i') and there are nc+6 columns, assume real/imag split.
    sample = rows[0]
    split_format = (len(sample) >= nc + 6 and 'i' not in sample[nc])

    x = np.array([float(r[0]) for r in rows])
    y = np.array([float(r[1]) for r in rows])
    z = np.array([float(r[2]) for r in rows]) if meta['has_z'] else None

    if split_format:
        # Columns: nc+0=Re(Ex), nc+1=Im(Ex), nc+2=Re(Ey), nc+3=Im(Ey), ...
        Ex = np.array([float(r[nc])   + 1j * float(r[nc+1]) for r in rows])
        Ey = np.array([float(r[nc+2]) + 1j * float(r[nc+3]) for r in rows])
        Ez = np.array([float(r[nc+4]) + 1j * float(r[nc+5]) for r in rows])
        print(f'  Format: real/imag split ({len(sample)} cols/row)')
    else:
        # Columns: nc+0=Ex(complex string), nc+1=Ey, nc+2=Ez
        Ex = np.array([_parse_complex(r[nc])     for r in rows])
        Ey = np.array([_parse_complex(r[nc + 1]) for r in rows])
        Ez = np.array([_parse_complex(r[nc + 2]) for r in rows])
        print(f'  Format: complex strings ({len(sample)} cols/row)')

    return {'x': x, 'y': y, 'z': z, 'Ex': Ex, 'Ey': Ey, 'Ez': Ez}, meta

# ╔══════════════════════════════════════════════════════════════╗
# ║                    GRID MATRIXING                            ║
# ╚══════════════════════════════════════════════════════════════╝

def matrix_field(data: dict, meta: dict):
    """
    Reshape flat field arrays into 2-D (Nx × Ny) complex grids.
    All coordinates are converted to µm.
    Supports arbitrary rectangular (N×M) grids.

    Returns
    -------
    Ex2d, Ey2d, Ez2d : (Nx, Ny) complex arrays
    x_um, y_um       : sorted unique coordinate vectors [µm]
    dx, dy           : grid spacing [µm]
    Nx, Ny           : grid dimensions
    """
    s = to_um(meta['length_unit'])
    x_um = data['x'] * s
    y_um = data['y'] * s

    x_uniq = np.unique(x_um)
    y_uniq = np.unique(y_um)

    def _expand_axis(axis_uniq, expected_n):
        if expected_n is None or expected_n <= len(axis_uniq) or len(axis_uniq) < 2:
            return axis_uniq

        step = float(np.median(np.diff(axis_uniq)))
        if not np.isfinite(step) or step <= 0:
            return axis_uniq

        # COMSOL exports are ordered monotonically; if the tail is truncated,
        # rebuild the missing coordinates using the observed start and spacing.
        rebuilt = axis_uniq[0] + step * np.arange(expected_n)
        return rebuilt

    expected_nodes = meta.get('nodes')
    expected_nx = expected_ny = None
    if expected_nodes:
        if len(x_uniq) > 0 and expected_nodes % len(x_uniq) == 0:
            expected_ny = expected_nodes // len(x_uniq)
        if len(y_uniq) > 0 and expected_nodes % len(y_uniq) == 0:
            expected_nx = expected_nodes // len(y_uniq)

    x_uniq = _expand_axis(x_uniq, expected_nx)
    y_uniq = _expand_axis(y_uniq, expected_ny)
    Nx, Ny = len(x_uniq), len(y_uniq)

    eps = np.finfo(float).eps
    xi = np.searchsorted(x_uniq, x_um)
    yi = np.searchsorted(y_uniq, y_um)

    def _fill(raw):
        arr = np.full((Nx, Ny), eps, dtype=complex)
        arr[xi, yi] = np.where(np.isnan(raw), eps, raw)
        return arr

    Ex2d = _fill(data['Ex'])
    Ey2d = _fill(data['Ey'])
    Ez2d = _fill(data['Ez'])

    dx = float(x_uniq[-1] - x_uniq[0]) / (Nx - 1) if Nx > 1 else 1.0
    dy = float(y_uniq[-1] - y_uniq[0]) / (Ny - 1) if Ny > 1 else 1.0

    return Ex2d, Ey2d, Ez2d, x_uniq, y_uniq, dx, dy, Nx, Ny

# ╔══════════════════════════════════════════════════════════════╗
# ║                    FFT PIPELINE                              ║
# ╚══════════════════════════════════════════════════════════════╝

def compute_fft(Ex2d, Ey2d, Ez2d, x_um, y_um, dx, dy, Nx, Ny,
                k0: float, h0_um: float, NA: float, expand: int = 1000):
    """
    Zero-pad → 2-D FFT → phase propagation + NA filter → iFFT.

    All units in µm / µm⁻¹.
    Array layout throughout: (dimEy rows, dimEx cols)  =  (ky axis, kx axis).

    Returns a result dict containing all intermediate and final arrays.
    """
    dimEx = Nx + 2 * expand
    dimEy = Ny + 2 * expand

    # ── extended real-space axes ──────────────────────────────────────
    x_ext = np.linspace(-dx * (dimEx - 1) / 2,  dx * (dimEx - 1) / 2, dimEx)
    y_ext = np.linspace(-dy * (dimEy - 1) / 2,  dy * (dimEy - 1) / 2, dimEy)

    # ── k-space axes (Nyquist: kmax = π/d) ───────────────────────────
    kmax_x = np.pi * (Nx - 1) / (x_um[-1] - x_um[0])
    kmax_y = np.pi * (Ny - 1) / (y_um[-1] - y_um[0])
    kx = np.linspace(-kmax_x, kmax_x, dimEx)
    ky = np.linspace(-kmax_y, kmax_y, dimEy)

    # ── zero-pad and transpose to (ky, kx) layout ────────────────────
    # Input Ex2d: (Nx, Ny) → after pad: (Nx+2e, Ny+2e) → .T: (Ny+2e, Nx+2e)
    #                                                         = (dimEy, dimEx)
    def _pad_T(A):
        return np.pad(A, expand, mode='constant').T

    Ex_p = _pad_T(Ex2d)
    Ey_p = _pad_T(Ey2d)
    Ez_p = _pad_T(Ez2d)

    normE0 = np.abs(Ex_p)**2 + np.abs(Ey_p)**2 + np.abs(Ez_p)**2

    # ── 2-D FFT ──────────────────────────────────────────────────────
    _fft2 = lambda A: np.fft.fftshift(np.fft.fft2(A))
    fEx = _fft2(Ex_p)
    fEy = _fft2(Ey_p)
    fEz = _fft2(Ez_p)
    normfE0 = np.abs(fEx)**2 + np.abs(fEy)**2 + np.abs(fEz)**2

    # ── NA mask + phase propagation ───────────────────────────────────
    KX, KY = np.meshgrid(kx, ky)           # shape (dimEy, dimEx)
    K_r2  = KX**2 + KY**2
    mask  = K_r2 <= (k0 * NA)**2
    KZ    = np.sqrt(np.maximum(k0**2 - K_r2, 0)) * mask
    phase = np.exp(1j * KZ * h0_um)

    ffEx = fEx * phase * mask
    ffEy = fEy * phase * mask
    ffEz = fEz * phase * mask

    # ── inverse FFT ───────────────────────────────────────────────────
    _ifft2 = lambda A: np.fft.ifft2(np.fft.ifftshift(A))
    normifE0 = (np.abs(_ifft2(ffEx))**2 +
                np.abs(_ifft2(ffEy))**2 +
                np.abs(_ifft2(ffEz))**2)

    return dict(
        kx=kx, ky=ky,
        x_ext=x_ext, y_ext=y_ext,
        x_um=x_um, y_um=y_um,
        normE0=normE0,
        fEx=fEx, fEy=fEy, fEz=fEz,
        ffEx=ffEx, ffEy=ffEy, ffEz=ffEz,
        normfE0=normfE0, normifE0=normifE0,
        normfEx=np.abs(fEx)**2,
        dimEx=dimEx, dimEy=dimEy,
    )

# ╔══════════════════════════════════════════════════════════════╗
# ║                    VISUALIZATION                             ║
# ╚══════════════════════════════════════════════════════════════╝

def visualize(res: dict, k0: float, NA: float, a_um: float,
              fname_base: str, output_dir: str,
              display_unit: str = 'um',
              plot_overview: bool = True,
              plot_polarization: bool = True,
              plot_cutlines: bool = True,
              plot_gauss: bool = True):
    """
    Generate and save static PNG figures in four groups:
      Fig A: 2×2 overview  (r-space, k-space log, k-space |Ex|², iFFT)
      Fig B: polarization polar plot
      Fig C: 3×2 cutlines  (k-cuts, θ-cuts, radiation patterns)
      Fig D: Gaussian fit on ky-cutline

    All figures are saved as PNG to output_dir.
    """
    kx = res['kx'];  ky = res['ky']
    x_ext = res['x_ext'];  y_ext = res['y_ext']
    x_um  = res['x_um'];   y_um  = res['y_um']
    normE0   = res['normE0']
    normfE0  = res['normfE0']
    normifE0 = res['normifE0']
    normfEx  = res['normfEx']
    ffEx = res['ffEx'];  ffEy = res['ffEy']
    dimEx = res['dimEx'];  dimEy = res['dimEy']

    s    = to_um(display_unit)
    xu   = display_unit
    ku   = f'1/{xu}'

    x_d      = x_ext / s
    y_d      = y_ext / s
    x_orig_d = x_um  / s
    y_orig_d = y_um  / s
    kx_d = kx * s
    ky_d = ky * s
    k0_d = k0 * s

    a_d   = a_um / s
    kD_d  = 4 * np.pi / (3 * a_d)
    G_d   = kD_d * np.sqrt(3)
    Gx_d  = G_d * np.cos(np.pi / 6)
    Gy_d  = G_d * np.sin(np.pi / 6)
    R_scl = 25 * a_d
    K_scl = 1.6 * G_d
    BZ_OFS = [(0, 0), (0, G_d), (Gx_d, Gy_d), (-Gx_d, Gy_d),
              (0, -G_d), (-Gx_d, -Gy_d), (Gx_d, -Gy_d)]

    tt    = np.linspace(0, 2 * np.pi, 200)
    bz_th = np.linspace(0, 2 * np.pi, 7)

    mid_x = dimEx // 2
    mid_y = dimEy // 2

    cut_ky = normfE0[:, mid_x]
    cut_kx = normfE0[mid_y, :]
    mask_ky = np.abs(ky) <= k0
    mask_kx = np.abs(kx) <= k0
    theta_ky = np.degrees(np.arcsin(np.clip(ky[mask_ky] / k0, -1, 1)))
    theta_kx = np.degrees(np.arcsin(np.clip(kx[mask_kx] / k0, -1, 1)))
    vals_ky  = (cut_ky / cut_ky.max())[mask_ky]
    vals_kx  = (cut_kx / cut_kx.max())[mask_kx]

    os.makedirs(output_dir, exist_ok=True)

    def out(suffix):
        return os.path.join(output_dir, fname_base + suffix + '.png')

    def _save(fig, path):
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved: {path}')

    # imshow extents: [xmin, xmax, ymin, ymax]
    r_ext = [x_d.min(), x_d.max(), y_d.min(), y_d.max()]
    k_ext = [kx_d.min(), kx_d.max(), ky_d.min(), ky_d.max()]

    # ══════════════════════════════════════════════════════════════════
    # Fig A: 2×2 overview
    # ══════════════════════════════════════════════════════════════════
    if plot_overview:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(fname_base, fontsize=13)

        # (0,0) r-space |E|²
        ax = axes[0, 0]
        im = ax.imshow(normE0 / normE0.max(),
                       origin='lower', extent=r_ext, aspect='equal',
                       cmap='hot', vmin=0, vmax=1, interpolation='bilinear')
        ax.set_xlim(x_orig_d.min(), x_orig_d.max())
        ax.set_ylim(y_orig_d.min(), y_orig_d.max())
        ax.set_xlabel(f'x ({xu})');  ax.set_ylabel(f'y ({xu})')
        ax.set_title('r-space |E|² (linear)')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # (0,1) iFFT near-field
        ax = axes[0, 1]
        im = ax.imshow(normifE0 / normifE0.max(),
                       origin='lower', extent=r_ext, aspect='equal',
                       cmap='hot', vmin=0, vmax=1, interpolation='bilinear')
        ax.set_xlim(x_orig_d.min(), x_orig_d.max())
        ax.set_ylim(y_orig_d.min(), y_orig_d.max())
        ax.set_xlabel(f'x ({xu})');  ax.set_ylabel(f'y ({xu})')
        ax.set_title(f'r-space |E|² – iFFT (NA={NA})')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # (1,0) k-space |Ex|²
        ax = axes[1, 0]
        im = ax.imshow(normfEx / normfE0.max(),
                       origin='lower', extent=k_ext, aspect='equal',
                       cmap='hot', interpolation='bilinear')
        ax.set_xlim(-k0_d * 1.1, k0_d * 1.1)
        ax.set_ylim(-k0_d * 1.1, k0_d * 1.1)
        ax.plot(k0_d * np.cos(tt), k0_d * np.sin(tt), 'w-', lw=1.5)
        ax.plot(NA * k0_d * np.cos(tt), NA * k0_d * np.sin(tt), 'w--', lw=1.5)
        ax.set_xlabel(f'kx ({ku})');  ax.set_ylabel(f'ky ({ku})')
        ax.set_title(f'k-space |Ex|² (NA={NA})')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # (1,1) k-space log
        ax = axes[1, 1]
        with np.errstate(divide='ignore'):
            klog = np.log10(normfE0 / normfE0.max())
        im = ax.imshow(klog, origin='lower', extent=k_ext, aspect='equal',
                       cmap='hot', vmin=-5, vmax=0, interpolation='bilinear')
        ax.set_xlim(-K_scl, K_scl);  ax.set_ylim(-K_scl, K_scl)
        for r_c, ls_c in [(k0_d, '-'), (NA * k0_d, '--')]:
            ax.plot(r_c * np.cos(tt), r_c * np.sin(tt), 'w', ls=ls_c, lw=1.5)
        for ox, oy in BZ_OFS:
            ax.plot(kD_d * np.cos(bz_th) + ox,
                    kD_d * np.sin(bz_th) + oy, color='cyan', lw=1)
        ax.set_xlabel(f'kx ({ku})');  ax.set_ylabel(f'ky ({ku})')
        ax.set_title('k-space (log)')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label='log₁₀(|fE|²/max)')

        fig.tight_layout()
        _save(fig, out('_A_overview'))

    # ══════════════════════════════════════════════════════════════════
    # Fig B: polarization polar plot
    # ══════════════════════════════════════════════════════════════════
    if plot_polarization:
        pol_th   = np.linspace(0, 2 * np.pi, 73)
        pol_Esum = np.array([
            np.sum(np.abs(ffEy * np.cos(t) - ffEx * np.sin(t))**2)
            for t in pol_th
        ])
        pol_dir = np.degrees(pol_th[np.argmax(pol_Esum)])
        pol_deg = ((pol_Esum.max() - pol_Esum.min()) /
                   abs(pol_Esum.max() + pol_Esum.min()))
        print(f'  Polarization degree = {pol_deg:.4f},  dominant θ = {pol_dir:.1f}°')

        fig, ax = plt.subplots(figsize=(5.5, 5.5),
                               subplot_kw={'projection': 'polar'})
        ax.plot(pol_th, pol_Esum / pol_Esum.max(), 'r-o', ms=3, lw=1.5)
        ax.set_theta_direction(-1)          # clockwise
        ax.set_theta_zero_location('N')     # 0° at top
        ax.set_title(
            f'{fname_base}\nPolarization  θ={pol_dir:.1f}°  degree={pol_deg:.3f}',
            fontsize=10, pad=20)
        fig.tight_layout()
        _save(fig, out('_B_polarization'))

    # ══════════════════════════════════════════════════════════════════
    # Fig C: 3×2 cutlines + radiation patterns
    # ══════════════════════════════════════════════════════════════════
    if plot_cutlines:
        fig, axes = plt.subplots(3, 2, figsize=(11, 12))
        fig.suptitle(fname_base, fontsize=12)
        lc = 'r'

        # (0,0) ky-cutline
        ax = axes[0, 0]
        ax.plot(ky_d, cut_ky / cut_ky.max(), lc, lw=1.5)
        ax.set_xlim(-k0_d * NA, k0_d * NA);  ax.set_ylim(0, 1)
        ax.set_xlabel(f'ky ({ku})');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('ky-cutline (kx=0)')

        # (0,1) kx-cutline
        ax = axes[0, 1]
        ax.plot(kx_d, cut_kx / cut_kx.max(), lc, lw=1.5)
        ax.set_xlim(-k0_d * NA, k0_d * NA);  ax.set_ylim(0, 1)
        ax.set_xlabel(f'kx ({ku})');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('kx-cutline (ky=0)')

        # (1,0) θ ky-cut
        ax = axes[1, 0]
        ax.plot(theta_ky, vals_ky, lc, lw=1.5)
        ax.set_xlim(-10, 10);  ax.set_ylim(0, 1)
        ax.set_xlabel('θ (deg)');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('ky-cut: θ vs |fE|²')

        # (1,1) θ kx-cut
        ax = axes[1, 1]
        ax.plot(theta_kx, vals_kx, lc, lw=1.5)
        ax.set_xlim(-10, 10);  ax.set_ylim(0, 1)
        ax.set_xlabel('θ (deg)');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('kx-cut: θ vs |fE|²')

        # (2,0) radiation pattern y-cut
        ax = axes[2, 0]
        ax.plot(theta_ky, vals_ky / vals_ky.max(), lc, lw=1.5)
        ax.set_xlim(-90, 90);  ax.set_ylim(0, 1)
        ax.set_xlabel('θ (deg)');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('Radiation pattern — y-cut')

        # (2,1) radiation pattern x-cut
        ax = axes[2, 1]
        ax.plot(theta_kx, vals_kx / vals_kx.max(), lc, lw=1.5)
        ax.set_xlim(-90, 90);  ax.set_ylim(0, 1)
        ax.set_xlabel('θ (deg)');  ax.set_ylabel('Intensity (norm.)')
        ax.set_title('Radiation pattern — x-cut')

        fig.tight_layout()
        _save(fig, out('_C_cutlines'))

    # ══════════════════════════════════════════════════════════════════
    # Fig D: Gaussian fit on ky-cutline
    # ══════════════════════════════════════════════════════════════════
    if plot_gauss:
        try:
            xdata = ky_d
            ydata = cut_ky / cut_ky.max()
            gauss = lambda x, a, b, c: a * np.exp(-((x - b) / c)**2)
            popt, _ = curve_fit(gauss, xdata, ydata, p0=[1, 0, 0.1 * s],
                                maxfev=10000)
            print(f'  Gaussian fit (ky): a={popt[0]:.4f}, '
                  f'b={popt[1]:.6f}, c={popt[2]:.4f} ({ku})')
            xfit = np.linspace(-0.5 * s, 0.5 * s, 1000)

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(xdata, ydata, 'b.', ms=2, label='data')
            ax.plot(xfit, gauss(xfit, *popt), 'r-', lw=2, label='Gaussian fit')
            ax.set_xlim(-0.4 * s, 0.4 * s);  ax.set_ylim(0, 1)
            ax.set_xlabel(f'ky ({ku})');  ax.set_ylabel('Intensity (norm.)')
            ax.set_title(f'{fname_base}  —  Gaussian fit (ky-cut)\n'
                         f'a={popt[0]:.4f},  b={popt[1]:.4e},  '
                         f'c={popt[2]:.4f} {ku}')
            ax.legend()
            fig.tight_layout()
            _save(fig, out('_D_gauss_fit'))
        except Exception as e:
            print(f'  Gaussian fit failed: {e}')

# ╔══════════════════════════════════════════════════════════════╗
# ║                    MAIN                                      ║
# ╚══════════════════════════════════════════════════════════════╝

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # ── collect files ─────────────────────────────────────────────────
    if SELECT_FILES is None:
        files = sorted(glob.glob(os.path.join(script_dir, '*.txt')))
    else:
        files = [os.path.join(script_dir, f) for f in SELECT_FILES]

    if not files:
        print('No .txt files found.')
        return

    # ── output directory ──────────────────────────────────────────────
    out_dir = os.path.join(script_dir, OUTPUT)
    os.makedirs(out_dir, exist_ok=True)

    # ── process each file ─────────────────────────────────────────────
    for fpath in files:
        fname_base = os.path.splitext(os.path.basename(fpath))[0]
        print(f'\n== {fname_base} ==')

        # 1. Load
        t0 = time.time()
        data, meta = load_comsol(fpath)
        print(f'  Loaded {len(data["x"])} pts  '
              f'unit={meta["length_unit"]}  has_z={meta["has_z"]}  '
              f'({time.time() - t0:.3f}s)')

        # ── resolve FREQ_THz and H0_UM ───────────────────────────────────
        if AUTO_PARAM:
            freq = meta['freq_thz']
            if freq is None:
                raise ValueError(f'Could not parse frequency from {fpath}. '
                                 f'Set AUTO_PARAM=False and provide FREQ_THz manually.')
            if data['z'] is not None:
                h0_um = float(np.mean(data['z'])) * to_um(meta['length_unit'])
            else:
                raise ValueError(f'No z column found in {fpath} to determine H0. '
                                 f'Set AUTO_PARAM=False and provide H0_UM manually.')
        else:
            freq  = FREQ_THz
            h0_um = H0_UM
        lda0_um = 299.792458 / freq
        k0_um   = 2 * np.pi / lda0_um

        print(f'  freq={freq:.4f} THz  lambda={lda0_um:.4f} um  '
              f'k0={k0_um:.4f} 1/um  h0={h0_um:.4f} um')

        # 2. Matrix
        t0 = time.time()
        Ex2d, Ey2d, Ez2d, x_um, y_um, dx, dy, Nx, Ny = matrix_field(data, meta)
        print(f'  Grid {Nx}x{Ny}  dx={dx:.4f} um  dy={dy:.4f} um  '
              f'({time.time() - t0:.3f}s)')

        # 3. FFT pipeline
        t0 = time.time()
        res = compute_fft(Ex2d, Ey2d, Ez2d, x_um, y_um, dx, dy, Nx, Ny,
                          k0_um, h0_um, NA, EXPAND)
        print(f'  FFT done ({time.time() - t0:.3f}s)')

        # 4. Visualize → save PNG
        visualize(res, k0_um, NA, A_UM, fname_base, out_dir,
                  display_unit='um',
                  plot_overview=PLOT_OVERVIEW,
                  plot_polarization=PLOT_POLARIZATION,
                  plot_cutlines=PLOT_CUTLINES,
                  plot_gauss=PLOT_GAUSS)


if __name__ == '__main__':
    main()
