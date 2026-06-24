import numpy as np
import matplotlib.pyplot as plt


def reciprocal_vectors(a1, a2):
    """Return b1, b2 satisfying ai · bj = 2π δij."""
    a = np.array([a1, a2], dtype=float).T
    b = 2 * np.pi * np.linalg.inv(a).T
    return b[:, 0], b[:, 1]


def bloch_hamiltonian(k, t0, t1, alpha, beta, gamma, mu):
    """
    6-site hexagonal unit-cell Hamiltonian with C2v modulation.

    Internal index (0-based) -> figure label (1-based):
        0 <-> 1  (orange, on-site +2mu)
        1 <-> 2  (blue,   on-site  -mu)
        2 <-> 3  (blue,   on-site  -mu)
        3 <-> 4  (orange, on-site +2mu)
        4 <-> 5  (blue,   on-site  -mu)
        5 <-> 6  (blue,   on-site  -mu)

    Hexagon layout (intracell):
          5 -- 4
         /      \\
        0        3
         \\      /
          1 -- 2

    Intracell:
        orange-blue  (0,1),(5,0),(2,3),(3,4)  ->  t0*beta
        blue-blue    (1,2),(4,5)              ->  t0*gamma

    Intercell:
        orange-orange  0 <-> 3  along  a1           ->  t1*alpha * exp(ik.a1)
        blue-blue      1 <-> 4  along  a2           ->  t1*gamma * exp(ik.a2)
        blue-blue      2 <-> 5  along  a3 = a2-a1  ->  t1*gamma * exp(ik.a3)

    Matches H_site in the reference figure exactly.
    """
    k = np.asarray(k, dtype=float)
    h = np.zeros((6, 6), dtype=complex)

    # On-site energies: orange +2mu, blue -mu
    h += np.diag(np.array([2*mu, -mu, -mu, 2*mu, -mu, -mu], dtype=float))

    def add_hop(i, j, amp, r=(0.0, 0.0)):
        phase = np.exp(1j * np.dot(k, r))
        h[i, j] += amp * phase
        h[j, i] += np.conj(amp * phase)

    a1 = np.array([ 1.0,             0.0])
    a2 = np.array([ 0.5, np.sqrt(3.0)/2.0])
    a3 = a2 - a1   # upper-left lattice vector

    # Intracell bonds
    for i, j in [(0, 1), (5, 0), (2, 3), (3, 4)]:   # orange-blue: t0*beta
        add_hop(i, j, t0 * beta)
    for i, j in [(1, 2), (4, 5)]:                    # blue-blue:   t0*gamma
        add_hop(i, j, t0 * gamma)

    # Intercell bonds
    add_hop(0, 3, t1 * alpha, a1)   # orange-orange: t1*alpha
    add_hop(1, 4, t1 * gamma, a2)   # blue-blue:     t1*gamma
    add_hop(2, 5, t1 * gamma, a3)   # blue-blue:     t1*gamma

    return h


def interpolate_k_path(points, n_per_segment):
    """Linearly interpolate a piecewise k-path through high-symmetry points."""
    k_list, x_list, ticks = [], [], [0.0]
    distance = 0.0
    for start, end in zip(points[:-1], points[1:]):
        start = np.asarray(start, dtype=float)
        end   = np.asarray(end,   dtype=float)
        seg     = end - start
        seg_len = np.linalg.norm(seg)
        for m in range(n_per_segment):
            u = m / n_per_segment
            k_list.append(start + u * seg)
            x_list.append(distance + u * seg_len)
        distance += seg_len
        ticks.append(distance)
    k_list.append(np.asarray(points[-1], dtype=float))
    x_list.append(distance)
    return np.array(k_list), np.array(x_list), np.array(ticks)


def calculate_bands(t0, t1, alpha, beta, gamma, mu, n_per_segment):
    """Compute band structure along Gamma-K-M-Gamma."""
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3.0) / 2.0])
    b1, b2 = reciprocal_vectors(a1, a2)

    G = np.array([0.0, 0.0])
    K = (2.0 * b1 + b2) / 3.0   # [4pi/3, 0]      — along kx
    M = b2 / 2.0                  # [0, 2pi/sqrt3]  — along ky

    k_path, x_path, ticks = interpolate_k_path(
        [M, G, K], n_per_segment=n_per_segment
    )
    bands = np.array([
        np.linalg.eigvalsh(
            bloch_hamiltonian(k, t0=t0, t1=t1,
                              alpha=alpha, beta=beta, gamma=gamma, mu=mu)
        )
        for k in k_path
    ])
    labels = ["M", r"$\Gamma$", "K"]
    return x_path, bands, ticks, labels


def plot_bands(t0, t1, alpha, beta, gamma, mu, n_per_segment):
    x_path, bands, ticks, labels = calculate_bands(
        t0=t0, t1=t1, alpha=alpha, beta=beta, gamma=gamma, mu=mu,
        n_per_segment=n_per_segment,
    )

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for n in range(bands.shape[1]):
        ax.plot(x_path, bands[:, n], color="black", lw=1.4)

    for x in ticks:
        ax.axvline(x, color="0.82", lw=0.8)

    ax.set_xlim(x_path[0], x_path[-1])
    ax.set_ylim(-2.0, 2.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Energy")
    ax.set_title(
        rf"t0={t0:.3g}, t1={t1:.3g},  "
        rf"$\alpha$={alpha:.3g}, $\beta$={beta:.3g}, $\gamma$={gamma:.3g},  "
        rf"$\mu$={mu:.3g}"
    )
    ax.grid(axis="y", color="0.9", lw=0.6)
    fig.tight_layout()
    return fig, ax


def _couplings_from(eta, zeta):
    """Derive all coupling parameters from the two physical control parameters."""
    t0    = 1-(eta-1)/2
    t1    = 1+(eta-1)/2
    alpha = zeta
    beta  = 1+(zeta-1)*0.6
    gamma = 1-(zeta-1)*0.5
    mu    = (zeta - 1)*0.2
    return t0, t1, alpha, beta, gamma, mu


def _band_title(t0, t1, alpha, beta, gamma, mu):
    return (
        rf"t0={t0:.3g}, t1={t1:.3g},  "
        rf"$\alpha$={alpha:.3g}, $\beta$={beta:.3g}, $\gamma$={gamma:.3g},  "
        rf"$\mu$={mu:.3g}"
    )


def interactive_plot(eta_init=1.0, zeta_init=1.0, n_per_segment=120,
                     eta_range=(0.80, 1.20), zeta_range=(0.50, 2.00)):
    """Band structure plot with real-time eta and zeta sliders."""
    from matplotlib.widgets import Slider

    t0, t1, alpha, beta, gamma, mu = _couplings_from(eta_init, zeta_init)
    x_path, bands, ticks, labels = calculate_bands(
        t0=t0, t1=t1, alpha=alpha, beta=beta, gamma=gamma, mu=mu,
        n_per_segment=n_per_segment,
    )

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    fig.subplots_adjust(bottom=0.22)

    lines = [ax.plot(x_path, bands[:, n], color="black", lw=1.4)[0]
             for n in range(bands.shape[1])]

    for x in ticks:
        ax.axvline(x, color="0.82", lw=0.8)

    ax.set_xlim(x_path[0], x_path[-1])
    ax.set_ylim(-2.0, 2.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Energy")
    title = ax.set_title(_band_title(t0, t1, alpha, beta, gamma, mu))
    ax.grid(axis="y", color="0.9", lw=0.6)

    # ── Sliders ──────────────────────────────────────────────────────────
    ax_eta  = fig.add_axes([0.15, 0.10, 0.70, 0.03])
    ax_zeta = fig.add_axes([0.15, 0.04, 0.70, 0.03])

    slider_eta  = Slider(ax_eta,  r'$\eta$',
                         eta_range[0],  eta_range[1],  valinit=eta_init,  valstep=0.005)
    slider_zeta = Slider(ax_zeta, r'$\zeta$',
                         zeta_range[0], zeta_range[1], valinit=zeta_init, valstep=0.005)

    def update(_):
        eta  = slider_eta.val
        zeta = slider_zeta.val
        t0, t1, alpha, beta, gamma, mu = _couplings_from(eta, zeta)
        _, new_bands, _, _ = calculate_bands(
            t0=t0, t1=t1, alpha=alpha, beta=beta, gamma=gamma, mu=mu,
            n_per_segment=n_per_segment,
        )
        for n, line in enumerate(lines):
            line.set_ydata(new_bands[:, n])
        title.set_text(_band_title(t0, t1, alpha, beta, gamma, mu))
        fig.canvas.draw_idle()

    slider_eta.on_changed(update)
    slider_zeta.on_changed(update)


def interactive_plot_full(t0_init=1.0, t1_init=1.0,
                           alpha_init=1.0, beta_init=1.0,
                           gamma_init=1.0, mu_init=0.0,
                           n_per_segment=120):
    """Band structure with independent sliders for all 6 coupling parameters."""
    from matplotlib.widgets import Slider

    x_path, bands, ticks, labels = calculate_bands(
        t0=t0_init, t1=t1_init, alpha=alpha_init, beta=beta_init,
        gamma=gamma_init, mu=mu_init, n_per_segment=n_per_segment,
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    fig.subplots_adjust(bottom=0.34)
    fig.canvas.manager.set_window_title("Full parameter control")

    lines = [ax.plot(x_path, bands[:, n], color="C0", lw=1.4)[0]
             for n in range(bands.shape[1])]

    for x in ticks:
        ax.axvline(x, color="0.82", lw=0.8)

    ax.set_xlim(x_path[0], x_path[-1])
    ax.set_ylim(-2.0, 2.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Energy")
    title = ax.set_title(_band_title(t0_init, t1_init, alpha_init,
                                     beta_init, gamma_init, mu_init))
    ax.grid(axis="y", color="0.9", lw=0.6)

    # Two-column slider layout
    #   Left  column: t0, alpha, beta
    #   Right column: t1, gamma, mu
    lx, rx, sw, sh = 0.08, 0.57, 0.36, 0.03
    rows = [0.05, 0.13, 0.22]   # bottom → top

    def make_slider(x, y, label, vmin, vmax, vinit):
        return Slider(fig.add_axes([x, y, sw, sh]),
                      label, vmin, vmax, valinit=vinit, valstep=0.01)

    sl_t0    = make_slider(lx, rows[2], r'$t_0$',      0.50, 2.0,  t0_init)
    sl_alpha = make_slider(lx, rows[1], r'$\alpha$',   0.30, 3.0,  alpha_init)
    sl_beta  = make_slider(lx, rows[0], r'$\beta$',    0.30, 3.0,  beta_init)
    sl_t1    = make_slider(rx, rows[2], r'$t_1$',      0.50, 2.0,  t1_init)
    sl_gamma = make_slider(rx, rows[1], r'$\gamma$',   0.30, 3.0,  gamma_init)
    sl_mu    = make_slider(rx, rows[0], r'$\mu$',     -1.00, 1.0,  mu_init)

    def update(_):
        t0    = sl_t0.val
        t1    = sl_t1.val
        alpha = sl_alpha.val
        beta  = sl_beta.val
        gamma = sl_gamma.val
        mu    = sl_mu.val
        _, new_bands, _, _ = calculate_bands(
            t0=t0, t1=t1, alpha=alpha, beta=beta, gamma=gamma, mu=mu,
            n_per_segment=n_per_segment,
        )
        for n, line in enumerate(lines):
            line.set_ydata(new_bands[:, n])
        title.set_text(_band_title(t0, t1, alpha, beta, gamma, mu))
        fig.canvas.draw_idle()

    for sl in (sl_t0, sl_t1, sl_alpha, sl_beta, sl_gamma, sl_mu):
        sl.on_changed(update)


if __name__ == "__main__":
    # ── Shared initial values ─────────────────────────────────────────────
    eta_init  = 1.0
    zeta_init = 1.0
    nk        = 120
    # ─────────────────────────────────────────────────────────────────────

    # Window 1: η / ζ control only (black bands)
    interactive_plot(eta_init=eta_init, zeta_init=zeta_init, n_per_segment=nk,
                     eta_range=(0.80, 1.20), zeta_range=(0.50, 2.00))

    # Window 2: full independent parameter control (blue bands)
    # Initial values are derived from the same η, ζ for easy comparison
    t0, t1, alpha, beta, gamma, mu = _couplings_from(eta_init, zeta_init)
    interactive_plot_full(t0_init=t0, t1_init=t1,
                          alpha_init=alpha, beta_init=beta,
                          gamma_init=gamma, mu_init=mu,
                          n_per_segment=nk)

    plt.show()
