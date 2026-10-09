import numpy as np
from matplotlib.colors import Normalize
from scripts.analysis import analyze_s5_cladding_scattering as analysis


def test_gs2b_row_scales_and_selective_output(monkeypatch, tmp_path):
    field = np.ones((22, 3, 3, 3), dtype=complex)
    field[..., 1] *= 20
    calls = []

    def inspect(fig, out, name):
        calls.append(name)
        axes = fig.axes[:12]
        limits = [ax.collections[0].get_clim() for ax in axes]
        assert len(fig.axes) == 14
        assert len(set(limits[:6])) == len(set(limits[6:])) == 1
        assert np.isclose(limits[6][1] / limits[0][1], 400)
        assert all(lo == 0 for lo, hi in limits)
        assert all(type(ax.collections[0].norm) is Normalize for ax in axes)
        assert all(ax.collections[0].cmap.name == 'magma' for ax in axes)
        analysis.plt.close(fig)

    monkeypatch.setattr(analysis, 'save', inspect)
    analysis.plot_angles(tmp_path, dict(F_V=field, q_m_inv=np.linspace(-1e5, 1e5, 3), frequency_Hz=2e14), include_cuts=False)
    assert calls == ['GS2b_cumulative_farfields']
