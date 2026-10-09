import matplotlib.pyplot as plt
import numpy as np

from scripts.analysis import plot_s5_volume_radiation as plotting
from scripts.analysis.validate_s5_bloch_radiation import pairs


def test_cavity_plots_preserve_shared_intensity_scale_and_scope(monkeypatch):
    saved = {}
    monkeypatch.setattr(plotting, 'save', lambda fig, name: saved.setdefault(name, fig))
    axis = np.array([-.08, 0., .08])
    contributions = np.zeros((3, 3, 19, 1, 2), complex)
    contributions[:, :, 0, 0, 0] = .4
    contributions[:, :, 1, 0, 0] = .1
    contributions[:, :, 0, 0, 1] = .2
    cube = np.zeros((19, 4, 2), complex)
    cube[0, 1] = [.4, .2]
    direct = np.array([0j, 1+0j])
    groups = {'Gamma_py': cube[0, 1], 'nonGamma_py': np.zeros(2),
              'other_bands': np.zeros(2), 'remainder': direct-cube[0, 1]}
    summary = {'basis_support': plotting.CAVITY_SUPPORT,
               'terms': pairs(cube), 'direct_air_F': pairs(direct), 'full_volume_F': pairs(direct),
               'basis_over_air': pairs(cube[0, 1]), 'full_vector_L2_residual': .15,
               'electric_L2_residual_inside_cavity': .15, 'volume_vs_air_relative_error': .1,
               'groups': {name: pairs(value) for name, value in groups.items()}}
    d = {'k19': np.column_stack([np.linspace(-.06, .06, 19), np.zeros(19)]), 'scale': 1., 'N': 19}
    try:
        plotting.plot(d, axis, contributions, cube, direct, direct, summary,
                      np.broadcast_to(direct, (3, 3, 2)))
        images = []
        for name in ['12_s5_volume_Ex_19_envelopes', '13_s5_volume_Ey_19_envelopes']:
            fig = saved[name]
            assert fig.axes[0].get_title() == 'Index map'
            assert 'cavity cells only' in fig._suptitle.get_text()
            assert not fig.axes[0].images
            panels = [ax.images[0] for ax in fig.axes if ax.images]
            assert len(panels) == 19
            assert all(im.get_cmap().name == 'magma' for im in panels)
            images.extend(panels)
        assert all(im.get_clim() == images[0].get_clim() for im in images)
        np.testing.assert_allclose(images[0].get_array().max(), .4**2)
        np.testing.assert_allclose(images[1].get_array().max(), .1**2)
        np.testing.assert_allclose(images[19].get_array().max(), .2**2)
        origin = saved['18_s5_nonGamma_py_interference']
        assert all('exact cavity-window null' in ax.get_title() for ax in origin.axes[1:])
        assert all(not ax.patches for ax in origin.axes[1:])
    finally:
        for fig in saved.values():
            plt.close(fig)
