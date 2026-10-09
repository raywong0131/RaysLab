"""S5 path and real saver checks without reading scientific data or solving."""
from unittest.mock import patch

import matplotlib.pyplot as plt
import pytest

from comsol_workflow.s5_paths import output_path
from scripts.analysis import visualize_s5_broadening as toy
from scripts.analysis import visualize_s5_broadening_2d as broadening
from scripts.analysis import plot_s5_volume_radiation as volume


def test_s5_layout_and_all_figure_savers(tmp_path):
    expected = {'field.npz': '01_results', 'figure.png': '10_overview',
                'figure.pdf': '11_pdf', 'report.md': '12_reports',
                'table.csv': '80_logs', 'summary.json': '80_logs',
                'run.log': '80_logs', 'old.zip': '90_history',
                's5_volume_preflight.json': '99_config', 'index.html': ''}
    for name, folder in expected.items():
        path = output_path(tmp_path, name)
        assert path == tmp_path / folder / name
        assert path.parent.is_dir()
    for name in ('../outside.npz', 'nested/file.png', ''):
        with pytest.raises(ValueError):
            output_path(tmp_path, name)
    savers = [lambda fig: toy.save_figure(fig, tmp_path, 'toy'),
              lambda fig: broadening.save(fig, tmp_path, 'broadening'),
              lambda fig: volume.save(fig, 'volume')]
    with patch.object(volume, 'OUTPUT', tmp_path):
        for save in savers:
            fig, ax = plt.subplots(figsize=(2, 2))
            ax.plot([0, 1], [0, 1])
            save(fig)
    for stem in ('toy', 'broadening', 'volume'):
        assert (tmp_path / '10_overview' / f'{stem}.png').read_bytes().startswith(b'\x89PNG')
        assert (tmp_path / '11_pdf' / f'{stem}.pdf').read_bytes().startswith(b'%PDF')
    assert not list(tmp_path.glob('*.png')) and not list(tmp_path.glob('*.pdf'))
