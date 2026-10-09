"""Integration contracts, using synthetic artifacts and no solver startup."""
import ast
from importlib import import_module
import json
import os
from pathlib import Path
import sys
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pytest

from comsol_workflow.dipolar_inputs import INPUT_MANIFEST_ENV, PROJECT_ROOT, load_dipolar_inputs
from comsol_workflow.figure_output import save_figure_formats
from comsol_workflow.output_paths import artifact_path, output_path, table_links, table_path
from comsol_workflow.result_io import read_json, write_csv, write_json
from scripts import workflow_cli


def test_all_legacy_modules_share_canonical_state():
    pairs = []
    for directory in ('comsol_workflow', 'scripts/analysis', 'scripts/run_main', 'scripts/run_sweep'):
        for file in (PROJECT_ROOT / directory).glob('*.py'):
            source = file.read_text(encoding='utf-8-sig')
            if 'redirect(__name__,' not in source:
                continue
            tree = ast.parse(source)
            call = next(n for n in ast.walk(tree) if isinstance(n, ast.Call)
                        and isinstance(n.func, ast.Name) and n.func.id == 'redirect')
            pairs.append((file.relative_to(PROJECT_ROOT).with_suffix('').as_posix().replace('/', '.'),
                          ast.literal_eval(call.args[1])))
    assert len(pairs) == 39
    for old_name, new_name in pairs:
        old, new = import_module(old_name), import_module(new_name)
        assert old is new
        marker = object()
        old._integration_probe = marker
        try:
            assert new._integration_probe is marker
        finally:
            del old._integration_probe


def test_input_defaults_and_explicit_relative_sources(tmp_path, monkeypatch):
    monkeypatch.delenv(INPUT_MANIFEST_ENV, raising=False)
    historical = load_dipolar_inputs()
    assert historical.mode_dir != historical.field_mode_dir
    assert historical.fourier_source.parent.parent == historical.mode_dir
    manifest = tmp_path / 'sources.json'
    manifest.write_text(json.dumps(dict(schema_version=1, series_dir='saved', mode_dir='overview/mode19',
        field_mode_dir='raw/mode19', source_config_dir='saved/config', finite_model='saved/model.mph',
        output_dir=str(historical.output_dir / 'synthetic'), run_id='new-run')), encoding='utf-8')
    before = set(tmp_path.rglob('*'))
    custom = load_dipolar_inputs(manifest)
    assert custom.series_dir == tmp_path / 'saved'
    assert custom.field_mode_dir == tmp_path / 'raw/mode19'
    assert custom.fourier_source == tmp_path / 'overview/mode19/13_lattice_fourier_Hz/finite_lattice_fourier_hz.npz'
    assert set(tmp_path.rglob('*')) == before
    assert not custom.output_dir.exists()
    monkeypatch.setenv(INPUT_MANIFEST_ENV, str(manifest))
    assert load_dipolar_inputs() == custom


def test_manifest_reaches_scientific_entrypoints_in_fresh_process(tmp_path):
    manifest = tmp_path / 'sources.json'
    manifest.write_text(json.dumps(dict(schema_version=1, series_dir='saved', mode_dir='overview',
        field_mode_dir='raw', source_config_dir='config', finite_model='source.mph',
        reference_dir='reference', boundary_case_dir='boundary', run_id='source-check')),
        encoding='utf-8')
    program = '''
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
from scripts.analysis import analyze_dipolar_singularity as analysis
from scripts.analysis import visualize_dipolar_broadening as fields
from scripts.analysis import visualize_dipolar_broadening_2d as fourier
from scripts.analysis import prepare_dipolar_gamma_origin as gamma
from scripts.run_main import run_dipolar_volume_export as volume
from scripts.run_main import run_dipolar_periodic_reference as periodic
s = load_dipolar_inputs()
assert analysis.MODE == s.mode_dir and fields.MODE == s.field_mode_dir
assert analysis.CONFIG == s.source_config_dir and fourier.SOURCE == s.fourier_source
assert volume.FINITE == s.finite_model and volume.CONFIG == s.source_config_dir / 'config.json'
assert periodic.OUT == s.reference_dir and periodic.S4 == s.boundary_case_dir
assert gamma.RUN == s.run_id and gamma.RESULT == s.output_dir
assert gamma.paths()['01_results'] == s.output_dir / '01_results' / s.gamma_group / s.run_id
assert gamma.input_path('01_results', 'A.npz') == s.output_dir / '01_results' / s.source_batch / 'A.npz'
import sys
assert 'mph' not in sys.modules
'''
    subprocess.run([sys.executable, '-B', '-c', program], cwd=PROJECT_ROOT,
                   env={**os.environ, INPUT_MANIFEST_ENV: str(manifest)}, check=True)
    assert list(tmp_path.iterdir()) == [manifest]


@pytest.mark.parametrize('data', [[], {'schema_version': 2}, {'schema_version': True},
    {'schema_version': 1, 'unknown': 2}, {'schema_version': 1, 'series_dir': 'other'},
    {'schema_version': 1, 'output_dir': '../outside'}, {'schema_version': 1, 'gamma_group': '../escape'},
    {'schema_version': 1, 'run_id': 'unsafe/run'}, {'schema_version': 1, 'mode_dir': None}])
def test_invalid_manifest_fails_without_writes(tmp_path, data):
    manifest = tmp_path / 'sources.json'
    manifest.write_text(json.dumps(data), encoding='utf-8')
    with pytest.raises((ValueError, TypeError)):
        load_dipolar_inputs(manifest)
    assert list(tmp_path.iterdir()) == [manifest]


def test_router_forwarding_describe_and_restoration(tmp_path, monkeypatch, capsys):
    saved_argv = sys.argv
    monkeypatch.setenv(INPUT_MANIFEST_ENV, 'prior-value')
    monkeypatch.setenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT', 'prior-layout')
    manifest = tmp_path / 'sources.json'
    manifest.write_text('{"schema_version": 1, "run_id": "synthetic"}', encoding='utf-8')
    calls = []
    def run(module, **options):
        calls.append((module, sys.argv[:], os.environ[INPUT_MANIFEST_ENV],
                      os.environ['COMSOL_WORKFLOW_OUTPUT_LAYOUT'], options))
        raise RuntimeError('simulated target failure')
    monkeypatch.setattr(workflow_cli, 'run_module', run)
    with pytest.raises(RuntimeError, match='simulated'):
        workflow_cli.dipolar_analysis_main(['broadening-2d', '--source-manifest', str(manifest),
                                           '--source', 'target-field.npz', '--check'])
    assert calls == [('scripts.analysis.visualize_dipolar_broadening_2d',
        ['scripts.analysis.visualize_dipolar_broadening_2d', '--source', 'target-field.npz', '--check'], str(manifest), 'unified',
        {'run_name': '__main__', 'alter_sys': True})]
    workflow_cli.dipolar_replay_main(['--source-manifest', str(manifest), '--describe-inputs'])
    assert json.loads(capsys.readouterr().out)['run_id'] == 'synthetic'
    assert len(calls) == 1
    assert sys.argv is saved_argv
    assert os.environ[INPUT_MANIFEST_ENV] == 'prior-value'
    assert os.environ['COMSOL_WORKFLOW_OUTPUT_LAYOUT'] == 'prior-layout'
    assert list(tmp_path.iterdir()) == [manifest]


def test_role_paths_tables_and_legacy_read_fallback(tmp_path, monkeypatch):
    monkeypatch.delenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT', raising=False)
    old = table_path(tmp_path / '80_logs/batch/run', 'weights.csv')
    old.write_text('legacy', encoding='utf-8')
    assert output_path(tmp_path, 'weights.csv').parent.name == '80_logs'
    monkeypatch.setenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT', 'unified')
    assert output_path(tmp_path, 'weights.csv').parent.name == '01_results'
    assert table_path(old.parent, 'weights.csv', existing=True) == old
    new = table_path(old.parent, 'weights.csv')
    assert new == tmp_path / '01_results/batch/run/weights.csv'
    new.write_text('new', encoding='utf-8')
    assert table_path(old.parent, 'weights.csv', existing=True) == new
    monkeypatch.delenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT')
    assert table_path(old.parent, 'weights.csv', existing=True) == old
    old.unlink()
    assert table_path(old.parent, 'weights.csv', existing=True) == new
    assert artifact_path(tmp_path, 'values.json', role='data').parent.name == '01_results'
    with pytest.raises(ValueError):
        artifact_path(tmp_path, 'values.json')
    for group in ('../escape', str(tmp_path.resolve())):
        with pytest.raises(ValueError):
            artifact_path(tmp_path, 'table.csv', group=group)
    read_root = tmp_path / 'not-created'
    table_path(read_root / '80_logs', 'absent.csv', existing=True)
    assert not read_root.exists()
    monkeypatch.setenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT', 'unified')
    assert table_links('`80_logs/values.csv` and (../80_logs/check.json)') == '`01_results/values.csv` and (../80_logs/check.json)'


def test_real_shared_figure_saver_preserves_display(tmp_path):
    fig, ax = plt.subplots(figsize=(3, 2))
    field = ax.imshow([[1, 2], [3, 4]], vmin=0, vmax=5)
    ax.set(xlim=(-.5, 1.5), ylim=(-.5, 1.5), title='Synthetic integration check')
    before = (field.get_clim(), ax.get_xlim(), ax.get_ylim())
    try:
        paths = save_figure_formats(fig, tmp_path, 'synthetic', pdf_group='batch/run', dpi=120,
                                    bbox_inches='tight', facecolor='white')
        assert paths == [tmp_path / '10_overview/synthetic.png', tmp_path / '11_pdf/batch/run/synthetic.pdf']
        assert paths[0].read_bytes().startswith(b'\x89PNG')
        assert paths[1].read_bytes().startswith(b'%PDF')
        assert plt.imread(paths[0]).shape[0] > 150
        assert (field.get_clim(), ax.get_xlim(), ax.get_ylim()) == before
    finally:
        plt.close(fig)


def test_serialization_preserves_complex_columns_and_atomic_json(tmp_path):
    csv_path = tmp_path / 'values.csv'
    write_csv(csv_path, [{'z': np.complex128(2+3j)}, {'label': 'other'}])
    assert csv_path.read_text(encoding='utf-8').splitlines() == ['z_re,z_im,label', '2.0,3.0,', ',,other']
    path = tmp_path / 'config.json'
    write_json(path, {'valid': True})
    with pytest.raises(ValueError):
        write_json(path, {'value': float('nan')})
    assert read_json(path) == {'valid': True}
    assert not path.with_name(path.name + '.tmp').exists()


def test_renamed_broadening_exports_self_contained_html(tmp_path, monkeypatch):
    from scripts.analysis import visualize_dipolar_broadening as program
    source = tmp_path / 'synthetic-source.txt'
    source.write_text('synthetic', encoding='utf-8')
    data = dict(source=str(source), source_sha256=program.digest(source),
                quadrant_sum_relative_error=0., fft_center_relative_error=0.,
                quadrants_raw=[[[1., 0.]]*3]*4, quadrants_normalized=[[[.25, 0.]]*3]*4,
                total_raw=[[4., 0.]]*3, total_normalized=[[1., 0.]]*3,
                abs_Ex_over_Ey=1., annotation='<synthetic>')
    out = tmp_path / 'results/generated'
    monkeypatch.setattr(program, 'ROOT', tmp_path)
    monkeypatch.setattr(program, 'real_diagnostics', lambda: data)
    def figure(directory, stem):
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.plot([0, 1], [0, 1])
        program.save_figure(fig, directory, stem)
    monkeypatch.setattr(program, 'toy_figure', lambda directory: figure(directory, '01_broadening_illustration'))
    monkeypatch.setattr(program, 'real_figure', lambda directory, _data: figure(directory, '02_mode19_normal_phasors'))
    monkeypatch.setenv('COMSOL_WORKFLOW_OUTPUT_LAYOUT', 'unified')
    monkeypatch.setattr(sys, 'argv', ['broadening', '--output', str(out)])
    program.main()
    html = (out / 'index.html').read_text(encoding='utf-8')
    assert '__S5_DATA__' not in html and '<synthetic>' not in html
    assert 'data:image/png;base64,' in html and 'data:application/pdf;base64,' in html
    assert 'download="normal_contributions.csv"' in html
    assert (out / '01_results/normal_contributions.csv').is_file()
    assert program.digest(source) == data['source_sha256']
