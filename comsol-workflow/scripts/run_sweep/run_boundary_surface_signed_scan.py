"""Authorized Gamma zeta scan, slab-surface Hz and signed 128-point boundary sums."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import os
from pathlib import Path
from types import SimpleNamespace
import traceback
import numpy as np
import pandas as pd
from scripts.run_main import run_boundary_analysis as s4
from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies, mesh_identity
from scripts.analysis.decompose_boundary_c6v import field_pair, weighted_parity, symmetry_coefficients
from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.boundary_integrals import hole_edges, line_rule, group_totals, signed_indicator, complex_columns

BASE = s4.ROOT / 'scripts/.out/unit_cell_2D'
ROOT = BASE / 'S4_s_eta0.96_zeta1-1.2_mesh5_edge128_20260918'
CONTROL = BASE / 'S4_surface_eta0.96_zeta1_vs1.156_20260918'
OLD = BASE / 'S4_surface_eta0.96_zeta_scan_20260918'
REFERENCE = OLD / 'zeta1.156'
ZETAS = sorted(set([i / 100 for i in range(100, 121)] + [i / 1000 for i in range(1150, 1161)]))


def validate_fields(directory, zeta):
    data = np.load(directory / 'surface_fields.npz')
    assert float(data['z_nm']) == 100
    holes = data['holes_m']
    np.testing.assert_allclose(np.linalg.norm(holes.mean(axis=1), axis=1), 820e-9*.96/3, rtol=1e-10)
    expected = np.full(6, 245e-9*np.sqrt((3-zeta*zeta)/2))
    expected[[0, 3]] = 245e-9*zeta
    np.testing.assert_allclose(np.linalg.norm(holes[:, 1]-holes[:, 0], axis=1), expected, rtol=1e-10)
    rows = []
    for j, edge in enumerate(hole_edges(holes)):
        xy, w = line_rule(edge, 128)
        np.testing.assert_allclose(xy, data['edge_xy_m_128'][j], atol=1e-20, rtol=1e-12)
        np.testing.assert_allclose(w, data['edge_weights_m_128'][j], atol=1e-25, rtol=1e-12)
        integral = w @ data['edge_Hz_raw_128'][j]
        rows.append({'group': edge.group, 'qx': integral*edge.normal[0], 'qy': integral*edge.normal[1]})
    table = pd.read_csv(directory / 'edge_integrals_raw.csv')
    raw = table[table.order == 128]
    np.testing.assert_allclose([r['qx'] for r in rows], raw.qx_re + 1j*raw.qx_im, rtol=1e-10, atol=1e-25)
    return data, group_totals(rows)


def prepare():
    path = ROOT / '99_config/config.json'
    if path.exists():
        manifest = s4.json_read(path)
        assert manifest['zetas'] == ZETAS
        for c in manifest['cases']:
            assert s4.digest(c['snapshot']) == c['snapshot_sha256']
            for p, h in c.get('source_hashes', {}).items():
                assert s4.digest(p) == h
        return manifest
    source = CONTROL / '99_config/parameter.json'
    cases = []
    for zeta in ZETAS:
        case_dir = ROOT / f'zeta{zeta:g}'
        snapshot = s4.json_read(source)
        snapshot['mesh_size'] = 5
        snapshot['unit_cell_2d'].update(b0_nm=245., eta=.96, zeta=zeta)
        pp = case_dir / '99_config/parameter.json'
        s4.write_json(pp, snapshot)
        params = s4.load_shared_parameters(pp)
        assert params.unit_cell_2d.eigenmode_pair_count == 2
        cached = CONTROL if zeta == 1 else OLD / f'zeta{zeta:g}' if 1.154 <= zeta <= 1.158 else None
        case = {'zeta': zeta, 'snapshot': str(pp), 'snapshot_sha256': s4.digest(pp), 'directory': str(case_dir)}
        if cached is not None:
            metadata = cached / ('execution.json' if zeta == 1 else 'export.json')
            audit = s4.json_read(metadata)
            assert audit['status'] == 'complete' and audit['mesh']['auto_size'] == 5
            assert audit['zeta'] == zeta
            validate_fields(cached, zeta)
            case.update(source_directory=str(cached), source_hashes={str(p): s4.digest(p) for p in
                (metadata, cached/'surface_fields.npz', cached/'edge_integrals_raw.csv', cached/'eigenfrequencies.csv')})
        cases.append(case)
    manifest = {'zetas': ZETAS, 'eta': .96, 'b0_nm': 245, 'a_nm': 820, 'H_nm': 200,
        'refractive_index': 3.3, 'mesh': 5, 'edge_order': 128, 'z_nm': 100, 'k_over_G': [0, 0],
        'user_mode': 'py', 'internal_mode': 'px', 'selection': 'Hz odd-x/even-y',
        'expected_eigenvalues': 4, 'shift': 'c_const/1.55[um]', 'new_solve_count': sum('source_directory' not in c for c in cases),
        'source_parameter_sha256': s4.digest(source), 'shared_parameter_sha256': s4.digest(s4.ROOT/'scripts/parameter.json'),
        'surface_side': 'slab_layer_domains_at_exact_z', 'cases': cases}
    s4.write_json(path, manifest)
    return manifest


def solve_case(client, band, config, case):
    directory = Path(case['directory'])
    done = directory / 'export.json'
    if done.exists() and s4.json_read(done)['status'] == 'complete':
        validate_fields(directory, case['zeta'])
        return
    params = s4.load_shared_parameters(case['snapshot'])
    hp = band.get_hole_params(params.unit_cell_2d.cell.to_fourier_params('surface_s_scan'))
    outer, holes, info = create_hexagon_design(band.A, hp)
    if min(i['min_dist'] for i in info) < band.D:
        raise ValueError('Geometry clearance failed')
    s4.write_json(directory/'99_config/realized_geometry.json', {'hole_params': hp.tolist(),
        'holes_um': np.asarray(holes).tolist(), 'outer_um': np.asarray(outer).tolist(), 'simulation_config': asdict(config)})
    runner = band.ReusableSimulationRun(config)
    model = runner.model
    audit = {'zeta': case['zeta'], 'status': 'running', 'new_eigensolves': 0}
    try:
        checkpoint = directory / 's4_gamma.mph'
        if checkpoint.exists():
            client.remove(model)
            model = client.load(str(checkpoint))
            runner.model = model
            runner.plane_datasets = {'center': 'cpl1', 'air': 'cpl2', 'yz': 'cpl3', 'xz': 'cpl4'}
        else:
            if (directory/'solve_attempt.json').exists():
                raise RuntimeError('Uncheckpointed attempt requires inspection before retry')
            s4.write_json(directory/'solve_attempt.json', {'zeta': case['zeta'], 'mesh': 5, 'Gamma': True})
            print(f"SOLVE zeta={case['zeta']}", flush=True)
            runner.build_and_run(.82, np.asarray(holes).tolist(), {'kx': 0., 'ky': 0.}, mesh_auto_size=5)
            model = runner.model
            model.save(str(checkpoint))
            audit['new_eigensolves'] = 1
        audit['mesh'] = mesh_identity(model)
        assert audit['mesh']['auto_size'] == 5
        assert int(model.java.study('std1').feature('eig').getInt('neigs')) == 2
        np.testing.assert_allclose([model.java.param().evaluate(k) for k in ('a', 'H', 'kx', 'ky')], [820e-9, 200e-9, 0, 0], atol=1e-18)
        modes = band.persist_k_point_solution(runner, directory/'01_results')
        for tag in ('s4interp', 's4surface_freq'):
            if tag in list(model.java.result().numerical().tags()):
                model.java.result().numerical().remove(tag)
        freq = frequencies(model)
        assert len(freq) == 4
        modes[['re', 'im', 'q']] = freq[['re', 'im', 'q']]
        modes.to_csv(directory/'eigenfrequencies.csv', index=False)
        indices = list(modes.index[(modes.p_weight > .9) & modes.is_valid])
        if len(indices) != 2:
            raise ValueError('Expected two p-dominant modes')
        sampler = s4.Sampler(SimpleNamespace(model=model))
        ref = np.load(REFERENCE/'surface_fields.npz')
        xy, area_xy, w = ref['center_probe_xy_m'], ref['area_xy_m'], ref['area_weights_m2']
        fields = [field_pair(sampler, indices, xy*r, 0., ['ewfd.Hz'])[:, 0] for r in ([1, 1], [-1, 1], [1, -1])]
        mode_audit = [{'mode_idx': int(idx), 'frequency_thz': float(modes.loc[idx, 're']), 'Q': float(modes.loc[idx, 'q']),
            **weighted_parity(*(h[:, j] for h in fields), np.ones(len(xy)))} for j, idx in enumerate(indices)]
        s4.write_csv(directory/'mode_audit.csv', mode_audit)
        good = [j for j, r in enumerate(mode_audit) if max(r['odd_x_error'], r['even_y_error']) < .05]
        if len(good) == 1:
            coeff = np.zeros(2, complex)
            coeff[good[0]] = 1
            audit.update(field_kind='single_eigenmode', mode_idx=int(indices[good[0]]), Q=mode_audit[good[0]]['Q'])
        else:
            if np.ptp(modes.loc[indices, 're']) / np.mean(modes.loc[indices, 're']) > 1e-3:
                raise ValueError('No unambiguous py mode; p pair is not near-degenerate')
            coeff = symmetry_coefficients(*fields, np.ones(len(xy)))
            audit['field_kind'] = 'constant_complex_p_subspace_combination_not_exact_single_eigenmode'
        audit['coefficients'] = [complex_columns({'mode_idx': int(i), 'coefficient': c}) for i, c in zip(indices, coeff)]
        audit['p_pair_frequency_thz'] = modes.loc[indices, 're'].tolist()
        audit['p_pair_Q'] = modes.loc[indices, 'q'].tolist()
        validation = [field_pair(sampler, indices, area_xy*r, 0., ['ewfd.Hz'])[:, 0] @ coeff for r in ([1, 1], [-1, 1], [1, -1])]
        audit['validation_parity'] = weighted_parity(*validation, w)
        if max(audit['validation_parity'].values()) > .06:
            raise ValueError('Independent py parity validation failed')
        surface = [field_pair(sampler, indices, area_xy*r, 100e-9, ['ewfd.Hz'])[:, 0] @ coeff for r in ([1, 1], [-1, 1], [1, -1])]
        audit['surface_parity'] = weighted_parity(*surface, w)
        if max(audit['surface_parity'].values()) > .1:
            raise ValueError('Surface py parity validation failed')
        xyz = sampler.plane(indices[0], area_xy[:5], 100e-9, ['x', 'y', 'z']).real
        np.testing.assert_allclose(xyz, np.column_stack([area_xy[:5], np.full(5, 100e-9)]), rtol=1e-9, atol=1e-14)
        edges = hole_edges(np.asarray(holes)*1e-6)
        rules = [line_rule(e, 128) for e in edges]
        points = np.concatenate([p for p, weights in rules])
        components = field_pair(sampler, indices, points, 100e-9, ['ewfd.Hz'])[:, 0].reshape(18, 128, 2)
        h = components @ coeff
        rows = []
        for j, e in enumerate(edges):
            integral = rules[j][1] @ h[j]
            rows.append({'order': 128, 'hole': e.hole, 'edge': e.edge, 'group': e.group,
                'length_m': e.length, 'nx': e.normal[0], 'ny': e.normal[1],
                'integral_hz': integral, 'qx': integral*e.normal[0], 'qy': integral*e.normal[1]})
        s4.write_csv(directory/'edge_integrals_raw.csv', rows)
        np.savez_compressed(directory/'surface_fields.npz', area_xy_m=area_xy, area_weights_m2=w,
            surface_Hz_raw=surface[0], center_validation_Hz_raw=validation[0], center_probe_xy_m=xy,
            center_probe_Hz_raw=fields[0] @ coeff, coefficients=coeff, mode_indices=indices,
            outer_m=np.asarray(outer)*1e-6, holes_m=np.asarray(holes)*1e-6, z_nm=100.,
            edge_xy_m_128=points.reshape(18, 128, 2), edge_weights_m_128=np.asarray([ww for p, ww in rules]),
            edge_Hz_raw_128=h, edge_Hz_components_128=components)
        validate_fields(directory, case['zeta'])
        audit['status'] = 'complete'
        s4.write_json(done, audit)
        print(f"EXPORTED zeta={case['zeta']} {audit['field_kind']}", flush=True)
    finally:
        client.remove(model)
        runner.model = None


def run(manifest):
    record = {'status': 'starting', 'completed_zetas': [], 'preflight': s4.runtime_preflight()}
    s4.write_json(ROOT/'execution.json', record)
    import mph
    from scripts.run_main import run_band_pair as band
    from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationConfig
    client = mph.start(version='6.3', cores=4)
    try:
        if not all(client.java.checkoutLicense(p) for p in ('COMSOL', 'WAVEOPTICS')):
            raise RuntimeError('Required license checkout failed')
        config = SimulationConfig(slab_height='200 [nm]', slab_refractive_index='3.3',
            eigenmode_count=2, mesh_auto_size=5, eigenfrequency_shift='c_const/1.55[um]')
        client.java.showProgress(str(ROOT/'comsol_progress.log'))
        for case in manifest['cases']:
            if 'source_directory' not in case:
                solve_case(client, band, config, case)
            else:
                print(f"REUSE zeta={case['zeta']}", flush=True)
            record['completed_zetas'].append(case['zeta'])
            s4.write_json(ROOT/'execution.json', record)
        from scripts.analysis.plot_boundary_signed_scan import analyze
        analyze(ROOT)
        assert s4.digest(s4.ROOT/'scripts/parameter.json') == manifest['shared_parameter_sha256']
        prepare()  # Check original cache hashes again.
        record.update(status='complete', analysis_completed=True, source_hashes_unchanged=True)
    except Exception as exc:
        record.update(status='failed', error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        client.clear()
        s4.write_json(ROOT/'execution.json', record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    manifest = prepare()
    print(f"Prepared {len(ZETAS)} points, {manifest['new_solve_count']} new solves, output={ROOT}", flush=True)
    if not args.prepare_only:
        run(manifest)
