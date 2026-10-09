"""Exact 19-k periodic reference fields for the S5 finite-mode validation.

Preparation is offline. --run reuses the saved Gamma mesh and solves only
missing non-Gamma points; --gamma-only exports the existing Gamma solution.
No shared parameters or source MPH files are modified.
"""
from __future__ import annotations

from comsol_workflow.dipolar_radiation import periodic_profiles, mode_assignment

import argparse
import json
import os
from pathlib import Path
import time
import traceback
from types import SimpleNamespace

import numpy as np
from scipy.optimize import linear_sum_assignment

from comsol_workflow.boundary_integrals import area_rule
from comsol_workflow.dipolar_lattice_analysis import load_data
from scripts.analysis.visualize_dipolar_broadening_2d import SOURCE, ROOT
from comsol_workflow.result_io import digest, read_json as json_read, write_json, write_csv

from comsol_workflow.dipolar_inputs import load_dipolar_inputs

INPUTS = load_dipolar_inputs()
S4 = INPUTS.boundary_case_dir
OUT = INPUTS.reference_dir






def prepare():
    d = load_data(SOURCE)
    geometry_path = S4 / '99_config/realized_geometry.json'
    geometry = json_read(geometry_path)
    model = S4 / '00_model/s4_gamma.mph'
    if geometry['simulation_config']['mesh_auto_size'] != 5:
        raise ValueError('Unexpected saved reference mesh')
    for sub in ['99_config', '01_results', '80_logs', '00_model']:
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    identity = {'version': 's5-periodic-reference-v1', 'source_model': str(model),
                'source_model_sha256': digest(model), 'source_Hz': str(SOURCE),
                'source_Hz_sha256': digest(SOURCE),
                'geometry_source': str(geometry_path), 'geometry_sha256': digest(geometry_path),
                'source_parameter_sha256': digest(S4 / '99_config/parameter.json'),
                'geometry': {'b0_nm': 245., 'eta': .96, 'zeta': 1.156},
                'mesh': 5, 'period_um': d['period'], 'expected_modes_per_point': 4,
                'phase_convention': 'native COMSOL exp(+i omega t); exp(-i k.r) Bloch phase',
                'normalization': 'one real Hz RMS at z=0 per eigenmode, applied to all fields',
                'k_um_inv': d['k19'].tolist(),
                'air_quadrature_order': 12, 'air_subdivisions': [2, 4],
                'air_heights_over_Hair': [.6, .75]}
    manifest = OUT / '99_config/manifest.json'
    if manifest.exists() and json_read(manifest) != identity:
        raise ValueError('Existing S5 reference identity differs; refusing cache reuse')
    write_json(manifest, identity)
    g = 4 * np.pi / np.sqrt(3) / d['period']
    points = [{'index': i, 'layer': 1 if i == 0 else 2 if i < 7 else 3,
               'kx_um_inv': float(k[0]), 'ky_um_inv': float(k[1]),
               'kx_over_G': float(k[0]/g), 'ky_over_G': float(k[1]/g),
               'action': 'reuse_Gamma' if i == 0 else 'solve_exact_k'}
              for i, k in enumerate(d['k19'])]
    write_csv(OUT / '99_config/k19_points.csv', points)
    return d, geometry, identity


def extract(model, d, geometry, index, modes, mesh, tracking_reference, tracking_index=0, stem=None):
    from scripts.run_main.run_boundary_analysis import Sampler
    j = model.java
    if 's4interp' in list(j.result().numerical().tags()):
        j.result().numerical().remove('s4interp')
    sampler = Sampler(SimpleNamespace(model=model))
    k = d['k19'][index]
    # Sampler uses the conjugated S4 convention; undo it consistently for S5.
    xy = d['rho_points'] * 1e-6
    hz = np.array([sampler.plane(m, xy, 0., ['ewfd.Hz'])[:, 0].conj() for m in range(4)])
    center_e = np.array([sampler.plane(m, xy, 0., ['ewfd.Ex', 'ewfd.Ey']).conj() for m in range(4)])
    norm = np.sqrt(np.mean(abs(hz)**2, axis=1))
    if not np.isfinite(norm).all() or np.any(norm <= 0):
        raise ValueError('Invalid center Hz normalization')
    u = periodic_profiles(hz/norm[:, None], d['rho_points'], k)
    if index == 0 and tracking_reference is None:
        assignment, overlap = np.arange(4), np.eye(4)
    else:
        assignment, overlap = mode_assignment(tracking_reference, u)
    target = int(assignment[1])
    target_overlap = float(overlap[1, target])
    alternative = float(np.max(np.delete(overlap[1], target)))
    tracking = json_read(S4/'99_config/parameter.json')['unit_cell_2d']
    minimum_overlap = float(tracking['minimum_overlap'])
    minimum_gap = float(tracking['minimum_ambiguity_gap'])
    if target_overlap < minimum_overlap or target_overlap - alternative < minimum_gap:
        raise ValueError(f'Ambiguous py tracking at index {index}: {overlap[1].tolist()}')
    outer = np.asarray(geometry['outer']) * 1e-6
    air_height = float(j.param().evaluate('H_air'))
    heights = np.array([.6, .75]) * air_height
    radiations = np.empty((2, 2, 4, 2), complex)
    sample_counts = []
    for refinement, subdiv in enumerate([2, 4]):
        air_xy, weights, _ = area_rule(outer, [], 12, subdiv)
        sample_counts.append(len(weights))
        phase = np.exp(1j*(air_xy/1e-6 @ k))
        for iz, z in enumerate(heights):
            for m in range(4):
                values = sampler.plane(m, air_xy, z, ['ewfd.Ex', 'ewfd.Ey']).conj()
                radiations[refinement, iz, m] = (weights*phase) @ values / weights.sum() / norm[m]
    omega = 2*np.pi*(modes.re.to_numpy()+1j*modes.im.to_numpy())*1e12
    kz = np.sqrt((omega/299792458.)**2-np.dot(k,k)*1e12)
    backpropagated = radiations * np.exp(1j*heights[None,:,None,None]*kz[None,None,:,None])
    stem = stem or f'k{index:02d}'
    filename = OUT / f'01_results/{stem}_reference.npz'
    np.savez_compressed(filename, k_um_inv=k, rho_points_um=d['rho_points'],
                        Hz_native=hz, E_xy_center_native=center_e, Hz_rms=norm,
                        u_Hz=u, air_z_m=heights, air_mean_E_over_Hz_rms=radiations,
                        air_mean_E_at_z0_over_Hz_rms=backpropagated,
                        frequencies_THz=modes[['re','im','q']].to_numpy(),
                        gamma_band_to_mode=assignment, overlap_to_tracking_reference=overlap)
    scale = float(np.max(abs(backpropagated[:,:,target])))
    audit = {'status': 'complete', 'index': index, 'k_um_inv': k.tolist(),
             'source_model_sha256': digest(S4/'00_model/s4_gamma.mph'),
             'mesh': mesh, 'py_mode_idx': target, 'tracking_reference_index': tracking_index,
             'py_overlap_to_tracking_reference': target_overlap,
             'py_overlap_gap': target_overlap-alternative,
             'tracking_minimum_overlap': minimum_overlap, 'tracking_minimum_gap': minimum_gap,
             'py_frequency_thz': float(modes.iloc[target].re), 'py_Q': float(modes.iloc[target].q),
             'quadrature_samples': sample_counts,
             'quadrature_change_over_max_c': float(np.max(abs(backpropagated[1,:,target]-backpropagated[0,:,target]))/max(scale,1e-300)),
             'height_change_over_max_c': float(np.max(abs(backpropagated[1,1,target]-backpropagated[1,0,target]))/max(scale,1e-300)),
             'data_sha256': digest(filename),
             'scope': 'periodic reference only; no finite-boundary radiation response inferred'}
    write_json(OUT/f'01_results/{stem}_audit.json', audit)
    print(f'k{index:02d}: saved; py mode={target}; overlap={target_overlap:.6f}; Q={audit["py_Q"]:.7g}', flush=True)
    return u[assignment]


def run(gamma_only=False, gamma_mesh=None):
    from scripts.run_main.run_boundary_analysis import runtime_preflight
    d, geometry, identity = prepare()
    preflight = runtime_preflight()
    import mph
    from jpype import JClass
    from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies, mesh_identity
    write_json(OUT/'99_config/preflight.json', preflight)
    client = mph.start(version='6.3', cores=4)
    model = None
    state = {'gamma_mesh': gamma_mesh, 'status': 'running', 'pid': os.getpid(), 'new_solves': 0, 'completed_indices': [],
             'start_time_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    execution_path = OUT / ('99_config/execution.json' if gamma_mesh is None else f'99_config/gamma_mesh{gamma_mesh}_execution.json')
    write_json(execution_path, state)
    try:
        util = JClass('com.comsol.model.util.ModelUtil')
        if not util.checkoutLicense('COMSOL','WAVEOPTICS'):
            raise RuntimeError('COMSOL and WAVEOPTICS license checkout failed')
        util.showProgress(str(OUT/'00_model/comsol_progress.log'))
        print('Loading saved Gamma model; no source writes.',flush=True)
        model = client.load(identity['source_model'])
        j = model.java
        np.testing.assert_allclose([float(j.param().evaluate('a')),float(j.param().evaluate('H'))],
                                   [d['period']*1e-6,200e-9],rtol=1e-12,atol=0)
        mesh = mesh_identity(model)
        if mesh['auto_size'] != 5 or int(j.study('std1').feature('eig').getInt('neigs')) != 2:
            raise ValueError('Unexpected mesh or eigenmode count')
        if any(abs(float(j.param().evaluate(k))) > 1e-12 for k in ['kx','ky']):
            raise ValueError('Source model is not a saved Gamma solution')
        print('Model identity checked: '+json.dumps(mesh),flush=True)
        history = {}
        if gamma_mesh is not None:
            with np.load(OUT/'01_results/k00_reference.npz',allow_pickle=False) as anchor:
                history[0] = anchor['u_Hz'][anchor['gamma_band_to_mode']]
            print(f'Gamma mesh control: rebuilding only the reference mesh at size {gamma_mesh}',flush=True)
            j.component('comp1').mesh('mesh1').autoMeshSize(gamma_mesh)
            j.component('comp1').mesh('mesh1').run()
            j.sol('sol1').clearSolutionData()
            j.sol('sol1').runAll()
            state['new_solves'] += 1
            mesh = mesh_identity(model)
        for index in range(1 if gamma_only or gamma_mesh is not None else 19):
            stem = f'k{index:02d}' if gamma_mesh is None else f'gamma_mesh{gamma_mesh}'
            data_path = OUT/f'01_results/{stem}_reference.npz'
            audit_path = OUT/f'01_results/{stem}_audit.json'
            if audit_path.exists():
                audit = json_read(audit_path)
                if audit['status'] != 'complete' or digest(data_path) != audit['data_sha256']:
                    raise ValueError('Invalid completed point cache')
                with np.load(data_path, allow_pickle=False) as cached:
                    history[index] = cached['u_Hz'][cached['gamma_band_to_mode']]
                print(f'k{index:02d}: reuse complete cache',flush=True)
            else:
                if index:
                    k = d['k19'][index]
                    for name, value in zip(['kx','ky'], k):
                        j.param().set(name,f'{value:.17g}[1/um]')
                    print(f'k{index:02d}: eigensolve at {k.tolist()} um^-1',flush=True)
                    j.sol('sol1').clearSolutionData()
                    j.sol('sol1').runAll()
                    state['new_solves'] += 1
                    if mesh_identity(model) != mesh:
                        raise ValueError('Mesh changed during a fixed-mesh solve')
                if 's4surface_freq' in list(j.result().numerical().tags()):
                    j.result().numerical().remove('s4surface_freq')
                modes = frequencies(model)
                if len(modes) != 4:
                    raise ValueError(f'Expected four modes, got {len(modes)}')
                if index == 0 and gamma_mesh is None and abs(float(modes.iloc[1].re)-198.38683823113476)>1e-7:
                    raise ValueError('Saved Gamma frequency does not match source identity')
                previous = min(history,key=lambda n: np.linalg.norm(d['k19'][n]-d['k19'][index])) if history else 0
                history[index] = extract(model,d,geometry,index,modes,mesh,history.get(previous),previous,stem)
            state['completed_indices'].append(index)
            write_json(execution_path,state)
        state.update(status='complete', source_model_unchanged=digest(S4/'00_model/s4_gamma.mph')==identity['source_model_sha256'])
        if not state['source_model_unchanged']:
            raise RuntimeError('Source model hash changed')
    except Exception as exc:
        state.update(status='failed', error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        write_json(execution_path,state)
        if model is not None:
            client.remove(model)
        client.clear()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--gamma-only', action='store_true')
    parser.add_argument('--gamma-mesh', type=int, choices=[3,4])
    args = parser.parse_args()
    if args.run:
        run(args.gamma_only,args.gamma_mesh)
    else:
        _,_,manifest = prepare()
        print(json.dumps(manifest,indent=2))
