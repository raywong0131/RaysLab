"""Approved four Gamma solves and five z=100 nm exports; original models stay intact."""
from __future__ import annotations
import argparse
import json
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.boundary_integrals import area_rule, complex_columns, group_totals, hole_edges, line_rule, parity_errors, phase_and_norm
from scripts.run_main import run_boundary_analysis as s4
from scripts.analysis.plot_boundary_cached_zeta import plot

ROOT = s4.ROOT / 'scripts/.out/unit_cell_2D/S4_surface_eta0.96_zeta_scan_20260918'


def frequencies(model):
    node = model.java.result().numerical().create('s4surface_freq', 'EvalGlobal')
    node.set('data', 'dset1')
    node.set('expr', ['ewfd.omega/(2*pi*1[THz])', 'ewfd.damp/(2*pi*1[THz])', 'ewfd.Qfactor'])
    node.set('unit', ['1', '1', '1'])
    result = np.asarray(node.getReal()).T[:, :3]
    # COMSOL reports the magnitude of Q; retain the signed damping in the returned table.
    np.testing.assert_allclose(np.abs(result[:, 0]/(2*result[:, 1])), result[:, 2], rtol=1e-7)
    return pd.DataFrame(result, columns=['re', 'im', 'q'])


def mesh_identity(model):
    mesh = model.java.component('comp1').mesh('mesh1')
    slab = list(model.java.component('comp1').selection('sel_slab_layer').entities(3))
    return {'auto_size': int(mesh.autoMeshSize()), 'vertices': int(mesh.getNumVertex()),
            'tetrahedra': int(mesh.getNumElem('tet')),
            'slab_tetrahedra': int(np.isin(np.asarray(mesh.getElemEntity('tet'), int), slab).sum())}


def choose_mode(sampler, modes, source_series, a, directory):
    point = source_series/'01_results/k_points/kx=0.000000_ky=0.000000'
    previous = pd.read_csv(point/'eigenfrequencies.csv')
    old_idx = int((previous.p_weight*previous.px_weight).idxmax())
    old = pd.read_parquet(point/f'eigenmodes/{old_idx:02d}_Hz_center.parquet')
    old = old.groupby(['x', 'y'], as_index=False)[['re', 'im']].mean()
    axis = np.linspace(-.3*a, .3*a, 17)
    xx, yy = np.meshgrid(axis, axis)
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    # Sampler exports the conjugated exp(-i wt) phasor; match that convention.
    ref = LinearNDInterpolator(old[['x','y']].to_numpy(), old.re.to_numpy()-1j*old.im.to_numpy())(xy/1e-6)
    rows, fields = [], {}
    for idx, m in modes.iterrows():
        h = sampler.plane(idx, xy, 0., ['ewfd.Hz'])[:, 0]
        hx = sampler.plane(idx, xy*[-1,1], 0., ['ewfd.Hz'])[:, 0]
        hy = sampler.plane(idx, xy*[1,-1], 0., ['ewfd.Hz'])[:, 0]
        overlap = abs(np.vdot(ref,h))/(np.linalg.norm(ref)*np.linalg.norm(h))
        errors = parity_errors(h,hx,hy)
        rows.append({'mode_idx':idx,'frequency_thz':m.re,'Q':m.q,
                     'overlap_cached_center':overlap,**errors})
        fields[idx] = h
    eligible = [r for r in rows if r['overlap_cached_center']>.95
                and max(r['odd_x_error'],r['even_y_error'])<.1
                and abs(r['frequency_thz']-previous.loc[old_idx,'re'])<.1]
    s4.write_csv(directory/'mode_audit.csv',rows)
    if len(eligible)!=1:
        raise ValueError('Target py mode ambiguous; inspect mode_audit.csv')
    selected = eligible[0]
    return selected, xy, fields[selected['mode_idx']]


def export_case(client, case, preparation):
    zeta = case['zeta']
    directory = ROOT/f'zeta{zeta:g}'
    directory.mkdir(exist_ok=True)
    done = directory/'export.json'
    if done.exists() and s4.json_read(done).get('status')=='complete':
        print(f'zeta={zeta}: reuse completed export',flush=True)
        return s4.json_read(done)
    original = Path(preparation['separate_Gamma_model']) if zeta==1.156 else Path(case['source_model'])
    digest = s4.digest(original)
    if zeta!=1.156 and digest!=case['source_mph_sha256']:
        raise ValueError('Source MPH identity changed')
    checkpoint = directory/'s4_gamma.mph'
    source = checkpoint if checkpoint.exists() else original
    audit={'zeta':zeta,'source':str(original),'source_sha256':digest,'new_eigensolves':0,
           'status':'loading','surface_z_nm':100,'surface_side':'slab_layer_domains_at_exact_z'}
    model = None
    start = time.monotonic()
    try:
        print(f'zeta={zeta}: loading {source.name}',flush=True)
        model = client.load(str(source))
        j = model.java
        a = float(j.param().evaluate('a'))
        height = float(j.param().evaluate('H'))
        np.testing.assert_allclose([a,height],[820e-9,200e-9],rtol=1e-12,atol=0)
        mesh_before=mesh_identity(model)
        audit['mesh']=mesh_before
        neigs = int(j.study('std1').feature('eig').getInt('neigs'))
        if neigs!=2 or mesh_before['auto_size']!=5:
            raise ValueError(f'Unexpected solver/mesh identity: {neigs}, {mesh_before}')
        audit['study_neigs']=neigs
        audit['expected_eigenvalue_count']=4
        audit['shift']=str(j.study('std1').feature('eig').getString('shift'))
        print(f'zeta={zeta}: mesh={mesh_before}; neigs={neigs}; shift={audit["shift"]}',flush=True)
        if zeta!=1.156 and not checkpoint.exists():
            if (directory/'solve_attempt.json').exists():
                raise RuntimeError('Prior uncheckpointed solve attempt; inspect before retry')
            j.param().set('kx','0*G')
            j.param().set('ky','0*G')
            s4.write_json(directory/'solve_attempt.json',{'started':True,'source_sha256':digest})
            print(f'zeta={zeta}: Gamma eigensolve started (existing mesh)',flush=True)
            solution=j.sol('sol1')
            solution.clearSolutionData()
            solution.runAll()
            audit['new_eigensolves']=1
            if mesh_identity(model)!=mesh_before:
                raise ValueError('Mesh changed during eigensolve')
            model.save(str(checkpoint))
            print(f'zeta={zeta}: solved; separate Gamma checkpoint saved',flush=True)
        if any(abs(float(j.param().evaluate(k)))>1e-12 for k in ('kx','ky')):
            raise ValueError('Loaded solution is not Gamma')
        for tag in ('s4interp','s4surface_freq'):
            if tag in list(j.result().numerical().tags()):
                j.result().numerical().remove(tag)
        modes=frequencies(model)
        if len(modes)!=4:
            raise ValueError(f'Expected four cached-mode counterparts, got {len(modes)}')
        modes.to_csv(directory/'eigenfrequencies.csv',index=False)
        sampler=s4.Sampler(SimpleNamespace(model=model))
        source_series=Path(case['source_model']).parents[1]
        selected,center_xy,center_h=choose_mode(sampler,modes,source_series,a,directory)
        mode=selected['mode_idx']
        audit.update(selected)
        outer,holes,_=create_hexagon_design(case['a_um'],case['hole_params'])
        outer=np.asarray(outer)*1e-6
        edges=hole_edges(np.asarray(holes)*1e-6)
        xy,weights,_=area_rule(outer,[],12,4)
        surface= sampler.plane(mode,xy,100e-9,['ewfd.Hz'])[:,0]
        # Explicit coordinate check guards against display-unit or cut-plane mistakes.
        xyz=sampler.plane(mode,xy[:5],100e-9,['x','y','z']).real
        np.testing.assert_allclose(xyz,np.column_stack([xy[:5],np.full(5,100e-9)]),rtol=1e-9,atol=1e-14)
        surface_parity=parity_errors(surface,
            sampler.plane(mode,xy*[-1,1],100e-9,['ewfd.Hz'])[:,0],
            sampler.plane(mode,xy*[1,-1],100e-9,['ewfd.Hz'])[:,0])
        audit['surface_parity']=surface_parity
        arrays={'area_xy_m':xy,'area_weights_m2':weights,'surface_Hz_raw':surface,
                'center_probe_xy_m':center_xy,'center_probe_Hz_raw':center_h,
                'outer_m':outer,'holes_m':np.asarray(holes)*1e-6,'z_nm':100.,'mode_idx':mode}
        raw_rows=[]
        for order in (128,256):
            rules=[line_rule(e,order) for e in edges]
            points=np.concatenate([p for p,w in rules])
            h=sampler.plane(mode,points,100e-9,['ewfd.Hz'])[:,0].reshape(18,order)
            arrays[f'edge_xy_m_{order}']=points.reshape(18,order,2)
            arrays[f'edge_weights_m_{order}']=np.asarray([w for p,w in rules])
            arrays[f'edge_Hz_raw_{order}']=h
            for i,e in enumerate(edges):
                integral=np.dot(rules[i][1],h[i])
                raw_rows.append({'order':order,'hole':e.hole,'edge':e.edge,'group':e.group,
                    'length_m':e.length,'nx':e.normal[0],'ny':e.normal[1],
                    'integral_hz':integral,'qx':integral*e.normal[0],'qy':integral*e.normal[1]})
        np.savez_compressed(directory/'surface_fields.npz',**arrays)
        s4.write_csv(directory/'edge_integrals_raw.csv',raw_rows)
        audit.update(status='complete',elapsed_seconds=time.monotonic()-start,
                     source_unchanged=s4.digest(original)==digest)
        if not audit['source_unchanged']:
            raise RuntimeError('Source MPH changed')
        s4.write_json(done,audit)
        print(f'zeta={zeta}: exported mode {mode}, Q={selected["Q"]:.8g}, overlap={selected["overlap_cached_center"]:.6f}',flush=True)
        return audit
    except Exception as exc:
        audit.update(status='failed',error=str(exc),traceback=traceback.format_exc())
        s4.write_json(directory/'export_failure.json',audit)
        raise
    finally:
        if model is not None:
            client.remove(model)


def analyze():
    preparation=s4.json_read(ROOT/'99_config/preparation.json')
    refdata=np.load(ROOT/'zeta1.156/surface_fields.npz')
    weights=refdata['area_weights_m2']
    reference=refdata['surface_Hz_raw']
    reference=reference*phase_and_norm(reference,weights)
    totals=[]
    for case in preparation['cases']:
        zeta=case['zeta']
        directory=ROOT/f'zeta{zeta:g}'
        data=np.load(directory/'surface_fields.npz')
        np.testing.assert_allclose(data['area_xy_m'],refdata['area_xy_m'],rtol=0,atol=1e-20)
        h=data['surface_Hz_raw']
        factor=phase_and_norm(h,weights,reference)/(case['a_um']*1e-6)
        overlap=abs(np.dot(weights,reference.conj()*h))/np.sqrt(np.dot(weights,abs(reference)**2)*np.dot(weights,abs(h)**2))
        if overlap<.9:
            raise ValueError('Surface mode overlap below threshold')
        audit=s4.json_read(directory/'export.json')
        raw=pd.read_csv(directory/'edge_integrals_raw.csv')
        for order in (128,256):
            rows=[]
            for r in raw[raw.order==order].to_dict('records'):
                rows.append({**r,'qx':complex(r['qx_re'],r['qx_im']),'qy':complex(r['qy_re'],r['qy_im'])})
            total=group_totals(rows)
            totals.append(complex_columns({'zeta':zeta,'eta':.96,'z_nm':100,'method':f'gauss{order}',
                'Q':audit['Q'],'frequency_thz':audit['frequency_thz'],'mode_idx':audit['mode_idx'],
                'Hz_surface_rms':np.sqrt(np.dot(weights,abs(h)**2)/sum(weights)),
                'phase_rotation_rad':np.angle(factor),'surface_overlap_reference':overlap,
                'f_A_raw':total['QAx'],'f_B_raw':total['QBx'],'f_raw':total['Qx'],
                'F_A':total['QAx']*factor,'F_B':total['QBx']*factor,'F':total['Qx']*factor,
                'cancellation_ratio':total['cancellation_ratio']}))
    table=pd.DataFrame(totals).sort_values(['zeta','method'])
    table.to_csv(ROOT/'boundary_sums.csv',index=False)
    plot(table,ROOT,methods=('gauss256','gauss128'),method_labels=('Gauss 256','Gauss 128'),
         section_label=r'$z=100$ nm, slab-side $H_z$',q_title=r'd  $Q$ at $\Gamma$')
    print(table[['zeta','method','F_re','F_im','cancellation_ratio','Q']].to_string(index=False),flush=True)
    return totals


def run():
    preparation=s4.json_read(ROOT/'99_config/preparation.json')
    record={'status':'starting','authorization':'User explicitly approved COMSOL and required exports',
            'planned_new_eigensolves':4,'z_nm':100,'phase_convention':'exp(-i omega t)',
            'normalization':'F=f/(a Hz_surface_RMS); common area grid; surface-Hz overlap reference zeta1.156',
            'edge_quadrature_orders':[128,256],'full_k_scan':False,'cases':[]}
    record['preflight']=s4.runtime_preflight()
    s4.write_json(ROOT/'execution.json',record)
    import mph
    client=mph.start(version='6.3',cores=4)
    try:
        # Export the existing reference before starting any new solve.
        for case in sorted(preparation['cases'],key=lambda c:(c['zeta']!=1.156,c['zeta'])):
            record['cases'].append(export_case(client,case,preparation))
            s4.write_json(ROOT/'execution.json',record)
        analyze()
        record['status']='complete'
        record['source_models_unchanged']=all(c['source_unchanged'] for c in record['cases'])
        record['new_eigensolves_this_invocation']=sum(c['new_eigensolves'] for c in record['cases'])
    except Exception as exc:
        record.update(status='failed',error=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        s4.write_json(ROOT/'execution.json',record)
        client.clear()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--run-approved',action='store_true')
    action.add_argument('--plot-only',action='store_true')
    args=parser.parse_args()
    if args.plot_only:
        analyze()
    else:
        run()
