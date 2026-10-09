"""Bounded px Gamma remeshing test; raw eigenmode Q, no symmetry mixing."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
from types import SimpleNamespace
import time
import traceback
import numpy as np
import pandas as pd
from scripts.run_main import run_boundary_analysis as s4
from scripts.run_sweep.run_boundary_mesh_convergence import explicit_frequencies, mesh_stats, snapshot_for_mesh, subspace_scores, scale_mesh_sizes
from scripts.run_sweep.run_boundary_dual_validation import assign_modes
from scripts.analysis.audit_s31_boundary import native_integral

ARCHIVE=s4.ROOT/'results/S4_boundary_analysis'
GROUP='10_px_mesh_test'
OUTPUT=s4.ROOT/'scripts/.out/unit_cell_2D/S4_px_b245_eta0.96_zeta0.8385-0.84_mesh4-3-2_20260929'
ZETAS=(.8385,.8393,.84)
LEVELS=(4,3,2)


def prepare(size_factors=None,zetas=ZETAS):
    controls=[(5,f) for f in size_factors] if size_factors else [(level,1.) for level in LEVELS]
    expected=sorted((z,level,factor) for z in zetas for level,factor in controls)
    config_path=OUTPUT/'99_config/config.json'
    if config_path.exists():
        config=s4.json_read(config_path)
        actual=sorted((c['zeta'],c['mesh'],c.get('size_factor',1.)) for c in config['cases'])
        if actual!=expected:raise ValueError('Existing config identity differs from requested cases')
        for p,h in config['source_sha256'].items():
            if s4.digest(p)!=h:raise ValueError(f'Source changed: {p}')
        return config
    df=pd.read_csv(ARCHIVE/'80_logs/09_px_py_validation/scan_results.csv')
    sources={};cases=[]
    for z in zetas:
        row=df[(df.user_mode=='px')&(df.zeta==z)].iloc[0]
        old_config=ARCHIVE/f'99_config/09_px_py_validation/zeta{z:g}'
        snapshot=old_config/'parameter.json'
        realized=old_config/'realized_geometry.json'
        fields=Path(row.field_file);model=Path(row.source_model)
        for p in (snapshot,realized,fields,model):sources[str(p)]=s4.digest(p)
        params=s4.json_read(snapshot)
        assert params['mesh_size']==5
        assert [params['unit_cell_2d'][k] for k in ('b0_nm','eta','zeta')]==[245.,.96,z]
        for level,factor in controls:
            folder=OUTPUT/(f'zeta{z:g}_h{factor:g}' if size_factors else f'zeta{z:g}_mesh{level}')
            p=folder/'parameter.json';s4.write_json(p,snapshot_for_mesh(params,level));s4.load_shared_parameters(p)
            cases.append(dict(zeta=z,mesh=level,size_factor=factor,directory=str(folder),snapshot=str(p),snapshot_sha256=s4.digest(p),source_model=str(model),source_fields=str(fields)))
    sources[str(s4.ROOT/'scripts/parameter.json')]=s4.digest(s4.ROOT/'scripts/parameter.json')
    config=dict(eta=.96,b0_nm=245,a_nm=820,H_nm=200,kx=0,ky=0,user_mode='px',internal_mode='py',expected_eigenvalues=4,shift='c_const/1.55[um]',cores=4,cases=cases,source_sha256=sources,Q_observation_target=1e8)
    s4.write_json(config_path,config)
    baseline=df[(df.user_mode=='px')&df.zeta.isin(zetas)].copy();baseline['mesh']=5;baseline['size_factor']=1.
    baseline.to_csv(OUTPUT/'99_config/baseline.csv',index=False)
    return config


def extract(model,case):
    data=np.load(case['source_fields']);freq=explicit_frequencies(model)
    assert len(freq)==4 and np.isfinite(freq).all()
    tags=list(model.java.result().numerical().tags())
    if 's4interp' in tags:model.java.result().numerical().remove('s4interp')
    sampler=s4.Sampler(SimpleNamespace(model=model))
    probe=data['center_probe_xy_m']
    fields=[np.stack([sampler.plane(i,probe*r,0.,['ewfd.Hz'])[:,0] for r in ([1,1],[-1,1],[1,-1])]) for i in range(len(freq))]
    candidates=np.column_stack([f[0] for f in fields])
    basis=np.column_stack([data['px_probe_Hz'],data['py_probe_Hz']])
    scores=subspace_scores(basis,candidates,np.ones(len(probe)))
    indices=np.argsort(scores)[-2:]
    if min(scores[indices])<.8 or min(scores[indices])-max(np.delete(scores,indices))<.3:raise ValueError(f'Ambiguous p subspace: {scores}')
    assigned,costs=assign_modes([fields[i] for i in indices]);slot=assigned[0];mode=int(indices[slot])
    h=fields[mode][0];ref=data['px_probe_Hz']
    overlap=float(abs(np.vdot(ref,h))/(np.linalg.norm(ref)*np.linalg.norm(h)))
    if overlap<.9:raise ValueError(f'px identity uncertain: overlap={overlap}')
    xy,w=data['area_xy_m'],data['area_weights_m2']
    surface=sampler.plane(mode,xy,100e-9,['ewfd.Hz'])[:,0]
    factor=s4.phase_and_norm(surface,w,data['px_surface_Hz'])
    air=[];omega=2*np.pi*1e12*complex(freq[mode,0],-freq[mode,1])
    for j,z in enumerate(data['air_z_m']):
        tag=f'pxmesh_air{j}';dataset=model.java.result().dataset().create(tag,'CutPlane');dataset.set('data','dset1');dataset.set('quickplane','xy');dataset.set('quickz',f'{z:.15g}[m]')
        v=native_integral(model,'pxmesh_int','IntSurface',None,['1','ewfd.Ex','ewfd.Ey'],['m^2','V*m','V*m'],mode,8,tag)
        air.append(factor*v[1:]/v[0]*np.exp(-1j*omega/299792458*(z-100e-9)))
    folder=Path(case['directory']);np.savez_compressed(folder/'fields.npz',center_probe_xy_m=probe,px_probe_Hz=h,px_surface_Hz=surface,area_xy_m=xy,area_weights_m2=w,factor=factor,air=np.asarray(air),candidate_probe_Hz=candidates)
    s4.write_csv(folder/'candidates.csv',[dict(mode_idx=i,frequency_THz=float(f[0]),damping_THz=float(f[1]),Q=float(f[2]),p_projection=float(scores[i]),selected_px=i==mode) for i,f in enumerate(freq)])
    return dict(zeta=case['zeta'],mesh=case['mesh'],Q=float(freq[mode,2]),frequency_thz=float(freq[mode,0]),damping_thz=float(freq[mode,1]),parity_error=float(costs[slot,0]),overlap_mesh5=overlap,mode_idx=mode,air_Ex=air[0][0],air_Ey=air[0][1],air2_Ex=air[1][0],air2_Ey=air[1][1],air_Ex_abs=float(abs(air[0][0])),air_Ey_abs=float(abs(air[0][1])))


def publish(config,rows,status):
    task=OUTPUT if OUTPUT.name.startswith('S4_px_') else OUTPUT.parent
    rows=[]
    for path in task.rglob('summary.json'):
        saved=s4.json_read(path)
        if 'result' in saved:
            row=dict(saved['result']);row.setdefault('size_factor',1.);rows.append(row)
    planned=sum(len(c['cases']) for p in task.rglob('config.json') if not (c:=s4.json_read(p)).get('superseded',False))
    baseline=pd.read_csv(ARCHIVE/'80_logs/09_px_py_validation/scan_results.csv')
    baseline=baseline[(baseline.user_mode=='px')&baseline.zeta.isin([r['zeta'] for r in rows])].copy()
    baseline['mesh']=5;baseline['size_factor']=1.
    baseline=baseline[['zeta','mesh','size_factor','Q','frequency_thz','parity_error']]
    table=pd.concat([baseline,pd.DataFrame(rows)],ignore_index=True).sort_values(['size_factor','mesh','zeta'],ascending=[False,False,True])
    path=ARCHIVE/'80_logs'/GROUP/'mesh_results.csv';path.parent.mkdir(parents=True,exist_ok=True);table.to_csv(path,index=False)
    best=table.loc[table.Q.idxmax()]
    summary=dict(status=status,completed_new_cases=len(rows),planned_new_cases=planned,best_sample=dict(zeta=float(best.zeta),mesh=int(best.mesh),size_factor=float(best.size_factor),Q=float(best.Q)),target_reached=bool(best.Q>=1e8),convergence_established=False,raw_eigenmode_Q=True,output=str(task))
    s4.write_json(ARCHIVE/'99_config'/GROUP/'summary.json',summary)
    report=ARCHIVE/'12_reports'/GROUP/'README.md';report.parent.mkdir(parents=True,exist_ok=True)
    lines=['# px 奇点附近网格加密测试','',f'状态：{status}；新求解 {len(rows)}/{planned}。固定η=0.96、b0=245 nm、a=820 nm、H=200 nm、Γ点。mesh数字越小越细。','', 'Q来自原始复本征频率；不混合本征场，不将拟合峰值作为结果。10^8是观察目标，不是保证值。','', '| ζ | mesh基础等级 | h/h0 | Q | Hz对称误差 |','|---|---|---|---|---|']
    lines += [f'| {r.zeta:g} | {int(r.mesh)} | {r.size_factor:g} | {r.Q:.6g} | {r.parity_error:.5g} |' for r in table.itertuples()]
    lines += ['', '自动mesh等级测试与真实尺寸缩放分开比较：h/h0=0.5、0.25明确缩小了全局及板层尺寸约束；基础mesh等级仍为5，但实际网格为自定义加密。', '空气复振幅按各ζ原mesh5场作相位参考；跨ζ对比默认使用模长。Q不依赖相位。', '',f'当前最高Q={best.Q:.6g}（mesh={int(best.mesh)}, h/h0={best.size_factor:g}, ζ={best.zeta:g}）。收敛性需结合相邻网格与参数共同判断。','',f'模型、日志、参数快照：{task}']
    report.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def run(config,allow_concurrent):
    preflight=s4.runtime_preflight(allow_concurrent=allow_concurrent)
    record=dict(status='starting',pid=os.getpid(),preflight=preflight,cores=4,completed=[],source_unchanged=None)
    record_path=OUTPUT/'99_config/execution.json';s4.write_json(record_path,record)
    client=model=None;rows=[]
    try:
        import mph
        client=mph.start(version='6.3',cores=4)
        if not all(client.java.checkoutLicense(p) for p in ('COMSOL','WAVEOPTICS')):raise RuntimeError('COMSOL/WAVEOPTICS license unavailable')
        record['license_checkout']=True
        for case in sorted(config['cases'],key=lambda c:(-c.get('size_factor',1.),-c['mesh'],abs(c['zeta']-(min(v['zeta'] for v in config['cases']) if OUTPUT.name.startswith('local_') else .8393)))):
            folder=Path(case['directory']);done=folder/'summary.json'
            if done.exists():
                saved=s4.json_read(done)
                for p,h in saved['output_hashes'].items():
                    if s4.digest(folder/p)!=h:raise ValueError('Cached output changed')
                rows.append(saved['result']);continue
            if s4.digest(case['snapshot'])!=case['snapshot_sha256']:raise ValueError('Snapshot changed')
            checkpoint=folder/'s4_gamma.mph'
            if (folder/'attempt.json').exists() and not checkpoint.exists():raise RuntimeError(f'Interrupted solve requires review: {folder}')
            record.update(status='loading',current_case=dict(zeta=case['zeta'],mesh=case['mesh'],size_factor=case.get('size_factor',1.)));s4.write_json(record_path,record)
            model=client.load(str(checkpoint if checkpoint.exists() else case['source_model']))
            np.testing.assert_allclose([model.java.param().evaluate(k) for k in ('a','H','kx','ky')],[820e-9,200e-9,0,0],rtol=1e-12,atol=1e-18)
            client.java.showProgress(str(folder/'comsol_progress.log'))
            mesh=model.java.component('comp1').mesh('mesh1')
            if not checkpoint.exists():
                assert mesh.isAutomatic() and int(mesh.autoMeshSize())==5
                before=mesh_stats(model);s4.write_json(folder/'attempt.json',dict(started=time.time(),source_model=case['source_model'],mesh=case['mesh']))
                model.java.sol('sol1').clearSolutionData()
                if case.get('size_factor',1.)<1:
                    changes=scale_mesh_sizes(mesh,case['size_factor']);mesh.clearMesh();s4.write_json(folder/'size_changes.json',changes)
                else:mesh.autoMeshSize(case['mesh'])
                record['status']='meshing';s4.write_json(record_path,record);print('MESH',case['mesh'],case['zeta'],flush=True)
                mesh.run();stats=mesh_stats(model);assert stats['elements']['tet']>before['elements']['tet']
                if case.get('size_factor',1.)<1:assert stats['slab_layer_tetrahedra']>before['slab_layer_tetrahedra']*1.25
                s4.write_json(folder/'mesh.json',stats)
                record['status']='solving';s4.write_json(record_path,record);print('SOLVE',case['mesh'],case['zeta'],stats,flush=True)
                started=time.monotonic();model.java.sol('sol1').runAll();record['last_solve_seconds']=time.monotonic()-started
                model.save(str(checkpoint))
            assert int(mesh.autoMeshSize())==case['mesh']
            record['status']='extracting';s4.write_json(record_path,record)
            result=s4.complex_columns(extract(model,case));result['tetrahedra']=mesh_stats(model)['elements']['tet'];result['slab_tetrahedra']=mesh_stats(model)['slab_layer_tetrahedra'];result['size_factor']=case.get('size_factor',1.)
            s4.write_json(done,dict(result=result,output_hashes={p:s4.digest(folder/p) for p in ('s4_gamma.mph','fields.npz','candidates.csv')}))
            rows.append(result);record['completed'].append(record['current_case']);publish(config,rows,'running')
            print('RESULT',result,flush=True);client.remove(model);model=None
        record['status']='complete';publish(config,rows,'complete')
    except Exception as exc:
        record.update(status='failed',error=str(exc),traceback=traceback.format_exc());raise
    finally:
        if model is not None:client.remove(model)
        if client is not None:
            client.clear()
            if client.port:client.disconnect()
        record['source_unchanged']=all(s4.digest(p)==h for p,h in config['source_sha256'].items());s4.write_json(record_path,record)
    if not record['source_unchanged']:raise RuntimeError('Protected source changed')


def main():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare-only',action='store_true');group.add_argument('--run',action='store_true');group.add_argument('--report-only',action='store_true')
    parser.add_argument('--allow-concurrent',action='store_true',help='Requires explicit user authorization')
    parser.add_argument('--size-factors',nargs='+',type=float)
    parser.add_argument('--zeta',nargs='+',type=float)
    args=parser.parse_args()
    if args.size_factors:
        if any(not 0<f<1 for f in args.size_factors):parser.error('Size factors must be between zero and one')
        OUTPUT=OUTPUT/('size_refinement' if not args.zeta else 'local_zeta'+'-'.join(f'{z:g}' for z in args.zeta)+'_h'+'-'.join(f'{f:g}' for f in args.size_factors))
    elif args.zeta:parser.error('Explicit zeta requires size factors')
    config=prepare(args.size_factors,args.zeta or ZETAS)
    if args.prepare_only:print(OUTPUT,len(config['cases']),'planned solves; no COMSOL',flush=True)
    elif args.report_only:publish(config,[],'complete' if all((Path(c['directory'])/'summary.json').exists() for c in config['cases']) else 'running')
    else:run(config,args.allow_concurrent)


if __name__=='__main__':main()
