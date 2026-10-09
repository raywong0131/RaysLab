"""One fixed exterior operator: FULL replay gate followed by GAMMA and REST."""
from __future__ import annotations
import argparse,json,subprocess,time,traceback
from pathlib import Path
import numpy as np
from scipy.constants import c
from scripts.analysis.prepare_dipolar_gamma_origin import paths, read, write, sha, input_path
from scripts.run_main.export_dipolar_gamma_boundary import progress, reflected_sample
from comsol_workflow.dipolar_radiation import plane_farfield
from scripts.analysis.analyze_dipolar_cladding_scattering import Z0


def contour_from_edges(segments):
    nxt={tuple(np.round(s[0],10)):tuple(np.round(s[1],10)) for s in segments}
    assert len(nxt)==len(segments)
    first=next(iter(nxt));current=first;vertices=[]
    while True:
        vertices.append(current);current=nxt[current]
        if current==first:break
        if len(vertices)>len(nxt):raise ValueError('Open cell union contour')
    assert len(vertices)==len(nxt)
    return np.array(vertices)


def planar_chart_boundaries(face_vertices,chart,tolerance):
    normal=np.array([0.,0.,1.]) if chart['kind']=='cap' else np.r_[chart['normal'],0.]
    offset=chart['zlo'] if chart['kind']=='cap' else chart['offset_um']
    return [bid for bid,xyz in face_vertices.items()
            if np.all(abs(xyz@normal-offset)<=tolerance)
            and xyz[:,2].min()>=chart['zlo']-tolerance
            and xyz[:,2].max()<=chart['zhi']+tolerance]


def create_tables():
    ps=paths();data=ps['01_results'];tables=data/'trace_tables';tables.mkdir(exist_ok=True)
    with np.load(data/'cap_boundary_inputs.npz') as z:
        xy=z['xy_um'];keep=(xy[:,0]>=-1e-10)&(xy[:,1]>=-1e-10)
        xy=xy[keep];f=z['FULL_EH'][keep,:3];g=z['GAMMA_EH'][keep,:3]
    def table(name,coords,f,g):
        arr=np.column_stack([coords,f.real,f.imag,g.real,g.imag]);path=tables/(name+'.txt')
        np.savetxt(path,arr,fmt='%.16e',header='coordinate1 coordinate2 FULL_Ex_re FULL_Ey_re FULL_Ez_re FULL_Ex_im FULL_Ey_im FULL_Ez_im GAMMA_Ex_re GAMMA_Ey_re GAMMA_Ez_re GAMMA_Ex_im GAMMA_Ey_im GAMMA_Ez_im',comments='% ')
        return str(path)
    charts=[dict(tag='gocap',kind='cap',path=table('cap',xy,f,g),arguments=['x/1[um]','y/1[um]'],zlo=.3,zhi=.3)]
    for chart in [0,1]:
        with np.load(data/f'side_boundary_inputs_chart{chart}.npz') as z:d={k:z[k] for k in z.files}
        xyz=d['xyz_um'];normals=d['normals'];groups={}
        for i in range(len(xyz)):
            normal=normals[i,:2];offset=float(normal@xyz[i,0,0,:2]);key=tuple(np.round(np.r_[normal,offset],9))
            groups.setdefault(key,[]).append(i)
        for key,ids in groups.items():
            n=np.array(key[:2]);t=np.array([-n[1],n[0]])
            pts=xyz[ids].reshape(-1,3);f=d['FULL_EH'][ids,...,:3].reshape(-1,3);g=d['GAMMA_EH'][ids,...,:3].reshape(-1,3)
            inside=(pts[:,0]>1e-8)&(pts[:,1]>1e-8)
            if not np.any(inside):continue
            seed=pts[np.flatnonzero(inside)[len(np.flatnonzero(inside))//2]].copy()
            lo=float(pts[:,2].min());hi=float(pts[:,2].max());seed[2]=(lo+hi)/2
            coords=np.column_stack([pts[:,:2]@t,pts[:,2]]);unique,inv,count=np.unique(np.round(coords,10),axis=0,return_inverse=True,return_counts=True)
            fv=np.zeros((len(unique),3),complex);gv=fv.copy();np.add.at(fv,inv,f);np.add.at(gv,inv,g);fv/=count[:,None];gv/=count[:,None]
            tag='gos'+str(len(charts));charts.append(dict(tag=tag,kind='side',chart=chart,path=table(tag,unique,fv,gv),arguments=[f'({t[0]:.17g}*x+{t[1]:.17g}*y)/1[um]','z/1[um]'],normal=n.tolist(),offset_um=key[2],seed_um=seed.tolist(),zlo=lo,zhi=hi))
    write(ps['99_config']/'interpolation_charts.json',charts);return charts


def install_interpolation(j,chart):
    from jpype.types import JArray,JString
    tag=chart['tag'];fn=j.func().create(tag,'Interpolation');fn.set('source','file');fn.set('filename',chart['path']);fn.set('nargs','2')
    names=[tag+'_'+group+axis+part for group in ['F','G'] for part in ['r','i'] for axis in 'xyz']
    fn.set('funcs',JArray(JString,2)([[name,str(i+1)] for i,name in enumerate(names)]));fn.set('argunit',['1','1']);fn.set('fununit',['V/m']*12);fn.set('interp','linear');fn.set('extrap','const');fn.importData()
    args=','.join(chart['arguments']);return [f'goFull*({tag}_F{a}r({args})+i*{tag}_F{a}i({args}))+goGamma*({tag}_G{a}r({args})+i*{tag}_G{a}i({args}))' for a in 'xyz']


def build_exterior(model,cfg,charts,util):
    from jpype.types import JArray,JString,JDouble,JInt
    from comsol_workflow.polygon_utils import clip_polygon_to_convex_region
    from scripts.run_main.run_dipolar_cladding_scattering import configure_study
    ps=paths();j=model.java;comp=j.component('comp1');geom=comp.geom('geom1')
    with np.load(ps['99_config']/'cell_mapping.npz') as z:contour=contour_from_edges(z['boundary_segments_um'])
    quarter=clip_polygon_to_convex_region(contour,np.array([[0.,0.],[100.,0.],[100.,100.],[0.,100.]]))
    assert len(quarter)>3
    if 'goprism' not in list(geom.feature().tags()):
        for tag in list(j.sol().tags()):j.sol().remove(tag)
        for tag in list(j.result().numerical().tags()):j.result().numerical().remove(tag)
        # Add the prism as an intersecting solid. Form Union retains its internal
        # boundaries; existing material selectors continue to use position/Boolean rules.
        geom.create('gowp','WorkPlane');wp=geom.feature('gowp');wp.set('quickplane','xy');wp.set('quickz','0[um]')
        wg=wp.geom();wg.create('gopoly','Polygon');wg.feature('gopoly').set('source','table');wg.feature('gopoly').set('table',JArray(JDouble,2)(quarter.tolist()));wg.feature('gopoly').set('type','solid')
        geom.create('goprism','Extrude');cut=geom.feature('goprism');cut.selection('input').set(['gowp']);cut.set('distance',['0.3[um]']);cut.set('selresult','on');cut.set('selresultshow','all')
        # Finalize uses the saved Form Union; internal partitions are retained.
        progress('Finalize geometry properties: '+str(list(geom.feature('fin').properties())),phase='geometry_setup')
        progress('G2 partitioning saved actual geometry at the complete-cell prism',status='running',phase='geometry')
        util.showProgress(str(ps['80_logs']/'FULL_COMSOL_progress.log'));geom.run()
        model.save(str(ps['00_model']/'exterior_partition.mph'))
        progress('Geometry partition checkpoint saved',phase='partition_complete')
    sel=comp.selection();tags=list(sel.tags());assert 'geom1_goprism_dom' in tags,tags[-20:]
    interior=list(map(int,comp.selection('geom1_goprism_dom').entities(3)));assert interior
    sel.create('goall','Explicit');comp.selection('goall').geom('geom1',3);comp.selection('goall').all()
    sel.create('gooutside','Difference');comp.selection('gooutside').set('entitydim',JInt(3));comp.selection('gooutside').set('add',['goall']);comp.selection('gooutside').set('subtract',['geom1_goprism_dom'])
    exterior=list(map(int,comp.selection('gooutside').entities(3)));assert exterior and not set(interior)&set(exterior)
    for name,inp in [('goinb','geom1_goprism_dom'),('gooutb','gooutside')]:
        sel.create(name,'Adjacent');comp.selection(name).set('entitydim',JInt(3));comp.selection(name).set('outputdim',JInt(2));comp.selection(name).set('input',[inp])
    sel.create('goboundary','Intersection');comp.selection('goboundary').set('entitydim',JInt(2));comp.selection('goboundary').set('input',['goinb','gooutb'])
    inner=list(map(int,comp.selection('goboundary').entities(2)));assert inner
    physics=comp.physics('ewfd');physics.selection().named('gooutside');physics.feature('s5current').active(False)
    j.param().set('goFull','1');j.param().set('goGamma','0')
    selection_tol=float(j.param().evaluate('selection_tol','um'))
    vertices=np.asarray(geom.getVertexCoord()).T
    face_vertices={bid:vertices[np.asarray(geom.getAdj(2,0,JInt(bid)))-1] for bid in inner}
    coverage=[];chart_records=[]
    for chart in charts:
        tag=chart['tag'];select=tag+'sel'
        ids=planar_chart_boundaries(face_vertices,chart,selection_tol)
        if not ids:raise ValueError('No inner boundaries for chart '+tag)
        sel.create(select,'Explicit');comp.selection(select).geom('geom1',JInt(2));comp.selection(select).set(JArray(JInt)(ids))
        if set(ids)&set(coverage):raise ValueError('Overlapping chart selections '+tag)
        coverage+=ids;expressions=install_interpolation(j,chart)
        if len(chart_records)%10==0:progress('Configuring boundary chart '+str(len(chart_records)+1)+'/'+str(len(charts)),phase='boundary_conditions')
        feature=physics.create(tag,'ElectricField',2);feature.selection().named(select);feature.set('E0',JArray(JString,1)(expressions))
        chart_records.append(dict(tag=tag,boundaries=ids,expressions=expressions))
    if set(coverage)!=set(inner):raise ValueError('Incomplete closed boundary selection: '+str(set(inner)-set(coverage)))
    # Reuse saved mesh settings, remove the prism interior from free tetrahedra.
    sel.create('gophysical','Difference');comp.selection('gophysical').set('entitydim',JInt(3));comp.selection('gophysical').set('add',['sel_physical_dom']);comp.selection('gophysical').set('subtract',['geom1_goprism_dom'])
    mesh=comp.mesh('mesh1');mesh.feature('ftet1').selection().named('gophysical')
    progress('G2 building one exterior mesh at the saved mesh5 settings',phase='mesh')
    mesh.run();configure_study(model)
    np.testing.assert_allclose(j.param().evaluate('s5freq'),cfg['frequency_Hz'],rtol=1e-12)
    write(ps['99_config']/'exterior_operator.json',dict(interior_domains=interior,exterior_domains=exterior,inner_boundaries=inner,charts=chart_records,original_current_disabled=True,cladding_geometry_unchanged=True,quarter_prism_polygon_um=quarter.tolist(),frequency_Hz=cfg['frequency_Hz'],mesh=5))
    model.save(str(ps['00_model']/'exterior_base.mph'))
    progress('G2 exterior mesh checkpoint saved',phase='mesh_complete')


def export_air(model,label):
    from scripts.run_main.run_dipolar_volume_export import NativeSampler
    ps=paths();sampler=NativeSampler(model)
    with np.load(input_path('01_results','A_air_EH.npz')) as z:axis=z['x_um'];valid=z['valid'];height=float(z['z_um'])
    xx,yy=np.meshgrid(axis,axis);points=np.column_stack([xx[valid],yy[valid],np.full(valid.sum(),height)])
    raw=reflected_sample(sampler,points);e=np.zeros((*valid.shape,3),complex);h=e.copy();e[valid]=raw[:,:3];h[valid]=raw[:,3:]
    np.savez_compressed(ps['01_results']/f'{label}_air_EH.npz',E_V_m=e,H_A_m=h,x_um=axis,y_um=axis,z_um=height,valid=valid)
    # Exterior-only slab response: never interpolate the excluded cavity interior.
    from matplotlib.path import Path as PolygonPath
    with np.load(ps['99_config']/'cell_mapping.npz') as z:contour=contour_from_edges(z['boundary_segments_um'])
    outside=valid & ~PolygonPath(contour).contains_points(np.column_stack([xx.ravel(),yy.ravel()]),radius=1e-7).reshape(valid.shape)
    slab=reflected_sample(sampler,np.column_stack([xx[outside],yy[outside],np.zeros(outside.sum())]))
    es=np.zeros((*valid.shape,3),complex);hs=es.copy();es[outside]=slab[:,:3];hs[outside]=slab[:,3:]
    np.savez_compressed(ps['01_results']/f'{label}_exterior_slab_EH.npz',E_V_m=es,H_A_m=hs,x_um=axis,y_um=axis,z_um=0.,valid=outside)
    return e,h,axis,height


def farfield(cfg,label):
    with np.load(paths()['01_results']/f'{label}_air_EH.npz') as z:e=z['E_V_m'];h=z['H_A_m'];axis=z['x_um'];height=float(z['z_um'])
    k0=2*np.pi*cfg['frequency_Hz']/c;q=np.linspace(-k0*np.sin(np.deg2rad(10)),k0*np.sin(np.deg2rad(10)),81)
    return plane_farfield(e,h,axis,height,q,k0),q


def replay_check(cfg):
    ps=paths()
    with np.load(input_path('01_results','A_air_EH.npz')) as z:ref={k:z[k] for k in z.files}
    with np.load(ps['01_results']/'FULL_air_EH.npz') as z:full={k:z[k] for k in z.files}
    f,q=farfield(cfg,'FULL');k0=2*np.pi*cfg['frequency_Hz']/c
    original=plane_farfield(ref['E_V_m'],ref['H_A_m'],ref['x_um'],float(ref['z_um']),q,k0)
    checks={key:float(np.linalg.norm(full[key]-ref[key])/np.linalg.norm(ref[key])) for key in ['E_V_m','H_A_m']}
    checks['center_complex_relative_error']=float(np.linalg.norm(f[40,40]-original[40,40])/np.linalg.norm(original[40,40]))
    checks['farfield_complex_relative_error']=float(np.linalg.norm(f-original)/np.linalg.norm(original))
    checks['pilot_limit']=.1;checks['passed']=all(v<.1 for k,v in checks.items() if k not in ['pilot_limit','passed'])
    checks['no_independent_phase_or_scale_fit']=True;checks['source_hashes']=cfg['source_hashes']
    write(ps['80_logs']/'FULL_replay_check.json',checks)
    np.savez_compressed(ps['01_results']/'FULL_reference_farfield.npz',FULL=f,reference_A=original,q_m_inv=q)
    return checks


def solve_cases(model,cfg,stage,client,util):
    ps=paths()
    labels=['FULL'] if stage=='full' else ['FULL','GAMMA','REST']
    for label in labels:
        done=ps['80_logs']/f'{label}_complete.json'
        if done.exists():
            assert read(done)['source_hashes']==cfg['source_hashes']
            if label=='FULL':
                check=replay_check(cfg)
                if not check['passed']:raise ValueError('Saved FULL replay failed; component solves gated')
            continue
        if label!='FULL':
            gate=read(ps['80_logs']/'FULL_replay_check.json');assert gate['passed'] and gate['source_hashes']==cfg['source_hashes']
        solved=ps['00_model']/f'{label}_solved.mph'
        if solved.exists():client.remove(model);model=client.load(str(solved))
        else:
            wf,wg={'FULL':(1,0),'GAMMA':(0,1),'REST':(1,-1)}[label]
            model.java.param().set('goFull',str(wf));model.java.param().set('goGamma',str(wg))
            util.showProgress(str(ps['80_logs']/f'{label}_COMSOL_progress.log'))
            progress(label+': solving fixed exterior frequency response',phase='solving',case=label)
            model.java.sol('sol1').runAll();model.save(str(solved))
        progress(label+': exporting common air-plane E/H',phase='exporting',case=label);export_air(model,label)
        write(done,dict(status='complete',model=str(solved),source_hashes=cfg['source_hashes']))
        if label=='FULL':
            check=replay_check(cfg)
            if not check['passed']:
                progress('FULL replay did not meet the recorded pilot criterion; component solves gated',status='implementation_unresolved',phase='FULL_replay_check');return
    if stage=='all':
        from scripts.analysis.analyze_dipolar_gamma_origin import report
        report(cfg)
    progress('Requested exterior cases completed',status='complete',phase='replay_complete')


def run(cfg,stage):
    from scripts.run_main.run_boundary_analysis import runtime_preflight
    ps=paths();record=read(ps['80_logs']/'surface_export_complete.json');assert record['status']=='complete'
    for path,meta in cfg['source_hashes'].items():
        if sha(path)!=meta['sha256']:raise ValueError('Source identity changed: '+path)
    write(ps['80_logs']/'replay_runtime_preflight.json',runtime_preflight())
    import mph
    from jpype import JClass
    client=mph.start(version='6.3',cores=16);util=JClass('com.comsol.model.util.ModelUtil')
    try:
        assert util.checkoutLicense('COMSOL','WAVEOPTICS')
        proc=subprocess.run(['powershell','-NoProfile','-Command',"Get-CimInstance Win32_Process -Filter \"Name='comsolmphserver.exe'\" | Select-Object ProcessId | ConvertTo-Json -Compress"],capture_output=True,text=True,check=True)
        progress('COMSOL started for exterior FULL replay',status='running',phase='loading_exterior',comsol_pid=int(json.loads(proc.stdout)['ProcessId']))
        base=ps['00_model']/'exterior_base.mph'
        if base.exists():model=client.load(str(base))
        else:
            partition=ps['00_model']/'exterior_partition.mph'
            model=client.load(str(partition) if partition.exists() else cfg['baseline']['model']);charts=read(ps['99_config']/'interpolation_charts.json') if (ps['99_config']/'interpolation_charts.json').exists() else create_tables();build_exterior(model,cfg,charts,util)
        solve_cases(model,cfg,stage,client,util)
    except Exception as exc:
        progress(str(exc),status='failed',exception_type=type(exc).__name__);traceback.print_exc();raise
    finally:
        util.showProgress(False);client.clear()
        if getattr(client,'server',None):client.server.stop()


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--stage',choices=['full','all'],default='full');p.add_argument('--execute',action='store_true');args=p.parse_args()
    if not args.execute:raise ValueError('Explicit --execute required')
    cfg=read(args.config);assert cfg['authorization']['execution_authorized'];run(cfg,args.stage)
if __name__=='__main__':main()
