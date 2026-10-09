"""Authorized saved-A surface export. Never starts COMSOL without --execute."""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, time, traceback
from pathlib import Path
import numpy as np
from scripts.analysis.prepare_dipolar_gamma_origin import paths, read, write, sha, input_path
from comsol_workflow.dipolar_decomposition import seam_metrics, shared_node_projection

EXPR=[f'ewfd.E{x}/(1[V/m])' for x in 'xyz']+[f'ewfd.H{x}/(1[A/m])' for x in 'xyz']
PARITIES=np.array([[-1,1,1,1,-1,-1],[-1,1,-1,1,-1,1],[1,1,-1,-1,-1,1]])

def progress(message,**fields):
    p=paths()['80_logs']/'execution.json';old=read(p) if p.exists() else {}
    if fields.get('status')=='running':
        old.pop('exception_type',None)
        old.pop('cells_exported',None)
    old.update(message=message,python_pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**fields)
    write(p,old);print(message,flush=True)

def reflected_sample(sampler,xyz):
    shape=xyz.shape[:-1];points=xyz.reshape(-1,3)
    positive=np.round(abs(points),12);unique,inverse=np.unique(positive,axis=0,return_inverse=True)
    raw=sampler(0,unique,EXPR)[inverse]
    for axis in range(3):
        raw*=np.where(points[:,axis,None]<-1e-12,PARITIES[axis],1)
        on=np.isclose(points[:,axis],0,atol=1e-12,rtol=0)
        raw[on]*=(1+PARITIES[axis])/2
    return raw.reshape(*shape,6)

def seam_export(model,cfg):
    from scripts.run_main.run_dipolar_volume_export import NativeSampler
    from comsol_workflow.simulation_utils import electric_reflection_matrix
    ps=paths();j=model.java
    np.testing.assert_array_equal(PARITIES[0,:3],np.diag(electric_reflection_matrix('y','PMC')))
    np.testing.assert_array_equal(PARITIES[1,:3],np.diag(electric_reflection_matrix('x','PEC')))
    for signs in PARITIES:np.testing.assert_array_equal(signs[3:],-signs[:3])
    assert 'dset1' in list(j.result().dataset().tags()) and 'sol1' in list(j.sol().tags())
    assert str(j.result().dataset('dset1').getString('solution'))=='sol1'
    study={str(tag):str(j.study('std1').feature(tag).getType()) for tag in j.study('std1').feature().tags()}
    assert study.get('freq')=='Frequency'
    sampler=NativeSampler(model);freq=sampler(0,np.array([[.4,.01,.05]]),['freq/1[Hz]'])[0,0]
    np.testing.assert_allclose(freq,cfg['frequency_Hz'],rtol=1e-12)
    height=float(j.param().evaluate('H_layer_0','um'));air=float(j.param().evaluate('H_air','um'));radius=float(j.param().evaluate('s5R','um'))
    assert abs(height-.2)<1e-12 and height/2<cfg['surface']['cap_z_um']<height/2+air
    pml=list(map(int,j.component('comp1').selection('sel_pml_dom').entities(3)))
    probe=sampler(0,np.array([[.4,.01,.3],[.4,.01,.05]]),['dom','ewfd.epsilonrxx'])
    assert int(probe[0,0].real) not in pml and abs(probe[0,1]-1)<1e-9
    identity=dict(frequency_Hz=float(freq.real),study=study,dataset='dset1',solution='sol1',solution_number=1,slab_height_um=height,air_height_um=air,source_radius_um=radius,pml_domains=pml,cap_in_air_outside_PML=True,source_feature_type=str(j.component('comp1').physics('ewfd').feature('s5current').getType()))
    write(ps['80_logs']/'loaded_model_identity.json',identity)
    with np.load(ps['99_config']/'cell_mapping.npz') as z:centers=z['centers_um'];corners=z['corners_um'];internal=z['internal_cell_edges']
    with np.load(input_path('01_results','A_cavity_EH.npz')) as z:
        loc=z['rho_xyz_um'][[20,500,2000]];chosen=[100,560,1050]
        expected=np.concatenate([z['E_V_m'][chosen][:,[20,500,2000]],z['H_A_m'][chosen][:,[20,500,2000]]],-1)
    xyz=np.pad(centers[chosen],((0,0),(0,1)))[:,None,:]+loc
    actual=reflected_sample(sampler,xyz)
    error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected));assert error<1e-7,error
    write(ps['80_logs']/'native_A_reproduction.json',dict(relative_error=error,independent_phase_alignment=False))
    u=np.linspace(0,1,61);edges=corners[:,None,:]*(1-u)[None,:,None]+np.roll(corners,-1,axis=0)[:,None,:]*u[None,:,None]
    local=np.concatenate([edges,np.full((6,61,1),cfg['surface']['cap_z_um'])],-1)
    values=np.empty((1141,6,61,6),complex)
    for start in range(0,len(centers),64):
        end=min(start+64,len(centers));xyz=local[None]+np.pad(centers[start:end],((0,0),(0,1)))[:,None,None,:]
        values[start:end]=reflected_sample(sampler,xyz)
        progress(f'G1 cap seam export: {end}/1141 complete cells',status='running',phase='cap_seam_export',cells_exported=end)
    gamma,records=seam_metrics(values,corners);rest=values-gamma
    lhs=values[internal[:,0,0],internal[:,0,1]];rhs=values[internal[:,1,0],internal[:,1,1],::-1]
    full_error=float(np.linalg.norm(lhs-rhs)/np.linalg.norm(lhs))
    closure=float(np.linalg.norm(values-gamma-rest)/np.linalg.norm(values))
    rest_mean=float(np.linalg.norm(rest.mean(0))/(np.linalg.norm(values)/np.sqrt(len(values))))
    assert full_error<1e-7 and closure<1e-12 and rest_mean<1e-12
    np.savez_compressed(ps['01_results']/'cap_edge_seam_EH.npz',FULL_EH=values,GAMMA_EH=gamma,REST_EH=rest,local_xyz_um=local,centers_um=centers,edge_tangent_xy=np.roll(corners,-1,axis=0)-corners,cap_z_um=cfg['surface']['cap_z_um'])
    result=dict(status='audit_ready',scope='top cap cell seams; not yet the complete closed surface',Gamma_edges=records,FULL_shared_edge_relative_error=full_error,FULL_GAMMA_REST_relative_error=closure,REST_cell_mean_relative_error=rest_mean,definition_changed=False,source_model_sha256=sha(cfg['baseline']['model']),FULL_replay_started=False)
    write(ps['80_logs']/'cap_seam_checks.json',result)
    progress('G1 cap seam audit exported; inspect compatibility before full trace construction',status='seam_audit_ready',phase='seam_gate')

def triangle_areas(vertices):
    a=vertices[:,1]-vertices[:,0];b=vertices[:,2]-vertices[:,0]
    return abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2

def export_common_grid(sampler,centers,local,path,label):
    meta=path.with_suffix('.progress.json');identity=hashlib.sha256(centers.tobytes()+local.tobytes()).hexdigest()
    if path.exists():
        state=read(meta);assert state['coordinate_sha256']==identity
        values=np.lib.format.open_memmap(path,mode='r+');assert values.shape==(len(centers),len(local),6)
        first=state['cells_complete']
    else:
        values=np.lib.format.open_memmap(path,mode='w+',dtype=np.complex128,shape=(len(centers),len(local),6));first=0
        write(meta,dict(coordinate_sha256=identity,cells_complete=0))
    for start in range(first,len(centers),8):
        end=min(start+8,len(centers));xyz=local[None]+np.pad(centers[start:end],((0,0),(0,1)))[:,None,:]
        values[start:end]=reflected_sample(sampler,xyz)
        values.flush();write(meta,dict(coordinate_sha256=identity,cells_complete=end))
        if start%64==0:progress(f'G1 {label}: {end}/1141 cells',status='running',phase=label,cells_exported=end)
    values.flush();return values


def complete_surface_export(model,cfg):
    from scipy.spatial import Delaunay,cKDTree
    from scripts.run_main.run_dipolar_volume_export import NativeSampler
    from comsol_workflow.dipolar_decomposition import complex_split
    ps=paths();sampler=NativeSampler(model)
    with np.load(ps['99_config']/'cell_mapping.npz') as z:centers=z['centers_um'];corners=z['corners_um'];boundary=z['boundary_cell_edge']
    u=np.linspace(0,1,61);edges=corners[:,None,:]*(1-u)[None,:,None]+np.roll(corners,-1,axis=0)[:,None,:]*u[None,:,None]
    tangent=np.roll(corners,-1,axis=0)-corners;normal=np.column_stack([tangent[:,1],-tangent[:,0]]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    axisx=np.linspace(-cfg['a_um']/2,cfg['a_um']/2,61);axisy=np.linspace(-cfg['a_um']/np.sqrt(3),cfg['a_um']/np.sqrt(3),61)
    xx,yy=np.meshgrid(axisx,axisy);xy=np.column_stack([xx.ravel(),yy.ravel()]);xy=xy[np.all(xy@normal.T<=cfg['a_um']/2+1e-12,axis=1)]
    xy=np.unique(np.round(np.concatenate([xy,edges.reshape(-1,2)]),12),axis=0)
    triangles=Delaunay(xy).simplices;verts=xy[triangles];areas=triangle_areas(verts)
    assert abs(areas.sum()-np.sqrt(3)*cfg['a_um']**2/2)<1e-10
    weights=np.zeros(len(xy));np.add.at(weights,triangles.ravel(),np.repeat(areas/3,3))
    local=np.column_stack([xy,np.full(len(xy),.3)])
    full=export_common_grid(sampler,centers,local,ps['01_results']/'all_cell_cap_EH.npy','cap_grid')
    gamma=np.asarray(full).mean(0);globalxy=(centers[:,None,:]+xy[None]).reshape(-1,2)
    unique,mapping,counts=np.unique(np.round(globalxy,10),axis=0,return_inverse=True,return_counts=True)
    projected=shared_node_projection(np.tile(gamma,(1141,1)),mapping,counts)
    actual=shared_node_projection(np.asarray(full).reshape(-1,6),mapping,counts)
    # Q_h is one fixed linear nodal trace map: identify duplicate global nodes
    # with equal weights, then common P1 face interpolation. Raw means stay saved.
    gamma_change=float(np.linalg.norm(projected[mapping]-np.tile(gamma,(1141,1)))/np.linalg.norm(np.tile(gamma,(1141,1))))
    full_change=float(np.linalg.norm(actual[mapping]-np.asarray(full).reshape(-1,6))/np.linalg.norm(full))
    assert full_change<1e-7
    cap_path=ps['01_results']/'cap_boundary_inputs.npz'
    np.savez_compressed(cap_path,xy_um=unique,FULL_EH=actual,GAMMA_EH=projected,REST_EH=actual-projected,raw_Gamma_local_EH=gamma,local_xy_um=xy,local_triangles=triangles,local_weights_um2=weights,cell_to_global_node=mapping.reshape(1141,-1),cap_z_um=.3)
    chart_data=[]
    for chart,(lo,hi,nz) in enumerate([(0.,.1,21),(.1,.3,41)]):
        z=np.linspace(lo,hi,nz);samplez=z.copy()
        if chart==0:samplez[-1]-=1e-7
        else:samplez[0]+=1e-7
        side_local=np.concatenate([np.broadcast_to(edges[:,None,:,:],(6,nz,61,2)),np.broadcast_to(samplez[None,:,None,None],(6,nz,61,1))],-1).reshape(-1,3)
        side=export_common_grid(sampler,centers,side_local,ps['01_results']/f'all_cell_sides_chart{chart}_EH.npy',f'side_chart{chart}')
        g=np.asarray(side).mean(0).reshape(6,nz,61,6)
        raw=np.asarray(side).reshape(1141,6,nz,61,6)
        face_full=np.array([raw[cell,edge] for cell,edge in boundary]);face_gamma=np.array([g[edge] for cell,edge in boundary])
        coords=np.array([edges[edge]+centers[cell] for cell,edge in boundary])
        # A boundary vertex can also belong to a cell with no exposed edge at it.
        # Include all incident cell-local edge views in Q_h, not only surface faces.
        allxy=np.round((centers[:,None,None,:]+edges[None]).reshape(-1,2),10)
        node_xy,node_map,node_counts=np.unique(allxy,axis=0,return_inverse=True,return_counts=True)
        lookup=cKDTree(node_xy);distance,face_map=lookup.query(coords.reshape(-1,2));assert distance.max()<1e-8
        qgamma=np.empty_like(face_gamma)
        for zi in range(nz):
            accum=shared_node_projection(np.tile(g[:,zi].reshape(-1,6),(1141,1)),node_map,node_counts)
            qgamma[:,zi]=accum[face_map].reshape(len(boundary),61,6)
        lower=np.zeros((len(boundary),nz,61,3));lower[...,:2]=coords[:,None];lower[...,2]=z[None,:,None]
        sw=np.ones(nz);sw[[0,-1]]=.5;uw=np.ones(61);uw[[0,-1]]=.5
        areaweights=np.linalg.norm(tangent[0])*(hi-lo)/(60*(nz-1))*sw[:,None]*uw[None,:]
        path=ps['01_results']/f'side_boundary_inputs_chart{chart}.npz'
        np.savez_compressed(path,xyz_um=lower,sampling_z_um=samplez,FULL_EH=face_full,GAMMA_EH=qgamma,REST_EH=face_full-qgamma,raw_Gamma_local_EH=g,normals=np.column_stack([normal[boundary[:,1]],np.zeros(len(boundary))]),weights_um2=areaweights,boundary_cell_edge=boundary,chart=chart,material_side='slab' if chart==0 else 'air')
        chart_data.append(dict(chart=chart,z_bounds_um=[lo,hi],logical_nodes=nz,one_sided_offset_um=1e-7,path=str(path)))
    definition=dict(raw_Gamma='Arithmetic complex mean across all 1141 cells at each common local sample',Q_h='Equal-weight identification of coincident global nodes, then one common linear nodal interpolation per planar/material chart; same mapping for FULL and Gamma; REST=Q_FULL-Q_Gamma',cap_raw_Gamma_to_Q_relative=gamma_change,cap_FULL_to_Q_relative=full_change,side_charts=chart_data,full_lower_half='Reflect z with E=(+,+,-), H=(-,-,+); coordinates and outward normals reflected. Stored upper-half coefficients plus this exact mapping represent both caps and all closed sides.',no_phase_alignment=True,no_independent_rescaling=True,source_frequency_Hz=cfg['frequency_Hz'])
    write(ps['99_config']/'gamma_definition.json',definition)
    write(ps['80_logs']/'surface_export_complete.json',dict(status='complete',source_sha256=sha(cfg['baseline']['model']),phase='G1',closed_surface_representation='upper-half E/H plus explicit exact zPMC reflection',Gamma_definition_preserved=True))
    progress('G1 closed boundary export and fixed linear nodal trace mapping complete',status='export_complete',phase='G1_complete')


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--stage',choices=['seam','export'],required=True);p.add_argument('--execute',action='store_true');args=p.parse_args()
    if not args.execute:raise ValueError('COMSOL export requires explicit --execute')
    cfg=read(args.config);assert cfg['authorization']['saved_field_export_authorized']
    done=paths()['80_logs']/'cap_seam_checks.json'
    if done.exists():
        assert read(done)['source_model_sha256']==cfg['source_hashes'][cfg['baseline']['model']]['sha256']
        if args.stage=='seam':return
    for path,meta in cfg['source_hashes'].items():
        if sha(path)!=meta['sha256']:raise ValueError('Input hash changed: '+path)
    from scripts.run_main.run_boundary_analysis import runtime_preflight
    write(paths()['80_logs']/'runtime_preflight.json',runtime_preflight())
    import mph
    from jpype import JClass
    client=mph.start(version='6.3',cores=16);util=JClass('com.comsol.model.util.ModelUtil')
    try:
        if not util.checkoutLicense('COMSOL','WAVEOPTICS'):raise RuntimeError('COMSOL/WAVEOPTICS checkout failed')
        proc=subprocess.run(['powershell','-NoProfile','-Command',"Get-CimInstance Win32_Process -Filter \"Name='comsolmphserver.exe'\" | Select-Object ProcessId | ConvertTo-Json -Compress"],capture_output=True,text=True,check=True)
        server=json.loads(proc.stdout)
        log=paths()['80_logs']/'COMSOL_export_progress.log';util.showProgress(str(log))
        progress('COMSOL started, licenses checked out; loading saved A (no new solve)',status='running',phase='loading_saved_A',comsol_pid=int(server['ProcessId']),progress_log=str(log))
        model=client.load(cfg['baseline']['model'])
        if not done.exists():seam_export(model,cfg)
        else:progress('Reusing identity-matched completed seam export',phase='resume_G1')
        if args.stage=='export':complete_surface_export(model,cfg)
    except Exception as exc:
        progress(str(exc),status='failed',phase='export',exception_type=type(exc).__name__);traceback.print_exc();raise
    finally:
        util.showProgress(False);client.clear()
        if getattr(client,'server',None):client.server.stop()
if __name__=='__main__':main()
