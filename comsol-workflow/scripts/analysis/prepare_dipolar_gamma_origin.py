"""Saved-data audit and exact global cell-label index-0 decomposition; no mph import."""
from __future__ import annotations

from comsol_workflow.dipolar_decomposition import complex_split, shared_node_projection, geometry_mapping, seam_metrics
import argparse, hashlib, json, shutil
from pathlib import Path
import numpy as np
from comsol_workflow.finite_lattice_fourier import direct_lattice_basis, hex_cyclic_indices, hex_xi_values, finite_lattice_fourier_components
from comsol_workflow.hex_lattice_utils import unit_cell_corners
from comsol_workflow.simulation_config_geometry import load_simulation_config_geometry

ROOT=Path(__file__).resolve().parents[2]
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
INPUTS=load_dipolar_inputs()
RESULT=INPUTS.output_dir
RUN=INPUTS.run_id
REQUEST=INPUTS.request_dir
CATEGORIES=['00_model','01_results','99_config','80_logs','12_reports','11_pdf']

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def paths():return {s:RESULT/s/INPUTS.gamma_group/RUN for s in CATEGORIES}
def input_path(category,name):return RESULT/category/INPUTS.source_batch/name





def audit():
    ps=paths()
    if any(p.exists() for p in ps.values()):raise FileExistsError('Run already exists; inspect it and use the saved config')
    snap=read(REQUEST/'parameter_snapshot.json');auth=read(REQUEST/'authorization.json')
    assert auth['implementation_authorized'] and auth['execution_authorized']
    assert snap['frequency_Hz']==197953000904452.1 and snap['mesh_label']==5
    geo_path=input_path('99_config','A_simulation_config.json');run_path=input_path('99_config','geometry_and_run_config.json')
    assert sha(geo_path)==snap['source_geometry_sha256'] and sha(run_path)==snap['source_run_sha256']
    geo=load_simulation_config_geometry(geo_path);hist=read(input_path('99_config','historical_config.json'));cfg=read(run_path)
    assert cfg['frequency_Hz']==snap['frequency_Hz'] and cfg['mesh']==5
    assert len(hist['bulk_points'])==1141 and len(hist['cladding_points'])==3540
    with np.load(input_path('01_results','A_cavity_EH.npz')) as d:
        centers=d['centers_um'];rho=d['rho_xyz_um'];assert d['E_V_m'].shape==d['H_A_m'].shape==(1141,3078,3)
    expected=np.array([[p['x'],p['y']] for p in hist['bulk_points']]);np.testing.assert_allclose(centers,expected,atol=1e-12)
    indices=np.rint(centers@np.linalg.inv(direct_lattice_basis(hist['a'])).T).astype(int)
    cyclic=hex_cyclic_indices(indices,19);assert len(set(cyclic))==1141
    corners,boundary,internal=geometry_mapping(centers,hist['a'])
    contour=np.array([corners[[edge,(edge+1)%6]]+centers[cell] for cell,edge in boundary])
    minimum_radius=float(np.min(np.linalg.norm(contour.mean(axis=1),axis=1)))
    assert cfg['source_radius_um']<minimum_radius and cfg['source_z_halfwidth_um']<.3
    files=[geo_path,run_path,input_path('99_config','historical_config.json'),input_path('00_model','A_solved.mph')]
    files += [input_path('01_results',name) for name in ['A_cavity_EH.npz','A_slab_EH.npz','A_air_EH.npz','A_air_check_EH.npz']]
    files += [REQUEST/n for n in ['source_plan.md','authorization.json','parameter_snapshot.json','request.json']]
    manifest={str(p):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files}
    for p in ps.values():p.mkdir(parents=True)
    for name in ['authorization.json','parameter_snapshot.json']:shutil.copyfile(REQUEST/name,ps['99_config']/name)
    np.savez_compressed(ps['99_config']/'cell_mapping.npz',centers_um=centers,indices=indices,cyclic_indices=cyclic,corners_um=corners,boundary_cell_edge=boundary,internal_cell_edges=internal,boundary_segments_um=contour)
    replay=dict(run_id=RUN,output_root=str(RESULT),authorization=auth,baseline=dict(model=str(input_path('00_model','A_solved.mph')),dataset='dset1',solution_number=1),frequency_Hz=cfg['frequency_Hz'],mesh=5,cell_count=1141,a_um=hist['a'],time_sign='+i omega t',source_hashes=manifest,geometry_config=str(geo_path),surface=dict(cap_z_um=.3,cap_xy_grid=[61,61],edge_samples=61,side_max_z_step_um=.005,footprint='exact complete cell union',boundary_segments=len(boundary),source_radius_um=cfg['source_radius_um'],material_interfaces_um=[-.1,.1]),decomposition=dict(method='global complete 1141-cell complex mean',cell_phase_alignment=False,independent_normalization=False),symmetry_mapping=dict(x_flip=[-1,1,1,1,-1,-1],y_flip=[-1,1,-1,1,-1,1],z_flip=[1,1,-1,-1,-1,1]),resource_budget=snap['resource_budget'],air_plane_z_um=1.1,observation_grid=dict(shape=[81,81],window_deg=[-10,10],exact_zero=True),completed_stages=['G0'],missing=['closed_surface_EH','compatible_global_Gamma_trace','exterior_model','FULL_replay'],trace_mapping='Common cell-edge samples before any averaging; opposite edges in reverse order. No seam averaging, smoothing or periodicization.')
    write(ps['99_config']/'replay_config.json',replay)
    write(ps['99_config']/'input_audit.json',dict(status='passed',source_hashes=manifest,quarter_holes=geo.hole_count,cell_count=1141,source_enclosed=True,cap_air_PML_check='pending saved model query',closed_surface='missing; export required',no_COMSOL_started=True))
    print('G0 passed: identities, 1141 cells, complete cyclic indexing, source enclosure; closed surface missing.',flush=True)
    decompose_volume(replay)
    return replay

def decompose_volume(cfg):
    ps=paths();checkpath=ps['80_logs']/'volume_decomposition_checks.json'
    if checkpath.exists():return read(checkpath)
    with np.load(input_path('01_results','A_cavity_EH.npz')) as z:d={k:z[k] for k in z.files}
    with np.load(ps['99_config']/'cell_mapping.npz') as z:order=np.argsort(z['cyclic_indices'])
    checks={};arrays={};xi=hex_xi_values(19)
    for kind in ['E_V_m','H_A_m']:
        field=d[kind];gamma,rest,check=complex_split(field);errors=[]
        for component in range(3):
            transformed=finite_lattice_fourier_components(field[order,:,component],d['rho_xyz_um'][:,:2],xi,cfg['a_um'])
            errors.append(float(np.linalg.norm(transformed['F'][0]/np.sqrt(1141)-gamma[:,component])/max(np.linalg.norm(gamma[:,component]),1e-300)))
            assert transformed['reconstruction_error']<1e-12 and transformed['parseval_error']<1e-12
        check['index0_mean_relative_error']=max(errors)
        assert max(check.values())<1e-12
        arrays['Gamma_'+kind]=gamma;arrays['REST_'+kind]=rest;checks[kind]=check
    np.savez_compressed(ps['01_results']/'volume_Gamma_REST.npz',**arrays,centers_um=d['centers_um'],rho_xyz_um=d['rho_xyz_um'],weights_half_m3=d['weights_half_m3'])
    write(checkpath,checks);print('Complete volume index-0 and REST checks passed.',flush=True);return checks


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['audit','decompose'],default='audit');p.add_argument('--config',type=Path);args=p.parse_args()
    if args.stage=='audit':audit()
    else:decompose_volume(read(args.config))
if __name__=='__main__':main()
