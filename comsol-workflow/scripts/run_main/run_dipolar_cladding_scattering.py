"""Authorized S5 single-source/single-frequency A/B/C driven COMSOL experiment."""
from __future__ import annotations
import argparse, hashlib, json, os, time, traceback
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scripts.analysis.analyze_dipolar_cladding_scattering import read, put, load, WORK, CONFIG, DATA, metric, project
from scripts.analysis.analyze_dipolar_singularity import digest

FREQ=197.9530009044521e12
MODELS=['B','A','C']
SOURCE_MODEL=WORK/'00_model/finite_mode19_source.mph'
PERIODIC_MODEL=WORK/'00_model/periodic_k00.mph'
STATE={}

def progress(out,message,**kw):
    STATE.update(kw,message=message,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    put(out/'80_logs/execution.json',STATE);print(message,flush=True)

def fingerprint(p):
    return dict(path=str(p.resolve()),bytes=p.stat().st_size,sha256=digest(p))

def geometry_snapshots(out,cfg):
    from comsol_workflow.geometry_utils import create_hexagon_design
    from comsol_workflow.polygon_utils import clip_polygon_to_convex_region,polygon_area
    from scripts.run_main.run_band_pair import get_hole_params
    from comsol_workflow.simulation_config_geometry import load_simulation_config_geometry
    shapes={k:create_hexagon_design(cfg['a'],get_hole_params(cfg[key]))[1] for k,key in [('cavity','bulk_params'),('cladding','cladding_params')]}
    boundary=np.asarray(cfg['quarter_boundary']);full={};quarter={}
    for label in MODELS:
        full[label]=[];quarter[label]=[]
        for kind,key in [('cavity','bulk_points'),('cladding','cladding_points')]:
            if label=='C' and kind=='cladding':continue
            hs=shapes['cavity' if label=='B' else kind]
            for pt in cfg[key]:
                center=np.array([pt['x'],pt['y']])
                for h in hs:
                    poly=np.asarray(h)+center;full[label].append(poly)
                    if np.any(poly.max(0)<-1e-10):continue
                    clipped=clip_polygon_to_convex_region(poly,boundary)
                    if len(clipped)>=3 and polygon_area(clipped)>1e-12:quarter[label].append(clipped)
        snap=dict(length_unit='um',footprint=boundary.tolist(),layers=[dict(height=cfg['simulation_common']['slab_height'],holes=[p.tolist() for p in quarter[label]])],model=label,frequency_Hz=FREQ,mesh=5,source_historical_config_sha256=digest(CONFIG/'config.json'))
        path=out/'99_config'/f'{label}_simulation_config.json';put(path,snap);parsed=load_simulation_config_geometry(path);assert parsed.hole_count==len(quarter[label])
    assert len(quarter['A'])==cfg['quarter_hole_count']
    # Every cavity hole and the slab footprint is identical; changes are only in peripheral cells.
    nc=6*len(cfg['bulk_points'])
    for label in ['B','C']:np.testing.assert_allclose(full[label][:nc],full['A'][:nc],atol=1e-13)
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from scripts.analysis.analyze_dipolar_cladding_scattering import save
    fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for ax,label in zip(axs,['A','B','C']):
        ax.add_collection(PolyCollection(full[label],facecolors='white',edgecolors='none'));ax.set_facecolor('#444c66');bound=np.asarray(cfg['finite_boundary']);ax.plot(*np.vstack([bound,bound[0]]).T,color='k',lw=.5);ax.set(xlim=(-33,33),ylim=(-29,29),aspect='equal',title={'A':'A: actual cavity-cladding','B':'B: cavity lattice continuation','C':'C: uniform peripheral slab'}[label],xlabel='x (um)',ylabel='y (um)')
    save(fig,out,'GS3a_geometry')
    return quarter

def prepare(out):
    if (out/'99_config/geometry_and_run_config.json').exists():raise FileExistsError('Run already prepared; resume without --prepare, never replace its configuration')
    cfg=read(CONFIG/'config.json');assert cfg['simulation_common']['mesh_auto_size']==5
    assert len(cfg['bulk_points'])==1141 and len(cfg['cladding_points'])==3540
    assert cfg['symmetry']==dict(id=1,x_boundary='PEC',y_boundary='PMC')
    identity=read(WORK/'80_logs/finite_model_identity.json');assert identity['parameters']['H_layer_0']=='200 [nm]'
    np.testing.assert_allclose(identity['frequency_THz']['re']*1e12,FREQ,rtol=1e-13)
    progress(out,'Preparing historical geometry and immutable input fingerprints',status='preparing')
    geometry_snapshots(out,cfg)
    inputs=[CONFIG/'config.json',SOURCE_MODEL,PERIODIC_MODEL,DATA/'periodic_k00_EH.npz',DATA/'finite_cavity_EH_o3.npz',DATA/'unitcell_clipped_integrals.npz']
    inputs += [DATA/f'finite_mesh_EH_{i:03d}.npz' for i in range(56)]
    manifest={str(p):fingerprint(p) for p in inputs}
    # Circle lies well inside the cavity union, with a zero slope at its edge.
    source_radius=.65*cfg['a']*cfg['bulk_radius']
    run=dict(user_authorization='2026-09-28: 按照当前spec的方案，直接执行。授权启动COMSOL等在内的高消耗资源计算',user_mode='py',internal_mode='px',mode_idx=19,frequency_Hz=FREQ,periodic_reference_frequency_Hz=198.3868382311348e12,mesh=5,solve_order=MODELS,number_of_initial_solves=3,source_radius_um=source_radius,source_z_halfwidth_um=.1,source_norm_target_A2_per_m=1.,source='Jamp * Gamma electric spatial template * cos(pi*r/(2*R))^2 * cos(pi*z/(0.2um))^2; zero outside r<R, |z|<0.1um',source_table_shape=[65,65,9],cores=16,estimated_memory_limit_GiB=768,walltime_budget_hours=24,concurrent_COMSOL_jobs=1,quarter_symmetry=cfg['symmetry'],z_symmetry='PMC at z=0; full 200nm slab reflected',source_norm_is_not_input_power=True,source_window_sensitivity_pending=True,frequency_sensitivity_pending=True,cache_policy='new output only; existing solved checkpoints may resume only with same input manifest',inputs=manifest,output=str(out),acceptance=dict(mode19_energy_overlap_min=.9,farfield_intensity_cosine_min=.9,source_norm_relative_tolerance=.005,small_projection_fraction=.01),geometry=dict(a_um=cfg['a'],slab_height=cfg['simulation_common']['slab_height'],refractive_index=cfg['simulation_common']['refractive_index'],cavity_cells=1141,cladding_cells=3540,shared_parameters_modified=False))
    put(out/'99_config/geometry_and_run_config.json',run);put(out/'99_config/historical_config.json',cfg)
    progress(out,'Preflight snapshots ready; COMSOL not started',status='prepared')


def source_table(client,out,cfg,run):
    path=out/'01_results/source_template.txt';meta=out/'80_logs/source_template.json'
    if path.exists() and meta.exists():return read(meta)
    from scripts.run_main.run_dipolar_volume_export import NativeSampler
    p=load('periodic_k00_EH.npz');model=client.load(str(PERIODIC_MODEL))
    try:
        sampler=NativeSampler(model);mode=int(p['mode_assignment'][1]);phase=p['phase_to_raw_reference'][1]
        old=sampler(mode,p['rho_xyz_um'],[f'ewfd.E{axis}/(1[V/m])' for axis in 'xyz'])*phase
        err=float(np.linalg.norm(old-p['E_V_m'][1])/np.linalg.norm(p['E_V_m'][1]));assert err<1e-7
        a=cfg['a'];a1=np.array([a,0]);a2=np.array([a/2,np.sqrt(3)*a/2]);u=np.linspace(-.5,.5,65);z=np.linspace(0,.1,9)
        uu,vv,zz=np.meshgrid(u,u,z,indexing='ij');xy=uu.ravel()[:,None]*a1+vv.ravel()[:,None]*a2
        shifts=np.array([i*a1+j*a2 for i in range(-2,3) for j in range(-2,3)]);local=xy-shifts[cKDTree(shifts).query(xy)[1]]
        xyz=np.column_stack([local,zz.ravel()]);xyz[:,2]=np.clip(xyz[:,2],1e-7,.1-1e-7)
        raw=sampler(mode,xyz,[f'ewfd.E{axis}/(1[V/m])' for axis in 'xyz'])*phase
        rms=float(np.sqrt(np.sum(abs(p['E_V_m'][1])**2*p['weights_half_m3'][:,None])/np.sum(p['weights_half_m3'])))
        raw=(raw/rms).reshape(65,65,9,3)
        # Match periodic endpoints; only the source interpolant is modified, never reference fields.
        raw[0]=raw[-1]=(raw[0]+raw[-1])/2;raw[:,0]=raw[:,-1]=(raw[:,0]+raw[:,-1])/2
        value=raw.reshape(-1,3);table=np.column_stack([uu.ravel(),vv.ravel(),zz.ravel(),value.real,value.imag]);np.savetxt(path,table,fmt='%.16e',header='u v z_um Ex_re Ey_re Ez_re Ex_im Ey_im Ez_im',comments='% ')
        info=dict(path=str(path),sha256=digest(path),shape=list(raw.shape),reference_reproduction_error=err,reference_mode_index=mode,reference_phase=[float(phase.real),float(phase.imag)],rms_V_m=rms,interpolation='periodic primitive coordinates, linear 3D; sampled interface endpoints 0.0001nm inside slab',source_normalization='Common saved full-cavity quadrature set to 1 A^2/m; identical amplitude retained in B/A/C. Actual COMSOL norm and source power checked independently after each solve.',source_is_not_exact_background=True)
        put(meta,info);return info
    finally:client.remove(model)

def attach_source(model,out,cfg,run):
    from jpype.types import JArray,JString
    j=model.java;j.param().set('s5freq',f'{FREQ:.17g}[Hz]');j.param().set('s5a',f"{cfg['a']:.17g}[um]");j.param().set('s5R',f"{run['source_radius_um']:.17g}[um]");j.param().set('s5Jamp','1[A/m^2]')
    func=j.func().create('s5Gamma','Interpolation');func.set('source','file');func.set('filename',str(out/'01_results/source_template.txt'));func.set('nargs','3')
    names=['s5Exr','s5Eyr','s5Ezr','s5Exi','s5Eyi','s5Ezi'];func.set('funcs',JArray(JString,2)([[name,str(i+1)] for i,name in enumerate(names)]));func.set('argunit',['1','1','1']);func.set('fununit',['1']*6);func.set('interp','linear');func.set('extrap','const');func.importData()
    var=j.component('comp1').variable().create('s5source');var.selection().geom('geom1',3);var.selection().all();var.set('s5u0','(x-y/sqrt(3))/s5a');var.set('s5v0','2*y/(sqrt(3)*s5a)');var.set('s5u','s5u0-floor(s5u0+0.5)');var.set('s5v','s5v0-floor(s5v0+0.5)');var.set('s5r','sqrt(x^2+y^2)');var.set('s5window','if(s5r<s5R,cos(pi*s5r/(2*s5R))^2,0)*if(abs(z)<0.1[um],cos(pi*z/(0.2[um]))^2,0)')
    for axis in 'xyz':var.set('s5J'+axis,f'if(s5r<s5R,s5Jamp*s5window*(s5E{axis}r(s5u,s5v,z/1[um])+i*s5E{axis}i(s5u,s5v,z/1[um])),0[A/m^2])')
    node=j.component('comp1').physics('ewfd').create('s5current','ExternalCurrentDensity',3);node.selection().named('sel_layer_0_extent');node.set('Je',JArray(JString,1)(['s5Jx','s5Jy','s5Jz']))

def integration(model,expressions,tag='s5int'):
    from jpype.types import JInt
    j=model.java
    if tag in list(j.result().numerical().tags()):j.result().numerical().remove(tag)
    n=j.result().numerical().create(tag,'IntVolume');n.set('data','dset1');n.selection().named('sel_layer_0_extent');n.set('intorder',JInt(6));n.set('expr',expressions);n.set('unit',['1']*len(expressions))
    return np.asarray(n.getReal()).ravel()+1j*np.asarray(n.getImag()).ravel()

def configure_study(model):
    j=model.java
    for collection in [j.result().numerical(),j.result().dataset(),j.sol(),j.study()]:
        for tag in list(collection.tags()):collection.remove(tag)
    j.study().create('std1');j.study('std1').create('freq','Frequency');j.study('std1').feature('freq').set('plist','s5freq');j.study('std1').createAutoSequences('all')
    sol=j.sol('sol1');features={str(t):str(sol.feature(t).getType()) for t in sol.feature().tags()}
    for tag,kind in features.items():
        if kind=='Stationary':
            st=sol.feature(tag)
            if 'd1' not in list(st.feature().tags()):st.create('d1','Direct')
            st.feature('d1').set('linsolver','pardiso')
            for child in st.feature().tags():
                if str(st.feature(child).getType())=='FullyCoupled':st.feature(child).set('linsolver','d1')
    if 'dset1' not in list(j.result().dataset().tags()):j.result().dataset().create('dset1','Solution')
    j.result().dataset('dset1').set('solution','sol1')
    return features

def canonical(poly):return tuple(sorted(map(tuple,np.round(np.asarray(poly),8))))

def verify_saved_geometry(model,out):
    from comsol_workflow.simulation_config_geometry import load_simulation_config_geometry
    expected=load_simulation_config_geometry(out/'99_config/A_simulation_config.json');wp=model.java.component('comp1').geom('geom1').feature('wp1').geom();polys=[]
    for tag in wp.feature().tags():
        node=wp.feature(tag)
        if str(node.getType())=='Polygon' and str(node.label()).startswith('Hole'):
            table=node.getStringMatrix('table');poly=np.array([[float(str(x).split('[')[0]) for x in row] for row in table]);polys.append(poly)
    if not polys:raise RuntimeError('Cannot verify saved model hole polygons')
    assert {canonical(p) for p in polys}=={canonical(p) for p in expected.holes},'Historical geometry reconstruction differs from saved MPH'
    j=model.java;assert abs(j.param().evaluate('H_layer_0')-200e-9)<1e-15
    put(out/'80_logs/saved_geometry_check.json',dict(matched=True,holes=len(polys),slab_height_m=float(j.param().evaluate('H_layer_0')),input_model_unchanged=True))

def export_model(model,out,label,cfg,run):
    from scripts.run_main.run_dipolar_volume_export import NativeSampler
    from matplotlib.path import Path as PolygonPath
    from comsol_workflow.simulation_utils import electric_reflection_matrix
    sampler=NativeSampler(model)
    actual_frequency=float(sampler(0,np.array([[.4,.01,.05]]),['freq/1[Hz]'])[0,0].real);np.testing.assert_allclose(actual_frequency,FREQ,rtol=1e-12)
    expr=[f'ewfd.E{x}/(1[V/m])' for x in 'xyz']+[f'ewfd.H{x}/(1[A/m])' for x in 'xyz'];sx=np.array([-1,1,1,1,-1,-1]);sy=np.array([-1,1,-1,1,-1,1])
    old=load('finite_cavity_EH_o3.npz');rho=old['rho_xyz_um'];centers=old['centers_um'];idx=np.flatnonzero((centers[:,0]>=-1e-10)&(centers[:,1]>=-1e-10));tree=cKDTree(centers);rtree=cKDTree(rho);values=np.zeros((len(centers),len(rho),6),complex);seen=np.zeros(len(centers),bool)
    for start in range(0,len(idx),8):
        batch=idx[start:start+8];xyz=rho[None]+np.pad(centers[batch],((0,0),(0,1)))[:,None];raw=sampler(0,np.abs(xyz).reshape(-1,3),expr).reshape(len(batch),len(rho),6);raw*=np.where(xyz[...,0,None]<0,sx,1);raw*=np.where(xyz[...,1,None]<0,sy,1)
        for fx,fy in [(1,1),(-1,1),(1,-1),(-1,-1)]:
            distance,mapping=rtree.query(rho*np.array([fx,fy,1]));assert distance.max()<1e-10
            distance,target=tree.query(centers[batch]*[fx,fy]);assert distance.max()<1e-10
            signs=(sx if fx<0 else np.ones(6))*(sy if fy<0 else np.ones(6));values[target]=raw[:,mapping]*signs;seen[target]=True
    assert seen.all();np.savez_compressed(out/'01_results'/f'{label}_cavity_EH.npz',E_V_m=values[...,:3],H_A_m=values[...,3:],rho_xyz_um=rho,centers_um=centers,weights_half_m3=old['weights_half_m3'])
    axis=np.linspace(-cfg['finite_Lx']/2,cfg['finite_Lx']/2,801);axis[len(axis)//2]=0;xx,yy=np.meshgrid(axis,axis);mask=PolygonPath(np.asarray(cfg['finite_boundary'])).contains_points(np.column_stack([xx.ravel(),yy.ravel()]),radius=-1e-8).reshape(xx.shape);refl=np.maximum(np.arange(len(axis)),len(axis)-1-np.arange(len(axis)))
    for z,tag in [(0.,'slab'),(1.1,'air'),(1.4,'air_check')]:
        quarter=mask&(xx>=0)&(yy>=0);xyz=np.column_stack([xx[quarter],yy[quarter],np.full(quarter.sum(),z)]);eh=np.zeros((*xx.shape,6),complex);eh[quarter]=sampler(0,xyz,expr);eh=eh[refl[:,None],refl[None,:]];eh*=np.where(xx[...,None]<0,sx,1);eh*=np.where(yy[...,None]<0,sy,1);eh[~mask]=0
        # Quarter-axis parities are exact model constraints, not independent checks.
        eh[:,len(axis)//2,:]*=(1+sx)/2;eh[len(axis)//2,:,:]*=(1+sy)/2
        np.savez_compressed(out/'01_results'/f'{label}_{tag}_EH.npz',E_V_m=eh[...,:3],H_A_m=eh[...,3:],x_um=axis,y_um=axis,z_um=z,valid=mask,frequency_Hz=FREQ)
    # Full reflected volume: 8 * (1/4) Re(E dot D* + H dot B*) for the fixed nondispersive materials.
    ints=integration(model,['8*(abs(comp1.s5Jx)^2+abs(comp1.s5Jy)^2+abs(comp1.s5Jz)^2)/(1[A^2/m^4])/(1[m^3])','-4*real(comp1.s5Jx*conj(ewfd.Ex)+comp1.s5Jy*conj(ewfd.Ey)+comp1.s5Jz*conj(ewfd.Ez))/(1[W/m^3])/(1[m^3])','2*real(ewfd.Ex*conj(ewfd.Dx)+ewfd.Hx*conj(ewfd.Bx)+ewfd.Ey*conj(ewfd.Dy)+ewfd.Hy*conj(ewfd.By)+ewfd.Ez*conj(ewfd.Dz)+ewfd.Hz*conj(ewfd.Bz))/(1[J/m^3])/(1[m^3])'])
    return dict(source_norm_A2_per_m=float(ints[0].real),source_delivered_power_W=float(ints[1].real),slab_energy_J=float(ints[2].real),source_power_definition='-0.5 Re integral J_ext dot E* dV; factor 8 full reflected volume')

def run(out):
    from scripts.run_main.run_boundary_analysis import runtime_preflight
    cfg=read(out/'99_config/historical_config.json');settings=read(out/'99_config/geometry_and_run_config.json')
    for path,identity in settings['inputs'].items():
        source=Path(path)
        if source.stat().st_size!=identity['bytes'] or digest(source)!=identity['sha256']:raise RuntimeError('Input identity changed: '+path)
    put(out/'80_logs/runtime_preflight.json',runtime_preflight())
    import mph
    from jpype import JClass
    client=mph.start(version='6.3',cores=settings['cores']);util=JClass('com.comsol.model.util.ModelUtil')
    if not util.checkoutLicense('COMSOL','WAVEOPTICS'):raise RuntimeError('COMSOL/WAVEOPTICS license checkout failed')
    progress(out,'COMSOL started; licenses checked out',status='running',COMSOL_started=True)
    try:
        progress(out,'Exporting periodic Gamma electric source template');source_table(client,out,cfg,settings)
        # Verify A before solving controls: cannot silently solve a different historical geometry.
        scale_path=out/'99_config/source_amplitude.json'
        if not (out/'00_model/A_prepared.mph').exists():
            progress(out,'Loading saved mode19 model for historical geometry verification');model=client.load(str(SOURCE_MODEL));verify_saved_geometry(model,out)
            attach_source(model,out,cfg,settings)
            if not scale_path.exists():
                from scipy.interpolate import RegularGridInterpolator
                table=np.loadtxt(out/'01_results/source_template.txt',comments='%');data=table[:,3:6]+1j*table[:,6:9]
                interp=RegularGridInterpolator((np.linspace(-.5,.5,65),np.linspace(-.5,.5,65),np.linspace(0,.1,9)),data.reshape(65,65,9,3),bounds_error=True)
                p=load('periodic_k00_EH.npz');f=load('finite_cavity_EH_o3.npz');rho=p['rho_xyz_um'];u=(rho[:,0]-rho[:,1]/np.sqrt(3))/cfg['a'];v=2*rho[:,1]/(np.sqrt(3)*cfg['a']);uvz=np.column_stack([u-np.floor(u+.5),v-np.floor(v+.5),rho[:,2]])
                template=interp(uvz);r=np.linalg.norm(f['centers_um'][:,None,:]+rho[None,:,:2],axis=-1);window=np.where(r<settings['source_radius_um'],np.cos(np.pi*r/(2*settings['source_radius_um']))**2,0)*np.cos(np.pi*rho[None,:,2]/.2)**2
                norm=float(np.sum(window**2*(2*p['weights_half_m3']*np.sum(abs(template)**2,-1))[None]))
                put(scale_path,dict(Jamp_A_m2=float(np.sqrt(settings['source_norm_target_A2_per_m']/norm)),unscaled_norm_A2_per_m=norm,reference='common saved full-cavity quadrature; actual COMSOL source norm checked after each solve',quadrature_norm_target_A2_per_m=1.))
            configure_study(model)
            model.save(str(out/'00_model/A_prepared.mph'))
            client.remove(model)
        for label in MODELS:
            completed=out/'80_logs'/f'{label}_completed.json'
            if completed.exists():progress(out,f'Reuse completed {label}',model=label);continue
            checkpoint=out/'00_model'/f'{label}_solved.mph';resume=checkpoint.exists();mesh_checkpoint=out/'00_model'/f'{label}_mesh_and_source.mph';mesh_resume=mesh_checkpoint.exists()
            progress(out,f'{label}: loading solved checkpoint' if resume else f'{label}: building driven model',model=label,phase='building')
            holder=None
            if resume:model=client.load(str(checkpoint))
            elif mesh_resume:model=client.load(str(mesh_checkpoint))
            elif label=='A':model=client.load(str(out/'00_model/A_prepared.mph'))
            else:
                from comsol_workflow.simulation_utils import SimulationRun,SimulationConfig,LayerSpec,BoundarySpec
                from comsol_workflow.simulation_config_geometry import load_simulation_config_geometry
                from comsol_workflow.mesh_constrcution import construct_manual_mesh
                geo=load_simulation_config_geometry(out/'99_config'/f'{label}_simulation_config.json');layers=[LayerSpec(height=cfg['simulation_common']['slab_height'],holes=list(geo.holes),refractive_index=cfg['simulation_common']['refractive_index'])]
                holder=SimulationRun(config=SimulationConfig(wavelength=cfg['simulation_common']['wavelength'],mesh_auto_size=5,mode_type='TE',symmetry_x_boundary='PEC',symmetry_y_boundary='PMC'));model=holder.model
                holder.build_geometry(BoundarySpec('quarter_hexagon_boundary',a=cfg['finite_hex_a']),layers,'finite_quarter');construct_manual_mesh(model.java,layers,mesh_auto_size=5)
                progress(out,f'{label}: meshing at mesh=5',phase='meshing');util.showProgress(str(out/'80_logs'/f'{label}_comsol_progress.log'));model.java.component('comp1').mesh('mesh1').run()
            try:
                j=model.java
                if not resume:
                    if label!='A' and not mesh_resume:attach_source(model,out,cfg,settings)
                    features=configure_study(model)
                    amp=read(scale_path)['Jamp_A_m2'];j.param().set('s5Jamp',f'{amp:.17g}[A/m^2]')
                    model.save(str(out/'00_model'/f'{label}_mesh_and_source.mph'))
                    put(out/'80_logs'/f'{label}_solver.json',dict(sequence=features,source_amplitude_A_m2=amp,frequency_Hz=FREQ,mesh_vertices=int(np.asarray(j.component('comp1').mesh('mesh1').getVertex()).shape[1])))
                    progress(out,f'{label}: solving same-source frequency-domain response',phase='solving');util.showProgress(str(out/'80_logs'/f'{label}_comsol_progress.log'));j.sol('sol1').runAll();model.save(str(checkpoint))
                progress(out,f'{label}: exporting common-coordinate complex E/H fields',phase='exporting');metrics=export_model(model,out,label,cfg,settings);metrics.update(model=label,status='complete',checkpoint=str(checkpoint));put(completed,metrics);progress(out,f'{label}: solve and exports complete',phase='complete')
            finally:
                if holder is not None:holder.clear()
                else:client.remove(model)
        progress(out,'Three driven models complete; offline comparison starting',phase='analysis')
        from scripts.analysis.analyze_dipolar_cladding_scattering import driven_comparison
        driven_comparison(out)
        progress(out,'Initial three-model experiment complete',status='complete',phase='complete')
    finally:
        util.showProgress(False);client.clear();
        if getattr(client,'server',None):client.server.stop()

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--prepare',action='store_true');args=p.parse_args();out=args.output.resolve()
    root=Path(__file__).resolve().parents[2]
    allowed=[root/'results/S5_dipolarSingularity_analysis',root/'scripts/.out/finite_cavity']
    if not any(out.is_relative_to(p.resolve()) and out!=p.resolve() for p in allowed) or not out.is_dir():raise ValueError('Output must be a pre-created run directory under results/S5_dipolarSingularity_analysis (legacy scripts/.out/finite_cavity also supported)')
    try:
        if args.prepare:prepare(out)
        else:run(out)
    except Exception as exc:
        progress(out,str(exc),status='failed',exception_type=type(exc).__name__);traceback.print_exc();raise

if __name__=='__main__':main()
