"""Export a new cavity plane from the saved A solution and reuse k analysis."""
from __future__ import annotations
import argparse,hashlib,os,time,traceback
import numpy as np
from scripts.analysis.prepare_dipolar_gamma_origin import paths, read, write, sha
from comsol_workflow.dipolar_decomposition import shared_node_projection
from scripts.run_main.export_dipolar_gamma_boundary import reflected_sample


def run(height):
    if not np.isfinite(height) or height<=0:raise ValueError('Expected positive plane height')
    base=paths();cfg=read(base['99_config']/'replay_config.json')
    tag=f"{height:g}".replace('.','p');name=base['99_config'].name+f'_z{tag}um_20261009'
    target={key:value.parent/name for key,value in base.items()}
    for path in target.values():path.mkdir(parents=True,exist_ok=True)
    def progress(message,**extra):
        state=dict(message=message,python_pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),height_um=height,**extra)
        write(target['80_logs']/'height_execution.json',state);print(message,flush=True)
    cfg.update(run_id=name,parent_run=read(base['99_config']/'replay_config.json')['run_id'],height_request='User requested same three figures at z=1.5 um on 2026-10-09',no_new_solve=True)
    cfg['surface']['cap_z_um']=height
    config=target['99_config']/'replay_config.json'
    if config.exists():assert read(config)==cfg,'Configuration mismatch'
    else:write(config,cfg)
    output=target['01_results']/'cap_boundary_inputs.npz'
    try:
        if not output.exists():
            from scripts.run_main.run_boundary_analysis import runtime_preflight
            from scripts.run_main.run_dipolar_volume_export import NativeSampler
            write(target['80_logs']/'runtime_preflight.json',runtime_preflight())
            model_hash=sha(cfg['baseline']['model'])
            with np.load(base['01_results']/'cap_boundary_inputs.npz') as z:
                xy=z['xy_um'];local=z['local_xy_um'];triangles=z['local_triangles'];weights=z['local_weights_um2'];mapping=z['cell_to_global_node']
                probe_index=np.linspace(0,len(xy)-1,9).astype(int);expected=z['FULL_EH'][probe_index]
            progress('Starting COMSOL to read saved A; no mesh or solve',status='running')
            import mph
            from jpype import JClass
            client=mph.start(version='6.3',cores=16);util=JClass('com.comsol.model.util.ModelUtil')
            try:
                if not util.checkoutLicense('COMSOL','WAVEOPTICS'):raise RuntimeError('License checkout failed')
                progress('Loading saved A_solved.mph',status='running')
                model=client.load(cfg['baseline']['model']);j=model.java;sampler=NativeSampler(model)
                assert str(j.result().dataset('dset1').getString('solution'))=='sol1'
                assert str(j.study('std1').feature('freq').getType())=='Frequency'
                slab=float(j.param().evaluate('H_layer_0','um'));air=float(j.param().evaluate('H_air','um'))
                assert slab/2<height<slab/2+air,'Plane not in untransformed air interval'
                with np.load(base['99_config']/'cell_mapping.npz') as z:centers=z['centers_um']
                probes=np.unique(np.column_stack([abs(centers),np.full(len(centers),height)]),axis=0)
                probes[:,:2]=np.maximum(probes[:,:2],1e-7)
                identity=sampler(0,probes,['freq/1[Hz]','dom','ewfd.epsilonrxx'])
                np.testing.assert_allclose(identity[:,0],cfg['frequency_Hz'],rtol=1e-12)
                pml=list(map(int,j.component('comp1').selection('sel_pml_dom').entities(3)))
                assert not np.isin(np.rint(identity[:,1].real).astype(int),pml).any()
                np.testing.assert_allclose(identity[:,2],1,atol=1e-9)
                old=reflected_sample(sampler,np.column_stack([xy[probe_index],np.full(len(probe_index),.3)]))
                error=float(np.linalg.norm(old-expected)/np.linalg.norm(expected));assert error<1e-7,error
                write(target['80_logs']/'loaded_identity.json',dict(model_sha256=model_hash,frequency_Hz=cfg['frequency_Hz'],slab_height_um=slab,air_height_um=air,plane_height_um=height,air_outside_PML=True,PML_domains=pml,previous_plane_probe_relative_error=error,dataset='dset1',solution_number=1))
                def checkpoint_sampler(mode,points,expressions):
                    cache=target['01_results']/'positive_cap_EH.npy';meta=target['80_logs']/'positive_export_progress.json'
                    identity=dict(source_sha256=model_hash,coordinates_sha256=hashlib.sha256(points.tobytes()).hexdigest(),expressions=expressions,shape=[len(points),6])
                    if cache.exists():
                        state=read(meta);assert state['identity']==identity
                        raw=np.load(cache,mmap_mode='r+');start=state['completed']
                    else:
                        raw=np.lib.format.open_memmap(cache,mode='w+',dtype=np.complex128,shape=(len(points),6));start=0
                        write(meta,dict(identity=identity,completed=0))
                    for lo in range(start,len(points),20000):
                        hi=min(lo+20000,len(points));raw[lo:hi]=sampler(mode,points[lo:hi],expressions);raw.flush()
                        write(meta,dict(identity=identity,completed=hi));progress(f'Plane export {hi}/{len(points)} symmetry-unique points',status='running',completed=hi,total=len(points))
                    return raw
                full=reflected_sample(checkpoint_sampler,np.column_stack([xy,np.full(len(xy),height)]))
            finally:
                client.clear()
                if getattr(client,'server',None):client.server.stop()
            progress('Constructing common-cell GAMMA and REST',status='running')
            counts=np.bincount(mapping.ravel(),minlength=len(xy));gamma=np.mean(full[mapping],axis=0)
            projected=shared_node_projection(np.tile(gamma,(len(mapping),1)),mapping.ravel(),counts)
            np.savez_compressed(output,xy_um=xy,FULL_EH=full,GAMMA_EH=projected,REST_EH=full-projected,raw_Gamma_local_EH=gamma,local_xy_um=local,local_triangles=triangles,local_weights_um2=weights,cell_to_global_node=mapping,cap_z_um=height)
            del full,projected,gamma
        progress('Analyzing complete k components and rendering Ex',status='running')
        from scripts.analysis.analyze_dipolar_gamma_origin import boundary_figure, combined_boundary_figure
        figures=[boundary_figure('Ex',target)]
        progress('Rendering Ey',status='running');figures.append(boundary_figure('Ey',target))
        progress('Rendering combined Ex/Ey',status='running');figures.append(combined_boundary_figure(target))
        progress('Plane export, decomposition and three figures completed',status='complete',figures=figures)
    except Exception as exc:
        progress(str(exc),status='failed',exception_type=type(exc).__name__);traceback.print_exc();raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height-um',type=float,required=True);p.add_argument('--execute',action='store_true');args=p.parse_args()
    if not args.execute:raise ValueError('Saved COMSOL field export requires --execute')
    run(args.height_um)
