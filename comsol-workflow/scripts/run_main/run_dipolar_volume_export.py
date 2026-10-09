"""Fixed-mesh S5 volume fields; existing finite solution is read, never solved."""
from __future__ import annotations

from comsol_workflow.output_paths import output_path

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial import cKDTree

from comsol_workflow.boundary_integrals import _vertical_span, polygon_area
from scripts.analysis.visualize_dipolar_broadening_2d import SERIES, OUTPUT, SOURCE
from comsol_workflow.dipolar_lattice_analysis import load_data
from comsol_workflow.result_io import digest, read_json as json_read, write_json
from scripts.run_main.run_boundary_analysis import runtime_preflight
from scripts.run_main.run_dipolar_periodic_reference import S4, OUT
from comsol_workflow.dipolar_radiation import mode_assignment

from comsol_workflow.dipolar_inputs import load_dipolar_inputs

INPUTS = load_dipolar_inputs()
FINITE = INPUTS.finite_model
CONFIG = INPUTS.source_config_dir / 'config.json'


def cell_rule(order=3):
    """One common rule resolving BOTH cavity and cladding material interfaces."""
    from comsol_workflow.geometry_utils import create_hexagon_design
    from scripts.run_main.run_band_pair import get_hole_params
    config = json_read(CONFIG)
    geometry = json_read(S4/'99_config/realized_geometry.json')
    outer = np.asarray(geometry['outer'])
    cavity = np.asarray(geometry['holes'])
    _, cladding, _ = create_hexagon_design(config['a'], get_hole_params(config['cladding_params']))
    holes = [*cavity, *cladding]
    # Rounded breakpoints remove 1e-16 rotations of nominally identical vertices.
    breaks = np.unique(np.round(np.concatenate([outer, *holes])[:, 0], 13))
    nodes, weights = leggauss(order)
    points, areas, masks = [], [], []
    for left, right in zip(breaks[:-1], breaks[1:]):
        for x, wx in zip((left+right)/2+(right-left)/2*nodes, (right-left)/2*weights):
            lo, hi = _vertical_span(outer, x)
            spans = [_vertical_span(h, x) for h in holes]
            ends = np.unique(np.round([lo, hi, *[v for s in spans if s for v in s]], 13))
            for low, high in zip(ends[:-1], ends[1:]):
                if high-low < 1e-12:
                    continue
                midpoint = (low+high)/2
                air = [any(s and s[0] < midpoint < s[1] for s in spans[i:i+6]) for i in [0,6]]
                for y, wy in zip(midpoint+(high-low)/2*nodes, (high-low)/2*weights):
                    points.append([x,y]); areas.append(wx*wy); masks.append(np.logical_not(air))
    xy, area, material = np.asarray(points), np.asarray(areas), np.asarray(masks).T
    np.testing.assert_allclose(area.sum(), abs(polygon_area(outer)), rtol=1e-10)
    for i, hs in enumerate([cavity,cladding]):
        expected = abs(polygon_area(outer))-sum(abs(polygon_area(h)) for h in hs)
        np.testing.assert_allclose(area@material[i], expected, rtol=1e-10)
    z = (nodes+1)*.05  # upper half of H=0.2 um, TE mirror at z=0
    xyz = np.column_stack([np.repeat(xy,len(z),axis=0),np.tile(z,len(xy))])
    w = np.repeat(area,len(z))*np.tile(weights*.05,len(xy))*1e-18
    return xyz, w, np.repeat(material,len(z),axis=1)


class NativeSampler:
    """COMSOL arbitrary-point interpolation with explicit SI expressions."""
    def __init__(self, model):
        from jpype.types import JArray,JDouble,JInt,JString
        self.JArray,self.JDouble,self.JInt,self.JString=JArray,JDouble,JInt,JString
        j=model.java
        if 's5volume' in list(j.result().numerical().tags()):
            j.result().numerical().remove('s5volume')
        self.node=j.result().numerical().create('s5volume','Interp')
        for key,value in [('data','dset1'),('ext',0.),('matherr','on'),('coorderr','on'),('recover','off')]:
            self.node.set(key,value)
        self.node.selection().all()

    def __call__(self,mode,xyz_um,expressions):
        self.node.set('expr',self.JArray(self.JString,1)(expressions))
        self.node.set('solnum',self.JArray(self.JInt,1)([mode+1]))
        result=[]
        for points in np.array_split(xyz_um,max(1,int(np.ceil(len(xyz_um)/20000)))):
            self.node.setInterpolationCoordinates(self.JArray(self.JDouble,2)(points.T.tolist()))
            re=np.asarray(self.node.getData(),float)
            im=np.asarray(self.node.getImagData(),float) if self.node.isComplex() else np.zeros_like(re)
            value=(re+1j*im).reshape(len(expressions),-1).T
            if value.shape!=(len(points),len(expressions)) or not np.isfinite(value).all():
                raise ValueError('Nonfinite/missing volume interpolation')
            result.append(value)
        return np.concatenate(result)


E_EXPR=[f'ewfd.E{c}/(1[V/m])' for c in 'xyz']


def finite_export(client, order=3):
    from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies
    from comsol_workflow.simulation_utils import electric_reflection_matrix
    from comsol_workflow.farfield_fft import load_air_field
    path=output_path(OUTPUT, f's5_finite_volume_o{order}.npz')
    if path.exists():
        print('Reuse finite volume export',flush=True); return
    config=json_read(CONFIG)
    if config['simulation_common']['mesh_auto_size']!=5 or config['symmetry']!={'id':1,'x_boundary':'PEC','y_boundary':'PMC'}:
        raise ValueError('Unexpected finite model configuration')
    before=FINITE.stat()
    print('Loading saved finite mode19 (no eigensolve)',flush=True)
    model=client.load(str(FINITE)); j=model.java
    try:
        print('Loaded finite; material selection: sel_layer_0_mat_dom',flush=True)
        print('Parameters: '+str(list(j.param().varnames())),flush=True)
        modes=frequencies(model)
        np.testing.assert_allclose(modes.iloc[19][['re','im']],[197.9530009044521,.0128930496899445],rtol=1e-9)
        omega=2*np.pi*complex(modes.iloc[19]['re'],modes.iloc[19]['im'])*1e12
        sampler=NativeSampler(model)
        print('Probe:',sampler(19,np.array([[.4,.01,.05]]),E_EXPR),flush=True)
        # The saved air field is a fully independent source and fixes its phase convention.
        emeta=json_read(output_path(OUTPUT, 'electric_19_summary.json'))
        air=load_air_field(Path(emeta['source']))
        ix=np.array([409,431,455]);iy=np.array([407,429,451])
        xyz=np.column_stack([air['x_axis'][ix],air['y_axis'][iy],np.full(3,emeta['plane_z_um'])])
        check=sampler(19,xyz,E_EXPR)[:,:2]
        old=np.column_stack([air['Ex'][iy,ix],air['Ey'][iy,ix]])
        native_error=float(np.linalg.norm(check-old)/np.linalg.norm(old))
        conjugate_error=float(np.linalg.norm(check.conj()-old)/np.linalg.norm(old))
        print(f'Air convention comparison native={native_error:.6g}, conjugate={conjugate_error:.6g}',flush=True)
        if min(native_error,conjugate_error)>1e-6:
            raise ValueError('Saved air field does not match the source eigenmode')
        # Convert old export to native before all direct comparisons.
        direct=np.array([complex(*emeta['C0_full_saved_air_V_m'][c]) for c in ['Ex','Ey']]) if 'C0_full_saved_air_V_m' in emeta else None
        if direct is None:
            from scripts.analysis.validate_dipolar_bloch_radiation import finite_air_check
            direct,_=finite_air_check(Path(emeta['source']))
        if conjugate_error<native_error:
            direct=direct.conj()
        # Native FEM integration covers the whole dielectric, including outer pad.
        node=j.result().numerical().create('s5fullintegral','IntVolume')
        node.set('data','dset1');node.set('solnum',sampler.JArray(sampler.JInt,1)([20]));node.selection().named('sel_layer_0_mat_dom')
        kval=omega/299792458.
        ks=f'({kval.real:.17g}+i*{kval.imag:.17g})[1/m]'
        expressions=[f'8*9.89*ewfd.Ey*cos({ks}*z)/1[V*m^2]',
                     '9.89/1[m^3]']
        node.set('expr',expressions);node.set('unit',['1','1']);node.set('intorder',sampler.JInt(6))
        integral=np.asarray(node.getReal()).ravel()+1j*np.asarray(node.getImag()).ravel()
        full=np.array([0j,integral[0]*kval**2/(4*np.pi)])
        print('Complete material integral F (V): '+str(full),flush=True)
        rho,w,material=cell_rule(order)
        tree=cKDTree(rho)
        flips=[np.array([sx,sy,1]) for sx,sy in [(1,1),(-1,1),(1,-1),(-1,-1)]]
        mappings=[]
        for flip in flips:
            distance,index=tree.query(rho*flip)
            if distance.max()>1e-10:
                raise ValueError('Quadrature is not mirror closed')
            mappings.append(index)
        ex=np.diag(electric_reflection_matrix('y','PMC'))
        ey=np.diag(electric_reflection_matrix('x','PEC'))
        signs=[np.ones(3),ex,ey,ex*ey]
        arrays={}
        for region,key in enumerate(['bulk_points','cladding_points']):
            centers=np.array([[p['x'],p['y']] for p in config[key]])
            center_tree=cKDTree(centers)
            ids=np.where((centers[:,0]>=-1e-10)&(centers[:,1]>=-1e-10))[0]
            values=np.empty((len(centers),len(rho),3),complex)
            seen=np.zeros(len(centers),bool)
            for step,start in enumerate(range(0,len(ids),8)):
                batch=ids[start:start+8]
                points=rho[None,:,:]+np.pad(centers[batch],((0,0),(0,1)))[:,None,:]
                if start==0:
                    domains=sampler(19,np.abs(points[0]),['dom'])[:,0].real.astype(int)
                    material_domains=np.asarray(j.component('comp1').selection('sel_layer_0_mat_dom').entities(3),int)
                    if not np.array_equal(np.isin(domains,material_domains),material[region]):
                        raise ValueError('Recorded geometry differs from actual COMSOL material at quadrature nodes')
                    print(f'{key}: actual material mask verified',flush=True)
                raw=sampler(19,np.abs(points).reshape(-1,3),E_EXPR).reshape(points.shape)
                raw*=np.where(points[:,:,0,None]<0,ex,1)
                raw*=np.where(points[:,:,1,None]<0,ey,1)
                for flip,mapping,sign in zip(flips,mappings,signs):
                    distance,target=center_tree.query(centers[batch]*flip[:2])
                    if distance.max()>1e-10: raise ValueError('Finite cells are not mirror closed')
                    values[target]=raw[:,mapping,:]*sign
                    seen[target]=True
                if step%10==0:
                    print(f'{key}: {min(start+8,len(ids))}/{len(ids)} quarter cells; {len(rho)} nodes/cell',flush=True)
            if not seen.all(): raise ValueError('Missing mirrored finite cells')
            arrays[f'E_{region}']=values
            arrays[f'centers_{region}_um']=centers
        after=FINITE.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise RuntimeError('Source finite MPH changed')
        metadata={'finite_model':str(FINITE),'finite_model_size':before.st_size,
                  'finite_model_mtime_ns':before.st_mtime_ns,'config_sha256':digest(CONFIG),
                  'mesh':5,'mode':19,'quadrature_order':order,'native_air_error':native_error,
                  'conjugate_air_error':conjugate_error,'time_convention':'native_COMSOL_exp(+iwt)',
                  'upper_half_slab_TE_mirror':True,'exported_regions':'actual finite fields in cavity and cladding cells; basis support is set by the analysis',
                  'complete_integral':'native FEM over all material, x/y and z mirror reconstruction'}
        np.savez_compressed(path,**arrays,rho_xyz_um=rho,weights_half_m3=w,material_masks=material,
                            omega_rad_s=omega,full_volume_F_V=full,direct_air_C_V_m=direct,
                            z_air_m=emeta['plane_z_um']*1e-6,metadata_json=json.dumps(metadata))
        write_json(output_path(OUTPUT, f's5_finite_volume_o{order}_audit.json'),metadata)
        print('Saved '+str(path),flush=True)
    finally:
        client.remove(model)


def periodic_export(client,order=3):
    from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies, mesh_identity
    d=load_data(SOURCE);rho,w,material=cell_rule(order)
    model=client.load(str(S4/'00_model/s4_gamma.mph'));j=model.java
    try:
        mesh=mesh_identity(model)
        if mesh['auto_size']!=5: raise ValueError('Mesh must remain 5')
        for index,k in enumerate(d['k19']):
            path=OUT/f'01_results/k{index:02d}_volume_o{order}.npz'
            if path.exists():
                print(f'k{index:02d}: reuse volume fields',flush=True);continue
            anchor=np.load(OUT/f'01_results/k{index:02d}_reference.npz')
            if index:
                for name,value in zip(['kx','ky'],k):j.param().set(name,f'{value:.17g}[1/um]')
                print(f'k{index:02d}: same-mesh solve for volume export',flush=True)
                j.sol('sol1').clearSolutionData();j.sol('sol1').runAll()
            if 's4surface_freq' in list(j.result().numerical().tags()):j.result().numerical().remove('s4surface_freq')
            modes=frequencies(model)
            sampler=NativeSampler(model)
            center=np.column_stack([d['rho_points'],np.zeros(len(d['rho_points']))])
            hz=np.array([sampler(m,center,['ewfd.Hz/(1[A/m])'])[:,0] for m in range(4)])
            target=anchor['Hz_native'][anchor['gamma_band_to_mode']]
            assignment,overlap=mode_assignment(target,hz)
            if np.min(overlap[np.arange(4),assignment])<.99:
                raise ValueError('Repeated mesh5 modes differ from exact-k references')
            norm=np.sqrt(np.mean(abs(hz)**2,axis=1))
            electric=np.array([sampler(int(m),rho,E_EXPR)/norm[m] for m in assignment])
            if mesh_identity(model)!=mesh:raise ValueError('Mesh changed')
            np.savez_compressed(path,E_native_over_Hz_rms=electric,rho_xyz_um=rho,
                                k_um_inv=k,gamma_band_to_mode=assignment,
                                frequencies_THz=modes[['re','im','q']].to_numpy()[assignment],
                                overlap=overlap,weights_half_m3=w,material_masks=material)
            print(f'k{index:02d}: volume saved',flush=True)
    finally:
        client.remove(model)


def run(stage,order):
    preflight=runtime_preflight()
    import mph
    from jpype import JClass
    client=mph.start(version='6.3',cores=4)
    write_json(output_path(OUTPUT, 's5_volume_preflight.json'),preflight)
    if not JClass('com.comsol.model.util.ModelUtil').checkoutLicense('COMSOL','WAVEOPTICS'):
        raise RuntimeError('License checkout failed')
    try:
        if stage=='finite':finite_export(client,order)
        else:periodic_export(client,order)
    finally:
        client.clear()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['finite','periodic'])
    p.add_argument('--order',type=int,default=3)
    a=p.parse_args();run(a.stage,a.order)
