"""Radiation quadrature on the existing mesh5, resolving both material patterns."""
from comsol_workflow.output_paths import output_path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial import ConvexHull
from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.boundary_integrals import polygon_area
from scripts.run_main.run_band_pair import get_hole_params
from scripts.run_main.run_dipolar_volume_export import NativeSampler, E_EXPR, CONFIG, S4, OUT, OUTPUT, SOURCE, load_data, json_read, mode_assignment, runtime_preflight, write_json


def tetra_rule(vertices,order=3):
    """Interior Duffy rule, exact for quadratic FEM E before the smooth phase."""
    nodes,weights=leggauss(order);nodes=(nodes+1)/2;weights=weights/2
    u,v,t=np.meshgrid(nodes,nodes,nodes,indexing='ij')
    wu,wv,wt=np.meshgrid(weights,weights,weights,indexing='ij')
    bary=np.stack([(1-u)*(1-v)*(1-t),u,(1-u)*v,(1-u)*(1-v)*t],axis=-1).reshape(-1,4)
    determinant=abs(np.linalg.det((vertices[1:]-vertices[0]).T))
    return bary@vertices,(wu*wv*wt*(1-u)**2*(1-v)).ravel()*determinant


def clipped_tetrahedra(vertices,triangle):
    """Intersect a tetrahedron with a vertical triangular air-hole prism."""
    points=np.asarray(vertices);triangle=np.asarray(triangle)
    if polygon_area(triangle)<0:triangle=triangle[::-1]
    for a,b in zip(triangle,np.roll(triangle,-1,axis=0)):
        n=np.array([-(b-a)[1],(b-a)[0],0.]);distance=(points-np.r_[a,0.])@n
        if np.min(distance)>=-1e-13:continue
        if np.max(distance)<=1e-13:return []
        keep=points[distance>=0].tolist()
        for i in range(len(points)):
            for j in range(i):
                if distance[i]*distance[j]<0:
                    keep.append((points[i]+distance[i]/(distance[i]-distance[j])*(points[j]-points[i])).tolist())
        points=np.unique(np.round(keep,14),axis=0)
        if len(points)<4:return []
        hull=ConvexHull(points)
        if hull.volume<1e-18:return []
        points=points[hull.vertices]
    if len(points)==4:return [points]
    hull=ConvexHull(points);center=points.mean(axis=0)
    return [np.vstack([center,points[face]]) for face in hull.simplices]


def mesh_rule(model):
    j=model.java;mesh=j.component('comp1').mesh('mesh1')
    vertices=np.asarray(mesh.getVertex(),float).T;elements=np.asarray(mesh.getElem('tet'),int).T
    entities=np.asarray(mesh.getElemEntity('tet'),int)
    slab=np.asarray(j.component('comp1').selection('sel_slab_layer').entities(3),int)
    silicon=np.asarray(j.component('comp1').selection('sel_slab_dom').entities(3),int)
    selected=np.flatnonzero(np.isin(entities,slab))
    if np.max(abs(vertices))>10 or np.max(abs(vertices))<.1:raise ValueError('Unexpected mesh-coordinate units')
    config=json_read(CONFIG);outer,holes,_=create_hexagon_design(config['a'],get_hole_params(config['cladding_params']))
    points=[];weights=[]
    for index in selected:
        tet=vertices[elements[index]];xyz,w=tetra_rule(tet)
        points.append(xyz);weights.append(np.vstack([w*np.isin(entities[index],silicon),w]))
        lower,upper=tet[:,:2].min(axis=0),tet[:,:2].max(axis=0)
        for hole in holes:
            if np.any(np.max(hole,axis=0)<lower) or np.any(np.min(hole,axis=0)>upper):continue
            for part in clipped_tetrahedra(tet,hole):
                xyz,w=tetra_rule(part)
                if w.sum()<1e-18:continue
                points.append(xyz);weights.append(np.vstack([np.zeros_like(w),-w]))
    xyz=np.concatenate(points);w=np.concatenate(weights,axis=1)*1e-18
    geometry=json_read(S4/'99_config/realized_geometry.json')
    for region,hs in enumerate([geometry['holes'],holes]):
        expected=(abs(polygon_area(outer))-sum(abs(polygon_area(h)) for h in hs))*.1e-18
        np.testing.assert_allclose(w[region].sum(),expected,rtol=1e-8,atol=0)
    print(f'Original mesh rule: {len(selected)} slab tets, {len(xyz)} interior points',flush=True)
    return xyz,w


def export(client):
    from scripts.run_sweep.run_boundary_surface_zeta_scan import mesh_identity
    d=load_data(SOURCE);model=client.load(str(S4/'00_model/s4_gamma.mph'));j=model.java
    try:
        mesh=mesh_identity(model)
        if mesh['auto_size']!=5:raise ValueError('Mesh must remain 5')
        rho,weights=mesh_rule(model)
        for index,k in enumerate(d['k19']):
            path=OUT/f'01_results/k{index:02d}_mesh_volume.npz'
            if path.exists():continue
            if index:
                for name,value in zip(['kx','ky'],k):j.param().set(name,f'{value:.17g}[1/um]')
                j.sol('sol1').clearSolutionData();j.sol('sol1').runAll()
            sampler=NativeSampler(model)
            center=np.column_stack([d['rho_points'],np.zeros(len(d['rho_points']))])
            hz=np.array([sampler(m,center,['ewfd.Hz/(1[A/m])'])[:,0] for m in range(4)])
            anchor=np.load(OUT/f'01_results/k{index:02d}_volume_o3.npz')
            old=np.load(OUT/f'01_results/k{index:02d}_reference.npz')
            assignment,overlap=mode_assignment(old['Hz_native'][old['gamma_band_to_mode']],hz)
            if np.min(overlap[np.arange(4),assignment])<.99:raise ValueError('Mode assignment changed')
            norm=np.sqrt(np.mean(abs(hz)**2,axis=1))
            small=np.array([sampler(int(m),anchor['rho_xyz_um'],E_EXPR)/norm[m] for m in assignment])
            previous=anchor['E_native_over_Hz_rms']
            phase=np.array([np.vdot(a,b)/np.vdot(a,a) for a,b in zip(small,previous)])
            relative=np.linalg.norm((small*phase[:,None,None]-previous).reshape(4,-1),axis=1)/np.linalg.norm(previous.reshape(4,-1),axis=1)
            if np.max(relative)>1e-7:raise ValueError('Eigenfield changed beyond its global phase')
            electric=np.array([sampler(int(m),rho,E_EXPR)/norm[m]*phase[b] for b,m in enumerate(assignment)])
            if mesh_identity(model)!=mesh:raise ValueError('Mesh changed')
            if index==0:
                from comsol_workflow.dipolar_radiation import pairs
                finite=np.load(output_path(OUTPUT, 's5_finite_volume_o3.npz'));k0=complex(finite['omega_rad_s'])/299792458.
                kval=f'({k0.real:.17g}+i*{k0.imag:.17g})[1/m]'
                node=j.result().numerical().create('s5meshcheck','IntVolume')
                node.set('data','dset1');node.selection().named('sel_slab_dom')
                node.set('intorder',sampler.JInt(8))
                node.set('expr',[f'2*9.89*ewfd.E{c}*cos({kval}*z)/1[V*m^2]' for c in 'xy'])
                native=[]
                for b,m in enumerate(assignment):
                    node.set('solnum',sampler.JArray(sampler.JInt,1)([int(m)+1]))
                    val=np.asarray(node.getReal()).ravel()+1j*np.asarray(node.getImag()).ravel()
                    native.append(val/norm[m]*phase[b]*k0*k0/(4*np.pi))
                native=np.asarray(native)
                numerical=k0*k0/(4*np.pi)*np.einsum('bnp,n->bp',electric[:,:,:2],2*9.89*weights[0]*np.cos(k0*rho[:,2]*1e-6))
                audit={'native_cell_F_over_Hrms':pairs(native),'mesh_quadrature_cell_F_over_Hrms':pairs(numerical),
                       'absolute_difference':abs(native-numerical).tolist(),'relative_norm_difference':float(np.linalg.norm(native-numerical)/np.linalg.norm(native))}
                write_json(OUT/'01_results/s5_mesh_quadrature_check.json',audit)
                print('Native FEM / mesh-rule check: '+str(audit),flush=True)
            np.savez_compressed(path,E_native_over_Hz_rms=electric,rho_xyz_um=rho,
                                region_weights_half_m3=weights,k_um_inv=k,phase_to_projection=phase)
            print(f'k{index:02d}: mesh-aligned radiation field saved',flush=True)
    finally:
        client.remove(model)


if __name__=='__main__':
    preflight=runtime_preflight()
    import mph
    from jpype import JClass
    client=mph.start(version='6.3',cores=4)
    write_json(output_path(OUTPUT, 's5_mesh_radiation_preflight.json'),preflight)
    if not JClass('com.comsol.model.util.ModelUtil').checkoutLicense('COMSOL','WAVEOPTICS'):raise RuntimeError('License unavailable')
    try:export(client)
    finally:client.clear()
