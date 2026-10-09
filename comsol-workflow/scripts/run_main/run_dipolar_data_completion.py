"""Complete mode19 E/H exports using saved COMSOL solutions; retain every MPH."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import traceback
import numpy as np
from matplotlib.path import Path as PolygonPath
from scipy.spatial import cKDTree
from scripts.analysis.analyze_dipolar_singularity import ROOT,SERIES,MODE,WORK,CONFIG,REF
from scripts.run_main.run_dipolar_volume_export import NativeSampler, cell_rule
from scripts.run_main.run_boundary_analysis import runtime_preflight
from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies
from comsol_workflow.simulation_utils import electric_reflection_matrix
from comsol_workflow.farfield_fft import load_air_field

MODEL=WORK/"00_model"
DATA=WORK/"01_results"
LOG=WORK/"80_logs"
EHS=[f"ewfd.E{c}/(1[V/m])" for c in "xyz"]+[f"ewfd.H{c}/(1[A/m])" for c in "xyz"]
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
SOURCE_FINITE=load_dipolar_inputs().finite_model
STATE={}


def put(name,data):
    (LOG/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")


def progress(message,**items):
    STATE.update(items,last_message=message,updated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()))
    put("comsol_execution.json",STATE)
    print(message,flush=True)


def fingerprint(path):
    path=Path(path)
    with path.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
    return dict(path=str(path.resolve()),bytes=path.stat().st_size,sha256=sha)


def snapshot(source,name):
    source=Path(source);target=(MODEL/name).resolve()
    assert target.parent==MODEL.resolve() and source.resolve()!=target
    record=fingerprint(source)
    if target.exists():
        if fingerprint(target)["sha256"]!=record["sha256"]:
            raise ValueError(f"Existing model snapshot differs: {target}")
    else:
        progress(f"Copying preserved model: {name}")
        shutil.copy2(source,target)
        assert fingerprint(target)["sha256"]==record["sha256"]
    records=json.loads((LOG/"model_manifest.json").read_text(encoding="utf8")) if (LOG/"model_manifest.json").exists() else {}
    records[name]={"source":record,"copy":fingerprint(target),"retain":True}
    put("model_manifest.json",records)
    return target


def reflect_sample(sampler,mode,xyz):
    sx=np.diag(electric_reflection_matrix("y","PMC"))
    sy=np.diag(electric_reflection_matrix("x","PEC"))
    rx=np.r_[sx,-sx];ry=np.r_[sy,-sy]
    values=sampler(mode,np.abs(xyz),EHS)
    values*=np.where(xyz[:,0,None]<0,rx,1)
    values*=np.where(xyz[:,1,None]<0,ry,1)
    # Saved finite model is upper-half TE; this branch only serves explicit negative-z queries.
    values*=np.where(xyz[:,2,None]<0,np.array([1,1,-1,-1,-1,1]),1)
    return values


def plane_export(model,sampler,cfg,name,z_um,axis,material_domains):
    path=DATA/f"{name}.npz"
    if path.exists():
        progress(f"Resume verified plane: {name}")
        with np.load(path) as d:
            assert float(d["z_um"])==z_um
            np.testing.assert_allclose(d["x_um"],axis)
        return
    xx,yy=np.meshgrid(axis,axis)
    mask=PolygonPath(np.asarray(cfg["finite_boundary"])).contains_points(
        np.column_stack([xx.ravel(),yy.ravel()]),radius=-1e-8).reshape(xx.shape)
    assert len(axis)%2==1
    np.testing.assert_allclose(axis,-axis[::-1],atol=1e-12)
    assert np.array_equal(mask,mask[::-1]) and np.array_equal(mask,mask[:,::-1])
    # Query the stored quarter solution once per distinct point, then reflect E/H.
    quarter=mask&(xx>=0)&(yy>=0)
    xyz=np.column_stack([xx[quarter],yy[quarter],np.full(quarter.sum(),z_um)])
    out=np.full((*xx.shape,6),np.nan+1j*np.nan)
    domains=np.full(xx.shape,-1,int)
    pieces=[]
    for start in range(0,len(xyz),30000):
        pieces.append(reflect_sample(sampler,19,xyz[start:start+30000]))
        if start%150000==0:progress(f"{name}: {min(start+30000,len(xyz))}/{len(xyz)} quarter points")
    out[quarter]=np.concatenate(pieces)
    domains[quarter]=sampler(19,xyz,["dom"])[:,0].real.astype(int)
    reflected=np.maximum(np.arange(len(axis)),len(axis)-1-np.arange(len(axis)))
    out=out[reflected[:,None],reflected[None,:]]
    domains=domains[reflected[:,None],reflected[None,:]]
    sx=np.diag(electric_reflection_matrix("y","PMC"));sy=np.diag(electric_reflection_matrix("x","PEC"))
    out*=np.where(xx[...,None]<0,np.r_[sx,-sx],1)
    out*=np.where(yy[...,None]<0,np.r_[sy,-sy],1)
    # Native E/H traces on FEM edges can depend on sub-ulp coordinate differences.
    # Evaluate symmetry axes at their EXACT requested coordinates, not a mirrored node.
    axis_mask=mask&((xx==0)|(yy==0))
    axis_xyz=np.column_stack([xx[axis_mask],yy[axis_mask],np.full(axis_mask.sum(),z_um)])
    out[axis_mask]=reflect_sample(sampler,19,axis_xyz)
    domains[axis_mask]=sampler(19,np.abs(axis_xyz),["dom"])[:,0].real.astype(int)
    eps=np.where(mask,np.where(np.isin(domains,material_domains),
        float(cfg["simulation_common"]["refractive_index"])**2,1.),np.nan)
    np.savez_compressed(path,E_V_m=out[...,:3],H_A_m=out[...,3:],x_um=axis,y_um=axis,
                        z_um=z_um,valid=mask,epsilon_r=eps,domain_id=domains,
                        mode_idx=19,time_convention="exp(+i omega t)")
    progress(f"Saved actual E/H plane: {name}")


def native_integrals(model,sampler,material_domains,slab_domains,k0):
    from scipy.constants import epsilon_0,mu_0
    j=model.java
    node=j.result().numerical().create("dipolar_integrals","IntVolume")
    node.set("data","dset1");node.set("solnum",sampler.JArray(sampler.JInt,1)([20]))
    node.set("intorder",sampler.JInt(6))
    eps="real(ewfd.epsilonrxx)"
    expressions=[f"8*epsilon0_const*{eps}*abs(ewfd.E{a})^2/4/1[J]" for a in "xyz"]
    expressions += [f"8*mu0_const*abs(ewfd.H{a})^2/4/1[J]" for a in "xyz"]
    node.selection().set(sampler.JArray(sampler.JInt,1)(list(map(int,slab_domains))))
    node.set("expr",expressions);node.set("unit",["1"]*6)
    energy=np.asarray(node.getReal()).ravel()
    node.selection().set(sampler.JArray(sampler.JInt,1)(list(map(int,material_domains))))
    sx=np.diag(electric_reflection_matrix("y","PMC"))
    sy=np.diag(electric_reflection_matrix("x","PEC"))
    mirror=2*(1+sx)*(1+sy)
    node.set("expr",[f"{mirror[i]:.17g}*(real(ewfd.epsilonrxx)-1)*ewfd.E{a}*cos({k0:.17g}[1/m]*z)/1[V*m^2]" for i,a in enumerate("xy")])
    node.set("unit",["1"]*2)
    raw=np.asarray(node.getReal()).ravel()+1j*np.asarray(node.getImag()).ravel()
    data={"slab_energy_components_J":energy.tolist(),"component_order":["Ex","Ey","Ez","Hx","Hy","Hz"],
          "U_ref_J":float(energy.sum()),"full_material_center_F_V_pairs":
          np.column_stack([(raw*k0*k0/(4*np.pi)).real,(raw*k0*k0/(4*np.pi)).imag]).tolist(),
          "quadrature_order":6,"mirror_factor":8,
          "slab_domains":list(map(int,slab_domains)),"material_domains":list(map(int,material_domains))}
    put("native_full_integrals.json",data)
    progress("Native complete-slab energy and material-source center integrals saved.")
    return data


def cells_export(sampler,cfg,material_domains,order=3):
    rho,w,material=cell_rule(order)
    np.savez_compressed(DATA/"cell_quadrature.npz",rho_xyz_um=rho,weights_half_m3=w,material_masks=material)
    flip_xy=[np.array([sx,sy]) for sx,sy in [(1,1),(-1,1),(1,-1),(-1,-1)]]
    tree=cKDTree(rho);mappings=[]
    sx=np.diag(electric_reflection_matrix("y","PMC"));sy=np.diag(electric_reflection_matrix("x","PEC"))
    signs=[np.ones(6),np.r_[sx,-sx],np.r_[sy,-sy],np.r_[sx*sy,sx*sy]]
    for f in flip_xy:
        distance,index=tree.query(rho*np.r_[f,1])
        assert distance.max()<1e-10
        mappings.append(index)
    for region,key in enumerate(["bulk_points","cladding_points"]):
        path=DATA/f"finite_{'cavity' if region==0 else 'cladding'}_EH_o{order}.npz"
        if path.exists():
            progress(f"Resume completed region: {key}")
            continue
        centers=np.array([[p["x"],p["y"]] for p in cfg[key]])
        center_tree=cKDTree(centers)
        ids=np.flatnonzero((centers[:,0]>=-1e-10)&(centers[:,1]>=-1e-10))
        values=np.empty((len(centers),len(rho),6),complex);seen=np.zeros(len(centers),bool)
        for start in range(0,len(ids),8):
            batch=ids[start:start+8]
            xyz=rho[None,:,:]+np.pad(centers[batch],((0,0),(0,1)))[:,None,:]
            if start==0:
                dom=sampler(19,np.abs(xyz[0]),["dom"])[:,0].real.astype(int)
                assert np.array_equal(np.isin(dom,material_domains),material[region])
            raw=reflect_sample(sampler,19,xyz.reshape(-1,3)).reshape(len(batch),len(rho),6)
            for flip,mapping,sign in zip(flip_xy,mappings,signs):
                distance,target=center_tree.query(centers[batch]*flip)
                assert distance.max()<1e-10
                values[target]=raw[:,mapping,:]*sign;seen[target]=True
            if start%80==0:progress(f"{key} E/H: {min(start+8,len(ids))}/{len(ids)} quarter cells")
        assert seen.all()
        np.savez_compressed(path,E_V_m=values[...,:3],H_A_m=values[...,3:],centers_um=centers,
                            rho_xyz_um=rho,weights_half_m3=w,material_mask=material[region],
                            time_convention="exp(+i omega t)")
        progress(f"Saved complete cell E/H: {key}")


def finite(client):
    from scipy.constants import c
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    assert cfg["simulation_common"]["mesh_auto_size"]==5
    assert cfg["symmetry"]=={"id":1,"x_boundary":"PEC","y_boundary":"PMC"}
    target=snapshot(SOURCE_FINITE,"finite_mode19_source.mph")
    progress("Loading preserved finite model; no finite eigensolve.")
    model=client.load(str(target));j=model.java
    try:
        f=frequencies(model)
        np.testing.assert_allclose(f.iloc[19][["re","im"]],[197.9530009044521,.0128930496899445],rtol=1e-9)
        sampler=NativeSampler(model)
        selections={}
        for tag in j.component("comp1").selection().tags():
            try:selections[str(tag)]=list(map(int,j.component("comp1").selection(str(tag)).entities(3)))
            except Exception:pass
        mesh=j.component("comp1").mesh("mesh1")
        params={str(k):str(j.param().get(str(k))) for k in j.param().varnames()}
        put("finite_model_identity.json",{"parameters":params,"selections":selections,
            "mesh_vertices":int(mesh.getNumVertex()),"mesh_tetrahedra":int(mesh.getNumElem("tet")),
            "frequency_THz":f.iloc[19].to_dict(),"model_path":str(target)})
        material=np.array(selections["sel_layer_0_mat_dom"],int)
        # Identify actual slab domains by native domain representatives, including hole air.
        vertices=np.asarray(mesh.getVertex(),float).T
        elements=np.asarray(mesh.getElem("tet"),int).T
        entities=np.asarray(mesh.getElemEntity("tet"),int)
        centroids=vertices[elements].mean(axis=1)
        slab_ids=np.unique(entities[(centroids[:,2]>1e-10)&(centroids[:,2]<.1-1e-10)])
        slab_ids=np.union1d(slab_ids,material)
        np.savez_compressed(DATA/"finite_mesh_geometry.npz",vertices_um=vertices,
                            tetrahedra=elements,domain_ids=entities,slab_domains=slab_ids,
                            material_domains=material)
        progress(f"Finite model loaded: {len(elements)} tetrahedra; slab domains {len(slab_ids)}")
        old=load_air_field(MODE/"11_simulation_exports/E_air.parquet")
        ix=np.array([409,431,455]);iy=np.array([407,429,451])
        xyz=np.column_stack([old["x_axis"][ix],old["y_axis"][iy],np.full(3,1.65)])
        check=reflect_sample(sampler,19,xyz)
        ref=np.column_stack([old[k][iy,ix] for k in ["Ex","Ey","Ez"]])
        error=float(np.linalg.norm(check[:,:3]-ref)/np.linalg.norm(ref))
        put("native_phase_check.json",{"relative_E_error":error,"native_phase_verified":error<1e-8})
        assert error<1e-8
        axis=np.linspace(-max(abs(np.asarray(cfg["finite_boundary"])[:,0])),max(abs(np.asarray(cfg["finite_boundary"])[:,0])),1601)
        plane_export(model,sampler,cfg,"finite_slab_z0_EH",0.,axis,material)
        for name,z in [("finite_surface_EH",.100001),("finite_air_EH",1.65),("finite_air_high_EH",1.95)]:
            plane_export(model,sampler,cfg,name,z,old["x_axis"],material)
        native_integrals(model,sampler,material,slab_ids,2*np.pi*float(f.iloc[19]["re"])*1e12/c)
        cells_export(sampler,cfg,material)
        progress("Finite E/H plane, cell and complete native integral exports finished.",finite_complete=True)
    finally:
        client.remove(model)



def periodic(client):
    from comsol_workflow.dipolar_radiation import mode_assignment
    manifest=json.loads((REF/"99_config/manifest.json").read_text(encoding="utf8"))
    source=Path(manifest["source_model"])
    target=snapshot(source,"periodic_Gamma_source.mph")
    quadrature=DATA/"cell_quadrature.npz"
    with np.load(quadrature) as d:
        rho=d["rho_xyz_um"];w=d["weights_half_m3"];masks=d["material_masks"]
    center_axis=np.linspace(-.82/2,.82/2,161)
    # Native source cell is a vertical hexagon; display on the known source sampling grid.
    with np.load(REF/"01_results/k00_reference.npz") as anchor:
        center_points=np.column_stack([anchor["rho_points_um"],np.zeros(len(anchor["rho_points_um"]))])
    for index,k in enumerate(manifest["k_um_inv"]):
        output=DATA/f"periodic_k{index:02d}_EH.npz"
        saved=MODEL/f"periodic_k{index:02d}.mph"
        if output.exists() and saved.exists():
            progress(f"Resume completed periodic k{index:02d}",periodic_last_completed=index)
            continue
        progress(f"Periodic k{index:02d}/18: load retained model")
        model=client.load(str(saved if saved.exists() else target));j=model.java
        try:
            if not saved.exists() and index:
                for name,value in zip(["kx","ky"],k):j.param().set(name,f"{value:.17g}[1/um]")
                progress(f"Periodic k{index:02d}: fixed mesh5 eigensolve, four modes")
                j.sol("sol1").clearSolutionData();j.sol("sol1").runAll()
            # Preserve the solved model before any export can fail.
            if not saved.exists():
                model.save(str(saved))
                progress(f"Retained MPH: periodic_k{index:02d}.mph")
            if "s4surface_freq" in list(j.result().numerical().tags()):
                j.result().numerical().remove("s4surface_freq")
            freq=frequencies(model)
            assert len(freq)==4
            sampler=NativeSampler(model)
            center=np.array([sampler(mode,center_points,EHS) for mode in range(4)])
            with np.load(REF/f"01_results/k{index:02d}_reference.npz") as old:
                reference=old["Hz_native"][old["gamma_band_to_mode"]]
            assignment,overlap=mode_assignment(reference,center[...,5])
            assert np.min(overlap[np.arange(4),assignment])>.99
            chosen=center[assignment]
            # Align only the global phase to the raw reference, then normalize full E/H together.
            phase=np.array([np.vdot(h,ref)/abs(np.vdot(h,ref)) for h,ref in zip(chosen[...,5],reference)])
            eh=np.array([sampler(int(mode),rho,EHS)*phase[b] for b,mode in enumerate(assignment)])
            chosen*=phase[:,None,None]
            from scipy.constants import epsilon_0,mu_0
            metric=2*w
            energy=np.sum(metric[None,:]*(epsilon_0/4*(1+9.89*masks[0])[None,:]*np.sum(abs(eh[...,:3])**2,-1)
                     +mu_0/4*np.sum(abs(eh[...,3:])**2,-1)),axis=1)
            assert np.all(energy>0)
            normalization=np.sqrt(energy)
            np.savez_compressed(output,E_V_m=eh[...,:3],H_A_m=eh[...,3:],E_center_V_m=chosen[...,:3],
                    H_center_A_m=chosen[...,3:],center_xyz_um=center_points,rho_xyz_um=rho,
                    weights_half_m3=w,epsilon_r=1+9.89*masks[0],cell_EM_energy_J=energy,
                    E_normalized=eh[...,:3]/normalization[:,None,None],
                    H_normalized=eh[...,3:]/normalization[:,None,None],
                    k_um_inv=k,mode_assignment=assignment,overlap=overlap,
                    phase_to_raw_reference=phase,frequencies_THz=freq[["re","im","q"]].to_numpy()[assignment],
                    time_convention="exp(+i omega t)")
            progress(f"Periodic k{index:02d}: complete E/H, z=0 and EM normalization saved.",periodic_last_completed=index)
        finally:
            client.remove(model)
    progress("All 19 periodic E/H references exported; all MPH files retained.",periodic_complete=True)



def mesh_export(client):
    """Independent native tetrahedral fields, including the actual outer pad."""
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    progress("Loading retained finite solution for air-domain and native-mesh audit.")
    model=client.load(str(MODEL/"finite_mode19_source.mph"));j=model.java
    try:
        sampler=NativeSampler(model)
        material=np.asarray(j.component("comp1").selection("sel_layer_0_mat_dom").entities(3),int)
        slab=np.asarray(j.component("comp1").selection("sel_layer_0_extent").entities(3),int)
        pml=np.asarray(j.component("comp1").selection("sel_pml_dom").entities(3),int)
        with np.load(DATA/"finite_air_EH.npz") as d:axis=d["x_um"]
        for name,z in [("finite_air_1100nm_EH",1.1),("finite_air_1400nm_EH",1.4)]:
            plane_export(model,sampler,cfg,name,z,axis,material)
            with np.load(DATA/f"{name}.npz") as d:
                assert not np.isin(d["domain_id"][d["valid"]],pml).any()
        put("air_domain_audit.json",{"excluded":["finite_air_high_EH.npz"],
            "reason":"z=1.95 um is inside PML, whose lower boundary is z=1.65 um",
            "valid_interior_air_planes_um":[1.1,1.4],"pml_domains":pml.tolist(),
            "original_interface_plane_um":1.65,"surface_air_side_um":.100001})
        with np.load(DATA/"finite_mesh_geometry.npz") as d:
            vertices=d["vertices_um"];tet=d["tetrahedra"];dom=d["domain_ids"]
        chosen=np.flatnonzero(np.isin(dom,slab))
        # Positive degree-2 rule; compare to independent order-6 native integrals.
        b=(5-np.sqrt(5))/20;a=1-3*b
        bary=np.full((4,4),b);np.fill_diagonal(bary,a)
        records=[]
        for part,start in enumerate(range(0,len(chosen),25000)):
            path=DATA/f"finite_mesh_EH_{part:03d}.npz"
            ids=chosen[start:start+25000]
            if not path.exists():
                v=vertices[tet[ids]]
                xyz=np.einsum("qv,tvc->tqc",bary,v).reshape(-1,3)
                w=np.repeat(abs(np.linalg.det(v[:,1:]-v[:,:1]))/24,4)*1e-18
                eh=sampler(19,xyz,EHS)
                np.savez_compressed(path,xyz_um=xyz,weights_octant_m3=w,
                    E_V_m=eh[:,:3],H_A_m=eh[:,3:],
                    material_mask=np.repeat(np.isin(dom[ids],material),4),tetrahedron_ids=ids,
                    epsilon_r=np.repeat(np.where(np.isin(dom[ids],material),
                        float(cfg["simulation_common"]["refractive_index"])**2,1.),4))
            records.append({"path":path.name,"tetrahedra":len(ids),"points":4*len(ids)})
            progress(f"Native slab tetrahedra: {min(start+len(ids),len(chosen))}/{len(chosen)}")
        put("finite_mesh_exports.json",{"parts":records,"rule":"positive symmetric degree 2, four points/tetrahedron",
            "octant":"x>=0,y>=0,z>=0; E/H vector and axial-vector reflection required",
            "slab_domain_count":len(slab),"material_domain_count":len(material)})
        progress("Actual full-slab tetrahedral E/H and interior-air exports complete.",mesh_complete=True)
    finally:
        client.remove(model)



def periodic_mesh_export(client):
    """Native tetrahedral radiation quadrature, on each retained periodic solution."""
    from scipy.constants import c
    b=(5-np.sqrt(5))/20
    bary=np.full((4,4),b);np.fill_diagonal(bary,1-3*b)
    for index in range(19):
        output=DATA/f"periodic_k{index:02d}_mesh_EH.npz"
        if output.exists():
            progress(f"Resume periodic native mesh k{index:02d}");continue
        model=client.load(str(MODEL/f"periodic_k{index:02d}.mph"));j=model.java
        try:
            sampler=NativeSampler(model)
            mesh=j.component("comp1").mesh("mesh1")
            vertices=np.asarray(mesh.getVertex(),float).T
            tet=np.asarray(mesh.getElem("tet"),int).T
            dom=np.asarray(mesh.getElemEntity("tet"),int)
            slab=np.asarray(j.component("comp1").selection("sel_slab_layer").entities(3),int)
            silicon=np.asarray(j.component("comp1").selection("sel_slab_dom").entities(3),int)
            selected=np.flatnonzero(np.isin(dom,slab));v=vertices[tet[selected]]
            xyz=np.einsum("qv,tvc->tqc",bary,v).reshape(-1,3)
            w=np.repeat(abs(np.linalg.det(v[:,1:]-v[:,:1]))/24,4)*1e-18
            material=np.repeat(np.isin(dom[selected],silicon),4)
            with np.load(DATA/f"periodic_k{index:02d}_EH.npz") as d:
                assignment=d["mode_assignment"];phase=d["phase_to_raw_reference"];normal=d["cell_EM_energy_J"]
            eh=np.array([sampler(int(mode),xyz,EHS)*phase[band] for band,mode in enumerate(assignment)])
            node=j.result().numerical().create("dipolar_native_center","IntVolume")
            node.set("data","dset1");node.selection().named("sel_slab_dom")
            node.set("intorder",sampler.JInt(6))
            kval=2*np.pi*197.9530009044521e12/c
            expr=[f"2*9.89*ewfd.E{a}*cos({kval:.17g}[1/m]*z)/1[V*m^2]" for a in "xy"]
            node.set("expr",expr);node.set("unit",["1","1"])
            centers=[]
            for band,mode in enumerate(assignment):
                node.set("solnum",sampler.JArray(sampler.JInt,1)([int(mode)+1]))
                value=np.asarray(node.getReal()).ravel()+1j*np.asarray(node.getImag()).ravel()
                centers.append(value*phase[band]/np.sqrt(normal[band])*kval*kval/(4*np.pi))
            np.savez_compressed(output,xyz_um=xyz,weights_half_m3=w,material_mask=material,
                E_normalized=eh[...,:3]/np.sqrt(normal)[:,None,None],
                H_normalized=eh[...,3:]/np.sqrt(normal)[:,None,None],
                source_center_native_F_V_per_sqrtJ=np.array(centers))
            progress(f"Periodic native FEM E/H k{index:02d}: {len(selected)} tets, {len(xyz)} points")
        finally:client.remove(model)
    progress("Native periodic source quadrature completed.",periodic_mesh_complete=True)



def convergence_export(client):
    from scipy.constants import c
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    model=client.load(str(MODEL/"finite_mode19_source.mph"));j=model.java
    try:
        sampler=NativeSampler(model)
        mat=np.asarray(j.component("comp1").selection("sel_layer_0_mat_dom").entities(3),int)
        with np.load(DATA/"finite_slab_z0_EH.npz") as d:axis=d["x_um"]
        plane_export(model,sampler,cfg,"finite_air_1100nm_fine_EH",1.1,axis,mat)
        physics={}
        for tag in j.component("comp1").physics("ewfd").feature().tags():
            f=j.component("comp1").physics("ewfd").feature(str(tag))
            try:
                physics[str(tag)]={"type":str(f.getType()),"active":bool(f.isActive()),
                    "selection":list(map(int,f.selection().entities()))}
            except Exception as exc:physics[str(tag)]={"metadata_error":str(exc)}
        put("finite_boundary_conditions.json",physics)
        node=j.result().numerical().create("dipolar_directional","IntVolume")
        node.set("data","dset1");node.set("solnum",sampler.JArray(sampler.JInt,1)([20]))
        node.selection().named("sel_layer_0_mat_dom");node.set("intorder",sampler.JInt(6))
        identity=json.loads((LOG/"finite_model_identity.json").read_text())
        k0=2*np.pi*identity["frequency_THz"]["re"]*1e12/c
        dx=(axis[1]-axis[0])*2e-6;dq=2*np.pi/(2001*dx)
        rows=[]
        for ix,iy in [(0,0),(20,0),(0,40),(50,40),(80,20),(100,20)]:
            qx,qy=ix*dq,iy*dq;kz=np.sqrt(k0*k0-qx*qx-qy*qy)
            n=np.array([qx,qy,kz])/k0
            expr=[]
            for pol,a in enumerate("xyz"):
                sx=(-1,1,1)[pol];sy=(-1,1,-1)[pol];sz=(1,1,-1)[pol]
                factors=[f"(exp(i*{v:.17g}[1/m]*{coord})+({sign})*exp(-i*{v:.17g}[1/m]*{coord}))"
                    for v,coord,sign in zip([qx,qy,kz],"xyz",[sx,sy,sz])]
                expr.append("9.89*ewfd.E"+a+"*"+"*".join(factors)+"/1[V*m^2]")
            node.set("expr",expr);node.set("unit",["1"]*3)
            value=np.asarray(node.getReal()).ravel()+1j*np.asarray(node.getImag()).ravel()
            far=(value-n*np.dot(n,value))*k0*k0/(4*np.pi)
            rows.append({"qx_index":ix,"qy_index":iy,"q_m_inv":[qx,qy],
                "F_V_pairs":np.column_stack([far.real,far.imag]).tolist()})
            progress(f"Native directional radiation check: q indices ({ix},{iy})")
        put("native_directional_radiation.json",{"intorder":6,"directions":rows})
    finally:client.remove(model)


def sampling_export(client):
    """Double spatial sampling from retained solutions, never re-solve or save MPH."""
    from comsol_workflow.hex_lattice_utils import unit_cell_corners
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    before={p.name:(p.stat().st_size,p.stat().st_mtime_ns) for p in MODEL.glob("*.mph")}
    checks=[]
    model=client.load(str(MODEL/"finite_mode19_source.mph"))
    try:
        freq=frequencies(model)
        np.testing.assert_allclose(freq.iloc[19][["re","im"]],[197.9530009044521,.0128930496899445],rtol=1e-9)
        sampler=NativeSampler(model)
        material=np.asarray(model.java.component("comp1").selection("sel_layer_0_mat_dom").entities(3),int)
        pml=np.asarray(model.java.component("comp1").selection("sel_pml_dom").entities(3),int)
        for oldname,newname,z in [
            ("finite_slab_z0_EH","finite_slab_z0_fine_EH",0.),
            ("finite_surface_EH","finite_surface_fine_EH",.100001),
            ("finite_air_1100nm_EH","finite_air_1100nm_fine_EH",1.1),
            ("finite_air_1400nm_EH","finite_air_1400nm_fine_EH",1.4),
            ("finite_air_EH","finite_air_fine_EH",1.65)]:
            with np.load(DATA/f"{oldname}.npz") as old:
                axis=old["x_um"]
                fineaxis=np.linspace(axis[0],axis[-1],2*len(axis)-1)
                plane_export(model,sampler,cfg,newname,z,fineaxis,material)
                if newname=="finite_slab_z0_fine_EH":
                    with np.load(DATA/f"{newname}.npz") as saved:
                        if "exact_axis_coordinates" not in saved.files:
                            record={k:saved[k] for k in saved.files}
                            xx,yy=np.meshgrid(fineaxis,fineaxis)
                            am=record["valid"]&((xx==0)|(yy==0))
                            xyz=np.column_stack([xx[am],yy[am],np.full(am.sum(),z)])
                            value=reflect_sample(sampler,19,xyz)
                            record["E_V_m"][am]=value[:,:3];record["H_A_m"][am]=value[:,3:]
                            np.savez_compressed(DATA/f"{newname}.npz",**record,exact_axis_coordinates=True)
                            progress("Re-evaluated exact axis coordinates in newly refined slab export")
                if newname=="finite_air_fine_EH":
                    with np.load(DATA/f"{newname}.npz") as saved:
                        if "exact_interface_coordinates" not in saved.files:
                            record={k:saved[k] for k in saved.files}
                            xx,yy=np.meshgrid(fineaxis,fineaxis);valid=record["valid"]
                            xyz=np.column_stack([xx[valid],yy[valid],np.full(valid.sum(),z)])
                            pieces=[]
                            for start in range(0,len(xyz),30000):
                                pieces.append(reflect_sample(sampler,19,xyz[start:start+30000]))
                                if start%150000==0:progress(f"Exact interface coordinates: {min(start+30000,len(xyz))}/{len(xyz)}")
                            value=np.concatenate(pieces)
                            record["E_V_m"][valid]=value[:,:3];record["H_A_m"][valid]=value[:,3:]
                            record["domain_id"][valid]=sampler(19,np.abs(xyz),["dom"])[:,0].real.astype(int)
                            record["epsilon_r"][valid]=np.where(np.isin(record["domain_id"][valid],material),float(cfg["simulation_common"]["refractive_index"])**2,1.)
                            np.savez_compressed(DATA/f"{newname}.npz",**record,exact_interface_coordinates=True)
                with np.load(DATA/f"{newname}.npz") as new:
                    np.testing.assert_allclose(new["x_um"][::2],axis,atol=1e-12)
                    valid=old["valid"]
                    assert np.array_equal(new["valid"][::2,::2],valid)
                    assert np.array_equal(new["domain_id"][::2,::2],old["domain_id"])
                    row={"source":oldname+".npz","refined":newname+".npz","old_shape":[len(axis)]*2,
                         "new_shape":[len(fineaxis)]*2,"old_step_um":float(axis[1]-axis[0]),
                         "new_step_um":float(fineaxis[1]-fineaxis[0]),"z_um":z}
                    for key in ["E_V_m","H_A_m"]:
                        a=old[key][valid];b=new[key][::2,::2][valid]
                        row[key+"_nested_relative_error"]=float(np.linalg.norm(a-b)/np.linalg.norm(a))
                        assert row[key+"_nested_relative_error"]<1e-8
                    if z in [1.1,1.4,.100001]:assert not np.isin(new["domain_id"][new["valid"]],pml).any()
                    checks.append(row)
                    put("sampling_export_validation.json",{"planes":checks,"status":"running"})
    finally:client.remove(model)
    periodic_checks=[]
    for index in range(19):
        output=DATA/f"periodic_k{index:02d}_z0_fine_EH.npz"
        with np.load(DATA/f"periodic_k{index:02d}_EH.npz") as old:
            points=old["center_xyz_um"]
            ax=[np.linspace(points[:,j].min(),points[:,j].max(),2*len(np.unique(points[:,j]))-1) for j in range(2)]
            xx,yy=np.meshgrid(*ax);xy=np.column_stack([xx.ravel(),yy.ravel()])
            keep=PolygonPath(unit_cell_corners(cfg["a"])).contains_points(xy,radius=1e-10)
            xyz=np.column_stack([xy[keep],np.zeros(keep.sum())])
            if not output.exists():
                progress(f"Periodic k{index:02d}: querying {len(xyz)} z=0 points for all four modes")
                model=client.load(str(MODEL/f"periodic_k{index:02d}.mph"))
                try:
                    sampler=NativeSampler(model);assignment=old["mode_assignment"]
                    np.testing.assert_allclose(frequencies(model)[["re","im","q"]].to_numpy()[assignment],old["frequencies_THz"],rtol=1e-9)
                    eh=np.array([sampler(int(mode),xyz,EHS)*old["phase_to_raw_reference"][b] for b,mode in enumerate(assignment)])
                    np.savez_compressed(output,center_xyz_um=xyz,E_center_V_m=eh[...,:3],H_center_A_m=eh[...,3:],
                        **{k:old[k] for k in ["k_um_inv","mode_assignment","phase_to_raw_reference","cell_EM_energy_J","frequencies_THz","time_convention"]})
                finally:client.remove(model)
            with np.load(output) as new:
                distance,ids=cKDTree(new["center_xyz_um"]).query(points)
                assert distance.max()<1e-12
                row={"k_index":index,"old_points":len(points),"new_points":len(xyz),"background_grid_before":[41,41],"background_grid_after":[81,81]}
                for key in ["E_center_V_m","H_center_A_m"]:
                    a=old[key];b=new[key][:,ids]
                    row[key+"_nested_relative_error"]=float(np.linalg.norm(a-b)/np.linalg.norm(a))
                    assert row[key+"_nested_relative_error"]<1e-8
                periodic_checks.append(row)
        progress(f"Verified periodic z=0 sampling k{index:02d}/18")
    after={p.name:(p.stat().st_size,p.stat().st_mtime_ns) for p in MODEL.glob("*.mph")}
    assert before==after
    put("sampling_export_validation.json",{"status":"complete","planes":checks,"periodic":periodic_checks,
        "model_files_unchanged":True,"model_count":len(after),"FEM_mesh_refined":False,"new_eigensolves":0})


def run(stage):
    for folder in [MODEL,DATA,LOG]:folder.mkdir(exist_ok=True)
    STATE.update(stage=stage,status="running",pid=os.getpid(),started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()))
    put("comsol_preflight.json",runtime_preflight())
    import mph
    from jpype import JClass
    client=mph.start(version="6.3",cores=4)
    try:
        util=JClass("com.comsol.model.util.ModelUtil")
        if not util.checkoutLicense("COMSOL","WAVEOPTICS"):raise RuntimeError("COMSOL/WAVEOPTICS license unavailable")
        util.showProgress(str(LOG/"comsol_progress.log"))
        progress("COMSOL 6.3 started; COMSOL and WAVEOPTICS licenses checked out.")
        if stage=="finite":finite(client)
        elif stage=="periodic":periodic(client)
        elif stage=="mesh":mesh_export(client)
        elif stage=="convergence":convergence_export(client)
        elif stage=="sampling":sampling_export(client)
        elif stage=="periodic_mesh":periodic_mesh_export(client)
        elif stage=="remaining":
            mesh_export(client)
            periodic(client)
        else:raise ValueError(stage)
        progress("Requested export stage completed.",status="complete")
    except Exception as exc:
        progress(str(exc),status="failed",traceback=traceback.format_exc())
        raise
    finally:
        client.clear()


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=["finite","periodic","mesh","remaining","periodic_mesh","convergence","sampling"])
    args=parser.parse_args()
    run(args.stage)
