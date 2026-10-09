"""Complete E/H energy and radiation analysis from this task's COMSOL exports."""
from __future__ import annotations
from comsol_workflow.output_paths import table_links, table_path
from comsol_workflow.dipolar_inputs import load_dipolar_inputs
import argparse,json,math,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.constants import epsilon_0,mu_0,c
from scipy.spatial import cKDTree
from scipy.interpolate import LinearNDInterpolator
import finufft
from scripts.analysis.analyze_dipolar_singularity import WORK, REPORTS, CONFIG, REF, Z0, plt, LogNorm, colorbar, save, write_json, log
from comsol_workflow.dipolar_radiation import forward, inverse, energy_group_gram, budget_rows
D=WORK/"01_results";L=WORK/"80_logs";O=WORK/"10_overview";R=REPORTS
NAMES=["Gamma-py","Around-Gamma py","Other cavity bands","Cavity residual","Actual cladding","Outer pad"]
FREQ=json.loads((L/"finite_model_identity.json").read_text(encoding="utf8"))["frequency_THz"]["re"]*1e12
K0=2*np.pi*FREQ/c
SZ=1e-7
FFT_SIZE=8001
FIGURE_DPI=380

def load(name):
    with np.load(D/name) as d:return {k:d[k] for k in d.files}

def native():
    return json.loads((L/"native_full_integrals.json").read_text(encoding="utf8"))

def axes():
    with np.load(D/"finite_air_1100nm_fine_EH.npz") as d:dx=np.diff(d["x_um"])[0]*1e-6
    q=(np.arange(FFT_SIZE)-FFT_SIZE//2)*2*np.pi/(FFT_SIZE*dx)
    h=int(np.ceil(K0/(q[1]-q[0])))
    small=q[FFT_SIZE//2-h:FFT_SIZE//2+h+1]
    qx,qy=np.meshgrid(small,small);mask=qx*qx+qy*qy<K0*K0
    kz=np.sqrt(np.maximum(K0*K0-qx*qx-qy*qy,0))
    n=np.stack([qx/K0,qy/K0,kz/K0],-1)
    return dx,q,small,mask,kz,n

def radiation(xyz_um,electric,weights,octant=False):
    """Positive spatial Fourier integral via NUFFT; z mirrors evaluated by moments."""
    dx,q,small,mask,kz,n=axes();dq=small[1]-small[0];side=len(small)
    x=np.ascontiguousarray(xyz_um[:,0]*1e-6*dq)
    y=np.ascontiguousarray(xyz_um[:,1]*1e-6*dq)
    # Finite slab upper half; t is bounded by 0.415 for this model.
    t=K0*xyz_um[:,2]*1e-6
    plan=finufft.Plan(1,(side,side),eps=1e-12,isign=1,nthreads=4)
    plan.setpts(x,y)
    out=np.zeros((side,side,3),complex)
    for pol in range(3):
        source=np.ascontiguousarray(electric[:,pol]*weights)
        for m in range(5):
            power=2*m+(pol==2)
            v=plan.execute(np.ascontiguousarray(source*t**power)).T
            if octant:
                # Exact vector reflections for this quarter model.
                sx=(-1,1,1)[pol];sy=(-1,1,-1)[pol]
                v=v+sx*v[:,::-1]+sy*v[::-1,:]+sx*sy*v[::-1,::-1]
            coefficient=2*(1j if pol==2 else 1)*(-1)**m/math.factorial(power)
            out[...,pol]+=coefficient*v*(kz/K0)**power
    transverse=out-n*np.sum(n*out,axis=-1)[...,None]
    far=transverse*K0*K0/(4*np.pi);far[~mask]=0
    return far

def fit():
    arrays=[load(f"periodic_k{i:02d}_EH.npz") for i in range(19)]
    B=np.concatenate([np.concatenate([a["E_normalized"],a["H_normalized"]],axis=-1) for a in arrays])
    k=np.repeat(np.array([a["k_um_inv"] for a in arrays]),4,axis=0)
    f=load("finite_cavity_EH_o3.npz");w=2*f["weights_half_m3"]
    eps=arrays[0]["epsilon_r"]
    metric=np.column_stack([np.repeat((w*epsilon_0*eps/4)[:,None],3,1),
                            np.repeat((w*mu_0/4)[:,None],3,1)])
    V=np.concatenate([f["E_V_m"],f["H_A_m"]],axis=-1)
    phase=np.exp(-1j*f["centers_um"]@k.T)
    bm=(B*np.sqrt(metric)[None]).reshape(len(B),-1)
    cellgram=bm.conj()@bm.T
    gram=cellgram*(phase.conj().T@phase)
    right=np.zeros(len(B),complex)
    for start in range(0,len(V),40):
        v=(V[start:start+40]*np.sqrt(metric)[None]).reshape(-1,bm.shape[1])
        overlaps=v@bm.conj().T
        right+=np.sum(phase[start:start+len(v)].conj()*overlaps,axis=0)
    coefficient=np.linalg.solve(gram,right)
    U=float(np.sum(abs(V)**2*metric))
    grouped=np.zeros((3,len(B)),complex)
    ids=np.arange(len(B))
    grouped[0,1]=coefficient[1]
    grouped[1,(ids%4==1)&(ids!=1)]=coefficient[(ids%4==1)&(ids!=1)]
    grouped[2,ids%4!=1]=coefficient[ids%4!=1]
    groupgram=energy_group_gram(gram,right,U,grouped)
    relative=float(groupgram[-1,-1].real/U)
    np.savez_compressed(D/"field_decomposition.npz",coefficients=coefficient,
        grouped_coefficients=grouped,gram_J=gram,right=right,group_gram_J=groupgram,
        k_um_inv=k,band_index=ids%4,basis_EH_normalized=B,phase=phase,
        cavity_energy_J=U,relative_residual_energy=relative)
    details={"basis_count":len(B),"condition_number":float(np.linalg.cond(gram)),
        "cavity_energy_J":U,"residual_energy_fraction":relative,
        "normal_equation_relative_error":float(np.linalg.norm(gram@coefficient-right)/np.linalg.norm(right)),
        "group_self_energy_J":np.diag(groupgram).real.tolist(),"group_total_energy_J":float(groupgram.sum().real)}
    write_json("fit_validation.json",details);log("Joint E/H fit: "+json.dumps(details))
    rows=[]
    for i,name in enumerate(NAMES[:4]):
        rows.append(dict(region=name,kind="self",second=name,energy_J=groupgram[i,i].real))
        for j in range(i):
            rows.append(dict(region=NAMES[j],kind="interference",second=name,energy_J=2*groupgram[j,i].real))
    pd.DataFrame(rows).to_csv(table_path(L, "energy_budget.csv"),index=False)

def padded_transform_roi(spectrum,dx,size,points,*,inverse_transform=False):
    """Exact padded 2D transform ROI via separability; retain physical phase and scale."""
    from scipy.fft import fft,ifft,fftshift,ifftshift
    transform=fft if inverse_transform else ifft
    value=spectrum
    assert points%2==size%2==value.shape[0]%2==value.shape[1]%2
    crop=np.arange(size//2-points//2,size//2+points//2+1)
    for axis in [1,0]:
        padding=(size-value.shape[axis])//2
        width=[(0,0),(0,0)];width[axis]=(padding,padding)
        value=fftshift(transform(ifftshift(np.pad(value,width),axes=axis),axis=axis),axes=axis)
        value=np.take(value,crop,axis=axis)
    return value/(size**2*dx**2) if inverse_transform else value*(size**2*dx**2)


def air():
    dx,q,small,mask,kz,n=axes();h=len(small)//2
    outputs={};details={}
    for name in ["finite_air_1100nm_fine_EH","finite_air_1400nm_fine_EH","finite_air_fine_EH"]:
        d=load(name+".npz")
        E=np.stack([padded_transform_roi(np.nan_to_num(d["E_V_m"][...,i]),dx,len(q),len(small)) for i in range(3)],-1)
        H=np.stack([padded_transform_roi(np.nan_to_num(d["H_A_m"][...,i]),dx,len(q),len(small)) for i in range(3)],-1)
        P=np.stack([H[...,1],-H[...,0]],-1);nt=n[...,:2]
        s=np.where(mask,kz/K0,1.)
        v=Z0/s[...,None]*(P-nt*np.sum(nt*P,-1)[...,None])
        up=np.zeros_like(E);down=np.zeros_like(E)
        up[...,:2]=(E[...,:2]+v)/2;down[...,:2]=(E[...,:2]-v)/2
        up[...,2]=-np.sum(nt*up[...,:2],-1)/s
        down[...,2]=np.sum(nt*down[...,:2],-1)/s
        up[~mask]=0;down[~mask]=0
        z=float(d["z_um"])*1e-6
        far=1j*kz[...,None]/(2*np.pi)*np.exp(1j*kz*z)[...,None]*up
        outputs[name]=far
        weight=(small[1]-small[0])**2/(2*Z0*(2*np.pi)**2)*kz/K0
        pu=float(np.sum(weight[...,None]*abs(up)**2))
        pdn=float(np.sum(weight[...,None]*abs(down)**2))
        details[name]={"z_um":float(d["z_um"]),"up_power_W":pu,
            "down_power_W":pdn,"down_up_power_ratio":pdn/pu,
            "center_F_V_pairs":np.column_stack([far[h,h].real,far[h,h].imag]).tolist()}
        if name=="finite_air_1100nm_fine_EH":
            raw=E*np.exp(1j*kz*(z-SZ))[...,None];raw[~mask]=0
            plus=up*np.exp(1j*kz*(z-SZ))[...,None]
            hs=np.cross(n,plus)/Z0
            fields={}
            for label,spec in [("E_raw_surface",raw),("E_up_surface",plus),("H_up_surface",hs)]:
                vals=[]
                for i in range(3):
                    vals.append(padded_transform_roi(spec[...,i],dx,len(q),len(d["x_um"]),inverse_transform=True))
                fields[label]=np.stack(vals,-1)
            np.savez_compressed(D/"direct_radiative_fields.npz",**fields,
                F_V=far,A_up_surface=plus,small_q_m_inv=small,lightcone_mask=mask,
                x_um=d["x_um"],kz_m_inv=kz,frequency_Hz=FREQ,z_surface_m=SZ,
                z_air_m=z,source_E_angular=E,source_H_angular=H,fft_size=FFT_SIZE)
    a=outputs["finite_air_1100nm_fine_EH"];b=outputs["finite_air_1400nm_fine_EH"]
    details["height_center_relative_error"]=float(np.linalg.norm(a[h,h]-b[h,h])/np.linalg.norm(a[h,h]))
    details["height_full_F_relative_error"]=float(np.linalg.norm((a-b)[mask])/np.linalg.norm(a[mask]))
    write_json("air_validation.json",details);log("New air E/H propagation: "+json.dumps(details))

def sampling_checks():
    """Isolate raw sampling effects at identical physical q, independent of padding."""
    dx,q,small,mask,kz,n=axes();h=len(small)//2;dq=small[1]-small[0]
    coarse=load("finite_air_1100nm_EH.npz");fine=load("direct_radiative_fields.npz")
    axis=coarse["x_um"]*1e-6;step=axis[1]-axis[0]
    xx,yy=np.meshgrid(axis,axis);valid=coarse["valid"]
    plan=finufft.Plan(1,(len(small),len(small)),eps=1e-12,isign=1,nthreads=4)
    plan.setpts(np.ascontiguousarray(xx[valid]*dq),np.ascontiguousarray(yy[valid]*dq))
    angular=[]
    for key in ["E_V_m","H_A_m"]:
        angular.append(np.stack([plan.execute(np.ascontiguousarray(coarse[key][...,i][valid]*step**2)).T for i in range(3)],-1))
    E,H=angular;P=np.stack([H[...,1],-H[...,0]],-1);nt=n[...,:2];s=np.where(mask,kz/K0,1.)
    up=np.zeros_like(E);up[...,:2]=(E[...,:2]+Z0/s[...,None]*(P-nt*np.sum(nt*P,-1)[...,None]))/2
    up[...,2]=-np.sum(nt*up[...,:2],-1)/s;up[~mask]=0
    far=1j*kz[...,None]/(2*np.pi)*np.exp(1j*kz*float(coarse["z_um"])*1e-6)[...,None]*up
    ref=fine["F_V"];weight=np.where(mask,dq*dq/(K0*np.where(mask,kz,1)),0.)
    power=lambda a:float(np.sum(abs(a)**2*weight[...,None])/(2*Z0))
    center=1j*K0/(2*np.pi)*np.exp(1j*K0*1.1e-6)*(E[h,h,:2]+Z0*np.array([H[h,h,1],-H[h,h,0]]))/2
    np.testing.assert_allclose(center,far[h,h,:2],rtol=1e-12,atol=1e-12)
    result={"comparison":"coarse 801 and fine 1601 fields evaluated at identical physical q; coarse direct NUFFT, fine padded FFT",
        "center_complex_relative_change":float(np.linalg.norm(far[h,h]-ref[h,h])/np.linalg.norm(ref[h,h])),
        "coarse_power_W_same_q":power(far),"fine_power_W_same_q":power(ref),
        "power_relative_change_same_q":abs(power(far)/power(ref)-1),
        "lightcone_weighted_complex_relative_change":float(np.sqrt(np.sum(abs(far-ref)**2*weight[...,None])/np.sum(abs(ref)**2*weight[...,None]))),
        "fine_FFT_center_vs_raw_sum_relative_error":float(np.linalg.norm(fine["source_E_angular"][h,h]-np.nansum(load("finite_air_1100nm_fine_EH.npz")["E_V_m"],axis=(0,1))*dx**2)/np.linalg.norm(fine["source_E_angular"][h,h]))}
    assert result["fine_FFT_center_vs_raw_sum_relative_error"]<1e-10
    path=L/"sampling_revision.json";revision=json.loads(path.read_text(encoding="utf8"))
    revision.update(same_q_sampling_check=result,q_step_after_um_inv=float(dq*1e-6),lightcone_points=int(mask.sum()),status="analysis_running")
    path.write_text(json.dumps(revision,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    log("Raw sampling convergence: "+json.dumps(result))


def source_fields():
    """Integrate fitted Bloch sources on native periodic meshes; actual regions on finite FEM."""
    dx,q,small,mask,kz,n=axes();side=len(small);dq=small[1]-small[0]
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    centers=np.array([[p["x"],p["y"]] for p in cfg["bulk_points"]])
    f=load("field_decomposition.npz");coeff=f["grouped_coefficients"]
    plan=finufft.Plan(1,(side,side),eps=1e-12,isign=1,nthreads=4)
    plan.setpts(np.ascontiguousarray(centers[:,0]*dq*1e-6),np.ascontiguousarray(centers[:,1]*dq*1e-6))
    grouped=np.zeros((6,side,side,3),complex);audit=[]
    native_reference=np.zeros((3,3),complex)
    for ik in range(19):
        p=load(f"periodic_k{ik:02d}_mesh_EH.npz");keep=p["material_mask"]
        window=plan.execute(np.ascontiguousarray(f["phase"][:,ik*4])).T
        for band in range(4):
            response=radiation(p["xyz_um"][keep],p["E_normalized"][band,keep],
                               p["weights_half_m3"][keep]*9.89)
            center=response[side//2,side//2,:2]
            exact=p["source_center_native_F_V_per_sqrtJ"][band]
            # Absolute error is assessed on the resulting finite-group center scale.
            error=(center-exact)*window[side//2,side//2]*f["coefficients"][4*ik+band]
            audit.append({"k":ik,"band":band,"finite_center_error_V":float(np.linalg.norm(error)),
                          "cell_center_F_norm_V_per_sqrtJ":float(np.linalg.norm(exact))})
            for group in range(3):
                grouped[group]+=coeff[group,4*ik+band]*window[...,None]*response
                native_reference[group,:2]+=coeff[group,4*ik+band]*window[side//2,side//2]*exact
        log(f"Native-mesh periodic radiation k{ik:02d}/18")
    for region,group in [("cavity",3),("cladding",4),("pad",5)]:
        data=load(f"finite_native_{region}_EH.npz");keep=data["material_mask"]
        field=radiation(data["xyz_um"][keep],data["E_V_m"][keep],
                        data["weights_octant_m3"][keep]*9.89,octant=True)
        grouped[group]=field-grouped[:3].sum(0) if group==3 else field
    np.savez_compressed(D/"group_radiation.npz",F_groups_V=grouped,names=NAMES,
        native_periodic_group_center_F_V=native_reference,fft_size=FFT_SIZE,small_q_m_inv=small)
    write_json("periodic_source_quadrature.json",{"per_basis":audit,
        "max_finite_center_error_V":max(row["finite_center_error_V"] for row in audit),
        "native_periodic_group_center_F_pairs":np.stack([native_reference.real,native_reference.imag],-1).tolist(),
        "rejected_source_files":[f"source_F_group{i}.npz" for i in range(6)],
        "reason":"geometrical cell Gauss rule inadequate for cancellation-sensitive material-source integral"})
    log("All six sources recomputed on independent native FEM quadratures.")


def mesh_audit():
    from comsol_workflow.hex_lattice_utils import unit_cell_corners
    from matplotlib.path import Path as Poly
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    centers=np.array([[p["x"],p["y"]] for key in ["bulk_points","cladding_points"] for p in cfg[key]])
    tree=cKDTree(centers);poly=Poly(unit_cell_corners(float(cfg["a"])))
    records=json.loads((L/"finite_mesh_exports.json").read_text(encoding="utf8"))["parts"]
    regions=[[],[],[]]
    energies=np.zeros((3,6));source=[];pad={k:[] for k in ["xyz_um","weights_octant_m3","E_V_m","H_A_m","material_mask","epsilon_r"]}
    center=np.zeros(3,complex);total_volume=0.
    for row in records:
        d=load(row["path"]);xyz=d["xyz_um"];w=d["weights_octant_m3"]
        _,ids=tree.query(xyz[:,:2]);inside=poly.contains_points(xyz[:,:2]-centers[ids],radius=1e-10)
        region=np.where(inside,np.where(ids<len(cfg["bulk_points"]),0,1),2)
        density=np.column_stack([epsilon_0/4*d["epsilon_r"][:,None]*abs(d["E_V_m"])**2,
                                mu_0/4*abs(d["H_A_m"])**2])
        for i in range(3):
            energies[i]+=np.sum(density[region==i]*w[region==i,None]*8,axis=0)
            regions[i].append({key:d[key][region==i] for key in pad})
        for key in pad:pad[key].append(d[key][region==2])
        use=d["material_mask"]
        source.append((xyz[use],d["E_V_m"][use],w[use]*9.89))
        center[1]+=np.sum(w[use]*9.89*d["E_V_m"][use,1]*np.cos(K0*xyz[use,2]*1e-6))*8*K0*K0/(4*np.pi)
        total_volume+=w.sum()*8
    np.savez_compressed(D/"finite_pad_EH.npz",**{k:np.concatenate(v) for k,v in pad.items()})
    for i,name in enumerate(["cavity","cladding","pad"]):
        np.savez_compressed(D/f"finite_native_{name}_EH.npz",**{k:np.concatenate([part[k] for part in regions[i]]) for k in pad})
    np.savez_compressed(D/"mesh_energy_regions.npz",energy_components_J=energies,names=["cavity","cladding","pad"])
    reference=native();expected=np.array(reference["full_material_center_F_V_pairs"]);expected=expected[:,0]+1j*expected[:,1]
    metrics={"mesh_energy_J":float(energies.sum()),"native_energy_J":reference["U_ref_J"],
        "energy_relative_error":float(abs(energies.sum()/reference["U_ref_J"]-1)),
        "center_source_relative_error":float(np.linalg.norm(center[:2]-expected)/np.linalg.norm(expected)),
        "source_center_V_pairs":np.column_stack([center.real,center.imag]).tolist(),
        "region_energy_components_J":energies.tolist(),"slab_volume_m3":float(total_volume)}
    write_json("mesh_validation.json",metrics);log("Independent full mesh audit: "+json.dumps(metrics))
    path=D/"full_mesh_radiation.npz"
    if not path.exists() or load(path.name).get("fft_size",0)!=FFT_SIZE:
        xyz=np.concatenate([s[0] for s in source]);E=np.concatenate([s[1] for s in source]);w=np.concatenate([s[2] for s in source])
        F=radiation(xyz,E,w,octant=True)
        np.savez_compressed(path,F_V=F,fft_size=FFT_SIZE,small_q_m_inv=axes()[2])

def self_check():
    # Check physical Fourier phase, area normalization, zero padding and ROI recovery.
    rng=np.random.default_rng(17)
    sample=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9));step=.031e-6;size=41
    spectrum=forward(sample,step,size)
    position=(np.arange(9)-4)*step
    x,y=np.meshgrid(position,position)
    for iy,ix in [(20,20),(22,17)]:
        qx=(ix-size//2)*2*np.pi/(size*step);qy=(iy-size//2)*2*np.pi/(size*step)
        expected=np.sum(sample*np.exp(1j*(qx*x+qy*y)))*step**2
        np.testing.assert_allclose(spectrum[iy,ix],expected,rtol=1e-12,atol=1e-26)
    np.testing.assert_allclose(inverse(spectrum,step)[16:25,16:25],sample,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(padded_transform_roi(sample,step,41,17),spectrum[12:29,12:29],rtol=1e-12,atol=1e-26)
    band=spectrum[16:25,16:25];full=np.zeros((41,41),complex);full[16:25,16:25]=band
    np.testing.assert_allclose(padded_transform_roi(band,step,41,17,inverse_transform=True),inverse(full,step)[12:29,12:29],rtol=1e-12,atol=1e-12)
    dx,q,small,mask,kz,n=axes();h=len(small)//2
    rng=np.random.default_rng(81)
    xyz=rng.uniform([0,0,0],[8,7,.1],(70,3));E=rng.normal(size=(70,3))+1j*rng.normal(size=(70,3))
    w=rng.uniform(.1,1,70)*1e-18
    F=radiation(xyz,E,w,octant=True)
    sx=np.array([-1,1,1]);sy=np.array([-1,1,-1]);sz=np.array([1,1,-1])
    for iy,ix in [(h,h),(h+3,h-9),(h+70,h+20)]:
        direction=n[iy,ix];out=np.zeros(3,complex)
        for xsign in [-1,1]:
            for ysign in [-1,1]:
                for zsign in [-1,1]:
                    point=xyz*np.array([xsign,ysign,zsign])
                    e=E*(sx if xsign<0 else 1)*(sy if ysign<0 else 1)*(sz if zsign<0 else 1)
                    out+=np.sum(e*(w*np.exp(1j*K0*1e-6*(point@direction)))[:,None],axis=0)
        direct=(out-direction*np.dot(direction,out))*K0*K0/(4*np.pi)
        np.testing.assert_allclose(F[iy,ix],direct,atol=1e-13,rtol=1e-9)
    log("PASS: NUFFT sign, vector parity, z moments and projection vs direct full-volume sum.")


def finalize_arrays():
    """Coherent budgets and full-lightcone surface fields of every source group."""
    dx,q,small,mask,kz,n=axes();h=len(small)//2
    direct=load("direct_radiative_fields.npz");groups=load("group_radiation.npz")["F_groups_V"]
    full=load("full_mesh_radiation.npz")["F_V"]
    summed=groups.sum(0);U=native()["U_ref_J"]
    dq=small[1]-small[0];dOmega=np.where(mask,dq*dq/(K0*np.where(mask,kz,1)),0.)
    power=lambda f:float(np.sum(abs(f)**2*dOmega[...,None])/(2*Z0))
    nearE=[];nearH=[]
    for f in groups:
        A=np.zeros_like(f)
        A[mask]=2*np.pi/(1j*kz[mask,None])*f[mask]*np.exp(-1j*kz[mask]*SZ)[:,None]
        H=np.cross(n,A)/Z0
        outputs=[]
        for spec in [A,H]:
            pieces=[]
            for i in range(3):
                pieces.append(padded_transform_roi(spec[...,i],dx,len(q),len(direct["x_um"]),inverse_transform=True))
            outputs.append(np.stack(pieces,-1))
        nearE.append(outputs[0]);nearH.append(outputs[1])
    np.savez_compressed(D/"radiative_fields.npz",F_groups_V=groups,F_direct_V=direct["F_V"],
        F_full_mesh_V=full,F_group_sum_V=summed,E_groups_surface_V_m=np.array(nearE),
        H_groups_surface_A_m=np.array(nearH),small_q_m_inv=small,x_um=direct["x_um"],
        lightcone_mask=mask,frequency_Hz=FREQ,z_surface_m=SZ,U_ref_J=U,names=NAMES,fft_size=FFT_SIZE)
    ref=direct["F_V"];den=np.linalg.norm(ref[h,h])
    metrics={"COMSOL_started":True,"legacy_derived_results_reused":False,
        "U_ref_J":U,"frequency_Hz":FREQ,"NA":1.,"slab_z_nm":0.,"surface_z_nm":SZ*1e9,
        "mesh_source_to_air_center_complex_error":float(np.linalg.norm(full[h,h]-ref[h,h])/den),
        "mesh_source_to_air_center_power_error":float(abs(np.sum(abs(full[h,h])**2)/np.sum(abs(ref[h,h])**2)-1)),
        "mesh_source_to_air_total_power_error":abs(power(full)/power(ref)-1),
        "mesh_source_to_air_lightcone_weighted_error":float(np.sqrt(np.sum(abs(full-ref)**2*dOmega[...,None])/np.sum(abs(ref)**2*dOmega[...,None]))),
        "cell_groups_to_mesh_center_complex_error":float(np.linalg.norm(summed[h,h]-full[h,h])/np.linalg.norm(full[h,h])),
        "cell_groups_to_mesh_total_power_error":abs(power(summed)/power(full)-1),
        "direct_power_W":power(ref),"source_power_W":power(full),"groups_power_W":power(summed),
        "center_Ex_over_Ey_power":float(abs(ref[h,h,0]/ref[h,h,1])**2),
        "direct_center_y_W_sr":float(abs(ref[h,h,1])**2/(2*Z0)),
        "source_center_y_W_sr":float(abs(full[h,h,1])**2/(2*Z0)),
        "power_per_U_s_inv":power(ref)/U,
        "source_group_sum_center_y_W_sr":float(abs(summed[h,h,1])**2/(2*Z0))}
    metrics["independent_physical_closure_passed"]=all([
        metrics["mesh_source_to_air_center_complex_error"]<=.02,
        metrics["mesh_source_to_air_center_power_error"]<=.05,
        metrics["mesh_source_to_air_total_power_error"]<=.05])
    metrics["mechanism"]="independent closure passed" if metrics["independent_physical_closure_passed"] else "independent closure not passed; attribution is provisional"
    write_json("validation.json",metrics)
    rows=budget_rows(NAMES,groups[:,h,h,:2],1/(2*Z0),"volume_source")
    for channel,index in [("Ex",0),("Ey",1)]:
        for name,f in [("direct_air",ref),("independent_full_mesh",full)]:
            rows.append(dict(scope=name,channel=channel,kind="reference",first=name,second=name,
                             value=float(abs(f[h,h,index])**2/(2*Z0))))
    pd.DataFrame(rows).assign(unit="W/sr").to_csv(table_path(L, "center_radiation_budget.csv"),index=False)
    cross=np.array([[np.sum(a.conj()*b*dOmega[...,None])/(2*Z0) for b in groups] for a in groups])
    np.savez_compressed(D/"radiation_interference.npz",power_gram_W=cross)
    r=[]
    for i in range(6):
        r.append(dict(first=NAMES[i],second=NAMES[i],kind="self",power_W=cross[i,i].real))
        for j in range(i):r.append(dict(first=NAMES[j],second=NAMES[i],kind="interference",power_W=2*cross[j,i].real))
    pd.DataFrame(r).to_csv(table_path(L, "radiation_power_budget.csv"),index=False)
    log("Independent closure: "+json.dumps(metrics))

def slab_groups():
    from comsol_workflow.hex_lattice_utils import unit_cell_corners
    from matplotlib.path import Path as Poly
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    actual=load("finite_slab_z0_fine_EH.npz");fitdata=load("field_decomposition.npz")
    axis=actual["x_um"][::2];xx,yy=np.meshgrid(axis,axis)
    xy=np.column_stack([xx.ravel(),yy.ravel()])
    centers=np.array([[p["x"],p["y"]] for key in ["bulk_points","cladding_points"] for p in cfg[key]])
    _,ids=cKDTree(centers).query(xy)
    local=xy-centers[ids]
    inside=Poly(unit_cell_corners(cfg["a"])).contains_points(local,radius=1e-10)
    reg=np.where(inside,np.where(ids<len(cfg["bulk_points"]),0,1),2)
    valid=actual["valid"][::2,::2].ravel();cavity=(reg==0)&valid
    V=np.concatenate([actual["E_V_m"][::2,::2],actual["H_A_m"][::2,::2]],-1).reshape(-1,6)
    fields=np.zeros((6,len(xy),6),complex)
    for ik in range(19):
        p=load(f"periodic_k{ik:02d}_z0_fine_EH.npz")
        normalized=np.concatenate([p["E_center_V_m"],p["H_center_A_m"]],-1)/np.sqrt(p["cell_EM_energy_J"])[:,None,None]
        interpolation=LinearNDInterpolator(p["center_xyz_um"][:,:2],normalized.transpose(1,0,2),fill_value=np.nan)
        local_field=interpolation(local[cavity])
        missing=~np.isfinite(local_field).all((1,2))
        if missing.any():
            # Only within half a source pixel of the hexagonal boundary: record this
            # nearest-node display extension; no interpolation enters volume fitting.
            _,near=cKDTree(p["center_xyz_um"][:,:2]).query(local[cavity][missing])
            local_field[missing]=normalized[:,near].transpose(1,0,2)
        phi=np.exp(-1j*centers[ids[cavity]]@p["k_um_inv"])
        for group in range(3):
            coeff=fitdata["grouped_coefficients"][group,4*ik:4*ik+4]
            fields[group,cavity]+=np.einsum("cbv,b,c->cv",local_field,coeff,phi)
    fields[3,cavity]=V[cavity]-fields[:3,cavity].sum(0)
    for region,group in [(1,4),(2,5)]:
        select=(reg==region)&valid;fields[group,select]=V[select]
    fields[:,~valid]=np.nan
    epsilon=actual["epsilon_r"][::2,::2].ravel()
    energies=epsilon_0/4*epsilon[None]*np.sum(abs(fields[...,:3])**2,-1)+mu_0/4*np.sum(abs(fields[...,3:])**2,-1)
    np.savez_compressed(D/"group_slab_z0.npz",energy_density_J_m3=energies.reshape(6,len(axis),len(axis)),
        x_um=axis,region_map=reg.reshape(len(axis),len(axis)),valid=valid.reshape(len(axis),len(axis)),
        display_method="actual finite E/H at z=0; periodic z=0 linear interpolation, boundary nearest source node only for display",
        reconstruction_relative_error=float(np.linalg.norm((fields.sum(0)-V)[valid])/np.linalg.norm(V[valid])))
    log("All six actual z=0 contribution field maps prepared.")

def figure_panel(fig,ax,values,axis,title,label,maximum,logarithmic=False,minimum=None):
    norm=LogNorm(maximum*1e-7 if minimum is None else minimum,maximum) if logarithmic else plt.Normalize(0,maximum)
    display=np.where(np.isfinite(values)&(values==0),norm.vmin,values) if logarithmic else values
    im=ax.imshow(display,origin="lower",extent=[axis[0],axis[-1]]*2,cmap="magma",norm=norm)
    ax.set(title=title,aspect="equal",xlabel=r"$x$ ($\mu$m)",ylabel=r"$y$ ($\mu$m)")
    colorbar(fig,ax,im,label)
    return im

def full_bz_surface_spectrum(view):
    """Complex upward angular spectrum, including decaying evanescent waves (exp(+iwt))."""
    dx,_,_,_,_,_=axes()
    size=2*FFT_SIZE-1
    q=(np.arange(size)-size//2)*2*np.pi/(size*dx)
    h=int(np.ceil(view*K0/(q[1]-q[0])))
    qs=q[len(q)//2-h:len(q)//2+h+1]
    d=load("finite_air_1100nm_fine_EH.npz")
    E=np.stack([padded_transform_roi(np.nan_to_num(d["E_V_m"][...,j]),dx,size,len(qs)) for j in range(3)],-1)
    H=np.stack([padded_transform_roi(np.nan_to_num(d["H_A_m"][...,j]),dx,size,len(qs)) for j in range(3)],-1)
    qx,qy=np.meshgrid(qs,qs)
    # Upgoing exp(-ikz*z) must decay outside the cone: Im(kz)<0.
    kz=np.sqrt(np.maximum(K0*K0-qx*qx-qy*qy,0))-1j*np.sqrt(np.maximum(qx*qx+qy*qy-K0*K0,0))
    nt=np.stack([qx,qy],-1)/K0;s=kz/K0
    assert np.min(abs(s))>1e-10, "Angular grid intersects the grazing singularity"
    P=np.stack([H[...,1],-H[...,0]],-1)
    v=Z0/s[...,None]*(P-nt*np.sum(nt*P,-1)[...,None])
    up=np.empty_like(E);up[...,:2]=(E[...,:2]+v)/2
    up[...,2]=-np.sum(nt*up[...,:2],-1)/s
    up*=np.exp(1j*kz*(float(d["z_um"])*1e-6-SZ))[...,None]
    direct=load("direct_radiative_fields.npz");old=direct["A_up_surface"]
    # Odd centered 16001 grid preserves zero exactly; other bins differ from 8001.
    np.testing.assert_allclose(up[len(qs)//2,len(qs)//2],old[len(old)//2,len(old)//2],rtol=1e-10,atol=1e-20)
    # Independent direct Fourier sums check off-center refined samples and phase.
    coord=(np.arange(d["E_V_m"].shape[0])-d["E_V_m"].shape[0]//2)*dx
    for iy,ix in [(len(qs)//2+7,len(qs)//2-11),(len(qs)//2-17,len(qs)//2+5)]:
        weight=np.exp(1j*(qs[iy]*coord[:,None]+qs[ix]*coord[None,:]))*dx**2
        for source,actual in [(d["E_V_m"],E),(d["H_A_m"],H)]:
            expected=np.sum(np.nan_to_num(source)*weight[...,None],axis=(0,1))
            np.testing.assert_allclose(actual[iy,ix],expected,rtol=1e-8,atol=1e-18)
    write_json("Fig1_angular_resolution.json",{"fft_size":size,"previous_fft_size":FFT_SIZE,"sampling_density_ratio":size/FFT_SIZE,"delta_k_m_inv":float(q[1]-q[0]),"center_angle_step_deg":float(np.rad2deg(np.arcsin((q[1]-q[0])/K0))),"checks":"center agrees with original upward spectrum; two refined off-center E/H bins checked against independent direct Fourier sums","source_sampling_changed":False})
    np.savez_compressed(D/"Fig1_full_BZ_angular_spectrum_2x.npz",q_m_inv=qs,E_up_surface_angular=up,
        kz_m_inv=kz,z_surface_m=SZ,z_input_m=float(d["z_um"])*1e-6,fft_size=size,
        source="finite_air_1100nm_fine_EH.npz",time_convention="exp(+i omega t)")
    return qs/K0,abs(up)**2


def smooth_material_display(values,mask,sigma_pixels):
    """Normalized Gaussian filter: only same-material values contribute."""
    from scipy.ndimage import gaussian_filter
    weight=gaussian_filter(mask.astype(float),sigma_pixels,mode="constant",cval=0.)
    numerator=gaussian_filter(np.where(mask[...,None],values,0.),(*sigma_pixels,0.),mode="constant",cval=0.)
    result=np.divide(numerator,weight[...,None],out=np.zeros_like(numerator),where=weight[...,None]>1e-12)
    assert np.isfinite(result[mask]).all()
    return result


def periodic_material_layers(points,fields,geometry):
    """Display-only interpolation within each material, with exact polygon clipping."""
    from scipy.interpolate import RBFInterpolator
    from matplotlib.path import Path as PolygonPath
    regions=np.zeros(len(points),int)
    for i,vertices in enumerate(geometry["holes"],1):
        regions[PolygonPath(vertices).contains_points(points)]=i
    layers=[];checks=[]
    packed=np.concatenate([fields.real,fields.imag],axis=-1)
    for i,vertices in enumerate([geometry["outer"],*geometry["holes"]]):
        vertices=np.asarray(vertices,float);inside=regions==i
        assert inside.sum()>=6
        interpolation=RBFInterpolator(points[inside],packed[inside],neighbors=min(32,int(inside.sum())),kernel="thin_plate_spline",smoothing=0.)
        error=float(np.max(abs(interpolation(points[inside])-packed[inside]))/np.max(abs(packed[inside])))
        assert error<1e-9
        size=641 if i==0 else 241
        x=np.linspace(vertices[:,0].min(),vertices[:,0].max(),size)
        y=np.linspace(vertices[:,1].min(),vertices[:,1].max(),size)
        xx,yy=np.meshgrid(x,y);query=np.column_stack([xx.ravel(),yy.ravel()])
        values=interpolation(query).reshape(size,size,-1)
        values=values[...,:fields.shape[1]]+1j*values[...,fields.shape[1]:]
        layers.append((vertices,[x[0],x[-1],y[0],y[-1]],values))
        checks.append({"region":i,"original_samples":int(inside.sum()),"display_grid":[size,size],"source_reconstruction_relative_error":error})
    return layers,checks


def figure2():
    """Existing near-field maps plus a spatially integrated linear-analyzer response."""
    plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    d=load("direct_radiative_fields.npz");axis=d["x_um"];E=d["E_raw_surface"]
    U=native()["U_ref_J"];nu=epsilon_0/4*abs(E)**2/U
    dx=float(axis[1]-axis[0])*1e-6
    ex,ey=E[...,0],E[...,1]
    xx=float(np.sum(abs(ex)**2)*dx**2);yy=float(np.sum(abs(ey)**2)*dx**2)
    xy=complex(np.sum(ex*np.conj(ey))*dx**2)
    angles=np.linspace(0.,2*np.pi,721)
    intensity=yy*np.cos(angles)**2+xx*np.sin(angles)**2-2*xy.real*np.sin(angles)*np.cos(angles)
    eigen=np.linalg.eigvalsh([[xx,xy.real],[xy.real,yy]])
    maximum=float(eigen[-1]);minimum=float(max(0.,eigen[0]));assert maximum>0
    contrast=(maximum-minimum)/(maximum+minimum)
    np.testing.assert_allclose(intensity[[0,180,360]],[yy,xx,yy],rtol=1e-12,atol=maximum*1e-14)
    test_angle=.731
    direct=float(np.sum(abs(ey*np.cos(test_angle)-ex*np.sin(test_angle))**2)*dx**2)
    calculated=yy*np.cos(test_angle)**2+xx*np.sin(test_angle)**2-2*xy.real*np.sin(test_angle)*np.cos(test_angle)
    np.testing.assert_allclose(direct,calculated,rtol=1e-12)
    pd.DataFrame({"angle_deg_from_y":np.rad2deg(angles),"integrated_squared_field_V2":intensity,
                  "normalized_intensity":intensity/maximum}).to_csv(table_path(L, "Fig2_polarization_radar.csv"),index=False)
    write_json("Fig2_polarization_radar.json",{"source":"direct_radiative_fields.npz:E_raw_surface",
        "plane_z_nm":100,"NA":1,"zero_angle":"+y; counterclockwise toward -x",
        "formula":"I(theta)=integral |Ey*cos(theta)-Ex*sin(theta)|^2 dxdy",
        "Ixx_V2":xx,"Iyy_V2":yy,"Ixy_real_V2":xy.real,"Ixy_imag_V2":xy.imag,
        "normalization_max_V2":maximum,"linear_analyzer_contrast":contrast,
        "interpretation":"spatially integrated surface-field linear-analyzer response; not full Stokes degree of polarization or calibrated far-field power",
        "map_order":["Full electric field","Polarization radar","Ey","Ex"],"layout":"2x2"})
    fig=plt.figure(figsize=(12,10),layout="constrained");grid=fig.add_gridspec(2,2)
    maps=[(fig.add_subplot(grid[0,0]),nu.sum(-1),"Full electric field"),
          (fig.add_subplot(grid[1,0]),nu[...,1],"Ey"),(fig.add_subplot(grid[1,1]),nu[...,0],"Ex")]
    for ax,value,label in maps:
        figure_panel(fig,ax,value,axis,label,r"$u_e/U_{\rm ref}$ (m$^{-3}$)",np.nanmax(nu.sum(-1)))
    polar=fig.add_subplot(grid[0,1],projection="polar")
    polar.plot(angles,intensity/maximum,color="#536b91",linewidth=1.7)
    polar.set_theta_zero_location("N");polar.set_ylim(0,1.05)
    polar.set_yticks([.25,.5,.75,1.]);polar.set_rlabel_position(45)
    polar.grid(color=".82",linewidth=.6)
    polar.set_title(f"Polarization radar | linear contrast={contrast:.3f}",pad=18)
    polar.set_xlabel("Analyzer angle (deg); 0 = y\nNormalized integrated intensity",labelpad=16)
    fig.suptitle("Radiative near field | NA=1, backpropagated to z=100 nm",fontsize=14)
    save(fig,"Fig.2_radiative_nearfield",dpi=FIGURE_DPI)


def figure4():
    """Cavity and cladding source contributions; coherent within each region."""
    from matplotlib.path import Path as PolygonPath
    from comsol_workflow.hex_lattice_utils import unit_cell_corners
    plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    slab=load("finite_slab_z0_fine_EH.npz");sa=slab["x_um"]
    energy=(epsilon_0*slab["epsilon_r"]*np.sum(abs(slab["E_V_m"])**2,-1)+mu_0*np.sum(abs(slab["H_A_m"])**2,-1))/4
    slab_peak=float(np.nanmax(energy));valid=slab["valid"]
    del slab
    centers=np.array([[p["x"],p["y"]] for key in ["bulk_points","cladding_points"] for p in cfg[key]])
    tree=cKDTree(centers);cell=PolygonPath(unit_cell_corners(cfg["a"]))
    region=np.full(energy.shape,2,dtype=np.int8)
    for start in range(0,len(sa),128):
        xx,yy=np.meshgrid(sa,sa[start:start+128]);xy=np.column_stack([xx.ravel(),yy.ravel()])
        _,ids=tree.query(xy);inside=cell.contains_points(xy-centers[ids],radius=1e-10)
        region[start:start+128]=np.where(inside,np.where(ids<len(cfg["bulk_points"]),0,1),2).reshape(xx.shape)
    slabmaps=np.array([np.where(valid,energy/slab_peak,np.nan),*[np.where(valid,np.where(region==i,energy/slab_peak,0.),np.nan) for i in (0,1)]])
    with np.load(D/"group_slab_z0.npz") as old:
        np.testing.assert_array_equal(region[::2,::2],old["region_map"])
    assert not np.any((region==0)&(region==1))
    source=load("group_radiation.npz");g=source["F_groups_V"]
    # Residual = actual cavity source minus fitted groups, so this is actual cavity radiation.
    far=np.array([g.sum(0),g[:4].sum(0),g[4]])
    with np.load(D/"full_mesh_radiation.npz") as full:
        closure=float(np.linalg.norm(far[0]-full["F_V"])/np.linalg.norm(full["F_V"]))
    assert closure<1e-10
    dx,q,small,mask,kz,n=axes();axis=load("direct_radiative_fields.npz")["x_um"]
    spectra=np.zeros_like(far)
    spectra[:,mask]=2*np.pi/(1j*kz[mask,None])*far[:,mask]*np.exp(-1j*kz[mask]*SZ)[None,:,None]
    near=[];near_total=[]
    for spectrum in spectra:
        magnetic=np.cross(n,spectrum)/Z0
        E=np.stack([padded_transform_roi(spectrum[...,j],dx,len(q),len(axis),inverse_transform=True) for j in range(3)],-1)
        H=np.stack([padded_transform_roi(magnetic[...,j],dx,len(q),len(axis),inverse_transform=True) for j in range(3)],-1)
        near.append(E)
        near_total.append((epsilon_0*np.sum(abs(E)**2,-1)+mu_0*np.sum(abs(H)**2,-1))/4)
    near=np.array(near);near_total=np.array(near_total);near_peak=float(np.max(near_total[1:]))
    near_electric=epsilon_0/4*abs(near[...,:2])**2/near_peak
    intensity=np.sum(abs(spectra)**2,-1);h=len(small)//2
    center_reference=float(np.max(intensity[1:,h,h]));angular=intensity/center_reference
    normalized_q=small/K0;step=float(normalized_q[1]-normalized_q[0])
    indices=np.flatnonzero(abs(normalized_q)<=np.sin(np.deg2rad(10))+step)
    selected=normalized_q[indices];angle_edges=np.rad2deg(np.arcsin(np.r_[selected-step/2,selected[-1]+step/2]))
    assert angle_edges[0]<-10 and angle_edges[-1]>10
    names=["Full field","Cavity only","Cladding only"]
    validation=json.loads((L/"validation.json").read_text(encoding="utf8"))
    np.savez_compressed(D/"Fig4_region_chain.npz",names=names,slab_x_um=sa,slab_normalized_energy=slabmaps,
        slab_region=region,near_x_um=axis,near_E_V_m=near,near_full_EM_J_m3=near_total,
        near_normalized_electric=near_electric,F_regions_V=far,A_regions_surface=spectra,
        q_m_inv=small,lightcone_mask=mask,angular_normalized_intensity=angular,
        slab_reference_J_m3=slab_peak,near_reference_J_m3=near_peak,angular_center_reference_V2_m2=center_reference)
    peaks=[]
    fig,axs=plt.subplots(3,4,figsize=(23,15),layout="constrained")
    boundary=np.asarray(cfg["inner_boundary"]);boundary=np.vstack([boundary,boundary[0]])
    for col,name in enumerate(names):
        images=[slabmaps[col],near_electric[col,...,0],near_electric[col,...,1]]
        maxima=[1.,.01,.03] if col==0 else [1.,.3,1.]
        if col==0:
            assert all(.1<=float(np.nanmax(v))/limit<=1. for v,limit in zip(images[1:],maxima[1:])), "Full-field near maps need a visible, non-saturating color range"
        labels=[r"$u_{EM}/u_{EM,slab,max}$ (1)",r"$u_{e,x}/u_{EM,near,max}$ (1)",r"$u_{e,y}/u_{EM,near,max}$ (1)"]
        titles=["Slab EM, z=0 nm","Ex radiative near field","Ey radiative near field"]
        for row,(value,limit,label,title) in enumerate(zip(images,maxima,labels,titles)):
            figure_panel(fig,axs[col,row],value,sa if row==0 else axis,f"{name} | {title}",label,limit)
            if row==0:axs[col,row].plot(boundary[:,0],boundary[:,1],color="white",linewidth=.8)
            peaks.append(dict(region=name,stage=title,maximum=float(np.nanmax(value)),color_max=limit,
                              fraction_above_color_max=float(np.mean(value[np.isfinite(value)]>limit))))
        value=angular[col][np.ix_(indices,indices)]
        im=axs[col,3].pcolormesh(angle_edges,angle_edges,value,shading="flat",cmap="magma",vmin=0,vmax=1,rasterized=True)
        axs[col,3].set(title=f"{name} | Far field intensity: Full E",aspect="equal",xlim=(-10,10),ylim=(-10,10),
            xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10],xlabel=r"$\theta_x$ (°)",ylabel=r"$\theta_y$ (°)")
        colorbar(fig,axs[col,3],im,r"$|\widetilde{\mathbf{E}}|^2/I_{center,ref}$ (1)")
        peaks.append(dict(region=name,stage="Angular spectrum",maximum=float(value.max()),color_max=1.,fraction_above_color_max=float(np.mean(value>1))))
    status="Independent source-to-air closure not passed; source attribution provisional" if not validation["independent_physical_closure_passed"] else "Independent source-to-air closure passed"
    fig.suptitle("Mode19 | Full field / cavity / cladding | Shared normalization\nFull field includes pad and coherent interference; regional rows show self terms\n"+status,fontsize=13)
    save(fig,"Fig.4_cavity_cladding_chain",dpi=FIGURE_DPI)
    write_json("Fig4_region_chain.json",{"layout":"3 rows x 4 columns","rows":names,
        "columns":["slab EM","Ex near field","Ey near field","angular spectrum (Fig.1 convention)"],
        "source_definition":"full=sum all six groups coherently, including pad; cavity=sum groups 0:4; cladding=group 4",
        "slab_reference_J_m3":slab_peak,"near_reference_J_m3":near_peak,
        "angular_center_reference_V2_m2":center_reference,"normalization":"shared across all three rows; previous cavity/cladding reference denominators retained",
        "regional_column_color_limits":[[0,1],[0,.3],[0,1],[0,1]],"full_field_column_color_limits":[[0,1],[0,.01],[0,.03],[0,1]],"angular_axes_deg":[-10,10],
        "source_sum_including_pad_relative_error":closure,"peaks":peaks,
        "independent_physical_closure_passed":validation["independent_physical_closure_passed"],
        "center_power_closure_error":validation["mesh_source_to_air_center_power_error"],
        "NA_power_closure_error":validation["mesh_source_to_air_total_power_error"],
        "COMSOL_started":False,"model_modified":False})


def figures(*,fig2_only=False,fig1_only=False,hz_part="imag"):
    if fig2_only:
        figure2()
        return
    assert hz_part in ("real","imag")
    plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    U=native()["U_ref_J"]
    cfg=json.loads((CONFIG/"config.json").read_text(encoding="utf8"))
    boundary=np.asarray(cfg["inner_boundary"],float)
    boundary=np.vstack([boundary,boundary[0]])
    from matplotlib.patheffects import Stroke,Normal
    boundary_style=dict(color="white",linewidth=.8,zorder=5,
        path_effects=[Stroke(linewidth=1.5,foreground="#333333"),Normal()],
        label="cavity-cladding boundary")
    slab=load("finite_slab_z0_fine_EH.npz");sa=slab["x_um"]
    ue=epsilon_0/4*slab["epsilon_r"][...,None]*abs(slab["E_V_m"])**2
    um=mu_0/4*np.sum(abs(slab["H_A_m"])**2,-1)
    total=ue.sum(-1)+um
    slab_peak=float(np.nanmax(total))
    assert np.isfinite(slab_peak) and slab_peak>0
    slabmaps=[ue[...,0]/slab_peak,ue[...,1]/slab_peak,total/slab_peak]
    slab_labels=[r"$u_{e,x}/u_{\mathrm{EM,max}}$ (1)",r"$u_{e,y}/u_{\mathrm{EM,max}}$ (1)",r"$u_{\mathrm{EM}}/u_{\mathrm{EM,max}}$ (1)"]
    slab_max=1.
    write_json("slab_normalization.json",{"definition":"u_EM,max = max_xy[(epsilon0 epsilon_r |E|^2 + mu0 |H|^2)/4] at z=0",
        "denominator_J_m3":slab_peak,"normalized_peaks":[float(np.nanmax(v)) for v in slabmaps],
        "source":"01_results/finite_slab_z0_fine_EH.npz","support":"full finite z=0 footprint including cavity, cladding, pad and hole air",
        "scope":["Fig.1 finite slab column"],"dimensionless":True,"color_limits":[1e-3,1.],
        "Fig1_slab":{"scale":"linear","color_limits":[0.,1.]},
        "local_energy_fraction":False,"COMSOL_started_this_revision":False})
    if not fig2_only:
        d=load("direct_radiative_fields.npz");axis=d["x_um"]
        nu=epsilon_0/4*abs(d["E_raw_surface"])**2/U
        upue=epsilon_0/4*np.sum(abs(d["E_up_surface"])**2,-1)/U
        upum=mu_0/4*np.sum(abs(d["H_up_surface"])**2,-1)/U
        near=[nu[...,0],nu[...,1],upue+upum]
        far=abs(d["F_V"])**2/(2*Z0*U);far[~d["lightcone_mask"]]=np.nan
        q=d["small_q_m_inv"]/K0
        fmap=[far[...,0],far[...,1],far.sum(-1)]
        p=load("periodic_k00_z0_fine_EH.npz");pe=p["center_xyz_um"][:,:2]
        # Material-resolved periodic shape map uses the same actual cavity geometry.
        manifest=json.loads((REF/"99_config/manifest.json").read_text(encoding="utf8"))
        geometry=json.loads(Path(manifest["geometry_source"]).read_text(encoding="utf8"))
        from matplotlib.path import Path as Poly
        hole=np.zeros(len(pe),bool)
        for hs in geometry["holes"]:hole|=Poly(hs).contains_points(pe[:,:2])
        eps=np.where(hole,1.,10.89);e=p["E_center_V_m"][1];hh=p["H_center_A_m"][1]
        pu=epsilon_0/4*eps[:,None]*abs(e)**2;pm=mu_0/4*np.sum(abs(hh)**2,-1)
        puvals=np.column_stack([pu[:,0],pu[:,1],pu.sum(-1)+pm]);pref=puvals[:,2].max()
        electric_peak=float(np.max(np.linalg.norm(e,axis=-1)))
        magnetic_peak=float(np.max(np.linalg.norm(hh,axis=-1)))
        phase_index=int(np.argmax(abs(e[:,1])))
        phase_angle=float(np.angle(e[phase_index,1]))
        aligned=e*np.exp(-1j*phase_angle)
        periodic_display=np.column_stack([aligned[:,0].real/electric_peak,aligned[:,1].real/electric_peak,puvals[:,2]/pref])
        assert np.max(abs(periodic_display[:,:2]))<=1.+1e-12
        assert aligned[phase_index,1].real>0 and abs(aligned[phase_index,1].imag)<electric_peak*1e-12
        write_json("Fig1_periodic_real_field.json",{"normalization":"max sqrt(|Ex|^2+|Ey|^2+|Ez|^2) on periodic z=0 cell",
            "denominator_V_m":electric_peak,"Hz_denominator_A_m":magnetic_peak,"Hz_display_part":hz_part,"phase_angle_rad":phase_angle,"rotation":"exp(-i phase_angle_rad)",
            "phase_reference":"largest |Ey| sample made positive real; same phase for Ex, Ey and Hz",
            "phase_reference_xyz_um":p["center_xyz_um"][phase_index].tolist(),"limits":[-1.,1.],"colormap":"RdBu_r",
            "Full_EM":"unchanged: u_EM / max(u_EM); magma 0..1","source":"periodic_k00_z0_fine_EH.npz"})
        from comsol_workflow.finite_lattice_fourier import first_bz_polygon
        cell=np.asarray(geometry["outer"],float)
        cell=np.vstack([cell,cell[0]])
        bz=first_bz_polygon(float(cfg["a"]))*1e6/K0
        bz=np.vstack([bz,bz[0]])
        view=1.08*max(1.,float(np.max(abs(bz))))
        layers,interpolation_checks=periodic_material_layers(pe,np.concatenate([e,hh],axis=-1),geometry)
        display_layers=[];saved={}
        sample_step=[float(np.median(np.diff(np.unique(np.round(pe[:,j],12))))) for j in (0,1)]
        smoothing_sigma_um=.75*float(np.sqrt(np.prod(sample_step)))
        for i,(vertices,extent,fields) in enumerate(layers):
            er=10.89 if i==0 else 1.
            field_e=fields[...,:3];field_h=fields[...,3:]
            energy=(epsilon_0*er*np.sum(abs(field_e)**2,-1)+mu_0*np.sum(abs(field_h)**2,-1))/4
            values=np.stack([(field_e[...,0]*np.exp(-1j*phase_angle)).real/electric_peak,
                             (field_e[...,1]*np.exp(-1j*phase_angle)).real/electric_peak,energy/pref,
                             getattr(field_h[...,2]*np.exp(-1j*phase_angle),hz_part)/magnetic_peak],-1)
            gx=np.linspace(extent[0],extent[1],values.shape[1]);gy=np.linspace(extent[2],extent[3],values.shape[0])
            xx,yy=np.meshgrid(gx,gy);xy=np.column_stack([xx.ravel(),yy.ravel()])
            mask=Poly(vertices).contains_points(xy,radius=1e-12)
            if i==0:
                for hs in geometry["holes"]:mask&=~Poly(hs).contains_points(xy)
            mask=mask.reshape(values.shape[:2])
            sigma=(smoothing_sigma_um/(gy[1]-gy[0]),smoothing_sigma_um/(gx[1]-gx[0]))
            before=values
            values=smooth_material_display(before,mask,sigma)
            interpolation_checks[i]["smoothing_sigma_um"]=smoothing_sigma_um
            interpolation_checks[i]["display_relative_L2_change"]=float(np.linalg.norm((values-before)[mask])/np.linalg.norm(before[mask]))
            display_layers.append((vertices,extent,values))
            saved[f"region_{i}_values"]=values;saved[f"region_{i}_extent_um"]=extent
        np.savez_compressed(D/"Fig1_periodic_material_display.npz",**saved)
        write_json("Fig1_periodic_interpolation.json",{"method":"per-material local thin-plate RBF, then material-masked normalized Gaussian (sigma=0.75 native grid step); exact polygon clips",
            "checks":interpolation_checks,"normalization_and_phase":"unchanged from original samples",
            "boundary_values":"one-sided extrapolation from same-material samples; no new FEM evaluation",
            "purpose":"visualization only; not used by energy or radiation calculations","COMSOL_started":False})
        near_peak=float(np.nanmax(upue+upum))
        near_display=[epsilon_0/4*abs(d["E_up_surface"][...,j])**2/U/near_peak for j in (0,1)]
        near_display.append((upue+upum)/near_peak)
        spectrum_q,spectrum=full_bz_surface_spectrum(view)
        angular_data=load("Fig1_full_BZ_angular_spectrum_2x.npz")
        np.testing.assert_allclose(angular_data["q_m_inv"]/K0,spectrum_q)
        sqx,sqy=np.meshgrid(spectrum_q,spectrum_q)
        angular_E=angular_data["E_up_surface_angular"]
        nz=np.sqrt(np.maximum(1-sqx**2-sqy**2,0))-1j*np.sqrt(np.maximum(sqx**2+sqy**2-1,0))
        angular_H=np.cross(np.stack([sqx,sqy,nz],-1),angular_E)/Z0
        angular_ue=epsilon_0/4*abs(angular_E)**2
        angular_uh=mu_0/4*abs(angular_H)**2
        angular_em=angular_ue.sum(-1)+angular_uh.sum(-1)
        center=len(spectrum_q)//2
        far_center=float(angular_em[center,center])
        assert np.isfinite(far_center) and far_center>0
        propagating=sqx**2+sqy**2<1
        np.testing.assert_allclose(angular_ue.sum(-1)[propagating],angular_uh.sum(-1)[propagating],rtol=1e-9,atol=1e-30)
        far_display=[angular_ue[...,0]/far_center,angular_ue[...,1]/far_center,angular_em/far_center,angular_uh[...,2]/far_center]
        angular_Hz=angular_H[...,2]
        np.savez_compressed(D/"Fig1_EM_angular_spectrum.npz",q_m_inv=angular_data["q_m_inv"],u_e_angular=angular_ue,u_h_angular=angular_uh,u_EM_angular=angular_em,normalization_EM_center=far_center)
        near_display.append(mu_0/4*abs(d["H_up_surface"][...,2])**2/U/near_peak)
        fig1_slabmaps=slabmaps+[mu_0/4*abs(slab["H_A_m"][...,2])**2/slab_peak]
        fig1_labels=slab_labels+[r"$u_{h,z}/u_{\mathrm{EM,max}}$ (1)"]
        assert abs(angular_Hz[center,center])==0
        write_json("Fig1_Hz.json",{"periodic_part":hz_part,"same_phase_as_Ex_Ey":True,
            "periodic_magnetic_peak_A_m":magnetic_peak,"periodic_Hz_real_peak_normalized":float(np.max(abs((hh[:,2]*np.exp(-1j*phase_angle)).real))/magnetic_peak),
            "slab_reference_J_m3":slab_peak,"near_reference_J_m3":near_peak*U,
            "angular_reference_EM_center":far_center,"angular_formula":"Hz=(kx/k0 * Ay - ky/k0 * Ax)/Z0; display=(mu0/4)|Hz|^2 / u_EM_angular(0)",
            "color_limits":[[-1,1],[0,.3],[0,.01],[0,.0003]],"notes":"Hz(center)=0 by transverse plane-wave condition; do not normalize by Hz(center)"})

        dq=float(spectrum_q[1]-spectrum_q[0])
        angular_indices=np.flatnonzero(abs(spectrum_q)<=np.sin(np.deg2rad(10.))+dq)
        angular_q=spectrum_q[angular_indices]
        angular_edges=np.rad2deg(np.arcsin(np.r_[angular_q-dq/2,angular_q[-1]+dq/2]))
        assert angular_edges[0]<-10 and angular_edges[-1]>10

        write_json("Fig1_display.json",{"slab_scale":[0.,1.],"near_denominator_J_m3":near_peak*U,
            "near_source":"E_up_surface and H_up_surface; z=100 nm; NA=1",
            "angular_spectrum_denominator_EM_center":far_center,"far_scale":[0.,1.],"far_scale_type":"linear",
            "far_reference":"Full electromagnetic energy angular spectrum at kx=ky=0; not dP/dOmega","far_global_max_over_center":float(np.nanmax(far_display[2])),
            "BZ_q_over_k0":bz.tolist(),"far_axis_limits":[-10.,10.],"far_axis_unit":"degree","far_angle_definition":"theta_x=asin(kx/k0), theta_y=asin(ky/k0)","far_angle_jacobian_applied":False,
            "Ex_row_color_limits":{"slab":[0.,.1],"surface":[0.,.3],"angular_spectrum":[0.,.03]},"Ex_row_saturation":"values above display maxima use top color; stored data unchanged",
            "outside_lightcone":"evanescent upward angular spectrum at z=100 nm; not radiated power","periodic_smoothing":"per-material local RBF plus normalized Gaussian; seven separately polygon-clipped layers"})
        fig,axs=plt.subplots(4,4,figsize=(22,17),layout="constrained")
        for row,(component,lab) in enumerate([(3,"Hz"),(0,"Ex"),(1,"Ey"),(2,"Full EM")]):
            from matplotlib.patches import Polygon
            for vertices,extent,values in display_layers:
                im=axs[row,0].imshow(values[...,component],interpolation="bilinear",interpolation_stage="data",origin="lower",extent=extent,
                    cmap="RdBu_r" if component!=2 else "magma",vmin=-1 if component!=2 else 0,vmax=1)
                im.set_clip_path(Polygon(vertices,closed=True,transform=axs[row,0].transData))
            axs[row,0].set(title=f"Unit cell Γ-point: {(('Im' if hz_part=='imag' else 'Re')+'(Hz)') if component==3 else ('Re('+lab+')' if component!=2 else lab)}",aspect="equal",xlabel=r"$x$ ($\mu$m)",ylabel=r"$y$ ($\mu$m)")
            axs[row,0].plot(cell[:,0],cell[:,1],color="#555555",linewidth=.8)
            for vertices in geometry["holes"]:
                edge=np.asarray(vertices,float)
                edge=np.vstack([edge,edge[0]])
                axs[row,0].plot(edge[:,0],edge[:,1],color="0.45",linewidth=.65)
            axs[row,0].set_xlim(cell[:,0].min()-.02,cell[:,0].max()+.02)
            axs[row,0].set_ylim(cell[:,1].min()-.02,cell[:,1].max()+.02)
            periodic_label=(rf"$\mathrm{{{'Im' if hz_part=='imag' else 'Re'}}}(H_z e^{{-i\phi}})/\max|\mathbf{{H}}|$ (1)" if component==3 else
                (rf"$\mathrm{{Re}}(E_{{{'x' if component==0 else 'y'}}}e^{{-i\phi}})/\max|\mathbf{{E}}|$ (1)" if component!=2 else r"$u/\max u_{\mathrm{EM,cell}}$ (1)"))
            bar=colorbar(fig,axs[row,0],im,periodic_label)
            if component!=2:bar.set_ticks(np.linspace(-1,1,9))
            figure_panel(fig,axs[row,1],fig1_slabmaps[component],sa,f"Field in slab: {lab}, z=0 nm",fig1_labels[component],.3 if component==3 else (.1 if component==0 else slab_max))
            axs[row,1].plot(boundary[:,0],boundary[:,1],**boundary_style)
            figure_panel(fig,axs[row,2],near_display[component],axis,f"Near field intensity: {lab}",fig1_labels[component],.01 if component==3 else (.3 if component==0 else 1.))
            im=axs[row,3].pcolormesh(angular_edges,angular_edges,far_display[component][np.ix_(angular_indices,angular_indices)],
                shading="flat",cmap="magma",vmin=0,vmax=.0003 if component==3 else (.03 if component==0 else 1.),rasterized=True)
            axs[row,3].set(title=f"Far field intensity: {lab}",aspect="equal",
                xlim=(-10,10),ylim=(-10,10),xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10],
                xlabel=r"$\theta_x$ (°)",ylabel=r"$\theta_y$ (°)")
            angular_label=(r"$\mathcal{U}_{h,z}/\mathcal{U}_{\mathrm{EM}}(0)$ (1)" if component==3 else
                (r"$\mathcal{U}_{e,j}/\mathcal{U}_{\mathrm{EM}}(0)$ (1)" if component!=2 else r"$\mathcal{U}_{\mathrm{EM}}/\mathcal{U}_{\mathrm{EM}}(0)$ (1)"))
            colorbar(fig,axs[row,3],im,angular_label)
        fig.suptitle("Mode19 / py | Periodic reference -> finite slab -> radiative surface -> angular spectrum (central angles)\n"
                     "z=0 nm inside slab; z=100 nm in air; NA=1; independent periodic shape scale; white outline: cavity-cladding boundary",fontsize=15)
        save(fig,"Fig.1_analysis_framework",dpi=FIGURE_DPI)
        if fig1_only:return
    figure2()
    fig,axs=plt.subplots(1,3,figsize=(16,5),layout="constrained")
    for ax,value,label in zip(axs,[upue,upum,upue+upum],["Electric energy","Magnetic energy","Full EM energy"]):
        figure_panel(fig,ax,value,axis,label,r"$u/U_{\rm ref}$ (m$^{-3}$)",np.nanmax(upue+upum))
    fig.suptitle("Upward radiative field | direction separated using independent air E/H | z=100 nm, NA=1")
    save(fig,"Fig.3_outgoing_EM_nearfield",dpi=FIGURE_DPI)
    figure4()
    log("Fig.1--Fig.4 updated; former six-group figure retired.")


def convergence_checks():
    dx,q,small,mask,kz,n=axes();h=len(small)//2
    full=load("full_mesh_radiation.npz")["F_V"]
    previous=json.loads((L/"resolution_revision.json").read_text(encoding="utf8"))
    checks=previous["native_directional_validation_before"]["native_directional_checks"]
    # Retain this unchanged quadrature audit at its ORIGINAL physical q.
    # Padding must not relabel the native COMSOL samples using refined-grid indices.
    fine=load("finite_air_1100nm_EH.npz");step=np.diff(fine["x_um"])[0]*1e-6
    e=np.nansum(fine["E_V_m"],axis=(0,1))*step**2
    hh=np.nansum(fine["H_A_m"],axis=(0,1))*step**2
    up=(e[:2]+Z0*np.array([hh[1],-hh[0]]))/2
    center=1j*K0/(2*np.pi)*np.exp(1j*K0*1.1e-6)*up
    direct=load("direct_radiative_fields.npz")["F_V"]
    air_error=float(np.linalg.norm(center-direct[h,h,:2])/np.linalg.norm(direct[h,h,:2]))
    table=[];dq=small[1]-small[0]
    domega=np.where(mask,dq*dq/(K0*np.where(mask,kz,1)),0.)
    for na in [.1,.3,.5,.8,.9,.95,1.]:
        select=mask&((n[...,0]**2+n[...,1]**2)<na*na)
        power=lambda v:float(np.sum(abs(v[select])**2*domega[select,None])/(2*Z0))
        table.append(dict(NA=na,air_power_W=power(direct),mesh_source_power_W=power(full),
                          relative_power_difference=abs(power(full)/power(direct)-1)))
    pd.DataFrame(table).to_csv(table_path(L, "NA_closure.csv"),index=False)
    result={"air_801_to_1601_center_complex_error":air_error,
        "coarse_air_center_F_V_pairs":np.column_stack([center.real,center.imag]).tolist(),
        "fine_air_center_F_V_pairs":np.column_stack([direct[h,h,:2].real,direct[h,h,:2].imag]).tolist(),
        "native_directional_checks":checks,
        "native_directional_audit_grid_fft_size":2001,
        "native_directional_audit_reused_unchanged_operator":True,
        "max_directional_quadrature_error":max(v["relative_complex_error"] for v in checks)}
    write_json("convergence_validation.json",result)
    energy=load("mesh_energy_regions.npz")["energy_components_J"];rows=[]
    for r,name in enumerate(["cavity","cladding","pad"]):
        for col,component in enumerate(["Ex","Ey","Ez","Hx","Hy","Hz"]):
            rows.append(dict(region=name,component=component,energy_J=energy[r,col],method="native tetrahedra, degree 2"))
    for component,value in zip(["Ex","Ey","Ez","Hx","Hy","Hz"],native()["slab_energy_components_J"]):
        rows.append(dict(region="whole slab",component=component,energy_J=value,method="COMSOL native intorder=6"))
    pd.DataFrame(rows).to_csv(table_path(L, "energy_components.csv"),index=False)
    f=load("field_decomposition.npz");cavity=float(energy[0].sum())
    a={"cell_fit_vs_native_mesh_cavity_energy_relative_error":abs(float(f["cavity_energy_J"])/cavity-1),
       "native_regions_energy_fraction":(energy.sum(1)/native()["U_ref_J"]).tolist(),
       "fit_self_fraction_Uref":(np.diag(f["group_gram_J"]).real/native()["U_ref_J"]).tolist()}
    write_json("energy_accuracy.json",a)

    budget=[]
    gram=f["group_gram_J"]
    for i,name in enumerate(NAMES[:4]):
        budget.append(dict(scope="cavity projection quadrature",region=name,kind="self",second=name,energy_J=gram[i,i].real))
        for j in range(i):
            budget.append(dict(scope="cavity projection quadrature",region=NAMES[j],kind="interference",second=name,energy_J=2*gram[j,i].real))
    budget.append(dict(scope="cavity projection quadrature",region="cavity",kind="total",second="all",energy_J=gram.sum().real))
    for name,value in zip(["cavity","cladding","pad"],energy.sum(1)):
        budget.append(dict(scope="native FEM quadrature",region=name,kind="regional total",second=name,energy_J=value))
    budget.extend([
        dict(scope="native FEM quadrature",region="whole slab",kind="total",second="all",energy_J=energy.sum()),
        dict(scope="COMSOL native order6",region="whole slab",kind="independent total",second="all",energy_J=native()["U_ref_J"]),
        dict(scope="integration difference; not a physical source",region="cavity",kind="native minus projection",
             second="quadrature discrepancy",energy_J=energy[0].sum()-gram.sum().real)])
    pd.DataFrame(budget).to_csv(table_path(L, "energy_budget.csv"),index=False)
    log("Convergence audit: "+json.dumps({**result,**a}))


def report():
    """Regenerate the complete per-figure calculation record from current numerical outputs."""
    v=json.loads((L/"validation.json").read_text());f=json.loads((L/"fit_validation.json").read_text())
    a=json.loads((L/"air_validation.json").read_text());mesh=json.loads((L/"mesh_validation.json").read_text())
    conv=json.loads((L/"convergence_validation.json").read_text());energy=json.loads((L/"energy_accuracy.json").read_text())
    g=load("group_radiation.npz")["F_groups_V"];h=g.shape[1]//2
    center_rows=[f"| {name} | {abs(value[0])**2/(2*Z0):.7g} | {abs(value[1])**2/(2*Z0):.7g} |" for name,value in zip(NAMES,g[:,h,h])]
    budget=pd.read_csv(table_path(L, "center_radiation_budget.csv", existing=True))
    y=budget[(budget.channel=="Ey")&(budget.scope=="volume_source")]
    interference=float(y[y.kind=="interference"].value.sum())
    U=v["U_ref_J"]
    slab_norm=json.loads((L/"slab_normalization.json").read_text(encoding="utf8"))
    text=table_links(f"""# mode19：周期 BIC 与有限腔中心辐射的数据、图示和计算报告

此前的数据补全阶段已启动 COMSOL 6.3，从有限腔保存解导出 mode19 的完整场，并完成 18 个非 Γ 周期点的新求解。
连同 Γ 点，共补齐 19 点、每点 4 模式的完整 E/H。有限腔未重新求本征值。
原模型保持不变；有限模型副本、周期 Γ 源副本、19 个周期解全部保存在 ../00_model/，没有删除 MPH。

**六图所需数据已经取得；物理闭合仍未通过。**
本报告严格区分数据可用、代数一致和独立物理验证。
全体源到独立空气面的中心复振幅误差为 **{v['mesh_source_to_air_center_complex_error']:.3%}**，
超过 spec 的 2% 目标；中心功率误差 **{v['mesh_source_to_air_center_power_error']:.3%}**，
全 NA 功率误差 **{v['mesh_source_to_air_total_power_error']:.3%}**，超过原定 5% 目标。
下面的分组结果是指定算子下的量化证据，不能宣称已经精确复现有限模型的全部辐射。

## 1. 当前结论与严格中心预算

1. 实际有限模式的法向辐射由 Ey 主导，Ex/Ey 中心功率比为 {v['center_Ex_over_Ey_power']:.3e}。
   Ex 的数值零部分来自 quarter 解按边界奇偶性镜像，不是新的独立对称性证明。
2. 新 E/H 联合投影在公共正定晶胞求积度量下解释 cavity 储能的
   {1-f['residual_energy_fraction']:.4%}；其余 {f['residual_energy_fraction']:.4%} 为 cavity 残差。
   所有 E/H 使用同一系数，不分别拟合 Ex 与 Ey。
3. around-Γ 离散 py 在严格 q=0 处几乎相消，符合完整离散晶格窗口的结构因子零。
   因而“有 around-Γ 储能”不能自动解释“严格中心亮点”。
4. Γ-py 的 Ey 中心项很小；显著中心项出现在实际 cladding 和 cavity 重建残差中，
   且存在强干涉。此结果支持进一步核查 cladding 与界面引起的非周期场修正。
5. cavity 残差只是实际场减去所选 76 基重建，不是已经独立求解或命名的物理模式。
   体源与空气尚有独立误差，因此不将当前分组提升为唯一且最终的微观机制。

| 体源分组 | 中心 Ex 自身项 (W/sr) | 中心 Ey 自身项 (W/sr) |
| --- | ---: | ---: |
{chr(10).join(center_rows)}

Ey 两两干涉合计 **{interference:.7g} W/sr**；体源相干总中心 Ey 为
**{v['source_group_sum_center_y_W_sr']:.7g} W/sr**；独立空气参考为
**{v['direct_center_y_W_sr']:.7g} W/sr**。
自身项不能直接相加当总功率，也不能被写成互斥贡献百分比。
逐项预算保存在 ../80_logs/center_radiation_budget.csv。

## 2. 身份、保存模型与原始数据

固定 cavity (245 nm,0.96,1.156)、cladding (242 nm,0.98,0.93)、20/20 层、mesh5。
用户 py 对应内部 px；实际 mode19 频率 {FREQ/1e12:.12f} THz。
晶格周期为 0.82 μm；有限模型参数 a=65.873333... μm 是计算包络宽度，不是晶格周期。

薄板实际范围 z∈[-100,100] nm，V_ref 包含 cavity、cladding、pad、孔内空气，不包含 PML。
quarter 解在 x≥0,y≥0,z≥0，按实际 PEC/PMC 和 TE z 镜像补全。
内部场图严格取 z=0 nm；可辐射近场取真实上表面 z=100 nm 空气侧。

| 分类目录 | 内容 |
| --- | --- |
| ../00_model/ | finite_mode19_source.mph、periodic_Gamma_source.mph、periodic_k00–k18.mph |
| ../01_results/ | 原始/派生复场、求积、坐标与掩膜、投影与辐射数组 |
| ../80_logs/ | CSV、JSON、模型/输入哈希、验证和执行日志、图的实际色标清单 |
| ../10_overview/ | PNG 图像预览 |
| ../11_pdf/ | PDF 图像 |
| 当前 12_reports/ | report.md、spec 快照及详细分析报告 |

有限源模型为 {load_dipolar_inputs().finite_model}。
完整绝对路径和源/副本 SHA-256 见 model_manifest.json；最终文件清单见 completion_manifest.json。
新周期解从保存 Γ 模型出发，保持几何和 mesh，每个 k 的解分别保存。
旧 Hz 参考仅用于模式身份、相位对齐和 k 点来源；没有使用旧投影系数或旧辐射响应填图。

关键数据（以下文件名相对 ../01_results/）：

| 数据 | 字段、物理量与用途 |
| --- | --- |
| finite_slab_z0_fine_EH.npz | 3201² 原始 E_V_m/H_A_m，x_um/y_um、epsilon_r、domain_id、valid；实际 slab 场图 |
| finite_surface_fine_EH.npz | z=100.001 nm 空气侧实际 E/H，仅作原始表面核查，不能冒充 NA 回溯场 |
| finite_air_1100nm_fine_EH.npz、finite_air_1400nm_fine_EH.npz | 两个真实空气层高度的 1601² E/H，用于方向分离与传播核对 |
| finite_air_1100nm_EH.npz | 保留原 801² E/H，与新主输入对照横向采样收敛 |
| finite_air_fine_EH.npz | z=1.65 μm 空气/PML 接口加密至 1601²，仅辅助核查 |
| finite_air_high_EH.npz | z=1.95 μm 在 PML 内，保留但排除自由空间传播 |
| finite_cavity_EH_o3.npz、finite_cladding_EH_o3.npz | 公共材料分辨规则求积复 E/H，用于完整能量投影 |
| finite_mesh_geometry.npz、finite_mesh_EH_*.npz | slab 实际四面体和 5,577,972 个正权重 quarter 求积点及 E/H |
| finite_native_cavity_EH.npz、finite_native_cladding_EH.npz、finite_native_pad_EH.npz | 同一实际 mesh 数据按 cell 支持分区；不是空气目标场的拟合结果 |
| periodic_kNN_EH.npz | 19 点四带公共体求积及原 z=0 数据、统一相位、归一化、模式匹配、频率 |
| periodic_kNN_z0_fine_EH.npz | 19 点四带新 81² 背景网格的六角胞内 z=0 E/H，沿用相同原生相位和能量尺度 |
| periodic_kNN_mesh_EH.npz | 每点原生 1035 slab 四面体、4140 节点上的源场和原生中心积分 |
| field_decomposition.npz | 共用 E/H 系数、Gram、分组系数、残差储能、k/band 索引 |
| direct_radiative_fields.npz | 独立空气路径的滤波 E、向上 E/H、完整复远场及坐标 |
| group_radiation.npz、radiative_fields.npz | 六组复远场和回溯 E/H、全体源与空气参考、共同 U_ref |
| group_slab_z0.npz | 六组 z=0 完整 EM 自身密度图，周期胞内插值仅用于显示 |
| radiation_interference.npz | 完整光锥功率的六组复 Gram，用于自身/干涉账目 |

## 3. 公共定义、单位和相位

COMSOL 原生时间约定 exp(+iωt)，向上平面波的空间因子为 exp[-i(q·r+kz z)]。
本分析定义 A(q)=∫E(r)exp(+iq·r)d²r，逆变换为负号。
数值上正号变换由 IFFT 实现，并乘 N²Δx²，不能从函数名猜测物理符号。
E 角谱 A 单位 V·m，H 角谱单位 A·m，远场振幅 F 单位 V。

板内采用无色散、无损 n=3.3 与空气、mu_r=1：
u_e=epsilon0*epsilon_r*(|Ex|²+|Ey|²+|Ez|²)/4；
u_m=mu0*(|Hx|²+|Hy|²+|Hz|²)/4。
全部 slab 域的 COMSOL intorder=6 积分给出共同 **U_ref={U:.12e} J**。
这是指定 slab 包络内的储能，不是无限开放空间全部准正规模能量。

Fig.1 的 Finite slab 列使用 z=0 完整截面电磁能量密度峰值 u_EM,max 作公共分母，色标无量纲。
其余有限近场/分组图仍用同一 U_ref：u/U_ref 单位 m^-3，(dP/dΩ)/U_ref 单位 s^-1 sr^-1。
周期参考另用该周期 z=0 完整 EM 峰值作形状归一化，并明确标注；不能与有限图按颜色比绝对能量。
原生本征模振幅任意，表中 J/W 仅在同一原生尺度下可比。
P/U_ref 是耦合速率量，不是无量纲能量份额，也不能直接作为完整 Q。

Ex/Ey 可各自进行相同线性 Fourier/逆变换并保留相位。
它们在严格法向构成两个正交出射通道；离轴必须保留 Ez 与完整 E/H，不能声称只有两项。

### 3.1 空气传播和可辐射近场

令 nt=q/k0、s=kz/k0、P=(Hy,-Hx)，独立空气 E/H 分离为：

E_t^±=[E_t ± Z0/s*(P-nt*(nt·P))]/2；
Ez^+=-q·E_t^+/kz，Ez^-=+q·E_t^-/kz。

只保留 |q|<k0，NA=1，排除 kz=0 的零测度边界。
从 z_m=1.1 μm 回溯到 z_s=0.1 μm，乘 **exp(+ikz*1.0 μm)**，再逆变换复振幅。
向上 H 由已经验证方向的角谱 H_up=(n×E_up)/Z0 重建。
原始 E 直接滤波与方向分离的纯向上 E/H 分别保存，Fig.2、Fig.3 不混用这两种量。

原始有效孔径外无场数据，FFT 中以零填充并显式承认有限窗口误差；
补零到 {FFT_SIZE}²（上一版为 4001²；本次空气间距减半，Δq 保持约 0.01907 μm^-1）只细化 q 采样，不增加原始空间分辨率或已知物理孔径。
逆变换先在完整补零网格上做，再取 ROI 画图，不把 ROI 重新截成 cavity 源。
近场滤波非局域，亮区不能直接等同局部独立辐射源。
能量密度不是局部穿面功率；后者为 Re(E×H*)_z/2。

### 3.2 远场、完整体源与数值积分

定义远区 E=F(n)exp(-ik0 r)/r：
F=i*kz/(2*pi)*exp(+ikz*z_m)*A_up，
dP/dΩ=|F|²/(2Z0)。
由 F 回溯时先换回平面角谱
A_s=2*pi/(i*kz)*exp(-ikz*z_s)*F，不能直接逆变换 F。

独立源路径采用空气背景 Green 算子：
F=k0²/(4*pi)*(I-n*n^T)*∫(epsilon_r-1)E(r)exp(+ik0*n·r)dV。
所有组在相同 Re(f) 评价；孔内 chi=0，实际 cavity/cladding/pad 均计入。
分组是线性源贡献，不是把截断场视为独立物理解。
此空气背景算子与有限域侧边界是否完全等价由闭合检查决定，不能先验宣布相同。

最终材料源在每个模型自身 mesh 上用对称正权重四点、degree-2 四面体求积；
与 COMSOL 原生 intorder=6 中心及离轴积分独立对照。
x/y 积分用 FINUFFT 2.5.1，正号、精度参数 1e-12。
z 镜像对 Ex/Ey 给出 2cos(kz z)，对 Ez 给出 2i sin(kz z)；
分别保留至 z^8/z^9。这里 |k0 z|≤0.415，已与完整三维直接求和核对。
未采用 Gaussian 核、拟合展宽或从空气目标反拟合传输算子。

全光锥传播功率为 Σ|F|²/(2Z0)*ΔqxΔqy/(k0*kz)，显式保留
d²q=k0² cos(theta)dΩ。像素模平方简单求和不能直接标为 W。

## 4. Fig.1：周期—slab—可辐射表面—远场总览

文件：Fig.1_analysis_framework.png 与../11_pdf/ 中的同名 PDF。四行依次为 Hz、Ex、Ey、Full EM。

Fig.1 四列显示名称统一为 Unit cell Γ-point、Field in slab、Near field intensity、Far field intensity，保留各行分量标识。此次仅修改标题；第四列的 Far field intensity 仍指前述归一化 EM 能量角谱，并未转换为 dP/dΩ。


**Hz 新增首行。** Fig.1 现为四行四列，依次 Hz、Ex、Ey、Full EM。Hz 周期场沿用由 Ey 参考点确定的同一个 phi；该相位下 Re(Hz)/max|H| 峰值约 0.000380，主要场位于虚部，因此首幅明确显示 Im(Hz exp(-i phi))/max|H|（max|H| 为原始单胞三分量复磁场矢量模长最大值），RdBu_r -1–1，不对 Hz 独立旋转相位。保留相同的材料分区插值、显示平滑和几何描边。slab 列为 mu0|Hz|²/4 除以原完整 slab EM 峰值，线性 0–0.3；归一化峰值约 0.318，超过 0.3 的值使用最高颜色，原始数据及归一化不变；surface 列为实际空气 E/H 方向分离并 NA=1 回溯所得 H_up,z 的 mu0|H_up,z|²/4，除以原完整表面 EM 峰值，线性 0–0.01。角谱使用已有复 A_E 计算 A_Hz=((kx/k0)A_Ey−(ky/k0)A_Ex)/Z0，绘制 (mu0|A_Hz|²/4)/U_EM(0)，与 Ex、Ey 和 Full EM 共用完整电磁能量中心角谱分母；严格中心 Hz=0，不以 Hz 自身中心作分母。角度范围仍为 ±10°，线性色标 0–0.0003。原 Ex/Ey/Full EM 的复场、相位与色标保持不变，第四列统一采用完整 EM 中心分母。参数记录在 ../80_logs/Fig1_Hz.json。此次不启动 COMSOL。

**第一列：周期结构。** 数据来自 periodic_k00_z0_fine_EH.npz 的 Γ-py 复 E/H（z=0）。在最大 |Ey| 采样点取 phi=arg(Ey)，对 Ex、Ey 共同乘 exp(-i phi)，使参考点 Ey 为正实数；不分别旋转两分量。前两行分别显示 Re(Ex exp(-i phi))/Emax、Re(Ey exp(-i phi))/Emax，Emax=max sqrt(|Ex|²+|Ey|²+|Ez|²)，在整个有效单胞取最大值。色带匹配参考 s4_edges_py.svg 的 RdBu_r（负蓝、零白、正红），范围 -1–1。这是归一化实部场，不是能量密度，展示符号与空间相位关系；复数模长归一化不保证实部图一定达到 ±1。第三行保留完整 EM 能量密度除以自身截面最大值，magma 0–1。三个面板均描出真实单胞 outer 六边形及 geometry_source 中全部 holes 的闭合边界（六个三角孔，共 18 条孔边）。为消除材料边界的网格锯齿，按真实几何将原始样本分成硅与六个独立空气孔；每组仅用本组复 E/H 数据作局部薄板样条 RBF 插值（最多 32 邻点，smoothing=0）。硅显示网格 641²、每孔 241²，随后计算统一相位的实部与完整电磁能量密度。随后对三个显示标量（两个有符号实部及非负 EM 能量密度）进行同材料掩膜内的归一化高斯平滑：v_smooth=G_sigma*(mask*v)/(G_sigma*mask)，sigma=0.75 sqrt(Δx_raw Δy_raw)，远小于晶格周期且随显示网格换算为像素尺度。硅掩膜排除六个孔，每个孔独立处理；不跨硅–空气界面混色，避免将区域外零值带入边缘。该步骤仅平滑显示标量，不修改复场，相位与原始峰值分母保持不变，也不将平滑后的峰值重新归一到 1。七层图像分别以精确矢量多边形裁剪并叠加，最后描出界面，保留材料界面的跃变。平滑尺度与各区域显示值相对 L2 变化保存在同一验证 JSON；平滑是可视化处理，不应据平滑图读取尖角局部定量峰值。新网格仅改善显示，不代表新增 COMSOL 采样或有限元精度；边界附近是同材料侧的插值/外推，尖角真实场仍受原采样限制。公共相位、Emax 和能量峰值均继续取原始样本，原始复场及辐射计算不变。平滑前插值的源样本重建相对误差要求低于 1e-9，见 ../80_logs/Fig1_periodic_interpolation.json；显示数据见 ../01_results/Fig1_periodic_material_display.npz。归一化分母、相位及参考点见 ../80_logs/Fig1_periodic_real_field.json。

**第二列：finite slab。** 数据来自 finite_slab_z0_fine_EH.npz，z=0 nm，含 cavity、cladding、pad 和孔内空气。计算 u_ex=epsilon0 epsilon_r |Ex|²/4、u_ey=epsilon0 epsilon_r |Ey|²/4，完整量另含 Ez 和全部 H 分量。三行统一除以完整截面的 max(u_EM)，线性色标 0–1。白色描边为实际 cavity–cladding 分界。

**第三列：可辐射表面。** 数据来自 direct_radiative_fields.npz 的 E_up_surface、H_up_surface：先从空气 E/H 角谱分离向上传播场，保留 NA=1 并带传播相位回溯至 z=100 nm，再逆变换。三行都使用同一套向上场。Ex/Ey 行为 epsilon0 |E_up,x/y|²/4，Full EM 为 epsilon0 |E_up|²/4 + mu0 |H_up|²/4；共同除以该表面 max(u_EM,up)，线性色标 0–1。采用峰值归一化定义，但分母取本表面自身的完整 EM 峰值，不与 slab 混用。它展示可耦合自由空间的表面场，不能等同于原始板内储能。

**第四列：完整布里渊区电磁能量角谱（含光锥外数据）。** 从 finite_air_1100nm_fine_EH.npz 的 z=1.1 μm 复数 E/H 出发，采用正号空间 FFT，各方向补零到 16001 点（原 8001 点；保留奇数网格及精确中心），角谱采样密度约提高一倍，二维像素数约为四倍，保留覆盖整个 BZ 的矩形 q 网格。沿用 E/H 方向分离算子得到向上角谱；令 kz=sqrt(k0²-q²)（光锥内）或 kz=-i sqrt(q²-k0²)（光锥外），与 COMSOL exp(+iωt)、向上 exp(-ikz z) 的衰减约定一致。乘 exp(+ikz*(z_input-z_surface))，回溯到 z_surface=100 nm，保留完整复振幅。

四行依次显示磁场 Hz 能量角谱、电场 Ex 能量角谱、电场 Ey 能量角谱及完整 EM 能量角谱。由完整复电场角谱按 A_H=(k/k0)×A_E/Z0 重建三个磁场分量，使用上述同一复 kz 分支。定义 Ue,j=epsilon0|A_E,j|²/4，Uh,j=mu0|A_H,j|²/4，U_EM=sum_j(Ue,j+Uh,j)。四行统一除以 U_EM(kx=ky=0)，完整 EM 中心为 1；所有色标上限沿用已有设置。相比旧电场或磁场各自归一化，Ex、Ey、Hz 数值减半；传播区电磁能量相等，因此完整 EM 归一化图与旧 Full E 图一致，光锥外不假定电磁能量相等。数值检查传播区电磁能量相等，rtol=1e-9。该量是平面场能量的角谱权重，不是局域实空间能量密度，也不是 dP/dΩ；需配合所用傅里叶积分测度理解，未加入功率或角度 Jacobian。派生数据保存于 ../01_results/Fig1_EM_angular_spectrum.npz，含电、磁各分量、完整 EM 及中心分母；原复场数据保持不变。此列用于追踪各电磁分量对中心角谱的贡献。

显示坐标改为 θx=arcsin(kx/k0)、θy=arcsin(ky/k0)，单位为度，横纵范围均为 -10°–10°；它们是方向余弦角，不是极角 θ 与方位角 φ，也不是 atan(kx/kz) 的投影角。仅改变坐标与显示窗口，不乘角度 Jacobian，不将角谱转换为每度功率密度。采用实际角度网格边界的 pcolormesh，避免把非均匀角度网格当作等间距像素。完整 BZ 角谱及光锥外数据仍保存在 NPZ 中，但当前中心视窗完全在光锥内，光锥和 BZ 边界均位于视野外，不画误导性的边界图例。倏逝回溯会放大高 q 误差，有限空气采样孔径也会造成谱泄漏，不能据其亮度直接归因自由空间功率。细化后的完整复数角谱保存在 ../01_results/Fig1_full_BZ_angular_spectrum_2x.npz，旧版文件保留。中心与既有 A_up_surface 核对一致（rtol=1e-10、atol=1e-20），另取两个细网格离轴点，将全部 E/H 分量与直接离散傅里叶求和核对（rtol=1e-8、atol=1e-18）。网格信息见 ../80_logs/Fig1_angular_resolution.json。增加补零细化的是角谱采样，不增加原始 COMSOL 采样精度或物理孔径；没有启动 COMSOL。

**Ex 行显示增强。** Ex 行第 2、3、4 图仍采用线性色标，下限均为 0，上限依次改为 0.1、0.3、0.03；完整场公共归一化分母不变。超出显示上限的数值显示最高颜色，不截断源数据。其余行仍为 0–1，因此调整后颜色亮度不能直接跨 Ex 与 Ey/Full 行比较，定量比较应看 colorbar 数值。

**所得数据与意义。** 四列现在均无量纲，但各列使用自己的物理参考量，颜色不能跨列解释为相同绝对能量。图用于比较周期储能、有限模式储能、可辐射表面场和角分布的形态；并不意味着对能量图片直接做 FFT，也不构成 Γ-BIC 中心亮点来源的因果证明。具体分母、中心与全局峰值关系、BZ 顶点和显示参数见 ../80_logs/Fig1_display.json。此次复用原始结果，未启动 COMSOL；Fig.1 调整时其余图保持不变，后续 Fig.4 重构见第 7 节。

坐标符号统一为 kx、ky，表示有限场傅里叶波数，未折叠到第一布里渊区；本次仅改显示命名，原始数组字段与计算保持不变。

## 5. Fig.2：NA=1 的可辐射近场

文件：[Fig.2_radiative_nearfield.png](../10_overview/Fig.2_radiative_nearfield.png)，../11_pdf/ 中的同名 PDF。

**数据起点。** finite_air_1100nm_fine_EH.npz 的原始 Ex/Ey/Ez，z=1.1 μm，
有效六角孔径之外补零。没有借用旧 NA=0.9 近场图或旧反演缓存。

**操作。** 三复分量分别做正号 Fourier，保留全 NA=1 光锥，
乘 exp[+ikz(z_m-z_s)] 并负号逆变换。回溯终点 z_s=100 nm，绝不延伸到 slab 内 z=0。
保留同一全复矢量，图示 Ex、Ey、全电能密度 epsilon0 Σ|E_rad|²/4，再除 U_ref。
按完整 {FFT_SIZE}² 补零逆变换定义得到 1601² ROI（利用可分离 FFT 逐轴变换后裁切，与完整二维变换后裁切等价）用于显示，不再截断源支持。

**所得数据。** direct_radiative_fields.npz 的 E_raw_surface，单位 V/m；
图像为 u_e/U_ref，单位 m^-3，三图共同用全电能密度峰值作线性色标上限。

**意义与用途。** 表示在当前孔径与传播模型下，能进入自由空间的复场回溯到表面的分布。
它是传播谱滤波，不是 slab 表面全场，不包括 evanescent。
中心亮点取决于可辐射 Ex/Ey 的复面积积分之模平方，不是其强度的面积积分。
亮区由于滤波非局域性不能单独证明辐射源就在该处。

**四面板布局更新。** Fig.2 采用 2×2：左上为原第三幅完整电场能量密度（含 Ex、Ey、Ez，但不含 H），右上为偏振雷达图，左下为 Ey，右下为 Ex。三张场图保持原数据、U_ref 归一化和共同色标。

雷达图使用同一 E_raw_surface 的复 Ex/Ey，在整个已采样近场平面上先计算每个位置的分析器透过强度再面积积分：I(θ)=∫|Ey cosθ−Ex sinθ|² dxdy。0° 指向 +y，逆时针朝 −x，与项目既有偏振图约定一致。计算保留同一点 Ex/Ey 的复相位交叉项，不先对不同位置的复场求和。曲线除以所有分析器角度中的最大积分强度，径向范围 0–1；显示 linear contrast=(Imax−Imin)/(Imax+Imin)，不是完整 Stokes 偏振度，也不是绝对远场透射功率。Ez 不参与平面线偏振分析器投影。该图用于比较整个可辐射近场的 x/y 偏振取向，不能代表每个位置具有相同偏振态。CSV 和 JSON 保存于 ../80_logs/Fig2_polarization_radar.*；已核对 0°/90° 分别恢复 Ey/Ex 积分及任意角度的直接复场投影。

## 6. Fig.3：方向分离后的完整电磁可辐射近场

文件：[Fig.3_outgoing_EM_nearfield.png](../10_overview/Fig.3_outgoing_EM_nearfield.png)，../11_pdf/ 中的同名 PDF。

**数据起点。** 与 Fig.2 同高、同孔径独立导出的 E 和 H。
先用第 3.1 节关系区分向上/向下，随后按向上平面波关系计算 H_up；
不是只凭 E 假定全部向上，也不是把推导空气 H 当作 slab H。

**操作。** 向上 E/H 角谱使用同一 NA 掩膜、同一相位回溯到 z=100 nm。
三图依次 epsilon0 Σ|E_up|²/(4U_ref)、mu0 Σ|H_up|²/(4U_ref)、二者和。
三图共享完整 EM 的线性色标上限；x/y 等比例。

**所得数据。** direct_radiative_fields.npz 的 E_up_surface 和 H_up_surface，
单位 V/m 与 A/m；图像单位 m^-3。
Fig.2 原始分量滤波场与本图的方向分离场分别保留，不用名称掩盖算子差异。

**意义与用途。** 在表面比较电能、磁能和完整传播场能量。
叠加场局部 u_e 与 u_m 不必处处相等，不能用局部不等否定各平面波阻抗关系。
z=1.1/1.4 μm 的向下/向上功率比分别
{a['finite_air_1100nm_fine_EH']['down_up_power_ratio']:.6%}、
{a['finite_air_1400nm_fine_EH']['down_up_power_ratio']:.6%}。
有限窗口亦可产生假向下项，不将全部小差异解释为真实反射。

## 7. Fig.4：cavity / cladding 两区域贡献链

文件：[Fig.4_cavity_cladding_chain.png](../10_overview/Fig.4_cavity_cladding_chain.png)，../11_pdf/ 中的同名 PDF。旧 Fig.4_farfield_channels 和 Fig.6_component_chain 的图文件由本图替代，科学数据与 MPH 不删除。

**完整场新增行。** 当前为三行四列：第一行 Full field，第二行 Cavity only，第三行 Cladding only。完整 slab 行取原始整个有效截面的实际 EM 能量（含 cavity、cladding、pad 与孔内空气）；完整辐射源按原六组复 F 相干求和，包含 pad 及所有组间干涉，再按同一传播算子计算 Ex/Ey 表面场与角谱。不是将下面两行的强度相加。完整行使用体源链路，与区域行一致，不混入 Fig.1 独立空气链路；因此仍受源—空气闭合未通过的限制。原 cavity/cladding 两行的复场、分母和色标均保持不变，完整行沿用相同分母（slab 完整峰值、原两个区域的 near EM 公共峰值、原两个区域中心角谱的较大值），不另行归一化。

**完整场近场显示修正。** 第一排第二、三图原色标 0–0.3、0–1 下，归一化峰值仅为 0.00630507、0.0173900（使用约 2.10%、1.74% 的颜色范围），导致视觉近黑。数据非零，未发现空数组或传播结果丢失。仅将两图色标分别调整为线性 0–0.01、0–0.03，保留三排公共归一化分母、原复场和积分值。调整后峰值占显示上限约 63.1%、58.0%，不发生峰值饱和。区域行继续使用 0–0.3、0–1；此后必须结合 colorbar 数字比较，不能直接按跨行颜色明暗比较强度。

**区域行与问题。** 第二排 Cavity only，第三排 Cladding only；从左到右为 slab EM、Ex radiative near field、Ey radiative near field、Far field intensity: Full E（用户所指 far-field 列，物理量严格沿当前 Fig.1）。分析的是同一个 mode19 的分区源贡献，不是重新求解删除另一种材料后的不同结构。outer pad 两排均排除。近场由非局域传播产生，不应再次按源区域裁剪：cavity 源的近场可以延伸到 cladding 位置，反之亦然。

**第一列输入与处理。** 读取 finite_slab_z0_fine_EH.npz 的 3201² 实际复 E/H 与 epsilon_r，在 z=0 计算 u_EM=(epsilon0 epsilon_r |E|²+mu0 |H|²)/4。根据实际 config 的 bulk_points、cladding_points 和单胞六边形逐点分类；与原 group_slab_z0.npz 在重合点核对一致。第二排只保留 cavity 的 u_EM，第三排只保留 cladding；其余有效位置置零，结构外保持空白。这里 cavity/cladding 均包含各自晶胞的孔内空气，孔内极化源则因 epsilon_r−1=0 自然不参与辐射积分。两排共同除以原始完整 z=0 截面的 max(u_EM)，线性 0–1，保留实际 cavity–cladding 白色边界。

**区域复数源。** 读取 group_radiation.npz 的六组 F_groups_V。F_cavity=F_Gamma-py+F_around-Gamma-py+F_other+F_residual；必须先加复振幅，不能相加原 Fig.6 四行的强度。残差定义保证这四项相加恢复实际 cavity 原生体源积分。F_cladding=原 Actual cladding 项。两区域之和加回 pad 应重建 full_mesh_radiation.npz，验证相对误差 <1e-10。该核对只是源分组代数一致性，不是独立 Maxwell 闭合证明。

源积分使用完整薄板厚度内的介质极化源，F=k0²/(4π)(I−nnᵀ)∫(epsilon_r−1)E exp(+ik0 n·r)dV。不是对第一行二维能量密度做 FFT。源场的空间裁剪发生在实际区域体源上，而非传播后裁剪近场强度。

**第二、三列：表面可辐射场。** 对每个区域取光锥 NA=1 内角谱 A_s=2π/(i kz) F exp(−i kz z_s)，z_s=100 nm，再通过与 Fig.1 相同的正号分析/负号逆变换约定回溯到表面。H_s=n×A_s/Z0，同样逆变换。得到两个区域各自的完整复 E/H，分别计算 u_ex=epsilon0 |Ex|²/4、u_ey=epsilon0 |Ey|²/4。共同分母取“两区域所有表面位置中完整 u_EM 的最大值”，不是对每排单独归一化；因此两排能直接比较。该定义沿用 Fig.1 的完整 EM 峰值归一化原则，但参考集合改成这两组区域源。Ex 列线性 0–0.3，Ey 列线性 0–1，与当前 Fig.1 一致。与 Fig.1 比绝对幅值时须查各自分母，不能只看颜色。

**第四列：中心角谱。** 计算 |A_s,x|²+|A_s,y|²+|A_s,z|²，取 cavity、cladding 两排的中心强度中较大者为公共分母；因此至少一排中心为 1，不抹掉两区域中心幅值差异。线性色标 0–1，超上限时显示顶色但原数据保留。θx=arcsin(kx/k0)、θy=arcsin(ky/k0)，两轴 ±10°，采用非均匀角度网格 pcolormesh。当前窗口完全位于光锥内，不需新增光锥外源积分。它是电场角谱强度，未乘 kz²/(8π² Z0)，因此不是 dP/dΩ；与当前 Fig.1 相同。所有列色标无量纲。

**结果文件。** ../01_results/Fig4_region_chain.npz 保存区域掩膜、slab 能量图、区域复远场/角谱、表面复 E、完整 EM 密度及归一化分母。../80_logs/Fig4_region_chain.json 保存色标、峰值、饱和比例、区域重建检查与独立闭合状态；figure_manifest.json 保存实际图面设置。

**物理意义与限制。** 本图比较 cavity 与 cladding 的储能位置和各自源所能产生的表面/角分布，显示区域自项。两排强度相加不等于完整模式强度；必须恢复区域之间的复数干涉，并计入 pad 才能还原完整体源结果。不能将两排看成独立求解的腔或把任一列的自项功率解释为独立损耗通道。原有体源与独立空气远场的闭合尚未通过（中心功率约 9.57%、全 NA 功率约 73.40% 差异），新图沿用并显式标注该限制，不能据它做已验证的定量因果结论。此次未启动 COMSOL，未改动或删除模型。

## 10. 验证结果、有效范围与后续研究

| 核对 | 结果 | 结论/用途 |
| --- | ---: | --- |
| 原生 mesh 储能 vs COMSOL intorder=6 | {mesh['energy_relative_error']:.6%} | 完整体求积与能量权重核对 |
| 原生 mesh 中心体源 vs COMSOL intorder=6 | {mesh['center_source_relative_error']:.6%} | 强相消中心积分核对 |
| 多方向原生积分最大复场差（既有固定物理 q） | {conv['max_directional_quadrature_error']:.6%} | 离轴源求积独立核对 |
| 两空气高度中心复场差 | {a['height_center_relative_error']:.6%} | 回溯相位/空气层传播核对 |
| 两高度完整 F 相对范数差 | {a['height_full_F_relative_error']:.6%} | 孔径与完整角度传播限制 |
| 801²→1601² 同高空气中心复场差 | {conv['air_801_to_1601_center_complex_error']:.6%} | 横向采样，不代表孔径收敛 |
| 投影规则 cavity 储能 vs 原生 mesh 分区 | {energy['cell_fit_vs_native_mesh_cavity_energy_relative_error']:.6%} | 投影能量比例仍有此求积度量误差 |
| 六来源 vs 完整体源中心复场差 | {v['cell_groups_to_mesh_center_complex_error']:.3e} | 代数/源分区检查，不是独立物理闭合 |
| 完整体源 vs 空气中心复场差 | {v['mesh_source_to_air_center_complex_error']:.6%} | **不通过 ≤2%** |
| 完整体源 vs 空气中心功率差 | {v['mesh_source_to_air_center_power_error']:.6%} | **不通过 ≤5%** |
| 完整体源 vs 空气全 NA 功率差 | {v['mesh_source_to_air_total_power_error']:.6%} | **不通过 ≤5%** |

实际 native mesh 区域能量约占 U_ref：cavity {energy['native_regions_energy_fraction'][0]:.5%}，
cladding {energy['native_regions_energy_fraction'][1]:.5%}，pad {energy['native_regions_energy_fraction'][2]:.5%}。
区域和各六分量数值见 energy_components.csv；与原生全域积分的微小差异不通过缩放抹去。

不同 NA 的空气/体源功率对照见 NA_closure.csv。
全 NA 特别受掠射角、有限孔径及侧边界影响；精确成因尚未分离，不能将其全部归于某一个误差源。
实际有限模型的 physics 边界特征已经直接读出，保存在 finite_boundary_conditions.json。
本次从保存解加密空间场采样并同步调整 Fourier 补零；体求积未改动。多方向 COMSOL 原生求积核对沿用此前固定物理 q 的验证，
没有把旧网格编号当作新网格的相同物理坐标。
镜像导致的 Ex 中心零与独立 E/H 数据检验在解释中分别标明。

明确排除但保留的数据：

- finite_air_high_EH.npz 在 PML 中，不用于自由空间传播。
- source_F_group0–5.npz 是公共规则求积的已拒绝试算；
  强相消辐射严重失准，最终 group_radiation.npz 使用原生 FEM 响应。
- 旧 S5 辐射、旧拟合系数与旧中心预算未用作本轮结论。
- 公共正定求积保留用于 E/H 投影，其储能精度与原生源积分分别记录，不伪称已完全收敛。

当前固定 mode19、mesh5 的原始场、19 个周期解以及六图所需数组已经补齐。
要进一步宣称机制完全闭合，还需要研究性计算：

1. 有限模型 mesh 收敛，核查独立空气/材料源的复相位差。
2. 增大物理空气孔径或使用完整包围面，分离高角辐射与有限窗口损失。
3. 核查侧边界、介质外缘与空气 Green 背景是否等价，必要时使用匹配实际边界的响应算子。
4. 提高联合 E/H 投影求积精度并扩展 k/带字典，进一步解释 cavity 残差。

上述未验证项没有被标成完成；没有删去残差、重新对齐目标远场相位、事后改归一化或放宽验收阈值。

## 11. 复现入口

COMSOL 导出/补算入口 scripts/run_main/run_dipolar_data_completion.py；
分析入口 scripts/analysis/analyze_dipolar_complete.py。
COMSOL 6.3 本机 4 核；FINUFFT 依赖记录在 pyproject.toml/uv.lock。
本次原始取点修订已启动 COMSOL，重新评价保存解的更密 E/H 空间点，不重新求本征值。
六图维持 {FIGURE_DPI} dpi（更早版本 190 dpi），PDF 同步更新。原始平面网格各方向间距减半；
mesh5 和体求积不变，不将此称为有限元网格加密。所有旧原始场和 MPH 保留。
实际阶段、版本、模型和来源哈希在 80_logs。
仅重绘使用分析入口，不重新启动本征求解。

建议先看 Fig.1 的实际场总览，再看新 Fig.4 与中心干涉预算。
每幅图均已列明数据起点、操作、输出量、单位、物理意义和解释限制。
""")
    revision_path=L/"resolution_revision.json"
    if revision_path.exists():
        revision=json.loads(revision_path.read_text(encoding="utf8"))
        if "images_after" in revision:
            text+="\n## 12. 上一轮图像和补零修订记录（历史）\n\n"
            text+="六图从 190 dpi 提高至 380 dpi；PNG 像素宽高约加倍，PDF 嵌入图同步更新。\n"
            text+="\n| 图 | 原 PNG 像素 | 新 PNG 像素 |\n| --- | --- | --- |\n"
            for name,item in revision["images_after"].items():
                before=revision["images_before"][name]["pixels"];after=item["pixels"]
                text+=f"| {name} | {before[0]}×{before[1]} | {after[0]}×{after[1]} |\n"
            text+="\nFourier 补零网格 2001²→4001²；Δq 从 0.0381340671 降为 0.0190717991 μm^-1。\n"
            text+="全 NA 光锥内采样由 37,177 增至 148,681；直接场、六来源及近场均按相同网格重新计算。\n"
            text+="原始物理采样未改变；补零没有创造新的物理细节或独立 k 模式。中心复振幅保持数值一致，\n"
            text+="精确误差及 PNG 尺寸验证见 ../80_logs/resolution_revision.json。没有启动 COMSOL，也没有修改 MPH。\n"
    sampling_path=L/"sampling_export_validation.json"
    if sampling_path.exists():
        sampling=json.loads(sampling_path.read_text(encoding="utf8"))
        text+="\n## 13. 本轮实际取点加密及图像数据更新\n\n"
        text+="有限和周期场均从保留的 COMSOL 解直接评价新点；旧 NPZ 保留。1.1 μm 空气面复用此前真实加密导出，其余新导出。\n"
        text+="\n| 截面 z (μm) | 原网格 | 新网格 | 新间距 (nm) | 重合点 E 相对差 | 重合点 H 相对差 |\n| --- | --- | --- | ---: | ---: | ---: |\n"
        for row in sampling["planes"]:
            text+=f"| {row['z_um']} | {row['old_shape'][0]}² | {row['new_shape'][0]}² | {row['new_step_um']*1000:.7f} | {row['E_V_m_nested_relative_error']:.3e} | {row['H_A_m_nested_relative_error']:.3e} |\n"
        text+="\n19 点每点四带的周期 z=0 截面，由 41² 背景网格改为 81²；仅保留真实六角胞内点。沿用原模式配对、相位与完整能量归一化。\n"
        text+="Fig.1 周期列读取新周期截面；Fig.1 有限板内列使用 3201² 原始截面。Fig.2/3/4 从 1601² 空气 E/H 重做 FFT、方向分离、NA=1 与带相位回溯。\n"
        text+="旧 Fig.6（已替代）板内列在新 3201² 数据中隔点读取 1601²（此前 801²），周期显示重建使用新单胞点；近场列也是 1601²。\n"
        text+="空气间距减半后补零 4001²→8001²，保持 Δq≈0.01907 μm^-1；回溯 ROI 1601²。补零不是新增独立 K 模式。\n"
        text+="体源仍用已审计的原生四面体求积，投影体求积与系数不变；不是有限元网格收敛试验。所有独立物理闭合误差仍按原 2%/5% 门限评价。\n"
        text+="实际导出及重合点检查见 ../80_logs/sampling_export_validation.json；加密前后功率、角谱和完整来源清单见 ../80_logs/sampling_revision.json。\n"
        revision=json.loads((L/"sampling_revision.json").read_text(encoding="utf8"))
        if "same_q_sampling_check" in revision:
            check=revision["same_q_sampling_check"]
            text+="\n为避免混淆补零与取点变化，将旧 801² 空气场用正号 NUFFT 直接评价于新 FFT 的完全相同物理 q，再进行相同 E/H 方向分离、回溯相位和 dΩ 功率积分。\n"
            text+=f"新旧采样的中心复振幅相对差为 **{check['center_complex_relative_change']:.6%}**；同 q 全 NA 功率由 **{check['coarse_power_W_same_q']:.9g} W** 变为 **{check['fine_power_W_same_q']:.9g} W**，相对差 **{check['power_relative_change_same_q']:.6%}**。\n"
            text+=f"光锥加权复场相对差 {check['lightcone_weighted_complex_relative_change']:.6%}；新 FFT 中心与原始场面积直接求和的相对差 {check['fine_FFT_center_vs_raw_sum_relative_error']:.3e}。\n"
        text+="补零计算采用可分离变换：每个坐标方向均补到 8001 点，完成该方向变换后裁切需要的输出；与完整 8001² 变换后裁切等价，并已用完整二维算子逐点验证。\n"
        text+="取样过程中发现对称轴与空气/PML 接口的有限元迹值对微小坐标舍入敏感，已在这些位置使用精确请求坐标重新原生取值；没有复制旧场值或放宽重合点容差。接口面仅辅助核查，主辐射输入始终为真实空气内部 z=1.1 μm。\n"
    spatial_report=R/"spatial_filter_report.md"
    if spatial_report.exists():text+="\n"+spatial_report.read_text(encoding="utf8")
    unitcell_report=R/"unitcell_synthesis_report.md"
    if unitcell_report.exists():text+="\n"+unitcell_report.read_text(encoding="utf8")
    cell_report=R/"unitcell_center_radiation_report.md"
    if cell_report.exists():text+="\n"+cell_report.read_text(encoding="utf8")
    sector_report=R/"sector_interference_report.md"
    if sector_report.exists():text+="\n"+sector_report.read_text(encoding="utf8")
    (R/"report.md").write_text(text,encoding="utf8")
    write_json("framework_status.json",{"COMSOL_started":True,"data_complete":True,
        "all_current_figures_populated":True,"current_figure_count":len(list(O.glob("Fig.*.png"))),"independent_physical_closure_passed":v["independent_physical_closure_passed"],
        "finite_slab_z_nm":0,"surface_z_nm":100,"report":"12_reports/report.md"})
    log("Detailed per-figure report updated; independent closure failure remains explicit.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("stage",choices=["fit","air","mesh","source","check","finish","report","refine"])
    args=p.parse_args()
    if args.stage=="refine":
        self_check();air();sampling_checks();mesh_audit();source_fields()
        finalize_arrays();convergence_checks();slab_groups();figures();report()
    elif args.stage=="finish":
        finalize_arrays();convergence_checks();figures();report()
    else:
        {"fit":fit,"air":air,"mesh":mesh_audit,"source":source_fields,"check":self_check,"report":report}[args.stage]()
