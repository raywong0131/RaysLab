"""S5 stage A and shared diagnostics. Offline; never starts COMSOL on import."""
from __future__ import annotations

from comsol_workflow.output_paths import table_links, table_path
from comsol_workflow.dipolar_radiation import plane_farfield
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.constants import c, epsilon_0, mu_0
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from comsol_workflow.hex_lattice_utils import unit_cell_corners
from comsol_workflow.figure_output import save_figure_formats
from scripts.analysis.analyze_dipolar_singularity import WORK, CONFIG, digest

Z0 = np.sqrt(mu_0 / epsilon_0)
DATA = WORK / '01_results'

def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))

def put(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf8')

def load(name):
    with np.load(DATA / name) as z:
        return {k:z[k] for k in z.files}

def save(fig, out, name):
    plt.rcParams.update({'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    for ax in fig.axes: ax.tick_params(direction='in')
    fig.canvas.draw()
    fig.set_layout_engine('none')
    for ax in fig.axes:
        info=getattr(ax,'_colorbar_info',{})
        parents=info.get('parents',[])
        if parents:
            positions=[p.get_position() for p in parents];bounds=ax.get_position()
            lo=min(p.y0 for p in positions);hi=max(p.y1 for p in positions)
            ax.set_axes_locator(None);ax.set_box_aspect(None);ax.set_position([bounds.x0,lo,bounds.width,hi-lo])
            ax.yaxis.label.set_rotation(270);ax.yaxis.labelpad=17
    save_figure_formats(fig, out, name, formats=('png',), dpi=240, bbox_inches='tight')
    save_figure_formats(fig, out, name, formats=('pdf',), bbox_inches='tight')
    plt.close(fig)

def metric(periodic):
    w=2*periodic['weights_half_m3']/4
    return np.column_stack([np.repeat((epsilon_0*periodic['epsilon_r']*w)[:,None],3,1),np.repeat((mu_0*w)[:,None],3,1)])

def project(fields, basis, weights):
    """One fixed complex reference; no independent per-cell phase rotations."""
    denominator=float(np.sum(abs(basis)**2*weights))
    coeff=np.sum(fields*basis.conj()*weights,axis=(-2,-1))/denominator
    energy=np.sum(abs(fields)**2*weights,axis=(-2,-1))
    residual=np.sum(abs(fields-coeff[...,None,None]*basis)**2*weights,axis=(-2,-1))
    return coeff, energy, residual

def local_projection(out):
    p=load('periodic_k00_EH.npz'); f=load('finite_cavity_EH_o3.npz')
    np.testing.assert_allclose(p['rho_xyz_um'],f['rho_xyz_um'],atol=1e-13)
    b=np.concatenate([p['E_normalized'][1],p['H_normalized'][1]],-1)
    v=np.concatenate([f['E_V_m'],f['H_A_m']],-1); m=metric(p)
    coeff,energy,residual=project(v,b,m)
    global_coeff=coeff.mean(); global_residual=float(np.sum(abs(v-global_coeff*b)**2*m)/energy.sum())
    cfg=read(CONFIG/'config.json'); shells=np.array([x['shell'] for x in cfg['bulk_points']])
    np.savez_compressed(out/'01_results/local_gamma_projection.npz',coefficient=coeff,energy_J=energy,residual_J=residual,centers_um=f['centers_um'],shells=shells,global_coefficient=global_coeff)
    pd.DataFrame(dict(cell_id=np.arange(len(coeff)),shell=shells,x_um=f['centers_um'][:,0],y_um=f['centers_um'][:,1],gamma_re_sqrtJ=coeff.real,gamma_im_sqrtJ=coeff.imag,energy_J=energy,residual_J=residual,residual_fraction=residual/energy)).to_csv(table_path(out / '80_logs', 'local_gamma_projection.csv'),index=False)
    rows=[]
    for s in np.unique(shells):
        mask=shells==s; aa=coeff[mask]
        rows.append(dict(shell=int(s),cells=int(mask.sum()),mean_gamma_re=float(aa.mean().real),mean_gamma_im=float(aa.mean().imag),phase_coherence=float(abs(aa.sum())/abs(aa).sum()),local_gamma_explained_fraction=float(1-residual[mask].sum()/energy[mask].sum())))
    pd.DataFrame(rows).to_csv(table_path(out / '80_logs', 'gamma_shell_projection.csv'),index=False)
    corners=f['centers_um'][:,None,:]+unit_cell_corners(cfg['a'])[None]
    fig,axs=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    values=[abs(coeff),np.angle(coeff*global_coeff.conjugate()),residual/energy]
    for ax,vals,title,cmap,lims,label in zip(axs,values,['Local Gamma-py amplitude','Phase relative to global Gamma','Local residual energy'],['magma','twilight_shifted','magma'],[(0,float(abs(coeff).max())),(-np.pi,np.pi),(0,1)],[r'$|a_j|$ (J$^{1/2}$)',r'Phase (rad)',r'$U_{res,j}/U_j$ (1)']):
        pc=PolyCollection(corners,array=vals,cmap=cmap,clim=lims,edgecolors='none');ax.add_collection(pc);ax.autoscale_view();ax.set(aspect='equal',xlabel='x (um)',ylabel='y (um)',title=title);fig.colorbar(pc,ax=ax,label=label,shrink=.8)
    save(fig,out,'GS1a_local_Gamma_projection')
    # Fixed z sample nearest slab center, not an independently rephased field.
    z=f['rho_xyz_um'][:,2];sel=np.isclose(z,z.min());xy=(f['centers_um'][:,None,:]+f['rho_xyz_um'][None,:,:2])[:,sel].reshape(-1,2)
    actual=v[:,:,:3][:,sel]; reconstruction=(coeff[:,None,None]*b[:,:3])[:,sel]
    cl=load('finite_cladding_EH_o3.npz'); ce=cl['E_V_m'][:,sel]; cxy=(cl['centers_um'][:,None,:]+cl['rho_xyz_um'][None,:,:2])[:,sel].reshape(-1,2)
    scale=max(float(abs(actual[...,:2]).max()),float(abs(ce[...,:2]).max()))
    fig,axs=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for pol in range(2):
        for j,(pos,field,title) in enumerate([(xy,actual,'Actual cavity'),(xy,reconstruction,'Local Gamma reconstruction'),(cxy,ce,'Actual cladding')]):
            artist=axs[pol,j].scatter(*pos.T,c=field[...,pol].real.ravel(),s=.1,linewidths=0,edgecolors='none',cmap='RdBu_r',vmin=-scale,vmax=scale,rasterized=True)
            axs[pol,j].set(aspect='equal',xlabel='x (um)',ylabel='y (um)',title=title+' | Re(E'+ 'xy'[pol]+')')
    fig.colorbar(artist,ax=axs.ravel().tolist(),label='Re(E) (V/m)',shrink=.8)
    fig.suptitle(f'Common phase and linear scale | z = {z.min()*1000:.3f} nm')
    save(fig,out,'GS1b_actual_and_Gamma_fields')
    return dict(local_Gamma_explained_fraction=float(1-residual.sum()/energy.sum()),global_single_Gamma_explained_fraction=1-global_residual,definition='energy inner product of full 3D E/H in cavity; scalar py projection per cell; not joint 76-basis coefficients',local_projection_not_pure_global_Bloch_mode=True)

def center_budget(out):
    d=load('unitcell_center_radiation.npz'); group=np.where(d['regions']=='cavity',0,np.where(d['regions']=='pad',21,d['shells']-19))
    f=np.array([d['F_center_V'][group==j].sum(0) for j in range(22)])
    names=['cavity']+[f'cladding_{j}' for j in range(1,21)]+['pad'];total=f.sum(0); cumulative=f.cumsum(0)
    rows=[]
    for j,name in enumerate(names):
        for pol,label in enumerate(['Ex','Ey']):
            rows.append(dict(region=name,channel=label,Re_F_V=f[j,pol].real,Im_F_V=f[j,pol].imag,self_W_sr=abs(f[j,pol])**2/(2*Z0),signed_W_sr=np.real(f[j,pol]*total[pol].conjugate())/(2*Z0),cumulative_Re_V=cumulative[j,pol].real,cumulative_Im_V=cumulative[j,pol].imag,cumulative_W_sr=abs(cumulative[j,pol])**2/(2*Z0)))
    pd.DataFrame(rows).to_csv(table_path(out / '80_logs', 'center_amplitudes.csv'),index=False)

    interface=[]
    for width in [1,2]:
        mask=(d['shells']>=20-width)&(d['shells']<20+width)
        fi=d['F_center_V'][mask].sum(0)
        for pol,ch in enumerate(['Ex','Ey']):
            interface.append(dict(width_each_side_cells=width,channel=ch,Re_F_V=fi[pol].real,Im_F_V=fi[pol].imag,self_W_sr=abs(fi[pol])**2/(2*Z0),signed_W_sr=np.real(fi[pol]*total[pol].conjugate())/(2*Z0),included_shells=','.join(map(str,range(20-width,20+width))),add_to_main_budget=False))
    pd.DataFrame(interface).to_csv(table_path(out / '80_logs', 'interface_band_budget.csv'),index=False)

    np.savez_compressed(out/'01_results/center_amplitudes.npz',F_V=f,cumulative_F_V=cumulative,names=names)
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for pol,color in [(0,'#2374ab'),(1,'#ba3030')]:
        chain=np.r_[0j,cumulative[:,pol]];axs[0].plot(chain.real,chain.imag,'o-',ms=3,color=color,label=['Ex','Ey'][pol]);axs[1].plot(np.arange(22),abs(cumulative[:,pol])**2/(2*Z0),color=color,label=['Ex','Ey'][pol]);axs[2].plot(np.arange(22),np.real(f[:,pol]*total[pol].conjugate())/(2*Z0),color=color,label=['Ex','Ey'][pol])
    axs[0].set(xlabel='Re(F) (V)',ylabel='Im(F) (V)',title='Cumulative normal phasor',aspect='equal');axs[0].legend()
    for ax,title,yl in [(axs[1],'Coherent cumulative center','dP/dOmega (W/sr)'),(axs[2],'Signed contribution','Contribution (W/sr)')]:
        ax.set(xlabel='Added cladding rings (21 = pad)',ylabel=yl,title=title);ax.axhline(0,color='.6',lw=.5)
    save(fig,out,'GS2a_coherent_center_buildup')
    return f

def radiation(xyz, weighted_e, q, k0):
    import finufft
    qx,qy=np.meshgrid(q,q);kz=np.sqrt(k0*k0-qx*qx-qy*qy);n=np.stack([qx,qy,kz],-1)/k0
    plan=finufft.Plan(1,(len(q),len(q)),eps=1e-11,isign=1,nthreads=4)
    plan.setpts(np.ascontiguousarray(xyz[:,0]*1e-6*(q[1]-q[0])),np.ascontiguousarray(xyz[:,1]*1e-6*(q[1]-q[0])))
    out=np.zeros((*qx.shape,3),complex);t=k0*xyz[:,2]*1e-6
    for pol in range(3):
        for m in range(5):
            power=2*m+(pol==2)
            v=plan.execute(np.ascontiguousarray(weighted_e[:,pol]*t**power)).T
            sx=(-1,1,1)[pol];sy=(-1,1,-1)[pol]
            v=v+sx*v[:,::-1]+sy*v[::-1]+sx*sy*v[::-1,::-1]
            out[...,pol]+=2*(1j if pol==2 else 1)*(-1)**m/math.factorial(power)*v*(kz/k0)**power
    return (out-n*np.sum(out*n,-1)[...,None])*k0*k0/(4*np.pi)

def angular_shells(out,expected):
    """Conservative geometric shell clipping; angular source affine at FEM nodes."""
    cache=out/'01_results/shell_farfields.npz'
    if cache.exists():
        with np.load(cache) as d:return {k:d[k] for k in d.files}
    from scripts.analysis.audit_dipolar_cell_quadrature import NORMALS,BARY_INV,clipped_tetrahedron_moment
    cfg=read(CONFIG/'config.json');points=cfg['bulk_points']+cfg['cladding_points'];centers=np.array([[p['x'],p['y']] for p in points]);shell=np.array([p['shell'] for p in points]);groups=np.where(shell<20,0,shell-19)
    selected=np.flatnonzero((centers[:,0]>=-1e-10)&(centers[:,1]>=-1e-10));tree=cKDTree(centers[selected]);geo=load('finite_mesh_geometry.npz');chosen=np.flatnonzero(np.isin(geo['domain_ids'],geo['slab_domains']));records=read(WORK/'80_logs/finite_mesh_exports.json')['parts']
    nodes=[[] for _ in range(22)];sources=[[] for _ in range(22)];closure=0.
    for part,record in enumerate(records):
        d=load(Path(record['path']).name);ids=d['tetrahedron_ids'];np.testing.assert_array_equal(ids,chosen[part*25000:part*25000+record['tetrahedra']])
        material=d['epsilon_r'].reshape(-1,4)[:,0]>1
        v=geo['vertices_um'][geo['tetrahedra'][ids[material]]];xyz=d['xyz_um'].reshape(-1,4,3)[material];w=d['weights_octant_m3'].reshape(-1,4)[material];e=d['E_V_m'].reshape(-1,4,3)[material]*(d['epsilon_r'].reshape(-1,4)[material]-1)[...,None]
        middle=v.mean(1);nearest=selected[tree.query(middle[:,:2])[1]];dist=np.einsum('tvc,nc->tvn',v,NORMALS)-(cfg['a']/2+centers[nearest]@NORMALS[:,:2].T)[:,None,:];whole=np.all(dist<=1e-12,axis=(1,2));gw=np.zeros((22,len(v),4))
        for g in range(21):
            mask=whole&(groups[nearest]==g);gw[g,mask]=w[mask]
        radius=np.linalg.norm(v[:,:,:2]-middle[:,None,:2],axis=-1).max(1)
        for ti in np.flatnonzero(~whole):
            for local in tree.query_ball_point(middle[ti,:2],radius[ti]+cfg['a']/np.sqrt(3)+1e-10):
                ci=selected[local];offs=cfg['a']/2+NORMALS[:,:2]@centers[ci]
                if np.any((v[ti]@NORMALS.T-offs).min(0)>1e-12):continue
                moments=clipped_tetrahedron_moment(v[ti],NORMALS,offs)
                if moments.sum()>1e-22:gw[groups[ci],ti]+=moments@BARY_INV*1e-18
        gw[21]=w-gw[:21].sum(0)
        closure=max(closure,float(np.max(abs(gw.sum(0)-w))))
        for g in range(22):
            use=np.max(abs(gw[g]),axis=1)>1e-34
            if np.any(use):nodes[g].append(xyz[use].reshape(-1,3));sources[g].append((e[use]*gw[g,use,:,None]).reshape(-1,3))
        print(f'shell clipping {part+1}/{len(records)}',flush=True)
    freq=read(WORK/'80_logs/analysis_config.json')['frequency_Hz'];k0=2*np.pi*freq/c;q=np.linspace(-k0*np.sin(np.deg2rad(10)),k0*np.sin(np.deg2rad(10)),81);fields=[]
    for g in range(22):
        print(f'angular shell {g}/21',flush=True);fields.append(radiation(np.concatenate(nodes[g]),np.concatenate(sources[g]),q,k0));nodes[g]=sources[g]=None
    fields=np.array(fields);error=float(np.linalg.norm(fields[:,40,40,:2]-expected)/np.linalg.norm(expected));assert error<1e-7,error
    np.savez_compressed(cache,F_V=fields,q_m_inv=q,frequency_Hz=freq,center_relative_error=error,weight_closure_m3=closure)
    return dict(F_V=fields,q_m_inv=q,frequency_Hz=freq,center_relative_error=error,weight_closure_m3=closure)

def plot_angles(out,d,*,include_cuts=True):
    f=d['F_V'].cumsum(0);angles=np.rad2deg(np.arcsin(d['q_m_inv']/(2*np.pi*float(d['frequency_Hz'])/c)));ids=[0,1,5,10,20,21];intensity=abs(f[...,:2])**2/(2*Z0)
    fig,axs=plt.subplots(2,6,figsize=(19,7),layout='constrained')
    for pol in range(2):
        limit=float(intensity[ids,...,pol].max())
        for col,g in enumerate(ids):
            im=axs[pol,col].pcolormesh(angles,angles,intensity[g,...,pol],vmin=0,vmax=limit,cmap='magma',shading='auto',rasterized=True);axs[pol,col].set(aspect='equal',xlabel='theta_x (deg)',ylabel='theta_y (deg)',title=('cavity' if g==0 else ('+ pad' if g==21 else f'+ {g} cladding rings'))+' | E'+'xy'[pol])
        fig.colorbar(im,ax=axs[pol].tolist(),label='E'+'xy'[pol]+' dP/dOmega (W/sr)',ticks=np.linspace(0,limit,5),shrink=1)
    save(fig,out,'GS2b_cumulative_farfields')
    if not include_cuts:return
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for ax,axis in zip(axs,['x','y']):
        for g in ids:ax.plot(angles,intensity[g,40,:,1] if axis=='x' else intensity[g,:,40,1],label=str(g))
        ax.set(xlabel='theta_'+axis+' (deg)',ylabel='Ey dP/dOmega (W/sr)');ax.legend(title='Added rings; 21=pad',fontsize=7)
    save(fig,out,'GS2c_angular_cuts')

def check():
    rng=np.random.default_rng(17);b=rng.normal(size=(7,6))+1j*rng.normal(size=(7,6));a=np.array([1+2j,-.4j]);m=np.ones((7,6));coef,u,r=project(a[:,None,None]*b,b,m);np.testing.assert_allclose(coef,a);assert np.max(r)<1e-27
    xyz=np.array([[.12,.2,.05],[.3,.4,.04]]);e=rng.normal(size=(2,3))+1j*rng.normal(size=(2,3));q=np.linspace(-.1,.1,9)*4e6;k0=4e6;f=radiation(xyz,e,q,k0);exact=8*np.sum(e[:,1]*np.cos(k0*xyz[:,2]*1e-6))*k0*k0/(4*np.pi);np.testing.assert_allclose(f[4,4,1],exact,rtol=1e-10);assert abs(f[4,4,0])<abs(exact)*1e-12
    qx,qy=np.meshgrid(q,q);kz=np.sqrt(k0*k0-qx*qx-qy*qy);directions=np.stack([qx,qy,kz],-1)/k0;direct=np.zeros((*qx.shape,3),complex)
    for sx in [-1,1]:
        for sy in [-1,1]:
            for sz in [-1,1]:
                for pt,field in zip(xyz,e):
                    phase=np.exp(1j*(qx*sx*pt[0]+qy*sy*pt[1]+kz*sz*pt[2])*1e-6)
                    direct+=phase[...,None]*field*np.array([sx*sy,1,sy*sz])
    direct=(direct-directions*np.sum(directions*direct,-1)[...,None])*k0*k0/(4*np.pi)
    np.testing.assert_allclose(f,direct,rtol=1e-9,atol=abs(exact)*1e-11)
    axis=np.linspace(-1,1,5);ef=np.zeros((5,5,3),complex);hf=ef.copy();ef[...,1]=1+2j;hf[...,0]=-(1+2j)/Z0
    ff=plane_farfield(ef,hf,axis,1.1,q,k0);area=(5*(axis[1]-axis[0])*1e-6)**2
    np.testing.assert_allclose(ff[4,4,1],1j*k0/(2*np.pi)*np.exp(1j*k0*1.1e-6)*(1+2j)*area,rtol=1e-12)
    np.testing.assert_allclose(plane_farfield(ef*.3,hf*.3,axis,1.1,q,k0),ff*.3,rtol=1e-12,atol=1e-20)
    print('PASS: complex projection; mirror/z-retarded radiation; upward E/H phase, units and linearity',flush=True)



def driven_comparison(out):
    cfg=read(out/'99_config/geometry_and_run_config.json');k0=2*np.pi*cfg['frequency_Hz']/c
    q=np.linspace(-k0*np.sin(np.deg2rad(10)),k0*np.sin(np.deg2rad(10)),81);theta=np.rad2deg(np.arcsin(q/k0))
    periodic=load('periodic_k00_EH.npz');basis=np.concatenate([periodic['E_normalized'][1],periodic['H_normalized'][1]],-1);weights=metric(periodic)
    target=load('finite_cavity_EH_o3.npz');target_v=np.concatenate([target['E_V_m'],target['H_A_m']],-1);target_u=float(np.sum(abs(target_v)**2*weights))
    ref_air=load('finite_air_1100nm_fine_EH.npz');ref_f=plane_farfield(np.nan_to_num(ref_air['E_V_m'][::2,::2]),np.nan_to_num(ref_air['H_A_m'][::2,::2]),ref_air['x_um'][::2],float(ref_air['z_um']),q,k0)
    far=[];slabs=[];airs=[];coeff=[];rows=[];checks={}
    for label in ['A','B','C']:
        with np.load(out/'01_results'/f'{label}_cavity_EH.npz') as z:v=np.concatenate([z['E_V_m'],z['H_A_m']],-1)
        aj,u,res=project(v,basis,weights);aa=aj.mean();total_u=float(u.sum());gamma_u=float(len(aj)*abs(aa)**2)
        overlap=float(abs(np.sum(v*target_v.conj()*weights))**2/(total_u*target_u))
        with np.load(out/'01_results'/f'{label}_air_EH.npz') as z:air={k:z[k] for k in z.files}
        with np.load(out/'01_results'/f'{label}_air_check_EH.npz') as z:high={k:z[k] for k in z.files}
        with np.load(out/'01_results'/f'{label}_slab_EH.npz') as z:slab={k:z[k] for k in z.files}
        f=plane_farfield(air['E_V_m'],air['H_A_m'],air['x_um'],float(air['z_um']),q,k0);fh=plane_farfield(high['E_V_m'],high['H_A_m'],high['x_um'],float(high['z_um']),q,k0)
        intensity=np.sum(abs(f)**2,-1);reference=np.sum(abs(ref_f)**2,-1);shape=float(np.vdot(intensity,reference).real/(np.linalg.norm(intensity)*np.linalg.norm(reference)))
        center=float(intensity[40,40]/(2*Z0));central_max=bool(intensity[40,40]>=intensity.max()*(1-1e-8))
        completion=read(out/'80_logs'/f'{label}_completed.json');rows.append(dict(model=label,frequency_Hz=cfg['frequency_Hz'],source_norm_A2_per_m=completion['source_norm_A2_per_m'],source_delivered_power_W=completion['source_delivered_power_W'],cavity_U_J=total_u,gamma_re_sqrtJ=float(aa.real),gamma_im_sqrtJ=float(aa.imag),Gamma_global_energy_J=gamma_u,Gamma_global_fraction=gamma_u/total_u,Gamma_local_fraction=float(1-res.sum()/total_u),mode19_energy_overlap=overlap,mode19_farfield_intensity_cosine=shape,center_W_sr=center,center_over_cavity_U_s_inv_sr=center/total_u,center_over_Gamma_U_s_inv_sr=center/gamma_u if gamma_u/total_u>=cfg['acceptance']['small_projection_fraction'] else None,center_is_global_max_in_10deg=central_max))
        checks[label]=dict(air_height_center_relative_error=float(np.linalg.norm(f[40,40]-fh[40,40])/np.linalg.norm(f[40,40])),air_height_map_relative_error=float(np.linalg.norm(f-fh)/np.linalg.norm(f)),source_norm_relative_error=float(abs(completion['source_norm_A2_per_m']/cfg['source_norm_target_A2_per_m']-1)),mode19_match=overlap>=cfg['acceptance']['mode19_energy_overlap_min'] and shape>=cfg['acceptance']['farfield_intensity_cosine_min'])
        far.append(f);airs.append(air);slabs.append(slab);coeff.append(aa)
    far=np.array(far);pairs=[(0,1),(2,1),(0,2)];difference=np.array([far[i]-far[j] for i,j in pairs]);linear=[]
    for i,j in pairs:
        diff=plane_farfield(airs[i]['E_V_m']-airs[j]['E_V_m'],airs[i]['H_A_m']-airs[j]['H_A_m'],airs[i]['x_um'],float(airs[i]['z_um']),q,k0)
        linear.append(float(np.linalg.norm(diff-(far[i]-far[j]))/max(np.linalg.norm(diff),np.finfo(float).tiny)))
    assert max(linear)<1e-10
    raw=abs(far[...,:2])**2/(2*Z0);diff_int=abs(difference[...,:2])**2/(2*Z0)
    np.savez_compressed(out/'01_results/driven_comparison.npz',models=['A','B','C'],F_V=far,difference_F_V=difference,differences=['A-B','C-B','A-C'],q_m_inv=q,frequency_Hz=cfg['frequency_Hz'],reference_mode19_F_V=ref_f,Gamma_coefficients=coeff,source_phase_common=True)
    pd.DataFrame(rows).to_csv(table_path(out / '80_logs', 'source_and_projection.csv'),index=False)
    center_rows=[]
    for label,f in zip(['A','B','C','A-B','C-B','A-C'],np.concatenate([far,difference])):
        for pol,ch in enumerate(['Ex','Ey']):center_rows.append(dict(model=label,channel=ch,Re_F_V=f[40,40,pol].real,Im_F_V=f[40,40,pol].imag,self_W_sr=abs(f[40,40,pol])**2/(2*Z0)))
    pd.DataFrame(center_rows).to_csv(table_path(out / '80_logs', 'driven_center_amplitudes.csv'),index=False)
    for stem,maps,names in [('GS3c_same_source_farfields',raw,['A','B','C']),('GS3d_difference_farfields',diff_int,['A-B','C-B','A-C'])]:
        fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained');limit=float(maps.max())
        for pol in range(2):
            for j,label in enumerate(names):
                im=axs[pol,j].pcolormesh(theta,theta,maps[j,...,pol],vmin=0,vmax=limit,cmap='magma',shading='auto',rasterized=True);axs[pol,j].set(aspect='equal',title=label+' | E'+'xy'[pol],xlabel='theta_x (deg)',ylabel='theta_y (deg)')
        fig.colorbar(im,ax=axs.ravel().tolist(),label='dP/dOmega (W/sr)',shrink=.8);save(fig,out,stem)
    fields=np.array([x['E_V_m'] for x in slabs]);delta=np.array([fields[i]-fields[j] for i,j in pairs]);axis=slabs[0]['x_um']
    for stem,values,names in [('GS3b_same_source_complex_fields',fields,['A','B','C']),('GS3e_difference_complex_fields',delta,['A-B','C-B','A-C'])]:
        fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained');limit=float(abs(values[...,:2].real).max())
        for pol in range(2):
            for j,label in enumerate(names):
                im=axs[pol,j].pcolormesh(axis,axis,values[j,...,pol].real,vmin=-limit,vmax=limit,cmap='RdBu_r',shading='auto',rasterized=True);axs[pol,j].set(aspect='equal',title=label+' | Re(E'+'xy'[pol]+')',xlabel='x (um)',ylabel='y (um)')
        fig.colorbar(im,ax=axs.ravel().tolist(),label='Re(E) (V/m)',shrink=.8);save(fig,out,stem)
    fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    for j,label in enumerate(['A','B','C']):
        ey=far[j,40,40,1];axs[0].quiver(0,0,ey.real,ey.imag,angles='xy',scale_units='xy',scale=1,color=['#b33b3b','#326aaf','#458f5c'][j]);axs[0].text(ey.real,ey.imag,label)
        axs[1].plot(theta,raw[j,40,:,1]/rows[j]['cavity_U_J'],label=label)
        if rows[j]['center_over_Gamma_U_s_inv_sr'] is not None:axs[2].plot(theta,raw[j,40,:,1]/rows[j]['Gamma_global_energy_J'],label=label)
    allcenter=far[:,40,40,1];bound=max(abs(allcenter.real).max(),abs(allcenter.imag).max())*1.15;axs[0].set(xlim=(-bound,bound),ylim=(-bound,bound),aspect='equal',xlabel='Re(Fy) (V)',ylabel='Im(Fy) (V)',title='Common-source center phasors')
    for ax,title in [(axs[1],'Per cavity stored energy'),(axs[2],'Per global Gamma projection energy')]:ax.set(xlabel='theta_x (deg)',ylabel='Ey dP/dOmega / U (1/s/sr)',title=title);ax.legend()
    save(fig,out,'GS3f_center_and_normalized_radiation')
    checks['source_norm_relative_spread']=float(np.ptp([r['source_norm_A2_per_m'] for r in rows])/rows[1]['source_norm_A2_per_m'])
    checks.update(linearity_errors=linear,source_phase_independent_alignment=False,normalization_is_not_equal_internal_field=True,frequency_scan_executed=False,window_sensitivity_executed=False,mesh_convergence_executed=False)
    put(out/'80_logs/driven_validation.json',checks)
    status="A \u7684\u5185\u90e8\u573a\u4e0e\u8fdc\u573a\u5f62\u72b6\u901a\u8fc7 mode19 \u5bf9\u5e94\u6027\u9608\u503c\uff1b\u8fd9\u4ec5\u786e\u8ba4\u6bd4\u8f83\u5bf9\u8c61\uff0c\u5c1a\u4e0d\u80fd\u5355\u72ec\u8bc1\u660e cladding \u673a\u5236\u3002" if checks['A']['mode19_match'] else "A \u7684\u53d7\u8feb\u573a\u5c1a\u672a\u540c\u65f6\u901a\u8fc7 mode19 \u5185\u90e8\u573a\u4e0e\u8fdc\u573a\u5f62\u72b6\u9608\u503c\uff1b\u5f53\u524d\u5bf9\u7167\u4e0d\u80fd\u76f4\u63a5\u89e3\u91ca\u76ee\u6807\u672c\u5f81\u6a21\u5f0f\uff0c\u9700\u68c0\u67e5\u6e90\u7a97\u53e3\u4e0e\u5171\u540c\u9891\u7387\u3002"
    report="# S5 \u540c\u6e90\u540c\u9891\u4e09\u7ed3\u6784\u5bf9\u7167\\n\\n".replace('\\n','\n')+status+'\n\n'
    report+="A \u4e3a\u5b9e\u9645 cavity\u2013cladding\uff0cB \u5c06 cavity \u6676\u683c\u5ef6\u62d3\u5230\u539f cladding \u533a\uff0cC \u5c06\u539f cladding \u533a\u6539\u4e3a\u65e0\u5b54\u5bbf\u4e3b\u8584\u677f\u3002\u5171\u540c\u9891\u7387\u4e3a "+str(cfg['frequency_Hz']/1e12)+" THz\uff0cmesh=5\u3002\u6e90\u7684\u5e45\u5ea6\u3001\u76f8\u4f4d\u3001\u7a97\u53e3\u548c\u652f\u6301\u533a\u57df\u76f8\u540c\uff1b\u76f8\u540c\u6e90\u8303\u6570\u4e0d\u4ee3\u8868\u76f8\u540c\u8f93\u5165\u529f\u7387\u3002"+'\n\n'
    report+=pd.DataFrame(rows).to_csv(index=False)+'\n\n'
    report+="\u539f\u59cb\u540c\u6e90\u5dee\u573a\u4fdd\u7559\u516c\u5171\u76f8\u4f4d\uff0c\u672a\u5bf9\u4e09\u7ec4\u5355\u72ec\u8c03\u76f8\u6216\u5f52\u4e00\u5316\u3002\u7a7a\u6c14\u9762 E/H \u4f7f\u7528\u540c\u4e00\u5411\u4e0a\u4f20\u64ad\u7b97\u5b50\uff0c\u65f6\u95f4\u7ea6\u5b9a exp(+i omega t)\uff0c\u53c2\u8003\u9762 z=0\uff1b\u5148\u51cf\u573a\u518d\u53d8\u6362\u4e0e\u5148\u53d8\u6362\u518d\u51cf\u573a\u7684\u8bef\u5dee\u89c1 driven_validation.json\u3002\u5dee\u573a\u5f3a\u5ea6\u4e0d\u7b49\u4e8e\u4e24\u7ec4\u5f3a\u5ea6\u4e4b\u5dee\u3002"+'\n\n'
    report+="\u6bcf\u5355\u4f4d cavity \u50a8\u80fd\u3001\u6bcf\u5355\u4f4d\u5168\u5c40 \u0393 \u6295\u5f71\u80fd\u91cf\u7684\u8f90\u5c04\u662f\u7b2c\u4e8c\u5957\u6bd4\u8f83\uff0c\u4e0d\u66ff\u4ee3\u539f\u59cb\u5dee\u573a\u3002\u6e90\u8303\u6570\u8bef\u5dee\u3001\u4e24\u4e2a\u7a7a\u6c14\u9762\u7684\u4e00\u81f4\u6027\u548c mode19 \u5bf9\u5e94\u6027\u987b\u5206\u522b\u6838\u9a8c\u3002\u5355\u9891\u5dee\u5f02\u53ef\u80fd\u5305\u542b\u5171\u632f\u4f4d\u7f6e\u6539\u53d8\uff0c\u6e90\u7a97\u53e3\u4e0e\u9891\u7387\u7a33\u5065\u6027\u5c1a\u672a\u5b8c\u6210\u3002Ex \u6cd5\u5411\u76f8\u6d88\u7531 quarter \u5bf9\u79f0\u6027\u7ea6\u675f\uff0c\u4e0d\u80fd\u89c6\u4e3a\u72ec\u7acb\u5168\u6a21\u578b\u9a8c\u8bc1\u3002"+'\n'
    (out/'12_reports/driven_comparison.md').write_text(report,encoding='utf8')
    validation_path=out/'80_logs/validation.json'
    validation=read(validation_path) if validation_path.exists() else {}
    validation.update(stage_B_complete=True,initial_driven_solves_completed=3,driven_checks=checks,scientific_mechanism_accepted=False,interpretation='Single-window single-frequency results require conditional interpretation and follow-up robustness checks')
    put(validation_path,validation)
    readme=out/'12_reports/README.md'
    with readme.open('a',encoding='utf8') as handle:handle.write(table_links('\n三组首轮计算完成：[比较报告](driven_comparison.md)、[源与投影指标](../80_logs/source_and_projection.csv)、[同源远场](../10_overview/GS3c_same_source_farfields.png)。稳健性与科学验收仍待核查。\n'))



def source_and_geometry_figures(out):
    from scipy.interpolate import RegularGridInterpolator
    from matplotlib.path import Path as PolygonPath
    cfg=read(out/'99_config/geometry_and_run_config.json');old=read(out/'99_config/historical_config.json');amp=read(out/'99_config/source_amplitude.json')['Jamp_A_m2']
    raw=np.loadtxt(out/'01_results/source_template.txt',comments='%');values=(raw[:,3:6]+1j*raw[:,6:9]).reshape(65,65,9,3)
    interp=RegularGridInterpolator((np.linspace(-.5,.5,65),np.linspace(-.5,.5,65),np.linspace(0,.1,9)),values)
    axis=np.linspace(-cfg['source_radius_um']*1.05,cfg['source_radius_um']*1.05,601);xx,yy=np.meshgrid(axis,axis);x=abs(xx);y=abs(yy);u=(x-y/np.sqrt(3))/old['a'];v=2*y/(np.sqrt(3)*old['a']);r=np.hypot(xx,yy);window=np.where(r<cfg['source_radius_um'],np.cos(np.pi*r/(2*cfg['source_radius_um']))**2,0)
    field=interp(np.column_stack([(u-np.floor(u+.5)).ravel(),(v-np.floor(v+.5)).ravel(),np.zeros(xx.size)])).reshape(*xx.shape,3);field*=np.where(xx[...,None]<0,[-1,1,1],1);field*=np.where(yy[...,None]<0,[-1,1,-1],1);field*=amp*window[...,None]
    fig,axs=plt.subplots(1,3,figsize=(13,4.5),layout='constrained');lim=float(abs(field[...,:2].real).max())
    for j,(value,title,cmap,low,high,label) in enumerate([(window,'Common xy source window','magma',0,1,'Window (1)'),(field[...,0].real,'Prescribed Re(Jx), z=0','RdBu_r',-lim,lim,'Current density (A/m2)'),(field[...,1].real,'Prescribed Re(Jy), z=0','RdBu_r',-lim,lim,'Current density (A/m2)')]):
        im=axs[j].pcolormesh(axis,axis,value,cmap=cmap,vmin=low,vmax=high,shading='auto',rasterized=True);axs[j].set(aspect='equal',xlabel='x (um)',ylabel='y (um)',title=title);fig.colorbar(im,ax=axs[j],label=label)
    save(fig,out,'GS3_source')
    axis=np.linspace(0,old['finite_Lx']/2,801);xx,yy=np.meshgrid(axis,axis);grid=np.column_stack([xx.ravel(),yy.ravel()]);arrays=[]
    for label in ['A','B','C']:
        geometry=read(out/'99_config'/f'{label}_simulation_config.json');valid=PolygonPath(geometry['footprint']).contains_points(grid,radius=1e-10).reshape(xx.shape);eps=np.where(valid,float(old['simulation_common']['refractive_index'])**2,1.)
        for poly in geometry['layers'][0]['holes']:
            poly=np.asarray(poly);ix=np.flatnonzero((axis>=poly[:,0].min())&(axis<=poly[:,0].max()));iy=np.flatnonzero((axis>=poly[:,1].min())&(axis<=poly[:,1].max()))
            if not len(ix) or not len(iy):continue
            gx,gy=np.meshgrid(axis[ix],axis[iy]);inside=PolygonPath(poly).contains_points(np.column_stack([gx.ravel(),gy.ravel()])).reshape(gx.shape);block=eps[np.ix_(iy,ix)];block[inside]=1.;eps[np.ix_(iy,ix)]=block
        arrays.append(eps)
    fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained');limit=float(old['simulation_common']['refractive_index'])**2-1
    for ax,delta,title in zip(axs,[arrays[0]-arrays[1],arrays[2]-arrays[1]],['epsilon(A) - epsilon(B)','epsilon(C) - epsilon(B)']):
        im=ax.pcolormesh(axis,axis,delta,vmin=-limit,vmax=limit,cmap='RdBu_r',shading='auto',rasterized=True);ax.set(aspect='equal',xlabel='x (um)',ylabel='y (um)',title=title)
        boundary=np.asarray(old['inner_boundary']);ax.plot(*np.vstack([boundary,boundary[0]]).T,color='k',lw=.5);ax.set_xlim(0,33);ax.set_ylim(0,29)
    fig.colorbar(im,ax=axs.tolist(),label='Relative permittivity difference (1)');save(fig,out,'GS3_geometry_difference')
    np.savez_compressed(out/'01_results/geometry_difference_preview.npz',axis_um=axis,epsilon_A=arrays[0],epsilon_B=arrays[1],epsilon_C=arrays[2],display_raster_only=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--check',action='store_true');p.add_argument('--skip-angular',action='store_true');args=p.parse_args()
    if args.check:check();return
    out=args.output.resolve();assert out.is_dir()
    state=local_projection(out);f=center_budget(out)
    if not args.skip_angular:
        d=angular_shells(out,f);plot_angles(out,d);state['shell_center_reproduction_relative_error']=float(d['center_relative_error'])
    state.update(status='complete' if not args.skip_angular else 'projection_and_center_only',COMSOL_started=False,physical_closure=read(WORK/'80_logs/unitcell_center_radiation_checks.json'),source_regions_are_disjoint=True,Ex_cancellation_imposed_by_quarter_symmetry=True)
    put(out/'80_logs/stage_A_validation.json',state)
    (out/'12_reports/stage_A.md').write_text('# S5 阶段 A\n\n仅使用保存的实际复场。cavity 局部 Γ-py 单基底投影可解释能量比例 '+f"{state['local_Gamma_explained_fraction']:.6%}"+'；全 cavity 单一 Γ 系数可解释比例 '+f"{state['global_single_Gamma_explained_fraction']:.6%}"+'。局部系数保留公共相位；逐胞包络不同，不能称为严格纯 Γ 波。\n\nGS1：三维 E/H 能量内积得到局部 Γ 投影与残差，场图位于最靠近 z=0 的已保存求积面。cladding 显示实际场，不复制周期参考。GS2：互不重叠区域的实际体源相干累加；中心积分复用经过几何裁切修正的结果，角分布采用相同四面体样本的仿射源重构和精确几何积分权重。色标为共同线性尺度，包含交叉项。\n\n这些结果描述同一个本征解中的贡献，并非删除区域后的解。Ex 法向相消由既有对称边界约束，不能作为独立全模型检验。体源与独立空气场仍保留原有物理闭合误差，几何求和精确不等于场已收敛。\n',encoding='utf8')
    print(json.dumps({k:v for k,v in state.items() if k!='physical_closure'},indent=2),flush=True)

if __name__=='__main__':main()
