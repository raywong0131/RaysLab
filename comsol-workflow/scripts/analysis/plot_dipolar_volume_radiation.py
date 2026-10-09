"""Independent electric-volume projection and radiation of the original 19 k points."""
from __future__ import annotations

from comsol_workflow.output_paths import output_path
from comsol_workflow.figure_output import save_figure_formats

import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scripts.analysis.visualize_dipolar_broadening_2d import OUTPUT, SOURCE
from comsol_workflow.dipolar_lattice_analysis import load_data
from comsol_workflow.dipolar_radiation import planar_K0_to_far, pairs, radiation_budget
from comsol_workflow.result_io import write_json, write_csv, digest, read_json as json_read
from scripts.run_main.run_dipolar_periodic_reference import OUT

CAVITY_SUPPORT = 'cavity cells only; basis zero in cladding and outer pad'


def fit_cells(electric, finite, k, rho, w, masks, regions=(0,)):
    """Factorized complex least squares, algebraically equal to the full dictionary."""
    if tuple(regions) != (0,):
        raise ValueError('Cavity periodic fields may only cover cavity cells (region 0)')
    gram=np.zeros((76,76),complex);right=np.zeros(76,complex);field_norm=0.
    basis=electric.reshape(76,len(rho),3)
    for region in regions:
        field=finite[f'E_{region}'];centers=finite[f'centers_{region}_um']
        phase=np.exp(-1j*centers@k.T)
        metric=2*w*(1+9.89*masks[region])*1e18
        cell_gram=np.einsum('bnp,cnp,n->bc',basis.conj(),basis,metric,optimize=True)
        lattice_gram=phase.conj().T@phase
        gram+=cell_gram*np.repeat(np.repeat(lattice_gram,4,axis=0),4,axis=1)
        fourier=(phase.conj().T@field.reshape(len(centers),-1)).reshape(19,len(rho),3)
        right+=np.einsum('kbnp,knp,n->kb',electric.conj(),fourier,metric,optimize=True).ravel()
        field_norm+=float(np.einsum('rnp,rnp,n->',field.conj(),field,metric,optimize=True).real)
    norms=np.sqrt(np.diag(gram).real)
    normalized=gram/norms[:,None]/norms[None,:]
    a=np.linalg.solve(normalized,right/norms)/norms
    residual=max(0.,field_norm-2*np.vdot(a,right).real+np.vdot(a,gram@a).real)
    return a.reshape(19,4),float(np.sqrt(residual/field_norm)),float(np.linalg.cond(normalized))


def responses(q,k,electric,rho,w,masks,centers,omega,region_weights=None):
    """Cavity-supported Bloch radiation; +iwt, including the slab z mirror.

    Cladding arrays are never accessed. Dropping zero cavity weights is exact.
    """
    k0=omega/299792458.;bands=electric.shape[1];count=19*bands
    material_w=9.89*(w*masks[0] if region_weights is None else region_weights[0])
    if not np.isfinite(material_w).all() or np.any(material_w<0):
        raise ValueError('Cavity quadrature weights must be finite and nonnegative')
    active=material_w>0
    result=np.zeros((len(q),19,bands,2),complex)
    if not active.any() or len(centers[0])==0:
        return result
    rho=rho[active];material_w=material_w[active]
    basis=electric[:,:,active,:].reshape(count,len(rho),3)
    components=[np.ascontiguousarray(basis[:,:,c].T) for c in range(3)]
    bloch_phase=np.exp(-1j*centers[0]@k.T)
    for start in range(0,len(q),128):
        sl=slice(start,start+128);qq=q[sl];kz=np.sqrt(k0*k0-np.sum(qq*qq,axis=1)*1e12)
        n=np.column_stack([qq*1e6,kz])/k0
        kernel=np.exp(1j*qq@rho[:,:2].T)*material_w[None,:]
        cosine=2*np.cos(kz[:,None]*rho[None,:,2]*1e-6)
        sine=2j*np.sin(kz[:,None]*rho[None,:,2]*1e-6)
        p=np.empty((len(qq),count,3),complex)
        for component in range(3):
            p[:,:,component]=(kernel*(sine if component==2 else cosine))@components[component]
        p-=n[:,None,:]*np.einsum('qc,qbc->qb',n,p)[:,:,None]
        lattice=np.exp(1j*qq@centers[0].T)@bloch_phase
        result[sl]=k0*k0/(4*np.pi)*p[:,:,:2].reshape(len(qq),19,bands,2)*lattice[:,:,None,None]
        if len(q)>128 and (start%1024==0 or start+128>=len(q)):
            print(f'Cavity radiation: {min(start+128,len(q))}/{len(q)} observation points',flush=True)
    return result


def save_cavity_support(a,t,error,condition,direct,full):
    terms=a[:,:,None]*t;total=terms.sum(axis=(0,1))
    summary={'basis_support':CAVITY_SUPPORT,
             'electric_L2_residual_inside_cavity':error,'Gram_condition':condition,
             'K0_total_F_V':pairs(total),'K0_total_over_air_Ey':pairs(total/direct[1]),
             'nonGamma_K0_total_F_V':pairs(terms[1:].sum(axis=(0,1))),
             'Gamma_py_F_V':pairs(terms[0,1]),
             'Gamma_other_bands_F_V':pairs(np.delete(terms[0],1,axis=0).sum(axis=0)),
             'radiation_residual_including_cladding_and_pad_F_V':pairs(full-total),
             'scope':'active cavity-only reference with original 19 k and mesh5 fields'}
    write_json(output_path(OUTPUT, 's5_cavity_support_control.json'),summary)
    np.savez_compressed(output_path(OUTPUT, 's5_cavity_support_control.npz'),A=a,T0=t,terms=terms,basis_support=CAVITY_SUPPORT)


def plot_support_control(d,summary):
    convert=lambda values: np.asarray(values)[:,0]+1j*np.asarray(values)[:,1]
    direct=convert(summary['direct_air_F']);full=convert(summary['full_volume_F'])
    cavity=convert(summary['basis_over_air'])
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    index_map(axes[0,0],d['k19']/d['scale'])
    axes[0,0].set(xlabel=r'$k_x/(2\pi/a)$',ylabel=r'$k_y/(2\pi/a)$')
    ax=axes[0,1];error=summary['electric_L2_residual_inside_cavity']
    ax.bar([0],[100*error],color='#28618a',width=.5)
    ax.text(0,100*error+2,f'{error:.2%}',ha='center')
    ax.set(title='Electric reconstruction inside cavity',xticks=[0],xticklabels=['Cavity 19 x 4'],
           ylabel='Relative L2 residual (%)',ylim=(0,100),xlim=(-.8,.8))
    ax.text(.5,.83,f'{d["N"]} cavity cells; zero support in cladding',transform=ax.transAxes,ha='center',fontsize=9)
    values=np.array([cavity,full/direct[1]-cavity,full/direct[1],direct/direct[1]])
    for pol,name in enumerate(['Ex','Ey']):
        ax=axes[1,pol];x=np.arange(4)
        ax.bar(x-.18,values[:,pol].real,.36,color='#28618a',label='Real')
        ax.bar(x+.18,values[:,pol].imag,.36,color='#b44c42',label='Imaginary')
        ax.set_xticks(x,['Cavity 19','Remainder','Full volume','Direct air'],rotation=20,ha='right')
        ax.set(title=f'{name}: center amplitude',ylabel=r'$F_\alpha(0)/F_y^{air}(0)$')
        ax.axhline(0,color='.6',lw=.6);ax.legend(frameon=False)
    fig.suptitle('Cavity-only reference | complete finite radiation is an independent comparison',fontsize=13)
    save(fig,'17_s5_support_control')


def plot_nonGamma_py(d,summary):
    raw=np.asarray(summary['terms']);cube=raw[...,0]+1j*raw[...,1]
    ref=complex(*summary['direct_air_F'][1]);values=cube[1:,1]/ref
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    index_map(axes[0],d['k19']/d['scale']);axes[0].set(xlabel=r'$k_x/(2\pi/a)$',ylabel=r'$k_y/(2\pi/a)$')
    for pol,name in enumerate(['Ex','Ey']):
        ax=axes[pol+1]
        ax.scatter(0,0,color='#28618a',s=25)
        ax.axhline(0,color='.8',lw=.6);ax.axvline(0,color='.8',lw=.6)
        ax.set(xlim=(-.01,.01),ylim=(-.01,.01),aspect='equal',xlabel=r'Re$(d/F_y^{air})$',ylabel=r'Im$(d/F_y^{air})$',
               title=f'{name}: exact cavity-window null at K=0')
        ax.set_xticks([-.01,0,.01]);ax.set_yticks([-.01,0,.01])
        ax.text(.04,.92,f'Max numerical |d_n/R| = {abs(values[:,pol]).max():.2g}',transform=ax.transAxes,fontsize=9)
        ax.text(.04,.82,'No physical phase or coherence ratio at a zero',transform=ax.transAxes,fontsize=8)
    fig.suptitle('18 non-Gamma py contributions | cavity cells only | same amplitude scale',fontsize=13)
    fig.supxlabel('Finite-precision residuals are recorded in the data; they are not enlarged into interference arrows.',fontsize=9)
    save(fig,'18_s5_nonGamma_py_interference')


def run(order=3):
    d=load_data(SOURCE);k=d['k19'];path=output_path(OUTPUT, f's5_finite_volume_o{order}.npz')
    finite=np.load(path);rho=finite['rho_xyz_um'];w=finite['weights_half_m3'];masks=finite['material_masks']
    omega=complex(finite['omega_rad_s']);centers=[finite[f'centers_{r}_um'] for r in range(2)]
    fields=[];provenance={str(path):digest(path)}
    for index in range(19):
        p=OUT/f'01_results/k{index:02d}_volume_o{order}.npz';data=np.load(p)
        np.testing.assert_allclose(data['rho_xyz_um'],rho,atol=1e-13,rtol=0)
        np.testing.assert_allclose(data['k_um_inv'],k[index],atol=1e-13,rtol=0)
        fields.append(data['E_native_over_Hz_rms']);provenance[str(p)]=digest(p)
    electric=np.asarray(fields)
    print('Projecting cavity electric field only on 19 x 4 Bloch fields',flush=True)
    a,error,condition=fit_cells(electric,finite,k,rho,w,masks)
    radiation_fields=[]
    for index in range(19):
        p=OUT/f'01_results/k{index:02d}_mesh_volume.npz';data=np.load(p)
        np.testing.assert_allclose(data['k_um_inv'],k[index],atol=1e-13,rtol=0)
        radiation_fields.append(data['E_native_over_Hz_rms']);provenance[str(p)]=digest(p)
        if index==0:
            radiation_rho=data['rho_xyz_um'];radiation_weights=data['region_weights_half_m3']
        else:
            np.testing.assert_allclose(data['rho_xyz_um'],radiation_rho,rtol=0,atol=1e-13)
            np.testing.assert_allclose(data['region_weights_half_m3'],radiation_weights,rtol=0,atol=1e-32)
    radiation_fields=np.asarray(radiation_fields)
    t=responses(np.zeros((1,2)),k,radiation_fields,radiation_rho,None,None,centers,omega,radiation_weights)[0]
    terms,total,inc,coh,cross=radiation_budget(a,t);cube=terms.reshape(19,4,2)
    full=finite['full_volume_F_V'];direct=planar_K0_to_far(finite['direct_air_C_V_m'],float(finite['z_air_m']),omega)
    scale=abs(direct[1])
    lattice0=np.exp(-1j*centers[0]@k.T).mean(axis=0)
    np.testing.assert_allclose(lattice0,np.eye(1,19)[0],atol=1e-12,rtol=0)
    metadata=json.loads(str(finite['metadata_json']))
    metadata['basis_support']=CAVITY_SUPPORT
    metadata['projection_region']='cavity'
    metadata['cavity_cell_count']=len(centers[0])
    groups={'Gamma_py':cube[0,1], 'nonGamma_py':cube[1:,1].sum(axis=0),
            'other_bands':np.delete(cube,1,axis=1).sum(axis=(0,1)),
            'unreconstructed_field_including_pad':full-total}
    np.testing.assert_allclose(sum(groups.values()),full,rtol=1e-10,atol=scale*1e-10)
    summary={'mesh':5,'zeta':1.156,'point_count':19,'bands_per_k':4,'units':'V; spherical outgoing F',
             'status':'fixed_mesh_calculation_complete; interpret with displayed residuals',
             'metadata':metadata,'basis_support':CAVITY_SUPPORT,
             'electric_L2_residual_inside_cavity':error,'L2_residual_domain':'cavity cells only',
             'full_vector_L2_residual':error,'normalized_Gram_condition':condition,
             'volume_vs_air_relative_error':float(np.linalg.norm(full-direct)/np.linalg.norm(direct)),
             'basis_vs_air_relative_error':float(np.linalg.norm(total-direct)/np.linalg.norm(direct)),
             'basis_over_air':pairs(total/direct[1]),'full_volume_F':pairs(full),
             'direct_air_F':pairs(direct),
             'radiation_quadrature':'original cavity FEM tetrahedra; 3x3x3 Duffy integration; zero cavity weights omitted',
             'radiation_quadrature_check':json_read(OUT/'01_results/s5_mesh_quadrature_check.json'),
             'groups':{key:pairs(value) for key,value in groups.items()},
             'terms':pairs(cube),'coherent':coh.tolist(),'incoherent_76_terms':inc.tolist(),
             'incoherent_19_points':np.sum(abs(cube.sum(axis=1))**2,axis=0).tolist(),
             'cross_19_points':(coh-np.sum(abs(cube.sum(axis=1))**2,axis=0)).tolist(),
             'nonGamma_py_coherence_ratio':[None,None],
             'nonGamma_origin_status':'exact cavity-window zero; computed residual is roundoff, not a physical phase',
             'cavity_window_origin_max_nonGamma':float(abs(lattice0[1:]).max()),
             'source_sha256':provenance,
             'scope':'19 original k; cavity-reference quartet only on cavity cells; cladding and outer pad belong to the complete-field remainder',
             'limitations':['fixed mesh5 does not prove an exact periodic zero',
                            'original 19-point set is not mirror closed',
                            'cavity reference is identically zero outside cavity cells',
                            'projection uses positive material-resolved quadrature; radiation uses original FEM tetrahedra',
                            'residual includes cavity reconstruction error, actual cladding field and outer pad; it is not uniquely attributed']}
    rows=[]
    for i in range(19):
        rows.append({'index':i,'kx_um_inv':k[i,0],'ky_um_inv':k[i,1],
                     **{f'{c}_{part}':value for c,v in zip(['Ex','Ey'],cube[i].sum(axis=0))
                        for part,value in [('re_V',v.real),('im_V',v.imag)]},
                     **{f'py_{c}_{part}':value for c,v in zip(['Ex','Ey'],cube[i,1])
                        for part,value in [('re_V',v.real),('im_V',v.imag)]}})
    write_csv(output_path(OUTPUT, 's5_volume_K0_terms.csv'),rows)
    write_json(output_path(OUTPUT, 's5_volume_summary.json'),summary)
    save_cavity_support(a,t,error,condition,direct,full)
    axis=np.linspace(-.08,.08,101);xx,yy=np.meshgrid(axis,axis)
    q=np.column_stack([xx.ravel(),yy.ravel()])*d['scale']
    print('Computing 2D radiation response from material volume',flush=True)
    weighted_fields=np.einsum('kb,kbnp->knp',a,radiation_fields)[:,None,:,:]
    contributions=responses(q,k,weighted_fields,radiation_rho,None,None,centers,omega,radiation_weights)
    direct_map=air_spectrum(axis*d['scale'],omega,float(finite['z_air_m']),summary['metadata'])
    np.testing.assert_allclose(direct_map[50,50],direct,rtol=1e-10,atol=scale*1e-12)
    np.savez_compressed(output_path(OUTPUT, 's5_volume_radiation.npz'),k_um_inv=k,A=a,T0=t,K0_terms=cube,
                        K_axis_over_2pi_a=axis,contributions=contributions.reshape(101,101,19,1,2),
                        direct_air_F=direct,direct_air_map_F=direct_map,full_volume_F=full,basis_support=CAVITY_SUPPORT,**groups)
    plot(d,axis,contributions.reshape(101,101,19,1,2),cube,direct,full,summary,direct_map)
    print(json.dumps({key:summary[key] for key in ['full_vector_L2_residual','normalized_Gram_condition',
          'volume_vs_air_relative_error','basis_vs_air_relative_error','basis_over_air']},indent=2),flush=True)


def air_spectrum(axis,omega,z_air,metadata):
    from comsol_workflow.farfield_fft import load_air_field
    from pathlib import Path
    air=load_air_field(Path(json_read(output_path(OUTPUT, 'electric_19_summary.json'))['source']))
    px=np.exp(1j*axis[:,None]*air['x_axis'][None,:])
    py=np.exp(1j*axis[:,None]*air['y_axis'][None,:])
    dx=np.mean(np.diff(air['x_axis']));dy=np.mean(np.diff(air['y_axis']))
    spectra=[]
    for component in ['Ex','Ey']:
        field=air[component]
        if metadata['conjugate_air_error']<metadata['native_air_error']:field=field.conj()
        spectra.append(py@field@px.T*dx*dy*1e-12)
    xx,yy=np.meshgrid(axis,axis);kz=np.sqrt((omega/299792458.)**2-(xx**2+yy**2)*1e12)
    return np.stack(spectra,axis=-1)*(1j*kz/(2*np.pi)*np.exp(1j*kz*z_air))[:,:,None]


def align_bar(fig,bar,axes):
    bar.set_label(bar.ax.get_ylabel(),rotation=270,labelpad=18)
    bar.ax.tick_params(which='both',direction='in')
    fig.canvas.draw();fig.set_layout_engine(None)
    boxes=[ax.get_position() for ax in np.asarray(axes).ravel()]
    box=bar.ax.get_position();bottom=min(b.y0 for b in boxes);top=max(b.y1 for b in boxes)
    bar.ax.set_position([box.x0,bottom,box.width,top-bottom])


def save(fig,name):
    save_figure_formats(fig, OUTPUT, name, dpi=190, bbox_inches='tight')
    plt.close(fig)


def index_map(ax,k):
    ax.scatter(*k.T,s=12,color='#46586b')
    for i,p in enumerate(k):ax.annotate(str(i),p,xytext=(3,3),textcoords='offset points',fontsize=7)
    ax.set(title='Index map',xlim=(-.08,.08),ylim=(-.08,.08),aspect='equal')


def plot(d,axis,contributions,cube,direct,full,summary,direct_map):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    ref=direct[1];k=d['k19']/d['scale'];terms=contributions.sum(axis=3)/ref
    vmax=float(np.max(abs(terms)**2));vmax=10**np.ceil(np.log10(vmax));vmin=vmax*1e-7
    extent=[axis[0],axis[-1],axis[0],axis[-1]]
    for pol,c in enumerate(['Ex','Ey']):
        fig,axes=plt.subplots(4,5,figsize=(14,11),layout='constrained')
        index_map(axes.flat[0],k)
        for i,ax in enumerate(axes.flat[1:]):
            im=ax.imshow(abs(terms[:,:,i,pol])**2,origin='lower',extent=extent,cmap='magma',norm=LogNorm(vmin,vmax))
            ax.set_title(f'{i}'+('  (Gamma)' if i==0 else ''))
            ax.plot(*k[i],marker='+',ms=5,color='#7bbbc0',mew=.8)
            ax.set_aspect('equal')
        for ax in axes.flat:
            ax.set_xticks([-.08,0,.08]);ax.set_yticks([-.08,0,.08])
        for ax in axes[-1]:ax.set_xlabel(r'$K_x/(2\pi/a)$')
        for ax in axes[:,0]:ax.set_ylabel(r'$K_y/(2\pi/a)$')
        bar=fig.colorbar(im,ax=axes,label=r'$|\sum_b A_{nb}T_{\alpha,nb}(K)|^2/|F_y^{air}(0)|^2$',shrink=.92,pad=.02)
        fig.suptitle(f'{c}: 19 radiation contributions | cavity cells only | mesh=5',fontsize=14)
        align_bar(fig,bar,axes)
        save(fig,f'{12+pol:02d}_s5_volume_{c}_19_envelopes')
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    common=max(float(np.max(abs(terms.sum(axis=2))**2)),float(np.max(np.sum(abs(terms)**2,axis=2))))
    common=10**np.ceil(np.log10(common))
    norm=LogNorm(max(common*1e-7,1e-16),common)
    for pol,c in enumerate(['Ex','Ey']):
        coherent=abs(terms[:,:,:,pol].sum(axis=2))**2;incoherent=(abs(terms[:,:,:,pol])**2).sum(axis=2)
        for j,(arr,title) in enumerate([(coherent,'Coherent sum'),(incoherent,'Incoherent reference')]):
            ax=axes[pol,j];im=ax.imshow(arr,origin='lower',extent=extent,norm=norm,cmap='magma')
            ax.set(title=f'{c}: {title}',xlabel=r'$K_x/(2\pi/a)$',ylabel=r'$K_y/(2\pi/a)$')
        ax=axes[pol,2];values=cube.sum(axis=1)[:,pol]/ref
        for i,v in enumerate(values):
            ax.annotate('',(v.real,v.imag),(0,0),arrowprops={'arrowstyle':'->','color':plt.cm.viridis(i/18),'lw':1})
            if abs(v)>0.025*max(np.max(abs(values)),1e-15):ax.text(v.real,v.imag,str(i),fontsize=7)
        total=values.sum();ax.scatter(total.real,total.imag,marker='*',s=90,color='#bf4e45',label='19-point sum')
        ax.scatter((direct[pol]/ref).real,(direct[pol]/ref).imag,marker='x',s=60,color='black',label='Direct air')
        ax.scatter((full[pol]/ref).real,(full[pol]/ref).imag,marker='+',s=70,color='#267588',label='Full volume')
        ax.axhline(0,color='.8',lw=.6);ax.axvline(0,color='.8',lw=.6)
        ax.set(title=f'{c}: complex contributions at K=0',xlabel=r'Re$(F_\alpha/F_y^{air})$',ylabel=r'Im$(F_\alpha/F_y^{air})$',aspect='equal')
        radius=1.18*max(np.max(abs(values)),abs(total),abs(direct[pol]/ref),abs(full[pol]/ref),1e-15)
        ax.set_xlim(-radius,radius);ax.set_ylim(-radius,radius)
        ax.legend(frameon=False,fontsize=7)
    bar=fig.colorbar(im,ax=axes[:,:2],label=r'Intensity / $|F_y^{air}(0)|^2$',shrink=.8,pad=.02)
    fig.suptitle('Cavity-only reference: coherent reconstruction and independent radiation',fontsize=14)
    for ax in axes[:,:2].flat:ax.set_xticks([-.08,-.04,0,.04,.08]);ax.set_yticks([-.08,-.04,0,.04,.08])
    align_bar(fig,bar,axes[:,:2])
    save(fig,'14_s5_volume_coherent_Ex_Ey')
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    labels=['Gamma py','non-Gamma py','other bands','remainder incl. cladding','full volume','direct air']
    vals=[np.asarray(v)[:,0]+1j*np.asarray(v)[:,1] for v in summary['groups'].values()]+[full,direct]
    vals=np.asarray(vals)/ref
    for pol,c in enumerate(['Ex','Ey']):
        ax=axes[pol];x=np.arange(len(labels))
        ax.bar(x-.18,vals[:,pol].real,.36,color='#28618a',label='Real')
        ax.bar(x+.18,vals[:,pol].imag,.36,color='#b44c42',label='Imaginary')
        ax.set_xticks(x,labels,rotation=30,ha='right');ax.axhline(0,color='.5',lw=.6)
        ax.set(title=f'{c}: K=0 complex-amplitude budget',ylabel=r'$F_\alpha/F_y^{air}(0)$')
        ax.legend(frameon=False)
    fig.suptitle(f'mesh=5 | cavity electric L2 residual: {summary["full_vector_L2_residual"]:.1%} | volume/air discrepancy: {summary["volume_vs_air_relative_error"]:.1%}',fontsize=12)
    save(fig,'15_s5_volume_K0_budget')
    prediction=contributions.sum(axis=(2,3))/ref;actual=direct_map/ref
    maps=[prediction,actual,actual-prediction]
    common=10**np.ceil(np.log10(max(np.max(abs(v)**2) for v in maps)))
    fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    for pol,c in enumerate(['Ex','Ey']):
        for j,(value,label) in enumerate(zip(maps,['19-point reconstruction','Direct air spectrum','Complex-field difference'])):
            ax=axes[pol,j]
            im=ax.imshow(abs(value[:,:,pol])**2,origin='lower',extent=extent,cmap='magma',norm=LogNorm(common*1e-7,common))
            ax.set(title=f'{c}: {label}',xlabel=r'$K_x/(2\pi/a)$',ylabel=r'$K_y/(2\pi/a)$')
    bar=fig.colorbar(im,ax=axes,label=r'Intensity / $|F_y^{air}(0)|^2$',shrink=.88,pad=.02)
    fig.suptitle('Cavity-only reference versus direct air | shared intensity scale',fontsize=14)
    for ax in axes.flat:ax.set_xticks([-.08,-.04,0,.04,.08]);ax.set_yticks([-.08,-.04,0,.04,.08])
    align_bar(fig,bar,axes)
    save(fig,'16_s5_volume_vs_direct_air')
    plot_support_control(d,summary)
    plot_nonGamma_py(d,summary)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--order',type=int,default=3)
    run(p.parse_args().order)
