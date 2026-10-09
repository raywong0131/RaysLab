"""Common complex-amplitude analysis; no COMSOL import or independent scaling."""
from __future__ import annotations
from comsol_workflow.output_paths import table_path
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.spatial import cKDTree
from scipy.constants import c
from scripts.analysis.prepare_dipolar_gamma_origin import paths, read, write, RESULT, RUN
from comsol_workflow.dipolar_radiation import plane_farfield
from comsol_workflow.figure_output import save_figure_formats
from scripts.analysis.analyze_dipolar_cladding_scattering import Z0


def save(fig,stem,result_paths=None):
    ps=result_paths or paths();name='s5_gamma_origin_'+(ps['99_config'].name if result_paths else RUN)+'_'+stem
    for ax in fig.axes:ax.tick_params(direction='in')
    fig.canvas.draw();fig.set_layout_engine('none')
    for ax in fig.axes:
        parents=getattr(ax,'_colorbar_info',{}).get('parents',[])
        if parents:
            pos=[x.get_position() for x in parents];box=ax.get_position();lo=min(p.y0 for p in pos);hi=max(p.y1 for p in pos)
            ax.set_axes_locator(None);ax.set_box_aspect(None);ax.set_position([box.x0,lo,box.width,hi-lo]);ax.yaxis.label.set_rotation(270);ax.yaxis.labelpad=16
    save_figure_formats(fig, RESULT, name, pdf_group=ps['11_pdf'].relative_to(RESULT/'11_pdf'),
                        dpi=600, bbox_inches='tight');plt.close(fig)
    return name+'.png'


def positive_zoom_fft2(field,axis_um,q):
    """FFT-based chirp-z transform, positive spatial sign and SI area factor."""
    from scipy.signal import ZoomFFT
    dx=float(axis_um[1]-axis_um[0])*1e-6
    transform=ZoomFFT(len(axis_um),[q[0],q[-1]],m=len(q),fs=2*np.pi/dx,endpoint=True)
    result=np.conj(transform(transform(np.conj(field),axis=-1),axis=0))*dx**2
    return result*np.exp(1j*(q[:,None]+q[None,:])*axis_um[0]*1e-6)


def cap_fft_spectra(xy,values,weights,frequency_Hz):
    # Area-conserving linear deposition onto a common grid; no cell phase fit.
    step=.02
    half=int(np.ceil(np.max(abs(xy))/step))+1
    axis=np.arange(-half,half+1)*step
    uv=(xy-axis[0])/step;ij=np.floor(uv).astype(int);fraction=uv-ij
    q=np.linspace(-2*np.pi*frequency_Hz/c*np.sin(np.deg2rad(10)),2*np.pi*frequency_Hz/c*np.sin(np.deg2rad(10)),201)
    spectra=[]
    for value in values:
        grid=np.zeros((len(axis),len(axis)),complex)
        for ox,oy in [(0,0),(0,1),(1,0),(1,1)]:
            weight=(fraction[:,0] if ox else 1-fraction[:,0])*(fraction[:,1] if oy else 1-fraction[:,1])
            np.add.at(grid,(ij[:,1]+oy,ij[:,0]+ox),value*weights*weight/(step*1e-6)**2)
        spectra.append(positive_zoom_fft2(grid,axis,q))
    spectra=np.array(spectra)
    center=np.sum(values*weights,axis=1)
    check_indices=[(i,j) for i in [0,100,200] for j in [0,100,200]]
    direct=np.array([np.sum(values*(weights*np.exp(1j*(q[j]*xy[:,0]+q[i]*xy[:,1])*1e-6)),axis=1) for i,j in check_indices]).T
    sample=np.array([spectra[:,i,j] for i,j in check_indices]).T
    errors=np.max(abs(sample-direct),axis=1)/np.max(abs(direct),axis=1)
    peaks=np.max(abs(spectra),axis=(1,2));scale=max(float(np.linalg.norm(peaks)),np.finfo(float).tiny)
    center_error=float(np.linalg.norm(spectra[:,100,100]-center))
    resolved=bool(abs(center[0])>1e-10*max(float(peaks[0]),np.finfo(float).tiny))
    checks=dict(grid_step_um=step,grid_shape=[len(axis)]*2,angular_samples=201,angular_limits_deg=[-10,10],spatial_sign='+i q dot r',method='Area-weighted linear deposition followed by separable ZoomFFT (chirp-z); cap window only',direct_quadrature_check_points=9,direct_quadrature_relative_max_errors=errors.tolist(),center_quadrature_relative_error=float(center_error/np.linalg.norm(center)) if np.linalg.norm(center)>1e-10*scale else None,center_quadrature_error_over_spectrum_peak=center_error/scale,center_abs_V_m=abs(center).tolist(),center_over_case_peak=(abs(center)/np.maximum(peaks,np.finfo(float).tiny)).tolist(),center_FULL_resolved=resolved,complex_linearity_relative=float(np.linalg.norm(spectra[0]-spectra[1]-spectra[2])/np.linalg.norm(spectra[0])),center_REST_over_FULL=float(abs(center[2]/center[0])) if resolved else None,no_independent_normalization=True)
    assert np.max(errors)<.01,checks
    assert checks['center_quadrature_error_over_spectrum_peak']<1e-8 and checks['complex_linearity_relative']<1e-10,checks
    return spectra,q,checks,center


def cap_lattice_samples(values,mapping,local_xy,local_weights,cyclic,k):
    """Exact cell-lattice FFT followed by the weighted intracell integral."""
    n=len(mapping)
    np.testing.assert_array_equal(np.sort(cyclic),np.arange(n))
    assert len(k)==n and mapping.shape[1]==len(local_xy)==len(local_weights)
    phase=np.exp(1j*(k@local_xy.T))*local_weights[None,:]*1e-12
    return np.array([np.sum(np.fft.ifft(v[mapping[np.argsort(cyclic)]],axis=0)*n*phase,axis=1) for v in values])


def cell_k_decomposition(cell_fields,local_weights_um2):
    """Unitary cell-label transform; input cells must be in cyclic order."""
    f=np.fft.ifft(cell_fields,axis=0,norm='ortho')
    reconstructed=np.fft.fft(f,axis=0,norm='ortho')
    weights=np.sum(abs(f)**2*local_weights_um2[None,:],axis=1)*1e-12
    norm=float(np.sum(abs(cell_fields)**2*local_weights_um2[None,:])*1e-12)
    checks=dict(reconstruction_relative=float(np.linalg.norm(reconstructed-cell_fields)/np.linalg.norm(cell_fields)),parseval_relative=float(abs(weights.sum()-norm)/norm))
    assert max(checks.values())<1e-12,checks
    return f,weights,reconstructed,checks


def boundary_figure(polarization='Ey',result_paths=None):
    if polarization not in ['Ex','Ey']:raise ValueError('Expected Ex or Ey')
    from scripts.analysis.visualize_dipolar_broadening_2d import SOURCE
    from comsol_workflow.finite_lattice_fourier import direct_lattice_basis,hex_spiral_indices,first_bz_polygon
    ps=result_paths or paths();labels=['FULL','GAMMA','REST'];cfg=read(ps['99_config']/'replay_config.json');channel={'Ex':0,'Ey':1}[polarization]
    with np.load(ps['01_results']/'cap_boundary_inputs.npz') as z:
        assert float(z['cap_z_um'])==cfg['surface']['cap_z_um']
        xy=z['xy_um'];values=np.stack([z[k][:,channel] for k in ['FULL_EH','GAMMA_EH','REST_EH']])
        local_xy=z['local_xy_um'];local_weights=z['local_weights_um2'];mapping=z['cell_to_global_node']
    with np.load(SOURCE) as z:
        basis=np.linalg.lstsq(z['cell_indices'],z['cell_centers'],rcond=None)[0].T
        np.testing.assert_allclose(basis,direct_lattice_basis(cfg['a_um']),atol=1e-12)
        distance,cell_index=cKDTree(z['cell_centers']).query(xy[mapping[:,0]]-local_xy[0])
        assert distance.max()<1e-8 and len(mapping)==cfg['cell_count']
        np.testing.assert_allclose(xy[mapping],z['cell_centers'][cell_index,None]+local_xy[None],atol=1e-8,rtol=0)
        cyclic=z['cyclic_indices'][cell_index];np.testing.assert_array_equal(np.sort(cyclic),np.arange(len(mapping)))
        np.testing.assert_array_equal(hex_spiral_indices(z['k_folded']),z['k_grid_indices'])
        selection=np.argsort(z['k_grid_indices']);k=z['k_folded'][selection]
    cyclic_mapping=mapping[np.argsort(cyclic)];counts=np.bincount(mapping.ravel(),minlength=len(xy))
    coefficients=[];weights=[];rebuilt=[];checks={}
    for label,value in zip(labels,values):
        f,w,cell_rebuilt,check=cell_k_decomposition(value[cyclic_mapping],local_weights)
        flat=cell_rebuilt.ravel();nodes=cyclic_mapping.ravel()
        global_rebuilt=(np.bincount(nodes,weights=flat.real,minlength=len(xy))+1j*np.bincount(nodes,weights=flat.imag,minlength=len(xy)))/counts
        check['global_reconstruction_relative']=float(np.linalg.norm(global_rebuilt-value)/np.linalg.norm(value))
        assert check['global_reconstruction_relative']<1e-12
        coefficients.append(f[selection]);weights.append(w[selection]);rebuilt.append(global_rebuilt);checks[label]=check
    coefficients=np.array(coefficients);weights=np.array(weights);rebuilt=np.array(rebuilt)
    checks['component_linearity_relative']=float(np.linalg.norm(coefficients[0]-coefficients[1]-coefficients[2])/np.linalg.norm(coefficients[0]))
    cross=2*np.real(np.sum(coefficients[1]*coefficients[2].conj()*local_weights[None,:],axis=1))*1e-12
    np.testing.assert_allclose(weights[0],weights[1]+weights[2]+cross,rtol=1e-9,atol=weights[0].sum()*1e-14)
    areaweights=np.bincount(mapping.ravel(),weights=np.tile(local_weights,len(mapping)),minlength=len(xy))*1e-12
    spectra,q,fftchecks,center=cap_fft_spectra(xy,rebuilt,areaweights,cfg['frequency_Hz'])
    reference_path=ps['01_results']/f'cap_{polarization}_fft.npz'
    if not reference_path.exists():
        direct,direct_q,_,direct_center=cap_fft_spectra(xy,values,areaweights,cfg['frequency_Hz'])
        np.savez_compressed(reference_path,**{f'A_{polarization[-1]}_V_m':direct,f'center_A_{polarization[-1]}_V_m':direct_center},q_m_inv=direct_q,cases=labels,frequency_Hz=cfg['frequency_Hz'],cap_z_um=cfg['surface']['cap_z_um'])
    with np.load(reference_path) as z:
        np.testing.assert_array_equal(z['q_m_inv'],q);reference=z[f'A_{polarization[-1]}_V_m']
    checks['reconstructed_spectrum_vs_original_relative']=float(np.linalg.norm(spectra-reference)/np.linalg.norm(reference))
    assert checks['reconstructed_spectrum_vs_original_relative']<1e-10
    denominator=float(weights[0].sum());fractions=weights/denominator
    checks.update(polarization=polarization,point_count=len(k),coefficient_convention='Unitary +i cell-label DFT; F in spiral display order; inverse uses cyclic_to_display',weight_definition='sum_local w*abs(F)**2, area weights in m2; electric L2 weight, not Joules',normalization='W_n / sum_n W_FULL; same denominator for all three cases',FULL_total_weight_V2=denominator,weight_sums_over_FULL=fractions.sum(axis=1).tolist(),Gamma_index0_fraction_of_own_weight=float(weights[1,0]/weights[1].sum()),coherent_reconstruction='Sum all 1141 complex cell-label components via inverse DFT, then cap Fourier transform; no per-component intensity sum',FFT=fftchecks)
    write(ps['80_logs']/f'cap_{polarization}_k_decomposition_checks.json',checks)
    np.savez_compressed(ps['01_results']/f'cap_{polarization}_k_decomposition.npz',F_V_m=coefficients,W_V2=weights,P_over_FULL=fractions,cross_weight_V2=cross,k_rad_um=k,indices=np.arange(len(k)),cyclic_to_display=np.argsort(selection),local_xy_um=local_xy,local_weights_um2=local_weights,cases=labels,frequency_Hz=cfg['frequency_Hz'])
    np.savez_compressed(ps['01_results']/f'cap_{polarization}_k_reconstruction.npz',A_V_m=spectra,q_m_inv=q,center_A_V_m=center,cases=labels,frequency_Hz=cfg['frequency_Hz'],cap_z_um=cfg['surface']['cap_z_um'])
    with (table_path(ps['80_logs'], f'cap_{polarization}_k_weights.csv')).open('w',newline='',encoding='utf-8') as file:
        writer=csv.writer(file);writer.writerow(['case','index','kx_rad_um','ky_rad_um','W_V2','W_over_FULL_total'])
        for i,label in enumerate(labels):
            for n,point in enumerate(k):writer.writerow([label,n,*point,weights[i,n],fractions[i,n]])
    axis=np.linspace(xy.min(),xy.max(),1501);xx,yy=np.meshgrid(axis,axis);distance,index=cKDTree(xy).query(np.column_stack([xx.ravel(),yy.ravel()]))
    maps=[np.ma.masked_where(distance.reshape(xx.shape)>.025,v.real[index].reshape(xx.shape)) for v in values]
    kscale=2*np.pi/cfg['a_um'];points=k/kscale;bz=first_bz_polygon(cfg['a_um'])/kscale;bz=np.vstack([bz,bz[0]])
    intensity=abs(spectra)**2;scales=(max(float(abs(v).max()) for v in maps),denominator,float(intensity.max()))
    return plot_cap_summary(axis,maps,q*1e-6/kscale,points[:19],points,weights,intensity,bz,scales,polarization,result_paths)


def plot_cap_summary(axis,maps,kaxis,points,all_points,all_intensity,intensity,bz,scales,polarization,result_paths=None):
    labels=['FULL','GAMMA','REST'];cap_limit,all_scale,fft_limit=scales
    combined=polarization=='ExEy';suffix=polarization[-1]
    field_label=r'$|E_{xy}|$' if combined else f'Re({polarization})'
    spectrum_label=r'$|A_x|^2+|A_y|^2$' if combined else rf'$|A_{suffix}|^2$'
    fig,axs=plt.subplots(3,3,figsize=(13,11.5),layout='constrained')
    if result_paths:fig.suptitle('z = '+str(read(result_paths['99_config']/'replay_config.json')['surface']['cap_z_um'])+' um')
    for i,label in enumerate(labels):
        im=axs[0,i].pcolormesh(axis,axis,maps[i]/(cap_limit or 1.),cmap='magma' if combined else 'RdBu_r',vmin=0 if combined else -1,vmax=1,shading='auto',rasterized=True);axs[0,i].set(aspect='equal',title=label+' | cap '+field_label,xlabel='x (um)',ylabel='y (um)')
        if combined:axs[0,i].set_facecolor(im.cmap(im.norm(0)))
        axs[1,i].set_facecolor('black')
        axs[1,i].plot(*bz.T,color='white',linewidth=.9,zorder=2)
        dots=axs[1,i].scatter(*all_points.T,c=all_intensity[i]/(all_scale or 1.),s=7,cmap='magma',norm=LogNorm(1e-4,1,clip=True),linewidths=0,zorder=3)
        axs[1,i].set(aspect='equal',xlim=(-.73,.73),ylim=(-.635,.635),xlabel=r'$K_x/(2\pi/a)$',ylabel=r'$K_y/(2\pi/a)$',title=label+f' | {len(all_points)} k-component weights')
        axs[1,i].set_xticks([-.6,-.3,0,.3,.6]);axs[1,i].set_yticks([-.6,-.3,0,.3,.6])
        fft_im=axs[2,i].pcolormesh(kaxis,kaxis,intensity[i]/(fft_limit or 1.),cmap='magma',vmin=0,vmax=1,shading='auto',rasterized=True)
        axs[2,i].scatter(*points.T,s=27,facecolors='none',edgecolors='#b7c1cf',linewidths=.55)
        axs[2,i].set_title(label+(' | cap angular spectrum' if combined else f' | cap FFT({polarization})'))
        for ax in axs[2:,i]:
            ax.set(aspect='equal',xlim=(kaxis[0],kaxis[-1]),ylim=(kaxis[0],kaxis[-1]),xlabel=r'$K_x/(2\pi/a)$',ylabel=r'$K_y/(2\pi/a)$')
            ax.set_xticks([-.08,-.04,0,.04,.08]);ax.set_yticks([-.08,-.04,0,.04,.08])
    fig.colorbar(im,ax=axs[0].tolist(),ticks=np.linspace(0,1,6) if combined else np.linspace(-1,1,9),format='%.2f',label=field_label+' normalized (1)')
    fig.colorbar(dots,ax=axs[1].tolist(),ticks=np.logspace(-4,0,5),format=LogFormatterMathtext(),label=r'$W_n/\sum_m W_{\mathrm{FULL},m}$ (1)')
    fig.colorbar(fft_im,ax=axs[2].tolist(),ticks=np.linspace(0,1,6),format='%.2f',label=spectrum_label+' normalized (1)')
    stem='01_boundary_inputs' if polarization=='Ey' else '01_boundary_inputs_'+polarization
    return save(fig,stem,result_paths) if result_paths else save(fig,stem)


def combined_boundary_figure(result_paths=None):
    """Combine saved physical Ex/Ey amplitudes before any display normalization."""
    from comsol_workflow.finite_lattice_fourier import first_bz_polygon
    ps=result_paths or paths();cfg=read(ps['99_config']/'replay_config.json');spectra=[];samples=[]
    for pol in ['Ex','Ey']:
        with np.load(ps['01_results']/f'cap_{pol}_k_reconstruction.npz') as z:
            if spectra:
                np.testing.assert_array_equal(z['q_m_inv'],q)
                np.testing.assert_array_equal(z['cases'],cases)
            q=z['q_m_inv'];cases=z['cases'];spectra.append(z['A_V_m'])
            assert float(z['frequency_Hz'])==cfg['frequency_Hz'] and float(z['cap_z_um'])==cfg['surface']['cap_z_um']
        with np.load(ps['01_results']/f'cap_{pol}_k_decomposition.npz') as z:
            points=z['k_rad_um']/(2*np.pi/cfg['a_um'])
            if samples:np.testing.assert_array_equal(points,all_points)
            np.testing.assert_array_equal(z['cases'],cases)
            np.testing.assert_array_equal(z['indices'],np.arange(cfg['cell_count']))
            all_points=points;samples.append(z['W_V2'])
    spectra=np.stack(spectra,axis=-1);samples=np.stack(samples,axis=-1)
    intensity=np.sum(abs(spectra)**2,axis=-1);all_intensity=np.sum(samples,axis=-1)
    with np.load(ps['01_results']/'cap_boundary_inputs.npz') as z:
        xy=z['xy_um'];axis=np.linspace(xy.min(),xy.max(),1501);xx,yy=np.meshgrid(axis,axis)
        distance,index=cKDTree(xy).query(np.column_stack([xx.ravel(),yy.ravel()]))
        maps=[np.ma.masked_where(distance.reshape(xx.shape)>.025,np.linalg.norm(z[k][index,:2],axis=-1).reshape(xx.shape)) for k in ['FULL_EH','GAMMA_EH','REST_EH']]
    scales=(max(float(v.max()) for v in maps),float(all_intensity[0].sum()),float(intensity.max()))
    kscale=2*np.pi/cfg['a_um'];kaxis=q*1e-6/kscale;bz=first_bz_polygon(cfg['a_um'])/kscale;bz=np.vstack([bz,bz[0]])
    checks=dict(components=['Ex','Ey'],includes_Ez=False,field_definition='sqrt(abs(Ex)**2+abs(Ey)**2)',spectrum_definition='abs(Ax)**2+abs(Ay)**2; transform complex components before taking intensity',normalization='Row2: Wxy_n/sum Wxy_FULL; rows1/3 shared peak per row; combine physical components first',field_scale_V_m=scales[0],FULL_total_weight_V2=scales[1],central_scale_V2_m2=scales[2],full_BZ_points=len(all_points),continuous_shape=list(intensity.shape[1:]),complex_linearity_relative=float(np.linalg.norm(spectra[0]-spectra[1]-spectra[2])/np.linalg.norm(spectra[0])),weight_sums_over_FULL=(all_intensity.sum(axis=1)/scales[1]).tolist(),cap_window_only=True)
    assert checks['complex_linearity_relative']<1e-10
    write(ps['80_logs']/'cap_ExEy_combined_checks.json',checks)
    np.savez_compressed(ps['01_results']/'cap_ExEy_combined.npz',cases=cases,components=['Ex','Ey'],A_V_m=spectra,component_W_by_polarization_V2=samples,intensity_V2_m2=intensity,component_W_V2=all_intensity,q_m_inv=q,k_normalized=all_points,cap_axis_um=axis,cap_magnitude_V_m=np.array([v.filled(np.nan) for v in maps]),display_scales=scales,frequency_Hz=cfg['frequency_Hz'],cap_z_um=cfg['surface']['cap_z_um'])
    return plot_cap_summary(axis,maps,kaxis,all_points[:19],all_points,all_intensity,intensity,bz,scales,'ExEy',result_paths)


def report(cfg):
    ps=paths();gate=read(ps['80_logs']/'FULL_replay_check.json')
    if not gate['passed']:raise ValueError('Cannot attribute components before FULL replay')
    labels=['FULL','GAMMA','REST'];k0=2*np.pi*cfg['frequency_Hz']/c;q=np.linspace(-k0*np.sin(np.deg2rad(10)),k0*np.sin(np.deg2rad(10)),81);theta=np.rad2deg(np.arcsin(q/k0));fields=[];air=[];slab=[]
    for label in labels:
        with np.load(ps['01_results']/f'{label}_air_EH.npz') as z:d={k:z[k] for k in z.files}
        fields.append(plane_farfield(d['E_V_m'],d['H_A_m'],d['x_um'],float(d['z_um']),q,k0));air.append(d)
        with np.load(ps['01_results']/f'{label}_exterior_slab_EH.npz') as z:slab.append({k:z[k] for k in z.files})
    f=np.array(fields);delta=f[0]-f[1]-f[2]
    linearity=float(np.linalg.norm(delta)/np.linalg.norm(f[0]));center=float(np.linalg.norm(delta[40,40])/np.linalg.norm(f[0,40,40]))
    r=float(abs(f[2,40,40,1])/abs(f[0,40,40,1]))
    fullI=np.sum(abs(f[0])**2,-1)/(2*Z0)
    # Fixed connected half-height region of the original A central lobe.
    with np.load(ps['01_results']/'FULL_reference_farfield.npz') as z:original=z['reference_A']
    referenceI=np.sum(abs(original)**2,-1)/(2*Z0)
    from scipy.ndimage import label as connected
    regions,_=connected(referenceI>=referenceI[40,40]/2);region=regions==regions[40,40]
    assert regions[40,40]!=0
    local_rest=float(np.linalg.norm(f[2][region])/np.linalg.norm(f[0][region]))
    checks=dict(output_linearity_relative=linearity,center_linearity_relative=center,REST_y_center_amplitude_ratio=r,REST_central_lobe_amplitude_norm_ratio=local_rest,central_lobe_definition='Connected original A intensity >= half normal intensity on fixed 81x81 grid',central_lobe_points=int(region.sum()),FULL_replay=gate,no_phase_alignment=True,no_independent_normalization=True)
    checks['verdict']='implementation-unresolved' if max(linearity,center)>1e-6 else ('Gamma-dominant with finite REST; report ratios, not strict exclusivity' if r<.1 else 'Gamma alone does not explain the center')
    write(ps['80_logs']/'linearity_checks.json',checks)
    np.savez_compressed(ps['01_results']/'farfield_components.npz',F_V=f,q_m_inv=q,cases=labels,central_lobe_mask=region,frequency_Hz=cfg['frequency_Hz'])
    rows=[]
    for index,case in enumerate(labels):
        for pol,channel in enumerate(['Ex','Ey']):
            v=f[index,40,40,pol];rows.append(dict(case=case,channel=channel,Re_V=float(v.real),Im_V=float(v.imag),abs_V=float(abs(v)),phase_rad=float(np.angle(v)),self_W_sr=float(abs(v)**2/(2*Z0))))
    cross=float(np.real(np.vdot(f[2,40,40],f[1,40,40]))/Z0);checks['center_interference_W_sr']=cross;write(ps['80_logs']/'linearity_checks.json',checks)
    with (table_path(ps['80_logs'], 'normal_amplitudes.csv')).open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    figures=[boundary_figure()]
    fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    replay_fields=[original,f[0],f[0]-original]
    for pol in range(2):
        limit=max(float((abs(v[...,pol])**2/(2*Z0)).max()) for v in replay_fields)
        for i,name in enumerate(['Original A','FULL replay','FULL minus A']):
            im=axs[pol,i].pcolormesh(theta,theta,abs(replay_fields[i][...,pol])**2/(2*Z0),vmin=0,vmax=limit,cmap='magma',shading='auto',rasterized=True)
            axs[pol,i].set(aspect='equal',title=name+' | E'+'xy'[pol],xlabel='theta_x (deg)',ylabel='theta_y (deg)')
        fig.colorbar(im,ax=axs[pol].tolist(),label='dP/dOmega (W/sr)')
    figures.insert(0,save(fig,'00_replay_validation'))
    intensity=abs(f[...,:2])**2/(2*Z0);fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    for pol in range(2):
        limit=float(intensity[...,pol].max())
        for i,name in enumerate(labels):
            im=axs[pol,i].pcolormesh(theta,theta,intensity[i,...,pol],vmin=0,vmax=limit,cmap='magma',shading='auto',rasterized=True);axs[pol,i].set(aspect='equal',title=name+' | E'+'xy'[pol],xlabel='theta_x (deg)',ylabel='theta_y (deg)')
        fig.colorbar(im,ax=axs[pol].tolist(),label='dP/dOmega (W/sr)')
    figures.append(save(fig,'02_common_farfields'))
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for pol in range(2):
        vals=f[:,40,40,pol];bound=max(abs(vals.real).max(),abs(vals.imag).max(),1e-20)*1.15
        for i,v in enumerate(vals):axs[pol].quiver(0,0,v.real,v.imag,angles='xy',scale_units='xy',scale=1,color=['black','#c43c39','#2c72b0'][i]);axs[pol].text(v.real,v.imag,labels[i])
        axs[pol].plot([vals[1].real,(vals[1]+vals[2]).real],[vals[1].imag,(vals[1]+vals[2]).imag],color='#2c72b0')
        axs[pol].set(xlim=(-bound,bound),ylim=(-bound,bound),aspect='equal',title='Normal E'+'xy'[pol],xlabel='Re(F) (V)',ylabel='Im(F) (V)')
    figures.append(save(fig,'03_center_phasors'))
    fig,axs=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for i,name in enumerate(labels):axs[0,0].plot(theta,intensity[i,40,:,1],label=name);axs[0,1].plot(theta,intensity[i,:,40,1],label=name)
    for i,direction in enumerate('xy'):axs[0,i].set(xlabel='theta_'+direction+' (deg)',ylabel='Ey dP/dOmega (W/sr)');axs[0,i].legend()
    axs[0,2].axis('off');axs[0,2].text(0,.8,f'REST center amplitude ratio: {r:.6g}\nREST central-lobe ratio: {local_rest:.6g}\nComplex linearity error: {linearity:.3g}')
    bound=max(float(abs(d['E_V_m'][...,1].real).max()) for d in slab)
    for i,d in enumerate(slab):
        im=axs[1,i].pcolormesh(d['x_um'],d['y_um'],np.ma.masked_where(~d['valid'],d['E_V_m'][...,1].real),vmin=-bound,vmax=bound,cmap='RdBu_r',shading='auto',rasterized=True);axs[1,i].set(aspect='equal',title=labels[i]+' | exterior Re(Ey)',xlabel='x (um)',ylabel='y (um)')
    fig.colorbar(im,ax=axs[1].tolist(),label='Re(Ey) (V/m)');figures.append(save(fig,'04_cuts_and_exterior'))
    text='# S5 全局 Γ 分量的辐射验证：综合结果\n\n'
    text+='本次以实际有限结构的 1141 个完整 cavity cell 复平均定义 Γ；REST 为完整余项。保持实际 cladding、mesh5、频率和边界映射一致，不逐 cell 调相，不对三组独立缩放。\n\n'
    text+='## 验证结果\n\n'
    text+='| FULL 重放指标 | 相对误差 |\\n| --- | ---: |\\n'.replace('\\n','\n')
    for key in ['E_V_m','H_A_m','center_complex_relative_error','farfield_complex_relative_error']:
        text+=f'| {key} | {gate[key]:.8g} |\n'
    text+=f'\n预先设置的 FULL 门槛为 {gate["pilot_limit"]:.0%}；实际误差见上表。\n\n'
    text+=f'- 远场复振幅线性闭合误差：{linearity:.8g}；中心误差：{center:.8g}。\n- REST/FULL 中心 Ey 振幅比：{r:.8g}。\n- REST/FULL 中心亮斑区域振幅范数比：{local_rest:.8g}。\n- 中心 Γ/REST 相干交叉项：{cross:.8g} W/sr。\n\n'
    text+='| 组别 | 分量 | 实部 (V) | 虚部 (V) | 模 (V) | 相位 (rad) |\n| --- | --- | ---: | ---: | ---: | ---: |\n'
    for row in rows:
        text+=f'| {row["case"]} | {row["channel"]} | {row["Re_V"]:.8g} | {row["Im_V"]:.8g} | {row["abs_V"]:.8g} | {row["phase_rad"]:.8g} |\n'
    text+='\n若分量接近数值零，其相位没有稳定物理意义。远场图每种偏振内三组共用线性色标；Ex 与 Ey 分行标注范围。边界图上排为 cap 实部、下排为 cap FFT 角谱，分别标注范围；cap FFT 是窗口角谱，不替代完整外域辐射验证，原始数据未独立归一化。\n\n'

    if max(linearity,center)>1e-6:
        text+='**实现验证未通过：输出线性关系不满足阈值，不作机制归因。**\n\n'
    elif r<.1:
        text+='本次离散边界分解下，Γ 对中心 Ey 振幅占主导，REST 仍具有上述有限贡献；不能表述为严格全部来自 Γ。中心亮斑邻域须结合区域比值单独判断。\n\n'
    else:
        text+='Γ 单独不足以解释中心 Ey 复振幅，REST 在中心具有实质有限贡献。不得据此强行支持预设的 Γ 独占解释。\n\n'
    text+='这些是振幅比，不是可加功率份额。强度包含 2Re(FΓ·FREST*)。本次 Γ 是有限 cell 标签 index-0，不等同于无限周期 Γ 本征态；本次 REST 是完整余项，不局限于历史的 19 个 k 点。A 为原 mode19 的实频驱动代表场。\n\n'
    text+='## 边界表示与适用范围\n\n原始 Γ 与共用节点映射 Qh 保留在 gamma_definition.json；边界检查保留在 boundary_input_checks.json。映射与有限采样引入离散误差，FULL 重放是独立门控。本实验不单独证明 cladding 的结构因果作用，也不证明严格的普适消去定理。\n\n'
    text+='## 图与原始结果\n\n'
    text+='\n\n'.join('!['+name+'](../../../10_overview/'+name+')' for name in figures)+'\n'
    (ps['12_reports']/'gamma_origin_report.md').write_text(text,encoding='utf-8')
    return checks


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);args=p.parse_args();report(read(args.config))
if __name__=='__main__':main()
