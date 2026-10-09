"""S5 Bloch/radiation validation using saved fields, without solver imports.

Default: exact-k Hz projection, periodic-Gamma controls, and independent finite
air-plane integrals. An optional --radiation-bundle supplies independently
computed full-vector A, T, and direct radiation for the actual closure test.
"""
from __future__ import annotations

from comsol_workflow.dipolar_radiation import weighted_projection, radiation_budget, volume_far_amplitude, planar_K0_to_far, pairs

from comsol_workflow.output_paths import output_path
from comsol_workflow.figure_output import save_figure_formats

import argparse
import json
from pathlib import Path

import numpy as np

from comsol_workflow.farfield_fft import load_air_field
from comsol_workflow.dipolar_lattice_analysis import load_data
from scripts.analysis.visualize_dipolar_broadening_2d import SOURCE, OUTPUT
from comsol_workflow.result_io import digest, read_json as json_read, write_json, write_csv
from scripts.run_main.run_dipolar_periodic_reference import OUT as REFERENCES












def finite_air_check(source):
    f = load_air_field(source)
    x,y = f['x_axis'],f['y_axis']
    dx,dy = float(np.mean(np.diff(x))),float(np.mean(np.diff(y)))
    np.testing.assert_allclose(np.diff(x),dx,atol=1e-12)
    np.testing.assert_allclose(np.diff(y),dy,atol=1e-12)
    e = np.stack([f['Ex'],f['Ey']],axis=-1)
    if not np.isfinite(e).all():
        raise ValueError('Invalid saved air field')
    total = e.sum(axis=(0,1))*dx*dy*1e-12
    reference = float(np.linalg.norm(total))
    if reference == 0:
        raise ValueError('Zero direct amplitude; a nonzero error reference is required')
    rows=[]
    for fraction in [.25,.5,.75,.9,1.]:
        mask=(abs(y[:,None])<=fraction*max(abs(y)))*(abs(x[None,:])<=fraction*max(abs(x)))
        c=(e*mask[:,:,None]).sum(axis=(0,1))*dx*dy*1e-12
        rows.append({'fraction_of_saved_box':fraction,'Cx_V_m':c[0],'Cy_V_m':c[1],
                     'difference_from_saved_full_over_reference':float(np.linalg.norm(c-total)/reference)})
    return total,rows


def consume_radiation_bundle(path,expected_k):
    """Accept real full-vector projection/response data, never inferred M*c."""
    with np.load(path,allow_pickle=False) as data:
        meta=json.loads(str(data['metadata_json']))
        required=['independent_response','full_vector_projection','same_observation_and_units']
        if not all(meta.get(key) is True for key in required) or not meta.get('response_provenance'):
            raise ValueError('Missing independently validated full-vector radiation provenance')
        np.testing.assert_allclose(data['k_um_inv'],expected_k,rtol=0,atol=1e-12)
        a,t,direct=data['A'],data['T'],data['C_direct']
        if a.shape!=(19,4) or direct.shape!=(2,) or not np.isfinite(direct).all() or np.linalg.norm(direct)==0:
            raise ValueError('Expected 19 x 4 coefficients and a finite nonzero two-component direct field')
        terms,total,inc,coh,cross=radiation_budget(a,t)
        error=float(np.linalg.norm(total-direct)/np.linalg.norm(direct))
        # Band axis follows the tracked Gamma quartet; index 1 is user py.
        cube=terms.reshape(19,4,2)
        groups={'Gamma_py':cube[0,1], 'nonGamma_py':cube[1:,1].sum(axis=0),
                'other_bands':np.delete(cube,1,axis=1).sum(axis=(0,1))}
        return {'metadata':meta,'terms':pairs(cube),'groups':{k:pairs(v) for k,v in groups.items()},
                'total':pairs(total),'direct':pairs(direct),'residual_direct_minus_19':pairs(direct-total),
                'relative_amplitude_closure_error':error,'incoherent':inc.tolist(),
                'coherent':coh.tolist(),'cross':cross.tolist(),
                'scope':'19-point full-vector test; residual includes omitted bands, k, regions and numerical error'}


def consume_volume_bundle(path, expected_k):
    """Construct T from 3D electric fields and material contrast, not far-field fitting."""
    with np.load(path,allow_pickle=False) as data:
        meta=json.loads(str(data['metadata_json']))
        if meta.get('time_convention')!='native_COMSOL_exp(+iwt)' or not meta.get('complete_material_volume'):
            raise ValueError('Native phasors and complete material-volume provenance required')
        if not meta.get('finite_field_source') or not meta.get('basis_source'):
            raise ValueError('Missing independently sourced finite field or Bloch dictionary')
        np.testing.assert_allclose(data['k_um_inv'],expected_k,rtol=0,atol=1e-12)
        basis,field=data['E_basis'],data['E_finite']
        xyz,weights,chi=data['xyz_m'],data['weights_m3'],data['epsilon_contrast']
        omega=complex(data['omega_rad_s'])
        if basis.shape!=(19,4,len(xyz),3) or field.shape!=(len(xyz),3):
            raise ValueError('Expected 19 x 4 Bloch vector fields on the finite volume nodes')
        metric=np.repeat(weights*(1+chi.real),3)
        # ponytail: dense fit for moderate dictionaries; use block Gram accumulation for large volume exports.
        a,fit,error,condition=weighted_projection(basis.reshape(76,-1).T,field.ravel(),metric)
        a=a.reshape(19,4)
        response=volume_far_amplitude(basis,chi,xyz,weights,omega)
        terms,total,inc,coh,cross=radiation_budget(a,response)
        complete=volume_far_amplitude(field,chi,xyz,weights,omega)
        residual=volume_far_amplitude(field-fit.reshape(field.shape),chi,xyz,weights,omega)
        direct=planar_K0_to_far(data['direct_air_amplitude_V_m'],float(data['z_air_m']),omega)
        scale=float(np.linalg.norm(direct))
        if direct.shape!=(2,) or not np.isfinite(direct).all() or scale==0:
            raise ValueError('Finite nonzero independent direct air amplitude required')
        np.testing.assert_allclose(total+residual,complete,atol=max(scale,np.linalg.norm(complete))*1e-10,rtol=1e-10)
        cube=terms.reshape(19,4,2)
        return {'units':'V; spherical outgoing amplitude F', 'metadata':meta,
                'full_vector_reconstruction_error':error,'normalized_Gram_condition':condition,
                'volume_vs_direct_air_relative_error':float(np.linalg.norm(complete-direct)/scale),
                'basis_vs_direct_air_relative_error':float(np.linalg.norm(total-direct)/scale),
                'groups':{'Gamma_py':pairs(cube[0,1]),'nonGamma_py':pairs(cube[1:,1].sum(axis=0)),
                          'other_bands':pairs(np.delete(cube,1,axis=1).sum(axis=(0,1))),
                          'field_residual':pairs(residual)},
                'terms':pairs(cube),'total_basis':pairs(total),'complete_volume':pairs(complete),
                'direct_air':pairs(direct),'incoherent':inc.tolist(),'coherent':coh.tolist(),'cross':cross.tolist()}


def plot_checks(d,rows,summary,apertures,k_union):
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    ax=axes[0,0]; k=d['k19']/d['scale'];extra=np.asarray(k_union)[19:]/d['scale']
    ax.scatter(*k.T,s=23,color='#28618a',label='Original 19 points')
    ax.scatter(*extra.T,s=22,marker='x',color='#999999',label='Mirror-closure controls (not solved)')
    for i,point in enumerate(k):
        ax.annotate(str(i),point,xytext=(3,3),textcoords='offset points',fontsize=8)
    ax.set(title='a  Index map and mirror-closure check',xlabel=r'$k_x/(2\pi/a)$',ylabel=r'$k_y/(2\pi/a)$',
           xlim=(-.08,.08),ylim=(-.08,.08),aspect='equal')
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    ax=axes[0,1]; x=np.arange(19)
    ax.plot(x,[100*r['single_py_relative_residual'] for r in rows],'.-',color='#b44c42',label='Single py band')
    ax.plot(x,[100*r['four_band_relative_residual'] for r in rows],'.-',color='#28618a',label='Four-band dictionary')
    ax.set(title='b  Internal Hz reconstruction at each k',xlabel='Original k index',ylabel='Relative L2 residual (%)',
           xticks=[0,3,6,9,12,15,18],ylim=(0,55))
    ax.legend(frameon=False)
    ax=axes[1,0]; g=summary['Gamma_controls'][0]
    raw=np.asarray(g['common_phase_c_by_quadrature_height']);c=raw[...,0]+1j*raw[...,1]
    x=np.arange(4)
    ax.semilogy(x,abs(c[:,:,0]).ravel(),'o-',color='#28618a',label='|cx| / Hz RMS')
    ax.semilogy(x,abs(c[:,:,1]).ravel(),'s-',color='#b44c42',label='|cy| / Hz RMS')
    ax.set(title=r'c  Periodic Gamma radiation | mesh=5',
           ylabel=r'Electric / magnetic amplitude ($\Omega$)',
           xticks=x,xticklabels=['Q1 low z','Q1 high z','Q2 low z','Q2 high z'])
    ax.legend(frameon=False,fontsize=8)
    ax=axes[1,1]; cy_full=complex(*summary['finite_direct_full_saved_air_V_m'][1])
    values=np.array([r['Cy_V_m'] for r in apertures])/cy_full
    ax.plot([r['fraction_of_saved_box'] for r in apertures],values.real,'o-',color='#28618a',label='Re(Cy / Cy full)')
    ax.plot([r['fraction_of_saved_box'] for r in apertures],values.imag,'s-',color='#b44c42',label='Im(Cy / Cy full)')
    cavity=complex(*summary['finite_direct_cavity_air_V_m'][1])
    ax.text(.04,.91,f'Cavity footprint: |Cy| / |Cy full| = {abs(cavity/cy_full):.4f}',transform=ax.transAxes,fontsize=9)
    ax.set(title='d  Finite air-plane aperture sensitivity',xlabel='Fraction of saved box side length',
           ylabel='Complex amplitude / common full-plane reference',ylim=(-.2,1.2))
    ax.legend(frameon=False,fontsize=8,loc='lower right')
    fig.suptitle('S5 verification | exact periodic references and finite-mode checks',fontsize=15)
    fig.supxlabel('Hz projection is a diagnostic. Fixed mesh=5. Independent electric-volume radiation results are shown in figures 12-15.',fontsize=9)
    save_figure_formats(fig, OUTPUT, '11_s5_validation_checks', dpi=180, bbox_inches='tight')
    plt.close(fig)


def run(response_path=None,volume_path=None):
    d=load_data(SOURCE)
    manifest_path=REFERENCES/'99_config/manifest.json'
    manifest=json_read(manifest_path)
    np.testing.assert_allclose(manifest['k_um_inv'],d['k19'],rtol=0,atol=1e-12)
    if digest(SOURCE)!=manifest['source_Hz_sha256']:
        raise ValueError('Finite source changed since reference preparation')
    provenance={str(SOURCE):digest(SOURCE),str(manifest_path):digest(manifest_path)}
    weights=np.full(len(d['rho_points']),1/len(d['rho_points']))
    coefficients=[]; references=[]; rows=[]; gauge_error=0.
    rng=np.random.default_rng(217)
    for n,k in enumerate(d['k19']):
        path=REFERENCES/f'01_results/k{n:02d}_reference.npz'
        audit=json_read(REFERENCES/f'01_results/k{n:02d}_audit.json')
        if digest(path)!=audit['data_sha256']:
            raise ValueError(f'Reference hash mismatch at index {n}')
        provenance[str(path)]=digest(path)
        with np.load(path,allow_pickle=False) as saved:
            q={key:saved[key] for key in saved.files}
        np.testing.assert_allclose(q['rho_points_um'],d['rho_points'],rtol=0,atol=1e-12)
        np.testing.assert_allclose(q['k_um_inv'],k,rtol=0,atol=1e-12)
        order=q['gamma_band_to_mode']
        u=q['u_Hz'][order].T
        field=d['F19'][n]*np.exp(1j*(d['rho_points']@k))
        a,_,residual,condition=weighted_projection(u,field,weights)
        single,_,single_residual,_=weighted_projection(u[:,1:2],field,weights)
        radiation=q['air_mean_E_at_z0_over_Hz_rms'][1,1,order]
        phase=np.exp(1j*rng.uniform(-np.pi,np.pi,4))
        rotated,_,_,_=weighted_projection(u*phase,field,weights)
        error=np.linalg.norm(rotated[:,None]*radiation*phase[:,None]-a[:,None]*radiation)
        scale=max(np.linalg.norm(a[:,None]*radiation),1e-300)
        gauge_error=max(gauge_error,float(error/scale))
        coefficients.append(a);references.append(q)
        rows.append({'index':n,'layer':1 if n==0 else 2 if n<7 else 3,
                     'kx_um_inv':float(k[0]),'ky_um_inv':float(k[1]),
                     'Hz_weight_in_all_k':float(d['P19'][n]),'single_py_relative_residual':single_residual,
                     'four_band_relative_residual':residual,'normalized_Gram_condition':condition,
                     'py_mode_idx':int(order[1]),'py_frequency_thz':float(q['frequencies_THz'][order[1],0]),
                     'py_Q':float(q['frequencies_THz'][order[1],2]),
                     'A_py':a[1],'single_py_A':single[0],
                     'periodic_cx_over_Hz_rms':radiation[1,0],'periodic_cy_over_Hz_rms':radiation[1,1]})
    if gauge_error>1e-10:
        raise ValueError('Gauge invariance failed')
    a=np.asarray(coefficients)
    # This explicitly tests the old hard-aperture/reference-channel hypothesis.
    # It is NOT the independently computed finite-structure response T.
    electric_meta=json_read(output_path(OUTPUT, 'electric_19_summary.json'))
    air=Path(electric_meta['source'])
    provenance[str(air)]=digest(air)
    if provenance[str(air)]!=electric_meta['sha256']:
        raise ValueError('Finite air source changed')
    direct,apertures=finite_air_check(air)
    cavity=np.array([complex(*electric_meta['C0_cavity_footprint_V_m'][c]) for c in ['Ex','Ey']])
    s0=np.exp(-1j*(d['cell_centers']@d['k19'].T)).mean(axis=0)
    zero_orthogonality=float(np.max(abs(s0-np.eye(1,19)[0])))
    channel_response=[]
    for n,q in enumerate(references):
        order=q['gamma_band_to_mode']
        f=q['frequencies_THz'][order]
        kz=np.sqrt((2*np.pi*(f[:,0]+1j*f[:,1])*1e12/299792458.)**2-np.dot(d['k19'][n],d['k19'][n])*1e12)
        c=q['air_mean_E_at_z0_over_Hz_rms'][1,1,order]
        c=c*np.exp(-1j*kz[:,None]*electric_meta['plane_z_um']*1e-6)
        channel_response.append(np.sqrt(d['N'])*d['area']*1e-12*s0[n]*c)
    terms,total,inc,coh,cross=radiation_budget(a,np.asarray(channel_response))
    controls=[]
    for stem in ['k00']:
        path=REFERENCES/f'01_results/{stem}_reference.npz'
        if not path.exists():
            continue
        provenance[str(path)]=digest(path)
        with np.load(path,allow_pickle=False) as q:
            target=int(q['gamma_band_to_mode'][1])
            u=q['u_Hz'][target]
            factor=np.exp(-1j*np.angle(np.vdot(references[0]['u_Hz'][1],u)))
            c=q['air_mean_E_at_z0_over_Hz_rms'][:,:,target]*factor
            qvalue=float(q['frequencies_THz'][target,2])
        controls.append({'reference':stem,'Q':qvalue,'common_phase_c_by_quadrature_height':pairs(c),
                         'cy_height_difference_abs':float(abs(c[1,1,1]-c[1,0,1])),
                         'cy_quadrature_difference_at_upper_height_abs':float(abs(c[1,1,1]-c[0,1,1]))})
    mirror=[]; k_union=d['k19'].tolist()
    for reflect in [np.array([-1,1]),np.array([1,-1])]:
        reflected=d['k19']*reflect
        distances=np.linalg.norm(reflected[:,None]-d['k19'][None,:],axis=-1).min(axis=1)
        mirror.append(float(distances.max()/d['scale']))
        for k in reflected:
            if np.min(np.linalg.norm(np.asarray(k_union)-k,axis=1))>1e-10:
                k_union.append(k.tolist())
    write_csv(output_path(OUTPUT, 's5_exact_bloch_projection.csv'),rows)
    write_csv(output_path(OUTPUT, 's5_air_aperture_check.csv'),apertures)
    write_csv(output_path(OUTPUT, 's5_mirror_control_points.csv'),[{'index':i,'kx_um_inv':k[0],'ky_um_inv':k[1],
              'is_original19':i<19} for i,k in enumerate(k_union)])
    np.savez_compressed(output_path(OUTPUT, 's5_bloch_validation.npz'),k_um_inv=d['k19'],A_Hz_projection=a,
                        hard_aperture_reference_terms_V_m=terms.reshape(19,4,2),
                        direct_full_saved_air_V_m=direct,direct_cavity_air_V_m=cavity)
    weighted_error=np.sqrt(sum(r['Hz_weight_in_all_k']*r['four_band_relative_residual']**2 for r in rows))
    summary={'status':'reference_and_Hz_projection_complete_mechanism_not_established',
             'periodic_references':str(REFERENCES),'point_count':19,
             'internal_Hz_retained_norm_fraction':float(d['P19'].sum()),
             'Hz_reconstruction_residual_relative_full_norm':float(np.sqrt(1-d['P19'].sum()+weighted_error**2)),
             'Gamma_Hz_single_py_residual':rows[0]['single_py_relative_residual'],
             'Gamma_Hz_four_band_residual':rows[0]['four_band_relative_residual'],
             'max_projection_Gram_condition':max(r['normalized_Gram_condition'] for r in rows),
             'gauge_product_relative_error':gauge_error,
             'hard_aperture_origin_orthogonality_error':zero_orthogonality,
             'original19_mirror_closure_max_distance_over_2pi_a':mirror,
             'mirror_closed_control_point_count':len(k_union),
             'finite_direct_full_saved_air_V_m':pairs(direct),'finite_direct_cavity_air_V_m':pairs(cavity),
             'hard_aperture_reference_hypothesis':{
                 'status':'diagnostic_only_not_finite_Maxwell_response',
                 'predicted_K0_V_m':pairs(total),
                 'relative_error_to_cavity_air':float(np.linalg.norm(total-cavity)/np.linalg.norm(cavity)),
                 'nonGamma_contribution_V_m':pairs(terms.reshape(19,4,2)[1:].sum(axis=(0,1))),
                 'limitations':['Hz-only z=0 coefficients','periodic modes have their own eigenfrequencies',
                                'no cladding or finite-boundary response','orthogonal S forbids offGamma origin mixing']},
             'Gamma_controls':controls,
             'scientific_gates':{'periodic_exact_zero':'not established by finite mesh and height checks',
                                  '19_exact_reference_points':'complete',
                                  'full_vector_internal_projection':'missing finite volume/surface vector fields',
                                  'independent_finite_radiation_response':'missing',
                                  'complete_farfield_closure':'not evaluated'},
             'source_sha256':provenance}
    if response_path is not None:
        summary['independent_response_test']=consume_radiation_bundle(response_path,d['k19'])
    if volume_path is not None:
        summary['volume_radiation_closure']=consume_volume_bundle(volume_path,d['k19'])
    if not all(digest(Path(p))==h for p,h in provenance.items()):
        raise RuntimeError('Input source changed during analysis')
    summary['inputs_unchanged']=True
    write_json(output_path(OUTPUT, 's5_validation_summary.json'),summary)
    plot_checks(d,rows,summary,apertures,k_union)
    print(json.dumps({k:summary[k] for k in ['status','internal_Hz_retained_norm_fraction',
          'Gamma_Hz_single_py_residual','Gamma_Hz_four_band_residual','gauge_product_relative_error',
          'mirror_closed_control_point_count','hard_aperture_reference_hypothesis']},indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radiation-bundle',type=Path)
    parser.add_argument('--volume-bundle',type=Path)
    args=parser.parse_args()
    run(args.radiation_bundle,args.volume_bundle)
