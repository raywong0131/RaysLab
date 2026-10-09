"""Single zeta=1 Gamma control, compared with the existing zeta=1.156 surface export."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import time
from pathlib import Path
from types import SimpleNamespace
import traceback

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scripts.run_main import run_boundary_analysis as s4
from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies, mesh_identity
from scripts.analysis.decompose_boundary_c6v import symmetry_coefficients, weighted_parity, field_pair
from comsol_workflow.geometry_utils import create_hexagon_design
from comsol_workflow.boundary_integrals import hole_edges, line_rule, group_totals, phase_and_norm, complex_columns

ROOT=s4.ROOT/'scripts/.out/unit_cell_2D/S4_surface_eta0.96_zeta1_vs1.156_20260918'
REFERENCE=s4.ROOT/'scripts/.out/unit_cell_2D/S4_surface_eta0.96_zeta_scan_20260918/zeta1.156'


def run():
    if (ROOT/'execution.json').exists() and s4.json_read(ROOT/'execution.json').get('status')=='complete':
        analyze()
        return
    record={'status':'starting','eta':.96,'zeta':1.,'z_nm':100,'new_eigensolves':0,
            'phase_convention':'exp(-i omega t)','edge_orders':[128,256],
            'source_hashes':{str(p):s4.digest(p) for p in REFERENCE.iterdir() if p.is_file()}}
    record['preflight']=s4.runtime_preflight()
    s4.write_json(ROOT/'execution.json',record)
    # Importing these capabilities is restricted to the explicitly requested solve.
    import mph
    from scripts.run_main import run_band_pair as band
    from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationConfig
    client=mph.start(version='6.3',cores=4)
    runner=None
    model=None
    try:
        if not all(client.java.checkoutLicense(p) for p in ('COMSOL','WAVEOPTICS')):
            raise RuntimeError('COMSOL/Wave Optics license checkout failed')
        record['product_checkout']='COMSOL and WAVEOPTICS passed'
        params=s4.load_shared_parameters(ROOT/'99_config/parameter.json')
        hp=band.get_hole_params(params.unit_cell_2d.cell.to_fourier_params('zeta1_surface'))
        outer,holes,info=create_hexagon_design(band.A,hp)
        np.testing.assert_allclose(hp[:,0],.82*.96/3,rtol=1e-12)
        np.testing.assert_allclose(hp[:,2],.245,rtol=1e-12)
        if min(i['min_dist'] for i in info)<band.D:
            raise ValueError('Geometry clearance failed')
        config=SimulationConfig(slab_height='200 [nm]',slab_refractive_index='3.3',
             eigenmode_count=2,mesh_auto_size=5,eigenfrequency_shift='c_const/1.55[um]')
        s4.write_json(ROOT/'99_config/realized_geometry.json',{'length_unit':'um','outer':outer.tolist(),
          'holes':np.asarray(holes).tolist(),'hole_params':hp.tolist(),'simulation_config':asdict(config),'k_over_G':[0,0]})
        checkpoint=ROOT/'s4_gamma.mph'
        if checkpoint.exists():
            runner=band.ReusableSimulationRun(config)
            client.remove(runner.model)
            model=client.load(str(checkpoint))
            runner.model=model
            runner.plane_datasets={'center':'cpl1','air':'cpl2','yz':'cpl3','xz':'cpl4'}
            print('Reusing this task Gamma checkpoint; no new solve',flush=True)
        else:
            if (ROOT/'solve_attempt.json').exists():
                raise RuntimeError('Previous uncheckpointed solve attempt; inspect before retry')
            s4.write_json(ROOT/'solve_attempt.json',{'eta':.96,'zeta':1,'mesh5':True,'Gamma':True})
            print('Building eta=.96 zeta=1 mesh5; one Gamma solve',flush=True)
            runner=band.ReusableSimulationRun(config)
            runner.model.java # initialize model before attaching progress
            client.java.showProgress(str(ROOT/'comsol_progress.log'))
            runner.build_and_run(.82,np.asarray(holes).tolist(),{'kx':0.,'ky':0.},mesh_auto_size=5)
            model=runner.model
            record['new_eigensolves']=1
            model.save(str(checkpoint))
            print('Solved; Gamma checkpoint saved',flush=True)
        record['mesh']=mesh_identity(model)
        modes=band.persist_k_point_solution(runner,ROOT/'01_results')
        for tag in ('s4interp','s4surface_freq'):
            if tag in list(model.java.result().numerical().tags()):
                model.java.result().numerical().remove(tag)
        freq=frequencies(model)
        modes[['re','im','q']]=freq[['re','im','q']]
        modes.to_csv(ROOT/'eigenfrequencies.csv',index=False)
        indices=list(modes.index[(modes.p_weight>.9)&modes.is_valid])
        if len(indices)!=2:
            raise ValueError('Expected two p-dominant modes')
        sampler=s4.Sampler(SimpleNamespace(model=model))
        ref=np.load(REFERENCE/'surface_fields.npz')
        xy=ref['center_probe_xy_m']
        fields=[field_pair(sampler,indices,xy*r,0.,['ewfd.Hz'])[:,0]
                for r in ([1,1],[-1,1],[1,-1])]
        audit=[]
        for j,idx in enumerate(indices):
            errors=weighted_parity(*(h[:,j] for h in fields),np.ones(len(xy)))
            audit.append({'mode_idx':idx,'p_weight':modes.loc[idx,'p_weight'],
                          'frequency_thz':modes.loc[idx,'re'],'Q':modes.loc[idx,'q'],**errors})
        s4.write_csv(ROOT/'mode_audit.csv',audit)
        good=[j for j,r in enumerate(audit) if max(r['odd_x_error'],r['even_y_error'])<.05]
        if len(good)==1:
            coeff=np.zeros(2,complex);coeff[good[0]]=1
            record['field_kind']='single_eigenmode'
        else:
            split=np.ptp(modes.loc[indices,'re'])/np.mean(modes.loc[indices,'re'])
            if split>1e-3:
                raise ValueError('p pair not near-degenerate; refusing a mixed-frequency target')
            coeff=symmetry_coefficients(*fields,np.ones(len(xy)))
            record['field_kind']='constant_complex_p_subspace_combination_not_exact_single_eigenmode'
        record['mode_indices']=[int(i) for i in indices]
        record['coefficients']=[complex_columns({'mode_idx':i,'coefficient':v}) for i,v in zip(indices,coeff)]
        record['training_parity']=weighted_parity(*(h@coeff for h in fields),np.ones(len(xy)))
        area_xy=ref['area_xy_m'];weights=ref['area_weights_m2']
        validation=[field_pair(sampler,indices,area_xy*r,0.,['ewfd.Hz'])[:,0]
                    for r in ([1,1],[-1,1],[1,-1])]
        record['validation_parity']=weighted_parity(*(h@coeff for h in validation),weights)
        if max(record['validation_parity'].values())>.06:
            raise ValueError('Target parity not resolved on independent points')
        center=validation[0]@coeff
        omega=modes.loc[indices,'re'].to_numpy()-1j*modes.loc[indices,'im'].to_numpy()
        omega_bar=np.dot(weights*center.conj(),validation[0]@(coeff*omega))/np.dot(weights,abs(center)**2)
        record['relative_frequency_spread']=float(np.sqrt(np.dot(weights,abs(validation[0]@(coeff*(omega-omega_bar)))**2)/np.dot(weights,abs(center)**2))/abs(omega_bar))
        record['p_pair_frequency_thz']=modes.loc[indices,'re'].tolist()
        record['p_pair_Q']=modes.loc[indices,'q'].tolist()
        record['frequency_split_GHz']=float(np.ptp(modes.loc[indices,'re'])*1000)
        surface_pair=field_pair(sampler,indices,area_xy,100e-9,['ewfd.Hz'])[:,0]
        surface=surface_pair@coeff
        hx=field_pair(sampler,indices,area_xy*[-1,1],100e-9,['ewfd.Hz'])[:,0]@coeff
        hy=field_pair(sampler,indices,area_xy*[1,-1],100e-9,['ewfd.Hz'])[:,0]@coeff
        record['surface_parity']=weighted_parity(surface,hx,hy,weights)
        actual=sampler.plane(indices[0],area_xy[:5],100e-9,['x','y','z']).real
        np.testing.assert_allclose(actual,np.column_stack([area_xy[:5],np.full(5,100e-9)]),rtol=1e-9,atol=1e-14)
        arrays={'area_xy_m':area_xy,'area_weights_m2':weights,'surface_Hz_raw':surface,
                'surface_Hz_components':surface_pair,'coefficients':coeff,'mode_indices':indices,
                'outer_m':np.asarray(outer)*1e-6,'holes_m':np.asarray(holes)*1e-6,'z_nm':100.,
                'center_training_fields':np.asarray(fields),'center_validation_fields':np.asarray(validation)}
        edges=hole_edges(np.asarray(holes)*1e-6)
        rows=[]
        for order in (128,256):
            rules=[line_rule(e,order) for e in edges]
            points=np.concatenate([p for p,w in rules])
            components=field_pair(sampler,indices,points,100e-9,['ewfd.Hz'])[:,0].reshape(18,order,2)
            h=components@coeff
            arrays[f'edge_xy_m_{order}']=points.reshape(18,order,2)
            arrays[f'edge_weights_m_{order}']=np.asarray([w for p,w in rules])
            arrays[f'edge_Hz_components_{order}']=components
            arrays[f'edge_Hz_raw_{order}']=h
            for j,e in enumerate(edges):
                integral=np.dot(rules[j][1],h[j])
                rows.append({'order':order,'hole':e.hole,'edge':e.edge,'group':e.group,
                    'length_m':e.length,'nx':e.normal[0],'ny':e.normal[1],
                    'integral_hz':integral,'qx':integral*e.normal[0],'qy':integral*e.normal[1]})
        np.savez_compressed(ROOT/'surface_fields.npz',**arrays)
        s4.write_csv(ROOT/'edge_integrals_raw.csv',rows)
        record['reference_unchanged']=all(s4.digest(p)==h for p,h in record['source_hashes'].items())
        if not record['reference_unchanged']:
            raise RuntimeError('Reference changed')
        record['status']='export_complete'
        s4.write_json(ROOT/'execution.json',record)
        analyze()
        record['status']='complete'
        print('Complete',record['field_kind'],record['p_pair_Q'],record['validation_parity'],flush=True)
    except Exception as exc:
        record.update(status='failed',error=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        s4.write_json(ROOT/'execution.json',record)
        if model is not None:
            client.remove(model)
        if runner is not None and hasattr(runner,'clear'):
            runner.model=None
        client.clear()


def analyze():
    ref=np.load(REFERENCE/'surface_fields.npz')
    w=ref['area_weights_m2']
    anchor=ref['surface_Hz_raw']*phase_and_norm(ref['surface_Hz_raw'],w)
    rows=[]
    for zeta,directory in ((1.,ROOT),(1.156,REFERENCE)):
        data=np.load(directory/'surface_fields.npz')
        h=data['surface_Hz_raw']
        factor=phase_and_norm(h,w,anchor)/820e-9
        overlap=abs(np.dot(w,anchor.conj()*h))/np.sqrt(np.dot(w,abs(anchor)**2)*np.dot(w,abs(h)**2))
        raw=pd.read_csv(directory/'edge_integrals_raw.csv')
        for order in (128,256):
            edge_rows=[{**r,'qx':complex(r['qx_re'],r['qx_im']),'qy':complex(r['qy_re'],r['qy_im'])}
                       for r in raw[raw.order==order].to_dict('records')]
            t=group_totals(edge_rows)
            rows.append(complex_columns({'zeta':zeta,'order':order,'Hz_surface_rms':np.sqrt(np.dot(w,abs(h)**2)/sum(w)),
                'phase_rotation_rad':np.angle(factor),'overlap_reference':overlap,
                'f_A_raw':t['QAx'],'f_B_raw':t['QBx'],'f_raw':t['Qx'],
                'F_A':t['QAx']*factor,'F_B':t['QBx']*factor,'F':t['Qx']*factor,'F_abs':abs(t['Qx']*factor),
                'cancellation_ratio':t['cancellation_ratio']}))
    table=pd.DataFrame(rows)
    table.to_csv(ROOT/'comparison.csv',index=False)
    e=s4.json_read(ROOT/'execution.json')
    qref=s4.json_read(REFERENCE/'export.json')['Q']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(9.5,6.5),layout='constrained')
    a,b,c,d=axes.flat
    main=table[table.order==256];alt=table[table.order==128];x=np.arange(2)
    for offset,key,label,color in ((-.24,'F_A_re',r'Re $F_A$','#c24a42'),(0,'F_B_re',r'Re $F_B$','#387ba8'),(.24,'F_re',r'Re $F$','#252525')):
        a.bar(x+offset,main[key],width=.22,label=label,color=color)
    a.set(title='a  A/B contributions',ylabel='Normalized boundary sum (1)');a.legend(frameon=False,fontsize=9)
    b.bar(x-.13,main.F_re,width=.24,label=r'Re $F$',color='#252525')
    b.bar(x+.13,main.F_im,width=.24,label=r'Im $F$',color='#b16d2e')
    b.set(title='b  Complex total',ylabel='Normalized boundary sum (1)');b.legend(frameon=False)
    c.bar(x-.13,main.cancellation_ratio,width=.24,label='Gauss 256',color='#252525')
    c.bar(x+.13,alt.cancellation_ratio,width=.24,label='Gauss 128',color='#a5adb3')
    c.set(title='c  Cancellation residual',ylabel=r'$|f|/(|f_A|+|f_B|)$ (1)');c.legend(frameon=False)
    d.scatter(np.full(len(e['p_pair_Q']),0),e['p_pair_Q'],c='#854f99',s=35,label=r'$\zeta=1$: original p pair')
    d.scatter([1],[qref],c='#252525',s=35,label=r'$\zeta=1.156$: selected $p_y$')
    d.set(title=r'd  Original eigenmode $Q$',ylabel=r'$Q$ (1)',yscale='log');d.legend(frameon=False,fontsize=8,loc='center right')
    for ax in (a,b): ax.axhline(0,color='.75',lw=.7)
    for ax in axes.flat:
        ax.set_xticks(x,['1','1.156']);ax.set_xlabel(r'$\zeta$ (1)');ax.set_xlim(-.6,1.6)
    fig.suptitle(r'$\eta=0.96$ | $z=100$ nm, slab-side $H_z$ | mesh 5',fontsize=13)
    fig.supxlabel(r'$F=f/(aH_{z,\mathrm{RMS}})$; common Hz-overlap phase. No Q assigned to a mixed eigenstate.',fontsize=9)
    fig.savefig(ROOT/'zeta1_vs1.156.png',dpi=200);fig.savefig(ROOT/'zeta1_vs1.156.pdf');plt.close(fig)
    print(table[['zeta','order','F_re','F_im','F_abs','cancellation_ratio','overlap_reference']].to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot-only',action='store_true')
    args=parser.parse_args()
    if args.plot_only: analyze()
    else: run()
