"""Read saved Gamma modes and integrate Ey over two actual horizontal surfaces."""
from __future__ import annotations
import argparse
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
from scripts.run_sweep.run_boundary_surface_signed_scan import ROOT as SOURCE, REFERENCE, prepare as scan_manifest
from scripts.analysis.audit_s31_boundary import native_integral
from comsol_workflow.boundary_integrals import phase_and_norm, polygon_area, complex_columns

OUT=SOURCE/'Ey_area_z0_z100'


def prepare():
    manifest=scan_manifest()
    cases=[]
    for c in manifest['cases']:
        directory=Path(c.get('source_directory',c['directory']))
        meta=directory/('execution.json' if c['zeta']==1 else 'export.json')
        audit=s4.json_read(meta)
        model=directory/'s4_gamma.mph'
        if not model.exists():
            model=Path(audit['source'])
        assert model.is_file() and audit['status']=='complete'
        cases.append({'zeta':c['zeta'],'directory':str(directory),'source_model':str(model),
            'source_metadata':str(meta),'source_model_sha256':s4.digest(model),
            'source_fields_sha256':s4.digest(directory/'surface_fields.npz'),
            'source_metadata_sha256':s4.digest(meta)})
    config={'cases':cases,'z_nm':[0,100],'native_orders':[4,8], 'new_eigensolves':0,
        'eta':.96,'b0_nm':245,'mesh':5,'user_mode':'py','k_over_G':[0,0],
        'integration':'whole cell including hole air; slab-layer side',
        'normalization':'surface Hz RMS = 1 A/m; same complex factor at both heights',
        'phase_reference_zeta':1.156,'phase_convention':'exp(-i omega t)'}
    path=OUT/'config.json'
    if path.exists():
        assert s4.json_read(path)==config,'Source/config identity changed'
    else:
        s4.write_json(path,config)
    return config


def export_case(client,case):
    target=OUT/f"zeta{case['zeta']:g}"
    done=target/'export.json'
    if done.exists() and s4.json_read(done)['status']=='complete':
        print('REUSE',case['zeta'],flush=True)
        return
    from jpype.types import JInt
    directory=Path(case['directory'])
    audit=s4.json_read(case['source_metadata'])
    data=np.load(directory/'surface_fields.npz')
    ref=np.load(REFERENCE/'surface_fields.npz')
    h,w=data['surface_Hz_raw'],data['area_weights_m2']
    np.testing.assert_allclose(data['area_xy_m'],ref['area_xy_m'],atol=1e-20,rtol=0)
    anchor=ref['surface_Hz_raw']*phase_and_norm(ref['surface_Hz_raw'],ref['area_weights_m2'])
    factor=phase_and_norm(h,w,anchor)  # Scale the surface Hz RMS to 1 A/m.
    expected_area=abs(polygon_area(data['outer_m']))
    if 'coefficients' in audit:
        indices=[r['mode_idx'] for r in audit['coefficients']]
        coeff=np.array([complex(r['coefficient_re'],r['coefficient_im']) for r in audit['coefficients']])
    else:
        indices=[audit['mode_idx']]
        coeff=np.ones(1,complex)
    pairs=[(i,c) for i,c in zip(indices,coeff) if c!=0]
    model=None
    try:
        print(f"LOAD zeta={case['zeta']}",flush=True)
        model=client.load(case['source_model'])
        jm=model.java
        comp=jm.component('comp1')
        np.testing.assert_allclose([jm.param().evaluate(k) for k in ('a','H','kx','ky')],[820e-9,200e-9,0,0],atol=1e-18)
        assert mesh_identity(model)==audit['mesh']
        for tag in ('s4interp','s4surface_freq'):
            if tag in list(jm.result().numerical().tags()):
                jm.result().numerical().remove(tag)
        actual=frequencies(model)
        expected=pd.read_csv(directory/'eigenfrequencies.csv')
        np.testing.assert_allclose(actual[['re','im','q']],expected[['re','im','q']],rtol=1e-8,atol=1e-12)
        sampler=s4.Sampler(SimpleNamespace(model=model))
        check=np.linspace(0,len(h)-1,80,dtype=int)
        reproduced=sum(c*sampler.plane(i,data['area_xy_m'][check],100e-9,['ewfd.Hz'])[:,0] for i,c in pairs)
        relative_error=float(np.linalg.norm(reproduced-h[check])/np.linalg.norm(h[check]))
        if relative_error>1e-8:
            raise ValueError(f'Saved field/mode coefficients do not reproduce Hz: {relative_error}')
        layer=set(int(i) for i in comp.selection('sel_slab_layer').entities(3))
        updown=np.asarray(comp.geom('geom1').getUpDown(),int)
        surface_rows=[]
        component_rows=[]
        selections={}
        for z_nm in (0,100):
            tag=f'eyarea_z{z_nm}'
            if tag in list(comp.selection().tags()):
                comp.selection().remove(tag)
            selection=comp.selection().create(tag,'Box')
            selection.set('entitydim',JInt(2))
            selection.set('condition','inside')
            for key,value in {'xmin':-1.,'xmax':1.,'ymin':-1.,'ymax':1.,
                              'zmin':z_nm/1000-1e-8,'zmax':z_nm/1000+1e-8}.items():
                selection.set(key,value)
            faces=[int(i) for i in selection.entities(2)]
            if not faces:
                raise ValueError(f'No horizontal faces at z={z_nm}')
            face_domains=[]
            for face in faces:
                domains=layer.intersection(updown[:,face-1])
                if len(domains)!=1:
                    raise ValueError(f'Expected one slab-layer side on face {face}: {domains}')
                face_domains.append((face,int(domains.pop())))
            selections[str(z_nm)]=face_domains
            for order in (4,8):
                total=0j
                integrated_area=0.
                integrated_z=0.
                for mode,c in pairs:
                    for face,domain in face_domains:
                        value=native_integral(model,'eyarea_int','IntSurface',[face],
                            ['1',f'side({domain},ewfd.Ey)','z'],['m^2','V*m','m^3'],mode,order)
                        total+=c*value[1]
                        if mode==pairs[0][0]:
                            integrated_area+=value[0].real
                            integrated_z+=value[2].real
                        component_rows.append(complex_columns({'zeta':case['zeta'],'z_nm':z_nm,'order':order,
                            'mode_idx':mode,'face':face,'slab_domain':domain,'coefficient':c,
                            'area_m2':value[0].real,'Ey_integral_raw_Vm':value[1]}))
                np.testing.assert_allclose(integrated_area,expected_area,rtol=1e-9,atol=0)
                np.testing.assert_allclose(integrated_z/integrated_area,z_nm*1e-9,atol=1e-15,rtol=1e-10)
                surface_rows.append(complex_columns({'zeta':case['zeta'],'z_nm':z_nm,'order':order,
                    'area_m2':integrated_area,'Ey_integral_raw_Vm':total,
                    'Ey_integral_normalized_Vm':total*factor,'factor':factor,
                    'Hz_surface_RMS_raw_Am':float(np.sqrt(w@abs(h)**2/w.sum())),
                    'field_kind':audit.get('field_kind','single_eigenmode')}))
        target.mkdir(parents=True,exist_ok=True)
        s4.write_csv(target/'native_components.csv',component_rows)
        s4.write_csv(target/'area_integrals.csv',surface_rows)
        assert s4.digest(case['source_model'])==case['source_model_sha256']
        s4.write_json(done,{'status':'complete','zeta':case['zeta'],'new_eigensolves':0,
            'surface_face_domains':selections,'source_model_unchanged':True,
            'Hz_reproduction_relative_error':relative_error})
        print(f"EXPORTED zeta={case['zeta']}",flush=True)
    finally:
        if model is not None:
            client.remove(model)


def plot(config):
    raw=pd.concat([pd.read_csv(OUT/f"zeta{c['zeta']:g}"/'area_integrals.csv') for c in config['cases']],ignore_index=True)
    value=lambda t:t.Ey_integral_normalized_Vm_re.to_numpy()+1j*t.Ey_integral_normalized_Vm_im.to_numpy()
    fine=raw[raw.order==8].sort_values(['zeta','z_nm']).copy()
    coarse=raw[raw.order==4].sort_values(['zeta','z_nm'])
    assert len(fine)==60 and len(coarse)==60
    np.testing.assert_allclose(fine[['zeta','z_nm']],coarse[['zeta','z_nm']])
    delta=abs(value(fine)-value(coarse))
    scale=abs(value(fine)).max()
    if not np.isfinite(value(fine)).all() or scale==0 or delta.max()>1e-7*scale:
        raise ValueError(f'Native quadrature not stable: maximum change/scan scale {delta.max()/scale}')
    fine['quadrature_change_abs_Vm']=delta
    fine.to_csv(OUT/'Ey_area_vs_zeta.csv',index=False)
    raw.to_csv(OUT/'Ey_area_quadrature.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for ax,key,title in zip(axes,('re','im'),('a  Real part','b  Imaginary part')):
        for z,color,marker in ((0,'#387ba8','o'),(100,'#c24a42','s')):
            t=fine[fine.z_nm==z]
            ax.plot(t.zeta,t[f'Ey_integral_normalized_Vm_{key}'],color=color,marker=marker,
                ms=3.5,lw=1.3,label=f'$z={z}$ nm')
        ax.set(title=title,xlabel=r'$\zeta$ (1)',ylabel=(r'Re $\int_S E_y\,dA$ (V m)' if key=='re' else r'Im $\int_S E_y\,dA$ (V m)'))
        ax.axhline(0,color='.8',lw=.7)
        ax.set_xticks([1,1.05,1.1,1.15,1.2])
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
        ax.margins(x=.025,y=.1)
        ax.legend(frameon=False)
    fig.suptitle(r'$\eta=0.96$ | $b=245$ nm | $p_y$ at $\Gamma$ | mesh 5 | Whole-cell area integral',fontsize=12)
    fig.supxlabel(r'Both heights share the same phase and amplitude: surface $H_{z,\mathrm{RMS}}=1$ A/m. Native surface quadrature.',fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'Ey_area_vs_zeta.{ext}',dpi=200)
    plt.close(fig)
    s4.write_json(OUT/'analysis_summary.json',{'status':'complete','case_count':30,'height_count':2,
        'maximum_quadrature_change_Vm':float(delta.max()),'maximum_change_over_scan_scale':float(delta.max()/scale),
        'area_m2':float(fine.area_m2.iloc[0]),'normalization':config['normalization'],
        'new_eigensolves':0,'mixed_state_zetas':[1.,1.01,1.02]})
    report=f'''# Ey截面面积分随ζ变化

使用当前mesh=5、η=0.96、b=245 nm、a=820 nm、H=200 nm、折射率3.3、Γ点用户py模式的30个已保存解。没有新建网格或重新求解，也没有保存源MPH。

S为完整六角晶胞横截面，包含介质和孔内空气，面积{fine.area_m2.iloc[0]:.12g} m²。在z=0与z=100 nm的实际几何面上按板层侧逐域积分，检查了总面积及积分z坐标。这里是∫Ey dA，不是面积平均，不是∫|Ey|² dA。

本征模振幅任意，直接比较未经归一化的原始积分没有固定幅度意义。因此沿用此前的Hz相位参考（ζ=1.156），每个ζ将z=100 nm截面的Hz RMS设为1 A/m，两个高度共用同一复数因子。图中的面积分单位为V·m，是在此明确幅度约定下的数值；原始积分和缩放因子同时保留于CSV。没有对每个Ey积分单独旋转相位。

沿用原来模式及系数，ζ=1、1.01、1.02为近简并模常数复系数组合，不是严格单一本征模。其余点为单一本征模。导出前用Hz复数采样验证了保存的模式/系数/幅度一致。

面积分采用COMSOL有限元面上的原生积分，阶数4和8最大差值为{delta.max():.6g} V·m，相对于整个扫描最大复幅度为{delta.max()/scale:.6g}。这只检查面积求积稳定性，不表示mesh=5已收敛。此前128点设置属于线积分，不能直接作为二维面积分规则。

Ey_area_vs_zeta.csv为阶数8的60行结果；Ey_area_quadrature.csv保留两个阶数；各ζ子目录保留逐面逐模原始积分。PNG/PDF左图为实部，右图为虚部，各比较两个高度；无拟合或强制零点。
'''
    (OUT/'report.md').write_text(report,encoding='utf-8')
    print(fine[fine.zeta.isin([1,1.155,1.156,1.2])][['zeta','z_nm','Ey_integral_normalized_Vm_re','Ey_integral_normalized_Vm_im']].to_string(index=False),flush=True)


def run(config):
    record={'status':'starting','new_eigensolves':0,'completed_zetas':[], 'preflight':s4.runtime_preflight()}
    s4.write_json(OUT/'execution.json',record)
    import mph
    client=mph.start(version='6.3',cores=4)
    try:
        if not all(client.java.checkoutLicense(p) for p in ('COMSOL','WAVEOPTICS')):
            raise RuntimeError('COMSOL/Wave Optics license unavailable')
        for case in config['cases']:
            export_case(client,case)
            record['completed_zetas'].append(case['zeta'])
            s4.write_json(OUT/'execution.json',record)
        plot(config)
        prepare()  # Recheck all saved source hashes.
        record.update(status='complete',analysis_completed=True,source_hashes_unchanged=True)
    except Exception as exc:
        record.update(status='failed',error=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        client.clear()
        s4.write_json(OUT/'execution.json',record)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--plot-only',action='store_true')
    args=parser.parse_args()
    config=prepare()
    print(f'Prepared {len(config["cases"])} saved models, zero solves; output={OUT}',flush=True)
    if args.plot_only:
        plot(config)
    elif not args.prepare_only:
        run(config)
