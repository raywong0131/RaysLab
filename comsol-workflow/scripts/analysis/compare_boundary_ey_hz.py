"""Read-only two-height comparison of Ey, weighted Hz gradient, and hole boundary sum."""
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
from scripts.analysis import export_boundary_ey_area_scan as ey
from scripts.analysis.audit_s31_boundary import native_integral
from scripts.analysis.decompose_boundary_c6v import scaled_terms
from scripts.analysis.plot_boundary_control import gradient_average
from comsol_workflow.boundary_integrals import hole_edges, polygon_edges, line_rule, polygon_area

OUT=ey.SOURCE/'Ey_Hz_boundary_comparison'
EPS_AIR,EPS_SLAB=1.,3.3**2


def cvalue(row,name):
    return complex(row[name+'_re'],row[name+'_im'])


def prepare():
    parent=ey.prepare()
    config={**parent,'boundary_order':128,'derivative':'elementwise d(laginterp(2,Hz),x)',
        'new_eigensolves':0,'source_Ey_config_sha256':s4.digest(ey.OUT/'config.json'),
        'source_Ey_hashes':{str(p):s4.digest(p) for case in parent['cases']
            for p in [ey.OUT/f"zeta{case['zeta']:g}"/'area_integrals.csv',ey.OUT/f"zeta{case['zeta']:g}"/'export.json']}}
    file=OUT/'config.json'
    if file.exists():
        assert s4.json_read(file)==config,'Source identity changed'
    else:
        s4.write_json(file,config)
    return config


def export_case(client,case):
    target=OUT/f"zeta{case['zeta']:g}"
    if (target/'export.json').exists() and s4.json_read(target/'export.json')['status']=='complete':
        print('REUSE',case['zeta'],flush=True)
        return
    source=Path(case['directory'])
    previous=ey.OUT/f"zeta{case['zeta']:g}"
    native_meta=s4.json_read(previous/'export.json')
    previous_Ey=pd.read_csv(previous/'area_integrals.csv')
    factor=cvalue(previous_Ey.iloc[0],'factor')
    audit=s4.json_read(case['source_metadata'])
    data=np.load(source/'surface_fields.npz')
    area=abs(polygon_area(data['outer_m']))
    if 'coefficients' in audit:
        pairs=[(r['mode_idx'],cvalue(r,'coefficient')) for r in audit['coefficients'] if cvalue(r,'coefficient')!=0]
    else:
        pairs=[(audit['mode_idx'],1+0j)]
    indices,coeff=zip(*pairs)
    coeff=np.array(coeff)
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
        freq=frequencies(model)
        old=pd.read_csv(source/'eigenfrequencies.csv')
        np.testing.assert_allclose(freq[['re','im','q']],old[['re','im','q']],rtol=1e-8,atol=1e-12)
        omegas=2*np.pi*(freq.loc[list(indices),'re'].to_numpy()-1j*freq.loc[list(indices),'im'].to_numpy())*1e12
        sampler=s4.Sampler(SimpleNamespace(model=model))
        check=np.linspace(0,len(data['surface_Hz_raw'])-1,80,dtype=int)
        h=sum(c*sampler.plane(i,data['area_xy_m'][check],100e-9,['ewfd.Hz'])[:,0] for i,c in pairs)
        href=data['surface_Hz_raw'][check]
        hz_error=float(np.linalg.norm(h-href)/np.linalg.norm(href))
        if hz_error>1e-8:
            raise ValueError('Hz field or mode coefficients changed')
        dielectric=set(int(i) for i in comp.selection('sel_slab_dom').entities(3))
        layer=set(int(i) for i in comp.selection('sel_slab_layer').entities(3))
        adjacency=np.asarray(comp.geom('geom1').getUpDown(),int)
        edges=[*hole_edges(data['holes_m']),*polygon_edges(data['outer_m'])]
        rules=[line_rule(e,128) for e in edges]
        points=np.concatenate([p for p,w in rules])
        weights=np.array([w for p,w in rules])
        arrays={'edge_xy_m':points.reshape(24,128,2),'edge_weights_m':weights,'mode_indices':indices,'coefficients':coeff}
        components,edge_rows,summary=[],[],[]
        for z in (0,100):
            faces=native_meta['surface_face_domains'][str(z)]
            assert all(domain in layer and domain in adjacency[:,face-1] for face,domain in faces)
            samples=np.stack([sampler.plane(i,points,z*1e-9,['ewfd.Hz'])[:,0].reshape(24,128) for i in indices],axis=-1)
            arrays[f'Hz_components_z{z}']=samples
            if z==100:
                reproduced=samples[:18]@coeff
                error=np.linalg.norm(reproduced-data['edge_Hz_raw_128'])/np.linalg.norm(data['edge_Hz_raw_128'])
                if error>1e-8:
                    raise ValueError(f'Previously exported 128-point boundary field changed: {error}')
            line=np.einsum('en,enm->em',weights,samples)
            normals=np.array([e.normal[0] for e in edges])
            qx=np.sum(line[:18]*normals[:18,None],axis=0)
            outer=np.sum(line[18:]*normals[18:,None],axis=0)/EPS_SLAB
            zeros=np.zeros(len(indices),complex)
            for j,(mode,c) in enumerate(pairs):
                for k,e in enumerate(edges):
                    edge_rows.append({'zeta':case['zeta'],'z_nm':z,'mode_idx':mode,'coefficient':c,
                        'omega_rad_s':omegas[j],'hole':e.hole,'edge':e.edge,'group':e.group,
                        'nx':e.normal[0],'length_m':e.length,'Hz_integral_raw_A':line[k,j],
                        'nxHz_integral_raw_A':line[k,j]*e.normal[0]})
            for order in (4,8):
                values=[]
                for mode,c in pairs:
                    total=np.zeros(4,complex)
                    selected_area=selected_z=0.
                    for face,domain in faces:
                        eps=EPS_SLAB if domain in dielectric else EPS_AIR
                        expressions=['1',f'side({domain},ewfd.Ey)',
                            f'side({domain},d(laginterp(2,ewfd.Hz),x))/{eps}',
                            f'side({domain},d(laginterp(2,ewfd.Hx),z))/{eps}',
                            f'side({domain},d(laginterp(2,ewfd.Hz),x))','z']
                        v=native_integral(model,'compare_native','IntSurface',[face],expressions,
                            ['m^2','V*m','A','A','A','m^3'],mode,order)
                        selected_area+=v[0].real
                        selected_z+=v[5].real
                        total+=v[1:5]
                        components.append({'zeta':case['zeta'],'z_nm':z,'order':order,'mode_idx':mode,
                            'coefficient':c,'face':face,'domain':domain,'epsilon_r':eps,'area_m2':v[0].real,
                            'Ey_raw_Vm':v[1],'weighted_dxHz_raw_A':v[2],'weighted_dzHx_raw_A':v[3],'dxHz_raw_A':v[4]})
                    np.testing.assert_allclose(selected_area,area,rtol=1e-9,atol=0)
                    np.testing.assert_allclose(selected_z/selected_area,z*1e-9,rtol=1e-10,atol=1e-15)
                    values.append(total)
                values=np.array(values)
                direct=factor*(coeff@values[:,0])
                gradient=factor*area*gradient_average(values[:,1],coeff,omegas,area)
                converted=scaled_terms(qx,zeros,outer,zeros,values[:,2],zeros,omegas,area,EPS_AIR,EPS_SLAB,coeff)
                boundary,outer_term,z_term=[factor*area*converted[k] for k in ('Ly','Oy','Zy')]
                oldrow=previous_Ey[(previous_Ey.z_nm==z)&(previous_Ey.order==order)].iloc[0]
                np.testing.assert_allclose(direct,cvalue(oldrow,'Ey_integral_normalized_Vm'),rtol=1e-8,atol=1e-24)
                summary.append({'zeta':case['zeta'],'z_nm':z,'order':order,'edge_order':128,
                    'area_m2':area,'Ey_Vm':direct,'Hz_gradient_Vm':gradient,'boundary_Vm':boundary,
                    'outer_Vm':outer_term,'thickness_Vm':z_term,'full_curl_Vm':gradient+z_term,
                    'boundary_full_Vm':boundary+outer_term+z_term,'curl_residual_Vm':direct-gradient-z_term,
                    'boundary_gradient_residual_Vm':gradient-boundary-outer_term,
                    'weighted_dxHz_raw_A':coeff@values[:,1],'dxHz_raw_A':coeff@values[:,3],
                    'weighted_dzHx_raw_A':coeff@values[:,2],'boundary_sum_raw_A':coeff@qx,
                    'weighted_outer_raw_A':coeff@outer,'factor':factor,
                    'field_kind':audit.get('field_kind','single_eigenmode')})
        s4.write_csv(target/'native_components.csv',components)
        s4.write_csv(target/'edge_components.csv',edge_rows)
        s4.write_csv(target/'comparison.csv',summary)
        np.savez_compressed(target/'edge_fields.npz',**arrays)
        assert s4.digest(case['source_model'])==case['source_model_sha256']
        s4.write_json(target/'export.json',{'status':'complete','zeta':case['zeta'],'source_model_unchanged':True,
            'new_eigensolves':0,'Hz_reproduction_error':hz_error,'surface_boundary_reproduction_error':float(error)})
        print(f"EXPORTED zeta={case['zeta']}",flush=True)
    finally:
        if model is not None:
            client.remove(model)


def plot(config):
    all_rows=pd.concat([pd.read_csv(OUT/f"zeta{c['zeta']:g}"/'comparison.csv') for c in config['cases']],ignore_index=True)
    fine=all_rows[all_rows.order==8].sort_values(['zeta','z_nm']).copy()
    coarse=all_rows[all_rows.order==4].sort_values(['zeta','z_nm'])
    assert len(fine)==len(coarse)==60
    columns=[c for c in fine if c.endswith(('_re','_im'))]
    assert np.isfinite(fine[columns]).all().all()
    changes={}
    for name in ('Ey_Vm','Hz_gradient_Vm','thickness_Vm'):
        v=lambda frame:frame[name+'_re'].to_numpy()+1j*frame[name+'_im'].to_numpy()
        changes[name]=float(np.max(abs(v(fine)-v(coarse)))/np.max(abs(v(fine))))
        if changes[name]>1e-7:
            raise ValueError(f'Native quadrature not stable: {name} {changes[name]}')
    fine.to_csv(OUT/'comparison_vs_zeta.csv',index=False)
    all_rows.to_csv(OUT/'quadrature_comparison.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    for name,curves in (
        ('three_way_comparison',[('Ey_Vm',r'$\int E_y\,dA$','#252525','o'),
            ('Hz_gradient_Vm',r'$H_z$ x-gradient term','#387ba8','s'),('boundary_Vm','Hole-boundary term (128)','#c24a42','^')]),
        ('full_balance_check',[('Ey_Vm',r'$\int E_y\,dA$','#252525','o'),
            ('full_curl_Vm','Hz gradient + thickness','#387ba8','s'),('boundary_full_Vm','Hole + outer + thickness','#c24a42','^')])):
        fig,axes=plt.subplots(2,2,figsize=(11,7.4),sharex=True,sharey='col',layout='constrained')
        for row,z in enumerate((0,100)):
            t=fine[fine.z_nm==z]
            for col,part in enumerate(('re','im')):
                ax=axes[row,col]
                for key,label,color,marker in curves:
                    ax.plot(t.zeta,t[f'{key}_{part}'],color=color,marker=marker,ms=3,lw=1.1,label=label)
                ax.axhline(0,color='.8',lw=.7)
                ax.set(title=f'{chr(97+2*row+col)}  z = {z} nm | '+('Real part' if part=='re' else 'Imaginary part'),
                    ylabel=('Re' if part=='re' else 'Im')+' contribution (V m)')
                ax.set_xticks([1,1.05,1.1,1.15,1.2])
                ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
                ax.margins(x=.025,y=.12)
                if row==1: ax.set_xlabel(r'$\zeta$ (1)')
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside upper center',ncol=3,frameon=False,fontsize=10)
        fig.supxlabel(r'$\eta=0.96$ | mesh 5 | Common phase and surface $H_{z,\mathrm{RMS}}=1$ A/m | Each mode retains its own complex frequency',fontsize=9)
        for ext in ('png','pdf'):
            fig.savefig(OUT/f'{name}.{ext}',dpi=200)
        plt.close(fig)
    value=lambda key:fine[key+'_re'].to_numpy()+1j*fine[key+'_im'].to_numpy()
    scale=float(max(abs(value('Ey_Vm'))))
    summary={'status':'complete','case_count':30,'height_count':2,'new_eigensolves':0,
        'max_quadrature_change_over_scan_scale':changes,
        'max_curl_residual_over_Ey_scan_scale':float(max(abs(value('curl_residual_Vm')))/scale),
        'max_boundary_gradient_residual_over_Ey_scan_scale':float(max(abs(value('boundary_gradient_residual_Vm')))/scale),
        'max_outer_over_Ey_scan_scale':float(max(abs(value('outer_Vm')))/scale),
        'max_thickness_over_Ey_scan_scale':float(max(abs(value('thickness_Vm')))/scale)}
    s4.write_json(OUT/'analysis_summary.json',summary)
    report=f'''# Ey、Hz导数与边界积分的统一比较

使用当前30个ζ点、z=0与100 nm、η=0.96、b=245 nm、mesh=5、Γ点用户py的现有解，无新求解。完整晶胞面积含孔内空气，采用明确板层侧。
沿用此前相位，表面Hz RMS统一为1 A/m；两个高度及所有项共用同一系数。主图的三条曲线均为V·m，表示面积积分而非平均值。

采用exp(-iωt)，εr为空气1、介质3.3²，法向从介质指向孔内空气。
A = ∫Ey dS。
D = -i/(ωε0) ∫[∂Hz/∂x]/εr dS。
B = +i/(ωε0)(1/εair-1/εslab) Σ孔边∫nx Hz dl。
Z = +i/(ωε0) ∫[∂Hx/∂z]/εr dS。
O = -i/(ωε0) ∮外边界 nx Hz/εr dl。

连续场满足A=D+Z；若界面连续性及分部积分条件满足，则D=B+O。三维模型中的Z不能预先丢弃，A、D、B也不应强制相等。原始导数积分和裸边界和单位A，保留于CSV；它们不能直接与单位V·m的Ey积分叠加。

导数沿用d(laginterp(2,H),x或z)逐单元重构，面积使用COMSOL原生积分阶数4/8。孔与外边每条128点高斯积分，保留与前次完全相同的界面采样定义。z=100 nm的边界场已与之前128点NPZ核验，直接Ey已与前次双高度导出核验。
ζ=1、1.01、1.02为近简并模常数复系数组合，各分量使用自己的复频率换算后才求和，不给组合场指定一个频率或Q。

4/8阶最大变化相对各量全扫描幅度：{changes}。
最大完整旋度残差|A-D-Z|/max|A|={summary['max_curl_residual_over_Ey_scan_scale']:.6g}。
最大导数与边界差|D-B-O|/max|A|={summary['max_boundary_gradient_residual_over_Ey_scan_scale']:.6g}。
最大外边界项/max|A|={summary['max_outer_over_Ey_scan_scale']:.6g}；最大厚度项/max|A|={summary['max_thickness_over_Ey_scan_scale']:.6g}。
这些分母均为全扫描Ey最大复幅度，避免接近零处相对误差发散。积分阶数稳定不等于mesh5收敛；重构导数、界面迹及单元间场不连续可能造成差异，不把差异自动归因于物理效应。

three_way_comparison.png/PDF：上排z=0、下排z=100 nm，左实部、右虚部，每幅A、D、B。
full_balance_check.png/PDF：相同排布，比较A、D+Z、B+O+Z，供判断遗漏项和数值差异。
comparison_vs_zeta.csv：60组最终复数结果；quadrature_comparison.csv：阶数对照。各ζ目录保留逐面逐模导数、逐边原始积分及128点Hz样本。所有线仅连接真实点，无拟合或强制零点。
'''
    (OUT/'report.md').write_text(report,encoding='utf-8')
    print(summary,flush=True)
    print(fine[fine.zeta.isin([1,1.156,1.2])][['zeta','z_nm','Ey_Vm_re','Ey_Vm_im','Hz_gradient_Vm_re','Hz_gradient_Vm_im','boundary_Vm_re','boundary_Vm_im']].to_string(index=False),flush=True)


def run(config):
    record={'status':'starting','new_eigensolves':0,'completed_zetas':[],'preflight':s4.runtime_preflight()}
    s4.write_json(OUT/'execution.json',record)
    import mph
    client=mph.start(version='6.3',cores=4)
    try:
        if not all(client.java.checkoutLicense(p) for p in ('COMSOL','WAVEOPTICS')):
            raise RuntimeError('COMSOL/Wave Optics checkout failed')
        for case in config['cases']:
            export_case(client,case)
            record['completed_zetas'].append(case['zeta'])
            s4.write_json(OUT/'execution.json',record)
        plot(config)
        prepare()
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
    print(f'Prepared 30 saved models, two heights, zero solves. Output: {OUT}',flush=True)
    if args.plot_only:
        plot(config)
    elif not args.prepare_only:
        run(config)
