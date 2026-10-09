"""Synthetic cavity-aperture air fields from actual nineteen-k unit-cell solutions."""
from comsol_workflow.output_paths import table_links, table_path
import json
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
from scipy.constants import epsilon_0,mu_0
from matplotlib.path import Path as Polygon
from comsol_workflow.hex_lattice_utils import unit_cell_corners
from scripts.analysis import analyze_dipolar_complete as a


def run():
    cfg=json.loads((a.CONFIG/"config.json").read_text(encoding="utf8"))
    centers=np.array([[p['x'],p['y']] for p in cfg['bulk_points']])
    shifts=centers[np.argsort(np.linalg.norm(centers,axis=1))[:7]]
    source=a.load('finite_air_1100nm_fine_EH.npz');x=source['x_um'];dx=float(x[1]-x[0])*1e-6
    xx,yy=np.meshgrid(x,x);xy=np.column_stack([xx.ravel(),yy.ravel()])
    _,ids=cKDTree(centers).query(xy);local=xy-centers[ids]
    support=Polygon(unit_cell_corners(cfg['a'])).contains_points(local,radius=1e-10)
    pos=xy[support];rho=local[support];cell_R=centers[ids[support]]
    fit=a.load('field_decomposition.npz');coeff=fit['coefficients'].reshape(19,4)[:,1]
    size=2*a.FFT_SIZE-1;dq=2*np.pi/(size*dx);half=int(np.ceil(np.sin(np.deg2rad(10))*a.K0/dq))+1
    q=np.arange(-half,half+1)*dq;qx,qy=np.meshgrid(q/a.K0,q/a.K0);nz=np.sqrt(1-qx*qx-qy*qy)
    nt=np.stack([qx,qy],-1);n=np.stack([qx,qy,nz],-1)
    synthesized=np.zeros((3,len(pos),6),complex);spectra=[];rows=[]
    cached_path=a.D/'Fig6_unitcell_synthesis.npz'
    start=0
    if cached_path.exists():
        cached=a.load(cached_path.name)
        np.testing.assert_array_equal(cached['support'],support.reshape(xx.shape))
        np.testing.assert_allclose(cached['q_m_inv'],q,rtol=0,atol=0)
        start=len(cached['k_E_up_air'])
        assert start in (7,19)
        np.testing.assert_allclose(cached['coefficients'],coeff[:start],rtol=0,atol=0)
        synthesized[:len(cached['weighted_EH_air_on_support'])]=cached['weighted_EH_air_on_support']
        spectra=list(cached['k_E_up_air'])
        rows=pd.read_csv(table_path(a.L, 'Fig6_unitcell_weights.csv', existing=True)).to_dict('records')[:start]
        print(f'Reuse verified synthesis for k00-k{start-1:02d}',flush=True)

    def transform(values):
        tangential=[]
        for j in [0,1,3,4]:
            plane=np.zeros(xx.size,complex);plane[support]=values[:,j]
            tangential.append(a.padded_transform_roi(plane.reshape(xx.shape),dx,size,len(q)))
        et=np.stack(tangential[:2],-1);ht=np.stack(tangential[2:],-1)
        p=np.stack([ht[...,1],-ht[...,0]],-1)
        v=a.Z0/nz[...,None]*(p-nt*np.sum(nt*p,-1)[...,None])
        out=np.empty((*nz.shape,3),complex);out[...,:2]=(et+v)/2
        out[...,2]=-np.sum(nt*out[...,:2],-1)/nz
        # Independent center integral, including outgoing electric/magnetic projection.
        integral=values.sum(0)*dx**2
        expected=(integral[:2]+a.Z0*np.array([integral[4],-integral[3]]))/2
        np.testing.assert_allclose(out[half,half,:2],expected,rtol=1e-9,atol=1e-16)
        return out

    for i in range(start,19):
        print(f'Unit cell k{i:02d}: interpolate periodic envelope, tile cavity, FFT',flush=True)
        d=a.load(f'periodic_k{i:02d}_air1100nm_EH.npz');k=d['k_um_inv']
        np.testing.assert_allclose(k,fit['k_um_inv'].reshape(19,4,2)[i,1],atol=1e-12)
        base=d['xyz_um'][:,:2]
        field=np.concatenate([d['E_V_m'],d['H_A_m']],-1)/np.sqrt(d['cell_EM_energy_J'])
        periodic=field*np.exp(1j*(base@k))[:,None]
        interior=Polygon(unit_cell_corners(cfg['a'])).contains_points(base,radius=-1e-7)
        raw_base=base;raw_periodic=periodic
        base=base[interior];periodic=periodic[interior]
        ghost_points=np.concatenate([base+shift for shift in shifts])
        ghost_values=np.tile(periodic,(7,1))
        # Shared boundary nodes are merged to avoid ambiguous triangulation vertices.
        unique,index,inverse=np.unique(np.round(ghost_points,12),axis=0,return_index=True,return_inverse=True)
        values=np.zeros((len(unique),6),complex);np.add.at(values,inverse,ghost_values)
        values/=np.bincount(inverse)[:,None]
        interp=LinearNDInterpolator(unique,values)
        native_interp=interp(base)
        reconstruction_error=float(np.linalg.norm(native_interp-periodic)/np.linalg.norm(periodic))
        assert reconstruction_error<1e-8
        boundary_error=float(np.linalg.norm(interp(raw_base[~interior])-raw_periodic[~interior])/np.linalg.norm(raw_periodic[~interior]))
        tiled=interp(rho)*np.exp(-1j*(pos@k))[:,None]*coeff[i]
        assert np.isfinite(tiled).all()
        # Algebraic check of the two equivalent Bloch conventions at actual global locations.
        check=(np.exp(-1j*(rho@k))*np.exp(-1j*(cell_R@k)))
        np.testing.assert_allclose(check,np.exp(-1j*(pos@k)),rtol=1e-12,atol=1e-12)
        synthesized[0 if i==0 else (1 if i<7 else 2)]+=tiled
        spectra.append(transform(tiled))
        rows.append(dict(index=i,kx_um_inv=k[0],ky_um_inv=k[1],weight_re=coeff[i].real,weight_im=coeff[i].imag,weight_abs=abs(coeff[i]),frequency_THz=d['frequency_THz'][0],Q=d['frequency_THz'][2],unit_cell_EM_J=float(d['cell_EM_energy_J']),periodic_node_reconstruction_relative_error=reconstruction_error,excluded_boundary_points=int((~interior).sum()),boundary_trace_relative_difference=boundary_error))
    spectra=np.asarray(spectra);combined=np.stack([spectra[0],spectra[1:7].sum(0),spectra[7:19].sum(0)])
    closures=[]
    for group in (1,2):
        direct_ring=transform(synthesized[group]);closure=float(np.linalg.norm(direct_ring-combined[group])/np.linalg.norm(direct_ring));assert closure<1e-10
        closures.append(closure)
    closure=closures[0]
    if start:
        np.testing.assert_allclose(combined[:2],cached['E_up_air'][:2],rtol=0,atol=0)
    h=np.cross(n,combined)/a.Z0
    near_small=(epsilon_0*np.sum(abs(synthesized[...,:3])**2,-1)+mu_0*np.sum(abs(synthesized[...,3:])**2,-1))/4
    near=np.zeros((3,xx.size));near[:,support]=near_small;near=near.reshape(3,*xx.shape)
    finite_near=(epsilon_0*np.sum(abs(source['E_V_m'])**2,-1)+mu_0*np.sum(abs(source['H_A_m'])**2,-1))/4
    near_ref=float(np.nanmax(finite_near))
    ref=a.load('Fig1_EM_angular_spectrum.npz');far_ref=float(ref['normalization_EM_center'])
    far=(epsilon_0*np.sum(abs(combined)**2,-1)+mu_0*np.sum(abs(h)**2,-1))/4/far_ref
    near/=near_ref
    np.savez_compressed(a.D/'Fig6_unitcell_synthesis.npz',x_um=x,support=support.reshape(xx.shape),supported_xy_um=pos,weighted_EH_air_on_support=synthesized,q_m_inv=q,k_E_up_air=spectra,E_up_air=combined,H_up_air=h,near_EM_normalized=near,far_EM_normalized=far,coefficients=coeff,near_reference_J_m3=near_ref,far_reference_EM_center=far_ref)
    pd.DataFrame(rows).to_csv(table_path(a.L, 'Fig6_unitcell_weights.csv'),index=False)
    checks=dict(cavity_cells=len(centers),support_pixels=int(support.sum()),air_z_um=1.1,fft_size=size,k_indices=[[0],list(range(1,7)),list(range(7,19))],py_band_index=1,normalization='same finite reference for all three rows; air full-EM maximum and Fig.1 far full-EM center',near_peaks=near.max((1,2)).tolist(),far_peaks=far.max((1,2)).tolist(),far_centers=far[:,half,half].tolist(),ring_coherent_reconstruction_error=closure,outer_ring_coherent_reconstruction_error=closures[1],finite_reference_frequency_THz=a.FREQ/1e12,mode_frequencies_THz=[float(r['frequency_THz']) for r in rows],limitations=['cavity aperture only, no cladding feedback','superposition of spatial eigenfields at differing eigenfrequencies is a frozen-field approximation','hard truncation of air plane produces aperture radiation; not a solved finite cavity','Gamma solution has finite numerical residual loss, not an exact zero-radiation field','native periodic-boundary normal-field traces differ; boundary nodes excluded from interpolation, discrepancy recorded'])
    a.write_json('Fig6_unitcell_synthesis_checks.json',checks)
    print(json.dumps(checks,indent=2),flush=True)




def plot_report():
    z=a.load('Fig6_unitcell_synthesis.npz')
    checks=json.loads((a.L/'Fig6_unitcell_synthesis_checks.json').read_text(encoding='utf8'))
    weights=pd.read_csv(table_path(a.L, 'Fig6_unitcell_weights.csv', existing=True))
    cfg=json.loads((a.CONFIG/'config.json').read_text(encoding='utf8'))
    centers=np.array([[p['x'],p['y']] for p in cfg['bulk_points']])
    k=weights[['kx_um_inv','ky_um_inv']].to_numpy()
    lattice_sum=np.exp(-1j*centers@k.T).sum(0)
    assert np.max(abs(lattice_sum[1:])/len(centers))<1e-11
    mid=len(z['q_m_inv'])//2
    per_k_center=epsilon_0/2*np.sum(abs(z['k_E_up_air'][:,mid,mid])**2,-1)/float(z['far_reference_EM_center'])
    weights['lattice_sum_abs_over_N']=abs(lattice_sum)/len(centers)
    weights['far_center_over_finite_center']=per_k_center
    weights.to_csv(table_path(a.L, 'Fig6_unitcell_weights.csv'),index=False)
    a.plt.rcParams.update({'font.size':10,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
    x=z['x_um'];q=z['q_m_inv'];dq=q[1]-q[0]
    edge=np.rad2deg(np.arcsin(np.r_[q-dq/2,q[-1]+dq/2]/a.K0))
    boundary=np.asarray(cfg['inner_boundary']);boundary=np.vstack([boundary,boundary[0]])
    outer_near_peak=float(z['near_EM_normalized'][2].max())
    outer_near_limit=float(np.ceil(outer_near_peak/10**np.floor(np.log10(outer_near_peak)))*10**np.floor(np.log10(outer_near_peak)))
    near_limits=[.003,.3,outer_near_limit]
    far_max=float(z['far_EM_normalized'].max())
    far_upper=float(10**np.ceil(np.log10(far_max)))
    far_lower=far_upper*1e-5
    far_norm=a.LogNorm(far_lower,far_upper)
    assert np.all(z['near_EM_normalized'].max((1,2))<=near_limits)
    assert np.all(z['far_EM_normalized'].max((1,2))<=far_upper)
    fit=a.load('field_decomposition.npz')
    kpoints=fit['k_um_inv'].reshape(19,4,2)[:,1]/(2*np.pi/cfg['a'])
    coeff19=fit['coefficients'].reshape(19,4)[:,1]
    coefficient_weight=abs(coeff19)**2/np.sum(abs(coeff19)**2)
    np.testing.assert_allclose(kpoints*(2*np.pi/cfg['a']),k,atol=1e-12)
    a.write_json('Fig6_kpoint_schematic.json',dict(k_over_2pi_a=kpoints.tolist(),coefficient_weight=coefficient_weight.tolist(),selected_rows=[[0],list(range(1,7)),list(range(7,19))],definition='|a_k,py|^2 / sum_19 |a_k,py|^2; coefficient weight, not additive physical energy',color_limits=[1e-5,1]))
    fig,axs=a.plt.subplots(3,3,figsize=(20,16),layout='constrained')
    for i,label in enumerate(['Gamma only (K0)','First ring (K1-K6), coherent','Outer ring (K7-K18), coherent']):
        selected=np.array([[0],list(range(1,7)),list(range(7,19))][i])
        inactive=np.setdiff1d(np.arange(19),selected)
        ax=axs[i,0];ax.set_facecolor('black')
        ax.scatter(*kpoints[inactive].T,s=55,facecolors='none',edgecolors='#54545d',linewidths=.8)
        dots=ax.scatter(*kpoints[selected].T,s=85,c=coefficient_weight[selected],cmap='magma',norm=a.LogNorm(1e-5,1),linewidths=0)
        for number,point in enumerate(kpoints):
            ax.annotate(str(number),point,xytext=(0,9),textcoords='offset points',ha='center',color='white' if number in selected else '#73737d',fontsize=9)
        ax.set(title='Selected K points: '+['Gamma (0)','first ring (1-6)','outer ring (7-18)'][i],xlabel=r'$k_x/(2\pi/a)$',ylabel=r'$k_y/(2\pi/a)$',aspect='equal',xlim=(-.095,.095),ylim=(-.095,.095),xticks=[-.08,-.04,0,.04,.08],yticks=[-.08,-.04,0,.04,.08])
        bar=a.colorbar(fig,ax,dots,r'$|a_{k,\mathrm{py}}|^2/\sum_{19}|a_{k,\mathrm{py}}|^2$ (1)')
        bar.set_ticks([1e-5,1e-3,1e-1,1])
        a.figure_panel(fig,axs[i,1],z['near_EM_normalized'][i],x,f'{label} | Air field intensity',r'$u_{\mathrm{EM,syn}}/u_{\mathrm{EM,finite,max}}$ (1)',near_limits[i])
        axs[i,1].plot(boundary[:,0],boundary[:,1],color='white',lw=.65)
        axs[i,1].set_xlim(-18,18);axs[i,1].set_ylim(-18,18)
        im=axs[i,2].pcolormesh(edge,edge,np.maximum(z['far_EM_normalized'][i],far_lower),cmap='magma',norm=far_norm,shading='flat',rasterized=True)
        axs[i,2].set(title=f'{label} | Far field intensity',xlabel=r'$\theta_x$ (°)',ylabel=r'$\theta_y$ (°)',xlim=(-10,10),ylim=(-10,10),aspect='equal',xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10])
        a.colorbar(fig,axs[i,2],im,r'$\mathcal{U}_{\mathrm{EM,syn}}/\mathcal{U}_{\mathrm{EM,finite}}(0)$ (1)')
    fig.suptitle('Unit-cell synthesis | Actual cavity aperture, 1141 cells | Air plane z = 1.1 um\nCommon finite-field references; far fields share one logarithmic scale | No cladding feedback',fontsize=13)
    a.save(fig,'Fig.6_unitcell_Gamma_first_ring',dpi=a.FIGURE_DPI)
    checks.update(near_display_limits=near_limits,far_display_limits=[far_lower,far_upper],far_display_scale="log",lattice_phase_sum_abs_over_N=(abs(lattice_sum)/len(centers)).tolist())
    a.write_json('Fig6_unitcell_synthesis_checks.json',checks)
    table=weights[['index','kx_um_inv','ky_um_inv','weight_re','weight_im','frequency_THz','lattice_sum_abs_over_N','far_center_over_finite_center']].to_csv(index=False)
    report=table_links(f"""
## Fig.6：unit-cell Γ、第一圈与外围 K 点的空气场拼接

文件：[Fig.6_unitcell_Gamma_first_ring.png](../10_overview/Fig.6_unitcell_Gamma_first_ring.png)，../11_pdf/ 中的同名 PDF。三行三列，依次为纯 Γ-py、K1–K6 py 相干叠加、外围 K7–K18 py 相干叠加；第一列为对应离散 K 点示意，第二列为空气场完整 EM 强度，第三列为出射远场完整 EM 角谱。拼接范围由用户确认：只取 cavity 的 {len(centers)} 个实际晶胞，cladding 和外部置零。已有 Fig.1–5 不变。

**新增第一列：K 点示意。** 参考用户指定 04_k_3layers_broadening_2d.png 的第一个面板，显示相同 1+6+12 的 19 个实际离散 K 点及编号，坐标为 kx/(2π/a)、ky/(2π/a)，黑底、magma 对数色带 1e-5–1。第一排仅 0 号 Γ 着色，第二排仅 1–6 着色，第三排仅 7–18 着色，其余点用暗灰空心轮廓表示未参与当前重构。颜色使用本次重构真实复系数的 |a_k,py|²/sum_19 |a_k,py|²，三排分母相同，不单独重归一，也不复制参考图旧 Hz 权重。由于基底可非正交，此处为系数权重，不应直接称作可相加的物理能量占比；相干计算仍使用完整复系数。元数据见 80_logs/Fig6_kpoint_schematic.json。

**输入及补导出。** K 点编号与用户提供的 S5_mode19_c245-0.96-1.156_cl242-0.98-0.93_mesh5_20260920 分立 K 点一致；使用本工作区保留的 periodic_k00–k18.mph（cavity b=245 nm、eta=0.96、zeta=1.156、mesh5）。从已经求解的 py 本征模中补导出 z=1.1 μm 的全复 E/H，各 4881 个横向点。启动了本机 COMSOL 6.3 做插值导出，没有运行本征求解；模型文件未保存、修改或删除，导出前后 SHA256 一致。模型身份、原始 z=0 样本回查误差及模态频率核对见 80_logs/Fig6_unitcell_air_export.json。空气原始数据为 01_results/periodic_kNN_air1100nm_EH.npz。此处直接展示 z=1.1 μm 空气场，不混称为 z=100 nm 的 NA 滤波近场。

**权重来源。** 使用 field_decomposition.npz 中完整 E/H 能量投影的 coefficients.reshape(19,4)[:,1]，而非远场强度或 |Hz|² 权重。每个 unit-cell 模态沿用 periodic_kNN_EH.npz 保存的 phase_to_raw_reference 和完整单胞储能 U_cell；其复 E/H 除以 sqrt(U_cell) 后乘对应复系数 a_k，不做独立相位重置。Γ 组乘 a_0；第一圈使用1–6号六个 a_k，外围使用7–18号十二个 a_k，分别组内相干相加，保留幅度及相位差。py_band_index=1 是当前基底的内部索引，不把 mode 编号当作 K 编号。

**相位与拼接公式。** COMSOL 约定 exp(+iωt)，空间 Bloch 因子 exp(-ik·r)。原始单胞场 F_k(rho) 已包含内部 Bloch 相位，因此先令 u_k(rho)=F_k(rho)exp(+ik·rho)。在 r=R+rho 处构造 a_k u_k(rho)exp(-ik·r)/sqrt(U_cell)，等价于 a_k F_k(rho)exp(-ik·R)/sqrt(U_cell)。两个公式已数值核对，不能再对 F_k 同时乘 exp(-ik·rho) 而重复计相位。使用实际 bulk_points 晶胞中心及 unit_cell_corners 几何确定支持区，共 {checks['support_pixels']} 个空气网格像素。图中白线为实际 cavity–cladding 边界。

**插值和边界核查。** 原生边界点存在法向场迹值不一致；不把左右边界点简单平均当作高精度物理解。每个单胞的 160 个边界点仅从拼接插值中排除，原始数据完整保留；采用其余内部节点的周期延拓（中心及相邻六胞）线性插值。内部样本重建误差最大 {weights['periodic_node_reconstruction_relative_error'].max():.3e}。相对于原始边界采样，周期插值在 Γ 边界点的相对差为 {weights.iloc[0]['boundary_trace_relative_difference']:.3%}，其他点见 CSV；这是弱 Γ 信号的精度限制，不能由内部插值回代准确而宣称 Γ 微弱辐射已收敛。需要进一步量化该微弱信号时，应做 mesh 与单胞边界迹值收敛，不应直接赋予其机制意义。

**空气到远场。** 每个 K 点的四个切向场 Ex、Ey、Hx、Hy 在同一 1601² 空气网格拼接，cavity 外补零；正号空间 FFT 等效补零至 16001²。使用与现有图相同的 E/H 向上传播分离算子，在参考有限腔频率 {a.FREQ/1e12:.12g} THz 下投影出 A_E^+，令 A_Ez^+=-(kx A_Ex^++ky A_Ey^+)/kz、A_H^+=n×A_E^+/Z0。完整 EM 角谱 U=(epsilon0 sum|A_E|²+mu0 sum|A_H|²)/4。远场范围 ±10°，横纵为 asin(kx/k0)、asin(ky/k0)，不增加角度 Jacobian。此图沿用 Far field intensity 的角谱定义，不是 dP/dΩ。仅比较强度无需把共同的空气平面传播相位回溯到表面；本数据中的复角谱参考面仍为 z=1.1 μm，不能未经相位修正就与表面复角谱做相干相加。

**合并检查。** 外围 K7–K18 的组内相干重构误差为 {checks['outer_ring_coherent_reconstruction_error']:.3e}。K1–K6 先相加复空气场再 FFT，与逐 K 的复角谱相加一致，L2 相对误差 {checks['ring_coherent_reconstruction_error']:.3e}。每个 K 点还将 FFT 中心的出射切向 E 与原始拼接复 E/H 的直接面积和独立对照。没有把各 K 点强度直接相加。

**归一化与显示。** 三行中图共用实际有限 mode19 在同一 z=1.1 μm 空气平面的 max(u_EM)={float(z['near_reference_J_m3']):.12g} J/m³；三行右图共用 Fig.1 完整有限场在中心的 EM 角谱分母 {float(z['far_reference_EM_center']):.12g}。因此数值具有共同参照，不是把各组单独归一到 1。空气场三排线性色标上限分别为 {near_limits}，不能只按空气图颜色亮度跨排比较幅度。第三列三幅远场统一 LogNorm [{far_lower:g}, {far_upper:g}]，恰跨五个数量级，分母及色标完全相同，可以直接比较颜色。低于显示下限或精确为零的像素显示最低颜色，原始强度数据不截断；对数仅用于颜色映射，角度坐标仍为线性。空气视窗 ±18 μm 覆盖完整 cavity，远场 ±10°。

**定量结果。** 外围 K7–K18 合并中心为 {checks['far_centers'][2]:.9g}，离轴峰值为 {checks['far_peaks'][2]:.9g}（同一有限场参照）。 相对实际有限腔中心，纯 Γ 组中心强度为 {checks['far_centers'][0]:.9g}；第一圈合并中心为 {checks['far_centers'][1]:.9g}，但离轴峰值为 {checks['far_peaks'][1]:.9g}。纯 Γ 的保存本征解 Q≈{weights.iloc[0]['Q']:.9g}，不是数学上精确的 BIC；本图没有人为将其空气场设为零。

**为什么第一圈没有重建中心亮点。** 对完整晶胞拼接，A_k(K)=a_k/sqrt(U_cell) × [单胞 F_k(rho) 的 K 空间积分] × sum_R exp(i(K-k)·R)。在 K=0，第一圈各 K 点的晶胞相位和除以 N 的绝对值最大仅 {np.max(abs(lattice_sum[1:])/len(centers)):.3e}。所以每个非 Γ 分量的中心就被离散晶格相位和压低，并不只是六个大中心振幅相互抵消。矩形网格栅格化和插值会产生小残差，本次每个非 Γ 分量的中心强度约为 1e-12–1e-10，相干合并中心仍很小。该简单“相同单胞复制＋离散 K 权重＋cavity 硬孔径”构造不能解释实际有限腔中心强度 1。

**物理适用范围。** 这是一组受控的空气场空间拼接模型，未包含有限界面引起的场形变、反馈和新的散射响应，不是重新求解得到的有限腔模。各 unit-cell 本征频率约 {weights['frequency_THz'].min():.6f}–{weights['frequency_THz'].max():.6f} THz，与有限腔频率不完全相同；把它们作为同一频率的空间基底叠加属于冻结空间场近似，不是多频时间平均辐射。硬截断空气场本身产生孔径衍射，不能将拼接后的 Γ 非零远场解释为无限周期 Γ-BIC 自发辐射，也不能将本次未解释的差额直接归因某一微观机制。

**文件与复现。** 01_results/Fig6_unitcell_synthesis.npz 保存三组拼接复 E/H、掩膜、逐 K 及合并出射复角谱、EM 图和分母。80_logs/Fig6_unitcell_weights.csv 保存 K 编号、复权重、频率、Q、归一化、边界差异、晶格相位和及各 K 中心强度；Fig6_unitcell_synthesis_checks.json 保存闭合与显示参数。执行 `.venv\\Scripts\\python.exe -B -m scripts.analysis.analyze_dipolar_unitcell_synthesis` 重构并绘图；`--plot-only` 仅从缓存绘图，不启动 COMSOL。新 unit-cell 空气采样导出入口为 scripts/run_main/export_dipolar_unitcell_air.py。

K 点与权重审计表（CSV 原始数值）：
```csv
{table}```
""")
    (a.R/'unitcell_synthesis_report.md').write_text(report,encoding='utf8')
    p=a.R/'report.md';text=p.read_text(encoding='utf8').split('\n## Fig.6：unit-cell ')[0]
    cell_report=a.R/'unitcell_center_radiation_report.md'
    appendix='\n'+cell_report.read_text(encoding='utf8') if cell_report.exists() else ''
    sector_report=a.R/'sector_interference_report.md'
    if sector_report.exists():appendix+='\n'+sector_report.read_text(encoding='utf8')
    p.write_text(text.rstrip()+'\n'+report+appendix,encoding='utf8')
    print('Fig.6 PNG/PDF and detailed report saved',flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--plot-only',action='store_true')
    if not parser.parse_args().plot_only:run()
    plot_report()
