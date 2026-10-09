"""Six spatial apertures on the actual upward complex near field; no COMSOL."""
from __future__ import annotations
from comsol_workflow.output_paths import table_links, table_path
import json
import numpy as np
import pandas as pd
from scipy.constants import epsilon_0, mu_0
from scripts.analysis import analyze_dipolar_complete as a

NAMES = ["Upper right", "Top", "Upper left", "Lower left", "Bottom", "Lower right"]


def run(*, redraw=False, figure5_only=False):
    d=a.load("direct_radiative_fields.npz")
    x=d["x_um"]; E=d["E_up_surface"]; H=d["H_up_surface"]
    assert np.isfinite(E).all() and np.isfinite(H).all()
    cfg=json.loads((a.CONFIG/"config.json").read_text(encoding="utf8"))
    boundary=np.asarray(cfg["inner_boundary"])
    origin=(boundary.max(0)+boundary.min(0))/2
    xx,yy=np.meshgrid(x-origin[0],x-origin[1])
    angle=np.arctan2(yy,xx)
    centers=np.deg2rad(np.arange(30,360,60))
    distance=np.abs(np.angle(np.exp(1j*(angle[None]-centers[:,None,None]))))
    masks=(distance<=np.pi/6+1e-12).astype(float)
    masks[:,np.hypot(xx,yy)<1e-12]=1
    masks/=masks.sum(0)
    np.testing.assert_allclose(masks.sum(0),1,rtol=0,atol=1e-15)
    # Half weights at shared rays preserve mirror symmetry; center belongs equally to six sectors.
    for i in range(3):np.testing.assert_allclose(masks[i],masks[i+3,::-1,::-1],atol=1e-12)
    dx=float(x[1]-x[0])*1e-6
    size=2*a.FFT_SIZE-1
    dq=2*np.pi/(size*dx)
    half=int(np.ceil(np.sin(np.deg2rad(10))*a.K0/dq))+1
    q=np.arange(-half,half+1)*dq
    qx,qy=np.meshgrid(q/a.K0,q/a.K0)
    nz=np.sqrt(1-qx*qx-qy*qy)
    n=np.stack([qx,qy,nz],-1)
    nt=n[...,:2]

    def transform(mask):
        # Fourier-transform tangential E/H, then linear outgoing-wave projection.
        et=np.stack([a.padded_transform_roi(mask*E[...,j],dx,size,len(q)) for j in (0,1)],-1)
        ht=np.stack([a.padded_transform_roi(mask*H[...,j],dx,size,len(q)) for j in (0,1)],-1)
        p=np.stack([ht[...,1],-ht[...,0]],-1)
        v=a.Z0/nz[...,None]*(p-nt*np.sum(nt*p,-1)[...,None])
        out=np.empty((*nz.shape,3),complex)
        out[...,:2]=(et+v)/2
        out[...,2]=-np.sum(nt*out[...,:2],-1)/nz
        return out

    if redraw:
        cached=a.load("spatial_filter_six_sectors.npz")
        np.testing.assert_allclose(cached["masks"],masks,rtol=0,atol=0)
        np.testing.assert_allclose(cached["q_m_inv"],q,rtol=0,atol=0)
        spectra=cached["E_up_angular"];full=cached["merged_E_up_angular"][0]
    else:
        spectra=[]
        for i in range(6):
            print(f"FFT sector {i+1}/6: {NAMES[i]}",flush=True)
            spectra.append(transform(masks[i]))
        spectra=np.asarray(spectra)
        full=transform(np.ones_like(xx))
    amplitude_error=float(np.linalg.norm(spectra.sum(0)-full)/np.linalg.norm(full))
    assert amplitude_error<1e-10
    merged_masks=np.stack([np.ones_like(xx),masks[[1,4]].sum(0),masks[[0,2,3,5]].sum(0)])
    merged_e=np.stack([full,spectra[[1,4]].sum(0),spectra[[0,2,3,5]].sum(0)])
    np.testing.assert_allclose((merged_e[1]+merged_e[2])/np.max(abs(full)),full/np.max(abs(full)),rtol=1e-10,atol=1e-12)
    sh=np.cross(n,spectra)/a.Z0
    mh=np.cross(n,merged_e)/a.Z0
    energy=lambda e,h:(epsilon_0*np.sum(abs(e)**2,-1)+mu_0*np.sum(abs(h)**2,-1))/4
    near=energy(E,H);near_peak=float(near.max())
    raw_far=energy(spectra,sh);merged_far=energy(merged_e,mh)
    far_center=float(merged_far[0,half,half]);assert far_center>0
    far=raw_far/far_center;mf=merged_far/far_center
    far_max=float(np.ceil(max(far.max(),mf.max())*10)/10)
    sector_max=float(np.ceil(far.max()*100)/100)
    center_e=spectra[:,half,half];center_h=sh[:,half,half]
    gram=(epsilon_0*(center_e.conj()@center_e.T)+mu_0*(center_h.conj()@center_h.T))/4/far_center
    center_self=float(np.trace(gram).real)
    center_interference=float(gram.sum().real-center_self)
    np.testing.assert_allclose(gram.sum().real,1,rtol=1e-10)
    full_cropped=a.load("Fig1_full_BZ_angular_spectrum_2x.npz")
    old=full_cropped["E_up_surface_angular"]; mid=len(old)//2
    old=old[mid-half:mid+half+1,mid-half:mid+half+1]
    np.testing.assert_allclose(full_cropped["q_m_inv"][mid-half:mid+half+1],q,rtol=1e-10,atol=1e-8)
    reference_error=float(np.linalg.norm(full-old)/np.linalg.norm(old))
    reference_center_error=float(np.linalg.norm(full[half,half]-old[half,half])/np.linalg.norm(old[half,half]))
    rows=[]
    all_names=NAMES+["Full field","Top + bottom","Left/right four"]
    all_masks=np.concatenate([masks,merged_masks]);all_far=np.concatenate([far,mf])
    for name,mask,ff in zip(all_names,all_masks,all_far):
        rows.append(dict(region=name,near_surface_energy_fraction=float(np.sum(near*mask**2)/np.sum(near)),center_EM_relative_to_full=float(ff[half,half]),peak_EM_relative_to_full_center=float(ff.max())))
    pd.DataFrame(rows).to_csv(table_path(a.L, "spatial_filter_regions.csv"),index=False)
    pd.DataFrame([dict(first=NAMES[i],second=NAMES[j],center_interference_relative_to_full=float(2*gram[i,j].real)) for i in range(6) for j in range(i+1,6)]).to_csv(table_path(a.L, "spatial_filter_interference.csv"),index=False)
    np.savez_compressed(a.D/"spatial_filter_six_sectors.npz",x_um=x,masks=masks,q_m_inv=q,E_up_angular=spectra,H_up_angular=sh,merged_E_up_angular=merged_e,merged_H_up_angular=mh,near_EM=near,far_EM_normalized=far,merged_far_EM_normalized=mf,center_gram=gram,near_denominator=near_peak,far_denominator=far_center)
    checks=dict(source="direct_radiative_fields.npz:E_up_surface,H_up_surface",z_nm=100,sector_center_deg=np.arange(30,360,60).tolist(),origin_um=origin.tolist(),fft_size=size,source_grid=list(E.shape[:2]),angular_limits_deg=[-10,10],far_color_limits=[0,far_max],sector_far_color_limits=[0,sector_max],near_color_limits=[0,1],partition_sum_max_error=float(abs(masks.sum(0)-1).max()),coherent_reconstruction_relative_error=amplitude_error,full_vs_Fig1_complex_relative_error=reference_error,full_vs_Fig1_center_complex_relative_error=reference_center_error,center_self_sum=center_self,center_pair_interference_sum=center_interference,center_top_bottom_interference=float(1-mf[1,half,half]-mf[2,half,half]),mask_edge_rule="equal weights on shared rays and at center",interpretation="spatial-aperture decomposition, not independent material-source radiation",COMSOL_started=False)
    a.write_json("spatial_filter_checks.json",checks)
    a.plt.rcParams.update({"font.size":10,"xtick.direction":"in","ytick.direction":"in","pdf.fonttype":42})
    edges=np.rad2deg(np.arcsin(np.r_[q-dq/2,q[-1]+dq/2]/a.K0))
    boundary=np.vstack([boundary,boundary[0]])

    def pair(fig,axes,mask,ff,label,maximum):
        a.figure_panel(fig,axes[0],near*mask**2/near_peak,x,f"{label} | Near field intensity",r"$u_{\mathrm{EM}}/u_{\mathrm{EM,full,max}}$ (1)",1)
        axes[0].plot(boundary[:,0],boundary[:,1],color="white",lw=.5,alpha=.65)
        for theta in np.deg2rad(np.arange(0,360,60)):
            axes[0].plot([origin[0],origin[0]+100*np.cos(theta)],[origin[1],origin[1]+100*np.sin(theta)],color="#5782ac",lw=.45,alpha=.65)
        axes[0].set_xlim(x[0],x[-1]);axes[0].set_ylim(x[0],x[-1])
        im=axes[1].pcolormesh(edges,edges,ff,cmap="magma",vmin=0,vmax=maximum,shading="flat",rasterized=True)
        axes[1].set(title=f"{label} | Far field intensity",xlabel=r"$\theta_x$ (°)",ylabel=r"$\theta_y$ (°)",xlim=(-10,10),ylim=(-10,10),aspect="equal",xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10])
        a.colorbar(fig,axes[1],im,r"$\mathcal{U}_{\mathrm{EM}}/\mathcal{U}_{\mathrm{EM,full}}(0)$ (1)")

    # COMSOL-style phase-zero real vector; one common scale across all sectors.
    step=max(1,int(round(np.sin(np.deg2rad(1.5))*a.K0/dq)))
    ids=half+np.arange(-half//step+1,half//step+1)*step
    ids=ids[(ids>=0)&(ids<len(q))&(abs(q[np.clip(ids,0,len(q)-1)]/a.K0)<np.sin(np.deg2rad(9.5)))]
    arrow_angle=np.rad2deg(np.arcsin(q[ids]/a.K0))
    axx,ayy=np.meshgrid(arrow_angle,arrow_angle)
    vectors=np.concatenate([spectra,merged_e])[:,ids[:,None],ids[None,:],:2].real
    vector_peak=float(np.max(np.linalg.norm(vectors,axis=-1)))
    sector_vector_peak=float(np.max(np.linalg.norm(vectors[:6],axis=-1)))
    assert vector_peak>0
    a.write_json("Fig5_farfield_vectors.json",dict(source="spatial_filter_six_sectors.npz:E_up_angular,merged_E_up_angular",components=["real(Ex)","real(Ey)"],phase_deg=0,time_convention="exp(+i omega t)",reference_plane_z_nm=100,common_vector_peak=vector_peak,sector_vector_peak=sector_vector_peak,vector_scale_groups="first two far-field columns: six sectors share sector peak; last far-field column: retains previous common peak",max_arrow_length_plot_deg=1.15,grid_angle_deg=arrow_angle.tolist(),note="Phase-zero real angular-spectrum vectors; not polarization ellipse axes; no per-region or per-pixel phase rotation. No custom COMSOL arrow settings were found."))
    def overlay(ax,index):
        v=vectors[index]/(sector_vector_peak if index<6 else vector_peak)*1.15
        visible=np.linalg.norm(v,axis=-1)>1.15e-3
        from matplotlib.patheffects import Stroke,Normal
        arrows=ax.quiver(axx[visible],ayy[visible],v[...,0][visible],v[...,1][visible],angles="xy",scale_units="xy",scale=1,pivot="mid",color="white",width=.004,headwidth=3.4,headlength=4.4,minlength=0,zorder=5)
        arrows.set_path_effects([Stroke(linewidth=.4,foreground="#252525"),Normal()])
        np.testing.assert_allclose(arrows.U,v[...,0][visible])
        np.testing.assert_allclose(arrows.V,v[...,1][visible])

    fig,axs=a.plt.subplots(3,6,figsize=(33,15),layout="constrained")
    for slot,i in enumerate([1,4,2,5,0,3]):
        pair_axes=axs[slot//2,2*(slot%2):2*(slot%2)+2]
        pair(fig,pair_axes,masks[i],far[i],NAMES[i],sector_max)
        overlay(pair_axes[1],i)
    for i,name in enumerate(all_names[6:]):
        pair(fig,axs[i,4:6],merged_masks[i],mf[i],name,far_max)
        overlay(axs[i,5],i+6)
    from matplotlib.quiver import Quiver
    assert sum(isinstance(c,Quiver) for ax in axs.flat for c in ax.collections)==9
    fig.suptitle("Spatial filtering | Six sectors and coherent sums | Arrows: Re(Ex, Ey), phase = 0 deg",fontsize=15)
    a.save(fig,"Fig.5_six_sector_spatial_filter",dpi=a.FIGURE_DPI)
    report=table_links(f"""
## 空间滤波：Fig.5 六扇区与相干合并

**目标和数据起点。** 本分析从 Fig.1 的 Near field intensity: Full EM 对应的同一套复场出发：`01_results/direct_radiative_fields.npz` 中 z=100 nm 的 E_up_surface、H_up_surface，1601×1601 网格，已完成空气 E/H 方向分离、NA=1 筛选及含相位回溯。Full EM 图只用于显示，绝不对能量强度图做 FFT。原始 COMSOL 数据和 Fig.1–4 不修改，不启动 COMSOL。

**切分。** 使用实际 config.json 的 inner_boundary 中心与晶格六边形方向，按示意图沿 0°、60°、120° 三条直线切成六个 60° 扇区：右上、上、左上、左下、下、右下。每个扇区覆盖整个已采样平面，并非只选择窄界面条带，也不是沿单个孔轮廓截断。白线描出实际 cavity–cladding 晶胞边界，蓝线为分区线。共享射线像素权重各 1/2，原点各 1/6，以保持镜像对称；处处 sum M_i=1。局部复场 E_i=M_i E，H_i=M_i H。显示近场为 M_i² u_EM/max(u_EM,full)，所有近场共用 0–1 色标。共享线上的各区能量不能直接加回完整能量，必须保留交叉项；这不影响复场精确相加。

**FFT 与向上辐射。** 四个切向分量分别按 A(kx,ky)=dxdy sum E/H(x,y) exp(+i(kx x+ky y)) 计算，等效补零至 {size}²，保留覆盖 ±10° 的中心网格；没有二次回溯相位，因为输入已在 z=100 nm。空间硬截断会使每个分区不再独立满足自由空间 Maxwell 约束，故不能只把截断 E 的三个 FFT 分量当成独立物理远场。采用与 Fig.1 相同的线性 E/H 向上传播分离：P=(A_Hy,-A_Hx)，n_t=(kx,ky)/k0，s=kz/k0；A_Et^+=(A_Et+Z0/s [P-n_t(n_t·P)])/2，A_Ez^+=-n_t·A_Et^+/s，再以 A_H^+=n×A_E^+/Z0 重建完整磁场。此处窗口完全位于光锥内。相当于对分区的切向场/等效面流作线性辐射分解，并不把各区视为重新求解后的独立腔。

**远场物理量与归一化。** U_EM=(epsilon0 sum|A_E^+|²+mu0 sum|A_H^+|²)/4。所有六区、两种合并及完整场均除以本次有限近场窗口的完整场中心 U_EM,full(0)，不是分别按各区峰值归一化。Far field intensity 沿用 Fig.1 的 EM 角谱定义，不是 dP/dΩ。横纵坐标为 asin(kx/k0)、asin(ky/k0)，±10°，等比例，不引入角度 Jacobian；Fig.5 六区远场统一线性色标 0–{sector_max:g}，Fig.5 最后两列的完整场及合并远场统一 0–{far_max:g}，各自覆盖全部面板峰值；两类面板使用相同归一化分母，但显示上限不同。补零细化采样，不增加原始信息或物理孔径。

**Fig.5 每组含义。** 整图为三行六列：每两列构成一个近场/远场组。前四列三排分别为上/下、左上/右下、右上/左下，每组左为保留该扇区后的完整 EM 近场强度，右为该扇区的出射 EM 角谱。用于比较各条界面邻域对应的角分布、中心振幅与定向性；不能仅凭某区近场亮就认为其独立远场功率更大。文件 `Fig.5_six_sector_spatial_filter.png/pdf`。

**Fig.5 远场矢量。** 叠加同一 E_up_angular 的 (Re Ex, Re Ey)，时间约定 exp(+iωt)，统一显示相位 0°，参考面 z=100 nm；这是 COMSOL 频域场实部箭头的表示方式，但未读取到本模型自定义绘图相位设置，不能声称已逐项匹配其 GUI 配置。前两组远场列的六个单独扇区共用其六区采样实部矢量最大模长作为长度分母；最后一组完整场及两组合并结果保留原九图公共峰值分母。两类箭头尺度分别统一，不应跨组按箭头长度比较绝对振幅，最大采样箭头长度对应图坐标 1.15°，非物理角偏转；横纵方向分别对应 Ex、Ey，不对箭头乘坐标 Jacobian。不逐像素/分区旋转相位，不使用 |Ex|、|Ey|，不把箭头解释为偏振椭圆长轴。箭头在约 1.5° 间隔的原角谱网格抽样，实部矢量模长低于共同峰值的 0.1% 不画，避免零场箭头。白色箭头带深色细边；背景仍为完整 EM 强度，色标不变。参数保存于 80_logs/Fig5_farfield_vectors.json。

**Fig.5 最后两列每排含义。** 第一排是同一有限近场窗口的完整场重算基准；第二排先将上、下两个扇区复场相加后计算；第三排合并另外四区。计算上利用 FFT 和出射算子的线性，直接相加各区复角谱，完全等价于先合并复近场再 FFT，绝非相加强度。这三组合并结果已纳入 Fig.5 最后两列；独立 Fig.6 图文件已移除。

**中心干涉与核对。** 六区复角谱和相对完整场的 L2 误差为 {amplitude_error:.3e}。中心六区自项之和为 {center_self:.9g}，两两干涉总和为 {center_interference:.9g}，二者和为 1。上下合并中心强度为 {mf[1,half,half]:.9g}，左右四区为 {mf[2,half,half]:.9g}；两组间干涉为 {checks['center_top_bottom_interference']:.9g}，三者和为 1。CSV 保存各区近场面积能量比例、中心强度、峰值及每对中心干涉项；这些比例是表面量，不是三维储能或全 NA 辐射功率。

**有限视窗与解释边界。** 当前近场仅保存了有限的 1601² 实空间窗口。重新 FFT 相当于另加一个有限矩形孔径；不能假定它严格恢复 Fig.1 直接由空气场得到的角谱。与 Fig.1 同网格比较，本次完整场中心复振幅相对差为 {reference_center_error:.6%}，±10°窗口复场 L2 相对差为 {reference_error:.6%}。上述分区闭合是本次同窗基准的线性核对，不能替代这个独立差异检查。硬分区本身产生衍射展宽；结果回答“哪些表面区域的相干场贡献于该角分布”，不单独证明辐射微观起源就在界面。若需定量消除有限窗误差，应从原角谱重建更大的复近场视窗并做窗口收敛；本次不新增 COMSOL 计算。

**可追踪数据。** `01_results/spatial_filter_six_sectors.npz` 保存分区掩膜、原近场 EM、每区及合并复 E/H 角谱、归一化分母、强度和中心 Gram 矩阵；复近场由掩膜乘上述源 NPZ 精确复现。`80_logs/spatial_filter_regions.csv`、`spatial_filter_interference.csv`、`spatial_filter_checks.json` 保存指标与验证。复现命令：`.venv\\Scripts\\python.exe -B -m scripts.analysis.analyze_dipolar_spatial_filter`。
""")
    (a.R/"spatial_filter_report.md").write_text(report,encoding="utf8")
    report_path=a.R/"report.md"
    text=report_path.read_text(encoding="utf8").split("\n## 空间滤波：Fig.5 六扇区")[0]
    unitcell_report=a.R/"unitcell_synthesis_report.md"
    appendix="\n"+unitcell_report.read_text(encoding="utf8") if unitcell_report.exists() else ""
    cell_report=a.R/"unitcell_center_radiation_report.md"
    if cell_report.exists():appendix+="\n"+cell_report.read_text(encoding="utf8")
    sector_report=a.R/"sector_interference_report.md"
    if sector_report.exists():appendix+="\n"+sector_report.read_text(encoding="utf8")
    report_path.write_text(text.rstrip()+"\n"+report+appendix,encoding="utf8")
    print(json.dumps(checks,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--redraw",action="store_true");parser.add_argument("--figure5-only",action="store_true")
    args=parser.parse_args();run(redraw=args.redraw,figure5_only=args.figure5_only)
