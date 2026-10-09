"""Offline four-Figure delivery for the approved S4 dual-mode scan."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scripts.run_sweep.run_boundary_dual_validation import ROOT, ARCHIVE, s4, analyze, prepare
from comsol_workflow.boundary_integrals import hole_edges

GROUP = "09_px_py_validation"
NAMES = ["figure1_boundary_validation", "figure2_radiation_zero", "figure3_Hz_and_boundaries", "figure4_boundary_cancellation"]
COLORS = {"px":"#2476b8", "py":"#cf3e43"}


def formula(mode, hole=None, edge=False):
    normal, sign = ("y", "-") if mode == "px" else ("x", "")
    integral = r"\sum_{j=1}^{6}\oint_{\partial\Omega_j}" if hole is None else rf"\oint_{{\partial\Omega_{{{hole}}}}}"
    if edge: integral = r"\int_{\partial\Omega_{j,e}}"
    return rf"$\mathrm{{Re}}\left[{sign}\frac{{i}}{{\omega\varepsilon_0}}\left(\frac{{1}}{{\varepsilon_{{\rm air}}}}-\frac{{1}}{{\varepsilon_{{\rm slab}}}}\right){integral}n_{normal}H_z\,\mathrm{{d}}l\right]$"


def crossings(x, y, valid):
    """Never bridge missing/mixed samples or invent exact zeros by tolerance."""
    out = []
    for i in range(len(x)):
        if not valid[i] or not np.isfinite(y[i]):
            continue
        if y[i] == 0:
            out.append({"left":float(x[i]), "right":float(x[i]), "linear_estimate":float(x[i])})
        elif i+1 < len(x) and valid[i+1] and np.isfinite(y[i+1]) and y[i]*y[i+1] < 0:
            out.append({"left":float(x[i]), "right":float(x[i+1]),
                        "linear_estimate":float(x[i]-y[i]*(x[i+1]-x[i])/(y[i+1]-y[i]))})
    return out


def summarize(df):
    summary = {}
    for mode in ("px", "py"):
        d = df[df.user_mode == mode].sort_values("zeta").reset_index(drop=True)
        ok = d.mode_status.eq("single_eigenmode").to_numpy()
        if not ok.any():
            raise ValueError(f"No identified {mode} eigenstates")
        candidates = d[ok]
        minimum = candidates.loc[candidates.inverse_Q.idxmin()]
        i = int(minimum.name)
        near = d.iloc[max(0,i-1):min(len(d),i+2)]
        record = {"sample_count":len(d), "mixed_zeta":d.loc[~ok,"zeta"].tolist(),
            "q_minimum_zeta":float(minimum.zeta), "Q_maximum":float(minimum.Q),
            "q_minimum_bracket":near.zeta.tolist(), "representative_zeta":pd.concat([candidates[candidates.zeta<minimum.zeta].tail(1),candidates[candidates.zeta==minimum.zeta],candidates[candidates.zeta>minimum.zeta].head(1)]).zeta.tolist(),
            "air_complex_minimum_zeta":float(candidates.iloc[np.argmin(np.hypot(candidates.air_amplitude_V_per_m_re,candidates.air_amplitude_V_per_m_im))].zeta)}
        for key, col in (("electric", "electric_display_Vm"), ("boundary", "boundary_display_Vm"), ("air", "air_display_V_per_m")):
            y = d[col].to_numpy()
            record[key+"_crossings"] = crossings(d.zeta.to_numpy(), y, ok)
            unresolved=[]
            reliable=np.flatnonzero(ok)
            for left,right in zip(reliable[:-1],reliable[1:]):
                if right>left+1 and y[left]*y[right]<0:
                    unresolved.append({"left":float(d.zeta.iloc[left]),"right":float(d.zeta.iloc[right]),
                        "ambiguous_interior_zeta":d.zeta.iloc[left+1:right].tolist()})
            record[key+"_opposite_endpoints_across_unresolved_points"]=unresolved
            delta = np.diff(y)[ok[:-1] & ok[1:]]
            record[key+"_increments"] = {"positive":int(sum(delta>0)),"negative":int(sum(delta<0)),"zero":int(sum(delta==0))}
        record["full_coverage"] = bool(d.zeta.min()==.8 and d.zeta.max()==1.2 and len(d)>=50)
        record["identity_complete"] = bool(ok.all())
        record["air_monotonic"] = bool(all(np.diff(d.air_display_V_per_m)>=0) or all(np.diff(d.air_display_V_per_m)<=0))
        record["existence_uniqueness_supported"] = bool(record["full_coverage"] and ok.all() and record["air_monotonic"] and len(record["air_crossings"])==1)
        qlo,qhi = min(record["q_minimum_bracket"]),max(record["q_minimum_bracket"])
        record["local_position_agreement"] = bool(len(record["electric_crossings"])==len(record["boundary_crossings"])==1 and
            max(record["electric_crossings"][0]["left"],record["boundary_crossings"][0]["left"],qlo) <=
            min(record["electric_crossings"][0]["right"],record["boundary_crossings"][0]["right"],qhi))
        record["air_plane_difference_at_Q_max_V_per_m"] = float(np.hypot(minimum.air_amplitude_V_per_m_re-minimum.air2_amplitude_V_per_m_re,minimum.air_amplitude_V_per_m_im-minimum.air2_amplitude_V_per_m_im))
        record["air_complex_magnitude_at_Q_max_V_per_m"] = float(np.hypot(minimum.air_amplitude_V_per_m_re,minimum.air_amplitude_V_per_m_im))
        summary[mode] = record
    return summary


def decorate(ax, title, ylabel, full=True, zero=True):
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel(r"$\zeta$")
    ax.set_ylabel(ylabel)
    ax.grid(False)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if zero: ax.axhline(0,color=".55",lw=.6,zorder=0)
    if full: ax.set_xlim(.8,1.2)
    ax.spines[["top","right"]].set_visible(True)


def curve(ax, d, column, **kwargs):
    y = d[column].to_numpy(dtype=float).copy()
    mixed = d.mode_status.ne("single_eigenmode").to_numpy()
    y[mixed] = np.nan
    ax.plot(d.zeta,y,**kwargs)
    # Unassigned original eigenstates remain visible but are never connected as pure modes.
    if mixed.any(): ax.scatter(d.zeta[mixed],d[column][mixed],marker="x",color=".55",s=22,zorder=5)


def reference(ax, summary):
    z = summary["q_minimum_bracket"]
    if len(z)>1: ax.axvspan(min(z),max(z),color="#d7e0e7",alpha=.55,zorder=-1)


def save(fig, name, target):
    for category,ext in (("10_overview","png"),("11_pdf","pdf")):
        path=(target/category if ext=="png" else target/category/GROUP)/f"{name}.{ext}"
        path.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(path,dpi=220)
    plt.close(fig)


def figure3(df, edges, summary, target):
    """Pair input Hz, spatial edge contributions and their exact signed sum."""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    selected=[]
    for mode in ("px","py"):
        z=summary[mode]["q_minimum_zeta"]
        row=df[(df.user_mode==mode)&(df.zeta==z)].iloc[0]
        if row.mode_status!="single_eigenmode":raise ValueError("Figure 3 requires identified eigenstates")
        with np.load(row.field_file) as source:data=dict(source)
        hz=(data[f"{mode}_surface_Hz"]*complex(row.factor_re,row.factor_im)).real
        values=edges[(edges.user_mode==mode)&(edges.zeta==z)].sort_values(["hole","edge"])
        if list(zip(values.hole,values.edge))!=[(h,e) for h in range(1,7) for e in range(1,4)]:
            raise ValueError("Figure 3 needs each of the 18 edges exactly once")
        np.testing.assert_allclose(values.display_Vm.sum(),row.boundary_display_Vm,rtol=1e-10,atol=1e-24)
        selected.append((mode,z,row,data,hz,values.display_Vm.to_numpy()))
    scale=max(float(np.max(abs(hz))) for _,_,_,_,hz,_ in selected)
    fig,axes=plt.subplots(2,3,figsize=(16,10),gridspec_kw={"width_ratios":[1,1,1.35]})
    fig.subplots_adjust(left=.055,right=.975,bottom=.13,top=.91,wspace=.62,hspace=.72)
    fig.suptitle(r"$\eta=0.96,\quad z=100\,\mathrm{nm}$",fontsize=12,y=.98)
    limits={}; choices={}
    for i,(mode,z,row,data,hz,values) in enumerate(selected):
        field_ax,spatial,bars=axes[i]
        xy=data["area_xy_m"]*1e9; outer=data["outer_m"]*1e9
        geometry=hole_edges(data["holes_m"])
        segments=np.array([[e.start,e.end] for e in geometry])*1e9
        extent=np.max(abs(outer))*1.04
        for ax in (field_ax,spatial):
            ax.set(xlabel=r"$x\ (\mathrm{nm})$",ylabel=r"$y\ (\mathrm{nm})$",xlim=(-extent,extent),ylim=(-extent,extent),aspect="equal")
            ax.plot(*np.vstack([outer,outer[0]]).T,color=".82" if ax is spatial else ".65",lw=.6)
        mesh=field_ax.tripcolor(Triangulation(xy[:,0],xy[:,1]),hz/scale if scale else hz,cmap="RdBu_r",vmin=-1,vmax=1,shading="gouraud",rasterized=True)
        field_ax.add_collection(LineCollection(segments,colors=".3",linewidths=.6))
        field_ax.set_title(rf"$p_{{{mode[-1]}}},\ \zeta={z:g}$ : $\mathrm{{Re}}(H_z)$",fontsize=10,pad=12)
        field_bar=fig.colorbar(mesh,cax=field_ax.inset_axes([1.035,0,.045,1]),ticks=np.arange(-1,1.01,.25),format="%.2f")
        field_bar.set_label(r"$\mathrm{Re}(H_z)$ (normalized, 1)",fontsize=8)
        field_bar.ax.tick_params(direction="in",labelsize=7)
        maximum=float(np.max(abs(values)))
        ticks=MaxNLocator(nbins=4,symmetric=True).tick_values(-1.08*(maximum or 1.),1.08*(maximum or 1.))
        limit=float(ticks[-1]); limits[mode]=limit
        norm=Normalize(-limit,limit)
        collection=LineCollection(segments,cmap="RdBu_r",norm=norm,linewidths=2.0,capstyle="round",joinstyle="round")
        collection.set_array(values); spatial.add_collection(collection)
        for j,hole in enumerate(data["holes_m"]*1e9,1):
            field_ax.text(*hole.mean(axis=0),str(j),ha="center",va="center",fontsize=9,fontweight="bold",color=".15")
            spatial.text(*hole.mean(axis=0),str(j),ha="center",va="center",fontsize=8,fontweight="bold",color=".3")
        for edge,segment in zip(geometry,segments):
            midpoint=segment.mean(axis=0)
            spatial.text(*(midpoint-edge.normal*22),str(edge.edge),ha="center",va="center",fontsize=7,color=".4")
        spatial.set_title(formula(mode,"j",edge=True),fontsize=8,pad=12)
        formatter=ScalarFormatter(useMathText=True); formatter.set_powerlimits((0,0))
        boundary_bar=fig.colorbar(collection,cax=spatial.inset_axes([1.035,0,.045,1]),ticks=ticks,format=formatter)
        boundary_bar.set_label(r"$\mathrm{V\,m}$",fontsize=9)
        boundary_bar.ax.tick_params(direction="in",labelsize=7)
        boundary_bar.ax.yaxis.get_offset_text().set_fontsize(8)
        bars.bar(np.arange(18),values,color=plt.get_cmap("RdBu_r")(norm(values)),width=.8)
        bars.axhline(0,color=".55",lw=.7)
        bars.set(xticks=np.arange(18),xticklabels=[f"{e.hole}:{e.edge}" for e in geometry],xlim=(-.7,17.7),ylim=(-limit,limit),xlabel=r"$\partial\Omega_j$: edge 1, 2, 3",ylabel=r"$\mathrm{V\,m}$")
        bars.set_title(formula(mode,"j",edge=True),fontsize=8,pad=12)
        bars.set_yticks(ticks); bars.tick_params(axis="x",labelrotation=65,labelsize=7)
        bars.ticklabel_format(axis="y",style="sci",scilimits=(0,0),useMathText=True)
        bars.set_box_aspect(1/1.35)
        comp=mode[-1]
        electric_text="0" if row.electric_display_Vm==0 else formatter.format_data(float(f"{row.electric_display_Vm:.4g}"))
        sum_text="0" if values.sum()==0 else formatter.format_data(float(f"{values.sum():.4g}"))
        field_ax.text(.5,-.29,rf"$\mathrm{{Re}}\int_S E_{comp}\,\mathrm{{d}}A={electric_text}\,\mathrm{{V\,m}}$",transform=field_ax.transAxes,ha="center",va="top",fontsize=9)
        bars.text(.5,-.29,formula(mode)+"\n"+rf"$={sum_text}\,\mathrm{{V\,m}}$",transform=bars.transAxes,ha="center",va="top",fontsize=8)
        for ax in axes[i]:
            ax.grid(False); ax.tick_params(which="both",direction="in",top=True,right=True)
        choices[mode]={"zeta":float(z),"field_file":str(row.field_file),"electric_integral_display_Vm":float(row.electric_display_Vm),"edge_sum_display_Vm":float(values.sum()),"edge_abs_max_Vm":maximum}
    save(fig,NAMES[2],target)
    return {"Hz_common_max_A_per_m":scale,"boundary_display_sign":-1,"field_colormap":"RdBu_r","field_clim":[-1,1],"field_ticks":np.arange(-1,1.01,.25).tolist(),"figure3":{"layout":"rows:px,py; columns:Re(Hz),colored boundaries,18 edge contributions","boundary_display":"uniform colored lines on the actual hole edges; no area filling","boundary_color_width_pt":2.0,"boundary_join":"round line caps at actual vertices","boundary_capstyle":"round","boundary_shadow":False,"boundary_backing":False,"boundary_outline":False,"edge_labels":"outside hole, 22 nm offset","normal_arrows":False,"selected":choices,"boundary_symmetric_limits_Vm":limits,"boundary_scale":"same color and bar scale within each mode; explicit V*m units","reference":"90_history/10_overview/03_control_eta095/s4_edges_py.png"}}

def plotted_samples(df):
    """User-excluded Gamma degeneracy sample; keep original tables intact."""
    return df.loc[~np.isclose(df.zeta,1.0,rtol=0,atol=1e-12)].copy()


def trend_fit(x,y, *, zero_at=None):
    """Display-only quadratic least squares, optionally constrained at a given zero."""
    x,y=np.asarray(x,float),np.asarray(y,float)
    if x.shape!=y.shape or len(x)<3 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Trend fit requires at least three finite paired samples")
    if len(np.unique(x))!=len(x):raise ValueError("Repeated zeta in trend fit")
    if zero_at is not None:
        if not np.isfinite(zero_at):raise ValueError("Fit zero must be finite")
        offset=x-zero_at
        coefficients=np.linalg.lstsq(np.c_[offset,offset**2],y,rcond=None)[0]
        return np.polynomial.Polynomial([0.,*coefficients])(np.polynomial.Polynomial([-zero_at,1.]))
    return np.polynomial.Polynomial.fit(x,y,2).convert()


def figure1(df, summary, target):
    df=plotted_samples(df)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.linewidth":.7,"lines.linewidth":1.3,"savefig.facecolor":"white"})
    fits={}; centers={}
    band_width=.002  # Common visual guide width, not a zero uncertainty interval.
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout="constrained")
    for j,mode in enumerate(("px","py")):
        d=df[df.user_mode==mode].sort_values("zeta")
        comp="x" if mode=="px" else "y"
        label=rf"$\mathrm{{Re}}\int_S E_{comp}\,\mathrm{{d}}A$"
        x=d.zeta.to_numpy(float)
        boundary=trend_fit(x,d.boundary_display_Vm.to_numpy(float))
        roots=[float(r.real) for r in boundary.roots() if abs(r.imag)<1e-10 and x.min()<=r.real<=x.max()]
        if not roots:raise ValueError(f"No in-range boundary fit zero for {mode}")
        center=min(roots,key=lambda r:abs(r-summary[mode]["q_minimum_zeta"]))
        centers[mode]=center
        fitted={"boundary_display_Vm":boundary,
                "electric_display_Vm":trend_fit(x,d.electric_display_Vm.to_numpy(float),zero_at=center)}
        for column,model in fitted.items():
            fits[mode+":"+column]={"sample_count":len(x),"coefficients_ascending_zeta":model.coef.tolist(),
                "zero_constraint":center if column=="electric_display_Vm" else None}
        # Include the exact shared crossing in both rendered lines.
        smooth_x=np.sort(np.r_[np.linspace(x.min(),x.max(),801),center])
        for i in range(2):
            ax=axes[i,j]
            for column,color,style,legend in (("electric_display_Vm","#242424","-",label),("boundary_display_Vm","#cc3737","--",formula(mode))):
                ax.plot(x,d[column],linestyle="none",marker="o",ms=2.8,color=color,zorder=3)
                ax.plot(smooth_x,fitted[column](smooth_x),color=color,ls=style,label=legend)
            decorate(ax,rf"$p_{comp}$ : "+label,r"$\mathrm{V\,m}$")
            ax.axvspan(center-band_width/2,center+band_width/2,color="#d7e0e7",alpha=.55,zorder=-1)
            if i==1:
                lo,hi=center-.004,center+.004
                ax.set_xlim(max(.8,lo),min(1.2,hi))
                subset=d[d.zeta.between(lo,hi)]
                vals=np.r_[subset[["electric_display_Vm","boundary_display_Vm"]].to_numpy().ravel(),*[f(smooth_x[(smooth_x>=lo)&(smooth_x<=hi)]) for f in fitted.values()]]
                if vals.size:
                    a,b=min(0,np.nanmin(vals)),max(0,np.nanmax(vals)); pad=max(b-a,abs(a),1e-22)*.12
                    ax.set_ylim(a-pad,b+pad)
            bottom,top=ax.get_ylim(); ax.set_ylim(bottom,top+.32*(top-bottom))
            ax.legend(loc="upper right",fontsize=8,frameon=False)
    save(fig,NAMES[0],target)
    return {"figure1":{"figure_aspect":"4:3","figure_size_inches":[12,9],"marker":"o","fit":"boundary: quadratic least squares; electric: quadratic least squares constrained to boundary fit zero","excluded_zeta":[1.0],"shared_display_zero":centers,"highlight_width_zeta":band_width,"highlight_is_uncertainty":False,"fit_scope":"one global fit per mode and calculation; identical fit in full range and zoom","display_only":True,"zero_acceptance_unchanged":True,"fits":fits}}


def positive_q_fit(x,y):
    """Positive display trend; fitted minimum is not a measured radiation zero."""
    from scipy.optimize import least_squares
    x,y=np.asarray(x,float),np.asarray(y,float)
    trend_fit(x,y)  # Reuse paired-sample validation.
    if np.any(y<=0):raise ValueError("Inverse Q must be positive for log fitting")
    center=x[np.argmin(y)]
    curvature=max(float(y.max()/np.ptp(x)**2),np.finfo(float).tiny)
    def residual(p):
        return np.log(np.exp(p[0])+np.exp(p[1])*(x-p[2])**2)-np.log(y)
    result=least_squares(residual,[np.log(y.min()),np.log(curvature),center],bounds=([-700,-700,x.min()],[700,700,x.max()]))
    if not result.success:raise RuntimeError(result.message)
    return np.array([np.exp(result.x[0]),np.exp(result.x[1]),result.x[2]])


def anchored_q_fit(x,y,split,peak_Q):
    """Fit each side with the same measured peak; allow a corner at the join."""
    from scipy.optimize import least_squares
    x,y=np.asarray(x,float),np.asarray(y,float)
    trend_fit(x,y)
    if peak_Q<=0 or np.any(y<=0):raise ValueError("Q and inverse Q must be positive")
    distance=np.abs(x-split); scale=float(y.max()); floor=1/peak_Q
    design=np.c_[distance,distance**2]
    initial=np.maximum(np.linalg.lstsq(design,(y-floor)/scale,rcond=None)[0],1e-10)
    result=least_squares(lambda p:np.log(floor/scale+design@p)-np.log(y/scale),initial,bounds=(0,np.inf))
    if not result.success:raise RuntimeError(result.message)
    return result.x*scale


def refined_px(df, mesh_data):
    """Keep one mesh family; convert its per-zeta reference to the scan phase."""
    d=mesh_data.loc[(mesh_data.mesh==5)&(mesh_data.size_factor==.25)].copy()
    d=d.merge(df.loc[df.user_mode=="px",["zeta","factor_re","factor_im"]],on="zeta",how="left",validate="one_to_one").sort_values("zeta")
    factor=d.factor_re.to_numpy()+1j*d.factor_im.to_numpy()
    air=d.air_Ex_re.to_numpy()+1j*d.air_Ex_im.to_numpy()
    if len(d)<3 or not np.all(np.isfinite(factor)&(np.abs(factor)>0)&np.isfinite(air)&np.isfinite(d.Q)&(d.Q>0)):
        raise ValueError("Refined px needs finite samples and a matching nonzero scan phase at every zeta")
    # Mesh exports already have Hz RMS=1 A/m: apply phase only, never renormalize.
    d["air_display_V_per_m"]=(air*factor/np.abs(factor)).real
    return d


def figure2(df, summary, target, *, mesh_comparison=True):
    df=plotted_samples(df)
    fits={}
    mesh_path=target/"80_logs/10_px_mesh_test/mesh_results.csv"
    fine=refined_px(df,plotted_samples(pd.read_csv(mesh_path))) if mesh_comparison and mesh_path.exists() else None
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.linewidth":.7,"lines.linewidth":1.3,"savefig.facecolor":"white"})
    fig,axes=plt.subplots(2,1,figsize=(9,9),layout="constrained")
    for mode in ("px","py"):
        d=df[df.user_mode==mode].sort_values("zeta"); comp="x" if mode=="px" else "y"
        x=d.zeta.to_numpy(float)
        split=float(summary[mode]["q_minimum_zeta"])
        for ax,column in zip(axes,("air_display_V_per_m","Q")):
            ax.plot(x,d[column],linestyle="none",marker="o",ms=2.8,color=COLORS[mode],zorder=3)
        anchor=d.loc[d.zeta==split].iloc[0]
        fits[mode]={"shared_peak_Q":float(anchor.Q),"shared_air_value":float(anchor.air_display_V_per_m),"split_zeta":split,"split_source":"sampled Q maximum near radiation zero","segments":{}}
        for side,mask in (("left",x<=split),("right",x>=split)):
            part=d.loc[mask]; segment_x=part.zeta.to_numpy(float)
            smooth_x=np.linspace(segment_x.min(),segment_x.max(),801)
            offset=segment_x-split
            coefficients=np.linalg.lstsq(np.c_[offset,offset**2],part.air_display_V_per_m-anchor.air_display_V_per_m,rcond=None)[0]
            air=np.polynomial.Polynomial([float(anchor.air_display_V_per_m),*coefficients])
            linear,quadratic=anchored_q_fit(segment_x,part.inverse_Q,split,float(anchor.Q))
            distance=np.abs(smooth_x-split)
            suffix=r", $h/h_0=1$" if fine is not None else ""
            axes[0].plot(smooth_x,air(smooth_x-split),color=COLORS[mode],label=rf"$p_{comp}$"+suffix if side=="left" else "_nolegend_")
            axes[1].plot(smooth_x,1/(1/anchor.Q+linear*distance+quadratic*distance**2),color=COLORS[mode],label=rf"$p_{comp}$"+suffix if side=="left" else "_nolegend_")
            fits[mode]["segments"][side]={"sample_count":len(part),"air_coefficients_ascending_zeta_offset":air.coef.tolist(),"inverse_Q_linear_quadratic_distance":[float(linear),float(quadratic)]}
        for ax in axes:ax.axvline(split,color=COLORS[mode],ls="--",lw=.8,alpha=.65,zorder=0)
    if fine is not None:
        anchor=fine.loc[fine.Q.idxmax()]
        # Local refinement has no measured left branch; do not extrapolate one.
        if anchor.zeta!=fine.zeta.min():raise ValueError("Review the fine-mesh fit when samples extend left of its Q maximum")
        x=fine.zeta.to_numpy(); smooth_x=np.linspace(x.min(),x.max(),801)
        air=trend_fit(x,fine.air_display_V_per_m)
        linear,quadratic=anchored_q_fit(x,1/fine.Q,float(anchor.zeta),float(anchor.Q))
        distance=smooth_x-anchor.zeta
        qfit=1/(1/anchor.Q+linear*distance+quadratic*distance**2)
        for ax,column,fit_y in zip(axes,("air_display_V_per_m","Q"),(air(smooth_x),qfit)):
            ax.plot(smooth_x,fit_y,color=COLORS["px"],ls="--",lw=1.6,label=r"$p_x$, $h/h_0=0.25$",zorder=4)
            ax.plot(x,fine[column],ls="none",marker="o",ms=4.5,mfc="white",mec=COLORS["px"],zorder=5)
        inset=axes[0].inset_axes([.32,.11,.31,.20])
        coarse=df[(df.user_mode=="px")&df.zeta.between(x.min(),x.max())]
        inset.axhline(0,color=".65",lw=.6)
        inset.plot(coarse.zeta,coarse.air_display_V_per_m,ls="none",marker="o",ms=2.4,color=COLORS["px"])
        inset.plot(smooth_x,air(smooth_x),ls="--",color=COLORS["px"])
        inset.plot(x,fine.air_display_V_per_m,ls="none",marker="o",ms=4,mfc="white",mec=COLORS["px"])
        inset.set(xlabel=r"$\zeta$",ylabel=r"$\mathrm{V/m}$",xticks=[x.min(),x.max()])
        inset.tick_params(which="both",direction="in",labelsize=8)
        inset.xaxis.label.set_size(8);inset.yaxis.label.set_size(8)
        fits["px_refined"]={"source":str(mesh_path),"mesh":5,"size_factor":.25,"sample_count":len(fine),"zeta_range":[float(x.min()),float(x.max())],"phase":"multiply air_Ex by scan factor / abs(scan factor); no amplitude rescaling","air_coefficients_ascending_zeta":air.convert().coef.tolist(),"inverse_Q_linear_quadratic_distance":[float(linear),float(quadratic)],"anchor_zeta":float(anchor.zeta),"anchor_Q":float(anchor.Q),"fit_scope":"sampled right side only; no extrapolation or blending with coarse mesh","above_Q_axis_limit":fine.loc[fine.Q>1e8,["zeta","Q"]].to_dict("records"),"convergence_established":False}
    decorate(axes[0],r"$\mathrm{Re}\left[\frac{e^{-ik_0(z_a-z_s)}}{|S|}\int_{S(z_a)}E_{x,y}\,\mathrm{d}A\right]$",r"$\mathrm{V/m}$")
    limit=max(abs(v) for v in axes[0].get_ylim()); axes[0].set_ylim(-limit,limit)
    axes[0].legend(loc="upper right",fontsize=9,frameon=False)
    decorate(axes[1],r"$Q$",r"$Q\ (1)$",zero=False)
    axes[1].set_yscale("log")
    bottom,top=axes[1].get_ylim(); axes[1].set_ylim(bottom,1e8)
    axes[1].legend(loc="upper right",frameon=False,fontsize=9,**({"bbox_to_anchor":(1,1.18),"borderaxespad":0} if fine is not None else {}))
    save(fig,NAMES[1],target)
    return {"figure2":{"layout":"2 rows, 1 column; px and py overlaid in each panel","lower_quantity":"Q","Q_display":"raw Q samples; reciprocal of fitted positive inverse_Q model","figure_aspect":"1:1","figure_size_inches":[9,9],"marker":"o","air_fit":"coarse mesh: independent left/right quadratic least squares; refined mesh: local quadratic on sampled range only","join":"shared measured endpoint at split; value continuity only, no slope matching","Q_axis_max":1e8,"inverse_Q_fit":"1/Q_peak + a*abs(zeta-split) + b*(zeta-split)^2; a,b nonnegative, independently fitted per side in log residuals","excluded_zeta":[1.0],"position_guides":"vertical dashed lines at sampled Q maxima","air_axis_zero_centered":True,"Q_peak_annotations":False,"display_only":True,"zero_acceptance_unchanged":True,"fits":fits}}


def figure4(df, edges, summary, target):
    df=plotted_samples(df)
    fig,axes=plt.subplots(2,2,figsize=(13.6,8.6),layout="constrained")
    holecolors=plt.get_cmap("tab10").colors[:6]
    for j,mode in enumerate(("px","py")):
        d=df[df.user_mode==mode].sort_values("zeta")
        ed=edges[edges.user_mode==mode]
        holes=ed.groupby(["zeta","hole"],as_index=False).display_Vm.sum()
        for hole in range(1,7):
            h=holes[holes.hole==hole].merge(d[["zeta","mode_status"]],on="zeta")
            curve(axes[0,j],h,"display_Vm",color=holecolors[hole-1],label=rf"$\partial\Omega_{{{hole}}}$")
        curve(axes[0,j],d,"boundary_display_Vm",color="black",lw=1.8,label=r"$\sum_{j=1}^{6}$")
        decorate(axes[0,j],rf"$p_{{{mode[-1]}}}$ : "+formula(mode,"j"),r"$\mathrm{V\,m}$")
        bottom,top=axes[0,j].get_ylim(); axes[0,j].set_ylim(bottom,top+.32*(top-bottom))
        axes[0,j].legend(loc="upper right",fontsize=8,ncol=2,frameon=False)
        reference(axes[0,j],summary[mode])
        z=summary[mode]["q_minimum_zeta"]; near=ed[ed.zeta==z].sort_values(["hole","edge"])
        axes[1,j].bar(np.arange(18),near.display_Vm,color=[holecolors[int(h)-1] for h in near.hole],width=.75)
        decorate(axes[1,j],rf"$p_{{{mode[-1]}}},\ \zeta={z:g}$ : "+formula(mode,"j",edge=True),r"$\mathrm{V\,m}$",full=False)
        axes[1,j].set_xticks(np.arange(18),[f"{int(h)}.{int(e)}" for h,e in zip(near.hole,near.edge)],rotation=90,fontsize=8)
        axes[1,j].set_xlabel(r"$\partial\Omega_j$: edge 1, 2, 3")
    save(fig,NAMES[3],target)
    return {"figure4":{"excluded_zeta":[1.0]}}


def figures(df, edges, summary, target, *, mesh_comparison=True):
    return {**figure1(df,summary,target),
            **figure2(df,summary,target,mesh_comparison=mesh_comparison),
            **figure3(df,edges,summary,target),
            **figure4(df,edges,summary,target),"filenames":NAMES}






def validate_saved(config,df,edges,root=ROOT):
    from scripts.run_sweep.run_boundary_dual_validation import boundary_values
    from comsol_workflow.boundary_integrals import line_rule, polygon_area
    supplemental=[]
    for case in config["cases"]:
        folder=Path(case["directory"]); export=s4.json_read(folder/"export.json")
        assert export["status"]=="complete" and export["mesh"]["auto_size"]==5
        for name,digest in export["output_hashes"].items():assert s4.digest(folder/name)==digest
        data=np.load(folder/"fields.npz")
        assert data["edge_xy_m"].shape==(18,128,2)
        for k,edge in enumerate(hole_edges(data["holes_m"])):
            xy,w=line_rule(edge,128)
            np.testing.assert_allclose(data["edge_xy_m"][k],xy,atol=1e-20)
            np.testing.assert_allclose(data["edge_weights_m"][k],w,atol=1e-24)
            np.testing.assert_allclose(data["normals"][k],edge.normal,atol=1e-14)
        raw=pd.read_csv(folder/"mode_results.csv")
        assert len(raw)==2 and raw.mode_idx.nunique()==2
        for mode in ("px","py"):
            row=df[(df.zeta==case["zeta"])&(df.user_mode==mode)].iloc[0]
            r=raw[raw.user_mode==mode].iloc[0]
            factor=complex(row.factor_re,row.factor_im)
            omega=2*np.pi*1e12*complex(r.frequency_thz,-r.damping_thz)
            recalculated=boundary_values(data[f"{mode}_edge_Hz"],data["edge_weights_m"],data["normals"],omega)
            np.testing.assert_allclose(recalculated,data[f"{mode}_edge_converted_Vm"],rtol=1e-12,atol=1e-24)
            w=data["area_weights_m2"]
            np.testing.assert_allclose(np.dot(w,abs(factor*data[f"{mode}_surface_Hz"])**2)/w.sum(),1,rtol=1e-12)
            np.testing.assert_allclose(r.area_m2,abs(polygon_area(data["outer_m"])),rtol=1e-12)
            ei=edges[(edges.zeta==case["zeta"])&(edges.user_mode==mode)]
            assert len(ei)==18 and list(ei.groupby("hole").size())==[3]*6
            np.testing.assert_allclose(ei.display_Vm.sum(),row.boundary_display_Vm,rtol=1e-10,atol=1e-24)
            values={"zeta":case["zeta"],"user_mode":mode,"mode_status":row.mode_status,"Q":row.Q}
            for component in ("Ex","Ey"):
                for key in (component+"_raw_Vm","boundary_"+component+"_raw_Vm","air_"+component+"_raw_V_per_m","air2_"+component+"_raw_V_per_m"):
                    value=complex(r[key+"_re"],r[key+"_im"])*factor
                    assert np.isfinite(value)
                    values[key.replace("_raw", "_common_phase")]=value
            supplemental.append(values)
    assert s4.digest(s4.ROOT/"scripts/parameter.json")==config["shared_parameter_sha256"]
    for path,digest in config["source_hashes"].items():assert s4.digest(path)==digest
    result={"status":"passed","geometries":len(config["cases"]),"mode_rows":len(df),"edges_per_mode":18,
            "gauss_points_per_edge":128,"raw_complex_fields_and_phase_recomputed":True,
            "shared_parameter_unchanged":True,"historical_source_hashes_verified":len(config["source_hashes"]),
            "mode_thresholds":{"center_parity":.1,"surface_parity":.15},
            "phase_reference":{"px_zeta":.9,"py_fields":config["reference_fields"]},
            "display_boundary_sign":-1,"not_a_scientific_acceptance_flag":True}
    s4.write_csv(root/"80_logs/supplemental_complex_components.csv",supplemental)
    s4.write_json(root/"80_logs/verification.json",result)
    return result

def deliver(config, df, edges, summary):
    """Archive only this run; source MPH files in groups 01-08 are never moved."""
    import shutil
    records=[]
    manifest=ARCHIVE/"80_logs"/GROUP/"delivery_manifest.json"
    if manifest.exists():
        for row in s4.json_read(manifest)["files"]:
            assert s4.digest(row["destination"])==row["sha256"]
        delivered=pd.read_csv(manifest.parent/"scan_results.csv")
        assert len(delivered)==len(df) and set(delivered.zeta)==set(df.zeta)
        return delivered
    def transfer(source,target,move=False):
        source,target=Path(source).resolve(),Path(target).resolve()
        if not source.is_relative_to(ROOT.resolve()) or not target.is_relative_to((ARCHIVE/target.relative_to(ARCHIVE).parts[0]/GROUP).resolve()):
            raise ValueError("Delivery path outside this task")
        digest=s4.digest(source); size=source.stat().st_size
        if target.exists() and s4.digest(target)!=digest:
            raise ValueError(f"Conflicting delivery file: {target}")
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():
            if move:source.rename(target)
            else:shutil.copy2(source,target)
        assert s4.digest(target)==digest and target.stat().st_size==size
        records.append({"source":str(source),"destination":str(target),"sha256":digest,"bytes":size,"operation":"move" if move else "copy"})
        return target
    for case in config["cases"]:
        folder=Path(case["directory"]); name=folder.name
        model=folder/"s4_gamma.mph"
        if model.exists():
            target=transfer(model,ARCHIVE/"00_model"/GROUP/name/"s4_gamma.mph",move=True)
            case["source_model"]=str(target); case["source_model_sha256"]=s4.digest(target)
            case["new_model_relocated"]=True
            df.loc[df.zeta==case["zeta"],"source_model"]=str(target)
        field=transfer(folder/"fields.npz",ARCHIVE/"01_results"/GROUP/name/"fields.npz")
        df.loc[df.zeta==case["zeta"],"field_file"]=str(field)
        for source in folder.glob("*.csv"):
            transfer(source,ARCHIVE/"80_logs"/GROUP/name/source.name)
        for source in (folder/"01_results/eigenmodes").glob("*.parquet"):
            transfer(source,ARCHIVE/"01_results"/GROUP/name/"eigenmodes"/source.name)
        for source in (folder/"99_config").glob("*.json"):
            transfer(source,ARCHIVE/"99_config"/GROUP/name/source.name)
        transfer(folder/"export.json",ARCHIVE/"99_config"/GROUP/name/"export.json")
    s4.write_json(ROOT/"99_config/config.json",config)
    for source in (ROOT/"99_config").glob("*.json"):
        target=ARCHIVE/"99_config"/GROUP/source.name
        # The current run config is an active index, unlike immutable source snapshots.
        s4.write_json(target,s4.json_read(source))
    logs=ARCHIVE/"80_logs"/GROUP; logs.mkdir(parents=True,exist_ok=True)
    for source in (ROOT/"80_logs/supplemental_complex_components.csv",ROOT/"80_logs/verification.json"):
        transfer(source,logs/source.name)
    df["boundary_display_sign"]=-1
    edges=edges.merge(df[["zeta","user_mode","mode_idx","mode_status","field_file","source_model","factor_re","factor_im"]],on=["zeta","user_mode"],validate="many_to_one")
    edges["boundary_display_sign"]=-1
    df.to_csv(logs/"scan_results.csv",index=False)
    edges.to_csv(logs/"edge_results.csv",index=False)
    edges.groupby(["zeta","user_mode","hole"],as_index=False)[["boundary_Vm_re","boundary_Vm_im","display_Vm"]].sum().to_csv(logs/"hole_results.csv",index=False)
    s4.write_json(logs/"scientific_summary.json",summary)
    s4.write_json(logs/"delivery_manifest.json",{"files":records,"file_count":len(records),"total_bytes":sum(r["bytes"] for r in records),"verified":True})
    s4.write_json(ARCHIVE/"99_config"/GROUP/"execution.json",s4.json_read(ROOT/"execution.json"))
    for source in (Path(__file__),s4.ROOT/"scripts/run_sweep/run_boundary_dual_validation.py"):
        target=ARCHIVE/"99_config"/GROUP/"code"/source.name
        target.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(source,target)
    return df



def write_report(config, df, summary, audit):
    def brackets(items,unresolved=()):
        values=[f"[{r['left']:.4f}, {r['right']:.4f}]" for r in items]
        values += [f"[{r['left']:.4f}, {r['right']:.4f}]（中间模式待核验）" for r in unresolved]
        return "; ".join(values) or "未检出异号的可靠端点"
    lines=["# S4：双模式孔洞边界提取与法向辐射", "",
        "四图显示，两支模式的电场积分与边界加和具有同向的整体变化趋势，并在预期参数侧出现过零现象；局部边界曲线及损耗存在起伏，精确位置一致性与全区间唯一性须分别审查。", "",
        "本报告区分计算交付完成与物理假设得到支持；以下判断均来自真实样本，不强制零点一致，也不以幅值或斜率大小吻合作为验收。", "",
        f"固定 η=0.96，ζ=0.8–1.2，mesh=5，z=100 nm（板层侧），每条孔边128点。共 {len(config['cases'])} 个不同几何、{len(df)} 条模式记录；复用30份旧模型，其余为本轮真实Γ点求解。", "",
        "| 模式 | 电场面积分实部变号区间 | 边界加和显示值变号区间 | 空气面振幅实部变号区间 | 最低1/Q的采样ζ |",
        "|---|---|---|---|---|"]
    for mode in ("px","py"):
        s=summary[mode]
        lines.append(f"| {mode} | {brackets(s['electric_crossings'],s['electric_opposite_endpoints_across_unresolved_points'])} | {brackets(s['boundary_crossings'],s['boundary_opposite_endpoints_across_unresolved_points'])} | {brackets(s['air_crossings'],s['air_opposite_endpoints_across_unresolved_points'])} | {s['q_minimum_zeta']:.4f} |")
    lines += ["", "## 如何读四张图", "",
        "1. **面积与边界比较**：整体画布4:3，左列px的Ex、右列py的Ey；上排全范围，下排零点附近。全部实际计算值统一用圆点显示，黑色为电场面积分实部，红色为孔边界换算结果整体取负后的实部。ζ=1从绘图和拟合中排除，原始记录保留；红线为剩余样本的二次最小二乘全区间拟合，黑线约束经过对应红线的拟合零点。共同零点是显示约束，不是独立计算零点吻合的证据；下排显示同一拟合曲线。拟合仅展示整体趋势，不用于零点存在、位置或唯一性的科学验收。原始数值及模式状态保留。淡色矩形统一宽度Δζ=0.002，中心位于共同拟合零点，仅作视觉定位。",
        "2. **实际模式辐射**：上排从均匀空气区独立积分并传播补偿后的主偏振振幅实部；下排同一个原始本征解的Q。功率、Q本身无正负，不能替代上排的有符号振幅。",
        "3. **Hz、边界位置与逐边贡献**：参照03_control_eta095/s4_edges_py.png的三列形式，px、py各一行，分别取当前Q最大采样ζ=0.8393、1.156。左列为统一相位后的Re(Hz)，下方列独立电场面积分；中列将每条边的实际换算显示贡献直接着色到孔边界，标注孔号、边号及介质指向空气的法向；右列用相同颜色画18条边的贡献，并列出实际加和值。两幅Hz场共用一次归一化，色标为[-1,1]；每行边界颜色与柱状图共用该模式的对称V·m范围，两个模式的量级分别标明。没有对逐边贡献除以近零总和。Re(Hz)颜色不能直接等同于边界贡献的正负，因为后者还包含法向和复系数i/ω。空间三角插值仅用于Hz场图显示，不参与积分或生成新的ζ样本。",
        "4. **边界贡献的相消**：上排分别画六孔完整贡献及总和，下排画损耗极小值采样处全部18条边。每孔三边各计一次，没有按期望结果合并孔组。", "",
        "图4中的灰色叉号保留未通过模式对称性检查的原始本征场数据；这些点不与纯px/py曲线连接，也不参与可靠过零计数。图2整体画布9×9英寸，上下两图分别叠加px/py；上下两图在各模式Q最大值对应ζ左右分段拟合，共享分界处实测值，仅函数值连续、不强制斜率平滑；上排辐射实部二次拟合，下排展示原始Q和正值1/Q分段拟合的倒数（对数坐标，上限10^8）；拟合只作趋势展示。图1和图2按用户要求统一圆点，ζ=1从绘图和拟合中排除；图2以竖直虚线指示实测Q峰，上图零轴居中。原始科学验收仍沿用模式状态筛选。若可靠端点异号但中间含待核验模式，表中单独保留该区间，不将其当作已经核验的连续过零，也不误报为没有异号证据。其Q属于原始本征场，不是非等频组合场的Q。", "", "## 科学验收", ""]
    for mode in ("px","py"):
        s=summary[mode]
        lines += [f"### {mode}", "",
            f"- 电场/边界/独立损耗是否已定位到唯一共同零点：{'是' if s['local_position_agreement'] else '否，当前数据未验证位置完全一致'}。区间重叠仅是采样分辨率内的位置证据，不等于严格复振幅为零。",
            f"- 全区间连续、单调且唯一变号的完整验收：{'得到数值支持' if s['existence_uniqueness_supported'] else '尚未通过'}。身份未通过的ζ：{s['mixed_zeta']}。",
            f"- 独立空气面完整复振幅最小采样ζ={s['air_complex_minimum_zeta']:.4f}；最大Q={s['Q_maximum']:.6g}，发生在ζ={s['q_minimum_zeta']:.4f}。",
            f"- 在该Q峰值处，主空气通道复振幅模长为{s['air_complex_magnitude_at_Q_max_V_per_m']:.6g} V/m；两个空气面传播补偿后差值模长为{s['air_plane_difference_at_Q_max_V_per_m']:.6g} V/m。这些残余如实保留，不把实部变号自动解释为严格零辐射。", ""]
        for key,label in (("electric","电场面积分"),("boundary","边界加和显示值"),("air","空气振幅")):
            inc=s[key+"_increments"]
            lines.append(f"{label}：相邻可靠样本增量中，上升{inc['positive']}段、下降{inc['negative']}段、相等{inc['zero']}段。该统计包含局部加密，不用来要求幅值吻合；原始折线用于判断整体趋势和局部起伏。")
        lines.append("")
    lines += ["## 公式、相位、单位与来源", "",
        r"px边界原式为 $-i(1/\varepsilon_{air}-1/\varepsilon_{slab})\sum_j\oint n_yH_z dl/(\omega\varepsilon_0)$；py为 $+i(1/\varepsilon_{air}-1/\varepsilon_{slab})\sum_j\oint n_xH_z dl/(\omega\varepsilon_0)$。法向从介质指向孔内空气，时间约定exp(-iωt)，ω使用对应本征解的完整复频率。",
        "",
        "每个模式的全部E、H、边界贡献和空气面数据共用同一个Hz参考相位及Hz RMS=1 A/m幅度因子。px固定参考为ζ=0.9的已核验本征场；py沿用旧ζ=1.156的Hz参考。没有用电场/边界各自的相位分别转实，没有逐点翻号。",
        "",
        "**图面符号约定**：边界计算的实部统一乘−1后绘制；依用户要求，图例直接省略Re之前的外部负号。没有将这个显示负号移进括号。px原公式自身的负号仍属于物理表达式。完整复数原值及显示值分别保存；图面不宣称原式带符号相等。",
        "",
        "电场面积分与边界换算均为V·m，裸Hz线积分为A，空气区面积平均振幅为V/m。原始NPZ中的px_air_Vm/py_air_Vm键名中的Vm不作为单位定义，其实际单位为V/m，以本说明、字段元数据和CSV的显式单位为准。",
        "",
        "电场原生积分阶数8，阶数4留作核验；积分覆盖整个六角晶胞，包括孔内空气，逐面显式指定板层侧。空气面位于z=0.6H_air与0.75H_air（H_air=1550 nm），处于均匀空气区，分别补偿到z=100 nm。1/Q由原始复本征频率获得；有限损耗和剩余复振幅不做人工归零。",
        "",
        "模式反射对称误差使用原始Hz场计算；中心面阈值0.1、表面阈值0.15。超限同时可能反映实际混合或mesh=5的离散对称误差，网格保持5，超限原始本征场保留待核验标记。用户px对应内部py，用户py对应内部px；Ex/Ey坐标不交换。",
        "", "## 数据与复现", "",
        "- [复数汇总与实际绘图值](../../80_logs/09_px_py_validation/scan_results.csv)",
        "- [逐边贡献](../../80_logs/09_px_py_validation/edge_results.csv)；[逐孔贡献](../../80_logs/09_px_py_validation/hole_results.csv)",
        "- [零点与验收记录](../../80_logs/09_px_py_validation/scientific_summary.json)",
        "- [补充：统一相位后的完整Ex/Ey及两个空气面](../../80_logs/09_px_py_validation/supplemental_complex_components.csv)；[数据复算与源保护核验](../../80_logs/09_px_py_validation/verification.json)",
        "- [逐文件交付与哈希清单](../../80_logs/09_px_py_validation/delivery_manifest.json)",
        "- [参数、模式映射与两轮实际清单](../../99_config/09_px_py_validation/config.json)",
        "- [执行文档](../../../../docs/spec/plan-execute_20260928_s4_boundary_validation_four_figures.md)", "",
        "离线重绘：在仓库根运行 `.venv/Scripts/python.exe -B -m scripts.analysis.plot_boundary_dual_validation`，不会启动COMSOL。旧01–08分组未改动；本轮新模型归入00_model/09_px_py_validation，计算恢复目录保留，新主展示入口为本报告。原始逐点表保留导出时的模型路径，移动后的位置由交付清单逐文件解析；活动汇总已使用当前位置。", ""]
    for i,name in enumerate(NAMES,1):
        lines.append(f"- Figure {i}：[PNG](../../10_overview/{name}.png) · [PDF](../../11_pdf/{GROUP}/{name}.pdf)")
    path=ARCHIVE/"12_reports"/GROUP/"README.md"; path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")

def main():
    global GROUP
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-only",action="store_true")
    parser.add_argument("--deliver",action="store_true")
    parser.add_argument("--figure",type=int,choices=[1,2,3,4],help="Redraw only the selected Figure from delivered data")
    parser.add_argument("--group",choices=["09_px_py_validation","11_px_py_mesh025","12_px_py_peak_refine"],default=GROUP)
    args=parser.parse_args()
    if args.figure is not None:
        GROUP=args.group
        if args.summary_only or args.deliver:parser.error("--figure cannot be combined with --summary-only/--deliver")
        logs=ARCHIVE/"80_logs"/GROUP
        df=pd.read_csv(logs/"scan_results.csv"); summary=s4.json_read(logs/"scientific_summary.json")
        if args.figure==1:audit=figure1(df,summary,ARCHIVE)
        elif args.figure==2:audit=figure2(df,summary,ARCHIVE,mesh_comparison=GROUP=="09_px_py_validation")
        else:audit={3:figure3,4:figure4}[args.figure](df,pd.read_csv(logs/"edge_results.csv"),summary,ARCHIVE)
        path=ARCHIVE/"99_config"/GROUP/"plot_contract.json"
        contract=s4.json_read(path); contract.update(audit); s4.write_json(path,contract)
        print(f"Updated Figure {args.figure} only; no COMSOL",flush=True)
        return
    config=prepare(); df=analyze(config)
    summary=summarize(df)
    s4.write_json(ROOT/"80_logs/scientific_summary.json",summary)
    print(summary,flush=True)
    if not args.summary_only:
        edges=pd.read_csv(ROOT/"80_logs/edge_results.csv")
        if args.deliver:
            if len(df)!=2*len(config["cases"]):raise ValueError("Incomplete scan")
            validate_saved(config,df,edges)
            df=deliver(config,df,edges,summary)
        audit=figures(df,edges,summary,ARCHIVE)
        s4.write_json(ARCHIVE/"99_config"/GROUP/"plot_contract.json",audit)
        if args.deliver:write_report(config,df,summary,audit)


if __name__=="__main__":main()