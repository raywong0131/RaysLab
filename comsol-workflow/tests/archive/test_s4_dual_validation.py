import sys
import numpy as np
from scripts.run_sweep.run_s4_dual_validation import assign_modes, boundary_values, parity
from comsol_workflow.s4_boundary import EPS0


def test_physical_modes_are_assigned_without_rotating_eigenfields():
    x=np.array([1.,2.,3.])
    # Candidate order deliberately opposite to physical px/py.
    fields=[(x,-x,x),(2j*x,2j*x,-2j*x)]
    before=[tuple(a.copy() for a in f) for f in fields]
    assignment,cost=assign_modes(fields)
    assert assignment==(1,0)
    assert cost[1,0]==cost[0,1]==0
    assert cost[0,0]>1 and cost[1,1]>1
    for f,b in zip(fields,before):
        for a,v in zip(f,b):np.testing.assert_array_equal(a,v)
    assert "mph" not in sys.modules


def test_boundary_coordinates_frequency_and_common_phase():
    h=np.array([[1+2j,3-1j],[2j,1-3j]])
    w=np.array([[.2,.3],[.1,.4]])
    normals=np.array([[1.,0.],[0.,-1.]])
    omega=2*np.pi*(200e12-1e8j)
    out=boundary_values(h,w,normals,omega)
    integral=(h*w).sum(axis=1)
    np.testing.assert_allclose(out[0],[0,1j*(1-1/3.3**2)*integral[0]/(omega*EPS0)])
    np.testing.assert_allclose(out[1],[1j*(1-1/3.3**2)*integral[1]/(omega*EPS0),0])
    phase=2.4*np.exp(.7j)
    np.testing.assert_allclose(boundary_values(h*phase,w,normals,omega),out*phase)


def test_mixed_eigenmode_is_reported_not_symmetrized():
    x=np.array([1.,2.,3.]); y=np.array([3.,2.,1.])
    h=x+y
    assert parity(h,-x+y,x-y,"py")>.1
    assert parity(h,-x+y,x-y,"px")>.1


def test_crossings_do_not_bridge_ambiguous_modes():
    from scripts.analysis.plot_s4_dual_validation import crossings, formula
    x=np.array([.8,.9,1.,1.1,1.2])
    y=np.array([-2.,-1.,0.,1.,2.])
    assert crossings(x,y,np.array([True,False,False,True,True]))==[]
    roots=crossings(x,y,np.ones(5,bool))
    assert len(roots)==1 and roots[0]["left"]==roots[0]["right"]==1.
    assert "n_y" in formula("px") and "n_x" in formula("py")
    assert not formula("px").startswith("$-\\mathrm{Re}")
    assert "mph" not in sys.modules

def test_four_figure_contract_and_common_Hz_scale(tmp_path,monkeypatch):
    import pandas as pd
    import matplotlib.pyplot as plt
    from scripts.analysis import plot_s4_dual_validation as plotting
    angle=np.arange(6)*np.pi/3
    outer=np.c_[np.cos(angle),np.sin(angle)]*400e-9
    holes=np.array([np.c_[np.cos(np.arange(3)*2*np.pi/3),np.sin(np.arange(3)*2*np.pi/3)]*25e-9+np.array([np.cos(a),np.sin(a)])*200e-9 for a in angle])
    xy=np.vstack([outer,[[0,0]],holes.mean(axis=1)])
    path=tmp_path/'fields.npz'
    np.savez(path,area_xy_m=xy,outer_m=outer,holes_m=holes,px_surface_Hz=np.arange(len(xy))+0j,py_surface_Hz=-2*np.arange(len(xy))+0j)
    rows=[]; edge_rows=[]
    for mode in ('px','py'):
        for z,y in ((.8,-1),(.9,-.5),(.98,0),(1.,20),(1.1,.5),(1.2,1)):
            rows.append(dict(zeta=z,user_mode=mode,mode_status='mixed_eigenmode' if z==.8 else 'single_eigenmode',electric_display_Vm=y,boundary_display_Vm=2*y,air_display_V_per_m=y,inverse_Q=.01+abs(y),Q=1/(.01+abs(y)),field_file=str(path),factor_re=1,factor_im=0))
            for h in range(1,7):
                for e in range(1,4):edge_rows.append(dict(zeta=z,user_mode=mode,hole=h,edge=e,display_Vm=y/9+(-1)**(h*3+e)*.1))
    summary={mode:dict(q_minimum_zeta=.98,q_minimum_bracket=[.9,.98,1.1],representative_zeta=[.8,.98,1.2],electric_crossings=[],boundary_crossings=[]) for mode in ('px','py')}
    captured=[]
    def check(fig,name,target):
        captured.append(name)
        for ax in fig.axes:
            assert not any(g.get_visible() for g in ax.get_xgridlines()+ax.get_ygridlines())
            legend=ax.get_legend()
            if legend is not None:assert legend._loc==1
            assert 'real part' not in ax.get_title() and not ax.get_title().startswith(('A:','B:','C:','D:'))
        if name=='figure1_boundary_validation':
            np.testing.assert_allclose(fig.get_size_inches(),[12,9])
            for ax in fig.axes:
                points=[line for line in ax.lines if line.get_marker()=="o"]
                assert len(points)==2 and all(len(line.get_xdata())==5 for line in points)
                assert not ax.collections
                assert not any(line.get_marker()=="x" for line in ax.lines)
                fitted=[line for line in ax.lines if len(line.get_xdata())==802]
                assert len(fitted)==2
                assert all(not np.any(np.isclose(line.get_xdata(),1.)) for line in points)
                patch=ax.patches[0]
                center=patch.get_x()+patch.get_width()/2
                np.testing.assert_allclose(patch.get_width(),.002)
                for line in fitted:
                    idx=np.argmin(abs(line.get_xdata()-center))
                    np.testing.assert_allclose(line.get_xdata()[idx],center,atol=1e-14)
                    np.testing.assert_allclose(line.get_ydata()[idx],0,atol=1e-12)
        if name=='figure2_radiation_zero':
            assert len(fig.axes)==2
            assert fig.axes[1].get_title()==r"$Q$"
            assert fig.axes[1].get_ylim()[1]==1e8
            np.testing.assert_allclose(sum(fig.axes[0].get_ylim()),0,atol=1e-12)
            assert not fig.axes[1].texts
            assert all(not ax.patches for ax in fig.axes)
            np.testing.assert_allclose(fig.get_size_inches(),[9,9])
            for ax in fig.axes:
                points=[line for line in ax.lines if line.get_marker()=="o"]
                fitted=[line for line in ax.lines if len(line.get_xdata())==801]
                assert len(points)==2 and len(fitted)==4
                for left,right in (fitted[:2],fitted[2:]):
                    np.testing.assert_allclose(left.get_ydata()[-1],right.get_ydata()[0],rtol=1e-12)
                    assert left.get_xdata()[-1]==right.get_xdata()[0]
                assert len(points[0].get_xdata())==5 and not ax.collections
                assert all(not np.any(np.isclose(line.get_xdata(),1.)) for line in points)
                guides=[line for line in ax.lines if line.get_linestyle()=="--"]
                assert len(guides)==2 and all(np.all(np.asarray(line.get_xdata())==.98) for line in guides)
                if ax.get_yscale()=="log":
                    assert all(np.all(line.get_ydata()>0) for line in fitted)
                    np.testing.assert_allclose(points[0].get_ydata(),1/(.01+np.array([1,.5,0,.5,1])))
        if name=='figure4_boundary_cancellation':
            for ax in fig.axes[:2]:
                assert all(not np.any(np.isclose(line.get_xdata(),1.)) for line in ax.lines if len(line.get_xdata())>2)
        if name=='figure3_Hz_and_boundaries':
            fig.canvas.draw()
            assert len(fig.axes)==6
            for start in (0,3):
                field,spatial,bars=fig.axes[start:start+3]
                assert field.collections[0].get_cmap().name=='RdBu_r'
                assert field.collections[0].get_clim()==(-1,1)
                assert len(bars.patches)==18
                np.testing.assert_allclose(spatial.collections[0].get_array(),[b.get_height() for b in bars.patches])
                np.testing.assert_allclose(spatial.collections[0].get_clim(),bars.get_ylim())
                assert len(spatial.collections[0].get_paths())==18
                assert all(len(path.vertices)==2 for path in spatial.collections[0].get_paths())
                assert len(spatial.collections)==1  # One colored stroke; no gray backing/shadow.
                assert spatial.collections[0].get_capstyle()=="round"
                np.testing.assert_allclose(spatial.collections[0].get_linewidths(),2.0)
                assert not spatial.patches  # No normal arrows covering the contributions.
                for ax in (field,spatial):
                    cax=ax.child_axes[0]
                    np.testing.assert_allclose([cax.get_position().y0,cax.get_position().y1],[ax.get_position().y0,ax.get_position().y1])
                np.testing.assert_allclose(field.child_axes[0].get_yticks(),np.arange(-1,1.01,.25))
        plt.close(fig)
    monkeypatch.setattr(plotting,'save',check)
    audit=plotting.figures(pd.DataFrame(rows),pd.DataFrame(edge_rows),summary,tmp_path)
    assert captured==plotting.NAMES
    assert audit['Hz_common_max_A_per_m']==24
    assert audit['boundary_display_sign']==-1

def test_opposite_endpoints_are_preserved_without_claiming_continuity():
    import pandas as pd
    from scripts.analysis.plot_s4_dual_validation import summarize
    rows=[]
    for mode in ('px','py'):
        for z,y,status in ((.8,-1.,'single_eigenmode'),(1.,0.,'mixed_eigenmode'),(1.2,1.,'single_eigenmode')):
            rows.append(dict(zeta=z,user_mode=mode,mode_status=status,Q=100.,inverse_Q=.01,
                electric_display_Vm=y,boundary_display_Vm=y,air_display_V_per_m=y,
                air_amplitude_V_per_m_re=y,air_amplitude_V_per_m_im=0.,air2_amplitude_V_per_m_re=y,air2_amplitude_V_per_m_im=0.))
    result=summarize(pd.DataFrame(rows))
    for mode in ('px','py'):
        assert result[mode]['electric_crossings']==[]
        assert result[mode]['electric_opposite_endpoints_across_unresolved_points']==[dict(left=.8,right=1.2,ambiguous_interior_zeta=[1.])]
        assert not result[mode]['existence_uniqueness_supported']

def test_display_fit_is_global_least_squares_not_interpolation():
    from scripts.analysis.plot_s4_dual_validation import trend_fit
    x=np.array([.8,.9,1.,1.1,1.2]);y=np.array([-2.,-.8,.5,.9,2.])
    model=trend_fit(x,y);residual=y-model(x)
    assert model.degree()==2
    assert np.max(abs(residual))>.01
    np.testing.assert_allclose([sum(residual),np.dot(x,residual),np.dot(x*x,residual)],0,atol=1e-12)


def test_positive_q_fit_recovers_positive_quadratic():
    from scripts.analysis.plot_s4_dual_validation import positive_q_fit
    x=np.linspace(.8,1.2,81)
    expected=np.array([2e-8,.03,1.156])
    y=expected[0]+expected[1]*(x-expected[2])**2
    np.testing.assert_allclose(positive_q_fit(x,y),expected,rtol=1e-5)


def test_figure2_refined_mesh_preserves_phase_scale_and_sample_range(tmp_path,monkeypatch):
    import pandas as pd
    import pytest
    import matplotlib.pyplot as plt
    from scripts.analysis import plot_s4_dual_validation as plotting
    x=np.array([.8,.85,.9,.98,1.,1.1,1.2])
    df=pd.DataFrame([dict(user_mode=mode,zeta=z,factor_re=0.,factor_im=2.,
        air_display_V_per_m=z-.98,Q=1/(.01+abs(z-.98)),inverse_Q=.01+abs(z-.98))
        for mode in ("px","py") for z in x])
    mesh=pd.DataFrame(dict(mesh=[5,5,5,2],size_factor=[.25,.25,.25,1.],
        zeta=[.8,.85,.9,.8],Q=[4.7e8,1e7,1e6,9e9],air_Ex_re=[1.,2.,3.,100.],air_Ex_im=[2.,3.,4.,100.]))
    corrected=plotting.refined_px(df,mesh)
    np.testing.assert_allclose(corrected.air_display_V_per_m,[-2.,-3.,-4.])
    np.testing.assert_allclose(corrected.Q,[4.7e8,1e7,1e6])
    with pytest.raises(ValueError,match="matching nonzero scan phase"):
        plotting.refined_px(df[df.zeta!=.85],mesh)
    path=tmp_path/"80_logs/10_px_mesh_test/mesh_results.csv"
    path.parent.mkdir(parents=True);mesh.to_csv(path,index=False)
    def check(fig,name,target):
        assert len(fig.axes)==2 and len(fig.axes[0].child_axes)==1
        assert fig.axes[1].get_ylim()[1]==1e8
        for ax in fig.axes:
            points=[line for line in ax.lines if line.get_marker()=="o"]
            assert len(points)==3
            np.testing.assert_allclose(points[-1].get_xdata(),[.8,.85,.9])
            refined=[line for line in ax.lines if line.get_linestyle()=="--" and len(line.get_xdata())>2]
            assert len(refined)==1
            np.testing.assert_allclose([refined[0].get_xdata()[0],refined[0].get_xdata()[-1]],[.8,.9])
        np.testing.assert_allclose(points[-1].get_ydata(),[4.7e8,1e7,1e6])
        assert not fig.axes[1].texts
        plt.close(fig)
    monkeypatch.setattr(plotting,"save",check)
    summary={mode:dict(q_minimum_zeta=.98,q_minimum_bracket=[.9,1.1]) for mode in ("px","py")}
    audit=plotting.figure2(df,summary,tmp_path)["figure2"]["fits"]["px_refined"]
    assert audit["sample_count"]==3 and audit["anchor_Q"]==4.7e8
    assert not audit["convergence_established"]


def test_save_keeps_current_overview_flat(tmp_path,monkeypatch):
    from scripts.analysis import plot_s4_dual_validation as plotting
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plotting,"GROUP","11_px_py_mesh025")
    for name in plotting.NAMES:
        fig=plt.figure(figsize=(1,1))
        plotting.save(fig,name,tmp_path)
    overview=tmp_path/"10_overview"
    assert sorted(p.name for p in overview.iterdir())==sorted(n+".png" for n in plotting.NAMES)
    assert all(p.is_file() for p in overview.iterdir())
    assert len(list((tmp_path/"11_pdf/11_px_py_mesh025").glob("*.pdf")))==4






def test_constrained_fit_preserves_samples_and_minimizes_residual():
    from scripts.analysis.plot_s4_dual_validation import trend_fit
    x=np.array([.8,.9,1.1,1.2]);y=np.array([2.,1.,-.5,-2.]);original=y.copy();zero=.97
    model=trend_fit(x,y,zero_at=zero)
    np.testing.assert_allclose(model(zero),0,atol=1e-12)
    residual=y-model(x);offset=x-zero
    np.testing.assert_allclose([np.dot(offset,residual),np.dot(offset**2,residual)],0,atol=1e-12)
    np.testing.assert_array_equal(y,original)
