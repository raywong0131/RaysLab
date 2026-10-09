import numpy as np
from scripts.analysis.prepare_s5_gamma_origin import complex_split,geometry_mapping,seam_metrics
from scripts.run_main.run_s5_gamma_origin import reflected_sample,PARITIES
from comsol_workflow.hex_lattice_utils import lattice_points,unit_cell_corners


def test_complete_average_and_seam_detects_jump():
    rng=np.random.default_rng(3);v=rng.normal(size=(7,6,61,6))+1j*rng.normal(size=(7,6,61,6))
    gamma,rest,checks=complex_split(v)
    assert max(checks.values())<1e-14
    np.testing.assert_allclose(gamma,np.mean(v,axis=0));np.testing.assert_allclose(rest+gamma,v)
    v[:]=0;v[...,0]=1;v[...,1]=2
    _,records=seam_metrics(v,unit_cell_corners(.82));assert max(r['relative_tangential_jump'] for r in records)==0
    v[:,0,:,1]+=1
    _,records=seam_metrics(v,unit_cell_corners(.82));assert records[0]['tangential_jump_rms_V_m']>0


def test_boundary_and_reflections():
    points=lattice_points(19,.82);centers=np.array([[p.x,p.y] for p in points]);corners,boundary,internal=geometry_mapping(centers,.82)
    assert len(centers)==1141 and len(boundary)==234
    assert len(boundary)+2*len(internal)==6*1141
    def sampler(mode,xyz,expr):
        assert np.all(xyz>=0)
        return np.ones((len(xyz),6),complex)
    xyz=np.array([[1,2,3],[-1,2,3],[1,-2,3],[1,2,-3]])
    got=reflected_sample(sampler,xyz)
    np.testing.assert_array_equal(got,np.vstack([np.ones(6),PARITIES]))


def test_triangle_area_2d_and_closed_contour():
    from scripts.run_main.run_s5_gamma_origin import triangle_areas
    from scripts.run_main.run_s5_gamma_replay import contour_from_edges
    triangles=np.array([[[0,0],[2,0],[0,3]],[[0,0],[0,3],[2,0]],[[0,0],[1,1],[2,2]]],float)
    np.testing.assert_allclose(triangle_areas(triangles),[3,3,0])
    corners=unit_cell_corners(.82)
    segments=np.stack([corners,np.roll(corners,-1,axis=0)],axis=1)
    contour=contour_from_edges(segments)
    np.testing.assert_allclose(contour,corners,atol=1e-10)


def test_shared_trace_projection_is_linear_and_preserves_full_seam():
    from scripts.analysis.prepare_s5_gamma_origin import shared_node_projection
    mapping=np.array([0,1,1,2]);counts=np.array([1,2,1])
    full=np.array([[1+2j],[3-1j],[3-1j],[4j]])
    gamma=np.array([[2j],[1j],[2j],[0j]])
    q=lambda v:shared_node_projection(v,mapping,counts)
    np.testing.assert_allclose(q(full),q(gamma)+q(full-gamma))
    np.testing.assert_allclose(q(full)[mapping],full)
    np.testing.assert_allclose(q(gamma)[1],1.5j)


def test_planar_chart_selection_excludes_touching_neighbor():
    from scripts.run_main.run_s5_gamma_replay import planar_chart_boundaries
    faces={1:np.array([[1,0,0],[1,1,0],[1,1,.1],[1,0,.1]]),
           2:np.array([[1,1,0],[0,1,0],[0,1,.1],[1,1,.1]]),
           3:np.array([[1,0,.1],[1,1,.1],[1,1,.3],[1,0,.3]])}
    chart=dict(kind='side',normal=[1,0],offset_um=1,zlo=0,zhi=.1)
    assert planar_chart_boundaries(faces,chart,.001)==[1]


def test_positive_zoom_fft_matches_direct_complex_integral():
    from scripts.analysis.analyze_s5_gamma_origin import positive_zoom_fft2
    axis=np.arange(-3,4)*.2;xx,yy=np.meshgrid(axis,axis)
    field=(1+.2*xx+1j*yy)*np.exp(-1j*(.7*xx-.4*yy))
    q=np.linspace(-2e6,2e6,9);phase=np.exp(1j*q[:,None]*axis[None,:]*1e-6)
    expected=phase@field@phase.T*(.2e-6)**2
    np.testing.assert_allclose(positive_zoom_fft2(field,axis,q),expected,rtol=1e-12,atol=1e-25)


def test_cap_fft_cancelling_center_avoids_unresolved_ratio():
    from scripts.analysis.analyze_s5_gamma_origin import cap_fft_spectra
    xy=np.array([[-.04,0.],[.04,0.]])
    full=np.array([-1.,1.],complex)
    values=np.array([full,.4*full,.6*full])
    spectrum,q,checks,center=cap_fft_spectra(xy,values,np.ones(2)*1e-12,2e14)
    assert not checks['center_FULL_resolved']
    assert checks['center_REST_over_FULL'] is None
    assert checks['center_quadrature_error_over_spectrum_peak']<1e-12
    np.testing.assert_allclose(spectrum[0],spectrum[1]+spectrum[2],atol=1e-25)


def test_cap_lattice_samples_matches_direct_integral_with_shared_nodes():
    from scripts.analysis.analyze_s5_gamma_origin import cap_lattice_samples
    from comsol_workflow.finite_lattice_fourier import direct_lattice_basis,hex_cyclic_indices,hex_xi_values,fold_xis_to_first_bz
    indices=np.array([[0,0],[1,0],[0,1],[-1,1],[-1,0],[0,-1],[1,-1]])[[3,0,6,1,5,2,4]]
    centers=indices@direct_lattice_basis(.82).T
    local=np.array([[0.,0.],[.41,0.],[-.41,0.]])
    weights=np.array([.2,.1,.1])
    xy,mapping=np.unique(np.round((centers[:,None]+local).reshape(-1,2),12),axis=0,return_inverse=True)
    mapping=mapping.reshape(7,3)
    rng=np.random.default_rng(34)
    full=rng.normal(size=len(xy))+1j*rng.normal(size=len(xy))
    gamma=rng.normal(size=len(xy))+1j*rng.normal(size=len(xy))
    values=np.array([full,gamma,full-gamma])
    k=fold_xis_to_first_bz(hex_xi_values(1),.82)
    actual=cap_lattice_samples(values,mapping,local,weights,hex_cyclic_indices(indices,1),k)
    global_weights=np.bincount(mapping.ravel(),weights=np.tile(weights,7))*1e-12
    expected=np.array([np.sum(values*(global_weights*np.exp(1j*xy@q)),axis=1) for q in k]).T
    np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-23)
    np.testing.assert_allclose(actual[0],actual[1]+actual[2],atol=1e-26)


def test_cap_summary_combined_modulus_and_shared_scales(monkeypatch):
    from scripts.analysis import analyze_s5_gamma_origin as m
    axis=np.linspace(-.08,.08,3);values=np.arange(9).reshape(3,3)/8
    maps=np.ma.array([values,.5*values,.25*values]);points=np.array([[0.,0.],[.04,0.]])
    discrete=np.array([[2.,1.],[1.,.5],[.5,.25]])
    bz=np.array([[-.5,0],[0,.5],[.5,0],[0,-.5],[-.5,0]])
    def inspect(fig,stem):
        assert len(fig.axes)==12 and not any(a.child_axes for a in fig.axes)
        assert stem=='01_boundary_inputs_ExEy'
        for ax in fig.axes[:3]+fig.axes[6:9]:assert ax.collections[0].get_clim()==(0,1)
        for ax in fig.axes[3:6]:
            assert ax.collections[0].get_clim()==(1e-4,1)
            assert isinstance(ax.collections[0].norm,m.LogNorm)
        np.testing.assert_allclose(fig.axes[10].get_yticks(),np.logspace(-4,0,5))
        assert 'E_{xy}' in fig.axes[0].get_title()
        assert fig.axes[0].collections[0].get_cmap().name=='magma'
        np.testing.assert_allclose(fig.axes[4].collections[0].get_array(),[.5,.25])
        np.testing.assert_allclose(fig.axes[8].collections[0].get_array(),.25*values)
        m.plt.close(fig)
        return 'checked'
    monkeypatch.setattr(m,'save',inspect)
    assert m.plot_cap_summary(axis,maps,axis,points,points,discrete,maps,bz,(1.,2.,1.),'ExEy')=='checked'


def test_cell_component_weights_parseval_and_coherent_reconstruction():
    from scripts.analysis.analyze_s5_gamma_origin import cell_k_decomposition
    n=7;r=np.arange(n)
    profile=np.array([1+2j,3-1j]);other=np.array([2-1j,-1j]);area=np.array([.1,.3])
    field=profile[None,:]+np.exp(-2j*np.pi*2*r[:,None]/n)*other[None,:]
    f,w,reconstructed,checks=cell_k_decomposition(field,area)
    np.testing.assert_allclose(f[0],np.sqrt(n)*profile,atol=1e-14)
    np.testing.assert_allclose(f[2],np.sqrt(n)*other,atol=1e-14)
    np.testing.assert_allclose(w[[0,2]],n*np.array([np.sum(abs(profile)**2*area),np.sum(abs(other)**2*area)])*1e-12)
    np.testing.assert_allclose(reconstructed,field,atol=1e-14)
    assert max(checks.values())<1e-12
