import json

import numpy as np
import pytest

from scripts.analysis.validate_s5_bloch_radiation import (
    consume_radiation_bundle, radiation_budget, weighted_projection,
)
from scripts.run_main.run_s5_periodic_reference import mode_assignment, periodic_profiles


def test_nonorthogonal_projection_and_radiation_gauge():
    rng = np.random.default_rng(19)
    basis = rng.normal(size=(31,4))+1j*rng.normal(size=(31,4))
    basis[:,1] += .6*basis[:,0]
    expected = np.array([1+.2j,-.8j,.3,-.1+.7j])
    field = basis@expected
    weights = rng.uniform(.2,2,31)
    a,fit,error,condition = weighted_projection(basis,field,weights)
    np.testing.assert_allclose(a,expected,atol=1e-13)
    np.testing.assert_allclose(fit,field,atol=1e-13)
    assert error < 1e-13 and condition > 1
    response = rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
    gauge = np.exp(1j*rng.uniform(-np.pi,np.pi,4))
    rotated,_,_,_ = weighted_projection(basis*gauge,field,weights)
    np.testing.assert_allclose(radiation_budget(a,response)[0],
                               radiation_budget(rotated,response*gauge[:,None])[0],atol=1e-13)


def test_constructive_destructive_and_quadrature_phases():
    for second,expected in [(1,4),(-1,0),(1j,2)]:
        terms,total,inc,coh,cross=radiation_budget(np.array([1,second]),np.ones((2,2)))
        np.testing.assert_allclose(inc,2)
        np.testing.assert_allclose(coh,expected,atol=1e-14)
        np.testing.assert_allclose(cross,expected-2,atol=1e-14)


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError,match='Rank-deficient'):
        weighted_projection(np.ones((3,2)),np.ones(3),np.ones(3))
    with pytest.raises(ValueError,match='positive'):
        weighted_projection(np.ones((3,1)),np.ones(3),np.array([1,0,1]))
    with pytest.raises(ValueError,match='matching'):
        radiation_budget(np.ones(3),np.ones((2,2)))


def test_bloch_sign_and_mode_matching():
    rng=np.random.default_rng(21)
    xy=rng.normal(size=(50,2)); k=np.array([.3,-.4])
    u=rng.normal(size=(4,50))+1j*rng.normal(size=(4,50))
    np.testing.assert_allclose(periodic_profiles(u*np.exp(-1j*(xy@k)),xy,k),u,atol=1e-14)
    permutation=np.array([2,0,3,1])
    assignment,_=mode_assignment(u,u[permutation]*np.exp(1j*rng.normal(size=4))[:,None])
    np.testing.assert_array_equal(permutation[assignment],np.arange(4))


def test_radiation_closure_requires_independent_full_vector_data(tmp_path):
    path=tmp_path/'response.npz'
    np.savez(path,metadata_json=json.dumps({'independent_response':False}))
    with pytest.raises(ValueError,match='provenance'):
        consume_radiation_bundle(path,np.zeros((19,2)))
    a=np.ones((19,4),complex); t=np.zeros((19,4,2),complex)
    t[1,1]=[1,1]; t[2,1]=[-1,1]
    metadata={'independent_response':True,'full_vector_projection':True,
              'same_observation_and_units':True,'response_provenance':'synthetic test'}
    np.savez(path,metadata_json=json.dumps(metadata),k_um_inv=np.zeros((19,2)),
             A=a,T=t,C_direct=np.array([0,2]))
    result=consume_radiation_bundle(path,np.zeros((19,2)))
    assert result['relative_amplitude_closure_error']==0
    np.testing.assert_allclose(result['cross'],[-2,2])


def test_volume_radiation_matches_dipole_and_air_spectrum():
    from scripts.analysis.validate_s5_bloch_radiation import volume_far_amplitude,planar_K0_to_far
    omega=2*np.pi*198e12; k=omega/299792458.
    volume=2e-24; z_source=30e-9; z_air=1.5e-6
    e=np.array([[2+1j,-3j,900.]])
    contrast=np.array([9.89]); xyz=np.array([[0,0,z_source]])
    f=volume_far_amplitude(e,contrast,xyz,np.array([volume]),omega)
    polarization=e[0,:2]*contrast[0]*volume
    expected=k*k/(4*np.pi)*polarization*np.exp(1j*k*z_source)
    np.testing.assert_allclose(f,expected,rtol=1e-14,atol=0)
    c_air=-1j*k/2*polarization*np.exp(-1j*k*(z_air-z_source))
    np.testing.assert_allclose(planar_K0_to_far(c_air,z_air,omega),expected,rtol=1e-14,atol=0)
    # A normal electric dipole has a radiation null on its axis.
    np.testing.assert_array_equal(volume_far_amplitude(np.array([[0,0,900]]),contrast,xyz,np.array([volume]),omega),[0,0])


def test_volume_bundle_keeps_radiation_of_projection_residual(tmp_path):
    from scripts.analysis.validate_s5_bloch_radiation import consume_volume_bundle,volume_far_amplitude
    rng=np.random.default_rng(302)
    nodes=32; xyz=rng.normal(size=(nodes,3))*1e-7
    weights=np.full(nodes,1e-24);chi=np.full(nodes,9.89)
    basis=rng.normal(size=(19,4,nodes,3))+1j*rng.normal(size=(19,4,nodes,3))
    a=rng.normal(size=(19,4))+1j*rng.normal(size=(19,4))
    field=np.einsum('kb,kbnp->np',a,basis)+rng.normal(size=(nodes,3))
    omega=2*np.pi*198e12;k=omega/299792458.;z_air=1.6e-6
    total=volume_far_amplitude(field,chi,xyz,weights,omega)
    direct_air=total*2*np.pi/(1j*k)*np.exp(-1j*k*z_air)
    path=tmp_path/'volume.npz'
    meta={'time_convention':'native_COMSOL_exp(+iwt)','complete_material_volume':True,
          'finite_field_source':'synthetic algebra test','basis_source':'synthetic algebra test'}
    np.savez(path,metadata_json=json.dumps(meta),k_um_inv=np.zeros((19,2)),E_basis=basis,E_finite=field,
             xyz_m=xyz,weights_m3=weights,epsilon_contrast=chi,omega_rad_s=omega,
             direct_air_amplitude_V_m=direct_air,z_air_m=z_air)
    result=consume_volume_bundle(path,np.zeros((19,2)))
    assert result['volume_vs_direct_air_relative_error']<1e-13
    assert result['full_vector_reconstruction_error']>0
    assert result['basis_vs_direct_air_relative_error']>1e-5


def test_factorized_volume_projection_and_radiation_match_dense_dictionary():
    from scripts.analysis.plot_s5_volume_radiation import fit_cells, responses
    from scripts.analysis.validate_s5_bloch_radiation import volume_far_amplitude
    rng = np.random.default_rng(21)
    k = rng.normal(size=(19,2))
    rho = rng.uniform(-.3,.3,(32,3)); rho[:,2] = abs(rho[:,2])*.2
    w = rng.uniform(.1,1,len(rho))*1e-22
    masks = rng.random((2,len(rho)))>.3
    electric = rng.normal(size=(19,4,len(rho),3))+1j*rng.normal(size=(19,4,len(rho),3))
    finite = {}; blocks=[]; fields=[]; metrics=[]; centers=[]
    for region in range(2):
        r = rng.normal(size=(7,2)); centers.append(r)
        field = rng.normal(size=(7,len(rho),3))+1j*rng.normal(size=(7,len(rho),3))
        finite[f'E_{region}']=field; finite[f'centers_{region}_um']=r
        b = np.einsum('rk,kbnp->kbrnp',np.exp(-1j*r@k.T),electric).reshape(76,-1).T
        blocks.append(b); fields.append(field.ravel())
        metrics.append(np.tile(np.repeat(2*w*(1+9.89*masks[region]),3),len(r)))
    a,error,cond=fit_cells(electric,finite,k,rho,w,masks)
    cavity,_,_=fit_cells(electric,finite,k,rho,w,masks,regions=(0,))
    cavity_expected,_,_,_=weighted_projection(blocks[0],fields[0],metrics[0])
    np.testing.assert_allclose(cavity.ravel(),cavity_expected,rtol=1e-8,atol=1e-9)
    expected,_,expected_error,_=weighted_projection(blocks[0],fields[0],metrics[0])
    np.testing.assert_allclose(a.ravel(),expected,rtol=1e-8,atol=1e-9)
    np.testing.assert_allclose(error,expected_error,rtol=1e-10)
    omega=2*np.pi*(198+1j*.01)*1e12
    response=responses(np.zeros((1,2)),k,electric,rho,w,masks,centers,omega)[0]
    combined=np.einsum('kb,kbnp->knp',a,electric)[:,None]
    combined_response=responses(np.zeros((1,2)),k,combined,rho,w,masks,centers,omega)[0,:,0]
    np.testing.assert_allclose(combined_response,(response*a[:,:,None]).sum(axis=1),rtol=1e-10,atol=1e-20)
    weighted_response=responses(np.zeros((1,2)),k,electric,rho,None,None,centers,omega,w[None,:]*masks)[0]
    np.testing.assert_allclose(weighted_response,response,rtol=1e-12,atol=1e-24)
    dense=np.zeros((19,4,2),complex)
    for region,r in enumerate(centers[:1]):
        pos=rho[None,:,:]+np.pad(r,((0,0),(0,1)))[:,None,:]
        b=np.einsum('rk,kbnp->kbrnp',np.exp(-1j*r@k.T),electric)
        for zsign in [1,-1]:
            dense+=volume_far_amplitude((b*np.array([1,1,zsign])).reshape(19,4,-1,3),
                np.tile(9.89*masks[region],len(r)),(pos*np.array([1,1,zsign])).reshape(-1,3)*1e-6,
                np.tile(w,len(r)),omega)
    np.testing.assert_allclose(response,dense,rtol=1e-11,atol=1e-22)
    # Actual cladding data must not enter either cavity coefficients or responses.
    poisoned={**finite,'E_1':np.full_like(finite['E_1'],np.nan),'centers_1_um':np.full_like(centers[1],np.nan)}
    cavity_again,_,_=fit_cells(electric,poisoned,k,rho,w,masks)
    np.testing.assert_allclose(cavity_again,a,rtol=1e-12,atol=1e-12)
    with pytest.raises(ValueError,match='only cover cavity'):
        fit_cells(electric,finite,k,rho,w,masks,regions=(0,1))
    cavity_response=responses(np.zeros((1,2)),k,electric,rho,w,masks,[centers[0],np.full_like(centers[1],np.nan)],omega)[0]
    np.testing.assert_allclose(cavity_response,response,rtol=1e-12,atol=1e-24)
    poisoned_weights=np.array([w*masks[0],np.full_like(w,np.nan)])
    np.testing.assert_allclose(responses(np.zeros((1,2)),k,electric,rho,None,None,centers,omega,poisoned_weights)[0],response,rtol=1e-12,atol=1e-24)
    # Check the full transverse projection at an off-axis observation against
    # an independent direct sum of both slab halves, not the factorized kernel.
    q=np.array([[.31,-.17]])
    k0=omega/299792458.;kz=np.sqrt(k0*k0-np.sum(q[0]**2)*1e12)
    direction=np.r_[q[0]*1e6,kz]/k0
    direct=np.zeros((19,4,2),complex)
    r=centers[0];pos=rho[None,:,:]+np.pad(r,((0,0),(0,1)))[:,None,:]
    b=np.einsum('rk,kbnp->kbrnp',np.exp(-1j*r@k.T),electric)
    for sign in [1,-1]:
        xyz=(pos*np.array([1,1,sign])).reshape(-1,3)*1e-6
        field=(b*np.array([1,1,sign])).reshape(19,4,-1,3)
        kernel=np.tile(w*9.89*masks[0],len(r))*np.exp(1j*k0*(xyz@direction))
        polarization=np.einsum('kbnc,n->kbc',field,kernel)
        transverse=polarization-np.einsum('kbc,c->kb',polarization,direction)[:,:,None]*direction
        direct+=k0*k0/(4*np.pi)*transverse[:,:,:2]
    np.testing.assert_allclose(responses(q,k,electric,rho,w,masks,centers,omega)[0],direct,rtol=1e-10,atol=1e-22)



def test_mesh_quadrature_and_clipped_material_volume():
    from scripts.run_main.run_s5_mesh_radiation import tetra_rule,clipped_tetrahedra
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]])
    points,weights=tetra_rule(vertices)
    np.testing.assert_allclose(weights.sum(),1/6)
    np.testing.assert_allclose(weights@points,np.ones(3)/24)
    np.testing.assert_allclose(weights@(points[:,0]**2),1/60)
    pieces=clipped_tetrahedra(vertices,np.array([[0,0],[.5,0],[0,.5]]))
    np.testing.assert_allclose(sum(tetra_rule(p)[1].sum() for p in pieces),1/12)
