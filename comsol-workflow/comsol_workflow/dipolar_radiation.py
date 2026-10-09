"""Shared complex projections and radiation operators; retain native phase conventions."""
import numpy as np
from scipy.constants import epsilon_0, mu_0
from scipy.fft import fft2, ifft2, fftshift, ifftshift
from scipy.optimize import linear_sum_assignment

from .farfield_fft import symmetric_pad

Z0 = np.sqrt(mu_0 / epsilon_0)


def weighted_projection(basis, field, weights):
    """Complex least squares; coefficients retain the supplied basis gauge."""
    u, f, w = np.asarray(basis,complex), np.asarray(field,complex), np.asarray(weights,float)
    if u.ndim != 2 or f.shape != (u.shape[0],) or w.shape != f.shape:
        raise ValueError('Incompatible basis, field, and quadrature shapes')
    if not all(np.isfinite(v).all() for v in (u,f,w)) or np.any(w <= 0):
        raise ValueError('Finite fields and strictly positive quadrature weights required')
    norms = np.sqrt(np.sum(w[:,None]*abs(u)**2,axis=0))
    fnorm = np.sqrt(np.sum(w*abs(f)**2))
    if fnorm == 0 or np.any(norms == 0):
        raise ValueError('Zero field or zero basis column')
    normalized = u/norms
    a,_,rank,singular = np.linalg.lstsq(np.sqrt(w)[:,None]*normalized,np.sqrt(w)*f,rcond=None)
    if rank != u.shape[1]:
        raise ValueError('Rank-deficient reference dictionary')
    a = a/norms
    reconstruction = u@a
    residual = float(np.sqrt(np.sum(w*abs(f-reconstruction)**2))/fnorm)
    return a,reconstruction,residual,float((singular[0]/singular[-1])**2)


def radiation_budget(coefficients, response):
    """T has final axis [Ex,Ey]; all leading axes enumerate basis terms."""
    a,t = np.asarray(coefficients,complex),np.asarray(response,complex)
    if t.shape != a.shape+(2,) or not np.isfinite(a).all() or not np.isfinite(t).all():
        raise ValueError('A and independent T must be finite and have matching basis axes')
    terms = (a[...,None]*t).reshape(-1,2)
    total = terms.sum(axis=0)
    incoherent = np.sum(abs(terms)**2,axis=0)
    coherent = abs(total)**2
    return terms,total,incoherent,coherent,coherent-incoherent


def volume_far_amplitude(electric, contrast, xyz_m, weights_m3, omega):
    """Normal F in E(r)=F exp(-ikr)/r, native +iwt, air background.

    Include the complete dielectric contrast. Electric final axes are (nodes,3).
    F has volts; the planar Fourier amplitude instead has V*m.
    """
    e,chi,r,w = (np.asarray(electric,complex),np.asarray(contrast,complex),
                 np.asarray(xyz_m,float),np.asarray(weights_m3,float))
    if e.ndim<2 or e.shape[-1]!=3 or r.shape!=(e.shape[-2],3) or chi.shape!=w.shape or w.shape!=(len(r),):
        raise ValueError('Volume fields, contrast, coordinates and weights must share nodes')
    if not all(np.isfinite(x).all() for x in [e,chi,r,w]) or np.any(w<=0) or not np.isfinite(omega) or abs(omega)==0:
        raise ValueError('Finite volume inputs, positive quadrature and nonzero omega required')
    k=omega/299792458.
    kernel=w*chi*np.exp(1j*k*r[:,2])
    return k*k/(4*np.pi)*np.einsum('...np,n->...p',e[...,:2],kernel)


def planar_K0_to_far(amplitude_V_m,z_m,omega):
    """Outgoing air angular spectrum -> normal spherical F (V)."""
    k=omega/299792458.
    return 1j*k/(2*np.pi)*np.asarray(amplitude_V_m)*np.exp(1j*k*z_m)


def pairs(values):
    return np.stack([np.real(values),np.imag(values)],axis=-1).tolist()


def forward(field, dx, size):
    """Physical +i q.r transform, input shape (y,x), SI spacing."""
    return fftshift(ifft2(ifftshift(symmetric_pad(field, size)))) * size**2 * dx**2


def inverse(amplitude, dx):
    """Inverse of forward: -i q.r, output on the full padded grid."""
    return fftshift(fft2(ifftshift(amplitude))) / (amplitude.shape[0]**2 * dx**2)


def backprop(amplitude, kz, distance, mask, dx):
    return inverse(amplitude * mask * np.exp(1j * kz * distance), dx)


def energy_group_gram(gram, right, total_energy, coefficients):
    """Groups use one joint-vector fit; last group is the exact field residual."""
    grouped = coefficients.conj() @ gram @ coefficients.T
    overlap = coefficients.conj() @ right
    result = np.empty((len(coefficients) + 1,) * 2, complex)
    result[:-1, :-1] = grouped
    result[:-1, -1] = overlap - grouped.sum(axis=1)
    result[-1, :-1] = result[:-1, -1].conj()
    result[-1, -1] = total_energy - 2 * overlap.sum().real + grouped.sum().real
    return result


def budget_rows(names, amplitudes, factor, scope):
    rows = []
    for pol, channel in enumerate(("Ex", "Ey")):
        for i, name in enumerate(names):
            rows.append(dict(scope=scope, channel=channel, kind="self", first=name,
                             second=name, value=float(factor * abs(amplitudes[i, pol])**2)))
            for j in range(i):
                rows.append(dict(scope=scope, channel=channel, kind="interference", first=names[j],
                                 second=name, value=float(2*factor*np.real(amplitudes[j, pol].conj()*amplitudes[i, pol]))))
        rows.append(dict(scope=scope, channel=channel, kind="total", first="all",
                         second="all", value=float(factor * abs(amplitudes[:, pol].sum())**2)))
    return rows


def periodic_profiles(hz, rho_um, k_um):
    return np.asarray(hz) * np.exp(1j * (np.asarray(rho_um) @ np.asarray(k_um)))


def mode_assignment(reference, current):
    a = reference / np.linalg.norm(reference, axis=1)[:, None]
    b = current / np.linalg.norm(current, axis=1)[:, None]
    overlap = abs(a.conj() @ b.T)
    rows, columns = linear_sum_assignment(-overlap)
    if not np.array_equal(rows, np.arange(len(reference))):
        raise ValueError('Incomplete mode assignment')
    return columns, overlap


def plane_farfield(electric,magnetic,axis_um,z_um,q,k0):
    """Common upward E/H angular operator, positive spatial transform, z=0 reference."""
    dx=float(axis_um[1]-axis_um[0])*1e-6
    transform=np.exp(1j*q[:,None]*axis_um[None,:]*1e-6)*dx
    ef=np.stack([transform@electric[...,i]@transform.T for i in range(3)],-1)
    hf=np.stack([transform@magnetic[...,i]@transform.T for i in range(3)],-1)
    qx,qy=np.meshgrid(q,q);kz=np.sqrt(k0*k0-qx*qx-qy*qy);nt=np.stack([qx,qy],-1)/k0;cosine=kz/k0
    hp=np.stack([hf[...,1],-hf[...,0]],-1)
    correction=Z0/cosine[...,None]*(hp-nt*np.sum(nt*hp,-1)[...,None])
    up=np.zeros_like(ef);up[...,:2]=(ef[...,:2]+correction)/2;up[...,2]=-np.sum(nt*up[...,:2],-1)/cosine
    return up*(1j*kz/(2*np.pi)*np.exp(1j*kz*z_um*1e-6))[...,None]
