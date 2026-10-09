"""Complete-cell complex decomposition and seam diagnostics, without solver imports."""
import numpy as np
from .hex_lattice_utils import unit_cell_corners


def complex_split(v):
    gamma=v.mean(axis=0);rest=v-gamma
    denom=max(float(np.linalg.norm(v)),np.finfo(float).tiny)
    checks=dict(reconstruction_relative=float(np.linalg.norm(v-gamma-rest)/denom),rest_mean_relative=float(np.linalg.norm(rest.mean(axis=0))/(denom/np.sqrt(len(v)))))
    return gamma,rest,checks


def shared_node_projection(values,mapping,counts):
    result=np.zeros((len(counts),values.shape[-1]),dtype=values.dtype)
    np.add.at(result,mapping,values)
    return result/counts[:,None]


def geometry_mapping(centers,a):
    corners=unit_cell_corners(a);edges={}
    for cell,center in enumerate(centers):
        for edge in range(6):
            xy=np.array([corners[edge],corners[(edge+1)%6]])+center
            key=tuple(sorted(map(tuple,np.round(xy,10))))
            edges.setdefault(key,[]).append([cell,edge])
    assert all(len(v) in (1,2) for v in edges.values())
    boundary=np.array([v[0] for v in edges.values() if len(v)==1],int)
    internal=np.array([v for v in edges.values() if len(v)==2],int)
    return corners,boundary,internal


def seam_metrics(values,corners):
    # Cap seams share the edge tangent. Check both this necessary trace condition
    # and the complete in-plane E jump; normal-to-seam components are kept separate.
    gamma=values.mean(axis=0);records=[]
    for edge in range(3):
        tangent=corners[(edge+1)%6]-corners[edge];tangent/=np.linalg.norm(tangent)
        left=gamma[edge,:,:2];right=gamma[edge+3,::-1,:2]
        l=left@tangent;r=right@tangent;delta=l-r
        scale=np.sqrt(np.mean((abs(l)**2+abs(r)**2)/2))
        records.append(dict(edge=edge,opposite_edge=edge+3,tangential_jump_rms_V_m=float(np.sqrt(np.mean(abs(delta)**2))),trace_rms_V_m=float(scale),relative_tangential_jump=float(np.linalg.norm(delta)/max(np.sqrt((np.linalg.norm(l)**2+np.linalg.norm(r)**2)/2),1e-300)),max_tangential_jump_V_m=float(abs(delta).max()),inplane_jump_relative=float(np.linalg.norm(left-right)/max(np.sqrt((np.linalg.norm(left)**2+np.linalg.norm(right)**2)/2),1e-300))))
    return gamma,records
