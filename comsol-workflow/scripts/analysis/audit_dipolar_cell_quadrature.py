"""Audit cell quadrature on saved quarter FEM fields; never starts COMSOL.

Compare point binning with exact cell clipping of a per-tetrahedron affine
source reconstructed from the same four native Gauss samples. This is an
integration diagnostic, not a claim of convergence of the original FEM field.
"""
from __future__ import annotations

from comsol_workflow.output_paths import table_path
import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from scripts.analysis import analyze_dipolar_complete as a
from scripts.analysis.analyze_dipolar_cell_radiation import grouped_sum
from scripts.analysis.analyze_dipolar_singularity import digest

BETA = (5 - np.sqrt(5)) / 20
BARY = np.full((4, 4), BETA)
np.fill_diagonal(BARY, 1 - 3 * BETA)
BARY_INV = np.linalg.inv(BARY)
NORMALS = np.column_stack([np.cos(np.arange(6)*np.pi/3),
                           np.sin(np.arange(6)*np.pi/3), np.zeros(6)])


def clipped_tetrahedron_moment(vertices, normals, offsets):
    """Return integral of the original tetrahedron's four barycentric coords.

    A convex clipped polyhedron is triangulated into pyramids about an interior
    point. Volume and affine moments are exact up to floating-point geometry.
    Coordinates and output volume use um and um^3 respectively.
    """
    identity = np.eye(4)
    faces = [identity[index] for index in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]]
    tol = 1e-12
    for normal, offset in zip(normals, offsets):
        vertex_distance = vertices @ normal - offset
        all_distances = np.concatenate(faces) @ vertex_distance
        if all_distances.max() <= tol:
            continue
        if all_distances.min() >= -tol:
            return np.zeros(4)
        new_faces, intersections = [], []
        for face in faces:
            distances = face @ vertex_distance
            output = []
            for i, current in enumerate(face):
                previous = face[i-1]
                dc, dp = distances[i], distances[i-1]
                inc, inp = dc <= tol, dp <= tol
                if inc != inp:
                    ratio = np.clip(dp / (dp-dc), 0., 1.)
                    crossing = previous + ratio * (current-previous)
                    output.append(crossing)
                    intersections.append(crossing)
                if inc:
                    output.append(current)
            if len(output) >= 3:
                new_faces.append(np.array(output))
        if intersections:
            cap = np.array(intersections)
            _, unique = np.unique(np.round(cap, 12), axis=0, return_index=True)
            cap = cap[np.sort(unique)]
            if len(cap) >= 3:
                xyz = cap @ vertices
                xyz -= xyz.mean(0)
                # All production clipping planes are vertical. The synthetic
                # check also exercises a general nonvertical clipping plane.
                auxiliary = np.eye(3)[np.argmin(abs(normal))]
                u = np.cross(normal, auxiliary)
                u /= np.linalg.norm(u)
                v = np.cross(normal, u)
                cap = cap[np.argsort(np.arctan2(xyz @ v, xyz @ u))]
                new_faces.append(cap)
        faces = new_faces
        if not faces:
            return np.zeros(4)
    middle = np.concatenate(faces).mean(0)
    result = np.zeros(4)
    for face in faces:
        for i in range(1,len(face)-1):
            tri = face[[0,i,i+1]]
            volume = abs(np.linalg.det((tri-middle) @ vertices)) / 6
            result += volume * (middle+tri.sum(0)) / 4
    return result


def check():
    tetra = np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    full = clipped_tetrahedron_moment(tetra,np.empty((0,3)),np.empty(0))
    np.testing.assert_allclose(full,np.ones(4)/24,rtol=1e-13)
    low = clipped_tetrahedron_moment(tetra,np.array([[1.,0.,0.]]),np.array([.5]))
    high = clipped_tetrahedron_moment(tetra,np.array([[-1.,0.,0.]]),np.array([-.5]))
    np.testing.assert_allclose(low+high,full,rtol=1e-13)
    np.testing.assert_allclose(high.sum(),1/48,rtol=1e-13)
    point_bin_volume=np.count_nonzero((BARY @ tetra)[:,0]>.5)/24
    # Regression: node membership assigns twice the correct clipped volume.
    assert abs(point_bin_volume/high.sum()-2)<1e-12
    np.testing.assert_allclose(high @ tetra, np.array([.625,.125,.125])/48,rtol=1e-13)
    # Known affine source, exactly recoverable from the original four samples.
    fields = (BARY @ tetra) @ np.array([1+2j,3-1j,2j]) + .7
    np.testing.assert_allclose(high @ BARY_INV @ fields,
        (np.array([.625,.125,.125]) @ np.array([1+2j,3-1j,2j])+.7)/48,rtol=1e-13)
    # Oblique plane partitions, translations, and reversing cell-processing order.
    normal = np.array([1.,2.,3.])
    first = clipped_tetrahedron_moment(tetra,[normal],[.6])
    second = clipped_tetrahedron_moment(tetra,[-normal],[-.6])
    np.testing.assert_allclose(first+second,full,rtol=1e-13)
    shift = np.array([14.,17.,.03])
    moved = clipped_tetrahedron_moment(tetra+shift,[normal],[.6+normal@shift])
    np.testing.assert_allclose(moved,first,rtol=2e-12,atol=1e-16)
    print("PASS: clipped volume, affine source, oblique partitions, translation.",flush=True)


def run(all_cells=False):
    cfg = json.loads((a.CONFIG/'config.json').read_text(encoding='utf8'))
    points = cfg['bulk_points'] + cfg['cladding_points']
    centers = np.array([[p['x'],p['y']] for p in points])
    nc = len(cfg['bulk_points'])
    count = len(points)
    tree = cKDTree(centers)
    if all_cells:
        selected = np.flatnonzero((centers[:,0]>=-1e-10)&(centers[:,1]>=-1e-10))
    else:
        # Adjacent anomalous cells wholly inside the directly solved quadrant.
        # Add every cavity cell in the local rectangular patch for comparison.
        selected = np.flatnonzero((np.arange(count)<nc)&(centers[:,0]>1.)&(centers[:,0]<6.)
                                 &(centers[:,1]>1.4)&(centers[:,1]<6.4))
    selected_tree = cKDTree(centers[selected])
    print(f"Audit {len(selected)} directly sampled cells",flush=True)
    geometry = a.load('finite_mesh_geometry.npz')
    chosen = np.flatnonzero(np.isin(geometry['domain_ids'],geometry['slab_domains']))
    records = json.loads((a.L/'finite_mesh_exports.json').read_text(encoding='utf8'))['parts']
    integrals = np.zeros((count,10),complex)
    binned = np.zeros((count,2),complex)
    volume = np.zeros(count)
    old_volume = np.zeros(count)
    counts = dict(full_tetrahedra=0,cut_cell_tetrahedron_pairs=0,loaded_parts=0)
    start_time = time.monotonic()
    for part,record in enumerate(records):
        tet_ids = chosen[part*25000:part*25000+record['tetrahedra']]
        vertices = geometry['vertices_um'][geometry['tetrahedra'][tet_ids]]
        midpoint = vertices.mean(1)
        radius = np.linalg.norm(vertices[:,:,:2]-midpoint[:,None,:2],axis=-1).max(1)
        distances,_ = selected_tree.query(midpoint[:,:2])
        candidate = np.flatnonzero(distances <= radius+cfg['a']/np.sqrt(3)+1e-10)
        if not len(candidate):
            continue
        data = a.load(record['path'])
        np.testing.assert_array_equal(tet_ids,data['tetrahedron_ids'])
        counts['loaded_parts'] += 1
        xyz = data['xyz_um'].reshape(-1,4,3)
        w = data['weights_octant_m3'].reshape(-1,4)
        eps = data['epsilon_r'].reshape(-1,4)
        source = data['E_V_m'].reshape(-1,4,3)[:,:,:2] * (2*(eps-1)*np.cos(a.K0*xyz[:,:,2]*1e-6))[:,:,None]
        electric=data['E_V_m'].reshape(-1,4,3)
        magnetic=data['H_A_m'].reshape(-1,4,3)
        density=np.concatenate([a.epsilon_0*eps[:,:,None]*abs(electric)**2,
                                a.mu_0*abs(magnetic)**2],axis=-1)/4
        l1=2*(eps-1)[:,:,None]*abs(electric[:,:,:2])
        samples=np.concatenate([source,2*density,l1],axis=-1)
        for t in candidate:
            candidates = selected_tree.query_ball_point(midpoint[t,:2],radius[t]+cfg['a']/np.sqrt(3)+1e-10)
            for local_id in candidates:
                j = selected[local_id]
                offsets = cfg['a']/2 + NORMALS[:,:2] @ centers[j]
                distance = vertices[t] @ NORMALS.T - offsets
                if np.any(distance.min(0) > 1e-12):
                    continue
                if np.all(distance <= 1e-12):
                    # The unchanged native whole-tetrahedron quadrature.
                    integrals[j] += np.sum(samples[t]*w[t,:,None],axis=0)
                    volume[j] += 2*w[t].sum()
                    counts['full_tetrahedra'] += 1
                else:
                    moment = clipped_tetrahedron_moment(vertices[t],NORMALS,offsets)
                    if moment.sum() < 1e-22:
                        continue
                    integrals[j] += (moment @ BARY_INV @ samples[t])*1e-18
                    volume[j] += 2*moment.sum()*1e-18
                    counts['cut_cell_tetrahedron_pairs'] += 1
        # Original point-bin path evaluated on exactly these same native data.
        ids = tree.query(xyz[candidate,:,:2].reshape(-1,2))[1]
        relative=xyz[candidate,:,:2].reshape(-1,2)-centers[ids]
        inside=np.all(relative @ NORMALS[:,:2].T <= cfg['a']/2+1e-10,axis=1)
        local_sources = (source[candidate] * w[candidate,:,None]).reshape(-1,2)
        binned += grouped_sum(ids[inside],local_sources[inside],count)
        old_volume += np.bincount(ids[inside],weights=(2*w[candidate].ravel())[inside],minlength=count)
        print(f"part {part+1}/{len(records)}; cut pairs {counts['cut_cell_tetrahedron_pairs']}; elapsed {time.monotonic()-start_time:.1f}s",flush=True)
    # Selected cells can touch the symmetry axes; assemble their full cells.
    factor = a.K0**2/(4*np.pi)
    total=integrals[:,:2].copy()
    energy=integrals[:,2:8].real.copy()
    source_l1=integrals[:,8:10].real.copy()
    for j in selected:
        mult = 2**int(abs(centers[j,0])<1e-10) * 2**int(abs(centers[j,1])<1e-10)
        total[j] *= factor*mult
        binned[j] *= factor*mult
        volume[j] *= mult
        old_volume[j] *= mult
        energy[j] *= mult
        source_l1[j] *= mult
        if mult>1:
            total[j,0] = 0
            binned[j,0] = 0
    exact_volume = np.sqrt(3)/2*cfg['a']**2*.2e-18
    baseline = a.load('unitcell_center_radiation_pointbin.npz')
    np.testing.assert_allclose(binned[selected],baseline['F_center_V'][selected],rtol=5e-10,atol=1e-11)
    np.testing.assert_allclose(old_volume[selected],baseline['volume_m3'][selected],rtol=1e-12)
    max_error = float(abs(volume[selected]/exact_volume-1).max())
    assert max_error < 1e-7,max_error
    suffix = 'all' if all_cells else 'patch'
    np.savez_compressed(a.D/f'cell_quadrature_audit_{suffix}.npz',cell_ids=selected,
        centers_um=centers[selected],F_pointbin_V=binned[selected],F_clipped_affine_V=total[selected],
        volume_pointbin_m3=old_volume[selected],volume_clipped_m3=volume[selected],
        energy_components_J=energy[selected],source_L1_V_m2=source_l1[selected],
        exact_cell_volume_m3=exact_volume,field_model='affine reconstruction from four FEM Gauss samples; no new COMSOL values')
    rows=[]
    for j in selected:
        rows.append(dict(cell_id=int(j),x_um=centers[j,0],y_um=centers[j,1],
            old_volume_relative_error=old_volume[j]/exact_volume-1,new_volume_relative_error=volume[j]/exact_volume-1,
            old_Fx_re=binned[j,0].real,old_Fx_im=binned[j,0].imag,new_Fx_re=total[j,0].real,new_Fx_im=total[j,0].imag,
            old_Fy_re=binned[j,1].real,old_Fy_im=binned[j,1].imag,new_Fy_re=total[j,1].real,new_Fy_im=total[j,1].imag,
            old_Ey_self_W_sr=abs(binned[j,1])**2/(2*a.Z0),new_Ey_self_W_sr=abs(total[j,1])**2/(2*a.Z0)))
    pd.DataFrame(rows).to_csv(table_path(a.L, f'cell_quadrature_audit_{suffix}.csv'),index=False)
    details=dict(method='exact hexagonal clipping with affine source interpolation on boundary tetrahedra',
        COMSOL_started=False,selected_cells=len(selected),old_max_cell_volume_error=float(abs(old_volume[selected]/exact_volume-1).max()),
        clipped_max_cell_volume_error=max_error,counts=counts,elapsed_s=time.monotonic()-start_time,
        per_cell_source_converged=False,interpretation='diagnostic of integration partition; not a new independently converged FEM evaluation')
    if all_cells:
        # Reconstruct complete cells, not four copies of a partial axis cell.
        full_f=np.zeros((count+1,2),complex)
        full_u=np.zeros((count+1,6))
        full_l1=np.zeros((count+1,2))
        full_v=np.zeros(count+1)
        seen=np.zeros(count,bool)
        for sx,sy in [(1,1),(-1,1),(1,-1),(-1,-1)]:
            distance,target=tree.query(centers[selected]*(sx,sy))
            assert distance.max()<1e-10
            full_f[target]=total[selected]*(sx*sy,1)
            full_u[target]=energy[selected]
            full_l1[target]=source_l1[selected]
            full_v[target]=volume[selected]
            seen[target]=True
        assert seen.all()
        # Conservative remainder of the same affine tetrahedron fields.
        # Partition identity, never calibration to the independent air field.
        full_f[-1]=baseline['F_center_V'].sum(0)-full_f[:-1].sum(0)
        full_u[-1]=baseline['energy_components_J'].sum(0)-full_u[:-1].sum(0)
        full_l1[-1]=baseline['source_L1_V_m2'].sum(0)-full_l1[:-1].sum(0)
        full_v[-1]=baseline['volume_m3'].sum()-full_v[:-1].sum()
        assert np.all(full_u>=0) and np.all(full_l1>=0) and np.all(full_v>0)
        details['config_sha256']=digest(a.CONFIG/'config.json')
        details['baseline_sha256']=digest(a.D/'unitcell_center_radiation_pointbin.npz')
        details['raw_mesh_inputs']={record['path']:digest(a.D/record['path']) for record in records}
        np.savez_compressed(a.D/'unitcell_clipped_integrals.npz',F_center_V=full_f,
            energy_components_J=full_u,source_L1_V_m2=full_l1,volume_m3=full_v,
            baseline_F_total_V=baseline['F_center_V'].sum(0),
            baseline_U_components_J=baseline['energy_components_J'].sum(0),
            config_sha256=details['config_sha256'],cell_centers_um=centers)
    a.write_json(f'cell_quadrature_audit_{suffix}.json',details)
    print(json.dumps(details,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--all-cells',action='store_true')
    args=parser.parse_args()
    check()
    if not args.check:
        run(args.all_cells)
