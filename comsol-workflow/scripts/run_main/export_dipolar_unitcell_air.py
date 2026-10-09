"""Read saved solved unit-cell models; export air E/H without solving or saving MPH."""
import json
import numpy as np
from scripts.analysis import analyze_dipolar_complete as a
from scripts.run_main.run_dipolar_data_completion import EHS, fingerprint
from scripts.run_main.run_dipolar_volume_export import NativeSampler
from scripts.run_main.run_boundary_analysis import runtime_preflight
from scripts.run_sweep.run_boundary_surface_zeta_scan import frequencies


def run():
    audit={"action":"read retained solutions and interpolate only", "preflight":runtime_preflight(),"models":[]}
    previous=a.L/"Fig6_unitcell_air_export.json"
    if previous.exists():audit["models"]=json.loads(previous.read_text(encoding="utf8")).get("models",[])
    import mph
    from jpype import JClass
    client=mph.start(version="6.3",cores=4)
    if not JClass("com.comsol.model.util.ModelUtil").checkoutLicense("COMSOL","WAVEOPTICS"):
        raise RuntimeError("COMSOL/WAVEOPTICS license unavailable")
    a.write_json("Fig6_unitcell_air_export.json",audit)
    xyz=a.load("periodic_k00_z0_fine_EH.npz")["center_xyz_um"].copy();xyz[:,2]=1.1
    try:
        for i in range(19):
            path=a.WORK/"00_model"/f"periodic_k{i:02d}.mph"
            output=a.D/f"periodic_k{i:02d}_air1100nm_EH.npz"
            before=fingerprint(path);old=a.load(f"periodic_k{i:02d}_EH.npz")
            if output.exists():
                with np.load(output) as z:assert str(z["model_sha256"])==before["sha256"]
                print(f"k{i:02d}: reuse verified air export",flush=True);continue
            print(f"k{i:02d}: loading retained py solution",flush=True)
            model=client.load(str(path))
            try:
                assigned=old["mode_assignment"]
                np.testing.assert_allclose(frequencies(model)[["re","im","q"]].to_numpy()[assigned],old["frequencies_THz"],rtol=1e-9)
                np.testing.assert_allclose([model.java.param().evaluate(v)*1e-6 for v in ['kx','ky']],old['k_um_inv'],atol=1e-12)
                sampler=NativeSampler(model);mode=int(assigned[1]);phase=old['phase_to_raw_reference'][1]
                check=sampler(mode,old['center_xyz_um'][::37],EHS)*phase
                expected=np.concatenate([old['E_center_V_m'][1,::37],old['H_center_A_m'][1,::37]],-1)
                error=float(np.linalg.norm(check-expected)/np.linalg.norm(expected));assert error<1e-9
                values=sampler(mode,xyz,EHS)*phase
                np.savez_compressed(output,xyz_um=xyz,E_V_m=values[:,:3],H_A_m=values[:,3:],k_um_inv=old['k_um_inv'],cell_EM_energy_J=old['cell_EM_energy_J'][1],frequency_THz=old['frequencies_THz'][1],phase_to_raw_reference=phase,model_sha256=before['sha256'],py_mode=mode)
                assert fingerprint(path)==before
                audit['models'].append(dict(index=i,identity=before,z0_check_relative_error=error,air_samples=len(xyz),output=str(output)))
                a.write_json("Fig6_unitcell_air_export.json",audit)
                print(f"k{i:02d}: exported {len(xyz)} air samples; z0 error={error:.3g}",flush=True)
            finally:client.remove(model)
        audit['status']='complete';a.write_json("Fig6_unitcell_air_export.json",audit)
    finally:client.clear()


if __name__=='__main__':run()
