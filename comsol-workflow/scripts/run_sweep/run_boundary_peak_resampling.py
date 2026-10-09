"""Add real same-mesh Gamma samples around the measured S4 radiation minima."""
from __future__ import annotations
import argparse
from copy import deepcopy
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scripts.run_sweep import run_boundary_full_refined_scan as full

s4=full.s4
ARCHIVE=full.ARCHIVE
SOURCE_ROOT=full.ROOT
GROUP="12_px_py_peak_refine"
ROOT=s4.ROOT/"scripts/.out/unit_cell_2D/S4_b245_eta0.96_peak-refine_mesh5_h0.25_edge128_px-py_20260930"


def sampling_plan(df):
    existing={round(float(z),5) for z in df.zeta.unique()}
    windows={};points=set()
    for mode in ("px","py"):
        d=df[(df.user_mode==mode)&df.mode_status.eq("single_eigenmode")].sort_values("zeta")
        if d.empty or not np.isfinite(d[["zeta","Q","air_display_V_per_m"]]).all().all() or (d.Q<=0).any():
            raise ValueError("Need finite identified source modes and positive Q")
        peak=float(d.loc[d.Q.idxmax(),"zeta"])
        crossings=[]
        for left,right in zip(d.itertuples(),list(d.itertuples())[1:]):
            if left.air_display_V_per_m*right.air_display_V_per_m<0:
                zero=left.zeta-left.air_display_V_per_m*(right.zeta-left.zeta)/(right.air_display_V_per_m-left.air_display_V_per_m)
                crossings.append((zero,left.zeta,right.zeta))
        if not crossings:raise ValueError(f"No measured sign bracket for {mode}")
        estimate,left,right=min(crossings,key=lambda item:abs(item[0]-peak))
        if abs(estimate-peak)>.0005:raise ValueError("Measured sign bracket needs a wider reviewed window")
        center=round(float(estimate),4)
        coarse={round(peak+i*.0001,5) for i in range(-5,6)}
        fine={round(center+i*.00002,5) for i in range(-5,6)}
        desired=coarse|fine
        if not all(.8<=z<=1.2 for z in desired):raise ValueError("Sampling outside authorized zeta interval")
        points.update(desired-existing)
        windows[mode]=dict(sampled_Q_peak=peak,measured_sign_bracket=[left,right],coordinate_only_zero_estimate=float(estimate),
            coarse_range=[min(coarse),max(coarse)],coarse_step=.0001,fine_range=[min(fine),max(fine)],fine_step=.00002,
            new_points=sorted(desired-existing))
    return windows,sorted(points)


def prepare():
    source_path=SOURCE_ROOT/"99_config/config.json"
    source=s4.json_read(source_path)
    table=ARCHIVE/"80_logs/11_px_py_mesh025/scan_results.csv"
    df=pd.read_csv(table);windows,zetas=sampling_plan(df)
    if len(zetas)>40:raise ValueError("Unexpected local sampling count")
    expected=sorted([c["zeta"] for c in source["cases"]]+zetas)
    config_path=ROOT/"99_config/config.json"
    if config_path.exists():
        config=s4.json_read(config_path)
        full.validate_config(config,expected)
        if config.get("sampling_windows")!=windows:raise ValueError("Sampling plan changed")
        return config
    from scripts.run_main import run_band_pair as band
    from comsol_workflow.geometry_utils import create_hexagon_design
    hashes={}
    def protect(path):
        path=Path(path);hashes[str(path)]=s4.digest(path)
    for path in [source_path,table,source["reference_fields"],s4.ROOT/"scripts/parameter.json"]:protect(path)
    cases=[]
    for original in source["cases"]:
        case=dict(original,reuse_export=True);cases.append(case)
        folder=Path(case["directory"])
        for name in ["export.json","fields.npz","mode_results.csv","native_components.csv","99_config/parameter.json","99_config/realized_geometry.json","99_config/mesh_verified.json"]:protect(folder/name)
    templates={mode:next(c for c in source["cases"] if c["zeta"]==w["sampled_Q_peak"]) for mode,w in windows.items()}
    for z in zetas:
        mode=min(windows,key=lambda m:abs(z-windows[m]["sampled_Q_peak"]))
        template=templates[mode];params=s4.json_read(template["snapshot"])
        params["unit_cell_2d"]["zeta"]=z
        folder=ROOT/f"zeta{z:g}";snapshot=folder/"99_config/parameter.json"
        s4.write_json(snapshot,params);parsed=s4.load_shared_parameters(snapshot)
        hp=band.get_hole_params(parsed.unit_cell_2d.cell.to_fourier_params("dual_boundary_scan"))
        outer,holes,info=create_hexagon_design(.82,hp)
        if min(v["min_dist"] for v in info)<band.D:raise ValueError("Geometry clearance failed")
        protect(template["coarse_model"])
        cases.append(dict(zeta=z,directory=str(folder),snapshot=str(snapshot),snapshot_sha256=s4.digest(snapshot),
            coarse_model=template["coarse_model"],geometry_holes_um=np.asarray(holes,float).tolist()))
    config={k:deepcopy(v) for k,v in source.items() if k not in ["cases","source_hashes"]}
    config.update(cases=sorted(cases,key=lambda c:c["zeta"]),source_hashes=hashes,sampling_windows=windows,
        previous_plot_group="11_px_py_mesh025",reused_models=len(source["cases"]),new_solves=len(zetas),
        shared_parameter_sha256=s4.digest(s4.ROOT/"scripts/parameter.json"))
    full.validate_config(config,expected)
    s4.write_json(config_path,config);s4.write_json(ARCHIVE/"99_config"/GROUP/"config.json",config)
    return config


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--prepare-only",action="store_true")
    action.add_argument("--run",action="store_true")
    action.add_argument("--plot-only",action="store_true")
    parser.add_argument("--allow-concurrent",action="store_true")
    args=parser.parse_args();config=prepare()
    print(f"{len(config['cases'])} total geometries; {config['new_solves']} new solves; {config['reused_models']} read-only exports; 4 cores; mesh5 h/h0=0.25",flush=True)
    print(config["sampling_windows"],flush=True)
    full.ROOT=ROOT;full.GROUP=GROUP
    code=ARCHIVE/"99_config"/GROUP/"code";code.mkdir(parents=True,exist_ok=True);shutil.copy2(__file__,code/Path(__file__).name)
    if args.run:full.run(config,args.allow_concurrent)
    elif args.plot_only:full.finish(config)


if __name__=="__main__":main()
