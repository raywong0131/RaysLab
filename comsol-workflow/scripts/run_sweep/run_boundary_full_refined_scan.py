"""Recompute the delivered S4 scan at one verified mesh size, then plot offline."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import shutil
import time
import traceback
import numpy as np
import pandas as pd
from scripts.run_sweep import run_boundary_dual_validation as dual
from scripts.run_sweep.run_boundary_mesh_convergence import mesh_stats, scale_mesh_sizes

s4=dual.s4
ARCHIVE=dual.ARCHIVE
GROUP="11_px_py_mesh025"
ROOT=s4.ROOT/"scripts/.out/unit_cell_2D/S4_b245_eta0.96_zeta0.8-1.2_mesh5_h0.25_edge128_px-py_20260929"
SIZE_FACTOR=.25


def validate_config(config, zetas):
    expected=dict(version=1,mesh=5,size_factor=SIZE_FACTOR,cores=4,eta=.96,b0_nm=245,a_nm=820,H_nm=200,n=3.3,edge_order=128,z_nm=100,Gamma=True,expected_eigenvalues=4,shift="c_const/1.55[um]",user_modes=["px","py"],internal_modes=["py","px"])
    if any(config.get(k)!=v for k,v in expected.items()) or [c["zeta"] for c in config["cases"]]!=list(zetas):
        raise ValueError("Existing config identity differs from the requested refined scan")
    for p,h in config["source_hashes"].items():
        if s4.digest(p)!=h:raise ValueError(f"Protected source changed: {p}")
    for c in config["cases"]:
        if s4.digest(c["snapshot"])!=c["snapshot_sha256"]:raise ValueError("Case snapshot changed")


def prepare():
    old_config=ARCHIVE/"99_config/09_px_py_validation/config.json"
    old=s4.json_read(old_config)
    table_path=ARCHIVE/"80_logs/09_px_py_validation/scan_results.csv"
    df=pd.read_csv(table_path)
    zetas=sorted(df.zeta.unique().tolist())
    if len(zetas)!=81 or zetas!=sorted(c["zeta"] for c in old["cases"]):raise ValueError("Original scan coverage changed")
    config_path=ROOT/"99_config/config.json"
    if config_path.exists():
        config=s4.json_read(config_path);validate_config(config,zetas);return config
    from scripts.run_sweep.run_boundary_px_mesh_test import OUTPUT as test_root
    cached={};sources={}
    def protect(p):
        p=Path(p);sources[str(p)]=s4.digest(p);return str(p)
    for p in test_root.rglob("summary.json"):
        record=s4.json_read(p);r=record.get("result",{})
        if r.get("mesh")==5 and r.get("size_factor")==SIZE_FACTOR:
            model=p.parent/"s4_gamma.mph"
            if s4.digest(model)!=record["output_hashes"][model.name]:raise ValueError("Fine model cache hash differs")
            if r["zeta"] in cached:raise ValueError("Duplicate fine-mesh cache")
            cached[r["zeta"]]=protect(model);protect(p)
    cases=[]
    for z in zetas:
        row=df[(df.user_mode=="px")&(df.zeta==z)].iloc[0]
        source_snapshot=ARCHIVE/f"99_config/09_px_py_validation/zeta{z:g}/parameter.json"
        params=s4.json_read(source_snapshot);protect(source_snapshot)
        if params["mesh_size"]!=5 or [params["unit_cell_2d"][k] for k in ("b0_nm","eta","zeta")]!=[245.,.96,z]:raise ValueError("Original parameter mismatch")
        folder=ROOT/f"zeta{z:g}";snapshot=folder/"99_config/parameter.json"
        s4.write_json(snapshot,params);s4.load_shared_parameters(snapshot)
        case=dict(zeta=z,directory=str(folder),snapshot=str(snapshot),snapshot_sha256=s4.digest(snapshot),coarse_model=protect(row.source_model),source_fields=protect(row.field_file))
        if z in cached:case["refined_model"]=cached[z]
        cases.append(case)
    protect(old_config);protect(table_path);protect(old["reference_fields"]);protect(s4.ROOT/"scripts/parameter.json")
    config=dict(version=1,mesh=5,size_factor=SIZE_FACTOR,cores=4,eta=.96,b0_nm=245,a_nm=820,H_nm=200,n=3.3,edge_order=128,z_nm=100,Gamma=True,expected_eigenvalues=4,shift="c_const/1.55[um]",user_modes=["px","py"],internal_modes=["py","px"],cases=cases,source_hashes=sources,reference_fields=old["reference_fields"],shared_parameter_sha256=s4.digest(s4.ROOT/"scripts/parameter.json"),reused_models=len(cached),new_solves=len(cases)-len(cached),display_boundary_sign=-1)
    validate_config(config,zetas);s4.write_json(config_path,config)
    s4.write_json(ARCHIVE/"99_config"/GROUP/"config.json",config)
    return config


def check_hole_volumes(actual,expected):
    """Allow one ppm in CAD volume measurement, while rejecting incorrect domains."""
    actual,expected=np.asarray(actual,float),np.asarray(expected,float)
    if actual.shape!=(6,) or expected.shape!=(6,) or not np.isfinite([actual,expected]).all() or np.any(expected<=0):
        raise ValueError("Need six finite positive hole volumes")
    np.testing.assert_allclose(actual,expected,rtol=1e-6,atol=1e-14)
    return dict(actual_um3=actual.tolist(),expected_um3=expected.tolist(),
                relative_errors=(abs(actual-expected)/expected).tolist(),rtol=1e-6,atol_um3=1e-14)


def refined_model(client,case,record,status_path):
    folder=Path(case["directory"]);checkpoint=folder/"s4_gamma.mph";done=folder/"solve.json"
    if "refined_model" in case:
        return Path(case["refined_model"])
    if done.exists():
        saved=s4.json_read(done)
        if saved["size_factor"]!=SIZE_FACTOR or s4.digest(checkpoint)!=saved["model_sha256"]:raise ValueError("Refined checkpoint identity differs")
        return checkpoint
    if checkpoint.exists() or (folder/"solve_attempt.json").exists():raise RuntimeError(f"Interrupted solve needs inspection: {folder}")
    model=client.load(case["coarse_model"])
    try:
        np.testing.assert_allclose([model.java.param().evaluate(k) for k in ("a","H","kx","ky")],[820e-9,200e-9,0,0],rtol=1e-12,atol=1e-18)
        mesh=model.java.component("comp1").mesh("mesh1")
        if not mesh.isAutomatic() or int(mesh.autoMeshSize())!=5:raise ValueError("Expected original automatic mesh5")
        before=mesh_stats(model)
        s4.write_json(folder/"solve_attempt.json",dict(started=time.time(),source_model=case["coarse_model"],size_factor=SIZE_FACTOR))
        model.java.sol("sol1").clearSolutionData()
        record.update(status="checking_geometry",current_zeta=case["zeta"]);s4.write_json(status_path,record)
        if "geometry_holes_um" in case:
            from jpype.types import JArray, JString, JInt
            from comsol_workflow.boundary_integrals import polygon_area
            comp=model.java.component("comp1");geom=comp.geom("geom1")
            holes=np.asarray(case["geometry_holes_um"],float)
            if holes.shape!=(6,3,2) or not np.isfinite(holes).all():raise ValueError("Invalid new hole geometry")
            for i,hole in enumerate(holes):
                geom.feature("wp1").geom().feature(f"pol{i+2}").set("table",JArray(JString,2)([[f"{x:.17g} [um]",f"{y:.17g} [um]"] for x,y in hole]))
                selection=comp.selection(f"sel_hole_dom_{i}")
                selection.set("posx",float(hole[:,0].mean()));selection.set("posy",float(hole[:,1].mean()))
            geom.run()
            measured=[];expected=[]
            for i,hole in enumerate(holes):
                table=geom.feature("wp1").geom().feature(f"pol{i+2}").getStringMatrix("table")
                vertices=np.array([[float(str(value).split()[0]) for value in row] for row in table])
                np.testing.assert_allclose(vertices,hole,rtol=1e-12,atol=1e-14)
                domains=list(comp.selection(f"sel_hole_dom_{i}").entities(3))
                if len(domains)!=1:raise ValueError("New hole selection is ambiguous")
                measure=geom.measureFinal();measure.selection().geom("geom1",3);measure.selection().set(JArray(JInt,1)(domains))
                measured.append(float(measure.getVolume()));expected.append(abs(polygon_area(hole))*.1)
            volumes=check_hole_volumes(measured,expected)
            s4.write_json(folder/"99_config/geometry_verified.json",dict(zeta=case["zeta"],hole_volumes_match=True,holes_um=holes.tolist(),
                polygon_vertices_match=True,volume_check=volumes,geometry_repair_tolerance_type=str(geom.repairTolType()),
                geometry_relative_repair_tolerance=float(geom.repairTol()),geometry_absolute_repair_tolerance_um=float(geom.absRepairTol())))
        changes=scale_mesh_sizes(mesh,SIZE_FACTOR);mesh.clearMesh()
        s4.write_json(folder/"99_config/size_changes.json",changes)
        record.update(status="meshing",current_zeta=case["zeta"]);s4.write_json(status_path,record)
        client.java.showProgress(str(folder/"comsol_progress.log"));print("MESH",case["zeta"],flush=True)
        mesh.run();after=mesh_stats(model)
        if after["automatic"] or after["slab_layer_tetrahedra"]<=1.25*before["slab_layer_tetrahedra"]:raise ValueError("Slab was not refined")
        s4.write_json(folder/"99_config/mesh.json",dict(before=before,after=after,size_factor=SIZE_FACTOR))
        record.update(status="solving",mesh=after);s4.write_json(status_path,record);print("SOLVE",case["zeta"],after,flush=True)
        start=time.monotonic();model.java.sol("sol1").runAll();model.save(str(checkpoint))
        s4.write_json(done,dict(size_factor=SIZE_FACTOR,model_sha256=s4.digest(checkpoint),elapsed_seconds=time.monotonic()-start))
        return checkpoint
    finally:client.remove(model)


def finish(config):
    from scripts.analysis import plot_boundary_dual_validation as plotting
    df=dual.analyze(config,root=ROOT)
    if len(df)!=2*len(config["cases"]) or df.duplicated(["zeta","user_mode"]).any():raise ValueError("Full dual-mode coverage missing")
    edges=pd.read_csv(ROOT/"80_logs/edge_results.csv")
    plotting.validate_saved(config,df,edges,root=ROOT)
    summary=plotting.summarize(df)
    for c in config["cases"]:
        folder=Path(c["directory"]);name=folder.name
        target=ARCHIVE/"01_results"/GROUP/name/"fields.npz";target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(folder/"fields.npz",target);df.loc[df.zeta==c["zeta"],"field_file"]=str(target)
        for p in [folder/"export.json",folder/"99_config/parameter.json",folder/"99_config/realized_geometry.json",folder/"99_config/mesh_verified.json"]:
            dest=ARCHIVE/"99_config"/GROUP/name/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        for p in folder.glob("*.csv"):
            dest=ARCHIVE/"80_logs"/GROUP/name/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    logs=ARCHIVE/"80_logs"/GROUP;logs.mkdir(parents=True,exist_ok=True)
    df.to_csv(logs/"scan_results.csv",index=False);edges.to_csv(logs/"edge_results.csv",index=False)
    edges.groupby(["zeta","user_mode","hole"],as_index=False).display_Vm.sum().to_csv(logs/"hole_results.csv",index=False)
    for name in ("verification.json","supplemental_complex_components.csv"):shutil.copy2(ROOT/"80_logs"/name,logs/name)
    s4.write_json(logs/"scientific_summary.json",summary)
    # This process renders one isolated delivery group; historical figures stay intact.
    previous=ARCHIVE/"90_history/10_overview"/config.get("previous_plot_group","09_px_py_validation")
    for name in plotting.NAMES:
        current=ARCHIVE/"10_overview"/f"{name}.png"
        if current.exists() and not (previous/current.name).exists():
            previous.mkdir(parents=True,exist_ok=True);shutil.copy2(current,previous/current.name)
    plotting.GROUP=GROUP
    audit=plotting.figures(df,edges,summary,ARCHIVE,mesh_comparison=False)
    s4.write_json(ARCHIVE/"99_config"/GROUP/"plot_contract.json",audit)
    lines=["# S4 全扫描统一加密网格", "", f"固定η=0.96、ζ=0.8–1.2，共{len(config['cases'])}个几何、{len(df)}条模式记录。全部使用基础mesh=5、h/h0=0.25的自定义网格；每孔边128点积分，板层侧z=100 nm。", "", "四图沿用原计算与显示约定。Fig.2全部来自同一网格设置，不叠加粗网格。Q为原始本征频率给出的值，图中纵轴上限为10^8，最高点实际数值见下表；图2以竖直虚线标记实测Q峰，上图零轴居中。ζ=1从绘图和拟合中排除，原始记录保留。图1黑线约束经过红线的拟合零点，等宽矩形居中于该交点；共同零点是绘图约束，不作为独立零点一致性的证据。", "", "统一Hz参考相位和Hz RMS=1 A/m归一化；用户px/py对应内部py/px，分别观察主偏振Ex/Ey。边界曲线取负但标题依既定要求省略外部负号。", "", "| 模式 | 采样最高Q | 对应ζ | 身份待核验点数 |", "|---|---|---|---|"]
    for mode in ("px","py"):
        r=summary[mode];lines.append(f"| {mode} | {r['Q_maximum']:.8g} | {r['q_minimum_zeta']:g} | {len(r['mixed_zeta'])} |")
    lines += ["", "相同网格的完整扫描已计算；更细网格收敛与严格辐射零点仍需独立验证。源数据与详细判据见scientific_summary.json。", "", f"原始模型与日志：{ROOT}", ""]
    for i,name in enumerate(plotting.NAMES,1):lines.append(f"- Figure {i}：[PNG](../../10_overview/{name}.png) · [PDF](../../11_pdf/{GROUP}/{name}.pdf)")
    report=ARCHIVE/"12_reports"/GROUP/"README.md";report.parent.mkdir(parents=True,exist_ok=True);report.write_text("\n".join(lines)+"\n",encoding="utf-8")
    code=ARCHIVE/"99_config"/GROUP/"code";code.mkdir(exist_ok=True)
    for p in (Path(__file__),Path(dual.__file__),Path(plotting.__file__)):shutil.copy2(p,code/p.name)
    outputs=[ARCHIVE/"10_overview"/f"{name}.png" for name in plotting.NAMES]+[ARCHIVE/"11_pdf"/GROUP/f"{name}.pdf" for name in plotting.NAMES]
    completion=dict(status="computed_and_plotted_pending_visual_review",geometries=len(config["cases"]),mode_rows=len(df),edge_rows=len(edges),size_factor=SIZE_FACTOR,source_hashes_verified=True,visual_inspection_passed=False,output_sha256={str(p):s4.digest(p) for p in outputs})
    s4.write_json(ARCHIVE/"99_config"/GROUP/"completion.json",completion)
    return completion


def run(config,allow_concurrent):
    status_path=ROOT/"99_config/execution.json"
    record=dict(status="starting",pid=os.getpid(),cores=4,completed=[],preflight=s4.runtime_preflight(allow_concurrent=allow_concurrent))
    s4.write_json(status_path,record);client=None
    try:
        import mph
        from scripts.run_main import run_band_pair as band
        from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationConfig
        client=mph.start(version="6.3",cores=4)
        if not all(client.java.checkoutLicense(p) for p in ("COMSOL","WAVEOPTICS")):raise RuntimeError("Required license unavailable")
        record["license_checkout"]=True;s4.write_json(status_path,record)
        simconfig=SimulationConfig(slab_height="200 [nm]",slab_refractive_index="3.3",eigenmode_count=2,mesh_auto_size=5,eigenfrequency_shift=config["shift"])
        for case in sorted(config["cases"],key=lambda c:(0 if c["zeta"]==.835 else 1 if c["zeta"]==.9 else 2,c["zeta"])):
            folder=Path(case["directory"])
            if case.get("reuse_export"):
                saved=s4.json_read(folder/"export.json")
                if saved["status"]!="complete":raise ValueError("Incomplete source export")
                for name,digest in saved["output_hashes"].items():
                    if s4.digest(folder/name)!=digest:raise ValueError("Source export changed")
                record["completed"].append(case["zeta"])
                s4.write_json(status_path,record)
                print("REUSE_POINT",case["zeta"],flush=True)
                continue
            source=refined_model(client,case,record,status_path)
            model=client.load(str(source))
            try:
                stats=mesh_stats(model)
                if stats["automatic"] or stats["auto_size"]!=5:raise ValueError("Expected custom refined mesh")
                s4.write_json(folder/"99_config/mesh_verified.json",dict(size_factor=SIZE_FACTOR,mesh=stats,source_model=str(source),source_sha256=s4.digest(source)))
            finally:client.remove(model)
            extract_case=dict(case,source_model=str(source),source_model_sha256=s4.digest(source))
            record.update(status="extracting",current_zeta=case["zeta"]);s4.write_json(status_path,record)
            dual.extract_case(client,band,simconfig,config,extract_case)
            record["completed"].append(case["zeta"]);s4.write_json(status_path,record)
            print("COMPLETE_POINT",case["zeta"],len(record["completed"]),"/",len(config["cases"]),flush=True)
        record["status"]="postprocessing";s4.write_json(status_path,record)
        finish(config);record["status"]="computed_and_plotted_pending_visual_review"
    except Exception as exc:
        record.update(status="failed",error=str(exc),traceback=traceback.format_exc());raise
    finally:
        if client is not None:
            client.clear()
            if client.port:client.disconnect()
        record["source_unchanged"]=all(s4.digest(p)==h for p,h in config["source_hashes"].items())
        s4.write_json(status_path,record)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--prepare-only",action="store_true");action.add_argument("--run",action="store_true");action.add_argument("--plot-only",action="store_true")
    parser.add_argument("--allow-concurrent",action="store_true")
    args=parser.parse_args();config=prepare()
    print(f"81 geometries, {config['new_solves']} new solves, {config['reused_models']} reused models; 4 cores; h/h0={SIZE_FACTOR}; output={ROOT}",flush=True)
    if args.run:run(config,args.allow_concurrent)
    elif args.plot_only:finish(config)


if __name__=="__main__":main()
