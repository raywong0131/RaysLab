import sys
import pytest
from scripts.run_sweep import run_s4_full_refined_scan as scan


def test_refined_scan_rejects_changed_identity_or_source(tmp_path):
    source=tmp_path/"source";source.write_text("immutable")
    config=dict(version=1,mesh=5,size_factor=.25,cores=4,eta=.96,b0_nm=245,a_nm=820,H_nm=200,n=3.3,edge_order=128,z_nm=100,Gamma=True,expected_eigenvalues=4,shift="c_const/1.55[um]",user_modes=["px","py"],internal_modes=["py","px"],cases=[dict(zeta=.9,snapshot=str(source),snapshot_sha256=scan.s4.digest(source))],source_hashes={str(source):scan.s4.digest(source)})
    scan.validate_config(config,[.9])
    for key,value in (("size_factor",.5),("cores",16),("edge_order",64)):
        with pytest.raises(ValueError,match="identity"):
            scan.validate_config(dict(config,**{key:value}),[.9])
    with pytest.raises(ValueError,match="identity"):scan.validate_config(config,[.8,.9])
    source.write_text("changed")
    with pytest.raises(ValueError,match="Protected source"):
        scan.validate_config(config,[.9])
    assert "mph" not in sys.modules


def test_completed_refined_model_is_reused_without_solver(tmp_path):
    model=tmp_path/"s4_gamma.mph";model.write_bytes(b"saved-solution")
    scan.s4.write_json(tmp_path/"solve.json",dict(size_factor=.25,model_sha256=scan.s4.digest(model)))
    case=dict(directory=str(tmp_path))
    assert scan.refined_model(None,case,{},tmp_path/"status.json")==model
    scan.s4.write_json(tmp_path/"solve.json",dict(size_factor=.5,model_sha256=scan.s4.digest(model)))
    with pytest.raises(ValueError,match="identity"):
        scan.refined_model(None,case,{},tmp_path/"status.json")


def test_local_plan_uses_measured_zero_bracket_and_deduplicates():
    import numpy as np
    import pandas as pd
    from scripts.run_sweep.run_s4_peak_resampling import sampling_plan
    rows=[]
    for mode,peak,zero in [("px",.835,.83513),("py",1.159,1.15905)]:
        for z in [peak-.005,peak,peak+.001]:
            rows.append(dict(user_mode=mode,zeta=z,Q=1/(1e-10+(z-zero)**2),air_display_V_per_m=zero-z,mode_status="single_eigenmode"))
    df=pd.DataFrame(rows);windows,points=sampling_plan(df)
    assert len(points)==36 and len(set(points))==36
    assert not set(points).intersection(df.zeta)
    for mode,zero in [("px",.83513),("py",1.15905)]:
        w=windows[mode]
        assert w["fine_range"][0]<zero<w["fine_range"][1]
        assert w["fine_step"]==.00002
        np.testing.assert_allclose(w["coordinate_only_zero_estimate"],zero)
    with pytest.raises(ValueError,match="sign bracket"):
        sampling_plan(df.assign(air_display_V_per_m=1.))
    assert "mph" not in sys.modules


def test_hole_volume_check_accepts_cad_roundoff_but_rejects_geometry_changes():
    import numpy as np
    expected=np.full(6,.0018117644302)
    actual=expected.copy();actual[0]+=2.08491565e-10
    result=scan.check_hole_volumes(actual,expected)
    assert result["relative_errors"][0]<result["rtol"]
    actual[0]=expected[0]*1.001
    with pytest.raises(AssertionError):scan.check_hole_volumes(actual,expected)
    with pytest.raises(ValueError):scan.check_hole_volumes(actual[:5],expected)
    with pytest.raises(ValueError):scan.check_hole_volumes(actual,np.zeros(6))
