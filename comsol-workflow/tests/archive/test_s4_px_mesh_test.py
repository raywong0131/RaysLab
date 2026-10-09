import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
import pytest
from scripts.run_main import run_s4_boundary as s4


def test_concurrent_comsol_requires_explicit_opt_in(monkeypatch):
    process=dict(ProcessId=123,Name='comsolmphserver.exe',CommandLine='existing-task')
    monkeypatch.setattr(subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=json.dumps([process])))
    monkeypatch.setattr(Path,'is_file',lambda p:True)
    monkeypatch.setattr(Path,'read_text',lambda *a,**k:'local test license')
    with pytest.raises(RuntimeError,match='do not share resources implicitly'):
        s4.runtime_preflight()
    assert s4.runtime_preflight(allow_concurrent=True)['concurrent_explicitly_allowed'] is True


def test_prepared_cases_cannot_be_reused_for_different_mesh_identity(tmp_path,monkeypatch):
    from scripts.run_sweep import run_s4_px_mesh_test as task
    monkeypatch.setattr(task,'OUTPUT',tmp_path)
    path=tmp_path/'99_config/config.json';path.parent.mkdir()
    path.write_text(json.dumps({'cases':[{'zeta':.836,'mesh':5,'size_factor':.5}],'source_sha256':{}}))
    with pytest.raises(ValueError,match='identity differs'):
        task.prepare([.25],[.836])
