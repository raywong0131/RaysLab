"""Path-only S4 checks; no COMSOL import, mesh, or eigensolve."""
import sys

import pytest

from scripts.run_main import run_s4_boundary as entry


def test_optional_session_root_is_scoped(tmp_path, monkeypatch):
    output_root = tmp_path / "unit_cell_2D"
    monkeypatch.setattr(entry, "UNIT_CELL_2D_OUTPUT_ROOT", output_root)
    monkeypatch.delenv("COMSOL_WORKFLOW_S4_SESSION_ROOT", raising=False)
    assert entry.configured_session_root() is None
    session = output_root / "S4_control"
    monkeypatch.setenv("COMSOL_WORKFLOW_S4_SESSION_ROOT", str(session))
    assert entry.configured_session_root() == session.resolve()
    for invalid in (output_root, tmp_path / "elsewhere", session / "nested"):
        monkeypatch.setenv("COMSOL_WORKFLOW_S4_SESSION_ROOT", str(invalid))
        with pytest.raises(ValueError, match="direct child"):
            entry.configured_session_root()


def test_session_prepare_and_identity_without_solver(tmp_path, monkeypatch):
    session = tmp_path / "unit_cell_2D/S4_control"
    monkeypatch.setattr(entry, "SESSION_ROOT", session)
    monkeypatch.setattr(entry, "TASK_ROOT", session / "tmp")
    source = entry.parameter_path_from_environment()
    source_sha = entry.digest(source)
    entry.main(["prepare", "--run-id", "session-check", "--parameter-file", str(source),
                "--eta", "0.95", "--zeta", "1", "--mesh", "5"])
    preparation = session / "tmp/s4_session-check"
    manifest = entry.prepared_manifest(preparation)
    output = entry.Path(manifest["output_directory"])
    assert output.parent == session / "results"
    assert not output.exists()
    assert manifest["cases"][0]["eta"] == 0.95
    assert manifest["cases"][0]["zeta"] == 1
    assert manifest["cases"][0]["mesh"] == 5
    assert manifest["execution"]["solve_count"] == 1
    assert entry.digest(source) == source_sha
    assert "mph" not in sys.modules
    with pytest.raises(FileExistsError):
        entry.main(["prepare", "--run-id", "session-check"])
    # The new path option does not waive manifest integrity.
    identity = preparation / "99_config/s4_manifest.json"
    manifest["output_directory"] = str(tmp_path / "invalid")
    entry.write_json(identity, manifest)
    with pytest.raises(RuntimeError, match="Manifest changed"):
        entry.prepared_manifest(preparation)
