import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import scripts.run_main.run_unit_cell_2d as unit_cell_2d
import scripts.run_sweep.run_unit_cell_2d_zeta_scan as scan


def _parameter_file(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    source = Path("scripts/parameter.json")
    data = json.loads(source.read_text(encoding="utf-8"))
    target = tmp_path / "parameter.json"
    target.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + '\n',
        encoding="utf-8",
    )
    return target, data


def _redirect_output_root(tmp_path: Path, monkeypatch) -> Path:
    output_root = tmp_path / "unit_cell_2D"
    monkeypatch.setattr(scan, "UNIT_CELL_2D_OUTPUT_ROOT", output_root)
    monkeypatch.setattr(
        unit_cell_2d,
        "UNIT_CELL_2D_OUTPUT_ROOT",
        output_root,
    )
    return output_root


def test_prepare_scan_changes_only_unit_cell_2d_zeta(tmp_path, monkeypatch):
    parameter_path, original = _parameter_file(tmp_path)
    output_root = _redirect_output_root(tmp_path, monkeypatch)

    prepared = scan.prepare_scan(parameter_path, (1.151, 1.161))

    assert prepared["scan_dir"].parent == output_root
    assert [case["zeta"] for case in prepared["cases"]] == pytest.approx(
        [1.151, 1.161]
    )
    assert json.loads(parameter_path.read_text(encoding="utf-8")) == original
    for case in prepared["cases"]:
        snapshot = json.loads(
            Path(case["parameter_path"]).read_text(encoding="utf-8")
        )
        expected = json.loads(json.dumps(original))
        expected["unit_cell_2d"]["zeta"] = case["zeta"]
        assert snapshot == expected
        assert snapshot["unit_cell_2d"]["isQuarter"] == 1
        assert f"-{case['zeta']:.3f})_" in Path(case["series_dir"]).name
        assert "_quarter_band-py" in Path(case["series_dir"]).name
    assert prepared["cases"][0]["series_dir"] != prepared["cases"][1][
        "series_dir"
    ]
    scan_config = json.loads(
        (prepared["scan_dir"] / "99_config" / scan.SCAN_CONFIG_FILENAME).read_text(
            encoding="utf-8"
        )
    )
    assert scan_config["isQuarter"] == 1


def test_run_scan_launches_cases_serially_in_requested_order(
    tmp_path, monkeypatch
):
    parameter_path, _original = _parameter_file(tmp_path)
    _redirect_output_root(tmp_path, monkeypatch)
    calls = []

    def fake_run(command, *, cwd, env, check):
        snapshot = json.loads(
            Path(env[scan.PARAMETER_PATH_ENV]).read_text(encoding="utf-8")
        )
        calls.append(
            {
                "command": command,
                "cwd": cwd,
                "check": check,
                "zeta": snapshot["unit_cell_2d"]["zeta"],
            }
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(scan.subprocess, "run", fake_run)
    summary = scan.run_scan(parameter_path, (1.151, 1.161))

    assert summary["status"] == "complete"
    assert [call["zeta"] for call in calls] == pytest.approx([1.151, 1.161])
    assert all(call["check"] is False for call in calls)
    assert all(
        call["command"][-1] == "scripts.run_main.run_unit_cell_2d"
        for call in calls
    )
    saved = json.loads(
        (
            Path(summary["scan_dir"])
            / "99_config"
            / scan.SCAN_SUMMARY_FILENAME
        ).read_text(encoding="utf-8")
    )
    assert saved["status"] == "complete"
    assert [case["status"] for case in saved["cases"]] == [
        "complete",
        "complete",
    ]


@pytest.mark.parametrize(
    "values",
    [(), (1.151, 1.151), (0.0,), (3.0**0.5,)],
)
def test_invalid_zeta_lists_are_rejected(values):
    with pytest.raises(ValueError):
        scan.normalize_zeta_values(values)


def test_pyproject_registers_zeta_scan_entry():
    source = Path("pyproject.toml").read_text(encoding="utf-8")
    assert (
        'comsol-unit-cell-2d-zeta-scan = '
        '"scripts.run_sweep.run_unit_cell_2d_zeta_scan:main"'
    ) in source
