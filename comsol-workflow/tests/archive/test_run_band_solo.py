from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_main" / "run_band_solo.py"


def load_solo():
    spec = importlib.util.spec_from_file_location("run_band_solo", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


solo = load_solo()


def test_solo_uses_exact_cavity_only_series_name():
    cavity = solo.ACTIVE_PARAMETERS.cavity
    expected = (
        f"unit_cell_({cavity.b0_nm:g}-{cavity.eta:g}-{cavity.zeta:g})_"
        f"mesh{solo.ACTIVE_PARAMETERS.mesh_size}"
    )

    assert solo.solo_series_name() == expected
    assert solo.OUT_DIR.parent == solo.UNIT_CELL_OUTPUT_ROOT
    assert solo.CAVITY_PARAMS == cavity.to_fourier_params("cavity_p_bic")
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    assert "ACTIVE_PARAMETERS.cladding" not in source
    assert "CLADDING_PARAMS" not in source


def test_run_band_solo_orchestrates_only_cavity(tmp_path, monkeypatch):
    calls = []
    cavity = pd.DataFrame({"cell": ["cavity"]})
    fits = pd.DataFrame({"status": ["ok"]})
    monkeypatch.setattr(
        solo.band,
        "case_k_points",
        lambda case_name: {"gamma_m": [case_name], "gamma_k": [case_name]},
    )
    monkeypatch.setattr(
        solo.band,
        "run_case",
        lambda case_name, params, output_root: (
            calls.append((case_name, params, Path(output_root))) or cavity
        ),
    )
    monkeypatch.setattr(solo.band, "_build_q_fits", lambda frame: fits)
    monkeypatch.setattr(
        solo.band,
        "plot_band_cases",
        lambda cases, output: calls.append(("band_plot", tuple(cases), Path(output))),
    )
    monkeypatch.setattr(
        solo.band,
        "plot_p_band_q_cases",
        lambda cases, output: calls.append(("q_plot", tuple(cases), Path(output))),
    )

    params = {"name": "cavity", "r_f0": 0.96, "b_square_f0": 1.0}
    result = solo.run_band_solo(tmp_path, params)

    assert calls[0] == ("cavity", params, tmp_path)
    assert [item[0] for item in calls[1:]] == ["band_plot", "q_plot"]
    assert result == {
        "output_dir": tmp_path,
        "cavity": cavity,
        "cavity_fits": fits,
    }
    assert (tmp_path / "10_overview" / "cavity" / "q_power_law_fits.csv").is_file()
    config = json.loads(
        (tmp_path / "99_config" / "config.json").read_text(encoding="utf-8")
    )
    assert config["workflow"] == "run_band_solo"
    assert set(config["case"]) == {"cavity"}
    assert set(config["sampled_points"]) == {"cavity"}
    assert set(config["parameters"]) == {
        "cavity",
        "mesh_size",
        "unit_cell_eigenmode_count",
    }
    assert not (tmp_path / "01_results" / "cladding").exists()
    assert not (tmp_path / "10_overview" / "cladding").exists()
    assert not (tmp_path / "99_config" / "cladding").exists()


def test_solo_center_update_ignores_cladding_identity_change(tmp_path):
    parameter_path = tmp_path / "parameter.json"
    source_data = json.loads(
        solo.ACTIVE_PARAMETERS.source_path.read_text(encoding="utf-8")
    )
    parameter_path.write_text(
        json.dumps(source_data, indent=2) + "\n",
        encoding="utf-8",
    )
    expected = solo.band.ACTIVE_PARAMETERS.__class__(
        **{
            **solo.band.ACTIVE_PARAMETERS.__dict__,
            "source_path": parameter_path,
        }
    )
    changed = json.loads(parameter_path.read_text(encoding="utf-8"))
    changed["cells"]["cladding"]["zeta"] = float(changed["cells"]["cladding"]["zeta"]) + 0.01
    parameter_path.write_text(
        json.dumps(changed, indent=2) + "\n",
        encoding="utf-8",
    )

    updated = solo.update_center_frequency_from_cavity_p2(
        198.4,
        expected_parameters=expected,
        path=parameter_path,
        workflow="run_band_solo",
        include_cladding=False,
    )
    saved = json.loads(parameter_path.read_text(encoding="utf-8"))

    assert updated.center_frequency_thz == 198.4
    assert saved["cells"]["cladding"]["zeta"] == changed["cells"]["cladding"]["zeta"]
    assert saved["structure"]["center_frequency_source"]["workflow"] == "run_band_solo"


def test_main_updates_center_frequency_as_solo_without_cladding(monkeypatch):
    frame = pd.DataFrame({"re": [198.4]})
    captured = {}
    monkeypatch.setattr(
        solo,
        "run_band_solo",
        lambda output_dir, params: {"cavity": frame},
    )
    monkeypatch.setattr(solo, "cavity_p2_gamma_frequency", lambda value: 198.4)

    class Updated:
        center_frequency_thz = 198.4
        source_path = Path("parameter.json")

    def fake_update(frequency, **kwargs):
        captured["frequency"] = frequency
        captured.update(kwargs)
        return Updated()

    monkeypatch.setattr(solo, "update_center_frequency_from_cavity_p2", fake_update)

    solo.main()

    assert captured["frequency"] == 198.4
    assert captured["expected_parameters"] is solo.ACTIVE_PARAMETERS
    assert captured["workflow"] == "run_band_solo"
    assert captured["include_cladding"] is False
