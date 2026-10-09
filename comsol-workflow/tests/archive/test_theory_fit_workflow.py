import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.fit_theory_band_parameters import (
    BAND_COMPARISON_FILENAME,
    FULL_BAND_COMPARISON_FILENAME,
    FULL_BAND_SECOND_ORDER_FILENAME,
    BAND_MODE_COLORS,
    BAND_PANEL_TITLES,
    gamma_runtime_mapping,
    labeled_ray,
)
from scripts.run_sweep.collect_theory_fit_bands import build_manifest, sampling_branches


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_band_fit_plot_contract_is_locked():
    assert BAND_COMPARISON_FILENAME == "band_fit_comparison.png"
    assert FULL_BAND_COMPARISON_FILENAME == "full_band_fit_comparison.png"
    assert FULL_BAND_SECOND_ORDER_FILENAME == "full_band_second_order_comparison.png"
    assert BAND_PANEL_TITLES == (r"$\Gamma$–K ($+k_x$)", r"$\Gamma$–M ($+k_y$)")
    assert BAND_MODE_COLORS == {
        "d_xy": "#b8b8b8",
        "d_x2_minus_y2": "#111111",
        "p_y": "#d62728",
        "p_x": "#1f77b4",
    }


def test_manifest_locks_eta_scan_sampling_and_solver_snapshot(tmp_path):
    fit_config = json.loads(
        (REPO_ROOT / "scripts" / "theory_band_fit.json").read_text(encoding="utf-8")
    )
    parameters = json.loads(
        (REPO_ROOT / "scripts" / "parameter.json").read_text(encoding="utf-8")
    )
    manifest = build_manifest(fit_config, parameters, tmp_path / "dataset")
    calibration = [item for item in manifest["structures"] if item["role"] == "calibration"]
    assert [item["cell"]["eta"] for item in calibration] == [
        0.94, 0.96, 0.98, 1.0, 1.02, 1.04
    ]
    assert manifest["status"] == "dry_run"
    assert manifest["will_start_comsol"] is False
    assert manifest["simulation"]["mesh_auto_size"] == parameters["mesh_size"]
    collector_source = (
        REPO_ROOT / "scripts" / "run_sweep" / "collect_theory_fit_bands.py"
    ).read_text(encoding="utf-8")
    pre_execute_source = collector_source.split("def execute_manifest", 1)[0]
    assert "from scripts.run_main import run_band_pair" not in pre_execute_source
    assert "import mph" not in pre_execute_source
    target = next(item for item in manifest["structures"] if item["role"] == "target")
    assert set(target["branches"]) == {
        "fit_x", "fit_y", "validation_1", "validation_2", "qa_1", "qa_2"
    }


def test_sampling_axes_use_approved_near_gamma_points():
    config = json.loads(
        (REPO_ROOT / "scripts" / "theory_band_fit.json").read_text(encoding="utf-8")
    )
    branches = sampling_branches(config, include_validation=False)
    assert [point["kx"] for point in branches["fit_x"]] == [
        0.0, 0.002, 0.004, 0.006, 0.01, 0.015, 0.02
    ]
    assert all(point["ky"] == 0.0 for point in branches["fit_x"])


def test_selected_band_reader_maps_internal_modes_to_s17_order():
    runtime = {"p1": "py", "p2": "px", "d1": "dy", "d2": "dx"}
    rows = []
    gamma_frequency = {"p1": 1.0, "p2": 2.0, "d1": 3.0, "d2": 4.0}
    for point_index, kx in enumerate([0.0, 0.01]):
        for label, dominant in runtime.items():
            rows.append({
                "direction": "fit_x", "point_index": point_index,
                "band_label": label, "kx": kx, "ky": 0.0,
                "re": gamma_frequency[label] + point_index * 0.1,
                "is_valid": True, "match_status": "gamma" if point_index == 0 else "matched",
                "dominant_mode": dominant,
            })
    frame = pd.DataFrame(rows)
    assert gamma_runtime_mapping(frame) == {
        "p1": "p_x", "p2": "p_y", "d1": "d_xy", "d2": "d_x2_minus_y2"
    }
    points, frequencies = labeled_ray(frame, "fit_x")
    assert np.allclose(points, [[0.0, 0.0], [0.01, 0.0]])
    assert np.allclose(frequencies[0], [3.0, 4.0, 2.0, 1.0])
