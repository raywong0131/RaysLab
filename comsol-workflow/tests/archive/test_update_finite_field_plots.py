from pathlib import Path

import pandas as pd
import pytest

from scripts.analysis import update_finite_field_plots as updater


def _make_mode_export(root: Path) -> Path:
    case_dir = root / "series" / "shift0.000"
    model_path = case_dir / "00_model" / "finite_quarter.mph"
    model_path.parent.mkdir(parents=True)
    model_path.write_bytes(b"saved model")
    mode_dir = case_dir / "01_results" / "mode0"
    export_dir = mode_dir / "11_simulation_exports"
    export_dir.mkdir(parents=True)
    pd.DataFrame({"x": [0.0], "y": [0.0], "re": [1.0], "im": [0.0]}).to_parquet(
        export_dir / "Hz_center.parquet",
        index=False,
    )
    for filename in ("Hz_Re_2d.png", "Hz_Im_2d.png", "Wem_2d.png"):
        (export_dir / filename).write_bytes(b"old")
    overview_dir = mode_dir / "10_overview"
    overview_dir.mkdir()
    pd.DataFrame(
        [
            {
                "re": 195.25,
                "im": 0.1,
                "q": 976.25,
                "mode_idx": 0,
                "symmetry_id": 1,
                "x_boundary": "PEC",
                "y_boundary": "PMC",
            }
        ]
    ).to_csv(overview_dir / "eigenfrequency.csv", index=False)
    return export_dir


def test_discover_mode_exports_finds_saved_field_inputs(tmp_path):
    export_dir = _make_mode_export(tmp_path)

    modes = updater.discover_mode_exports(tmp_path)

    assert len(modes) == 1
    assert modes[0].export_dir == export_dir.resolve()
    assert modes[0].mode_idx == 0
    assert modes[0].frequency_thz == 195.25
    assert modes[0].quality_factor == 976.25
    assert modes[0].is_quarter is True
    assert modes[0].model_path.name == "finite_quarter.mph"


def test_transactional_replace_rolls_back_all_images_on_failure(
    monkeypatch,
    tmp_path,
):
    staging_dir = tmp_path / "stage"
    first_staged = staging_dir / "new" / "first.png"
    second_staged = staging_dir / "new" / "second.png"
    first_staged.parent.mkdir(parents=True)
    first_staged.write_bytes(b"new first")
    second_staged.write_bytes(b"new second")
    first_destination = tmp_path / "first.png"
    second_destination = tmp_path / "second.png"
    first_destination.write_bytes(b"old first")
    second_destination.write_bytes(b"old second")
    real_replace = updater.os.replace
    calls = 0

    def fail_second_replace(source, destination):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected replacement failure")
        return real_replace(source, destination)

    monkeypatch.setattr(updater.os, "replace", fail_second_replace)

    with pytest.raises(OSError, match="injected replacement failure"):
        updater.replace_images_transactionally(
            [
                (first_staged, first_destination),
                (second_staged, second_destination),
            ],
            staging_dir,
        )

    assert first_destination.read_bytes() == b"old first"
    assert second_destination.read_bytes() == b"old second"
    assert (staging_dir / "backup").is_dir()
