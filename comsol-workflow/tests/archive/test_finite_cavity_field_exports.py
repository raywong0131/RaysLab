import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts.run_main import run_finite, run_finite_quarter, run_strip_1d


class FakeFiniteSimulation:
    def __init__(self):
        self.image_calls = []
        self.three_d_calls = []

    def get_eigenfrequencies(self):
        return [[194.9, 0.1, 1000.0]]

    def get_2d_fields(self, mode_idx, expr, plane):
        points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        return points, np.array([1.0, 2.0, 3.0])

    def export_2d_fields(self, *args, **kwargs):
        self.image_calls.append((args, kwargs))

    def export_3d_fields(self, *args, **kwargs):
        self.three_d_calls.append((args, kwargs))


def test_full_finite_export_uses_requested_images_and_keeps_hz_data(monkeypatch, tmp_path):
    sim = FakeFiniteSimulation()
    plot_calls = []
    monkeypatch.setattr(
        run_finite,
        "save_field_plot",
        lambda path, _coordinates, _values, **kwargs: plot_calls.append(
            (Path(path).name, kwargs)
        ),
    )
    monkeypatch.setattr(
        run_finite,
        "mode_validity_for_boundary",
        lambda *_args, **_kwargs: [{"is_valid": True}],
    )
    monkeypatch.setattr(
        "comsol_workflow.energy_recovery.interpolate_field",
        lambda _source_points, values, _query_points: values,
    )

    run_finite.export_full_finite_results(
        sim,
        tmp_path,
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        label="finite_cavity",
    )

    assert [
        (
            name,
            kwargs["frequency_thz"],
            kwargs["quality_factor"],
            kwargs["quantity_label"],
            kwargs["unit_label"],
            kwargs["color_scale_mode"],
        )
        for name, kwargs in plot_calls
    ] == [
        ("00_Hz_Re_2d.png", 194.9, 1000.0, "Re(Hz)", "A/m", "linearsymmetric"),
        ("00_Hz_Im_2d.png", 194.9, 1000.0, "Im(Hz)", "A/m", "linearsymmetric"),
        ("00_Wem_2d.png", 194.9, 1000.0, "Wem", "J/m³", "linear"),
    ]
    assert sim.image_calls == []
    assert sim.three_d_calls == []
    field = pd.read_parquet(tmp_path / "00_Hz_center.parquet")
    assert list(field.columns) == ["x", "y", "re", "im"]


def test_quarter_export_uses_same_annotated_field_plot_contract(monkeypatch, tmp_path):
    sim = FakeFiniteSimulation()
    plot_calls = []
    triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    monkeypatch.setattr(
        run_finite,
        "mode_validity_for_boundary",
        lambda *_args, **_kwargs: [{"is_valid": True}],
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "_reconstruct_hz_field",
        lambda *_args, **_kwargs: (triangle, np.array([1 + 2j, 2 + 3j, 3 + 4j])),
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "_reconstruct_wem_field",
        lambda *_args, **_kwargs: (triangle, np.array([1.0, 2.0, 3.0])),
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "save_field_plot",
        lambda path, _coordinates, _values, **kwargs: plot_calls.append(
            (Path(path).name, kwargs)
        ),
    )
    monkeypatch.setattr(run_finite, "RUN_FARFIELD_FFT", False)

    run_finite_quarter.export_reconstructed_results(
        sim,
        tmp_path,
        triangle,
        finite_hex_a=2.0,
        case={"id": 1, "x_boundary": "PEC", "y_boundary": "PMC"},
    )

    assert [
        (name, kwargs["quantity_label"], kwargs["unit_label"])
        for name, kwargs in plot_calls
    ] == [
        ("00_Hz_Re_2d.png", "Re(Hz)", "A/m"),
        ("00_Hz_Im_2d.png", "Im(Hz)", "A/m"),
        ("00_Wem_2d.png", "Wem", "J/m³"),
    ]
    assert all(kwargs["frequency_thz"] == 194.9 for _name, kwargs in plot_calls)
    assert all(kwargs["quality_factor"] == 1000.0 for _name, kwargs in plot_calls)


def test_quarter_export_allows_no_valid_modes(monkeypatch, tmp_path):
    sim = FakeFiniteSimulation()
    triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    monkeypatch.setattr(
        run_finite,
        "mode_validity_for_boundary",
        lambda *_args, **_kwargs: [{"is_valid": False}],
    )
    monkeypatch.setattr(run_finite, "RUN_FARFIELD_FFT", False)

    frame = run_finite_quarter.export_reconstructed_results(
        sim,
        tmp_path,
        triangle,
        finite_hex_a=2.0,
        case={"id": 4, "x_boundary": "PMC", "y_boundary": "PMC"},
    )

    assert frame["is_valid"].tolist() == [False]
    assert (tmp_path / "eigenfrequencies.csv").is_file()
    assert not list(tmp_path.glob("*_Hz_center.parquet"))


def test_quarter_candidate_export_is_lightweight_and_full_export_is_selected_only(
    monkeypatch,
    tmp_path,
):
    sim = FakeFiniteSimulation()
    sim.get_eigenfrequencies = lambda: [
        [194.9, 0.1, 1000.0],
        [195.1, 0.1, 900.0],
    ]
    triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    hz_calls = []
    wem_calls = []
    plot_calls = []
    monkeypatch.setattr(
        run_finite,
        "mode_validity_for_boundary",
        lambda *_args, **_kwargs: [
            {"is_valid": True},
            {"is_valid": True},
        ],
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "_reconstruct_hz_field",
        lambda _sim, mode_idx, _case: (
            hz_calls.append(mode_idx)
            or (
                triangle,
                np.full(3, mode_idx + 1j * (mode_idx + 1)),
            )
        ),
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "_reconstruct_wem_field",
        lambda _sim, mode_idx: (
            wem_calls.append(mode_idx)
            or (triangle, np.full(3, mode_idx + 1.0))
        ),
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "save_field_plot",
        lambda path, *_args, **_kwargs: plot_calls.append(Path(path).name),
    )
    monkeypatch.setattr(run_finite, "RUN_FARFIELD_FFT", False)

    eigen_df = run_finite_quarter.export_reconstructed_candidate_hz_results(
        sim,
        tmp_path,
        triangle,
        case={"id": 1, "x_boundary": "PEC", "y_boundary": "PMC"},
    )

    assert hz_calls == [0, 1]
    assert plot_calls == []
    assert wem_calls == []
    assert (tmp_path / "00_Hz_center.parquet").is_file()
    assert (tmp_path / "01_Hz_center.parquet").is_file()

    run_finite_quarter.export_reconstructed_selected_results(
        sim,
        tmp_path,
        eigen_df,
        [1],
        finite_hex_a=2.0,
        case={"id": 1, "x_boundary": "PEC", "y_boundary": "PMC"},
    )

    assert plot_calls == [
        "01_Hz_Re_2d.png",
        "01_Hz_Im_2d.png",
        "01_Wem_2d.png",
    ]
    assert wem_calls == [1]


def test_strip_export_matches_finite_field_image_contract(monkeypatch, tmp_path):
    sim = FakeFiniteSimulation()
    plot_calls = []
    monkeypatch.setattr(
        run_strip_1d,
        "save_field_plot",
        lambda path, _coordinates, _values, **kwargs: plot_calls.append(
            (Path(path).name, kwargs)
        ),
    )
    monkeypatch.setattr(
        run_strip_1d,
        "mode_validity",
        lambda *_args, **_kwargs: [{"is_valid": True}],
    )
    monkeypatch.setattr(
        "comsol_workflow.energy_recovery.interpolate_field",
        lambda _source_points, values, query_points: np.resize(
            np.asarray(values),
            len(query_points),
        ),
    )

    run_strip_1d.export_simulation_results(
        sim,
        tmp_path,
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        label="strip_1d",
    )

    image_contract = [
        (
            name,
            kwargs["frequency_thz"],
            kwargs["quality_factor"],
            kwargs["quantity_label"],
            kwargs["unit_label"],
            kwargs["color_scale_mode"],
        )
        for name, kwargs in plot_calls
    ]
    assert image_contract == [
        ("00_Hz_Re_2d.png", 194.9, 1000.0, "Re(Hz)", "A/m", "linearsymmetric"),
        ("00_Hz_Im_2d.png", 194.9, 1000.0, "Im(Hz)", "A/m", "linearsymmetric"),
        ("00_Wem_2d.png", 194.9, 1000.0, "Wem", "J/m³", "linear"),
    ]
    assert sim.image_calls == []
    assert sim.three_d_calls == []
    field = pd.read_parquet(tmp_path / "00_Hz_center.parquet")
    assert list(field.columns) == ["x", "y", "re", "im"]


def test_finite_quarter_defaults_to_user_py_symmetry_and_parameterized_name():
    assert run_finite_quarter.OUTPUT_ROOT == run_finite.OUTPUT_ROOT
    assert (
        run_finite.EIGENMODE_COUNT
        == run_finite.ACTIVE_PARAMETERS.finite_eigenmode_count
    )
    assert (
        run_finite_quarter.quarter_simulation_config(
            run_finite_quarter.symmetry_case(1)
        ).eigenmode_count
        == run_finite.EIGENMODE_COUNT
    )
    assert run_finite_quarter.validate_symmetry_ids(
        list(run_finite_quarter.SYMMETRY_IDS)
    ) == list(run_finite_quarter.SYMMETRY_IDS)
    assert run_finite_quarter.SYMMETRY_CASES[1] == {
        "x_boundary": "PEC",
        "y_boundary": "PMC",
    }
    assert run_finite_quarter.symmetry_dirname(1) == "symmetry_1_xPEC_yPMC"
    assert run_finite_quarter.OUT_DIR == (
        run_finite_quarter.OUTPUT_ROOT
        / run_finite.ACTIVE_PARAMETERS.structure_series_label("finite_quarter")
    )
    if run_finite.ACTIVE_PARAMETERS.cladding_shift_geometry.kind == "ellipse":
        assert run_finite_quarter.OUT_DIR.name.endswith("_ellipse-shift")
    else:
        assert not run_finite_quarter.OUT_DIR.name.endswith("_ellipse-shift")
    assert run_finite_quarter.quarter_case_stem(0.03) == "shift0.030"
    assert run_finite_quarter.quarter_case_stem(0.05, 0.12) == (
        "shiftx0.050_shifty0.120"
    )
    assert run_finite_quarter.quarter_output_case_dir(0.03, 1, [1]) == (
        run_finite_quarter.OUT_DIR / "shift0.030"
    )
    assert run_finite_quarter.quarter_output_case_dir(0.03, 1, [1, 2]) == (
        run_finite_quarter.OUT_DIR
        / "shift0.030"
        / "symmetry_1_xPEC_yPMC"
    )
    assert run_finite_quarter.quarter_output_case_dir(
        0.05,
        1,
        [1],
        y_shift_factor=0.12,
    ) == (
        run_finite_quarter.OUT_DIR / "shiftx0.050_shifty0.120"
    )


def test_finite_main_scans_every_ratio_derived_pair(monkeypatch, tmp_path):
    pairs = ((0.0, 0.0), (0.05, 0.08), (0.09, 0.144))
    calls = []
    monkeypatch.setattr(run_finite, "OUT_DIR", tmp_path / "finite")
    monkeypatch.setattr(
        run_finite,
        "ACTIVE_PARAMETERS",
        SimpleNamespace(finite_cladding_shift_factor_pairs=pairs),
    )
    monkeypatch.setattr(
        run_finite,
        "run_finite_cavity_shift",
        lambda x_factor, y_factor: (
            calls.append((x_factor, y_factor))
            or {"out_dir": f"case-{x_factor:.3f}-{y_factor:.3f}"}
        ),
    )

    run_finite.main()

    assert calls == list(pairs)
    summary = json.loads((run_finite.OUT_DIR / "run_summary.json").read_text())
    assert len(summary["runs"]) == len(pairs)


def test_finite_quarter_main_scans_every_resolved_pair(monkeypatch, tmp_path):
    pairs = ((0.0, 0.08), (0.05, 0.0))
    calls = []
    monkeypatch.setattr(run_finite_quarter, "OUT_DIR", tmp_path / "quarter")
    monkeypatch.setattr(
        run_finite,
        "ACTIVE_PARAMETERS",
        SimpleNamespace(finite_cladding_shift_factor_pairs=pairs),
    )
    monkeypatch.setattr(
        run_finite_quarter,
        "run_finite_quarter_shift",
        lambda x_factor, y_factor: (
            calls.append((x_factor, y_factor))
            or [{"out_dir": f"case-{x_factor:.3f}-{y_factor:.3f}"}]
        ),
    )

    run_finite_quarter.main()

    assert calls == list(pairs)
    summary = json.loads(
        (run_finite_quarter.OUT_DIR / "run_summary.json").read_text()
    )
    assert len(summary["runs"]) == len(pairs)


def test_finite_model_geometry_figure_set_is_compact_and_model_local(
    monkeypatch,
    tmp_path,
):
    calls = []
    monkeypatch.setattr(
        run_finite,
        "save_finite_simulation_plot",
        lambda path, *_args, **kwargs: calls.append(
            (Path(path).name, kwargs)
        ),
    )
    monkeypatch.setattr(
        run_finite,
        "save_full_lattice_shift_plot",
        lambda path, *_args, **_kwargs: calls.append((Path(path).name, {})),
    )
    monkeypatch.setattr(
        run_finite,
        "save_full_lattice_normal_fan_plot",
        lambda path, *_args, **_kwargs: calls.append((Path(path).name, {})),
    )

    run_finite.save_finite_model_geometry_figures(
        tmp_path / "00_model",
        [],
        [],
        np.empty((0, 2)),
        np.empty((0, 2)),
        {},
    )

    assert [name for name, _kwargs in calls] == [
        "full_finite_simulation.png",
        "full_lattice_modulation_vectors.png",
        "full_lattice_normal_fan_regions.png",
    ]
    assert calls[0][1] == {"annotated": False}
    assert run_finite.SHIFT_ARROW_PLOT_SCALE == 10.0


def test_finite_model_geometry_figure_uses_ellipse_diagnostic_for_ellipse_mode(
    monkeypatch,
    tmp_path,
):
    calls = []
    monkeypatch.setattr(
        run_finite,
        "save_finite_simulation_plot",
        lambda path, *_args, **_kwargs: calls.append(Path(path).name),
    )
    monkeypatch.setattr(
        run_finite,
        "save_full_lattice_shift_plot",
        lambda path, *_args, **_kwargs: calls.append(Path(path).name),
    )
    monkeypatch.setattr(
        run_finite,
        "save_full_lattice_ellipse_shift_plot",
        lambda path, *_args, **_kwargs: calls.append(Path(path).name),
    )
    monkeypatch.setattr(
        run_finite,
        "save_full_lattice_normal_fan_plot",
        lambda *_args, **_kwargs: pytest.fail("groups diagnostic was selected"),
    )

    run_finite.save_finite_model_geometry_figures(
        tmp_path / "00_model",
        [],
        [],
        np.empty((0, 2)),
        np.empty((0, 2)),
        {"cladding_shift_geometry": {"kind": "ellipse"}},
    )

    assert calls == [
        "full_finite_simulation.png",
        "full_lattice_modulation_vectors.png",
        "full_lattice_ellipse_shift.png",
    ]


def test_finite_output_layout_is_flat_mode_centric_and_cleans_staging(tmp_path):
    case_dir = tmp_path / "shift0.050"
    paths = run_finite.prepare_finite_case_output(
        case_dir,
        model_filename="finite_quarter.mph",
        resume=False,
    )
    export_dir = paths["export_dir"]
    export_dir.mkdir(parents=True)
    eigen_df = pd.DataFrame(
        [
            {"re": 196.1, "im": -0.01, "q": 100.0, "is_valid": True, "mode_idx": 0},
            {"re": 198.4, "im": -0.001, "q": 1000.0, "is_valid": True, "mode_idx": 11},
        ]
    )
    eigen_df.to_csv(export_dir / "eigenfrequencies.csv", index=False)
    (export_dir / "air_field_metadata.json").write_text("{}", encoding="utf-8")
    for mode_idx in (0, 11):
        for suffix in (
            "Hz_center.parquet",
            "Hz_Re_2d.png",
            "Hz_Im_2d.png",
            "Wem_2d.png",
            "E_air.parquet",
        ):
            (export_dir / f"{mode_idx:02d}_{suffix}").write_bytes(b"result")

    farfield_dir = paths["staging_dir"] / "farfield_fft"
    farfield_dir.mkdir()
    (farfield_dir / "config.json").write_text("{}", encoding="utf-8")
    farfield_rows = []
    for mode_idx in (0, 11):
        mode_dir = farfield_dir / f"mode_{mode_idx:02d}"
        mode_dir.mkdir()
        plots = []
        for name in (
            "A_overview.png",
            "B_polarization.png",
            "C_cutlines.png",
            "D_gauss_fit.png",
        ):
            path = mode_dir / name
            path.write_bytes(b"plot")
            plots.append(str(path))
        payload = {
            "status": "ok",
            "mode_idx": mode_idx,
            "source_file": str(export_dir / f"{mode_idx:02d}_E_air.parquet"),
            "plots": plots,
        }
        (mode_dir / "summary.json").write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        farfield_rows.append(payload)
    pd.DataFrame(farfield_rows).to_csv(
        farfield_dir / "farfield_summary.csv",
        index=False,
    )

    lattice_dir = paths["staging_dir"] / "finite_lattice_fourier_hz"
    lattice_dir.mkdir()
    (lattice_dir / "index_r_cavity.png").write_bytes(b"bulk")
    (lattice_dir / "index_k_1stBZ.png").write_bytes(b"grid")
    lattice_mode = lattice_dir / "mode_11"
    lattice_mode.mkdir()
    for name in (
        "finite_lattice_fourier_hz.npz",
        "intensity_maps_normH.png",
        "k_weight_Hz.png",
        "k_weight_tops_norm.png",
        "k_weight_tops_phase.png",
        "p_subspace_mode_decomposition.csv",
        "top_k_peaks.csv",
    ):
        (lattice_mode / name).write_bytes(b"lattice")
    pd.DataFrame(
        [{"mode_idx": 11, "output_dir": str(lattice_mode)}]
    ).to_csv(
        lattice_dir / "finite_lattice_fourier_summary.csv",
        index=False,
    )

    score_df = pd.DataFrame(
        [
            {"mode_idx": 0, "mode_score": 0.1},
            {"mode_idx": 11, "mode_score": 0.9},
        ]
    )
    score_df.to_csv(paths["staging_dir"] / "mode_scores.csv", index=False)
    pd.DataFrame(
        [
            {"mode_idx": 0, "analysis_selected": False},
            {"mode_idx": 11, "analysis_selected": True},
        ]
    ).to_csv(
        paths["staging_dir"] / "analysis_mode_selection.csv",
        index=False,
    )
    (paths["staging_dir"] / "objective.json").write_text(
        "{}",
        encoding="utf-8",
    )

    run_finite.finalize_finite_case_output(
        case_dir,
        model_filename="finite_quarter.mph",
    )

    assert not paths["staging_dir"].exists()
    assert (
        paths["overview_dir"] / "analysis_mode_selection.csv"
    ).is_file()
    assert (paths["results_dir"] / "mode0").is_dir()
    assert (paths["results_dir"] / "mode11").is_dir()
    assert not (paths["results_dir"] / "mode00").exists()
    assert (
        paths["results_dir"]
        / "mode11"
        / "11_simulation_exports"
        / "Hz_center.parquet"
    ).is_file()
    assert (
        paths["results_dir"]
        / "mode0"
        / "12_farfield_FFT"
        / "farfield_summary.json"
    ).is_file()
    assert not (
        paths["results_dir"] / "mode0" / "13_lattice_fourier_Hz"
    ).exists()
    assert (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "finite_lattice_fourier_hz.npz"
    ).is_file()
    assert (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "intensity_maps_normH.png"
    ).is_file()
    assert (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "k_weight_Hz.png"
    ).is_file()
    assert not (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "k_weight_first_bz.png"
    ).exists()
    assert (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "k_weight_tops_norm.png"
    ).is_file()
    assert (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "k_weight_tops_phase.png"
    ).is_file()
    assert not (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "top_p_abs.png"
    ).exists()
    assert not (
        paths["results_dir"]
        / "mode11"
        / "13_lattice_fourier_Hz"
        / "top_p_phase.png"
    ).exists()
    assert (
        paths["overview_dir"] / "index_r_cavity.png"
    ).is_file()
    assert (
        paths["overview_dir"] / "index_k_1stBZ.png"
    ).is_file()
    assert not (paths["overview_dir"] / "lattice_bulk_cyclic_index.png").exists()
    assert not (paths["overview_dir"] / "lattice_finite_k_grid_first_bz.png").exists()
    assert not list(
        paths["results_dir"].rglob("index_r_cavity.png")
    )
    assert not list(
        paths["results_dir"].rglob("index_k_1stBZ.png")
    )
    assert (paths["config_dir"] / "objective.json").is_file()

    mode_root_files = [
        path
        for path in paths["results_dir"].glob("mode*")
        for path in path.iterdir()
        if path.is_file()
    ]
    assert mode_root_files == []
    nested_categories = [
        path
        for path in paths["results_dir"].rglob("*")
        if path.is_dir() and path.name in {"data", "metrics", "figures"}
    ]
    assert nested_categories == []

    farfield_summary = json.loads(
        (
            paths["results_dir"]
            / "mode11"
            / "12_farfield_FFT"
            / "farfield_summary.json"
        ).read_text(encoding="utf-8")
    )
    assert Path(farfield_summary["source_file"]).is_file()
    assert all(Path(path).is_file() for path in farfield_summary["plots"])
    lattice_summary = pd.read_csv(
        paths["overview_dir"] / "finite_lattice_fourier_summary.csv"
    )
    assert Path(lattice_summary.loc[0, "output_dir"]) == (
        paths["results_dir"] / "mode11" / "13_lattice_fourier_Hz"
    )


def test_quarter_geometry_only_uses_formal_default_case_layout(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setattr(run_finite_quarter, "OUT_DIR", tmp_path)
    monkeypatch.setattr(run_finite_quarter, "RUN_COMSOL", False)
    monkeypatch.setattr(run_finite_quarter, "RESUME_FROM_EXISTING_MPH", False)
    monkeypatch.setattr(run_finite_quarter, "SYMMETRY_IDS", [1])
    monkeypatch.setattr(
        run_finite,
        "save_finite_shift_profile_diagnostics",
        lambda *_args, **_kwargs: None,
    )
    geometry_calls = []
    monkeypatch.setattr(
        run_finite,
        "save_finite_model_geometry_figures",
        lambda model_dir, *_args, **_kwargs: geometry_calls.append(
            Path(model_dir)
        ),
    )
    monkeypatch.setattr(run_finite, "full_cell_records", lambda _metadata: [])

    summary = run_finite_quarter.run_symmetry_case(
        0.05,
        1,
        [],
        [],
        {
            "quarter_boundary": np.array(
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
            ),
            "finite_boundary": np.array(
                [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
            ),
        },
        2.0,
        0.12,
    )

    case_dir = tmp_path / "shiftx0.050_shifty0.120"
    assert summary["out_dir"] == str(case_dir)
    assert summary["output_layout"] == run_finite.OUTPUT_LAYOUT_VERSION
    assert (case_dir / "00_model").is_dir()
    assert (case_dir / "01_results").is_dir()
    assert (case_dir / "10_overview").is_dir()
    assert (case_dir / "80_logs").is_dir()
    assert (case_dir / "99_config" / "config.json").is_file()
    assert not (case_dir / run_finite.STAGING_DIRNAME).exists()
    assert not (case_dir / "symmetry_1_xPEC_yPMC").exists()
    assert geometry_calls == [case_dir / "00_model"]


def test_finite_output_layout_rejects_unknown_staging_before_moving(tmp_path):
    case_dir = tmp_path / "shift0.070"
    paths = run_finite.prepare_finite_case_output(
        case_dir,
        model_filename="finite_cavity.mph",
        resume=False,
    )
    paths["export_dir"].mkdir(parents=True)
    pd.DataFrame(
        [{"re": 198.4, "im": -0.001, "q": 1000.0, "mode_idx": 0}]
    ).to_csv(paths["export_dir"] / "eigenfrequencies.csv", index=False)
    unknown = paths["staging_dir"] / "unexpected.bin"
    unknown.write_bytes(b"unknown")

    with pytest.raises(ValueError, match="Unknown staged finite output"):
        run_finite.finalize_finite_case_output(
            case_dir,
            model_filename="finite_cavity.mph",
        )

    assert (paths["export_dir"] / "eigenfrequencies.csv").is_file()
    assert unknown.is_file()
    assert not (paths["overview_dir"] / "eigenfrequencies.csv").exists()
