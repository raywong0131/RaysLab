"""Focused validity checks for strip and finite Fourier scoring."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import comsol_workflow.lattice_fourier_postprocess as lattice_postprocess
from comsol_workflow.lattice_fourier_postprocess import (
    finite_mode_envelope_metrics,
    score_finite_modes,
    score_strip_modes,
    select_finite_fundamental_mode,
    select_finite_mode_closest_to_frequency,
    selected_mode_indices,
)


def _write_decomposition(
    path: Path,
    *,
    component_name: str,
    gamma_subspace_p: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                component_name: 0,
                "k_weight_fraction": 0.8,
                "subspace_p": gamma_subspace_p,
                "mode_px": 0.95,
                "dominant_subspace": "p",
                "dominant_mode": "px",
                "subspace_matches_mode": True,
            }
        ]
    ).to_csv(path, index=False)


def test_finite_single_mode_selects_closest_frequency_and_records_audit(tmp_path):
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame(
        {
            "mode_idx": [0, 1],
            "re": [197.0, 198.2],
            "q": [800.0, 900.0],
            "is_valid": [True, True],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    output = tmp_path / "analysis_mode_selection.csv"

    selected, selection = select_finite_mode_closest_to_frequency(
        export_dir,
        output,
        target_frequency_thz=198.0,
    )

    assert selected == 1
    assert selection.loc[
        selection["mode_idx"].eq(1), "analysis_selected"
    ].item()
    assert abs(
        selection.set_index("mode_idx").loc[
            1, "frequency_distance_to_target_thz"
        ]
        - 0.2
    ) < 1e-12
    assert output.is_file()
    eigen = pd.read_csv(export_dir / "eigenfrequencies.csv")
    assert eigen["analysis_selected"].tolist() == [False, True]


def test_finite_single_mode_candidate_subset_is_respected(tmp_path):
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame(
        {
            "mode_idx": [0, 1],
            "re": [197.0, 198.2],
            "q": [800.0, 900.0],
            "is_valid": [True, True],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    output = tmp_path / "analysis_mode_selection.csv"

    selected, _selection = select_finite_mode_closest_to_frequency(
        export_dir,
        output,
        target_frequency_thz=198.0,
        candidate_mode_indices=[0],
    )

    assert selected == 0
    assert output.is_file()


@pytest.mark.parametrize(
    ("target_internal_mode", "expected_mode"),
    [("px", 0), ("py", 1), ("dx", 2), ("dy", 3)],
)
def test_finite_gamma_component_selection_uses_requested_mode_weight(
    monkeypatch,
    tmp_path,
    target_internal_mode,
    expected_mode,
):
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame(
        {
            "mode_idx": [0, 1, 2, 3],
            "re": [200.01, 200.02, 200.03, 200.04],
            "q": [1000.0, 900.0, 800.0, 700.0],
            "is_valid": [True, True, True, True],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    component_weights = {
        0: {"mode_px": 0.91, "mode_py": 0.02, "mode_dx": 0.03, "mode_dy": 0.01},
        1: {"mode_px": 0.03, "mode_py": 0.92, "mode_dx": 0.01, "mode_dy": 0.02},
        2: {"mode_px": 0.02, "mode_py": 0.01, "mode_dx": 0.93, "mode_dy": 0.03},
        3: {"mode_px": 0.01, "mode_py": 0.03, "mode_dx": 0.02, "mode_dy": 0.94},
    }

    monkeypatch.setattr(
        lattice_postprocess,
        "sample_finite_cell_fields",
        lambda _export_dir, mode_idx, rho_grid, cells: np.full(
            (len(cells), len(rho_grid["points"])),
            mode_idx + 1,
            dtype=complex,
        ),
    )

    def fake_decompose(_rho, profiles, _weights, _fractions, _index_name):
        mode_idx = int(round(float(np.real(profiles[0, 0])))) - 1
        weights = component_weights[mode_idx]
        dominant = max(weights, key=weights.get).removeprefix("mode_")
        return pd.DataFrame(
            [
                {
                    "subspace_s": 0.01,
                    "subspace_p": weights["mode_px"] + weights["mode_py"],
                    "subspace_d": weights["mode_dx"] + weights["mode_dy"],
                    "subspace_f": 0.01,
                    "mode_s": 0.01,
                    **weights,
                    "mode_f": 0.01,
                    "dominant_subspace": "p" if dominant.startswith("p") else "d",
                    "dominant_subspace_fraction": 0.9,
                    "dominant_mode": dominant,
                    "dominant_mode_fraction": weights[f"mode_{dominant}"],
                    "subspace_matches_mode": True,
                }
            ]
        )

    monkeypatch.setattr(
        lattice_postprocess,
        "decompose_profiles",
        fake_decompose,
    )
    output = tmp_path / "analysis_mode_selection.csv"
    selected, selection = select_finite_fundamental_mode(
        export_dir,
        output,
        cells=[SimpleNamespace(x=0.0, y=0.0, shell=0)],
        period=0.82,
        rho_grid_size=3,
        target_internal_mode=target_internal_mode,
        target_frequency_thz=200.0,
    )

    assert selected == expected_mode
    assert output.is_file()
    selected_row = selection.loc[selection["analysis_selected"]].iloc[0]
    assert selected_row["target_internal_mode"] == target_internal_mode
    assert selected_row["target_component_column"] == f"mode_{target_internal_mode}"
    assert selected_row["selection_reason"] == (
        "selected_centered_nodeless_fundamental_envelope"
    )
    eigen = pd.read_csv(export_dir / "eigenfrequencies.csv")
    assert eigen.loc[eigen["analysis_selected"], "mode_idx"].tolist() == [
        expected_mode
    ]


def test_finite_gamma_component_selection_uses_frequency_only_as_tie_break(
    monkeypatch,
    tmp_path,
):
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame(
        {
            "mode_idx": [0, 1],
            "re": [199.0, 200.1],
            "q": [1000.0, 900.0],
            "is_valid": [True, True],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    monkeypatch.setattr(
        lattice_postprocess,
        "sample_finite_cell_fields",
        lambda _export_dir, mode_idx, rho_grid, cells: np.full(
            (len(cells), len(rho_grid["points"])),
            mode_idx + 1,
            dtype=complex,
        ),
    )
    monkeypatch.setattr(
        lattice_postprocess,
        "decompose_profiles",
        lambda *_args, **_kwargs: pd.DataFrame(
            [
                {
                    "subspace_s": 0.0,
                    "subspace_p": 1.0,
                    "subspace_d": 0.0,
                    "subspace_f": 0.0,
                    "mode_s": 0.0,
                    "mode_px": 0.8,
                    "mode_py": 0.2,
                    "mode_dx": 0.0,
                    "mode_dy": 0.0,
                    "mode_f": 0.0,
                    "dominant_subspace": "p",
                    "dominant_subspace_fraction": 1.0,
                    "dominant_mode": "px",
                    "dominant_mode_fraction": 0.8,
                    "subspace_matches_mode": True,
                }
            ]
        ),
    )

    selected, _selection = select_finite_fundamental_mode(
        export_dir,
        tmp_path / "analysis_mode_selection.csv",
        cells=[SimpleNamespace(x=0.0, y=0.0, shell=0)],
        period=0.82,
        rho_grid_size=3,
        target_internal_mode="px",
        target_frequency_thz=200.0,
    )

    assert selected == 1


def test_finite_envelope_metrics_detect_radial_nodes():
    cells = [
        SimpleNamespace(x=0.0, y=0.0, shell=0),
        SimpleNamespace(x=1.0, y=0.0, shell=1),
        SimpleNamespace(x=2.0, y=0.0, shell=2),
    ]
    carrier = np.array([1.0, 2.0, -1.0, -2.0], dtype=complex)
    nodeless = np.array([1.0, 0.6, 0.2])[:, None] * carrier
    two_nodes = np.array([1.0, -0.6, 0.2])[:, None] * carrier

    ground = finite_mode_envelope_metrics(
        nodeless,
        cells,
        period=1.0,
    )
    excited = finite_mode_envelope_metrics(
        two_nodes,
        cells,
        period=1.0,
    )

    assert ground["peak_shell"] == 0
    assert ground["radial_sign_changes"] == 0
    assert ground["outward_resurgence"] == pytest.approx(0.0)
    assert excited["radial_sign_changes"] == 2


def test_finite_fundamental_selection_rejects_purer_nodal_mode(
    monkeypatch,
    tmp_path,
):
    export_dir = tmp_path / "exports"
    export_dir.mkdir()
    pd.DataFrame(
        {
            "mode_idx": [0, 1],
            "re": [198.0, 198.2],
            "q": [1000.0, 900.0],
            "is_valid": [True, True],
        }
    ).to_csv(export_dir / "eigenfrequencies.csv", index=False)
    cells = [
        SimpleNamespace(x=0.0, y=0.0, shell=0),
        SimpleNamespace(x=1.0, y=0.0, shell=1),
        SimpleNamespace(x=2.0, y=0.0, shell=2),
    ]

    def fake_fields(_export_dir, mode_idx, rho_grid, _cells):
        amplitudes = (
            np.array([1.0, 0.6, 0.2])
            if mode_idx == 0
            else np.array([1.0, -0.6, 0.2])
        )
        return np.repeat(
            amplitudes[:, None],
            len(rho_grid["points"]),
            axis=1,
        ).astype(complex)

    decomposition_weights = iter([0.8, 0.99])

    def fake_decompose(*_args, **_kwargs):
        px = next(decomposition_weights)
        return pd.DataFrame(
            [
                {
                    "subspace_s": 0.0,
                    "subspace_p": 1.0,
                    "subspace_d": 0.0,
                    "subspace_f": 0.0,
                    "mode_s": 0.0,
                    "mode_px": px,
                    "mode_py": 1.0 - px,
                    "mode_dx": 0.0,
                    "mode_dy": 0.0,
                    "mode_f": 0.0,
                    "dominant_subspace": "p",
                    "dominant_subspace_fraction": 1.0,
                    "dominant_mode": "px",
                    "dominant_mode_fraction": px,
                    "subspace_matches_mode": True,
                }
            ]
        )

    monkeypatch.setattr(
        lattice_postprocess,
        "sample_finite_cell_fields",
        fake_fields,
    )
    monkeypatch.setattr(
        lattice_postprocess,
        "decompose_profiles",
        fake_decompose,
    )

    selected, selection = select_finite_fundamental_mode(
        export_dir,
        tmp_path / "analysis_mode_selection.csv",
        cells=cells,
        period=1.0,
        rho_grid_size=3,
        target_internal_mode="px",
        target_frequency_thz=198.1,
    )

    assert selected == 0
    assert selection.set_index("mode_idx").loc[0, "fundamental_eligible"]
    assert not selection.set_index("mode_idx").loc[
        1, "fundamental_eligible"
    ]
    assert selection.set_index("mode_idx").loc[
        1, "target_component_weight"
    ] > selection.set_index("mode_idx").loc[0, "target_component_weight"]


def test_strip_mode_selection_can_ignore_old_geometry_validity(tmp_path):
    eigen_path = tmp_path / "eigenfrequencies.csv"
    eigen_df = pd.DataFrame(
        [
            {"mode_idx": 0, "is_valid": False},
            {"mode_idx": 1, "is_valid": True},
        ]
    )

    assert selected_mode_indices(
        eigen_df,
        None,
        eigen_path,
        respect_validity=False,
    ) == [0, 1]


def test_strip_scoring_uses_strict_gamma_subspace_threshold(tmp_path):
    fourier_dir = tmp_path / "strip_bulk_fourier_hz"
    fourier_dir.mkdir()
    pd.DataFrame(
        [
            {"mode_idx": 0, "frequency": 198.0, "q": 1000.0, "is_valid": True},
            {"mode_idx": 1, "frequency": 199.0, "q": 900.0, "is_valid": False},
        ]
    ).to_csv(fourier_dir / "strip_bulk_fourier_summary.csv", index=False)
    _write_decomposition(
        fourier_dir / "mode_00/p_subspace_mode_decomposition.csv",
        component_name="m",
        gamma_subspace_p=0.9,
    )
    _write_decomposition(
        fourier_dir / "mode_01/p_subspace_mode_decomposition.csv",
        component_name="m",
        gamma_subspace_p=0.9001,
    )

    _, score_df = score_strip_modes(tmp_path)

    assert score_df["mode_idx"].tolist() == [0, 1]
    assert score_df["is_valid"].tolist() == [False, True]
    objective = json.loads((tmp_path / "objective.json").read_text(encoding="utf-8"))
    assert objective["best_mode"]["mode_idx"] == 1
    assert objective["attrs"]["validity_operator"] == ">"
    assert objective["attrs"]["validity_threshold"] == 0.9


def test_finite_scoring_keeps_existing_summary_validity_filter(tmp_path):
    fourier_dir = tmp_path / "finite_lattice_fourier_hz"
    fourier_dir.mkdir()
    pd.DataFrame(
        [
            {"mode_idx": 0, "frequency": 198.0, "q": 1000.0, "is_valid": False},
            {"mode_idx": 1, "frequency": 199.0, "q": 900.0, "is_valid": True},
        ]
    ).to_csv(fourier_dir / "finite_lattice_fourier_summary.csv", index=False)
    _write_decomposition(
        fourier_dir / "mode_01/p_subspace_mode_decomposition.csv",
        component_name="t",
        gamma_subspace_p=0.2,
    )

    _, score_df = score_finite_modes(tmp_path)

    assert score_df["mode_idx"].tolist() == [1]
    assert "is_valid" not in score_df.columns
