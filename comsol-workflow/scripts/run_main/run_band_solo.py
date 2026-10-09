#!/usr/bin/env python3
"""Calculate the cavity unit-cell band without a cladding comparison."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main import run_band_pair as band  # noqa: E402
from scripts.run_main.parameter_config import (  # noqa: E402
    ACTIVE_PARAMETERS,
    UNIT_CELL_OUTPUT_ROOT,
    update_center_frequency_from_cavity_p2,
)


def solo_series_name() -> str:
    """Return the cavity-only series identity derived from shared parameters."""
    cavity_label = ACTIVE_PARAMETERS.cavity.label("")
    return f"unit_cell_{cavity_label}_mesh{ACTIVE_PARAMETERS.mesh_size}"


OUT_DIR = UNIT_CELL_OUTPUT_ROOT / solo_series_name()
CAVITY_PARAMS = ACTIVE_PARAMETERS.cavity.to_fourier_params("cavity_p_bic")
EIGENMODE_COUNT = ACTIVE_PARAMETERS.unit_cell_eigenmode_count
MESH_AUTO_SIZE = ACTIVE_PARAMETERS.mesh_size
UPDATE_CENTER_FREQUENCY_AFTER_RUN = True


def _solo_parameter_metadata() -> dict[str, object]:
    return {
        "cavity": ACTIVE_PARAMETERS.cavity.to_dict(),
        "mesh_size": MESH_AUTO_SIZE,
        "unit_cell_eigenmode_count": EIGENMODE_COUNT,
    }


def run_band_solo(
    output_dir: Path,
    cavity_params: dict[str, float],
) -> dict[str, object]:
    """Run only the cavity unit cell and write single-cell band/Q outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    branches = band.case_k_points("cavity")
    band.write_json_atomic(
        output_dir / band.CONFIG_DIRNAME / "config.json",
        {
            "workflow": "run_band_solo",
            "execution_engine": "reused_model",
            "model_reuse": {
                "scope": "one COMSOL Model for the cavity unit-cell geometry",
                "same_mesh": "reuse without mesh rebuild",
                "changed_mesh": "rebuild mesh1 only",
            },
            "geometry": {
                "a": band.A,
                "r_0": band.R_0,
                "b_0": band.B_0,
                "d": band.D,
            },
            "simulation": {
                "slab_height": band.SLAB_HEIGHT,
                "refractive_index": band.REFRACTIVE_INDEX,
                "eigenmode_count": EIGENMODE_COUNT,
                "mesh_auto_size": MESH_AUTO_SIZE,
                "mode_type": band.MODE_TYPE,
                "wavelength": band.WAVELENGTH,
                "eigenfrequency_shift": band.EIGENFREQUENCY_SHIFT,
            },
            "frequency_window_thz": [
                band.FREQUENCY_MIN_THZ,
                band.FREQUENCY_MAX_THZ,
            ],
            "band_points_per_arm": band.BAND_POINTS_PER_ARM,
            "q_fit_k_magnitudes": band.Q_FIT_K_MAGNITUDES,
            "band_frequency_tolerance_thz": band.BAND_FREQUENCY_TOLERANCE_THZ,
            "band_minimum_overlap": band.BAND_MINIMUM_OVERLAP,
            "case": {"cavity": cavity_params},
            "parameters": _solo_parameter_metadata(),
            "sampled_points": {"cavity": branches},
        },
    )
    cavity = band.run_case("cavity", cavity_params, output_dir)
    cavity_fits = band._build_q_fits(cavity)
    cavity_overview = band.unit_cell_case_paths(output_dir, "cavity")["overview"]
    overview_dir = output_dir / band.OVERVIEW_DIRNAME
    band.write_csv_atomic(cavity_overview / "q_power_law_fits.csv", cavity_fits)
    band.plot_band_cases(
        (("Cavity cell", cavity),),
        overview_dir / "unit_cell_band.png",
    )
    band.plot_p_band_q_cases(
        (("Cavity cell", cavity, cavity_fits),),
        overview_dir / "p_bands_q_vs_k.png",
    )
    return {
        "output_dir": output_dir,
        "cavity": cavity,
        "cavity_fits": cavity_fits,
    }


def cavity_p2_gamma_frequency(cavity: pd.DataFrame) -> float:
    """Expose the shared strict p2@Gamma extraction for the solo workflow."""
    return band.cavity_p2_gamma_frequency(cavity)


def main() -> None:
    result = run_band_solo(OUT_DIR, CAVITY_PARAMS)
    if UPDATE_CENTER_FREQUENCY_AFTER_RUN:
        frequency_thz = cavity_p2_gamma_frequency(result["cavity"])
        updated = update_center_frequency_from_cavity_p2(
            frequency_thz,
            expected_parameters=ACTIVE_PARAMETERS,
            workflow="run_band_solo",
            include_cladding=False,
        )
        print(
            "Updated shared center frequency from cavity p2@Gamma: "
            f"{updated.center_frequency_thz:.12g} THz -> {updated.source_path}"
        )


if __name__ == "__main__":
    main()
