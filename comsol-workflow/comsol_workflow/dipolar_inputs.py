"""Explicit saved-data sources, with the historical dipolar task as a preset."""
from dataclasses import asdict, dataclass, fields
import os
from pathlib import Path
import re

from .result_io import read_json


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_MANIFEST_ENV = "COMSOL_WORKFLOW_DIPOLAR_INPUTS"
HISTORICAL_CASE = "f197.980_e15_20260903/symmetry_1_xPEC_yPMC"
HISTORICAL_RUN_ID = "DS-20261006T073408Z-s5-gamma-origin-6c7bd2"


@dataclass(frozen=True)
class DipolarInputs:
    series_dir: Path
    mode_dir: Path
    field_mode_dir: Path
    fourier_source: Path
    output_dir: Path
    reference_dir: Path
    source_config_dir: Path
    finite_model: Path
    boundary_case_dir: Path
    run_id: str
    request_dir: Path
    source_batch: str = "02_gamma_cladding"
    gamma_group: str = "03_gamma_origin"

    def to_metadata(self):
        return {name: str(value) if isinstance(value, Path) else value
                for name, value in asdict(self).items()}


def load_dipolar_inputs(manifest=None):
    """Resolve paths relative to a supplied manifest; do not open scientific data."""
    value = manifest if manifest is not None else os.environ.get(INPUT_MANIFEST_ENV)
    path = Path(value).resolve() if value else None
    data = read_json(path) if path else {}
    if not isinstance(data, dict):
        raise TypeError("Dipolar source manifest must be an object")
    unknown = set(data) - {field.name for field in fields(DipolarInputs)} - {"schema_version"}
    if unknown or (path and (type(data.get("schema_version")) is not int or data["schema_version"] != 1)):
        raise ValueError("Expected source schema_version=1 and supported DipolarInputs fields")
    if "series_dir" in data and not {"mode_dir", "field_mode_dir", "source_config_dir", "finite_model"} <= data.keys():
        raise ValueError("Changing series_dir requires explicit mode_dir, field_mode_dir, source_config_dir and finite_model")
    base = path.parent if path else PROJECT_ROOT

    def resolved(name, default):
        item = data.get(name, str(default))
        if not isinstance(item, str) or not item:
            raise TypeError(f"{name} must be a nonempty path string")
        candidate = Path(item)
        return (candidate if candidate.is_absolute() else base / candidate).resolve()

    series = resolved("series_dir", PROJECT_ROOT / "scripts/.out/finite_cavity/finite_quarter_all_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5")
    mode = resolved("mode_dir", series / "10_overview/symmetry_1_xPEC_yPMC_mode19")
    run_id = data.get("run_id", HISTORICAL_RUN_ID)
    if not isinstance(run_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", run_id):
        raise ValueError("Invalid dipolar run_id")
    output = resolved("output_dir", PROJECT_ROOT / "results/S5_dipolarSingularity_analysis")
    allowed = (PROJECT_ROOT / "results/S5_dipolarSingularity_analysis").resolve()
    if not output.is_relative_to(allowed):
        raise ValueError("Dipolar output_dir must stay under results/S5_dipolarSingularity_analysis")
    groups = {}
    for name, default in (("source_batch", "02_gamma_cladding"), ("gamma_group", "03_gamma_origin")):
        group = data.get(name, default)
        if not isinstance(group, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,96}", group):
            raise ValueError(f"Invalid {name}")
        groups[name] = group
    return DipolarInputs(
        series, mode,
        resolved("field_mode_dir", series / "01_results" / HISTORICAL_CASE / "mode19"),
        resolved("fourier_source", mode / "13_lattice_fourier_Hz/finite_lattice_fourier_hz.npz"),
        output,
        resolved("reference_dir", PROJECT_ROOT / "scripts/.out/unit_cell_2D/S5_reference_b245_eta0.96_zeta1.156_mesh5_k19_20260921"),
        resolved("source_config_dir", series / "99_config/source_runs" / HISTORICAL_CASE),
        resolved("finite_model", series / "00_model" / HISTORICAL_CASE / "finite_quarter.mph"),
        resolved("boundary_case_dir", PROJECT_ROOT / "scripts/.out/unit_cell_2D/S4_b245_eta0.96_zeta1.156_mesh5_Gamma_1case_20260916-gamma-v2/case000_eta0.96_zeta1.156_mesh5"),
        run_id,
        resolved("request_dir", PROJECT_ROOT.parent / "dipolar-exchange/requests" / run_id),
        **groups,
    )
