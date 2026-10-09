#!/usr/bin/env python3
"""Plan or execute the COMSOL dataset for S17 parameter fitting.

The default action is a dry run that writes only a manifest. COMSOL is
imported and started only after an explicit ``--execute`` flag.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DEFAULT_FIT_CONFIG = SCRIPT_DIR / "theory_band_fit.json"
DEFAULT_PARAMETER_CONFIG = SCRIPT_DIR / "parameter.json"
OUTPUT_ROOT = SCRIPT_DIR / ".out" / "unit_cell_band"
DATASET_NAME = "theoryfit_dataset_b245_eta0.94-1.04_target_eta0.96_zeta1.156_mesh5"


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"JSON root must be an object: {path}")
    return value


def make_ray(name: str, coordinates: list[list[float]]) -> list[dict[str, object]]:
    """Build the point dictionaries consumed by ``run_band_pair.run_case``."""

    if not coordinates or coordinates[0] != [0.0, 0.0]:
        raise ValueError(f"Ray {name!r} must start at Gamma")
    norms = [math.hypot(float(x), float(y)) for x, y in coordinates]
    endpoint = max(norms) or 1.0
    points = []
    for x, y in coordinates:
        kx, ky = float(x), float(y)
        norm = math.hypot(kx, ky)
        points.append({
            "direction": name,
            "kx": kx,
            "ky": ky,
            "kx_str": f"{kx:.6f}",
            "ky_str": f"{ky:.6f}",
            "k_norm": norm,
            "path_fraction": norm / endpoint,
        })
    return points


def sampling_branches(config: dict[str, Any], *, include_validation: bool) -> dict[str, list[dict[str, object]]]:
    sampling = config["sampling"]
    magnitudes = [float(value) for value in sampling["fit_axis_magnitudes_over_G"]]
    branches = {
        "fit_x": make_ray("fit_x", [[value, 0.0] for value in magnitudes]),
        "fit_y": make_ray("fit_y", [[0.0, value] for value in magnitudes]),
    }
    if include_validation:
        for index, point in enumerate(sampling["validation_points_over_G"], start=1):
            branches[f"validation_{index}"] = make_ray(
                f"validation_{index}", [[0.0, 0.0], list(point)]
            )
        for index, point in enumerate(sampling["qa_points_over_G"], start=1):
            branches[f"qa_{index}"] = make_ray(
                f"qa_{index}", [[0.0, 0.0], list(point)]
            )
    return branches


def build_manifest(
    fit_config: dict[str, Any], parameter_config: dict[str, Any], output_dir: Path,
) -> dict[str, Any]:
    from scripts.run_main.parameter_config import parameter_data_to_flat

    resolved_parameters = parameter_data_to_flat(parameter_config)
    calibration = fit_config["calibration"]
    target = fit_config["target"]
    structures = []
    for eta in calibration["eta_values"]:
        structures.append({
            "structure_id": f"calibration_eta{float(eta):.3f}",
            "role": "calibration",
            "cell": {
                "b0_nm": float(calibration["b0_nm"]),
                "eta": float(eta),
                "zeta": float(calibration["zeta"]),
            },
            "relative_output": f"calibration_zeta1/eta{float(eta):.3f}",
            "branches": sampling_branches(fit_config, include_validation=False),
        })
    structures.append({
        "structure_id": "target_eta0.960_zeta1.156",
        "role": "target",
        "cell": {key: float(value) for key, value in target.items()},
        "relative_output": "target_zeta1.156",
        "branches": sampling_branches(fit_config, include_validation=True),
    })
    return {
        "schema_version": 1,
        "status": "dry_run",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "will_start_comsol": False,
        "output_dir": str(output_dir.resolve()),
        "model": fit_config["model"],
        "basis_mapping": fit_config["basis_mapping"],
        "units": {
            "frequency": fit_config["frequency_units"],
            "velocity": fit_config["velocity_units"],
            "input_wavevector": fit_config["wavevector_input_units"],
            "physical_wavevector": "1/um",
            "G_inv_um": 4.0 * math.pi / (
                math.sqrt(3.0) * float(fit_config["lattice_constant_um"])
            ),
        },
        "geometry": {
            "lattice_constant_um": float(fit_config["lattice_constant_um"]),
            "minimum_clearance_um": 0.02,
        },
        "simulation": {
            "mesh_auto_size": int(resolved_parameters["mesh_size"]),
            "eigenmode_count": int(resolved_parameters["unit_cell_eigenmode_count"]),
            "eigenfrequency_shift": "c_const/1.55[um]",
            "center_frequency_thz_snapshot": float(resolved_parameters["center_frequency_thz"]),
            "model_reuse": "one COMSOL Model per fixed eta geometry",
        },
        "fit_config_snapshot": fit_config,
        "parameter_config_snapshot": parameter_config,
        "structures": structures,
    }


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def execute_manifest(manifest: dict[str, Any], manifest_path: Path) -> None:
    """Execute all structures; this is the only function importing COMSOL workflow code."""

    from scripts.run_main import run_band_pair as band
    from scripts.run_main.parameter_config import CellParameters

    output_dir = Path(manifest["output_dir"])
    manifest["status"] = "running"
    manifest["will_start_comsol"] = True
    write_manifest(manifest_path, manifest)
    for structure in manifest["structures"]:
        cell = CellParameters.from_mapping(structure["cell"], structure["structure_id"])
        case_output = output_dir / structure["relative_output"]
        band.run_case(
            "cavity",
            cell.to_fourier_params(structure["structure_id"]),
            case_output,
            branches=structure["branches"],
        )
    manifest["status"] = "complete"
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    write_manifest(manifest_path, manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_FIT_CONFIG)
    parser.add_argument("--parameters", type=Path, default=DEFAULT_PARAMETER_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT / DATASET_NAME)
    parser.add_argument(
        "--execute", action="store_true",
        help="Explicitly start COMSOL after writing and printing the manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fit_config = _load_json(args.config)
    parameter_config = _load_json(args.parameters)
    manifest_path = args.output_dir / "99_config" / "dry_run_manifest.json"
    if args.execute and manifest_path.exists():
        manifest = _load_json(manifest_path)
        if manifest.get("status") not in {"dry_run", "running"}:
            raise FileExistsError(
                f"Refusing to execute an already finalized dataset: {manifest_path}"
            )
    else:
        manifest = build_manifest(fit_config, parameter_config, args.output_dir)
    write_manifest(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"Manifest: {manifest_path.resolve()}")
    if args.execute:
        execute_manifest(manifest, manifest_path)
    else:
        print("Dry run only: COMSOL was not imported or started. Use --execute after review.")


if __name__ == "__main__":
    main()
