#!/usr/bin/env python3
'Run a serial, resumable unit-cell 2D zeta scan.'

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.unit_cell_2d_analysis import (  # noqa: E402
    formal_mode_label,
    unit_cell_2d_series_name,
)
from scripts.run_main.parameter_config import (  # noqa: E402
    PARAMETER_PATH_ENV,
    UNIT_CELL_2D_OUTPUT_ROOT,
    load_shared_parameters,
    parameter_path_from_environment,
)
from scripts.run_main.run_unit_cell_2d import resolved_series_dir  # noqa: E402


DEFAULT_ZETA_VALUES = (1.151, 1.161)
CASE_PARAMETER_FILENAME = "scan_parameter.json"
SCAN_CONFIG_FILENAME = "scan_config.json"
SCAN_SUMMARY_FILENAME = "scan_summary.json"


def _format_decimal(value: float) -> str:
    return format(float(value), ".12g")


def normalize_zeta_values(values: Sequence[float]) -> tuple[float, ...]:
    normalized = tuple(float(value) for value in values)
    if not normalized:
        raise ValueError("At least one zeta value is required.")
    if any(not 0.0 < value < 3.0**0.5 for value in normalized):
        raise ValueError("Each zeta must satisfy 0 < zeta < sqrt(3).")
    if len(set(normalized)) != len(normalized):
        raise ValueError("zeta values must be unique.")
    return normalized


def _read_parameter_data(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("The parameter file must contain one JSON object.")
    settings = data.get("unit_cell_2d")
    if not isinstance(settings, dict):
        raise ValueError("The parameter file has no unit_cell_2d object.")
    return data


def _write_json_atomic(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
        encoding="utf-8",
    )
    temporary.replace(path)


def build_case_parameter_data(
    base_data: dict[str, object],
    zeta: float,
) -> dict[str, object]:
    case_data = deepcopy(base_data)
    settings = case_data["unit_cell_2d"]
    if not isinstance(settings, dict):
        raise ValueError("unit_cell_2d must be a JSON object.")
    settings["zeta"] = float(zeta)
    return case_data


def case_series_dir(case_data: dict[str, object]) -> Path:
    settings = case_data["unit_cell_2d"]
    if not isinstance(settings, dict):
        raise ValueError("unit_cell_2d must be a JSON object.")
    name = unit_cell_2d_series_name(
        b0_nm=float(settings["b0_nm"]),
        eta=float(settings["eta"]),
        zeta=float(settings["zeta"]),
        mesh_size=int(case_data["mesh_size"]),
        q_max_over_g=float(settings["q_max_over_G"]),
        q_points_per_axis=int(settings["q_points_per_axis"]),
        target_bands=tuple(str(item) for item in settings["target_bands"]),
        is_quarter=int(settings.get("isQuarter", 0)),
    )
    return UNIT_CELL_2D_OUTPUT_ROOT / name


def scan_series_dir(
    base_data: dict[str, object],
    zeta_values: Sequence[float],
) -> Path:
    settings = base_data["unit_cell_2d"]
    if not isinstance(settings, dict):
        raise ValueError("unit_cell_2d must be a JSON object.")
    bands = "-".join(formal_mode_label(item) for item in settings["target_bands"])
    zeta_label = "-".join(_format_decimal(value) for value in zeta_values)
    is_quarter = int(settings.get("isQuarter", 0))
    if is_quarter not in {0, 1}:
        raise ValueError("unit_cell_2d.isQuarter must be 0 or 1")
    domain = "quarter" if is_quarter else "full"
    name = (
        "unit_cell_2D_zeta_scan_"
        f"({_format_decimal(settings['b0_nm'])}-"
        f"{_format_decimal(settings['eta'])})"
        f"_mesh{int(base_data['mesh_size'])}"
        f"_kmax{_format_decimal(settings['q_max_over_G'])}"
        f"_uniform{int(settings['q_points_per_axis'])}"
        f"_{domain}_band-{bands}_zeta{zeta_label}"
    )
    return UNIT_CELL_2D_OUTPUT_ROOT / name


def prepare_scan(
    parameter_path: Path,
    zeta_values: Sequence[float],
) -> dict[str, object]:
    parameter_path = Path(parameter_path).resolve()
    values = normalize_zeta_values(zeta_values)
    load_shared_parameters(parameter_path)
    base_data = _read_parameter_data(parameter_path)
    scan_dir = scan_series_dir(base_data, values)
    (scan_dir / "80_logs").mkdir(parents=True, exist_ok=True)
    (scan_dir / "99_config").mkdir(parents=True, exist_ok=True)

    cases = []
    for zeta in values:
        case_data = build_case_parameter_data(base_data, zeta)
        series_dir = case_series_dir(case_data)
        snapshot_path = series_dir / "99_config" / CASE_PARAMETER_FILENAME
        if snapshot_path.exists():
            existing = json.loads(snapshot_path.read_text(encoding="utf-8"))
            if existing != case_data:
                raise RuntimeError(
                    f"Existing case parameter snapshot conflicts: {snapshot_path}"
                )
        else:
            _write_json_atomic(snapshot_path, case_data)
        parsed = load_shared_parameters(snapshot_path)
        if resolved_series_dir(parsed).resolve() != series_dir.resolve():
            raise RuntimeError(f"Case path mismatch for zeta={zeta:g}")
        cases.append(
            {
                "zeta": float(zeta),
                "series_dir": series_dir,
                "parameter_path": snapshot_path,
            }
        )

    prepared = {
        "workflow": "run_unit_cell_2d_zeta_scan",
        "base_parameter_path": parameter_path,
        "scan_dir": scan_dir,
        "zeta_values": values,
        "cases": cases,
    }
    _write_json_atomic(
        scan_dir / "99_config" / SCAN_CONFIG_FILENAME,
        {
            "workflow": prepared["workflow"],
            "base_parameter_path": str(parameter_path),
            "zeta_values": list(values),
            "isQuarter": int(base_data["unit_cell_2d"].get("isQuarter", 0)),
            "cases": [
                {
                    "zeta": case["zeta"],
                    "series_dir": str(case["series_dir"]),
                    "parameter_path": str(case["parameter_path"]),
                }
                for case in cases
            ],
        },
    )
    return prepared


def _write_scan_summary(
    prepared: dict[str, object],
    status: str,
    cases: list[dict[str, object]],
) -> None:
    _write_json_atomic(
        Path(prepared["scan_dir"]) / "99_config" / SCAN_SUMMARY_FILENAME,
        {
            "workflow": prepared["workflow"],
            "status": status,
            "base_parameter_path": str(prepared["base_parameter_path"]),
            "zeta_values": list(prepared["zeta_values"]),
            "cases": cases,
        },
    )


def run_scan(
    parameter_path: Path,
    zeta_values: Sequence[float],
) -> dict[str, object]:
    prepared = prepare_scan(parameter_path, zeta_values)
    case_states = [
        {
            "zeta": case["zeta"],
            "series_dir": str(case["series_dir"]),
            "parameter_path": str(case["parameter_path"]),
            "status": "pending",
        }
        for case in prepared["cases"]
    ]
    _write_scan_summary(prepared, "running", case_states)
    print(f"Unit-cell 2D zeta scan: {prepared['scan_dir']}", flush=True)

    command = [
        sys.executable,
        "-m",
        "scripts.run_main.run_unit_cell_2d",
    ]
    for index, (case, state) in enumerate(
        zip(prepared["cases"], case_states),
        start=1,
    ):
        state["status"] = "running"
        _write_scan_summary(prepared, "running", case_states)
        print(
            f"[{index}/{len(case_states)}] zeta={case['zeta']:.3f} -> "
            f"{case['series_dir']}",
            flush=True,
        )
        environment = os.environ.copy()
        environment[PARAMETER_PATH_ENV] = str(case["parameter_path"])
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            check=False,
        )
        state["returncode"] = int(completed.returncode)
        if completed.returncode != 0:
            state["status"] = "failed"
            _write_scan_summary(prepared, "failed", case_states)
            raise subprocess.CalledProcessError(completed.returncode, command)
        state["status"] = "complete"
        state["run_summary"] = str(
            Path(case["series_dir"]) / "99_config" / "run_summary.json"
        )
        _write_scan_summary(prepared, "running", case_states)

    _write_scan_summary(prepared, "complete", case_states)
    print("Unit-cell 2D zeta scan status: complete", flush=True)
    return {
        "status": "complete",
        "scan_dir": str(prepared["scan_dir"]),
        "cases": case_states,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a serial unit-cell 2D zeta scan."
    )
    parser.add_argument(
        "--zeta",
        nargs="+",
        type=float,
        default=DEFAULT_ZETA_VALUES,
        help="Ordered zeta values; defaults to 1.151 1.161.",
    )
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=parameter_path_from_environment(),
        help="Read all non-zeta settings from this JSON file.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Write parameter snapshots and the scan manifest without COMSOL.",
    )
    return parser


def main(argv: list[str] | None = None) -> dict[str, object]:
    args = build_parser().parse_args(argv)
    if args.prepare_only:
        prepared = prepare_scan(args.parameter_file, args.zeta)
        print(f"Scan output: {prepared['scan_dir']}")
        for case in prepared["cases"]:
            print(
                f"zeta={case['zeta']:.3f}: {case['series_dir']} "
                f"[{case['parameter_path']}]"
            )
        return prepared
    return run_scan(args.parameter_file, args.zeta)


if __name__ == "__main__":
    main()
