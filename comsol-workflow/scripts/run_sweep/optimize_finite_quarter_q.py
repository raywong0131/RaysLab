#!/usr/bin/env python3
"""Optimize the user-visible finite-quarter py-mode Q without field plotting."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Callable, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_PATH = Path(__file__).resolve()
SCRIPTS_DIR = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPTS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_main.parameter_config import (  # noqa: E402
    FINITE_CAVITY_OUTPUT_ROOT,
    PARAMETER_PATH,
    PARAMETER_PATH_ENV,
    SharedParameters,
    load_shared_parameters,
)


EXPECTED_CAVITY = {"b0_nm": 245.0, "eta": 0.960, "zeta": 1.155}
START_CLADDING = {"b0_nm": 243.7, "eta": 0.980, "zeta": 0.928}
EXPECTED_CAVITY_LAYERS = 20
EXPECTED_CLADDING_LAYERS = 20
EXPECTED_MESH_SIZE = 9
FIXED_ETA = 0.980
FIXED_SHIFT = 0.0
SYMMETRY_ID = 1
TARGET_Q = 6700.0
STRETCH_Q = 7000.0
ZETA_BOUNDS = (0.912, 0.944)
ZETA_COARSE_STEP = 0.004
ZETA_FINE_STEP = 0.001
B0_BOUNDS_NM = (241.7, 245.7)
B0_COARSE_STEP_NM = 0.5
B0_FINE_STEP_NM = 0.1
DEFAULT_MAX_NEW_CANDIDATES = 20
Q_IMPROVEMENT_TOLERANCE = 1e-9

CSV_COLUMNS = [
    "iteration",
    "stage",
    "source",
    "parameter_key",
    "b0_nm",
    "eta",
    "zeta",
    "symmetry_id",
    "mode0_frequency_thz",
    "mode0_q",
    "mode0_is_valid",
    "mode1_frequency_thz",
    "mode1_q",
    "mode1_is_valid",
    "selected_mode_idx",
    "py_frequency_thz",
    "py_q",
    "is_new_best",
    "best_q_so_far",
    "elapsed_seconds",
]


@dataclass(frozen=True)
class Candidate:
    b0_nm: float
    eta: float
    zeta: float

    @property
    def key(self) -> str:
        return candidate_key(self.b0_nm, self.eta, self.zeta)


class SearchBudgetExhausted(RuntimeError):
    pass


def _format_decimal(value: float, places: int = 3) -> str:
    return f"{float(value):.{places}f}".rstrip("0").rstrip(".")


def candidate_key(b0_nm: float, eta: float, zeta: float) -> str:
    return f"b0={float(b0_nm):.3f}|eta={float(eta):.3f}|zeta={float(zeta):.3f}"


def optimization_series_name(parameters: SharedParameters) -> str:
    return (
        "finite_quarter_q_optimization_"
        f"{parameters.cavity_layers}-{parameters.cladding_layers}_"
        f"{parameters.cavity.label('cav')}_"
        f"clad_start({_format_decimal(parameters.cladding.b0_nm)}"
        f"-{_format_decimal(parameters.cladding.eta)}"
        f"-{_format_decimal(parameters.cladding.zeta)})_"
        f"mesh{parameters.mesh_size}_sym{SYMMETRY_ID}"
    )


def validate_approved_parameters(parameters: SharedParameters) -> None:
    actual_cavity = parameters.cavity.to_dict()
    actual_cladding = parameters.cladding.to_dict()
    for field, expected in EXPECTED_CAVITY.items():
        if not math.isclose(actual_cavity[field], expected, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                f"Approved cavity {field}={expected} but parameter file has "
                f"{actual_cavity[field]}"
            )
    for field, expected in START_CLADDING.items():
        if not math.isclose(actual_cladding[field], expected, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                f"Approved starting cladding {field}={expected} but parameter file has "
                f"{actual_cladding[field]}"
            )
    if parameters.cavity_layers != EXPECTED_CAVITY_LAYERS:
        raise ValueError("This optimization requires cavity_layers=20")
    if parameters.cladding_layers != EXPECTED_CLADDING_LAYERS:
        raise ValueError("This optimization requires cladding_layers=20")
    if parameters.mesh_size != EXPECTED_MESH_SIZE:
        raise ValueError("This optimization requires mesh_size=9")
    if tuple(parameters.cladding_shift_factors) != (FIXED_SHIFT,):
        raise ValueError("This optimization requires cladding_shift_factors=[0.0]")
    if parameters.cladding_shift_profile.kind != "uniform":
        raise ValueError("This optimization requires a uniform cladding shift profile")


def build_candidate_parameter_data(
    baseline_data: dict[str, object],
    candidate: Candidate,
) -> dict[str, object]:
    data = json.loads(json.dumps(baseline_data))
    cells = data if data["schema_version"] == 1 else data["cells"]
    structure = data if data["schema_version"] == 1 else data["structure"]
    cells["cladding"] = {
        "b0_nm": float(candidate.b0_nm),
        "eta": float(candidate.eta),
        "zeta": float(candidate.zeta),
    }
    structure["cladding_shift_factors"] = [FIXED_SHIFT]
    structure["cladding_shift_profile"] = {"kind": "uniform"}
    return data


def _csv_bool(value: object) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def normalize_modes(raw_modes: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    modes = []
    for index, raw in enumerate(raw_modes):
        mode_idx = int(raw.get("mode_idx", index))
        modes.append(
            {
                "mode_idx": mode_idx,
                "frequency_thz": float(raw.get("frequency_thz", raw.get("re"))),
                "q": float(raw["q"]),
                "is_valid": _csv_bool(raw.get("is_valid", True)),
            }
        )
    modes.sort(key=lambda mode: int(mode["mode_idx"]))
    if len(modes) != 2:
        raise ValueError(
            f"SYMMETRY_ID={SYMMETRY_ID} must return exactly two modes; got {len(modes)}"
        )
    if [int(mode["mode_idx"]) for mode in modes] != [0, 1]:
        raise ValueError("Expected the two COMSOL solutions to be mode0 and mode1")
    return modes


def select_high_q_py_mode(modes: Iterable[dict[str, object]]) -> dict[str, object]:
    normalized = normalize_modes(modes)
    valid = [mode for mode in normalized if bool(mode["is_valid"])]
    if not valid:
        raise ValueError("Neither SYMMETRY_ID=1 solution passed the validity check")
    return max(valid, key=lambda mode: float(mode["q"]))


def read_eigenfrequency_csv(path: Path) -> list[dict[str, object]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return normalize_modes(rows)


def _numeric_row(row: dict[str, object]) -> dict[str, object]:
    converted = dict(row)
    for field in (
        "iteration",
        "symmetry_id",
        "selected_mode_idx",
    ):
        converted[field] = int(float(converted[field]))
    for field in (
        "b0_nm",
        "eta",
        "zeta",
        "mode0_frequency_thz",
        "mode0_q",
        "mode1_frequency_thz",
        "mode1_q",
        "py_frequency_thz",
        "py_q",
        "best_q_so_far",
        "elapsed_seconds",
    ):
        converted[field] = float(converted[field])
    for field in ("mode0_is_valid", "mode1_is_valid", "is_new_best"):
        converted[field] = _csv_bool(converted[field])
    return converted


def read_history(path: Path) -> list[dict[str, object]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return [_numeric_row(row) for row in csv.DictReader(stream)]


def write_history(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def result_row(
    *,
    iteration: int,
    stage: str,
    source: str,
    candidate: Candidate,
    modes: list[dict[str, object]],
    previous_best_q: float,
    elapsed_seconds: float,
) -> dict[str, object]:
    modes = normalize_modes(modes)
    selected = select_high_q_py_mode(modes)
    py_q = float(selected["q"])
    is_new_best = py_q > previous_best_q + Q_IMPROVEMENT_TOLERANCE
    return {
        "iteration": int(iteration),
        "stage": stage,
        "source": source,
        "parameter_key": candidate.key,
        "b0_nm": float(candidate.b0_nm),
        "eta": float(candidate.eta),
        "zeta": float(candidate.zeta),
        "symmetry_id": SYMMETRY_ID,
        "mode0_frequency_thz": float(modes[0]["frequency_thz"]),
        "mode0_q": float(modes[0]["q"]),
        "mode0_is_valid": bool(modes[0]["is_valid"]),
        "mode1_frequency_thz": float(modes[1]["frequency_thz"]),
        "mode1_q": float(modes[1]["q"]),
        "mode1_is_valid": bool(modes[1]["is_valid"]),
        "selected_mode_idx": int(selected["mode_idx"]),
        "py_frequency_thz": float(selected["frequency_thz"]),
        "py_q": py_q,
        "is_new_best": is_new_best,
        "best_q_so_far": max(previous_best_q, py_q),
        "elapsed_seconds": float(elapsed_seconds),
    }


def best_row(rows: Iterable[dict[str, object]]) -> dict[str, object]:
    rows = list(rows)
    if not rows:
        raise ValueError("No completed Q evaluations are available")
    return max(rows, key=lambda row: float(row["py_q"]))


def directional_axis_search(
    *,
    start_value: float,
    coarse_step: float,
    fine_step: float,
    bounds: tuple[float, float],
    evaluate_value: Callable[[float, str], dict[str, object]],
    stage_prefix: str,
) -> dict[str, object]:
    """Directionally search one scalar axis, then refine around its best value."""

    def evaluate(value: float, stage: str) -> dict[str, object]:
        value = round(float(value), 9)
        if not bounds[0] - 1e-12 <= value <= bounds[1] + 1e-12:
            raise ValueError(f"Search value {value} lies outside {bounds}")
        return evaluate_value(value, stage)

    center = evaluate(start_value, f"{stage_prefix}_center")
    minus_value = max(bounds[0], start_value - coarse_step)
    plus_value = min(bounds[1], start_value + coarse_step)
    coarse_rows = [center]
    if not math.isclose(minus_value, start_value, abs_tol=1e-12):
        coarse_rows.append(evaluate(minus_value, f"{stage_prefix}_coarse"))
    if not math.isclose(plus_value, start_value, abs_tol=1e-12):
        coarse_rows.append(evaluate(plus_value, f"{stage_prefix}_coarse"))
    current = best_row(coarse_rows)
    current_value = float(
        current["zeta"] if stage_prefix.startswith("zeta") else current["b0_nm"]
    )
    if current is not center:
        direction = 1.0 if current_value > start_value else -1.0
        while True:
            next_value = round(current_value + direction * coarse_step, 9)
            if next_value < bounds[0] - 1e-12 or next_value > bounds[1] + 1e-12:
                break
            trial = evaluate(next_value, f"{stage_prefix}_coarse")
            if float(trial["py_q"]) <= float(current["py_q"]) + Q_IMPROVEMENT_TOLERANCE:
                break
            current = trial
            current_value = next_value

    fine_rows = [current]
    for sign in (-1.0, 1.0):
        value = round(current_value + sign * fine_step, 9)
        if bounds[0] - 1e-12 <= value <= bounds[1] + 1e-12:
            fine_rows.append(evaluate(value, f"{stage_prefix}_fine"))
    refined = best_row(fine_rows)
    refined_value = float(
        refined["zeta"] if stage_prefix.startswith("zeta") else refined["b0_nm"]
    )
    if refined is not current:
        direction = 1.0 if refined_value > current_value else -1.0
        fine_limit = coarse_step + 1e-12
        while abs(refined_value - current_value) < fine_limit:
            next_value = round(refined_value + direction * fine_step, 9)
            if next_value < bounds[0] - 1e-12 or next_value > bounds[1] + 1e-12:
                break
            trial = evaluate(next_value, f"{stage_prefix}_fine")
            if float(trial["py_q"]) <= float(refined["py_q"]) + Q_IMPROVEMENT_TOLERANCE:
                break
            refined = trial
            refined_value = next_value
    return refined


def plot_q_history(rows: list[dict[str, object]], path: Path) -> None:
    ordered = sorted(rows, key=lambda row: int(row["iteration"]))
    iterations = [int(row["iteration"]) for row in ordered]
    q_values = [float(row["py_q"]) for row in ordered]
    best = best_row(ordered)
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=True)
    ax.plot(iterations, q_values, linestyle="none", marker="o", color="#2864dc")
    ax.plot(iterations, q_values, linewidth=1.1, alpha=0.35, color="#2864dc")
    ax.axhline(TARGET_Q, color="#d97706", linestyle="--", label="Q = 6700")
    ax.axhline(STRETCH_Q, color="#b91c1c", linestyle=":", label="Q = 7000")
    ax.scatter(
        [int(best["iteration"])],
        [float(best["py_q"])],
        marker="*",
        s=150,
        color="#b91c1c",
        zorder=4,
        label="best",
    )
    ax.set_title("Finite-Quarter Py-Mode Q Optimization")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Py-mode Q")
    ax.grid(True, alpha=0.25)
    ax.legend()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def worker_main(args: argparse.Namespace) -> int:
    parameter_file = Path(args.parameter_file).resolve()
    os.environ[PARAMETER_PATH_ENV] = str(parameter_file)

    # These imports must happen after the task-local parameter path is installed.
    from comsol_workflow.simulation_utils import BoundarySpec, SimulationRun
    from scripts.run_main import run_finite as finite
    from scripts.run_main import run_finite_quarter as quarter

    parameters = load_shared_parameters(parameter_file)
    if parameters.cavity_layers != EXPECTED_CAVITY_LAYERS:
        raise ValueError("Worker requires cavity_layers=20")
    if parameters.cladding_layers != EXPECTED_CLADDING_LAYERS:
        raise ValueError("Worker requires cladding_layers=20")
    if parameters.mesh_size != EXPECTED_MESH_SIZE:
        raise ValueError("Worker requires mesh9")
    if not math.isclose(parameters.cladding.eta, FIXED_ETA, abs_tol=1e-12):
        raise ValueError("Worker requires fixed cladding eta=0.980")
    if tuple(parameters.cladding_shift_factors) != (FIXED_SHIFT,):
        raise ValueError("Worker requires cladding shift 0.000")
    if parameters.cladding_shift_profile.kind != "uniform":
        raise ValueError("Worker requires uniform cladding shift")

    work_dir = Path(args.work_dir).resolve()
    result_path = Path(args.result_path).resolve()
    progress_log = Path(args.progress_log).resolve()
    model_path = work_dir / "finite_quarter.mph"
    work_dir.mkdir(parents=True, exist_ok=True)
    progress_log.parent.mkdir(parents=True, exist_ok=True)

    finite.RUN_FARFIELD_FFT = False
    finite.RUN_FINITE_LATTICE_FOURIER_POSTPROCESS = False
    finite.set_cladding_inward_shift_factor(FIXED_SHIFT)
    case = quarter.symmetry_case(SYMMETRY_ID)
    start_time = time.perf_counter()
    (
        _full_records,
        quarter_records,
        _geometry_metadata,
        _quarter_boundary,
        finite_hex_a,
    ) = quarter.build_quarter_case_geometry()
    layers = finite.simulation_layers(finite.polygons_from_records(quarter_records))
    validity_boundary = finite.strip_validity_boundary()

    with SimulationRun(
        config=quarter.quarter_simulation_config(case),
        progress_log_path=progress_log,
    ) as simulation:
        if model_path.is_file():
            finite.attach_saved_model(simulation, model_path)
            finite.validate_geometry_checkpoint_model(simulation)
        else:
            simulation.build_geometry(
                BoundarySpec("quarter_hexagon_boundary", a=finite_hex_a),
                layers,
                simulation_mode="finite_quarter",
                k=None,
            )
            simulation.model.save(str(model_path))
        simulation.run_simulation(
            mesh_save_path=model_path,
            mesh_builder=finite.finite_mesh_builder(layers),
        )
        simulation.model.save(str(model_path))
        eigenfrequencies = simulation.get_eigenfrequencies()
        validity = finite.mode_validity_for_boundary(
            simulation,
            len(eigenfrequencies),
            validity_boundary,
        )

    modes = normalize_modes(
        [
            {
                "mode_idx": idx,
                "frequency_thz": values[0],
                "q": values[2],
                "is_valid": validity[idx]["is_valid"],
            }
            for idx, values in enumerate(eigenfrequencies)
        ]
    )
    selected = select_high_q_py_mode(modes)
    _write_json(
        result_path,
        {
            "candidate": parameters.cladding.to_dict(),
            "symmetry_id": SYMMETRY_ID,
            "modes": modes,
            "selected_mode": selected,
            "elapsed_seconds": time.perf_counter() - start_time,
            "model_path": str(model_path),
            "progress_log": str(progress_log),
        },
    )
    print(
        f"worker complete: b0={parameters.cladding.b0_nm:.3f}, "
        f"eta={parameters.cladding.eta:.3f}, zeta={parameters.cladding.zeta:.3f}, "
        f"selected mode={selected['mode_idx']}, Q={selected['q']:.6f}",
        flush=True,
    )
    return 0


class OptimizationCoordinator:
    def __init__(
        self,
        parameters: SharedParameters,
        baseline_data: dict[str, object],
        output_root: Path,
        *,
        max_new_candidates: int,
    ) -> None:
        self.parameters = parameters
        self.baseline_data = baseline_data
        self.output_root = output_root
        self.max_new_candidates = int(max_new_candidates)
        self.overview_dir = output_root / "10_overview"
        self.model_dir = output_root / "00_model"
        self.logs_dir = output_root / "80_logs"
        self.config_dir = output_root / "99_config"
        self.staging_dir = output_root / ".staging"
        self.history_path = self.overview_dir / "q_optimization.csv"
        self.plot_path = self.overview_dir / "q_vs_iteration.png"
        self.rows = read_history(self.history_path)
        self._row_by_key = {str(row["parameter_key"]): row for row in self.rows}
        self.new_candidate_count = sum(
            1 for row in self.rows if str(row["source"]) == "worker"
        )

    def prepare_output(self) -> None:
        for directory in (
            self.overview_dir,
            self.model_dir,
            self.logs_dir,
            self.config_dir,
            self.staging_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)
        config_path = self.config_dir / "optimization_config.json"
        config = {
            "workflow": "finite_quarter_py_q_optimization",
            "approved_baseline": self.parameters.to_metadata(),
            "fixed": {
                "cladding_eta": FIXED_ETA,
                "cladding_shift": FIXED_SHIFT,
                "cladding_shift_profile": {"kind": "uniform"},
                "symmetry_id": SYMMETRY_ID,
                "mesh_size": EXPECTED_MESH_SIZE,
                "cavity_layers": EXPECTED_CAVITY_LAYERS,
                "cladding_layers": EXPECTED_CLADDING_LAYERS,
            },
            "search": {
                "priority": ["zeta", "b0_nm"],
                "zeta_bounds": list(ZETA_BOUNDS),
                "zeta_coarse_step": ZETA_COARSE_STEP,
                "zeta_fine_step": ZETA_FINE_STEP,
                "b0_bounds_nm": list(B0_BOUNDS_NM),
                "b0_coarse_step_nm": B0_COARSE_STEP_NM,
                "b0_fine_step_nm": B0_FINE_STEP_NM,
                "max_new_candidates": self.max_new_candidates,
                "target_q": TARGET_Q,
                "stretch_q": STRETCH_Q,
            },
            "mode_selection": (
                "exactly two SYMMETRY_ID=1 solutions; select the valid higher-Q mode"
            ),
            "plotting_during_worker": False,
        }
        if config_path.is_file():
            existing = json.loads(config_path.read_text(encoding="utf-8"))
            if existing != config:
                raise RuntimeError(
                    f"Existing optimization configuration differs: {config_path}"
                )
        else:
            _write_json(config_path, config)

    def existing(self, candidate: Candidate) -> dict[str, object] | None:
        return self._row_by_key.get(candidate.key)

    def _previous_best_q(self) -> float:
        if not self.rows:
            return -math.inf
        return float(best_row(self.rows)["py_q"])

    def add_baseline(self, modes: list[dict[str, object]]) -> dict[str, object]:
        candidate = Candidate(**START_CLADDING)
        existing = self.existing(candidate)
        if existing is not None:
            return existing
        row = result_row(
            iteration=0,
            stage="baseline",
            source="existing_baseline",
            candidate=candidate,
            modes=modes,
            previous_best_q=self._previous_best_q(),
            elapsed_seconds=0.0,
        )
        self._append_row(row)
        return row

    def _append_row(self, row: dict[str, object]) -> None:
        self.rows.append(row)
        self.rows.sort(key=lambda item: int(item["iteration"]))
        self._row_by_key[str(row["parameter_key"])] = row
        write_history(self.history_path, self.rows)
        plot_q_history(self.rows, self.plot_path)

    def _promote_candidate(
        self,
        *,
        row: dict[str, object],
        candidate: Candidate,
        model_path: Path,
        progress_log: Path,
        parameter_path: Path,
    ) -> None:
        best_model_path = self.model_dir / "finite_quarter.mph"
        temporary_model = self.model_dir / "finite_quarter.mph.tmp"
        temporary_model.unlink(missing_ok=True)
        shutil.move(str(model_path), str(temporary_model))
        os.replace(temporary_model, best_model_path)
        if progress_log.is_file():
            shutil.copy2(progress_log, self.model_dir / "comsol_progress.log")
        shutil.copy2(parameter_path, self.config_dir / "best_parameter.json")
        _write_json(
            self.config_dir / "best_parameters.json",
            {
                "iteration": int(row["iteration"]),
                "cladding": {
                    "b0_nm": candidate.b0_nm,
                    "eta": candidate.eta,
                    "zeta": candidate.zeta,
                },
                "symmetry_id": SYMMETRY_ID,
                "selected_mode_idx": int(row["selected_mode_idx"]),
                "frequency_thz": float(row["py_frequency_thz"]),
                "q": float(row["py_q"]),
                "model_path": str(best_model_path),
            },
        )

    def evaluate(self, candidate: Candidate, stage: str) -> dict[str, object]:
        if not math.isclose(candidate.eta, FIXED_ETA, abs_tol=1e-12):
            raise ValueError("Optimization attempted to change fixed eta")
        existing = self.existing(candidate)
        if existing is not None:
            print(f"Reusing {candidate.key}: Q={float(existing['py_q']):.6f}", flush=True)
            return existing
        if self.new_candidate_count >= self.max_new_candidates:
            raise SearchBudgetExhausted(
                f"Reached maximum of {self.max_new_candidates} new candidates"
            )

        iteration = max((int(row["iteration"]) for row in self.rows), default=0) + 1
        candidate_dir = self.staging_dir / f"iteration{iteration:03d}"
        if candidate_dir.exists():
            shutil.rmtree(candidate_dir)
        candidate_dir.mkdir(parents=True)
        parameter_path = candidate_dir / "parameter.json"
        result_path = candidate_dir / "result.json"
        progress_log = self.logs_dir / f"iteration{iteration:03d}_comsol_progress.log"
        stdout_log = self.logs_dir / f"iteration{iteration:03d}_stdout.log"
        stderr_log = self.logs_dir / f"iteration{iteration:03d}_stderr.log"
        candidate_data = build_candidate_parameter_data(self.baseline_data, candidate)
        _write_json(parameter_path, candidate_data)

        command = [
            sys.executable,
            str(SCRIPT_PATH),
            "--worker",
            "--parameter-file",
            str(parameter_path),
            "--work-dir",
            str(candidate_dir),
            "--result-path",
            str(result_path),
            "--progress-log",
            str(progress_log),
        ]
        print(
            f"Iteration {iteration}: stage={stage}, b0={candidate.b0_nm:.3f}, "
            f"eta={candidate.eta:.3f}, zeta={candidate.zeta:.3f}",
            flush=True,
        )
        start_time = time.perf_counter()
        environment = os.environ.copy()
        environment[PARAMETER_PATH_ENV] = str(parameter_path.resolve())
        with stdout_log.open("w", encoding="utf-8") as stdout_stream, stderr_log.open(
            "w", encoding="utf-8"
        ) as stderr_stream:
            process = subprocess.run(
                command,
                cwd=REPO_ROOT,
                env=environment,
                stdout=stdout_stream,
                stderr=stderr_stream,
                check=False,
            )
        elapsed = time.perf_counter() - start_time
        if process.returncode != 0 or not result_path.is_file():
            raise RuntimeError(
                f"Iteration {iteration} worker failed with code {process.returncode}; "
                f"see {stderr_log}"
            )
        result = json.loads(result_path.read_text(encoding="utf-8"))
        previous_best_q = self._previous_best_q()
        row = result_row(
            iteration=iteration,
            stage=stage,
            source="worker",
            candidate=candidate,
            modes=result["modes"],
            previous_best_q=previous_best_q,
            elapsed_seconds=elapsed,
        )
        self._append_row(row)
        self.new_candidate_count += 1
        model_path = candidate_dir / "finite_quarter.mph"
        if bool(row["is_new_best"]):
            self._promote_candidate(
                row=row,
                candidate=candidate,
                model_path=model_path,
                progress_log=progress_log,
                parameter_path=parameter_path,
            )
        if model_path.is_file():
            model_path.unlink()
        shutil.rmtree(candidate_dir, ignore_errors=True)
        print(
            f"Iteration {iteration} complete: mode{row['selected_mode_idx']}, "
            f"f={float(row['py_frequency_thz']):.6f} THz, "
            f"Q={float(row['py_q']):.6f}, best={float(row['best_q_so_far']):.6f}",
            flush=True,
        )
        return row

    def promote_existing_baseline_if_needed(self, baseline_model: Path) -> None:
        best = best_row(self.rows)
        if str(best["source"]) != "existing_baseline":
            return
        destination = self.model_dir / "finite_quarter.mph"
        if not destination.is_file():
            shutil.copy2(baseline_model, destination)
        baseline_progress = baseline_model.parent / "comsol_progress.log"
        if baseline_progress.is_file():
            shutil.copy2(baseline_progress, self.model_dir / "comsol_progress.log")
        _write_json(
            self.config_dir / "best_parameters.json",
            {
                "iteration": int(best["iteration"]),
                "cladding": START_CLADDING,
                "symmetry_id": SYMMETRY_ID,
                "selected_mode_idx": int(best["selected_mode_idx"]),
                "frequency_thz": float(best["py_frequency_thz"]),
                "q": float(best["py_q"]),
                "model_path": str(destination),
                "source": str(baseline_model),
            },
        )


def baseline_result_paths(parameters: SharedParameters) -> tuple[Path, Path]:
    series = FINITE_CAVITY_OUTPUT_ROOT / parameters.structure_series_label(
        "finite_quarter"
    )
    case = series / "shift0.000"
    return (
        case / "10_overview" / "eigenfrequencies.csv",
        case / "00_model" / "finite_quarter.mph",
    )


def coordinator_main(args: argparse.Namespace) -> int:
    parameters = load_shared_parameters(PARAMETER_PATH)
    validate_approved_parameters(parameters)
    baseline_data = json.loads(Path(PARAMETER_PATH).read_text(encoding="utf-8"))
    output_root = (
        Path(args.output_root).resolve()
        if args.output_root is not None
        else FINITE_CAVITY_OUTPUT_ROOT / optimization_series_name(parameters)
    )
    coordinator = OptimizationCoordinator(
        parameters,
        baseline_data,
        output_root,
        max_new_candidates=args.max_new_candidates,
    )
    coordinator.prepare_output()
    baseline_eigen_path, baseline_model_path = baseline_result_paths(parameters)
    baseline_modes = read_eigenfrequency_csv(baseline_eigen_path)
    baseline = coordinator.add_baseline(baseline_modes)

    print("Approved finite-quarter Q optimization preflight:", flush=True)
    print(f"  output: {output_root}", flush=True)
    print("  cavity: 245.000-0.960-1.155; layers: 20-20; mesh9", flush=True)
    print("  cladding start: 243.700-0.980-0.928; eta fixed", flush=True)
    print("  shift=0.000 uniform; symmetry ID=1; exactly two modes", flush=True)
    print(f"  baseline selected mode{baseline['selected_mode_idx']}: Q={baseline['py_q']}", flush=True)
    print("  first new candidates: zeta=0.924 and zeta=0.932", flush=True)
    if args.dry_run:
        return 0

    try:
        zeta_best = directional_axis_search(
            start_value=START_CLADDING["zeta"],
            coarse_step=ZETA_COARSE_STEP,
            fine_step=ZETA_FINE_STEP,
            bounds=ZETA_BOUNDS,
            stage_prefix="zeta",
            evaluate_value=lambda value, stage: coordinator.evaluate(
                Candidate(START_CLADDING["b0_nm"], FIXED_ETA, value),
                stage,
            ),
        )
        best_zeta = float(zeta_best["zeta"])

        b0_best = directional_axis_search(
            start_value=START_CLADDING["b0_nm"],
            coarse_step=B0_COARSE_STEP_NM,
            fine_step=B0_FINE_STEP_NM,
            bounds=B0_BOUNDS_NM,
            stage_prefix="b0",
            evaluate_value=lambda value, stage: coordinator.evaluate(
                Candidate(value, FIXED_ETA, best_zeta),
                stage,
            ),
        )
        best_b0 = float(b0_best["b0_nm"])

        final_rows = [b0_best]
        for delta in (-ZETA_FINE_STEP, ZETA_FINE_STEP):
            final_zeta = round(best_zeta + delta, 9)
            if ZETA_BOUNDS[0] <= final_zeta <= ZETA_BOUNDS[1]:
                final_rows.append(
                    coordinator.evaluate(
                        Candidate(best_b0, FIXED_ETA, final_zeta),
                        "zeta_joint_recheck",
                    )
                )
        final_best = best_row(final_rows)
        print(
            f"Joint local best: b0={float(final_best['b0_nm']):.3f}, "
            f"eta={float(final_best['eta']):.3f}, "
            f"zeta={float(final_best['zeta']):.3f}, "
            f"Q={float(final_best['py_q']):.6f}",
            flush=True,
        )
    except SearchBudgetExhausted as exc:
        print(f"Search budget reached: {exc}", flush=True)

    plot_q_history(coordinator.rows, coordinator.plot_path)
    coordinator.promote_existing_baseline_if_needed(baseline_model_path)
    best = best_row(coordinator.rows)
    print(
        f"Optimization best after {len(coordinator.rows) - 1} new iterations: "
        f"b0={float(best['b0_nm']):.3f}, eta={float(best['eta']):.3f}, "
        f"zeta={float(best['zeta']):.3f}, Q={float(best['py_q']):.6f}",
        flush=True,
    )
    print(f"CSV: {coordinator.history_path}", flush=True)
    print(f"Plot: {coordinator.plot_path}", flush=True)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-new-candidates", type=int, default=DEFAULT_MAX_NEW_CANDIDATES)
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--parameter-file", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--work-dir", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--result-path", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--progress-log", type=Path, default=None, help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.worker:
        required = (
            args.parameter_file,
            args.work_dir,
            args.result_path,
            args.progress_log,
        )
        if any(value is None for value in required):
            raise ValueError("Worker mode requires all worker paths")
        return worker_main(args)
    if args.max_new_candidates <= 0:
        raise ValueError("--max-new-candidates must be positive")
    return coordinator_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
