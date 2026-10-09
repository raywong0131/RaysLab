#!/usr/bin/env python3
"""Transactionally redraw existing finite-cavity Hz and Wem field images."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import shutil

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from comsol_workflow.field_plotting import save_field_plot
from comsol_workflow.simulation_utils import (
    SimulationRun,
    reconstruct_quarter_scalar_field,
)


DEFAULT_ROOT = Path(__file__).resolve().parents[1] / ".out" / "finite_cavity"
STAGING_DIRNAME = ".field_plot_update_staging"


@dataclass(frozen=True)
class ModeExport:
    export_dir: Path
    mode_dir: Path
    case_dir: Path
    model_path: Path
    mode_idx: int
    frequency_thz: float
    quality_factor: float
    is_quarter: bool

    @property
    def hz_path(self) -> Path:
        return self.export_dir / "Hz_center.parquet"


def discover_mode_exports(root: Path) -> list[ModeExport]:
    root = Path(root).resolve()
    exports: list[ModeExport] = []
    for hz_path in sorted(root.rglob("Hz_center.parquet")):
        export_dir = hz_path.parent
        if export_dir.name != "11_simulation_exports":
            continue
        mode_dir = export_dir.parent
        case_dir = export_dir.parents[2]
        eigen_path = mode_dir / "10_overview" / "eigenfrequency.csv"
        model_path = case_dir / "00_model" / "finite_quarter.mph"
        if not eigen_path.is_file():
            raise FileNotFoundError(f"Missing per-mode eigenfrequency table: {eigen_path}")
        eigen_df = pd.read_csv(eigen_path)
        if len(eigen_df) != 1:
            raise ValueError(f"Expected one eigenfrequency row in {eigen_path}")
        row = eigen_df.iloc[0]
        mode_idx = int(row["mode_idx"])
        is_quarter = all(
            column in eigen_df.columns
            for column in ("symmetry_id", "x_boundary", "y_boundary")
        )
        if not model_path.is_file():
            fallback_models = sorted((case_dir / "00_model").glob("*.mph"))
            if len(fallback_models) != 1:
                raise FileNotFoundError(
                    f"Expected one saved MPH model under {case_dir / '00_model'}"
                )
            model_path = fallback_models[0]
        required_images = (
            export_dir / "Hz_Re_2d.png",
            export_dir / "Hz_Im_2d.png",
            export_dir / "Wem_2d.png",
        )
        missing_images = [path for path in required_images if not path.is_file()]
        if missing_images:
            raise FileNotFoundError(f"Missing field images: {missing_images}")
        exports.append(
            ModeExport(
                export_dir=export_dir,
                mode_dir=mode_dir,
                case_dir=case_dir,
                model_path=model_path,
                mode_idx=mode_idx,
                frequency_thz=float(row["re"]),
                quality_factor=float(row["q"]),
                is_quarter=is_quarter,
            )
        )
    if not exports:
        raise FileNotFoundError(f"No finite-cavity mode exports found under {root}")
    return exports


def _stage_path(staging_dir: Path, export_index: int, filename: str) -> Path:
    return staging_dir / "new" / f"mode_{export_index:04d}" / filename


class _SimulationView:
    """Borrow field extraction without owning or clearing the loaded model."""

    get_2d_fields = SimulationRun.get_2d_fields

    def __init__(self, model):
        self.model = model
        self.plane_datasets = {"center": "cpl1"}


def stage_hz_images(
    exports: list[ModeExport],
    staging_dir: Path,
    *,
    dpi: int,
) -> list[tuple[Path, Path]]:
    replacements: list[tuple[Path, Path]] = []
    for export_index, mode in enumerate(exports):
        field = pd.read_parquet(mode.hz_path)
        coordinates = field[["x", "y"]].to_numpy(dtype=float)
        for filename, column, quantity_label in (
            ("Hz_Re_2d.png", "re", "Re(Hz)"),
            ("Hz_Im_2d.png", "im", "Im(Hz)"),
        ):
            staged_path = _stage_path(staging_dir, export_index, filename)
            save_field_plot(
                staged_path,
                coordinates,
                field[column].to_numpy(dtype=float),
                frequency_thz=mode.frequency_thz,
                quality_factor=mode.quality_factor,
                quantity_label=quantity_label,
                unit_label="A/m",
                color_scale_mode="linearsymmetric",
                dpi=dpi,
            )
            replacements.append((staged_path, mode.export_dir / filename))
        print(
            f"Staged Hz {export_index + 1}/{len(exports)}: {mode.mode_dir}",
            flush=True,
        )
    return replacements


def stage_wem_images(
    client,
    exports: list[ModeExport],
    staging_dir: Path,
    *,
    dpi: int,
) -> list[tuple[Path, Path]]:
    replacements: list[tuple[Path, Path]] = []
    grouped: dict[Path, list[tuple[int, ModeExport]]] = {}
    for export_index, mode in enumerate(exports):
        grouped.setdefault(mode.model_path, []).append((export_index, mode))

    for model_number, (model_path, modes) in enumerate(grouped.items(), start=1):
        print(
            f"Loading saved model {model_number}/{len(grouped)}: {model_path}",
            flush=True,
        )
        model = client.load(str(model_path))
        try:
            sim_run = _SimulationView(model)
            for export_index, mode in modes:
                coordinates, values = sim_run.get_2d_fields(
                    mode.mode_idx,
                    "ewfd.Wav",
                    "center",
                )
                coordinates = np.asarray(coordinates, dtype=float)[:, :2]
                values = np.real(np.asarray(values))
                if mode.is_quarter:
                    coordinates, values = reconstruct_quarter_scalar_field(
                        coordinates,
                        values,
                        x_reflection_sign=1.0,
                        y_reflection_sign=1.0,
                    )
                staged_path = _stage_path(staging_dir, export_index, "Wem_2d.png")
                save_field_plot(
                    staged_path,
                    coordinates,
                    values,
                    frequency_thz=mode.frequency_thz,
                    quality_factor=mode.quality_factor,
                    quantity_label="Wem",
                    unit_label="J/m³",
                    color_scale_mode="linear",
                    dpi=dpi,
                )
                replacements.append((staged_path, mode.export_dir / "Wem_2d.png"))
                print(f"  Staged Wem: {mode.mode_dir.name}", flush=True)
        finally:
            client.remove(model)
    return replacements


def validate_staged_images(
    replacements: list[tuple[Path, Path]],
    *,
    expected_count: int,
) -> None:
    if len(replacements) != expected_count:
        raise ValueError(
            f"Expected {expected_count} staged images, got {len(replacements)}"
        )
    staged_paths = [staged for staged, _destination in replacements]
    destination_paths = [destination for _staged, destination in replacements]
    if len(set(staged_paths)) != expected_count:
        raise ValueError("Duplicate staged image path detected")
    if len(set(destination_paths)) != expected_count:
        raise ValueError("Duplicate destination image path detected")
    invalid = [path for path in staged_paths if not path.is_file() or path.stat().st_size == 0]
    if invalid:
        raise ValueError(f"Missing or empty staged images: {invalid}")


def replace_images_transactionally(
    replacements: list[tuple[Path, Path]],
    staging_dir: Path,
) -> None:
    backup_dir = staging_dir / "backup"
    backup_dir.mkdir(parents=True, exist_ok=False)
    backups: list[tuple[Path, Path]] = []
    for index, (_staged, destination) in enumerate(replacements):
        backup = backup_dir / f"{index:04d}_{destination.name}"
        shutil.copy2(destination, backup)
        backups.append((backup, destination))

    try:
        for staged, destination in replacements:
            os.replace(staged, destination)
    except Exception:
        for backup, destination in backups:
            shutil.copy2(backup, destination)
        raise


def update_existing_field_plots(root: Path, *, dpi: int) -> int:
    root = Path(root).resolve()
    staging_dir = root / STAGING_DIRNAME
    if staging_dir.exists():
        raise FileExistsError(
            f"Staging directory already exists; inspect or remove it first: {staging_dir}"
        )
    exports = discover_mode_exports(root)
    print(
        f"Discovered {len(exports)} modes and "
        f"{len({mode.model_path for mode in exports})} saved models.",
        flush=True,
    )
    staging_dir.mkdir(parents=True)
    try:
        replacements = stage_hz_images(exports, staging_dir, dpi=dpi)
        import mph

        client = mph.start(cores=1)
        try:
            replacements.extend(
                stage_wem_images(client, exports, staging_dir, dpi=dpi)
            )
        finally:
            client.disconnect()
        validate_staged_images(replacements, expected_count=3 * len(exports))
        replace_images_transactionally(replacements, staging_dir)
        print(f"Updated {len(replacements)} field images.", flush=True)
        return len(replacements)
    except Exception:
        print(f"Update failed; staging retained at {staging_dir}", flush=True)
        raise
    finally:
        # A successful replacement leaves only backups and empty new directories.
        if (staging_dir / "backup").is_dir() and not any(
            staged.is_file() for staged in (staging_dir / "new").rglob("*.png")
        ):
            shutil.rmtree(staging_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="perform the update; without this flag only inventory is printed",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    exports = discover_mode_exports(args.root)
    model_count = len({mode.model_path for mode in exports})
    print(
        f"finite-cavity inventory: {len(exports)} modes, "
        f"{3 * len(exports)} images, {model_count} saved models"
    )
    if not args.apply:
        print("Dry run only; pass --apply to update the existing images.")
        return
    update_existing_field_plots(args.root, dpi=args.dpi)


if __name__ == "__main__":
    main()
