#!/usr/bin/env python3
"""Export one existing finite-cavity mode as fully vector SVG field plots."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-comsol-workflow")

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from comsol_workflow.field_plotting import (
    save_field_plot_hybrid_svg,
    save_field_plot_svg,
)
from comsol_workflow.simulation_utils import reconstruct_quarter_scalar_field
from scripts.analysis.update_finite_field_plots import (
    _SimulationView,
    discover_mode_exports,
)


STAGING_DIRNAME = ".field_svg_staging"


def _raster_image_count(path: Path) -> int:
    marker = b"<image"
    overlap = b""
    count = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            data = overlap + chunk
            count += data.count(marker)
            overlap = data[-(len(marker) - 1) :]
    return count


def _validate_svg(path: Path, *, expected_raster_images: int) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty SVG: {path}")
    with path.open("rb") as stream:
        header = stream.read(4096)
    if b"<svg" not in header:
        raise ValueError(f"Not an SVG document: {path}")
    raster_images = _raster_image_count(path)
    if raster_images != expected_raster_images:
        raise ValueError(
            f"Expected {expected_raster_images} raster image nodes in {path}, "
            f"found {raster_images}"
        )


def export_mode_svg(
    path: Path,
    *,
    hybrid: bool = False,
    dpi: int = 600,
) -> list[Path]:
    modes = discover_mode_exports(path)
    if len(modes) != 1:
        raise ValueError(f"Expected exactly one mode, found {len(modes)}")
    mode = modes[0]
    staging_dir = mode.export_dir / STAGING_DIRNAME
    if staging_dir.exists():
        raise FileExistsError(f"SVG staging already exists: {staging_dir}")
    staging_dir.mkdir()
    staged_paths: list[Path] = []
    suffix = "_hybrid.svg" if hybrid else ".svg"
    save_svg = save_field_plot_hybrid_svg if hybrid else save_field_plot_svg
    expected_raster_images = 1 if hybrid else 0
    try:
        field = pd.read_parquet(mode.hz_path)
        coordinates = field[["x", "y"]].to_numpy(dtype=float)
        for filename, column, quantity_label in (
            (f"Hz_Re_2d{suffix}", "re", "Re(Hz)"),
            (f"Hz_Im_2d{suffix}", "im", "Im(Hz)"),
        ):
            staged_path = staging_dir / filename
            kind = "hybrid" if hybrid else "full-vector"
            print(f"Rendering {kind} {filename} ...", flush=True)
            save_kwargs = {
                "frequency_thz": mode.frequency_thz,
                "quality_factor": mode.quality_factor,
                "quantity_label": quantity_label,
                "unit_label": "A/m",
                "color_scale_mode": "linearsymmetric",
            }
            if hybrid:
                save_kwargs["dpi"] = dpi
            save_svg(
                staged_path,
                coordinates,
                field[column].to_numpy(dtype=float),
                **save_kwargs,
            )
            _validate_svg(
                staged_path,
                expected_raster_images=expected_raster_images,
            )
            print(
                f"Finished {filename}: {staged_path.stat().st_size / 1024**2:.1f} MiB",
                flush=True,
            )
            staged_paths.append(staged_path)

        import mph

        client = mph.start(cores=1)
        try:
            model = client.load(str(mode.model_path))
            try:
                sim_run = _SimulationView(model)
                wem_coordinates, wem = sim_run.get_2d_fields(
                    mode.mode_idx,
                    "ewfd.Wav",
                    "center",
                )
                wem_coordinates = np.asarray(wem_coordinates, dtype=float)[:, :2]
                wem = np.real(np.asarray(wem))
                if mode.is_quarter:
                    wem_coordinates, wem = reconstruct_quarter_scalar_field(
                        wem_coordinates,
                        wem,
                        x_reflection_sign=1.0,
                        y_reflection_sign=1.0,
                    )
                filename = f"Wem_2d{suffix}"
                staged_path = staging_dir / filename
                kind = "hybrid" if hybrid else "full-vector"
                print(f"Rendering {kind} {filename} ...", flush=True)
                save_kwargs = {
                    "frequency_thz": mode.frequency_thz,
                    "quality_factor": mode.quality_factor,
                    "quantity_label": "Wem",
                    "unit_label": "J/m³",
                    "color_scale_mode": "linear",
                }
                if hybrid:
                    save_kwargs["dpi"] = dpi
                save_svg(
                    staged_path,
                    wem_coordinates,
                    wem,
                    **save_kwargs,
                )
                _validate_svg(
                    staged_path,
                    expected_raster_images=expected_raster_images,
                )
                print(
                    f"Finished {filename}: "
                    f"{staged_path.stat().st_size / 1024**2:.1f} MiB",
                    flush=True,
                )
                staged_paths.append(staged_path)
            finally:
                client.remove(model)
        finally:
            client.disconnect()

        outputs = []
        for staged_path in staged_paths:
            output = mode.export_dir / staged_path.name
            os.replace(staged_path, output)
            outputs.append(output)
        return outputs
    except Exception:
        print(f"SVG export failed; staging retained at {staging_dir}", flush=True)
        raise
    finally:
        if staging_dir.is_dir() and not any(staging_dir.iterdir()):
            shutil.rmtree(staging_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="one finite-cavity mode directory")
    parser.add_argument(
        "--hybrid",
        action="store_true",
        help="rasterize only the field at --dpi while retaining vector annotations",
    )
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for output in export_mode_svg(args.path, hybrid=args.hybrid, dpi=args.dpi):
        print(output)


if __name__ == "__main__":
    main()
