#!/usr/bin/env python3
"""Regenerate unit-cell 2D polarization and winding figures without COMSOL."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.unit_cell_2d_analysis import (  # noqa: E402
    run_analysis,
    unit_cell_2d_series_name,
)
from scripts.run_main.parameter_config import (  # noqa: E402
    UNIT_CELL_2D_OUTPUT_ROOT,
    load_shared_parameters,
    parameter_path_from_environment,
)


def default_series_dir() -> Path:
    parameters = load_shared_parameters(parameter_path_from_environment())
    settings = parameters.unit_cell_2d
    name = unit_cell_2d_series_name(
        b0_nm=settings.cell.b0_nm,
        eta=settings.cell.eta,
        zeta=settings.cell.zeta,
        mesh_size=parameters.mesh_size,
        q_max_over_g=settings.q_max_over_g,
        q_points_per_axis=settings.q_points_per_axis,
        target_bands=settings.target_bands,
        is_quarter=settings.is_quarter,
    )
    return UNIT_CELL_2D_OUTPUT_ROOT / name


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Redraw unit-cell 2D c(k), polarization, and winding outputs."
    )
    parser.add_argument(
        "series_dir",
        nargs="?",
        type=Path,
        help="Existing series directory; defaults to the active parameter.json series.",
    )
    return parser


def main(argv: list[str] | None = None) -> dict[str, object]:
    args = build_parser().parse_args(argv)
    parameters = load_shared_parameters(parameter_path_from_environment())
    series_dir = args.series_dir or default_series_dir()
    summary = run_analysis(
        series_dir,
        angle_methods=parameters.unit_cell_2d.angle_methods,
        mask_relative_threshold=(
            parameters.unit_cell_2d.mask_relative_threshold
        ),
    )
    print(f"Unit-cell 2D analysis complete: {series_dir}")
    print("Winding summary:", series_dir / "99_config" / "winding_summary.json")
    return summary


if __name__ == "__main__":
    main()
