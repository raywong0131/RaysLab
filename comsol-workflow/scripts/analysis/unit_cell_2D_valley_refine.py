#!/usr/bin/env python3
"""Prepare, solve, resume, and plot COMSOL-local refinement of the py valley."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from comsol_workflow.unit_cell_2d_valley_refinement import (  # noqa: E402
    load_source_scan,
    prepare_valley_refinement,
    refined_series_name,
)
from scripts.run_main.parameter_config import (  # noqa: E402
    UNIT_CELL_2D_OUTPUT_ROOT,
    load_shared_parameters,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Use a completed uniform unit-cell scan to plan and execute "
            "real-COMSOL local refinement of the py low-value valley."
        )
    )
    parser.add_argument("source_series", type=Path)
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=REPO_ROOT / "scripts" / "parameter.json",
        help="Task-local parameter snapshot matching the source series.",
    )
    stage = parser.add_mutually_exclusive_group()
    stage.add_argument(
        "--prepare-only",
        action="store_true",
        help="Create and inspect the round-1 sampling plan without importing mph.",
    )
    stage.add_argument(
        "--plot-only",
        action="store_true",
        help="Merge completed real COMSOL points and redraw only the refined B4.",
    )
    stage.add_argument(
        "--prepare-repair-only",
        action="store_true",
        help=(
            "Audit the completed run and create a deduplicated convergence "
            "correction plan without starting COMSOL."
        ),
    )
    stage.add_argument(
        "--run-approved-repair",
        action="store_true",
        help=(
            "Verify and solve only the separately reviewed convergence-repair "
            "CSV, preserving all existing refined caches."
        ),
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Reject an existing point cache instead of resuming it.",
    )
    parser.add_argument(
        "--allow-concurrent-comsol",
        action="store_true",
        help="Allow solving while another project COMSOL process is detected.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    parameters = load_shared_parameters(args.parameter_file)

    if args.plot_only or args.prepare_repair_only or args.run_approved_repair:
        # Plot-only must preserve the approved, dynamically extended sampling
        # plan.  Calling prepare_valley_refinement() here would replace it with
        # the round-1 plan before the completed 709-point table is read.
        source = load_source_scan(args.source_series)
        series_dir = (
            Path(UNIT_CELL_2D_OUTPUT_ROOT).resolve()
            / refined_series_name(parameters, source)
        )
        if args.prepare_repair_only:
            from scripts.run_sweep.unit_cell_2d_valley_solver import (  # noqa: PLC0415
                prepare_convergence_repair,
            )

            repair = prepare_convergence_repair(series_dir, parameters)
            print(
                "Convergence repair prepared: "
                f"new={repair['proposed_new_point_count']}, "
                f"projected={repair['projected_refined_point_count']}, "
                f"figure={repair['sampling_plan_figure']}"
            )
        elif args.run_approved_repair:
            from scripts.run_sweep.unit_cell_2d_valley_solver import (  # noqa: PLC0415
                plot_completed_refinement,
                solve_approved_convergence_repair,
            )

            repair = solve_approved_convergence_repair(
                series_dir,
                parameters,
                resume=not args.no_resume,
                allow_concurrent_comsol=args.allow_concurrent_comsol,
            )
            output = plot_completed_refinement(series_dir, parameters)
            print(
                "Approved convergence repair finished: "
                f"status={repair['status']}, "
                f"solves={repair['solve_count']}, figure={output}"
            )
        else:
            from scripts.run_sweep.unit_cell_2d_valley_solver import (  # noqa: PLC0415
                plot_completed_refinement,
            )

            output = plot_completed_refinement(series_dir, parameters)
            print(f"Refined B4: {output}")
        return 0

    summary = prepare_valley_refinement(
        args.source_series,
        parameters,
        UNIT_CELL_2D_OUTPUT_ROOT,
    )
    print(f"Valley-refinement output: {summary['output_dir']}")
    print(
        "Prepared: "
        f"candidates={summary['candidate_count']}, "
        f"segments={summary['segment_count']}, "
        f"junctions={summary['junction_count']}, "
        f"round-1 new points={summary['round_1_unique_new_point_count']}"
    )
    print(f"Sampling plan: {summary['sampling_plan_figure']}")
    if args.prepare_only:
        print("COMSOL was not started (--prepare-only).")
        return 0

    # The expensive stage is imported only after preparation.  This invariant
    # is tested so prepare-only can be reviewed safely on a busy workstation.
    from scripts.run_sweep.unit_cell_2d_valley_solver import (  # noqa: PLC0415
        plot_completed_refinement,
        solve_valley_refinement,
    )

    solve_valley_refinement(
        Path(summary["output_dir"]),
        parameters,
        resume=not args.no_resume,
        allow_concurrent_comsol=args.allow_concurrent_comsol,
    )
    plot_completed_refinement(Path(summary["output_dir"]), parameters)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
