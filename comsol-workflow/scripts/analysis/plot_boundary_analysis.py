"""Render S4 boundary contributions from one saved case; no COMSOL required."""
import argparse

from comsol_workflow.boundary_plotting import plot_s4_case


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", help="S4 case with 80_logs (or legacy 01_results) and 99_config")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--edges-only", action="store_true", help="Redraw edge figures; compact cases retain only the py main figure")
    selection.add_argument("--validation-only", action="store_true", help="Render only the py S31 diagnostic figure")
    args = parser.parse_args(argv)
    for path in plot_s4_case(args.case_dir, edges_only=args.edges_only, validation_only=args.validation_only):
        print(path)


if __name__ == "__main__":
    main()
