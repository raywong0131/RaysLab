"""Small public routers; each scientific entrypoint retains its own parser."""
import argparse
import json
import os
from runpy import run_module
import sys

from comsol_workflow.dipolar_inputs import INPUT_MANIFEST_ENV, load_dipolar_inputs


DIPOLAR_ANALYSES = {
    "singularity": "analyze_dipolar_singularity", "complete": "analyze_dipolar_complete",
    "broadening": "visualize_dipolar_broadening", "broadening-2d": "visualize_dipolar_broadening_2d",
    "bloch-validation": "validate_dipolar_bloch_radiation", "gamma-origin": "analyze_dipolar_gamma_origin",
    "cladding-scattering": "analyze_dipolar_cladding_scattering", "angle-filter": "plot_dipolar_angle_filter",
    "cell-radiation": "analyze_dipolar_cell_radiation", "sector-interference": "analyze_dipolar_sector_interference",
    "spatial-filter": "analyze_dipolar_spatial_filter", "unitcell-synthesis": "analyze_dipolar_unitcell_synthesis",
    "volume-radiation": "plot_dipolar_volume_radiation", "gamma-mechanism": "plot_dipolar_gamma_mechanism",
    "cell-quadrature": "audit_dipolar_cell_quadrature",
}
BOUNDARY_SCANS = {
    "surface-zeta": "run_boundary_surface_zeta_scan", "mesh-convergence": "run_boundary_mesh_convergence",
    "dual-validation": "run_boundary_dual_validation", "full-refined": "run_boundary_full_refined_scan",
    "peak-resampling": "run_boundary_peak_resampling", "px-mesh": "run_boundary_px_mesh_test",
    "signed-surface": "run_boundary_surface_signed_scan",
}


def _dispatch(package, entries, argv=None, *, dipolar=False):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    if isinstance(entries, dict):
        parser.add_argument("stage", choices=entries)
    if dipolar:
        parser.add_argument("--source-manifest", help="schema v1 source JSON; paths relative to that file")
        parser.add_argument("--describe-inputs", action="store_true", help="show sources without starting COMSOL or writing outputs")
    args, remaining = parser.parse_known_args(argv)
    saved_env = {key: os.environ.get(key) for key in (INPUT_MANIFEST_ENV, "COMSOL_WORKFLOW_OUTPUT_LAYOUT")}
    saved_argv = sys.argv
    try:
        if dipolar and args.source_manifest:
            os.environ[INPUT_MANIFEST_ENV] = args.source_manifest
        os.environ["COMSOL_WORKFLOW_OUTPUT_LAYOUT"] = "unified"
        if dipolar:
            inputs = load_dipolar_inputs()
            if args.describe_inputs:
                print(json.dumps(inputs.to_metadata(), indent=2, ensure_ascii=False))
                return
        module = package + "." + (entries[args.stage] if isinstance(entries, dict) else entries)
        sys.argv = [module, *(remaining[1:] if remaining[:1] == ["--"] else remaining)]
        run_module(module, run_name="__main__", alter_sys=True)
    finally:
        sys.argv = saved_argv
        for key, value in saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def dipolar_analysis_main(argv=None):
    return _dispatch("scripts.analysis", DIPOLAR_ANALYSES, argv, dipolar=True)


def boundary_scan_main(argv=None):
    return _dispatch("scripts.run_sweep", BOUNDARY_SCANS, argv)


def dipolar_reference_main(argv=None):
    return _dispatch("scripts.run_main", "run_dipolar_periodic_reference", argv, dipolar=True)


def dipolar_volume_main(argv=None):
    return _dispatch("scripts.run_main", "run_dipolar_volume_export", argv, dipolar=True)


def dipolar_replay_main(argv=None):
    return _dispatch("scripts.run_main", "run_dipolar_boundary_replay", argv, dipolar=True)
