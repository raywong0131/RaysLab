"""Re-extract S4 sections from a saved single-case MPH without an eigensolve."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from types import SimpleNamespace
import uuid

import numpy as np
import pandas as pd

from comsol_workflow.boundary_plotting import plot_s4_case, _load_area_field
from scripts.run_main import run_boundary_analysis as s4


def reprocess(case_dir, replace_results=False, area_field_only=False):
    case_dir = Path(case_dir).resolve()
    manifest = s4.json_read(case_dir.parent / "99_config/s4_manifest.json")
    if len(manifest["cases"]) != 1 or manifest["cases"][0]["case_id"] != case_dir.name:
        raise ValueError("Reprocessing requires a recorded single-case S4 run")
    case = manifest["cases"][0]
    if s4.digest(case_dir.parent / case["snapshot"]) != case["snapshot_sha256"]:
        raise ValueError("Saved parameter identity changed")
    source = case_dir / "00_model/s4_gamma.mph"
    geometry = s4.json_read(case_dir / "99_config/realized_geometry.json")
    old_summary = s4.json_read(case_dir / "99_config/s4_summary.json")
    task = s4.TASK_ROOT / ("s4_postprocess_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:6])
    staged = task / "case"
    data_dir = s4.data_directory(case_dir)
    output = staged / data_dir.name
    output.mkdir(parents=True)
    settings = {**manifest["s4"], "z_nm": [0.0, 50.0, 100.0], "primary_z_nm": 0.0,
                "surface_side": "slab_layer_domains_at_exact_z", "normalization": "cell_Hz_RMS_at_z0"}
    settings.pop("z_over_H", None)
    audit = {"source_case": str(case_dir), "source_mph_sha256": s4.digest(source),
             "eigensolves": 0, "settings": settings, "source_files": s4.source_hashes(),
             "status": "starting"}
    s4.write_json(task / "audit.json", audit)
    model, client = None, None
    try:
        audit["preflight"] = s4.runtime_preflight()
        import mph
        client = mph.start(version="6.3")
        model = client.load(str(source))
        parameters = model.java.param()
        if any(abs(float(parameters.evaluate(k))) > 1e-12 for k in ("kx", "ky")):
            raise ValueError("MPH is not a Gamma solution")
        if not np.isclose(float(parameters.evaluate("a")), manifest["constants"]["A"] * 1e-6, rtol=1e-12, atol=0):
            raise ValueError("MPH lattice identity mismatch")
        height = float(parameters.evaluate("H"))
        s4.section_coordinates(settings, height)
        parameters.set("s4_expected_H", geometry["simulation_config"]["slab_height"])
        if not np.isclose(height, float(parameters.evaluate("s4_expected_H")), rtol=1e-12, atol=0):
            raise ValueError("MPH slab identity mismatch")
        parameters.remove("s4_expected_H")
        modes = pd.read_csv(data_dir / "eigenfrequencies.csv")
        # Global evaluation only; numerical.run() does not run the study/solver.
        frequencies = model.java.result().numerical("gev1")
        frequencies.run()
        actual = np.asarray(frequencies.getReal()).T
        if actual.shape != (len(modes), 3) or not np.allclose(actual[:, :2], modes[["re", "im"]].to_numpy(), rtol=1e-10, atol=1e-12):
            raise ValueError("MPH eigenfrequencies do not match the saved candidate table")
        runner = SimpleNamespace(model=model, config=SimpleNamespace(**geometry["simulation_config"]))
        sampler = s4.Sampler(runner)
        if area_field_only:
            # Reuse the exact saved quadrature and phase; do not recompute other products.
            settings = s4.json_read(case_dir / "99_config/s4_postprocess.json")
            area = pd.read_csv(data_dir / "area_comparison.csv")
            area = area[area.z_nm == 0]
            refinement = int(area.refinement.max())
            area = area[area.refinement == refinement]
            if len(area) != 1:
                raise ValueError("Expected one saved primary area integral")
            mode = int(old_summary["mode_idx"])
            if not np.isclose(modes.loc[mode, "re"], old_summary["frequency_thz"], rtol=1e-12):
                raise ValueError("Saved area mode identity mismatch")
            outer = np.asarray(geometry["outer"]) * 1e-6
            holes = [np.asarray(h) * 1e-6 for h in geometry["holes"]]
            xy, weights, _ = s4.area_rule(outer, holes, settings["area_order"], settings["area_subdivisions"] * refinement)
            factor = complex(old_summary["normalization_factor_re"], old_summary["normalization_factor_im"])
            values = factor * sampler.plane(mode, xy, 0.0, ["ewfd.Ex", "ewfd.Ey"])
            field_path = s4.export_area_field(output, sampler, mode, outer, xy, weights, values, factor, refinement)
            _load_area_field(staged, area, old_summary, refinement)
            if s4.digest(source) != audit["source_mph_sha256"]:
                raise RuntimeError("Source MPH changed; refusing to publish")
            audit.update(status="area_field_staged", area_field_only=True, quadrature_samples=len(xy),
                         selected_mode=mode, integral_matches_saved_csv=True)
            if replace_results:
                target = data_dir / field_path.name
                if target.exists():
                    backup = task / "backup" / target.relative_to(case_dir.parent)
                    backup.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(target, backup)
                shutil.copy2(field_path, target)
                audit.update(status="area_field_published", published_files=[str(target)])
            return task
        mode, chosen, _, phase = s4.select_mode(runner, sampler, modes, manifest["constants"]["A"] * 1e-6,
                                               0.0, output, settings)
        if mode != old_summary["mode_idx"]:
            raise ValueError("z=0 mode identification differs from the original target")
        omega = 2 * np.pi * complex(modes.loc[mode, "re"], -modes.loc[mode, "im"]) * 1e12
        holes = [np.asarray(h) * 1e-6 for h in geometry["holes"]]
        summary = s4.evaluate_case(runner, sampler, mode, omega, np.asarray(geometry["outer"]) * 1e-6,
                                   holes, settings, phase, output)
        result_summary = {**case, **chosen, **summary, "status": "computed_pending_review",
                          "postprocess_audit": str(task / "audit.json")}
        s4.write_json(staged / "99_config/s4_summary.json", result_summary)
        s4.write_json(staged / "99_config/s4_postprocess.json", settings)
        factor = complex(summary["normalization_factor_re"], summary["normalization_factor_im"])
        predictions = s4.frozen_predictions(sampler, mode, holes, settings, factor, height,
                                            [(case["case_id"], holes)], case["case_id"])
        s4.write_csv(task / "frozen_field_predictions.csv", predictions)
        client.remove(model)
        model = None
        if s4.digest(source) != audit["source_mph_sha256"]:
            raise RuntimeError("Source MPH changed; refusing to publish")
        plot_s4_case(staged)
        audit.update(status="staged", source_mph_unchanged=True, selected_mode=mode,
                     sections_nm=summary["sections_nm"])
        if replace_results:
            # All new products are ready before any original is replaced. Back up
            # every replacement first so a Windows file-lock failure is recoverable.
            replacements = [(p, case_dir / p.relative_to(staged)) for p in staged.rglob("*") if p.is_file()]
            replacements.append((task / "frozen_field_predictions.csv", case_dir.parent / "80_logs/frozen_field_predictions.csv"))
            root_updates = [case_dir.parent / name for name in (
                "80_logs/case_table.csv", "80_logs/zero_estimates.csv", "80_logs/acceptance.md", "99_config/execution.json")]
            targets = [target for _, target in replacements] + root_updates
            for target in targets:
                if not target.resolve().is_relative_to(case_dir.parent):
                    raise ValueError("Replacement escaped the S4 run directory")
                if target.exists():
                    backup = task / "backup" / target.relative_to(case_dir.parent)
                    backup.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(target, backup)
            audit["replaced_files"] = [str(p) for p in targets]
            s4.write_json(task / "audit.json", audit)
            for source_file, target in replacements:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_file, target)
            execution = s4.json_read(case_dir.parent / "99_config/execution.json")
            execution.update(cases=[result_summary], status="complete_pending_scientific_review",
                             postprocess_audit=str(task / "audit.json"))
            s4.write_json(case_dir.parent / "99_config/execution.json", execution)
            s4.report(case_dir.parent)
            audit["status"] = "published_pending_scientific_review"
    except Exception as exc:
        audit.update(status="failed", error=str(exc))
        raise
    finally:
        if model is not None:
            client.remove(model)
        audit["source_mph_unchanged"] = s4.digest(source) == audit["source_mph_sha256"]
        s4.write_json(task / "audit.json", audit)
        print(task, flush=True)
    return task


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir")
    parser.add_argument("--replace-results", action="store_true", help="publish staged products with task-local backups")
    parser.add_argument("--area-field-only", action="store_true", help="extract only z=0 E samples; preserve existing integrals and figures")
    args = parser.parse_args(argv)
    reprocess(args.case_dir, args.replace_results, args.area_field_only)


if __name__ == "__main__":
    main()
