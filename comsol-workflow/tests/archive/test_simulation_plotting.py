from tests._paths import output_root
import json
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis
from comsol_workflow.geometry_utils import create_hexagon_design, visualize_hexagon_design
from comsol_workflow.simulation_spatial.hexagon_unit_cell import (
    SimulationConfig as HexagonSimulationConfig,
)
from comsol_workflow.simulation_utils import (
    BoundarySpec,
    LayerSpec,
    SimulationConfig as UnifiedSimulationConfig,
)
from comsol_workflow.simulation_utils import SimulationRun as UnifiedSimulationRun


def log(message):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {message}", flush=True)


def unified_config_from_spatial(spatial_config):
    return UnifiedSimulationConfig(
        wavelength=spatial_config.wavelength,
        air_height=spatial_config.air_height,
        air_layer_top_factor=spatial_config.air_layer_top_factor,
        pml_layer_top_factor=spatial_config.pml_layer_top_factor,
        air_refractive_index=spatial_config.air_refractive_index,
        air_extinction_coefficient=spatial_config.air_extinction_coefficient,
        eigenfrequency_shift=spatial_config.eigenfrequency_shift,
        eigenmode_count=spatial_config.eigenmode_count,
        air_cutplane_z=spatial_config.air_cutplane_z,
        selection_tolerance=spatial_config.selection_tolerance,
        mode_type=spatial_config.mode_type,
    )


def layer_from_spatial(spatial_config, holes):
    return LayerSpec(
        height=spatial_config.slab_height,
        holes=[np.asarray(hole, dtype=float) for hole in holes],
        refractive_index=spatial_config.slab_refractive_index,
        extinction_coefficient=spatial_config.slab_extinction_coefficient,
        label="slab",
    )


def build_hexagon_case(out_root):
    a = 0.82
    r_0 = 0.82 / 3
    b_0 = 0.23
    d = 0.02

    symmetry_config = {
        "r": 5,
        "theta": None,
        "b": 3,
        "phi": None,
    }

    r_f0 = 0.97
    r_fs = 0.0
    theta_fs = 0.0
    b_square_f0 = 1.0
    b_square_fs = 0.1
    phi_fs = 0.0

    fourier_min = get_fourier_basis_min_nonzero_values(6)

    r_f = np.zeros(6)
    theta_f = np.zeros(6)
    b_square_f = np.zeros(6)
    phi_f = np.zeros(6)

    r_f[0] = 1
    b_square_f[0] = 1

    if symmetry_config["r"] is not None:
        r_f[int(symmetry_config["r"])] = r_fs
    r_f *= r_f0

    if symmetry_config["theta"] is not None:
        theta_f[int(symmetry_config["theta"])] = theta_fs

    if symmetry_config["b"] is not None:
        b_square_f[int(symmetry_config["b"])] = b_square_fs
    b_square_f *= b_square_f0

    if symmetry_config["phi"] is not None:
        phi_f[int(symmetry_config["phi"])] = phi_fs

    r_f = r_0 * r_f
    b_square_f = b_0 * b_0 * b_square_f

    r = fourier_to_standard_basis(r_f / fourier_min, 6)
    theta = fourier_to_standard_basis(theta_f / fourier_min, 6)
    b_square = fourier_to_standard_basis(b_square_f / fourier_min, 6)
    b = np.sqrt(np.maximum(b_square, 0.0))
    phi = fourier_to_standard_basis(phi_f / fourier_min, 6)

    hole_params = np.stack([r, theta, b, phi], axis=0)
    hexagon, holes, info = create_hexagon_design(a, hole_params.T)
    triangles = [h.tolist() for h in holes]

    case_root = os.path.join(out_root, "hexagon_unit_cell")
    os.makedirs(case_root, exist_ok=True)
    with open(os.path.join(case_root, "config.json"), "w") as f:
        json.dump(
            {
                "a": a,
                "r_0": r_0,
                "b_0": b_0,
                "d": d,
                "symmetry_config": symmetry_config,
                "params": {
                    "r_f0": r_f0,
                    "r_fs": r_fs,
                    "theta_fs": theta_fs,
                    "b_square_f0": b_square_f0,
                    "b_square_fs": b_square_fs,
                    "phi_fs": phi_fs,
                },
                "hole_params": hole_params.tolist(),
                "hexagon": hexagon.tolist(),
                "triangles": triangles,
                "info": info,
            },
            f,
            indent=2,
        )

    visualize_hexagon_design(
        hexagon,
        holes,
        info,
        filename=os.path.join(case_root, "design.png"),
        annotate=False,
    )
    visualize_hexagon_design(
        hexagon,
        holes,
        info,
        filename=os.path.join(case_root, "design_annotated.png"),
        annotate=True,
    )

    spatial_config = HexagonSimulationConfig(eigenmode_count=4, mode_type="TE")
    return {
        "a": a,
        "root": case_root,
        "config": unified_config_from_spatial(spatial_config),
        "layers": [layer_from_spatial(spatial_config, triangles)],
    }


def write_status(status_path, rows, **row):
    rows.append(row)
    pd.DataFrame(rows).to_csv(status_path, index=False)


def export_plot(sim_run, status_path, rows, mode_idx, dimension, expr, expr_name, save_path):
    png_name = (
        f"{mode_idx:02d}_{expr_name}_3d.png"
        if dimension == "3d"
        else f"{mode_idx:02d}_{expr_name}_center_2d.png"
    )
    png_path = os.path.join(save_path, png_name)
    start = time.monotonic()
    log(f"export start mode={mode_idx} dim={dimension} expr={expr} path={png_path}")
    write_status(
        status_path,
        rows,
        mode_idx=mode_idx,
        dimension=dimension,
        expr=expr,
        expr_name=expr_name,
        png_path=png_path,
        status="started",
        elapsed_s=0.0,
        message="",
    )
    try:
        if dimension == "3d":
            sim_run.export_3d_fields(mode_idx, expr, expr_name, save_path)
        else:
            sim_run.export_2d_fields(
                mode_idx,
                expr,
                expr_name,
                save_path,
                plane="center",
                export_data=False,
                export_image=True,
            )
    except Exception as exc:
        write_status(
            status_path,
            rows,
            mode_idx=mode_idx,
            dimension=dimension,
            expr=expr,
            expr_name=expr_name,
            png_path=png_path,
            status="error",
            elapsed_s=time.monotonic() - start,
            message=str(exc),
        )
        raise
    write_status(
        status_path,
        rows,
        mode_idx=mode_idx,
        dimension=dimension,
        expr=expr,
        expr_name=expr_name,
        png_path=png_path,
        status="ok",
        elapsed_s=time.monotonic() - start,
        message="",
    )
    log(f"export done mode={mode_idx} dim={dimension} expr={expr}")


def run_plotting_case(case):
    save_path = os.path.join(case["root"], "simulation_utils")
    plot_path = os.path.join(save_path, "plots")
    os.makedirs(plot_path, exist_ok=True)
    status_path = os.path.join(save_path, "plot_export_status.csv")
    status_rows = []

    log("simulation start")
    simulation_start = time.monotonic()
    with UnifiedSimulationRun(config=case["config"]) as sim_run:
        sim_run.build_and_run(
            BoundarySpec("hexagon", a=case["a"]),
            case["layers"],
            simulation_mode="unit_cell",
        )
        log(f"simulation build_and_run done elapsed={time.monotonic() - simulation_start:.2f}s")

        sim_run.model.save(os.path.join(save_path, "design.mph"))
        log(f"model save done elapsed={time.monotonic() - simulation_start:.2f}s")

        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])
        df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
        log(f"eigenfrequencies save done elapsed={time.monotonic() - simulation_start:.2f}s")

        log("plot export start")
        for mode_idx in range(len(df)):
            export_plot(sim_run, status_path, status_rows, mode_idx, "3d", "ewfd.normE", "normE", plot_path)
            export_plot(sim_run, status_path, status_rows, mode_idx, "3d", "ewfd.normH", "normH", plot_path)
            export_plot(sim_run, status_path, status_rows, mode_idx, "2d", "ewfd.Hz", "ReHz", plot_path)
            export_plot(sim_run, status_path, status_rows, mode_idx, "2d", "ewfd.Hz*(-i)", "ImHz", plot_path)
            export_plot(sim_run, status_path, status_rows, mode_idx, "2d", "ewfd.normH", "normH", plot_path)
        log("plot export done")

    return df


if __name__ == "__main__":
    out_root = str(output_root() / 'test_simulation_plotting')
    if os.path.exists(out_root):
        shutil.rmtree(out_root)
    os.makedirs(out_root, exist_ok=True)

    case = build_hexagon_case(out_root)
    result = run_plotting_case(case)
    assert len(result) > 0, "No eigenmodes were computed"
