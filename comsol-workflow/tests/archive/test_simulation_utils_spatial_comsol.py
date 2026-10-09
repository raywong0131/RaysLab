from tests._paths import output_root
import json
import os
import shutil
import sys
import traceback

import numpy as np
import pandas as pd

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis
from comsol_workflow.energy_recovery import integrate_field, interpolate_field
from comsol_workflow.geometry_utils import create_hexagon_design, visualize_hexagon_design
from comsol_workflow.simulation_spatial.hexagon_unit_cell import (
    SimulationConfig as HexagonSimulationConfig,
)
from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationRun as HexagonSimulationRun
from comsol_workflow.simulation_spatial.rectangle_finite_size import (
    SimulationConfig as RectangleSimulationConfig,
)
from comsol_workflow.simulation_spatial.rectangle_finite_size import SimulationRun as RectangleSimulationRun
from comsol_workflow.simulation_spatial.square_unit_cell import (
    SimulationConfig as SquareSimulationConfig,
)
from comsol_workflow.simulation_spatial.square_unit_cell import SimulationRun as SquareSimulationRun
from comsol_workflow.simulation_utils import (
    BoundarySpec,
    LayerSpec,
    SimulationConfig as UnifiedSimulationConfig,
)
from comsol_workflow.simulation_utils import SimulationRun as UnifiedSimulationRun


def field_expr_for_mode(mode_type):
    mode_type = mode_type.upper()
    if mode_type == "TE":
        return "ewfd.Hz", "Hz"
    if mode_type == "TM":
        return "ewfd.Ez", "Ez"
    raise ValueError("mode_type must be 'TE' or 'TM'")


def circle_polygon(x0, y0, radius, vertex_count=16):
    angles = np.linspace(0, 2 * np.pi, vertex_count, endpoint=False)
    return np.column_stack([
        x0 + radius * np.cos(angles),
        y0 + radius * np.sin(angles),
    ])


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


def two_layers_from_spatial(spatial_config, holes):
    # Unified geometry models layer 0 with height/2. These input heights
    # give modeled layer thicknesses H/6 and H/3, which sum to H/2.
    layer_heights = [
        f"({spatial_config.slab_height})/3",
        f"({spatial_config.slab_height})/3",
    ]
    return [
        LayerSpec(
            height=height,
            holes=[np.asarray(hole, dtype=float) for hole in holes],
            refractive_index=spatial_config.slab_refractive_index,
            extinction_coefficient=spatial_config.slab_extinction_coefficient,
            label=f"slab_split_{idx}",
        )
        for idx, height in enumerate(layer_heights)
    ]


def mode_validity(sim_run, mode_count):
    is_valid_mode = []
    for mode_idx in range(mode_count):
        coords_normH_xz, normH_xz = sim_run.get_2d_fields(mode_idx, "ewfd.normH", "xz")
        coords_normE_xz, normE_xz = sim_run.get_2d_fields(mode_idx, "ewfd.normE", "xz")
        coords_normH_yz, normH_yz = sim_run.get_2d_fields(mode_idx, "ewfd.normH", "yz")
        coords_normE_yz, normE_yz = sim_run.get_2d_fields(mode_idx, "ewfd.normE", "yz")

        z_peakH_xz = float(coords_normH_xz[np.argmax(normH_xz), 1])
        z_peakE_xz = float(coords_normE_xz[np.argmax(normE_xz), 1])
        z_peakH_yz = float(coords_normH_yz[np.argmax(normH_yz), 1])
        z_peakE_yz = float(coords_normE_yz[np.argmax(normE_yz), 1])
        valid = (
            (z_peakH_xz <= 1.0) and
            (z_peakE_xz <= 1.0) and
            (z_peakH_yz <= 1.0) and
            (z_peakE_yz <= 1.0)
        )
        is_valid_mode.append(valid)
    return is_valid_mode


def export_center_fields(sim_run, df, field_expr, component_name, save_path):
    for mode_idx in range(len(df)):
        if not df.at[mode_idx, "is_valid"]:
            continue
        coords_re, field_re = sim_run.get_2d_fields(mode_idx, field_expr, "center")
        coords_im, field_im = sim_run.get_2d_fields(mode_idx, f"{field_expr}*(-i)", "center")
        field_im_aligned = interpolate_field(coords_im, field_im, coords_re)
        pd.DataFrame({
            "x": coords_re[:, 0],
            "y": coords_re[:, 1],
            "re": field_re,
            "im": field_im_aligned,
        }).to_parquet(os.path.join(save_path, f"{mode_idx:02d}_{component_name}_center.parquet"))


def run_unified_simulation(boundary, layers, simulation_mode, config, field_expr, component_name, save_path):
    os.makedirs(save_path, exist_ok=True)
    with UnifiedSimulationRun(config=config) as sim_run:
        sim_run.build_and_run(boundary, layers, simulation_mode=simulation_mode)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        df["is_valid"] = mode_validity(sim_run, len(df))

        # Export selected center-plane field for valid modes
        export_center_fields(sim_run, df, field_expr, component_name, save_path)

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def run_hexagon_spatial(a, holes, config, field_expr, component_name, save_path):
    os.makedirs(save_path, exist_ok=True)
    with HexagonSimulationRun(config=config) as sim_run:
        sim_run.build_and_run(a, holes)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        df["is_valid"] = mode_validity(sim_run, len(df))

        # Export selected center-plane field for valid modes
        export_center_fields(sim_run, df, field_expr, component_name, save_path)

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def run_square_spatial(a, holes, config, field_expr, component_name, save_path):
    os.makedirs(save_path, exist_ok=True)
    with SquareSimulationRun(config=config) as sim_run:
        sim_run.build_and_run(a, holes)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        df["is_valid"] = mode_validity(sim_run, len(df))

        # Export selected center-plane field for valid modes
        export_center_fields(sim_run, df, field_expr, component_name, save_path)

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def run_rectangle_spatial(Lx, Ly, holes, config, field_expr, component_name, save_path):
    os.makedirs(save_path, exist_ok=True)
    with RectangleSimulationRun(config=config) as sim_run:
        sim_run.build_and_run(Lx, Ly, holes)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        df["is_valid"] = mode_validity(sim_run, len(df))

        # Export selected center-plane field for valid modes
        export_center_fields(sim_run, df, field_expr, component_name, save_path)

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def compare_simulations(save_path_a, save_path_b, component_name, freq_rtol=1e-4, overlap_atol=0.99, q_rtol=None):
    """Compare valid modes from two simulations.

    Loads eigenfrequencies from both paths, filters to is_valid==True rows
    (sorted by Re(freq)), and compares pairs one-by-one until either side
    runs out of valid modes.

    For each pair:
    - Frequency check: Re(freq), Im(freq), Q within freq_rtol / absolute
    - Complex field overlap: |<field_a | field_b>| / sqrt(||field_a||^2 ||field_b||^2)
      must be >= overlap_atol

    Parameters
    ----------
    save_path_a, save_path_b : str
        Directories produced by the unified and spatial simulations.
    component_name : str
        Field component saved to parquet, for example ``Hz`` or ``Ez``.
    freq_rtol : float
        Relative tolerance on Re(freq) for frequency check.
    overlap_atol : float
        Minimum acceptable |overlap| for field check.
    q_rtol : float, optional
        Relative tolerance on Q. Defaults to ``freq_rtol``.

    Returns
    -------
    dict with keys ``passed`` and ``results_df``.
    """
    if q_rtol is None:
        q_rtol = freq_rtol

    print(f"\n{'='*70}")
    print(f"COMPARING TWO SIMULATIONS ({component_name})")
    print(f"  A: {save_path_a}")
    print(f"  B: {save_path_b}")
    print(f"{'='*70}")

    df_a = pd.read_csv(os.path.join(save_path_a, "eigenfrequencies.csv"))
    df_b = pd.read_csv(os.path.join(save_path_b, "eigenfrequencies.csv"))

    # Filter to valid modes, sort by Re(freq)
    valid_a = (
        df_a[df_a["is_valid"]]
        .reset_index()
        .rename(columns={"index": "mode_idx"})
        .sort_values("re")
        .reset_index(drop=True)
    )
    valid_b = (
        df_b[df_b["is_valid"]]
        .reset_index()
        .rename(columns={"index": "mode_idx"})
        .sort_values("re")
        .reset_index(drop=True)
    )

    n_pairs = min(len(valid_a), len(valid_b))
    print(f"Valid modes: A={len(valid_a)}, B={len(valid_b)}  ->  comparing {n_pairs} pairs")
    assert n_pairs >= 1, "No valid modes available for comparison"

    rows = []
    all_passed = True

    print(f"\n{'Pair':>4}  {'idx_a':>5}  {'idx_b':>5}  "
          f"{'re_a':>10}  {'re_b':>10}  {'re_err':>9}  "
          f"{'im_a':>10}  {'im_b':>10}  {'im_err':>9}  "
          f"{'q_a':>9}  {'q_b':>9}  {'q_err':>9}  "
          f"{'|overlap|':>10}  {'eig_ok':>7}  {'ovlp_ok':>7}")

    for pair_idx in range(n_pairs):
        row_a = valid_a.iloc[pair_idx]
        row_b = valid_b.iloc[pair_idx]
        orig_idx_a = int(row_a["mode_idx"])
        orig_idx_b = int(row_b["mode_idx"])

        re_a, re_b = float(row_a["re"]), float(row_b["re"])
        im_a, im_b = float(row_a["im"]), float(row_b["im"])
        q_a, q_b = float(row_a["q"]), float(row_b["q"])

        re_err = abs(re_a - re_b) / (abs(re_a) + 1e-30)
        im_err = abs(im_a - im_b) / (abs(im_a) + 1e-30)
        q_err = abs(q_a - q_b) / (abs(q_a) + 1e-30)
        re_ok = re_err < freq_rtol
        im_ok = im_err < freq_rtol
        q_ok = q_err < q_rtol
        eig_ok = re_ok and im_ok and q_ok

        # Field overlap
        parquet_a = os.path.join(save_path_a, f"{orig_idx_a:02d}_{component_name}_center.parquet")
        parquet_b = os.path.join(save_path_b, f"{orig_idx_b:02d}_{component_name}_center.parquet")

        overlap = float("nan")
        ovlp_ok = False
        if os.path.exists(parquet_a) and os.path.exists(parquet_b):
            fa = pd.read_parquet(parquet_a)
            fb = pd.read_parquet(parquet_b)

            pts_a = fa[["x", "y"]].to_numpy()
            field_a = fa["re"].to_numpy() + 1j * fa["im"].to_numpy()
            pts_b = fb[["x", "y"]].to_numpy()
            field_b = fb["re"].to_numpy() + 1j * fb["im"].to_numpy()

            # Interpolate B onto A's mesh
            field_b_at_a = (
                interpolate_field(pts_b, np.real(field_b), pts_a)
                + 1j * interpolate_field(pts_b, np.imag(field_b), pts_a)
            )

            integrand = np.conj(field_a) * field_b_at_a
            inner = (
                integrate_field(pts_a, np.real(integrand))
                + 1j * integrate_field(pts_a, np.imag(integrand))
            )
            norm_a = integrate_field(pts_a, np.abs(field_a) ** 2)
            norm_b = integrate_field(pts_a, np.abs(field_b_at_a) ** 2)
            overlap = float(np.abs(inner) / np.sqrt(norm_a * norm_b + 1e-30))
            ovlp_ok = overlap >= overlap_atol
        else:
            print(f"  pair {pair_idx}: {component_name} parquet missing, skipping overlap")

        pair_passed = eig_ok and ovlp_ok
        all_passed = all_passed and pair_passed

        print(f"{pair_idx:4d}  {orig_idx_a:5d}  {orig_idx_b:5d}  "
              f"{re_a:10.4f}  {re_b:10.4f}  {re_err:9.2e}  "
              f"{im_a:10.4e}  {im_b:10.4e}  {im_err:9.2e}  "
              f"{q_a:9.2e}  {q_b:9.2e}  {q_err:9.2e}  "
              f"{overlap:10.6f}  {'PASS' if eig_ok else 'FAIL':>7}  {'PASS' if ovlp_ok else 'FAIL':>7}")

        rows.append({
            "pair": pair_idx,
            "idx_a": orig_idx_a,
            "idx_b": orig_idx_b,
            "re_a": re_a,
            "re_b": re_b,
            "re_err": re_err,
            "im_a": im_a,
            "im_b": im_b,
            "im_err": im_err,
            "q_a": q_a,
            "q_b": q_b,
            "q_err": q_err,
            "overlap": overlap,
            "re_ok": re_ok,
            "im_ok": im_ok,
            "q_ok": q_ok,
            "eig_ok": eig_ok,
            "ovlp_ok": ovlp_ok,
        })

    results_df = pd.DataFrame(rows)
    summary_csv = os.path.join(save_path_a, "comparison_summary.csv")
    results_df.to_csv(summary_csv, index=False)
    print(f"\nSaved comparison summary to {summary_csv}")

    print(f"\n{'='*70}")
    print(f"OVERALL: {'PASS' if all_passed else 'FAIL'}")
    if len(valid_a) != len(valid_b):
        print(f"  WARNING: unequal valid mode counts (A={len(valid_a)}, B={len(valid_b)})")
    print(f"{'='*70}\n")

    return {"passed": all_passed, "results_df": results_df}


def build_hexagon_case(out_root, multilayer=False):
    # --- Parameters copied from tests/test_simulation_spatial.py ---
    a = 0.82
    r_0 = 0.82 / 3  # a/3: double dirac cone
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

    # --- Convert to standard basis coordinates ---
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

    # --- Build geometry ---
    hexagon, holes, info = create_hexagon_design(a, hole_params.T)
    triangles = [h.tolist() for h in holes]

    config = {
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
    }

    case_name = "hexagon_unit_cell_multilayer" if multilayer else "hexagon_unit_cell"
    case_root = os.path.join(out_root, case_name)
    os.makedirs(case_root, exist_ok=True)
    config["layer_mode"] = "two_layer" if multilayer else "single_layer"
    with open(os.path.join(case_root, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    visualize_hexagon_design(
        hexagon, holes, info,
        filename=os.path.join(case_root, "design.png"),
        annotate=False,
    )
    visualize_hexagon_design(
        hexagon, holes, info,
        filename=os.path.join(case_root, "design_annotated.png"),
        annotate=True,
    )

    spatial_config = HexagonSimulationConfig(eigenmode_count=4, mode_type="TE")
    unified_config = unified_config_from_spatial(spatial_config)
    field_expr, component_name = field_expr_for_mode(spatial_config.mode_type)
    layers = (
        two_layers_from_spatial(spatial_config, triangles)
        if multilayer else
        [layer_from_spatial(spatial_config, triangles)]
    )

    return {
        "name": case_name,
        "root": case_root,
        "field_expr": field_expr,
        "component_name": component_name,
        "run_new": lambda path: run_unified_simulation(
            BoundarySpec("hexagon", a=a),
            layers,
            "unit_cell",
            unified_config,
            field_expr,
            component_name,
            path,
        ),
        "run_old": lambda path: run_hexagon_spatial(
            a,
            triangles,
            spatial_config,
            field_expr,
            component_name,
            path,
        ),
    }


def build_square_case(out_root, multilayer=False):
    a = 0.82
    holes = [circle_polygon(0.0, 0.0, 0.18, vertex_count=16)]
    spatial_config = SquareSimulationConfig(eigenmode_count=4, mode_type="TM")
    unified_config = unified_config_from_spatial(spatial_config)
    field_expr, component_name = field_expr_for_mode(spatial_config.mode_type)
    layers = (
        two_layers_from_spatial(spatial_config, holes)
        if multilayer else
        [layer_from_spatial(spatial_config, holes)]
    )

    case_name = "square_unit_cell_multilayer" if multilayer else "square_unit_cell"
    case_root = os.path.join(out_root, case_name)
    os.makedirs(case_root, exist_ok=True)
    with open(os.path.join(case_root, "config.json"), "w") as f:
        json.dump({
            "a": a,
            "holes": [hole.tolist() for hole in holes],
            "mode_type": spatial_config.mode_type,
            "layer_mode": "two_layer" if multilayer else "single_layer",
        }, f, indent=2)

    return {
        "name": case_name,
        "root": case_root,
        "field_expr": field_expr,
        "component_name": component_name,
        "run_new": lambda path: run_unified_simulation(
            BoundarySpec("square", a=a),
            layers,
            "unit_cell",
            unified_config,
            field_expr,
            component_name,
            path,
        ),
        "run_old": lambda path: run_square_spatial(
            a,
            holes,
            spatial_config,
            field_expr,
            component_name,
            path,
        ),
    }


def build_rectangle_case(out_root, multilayer=False):
    Lx = 1.6
    Ly = 1.2
    holes = [circle_polygon(0.0, 0.0, 0.14, vertex_count=16)]
    spatial_config = RectangleSimulationConfig(eigenmode_count=4, mode_type="TM")
    unified_config = unified_config_from_spatial(spatial_config)
    field_expr, component_name = field_expr_for_mode(spatial_config.mode_type)
    layers = (
        two_layers_from_spatial(spatial_config, holes)
        if multilayer else
        [layer_from_spatial(spatial_config, holes)]
    )

    case_name = "rectangle_finite_size_multilayer" if multilayer else "rectangle_finite_size"
    case_root = os.path.join(out_root, case_name)
    os.makedirs(case_root, exist_ok=True)
    with open(os.path.join(case_root, "config.json"), "w") as f:
        json.dump({
            "Lx": Lx,
            "Ly": Ly,
            "holes": [hole.tolist() for hole in holes],
            "mode_type": spatial_config.mode_type,
            "layer_mode": "two_layer" if multilayer else "single_layer",
        }, f, indent=2)

    return {
        "name": case_name,
        "root": case_root,
        "field_expr": field_expr,
        "component_name": component_name,
        "run_new": lambda path: run_unified_simulation(
            BoundarySpec("rectangle", a=Lx, b=Ly),
            layers,
            "finite_size",
            unified_config,
            field_expr,
            component_name,
            path,
        ),
        "run_old": lambda path: run_rectangle_spatial(
            Lx,
            Ly,
            holes,
            spatial_config,
            field_expr,
            component_name,
            path,
        ),
    }


def run_case(case):
    print(f"\n{'#'*70}")
    print(f"CASE: {case['name']}")
    print(f"{'#'*70}")

    new_save_path = os.path.join(case["root"], "new")
    old_save_path = os.path.join(case["root"], "old")

    errors = []

    print(f"\n--- Run unified simulation_utils for {case['name']} ---")
    try:
        df_new = case["run_new"](new_save_path)
        print("unified simulation_utils eigenfrequencies:")
        print(df_new.to_string())
    except Exception as exc:
        errors.append(f"new simulation failed: {exc}")
        traceback.print_exc()

    print(f"\n--- Run spatial baseline for {case['name']} ---")
    try:
        df_old = case["run_old"](old_save_path)
        print("simulation_spatial eigenfrequencies:")
        print(df_old.to_string())
    except Exception as exc:
        errors.append(f"old simulation failed: {exc}")
        traceback.print_exc()

    if errors:
        return {
            "passed": False,
            "results_df": pd.DataFrame({
                "case": [case["name"]],
                "error": ["; ".join(errors)],
            }),
        }

    print(f"\n--- Compare {case['name']} ---")
    try:
        return compare_simulations(
            new_save_path,
            old_save_path,
            case["component_name"],
            freq_rtol=1e-4,
            q_rtol=1e-4,
            overlap_atol=0.99,
        )
    except Exception as exc:
        traceback.print_exc()
        return {
            "passed": False,
            "results_df": pd.DataFrame({
                "case": [case["name"]],
                "error": [f"comparison failed: {exc}"],
            }),
        }


if __name__ == "__main__":
    out_root = str(output_root() / 'test_simulation_utils_spatial_comsol')
    if os.path.exists(out_root):
        shutil.rmtree(out_root)
    os.makedirs(out_root, exist_ok=True)

    cases = [
        build_hexagon_case(out_root),
        build_hexagon_case(out_root, multilayer=True),
        build_square_case(out_root),
        build_square_case(out_root, multilayer=True),
        build_rectangle_case(out_root),
        build_rectangle_case(out_root, multilayer=True),
    ]

    all_results = {}
    failures = []
    for case in cases:
        result = run_case(case)
        all_results[case["name"]] = result["results_df"]
        if not result["passed"]:
            failures.append(case["name"])

    summary = {
        case_name: results_df.to_dict(orient="records")
        for case_name, results_df in all_results.items()
    }
    with open(os.path.join(out_root, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    assert not failures, "Simulation comparison FAILED for: " + ", ".join(failures)
