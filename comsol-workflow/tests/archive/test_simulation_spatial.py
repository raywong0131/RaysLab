from tests._paths import output_root
import os
import json
import numpy as np
import pandas as pd

from comsol_workflow.basis_utils import get_fourier_basis_min_nonzero_values, fourier_to_standard_basis
from comsol_workflow.geometry_utils import create_hexagon_design, visualize_hexagon_design
from comsol_workflow.simulation_utils import SimulationRun
from comsol_workflow.simulation_spatial.hexagon_unit_cell import SimulationRun as SimulationRunSpatial
from comsol_workflow.energy_recovery import interpolate_field, integrate_field


def run_simulation(a, triangles, save_path):
    os.makedirs(save_path, exist_ok=True)
    with SimulationRun() as sim_run:
        sim_run.build_and_run(a, triangles)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        is_valid_mode = []
        for mode_idx in range(len(df)):
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
        df["is_valid"] = is_valid_mode

        # Export Hz fields for valid modes
        for mode_idx in range(len(df)):
            if not df.at[mode_idx, "is_valid"]:
                continue
            coords_re, reHz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz", "center")
            coords_im, imHz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz*(-i)", "center")
            imHz_aligned = interpolate_field(coords_im, imHz, coords_re)
            pd.DataFrame({
                "x": coords_re[:, 0], "y": coords_re[:, 1],
                "re": reHz, "im": imHz_aligned,
            }).to_parquet(os.path.join(save_path, f"{mode_idx:02d}_Hz_center.parquet"))

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def run_simulation_spatial(a, triangles, save_path):
    os.makedirs(save_path, exist_ok=True)
    with SimulationRunSpatial() as sim_run:
        sim_run.build_and_run(a, triangles)
        sim_run.model.save(os.path.join(save_path, "design.mph"))
        eigenfrequencies = sim_run.get_eigenfrequencies()
        df = pd.DataFrame(eigenfrequencies, columns=["re", "im", "q"])

        # Step 1: check validity for each mode
        is_valid_mode = []
        for mode_idx in range(len(df)):
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
        df["is_valid"] = is_valid_mode

        # Export Hz fields for valid modes
        for mode_idx in range(len(df)):
            if not df.at[mode_idx, "is_valid"]:
                continue
            coords_re, reHz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz", "center")
            coords_im, imHz = sim_run.get_2d_fields(mode_idx, "ewfd.Hz*(-i)", "center")
            imHz_aligned = interpolate_field(coords_im, imHz, coords_re)
            pd.DataFrame({
                "x": coords_re[:, 0], "y": coords_re[:, 1],
                "re": reHz, "im": imHz_aligned,
            }).to_parquet(os.path.join(save_path, f"{mode_idx:02d}_Hz_center.parquet"))

    df.to_csv(os.path.join(save_path, "eigenfrequencies.csv"), index=False)
    return df


def compare_simulations(save_path_a, save_path_b, freq_rtol=1e-4, overlap_atol=0.99, q_rtol=None):
    """Compare valid modes from two simulations.

    Loads eigenfrequencies from both paths, filters to is_valid==True rows
    (sorted by Re(freq)), and compares pairs one-by-one until either side
    runs out of valid modes.

    For each pair:
    - Frequency check: Re(freq), Im(freq), Q within freq_rtol / absolute
    - Complex Hz field overlap: |<hz_a | hz_b>| / sqrt(||hz_a||^2 ||hz_b||^2)
      must be >= overlap_atol

    Parameters
    ----------
    save_path_a, save_path_b : str
        Directories produced by run_simulation / run_simulation_spatial.
    freq_rtol : float
        Relative tolerance on Re(freq) for frequency check.
    overlap_atol : float
        Minimum acceptable |overlap| for Hz field check.
    q_rtol : float, optional
        Relative tolerance on Q. Defaults to ``freq_rtol``.

    Returns
    -------
    dict with keys ``passed`` and ``results_df``.
    """
    if q_rtol is None:
        q_rtol = freq_rtol

    print(f"\n{'='*70}")
    print("COMPARING TWO SIMULATIONS")
    print(f"  A: {save_path_a}")
    print(f"  B: {save_path_b}")
    print(f"{'='*70}")

    df_a = pd.read_csv(os.path.join(save_path_a, "eigenfrequencies.csv"))
    df_b = pd.read_csv(os.path.join(save_path_b, "eigenfrequencies.csv"))

    # Filter to valid modes, sort by Re(freq)
    valid_a = df_a[df_a["is_valid"]].sort_values("re").reset_index(drop=True)
    valid_b = df_b[df_b["is_valid"]].sort_values("re").reset_index(drop=True)

    n_pairs = min(len(valid_a), len(valid_b))
    print(f"Valid modes: A={len(valid_a)}, B={len(valid_b)}  ->  comparing {n_pairs} pairs")

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
        orig_a = int(row_a.name) if "index" not in row_a else int(row_a["index"])
        orig_b = int(row_b.name) if "index" not in row_b else int(row_b["index"])

        # Recover original mode indices (positions in the unsorted df)
        orig_idx_a = df_a.index[df_a["is_valid"] & (df_a["re"] == row_a["re"])][0]
        orig_idx_b = df_b.index[df_b["is_valid"] & (df_b["re"] == row_b["re"])][0]

        re_a, re_b = float(row_a["re"]), float(row_b["re"])
        im_a, im_b = float(row_a["im"]), float(row_b["im"])
        q_a,  q_b  = float(row_a["q"]),  float(row_b["q"])

        re_err = abs(re_a - re_b) / (abs(re_a) + 1e-30)
        im_err = abs(im_a - im_b) / (abs(im_a) + 1e-30)
        q_err = abs(q_a - q_b) / (abs(q_a) + 1e-30)
        re_ok = re_err < freq_rtol
        im_ok = im_err < freq_rtol
        q_ok = q_err < q_rtol
        eig_ok = re_ok and im_ok and q_ok

        # Field overlap
        parquet_a = os.path.join(save_path_a, f"{orig_idx_a:02d}_Hz_center.parquet")
        parquet_b = os.path.join(save_path_b, f"{orig_idx_b:02d}_Hz_center.parquet")

        overlap = float("nan")
        ovlp_ok = False
        if os.path.exists(parquet_a) and os.path.exists(parquet_b):
            fa = pd.read_parquet(parquet_a)
            fb = pd.read_parquet(parquet_b)

            pts_a = fa[["x", "y"]].to_numpy()
            hz_a  = fa["re"].to_numpy() + 1j * fa["im"].to_numpy()
            pts_b = fb[["x", "y"]].to_numpy()
            hz_b  = fb["re"].to_numpy() + 1j * fb["im"].to_numpy()

            # Interpolate B onto A's mesh
            hz_b_at_a = (
                interpolate_field(pts_b, np.real(hz_b), pts_a)
                + 1j * interpolate_field(pts_b, np.imag(hz_b), pts_a)
            )

            integrand = np.conj(hz_a) * hz_b_at_a
            inner = (
                integrate_field(pts_a, np.real(integrand))
                + 1j * integrate_field(pts_a, np.imag(integrand))
            )
            norm_a = integrate_field(pts_a, np.abs(hz_a) ** 2)
            norm_b = integrate_field(pts_a, np.abs(hz_b_at_a) ** 2)
            overlap = float(np.abs(inner) / np.sqrt(norm_a * norm_b + 1e-30))
            ovlp_ok = overlap >= overlap_atol
        else:
            print(f"  pair {pair_idx}: Hz parquet missing, skipping overlap")

        pair_passed = eig_ok and ovlp_ok
        all_passed = all_passed and pair_passed

        print(f"{pair_idx:4d}  {orig_idx_a:5d}  {orig_idx_b:5d}  "
              f"{re_a:10.4f}  {re_b:10.4f}  {re_err:9.2e}  "
              f"{im_a:10.4e}  {im_b:10.4e}  {im_err:9.2e}  "
              f"{q_a:9.2e}  {q_b:9.2e}  {q_err:9.2e}  "
              f"{overlap:10.6f}  {'PASS' if eig_ok else 'FAIL':>7}  {'PASS' if ovlp_ok else 'FAIL':>7}")

        rows.append({
            "pair": pair_idx,
            "idx_a": orig_idx_a, "idx_b": orig_idx_b,
            "re_a": re_a, "re_b": re_b, "re_err": re_err,
            "im_a": im_a, "im_b": im_b, "im_err": im_err,
            "q_a": q_a, "q_b": q_b, "q_err": q_err,
            "overlap": overlap,
            "re_ok": re_ok, "im_ok": im_ok, "q_ok": q_ok,
            "eig_ok": eig_ok, "ovlp_ok": ovlp_ok,
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


if __name__ == "__main__":
    out_root = str(output_root() / 'test_simulation_spatial')

    # --- Parameters ---
    a = 0.82
    r_0 = 0.82 / 3  # a/3: double dirac cone
    b_0 = 0.23
    d = 0.02

    symmetry_config = {
        'r': 5,
        'theta': None,
        'b': 1,
        'phi': None,
    }

    r_f0 = 1.0          # relative scale for r
    r_fs = 0.0          # relative modulation for r
    theta_fs = 0.0      # absolute modulation for theta
    b_square_f0 = 1.0   # relative scale for b^2
    b_square_fs = 0.0   # relative modulation for b^2
    phi_fs = 0.0        # absolute modulation for phi

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

    hole_params = np.stack([r, theta, b, phi], axis=0)  # (4, 6)

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

    os.makedirs(out_root, exist_ok=True)
    with open(os.path.join(out_root, "config.json"), "w") as f:
        json.dump(config, f, indent=2)

    visualize_hexagon_design(
        hexagon, holes, info,
        filename=os.path.join(out_root, "design.png"),
        annotate=False,
    )
    visualize_hexagon_design(
        hexagon, holes, info,
        filename=os.path.join(out_root, "design_annotated.png"),
        annotate=True,
    )

    # --- Run simulation (simulation_utils.SimulationRun) ---
    sim_save_path = os.path.join(out_root, "simulation")
    df_sim = run_simulation(a, triangles, sim_save_path)
    print("simulation_utils eigenfrequencies:")
    print(df_sim.to_string())

    # --- Run simulation (simulation_spatial.hexagon_unit_cell.SimulationRun) ---
    sim_spatial_save_path = os.path.join(out_root, "simulation_spatial")
    df_sim_spatial = run_simulation_spatial(a, triangles, sim_spatial_save_path)
    print("simulation_spatial eigenfrequencies:")
    print(df_sim_spatial.to_string())

    # --- Compare the two simulations ---
    results = compare_simulations(sim_save_path, sim_spatial_save_path)
    assert results["passed"], "Simulation comparison FAILED"

