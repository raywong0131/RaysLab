import os
import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from loguru import logger
from scipy.optimize import linear_sum_assignment

from .energy_recovery import (
    field_to_vector,
    fourier_subspace_energies_field,
    integrate_field,
    interpolate_field,
)
from .basis_utils import standard_to_fourier_basis


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_hz_center(k_dir, mode_idx):
    path = os.path.join(k_dir, "eigenmodes", f"{mode_idx:02d}_Hz_center.parquet")
    df = pd.read_parquet(path)
    pts = df[["x", "y"]].values
    hz = df["re"].values + 1j * df["im"].values
    return pts, hz


def load_hz_center(k_dir, mode_idx):
    """Public wrapper for loading a saved complex center-plane Hz field."""
    return _load_hz_center(k_dir, mode_idx)


def mode_composition_from_frames(df_re, df_im):
    """Return normalized rotational-subspace and Fourier-mode field weights."""
    res = fourier_subspace_energies_field(df_re, df_im, N=6)
    total_energy = float(res["total_energy"])
    if not np.isfinite(total_energy) or total_energy <= 0.0:
        raise ValueError("Mode field has non-positive total energy")

    subspace_values = [res["E_dc"], *res["E_pairs"]]
    subspace_values.append(0.0 if res["E_nyquist"] is None else res["E_nyquist"])
    subspace_weights = np.real(np.asarray(subspace_values, dtype=complex)) / total_energy
    subspace_names = ("s", "p", "d", "f")

    standard_coeffs = field_to_vector(df_re, df_im, N=6)
    fourier_coeffs = standard_to_fourier_basis(standard_coeffs, N=6)
    mode_energies = np.abs(fourier_coeffs) ** 2
    mode_total = float(np.sum(mode_energies))
    if not np.isfinite(mode_total) or mode_total <= 0.0:
        raise ValueError("Mode field has non-positive Fourier coefficient energy")
    mode_weights = mode_energies / mode_total
    mode_names = ("s", "px", "py", "dx", "dy", "f")

    result = {
        f"{name}_weight": float(weight)
        for name, weight in zip(subspace_names, subspace_weights)
    }
    for name, weight in zip(mode_names, mode_weights):
        key = f"{name}_mode_weight" if name in {"s", "f"} else f"{name}_weight"
        result[key] = float(weight)
    result["dominant_subspace"] = subspace_names[int(np.argmax(subspace_weights))]
    result["dominant_mode"] = mode_names[int(np.argmax(mode_weights))]
    return result


def load_mode_composition(k_dir, mode_idx):
    """Load a saved complex center-plane Hz field and analyze its composition."""
    path = os.path.join(k_dir, "eigenmodes", f"{mode_idx:02d}_Hz_center.parquet")
    df = pd.read_parquet(path)
    df_re = pd.DataFrame({0: df["x"], 1: df["y"], 2: df["re"]})
    df_im = pd.DataFrame({0: df["x"], 1: df["y"], 2: df["im"]})
    return mode_composition_from_frames(df_re, df_im)


def assign_mode_candidates(
    previous_frequencies,
    candidate_frequencies,
    overlap_matrix,
    frequency_tolerance,
    minimum_overlap,
):
    """Globally assign unique candidates under frequency and overlap gates."""
    previous = np.asarray(previous_frequencies, dtype=float)
    candidates = np.asarray(candidate_frequencies, dtype=float)
    overlaps = np.asarray(overlap_matrix, dtype=float)
    expected_shape = (len(previous), len(candidates))
    if overlaps.shape != expected_shape:
        raise ValueError(f"overlap_matrix shape {overlaps.shape} != {expected_shape}")
    if frequency_tolerance <= 0.0:
        raise ValueError("frequency_tolerance must be positive")

    frequency_difference = np.abs(previous[:, None] - candidates[None, :])
    allowed = (
        (frequency_difference <= frequency_tolerance)
        & np.isfinite(overlaps)
        & (overlaps >= minimum_overlap)
    )
    score = overlaps - 0.01 * frequency_difference / frequency_tolerance
    real_cost = np.where(allowed, -score, 1.0e6)
    dummy_cost = np.zeros((len(previous), len(previous)), dtype=float)
    cost = np.concatenate([real_cost, dummy_cost], axis=1)
    row_indices, column_indices = linear_sum_assignment(cost)
    assigned_columns = dict(zip(row_indices.tolist(), column_indices.tolist()))

    results = []
    for previous_idx in range(len(previous)):
        candidate_idx = assigned_columns[previous_idx]
        matched = candidate_idx < len(candidates) and allowed[previous_idx, candidate_idx]
        results.append(
            {
                "matched": bool(matched),
                "candidate_index": int(candidate_idx) if matched else None,
                "frequency_difference": (
                    float(frequency_difference[previous_idx, candidate_idx]) if matched else None
                ),
                "overlap": float(overlaps[previous_idx, candidate_idx]) if matched else None,
            }
        )
    return results


def _field_overlap(field1, field2):
    """Normalised |<hz_a | hz_b>| with hz_b interpolated onto the mesh of hz_a.

    Returns a scalar in [0, 1], or 0.0 on any numerical degeneracy.
    """
    pts_a, hz_a = field1
    pts_b, hz_b = field2

    norm_a = float(np.real(integrate_field(pts_a, np.conj(hz_a) * hz_a)))

    hz_b_on_a = (
        interpolate_field(pts_b, np.real(hz_b), pts_a)
        + 1j * interpolate_field(pts_b, np.imag(hz_b), pts_a)
    )
    norm_b = float(np.real(integrate_field(pts_a, np.conj(hz_b_on_a) * hz_b_on_a)))

    inner = integrate_field(pts_a, np.conj(hz_a) * hz_b_on_a)
    return float(np.abs(inner) / np.sqrt(norm_a * norm_b))


def field_overlap(field1, field2):
    """Public wrapper for normalized complex-Hz field overlap."""
    return _field_overlap(field1, field2)


def _field_decomposition(k_dir, mode_idx):
    composition = load_mode_composition(k_dir, mode_idx)
    subspace_names = ("s", "p", "d", "f")
    mode_names = ("s", "px", "py", "dx", "dy", "f")
    dominant_subspace = subspace_names.index(composition["dominant_subspace"])
    dominant_mode = mode_names.index(composition["dominant_mode"])
    dominant_mode_key = (
        f"{composition['dominant_mode']}_mode_weight"
        if composition["dominant_mode"] in {"s", "f"}
        else f"{composition['dominant_mode']}_weight"
    )
    return (
        dominant_subspace,
        composition[f"{composition['dominant_subspace']}_weight"],
        dominant_mode,
        composition[dominant_mode_key],
    )


# ---------------------------------------------------------------------------
# BandConnector
# ---------------------------------------------------------------------------


class BandConnector:
    """Discover k-point directories, load eigenfrequency data, and connect
    eigenmodes into continuous bands across the Brillouin zone.

    Parameters
    ----------
    base_path : str
        Directory that contains the per-k-point sub-directories, each named
        ``kx=<kx_str>_ky=<ky_str>``.
    """

    def __init__(self, base_path):
        self.base_path = base_path

    # ------------------------------------------------------------------
    # Discovery
    # ------------------------------------------------------------------

    def discover_k_points(self):
        """Scan *base_path* for ``kx=..._ky=...`` directories.

        Returns
        -------
        list of (kx_str, ky_str) tuples, sorted lexicographically.
        """
        pattern = re.compile(r"^kx=(?P<kx>[^_]+)_ky=(?P<ky>.+)$")
        k_points = []
        for name in os.listdir(self.base_path):
            m = pattern.match(name)
            if m and os.path.isdir(os.path.join(self.base_path, name)):
                k_points.append((m.group("kx"), m.group("ky")))
        return sorted(k_points)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def load_band_data(self, k_points=None):
        """Load per-k-point ``eigenfrequencies.csv`` files into dicts.

        Parameters
        ----------
        k_points : list of (kx_str, ky_str), optional
            If *None*, calls :meth:`discover_k_points` to find them.

        Returns
        -------
        list of dicts with keys ``kx_str``, ``ky_str``, ``modes`` (DataFrame).
        """
        if k_points is None:
            k_points = self.discover_k_points()

        results = []
        for kx_str, ky_str in k_points:
            k_save_path = os.path.join(self.base_path, f"kx={kx_str}_ky={ky_str}")
            csv_path = os.path.join(k_save_path, "eigenfrequencies.csv")
            if not os.path.exists(csv_path):
                logger.warning(f"Missing eigenfrequencies.csv for kx={kx_str} ky={ky_str}, skipping")
                continue
            df = pd.read_csv(csv_path)
            results.append({
                "kx_str": kx_str,
                "ky_str": ky_str,
                "modes": df,
            })
        return results

    # ------------------------------------------------------------------
    # Band connection
    # ------------------------------------------------------------------

    def connect(
        self,
        band_data,
        freq_tolerance,
        overlap_threshold,
    ):
        """Connect eigenmode results across k-points into continuous bands.

        Parameters
        ----------
        band_data :
            List returned by :meth:`load_band_data`.
        freq_tolerance :
            Maximum frequency difference (same units as the eigenfrequency
            column) to consider two modes as candidates for the same band.
        overlap_threshold :
            Minimum field overlap score required to match two modes.

        Returns
        -------
        List of band dicts, each with keys ``band_id``, ``k_idx``,
        ``mode_idx``, ``freq``, ``kx``, ``ky``, ``q``, ``first_mode_fourier_argmax``,
        ``first_mode_fourier_max``.
        """
        save_path = self.base_path
        n_k = len(band_data)
        kx_arr = np.array([float(e["kx_str"]) for e in band_data])
        ky_arr = np.array([float(e["ky_str"]) for e in band_data])
        freq_tol = freq_tolerance

        k_norm = np.sqrt(kx_arr**2 + ky_arr**2)
        k_angle = np.arctan2(ky_arr, kx_arr) % (2 * np.pi)
        process_order = sorted(range(n_k), key=lambda i: (k_norm[i], k_angle[i]))

        bands = []
        processed = set()

        for pos in process_order:
            df = band_data[pos]["modes"]
            df["band_id"] = -1
            valid_idxs = np.where(df["is_valid"].values)[0]

            if not processed:
                for mode_idx in valid_idxs:
                    bid = len(bands)
                    df.iat[mode_idx, df.columns.get_loc("band_id")] = bid
                    q_val = float(df.iat[mode_idx, df.columns.get_loc("q")]) if "q" in df.columns else float("NaN")
                    
                    # Compute decomposition for first mode
                    k_dir = os.path.join(
                        save_path, 
                        f"kx={band_data[pos]['kx_str']}_ky={band_data[pos]['ky_str']}"
                    )
                    dominant_subspace, dominant_subspace_energy, dominant_mode, dominant_mode_energy = _field_decomposition(k_dir, int(mode_idx))
                    
                    bands.append({
                        "band_id":  bid,
                        "dominant_subspace": dominant_subspace,
                        "dominant_subspace_energy": dominant_subspace_energy,
                        "dominant_mode": dominant_mode,
                        "dominant_mode_energy": dominant_mode_energy,
                        "kx":       [band_data[pos]["kx_str"]],
                        "ky":       [band_data[pos]["ky_str"]],
                        "freq":     [float(df.iat[mode_idx, df.columns.get_loc("re")])],
                        "q":        [q_val],
                        "mode_idx": [int(mode_idx)],
                    })
                processed.add(pos)
                continue

            pred_pos = min(
                processed,
                key=lambda p: (kx_arr[p] - kx_arr[pos])**2 + (ky_arr[p] - ky_arr[pos])**2,
            )
            df_pred = band_data[pred_pos]["modes"]
            df_pred_valid = df_pred[df_pred["is_valid"]]
            taken = {}

            if valid_idxs.size and not df_pred_valid.empty:
                pred_kdir = os.path.join(
                    save_path,
                    f"kx={band_data[pred_pos]['kx_str']}_ky={band_data[pred_pos]['ky_str']}",
                )
                curr_kdir = os.path.join(
                    save_path,
                    f"kx={band_data[pos]['kx_str']}_ky={band_data[pos]['ky_str']}",
                )
                ref_fields = {idx: _load_hz_center(pred_kdir, idx) for idx in df_pred_valid.index}
                cand_fields = {idx: _load_hz_center(curr_kdir, idx) for idx in valid_idxs}

                for bid in sorted(df_pred_valid["band_id"].unique()):
                    if bid < 0:
                        continue
                    pred_mode_idx = int(df_pred.index[df_pred["band_id"] == bid][0])
                    pred_freq = df_pred.iat[pred_mode_idx, df_pred.columns.get_loc("re")]
                    cands = [
                        int(mode_idx)
                        for mode_idx in valid_idxs
                        if abs(df.iat[mode_idx, df.columns.get_loc("re")] - pred_freq) <= freq_tol
                    ]

                    rf = ref_fields[pred_mode_idx]
                    scores = {m: _field_overlap(rf, cand_fields[m]) for m in cands}
                    ranked = sorted(cands, key=lambda m: scores[m], reverse=True)
                    chosen_mode_idx = None

                    for best_mode_idx in ranked:
                        if scores[best_mode_idx] < overlap_threshold:
                            break
                        if best_mode_idx not in taken:
                            chosen_mode_idx = best_mode_idx
                            break

                    if chosen_mode_idx is None:
                        continue

                    freq = float(df.iat[chosen_mode_idx, df.columns.get_loc("re")])
                    q_val = float(df.iat[chosen_mode_idx, df.columns.get_loc("q")])

                    bands[bid]["kx"].append(band_data[pos]["kx_str"])
                    bands[bid]["ky"].append(band_data[pos]["ky_str"])
                    bands[bid]["freq"].append(freq)
                    bands[bid]["q"].append(q_val)
                    bands[bid]["mode_idx"].append(int(chosen_mode_idx))

                    df.iat[chosen_mode_idx, df.columns.get_loc("band_id")] = bid
                    taken[chosen_mode_idx] = bid

            for mode_idx in valid_idxs:
                if mode_idx not in taken:
                    bid = len(bands)
                    df.iat[mode_idx, df.columns.get_loc("band_id")] = bid
                    q_val = float(df.iat[mode_idx, df.columns.get_loc("q")]) if "q" in df.columns else float("NaN")
                    
                    # Compute decomposition for first mode
                    k_dir = os.path.join(save_path, f"kx={band_data[pos]['kx_str']}_ky={band_data[pos]['ky_str']}")
                    dominant_subspace, dominant_subspace_energy, dominant_mode, dominant_mode_energy = _field_decomposition(k_dir, int(mode_idx))
                    
                    bands.append({
                        "band_id":  bid,
                        "dominant_subspace": dominant_subspace,
                        "dominant_subspace_energy": dominant_subspace_energy,
                        "dominant_mode": dominant_mode,
                        "dominant_mode_energy": dominant_mode_energy,
                        "kx":       [band_data[pos]["kx_str"]],
                        "ky":       [band_data[pos]["ky_str"]],
                        "freq":     [float(df.iat[mode_idx, df.columns.get_loc("re")])],
                        "q":        [q_val],
                        "mode_idx": [int(mode_idx)],
                    })

            processed.add(pos)

        return bands

    # ------------------------------------------------------------------
    # Save / load
    # ------------------------------------------------------------------

    @staticmethod
    def save(bands, csv_path):
        """Write per-band per-k-point data to ``<save_path>/<name>.csv``.

        Returns the DataFrame that was written.
        """
        dfs = []
        for band in bands:
            n = len(band["freq"])
            df = pd.DataFrame({
                "band_id": [band["band_id"]] * n,
                "dominant_subspace": [band["dominant_subspace"]] * n,
                "dominant_subspace_energy": [band["dominant_subspace_energy"]] * n,
                "dominant_mode": [band["dominant_mode"]] * n,
                "dominant_mode_energy": [band["dominant_mode_energy"]] * n,
                "kx": band["kx"],
                "ky": band["ky"],
                "freq": band["freq"],
                "q": band["q"],
                "mode_idx": band["mode_idx"],
            })
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        tmp_path = csv_path + ".tmp"
        df.to_csv(tmp_path, index=False)
        os.replace(tmp_path, csv_path)
        return df

    @staticmethod
    def load(csv_path):
        """Reconstruct band dicts from a CSV saved by :meth:`save`."""
        df = pd.read_csv(csv_path, dtype={"kx": str, "ky": str})
        bands = []
        for band_id, grp in df.groupby("band_id"):
            bands.append({
                "band_id":  int(band_id),
                "dominant_subspace": int(grp.iloc[0]["dominant_subspace"]),
                "dominant_subspace_energy": grp.iloc[0]["dominant_subspace_energy"],
                "dominant_mode": int(grp.iloc[0]["dominant_mode"]),
                "dominant_mode_energy": grp.iloc[0]["dominant_mode_energy"],
                "kx":       grp["kx"].tolist(),
                "ky":       grp["ky"].tolist(),
                "freq":     grp["freq"].tolist(),
                "q":        grp["q"].tolist(),
                "mode_idx": grp["mode_idx"].tolist(),
                
            })
        return bands

    # ------------------------------------------------------------------
    # High-level interface
    # ------------------------------------------------------------------

    def get_connected_bands(
        self,
        k_points = None,
        force_recompute = False,
        freq_tolerance = 5,
        overlap_threshold = 0.0,
    ):
        csv_path = os.path.join(self.base_path, "connected_bands.csv")
        if not force_recompute and os.path.exists(csv_path):
            return self.load(csv_path)

        band_data = self.load_band_data(k_points)
        bands = self.connect(band_data, freq_tolerance=freq_tolerance, overlap_threshold=overlap_threshold)
        self.save(bands, csv_path)
        return bands

