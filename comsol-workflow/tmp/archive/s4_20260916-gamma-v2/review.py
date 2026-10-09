"""Reproduce the S4 single-point review from saved CSVs; no COMSOL import."""
import csv
import json
from pathlib import Path

import numpy as np

task = Path(__file__).resolve().parent
manifest = json.loads((task / "99_config/s4_manifest.json").read_text(encoding="utf8"))
output = Path(manifest["output_directory"])
results = output / manifest["cases"][0]["case_id"] / "01_results"


def read(name):
    with (results / name).open(encoding="utf8", newline="") as stream:
        return list(csv.DictReader(stream))


def value(row, key):
    return complex(float(row[key + "_re"]), float(row[key + "_im"]))


groups, area = read("group_totals.csv"), read("area_comparison.csv")
air, traces = read("air_radiation.csv"), read("interface_traces.csv")
mode = next(row for row in read("mode_identification.csv") if row["mode_idx"] == "1")
assert len(read("edge_geometry.csv")) == 18
assert len(read("edge_integrals.csv")) == 108
assert len(groups) == len(area) == 6 and len(traces) == 108 and len(air) == 2
changes = []
for fine in (row for row in groups if row["refinement"] == "2"):
    coarse = next(row for row in groups if row["refinement"] == "1" and row["z_over_H"] == fine["z_over_H"])
    changes.append({"z_over_H": float(fine["z_over_H"]), **{
        key: abs(value(fine, key) - value(coarse, key)) / abs(value(fine, key))
        for key in ("QAx", "QBx")}})
air_vectors = [np.array([value(row, "cx_at_z0"), value(row, "cy_at_z0")]) for row in air]
air_difference = float(np.linalg.norm(air_vectors[1] - air_vectors[0]))
fine_area = [row for row in area if row["refinement"] == "2"]
scale = min(float(row["amplitude_scale"]) for row in fine_area)
review = {
    "status": "needs_improvement",
    "parameter_identity": manifest["cases"][0],
    "source_output": str(output),
    "mode_idx": 1,
    "frequency_thz": float(mode["frequency_thz"]),
    "parity_max_error": max(float(mode[k]) for k in ("odd_x_error", "even_y_error")),
    "identity_max_error": max(float(row[k]) for row in fine_area for k in ("Ex_identity_error", "Ey_identity_error")),
    "outer_max_error": max(float(row["outer_error"]) for row in fine_area),
    "group_refinement_relative_changes": changes,
    "group_refinement_max_change": max(row[key] for row in changes for key in ("QAx", "QBx")),
    "identity_convergence": "not monotone across all sections/components; two resolutions do not establish convergence",
    "trace_max_relative_difference": max(float(row["Hz_trace_relative_difference"]) for row in traces),
    "trace_median_by_offset_m": {
        offset: float(np.median([float(row["Hz_trace_relative_difference"]) for row in traces if row["offset_m"] == offset]))
        for offset in sorted({row["offset_m"] for row in traces})},
    "all_trace_domains_distinct": all(row["distinct_domain_each_point"] == "True" for row in traces),
    "air_height_absolute_difference_V_per_m": air_difference,
    "air_height_relative_difference": air_difference / float(np.linalg.norm(air_vectors[0])),
    "air_height_difference_over_finite_section_scale": air_difference / scale,
    "finite_section_scale_V_per_m": scale,
    "air_cy_absolute_difference_V_per_m": float(abs(air_vectors[1][1] - air_vectors[0][1])),
    "air_cy_at_z0_magnitudes_V_per_m": [float(abs(vector[1])) for vector in air_vectors],
    "gates": {"identity_target": 1e-3, "parity_target": 1e-3, "group_change_target": .01},
    "conclusion": "P2 numerical acceptance fails parity and group convergence; no P3 scan or radiation-zero claim",
    "limitations": ["Single geometry and mesh cannot establish zeta trends or a radiation zero",
                    "Air relative error is ill-conditioned near a zero; absolute and finite-scale errors are also reported",
                    "Nonzero unwanted polarization and interface differences require mesh/integration review"],
}
assert review["identity_max_error"] < 1e-3
assert review["parity_max_error"] > 1e-3
assert review["group_refinement_max_change"] > .01
(task / "scientific_review.json").write_text(json.dumps(review, indent=2) + "\n", encoding="utf8")
print(json.dumps(review, indent=2))
