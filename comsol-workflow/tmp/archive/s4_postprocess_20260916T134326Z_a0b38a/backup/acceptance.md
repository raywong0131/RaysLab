# S4 execution report

Run: 20260916-gamma-v2

Computed cases: 1 / 1.

Scientific acceptance remains pending. A completed solve is not proof of a radiation zero.

Inspect mode_identification.csv, edge_integrals.csv, group_totals.csv, area_comparison.csv,
interface_traces.csv and air_radiation.csv in each case's 01_results directory.

The 3D correction terms and the outgoing air-plane amplitude are distinct outputs.

frozen_field_predictions.csv compares fixed edge means and the fixed spatial Hz field.

zero_estimates.csv contains complex linear-interpolation diagnostics, not solved zeros.

Mesh convergence and A/B parameter trends still require scientific review.

## Reviewed 2026-09-16: needs improvement

P2 numerical acceptance did not pass. Maximum parity error: 0.0220874 (target 0.001); maximum A/B quadrature change: 0.0137704 (target 0.01). The area/boundary identity error is 0.000765117, but convergence is not established. No P3 scan was launched.

The reproducible review and absolute air-field errors are in `tmp/s4_20260916-gamma-v2/scientific_review.json` and `review.py` at repository root. The completed execution record remains a record of computation, not scientific acceptance.
