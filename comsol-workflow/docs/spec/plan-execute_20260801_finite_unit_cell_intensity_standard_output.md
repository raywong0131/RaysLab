# Finite unit-cell intensity standard output: plan and execution record

Date: 2026-08-01
Status: complete

## Objective

Make the approved center-cell-normalized unit-cell intensity map a standard
output for every valid full-finite and quarter-finite mode. The output must not
depend on running the later trend analysis, enabling Fourier decomposition, or
including that mode in an optional Fourier mode subset.

## Plot contract

- Output: `01_results/modeN/13_lattice_fourier_Hz/intensity_maps_normH.png`.
- Source: each mode's reconstructed `Hz_center.parquet` plus the complete
  cavity and cladding cell-center lists.
- Cell intensity: `mean(|Hz|^2)` over the common within-cell sampling grid.
  Since every cell has the same area and sampling contract, the center-normalized
  ratio is the same as the ratio of complete-cell integrals.
- Normalization: use the cavity center cell `(i,j)=(0,0)` independently for
  each mode and plot `I_R/I_0`.
- Region: complete cavity, including its outermost layer, and all cladding
  cells.
- Colormap: `RdYlBu_r`.
- Color scale: linear `0--2`, with `1` at the colorbar midpoint; ticks are
  `0, 0.5, 1.0, 1.5, 2+`, and values at or above `2` use the endpoint color.
- Geometry: unit-cell polygons, equal x-y aspect, with no special center-cell
  outline or marker.

## Implementation

1. The shared finite-lattice postprocessor samples cavity and cladding cells in
   one interpolation pass for modes that undergo Fourier decomposition.
2. Modes omitted from an explicit Fourier subset receive the map through a
   map-only path.
3. When Fourier decomposition is disabled, full finite and quarter finite call
   the map-only path for every valid mode.
4. Quarter finite uses its reconstructed full-field parquet, so both workflows
   share the same plotting and normalization implementation.
5. The existing finalizer moves staged mode maps into the formal mode-centric
   result layout without changing existing Fourier NPZ, CSV, score, or trend
   products.

## Verification record

- Focused tests: `42 passed` for `test_finite_flatness_trend.py`,
  `test_finite_cavity_field_exports.py`, and `test_run_strip_shift_scan.py`
  using a controlled workspace basetemp.
- Syntax: `py_compile` passed for the shared postprocessor, finite, quarter,
  strip, and trend modules.
- Static diff QA: `git diff --check` passed (line-ending notices only).
- Visual QA: rendered 1,261 real cells from the saved shift `0.060` trend data
  and inspected the PNG at original resolution. Equal aspect, cell geometry,
  linear colorbar, center midpoint, and cavity-to-cladding contrast matched the
  agreed contract.

## Rollback

Remove the shared map helpers and their calls from full/quarter finite, remove
the added cladding-cell argument, and remove the focused tests and documentation.
Existing numerical result formats remain otherwise unchanged.
