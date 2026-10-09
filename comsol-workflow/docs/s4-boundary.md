# S4 boundary-contribution calculation

Run commands from `I:/codeXproject/comsol-workflow`. This entry is independent
of the standard k scan and does not change `scripts/parameter.json`.


## 当前汇总图目录

`results/S4_boundary_analysis/10_overview/` 仅平铺当前四张主图PNG，不建立来源或批次子文件夹。
被替代及辅助图归入 `90_history/10_overview/`，保留来源分组；PDF继续使用 `11_pdf/`。
新的完整计算结果更新四个固定文件名，并保留上一版历史图。原始模型和科学数据不随图目录整理迁移。

## Single-point preparation

For an explicitly authorized self-contained task, set the same session root
before both preparation and execution:

```powershell
$env:COMSOL_WORKFLOW_S4_SESSION_ROOT = "I:/codeXproject/comsol-workflow/scripts/.out/unit_cell_2D/S4_control_b245_eta0.95_zeta1_mesh5_20260917"
```

This optional directory must be a direct child of `scripts/.out/unit_cell_2D`.
Preparation, backend checks and audit staging then use its `tmp/` child;
formal runs use its `results/` child. The parameter-derived run names,
manifest hashes and no-overwrite checks remain unchanged. Put launcher logs
under the same container and set `TMP`, `TEMP`, and `TMPDIR` to its runtime
temporary directory when starting a process. This changes paths only, not
the solver, mode selection or integration. Without this variable the legacy
paths below remain unchanged.

```powershell
uv run python -m scripts.run_main.run_boundary_analysis prepare --run-id gamma-review --zeta 1.156
```

This reads `unit_cell_2d` and the root mesh setting. It creates a snapshot,
code hashes, and `99_config/s4_manifest.json` under `tmp/s4_<run_id>/`.
Preparation and backend checks never create a formal result directory.
The S4 baseline defaults to the approved zeta=1.156 even if the shared value
changes; the source value is recorded separately. Other defaults inspected
on 2026-09-16 were b0=245 nm, eta=0.96,
mesh=5, two requested eigenmode pairs, and k=(0,0). Other physical defaults
come from the existing band module. Preparation prints the actual values.
The default sections are `--z-nm 0 50 100`. The primary section and Hz RMS
normalization plane are z=0; the other two sections are numerical references.
Sections must start at zero and lie within the closed simulated half slab.
For H=200 nm, z=100 nm is its top interface: interpolation selects the
`sel_slab_layer` domains at the exact requested coordinate, including hole air,
so the top-plane data are the slab-side trace. No small z offset is substituted.

## Explicit execution

```powershell
uv run python -m scripts.run_main.run_boundary_analysis run tmp/s4_gamma-review
```

This starts one Gamma calculation, using the existing model builder. It saves
an MPH checkpoint before extraction. It does not call the complete k scan,
the standard plotting pipeline, or a shared-parameter update function.
Existing COMSOL processes block execution; sharing resources is not implicit.
An attempted run is never overwritten. Code or snapshot changes require a
new preparation. Manifest integrity and resolved snapshot/output paths are
checked before execution. Process/PATH and local or network license checks,
then COMSOL/Wave Optics product checkout, precede formal directory creation.
Partial outputs are retained for diagnosis. Runtime diagnostics before model
creation stay in the preparation directory; model progress is recorded in
each case's `00_model/comsol_progress.log`.

The target uses the manuscript convention: Hz is odd under x reflection and
even under y reflection. The selection also checks the existing p weights.
The program retains internal mode indices and component names. It does not
equate p2 with the manuscript target without checking its field.
Mixed or ambiguous states, including unresolved degeneracy at eta=zeta=1,
stop the case and preserve the candidate table; automatic rotation inside a
degenerate subspace is not implemented in this version.

## Explicit local scan and mesh checks

```powershell
uv run python -m scripts.run_main.run_boundary_analysis prepare --run-id s4-zeta-five --zeta 1.154 1.155 1.156 1.157 1.158 --max-cases 5
uv run python -m scripts.run_main.run_boundary_analysis prepare --run-id s4-mesh-three --mesh 5 4 3 --max-cases 3
```

Preparation alone does not authorize or launch these solves. `--eta` provides
a task-local position change. `--parameter-file` accepts a complete snapshot
through the existing shared loader. The program uses one reference field per
mesh for the length-only and fixed-Hz predictions. These approximations use
the actual reference geometry, not an assumed zeta=1 reference.

## Outputs

The manifest records the exact formal output directory under
`scripts/.out/unit_cell_2D/S4_b<b0>_eta<eta>_zeta<zeta>_mesh<M>_Gamma_<N>case_<run_id>/`.
The geometry label identifies the first case; all case identities and the
explicit total solve budget are stored in the manifest.
Each case contains `00_model`, `01_results`, and `99_config`.

| File | Meaning |
|---|---|
| `mode_identification.csv` | Candidate frequency, p weights, parity errors and overlap |
| `edge_geometry.csv` | Actual 18 edges, vertices, lengths and dielectric-to-air normals |
| `edge_integrals.csv` | Complex Hz line integrals and x/y contributions at two resolutions |
| `group_totals.csv` | Six A edges and twelve B edges, complex sums and relative phase |
| `area_comparison.csv` | Raw Ex/Ey area integrals (V m), averages (V/m), boundary, outer and z-derivative terms; display/reference sections explicitly marked |
| `area_field_z0.npz` | Complex Ex/Ey at the actual quadrature nodes and weights, plus a separate masked display raster at z=0; SI units |
| `interface_traces.csv` | Two material traces, domain IDs and two distances from each interface |
| `air_radiation.csv` | Independent G=0 amplitudes at two air heights, including propagation correction |
| `frozen_field_predictions.csv` | Length-only and fixed-spatial-Hz predictions for requested geometries |
| `zero_estimates.csv` | Complex linear-interpolation estimates, only when multiple zeta values exist |
| `s4_edges_py.png/.svg` | py equation comparison at z=0: Ey area map/integral, numbered boundary geometry, and signed edge contributions with their complex sum |
| `s4_edges_py_Ex_reference.png/.svg` | Ex transverse reference of the same py mode; not a px-mode result |
| `s4_group_contributions.png/.svg` | A, B and A+B real/imaginary bars at z=0 for x/y components |
| `s4_area_comparison.png/.svg` | Direct electric-field area average versus hole, outer, 3D and reconstructed terms at z=0, in V/m |
| `s4_validation_py.png/.svg` | py / Ey S31 diagnostic: area versus boundary, integration sensitivity, and outer/thickness contribution balance at z=0 |

Complex CSV quantities have `_re` and `_im` columns. Coordinates, lengths,
areas, field values and derivatives use SI units. E, H and complex frequency
are conjugated together to express the results with exp(-i omega t).
One Hz RMS normalization factor multiplies all field components.

The current model is a three-dimensional half slab. The program retains the
z-derivative and outer-boundary terms when comparing the area and boundary
integrals. It does not identify the bare boundary sum with the field in air.
Two numerical resolutions estimate integration error; this is not a substitute
for mesh convergence. A real-part sign change alone is not reported as a zero.
Interpolated estimates require actual solved-point and convergence checks.
Automated root refinement, final figures, and a scientific acceptance decision
are not included in this version.

## Boundary-contribution figures

For a completed C6v p-pair decomposition, use:

```powershell
uv run python -B -m scripts.analysis.plot_boundary_control <session-root>
```

The session must contain `results/symmetry_decomposition_v1/`. Add
`--extract-field` once if its `80_logs/area_field_z0.npz` cache is missing:
this explicitly reads the saved MPH without meshing, solving or saving it,
reuses the accepted constant complex coefficients, and checks the cached
quadrature against the saved Ey/Ex integrals. Later plotting needs no COMSOL.
The current control main figure goes to session-level
`10_overview/s4_edges_py.png/.svg`. It uses dimensional V/m values: the left
formula is `mean(Ey) = integral Ey(r) d²r / A_cell`; the middle formula
evaluates the weighted x-derivative of Hz; the right formula evaluates the
existing hole-boundary sum. All three captions start with `mean(Ey)` and
retain the complex results in one common phase. The field and edge
colors/bars show real parts in V/m after multiplication by the same unit
phasor, chosen to make the direct Ey average positive real. This changes no
amplitude, relative phase or complex relative error; there is no area-amplitude
normalization. Original and displayed values are both saved. No overall formula
or footer is shown. The two spatial frames stay square and the bar frame
stays 4:3, with equal gaps and at most five bar-axis ticks.

Add `--extract-gradient` once to integrate
`d(laginterp(2,ewfd.Hz),x)/epsilon` on the saved z=0 solution: two Python
quadrature levels and native orders 4/8/12 with explicit material sides.
This never runs the mesh/solver or saves the model. The middle caption uses
the fine Python level, matching the direct-area and boundary calculations;
the native values are independent integration checks. Each original mode
keeps its complex frequency in the Maxwell conversion; the p-pair is not
labeled as an exact single-frequency eigenstate. In 3D, the gradient-only
expression omits the thickness term and need not equal the direct Ey average.

`80_logs/method_comparison.csv` and `edge_contributions_vm.csv` contain
dimensional results. Gradient components, extraction identity, field cache,
normalized reference tables and Chinese analysis also stay in `80_logs/`.
Existing source results and the previously generated normalized
`s31_diagnostic_py` figure are preserved. Subsequent plotting uses the caches
without COMSOL.

The reviewed `gamma-v2/case000_eta0.96_zeta1.156_mesh5` case uses a compact
layout: `10_overview/` contains only `s4_edges_py` and `s31_diagnostic_py`
(PNG/SVG); `80_logs/` contains CSV tables, the area NPZ, eigenmode Parquet
files, the S31 report/audit and a ZIP archive of the four supplementary figure
pairs. The MPH/progress and configuration files stay in `00_model/` and
`99_config/`. No numerical data are discarded. The case's
`80_logs/output_layout.json` records every original file, its final location
or ZIP member, size and SHA-256.

Run-level tables and reports belong in the run's own `80_logs/`:
`case_table.csv`, `frozen_field_predictions.csv`, `acceptance.md`, and
`zero_estimates.csv` when interpolation estimates exist. They summarize
cases and are separate from each case's `80_logs/` numerical data.
The run/report/re-extraction entry points use these paths, while run
configuration stays in the run's `99_config/`.

Readers support both this compact layout and legacy `01_results/` cases.
For a compact case the default plot command (including `--edges-only`)
updates only `10_overview/s4_edges_py.png/.svg`; it preserves the independent
S31 audit figure. `--validation-only` remains an explicit opt-in for the older
integration diagnostic. Re-extraction stages and publishes the matching
layout. Existing audit records retain their historical source paths.

Each newly computed case automatically generates the five PNG figures and
their SVG versions (rasterized E map, vector axes/annotations) in the existing `01_results/` directory. To generate
only these figures from a saved case, without COMSOL or changes to CSV/MPH:

```powershell
uv run python -m scripts.analysis.plot_boundary_analysis <run-directory>/<case-directory>
```

The plotter uses the highest saved integration refinement. It verifies all
18 numbered edges at each section and checks their complex sums against
`group_totals.csv`. Missing edges, nonfinite contributions or inconsistent
group totals stop plotting.

The two three-panel figures use a single complex reference `I_y = integral Ey dS`
at z=0 and the highest saved refinement. The Ey integral becomes exactly 1;
the Ex reference uses the same denominator, so its integral remains `I_x/I_y`.
All three amplitude panels are dimensionless. The boundary geometry colorbar
and edge bars share a symmetric range fitted to the largest x/y edge contribution;
the field map retains its wider symmetric range. This changes only the display
limits, not the common area normalization, and is stated on the figure.
Coordinates remain in micrometres. No panel is divided by its
own maximum. Zero or numerically unresolved I_y is rejected, not replaced by
an arbitrary divisor. Raw CSV/NPZ data remain unchanged.
The edge geometry shows `Re(c_e)` with RdBu_r. Labels A1/A4 identify group A;
B2/B3/B5/B6 identify group B. Small edge labels correspond to the CSV edge IDs;
arrows point from dielectric into each hole. Bar colors identify the groups.
Each edge figure has three panels: the 2D electric-field map on the left,
the numbered boundary geometry in the middle, and compact `Re(c_e)` bars on
the right. The bars share a symmetric y range for x/y with at most five
numbered y ticks. The right plotting frame
has width:height = 4:3; the first two frames stay square, with aligned top/bottom
edges and equal horizontal gaps between the three plotting frames.
Complex division by I_y fixes a common phase for the signed projection.
No per-edge phase rotation or signed-magnitude substitution is applied. Bars
sum to the projection of the complex total; the omitted quadrature remains in
CSV and is required to assess a full complex zero. Add `--edges-only` to update only
these two PNG/SVG pairs. All figures contain only z=0; z=50/100 nm are retained
in CSV rows with `section_role=reference`, alongside `z_nm` and `z_over_H`.

Mode labels and figure filenames use the user's px/py convention. The primary
integral for py is `integral Ey dS`: `s4_edges_py` is explicitly titled
`py: Ey integral`. Its dimensionless boundary contribution is
`c_e = i*Delta(1/epsilon)/(omega*epsilon0) * nx*integral Hz ds / I_y`.
The Ex figure uses the negative of this coefficient with ny and is only the transverse
reference of this same py mode, titled
`py: Ex reference`; it is not a px-mode result. The Maxwell coefficient first
converts each Hz line integral into E area-integral units (V m); division by I_y
then makes the contribution dimensionless. E maps show `Re(F)` where `F=S*E/I_y`
and S is the full cell area. Thus `integral F dS/S = I/I_y`, equal to 1 for Ey.
Multiplication by S is necessary: E/I_y alone has inverse-area units. Hole outlines and the full
cell boundary mark the integration region: both dielectric and hole air are
included. Below the map, `integral F_y dS/S = integral Ey dS/I_y = 1` gives
the area reference. Below the edge bars, the Maxwell-scaled contour
integral equals `C_boundary = sum(c_e)`, with its full complex numerical value.
The footer compares it directly with the target 1 using `abs(C_boundary-1)`;
the Ex reference instead compares against `I_x/I_y`. Headings and captions
focus on S31; raw dimensional integrals remain in CSV. Positive/negative contributions are
preserved; the boundary sum is not forced to 1, since outer/3D terms and
numerical residuals remain separate. The actual quadrature samples in
`area_field_z0.npz` must reproduce `area_comparison.csv` before plotting;
the separate display raster is never used to calculate the integral.

The separate group figure retains `C_ref=max_edge|integral Hz ds|` normalization
and plots real and imaginary parts separately, with the complex
sum A+B in black. The area figure shows direct area averages and the individual
terms in the boundary reconstruction on separate real/imaginary axes.
The plotted group quantities are boundary integrals,
not the outgoing air amplitude. The plots do not change the numerical or
scientific acceptance decision.

## S31 diagnostic figure

```powershell
uv run python -m scripts.analysis.plot_boundary_analysis <run-directory>/<case-directory> --validation-only
```

This writes only `s4_validation_py.png/.svg`, reading `area_comparison.csv`
and `s4_summary.json`; no COMSOL, field raster or changes to existing data are
needed. The original four figures remain supplementary outputs. The three
equally spaced 4:3 panels show:

1. Direct area average A versus the Maxwell-scaled hole boundary sum B.
2. A and B at the saved integration levels of the same FEM solution.
3. A, B, outer-boundary term O, thickness-derivative term Z, and B+O+Z.

Every panel and integration level uses the same complex finest-level area
average A_f as its divisor. Signed real projections are plotted, while the
reported residuals use complete complex differences: |A-B|/|A_f| and
|A-B-O-Z|/|A_f|. Area change is |A_f-A_previous|/|A_f|. The coarse area
is not separately normalized to one, so changes in sign remain visible.
An exactly zero A_f switches the figure to absolute V/m values. Ratios near
cancellation are sensitive and are not sufficient evidence of accuracy.
Two integration levels do not establish convergence and are not a mesh
convergence study. An exported Z=0 is explicitly marked as unvalidated;
the figure does not infer a zero thickness derivative from symmetry.

## Re-extract sections from an existing Gamma model

```powershell
uv run python -m scripts.analysis.reprocess_boundary_analysis <run-directory>/<case-directory> --replace-results
```

This single-case postprocessor checks the recorded parameter snapshot and the
MPH's Gamma wavevector, lattice, slab height and eigenfrequencies, identifies
the same target at z=0, then reuses the shared S4 integration and plotting code.
It starts COMSOL only to read/evaluate the saved solution: no mesh/eigensolve
and no MPH save. It first stages all new products under `tmp/s4_postprocess_*`.
With `--replace-results`, it backs up the replaced files there and publishes
the S4 CSVs, summary, figures and derived single-case report in place. Without
that flag, it keeps only staged results. Original eigenmode field data and the
MPH stay unchanged. A new `99_config/s4_postprocess.json` records the current
section settings; the original solve manifest remains historical metadata.

For an existing z=0 case missing only the electric-field samples, extract those
alone and then redraw the two three-panel figures:

```powershell
uv run python -m scripts.analysis.reprocess_boundary_analysis <run-directory>/<case-directory> --area-field-only --replace-results
uv run python -m scripts.analysis.plot_boundary_analysis <run-directory>/<case-directory> --edges-only
```

`--area-field-only` preserves all existing integrals, summaries and figures,
checks the sampled integrals against the saved z=0 CSV, and publishes only
`area_field_z0.npz`. This read-only-model operation starts COMSOL for field
evaluation, with no study solve or model save.

## Verification without a new eigensolve

For a saved single-case Gamma solution, audit the independent native FEM
integrals and field derivatives without meshing, solving or saving the MPH:

```powershell
uv run python -m scripts.analysis.audit_s31_boundary <run-directory>/<case-directory> --output tmp/s31_audit/<new-task>
```

The task directory must be new and inside `tmp`. The audit checks the snapshot,
Gamma vector, lattice, height, frequencies, license and running processes. It
records SHA-256 hashes of all existing case files before and after evaluation.
It compares native surface and line integration orders 4/8/12, explicit air
and dielectric boundary traces, translated periodic pairs, local finite
differences, and three Python quadrature resolutions. All use the original
complex field normalization and z=0 as the primary section.

For vector-element magnetic fields, direct `d(ewfd.Hx,z)` can export zero
despite a nonconstant field. The audit independently evaluates
`d(laginterp(k,ewfd.Hx),z)` for k=1/2/3 and compares with finite differences.
The shared S4 extraction now uses k=2 for both Hx and Hy thickness derivatives
and records `derivative_evaluation=elementwise_laginterp2` in its CSV and
summary. Existing saved S4 results are not rewritten by this code change.
Material weights at real boundary faces use explicit domain-sided integrals;
the boundary `dom` entity number is not interpreted as a 3D material ID.
Positive-z cut planes must reproduce both the full and dielectric areas.
The audit does not silently replace any original S4 data or adopt a boundary
trace solely because it makes the equality agree.

`--report-only` regenerates `diagnostic_py.png/.svg`, `report.md`, and the
balance-sensitivity table from a completed task without COMSOL. Its three
equally spaced 4:3 panels show area quadrature sensitivity, signed boundary
trace dependence, and the independently reconstructed thickness term. Every
panel uses the same complex native area reference. Integration stability on
one solution does not establish FEM mesh convergence.

Reviewed audit products can be published alongside the original case as
`s31_*.csv`, `s31_audit.json`, `s31_report.md`, and
`s31_diagnostic_py.png/.svg`; the historical `s4_*` products remain intact.

```powershell
uv run python -m pytest tests/test_s4_boundary.py -q
uv run python -m scripts.run_main.run_boundary_analysis check-model <existing-model.mph>
```

The first command does not import mph. The second checks the installed COMSOL
API against an existing solution. It loads the model read-only, reads complex
fields and derivatives, checks coordinates, then releases it. It neither
solves nor saves the model. Its report explicitly records the source k value.
A saved model from the end of a k scan is not assumed to contain Gamma.
Backend audit files stay under `tmp/s4_backend_check_<timestamp>-<id>/`.
The removed historical pilot directories must not be reused.

API references: [COMSOL Interp](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_results.52.075.html)
and [COMSOL sign convention](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.11.html).
