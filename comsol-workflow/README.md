# COMSOL Workflow

A Windows-local Python controller for photonic-crystal unit-cell, strip, and
finite-cavity eigenfrequency simulations with COMSOL 6.3. It builds geometry
and mesh, tracks modes, and exports frequency, Q, fields, Fourier, and far-field
results. No remote server or cluster scheduler is required.

## Development workflow

新增需求先写 [工程 Spec](docs/spec/)，实现、必要验证及按需晋升连续完成：

```text
新可复用算法：tmp/<topic>/<algorithm>.py → comsol_workflow/<owner>.py
新计算脚本：  tests/<topic>/<verb>_<capability>.py → scripts/{run_main,run_sweep,analysis}/
离线验证：    tests/<topic>/checks/test_<capability>.py（随晋升改为验证正式实现）
```

算法与脚本分别负责计算实现和任务编排；试验入口可调用核心包及明确的试验算法，
正式代码不反向依赖 tmp/tests。先复用既有算法，不复制已有实现。
从仓库根按包调用候选，例如 `uv run python -B -m tests.boundary_sampling.analyze_boundary_sampling`；
这是命名示例，不是现有命令，也不授权实际求解。

试验源码保留于主题目录，合成数据/验证缓存进入 tmp/.work 或 tests/.work 的任务目录。
真正 COMSOL 输出仍进入当前已授权结果根。达到来源身份、数值合同、输出和 PNG 检查门槛后，
适合复用的算法合入核心包，入口按求解/扫描/离线分析晋升，更新调用者、验证和文档并消除重复实现。
完整强制流程与命名见 [AGENTS](AGENTS.md#新算法与新脚本的两条晋升通道)。

归档及冗余清理后保留的旧内容位于 [tmp/archive](tmp/archive/) 和 [tests/archive](tests/archive/)，
分类与待决定的晋升项见 [清理报告](docs/staging_cleanup.md)。历史回归仍通过 check_offline.py 显式运行，默认 pytest 不递归 archive。
`test_*.py` 仅用于纯检查，候选 run/export/analyze/plot 脚本必须显式调用。
详见 [算法试验区](tmp/README.md)、[候选脚本与验证区](tests/README.md) 和
[本轮规格及归档回执](docs/spec/plan-execute_20261009_staging_workflow.md)。

工程计划统一命名 `docs/spec/plan-execute_YYYYMMDD_<topic>.md`，能力名用 lower_snake_case；
任务报告写入对应结果的 12_reports。绘图统一执行 [科研绘图合同](docs/agent-plotting.md)，
复用共享路径和图保存函数；图的不同格式共用 stem，PNG/PDF/报告按职责分类。
历史科学命名、准备清单和精确输出合同保留，不自动改写或迁移。

## Quick start

Requirements: Python 3.12+, [uv](https://docs.astral.sh/uv/), COMSOL 6.3, and a
valid local license. From the repository root in PowerShell:

```powershell
uv sync
$env:PATH = "C:\Program Files\COMSOL\COMSOL63\Multiphysics\bin\win64;$env:PATH"
```

Edit [scripts/parameter.json](scripts/parameter.json), then choose a workflow:

| Workflow | Command |
| --- | --- |
| Cavity/cladding bands, mode composition, and Q-k | `uv run comsol-band-pair` |
| Cavity-only bands and Q-k | `uv run comsol-band-solo` |
| Cartesian unit-cell 2D scan with A–D analysis | `uv run comsol-unit-cell-2d` |
| 1D strip and cladding-shift scan | `uv run comsol-strip-1d` |
| Full finite cavity | `uv run comsol-finite` |
| Symmetry-reduced finite cavity | `uv run comsol-finite-quarter` |
| All four quarter symmetry sectors | `uv run comsol-finite-quarter-all --preflight-only` |

These commands start COMSOL calculations. Check the resolved geometry, mesh,
mode count, frequency or k range, output path, and resume settings before
launching. If COMSOL cannot be found, check `where.exe comsol` and `PATH` in
the same PowerShell session. File entry points remain available under
[scripts/run_main](scripts/run_main); registered commands are listed in
[pyproject.toml](pyproject.toml).

## Configuration

Boundary and dipolar programs now use capability names. Historical `s4`/`s5`
modules and the four-modes entry remain aliases to the same implementation.
The integration spec records the complete mapping and preservation checks:
[workflow integration](docs/spec/plan-execute_20261009_workflow_integration.md).

```powershell
uv run comsol-boundary-analysis --help
uv run comsol-boundary-scan --help
uv run comsol-dipolar-analyze broadening-2d --describe-inputs
uv run comsol-dipolar-analyze broadening-2d -- --help
uv run comsol-dipolar-reference --describe-inputs
uv run comsol-dipolar-volume-export --describe-inputs
uv run comsol-dipolar-boundary-replay --describe-inputs
```

These examples only show help or resolved sources. Analysis stages use saved
data; reference, export and replay stages can start COMSOL when their target
options request execution. Each target retains its existing parser and scientific
acceptance checks. `--describe-inputs` opens no scientific data and creates no outputs.

Dipolar commands accept `--source-manifest <JSON>` (or
`COMSOL_WORKFLOW_DIPOLAR_INPUTS`). Schema version 1 resolves relative source
paths against the manifest, with historical locations as the default preset.
For a different series, supply `series_dir`, `mode_dir`, `field_mode_dir`,
`source_config_dir` and `finite_model` together. Other available fields are
`fourier_source`, `reference_dir`, `boundary_case_dir`, `output_dir`, `run_id`,
`request_dir`, `source_batch` and `gamma_group`. Output stays under
`results/S5_dipolarSingularity_analysis`; changing locations does not generalize
the task's accepted geometry, frequency or decomposition contract.

Public dipolar commands classify numerical CSV/NPZ as `01_results`, PNG as
`10_overview`, PDF as `11_pdf`, reports as `12_reports`, checks/logs as `80_logs`,
and inputs as `99_config`. Batch groups stay inside format directories, with
main PNG files flat. Direct historical calls retain their layout; numeric-table
readers support both `01_results` and `80_logs`. Existing finite and boundary
output contracts and previously saved results are retained.

Run the offline regression checks without starting COMSOL:

```powershell
uv run python -B scripts/check_offline.py
```

The runner isolates solver-import tests from scientific tests and rejects real
`mph.start()` calls. It excludes the three real COMSOL example suites.

[scripts/parameter.json](scripts/parameter.json) is the shared, comment-free
JSON configuration (schema v2). See the [parameter guide](docs/parameters.md)
for field ownership, examples, and compatibility. Execution, recovery and
post-processing switches remain in the entry scripts; each result records
its resolved configuration under `99_config/`.

| Section | Parameter source and behavior |
| --- | --- |
| `mesh_size` | Root-level mesh selector for every primary workflow. |
| `cells` | Shared cavity/cladding geometry for band, strip and finite; band solo uses only cavity. |
| `structure` | Layers, center frequency and its provenance, shift factors and layer profile for strip/finite. |
| `unit_cell_band` | Band pair/solo eigenmode count; solver shift remains `c_const/1.55[um]`. |
| `unit_cell_2d` | Independent geometry, mode-pair count, k-grid, tracking, analysis and refinement. |
| `strip_1d` | Strip eigenmode count, cut angle in degrees and kx. |
| `finite_cavity` | Footprint, eigenmode count, directional shifts and ordinary quarter symmetry IDs. |

A successful complete band pair or solo run updates
`structure.center_frequency_thz` from cavity `p2@Gamma`. Sweeps, reusable
functions and redraws do not update the shared file. Legacy schema v1
snapshots remain readable and retain their original layout when updated.

Strip scans `structure.cladding_shift_factors`. Finite/quarter use the same
list with x=factor and y=factor times `finite_cavity.cladding_y_over_x_shift_ratio`;
alternatively, remove that ratio and supply both `cladding_x_shift_factor`
and `cladding_y_shift_factor` inside `finite_cavity`. The inputs are mutually
exclusive. Layer shifts are absolute offsets from the ideal lattice.

`finite_cavity.geometry` selects `hex` (also the missing-field default) or
`square`; `finite_cavity.cladding_shift_geometry` independently selects
`groups` or `ellipse`. These selectors affect finite/quarter only.
Geometry, profile and Fourier rules are documented in
[the simulation reference](docs/agent-simulation.md).

### Cell parameter convention

Compact `(b0, eta, zeta)` inputs map to the existing Fourier geometry:

```text
r_f0        = eta
b_square_f0 = (b0 / 230 nm)^2
b_square_f3 = (zeta^2 - 1) / 2
```

The left/right triangle side is `b0*zeta`; the other four sides are
`b0*sqrt((3-zeta^2)/2)`, preserving fill factor. User-facing p-mode names differ
from stored internal components; see [the mode mapping](CONTEXT.md#language).

## Workflow details

### Unit-cell 2D and refinement

`isQuarter=1` solves an `N x N` grid on `[0,K]^2`, where
`N=q_points_per_axis` and `K=q_max_over_G`. A read-only C2v view produces the
`(2N-1) x (2N-1)` display; mirrored points are never stored as COMSOL results.
`isQuarter=0` independently solves `N x N` points on `[-K,K]^2`.

The main run includes A–D analysis. To redraw an existing complete series:

```powershell
uv run python scripts/analysis/unit_cell_2D.py <series-dir>
```

This command redraws the full analysis; use the shared plotting function when
only one figure is requested. Matching task snapshots, ζ scans, recovery,
completion checks, and local refinement are covered in the
[Unit-cell 2D reference](docs/agent-unit-cell-2d.md).

Valley refinement adds real COMSOL solves. Start with
`comsol-unit-cell-2d-valley-refine ... --prepare-only` and review the sampling
plan before approving a solve. Convergence repairs require their own approved
plan; neither process replaces amplitudes with fitted values.

### Finite cavity and symmetry sectors

Ordinary full/quarter runs enable `FINITE_SINGLE_MODE_ANALYSIS_ENABLED` by
default, selecting the valid mode closest to the target frequency for downstream
analysis. Quarter runs solve the
first quadrant, reconstruct full fields, and reuse the finite analysis.
`SYMMETRY_IDS` selects the boundary pair along the x/y axes:

| ID | x-axis boundary | y-axis boundary |
| --- | --- | --- |
| 1 (default) | PEC | PMC |
| 2 | PMC | PEC |
| 3 | PEC | PEC |
| 4 | PMC | PMC |

For the **complete four-sector spectrum**, use the dedicated batch entry:

```powershell
uv run python scripts/run_main/run_finite_quarter_all_symmetries.py --preflight-only
# After reviewing the preflight:
uv run python scripts/run_main/run_finite_quarter_all_symmetries.py
```

For each shift pair, this batch builds one model and mesh, solves IDs 1–4 at
the shared center frequency, and retains each sector in `finite_quarter_all.mph`
as `solsym1`–`solsym4` / `dsetsym1`–`dsetsym4`. It preserves all eigenpairs and
analyzes every valid mode, independently of the ordinary single-mode switch;
it does not read a unit-cell target table.

The `finite_quarter_all_<cavity>-<cladding>...` series stores the full `4N`
spectrum (`N=finite_cavity.eigenmode_count`) in `10_overview/eigenfrequencies.csv` and
an ungrouped, linear-axis frequency-Q scatter in `10_overview/f_Q.png`.
Mode identities such as `mode7_1_xPEC_yPMC` retain local COMSOL indices and
boundary metadata without using sector identity in scientific plot styling.
This batch requires fresh output; `--output-dir` may select a non-conflicting
series under `scripts/.out/finite_cavity/`.

## Results

Generated data is ignored by Git and belongs to four workflow roots:

```text
scripts/.out/
  unit_cell_band/   # bands and path sweeps
  unit_cell_2D/     # Cartesian scans and refinements
  strip_1d/         # strips and trend analysis
  finite_cavity/    # full, quarter, and four-sector finite calculations
```

Series names derive from actual geometry, layer counts, mesh, and other
result-changing settings. Finite directional cases use
`shiftxX.XXX_shiftyY.YYY`; strip uses `cut_00/` or `cut_60/` followed by
`shiftX.XXX`. Square finite series add `_square`, ellipse shifts add
`_ellipse-shift`, and non-uniform profiles add `shiftprof-...`.

Within a case, `00_model` holds retained models/progress, `01_results` holds
numerical or per-mode data, `10_overview` holds PNG figures (and exact legacy overview CSV contracts), `80_logs`
holds other logs, and `99_config` holds configuration. Optional directories
are created only when needed. Per-mode fields use `11_simulation_exports`,
finite far-fields use `12_farfield_FFT`, and Fourier results use workflow-specific
`13_*` directories. Strip retains only modes with `gamma_subspace_p > 0.9`,
while case-level summaries retain all screening evidence.

Existing results are protected. Resume only where the workflow supports it
and configuration identity matches. A running process or saved MPH alone is
not proof of completion: check exit status, logs, summary, and expected outputs.
See [output contracts](docs/agent-simulation.md#输出与布局) for exact layouts.

## Code and further reading

Promoted reusable geometry, simulation, and analysis live in
[comsol_workflow](comsol_workflow). Orchestration lives in
[scripts/run_main](scripts/run_main), COMSOL-backed scans in
[scripts/run_sweep](scripts/run_sweep), and saved-result analysis/refinement
tools in [scripts/analysis](scripts/analysis). Use package-qualified imports;
importing the package root does not start COMSOL.

- [Script catalog](scripts/README.md#活动-python-程序): cavity/cladding searches, scans, and trend tools.
- [Simulation reference](docs/agent-simulation.md): geometry, mesh, finite/strip processing, and explicit-hole imports.
- [Plotting reference](docs/agent-plotting.md): figure scope, scales, preview approval, and visual checks.
- [Agent instructions](AGENTS.md): repository editing, computation, and result boundaries.

The independent [blueprints](../blueprints/) and
[comsol-data-analysis](../comsol-data-analysis/) projects are sibling directories.
The development environment uses `../blueprints` for the geometry-parity tests;
keep that sibling checkout when installing the `dev` dependency group.

## Validation

Choose tests for the changed workflow. For example, these checks do not launch
COMSOL:

```powershell
uv run python -B -m pytest -q tests/archive/test_run_band_pair.py tests/archive/test_farfield_fft.py
git diff --check
```

For documentation-only changes, check links and the diff. Actual COMSOL
integration checks require the same parameter, license, output, and resource
preflight as a calculation.
