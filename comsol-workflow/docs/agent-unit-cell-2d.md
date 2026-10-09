# Unit-cell 2D 操作与科学合同

仅在运行、修改或重绘 Unit-cell 2D、ζ 扫描和 valley refinement 时读取。通用授权见 [AGENTS](../AGENTS.md)。

### 参数权威源与网格语义

- 主计算只从当前参数文件的 `unit_cell_2d` 对象读取紧凑几何 `b0_nm/eta/zeta`、
  `eigenmode_pair_count`、`q_max_over_G`、`q_points_per_axis`、`isQuarter`、
  `target_bands`、模式跟踪阈值和分析参数；mesh 单独读取根级 `mesh_size`。不要误用
  `cells.cavity` 几何或 `unit_cell_band.eigenmode_count`。默认参数文件为
  `scripts/parameter.json`，任务级快照可通过环境变量
  `COMSOL_WORKFLOW_PARAMETER_PATH` 指定。
- `q_axis_over_G` 已废弃且会被配置加载器拒绝。只给 `q_max_over_G`（$K$）和
  `q_points_per_axis`（单轴点数 $N$）；两个方向强制使用同一范围、点数和间隔。
- `isQuarter=1` 时，COMSOL 真正求解第一象限
  $(k_x/G,k_y/G)\in[0,K]^2$ 的 $N\times N$ 点，单轴间隔为 $K/(N-1)$；绘图时通过
  只读 C2v 访问器覆盖 $(2N-1)\times(2N-1)$ 坐标。镜像坐标不得写入 CSV、Parquet、
  point cache 或 MPH 数据目录。
- `isQuarter=0` 时，COMSOL 独立求解 $[-K,K]^2$ 的 $N\times N$ 点，$N$ 必须为不小于
  3 的奇数，单轴间隔为 $2K/(N-1)$。Quarter 的 $N$ 可为任意不小于 2 的整数。
- `eigenmode_pair_count=M` 表示每个 k 点请求 $M$ 对、预期 $2M$ 个本征值；它与
  unit-cell band 使用的 `unit_cell_band.eigenmode_count` 无关。求解 shift 固定为
  `c_const/1.55[um]`。
- 当前只接受 `normalization_kind="planar_l2"`、`field_plane="air"`，且
  `analysis_enabled` 必须保持 `true`。正式名称由内部 `p1/p2` 通过当前映射显示为
  `px/py`；设置和 cache 中继续保留内部 band 名，不手工改写存储字段。

### 昂贵计算启动前检查

启动前必须从实际参数文件重新解析并报告：几何三元组、mesh、Quarter/full、$K$、$N$、
真实 COMSOL 点数、绘图坐标数、目标 band、每点预期本征值数、输出 series 路径和是否运行
默认 analysis。还要检查当前 Python/COMSOL 进程、COMSOL 6.3 PATH、license 和目标目录
的 `99_config/config.json`；不得仅凭目录名或对话记忆启动。

series 名由实际参数自动生成，形式为：

```text
unit_cell_2D_(b0-eta-zeta)_meshM_kmaxK_uniformN_quarter|full_band-<formal-bands>
```

同名目录已有 `config.json` 时，程序只允许完全相同的 cache identity；身份不同会拒绝
运行，不得为了绕过检查覆盖配置或混用旧 point cache。

### 主计算、恢复与默认分析

在仓库根目录、已配置 COMSOL PATH 的同一 PowerShell 会话中运行：

```powershell
uv run comsol-unit-cell-2d
```

等价的文件入口为：

```powershell
uv run python scripts\run_main\run_unit_cell_2d.py
```

主入口按从 Gamma 向外的顺序跟踪目标 band，在同一个 reusable COMSOL model 中逐点更新
k、重求并保存结果；首轮结束后仅对 missing/ambiguous 点做一次基于真实重求的确定性修复，
绝不插值生成 c(k)。中断后以完全相同参数重跑同一命令，会逐点复用通过 identity 和文件
完整性检查的 cache。不要另写恢复脚本，也不要删除已完成点。

正式运行默认在求解后自动调用 `run_analysis()`，生成标准 A--D 的 $|c(k)|$、分量、偏振
椭圆/相位与 winding 输出。`--skip-analysis` 仅用于诊断：

```powershell
uv run comsol-unit-cell-2d --skip-analysis
```

此时成功求解的 `run_summary.json` 状态保持 `solved` 且
`analysis_completed=false`，不能当作正式完整结果。参数中的 `analysis_enabled=false` 不是
替代方式，配置加载器会直接拒绝。

如需使用任务级参数快照，不修改共享 `parameter.json`：

```powershell
$env:COMSOL_WORKFLOW_PARAMETER_PATH = "I:\path\to\task_parameter.json"
uv run comsol-unit-cell-2d
Remove-Item Env:COMSOL_WORKFLOW_PARAMETER_PATH
```

### 不启动 COMSOL 的标准重绘

已有完整 `polarization_grid.csv`/Parquet 时，使用分析入口而不是重跑求解：

```powershell
uv run python scripts\analysis\unit_cell_2D.py `
  scripts\.out\unit_cell_2D\<series-dir>
```

省略 `<series-dir>` 时，入口按当前参数文件推导默认 series。显式给定历史 series 时，
Quarter/full 解释以该结果的 `99_config/config.json` 为准；`angle_methods` 和
`mask_relative_threshold` 仍来自当前参数文件，重绘前必须核对。只要求某一张图时，应调用
共享绘图函数窄范围原位重绘，不得用完整 analysis 顺带覆盖其余图。

`unit_cell_2D.py` 不再提供 minimum-curvature、拟合 floor 或 topology reconstruction 的
历史 CLI；已有比较文件只作为历史记录保留，不得重建为标准 B4 或物理结果。真实 valley
改善只允许走 COMSOL 局部加密流程。

### ζ 参数扫描

ζ 扫描串行调用主入口，每个 case 写独立参数快照并沿用相同的恢复规则。始终显式给出本次
ζ 列表，先用 `--prepare-only` 审核目录和快照：

```powershell
uv run comsol-unit-cell-2d-zeta-scan --zeta 1.155 1.156 1.157 --prepare-only
uv run comsol-unit-cell-2d-zeta-scan --zeta 1.155 1.156 1.157
```

其他参数全部来自 `scripts/parameter.json`，也可增加
`--parameter-file <task-parameter.json>`。扫描 wrapper 位于 `unit_cell_2D_zeta_scan_...`，
实际 case 仍分别写入各自 `unit_cell_2D_(...)` series；任一 case 失败时串行扫描停止并在
`99_config/scan_summary.json` 保留状态。

### COMSOL valley 局部加密

局部加密必须消费已完成的标准 Cartesian series，并使用与 source 完全匹配的参数快照。
第一步永远只生成采样计划，不导入 `mph`：

```powershell
uv run comsol-unit-cell-2d-valley-refine `
  <source-series-dir> --parameter-file <matching-parameter-snapshot> `
  --prepare-only
```

检查 sampling-plan PNG、去重后的新增点数、法向/纵向覆盖和 `maximum_new_points` 后，只有
用户明确批准才可去掉 `--prepare-only` 启动 COMSOL。默认恢复既有 refined point cache；
`--no-resume` 会拒绝已有 cache，而不是清空它。检测到其他项目 COMSOL 进程时默认拒绝，
只有用户明确允许共享资源才可加 `--allow-concurrent-comsol`。完成后只重绘 refined B4：

```powershell
uv run comsol-unit-cell-2d-valley-refine `
  <source-series-dir> --parameter-file <matching-parameter-snapshot> --plot-only
```

若已完成结果存在数值不收敛点，先生成独立修复计划并审核，不得直接追加求解：

```powershell
uv run comsol-unit-cell-2d-valley-refine `
  <source-series-dir> --parameter-file <matching-parameter-snapshot> `
  --prepare-repair-only
```

只有用户批准该修复计划后，才运行完全相同 source 和参数快照的
`--run-approved-repair`；该步骤只求解已审核、去重后的修复点，并保留现有 refined cache：

```powershell
uv run comsol-unit-cell-2d-valley-refine `
  <source-series-dir> --parameter-file <matching-parameter-snapshot> `
  --run-approved-repair
```

Quarter source 的候选坐标在交给 COMSOL 前规范到第一象限；其他象限只在最终 Delaunay
线性显示时临时镜像。不得用拟合支线、minimum-curvature、人工 floor 或中心恢复替代真实
COMSOL 振幅。

### 输出布局与完成判据

每个标准 series 固定使用：

```text
scripts/.out/unit_cell_2D/<series>/
  00_model/       geometry.png, unit_cell_2D.mph
  01_results/     k_points/ 与 polarization_grid.parquet
  10_overview/    标准 A--D PNG
  80_logs/        polarization_grid.csv, band_tracking_2d.csv, winding scan CSV
  99_config/      config.json, run_summary.json, winding/symmetry summaries
```

只有同时满足以下条件才报告正式完成：进程已退出且无未解释错误；
`99_config/run_summary.json` 为 `status="complete"`、`analysis_completed=true`；真实点数与
$N^2$ 一致；源表行数等于真实点数乘目标 band 数；`match_status_counts` 中 missing/ambiguous
已被理解；CSV、Parquet、MPH、标准 A--D PNG 和 analysis summary 均存在。Quarter 的
`plot_coordinate_count` 只是只读显示坐标数，不是 COMSOL 数据点数。主计算、ζ 扫描、
底层重绘和 refinement 都不得隐式改写 `scripts/parameter.json`。

## 已批准的 A–D 专用显示与 refinement 合同

以下精确数值仅用于其明确指定的既有输出；修改时核对当前实现、测试和本次要求，不推广为其他图默认值。

- Unit-cell 2D B4 follows exactly the same component-map contract as B1--B3:
  use the saved raw `cx_raw/cy_raw` coefficients, divide both component
  magnitudes by the single full-grid `max sqrt(|cx_raw|^2+|cy_raw|^2)`, densify
  with the shared `_smooth_grid()` cubic display interpolation, and render B4
  with `magma` plus the shared logarithmic colorbar rules. B4 must not have a
  separate coherency-matrix, Gaussian, nodal-line, curve-fitting, center-pixel,
  or smoothing-only path. A retained `*.pre_regenerate_*.png` is a protected
  historical reference and must not be overwritten by a standard redraw.

- Quarter-series B4 keeps that data and interpolation contract but uses the
  fixed logarithmic colorbar interval `10^-4` through `10^-1`, with both
  endpoints labeled as powers of ten. Full-series B4 retains its data-derived
  logarithmic interval.

- Quarter-series B2 likewise keeps the shared component-map data contract but
  uses the fixed logarithmic colorbar interval `10^-4` through `10^0`, with
  both endpoints labeled as powers of ten. Full-series B2 retains its
  data-derived logarithmic interval.

- Quarter-series B1 keeps the shared component-map data contract but uses the
  fixed linear colorbar interval `0` through `1`. Full-series B1 retains its
  data-derived linear interval.

- Quarter-series A1 and A2 display normalized
  `|c(k)| = sqrt(|cx(k)|^2 + |cy(k)|^2)` rather than normalized intensity.
  Their retained filenames remain `A1_ck_intensity_linear.png` and
  `A2_ck_intensity_log.png` for output compatibility. Full-series A1/A2 retain
  the historical normalized-intensity definition.

- Quarter-series D contains only one row: the unscaled winding-loop
  polarization ellipses and its winding-number curve. Omit both historical
  `cy * 10` panels. Full-series D retains the two-row compatibility layout.

- Quarter-series C2 reconstructs the pi-periodic polarization director angle
  from first-quadrant Stokes data. Interior points change director-angle sign
  under one mirror and recover it under two mirrors; the negative halves of
  `kx=0` and `ky=0` retain the same director branch as their solved positive
  half-axes before doubled-angle interpolation. This C2-only display convention
  must not alter Jones vectors, C1, D, winding, or stored source rows.

- Quarter-series C1 samples `N x N` evenly spaced display positions from the
  transient `(2N-1) x (2N-1)` C2v view, where `N=q_points_per_axis`; it does not
  plot all mirrored positions or persist sampled/mirrored rows. Each position
  retains its arrowed Jones ellipse. Use red for positive
  `S3=-2 Im(cx*conj(cy))` (counter-clockwise), blue for negative S3 (clockwise),
  and neutral grey only for exact zero. Full-series C1 remains unchanged.

- Every standard quarter-series kx-ky panel in A--D explicitly labels both
  sampled coordinate endpoints and uses `0.05 G` major-tick spacing. This
  applies to A1/A2, B1--B4, C1/C2, and the left winding-loop panel in D; it
  does not alter non-k-space axes. Full-series tick placement remains
  data-derived.

- Quarter-series kx-ky tick labels omit insignificant trailing zeros on both
  positive and negative values: use `0`, `0.1`, `-0.1`, `0.05`, and `-0.05`
  rather than `0.00`, `0.10`, or `-0.10`.

- Unit-cell 2D supports explicit `isQuarter=0/1`. In quarter mode,
  `q_points_per_axis=N` means an `N x N` real COMSOL grid on
  `kx>=0, ky>=0`; C2v is applied only by a read-only plotting/winding accessor
  to cover a `(2N-1) x (2N-1)` display. Never materialize mirror points in the
  polarization table, point cache, parquet, or COMSOL directories. Missing
  `isQuarter` in a historical snapshot means legacy full scan. New series names
  carry `_quarter` or `_full`, and local refinement canonicalizes all requested
  COMSOL coordinates to the first quadrant when its source is quarter data.

- Minimum-curvature, fitted-floor, topology-adaptive branch reconstruction, and
  B5 curve-overlay routes are rejected as physical-result generators. Retain
  their existing files only as historical comparisons; do not extend or use
  them to replace standard B4.

- The approved physical refinement route uses the completed Cartesian scan only
  to choose additional k coordinates. Detect every strict row/column local
  minimum without an amplitude threshold; permit a three-point parabola only
  for coordinate refinement; connect adjacent candidates by mutual-nearest
  piecewise-linear segments; and use a real two-dimensional COMSOL stencil when
  slice counts or directions are ambiguous. Every added amplitude must come
  from the same-mesh COMSOL eigensolve with Hz-overlap mode tracking. Later
  rounds and recentering are selected only from actual computed minima. Merge
  the saved coarse/refined raw `cx_raw/cy_raw`, divide both components by one
  `max sqrt(|cx_raw|^2+|cy_raw|^2)` over the full real-node set, and use only
  Delaunay barycentric linear interpolation for display. Planar-L2 fields may
  remain as diagnostics and mode-tracking metadata but do not define refined B4.
  Never add fitted branches, topology constraints, floors, Gaussian/minimum-
  curvature fields, center restoration, or amplitude corrections. Prepare-only
  sampling plans require user review before COMSOL starts, and the new result
  never overwrites standard or historical B4.

- A convergence-repair point CSV is not itself a complete audit manifest:
  logical windows whose required coordinates already exist may add zero rows.
  Persist and hash a separate target manifest covering every normal/junction
  window, and require every target to pass completeness and interior-minimum
  gates before declaring the refined result complete.

- Sparse longitudinal anchor plans also audit every extracted piecewise
  segment for uncovered midpoint gaps. A gap-fill coordinate may be added only
  when it is farther than the configured longitudinal spacing from planned
  nodes and differs from both coarse and new coordinates at solver precision.
  Gap-fill points never supply fitted amplitudes or define an adaptive
  minimum-search stencil.
