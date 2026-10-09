# 仿真、参数与结果合同

按任务读取相应节；Unit-cell 2D 的独立合同见 [2D 文档](agent-unit-cell-2d.md)。
模式命名先查 [Memory](../MEMORY.md#scientific-vocabulary-and-invariant-mappings)。本页不会授权启动额外计算。

## 核心模块与数据流

```text
scripts/parameter.json 共享参数 + 各脚本专属配置
  -> geometry/lattice 模块生成孔与边界
  -> simulation runner + mesh builder 构建并求解 COMSOL 模型
  -> 导出频率、Q、场、CSV/parquet/JSON
  -> band/Fourier/far-field 模块分析与绘图
```

主要底层模块：

- `geometry_utils.py`、`hex_lattice_utils.py`、`polygon_utils.py`：几何与裁剪；
- `mesh_constrcution.py`：finite cavity 手动 mesh；
- `simulation_utils.py`：strip/finite 通用 COMSOL runner；
- `simulation_spatial/hexagon_unit_cell.py`：六角 unit-cell runner；
- `simulation_spatial/square_unit_cell.py`、`rectangle_finite_size.py`：保留的可扩展
  本机 runner；
- `basis_utils.py`、`energy_recovery.py`、`band_connector.py`：模式与能带；
- `finite_lattice_fourier.py`、`lattice_fourier_postprocess.py`、
  `farfield_fft.py`：有限结构后处理。

## COMSOL 修改规则

- 优先使用 named/position-based selections，避免依赖固定 entity ID。
- physics、materials 和 mesh 绑定到稳定 named selections。
- 外边界可用小 `Ball` seed，并用 `groupcontang = "on"` 扩展连续边界。
- `SimulationConfig.selection_tolerance` 是 selection seed 和 z-range padding 的共享
  容差；除非有明确原因，不要局部改成不一致的值。
- silicon slab 应构造为 slab layer 减去 hole domains，不要假定固定 domain 编号。
- COMSOL 6.3 中 `Adjacent` 不接受 `exterior = "outside"`；component
  `Difference` selection 使用 `add` 和 `subtract`。

## 共享参数、几何与命名

当前参数格式为 schema v2，完整字段归属见 [参数指南](parameters.md)。
`cells` 存放共享 cavity/cladding 几何，`structure` 存放层数、中心频率及其来源、
shift 列表和逐层 profile；各入口专属参数分属 `unit_cell_band`、`unit_cell_2d`、
`strip_1d` 和 `finite_cavity`，mesh 仍在根级。下文物理规则保持不变；Python 属性
和结果 metadata 沿用原字段名。v1 历史快照继续由统一加载器兼容，不改写已有结果。

- 主入口的 cavity/cladding 紧凑几何与 mesh 只从 `scripts/parameter.json` 读取；
  full finite 与 finite-quarter 还从该文件读取 `finite_cavity.geometry`，只允许精确小写
  `hex` 或 `square`，缺失时回退 `hex`。该选择与 `finite_cavity.cladding_shift_geometry` 独立，
  不影响 unit-cell 或 strip。`square` 系列增加 `_square`，hex 历史目录不变；square
  的 full/quarter 必须复用同一 `FiniteGeometryPlan`，只用完整 cavity parent cells 做
  arbitrary finite DFT，不得调用 hex cyclic quotient。
  strip、finite 与 finite-quarter 的 cavity/cladding 层数和逐层 profile 也只从该文件
  读取。Hex finite 的层数采用 1-based 语义：`cavity_layers=N` 对应 shell `0..N-1`，
  中心 cell 为第 1 层；cladding 对应 shell `N..N+M-1`，其局部 profile layer 仍从 1
  开始。strip 始终扫描 `structure.cladding_shift_factors`。finite 与 finite-quarter 默认也扫描
  该列表，并以 `finite_cavity.cladding_y_over_x_shift_ratio` 派生 `x=factor`、`y=factor*ratio`；
  非默认输入可移除 ratio 并同时提供 `finite_cavity.cladding_x_shift_factor` 与
  `finite_cavity.cladding_y_shift_factor`，此时只运行该显式 pair。ratio 与显式 x/y 不得同时出现，
  x/y 也不得只给一项。`finite_cavity.cladding_shift_geometry` 的 `{"kind":"groups"}` 保留现有
  normal-fan 分组：x 组固定为
  side `{0,2,3,5}`/corner `{0,3}`，y 组固定为 side `{1,4}`/corner
  `{1,2,4,5}`；`{"kind":"ellipse"}` 则让每个 cladding cell 沿指向原点的径向移动，
  幅值为 `A*g(L)*(fx*cos(theta)^2+fy*sin(theta)^2)`，不使用 corner 补偿，并允许
  任一 factor 为 0。geometry 字段缺失时必须回退到 `groups`。finite 与
  finite-quarter 的本征模式求解数统一
  读取 `finite_cavity.eigenmode_count`。逐层 shift 是相对理想晶格的绝对位移，不得按
  shell 累加；finite 与 strip 必须复用 `comsol_workflow` 中的同一 profile 实现。
  `parameter.json` 必须保持为不含注释的标准 JSON；不得写入 `//`、`#`、块注释或
  profile 说明字段。需要恢复历史 uniform cladding 时，显式设置
  `"cladding_shift_profile": {"kind": "uniform"}`，并移除 `scale_layers` 与
  `power`；相关解释写入文档而不是参数文件。
  Unit-cell 不使用共享结构层数和 cladding shift，且 eigensolver 基准始终为
  `c_const/1.55[um]`；完整 unit-cell 主流程成功后默认以 cavity `p2@Gamma`
  更新共享中心频率。扫描、底层复用函数和重绘不得隐式改写该文件。

- Finite 与 strip 的系列目录必须从共享配置派生层数，并分别命名为
  `finite_<cavity_layers>-<cladding_layers>_<cell-and-mesh-label>` 和
  `strip1d_<cavity_layers>-<cladding_layers>_<cell-and-mesh-label>`；非 uniform
  cladding-shift profile 必须增加稳定 `shiftprof-...` 系列后缀。Strip 远端目标保留在
  `shiftX.XXX` case；finite/quarter directional case 必须命名为
  `shiftxX.XXX_shiftyY.YYY`；ellipse finite/quarter 系列必须增加 `_ellipse-shift`
  后缀，且该后缀不得影响 strip；strip cut angle 单独分组为 `cut_60/` 或 `cut_00/`。

- Full finite and quarter finite select one mutually exclusive footprint with
  `finite_cavity.geometry`: missing or `hex` preserves the historical full-cell
  hexagon, while `square` uses the approved left/right half-cell rectangle and
  adds `_square` to the finite series. This setting is independent of
  `finite_cavity.cladding_shift_geometry` and does not affect unit-cell or strip. Classify
  square fragments on the ideal lattice before shifting; keep both fragments
  of a full cladding parent rigid, move only the cladding fragment of a mixed
  parent, and use `square_axis_components_v1` for groups without hex corner
  compensation. Square Fourier uses direct arbitrary finite-point DFT over
  complete cavity parents only, never the hex cyclic quotient.

- `STRIP_CUT_ANGLE` selects a physical strip-boundary orientation, not an
  internal cavity-hole rotation. `0` degrees represents the two horizontal
  edges; `60` degrees represents the four symmetry-equivalent slanted edges.
  Build either case by rigidly transforming the complete cavity/cladding
  lattice, displacement vectors, cell centers, and relevant boundaries into a
  cut-local frame before clipping. It is a strip-only, one-angle-per-run
  control and must not alter finite-cavity or unit-cell workflows.

Directional finite metadata and resume checks preserve both x/y factors and both resolved layer profiles.
For hex `groups`, retain the established 12-region partition; square uses its own documented axis-component rule.

正式 strip 保持 `RUN_DIAGNOSTIC_PLOTS = False`；per-mode 场由 Python `save_field_plot`
生成，不为 `Hz_Re`、`Hz_Im`、`Wem` 创建或运行 COMSOL `PlotGroup2D` / `Image Export`。

## 主入口工作流

保留各自入口的合同；不要把 full/quarter 单模式默认套用到明确选择的 all-modes 专用入口。

- The unit-cell stage exposes `run_main/run_band_pair.py` for the cavity/cladding
  comparison and `run_main/run_band_solo.py` for a cavity-only calculation.
  `run_main/run_strip_1d.py` handles 1D-cut validation, and
  `run_main/run_finite.py` handles the final finite-size cavity calculation. Put
  COMSOL-backed parameter scans in `scripts/run_sweep/`, and analysis or
  plotting scripts that do not launch COMSOL in `scripts/analysis/`. New
  calculation scripts should support one of these three stages and reuse the
  core modules under `comsol_workflow/`.
- The primary runners share the relevant compact geometry and mesh from
  `scripts/parameter.json`. Strip, finite, and quarter finite also consume its
  cavity/cladding layer counts, center frequency, cladding-shift list, and
  shift profile. Full finite and quarter finite share `finite_cavity.eigenmode_count`;
  unit-cell keeps its separate `unit_cell_band.eigenmode_count`. Non-uniform profiles must be encoded in finite/strip series
  identity and resume metadata so they cannot reuse uniform geometry.
  Unit-cell keeps `c_const/1.55[um]` as its solver shift and ignores shared
  structure layer counts and cladding shift; only
  a successful complete unit-cell `main()` may update the shared center
  frequency, using the cavity `p2@Gamma` result.
- A complete pair band calculation compares cavity and cladding cells on the
  same Gamma-M and Gamma-K sampling. A solo band calculation uses the same
  sampling and tracking but reads and solves only the cavity case. Both select
  two p-dominant and two d-dominant Gamma bands and track their identities away
  from Gamma.
- Unit-cell Q-k analysis focuses on the p bands for both cells. Fit Gamma-M and
  Gamma-K separately; present Gamma-M on the negative signed-k half-axis and
  Gamma-K on the positive half-axis.
- For cladding alignment, first align the explicitly translated target p mode
  at Gamma. Treat placement near the cladding midgap as a secondary soft
  objective. Define the Gamma p-d gap from the highest-frequency selected p
  mode and the lowest-frequency selected d mode rather than fixed band labels.
  Tune `b0` first, then `zeta`, and introduce `eta` only if needed.
- `CLADDING_INWARD_SHIFT_FACTORS` is always a non-empty sequence, including for
  a one-point strip calculation; use a one-element sequence rather than a
  scalar. Layer-resolved values are absolute displacements from the ideal
  lattice and must never be accumulated from the previously shifted shell.
  Gamma-trend selection uses the mode minimizing the absolute frequency
  difference from the target, not only modes below the target frequency. By
  default, the strip runner automatically generates the trend after a scan
  containing more than three completed cladding-shift values; an explicit
  per-run `True`/`False` override remains available. Supplemental scans must
  merge all completed `shiftX.XXX` cases in the same cut directory into
  `run_summary.json` instead of replacing the earlier scan summary.
- The finite-cavity workflow should keep simulation, validity filtering,
  center-plane field export, far-field analysis, finite-lattice Fourier
  analysis, and mode scoring as composable stages. Far-field work processes
  only modes marked valid.
- Full finite and finite-quarter default to selecting exactly one valid mode by
  the smallest absolute difference from the configured target frequency before
  expensive downstream analysis. Do not infer the fundamental envelope from
  spatial intensity. Retain candidate raw fields and a case-level frequency
  selection audit, but run far-field, the center-normalized unit-cell map,
  finite-lattice Fourier, and scoring only for the selected mode. Keep this
  behavior optional so a deliberate all-valid-mode analysis remains possible.
- Finite and finite-quarter calculations always use the manual mesh builder;
  `MESH_AUTO_SIZE` remains the single level selector. Its levels map to
  `(maximum element growth rate, curvature factor, narrow-region resolution)`
  as follows: `1=(1.3,0.2,1.0)`, `2=(1.35,0.3,0.85)`,
  `3=(1.4,0.4,0.7)`, `4=(1.45,0.5,0.6)`, `5=(1.5,0.6,0.5)`,
  `6=(1.6,0.7,0.4)`, `7=(1.7,0.8,0.3)`, `8=(1.85,0.9,0.2)`, and
  `9=(2.0,1.0,0.1)`. Set the supported COMSOL size values directly; do not
  reintroduce unsupported COMSOL 6.3 activation properties such as
  `hgradactive`.
- Finite-quarter solves only the first quadrant (`x >= 0`, `y >= 0`) and
  reconstructs the other three quadrants from field parity. The symmetry ID
  mapping along the x/y coordinate-axis boundaries is
  `1={PEC,PMC}`, `2={PMC,PEC}`, `3={PEC,PEC}`, and `4={PMC,PMC}`;
  ID 1 is the default user-targeted fundamental case. Keep
  `finite_quarter.mph`, but do not retain raw quarter-field exports or quarter
  figures once the reconstructed full-field products exist.
- Finite-cavity and strip field images retain center-plane real and imaginary
  Hz with a centered symmetric scale and the `Wave` color table, and
  time-averaged electromagnetic energy density with a linear scale and
  `HeatCamera`. Finite-cavity `11_simulation_exports` use transparent PNG and,
  when an SVG is requested, the compact hybrid form: rasterize only the field
  distribution at 600 DPI while keeping the frame, geometry ticks and labels,
  metadata, and color bar as vectors. Use the `*_hybrid.svg` suffix; do not use
  the very large full-vector field SVG as the routine output. Keep frequency
  (two decimal places), Q, quantity, and unit on one line close to the frame,
  and save PNG and SVG with a tight transparent bounding box and only a small
  safety margin. Include color bars; do not restore the removed
  three-dimensional E/H image exports. Strip runs default to no geometry or
  Fourier diagnostic PNG output; enable such images only through an explicit
  diagnostic switch.
- Far-field sampling uses a physical square enclosing the centered finite
  footprint, represents the origin exactly, and zero-fills points outside the
  selected hexagon or rectangle. Grid and padded FFT sizes remain compatible odd sizes. Read
  the actual xy-air plane position from COMSOL unless an explicit propagation
  distance override is supplied; keep numerical aperture configurable.
- Analysis and plotting stages should be reproducible from saved numerical
  outputs without rerunning a completed eigensolve.

## 输出与布局

- `scripts/.out/` uses four established first-level workflow directories:
  `unit_cell_band/` for unit-cell band/path runs, `unit_cell_2D/` for Cartesian
  unit-cell 2D scans and their local refinements, `strip_1d/` for strip runs and
  analyses, and `finite_cavity/` for full and quarter finite runs. Use
  parameter-bearing, non-colliding directories and derive names and numeric
  precision from the actual configuration rather than manually embedding stale
  labels.

- Group finite and strip results under parameterized series directories named
  `finite_<cavity_layers>-<cladding_layers>_<cell-and-mesh-label>` and
  `strip1d_<cavity_layers>-<cladding_layers>_<cell-and-mesh-label>`.
  Store shift cases below those directories. Group strip cases first by the
  two-digit cut directory `cut_60/` or `cut_00/`, then use shift-only case names
  such as `shift0.050/`. Keep strip trend products inside the cut directory's
  fixed `strip1d_trend/` subdirectory, generated by the deliberately named analysis entry
  `scripts/analysis/srip1d_trend.py`. That directory contains only
  `strip1d_trend.csv` and `strip1d_trend.png`; report selection problems in the
  terminal and never generate an issues TXT. Preserve existing series and case
  directories.

- Name complete pair series `unitcell_band_<cell-and-mesh-label>`. Name a solo
  cavity series `unit_cell_(<b0>-<eta>-<zeta>)_mesh<mesh>`, without `cav` or
  `clad` in the label. Pair runs use `01_results/cavity/k_points/` and
  `01_results/cladding/k_points/` for per-k numerical results, the matching
  `10_overview/<case>/` directories for derived outputs, and root overview
  comparison figures. Solo runs retain only the corresponding cavity paths and
  use `unit_cell_band.png` and `p_bands_q_vs_k.png` instead of comparison
  figures. Store run and per-cell configuration under `99_config/`; pair-only
  validation also belongs there. Do not create empty `00_model` or `80_logs`
  directories when a unit-cell result retained no MPH model or logs.

- Use one stable two-digit output-prefix vocabulary; a prefix has exactly one
  meaning and must not be reused for another function. For finite and
  finite-quarter results, the durable mapping is: `00_model` for the MPH model,
  its `comsol_progress.log`, and the compact shared full-geometry figure set;
  `01_results` for the collection of per-mode
  directories, `10_overview` for CSV and overview-PNG summaries,
  `11_simulation_exports` for per-mode field exports, `12_farfield_FFT` for
  per-mode far-field results, `13_lattice_fourier_Hz` for per-mode
  finite-lattice Fourier results, `80_logs` for other logs, and `99_config` for
  JSON configuration/metadata, including `objective.json`. Keep the series-level `run_summary.json`
  directly in the series root. Use `shiftxX.XXX_shiftyY.YYY` for finite/quarter cases;
  only strip uses `shiftX.XXX`. Name per-mode directories as `mode<index>` without
  zero-padding: `mode0`, `mode1`, ..., `mode9`, `mode10`, `mode11`, etc. A
  per-mode root contains directories only, never files. Keep the parent name
  `01_results`; do not rename it to `modes`, because `modes` and `model` are
  visually easy to confuse.

- Strip cases use the `cut_XX/shiftX.XXX` hierarchy and the same mode-centric vocabulary.
  Keep MPH/progress in `00_model`, aggregates in case-level `10_overview`, JSON
  configuration/objective metadata in `99_config`, per-mode fields in
  `11_simulation_exports`, and per-mode strip quotient/bulk Fourier in
  `13_strip_bulk_fourier_Hz`. Strictly use `gamma_subspace_p > 0.9` as strip
  validity: retain all modes and screening evidence in case-level summaries,
  but keep only passing modes in `01_results`. Each passing mode also contains
  `11_simulation_exports/cutline_Wem.png`, normalized per mode and sampled along
  `x=0` with y on the horizontal axis. Natural mode numbers are never zero-padded, mode
  roots contain directories only, and missing analyses do not get empty
  directories. The strip runner uses case-local `.staging`, preflights unknown
  files and conflicts, rewrites summary paths, and removes staging only after
  success. It does not run full-finite calculations.

- Formal strip calculations use Matplotlib `Agg` and do not open plot windows.
  Keep diagnostic geometry, shift-profile, and Fourier plots disabled; retain
  only the per-mode `Hz_Re`, `Hz_Im`, `Wem`, and `cutline_Wem` result images.
  Generate the first three from already-read field data with Python rather
  than COMSOL `PlotGroup2D` or `Image Export`.

- Full and quarter finite runs keep exactly three model-level geometry figures
  in `00_model`: `full_finite_simulation.png`,
  `full_lattice_modulation_vectors.png`, and one shift-geometry diagnostic:
  `full_lattice_normal_fan_regions.png` for `groups` or
  `full_lattice_ellipse_shift.png` for `ellipse`. The first retains the
  unannotated full simulation view. The vector plot uses `10x` cladding-shift
  arrows, and the shift diagnostics contain no strip-cut guides. Square-cavity
  overview plots must show the actual triangular-lattice hexagonal cells, never
  the rectangular supports used only for half-cell classification; merge rigid
  compatible halves into one hexagon and retain true hexagonal halves for mixed
  or outer parents. Do not generate the annotated simulation or modulated-hole
  overview variants from these runners.

- Keep files directly in the per-mode result directories
  `11_simulation_exports`, `12_farfield_FFT`, and `13_lattice_fourier_Hz`;
  do not add `data`, `metrics`, or `figures` subdirectories, because the extra
  nesting makes result inspection harder. Each directory contains only results
  belonging to that analysis and mode. Apart from the three model-level figures
  assigned to `00_model`, shared Brillouin-zone, basis, sampling-grid, and
  shift-profile information belongs once in the shift-level `10_overview`.
  Do not create an empty analysis directory when an
  analysis produced no result. A result summary such as
  `farfield_summary.json` remains directly in its analysis directory; it is not
  configuration and must not be moved to `99_config`.

- `run_finite.py` and `run_finite_quarter.py` must generate this layout through
  a case-local `.staging` directory so existing post-processors can keep their
  reusable interfaces. Finalization must preflight unknown files and target
  conflicts, rewrite result paths to their final locations, and remove staging
  after success. MPH checkpoints and `comsol_progress.log` bypass staging and
  remain together in `00_model` throughout the run. Default quarter
  `SYMMETRY_IDS=[1]` uses the shift directory directly; multiple IDs or a
  non-default ID use separate symmetry-case directories below the shift.

- Finite-cavity and unit-cell output directory names must distinguish the cell
  parameters and other result-changing controls sufficiently to prevent one
  geometry from overwriting another.

- At launch and completion, report the resolved output directory and provide
  direct paths to the main configuration, logs, machine-readable summaries, and
  figures.

## 显式孔洞 simulation config 几何读取

对于 *_simulation_config.json 这类已经保存显式孔洞多边形的文件，使用
comsol_workflow.simulation_config_geometry.load_simulation_config_geometry() 读取
二维几何。该接口只读取以下字段：

- length_unit，当前必须为 "um"；
- footprint，计算区域的二维凸多边形；
- layers[layer_index].holes，指定层中的孔洞多边形，默认
  layer_index=0。

从仓库根目录调用示例：

```python
from pathlib import Path

from comsol_workflow.simulation_config_geometry import (
    load_simulation_config_geometry,
)

config_path = (
    Path("results")
    / "finite_cavity_reference_results"
    / "finite_cavity_xy_fourier26_bulk5_cladding6_wide_bounds_mma_move0p01_eval0031"
    / "quarter_simulation_config.json"
)
geometry = load_simulation_config_geometry(config_path)

# 传给几何构建器的二维数据，单位为 um。
footprint = geometry.footprint       # shape: (vertex_count, 2)
holes = geometry.holes                # tuple[np.ndarray], each shape: (n, 2)

print(geometry.hole_count)
print(geometry.footprint_bounds)       # (xmin, xmax, ymin, ymax)
```

返回对象 SimulationConfigGeometry 的数组是只读的，并包含
source_path、source_sha256、layer_index、footprint 和 holes，同时提供
hole_count、footprint_width、footprint_height 与 summary() 等便捷属性。
读取器会拒绝非法 JSON、非 um 单位、退化或非逆时针凸多边形、重复顶点，以及超出
footprint 的孔洞。

该接口不读取或覆盖源文件中的 materials、physics、boundary_conditions、
shift_frequency、mode_count 等求解配置。接入 finite/quarter runner 时，仅使用
footprint 和 holes 构建显式横向几何；层数、层厚、材料、mesh、中心频率、本征模式
数、输出目录和后处理开关仍必须按照 scripts/parameter.json 及现有 runner 规则处理。
不得仅凭 config 文件名或其历史求解字段替换 parameter.json，也不得把缺失的
cavity/cladding cell 拓扑元数据静默猜测出来。
