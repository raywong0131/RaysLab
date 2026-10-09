# unit-cell-band v2 COMSOL Model 复用重构

日期：2026-07-29

## 实施前方案概述

### 目标

新增 `scripts/run_main/run_unit_cell_band_v2.py`，与现有
`run_unit_cell_band.py` 并存。在 cavity 或 cladding 的几何参数保持不变时，
一个 case 只创建一个 COMSOL client 中的 Model、只构建一次几何，并在全部 k 点
之间复用；每个新 k 点只更新 Floquet `kx/ky`、清除旧 solution data 并重新执行
eigenfrequency solver。mesh 设置变化时继续复用同一个 Model 和几何，仅更新并
重跑 mesh 后再求解。

重构不得改变现有程序的几何、物理场、k 路径、缓存、模式有效性、模式成分、
Gamma 选模、band tracking、Q 拟合、绘图、输出数据字段和中心频率更新功能。

### 共存与兼容策略

- 保留 `scripts/run_main/run_unit_cell_band.py` 及现有
  `comsol-unit-cell-band` CLI，不改变 v1 的逐 k 点新建 Model 行为。
- 新增 `comsol-unit-cell-band-v2` CLI，指向新入口。
- v2 使用独立系列目录
  `scripts/.out/unit_cell_band/unitcell_band_v2_<cell-parameters>_mesh<...>/`，
  不读取或覆盖 v1 的正式结果目录。
- v2 不导入 `run_unit_cell_band.py`，以独立源码保存当前完整计算和后处理流程；
  两个入口互不调用，后续可以分别验证、修改和回滚。
- 不修改现有 `run_unit_cell_band.py`，也不修改其当前使用的
  `comsol_workflow.simulation_spatial.hexagon_unit_cell.SimulationRun`。v2 在自己的
  模块内定义复用 runner：首次调用委托现有底层完成完整构建，后续同几何调用由
  v2 自己更新 k、按需重跑 mesh 并求解。几何签名变化时明确拒绝，要求创建新
  runner。

### COMSOL 生命周期

```text
一个 case（cavity 或 cladding）
  -> 创建一次 v2 ReusableSimulationRun / COMSOL Model
  -> 第一个未缓存 k：构建几何、named selections、材料、physics、study、mesh、solve
  -> 后续未缓存 k：更新 kx/ky -> clearSolutionData -> runAll
  -> mesh size 改变：更新 autoMeshSize -> mesh.run -> clearSolutionData -> runAll
  -> 已缓存 k：直接读 CSV/parquet，不触发 COMSOL solve
  -> case 完成后移除 Model
```

cavity 和 cladding 几何不同，因此各自使用一个 Model；两者仍位于同一 Python
进程和同一 `mph` client/session 中串行执行。

### 几何与 mesh 身份

- 以晶格常数 `a` 和全部 hole vertices 的规范化数值元组作为几何签名。
- 同一 runner 后续调用必须具有完全相同的几何签名；不允许在已有 Model 上静默
  替换几何。
- 记录当前已经应用的 `mesh_auto_size`。值未变化时不重跑 mesh；值变化时只调用
  mesh size 更新和 `mesh1.run()`，不重建 geometry、materials、physics 或 study。
- runner 暴露只读统计：Model 构建次数、mesh 构建次数、solve 次数和当前 k，
  用于测试、metadata 和运行诊断。

### v2 数据流与功能保持

1. 与 v1 一样在模块导入时读取 `parameter.json`，生成 cavity/cladding compact
   geometry 和非冲突输出标签。
2. 与 v1 一样生成 Gamma-M、Gamma-K 及 Q-fit dense k 点，并去重。
3. 每个 case 先检查几何间距并输出 geometry 图。
4. 对 k 点逐项检查现有 v2 缓存；首个 cache miss 时创建 runner，以后复用到该
   case 结束。若全命中，则不启动 COMSOL。
5. 每个新解继续执行与 v1 相同、但保存在 v2 文件内的有效模式检查、复数中心面
   Hz parquet 和模式成分流程。
6. 完成 Gamma `p1/p2/d1/d2` 识别、两个方向 tracking、mode-composition 汇总、
   unmatched 校验、p-band Q 拟合和对比图。
7. 完整成功后，与 v1 一样以 cavity `p2@Gamma` 更新共享中心频率；参数身份在
   写入前继续由 `parameter_config` 校验。

### 输出与 metadata

- 保持 `01_results/`、`10_overview/`、`99_config/` 的既有布局和所有 CSV/PNG。
- v2 根配置增加 `execution_engine: "reused_model_v2"`。
- 每个 case 配置增加 Model/mesh/solve 统计、cache hit/miss 数，以及几何复用说明。
- 原子写入继续复用 v1 的 CSV/JSON helpers。

### 验证计划

#### 静态与单元测试

- 新 v2 模块可导入且具有独立 OUT_DIR/CLI。
- 同一 case 多个 cache miss 只创建一个 `SimulationRun`，并对每个 miss 调用一次
  solve；cache hit 不调用 COMSOL。
- cavity 与 cladding 各自创建独立 runner。
- 全 cache hit 时不创建 runner。
- v2 自身的 Gamma 选模、tracking、Q 拟合、绘图和中心频率更新保持可用。
- v2 runner 的几何签名稳定；不同几何被拒绝；现有底层 runner 文件无差异。
- mesh size 未变化时不重跑 mesh；变化时 Model build count 保持 1、mesh build
  count 增加、solve count 正常增加。
- 运行 `py_compile`、聚焦 pytest、相关 v1 回归和 `git diff --check`。

#### 本机 COMSOL 聚焦验证

启动前再次核对：固定 unit-cell 几何、mesh5、固定 solver shift、少量 k 点、独立
任务级输出/无正式结果覆盖。实机验证至少包括：

1. 同一 runner 连续求解 Gamma 与一个非零 k，确认 Model build count 为 1、mesh
   build count 为 1、solve count 为 2；
2. 同一 runner 将 mesh5 改为 mesh6，确认只重跑 mesh，Model build count 仍为 1；
3. 用一个全新 Model 单独求解相同非零 k，与复用 Model 的频率/Q 数值比较；
4. 确认 repeated result evaluation、场读取和 mode index 均来自最新 solution，而非
   前一个 k 的残留结果。

### 风险与防护

- COMSOL solver sequence 可能缓存旧参数：每次复用求解前显式调用
  `sol1.clearSolutionData()`，再 `runAll()`；结果读取前重新运行 numerical feature。
- mesh 改变会使旧 solution 失效：先重跑 mesh，再清 solution 并求解。
- 结果 dataset/table 可能保留旧 plot data：v2 不依赖 table 累积内容，
  `get_eigenfrequencies()` 和 field interpolation 都重新运行当前 numerical feature；
  实机测试验证这一点。
- k 点中途失败时保持已有 per-k 缓存；runner 在 context exit 时移除 Model，后续可
  从已完成 k 点恢复。
- v2 尚未验证通过前不替换 v1，不修改现有 CLI 指向。

### 回滚方案

- 删除 v2 入口、v2 CLI 和聚焦测试即可恢复；现有 v1 与底层 runner 从始至终不在
  本次任务中修改。
- v2 结果位于独立目录，可单独保留或清理，不影响 v1 正式结果。
- 不修改或删除任何现有计算结果，也不自动提交或 push。

## 实施记录与结果

### 实际改动

- 新增 `scripts/run_main/run_unit_cell_band_v2.py`。该文件不导入 v1，完整保存几何、
  有效模式、场、模式成分、Gamma 选模、band tracking、Q 拟合、绘图和中心频率
  更新流程。
- v2 内新增 `ReusableSimulationRun`，以 `a + hole vertices` 作为几何签名。第一次
  调用沿用现有底层 runner 完整构建；后续同几何只更新 `kx/ky`、调用
  `clearSolutionData()` 和 `runAll()`。mesh size 改变时额外执行
  `autoMeshSize()` 与 `mesh1.run()`，不重建 Model/geometry/physics/study。
- v2 的每个 cavity/cladding case 延迟创建一个 runner：全 cache hit 时不启动
  COMSOL；存在 miss 时从第一个 miss 复用到 case 结束，并在结束或异常时移除
  Model。
- v2 使用独立 `unitcell_band_v2_<parameter-label>` 输出目录；根配置和 case 配置
  记录执行引擎、Model/mesh/solve 次数与 cache hit/miss。
- 新增 `comsol-unit-cell-band-v2` CLI，并在 `scripts/README.md` 记录共存方式。
- 新增 `tests/test_run_unit_cell_band_v2.py`，覆盖独立入口、几何保护、同 mesh
  复用、mesh-only 重构、调用方持有 runner、每 case 至多一个 runner 和全缓存不
  创建 runner。

### 原程序恢复记录

实施早期曾短暂尝试扩展现有底层 runner，并抽取 v1 的结果持久化函数；用户明确
要求完全独立后，这些改动已逐项撤回。现有
`comsol_workflow/simulation_spatial/hexagon_unit_cell.py` 无本任务内容差异，v1
也未保留复用逻辑。排查原程序运行异常时确认 `EIGENMODE_COUNT=4` 导致 Gamma-M
末段有效模式不足；按用户确认恢复为仓库基线的 12。现有 4-mode 正式缓存未删除、
移动或覆盖，并且原程序当前仍会把它们识别为完整缓存；后续若继续该旧目录，必须
由用户明确选择清理/覆盖或改用新目录重算。

### 验证结果

- `python -m py_compile scripts/run_main/run_unit_cell_band.py
  scripts/run_main/run_unit_cell_band_v2.py`：通过。
- `pytest -q tests/test_run_unit_cell_band_v2.py`：`7 passed`。
- v1 聚焦回归：`85 passed, 1 deselected`。排除项仍硬编码期望 mesh9，而当前共享
  配置是 mesh5；不是 v1/v2 执行逻辑失败。
- `git diff --check`：通过。
- `uv sync`：通过；本地 console scripts 同时包含原
  `comsol-unit-cell-band` 与新 `comsol-unit-cell-band-v2`。
- 本机 COMSOL 6.3 聚焦验证（cavity、12 requested modes）：
  - Gamma 后复用求解 `k=(0.002, 0)`：Model/mesh/solve 为 `1/1/2`；
  - 同一 Model 改 mesh5 为 mesh6 并求解 `k=(0.004, 0)`：`1/2/3`；
  - 全新 Model、mesh5 单独求解 `k=(0.002, 0)`：`1/1/1`；
  - 复用与全新 Model 均返回 24 个结果，最大频率差
    `4.2995e-7 THz`，Q 最大相对差 `2.1767e-5`；
  - Gamma、复用非零 k、mesh6 以及 fresh-model 四次场读取均成功且各有 3120 个
    中心面采样点。

### 遗留问题与最终状态

- v2 已实现并通过静态、单元、v1 回归和本机 COMSOL 数值验证，可以与 v1 共存。
- 尚未替换或删除 v1；是否正式替换需用户在 v2 完整工作流结果验收后另行确认。
- 本任务生成的 pytest/COMSOL 临时验证目录均已清理，没有修改既有正式结果。
