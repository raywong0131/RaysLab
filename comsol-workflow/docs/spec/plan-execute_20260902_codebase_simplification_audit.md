# 当前 COMSOL workflow 代码库简化审核与执行方案

日期：2026-09-02
状态：审核完成；Batch 1--2 已实施并验证，候选 3--6 尚未实施
方法：`simplify-codebase` Broad audit

## 一、目标

本方案的目标是减少当前 COMSOL 本机控制器需要长期保持一致的概念、旧入口、重复实现和
无主扩展面，而不是为了缩短文件或统一代码风格。优先选择已经被现行合同取代、在生产入口
不可达、并且具有明确验证边界的纵向删除。

审计阶段只形成 spec；实施阶段当前仅完成候选 1。全过程不移动或删除正式计算结果、
不启动 COMSOL。后续候选仍须按本文的独立批次逐项批准、逐项验证，并继续在本文回填实际改动和测试结果。

## 二、范围、排除项与基线

### 2.1 纳入范围

- 根级项目合同：`AGENTS.md`、`README.md`、`pyproject.toml`、`.gitignore`；
- `comsol_workflow/` 的几何、数值、COMSOL runner、后处理模块；
- `scripts/run_main/`、`scripts/run_sweep/`、`scripts/analysis/`；
- `tests/`；
- `docs/` 和 `docs/spec/` 中与当前 COMSOL workflow 有关的记录；
- `scripts/.out/` 只读取目录、JSON 和文件名以识别持久化兼容合同，不读取或改写大型场、
  MPH、parquet、PNG 等正式数据。

### 2.2 明确排除

- `blueprints/` 整个目录；
- `comsol-data-analysis/` 整个目录；
- `.git/`、`.venv/`、`.understand-anything/` 和测试/工具缓存内容；
- 与被排除目录绑定的实现选择和行为判断。`tests/test_blueprint_comsol_geometry_parity.py`
  虽位于根 `tests/`，但跨越被排除的 blueprint 边界，本轮不据此提出接口删除。

### 2.3 工作树基线

审核时分支为 `raywong/finite-cavity`，跟踪
`origin/raywong/finite-cavity`。工作树在本轮开始前已有用户修改，包括
`AGENTS.md`、`MEMORY.md`、`scripts/parameter.json`、删除的 `scripts/__init__.py`、
两个排除目录及若干未跟踪 handoff/临时文件。本轮不接管、不改写这些内容。

`.codex_process_probe.py`、`.task-mplconfig/`、`.vscode/` 等未跟踪项属于工作区状态，
不是可自动删除的产品死代码；需由原作者确认用途后另行清理。

## 三、当前合同覆盖图

| 域 | 已检查的当前所有者与入口 | 结论或盲点 |
| --- | --- | --- |
| 安装命令 | `pyproject.toml:[project.scripts]` 的 8 个命令 | 是安装后 CLI 的当前所有者；README 未完整列出全部命令 |
| 直接主入口 | band pair/solo、unit-cell 2D、strip、finite、finite-quarter、four-modes | `run_finite_quarter_four_modes.py` 是有文档和测试的直接入口，不能因未注册 console script 判死 |
| 扫描与分析 | `scripts/run_sweep/*.py`、`scripts/analysis/*.py` | 多数是有文档的直接文件入口；一次性迁移/恢复脚本需单独分类 |
| 参数 | `scripts/parameter.json` + `scripts/run_main/parameter_config.py` | 共享参数、环境变量快照和历史参数解释是持久化边界，默认保留 |
| Unit-cell 求解 | `run_band_pair.py`、`run_band_solo.py`、`hexagon_unit_cell.py` | solo 是复用 pair 实现的薄入口，边界合理 |
| Unit-cell 2D | `run_unit_cell_2d.py`、`unit_cell_2d_analysis.py`、valley refinement | 标准 A--D 与真实 COMSOL refinement 活跃；三条拟合重建路线已被正式否决但仍可调用 |
| Strip | `run_strip_1d.py` + layout v3 staging/finalize | 当前本地发现的 strip config/run summary 均为 layout v3；脚本仍含被入口明确拒绝的 full-finite 和永久禁用诊断图 |
| Full finite | `run_finite.py` 当前 shift -> `run_finite_full_case()` 流程 | 活跃路径使用共享 `lattice_fourier_postprocess.py`；同文件仍保留旧 gap、旧 unit-cell 和旧本地 Fourier 流水线 |
| Quarter finite | `run_finite_quarter.py` 复用 `run_finite.py` 的几何、输出和后处理 | full/quarter 分离承载不同边界条件、场重建和生命周期，不能仅因相似而合并 |
| COMSOL 生命周期 | `SimulationRun`、checkpoint MPH、resume metadata、staging/finalize | 进程、恢复、原子发布和数据防丢是强边界；简化不得削弱 |
| 持久化结果 | 四个 `scripts/.out/` 一级目录及 case-local config | 存在 legacy full unit-cell 2D、scalar/directional finite 等历史结果；读取兼容仍有真实消费者 |
| 底层 package | `comsol_workflow/` 21 个活动/保留模块 | 三个模块只有测试/README 消费，且 README 明示“为未来保留” |
| 测试 | 根 `tests/` | 大量测试是活跃合同；少数测试只维持不可达或已否决能力 |
| 文档/决策 | AGENTS、三层 README、约 40 个历史 plan/spec、handoff | 历史 spec 保存设计原因；当前入口清单和数量出现可证实漂移 |
| 外部消费者 | setuptools 会发布 `comsol_workflow*` 和 `scripts*` | 仓库外直接 import 无法从本仓库证明不存在；删除公开模块时需降低置信度 |

## 四、候选排序摘要

收益和置信度分别排序。预计行数仅用于显示维护负担，不是批准删除的依据。

| 顺序 | 候选 | 置信度 | 预期净效果 | 风险 |
| ---: | --- | --- | --- | --- |
| 1 | 删除 `run_finite.py` 中被取代的旧 finite/unit-cell/Fourier 流水线 | 高 | 约 900 行生产代码；移除三套旧输出/模式概念 | 中低 |
| 2 | 删除已正式否决的 unit-cell 2D 拟合重建能力 | 高（仓库内） | 约 2,000 行生产代码、约 400 行专属测试、3 个 CLI 开关 | 中 |
| 3 | 删除 strip 中不可达 full-finite 与永久禁用诊断图分支 | 高 | 约 600--800 行生产代码及专属测试 | 低至中 |
| 4 | 收敛入口与运行文档的多重事实源 | 高 | 修复已发生的漂移，减少每次变更的同步面 | 低 |
| 5 | 在外部历史数据盘点后删除一次性 layout 迁移能力 | 中 | 约 570 行生产代码及迁移专属测试/说明 | 中 |
| 6 | 删除仅“为未来保留”的三个模板模块 | 中低 | 1,413 行生产代码及约 500 行支持测试 | 中高，需产品决定 |

推荐先实施候选 1；候选 2、3 各自作为后续独立批次。候选 5、6 不应随前三项顺带删除。

## 五、证据记录

### 5.1 候选 1：删除 `run_finite.py` 的旧流水线

**Candidate**
删除当前正式 finite shift 流程之外的三组旧实现：

1. `RUN_INFINITE_UNIT_CELL`、`INFINITE_UNIT_CELL_OUT_DIR`、
   `INFINITE_UNIT_CELL_RUNS` 和 `run_infinite_unit_cell_cases()`；
2. `build_finite_cavity(gap)`、旧 `draw_geometry()`/
   `draw_annotated_geometry()`、`run_finite_simulation()`、`run_gap()`；
3. 本地旧 lattice-Fourier 链：`unit_cell_rho_grid()`、
   `read_hz_center_field()`、`interpolate_complex_field()`、
   `finite_lattice_fourier_dir()`、`decompose_profile_modes()`、旧绘图函数、
   `finite_lattice_fourier_mode()` 和本地 `run_finite_lattice_fourier_postprocess()`。

实施时必须再次逐符号搜索，不能仅按行号整段删除。

**Burden**
同一个 4,000 行主入口同时描述正式 shift 流程、旧 gap preview、旧 unit-cell 求解和旧
per-mode Fourier 布局。旧 Fourier 又与当前共享模块使用相同函数名，迫使正式 import 使用
`run_exported_*` 等别名区分两个所有者。

**Reachability**

- `main()` 的正式路径是参数 pair -> `run_finite_cavity_shift()` ->
  `run_finite_full_case()`；
- `run_gap()` 没有仓库内调用者；
- 旧本地 Fourier 只在自己的旧调用链中互相调用，正式 full/quarter 路径调用的是
  `comsol_workflow.lattice_fourier_postprocess.run_finite_lattice_fourier_postprocess()`；
- `RUN_INFINITE_UNIT_CELL=False`，代码注释也称其为 old branch；当前独立 unit-cell
  入口已由 band pair/solo 拥有；
- 测试仅对旧符号做 `hasattr` 兼容 patch 或检查旧输出根的父目录，没有行为消费者；
- 本地结果中没有发现 `unit_cell_band/single_hexagon_unit_cell/`。

**Rationale**
Git 历史显示这些代码早于后续“standardize workflow parameters/result layouts”和
“expand finite mode selection”提交。当前共享后处理模块和正式 mode-centric 输出已经接管
其职责，原始理由不再成立。

**Cut**
删除以上符号、只被它们使用的 import/常量/测试断言和注释；保留当前
`run_finite_full_case()`、checkpoint/resume、staging/finalize、single-mode selection、
far-field 和共享 lattice-Fourier 调用。

**Consequence**
不再允许通过手工改常量恢复旧的 `single_hexagon_unit_cell`、gap preview、两位补零 mode
和根级 Fourier 目录。对应能力分别由 band 主入口和当前 finite 流程替代。

**Confidence / risk**
仓库内置信度高；风险中低。主要不确定性是有人在仓库外直接 import 这些未声明为 public
但可导入的函数。

**Proof**

- residue search 对以上符号为零；
- `tests/test_run_finite_farfield.py`、`tests/test_finite_cavity_field_exports.py`、
  `tests/test_run_strip_shift_scan.py`、`tests/test_run_finite_quarter_four_modes.py` 通过；
- full 和 quarter 的 preflight/config/output-path 单元测试保持一致；
- 不启动 COMSOL的情况下，正式 flow 的 monkeypatched 测试仍只调用共享 Fourier owner；
- `py_compile`、全仓相关 pytest、`git diff --check` 通过。

**Net effect**
删除旧概念和重复实现，不新增适配层；这是当前收益/置信度最佳的第一批。

### 5.2 候选 2：删除被否决的 unit-cell 2D 拟合重建

**Candidate**
删除 topology-adaptive valley reconstruction、minimum-curvature、quadratic B5 fit 三套
实验物理结果生成器，以及对应 CLI：

- `--valley-reconstruction-b4-only`；
- `--minimum-curvature-b4-only`；
- `--b5-quadratic-valley-fit-only`。

**Burden**
`unit_cell_2d_analysis.py` 为这些路线保留约 2,000 行优化、拟合、拓扑分类、诊断表和绘图
逻辑；`test_unit_cell_2d_analysis.py` 还维护约 400 行专属测试，CLI 和长篇文档继续把它们
表现为可用能力。

**Reachability**

- 标准 `run_analysis()` 不调用这些函数；
- 生产消费者只有 `scripts/analysis/unit_cell_2D.py` 的三个显式可选开关；
- 其余仓库消费者是专属测试和历史文档；
- `docs/unit_cell_2d_b4_nodal_line_reconstruction.md` 顶部已明确写明所有这类路线均被
  **rejected as physical-result generators**；
- 当前权威替代是 `unit_cell_2D_valley_refine.py` + 真实 COMSOL 新点 + Delaunay 线性显示；
- 本地仍存在 legacy full unit-cell 2D 结果，因此应保留 full-grid 标准读取/绘图语义，
  但这不要求保留拟合重建器。

**Rationale**
这些路线曾用于改善粗网格 valley 的视觉连续性。后续设计记录认定拟合 floor、人工中心
恢复和拓扑假设不能作为物理振幅来源，原动机已经由真实 COMSOL 局部加密取代。

**Cut**
按消费者图删除专属 dataclass、检测/拟合/优化函数、三类 redraw/plot 函数、CLI 分支和
专属测试。保留标准 A--D、`FullFieldView`、C2v quarter view、winding、标准 B4、source
hash 和真实 COMSOL refinement。历史说明压缩为短的“拒绝方案”记录，保留假设、失败原因、
被放弃能力和重新引入条件，不抹去教训。

**Consequence**
以后不能从旧 full-grid CSV 重新生成三类历史比较图；已有 PNG/CSV/JSON 结果保持原位，
Git 历史仍能恢复实现。该能力删除属于可观察 CLI 收缩，实施前需用户明确批准。

**Confidence / risk**
仓库内置信度高；风险中等，来自可能存在的仓库外手工重绘流程。

**Proof**

- 三个 CLI 开关和所有专属符号 residue 为零；
- `tests/test_unit_cell_2d_analysis.py` 删除专属测试后，其余标准 A--D、C2v、winding 测试通过；
- `tests/test_unit_cell_2d_valley_refinement.py`、`tests/test_run_unit_cell_2d.py`、
  `tests/test_run_unit_cell_2d_zeta_scan.py` 通过；
- 对固定 fixture 比较标准 A--D 文件集合、summary schema、clim/坐标合同；
- 不改写任何 `scripts/.out/unit_cell_2D/` 文件。

**Net effect**
真正退休三套已否决模型，不以新 wrapper 或 archive package 继续维护同一义务。

### 5.3 候选 3：删除 strip 的不可达 full-finite 和诊断模式

**Candidate**

1. 删除 `RUN_FINITE_FULL`、`run_finite_full_case()` 及仅由它消费的 full-finite helper；
2. 删除 `RUN_DIAGNOSTIC_PLOTS` 和其恒假分支下的几何/shift/Fourier 诊断图生成函数，正式
   Fourier 调用直接固定 `generate_plots=False`。

**Burden**
strip 主入口为已经禁止的第二套 finite solver 和不允许生成的诊断输出保留约
600--800 行代码、绘图概念和 fake-COMSOL 测试。

**Reachability**

- `run_case()` 在 `RUN_FINITE_FULL=True` 时立即抛错并要求改用 `run_finite.py`；
- `AGENTS.md` 明确规定 strip 不得运行 full finite；
- full-finite 函数只有专属测试调用，没有生产调用；
- `RUN_DIAGNOSTIC_PLOTS=False`，AGENTS 要求保持关闭；活动调用仅把该值传给共享
  Fourier 后处理；
- 当前本地所有发现的 strip 配置和 run summary 都是 layout v3。

**Rationale**
这些分支来自 strip/finite 尚未分离和输出布局仍在调试的阶段；当前输出根、入口与正式图
合同已经分离。

**Cut**
先做 3A full-finite，再做 3B diagnostic，各自独立提交和验证。保留 strip 几何构建、
cut 坐标变换、正式 `Hz_Re/Hz_Im/Wem`、`cutline_Wem.png`、评分、staging/finalize 和
progress log。

**Consequence**
不能再通过修改源码常量恢复旧调试图；需要调试时应增加一次性、明确授权的纯分析工具，
而不是把不可达分支重新放回正式入口。full finite 继续由 `run_finite.py` 提供。

**Confidence / risk**
高置信度。风险低至中，主要是未来几何调试便利性下降。

**Proof**

- `RUN_FINITE_FULL`、`RUN_DIAGNOSTIC_PLOTS` 和专属 helper residue 为零；
- `tests/test_run_strip_output_layout.py`、`tests/test_run_strip_shift_scan.py`、
  `tests/test_run_strip_bulk_rotation.py`、`tests/test_strip_1d_fourier_quotient.py` 通过；
- monkeypatched `run_strip_cut_case()` 仍生成相同 config、MPH 路径、正式场和评分调用；
- 正式测试断言不再维护“恒为 False 的开关”，改为断言输出集合不含诊断文件。

**Net effect**
删除两个不可达状态，不引入新配置。

### 5.4 候选 4：收敛文档事实源

**Candidate**
把入口、输出根、模块清单和安全规则分别交给唯一文档所有者，并删去重复清单：

- `pyproject.toml`：安装后 CLI；
- `scripts/README.md`：直接脚本入口和最小使用说明；
- `README.md`：项目定位和主数据流；
- `AGENTS.md`：agent 安全、昂贵计算审批、输出/兼容约束；
- `comsol_workflow/README.md`：底层模块所有权。

**Burden / evidence**

- `AGENTS.md` 仍引用不存在的 `run_unit_cell_band.py` 和对应不存在测试；
- `scripts/README.md` 前文承认四个主入口，却又写“上述三个工作流根目录”和“三个活动
  子目录”；
- `comsol_workflow/README.md` 写“15 个 Python 模块”，当前实际清单为 21 个模块；
- 根 README 写“three primary entry points”，而安装命令已有 8 个；
- 根 README、scripts README、AGENTS 和 handoff 合计超过 1,200 个非空行，重复描述
  相同参数、路径和完成判据。

**Cut**
先修正错误引用，再把长参数语义和操作合同保留在 AGENTS/对应 spec，README 只链接；不
生成新的文档同步器。当前未跟踪的 handoff 是用户已有工作，第一批不删除；待其交接目的
结束后再由用户决定是否吸收唯一信息并移除。

**Consequence**
README 不再独立复制全部 agent 运行规约，读者需沿明确链接查看详细合同。

**Confidence / risk / proof**
置信度高、风险低。通过精确路径搜索、Markdown link 检查、命令清单与 pyproject 人工
对照验证。文档批次不得顺带改程序行为。

**Net effect**
减少同步义务并修复已出现的 split truth；不新增生成系统。

### 5.5 候选 5：退休一次性迁移能力（需外部盘点）

**Candidate**

- `scripts/analysis/migrate_strip_output_layout_v3.py`；
- `run_band_pair.py:migrate_legacy_unit_cell_layout()`、`LEGACY_CASE_OVERVIEW_FILES`、
  `_tree_file_stats()` 和专属测试。

**Burden**
约 570 行迁移代码拥有文件移动、staging、MPH 加载、旧路径识别和校验责任；它们不参与
新计算。

**Reachability / evidence**

- 当前本地 strip 结果全部为 layout v3；未发现 `.layout_v3_migration` 证据；
- 当前本地 unit-cell band 结果使用 `01_results/10_overview/99_config` 正式布局；
- unit-cell 迁移函数没有生产调用者，只有专属测试和历史 spec；
- strip 迁移器是直接手工入口，scripts README 仍公开说明；
- 无法证明其他磁盘、备份或外部路径不存在待迁移结果。

**Consequence**
删除后，未迁移的外部旧结果不能再用现成工具升级。现有正式结果不受影响，代码可从 Git
历史恢复。

**Confidence / risk**
本仓库证据中等偏高，整体置信度中等；外部数据盲点使其不能自动批准。

**批准前所需事实**

1. 用户确认没有仍需迁移的仓库外 layout v2 strip/unit-cell 结果；或
2. 先对明确的数据根运行 dry-run 并完成最后一次迁移；
3. 记录迁移工具最后支持的 layout 和恢复 commit。

### 5.6 候选 6：删除无主“未来模板”（需产品决定）

**Candidate**

- `comsol_workflow/finite_patch_lattice.py`；
- `comsol_workflow/simulation_spatial/square_unit_cell.py`；
- `comsol_workflow/simulation_spatial/rectangle_finite_size.py`；
- 只维护这些模块的测试和 README 条目。

**Burden**
三个模块共 1,413 个非空生产代码行，另有约 500 行专属或混合测试。它们复制晶格生成、
COMSOL selection/材料/求解/export 等责任。

**Reachability**

- `finite_patch_lattice.py` 只有 `test_bulk_side_corner_lattice_generation.py` 和 package
  README 消费；
- square/rectangle runner 只有 `test_simulation_utils_spatial_comsol.py` 和文档消费；
- package README 明确写“未被当前三入口调用”“作为今后模板保留”；
- Git 历史显示它们来自早期通用 spatial/finite 尝试；
- setuptools 会发布整个 package，因此仓库外 import 仍是未知消费者。

**Consequence**
删除会放弃未来正方形 unit-cell、独立 rectangle finite runner 和通用 side/corner patch
构建模板。当前 `finite_geometry="square"` 不依赖这些模块，而是复用统一
`FiniteGeometryPlan` + `simulation_utils.py`，因此当前正式 square finite 能力不受影响。

**Confidence / risk**
仓库内消费者图置信度高，但产品决策和外部 import 不确定，整体为中低。只有用户明确表示
不再维护这些未来能力时才执行；不要用“先移动到 archive/”替代删除，那会保留同一维护义务。

**Proof**
删除后搜索模块名和公开类，移除专属测试；运行 finite geometry、simulation utils 和全部
非 blueprint 测试，确认当前 hex/square full/quarter 合同仍在统一 runner 上通过。

## 六、明确保留或暂不处理的高价值线索

### 6.1 保留：finite full 与 quarter 两个入口

quarter 已复用 full 的参数、几何 plan、输出 finalizer 和后处理；它自身仍拥有 PEC/PMC
组合、第一象限裁剪、场奇偶恢复和独立 COMSOL 生命周期。进一步强行合并只会把条件分支
塞回一个更大的状态机，不能减少合同。

### 6.2 保留：`run_band_solo.py`

该文件只有约 114 个非空行，复用 pair 的求解、tracking 和绘图实现，同时提供明确的
“只算 cavity、不创建 cladding/comparison”命令合同。删除会移除真实能力，合并参数开关
反而增加 pair 入口状态。

### 6.3 暂保留：`recover_finite_quarter_selected_mode.py`

它是 solved MPH 到完整后处理的恢复边界，包含 dry-run、staging、输入校验和 backup。
当前仍有 finite-quarter 计算活动迹象，且恢复逻辑属于数据防丢能力；在所有相关进程退出、
当前 four-mode case 完成并确认不再有中断 case 前，不应删除。之后可重新审计其真实使用。

### 6.4 保留：历史参数和结果读取兼容

本地确有 legacy full unit-cell 2D、scalar/directional finite 等历史结果。缺失
`isQuarter` 的 full 解释、旧 config 默认值和 trend 对历史路径的只读兼容仍有持久化消费者。
这些不是“看起来旧”就能删除的分支。

### 6.5 暂保留：`intensity_maps_normH.png` 兼容名

该名字广泛存在于正式结果和当前文档中，且可能被外部分析引用。虽然 flatness trend 会把
canonical PNG 复制为同字节别名，但没有外部消费者清单或弃用窗口，当前不具备删除权限。

### 6.6 不做：把原子 JSON writer 全部抽成新 utility

多处函数文本相似，但分别位于参数快照、求解状态、分析 summary 等不同持久化边界。仅为
消除几行重复而引入新的共享模块会搬运而非消除复杂度。只有在统一 durable-write 合同的
独立目标下才应处理。

### 6.7 保留：`srip1d_trend.py` 当前拼写

这是已记录且被直接命令使用的历史文件名。简单改名需要兼容 wrapper，净概念数不降；若
未来决定破坏兼容，可与一次明确的 CLI 收缩一起处理。

### 6.8 不做：批量删除历史 `docs/spec/`

大量 plan 已完成并记录了输出布局、模式命名和失败原因。它们不是运行时代码，但保留了
防止旧错误重现的设计证据。只应在当前 owner 已吸收唯一理由、入链修复完成后逐份判定，
不能按日期批量清理。

### 6.9 未决：`SI_thoery.md` 与 `SI_thoery_v2.md`

两份文档有大段重合，但 v2 声明旧名为 alias，当前数值代码又明确引用 v2 的 Appendix C；
历史拟合 spec 同时引用两者。它们也存在实质内容差异，不能仅凭相似判重复。需要论文作者
确认哪份是 manuscript canonical、另一份是否仍含唯一推导，再决定合并；本轮不排名删除。
`inner_product_decomposition.md`/`_complex.md` 和两份 parametrization 文档分别处理不同
数学域/约束，也不判为重复。

## 七、建议实施顺序

### Batch 0：确认授权与静止状态

1. 等待当前 COMSOL/Python 任务退出，确认没有源码仍在被活动任务读取；
2. 重新检查 `git status -sb`，保护本轮之前的用户修改；
3. 用户分别确认：
   - 是否放弃 unit-cell 2D 三个实验重绘 CLI；
   - 是否放弃 strip 诊断图的源码开关；
   - 是否存在仓库外待迁移结果；
   - 是否放弃三个未来模板模块；
4. 对第一批相关测试取得不启动 COMSOL的 baseline。

### Batch 1：旧 finite 流水线

- 实施候选 1；
- 一个 ownership boundary、一个可审查 diff；
- 验证后在本文回填实际删除符号、行数、测试和残余风险。

### Batch 2：被否决的 unit-cell 2D 重建

- 仅在用户批准 CLI 收缩后实施候选 2；
- 先保留简短的 rejected-alternative 理由，再删除实现、CLI 和专属测试；
- 不改标准 A--D，不改真实 valley refinement，不触碰历史结果。

### Batch 3：strip 不可达状态

- Batch 3A 删除 full-finite fossil；
- Batch 3B 删除 disabled diagnostics；
- 每批分别做 residue、strip 聚焦测试和输出合同检查。

### Batch 4：文档事实源

- 修复已确认的入口/数量错误；
- 删除重复说明并改为唯一 owner 链接；
- 不改写仍在使用的用户 handoff，除非用户单独批准。

### Batch 5：条件性删除

- 完成外部旧结果盘点后，才处理候选 5；
- 完成产品能力确认后，才处理候选 6；
- 两者必须分开实施，不能作为“顺手清理”。

## 八、统一验收合同

每个批次至少完成以下适用层级，并逐项报告，不能用一个窄测试代表全部通过：

1. **Residue**：搜索已删符号、flag、路径、输出文件名、README 和测试引用；
2. **Decisive test**：运行最能暴露该删除错误的聚焦测试；
3. **Import/compile**：对受影响文件运行 `py_compile`；
4. **Local gates**：运行对应 workflow 的纯 Python/mock pytest；
5. **Repository gates**：资源允许时运行除 blueprint 边界测试外的根项目 pytest；
6. **Boundary comparison**：比较正式 CLI、config schema、case 路径、mode 编号、
   staging/finalize、resume 和 summary；
7. **COMSOL safety**：验证命令不得创建 `SimulationRun` 或 eigensolve；如最终需要 COMSOL
   smoke，必须另行做昂贵计算 preflight 和用户授权；
8. **Data safety**：`scripts/.out/` 的既有文件数、mtime 和 hash 不因源码简化改变；
9. **Diff audit**：`git diff --check`、逐文件 diff 和工作树范围检查。

建议命令集合（按批次取子集）：

```powershell
uv run python -m pytest -q `
  tests\test_finite_cavity_field_exports.py `
  tests\test_run_finite_farfield.py `
  tests\test_run_finite_quarter_four_modes.py `
  tests\test_run_strip_shift_scan.py

uv run python -m pytest -q `
  tests\test_run_unit_cell_2d.py `
  tests\test_unit_cell_2d_analysis.py `
  tests\test_unit_cell_2d_valley_refinement.py

uv run python -m pytest -q `
  tests\test_run_strip_output_layout.py `
  tests\test_run_strip_bulk_rotation.py `
  tests\test_strip_1d_fourier_quotient.py

uv run python -m py_compile scripts\run_main\run_finite.py
uv run python -m py_compile scripts\run_main\run_strip_1d.py
uv run python -m py_compile comsol_workflow\unit_cell_2d_analysis.py
git diff --check
```

## 九、回滚策略

- 每个候选独立提交；不把多个 ownership boundary 混入同一提交；
- 所有建议均为源码/文档删除，不迁移持久化结果，正常回滚为 revert 对应提交；
- 若 decisive test 暴露仓库外/动态消费者，立即恢复当前批次并把该消费者写入本方案；
- 迁移工具只有在外部结果盘点完成后删除；必要时可从记录的 Git commit 恢复，而不是长期
  保留 archive 副本；
- 不回滚或覆盖本轮之前的用户工作树修改。

## 十、本轮实际执行回执

### 10.1 审计阶段

**Scope**
除 `blueprints/`、`comsol-data-analysis/` 外的 Broad simplification audit；只新增本文。

**Baseline**
已读取项目指令、入口清单、共享参数、主要 README、相关历史 spec、生产/测试符号引用和
本地结果 metadata；已记录 dirty worktree。Git 沙箱内一度因 Windows `0xc0000142` 无法
启动，改为经批准在沙箱外完成只读 `git status`/`git log`，未写 Git 元数据。

**Artifacts**
仅新增 `docs/spec/plan-execute_20260902_codebase_simplification_audit.md`。

**Behavior**
没有修改运行行为、配置、测试、计算结果或 COMSOL 进程。

**Verification**
本轮执行的是只读消费者搜索、当前结果布局检查、Git history/status 和最终 Markdown/diff
检查；没有运行 pytest，因为尚未改代码且当前工作站有活动/近期 COMSOL 资源压力迹象。

**Residual risk**
仓库外 Python imports、外部历史结果目录和作者个人重绘命令不可从当前仓库穷尽；候选 2、
5、6 因此保留明确批准门槛。

**Undo**
若不保留本审核，只删除本 spec 文件即可；本轮没有其他副作用。

### 10.2 Batch 1：删除旧 finite/unit-cell/Fourier 流水线

**Scope**
只处理候选 1 的单一 ownership boundary：`scripts/run_main/run_finite.py` 中已被正式
finite shift 流程、band 入口和共享 `lattice_fourier_postprocess.py` 取代的旧实现。

**Baseline**
实施前确认项目相关 Python/COMSOL 进程已经退出；目标文件存在用户先行完成的 finite
1-based layer 语义修改，包含 `CAVITY_LAYERS`、`BULK_RADIUS=CAVITY_LAYERS-1` 和
parent-layer metadata。本批把这些改动视为必须保留的基线。聚焦测试基线为 162 passed。

**Retired obligation**
删除以下不再受正式入口支持的维护义务：

- `RUN_INFINITE_UNIT_CELL`、`INFINITE_UNIT_CELL_OUT_DIR`、`INFINITE_UNIT_CELL_RUNS`
  及内嵌的旧 periodic unit-cell 求解/模式分析；
- gap-based `build_finite_cavity()`、旧裁剪预览、`run_finite_simulation()` 和 `run_gap()`；
- `run_finite.py` 本地的 finite-lattice Fourier 采样、分解、绘图和输出流水线；
- 上述路径专属 imports、配置、`main()` dispatch 和两处只维护旧表面的测试断言。

**Artifacts**
修改 `scripts/run_main/run_finite.py`、`tests/test_finite_cavity_field_exports.py` 和
`tests/test_run_band_pair.py`；更新本文。未创建新源码目录，未修改 `scripts/.out/`。

**Realized net effect**
`run_finite.py` 从实施前 4,055 行降至 3,097 行，删除 958 行；测试删除 4 行旧兼容
维护，共净删除 962 行，没有新增 wrapper、迁移层或依赖。共享 Fourier owner 和正式
full/quarter 入口保持唯一实现。

**Behavior**
保留当前 full finite、quarter finite、checkpoint/resume、mode selection、场导出、
共享 lattice-Fourier、intensity map、far-field、mode-centric layout，以及用户已有的
1-based finite layer 语义。主动放弃从 `run_finite.py` 手工恢复旧 single-unit-cell、gap
preview、旧根级 Fourier 输出的能力。没有启动 COMSOL，没有覆盖或删除正式结果。

**Verification**

- residue：`run_finite.py` 和正式调用路径中旧常量、函数和本地 Fourier owner 均为零；
  搜索剩余命中仅为共享模块权威实现、独立 polygon fixture 和本文历史证据；
- compile：`py_compile` 对修改的源码及两份测试通过；
- decisive baseline/post-check：同一组 finite/band/quarter/strip 聚焦测试实施前后均为
  162 passed；
- expanded finite/cavity/far-field：组合运行曾出现 1 个 Matplotlib 跨文件 figure
  污染失败；失败测试单独运行通过，随后按文件隔离为 far-field 16 passed、其余 70 passed；
- repository gate：`pytest -q tests --ignore=tests\test_blueprint_comsol_geometry_parity.py`
  为 524 passed、3 个既存 `PytestReturnNotNoneWarning`；
- lint：当前环境未安装 `ruff`，该层不可用；以 compile、pytest、残留搜索和 diff audit
  代替，但不把这些检查描述为 lint 通过；
- data/process：写入前确认无项目 Python/COMSOL 进程；本批没有调用 `SimulationRun`、
  eigensolve 或结果生成入口。

**Residual risk**
仓库外代码仍可能直接 import 已删除的未声明 public 函数；当前包发布边界使该风险无法在
仓库内完全排除。独立 `test_cavity_hole_polygon_clipping.py` 继续测试通用 polygon API，
不再被解释为旧 `run_finite.py` 流程的行为消费者。

**Retained candidates**
候选 2 涉及三个历史 unit-cell 2D CLI，候选 3 涉及 strip 诊断能力；二者均保留并等待
独立批准。候选 5、6 继续受外部数据盘点和产品能力决策门槛约束。

**Undo**
本批最终随本次简化提交记录；`run_finite.py` 与 `test_run_band_pair.py` 在本批前已有用户修改，不得用
整文件 `git checkout` 回滚。若撤销 Batch 1，应只反向应用本批删除 hunks，恢复旧 imports、
常量、函数、dispatch 和两处测试断言，同时保留用户的 finite layer 语义修改。无数据或
配置恢复步骤。

### 10.3 Batch 2：退役三套 unit-cell 2D 拟合重建路径

**Scope**
按用户明确批准，移除 topology-adaptive valley reconstruction、minimum-curvature B4 和
quadratic B5 fit 的离线 CLI、专属实现、测试和当前可调用文档；不处理标准 A--D 或真实
COMSOL valley refinement。

**Baseline**
写入前确认没有项目相关 Python/COMSOL 进程。`scripts/analysis/unit_cell_2D.py`、
`comsol_workflow/unit_cell_2d_analysis.py`、专属测试和退役长文没有用户已有改动；
`AGENTS.md` 有大量用户改动，因此只精确替换其中三行旧 CLI 陈述。聚焦基线为 99 passed。

**Retired obligation**

- 删除 `--valley-reconstruction-b4-only`、`--minimum-curvature-b4-only` 和
  `--b5-quadratic-valley-fit-only` 三个互斥开关及其 dispatch/import；
- 删除低值检测、中心不定二次型、bootstrap topology、X-like 多项式分支、人工 floor、
  low-value relaxation、薄板 minimum-curvature 优化和中心二次拟合实现；
- 删除三套重建 PNG/CSV/JSON writer 与专属 dataclass、SciPy imports；
- 删除 12 个专属测试实例和只服务它们的 fixtures/imports；
- 将原 reconstruction 长文从当前可复现说明收敛为 79 行退役方案记录。

**Artifacts**
修改 `comsol_workflow/unit_cell_2d_analysis.py`、`scripts/analysis/unit_cell_2D.py`、
`tests/test_unit_cell_2d_analysis.py`、`docs/unit_cell_2d_b4_nodal_line_reconstruction.md`、
`AGENTS.md` 和本文。未修改 `scripts/.out/`，未删除既有历史 PNG/CSV/JSON。

**Realized net effect**
核心分析模块删除 2,043 行，离线入口删除 62 行，活动测试删除 509 行；活动代码和测试合计
净删除 2,614 行。历史长文 diff 为 54 行新增、2,025 行删除，最终保留 79 行。没有新增
wrapper、兼容开关、模型实现或依赖。

**Behavior**
标准 `run_analysis()`、A--D、标准 B4、`_sha256_file` source audit、legacy full-grid
正常读取、Quarter C2v 只读视图、winding 和 COMSOL valley refinement 保持。离线入口帮助
现在只接受可选 `series_dir`；主动放弃重新生成三套历史比较输出的能力。已有结果保持原位，
但其文件名不再代表当前支持接口。

**Verification**

- CLI：`unit_cell_2D.py --help` 只显示 `series_dir` 和 `--help`；
- residue：活动 scripts/core/tests/AGENTS/当前退役文档中三个 flag、专属函数和 diagnostics
  schema 均为零；全仓剩余 flag 命中只在历史 specs 和本文；
- compile：核心模块、离线入口和测试 `py_compile` 通过；
- focused：基线 99 passed；第一次修改后检查因误删标准 winding parametrization 装饰器为
  83 passed、1 error，随后只恢复该装饰器，最终为 87 passed；减少的 12 个实例均属于
  被退役能力；
- repository gate：`pytest -q tests --ignore=tests\test_blueprint_comsol_geometry_parity.py`
  为 512 passed、3 个既存 `PytestReturnNotNoneWarning`；
- diff：目标文件 `git diff --check` 通过；删除接缝和标准函数清单已人工复核；
- process/data：没有调用 `mph`、`SimulationRun`、eigensolve、重绘或结果 writer。

**Residual risk**
仓库外代码可能直接 import 被删除的未声明 public 函数；已有历史输出只能读取，不能由当前
代码重建。两份历史 spec 继续记录当时的 CLI 和输出，不应被误读为当前操作说明；79 行
退役记录已把当前 owner、替代方案和重新引入条件置于文首。

**Decision-record handling**
旧长文分类为 fully displaced implementation record，但其失败教训仍有价值；已保留三种
先验、为何不能代表 Maxwell 解、主动放弃能力、当前 COMSOL refinement owner 和恢复条件。
历史执行 spec 不重写，以免抹除当时事实。

**Retained candidates**
候选 3 的 strip full-finite fossil 与 disabled diagnostics、候选 4 的文档事实源收敛仍未
实施；候选 5、6 继续等待外部数据盘点或产品决策。

**Undo**
本批最终随本次简化提交记录。撤销时应只反向应用 Batch 2 的六个文件 hunks；尤其 `AGENTS.md` 和本文含有
用户/前一批改动，不能整文件 checkout。恢复旧 CLI 还必须同时恢复专属实现、测试和长文，
避免留下可调用但无维护证据的半恢复状态。无数据或配置恢复步骤。

## 十一、建议决策

候选 1 和 2 已完成。后续建议顺序为：候选 3 -> 候选 4；候选 5 和 6 暂缓。

下一批候选 3 将删除 strip 中入口已拒绝的 full-finite fossil，并评估移除永久关闭的
diagnostic plot 分支；仍应作为独立批次获得用户确认后实施。
