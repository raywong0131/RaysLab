# COMSOL workflow 与 blueprints 对话交接

日期：2026-08-31

## 接手说明

本文件用于把当前对话交给新的 Codex 对话继续处理。接手后先从仓库根目录
`I:\codeXproject\comsol-workflow` 工作，并完整读取：

1. `AGENTS.md`
2. `MEMORY.md`
3. 涉及 blueprints 时再读取 `blueprints/AGENTS.md`
4. 与当前具体任务对应的 `docs/spec/plan-execute_*.md`

不要仅凭本 handoff 启动 COMSOL、覆盖正式输出或修改加工版图。本文件记录的是交接时刻
状态；准确参数和行为仍以当前配置、代码、测试及结果 metadata 为准。

## 项目定位

- 外层 `comsol-workflow` 是 Windows 本机 COMSOL 6.3 控制器，负责 unit cell、strip 1D、
  finite cavity、finite-quarter、模式分析、Fourier 和 far-field 后处理。
- `blueprints/` 是可独立迁出的微纳加工版图项目，不得依赖 `mph`、`scripts` 或
  `comsol_workflow`，也不得启动 COMSOL。
- 所有正式 COMSOL 输出只允许位于 `scripts/.out/unit_cell_band/`、
  `scripts/.out/strip_1d/` 和 `scripts/.out/finite_cavity/`。
- blueprints 正式运行输出位于 `blueprints/main/.out/`，不纳入版本控制。

## Git 状态

交接前检查结果：

- 仓库：`https://github.com/jia-yli/comsol-workflow.git`
- 分支：`raywong/finite-cavity`
- `HEAD` 与 `origin/raywong/finite-cavity` 同步指向 `63ee7e3`
- 最近三个关键提交：
  - `63ee7e3 parameterize-blueprint-window-size-and-output-names`
  - `b0fc74c add-unit-cell-2d-and-blueprint-fabrication-workflows`
  - `d669e7b expand-finite-mode-selection-and-band-analysis-workflows`
- 创建本 handoff 前没有已暂存或未暂存的受版本控制文件。
- `.vscode/` 未跟踪，未纳入提交。
- `comsol-data-analysis/` 是带独立 `.git`、环境和结果的嵌套仓库，必须继续排除，禁止当作
  外层仓库普通目录暂存。

创建本文件后，工作树会新增这一份 handoff；除非用户明确要求，不自动 commit 或 push。

### Git 凭据注意事项

Codex 的非交互进程无法访问本机 `wincredman`，会出现
`Unable to persist credentials with the 'wincredman' credential store`。用户在 VS Code
交互式终端执行原生 `git push` 可以成功。未来若需要 push：

1. 先按 `AGENTS.md` 运行状态、diff、测试和提交检查；
2. Codex 可以完成选择性暂存与 commit；
3. 若工具进程仍被 wincredman 阻塞，请用户在已登录的 VS Code 终端执行
   `git push origin raywong/finite-cavity`；
4. 随后以 `git status -sb` 和 `git log --decorate` 核对本地跟踪引用。

## 当前共享 COMSOL 参数

权威来源：`scripts/parameter.json`。

交接时主要配置为：

- cavity：`b0=245 nm, eta=0.96, zeta=1.156`
- cladding：`b0=242 nm, eta=0.98, zeta=0.93`
- finite 层数：`20/20`
- mesh：`5`
- finite eigensolver 候选数：`10`
- finite geometry：`hex`
- cladding shift：`0`，`groups + uniform`
- 共享中心频率：cavity `p2@Gamma` 来源
- unit-cell 2D 使用独立子配置，目前 cavity `zeta=1.158`、`isQuarter=1`、目标 band 为
  `p2`，并启用 valley refinement。

昂贵计算前必须重新读取 JSON，并执行 `AGENTS.md` 中的 preflight；不要把上述快照当作
永久默认值。

## 已完成的 COMSOL/分析工作

### finite-quarter 四个 Gamma 基准模式

主设计与执行记录：

`docs/spec/plan-execute_20260804_finite_quarter_four_modes.md`

关键长期结论：

- 顺序为 `p2 -> p1 -> d1 -> d2`；内部映射分别使用 symmetry `1/2/4/3`。
- Gamma 分量只用于模式身份确认；基模主判据是中心峰值、无显著径向符号翻转、峰后无
  明显回升的空间包络。
- 不得回退到“频率最近”或“Gamma 分量最高”。无合格候选时必须保留审计证据并失败。
- p2 已确认 `mode14` 为基模；具体频率、Q、审计图和恢复路径见上述 plan-execute 文档及
  正式结果目录。
- `scripts/analysis/recover_finite_quarter_selected_mode.py` 用于从 solved MPH 恢复单个已批准
  模式的后处理，不得重新建模、mesh 或 eigensolve。

### unit-cell 2D

当前已加入完整 2D q-grid 分析、zeta scan、COMSOL valley refinement 和 C2v quarter
支持。接手时优先阅读：

- `docs/spec/plan-execute_20260821_unit_cell_2D.md`
- `docs/spec/plan-execute_20260822_unit_cell_2D_zeta_scan.md`
- `docs/spec/plan-execute_20260825_unit_cell_2D_comsol_valley_refinement.md`
- `docs/spec/plan-execute_20260826_unit_cell_2D_C2v_quarter.md`
- `docs/unit_cell_2d_b4_nodal_line_reconstruction.md`

主要入口：

- `scripts/run_main/run_unit_cell_2d.py`
- `scripts/analysis/unit_cell_2D.py`
- `scripts/analysis/unit_cell_2D_valley_refine.py`
- `scripts/run_sweep/run_unit_cell_2d_zeta_scan.py`
- `scripts/run_sweep/unit_cell_2d_valley_solver.py`

## blueprints 当前状态

### 结构与入口

- 唯一正式参数：`blueprints/main/parameter.json`
- 运行入口：`blueprints/main/run_blueprint.py`
- 生产逻辑：`blueprints/workflow/`
- SEM 标定数据：`blueprints/data/`
- 独立测试：`blueprints/tests/`

不要恢复旧的 `blueprints/src/`、根目录 `parameter.json` 或根目录 `run_blueprint.py`。

### 当前参数与命名

交接时 blueprints 参数为：

- finite 层数：cavity/cladding `15/15`
- cavity/cladding 紧凑参数仍为 `245-0.96-1.156` / `242-0.98-0.93`
- triangle side correction：`global_delta_nm=-10`、
  `envelope_modulation_strength=1`
- SEM profile：`sem-air-hole-220nm-radial-v1`
- outer windows：`segmented_hex_ring_windows_v2`
- outer-window radial width 不再是固定值，而是
  `0.4 * (cavity_layers + cladding_layers) - 3.0 um`；15/15 时解析为 `9.0 um`
- 当前仍是 `engineering_preview`；dose、polarity 与生产 DRC 阈值尚未标定，不得描述为
  fabrication-ready。

case 与 DXF 名称显式带 `(global_delta_nm,envelope_modulation_strength)`，不再显示 identity
短哈希；完整 identity 仍保存在 `layout_report.json`。同名但 identity 不一致时必须拒绝
覆盖。

最新相关设计记录：

- `docs/spec/plan-execute_20260827_blueprints_project_restructure.md`
- `docs/spec/plan-execute_20260827_blueprints_workflow_flatten.md`
- `docs/spec/plan-execute_20260827_blueprints_workflow_function_groups.md`
- `docs/spec/plan-execute_20260828_blueprints_SEM_correction_apply.md`
- `docs/spec/plan-execute_20260828_blueprints_global_hole_shrink_outer_window_dimensions.md`
- `docs/spec/plan-execute_20260828_blueprints_hex_ring_windows_v2.md`
- `docs/spec/plan-execute_20260828_blueprints_preview_three_positions_correction_function.md`
- `docs/spec/plan-execute_20260828_blueprints_radial_preview_strip.md`

### 最新验证

最后一轮 blueprints 更新已完成：

- 独立测试：`99 passed`
- `py_compile`：通过
- CLI preflight：
  - case：`dxf_15-15_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_(-10,1)`
  - features：`16868`
  - DRC：`passed=True, errors=0, warnings=0`
  - preflight 未写输出
- `git diff --check`：通过

这些验证结果对应提交 `63ee7e3`；修改后必须重新运行，不能沿用为新改动的证明。

## 新对话建议的接手顺序

1. 读取三份 agent/memory 文档和用户最新请求。
2. 执行 `git status -sb`，确认除了本 handoff、`.vscode/` 和独立
   `comsol-data-analysis/` 外是否出现新改动。
3. 根据用户的新任务，只读取对应 plan-execute 和实现文件，不要无目的遍历 `.out`。
4. 若是 blueprints 修改，先用现有结果/预览进行设计确认；跨模块变化先更新中文
   plan-execute，再实施和回填。
5. 若是 COMSOL 昂贵计算，完成目标模式映射、几何、mesh、中心频率、恢复策略、输出目录
   和后处理开关的 preflight 后再启动。
6. 修改后按风险运行聚焦测试、完整测试、`py_compile`、preflight 和
   `git diff --check`；只暂存本次任务文件。

当前没有需要自动恢复或继续监控的 COMSOL 进程记录。若用户要求检查运行状态，必须重新
查询本机进程与对应 progress/stderr，而不是依据旧对话判断。
