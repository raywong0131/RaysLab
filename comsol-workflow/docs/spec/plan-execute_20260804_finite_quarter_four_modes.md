# finite-quarter 四个 Gamma 基准模式批处理：方案与执行记录

日期：2026-08-04

## 方案概述

新增一个独立的 finite-quarter 四模批处理入口。入口从当前或用户显式指定的
unit-cell `selected_bands.csv` 中读取 cavity 的 `p2/p1/d1/d2` Gamma 点频率与
project-internal `px/py/dx/dy` 分量，按 `p2 -> p1 -> d1 -> d2` 顺序调用现有
`run_finite_quarter.py` 求解路径。每次求解使用对应 Gamma 频率作为 eigensolver
shift，并按内部模式分量选择一组 x/y 轴 PEC/PMC 边界。

本任务不做空间包络基模判断，也不使用 Fourier 分量再次确认有限腔模式身份；保留
现有 `is_valid` 物理有效性检查，并在有效候选中选择与目标频率绝对差最小的 mode。

## 范围

- 新增 `scripts/run_main/run_finite_quarter_four_modes.py`。
- 参数化 `scripts/run_main/run_finite_quarter.py` 的目标频率、目标 band、批次
  symmetry IDs 和 unit-cell 来源元数据；原 `SYMMETRY_IDS=[1]` 入口保持兼容。
- 让 `scripts/run_main/run_finite.py` 的最近频率选择辅助函数接受可选的显式目标频率，
  不改变 full finite 默认目标。
- 在 `tests/` 增加 unit-cell Gamma 数据解析、分量到边界映射、固定运行顺序、目标频率
  传递与无 COMSOL orchestration 测试。
- 更新 README 使用说明。
- 不改 `scripts/parameter.json`，不启动 COMSOL，不移动或覆盖现有正式结果。

## 数据流

```text
scripts/parameter.json
  + unit-cell series/99_config/config.json
  + unit-cell series/10_overview/cavity/selected_bands.csv
  -> 校验 cavity 参数、mesh、Gamma 双分支频率与 dominant_mode
  -> p2/p1/d1/d2 任务列表
  -> internal px/py/dx/dy 映射到 symmetry ID 1/2/3/4
  -> 对每个 finite shift pair 构建一次 quarter geometry
  -> 按 p2, p1, d1, d2 顺序运行四个独立 SimulationRun
  -> 每个 symmetry case 写入原 mode-centric finite 输出
  -> four_mode_run_summary.json 记录来源、顺序、目标与状态
```

分量到边界映射：

| internal mode | symmetry ID | x-axis boundary | y-axis boundary |
| --- | ---: | --- | --- |
| `px` | 1 | PEC | PMC |
| `py` | 2 | PMC | PEC |
| `dx` | 3 | PEC | PEC |
| `dy` | 4 | PMC | PMC |

## 输入与兼容策略

- 默认 unit-cell series 是当前 `ACTIVE_PARAMETERS.case_parameter_label` 派生的标准目录。
- 允许通过命令行显式指定既有 unit-cell series；此时仍校验 cavity 紧凑参数和 mesh，
  cladding 参数不作为 cavity Gamma 数据兼容条件。
- Gamma-M 与 Gamma-K 的 Gamma 行必须同时存在，且同一 band 的频率和
  `dominant_mode` 必须一致。
- 四个 band 必须恰好映射到四个不同 symmetry ID；歧义或重复在启动 COMSOL 前失败。
- `finite_eigenmode_count` 继续只从共享参数读取。当前值为 1；若将来大于 1，则在有效
  候选中选目标频率最近者。
- 原 `run_finite_quarter.py` 单 symmetry 用法、输出命名、恢复开关和后处理开关不变。
- 四模批处理使用既有 `symmetry_<id>_x..._y...` 子目录，不引入第四个一级输出目录。

## 验证方案

- unit test：合法 Gamma 双分支解析并保持 `p2,p1,d1,d2` 顺序。
- unit test：缺行、频率不一致、分量不一致、cavity/mesh 不匹配、symmetry 重复均失败。
- unit test：显式目标频率进入 `SimulationConfig.eigenfrequency_shift` 和最近频率选择。
- orchestration test：mock quarter runner，确认每个 shift 内严格执行四次且不启动 COMSOL。
- 运行相关 pytest、两个新入口的 `py_compile` 和 `git diff --check`。

## 回滚方案

删除新增四模入口和测试，撤回三个现有入口中的可选参数即可。由于默认参数保持原值、
不修改共享 JSON 且不执行正式计算，回滚不涉及结果迁移。

## 执行记录

### 实际改动

- 新增 `run_finite_quarter_four_modes.py`：
  - 默认从当前共享参数派生 unit-cell series，也接受 `--unit-cell-series`；
  - 校验保存的 cavity `b0_nm/eta/zeta` 与 mesh，明确忽略不参与 cavity Gamma 数据的
    cladding 参数；
  - 从 Gamma-M/Gamma-K 两个 Gamma 行解析四个 band，校验频率、有效性和内部模式一致；
  - 动态执行 `px->1`、`py->2`、`dx->3`、`dy->4` 映射，并要求四个 ID 不重复；
  - 每个 shift 只构建一次 Python quarter geometry，随后严格按
    `p2,p1,d1,d2` 调用四次现有 quarter runner；
  - 增加 `--preflight-only`、输出冲突检查、增量 batch summary 和首错停止记录。
- `run_finite_quarter.py` 新增可选目标频率、band、内部模式、unit-cell 来源和 batch
  symmetry IDs；目标频率同时进入 eigensolver shift、case metadata 和下游最近频率选择。
  未传入时继续使用原共享中心频率和 `SYMMETRY_IDS=[1]`。
- `run_finite.py::resolve_finite_analysis_mode_indices` 接受可选显式目标，full finite 默认
  调用不变。
- 新增 8 个四模聚焦测试；同时把两个既有 active-config 测试从硬编码 ellipse 改为按
  当前 `cladding_shift_geometry` 派生预期，没有修改参数文件。
- 更新根 README 与 scripts README。

### 实际预检

使用既有、cavity 参数与 mesh 匹配但 cladding 不同的 unit-cell series 执行：

```powershell
uv run python scripts\run_main\run_finite_quarter_four_modes.py `
  --preflight-only `
  --unit-cell-series "scripts\.out\unit_cell_band\unitcell_band_cav(245-0.96-1.156)_clad(242-0.986-0.917)_mesh5"
```

结果解析为 `p2/p1/d1/d2 -> symmetry 1/2/4/3`，频率分别为
`198.386838231/196.564396032/200.461469059/203.242606229 THz`；当前
`x=0.05,y=0.08` hex quarter geometry 为 1912 个孔。预检未启动 COMSOL，也未写正式
结果。

用户随后要求按最新 `parameter.json` 且不引入 cladding shift。复核时该文件已经由
用户设为 `cavity_layers=20`、`cladding_layers=20`、
`cladding_shift_factors=[0.0]`、`cladding_y_over_x_shift_ratio=1.0`，因此本轮没有再次
修改参数。以相同 unit-cell cavity 来源重新预检后，finite pair 为 `x=0,y=0`，hex
quarter geometry 为 7422 个孔，四模顺序、频率和 symmetry 映射不变；仍未启动
COMSOL 或写正式结果。

### 测试结果

- 聚焦与关联回归：`86 passed`。
- 更广 finite/quarter/profile/geometry 回归：`62 passed`。
- 三个相关入口 `py_compile`：通过。
- `git diff --check`：通过；仅报告工作树现有 LF/CRLF 转换警告。

### 遗留问题与最终状态

- 当前共享参数精确派生的 unit-cell series 不存在，因此默认命令会在 COMSOL 前明确
  失败；可先运行当前 unit-cell workflow，或像上述预检一样显式指定 cavity 参数与 mesh
  匹配的既有 series。
- 本任务没有启动 COMSOL、没有创建或覆盖正式结果；最新无 shift 参数由用户预先写入，
  本轮没有再次修改 `parameter.json`。
- 实施完成，等待用户决定是否进行实际四模计算。

## 2026-08-05 选择规则变更

用户将有限腔模式选择规则从“离目标 Γ 频率最近”改为“边界条件对应的
`px/py/dx/dy` Γ 分量权重最高”。四组计算仍使用 unit-cell cavity band 的 Γ
频率作为 eigensolver shift，但该频率只用于求解定位，并在分量权重完全并列时作为
次级排序条件，不再是首要选择条件。

更新后的数据流如下：

```text
每个 symmetry 求解 parameter.json 中 finite_eigenmode_count 个候选
  -> 导出有效候选的中心面 Hz（仅用于模式判别）
  -> 在 cavity bulk cells 上计算 Γ profile
  -> symmetry 1/2/4/3 分别比较 mode_px/mode_py/mode_dy/mode_dx
  -> 权重降序、Γ 频率距离升序、mode_idx 升序确定唯一模式
  -> analysis_mode_selection.csv 记录全部候选权重与选择理由
  -> 只为选中模式生成 Hz/Wem 正式图、air-plane、far-field、完整 lattice Fourier 和评分
```

兼容与输出策略：

- 规则只作用于传入 `target_internal_mode` 的 four-mode quarter 调用；full finite
  及普通 quarter 的最近频率兼容入口保持不变。
- 候选模式的 `Hz_center.parquet` 是完成 Γ 分量筛选所需的最小数据，不对未选中模式
  生成正式场图、Wem、air-plane、far-field、完整 Fourier 或 score。
- 不修改、移动或覆盖旧的最近频率规则结果。已有非空 case 继续由输出冲突检查拒绝
  覆盖；重新正式运行时需使用未占用的输出目标或由用户明确处置旧结果。
- 当前 `parameter.json` 已由用户设置 `finite_eigenmode_count=10`、shift=0；本次代码
  修改不再更改这些参数，也不启动 COMSOL。

验证补充：

- 单元测试覆盖四个目标分量、权重优先、频率并列判据、无效候选和审计 CSV。
- quarter 导出测试确认候选只生成中心面 Hz，正式场图/Wem/air-plane 仅属于选中模式。
- orchestration 测试确认 `target_internal_mode` 传入新的选择器，并且下游所有完整处理
  只收到一个选中 mode index。

### 2026-08-05 实际改动与验证

- 新增通用的 finite Γ 分量选择器：只计算候选的 Γ profile，输出完整
  `mode_px/mode_py/mode_dx/mode_dy` 证据，并按“目标分量权重降序、目标频率距离升序、
  mode index 升序”选出唯一模式。
- finite-quarter 导出拆成候选 `Hz_center` 阶段和选中模式完整导出阶段；未选中模式不再
  生成 Hz/Wem 正式图或 air-plane 数据，因此也不会进入 far-field、完整 lattice
  Fourier、intensity map 或 score。
- 四模入口要求 `finite_eigenmode_count >= 2`，批处理摘要记录选择规则和四种 symmetry
  对应分量。当前共享参数为 10 个候选、无 cladding shift。
- 更新根 README、scripts README、审计 CSV 和 case metadata 说明。
- 聚焦测试：35 passed；finite/strip/profile 扩展回归：87 passed；完整仓库回归：
  404 passed，只有 3 条既有 `PytestReturnNotNoneWarning`。
- 四个相关 Python 文件 `py_compile` 通过，`git diff --check` 通过（仅 Windows
  LF/CRLF 提示）。本次未启动 COMSOL，也未移动、删除或覆盖既有正式结果。

最终状态：代码实现和本地无 COMSOL 验证完成。旧的最近频率规则结果仍保留原位；其
输出目录会触发现有冲突保护，因此再次正式计算前需要用户明确选择新的输出位置，或
明确授权如何处理旧结果。

## 2026-08-05 p2 基模审批与选择规则二次修正

### 审批结论

用户检查后明确否决“目标 `px/py/dx/dy` Γ 分量权重最高即为基模”的规则，并确认
当前 p2 symmetry sector 的 `mode 14` 是 p2 基模。离线审计证据为：

- mode 14 的 cell-averaged intensity 在中心壳层取峰值，并从 shell 0 到 shell 20
  严格向外衰减；
- 以中心 cell 场型为载波参考得到的复包络在显著壳层内无符号翻转，最大相位偏差约
  3.14 度；
- mode 9 虽然 `mode_px=0.99945`，但复包络有两次径向符号翻转，并呈明显多瓣高阶
  结构；
- mode 14 频率为 197.982801 THz，Q 为 8648.38，距 p2@Gamma 目标约
  0.404037 THz。

审批样例保存在 p2 case 的 `10_overview/p2_fundamental_mode_screen.csv` 和
`p2_fundamental_envelope_screen.png`。用户已明确确认 mode 14，因此本节形式可以固化
到共享选择器。

### 修正后的正式选择合同

每个 symmetry 仍只在有效候选内计算目标 `px/py/dx/dy` Γ 分量，目标分量用于确认
内部模式身份；基模主判据改为空间包络：

1. 在 cavity bulk cells 上采样复数 `Hz`，以中心 cell profile 为载波参考，计算每个
   cell 的复包络投影；
2. 对 cell intensity 和复包络分别做 shell 平均；
3. 基模候选必须目标分量占优、强度峰值位于中心、中心/峰值壳层不低于 0.9、显著
   壳层无径向符号翻转，且峰后归一化强度累计回升不超过 0.05；
4. 多个候选同时通过时，依次按径向节点数、峰值壳层、向外回升、中心占优、目标分量
   权重、Γ 频率距离和 mode index 排序；
5. 无候选通过时保留完整审计 CSV 并终止，不静默回退到纯度最高或频率最近模式。

`analysis_mode_selection.csv` 必须同时保存模式身份和包络指标。未选中模式仍只保留完成
审计所需的 `Hz_center`，完整处理继续限定为唯一基模。

### p2 恢复后处理方案

当前 p2 eigensolve 已完成并保存 solved MPH，不重新求解。增加一个显式的已求解
quarter 单模恢复入口：加载已有 MPH，核对 case/config/eigen table/mode index 后，只
导出 mode 14 缺少的 Wem 与 air-plane，复用现有 Hz 生成正式场图、far-field、完整
lattice Fourier、intensity map 和 score。原 mode 9 的完整处理结果移动到 series 内的
任务备份位置，不删除；随后把正式 case 的 selection、overview summary 和 objective
更新为 mode 14。

恢复入口不得调用 geometry build、mesh builder 或 eigensolve，且必须对目标冲突、
未知文件、MPH solution 和最终输出集合做显式检查。完成后视觉检查 mode 14 正式图，
运行聚焦测试、完整 pytest、相关 `py_compile` 与 `git diff --check`。
