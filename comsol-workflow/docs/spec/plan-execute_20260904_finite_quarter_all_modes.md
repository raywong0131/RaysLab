# finite-quarter 全对称性单模型批处理：方案与执行记录

日期：2026-09-04

## 2026-09-06 正式运行补充规格

本次正式运行使用 `center_frequency_thz=203.6`、`finite_eigenmode_count=6`，并只读加载
以下已完成模型作为 geometry/mesh 来源：

```text
scripts/.out/finite_cavity/
  fq4m_retry_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f197.980_e15_20260903/
  shiftx0.000_shifty0.000/symmetry_1_xPEC_yPMC/00_model/finite_quarter.mph
```

源配置已确认是 hex、20--20、cavity `245/0.96/1.156`、cladding
`242/0.98/0.93`、mesh 5、uniform groups shift `(0,0)`，quarter hole 数为 7061。
旧的 `197.98 THz` shift 和每 sector 15 个本征模只属于旧 solution，不作为新求解参数。

复用路径的执行合同为：

1. 不覆盖或改写源 MPH、源配置和源结果；目标目录仍以
   `finite_quarter_all_20-20...` 开头，并使用独立的 `finite_quarter_all.mph`。
2. 加载后校验 geometry、层数、compact cell 参数、shift、mesh 和 COMSOL 必需 tag；
   将旧模型的 generic symmetry boundary feature 升级为可切换的 PEC/PMC feature。
3. 同时更新 study 与既有 eigenvalue solver 的 shift/neigs，然后清理旧 working solution
   data；不重建 geometry，不运行 mesh。
4. 依次执行 symmetry ID `1 -> 2 -> 3 -> 4`。每轮 solve 后先复制 solution 副本并验证，
   再保存目标 MPH、切换边界并开启下一轮。
5. 正式完成计数改为 `model_build_count=0`、`mesh_build_count=0`、
   `reused_model_count=1`、`reused_mesh_count=1`、`solve_count=4`、
   `solution_copy_count=4`；聚合表预期 24 行。
6. 启动前重新检查 Python/COMSOL 进程、COMSOL 6.3 PATH、目标冲突和源文件状态；若存在
   来源不明的 COMSOL 进程，不并发启动，也不擅自终止该进程。

参数覆盖以本次用户指令为准；保留当前 `parameter.json` 中的 203.6 THz，并将
`finite_eigenmode_count` 更新为 6，使入口、配置快照、目录身份和后处理读取同一权威值。

## 当前状态

规格已实施，当前状态为 `implemented_pending_comsol_smoke`。源码、测试和文档重构已经
完成；未启动新的 COMSOL 进程，未执行 20-20 正式求解，也未创建、移动或覆盖正式计算
结果。COMSOL solution-copy 的真实保存/重载 smoke 仍需在本机 COMSOL 空闲且 PATH 配置
完成后执行。

## 背景与现状

`scripts/run_main/run_finite_quarter_four_modes.py` 当前从 unit-cell Gamma 结果构建
`p2/p1/d1/d2` 四个目标，并依次调用 `run_finite_quarter.py::run_symmetry_case()`。
它已经能在每个 shift pair 内复用一次 Python 几何构造，但每个 symmetry case 仍会：

- 创建独立 `SimulationRun` 和 COMSOL Model；
- 独立构建或加载 MPH、mesh 和 solution；
- 写入独立 `symmetry_<id>_x..._y...` 结果目录；
- 使用各自的 unit-cell Gamma 目标频率和单模选择规则。

因此，现有实现不是“只改变边界条件的四次同模型求解”，也没有形成可替代 full finite
本征谱的一套统一数据集。

## 目标

1. 用同一个 quarter geometry 覆盖四种 x/y mirror symmetry sector，并把四个 sector 的
   本征谱并集作为 full-size model 对应的完整对称性数据集。
2. 对每个 cladding shift pair 只创建一个 COMSOL Model、构建一次 geometry、生成一次
   mesh，并在这个 Model 中执行四次 eigensolve。
3. 每轮求解完成后，在同一个 MPH 内复制并保留该轮 solution 副本；确认副本可读后，才
   清理工作 solution、切换 symmetry boundary 并开始下一轮。
4. 四轮结果最终化到同一个 `01_results/`，不再生成四个 symmetry case 目录。
5. 在统一 `10_overview/` 中增加频率-Q 散点图，横纵坐标都使用线性坐标。
6. series 根目录使用 `finite_quarter_all_<cavity_layers>-<cladding_layers>_...` 命名。

## 非目标

- 不修改 `scripts/parameter.json`。
- 不改变 cavity/cladding geometry、材料、层厚、mesh、中心频率、本征模式数、far-field、
  Fourier 或其他数值参数。
- 不改变 `run_finite.py` 的 full finite 默认行为。
- 不改变 `run_finite_quarter.py` 的既有单 symmetry 入口及其兼容输出。
- 不迁移、重命名、删除或覆盖既有正式结果。
- 不在科研图的颜色、图例、标题或标注中展示 symmetry ID、PEC 或 PMC。

## 参数权威源与求解合同

四轮求解只允许改变 symmetry boundary。其余设置均从本次运行实际解析的
`scripts/parameter.json` 和现有 finite 配置读取，并在四轮之间保持相同：

- finite geometry 与 cavity/cladding compact geometry；
- cavity/cladding layers 和 shift pair；
- slab、air、PML、材料与 mode type；
- mesh size 和 mesh construction；
- `center_frequency_thz` 对应的 eigensolver shift；
- `finite_eigenmode_count`；
- far-field、finite lattice Fourier、评分及绘图数值参数。

four-modes 入口不再从 unit-cell `selected_bands.csv` 为四轮提供不同 Gamma shift，也不再
以 `p2/p1/d1/d2` 作为四次求解的目标。当前参数为 `cavity_layers=20`、
`cladding_layers=20`、`mesh_size=5`、`finite_eigenmode_count=15`、
`center_frequency_thz=197.98`、hex geometry 和一个 `(x=0, y=0)` shift pair；实际运行前
必须重新解析参数文件，不能把这些当前值硬编码进程序。

## Symmetry 内部映射

执行顺序固定为 symmetry ID `1 -> 2 -> 3 -> 4`：

| symmetry ID | x-axis | y-axis | 内部存储后缀 |
| ---: | --- | --- | --- |
| 1 | PEC | PMC | `_1_xPEC_yPMC` |
| 2 | PMC | PEC | `_2_xPMC_yPEC` |
| 3 | PEC | PEC | `_3_xPEC_yPEC` |
| 4 | PMC | PMC | `_4_xPMC_yPMC` |

ID 和 PEC/PMC 组合属于内部计算身份：用于边界切换、solution tag、mode 文件夹消歧、
缓存身份和审计元数据。它们不得进入 f-Q 图或其他科研结果图的颜色编码、图例、标题和
点标注。mode 文件夹后缀是用户指定的技术存储键，不作为物理模式名称展示。

## Series 与输出布局

series 根目录由当前 finite-quarter 系列标签派生，只把前缀改为
`finite_quarter_all`，从而原样保留 cell/mesh、square、非 uniform shift profile 和
ellipse-shift 等稳定后缀。当前参数应得到：

```text
scripts/.out/finite_cavity/
  finite_quarter_all_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
```

每个实际 shift pair 继续使用现有 directional case 命名。当前示例布局为：

```text
finite_quarter_all_20-20_.../
  shiftx0.000_shifty0.000/
    00_model/
      finite_quarter_all.mph
      comsol_progress.log
      <existing geometry figures>
    01_results/
      mode0_1_xPEC_yPMC/
      ...
      mode14_1_xPEC_yPMC/
      mode0_2_xPMC_yPEC/
      ...
      mode14_4_xPMC_yPMC/
    10_overview/
      eigenfrequencies.csv
      f_Q.png
      <existing merged overview tables and figures>
    80_logs/
    99_config/
      config.json
      run_summary.json
      <existing postprocess configs>
```

同一 symmetry sector 内的 `mode_idx` 继续对应 COMSOL 本轮 solution 的 looplevel，不做
全局重编号。跨 sector 唯一键为完整 mode 目录名，例如
`mode7_3_xPEC_yPEC`。对当前 `finite_eigenmode_count=15`，统一本征频率表的预期行数为
`4 * 15 = 60`。

## 单 Model、单 mesh、四 solution 副本生命周期

每个 shift pair 使用如下生命周期：

```text
创建一个 SimulationRun / COMSOL Model
  -> 构建一次 quarter geometry、materials、physics、study
  -> 构建并运行一次 mesh
  -> ID 1: 激活 PEC/PMC -> solve working solution
       -> 复制为持久 solution 副本 -> 校验副本 -> 导出本轮结果
  -> 清理 working solution data
  -> ID 2: 激活 PMC/PEC -> solve
       -> 复制副本 -> 校验 -> 导出
  -> 清理 working solution data
  -> ID 3: 激活 PEC/PEC -> solve
       -> 复制副本 -> 校验 -> 导出
  -> 清理 working solution data
  -> ID 4: 激活 PMC/PMC -> solve
       -> 复制副本 -> 校验 -> 导出
  -> 保存包含四个 solution 副本的同一个 finite_quarter_all.mph
```

实现时为每条 symmetry axis 预建 PEC 和 PMC physics feature，并通过 active state 切换；
两者始终绑定现有稳定 named selection。切换边界不得重建 geometry、materials 或 mesh。

工作 solution 保留现有 study/solver 连接，用于四次顺序求解。每轮成功后复制为稳定、
互不覆盖的 solution tag；tag 与人类可读 label 可含内部 symmetry 标识。复制完成后必须：

1. 确认副本 tag 存在且不覆盖前一轮；
2. 从副本读取本征频率/Q，逐项与刚完成的工作 solution 比较；
3. 至少抽查一个有效 mode 的场值可从副本 dataset 读取；
4. 保存同一个 MPH 路径；
5. 只有上述检查通过，才允许清理 working solution data 并进入下一轮。

最终 MPH 必须包含四个可独立访问的持久 solution 副本，而不是只保留第四轮解。配置和
run summary 记录内部 solution tag、dataset tag、边界组合、完成时间与验证结果。若 COMSOL
6.3 的 Java API 对 solution copy 的实际调用形式与预期不同，先用最小 geometry smoke
model 验证 copy、保存、重新加载和读取；验证通过前不得启动正式 20-20 计算。

执行统计写入 `run_summary.json`，正式完成时必须满足：

```text
model_build_count = 1
mesh_build_count = 1
solve_count = 4
solution_copy_count = 4
```

## 统一数据集与最终化

每轮仍在 case-local `.staging/` 的内部 symmetry 分区写出数据，避免四轮重复
`mode_idx=0..N-1` 相互覆盖。四轮全部完成后才统一最终化：

1. 合并四份 `eigenfrequencies.csv`，保留全部本征值和有效性判断；
2. 内部增加 `mode_uid`、`local_mode_idx`、`symmetry_id`、`x_boundary`、
   `y_boundary`、`solution_tag` 和 `dataset_tag`；
3. 使用 `mode_uid` 将场、air-plane、far-field、lattice Fourier 和 score 路由到带后缀的
   mode 目录；
4. 合并 overview 表，并重写其中的 mode 输出路径；
5. 检查未知 staging 文件、重复目的路径、缺失模式和行数；
6. 无冲突后一次性移动到正式目录并删除空 staging。

完整数据集必须保留每个 solver eigenpair。物理有效性规则保持现状：所有模式都有本征频率
和有效性证据；有效模式执行完整场、far-field、Fourier 和评分处理；无效模式不伪造场或
后处理结果。four-modes 工作流不再执行“每个 sector 只选一个 fundamental mode”的单模
筛选。

为避免维护两套搬运逻辑，应把现有 finite finalizer 的 mode 目录解析和 move-plan 部分
最小参数化，使普通 finite/quarter 仍以原默认命名运行，而 all-symmetry workflow 传入
mode directory resolver 和四个 staging partition。

## f-Q 图合同

- 输出文件：`10_overview/f_Q.png`。
- 唯一数据源：统一 `10_overview/eigenfrequencies.csv`。
- 横轴：`Frequency (THz)`，取本征频率实部 `re`。
- 纵轴：`Q`，取 `q`。
- 图形：不连线的散点图。
- x/y 均显式设置为 linear scale，不使用对数变换。
- 不归一化；坐标范围根据本次有限数据自动留出小幅边距。
- 所有具有有限 frequency 和 Q 的模式使用完全相同的颜色、marker 和尺寸。
- 不按 symmetry ID、PEC/PMC、有效性或 mode 类型分色、分组或分面。
- 不显示 symmetry/PEC/PMC 图例，不在点旁标注内部 mode 身份。
- 无冗余标题；刻度向内，使用克制的科研绘图样式。
- frequency 或 Q 非有限的行仍保留在 CSV，并在内部 summary 记录未绘制计数。

该图只表达完整离散本征谱中的 frequency-Q 分布，不承担 symmetry 身份识别职能。

## 程序修改范围

预计修改以下现有文件，最终以最小可行 diff 为准：

- `scripts/run_main/run_finite_quarter_four_modes.py`
  - 改为 all-symmetry batch orchestration；
  - 派生 `finite_quarter_all_...` series；
  - 管理单 Model、四轮状态、统一 staging 和 f-Q 图。
- `scripts/run_main/run_finite_quarter.py`
  - 拆出 caller-owned `SimulationRun` 的单 sector 求解与导出步骤；
  - 保留原单 symmetry 包装入口。
- `comsol_workflow/simulation_utils.py`
  - 增加有限 quarter 边界 active-state 切换、working solution 重求和 solution copy/readback
    所需的最小公共能力。
- `scripts/run_main/run_finite.py`
  - 最小参数化现有 finalizer，以支持复合 mode 目录名和多 staging partition；默认行为不变。
- `tests/test_run_finite_quarter_four_modes.py` 及相关 focused tests
  - 更新 orchestration、命名、合并与绘图合同测试。
- README/`scripts/README.md`
  - 更新入口目的、命令、单 MPH 与统一输出布局说明。

不新增第二套 geometry、field reconstruction、far-field 或 Fourier 实现。

## 昂贵计算前预检

正式运行前必须从实际参数文件重新报告并检查：

- cavity/cladding compact geometry；
- cavity/cladding layers；
- finite geometry 与 shift geometry/profile；
- x/y shift pair；
- mesh size、中心频率和每轮 eigenmode count；
- symmetry 执行顺序与预计总模式行数；
- series/case 路径和唯一 MPH 路径；
- 预期 model/mesh/solve/solution-copy 次数；
- 当前 Python/COMSOL 进程；
- COMSOL 6.3 PATH、license；
- 目标目录、`99_config/config.json` 和 partial staging 冲突。

`--preflight-only` 不导入 `mph`、不启动 COMSOL，也不写正式结果。

## 失败、恢复与完成判据

- 任一轮 solve、solution copy、readback、导出或 finalization 失败时立即停止后续工作，并在
  `run_summary.json` 记录失败阶段和内部 symmetry 身份。
- 未通过副本校验时不得清理当前工作解。
- 未完成四个副本时不得把 batch 标记为 complete。
- 正式目录已有结果时默认拒绝覆盖；不通过删除旧结果绕过身份检查。
- `.staging` 中的 partial 数据不等于正式结果；恢复行为只能在配置身份完全相同且已完成
  solution 副本可重新加载验证时复用。

正式完成必须同时满足：

- COMSOL 进程正常退出且无未解释错误；
- 一个 MPH 可重新加载，并能枚举和读取四个 solution 副本；
- model/mesh/solve/copy 计数分别为 `1/1/4/4`；
- 聚合本征频率表行数为 `4 * finite_eigenmode_count`；
- 每行 `mode_uid` 唯一且与 mode 目录一致；
- `01_results` 不含四个 symmetry case 子目录；
- 必需 CSV、Parquet、图像、MPH、config、summary 和 `f_Q.png` 完整；
- f-Q 图为双线性坐标，点样式不显化 symmetry 或 PEC/PMC；
- `.staging` 已完成清理。

## 验证方案

### 无 COMSOL 自动测试

- series 名精确以 `finite_quarter_all_<layers>...` 开头并保留所有现有稳定后缀；
- symmetry 顺序固定为 `1,2,3,4`，内部边界映射正确；
- mock runner 验证 geometry build 1 次、mesh 1 次、solve 4 次、solution copy 4 次；
- 前一 solution copy 校验完成前不能清理 working solution 或进入下一轮；
- 四轮接收相同的 eigensolver shift、eigenmode count、mesh 和后处理参数；
- `4N` 行聚合、`mode_uid` 唯一、mode 目录后缀和 overview 路径重写正确；
- 统一结果中不存在 `symmetry_*` case 子目录；
- f-Q 图 x/y scale 都为 `linear`，点数等于有限 frequency/Q 行数，只有一种点样式，且
  无 symmetry/PEC/PMC 图例或文本；
- 普通 full finite 和单 symmetry quarter 的既有 finalizer 行为保持不变。

### COMSOL 6.3 最小验证

先用小 geometry 和小 eigenmode count 验证：

1. 预建的 PEC/PMC feature 能通过 active state 正确切换；
2. 第二至第四轮不重建 geometry 和 mesh；
3. solution copy 在清理 working solution 后仍可读取；
4. MPH 保存并重新加载后四个副本仍存在，频率/Q 和抽查场值不变；
5. 四个副本绑定的边界身份和数据不串轮。

该 smoke 验证通过并获得用户对正式计算的明确授权后，才允许运行 20-20 正式任务。

### 静态与回归检查

运行相关 finite/quarter/finalizer/field/far-field 测试、相关入口 `py_compile`、
`git diff --check`，并根据实际改动范围运行完整 pytest。导出后检查实际 `f_Q.png`，确认
线性坐标、刻度、边距和点样式符合合同。

## 回滚方案

- 恢复 four-modes 入口原 target-driven orchestration；
- 撤回 SimulationRun 的 symmetry toggle/solution-copy 辅助能力；
- 撤回 finalizer 的可选复合命名参数，普通默认调用不受影响；
- 删除本任务新增测试和文档更新。

由于实施阶段默认不覆盖既有正式结果，代码回滚不需要迁移或删除历史数据。

## 执行记录

### 实际改动

- `run_finite_quarter_four_modes.py` 已改为固定 symmetry ID `1 -> 2 -> 3 -> 4` 的
  all-symmetry 编排，不再读取 unit-cell Gamma target，也不再执行每 sector 单模筛选。
- 每个 shift pair 只创建一个 `SimulationRun`，只构建一次 COMSOL geometry 和 mesh；
  第一轮通过原 `run_simulation()` 求解，后三轮只切换预建 PEC/PMC feature 并重求
  `sol1`。
- 每轮将工作解复制为 `solsym1..4`，建立 `dsetsym1..4`，完成频率逐项比对和有效模式场
  抽查后保存同一个 `00_model/finite_quarter_all.mph`，随后才进入下一轮。
- 四个内部 symmetry partition 继续复用现有 quarter 全模式导出、far-field、finite
  lattice Fourier、评分和 finite finalizer；最终统一合并到一个 `01_results`，mode
  目录使用 `modeN_<id>_xPEC|PMC_yPEC|PMC` 存储身份。合并前检查 `4N` 唯一身份、未知
  文件、目标冲突和共享索引图字节一致性，内部配置归档到
  `99_config/internal_symmetry/<id>/`。
- 新增 `10_overview/f_Q.png`：从聚合 `eigenfrequencies.csv` 读取 frequency/Q，使用
  单一 scatter 样式和显式 linear x/y scale；图中无 symmetry/PEC/PMC 分色、图例、
  标题或点标注。非有限行保留在 CSV，并在 run summary 记录未绘制数量。
- series 命名改为 `finite_quarter_all_<cavity_layers>-<cladding_layers>...`，保留 cell、
  mesh、shift profile、square 和 ellipse 等现有稳定后缀。`--output-dir` 仍被限制在
  `scripts/.out/finite_cavity/` 下。
- `SimulationRun` 新增 finite-quarter boundary active-state 切换、只重求工作 solution、
  solution copy/dataset 建立能力；普通 single-quarter 会同时预建两套边界 feature，但
  仍按原 config 只激活一套，原输出命名和 finalizer 行为不变。
- README 与 `scripts/README.md` 已同步更新。实际最小实现未修改 `run_finite.py` 或
  `run_finite_quarter.py` 的公共 API，而是直接复用其现有接口。

### 测试与 COMSOL 验证

- `py_compile`：核心入口、`simulation_utils.py` 及相关测试文件通过。
- 聚焦测试：`22 passed`，覆盖 series/mode 命名、solution-copy 校验顺序、四 partition
  合并与路径重写、linear f-Q 绘图合同、单 model/mesh + 四 solve/copy mock 时序、普通
  single-quarter mode 命名兼容，以及底层边界切换/重求/复制 API。
- finite/quarter/SimulationRun 回归：`60 passed`。
- 全仓 pytest：`522 passed, 3 warnings`；3 个 warning 均来自既有 strip 测试返回非
  `None`，与本次改动无关。
- 使用当前实际参数运行 `--preflight-only` 成功：hex、20-20、mesh5、每 sector 15
  modes、预计总计 60 modes；目标 series 为
  `finite_quarter_all_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5`，quarter
  hole 数为 7061。该命令未启动 COMSOL，也未写正式结果。
- preflight 检测到当前 PATH 中没有 `comsol`，且已有 `comsol.exe` PID 22624；因此
  未执行真实 COMSOL smoke 或正式计算。license 也未在不启动 COMSOL 的 preflight 中
  消耗或验证。
- 用临时样例实际生成过 `f_Q.png`，文件写出成功；本轮 `view_image` 因 Windows sandbox
  helper 的 `0xc0000142` 故障无法载入，临时目录已验证路径后删除。

### 遗留问题与最终状态

代码与无 COMSOL 验证已完成，最终实现状态为 `implemented_pending_comsol_smoke`。
下一步应先在没有其他项目 COMSOL 进程的会话中配置 COMSOL 6.3 PATH，用小 geometry/
小 eigenmode count 验证四个 solution copy 在清理 `sol1`、保存并重新加载 MPH 后仍可
读取频率/Q和场。只有该 smoke 通过且用户明确批准后，才运行正式 20-20 计算。

### 2026-09-06 复用模型实施与启动状态

- `parameter.json` 的权威求解字段已按本次指令更新为
  `center_frequency_thz=203.6`、`finite_eigenmode_count=6`；
  `center_frequency_source.frequency_thz` 同为 203.6。
- all-symmetry 入口新增 `--reuse-meshed-model`。它校验标准结果配置、compact geometry、
  20--20 层数、mesh 5、shift geometry/profile、shift pair、quarter hole 数、源/目标路径隔离，
  然后只读加载源 MPH。
- `SimulationRun.prepare_reused_finite_quarter_model()` 校验既有 study、solution、dataset 和
  numerical tag，禁用旧 generic x/y symmetry feature，补齐四个标准 PEC/PMC feature，并将
  study `eig` 与 solver `sol1/e1` 同时更新为 203.6 THz、6 modes。
- 复用分支不调用 geometry builder 或 mesh builder；第一轮到第四轮都只重求 `sol1`，每轮
  复制和验证完成后保存新的目标 MPH。完成计数合同为 `0/0/1/1/4/4`。
- 新增并发安全门：发现已有 COMSOL/mphserver 时默认拒绝启动；只有显式
  `--allow-concurrent-comsol` 才允许共享本机资源。
- 实际 preflight 已通过：hex、20--20、mesh 5、203.6 THz、每 sector 6 modes、预计 24 行、
  quarter holes 7061、shift `(0,0)`、COMSOL 6.3 可在 PATH 找到。源配置 SHA-256 为
  `8b25ce0786b8dc889f8e3852a0267ee6dd15c952037ba5acd88963531ca38ed2`；四个旧 sector 的
  geometry identity hash 全部为
  `6508f3706559be269121c609158e4bd3e0a66f3542bf91998764ed85f8747f33`。
- 目标目录为
  `finite_quarter_all_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f203.600_e6_20260906`；
  preflight 没有创建它，也没有改写源结果。
- 验证结果：复用聚焦测试 `24 passed`，finite/quarter 回归 `55 passed`，全仓
  `524 passed, 3 warnings`；3 个 warning 是既有 strip 测试返回非 `None`，与本任务无关；
  `py_compile` 与 `git diff --check` 通过。
- 正式启动当前暂停：检测到 PID 22624 的既有 `comsol.exe`，自 2026-08-18 启动，命令行
  正在打开 `E:\MA\240127 magic-angle_3D single cell 4.41 cellwaiGaAs D=123 R0=267 without Au without SiO2.mph`，
  工作集约 2.7 GB。该进程来源于其他模型，未擅自终止，也未并发启动新求解。

用户已明确允许与既有 COMSOL 会话并发。正式任务于 2026-09-06 18:49:55 启动；本任务
Python PID 为 22456，子进程 `comsolmphserver.exe` PID 为 35328。COMSOL progress 已进入
symmetry 1 的“编译方程：特征频率”，未出现 geometry 或 computational mesh 重建阶段，
确认当前正在使用源 MPH 的既有 mesh。启动时目标 `run_summary.json` 为 `status=running`。

当前状态：`running_first_symmetry_solve`。

### 2026-09-06 solution-copy 故障与修复

- 首次正式运行的 symmetry 1 在 2026-09-06 21:42:43 完成，solver 用时 9138 s；随后在
  `jmodel.sol().copy("solsym1", "sol1")` 处收到 COMSOL 6.3“unsupported operation”。
  运行以 `status=failed` 退出，计数为 reused model/mesh `1/1`、solve `1`、copy `0`；
  未生成目标 MPH 或 mode results，源 MPH 未改变。
- 本机 COMSOL 6.3 Java API 文档明确规定 solution data 副本接口为
  `model.sol("sol1").copySolution("solsym1")`；solver-list `copy()` 不支持复制 solution。
- `SimulationRun.copy_solution()` 已改为调用 working solver 的 `copySolution()`，然后从
  新 solver tag 设置 label 并创建对应 Solution dataset。其余编排和数据格式不变。
- 真实 COMSOL 小模型 smoke 已通过：`solcopy` 与 `dsetcopy` 创建成功；清理并重求 `sol1`
  后副本数值不变；保存 5.2 MB MPH 并重载后，副本仍可读。临时模型和临时脚本均已清理。
- 修复后验证：聚焦测试 `24 passed`，全仓 `524 passed, 3 warnings`，`py_compile` 和
  `git diff --check` 通过。warning 仍为既有 strip 测试返回非 `None`。

当前状态：`solution_copy_fixed_and_smoke_verified_ready_to_retry`。

修复后的正式任务已于 2026-09-06 22:51:51 使用非冲突目录
`finite_quarter_all_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f203.600_e6_20260906_retry1`
重新启动。失败目录原样保留，不覆盖也不删除；retry1 的 Python PID 为 36468，
`comsolmphserver.exe` PID 为 21188，当前正在加载源 MPH。

当前状态：`retry1_running`。

### 2026-09-07 正式运行完成记录

- 203.6 THz、`finite_eigenmode_count=6` 的四轮 COMSOL 求解全部完成；每轮实际返回 12 个本征解，
  `solve_count=4`、`solution_copy_count=4`。每轮均先复制为 `solsym1..4` / `dsetsym1..4` 并保存同一个
  `finite_quarter_all.mph`，随后才进入下一轮。最终 MPH 为 21,117,055,898 bytes。
- 恢复后逐项验证四个保存解的本征频率，`pending=[]`，没有再次触发 eigensolve。symmetry 1/2/3 分别有
  6/2/3 个有效模式；symmetry 4 的 12 个本征解均未通过有效性判据，因此保留全部模式和频率/Q，跳过
  不可定义的 lattice Fourier 与评分。
- 后处理首次失败的根因是 Windows `MAX_PATH`：内部 staging 的
  `p_subspace_mode_decomposition.csv` 绝对路径为 267 字符。恢复时通过临时短盘符完成 11 个有效模式的
  并行后处理；代码随后把内部 partition 缩短为 `s1..s4`，最终 mode 后缀和公开布局不变。
- 最终统一数据集包含 48 个目录，命名均为 `modeN_<id>_xPEC|PMC_yPEC|PMC`；聚合
  `eigenfrequencies.csv` 为 48 行、48 个唯一 `mode_uid`。far-field、finite-lattice Fourier 和 mode-score
  聚合表各有 11 行。
- `10_overview/f_Q.png` 使用 48 个有限 frequency/Q 点，x/y 均为线性坐标，使用单一散点样式，
  不显示 symmetry ID 或 PEC/PMC 图例。PNG 可完整解码，尺寸为 1144×880。
- `run_summary.json` 已为 `status=complete`，`.staging` 已清理，本任务 Python/mphserver 均已退出；
  用户原有 COMSOL PID 22624 未改动。
- 收尾验证：聚焦测试 `25 passed`，相关入口 `py_compile` 与 `git diff --check` 通过。`view_image` 因
  本机 sandbox helper 的 `0xC0000142` 无法加载，改用 Pillow 完成 PNG 解码、尺寸和像素范围检查。

当前状态：`complete`。
