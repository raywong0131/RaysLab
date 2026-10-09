# Unit-cell 二维 Γ 点扫描、偏振与 winding 分析实施计划
- 文档日期：2026-08-21
- 当前状态：**功能与小网格实机验证已完成，等待正式 19x19 运行前确认**
- 目标入口：`scripts/run_main/run_unit_cell_2d.py`、`scripts/analysis/unit_cell_2D.py`
- 正式输出根：`scripts/.out/unit_cell_band/unit_cell_2D/`
- 外部参考：Workspace 中的独立项目 `comsol-data-analysis`

## 1. 背景与目标

在现有六角 unit-cell 本机 COMSOL 工作流上，新增 Γ 点附近的二维 Cartesian k 网格扫描，持续追踪指定的 Γ 点模态（首个正式目标为 `p2@Gamma`），计算其远场零级偏振系数

\[
c_x(\mathbf{k})=\frac{1}{\sqrt A}\int_A E_x(\mathbf r) e^{i\mathbf k\cdot\mathbf r}\,\mathrm dA,
\qquad
c_y(\mathbf{k})=\frac{1}{\sqrt A}\int_A E_y(\mathbf r) e^{i\mathbf k\cdot\mathbf r}\,\mathrm dA,
\]

并从同一份结构化结果生成：

1. `c(k)` 强度二维 Color Map；
2. 偏振方向二维 Color Map；
3. 全二维网格偏振矢量图；
4. winding loop 与展开偏振角的两面板图；
5. shrinking-square-loop winding 稳定性图和数值摘要；
6. `comsol-data-analysis` 中 period-2D polarization/winding 工作流现有角度定义与绘图能力的兼容输出。

本计划只定义实施方案。本阶段不修改计算代码、不启动 COMSOL、不生成正式结果。

## 2. 外部参考边界

`comsol-data-analysis` 只作为独立参考文献和数值回归对照，保持以下硬边界：

- 不复制整个项目或其 package 到本仓库；
- 不在 `pyproject.toml` 中添加对它的依赖；
- 不通过 `sys.path`、editable install、相对路径或运行时 import 调用它；
- 不改动其源码、配置和既有结果；
- 可在开发验证时只读其已保存 CSV，比较表结构、角度定义和 winding 结果；
- 正式运行与测试不得要求参考项目存在。

允许按本项目结构重新实现并测试的参考思想包括：批量读取 `Hz/Ex/Ey`、`cx/cy` 投影、`real`/`imag`/`complex-axis`/`stokes-axis` 角度、square-loop winding 和 shrinking-loop 扫描。实现需写成当前项目自己的小型、可测试模块，不保留跨项目耦合。

## 3. 需求与仓库约束对照

### 3.1 入口拆分

用户指定的 `scripts/analysis/unit_cell_2D.py` 必须存在，但 `scripts/analysis/` 按仓库约定不得启动 COMSOL。因此采用两个入口：

- `scripts/run_main/run_unit_cell_2d.py`：从 `scripts/parameter.json` 读取全部输入，其中几何参数严格取 `unit_cell_2d.b0_nm/eta/zeta`，然后构建/复用 COMSOL 模型、执行二维 k 扫描并写原始数据；求解和数据聚合成功后，默认调用同一个纯 Python 分析模块生成全部图与 winding 摘要；
- `scripts/analysis/unit_cell_2D.py`：只读取一个既有 series 目录中的表格和 metadata，重算 winding、重绘全部图，不构造 `SimulationRun`，也不 import/启动 `mph`。

两个入口不得各自维护一套 winding 或绘图实现。

### 3.2 输出根

仓库规定 `scripts/.out/` 的一级目录只能是 `unit_cell_band/`、`strip_1d/`、`finite_cavity/`。因此用户要求的新 `unit_cell_2D` 文件夹建立为：

```text
scripts/.out/unit_cell_band/unit_cell_2D/
```

不创建 `scripts/.out/unit_cell_2D/` 这一第四个一级输出目录。

### 3.3 参数来源

- `run_unit_cell_2d.py` 的运行参数统一从 `scripts/parameter.json` 读取，不在 runner 中保留另一套会漂移的数值默认值；
- 二维扫描的 unit-cell 几何参数严格读取 `scripts/parameter.json` 的 `unit_cell_2d.b0_nm/eta/zeta`，不读取顶层 `cavity` 或 `cladding`，也不接受代码内另行覆盖；顶层 `cavity` 继续服务既有工作流。为兼容不运行本二维 workflow 的旧版参数文件，仅当整个 `unit_cell_2d` 对象缺失时，解析层以顶层 `cavity` 构造默认二维配置；对象一旦存在，三个几何字段必须完整；
- mesh 等级、unit-cell eigenmode count 和既有求解器参数继续读取 `scripts/parameter.json`；
- 二维 q 轴、目标 band、追踪门限、归一化版本和绘图开关放入 `scripts/parameter.json` 的专属 `unit_cell_2d` 对象，并由参数解析层完成类型、范围和互斥校验；
- unit-cell eigensolver 基准继续固定为 `c_const/1.55[um]`；
- runner 只允许非数值型的执行控制参数（例如显式 `--skip-analysis`）；所有影响结果的参数仍以 `parameter.json` 和运行时快照 `99_config/config.json` 为准；
- 本 workflow 不在运行中隐式改写 `scripts/parameter.json`，也不改变 unit-cell 主流程成功后更新共享中心频率的既有规则；新二维扫描本身不更新共享中心频率。

### 3.4 默认 Analysis 行为

- `run_unit_cell_2d.py` 在全部要求的 k 点求解/恢复并成功写出聚合表后，默认直接执行 `unit_cell_2d_analysis`，生成本计划列出的全部 Color Map、偏振矢量图、winding 图和摘要；
- 自动 Analysis 与 `scripts/analysis/unit_cell_2D.py` 调用同一实现和同一绘图合同；
- Analysis 失败时 runner 返回非零状态并保留已经完成的求解缓存，不能把“COMSOL 求解成功但分析失败”报告成完整成功；
- 仅为排查或分阶段恢复提供显式 `--skip-analysis`，其默认值必须为 `False`，并记录到 `run_summary.json`。

## 4. 范围与非范围

### 4.1 本轮实施范围

- cavity unit-cell、当前 hex unit-cell runner；
- Γ 附近二维 k 网格 eigensolve；
- 可恢复的单 Model 重用扫描；
- 目标模态的二维 band tracking；
- `Hz` 追踪场和 `Ex/Ey` 偏振投影；
- 参考兼容表、聚合表、Color Map、矢量图与 winding 分析；
- 聚焦单元测试、纯 Python 集成测试、小网格实机验证和正式运行前检查。

### 4.2 非范围

- 不把 `comsol-data-analysis` 合并为本仓库子项目；
- 不改 strip、finite、finite-quarter 工作流；
- 不重新引入 Slurm、SSH、远程 scratch 或后台 launcher；
- 不删除、移动或覆盖任何既有正式结果；
- 不在计划获确认前运行 19×19 正式 COMSOL 扫描；
- 不把逐点未归一化的 COMSOL 本征矢量幅值直接宣称为可跨 k 比较的物理强度。

## 5. 可复用基础与拟修改文件

### 5.1 直接复用

- `scripts/run_main/run_band_pair.py`
  - `ReusableSimulationRun` 的固定几何单 Model 重用；
  - k 参数更新、重新求解、频率和中心 `Hz` 缓存；
  - 频率门限与复数 `Hz` overlap 的 branch selection 思路。
- `comsol_workflow/band_connector.py`
  - `assign_mode_candidates` 的 Hungarian 全局唯一指派。
- `comsol_workflow/simulation_spatial/hexagon_unit_cell.py`
  - 现有六角 unit-cell 构建、2D 场读取与 `compute_polarization()` 数值定义。
- `scripts/run_main/run_band_solo.py`、`scripts/run_main/run_unit_cell_band.py`
  - cavity-only unit-cell 配置、Γ 点模式分类和 `p2` 命名选择。

### 5.2 计划新增

- `scripts/run_main/run_unit_cell_2d.py`
  - COMSOL 二维扫描主入口。
- `scripts/analysis/unit_cell_2D.py`
  - 用户指定的纯分析/重绘入口。
- `comsol_workflow/unit_cell_2d_analysis.py`
  - 表 schema、Stokes 量、角度、loop、winding、图形生成和结果 QA；禁止依赖 `mph`。
- `tests/test_unit_cell_2d_analysis.py`
  - 解析场、相位不变性、winding、非均匀网格和绘图合同测试。
- `tests/test_run_unit_cell_2d.py`
  - 配置、命名、扫描调度、缓存和不启动 COMSOL 的入口测试。

### 5.3 计划修改

- `comsol_workflow/simulation_spatial/hexagon_unit_cell.py`
  - 增加一次 `Interp` 批量获取复数 `Hz/Ex/Ey` 及归一化积分所需量的接口；保留旧接口兼容。
- `comsol_workflow/band_connector.py` 或新 runner 的纯函数区
  - 在现有 Hungarian 指派之上增加二维多前驱 score 聚合；若逻辑超过 runner 可读性，则落到独立的纯 Python 模块，避免污染 COMSOL runner。
- 相关现有测试
  - 为批量 field extraction 的表达式顺序、实虚部重组和旧 `compute_polarization()` 对照补测试。
- 项目输出说明文档（若当前仓库已有统一结果说明入口）
  - 登记 `unit_cell_2D` 的嵌套输出位置和重绘命令，不另造重复文档体系。

不在 package 根 `__init__.py` 中引入 `mph` 或重型模块。

## 6. 正式数据流

```text
scripts/parameter.json
  -> 读取 unit_cell_2d 专属几何与扫描配置 + mesh / eigenmode count
  -> 构建一次 ReusableSimulationRun
  -> 求解 Gamma 并按现有模式成分选择 p2
  -> 按由内向外的二维 shell 顺序更新 kx/ky 并 eigensolve
  -> 频率门限 + 多内侧邻点复数 Hz overlap + Hungarian 唯一指派
  -> 批量读取目标 mode 的 Hz/Ex/Ey 和归一化量
  -> 写 point cache、tracking evidence、polarization table
  -> runner 默认调用纯 Python unit_cell_2d_analysis
  -> Color Maps / vector maps / winding figures / CSV / JSON
```

COMSOL 模型在一次进程中只构建一次。恢复运行时仍先按当前配置构建一次模型，然后跳过 cache identity 完整且校验通过的 k 点；不得仅凭文件存在即视为有效缓存。

## 7. 二维 k 网格

### 7.1 默认正式网格

首个正式配置复用参考数据的 Γ 邻域多尺度轴：

```text
q = k/G = [
  -0.10, -0.08, -0.06, -0.04, -0.02, -0.01,
  -0.001, -0.0001, -0.00001,
   0,
   0.00001, 0.0001, 0.001,
   0.01, 0.02, 0.04, 0.06, 0.08, 0.10
]
```

取 `qx × qy`，共 `19×19=361` 点。轴在文件中严格递增；`G` 使用当前 unit-cell reciprocal-lattice 定义，不从参考项目复制常数。每行同时保存：

- `qx_over_G`、`qy_over_G`；
- 传给 COMSOL 的物理 `kx`、`ky` 及单位；
- reciprocal-scale `G` 及其来源版本。

正式实现须用当前 reciprocal helper 计算并测试 `q -> k`，禁止散落硬编码。

### 7.2 求解顺序

求解顺序与最终表格顺序分离：

1. Γ 点先求解并确定目标 band；
2. 其余点按到 Γ 的离散 shell/半径由小到大；
3. 同一 shell 使用确定性排序，确保恢复运行和测试可重复；
4. 最终表按 `qy`、`qx` 稳定排序，绘图不依赖求解顺序。

小型实机验证使用独立的 3×3 或 5×5 临时网格和非正式 series label，不覆盖正式 19×19 目录。

## 8. Γ 点选择与二维 band tracking

### 8.1 Γ 点锚定

- 使用二维 workflow 专属的 `unit_cell_2d.eigenmode_pair_count` 作为请求本征解对数；当前默认值为 `2`，原样写入 COMSOL `neigs=2`，按成对返回语义预期得到 4 个本征解，不在 Python 侧再次乘 2；顶层 `unit_cell_eigenmode_count` 继续保持原值并仅服务 `unit_cell_band` 系列；
- 复用现有 Γ 点模式成分/命名逻辑选择 `p2`，不得把固定 `solnum` 当作 `p2`；
- 保存 Γ 点所有候选的频率、Q、模式成分、内部 mode index、目标 label 和选择证据；
- 若 `p2` 缺失、并列或不满足既有分类门限，正式扫描停止并给出可诊断错误，不猜测 band。

实现支持目标 label 列表，但首个正式 series 默认只追踪 `p2`。扩展到 `p1/d1/d2` 时复用同一数据流和每-band 输出，不复制脚本。

### 8.2 多前驱二维追踪

参考项目的单最近前驱贪心策略不直接采用。对每个待求点：

1. 从已成功的、更靠近 Γ 的轴向/对角相邻点选择最多若干前驱；
2. 对当前全部 eigenmode candidates，先应用频率差门限；
3. 在相同采样网格上归一化复数 `Hz`，计算相位不敏感 overlap `|<H_a,H_b>|`；
4. 将多个前驱的 overlap 和频率连续性聚合为 candidate score；必须记录各前驱分数，而非只保存最终值；
5. 多 band 配置使用 Hungarian 指派，确保一个 candidate 不被两个 band 占用；
6. 保存 selected score、次优 score、ambiguity gap、频率差、candidate count、前驱坐标和状态；
7. 低于 overlap 门限、gap 过小或所有前驱无效时标记 `missing/ambiguous`，不静默换带；
8. 第一遍结束后允许使用另一方向已完成邻点做一次确定性的 repair pass，修复仍须通过相同门限并保留前后证据。

无效/歧义点在聚合表中保留，图中 mask；不得通过插值伪造 `c(k)` 或 winding 输入。winding loop 只使用所有边界点均有效的闭合 loop。

## 9. 场提取、`c(k)` 与归一化

### 9.1 批量提取

为减少 COMSOL/Java 往返，在目标 mode 上用一次批量 `Interp` 获取并重组：

- `Hz` 的 real/imag：用于 band tracking；
- `Ex*exp(i*(kx*x+ky*y))` 的 real/imag；
- `Ey*exp(i*(kx*x+ky*y))` 的 real/imag；
- 归一化积分所需的场强/能量表达式。

表达式符号、采样平面、采样范围、积分权重和 area 与旧 `compute_polarization()` 保持一致，并通过同一小模型/同一 solnum 的数值对照锁定。批量接口失败时应明确报错，不在正式扫描中悄悄退回每表达式多次查询造成不可预期耗时。

### 9.2 两套系数

同时保存：

- `cx_raw`、`cy_raw`：与参考程序一致，仅含 `1/sqrt(A)` area normalization，作为主 `c(k)` 强度 Color Map 的数据；该图用于复刻参考合同，不把绝对幅值宣称为模式能量归一化后的物理强度；
- `cx_norm`、`cy_norm`：除以同一 mode 的正实数场范数，作为方向、Stokes、winding 和跨独立求解稳定性审计数据。正实数缩放不改变偏振方向或 winding。

第一阶段采用与投影同一采样面上的

\[
N_E=\sqrt{\int_A (|E_x|^2+|E_y|^2)\,\mathrm dA},
\qquad
\tilde c_{x,y}=c_{x,y}^{raw}/N_E.
\]

同时把 `N_E` 写入表中。若实机验证证明当前 COMSOL 模型已有稳定、可读取的全 cell electromagnetic-energy normalization，则在正式 19×19 运行前通过版本化配置切换为该物理能量范数，并保留 planar norm 作为 QA 列；不得无记录地改变定义。

所有 complex 系数在 CSV 中拆为 real/imag 列，在 parquet 中可额外保留 complex-friendly 表示。表中必须有 `normalization_kind` 和 `polarization_definition_version`。

### 9.3 派生量

主强度定义为：

\[
C(\mathbf k)=|c_x^{raw}|^2+|c_y^{raw}|^2.
\]

Stokes 定义为：

\[
S_0=|\tilde c_x|^2+|\tilde c_y|^2,\quad
S_1=|\tilde c_x|^2-|\tilde c_y|^2,\quad
S_2=2\operatorname{Re}(\tilde c_x\tilde c_y^*),
\]

\[
S_3=-2\operatorname{Im}(\tilde c_x\tilde c_y^*),\qquad
\psi=\tfrac12\operatorname{atan2}(S_2,S_1)\pmod\pi.
\]

`S3` 的符号 convention 必须写入 metadata 并由测试锁定。主偏振方向使用 `stokes-axis`，因其对 COMSOL 本征模的任意整体复相位不敏感。

兼容实现还提供参考程序的 `real`、`imag`、`complex-axis` 和 `stokes-axis` 四种 `field_angles`；`real/imag` 输出明确标记为 gauge-dependent diagnostic，不用于主结论。

## 10. winding 定义与质量控制

### 10.1 Square loop

- 以 Γ 为中心，在非均匀 Cartesian 网格上按各正半宽构造逆时针 square loop；
- 每个角点只出现一次，闭合差分单独补首点；
- `stokes-axis` 按其 `π` 周期进行 period-aware unwrap；
- winding 按完整 loop 的累计角变化除以 `2π`；
- 保存 raw winding、nearest integer/half-integer diagnostic、residual、point count、最小 `C`、最小有效系数范数和 loop 完整性。

实现以解析场 `c=(kx,ky)` 的 `+1` 和反向场的 `-1` 作为符号/方向基准，避免因 loop 顺序或坐标轴翻转改变符号。

### 10.2 Shrinking-loop scan

对网格中所有可用 half-width 从外向内扫描。仅当 loop 边界所有追踪点有效、且最小系数范数高于配置阈值时计为有效。主 winding 报告需同时满足：

- 至少若干相邻尺度给出一致 winding；
- residual 在门限内；
- 不依赖单个接近零或 ambiguous 的边界点。

若不同尺度不一致，输出 `unstable`，保留全部曲线和证据，不选择“看起来正确”的单个 loop。

## 11. 绘图合同

所有图从聚合表生成；可比较面板使用一次完成的共享 scale；`kx-ky` 始终等比例。正式实现前测试中锁定文件名、标题、轴范围、colormap、mask、colorbar 和 overlay 数量，导出后必须查看实际 PNG。

### 11.1 `c(k)` 强度 Color Map

- 文件：`ck_intensity_color_map_<band>.png`；
- 数据：`C=|cx_raw|^2+|cy_raw|^2`，严格复用参考程序的 projected-component 定义；
- 非均匀 19×19 网格使用真实 cell edges 的 `pcolormesh`，不得用等间距 `imshow`；
- 主图使用全网格统一 linear scale，附加 `ck_intensity_log_color_map_<band>.png` 展示跨数量级结构；
- log 色标以 `10^n` 表示，floor、clim 和 ticks 从当前有效数据确定并写入 metadata，不继承历史图数值；
- `missing/ambiguous` 点 mask；图题和 colorbar 明确写明 raw projection intensity，不能把它称为模式能量归一化后的绝对物理强度。

### 11.2 偏振方向 Color Map

- 文件：`polarization_angle_color_map_<band>_stokes_axis.png`；
- 数据：`psi`，周期 `π`；
- 使用 cyclic colormap 和固定物理周期的 colorbar；
- `C` 低于相对阈值或追踪无效处 mask，避免在 BIC/近零点显示任意角度；
- 不对角度做普通线性插值。

### 11.3 偏振矢量图

- 文件：`polarization_vector_map_<band>_<angle_method>.png`；
- 全有效网格为克制的灰色 headless vectors；
- 选中 winding loop 为蓝色 vectors，loop 为红线；不额外强调 Gamma 点；
- 向量只表示方向，长度不伪装为强度；强度由 Color Map 单独表达；
- 低强度点不画箭头；坐标等比例，不拉伸拓扑纹理；
- 正式输出为四种 angle method 生成兼容版本，`stokes_axis` 为主版本，其他版本在标题/metadata 中标注 diagnostic。

### 11.4 Winding 两面板图

- 文件：`winding_loop_and_angle_<band>_<angle_method>.png`；
- 左面板：全网格 vectors、选中 loop、loop vectors、Γ；
- 右面板：沿逆时针 loop 的 unwrapped polarization angle、线性趋势和计算 winding；
- 标注 loop half-width、winding、residual 和最小 coefficient norm；
- 不在图中重复堆叠原始点百分数或无信息图例。

### 11.5 Shrinking-loop 稳定性图

- 文件：`winding_loop_scan_<band>_<angle_method>.png`；
- 横轴为 loop half-width `k/G`，纵轴显示 raw winding；
- 同图或紧凑副面板显示 residual/min norm，不能让不同量纲共用无说明的轴；
- 有效、无效、缺点 loop 有清晰但克制的区别；
- 数值源同时写入 `winding_loop_scan_<band>_<angle_method>.csv`。

图中 frequency、Q、band label、k 范围和 normalization 必须从当前 metadata 读取，不硬编码参考结果中的 mode 4、`+1` 或标题数值。

## 12. 输出布局与命名

正式路径：

```text
scripts/.out/unit_cell_band/unit_cell_2D/
  unit_cell_2D_<canonical-cell-and-mesh-label>_kmax0.1_uniform19_band-p2/
    01_results/
      k_points/
        kx<signed-qx>_ky<signed-qy>/
          eigenfrequencies.csv
          tracked_fields.parquet
          point_metadata.json
      polarization_grid.csv
      polarization_grid.parquet
      band_tracking_2d.csv
    10_overview/
      ck_intensity_color_map_p2.png
      ck_intensity_log_color_map_p2.png
      polarization_angle_color_map_p2_stokes_axis.png
      polarization_vector_map_p2_<angle_method>.png
      winding_loop_and_angle_p2_<angle_method>.png
      winding_loop_scan_p2_<angle_method>.png
      winding_loop_scan_p2_<angle_method>.csv
      winding_summary.csv
      winding_summary.json
    99_config/
      config.json
      run_summary.json
```

命名细则：

- `<canonical-cell-and-mesh-label>` 必须调用/复用现有 unit-cell 参数 formatter，包含实际 `unit_cell_2d.(b0_nm, eta, zeta)` 与 mesh，不手拼另一套浮点格式；
- 例示 series label 只说明字段，不作为绕过现有 formatter 的硬编码；
- k 点目录使用固定 signed-decimal formatter，`+0.00000` 和 `-0.00000` 统一为 `0.00000`，避免同一点出现两个目录；
- target band label 使用物理名称 `p2`，内部 solnum 只作为表字段；
- case 根只放上述目录，聚合图只在 `10_overview`；
- 不保存 MPH，除非后续用户明确要求；因此不创建空 `00_model`；
- 已存在同名正式 series 时默认拒绝覆盖；resume 只补齐同一 cache identity 的缺点。

为兼容参考表，可按每个目标 band 额外导出 `polarization_mode_<display-index>.csv`，但它必须由 `polarization_grid` 生成，不能成为第二份权威数据。

## 13. 表结构与 provenance

`polarization_grid` 每个 k 点/target band 一行，至少包含：

- 网格与求解：`qx_over_G, qy_over_G, kx, ky, G, solve_order`；
- identity：`band_label, gamma_mode_label, solnum, mode_index`；
- 频率/Q：`frequency_hz, frequency_thz, q_factor`；
- 原始系数：`cx_raw_re/im, cy_raw_re/im`；
- 归一化：`normalization_kind, normalization_value, cx_norm_re/im, cy_norm_re/im`；
- 派生量：`C, S0, S1, S2, S3, psi_stokes`；
- 追踪 QA：`match_status, selected_score, runner_up_score, ambiguity_gap, frequency_delta, candidate_count`；
- provenance：`polarization_definition_version, field_plane, grid_shape, cache_identity`。

`band_tracking_2d.csv` 另保存每个候选/前驱的细粒度证据，使换带问题可审计。CSV/JSON 原子写入；parquet 是性能副本，CSV/JSON 为长期可读合同。

## 14. 缓存、恢复与失败处理

cache identity 至少包含：

- cell 参数与 geometry 版本；
- mesh 等级/实际 mesh 设置；
- COMSOL 版本、eigensolver 基准和 eigenmode count；
- 完整 q 轴、`G` 定义和 solve-order 版本；
- target band label 与 Γ 分类版本；
- field plane、采样分辨率、表达式和 polarization definition version；
- normalization kind；
- tracking frequency/overlap/ambiguity 门限。

每个点只有在频率表、目标场、metadata 和 checksum/schema 都完整时才标记完成。进程中断留下的临时文件不进入正式聚合表。聚合阶段列出 solved/cached/missing/ambiguous 数量；只要正式要求的 loop 不完整，winding summary 就不能标记成功。

COMSOL license、内存或 Java 异常应保留最后成功点和 progress log，并允许同 identity 恢复，不删除既有点。改变定义版本或网格时必须创建新的非冲突 series，不混用旧 cache。

## 15. 实施阶段与验收门

### 阶段 A：纯 Python 合同和解析测试

1. 实现网格、路径 formatter、表 schema 和配置快照；
2. 实现 Stokes/四种角度、period-aware unwrap、square-loop 与 shrinking-loop；
3. 用解析向量场完成 `+1/-1/0` winding 测试；
4. 实现全部绘图并用合成非均匀网格锁定合同；
5. 确认 `scripts/analysis/unit_cell_2D.py` 的 import/执行路径不加载 `mph`。

验收：不需要 COMSOL 即可从 fixture table 生成所有预期表和 PNG。

### 阶段 B：COMSOL 批量提取与单点对照

1. 增加 batched `Hz/Ex/Ey` extractor；
2. 在同一模型、同一 k、同一 solnum 对照旧 `get_2d_fields()` 和 `compute_polarization()`；
3. 核对复数符号、相位因子、积分权重、cut plane、area 和 normalization；
4. 检查不同 COMSOL 本征矢量整体复缩放下 `cx_norm/cy_norm` 与 Stokes 方向不变。

验收：旧/新 raw polarization 在定义的数值容差内一致，批量字段形状和表达式顺序有测试保护。

### 阶段 C：二维追踪与小网格实机验证

1. 构建 3×3 或 5×5 Γ 邻域临时 series；
2. 检查 Γ `p2`、每点 candidates、多个前驱分数和 Hungarian 唯一性；
3. 人工抽查轴向/对角点的频率连续性和 `Hz` overlap；
4. 从已保存表运行 `scripts/analysis/unit_cell_2D.py`，确认重绘不触发 COMSOL；
5. 查看实际 PNG 的标题、等比例、mask、箭头、loop、colorbar 与对齐；
6. 清理 smoke test 临时输出，不触碰正式结果。

验收：小网格可中断恢复，没有静默换带，所有绘图能从离线表复现。

### 阶段 D：正式 19×19 运行前审批

启动昂贵任务前再次向用户报告并核对：

- 当前 `unit_cell_2d.b0_nm/eta/zeta` 参数与 canonical series label；
- mesh5、2 对本征解（预期 4 个解）、`c_const/1.55[um]`；
- target=`p2@Gamma`；
- 19×19 多尺度 q 轴和实际 `G`；
- normalization kind 与 field plane；
- tracking 门限、缓存/恢复设置；
- 正式输出路径不存在冲突；
- COMSOL 6.3 PATH、license、核数和预计资源；
- 全部后处理开关和预期图清单。

只有用户确认后才运行正式 361 点扫描。

### 阶段 E：正式结果 QA

1. 监控到 mesh 完成和首个完整非 Γ 点，检查 progress/stderr；
2. 完成后检查 361 点完整性、missing/ambiguous、频率连续性和缓存摘要；
3. 检查主 winding 的多尺度稳定性，不以参考项目的 `+1` 作为预设答案；
4. 逐张查看实际 PNG；
5. 运行离线重绘并比较结构化数据/图片 metadata；
6. 回填本文档的实际改动、测试、结果路径、遗留问题和最终状态。

## 16. 测试清单

### 16.1 纯 Python 单元测试

- 默认 q 轴严格递增、包含唯一 Γ、Cartesian product 正好 361 点；
- `q -> k` 使用当前 reciprocal 定义，x/y 轴和单位正确；
- solve order 从 Γ 向外且确定性可复现；
- signed-zero 和 series/k-point 路径 formatter；
- cache identity 对 cell、mesh、网格、band、归一化和定义版本敏感；
- Γ `p2` 选择不依赖固定 solnum；
- 多前驱 score、Hungarian 唯一指派、runner-up gap、missing/ambiguous/repair；
- complex field overlap 对整体复相位不变；
- Stokes angle 对 `c -> c*exp(i phi)` 不变；
- `real/imag/complex-axis/stokes-axis` 与参考定义的只读 fixture 对照；
- 解析 `c=(kx,ky)` 为 `+1`，反向 loop/场符号测试，常向量场为 `0`；
- shrinking loops 的稳定、缺点、近零边界和 residual 案例；
- 非均匀 `pcolormesh` cell edges、等比例坐标、cyclic colormap、mask 和 log ticks；
- 精确文件名、标题、colorbar label/ticks、overlay 数量；
- analysis 入口不实例化 `ReusableSimulationRun/SimulationRun`，不 import `mph`。

### 16.2 COMSOL 聚焦验证

- batched `Hz` 与旧单表达式读取一致；
- batched `cx_raw/cy_raw` 与旧 `compute_polarization()` 一致；
- 同一 mode 的 normalization value 有限且非零；
- 3×3/5×5 model reuse 不重建 geometry；
- resume 会跳过完整点并拒绝 identity 不匹配点。

### 16.3 仓库检查

计划实施后至少运行：

```powershell
uv run python -m pytest -q tests\test_unit_cell_2d_analysis.py tests\test_run_unit_cell_2d.py
uv run python -m pytest -q <受影响的 unit-cell / band-connector 测试>
uv run python -m py_compile scripts\run_main\run_unit_cell_2d.py
uv run python -m py_compile scripts\analysis\unit_cell_2D.py
git diff --check
```

实际 COMSOL smoke test 与 pytest 分开报告，不把环境/license 失败描述成代码单测失败。

## 17. 风险与应对

| 风险 | 影响 | 应对 |
|---|---|---|
| Γ 附近简并/近简并导致换带 | winding 和 Color Map 失真 | Γ 物理 label 锚定、多前驱 overlap、Hungarian、ambiguity gap、repair 和全证据表 |
| COMSOL 本征矢量任意幅值/相位 | raw 强度不能解释为能量归一化绝对强度，`Re(c)` 方向不稳定 | raw 与 normalized 分列；主强度图明确复刻参考 raw projection；主角度/winding 用 normalized Stokes；相位不变性测试 |
| BIC 附近 `|c|≈0` | 角度噪声和伪箭头 | 低强度 mask；loop 记录 min norm；多尺度稳定性判定 |
| 非均匀多尺度网格被按等间距绘制 | 几何关系错误 | 使用真实坐标和 cell edges 的 `pcolormesh`，强制 equal aspect |
| 361 点耗时或中断 | 结果丢失 | 单 Model 重用、point-level 原子缓存、严格 identity、可恢复运行 |
| 外部参考项目变化/缺失 | 正式流程不可复现 | 参考只用于开发期只读对照；本仓库实现和 fixture 自包含 |
| 输出路径违反三级根规范 | 污染 `.out` 布局 | 固定放在 `unit_cell_band/unit_cell_2D`，路径测试锁定 |
| 批量 COMSOL 表达式顺序/复数重组错误 | 系数系统性错误 | 与旧接口逐字段实机对照，表达式 schema 版本化 |

## 18. 回滚方案

- 新 workflow 使用独立入口和独立嵌套输出根，不改变现有主入口行为；
- `hexagon_unit_cell.py` 的旧字段读取和 `compute_polarization()` 接口保留，批量接口为增量能力；
- 若二维追踪未通过验收，可移除新入口/新纯模块和相应测试，不需要修改或迁移既有结果；
- 不通过回滚删除任何 `scripts/.out` 正式数据；测试/smoke 临时目录按任务级精确路径清理；
- 不使用 `git reset --hard` 或覆盖用户已有修改。

## 19. 实施 Checklist

- [x] 用户确认本 Plan 并授权开始实施
- [x] 在 `scripts/parameter.json` 增加并校验 `unit_cell_2d` 专属配置，几何读取锁定为其内部 `b0_nm/eta/zeta`
- [x] 实现纯 Python schema、角度、winding 和绘图合同
- [x] 实现离线 `scripts/analysis/unit_cell_2D.py`
- [x] 实现批量 COMSOL field/polarization extraction
- [x] 实现 `run_unit_cell_2d.py`、二维调度、追踪和恢复
- [x] 验证 runner 成功后默认执行 Analysis，`--skip-analysis` 仅作显式诊断开关
- [x] 完成纯 Python 与聚焦单元测试
- [x] 完成单点 extractor 数值对照
- [x] 完成 3x3 实机 smoke test、离线重绘与实际 PNG 检查
- [ ] 向用户提交正式 19×19 运行前参数核对
- [ ] 经用户确认后运行正式扫描
- [ ] 完成正式结果 QA 和离线重绘验证
- [x] 回填当前阶段执行记录与状态

## 20. 执行记录（实施后回填）

### 实际改动

- 新增 `unit_cell_2d` 参数对象并由 `parameter_config.py` 严格校验；按用户后续要求，几何改为只读取该对象内部的 `b0_nm/eta/zeta`，不再读取顶层 `cavity`。
- 新增 `run_unit_cell_2d.py`，实现 2D Cartesian 调度、Gamma 锚定、多前驱 overlap、Hungarian、ambiguity gate、point cache 和一次 repair pass；求解完成后默认 Analysis。
- 新增纯 Python `unit_cell_2d_analysis.py` 与 `scripts/analysis/unit_cell_2D.py`，生成 raw `c(k)` linear/log Color Map、Stokes 角度图、四种偏振矢量图、winding 双面板图及 shrinking-loop CSV/PNG。
- `hexagon_unit_cell.py` 增加一次 Interp 的批量复数字段读取和 polarization 数据接口，旧接口保留。
- `comsol-data-analysis/` 始终只读参考，未 import、未修改、未并入本项目。

### 测试与实机验证

- 纯 Python 和既有 unit-cell 聚焦回归共 147 项通过；相关模块 `py_compile` 与 `git diff --check` 通过。
- 合成解析场锁定四种角度定义、`+1/-1/0` winding、全局复相位不变性、非均匀 cell edges、magma/cyclic 色图、等比例坐标及 `10^n` log ticks。
- COMSOL Gamma 单点对照：新批量接口的 `cx_raw/cy_raw` 与旧 `compute_polarization()` 逐复数完全相同，比值均为 `1+0j`。
- 3x3 smoke：9/9 点完成，Gamma 加 8 个 matched，频率跨度 0.005592 THz，overlap 约 0.999，无 missing/ambiguous/repair；默认 Analysis 状态 `complete`。
- 使用保存的 CSV 再次运行离线入口，确认不加载/启动 `mph`；已检查实际 log Color Map、矢量图和 winding 双面板 PNG，并据此修正非均匀局部箭头长度与 log ticks。

### 正式输出

尚未启动或生成正式 19x19 结果。smoke 与绘图 QA 仅位于系统临时目录，并在交付前清理。

### 遗留问题

- 历史 3x3 smoke 当时把 `neigs=10` 写入求解器，`Global Evaluation` 返回 20 行、其中 14 行通过既有有效性筛选，验证了当前求解器的成对返回行为。按用户后续要求，二维 workflow 现默认使用 `unit_cell_2d.eigenmode_pair_count=2`，即请求 2 对、预期 4 个本征解；`unit_cell_band` 系列继续使用顶层 `unit_cell_eigenmode_count=10`，不受本次修改影响。
- 正式 19x19 仍需用户再次确认当前 `unit_cell_2d` 几何、raw 主强度合同、mesh/资源和输出路径后才启动。

### 最终状态

功能实现、纯 Python 验证、单点 extractor 对照和 3x3 COMSOL smoke 均已完成；当前停在正式 19x19 启动审批门前。

## 2026-08-22 图 D 四面板标准输出

- 用户已批准将样例图 D 的 2×2 形式固化为共享标准输出，文件名仍为
  `D_BIC_winding.png`，不新增版本目录或兼容图片。
- 上排保持原始 `complex-axis` 偏振椭圆与逆时针最外层有效 square-loop winding：
  左侧为全网格椭圆、外层 loop 和 loop 椭圆，右侧为展开角及其线性参考。
- 下排对每个点的 `cy_norm_re/cy_norm_im` 同时乘以 10，重新派生偏振量、扫描稳定
  偏振量，并在与上排相同的最外层 loop 上重新计算 winding；左侧标题明确标出
  `$c_y \\times 10$`。
- 右侧两幅图只显示整数标题 `Winding number = X`，不在图中显示 residual；
  residual 和 shrinking-loop 多尺度稳定性仍保留在数值日志中用于质量审计。
- 左侧两幅图继续使用真实 k 坐标、等比例坐标轴，并在数据边界外延一个实际采样间隔。
- 本次批量重绘只覆盖现有单组结果的 `10_overview/D_BIC_winding.png`，不运行
  COMSOL，不覆盖 A–C，也不处理不含 `polarization_grid.csv` 的 sweep 汇总目录。

### 执行结果

- 共享 `plot_winding_detail()` 已改为正式 2×2 生成器；`run_analysis()` 即使
  complex-axis shrinking-loop 没有多尺度稳定结果，也会使用最外层有效 loop 生成
  图 D，同时继续在 winding summary 中如实记录 `no_stable_loop`。
- 已原位重绘 ζ=1.151、1.154、1.155、1.156、1.157、1.158、1.161 共 7 组结果。
  上排 winding 依次为 −1、−1、−1、−1、−1、−1、0；下排
  `c_y × 10` 后依次为 0、0、−2、−2、−2、−2、−1。
- 实际 PNG 均为 3120×2640；已抽查 ζ=1.156 和 ζ=1.161，标题、四面板布局、
  等比例坐标、一个采样间隔的外边界和裁切均正常。
- 聚焦测试结果为 `18 passed`（analysis）和 `23 passed`（runner）；相关
  `py_compile` 与 `git diff --check` 通过。
- 本次没有启动 COMSOL，没有重写 A–C，没有新建结果目录。
# 2026-08-21 输出合同修订与执行记录

- 正式输出根改为 `scripts/.out/unit_cell_2D/`，不再位于 `unit_cell_band` 下。
- `00_model` 保存几何图和后续运行保存的 MPH；`80_logs` 保存汇总 CSV；
  `10_overview` 保留原 9 张正式图；2026-08-25 起额外保留一张明确标为
  smoothing-only 的 B4 方法参照图，不把它误作第二张正式物理解。
- 正式图合同为 A1/A2 总强度 linear/log，B1/B2 为归一化
  `|cx(k)|` linear/log，B3/B4 为归一化 `|cy(k)|` linear/log，
  C1/C2 为 complex-axis 矢量图和 stokes-axis 相位图，D 为 complex-axis
  BIC winding 双面板图。
- 当前 361 点结果已迁移到新一级目录并离线重绘；没有重跑 COMSOL。
- D 图使用逆时针闭合 square loop。当前数据在 half-width=0.1 上得到
  total angle=-4pi、winding=-2、residual=0，第二面板与数值定义自洽。
- 聚焦测试：`33 passed`。最终完整检查与 MPH 保存能力仍按本文件后续执行记录更新。
# 2026-08-21 q 网格参数修订

- 删除输入参数 `unit_cell_2d.q_axis_over_G`，不再接受人工坐标列表。
- 新输入为 `q_max_over_G` 和 `q_points_per_axis`；当前值分别为
  `0.1` 与 `19`。
- Kx、Ky 强制复用
  `linspace(-q_max_over_G, q_max_over_G, q_points_per_axis)`，因此两方向
  边界、点数和间隔完全一致。
- `q_points_per_axis` 必须是大于等于 3 的奇数，以保证 Gamma 点唯一存在。
- 新等间距网格使用 `_uniform19` 系列后缀，不复用旧
  `_multiscale19` 非均匀网格缓存；旧正式结果保持原位。
- 本次只修改参数、派生逻辑、缓存身份、命名和测试，不启动 COMSOL。

# 2026-08-25 uniform25 B4 无拓扑先验重绘

> **2026-08-25 最新回滚：** 用户否决本节此前采用的 coherency-matrix、Gaussian
> 与 smoothing-only 路径。当前正式 B4 已恢复为与 B1--B3 完全相同的共享绘图合同：
> 使用 `cx_raw/cy_raw`，共同除以全网格 `max|c(k)|`，经共享 `_smooth_grid()` 三次
> 显示插值后用 `magma + LogNorm` 输出。正式输出仅为
> `B4_cy_magnitude_log.png`；不再生成 smoothing-only 图，也不做 nodal-line 优化。
> `B4_cy_magnitude_log.pre_regenerate_20260824.png` 作为受保护参考，不得覆盖。
> 下方 coherency/Gaussian 内容仅保留为已否决的历史执行记录，不代表当前实现。

- 已按共享 B 图路径原位生成 ζ=1.154、1.155、1.156、1.157 四张
  `B4_cy_magnitude_log.png`；没有改写 A、B1--B3、C、D 或数值数据。
- ζ=1.156 的参考图 SHA256 在重绘前后均为
  `7D17293ACA90AEF6C21AA4160B5EE6AD0D6EAEBC48B6EAD30BE7E89BE6F970F1`。
- 四组均未生成 `B4_cy_magnitude_log_smoothing_only.png`，实际 PNG 视觉检查通过；
  `tests/test_unit_cell_2d_analysis.py` 为 `22 passed`，相关模块 `py_compile` 通过。

- 用户否决 topology-constrained reconstruction。branch extraction、Gamma-constrained
  fit、distance suppression、line/Gamma floor 和中心恢复均把 nodal-line 先验写入
  显示图，不能构成真实连续场或精确 nodal lines 的独立证据。
- 数据审计进一步发现 raw COMSOL 本征矢量在 k 网格上带有任意逐点尺度；四组
  raw/planar-L2-normalized 比例的变异系数为 0.564–0.569，旧 raw B4 不能用于比较
  参数演化。正式 B4 改用保存的 `cx_norm/cy_norm`，并继续使用全网格共同的
  `max(sqrt(|c_x|^2+|c_y|^2))` 归一化。
- 正式连续重建对象为 gauge-invariant coherency matrix `J = c c^dagger`。对
  `Jxx/Jyy/Jxy_re/Jxy_im` 使用相同的分片双线性插值和各向同性 Gaussian，输出
  `sqrt(Jyy)`；双线性和 Gaussian 都是正权平均，因而保持半正定和 trace 上界，且
  不识别暗谷数量、位置或方向。
- 四组共同比较 `sigma=0.25/0.50/0.75 Delta q` 后固定为 `0.50 Delta q`：其 Nyquist
  响应约 0.291，低值区域单格 Laplacian 起伏相对 0.25 再下降 36%–41%，而 0.75 的
  全场改动已偏强。该尺度由采样间隔决定，不按 zeta 或曲线形状调参。
- 每组额外生成 `B4_cy_magnitude_log_smoothing_only.png`，直接对归一化 `|c_y|`
  做同尺度双线性和 Gaussian，作为与正式“强度先平均、再开平方”路径的参照。
- 已预检 ζ=1.154、1.155、1.156、1.157 四组 `kmax0.15_uniform25`：每组均为
  625 行、25×25、band=`p2`，状态只含 `gamma/matched`，无重复、缺点或无效点，最低
  selected score 为 0.920857。
- 数学推导、失败候选、尺度选择、数据审计、精确解边界、代码入口和绘图合同已整理到
  `docs/unit_cell_2d_b4_nodal_line_reconstruction.md`。严格的 exact nodal line 仍需
  新增自适应 COMSOL eigensolve/root refinement；577×577 显示像素不宣称为新解。
- 已原位覆盖四组 `10_overview/B4_cy_magnitude_log.png`，并各新增一张
  `B4_cy_magnitude_log_smoothing_only.png`。没有覆盖 A、B1–B3、C、D、CSV、JSON、
  MPH 或参数文件，也没有启动 COMSOL。
- 8 张实际 PNG 已逐张检查：标题无裁切，0.05 坐标刻度、equal aspect、magma、
  colorbar 等高、真实最大值和刻度防重叠均正常；正式图与仅平滑参照图显示一致的参数
  演化，正式图对孤立单点低坑抑制更强。
- analysis 聚焦测试为 `27 passed`；analysis + runner + zeta-scan 联合测试为
  `57 passed, 1 failed`。唯一失败是既有 runner 测试硬编码期望 zeta=1.156，而用户
  当前 `parameter.json` 为 zeta=1.154，与本次 B4 代码无关，未擅自修改参数或该测试。
  相关模块 `py_compile` 与 `git diff --check` 通过；方法候选临时目录已清理。

# 2026-08-25 B4 保值约束最小曲率连续场

## 方案概述

- 保留既有 `B4_cy_magnitude_log.png`，新增独立输出
  `B4_cy_magnitude_log_minimum_curvature.png`。
- 数据仍为标准 B 图使用的 `cx_raw/cy_raw`，并以全网格共同的
  `max sqrt(|cx_raw|^2+|cy_raw|^2)` 归一化。
- 在包含全部原始采样节点的稠密 Cartesian 网格上，求解离散薄板能量
  `sum(Dxx^2 + 2 Dxy^2 + Dyy^2)` 的凸最小化；625 个原始节点作为硬等式约束，
  同时施加 `0 <= M <= 1` 的物理幅值约束。
- 不输入 valley/nodal-line 的条数、方向、交点、中心位置、连通性或曲线模型；
  是否出现支线、相交、断裂或闭合结构完全由样本与上述通用约束决定。
- 稠密节点之间采用分片双线性延拓，因此连续且保持区间约束。log 图只在显示层
  使用统一的正值下限，不改变求解场和任何原始采样节点。

## 绘图合同

- 静态 Matplotlib heatmap，300 DPI；每组 uniform25 生成一张新图。
- 标题保持 `Normalized |c_y(k)| (py)`；`magma`、equal aspect、0.05 坐标刻度、
  colorbar 等高、最大值显式标注和重叠刻度删除规则均复用标准 B4。
- 不覆盖标准 B4、参考图、A/B1--B3/C/D、CSV、JSON 或 MPH。

## 验证与回滚

- 数值测试锁定原始节点误差、`[0,1]` 约束、目标能量下降、文件名和标准 B4 不变。
- 对四组 ζ=1.154、1.155、1.156、1.157 的实际 PNG 做视觉 QA，并记录求解收敛。
- 回滚只需删除新增入口与新文件；标准 B4 始终保留，不需要恢复备份。

## 执行结果

- 已在 `comsol_workflow/unit_cell_2d_analysis.py` 实现解析离散薄板能量/梯度、
  sample-exact box-constrained 优化、独立绘图函数和批量重绘入口；分析 CLI 新增
  `--minimum-curvature-b4-only`，不会触发标准图重绘或 COMSOL。
- ζ=1.154、1.155、1.156、1.157 四组均生成
  `B4_cy_magnitude_log_minimum_curvature.png`，并在 `99_config` 写入独立诊断 JSON。
- 四组 625 个样本的最大误差均为 0，迭代数为 186--193，最终弯曲能量相对双线性
  初值下降；重建场均位于 `[0,1]`，全场最小值严格为正，没有人为产生 exact zero line。
- ζ=1.156 完成 4/8/12 细分收敛检查；8 与 12 在共同节点上的 RMS/最大差为
  `9.15e-6/5.83e-5`，故固定 8 细分的 193×193 优化网格。
- 四张实际 PNG 已检查，标题、0.05 坐标刻度、equal aspect、magma、colorbar 等高、
  最大值刻度及画布布局正常。标准 B4 与 ζ=1.156 受保护参考图哈希未改变。
- 聚焦测试为 `26 passed`；相关模块 `py_compile` 和 `git diff --check` 在最终交付前
  再次执行。

# 2026-08-25 B4 low-value 局部放宽与 valley-floor 平滑

## 方案概述

- 标准 `B4_cy_magnitude_log.png`、历史参考图以及 A/B1--B3/C/D 均保持不变；只覆盖
  已独立命名的 `B4_cy_magnitude_log_minimum_curvature.png`。
- 在原始 25×25 `log10(|c_y|/Cmax)` 样本上，以一维 Otsu 最大类间方差自动分出
  low-value 集合，不手工输入支线条数、方向、交点或位置。
- 只在 low-value 点之间按原始网格两步邻域建立无向图；逐点暂时移除中心样本，使用
  其余 low-value 邻居在 `log10` 幅值上做空间加权局部线性拟合。邻域退化时回退为
  加权局部平均，预测值限制在邻居范围内以避免外推。
- high-value 样本继续作为最小曲率场的硬约束；low-value 样本改用图平滑后的局部值，
  并记录原值、修正值、点数、阈值、局部拟合/回退次数、RMS/最大修正和连通分量数。
- 修正后的完整 25×25 网格继续进入 `[0,1]` box-constrained minimum-curvature 求解；
  新增的先验只有“valley floor 在局部邻域内平滑”，不加入曲线拓扑。

## 验证与回滚

- 测试锁定 Otsu 数据分割、图平滑常量保持、非 low-value 点不变、修改诊断和输出合同。
- 先以 ζ=1.156 评估实际修改幅度和视觉效果，再以完全相同的自动规则覆盖四组新图。
- 回滚到上一版只需禁用 low-value 预处理；标准 B4 无需恢复。

## 执行结果

- 已实现一维 Otsu low-value 自动分割、两步邻域 leave-one-out 加权局部线性预测、
  `median + 1.4826 MAD` 稳健残差估计，以及只对超过 `3×MAD` 的 low-value 异常点
  进行有限修正。预测值始终限制在邻域值域内；未输入暗线条数、方向、交点或连通性。
- Gamma 原始采样点参与邻域拟合但被硬保护，四组审计表中的 Gamma 修正量均为 0。
  ζ=1.156 的 Gamma 值在处理前后均为 `1.7540423187208764e-05`，也是重建场全局
  最小值；实际 PNG 中中心奇点位于两条低谷交汇处清晰可见，且没有额外中心 marker、
  黑点恢复或绘图后处理。
- ζ=1.154、1.155、1.156、1.157 的 low-value 点数分别为 181、203、187、205；
  实际修改点数分别为 17、16、13、13。修正后的 25×25 节点进入 193×193、
  `[0,1]` box-constrained minimum-curvature 求解，四组均收敛并降低弯曲能量。
- 四组 `B4_cy_magnitude_log_minimum_curvature.png` 已原位更新；每组同时写入
  `80_logs/B4_minimum_curvature_low_value_relaxation.csv` 和
  `99_config/B4_minimum_curvature_diagnostics.json`。标准 B4、ζ=1.156 历史参考图、
  A/B1--B3/C/D、CSV 原始数据和 MPH 均未改动，未启动 COMSOL。
- 聚焦测试在 Gamma 保护完成后为 `27 passed`；最终交付前再次执行 `py_compile`、
  聚焦 pytest 与 `git diff --check`。

# 2026-08-25 ζ=1.156 B5 二次型 valley 拟合

## 方案与绘图合同

- 只处理现有 ζ=1.156、`kmax0.15_uniform25` 结果，不启动 COMSOL，不覆盖标准 B4、
  minimum-curvature B4 或其他 A--D 图片。
- 拟合窗口固定为 `|kx/G| <= 0.05` 且 `|ky/G| <= 0.05`；只使用窗口内原始
  `cx_raw/cy_raw` 经全网格共同 `max|c(k)|` 归一化后的 `|cy|` 采样值。
- 在该窗口内部独立执行一维 log-amplitude Otsu，自动选取 low-value 点。拟合模型严格为
  `A(kx/G)^2 + B(ky/G)^2 = C`；本次按用户给定物理前提固定 `C=0`，并以
  `A^2+B^2=1` 消除整体尺度简并、`A>=0` 固定符号。
- 使用低值点坐标的正交最小二乘/SVD 求解最小代数残差方向；不修改原始幅值，也不使用
  B4 图像像素反推曲线。若 `A*B>=0`，模型不存在两条实分支，必须报告而不是强画。
- 新图命名为 `10_overview/B5_cy_magnitude_log_quadratic_valley_fit.png`。背景复用
  B4 minimum-curvature 的完整 field、`magma + LogNorm`、标题、坐标、equal aspect、
  colorbar 等高与最大值规则；只在拟合窗口内叠加两条二次型零等值分支，不外推。
- 图例给出归一化后的拟合方程和窗口范围；拟合系数、选点数、Otsu 阈值、残差与分支
  斜率写入 `99_config/B5_quadratic_valley_fit.json`，入选点写入
  `80_logs/B5_quadratic_valley_fit_points.csv`。

## 验证与回滚

- 测试锁定 `C=0`、`A^2+B^2=1`、符号约定、选点窗口、两条实分支、曲线不越出
  ±0.05、输出文件名、magma/log 色标和 B4 文件不变。
- 导出后检查 ζ=1.156 实际 PNG，确认中心奇点仍可见、两条拟合线贴合 low-value
  valley、图例不遮挡主体且 colorbar/画布无裁切。
- 回滚只需删除 B5 入口、B5 PNG 与对应独立 JSON/CSV；原始 B4 无需恢复。

## 执行结果

- 已在共享分析模块实现窗口内 Otsu 选点、固定 `C=0` 的二次齐次 TLS/SVD 拟合、
  `A^2+B^2=1`/`A>=0` 唯一化、两条实分支检查、B5 绘图和 JSON/CSV 审计输出；
  离线 CLI 新增 `--b5-quadratic-valley-fit-only`。
- ζ=1.156 的 ±0.05 窗口包含 81 个原始点，自动选中 23 个 low-value 点。拟合结果为
  `A=0.9965552185`、`B=-0.0829318783`、`C=0`，等价于
  `ky=±3.4664898502 kx`；RMS/最大代数残差为 `9.22305e-05/1.55712e-04`。
- 已生成 `B5_cy_magnitude_log_quadratic_valley_fit.png`、
  `B5_quadratic_valley_fit.json` 和 `B5_quadratic_valley_fit_points.csv`。拟合线仅覆盖
  ±0.05 窗口；为避免遮住数据，在 Γ 周围只对 overlay 留出四分之一个原始采样间隔
  (`0.003125G`) 的透明缺口，数学曲线仍由 `C=0` 在 Γ 相交；最终拟合线按用户要求
  使用灰色虚线。
- 实际 PNG 已检查：中心奇点可见，拟合线贴合局部 low-value valley，图例、标题、
  equal aspect、colorbar 等高、最大值刻度和画布裁切均正常。
- 标准 B4、minimum-curvature B4 与历史参考图 SHA256 未改变；没有重绘其他图、修改
  原始数值数据或启动 COMSOL。
- 新增 B5 聚焦测试为 `2 passed`；完整 analysis 测试、`py_compile` 和
  `git diff --check` 在最终交付前再次执行。

# 2026-08-25 uniform25 B4 拓扑自适应低谷重建

## 方案概述

- 原始 `10_overview/B4_cy_magnitude_log.png` 是标准 B 图合同的基准输出，先记录
  SHA256，整个流程不得覆盖。新图固定命名为
  `B4_cy_magnitude_log_valley_reconstruction.png`。
- 数据源继续使用 `cx_raw/cy_raw`，`|c_x|` 与 `|c_y|` 共同除以完整网格
  `max sqrt(|c_x|^2+|c_y|^2)`；不启动 COMSOL，不改原始 CSV。
- 在 25×25 原始 log-amplitude 网格上自动检测逐行、逐列局部低谷，以 Otsu 低值类
  限制候选，并用局部抛物线估计亚网格低谷位置。每个截面最多保留两个数据支持的候选，
  不为 Γ 人工补点或补第二支线。
- 允许的显式先验是 low-value valley 在 Γ 邻域具有 X/类 X 或中心双曲线拓扑。用
  归一化的不定二次型
  `q^T H q = C, det(H) < 0` 做稳健拟合；`C` 自由拟合，因此支线不被强制经过 Γ。
  该模型统一容纳 X-like、左右双谷和上下双谷。候选按实际双曲线的实分支配对，不按
  ζ 或预期方向预先分组。
- 用候选重采样得到拓扑置信诊断：`C=0` 落入区间或分裂低于半网格分辨率时，数据选择
  嵌套 `C=0` 临界模型；否则保留自由拟合的非零 `C`，并按开口主轴标记左右或上下双谷。
  Γ 候选不参与自由 `C` 拟合，因此临界连接不是由中心单点强迫出来的。
- 二次型只负责中心拓扑、特征轴与初始配对。X-like 的两条远场支线分别在特征坐标中做
  经过 Γ 的稳健奇次五阶拟合，以恢复参考图中由数据支持的弯曲；非临界双曲线保留二次型
  实分支。
- 在分片线性连续背景上使用两支线距离乘积重建均匀 valley，宽度为 `0.65Δq`，普通线底
  为分量最大值的 `1e-3`。仅对数据选出的 X-like 状态，线底在 `0.5Δq` 内平滑过渡到
  真实 Γ 值；最后再逐位恢复 Γ 样本与原始最大值，不画中心 marker 或圆形暗斑。

## 绘图合同

- 单面板静态 Matplotlib PNG，300 DPI；标题保持
  `Normalized |c_y(k)| (py)`，使用 `magma + LogNorm`。
- 坐标范围来自原始 uniform25 网格，主刻度间隔 0.05，`kx-ky` 等比例；colorbar 与
  图框等高，显式显示实际最大值，并删除与最大值标签重叠的最近 decade。
- 不叠加拟合线、中心 marker 或拓扑文字；拟合证据只写入 CSV/JSON，避免把诊断线混入
  标量场图。

## 输出、验证与回滚

- 每组在 `80_logs` 写候选/支线配对及重建采样审计 CSV，在 `99_config` 写拟合矩阵、
  特征轴、`C`、bootstrap 区间、拓扑分类、残差、支线支持数、Γ/最大值保真与输入输出
  哈希。
- 合成测试必须覆盖 X-like、左右双谷、上下双谷，另锁定 Γ 严格保值、`det(H)<0`、
  原始 B4 不变、输出命名、magma/log 色标、0.05 刻度和 colorbar 最大值。
- 先验证 ζ=1.156，再以同一算法处理 ζ=1.154、1.155、1.156、1.157；逐张检查实际
  PNG，并核对四张标准 B4 的前后 SHA256。
- 回滚只需移除新 CLI 入口、代码路径及独立新 PNG/CSV/JSON；标准 B4 无需恢复备份。

## 执行结果

- 已在 `comsol_workflow/unit_cell_2d_analysis.py` 实现 Otsu 低谷筛选、逐行/逐列亚网格
  极小点、候选合并、不定二次型多初值稳健拟合、128 次 bootstrap、三类拓扑选择、支线
  配对、X-like 奇次五阶远场拟合、距离乘积重建和 Γ 逐位保值；离线 CLI 新增
  `--valley-reconstruction-b4-only`，没有启动 COMSOL。
- 四组实际分类为：1.154/1.155 `upper-lower`，1.156 `x-like`，1.157
  `left-right`。1.156 的自由 `C=-0.0007692373`，95% 区间
  `[-0.00295683, 0.000406954]` 包含 0，故选择临界 `C=0`；其他三组区间均排除 0。
- ζ=1.156 的两个中间错误版本已否决：非零 `C` 造成中心场分裂和孤立暗点；全局二次
  `C=0` 又把远处弯曲低谷拉成直线。最终版用二次型判别中心拓扑，用配对奇次五阶曲线
  恢复远场弧度，并用半步线底过渡凸显真实中心奇点。
- 四组均已生成 `B4_cy_magnitude_log_valley_reconstruction.png`、候选 CSV、625 点重建
  审计 CSV 和诊断 JSON。四组 `gamma_bitwise_equal=true`、`det(H)<0`，标准 B4 的
  SHA256 前后完全一致；原始 B4、受保护参考图、A/B1--B3/C/D、CSV 和 MPH 未改动。
- 四张实际 PNG 已检查：远场低谷连续且无离散深坑，1.156 的弧度和中心奇点接近受保护
  参考图效果；0.05 坐标刻度、equal aspect、magma、最终场决定的 log 下限、统一
  `10^-1` colorbar 上限、
  等高 colorbar 和画布裁切均正常。
- 新增拓扑/Γ/绘图聚焦测试为 `5 passed`；最终完整 analysis 测试为 `34 passed`，相关
  `py_compile` 与 `git diff --check` 均通过。四张最终 PNG 已逐一视觉检查，colorbar
  顶端均明确显示 `10^-1`。

# 2026-08-25 valley-reconstruction B4 中心 3×3 保值与统一色标

## 方案与绘图合同

- 只原位更新四组 uniform25 的
  `10_overview/B4_cy_magnitude_log_valley_reconstruction.png` 及其对应审计 CSV/JSON；
  标准 `B4_cy_magnitude_log.png`、历史参考图和其他 A--D 图保持不变，不启动 COMSOL。
- 中心保护区定义为原始 25×25 采样网格中以 Γ 为中心、沿两个索引方向各延伸一格的
  3×3 方块，共 9 个不可变源样本；处理过程不得修改或回写这些 `float64` 输入值。
- 实测 Γ 与最近邻相差约两阶；若把 3×3 九点同时作为连续图像的硬插值约束，会形成
  宽暗团并使两条 valley 起伏。故连续场与不可变源数据分层：Γ 是唯一硬显示约束，
  其余八点保留在审计表中作为原始证据，但不强迫正则化场逐点穿过。
- 四张图统一使用 `magma + LogNorm`，固定 `vmin=1e-5`、`vmax=1e-1`；colorbar 显式
  显示 `$10^{-5}$` 与 `$10^{-1}$`，中间刻度为完整十进制幂。

## 验证与回滚

- 测试锁定 9 个源点逐位不被修改、Γ 重建值逐位相等、源点审计标记、固定 clim、端点 ticks/labels、
  原始 B4 哈希不变及既有拓扑分类。
- 重绘后逐张检查中心区域连续性、ζ=1.156 奇点、colorbar 等高与端点无裁切；运行完整
  analysis pytest、相关 `py_compile` 和 `git diff --check`。
- 回滚仅需恢复本节涉及的共享重建/绘图代码并重新离线绘制；无需恢复或重算 COMSOL。

## 执行结果（非中心计算最小值修订）

- 否决“每条支线分别拟合 log-amplitude floor”的中间版本。该版本在 ζ=1.156 中把
  远端 floor 抬高到约 `5.4e-3`，使参考图中连续可见的低谷局部消失；离散网格到
  精确谷底的横向距离残差不应被解释为沿 valley 的本征幅值变化。
- 当前实现继续使用数据选择的拓扑、支线配对和远场曲率，但两条支线共享中心邻域之外
  所有有效低谷候选的全局计算最小值。ζ=1.156 的锚点为
  `(kx/G, ky/G)=(-0.02480924, -0.07496036)`，归一化幅值为 `3.5628723e-4`。
  Γ 仍逐位保持为 `1.7540423187208764e-5`，并只在半个原始采样间隔内过渡。
- 四组 line floor 分别为：ζ=1.154 `1.2651644e-3`、ζ=1.155 `2.0764816e-4`、
  ζ=1.156 `3.5628723e-4`、ζ=1.157 `4.3299046e-4`。对应拓扑仍为
  `upper-lower / upper-lower / x-like / left-right`。
- 四张 `B4_cy_magnitude_log_valley_reconstruction.png` 已离线重绘并逐张检查：低谷连续、
  远场弧度保留，ζ=1.156 的 Γ 奇点清晰；所有图继续使用固定
  `LogNorm(1e-5, 1e-1)`。四组 `gamma_bitwise_equal=true`、中心 3×3 源数据不可变，
  标准 B4 前后 SHA256 一致；未启动 COMSOL，也未覆盖标准 B4 或历史参考图。
- 最终验证为 `34 passed`；相关 `py_compile` 与 `git diff --check` 均通过。
