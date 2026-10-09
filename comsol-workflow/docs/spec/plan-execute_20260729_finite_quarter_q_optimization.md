# plan-execute：20–20 finite-quarter py-mode Q 优化

日期：2026-07-29

## 实施前方案概述

### 目标

在固定 cavity `245-0.960-1.155`、20 层 cavity、20 层 cladding 和 mesh9
条件下，从 cladding `243.7-0.980-0.928` 出发，仅微调 cladding 的 `zeta`
与 `b0`，最大化 finite-quarter 中用户定义 py-mode 的 Q。目标首先达到
`Q > 6700`，并尽量达到 `Q > 7000`。

### 已确认约束

- 参数优先级为 `zeta > b0`；本次不改变 `eta=0.980`；
- `SYMMETRY_ID=1` 固定不变；项目内部该扇区对应用户定义的 py-mode；
- COMSOL 在该设置下实际返回两个解。低 Q 解通常为 `mode0`，高 Q 解通常为
  `mode1`；程序不得硬编码 mode 序号，而应从两个有效解中选择较高 Q 解；
- cladding shift 固定为 `0.000`，profile 固定为 `uniform`；
- eigensolver 中心频率固定为当前 cavity `p2@Gamma` 的共享中心频率；
- 除最终的“迭代次数–py-mode Q”点图外，不运行或生成几何、场、far-field、
  finite-lattice Fourier 等图像与后处理；
- 每次候选完成后立即更新 CSV，以支持中断续算；
- 最终保留迭代 CSV、Q 点图、最优参数记录和最优候选 MPH；其余候选 MPH
  作为过程文件清理。

### 搜索策略

1. 复用已有起点结果：`b0=243.7 nm, eta=0.980, zeta=0.928`，
   `Q=5070.723304628295`，不重复求解。
2. `zeta` 方向搜索：先计算 `0.924` 与 `0.932`，从改善方向以 `0.004`
   前进，首次下降后在局部最优点周围以 `0.001` 细化；硬边界
   `0.912 <= zeta <= 0.944`。
3. 固定最优 `zeta`，对 `b0` 先计算 `243.2 nm` 与 `244.2 nm`，从改善方向
   以 `0.5 nm` 前进，首次下降后以 `0.1 nm` 细化；硬边界
   `241.7 nm <= b0 <= 245.7 nm`。
4. 在联合最优 `b0` 下复核 `zeta_best +/- 0.001`。
5. 达到 `Q > 6700` 后仍继续完成局部细化，以争取 `Q > 7000`；默认最大
   20 个新候选，防止搜索失控。

### 数据流与输出

- 新增 `scripts/run_sweep/optimize_finite_quarter_q.py`；协调进程只管理候选、
  恢复、CSV 与最终点图。
- 每个候选通过独立子进程加载专用临时参数 JSON，避免改写
  `scripts/parameter.json`，并复用现有 finite/quarter 几何、手动 mesh9、物理场和
  solver 构建函数。
- worker 只保存并读取 MPH、本征频率、Q 和有效性，不调用正式 main 中的绘图及
  Fourier/far-field 后处理。
- 输出归入唯一允许的一级结果目录 `scripts/.out/finite_cavity/`，使用独立优化目录；
  临时参数和候选 MPH 位于该任务自己的 `.staging/` 中。
- CSV 至少包含：迭代编号、阶段、`b0/eta/zeta`、两个 mode 的频率/Q/有效性、
  选中 mode、py-mode 频率/Q、是否刷新历史最优和累计最优 Q。

### 兼容性与安全

- 为独立 worker 增加可选的参数文件环境变量；未设置时，四个主程序仍严格读取
  原有 `scripts/parameter.json`，默认行为不变。
- worker 对 cavity、层数、mesh、eta、shift/profile 和 `SYMMETRY_ID` 做前置检查；
  不满足批准配置时拒绝启动昂贵计算。
- 若返回解数量不是两个，或没有有效解，候选标记失败并停止自动搜索，避免误选模式。
- 已完成候选由参数键恢复，不重复运行；最优模型复制成功后才清理对应候选 MPH。

### 验证计划

- 单元测试：参数覆盖文件、候选去重/恢复、高 Q 双解选择、方向搜索、细化、CSV
  与最终点图数据；
- `py_compile` 新脚本和受影响配置模块；
- 聚焦 pytest；
- `git diff --check`；
- 实际运行前打印并人工核对起点、边界、固定参数、输出目录和首批候选；
- 实际运行后逐行核对 CSV、最优参数、最优 MPH 与点图。

### 回滚方案

- 新脚本可独立删除，不影响四个主入口；
- 参数文件覆盖机制删除后默认主流程仍使用 `scripts/parameter.json`；
- 优化结果全部位于独立任务目录，不改动已有正式 finite-quarter 结果。

## 实施记录与结果

### 已完成实现

- 新增 `scripts/run_sweep/optimize_finite_quarter_q.py`，实现断点续算、候选去重、
  `zeta -> b0` 分阶段搜索、固定 `eta=0.980`、双解有效性检查与高 Q 模式选择；
- 为 worker 增加任务局部参数文件覆盖能力，未设置覆盖时仍读取正式的
  `scripts/parameter.json`，四个主入口的默认行为不变；
- 候选计算只执行几何、mesh9、本征求解和有效性判定，不调用场图、几何图、
  far-field FFT 或 finite-lattice Fourier 后处理；
- 每个候选完成后更新迭代 CSV 与 Q 点图，仅提升并保留当前最优 MPH。

### 验证结果

- 新增优化器聚焦测试，覆盖参数覆盖、严格双解检查、非硬编码 mode 选择、候选生成、
  恢复、CSV 与点图；
- 聚焦回归测试结果：`101 passed`；
- `py_compile`、`git diff --check` 和正式运行前 dry-run 均通过；
- dry-run 已核对起点 Q=`5070.723304628295`，首批候选为
  `zeta=0.924` 与 `zeta=0.932`。

### 正式运行状态

- 已于 2026-07-29 启动正式优化；
- 当前正在计算第 1 个新候选：`b0=243.700 nm, eta=0.980, zeta=0.924`；
- 最优参数、最终 Q、最终 CSV/点图和最优 MPH 待搜索结束后补充。
