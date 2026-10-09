# Finite 基模包络筛选：计划与执行记录

日期：2026-08-02
状态：已撤回（由最近目标频率单模选择取代）

> 2026-08-02 用户明确撤回本文件描述的空间包络判据。当前正式行为见
> `plan-execute_20260802_finite_nearest_frequency_single_mode.md`；中心峰值、峰值壳层和
> 壳层向外回升不再参与 mode 选择。本文件仅保留为决策历史。

## 目标

为 full finite 与 finite-quarter 增加共用、可关闭的基模包络筛选。开关开启时，
先在所有候选有效 mode 中识别有限腔 0 阶包络，只允许选中的一个 mode 进入
far-field、完整 cavity+cladding cell map、finite-lattice Fourier、胞内对称性分解和
mode scoring。其余 mode 仅保留筛选所需的原始导出与审计信息，不执行昂贵分析。

## 判据

- 使用 cavity 内每个完整 unit cell 的粗网格 `mean(|Hz|^2)`，不使用单点场值。
- 由于先对完整 cell 积分，p/d 胞内的符号节点不会被误认为有限腔包络节点。
- 对 cell intensity 计算中心 cell/全局峰值、壳层平均强度峰值位置、壳层向外回升量。
- 基模候选要求中心 cell/全局峰值不低于 `0.5`、壳层平均峰值位于中心区
  `max(前 2 壳层, cavity 半径的前 15%)`，且归一化壳层均值的累计向外回升量不超过
  `0.15`；在通过者中按
  中心占优且向外回升更少的综合分数选出唯一 mode。第三项使“中心仍最强但外侧再次
  出现明显亮环”的有节点包络也不能通过。
- 若没有 mode 通过绝对判据，保留审计 CSV 并中止后续分析，不静默选择高阶包络。

## 配置与兼容

- 在 `run_finite.py` 提供布尔开关，finite-quarter 复用同一开关和实现；本次按用户要求
  默认开启。
- 开关关闭时保持现有行为：所有有效 mode 进入已启用的后处理阶段，标准 intensity
  map 仍覆盖所有有效 mode。
- 显式 Fourier mode 列表若存在，先限定筛选候选集，再在其中选择唯一基模。
- 筛选结果写入 case 的 `10_overview/fundamental_mode_screen.csv`，并在
  `eigenfrequencies.csv` 中记录 `analysis_selected`。

## 数据流

1. COMSOL eigensolve 与现有有效性检查保持不变，并导出候选 mode 的重建完整 `Hz` 场。
2. 默认使用与正式 finite Fourier 相同的 41 点边长 unit-cell rho 网格对 cavity cells
   采样，得到候选包络指标。
3. 选择唯一 0 阶包络 mode，并把其 mode index 显式传给 far-field 与 lattice Fourier。
4. lattice Fourier 不再为未选 mode 补做 map-only 处理。
5. finalizer 把筛选 CSV 移至正式 `10_overview`，保留候选 mode 的原始字段用于审计。

## 验证

- 纯数值单元测试覆盖基模/有节点包络区分、唯一选择、无通过候选时报错与关闭开关兼容。
- far-field 和 lattice Fourier 测试锁定显式 mode 子集只处理选中 mode。
- 使用已完成的 20+20 finite-quarter mode 0/1 字段离线筛选，确认选择视觉上无包络节点的
  mode 1，不启动 COMSOL、不覆盖正式结果。
- 运行相关 pytest、`py_compile` 与 `git diff --check`。

## 回滚

关闭配置开关即可恢复所有有效 mode 的分析。代码级回滚则移除筛选函数、far-field 的
可选 mode 参数、full/quarter 接线、审计输出映射与相关测试；现有求解及原始字段格式不变。

## 执行记录

- 在共享 `lattice_fourier_postprocess` 中实现完整 cell `mean(|Hz|^2)` 包络指标、
  绝对门槛、唯一选择、CSV 审计及 eigenfrequency 选择标记。
- Full finite 与 finite-quarter 均在导出重建 `Hz` 后、far-field/Fourier 前执行筛选；
  开关默认开启，关闭后恢复所有有效 mode 的原行为。far-field 与 finite Fourier 都只
  接收选中的 mode index，未选 mode 不再生成标准 map 或执行分解。
- finalizer 将 `fundamental_mode_screen.csv` 放入 `10_overview`；候选 mode 原始场仍保留，
  用于筛选复核和避免不可审计的静默丢弃。
- 对保存的 20+20 finite-quarter `shiftx0.050_shifty0.080` 做了两次离线校验，未启动
  COMSOL。41 网格下 mode 0 的中心/峰值为 `0.4135`、外向回升 `0.5261`、峰值壳层为
  最外层 20，判为高阶包络；mode 1 对应 `0.9710`、`0.0304`、壳层 1，唯一入选。
  11 网格复核仍选择 mode 1；正式默认保留 41 以提高胞内积分稳定性。
- 聚焦回归为 `68 passed`；四个修改模块均通过 `py_compile`，`git diff --check`
  通过。完整测试集为 `369 passed, 1 failed`；唯一失败是当前用户参数
  `cladding_y_over_x_shift_ratio=1.0` 与既有测试固定期待 `1.6` 不一致，和本次筛选无关，
  因此未覆盖当前参数文件。

### 2026-08-02 小 cavity 修正

10 层 cavity 的 `shiftx0.050_shifty0.080` 实际基模 mode 1 的壳层均值峰值位于第 2
壳层，中心/全局峰值为 `0.9762`、向外回升仅 `0.0084`。旧的纯比例门槛
`peak_shell / radius <= 0.15` 只能允许到第 1 壳层，因离散壳层取整误拒绝该模式。
门槛修正为至少允许中心前 2 壳层；mode 0 的峰值仍在最外第 10 壳层且向外回升
`0.7492`，继续被拒绝。新增半径 10 的第 2 壳层通过、第 3 壳层拒绝回归测试。
