# Finite 最近目标频率单模选择：计划与执行记录

日期：2026-08-02
状态：已完成

## 目标

撤回 finite/full-quarter 的空间基模包络判断。默认只保留一个与目标频率最接近的有效
mode 进入完整后处理，其余 mode 仅保留原始导出和选择审计。

## 选择规则

1. 从 `eigenfrequencies.csv` 读取所有 `is_valid=True` 的候选 mode；若配置了显式 mode
   子集，先限制候选集。
2. 对每个候选计算 `abs(frequency - target_frequency)`。
3. 选择绝对差最小的唯一 mode；完全相同时选择较小的 `mode_idx`，保证可复现。
4. 不读取 Hz 空间分布，不计算中心峰值、壳层位置或向外回升，也不存在空间门槛失败。

## 输出与兼容

- 开关：`FINITE_SINGLE_MODE_ANALYSIS_ENABLED=True`。
- 目标：`FINITE_SINGLE_MODE_TARGET_FREQUENCY_THZ`，默认读取 cavity 的目标 Gamma 频率。
- 审计：`10_overview/analysis_mode_selection.csv`，记录每个候选的频率、目标、绝对差、
  `analysis_selected` 和选择原因。
- 关闭开关时恢复所有有效 mode 的既有后处理。
- far-field、标准 cell map、finite Fourier 和 scoring 仍只接收选中的 mode index。

## 验证与回滚

- 单元测试覆盖最近频率选择、显式候选子集和关闭开关。
- 用失败现场 `shiftx0.050_shifty0.080` 的现有 `.staging` 数据验证，不重新运行 COMSOL。
- 运行聚焦 pytest、`py_compile` 和 `git diff --check`。
- 回滚时恢复全有效 mode 分析；不恢复已由用户撤回的空间包络判据。

## 执行记录

- 已删除空间包络指标、门槛和对应选择路径；full finite 与 finite-quarter 共用最近目标
  频率选择器。
- 失败现场目标为 `198.38683823113445 THz`：mode 0 距离 `1.187373 THz`，mode 1
  距离 `0.239334 THz`，因此选择 mode 1。
- 复用现场已有 COMSOL 导出和已完成 far-field，仅重新执行 mode 1 的 finite Fourier、
  scoring 与 finalizer；未重新启动 COMSOL。正式结果已归档，`.staging` 已移除，mode 0
  没有 far-field/Fourier 目录。
- 现场旧 `fundamental_mode_screen.csv` 已移除，正式审计为
  `10_overview/analysis_mode_selection.csv`；case config 已同步为新选择合同。
- 聚焦测试 `54 passed`，四个相关模块通过 `py_compile`，`git diff --check` 通过。
- 全仓测试为 `368 passed, 1 failed`；唯一失败是当前参数扫描只保留一个 `0.05`
  factor，而既有测试固定读取第 6 个 factor，属于当前参数文件与旧测试假设不一致，和本次
  单模选择无关，未覆盖用户当前扫描参数。
