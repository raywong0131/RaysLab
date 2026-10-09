# 2026-09-01 finite-quarter 计划任务

## 实施前方案

### 目标

在不覆盖既有
`finite_quarter_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5`
结果的前提下，创建一个今天 20:00 启动的一次性 Windows 计划任务。任务运行
finite-quarter 第一象限求解，统一使用基准频率 `197.98 THz`，依次计算 symmetry
ID `1, 2, 3, 4`。

### 固定参数

| 项目 | 值 |
| --- | --- |
| cavity | `b0_nm=245`, `eta=0.96`, `zeta=1.156` |
| cladding | `b0_nm=242`, `eta=0.98`, `zeta=0.93` |
| cavity/cladding layers | `20 / 20` |
| finite geometry | `hex` |
| mesh | `5`，手动 finite mesh |
| shift | `x=0`, `y=0`，uniform/groups |
| eigenmode count | `20` 对；每个 symmetry ID 预期约 `40` 个本征值 |
| frequency shift | `197.98 THz` |
| symmetry IDs | `[1, 2, 3, 4]` |
| resume | `false`，新结果根目录，不复用旧 MPH |
| 后处理 | 单模式选择、far-field FFT、finite-lattice Fourier、field/geometry plots 开启；infinite unit-cell 分支关闭 |

### 数据流与输出

任务级参数快照为
`comsol_task_20260901_finite_quarter_f197_98_eigen20_parameter.json`。
启动器调用现有 `scripts.run_main.run_finite_quarter`，只在进程内设置
`SYMMETRY_IDS=[1,2,3,4]` 和新的输出根，不复制求解实现。

正式结果根目录为：

```text
scripts/.out/finite_cavity/
finite_quarter_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5_f197.980_eigen20_20260901/
```

四个 symmetry case 将位于 `shift0.000/symmetry_*` 下；启动器日志位于结果根的
`task_launcher.log`。旧目录保持原样。

## 实际执行记录

- 参数快照已创建，且通过 `python -m json.tool` 校验。
- 新结果根目录已创建；当前仅有预检日志，没有启动 COMSOL 求解。
- sandbox 外的 Python/几何预检已通过：几何参数、层数、mesh、197.98 THz、20
  对、四个 symmetry ID、`shift_pairs=((0.0, 0.0),)` 和 7422 个第一象限孔均已解析。
- Windows 计划任务已注册并核验为 `Ready`，下一次运行时间为
  `2026/9/1 20:00:00`；20:00 前未启动正式 COMSOL 求解。

## 验证与回滚

计划任务已核验任务名、触发时间、交互式账户、状态和启动器路径。任务启动后以
`task_launcher.log`、各 case 的 `00_model/comsol_progress.log`、
`99_config/config.json` 和 `run_summary.json` 追踪；不把任务创建误报成计算完成。

如需取消本次尚未运行的任务，仅删除对应计划任务即可；不删除新结果根或旧结果。
