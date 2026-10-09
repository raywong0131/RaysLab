# unit-cell-band V2 正式替换计划与执行记录

## 方案概述

将已经完成验证的 reusable-model V2 实现提升为唯一的
`scripts/run_main/run_unit_cell_band.py`。同一组几何参数只创建一个 COMSOL Model；
k 点变化时复用该 Model，mesh 参数变化时仍复用 Model、仅重建 mesh。

## 范围

- 用 V2 实现替换现有 `scripts/run_main/run_unit_cell_band.py`。
- 删除独立的 `scripts/run_main/run_unit_cell_band_v2.py`。
- 保留唯一命令 `comsol-unit-cell-band`，移除 `comsol-unit-cell-band-v2`。
- 将 V2 专属测试改为唯一主入口的 model-reuse 测试。
- 更新 `scripts/README.md` 中的入口与执行逻辑说明。
- 保持 `scripts/parameter.json` 中 `unit_cell_eigenmode_count = 6` 及其读取方式不变。

## 数据流

```text
scripts/parameter.json
  -> scripts/run_main/run_unit_cell_band.py
  -> 每组几何参数创建一个 reusable SimulationRun / COMSOL Model
  -> 首个 k 点构建 geometry + mesh，后续 k 点更新 Floquet 参数并求解
  -> mesh 配置变化时清理并重建 mesh，继续复用同一 Model
  -> scripts/.out/unit_cell_band/unitcell_band_<参数标签>/
```

## 兼容策略

- 原命令 `comsol-unit-cell-band` 和原 Python 路径保持不变。
- 标准输出前缀恢复为 `unitcell_band_`，不再生成带 V2 前缀的新目录。
- 已存在的 `unitcell_band_v2_*` 结果不移动、不覆盖、不删除。
- 共享参数、模式命名、缓存文件格式和后处理功能保持不变。
- metadata 中执行引擎统一记为 `reused_model`，不再暴露版本后缀。

## 验证方案

- 编译检查唯一主入口。
- 运行 unit-cell 主入口与 reusable-model 聚焦测试。
- 检查 `pyproject.toml` 只保留标准 CLI。
- 搜索活动代码、测试和 README，确认无 V2 入口或 V2 输出前缀残留。
- 执行 `git diff --check`。

本次替换不重复启动昂贵 COMSOL 求解；被提升的 V2 实现此前已经完成实际 COMSOL
验证，本次重点验证文件切换、入口兼容和回归测试。

## 回滚方案

若回归测试失败，可依据本次 diff 恢复原主入口和独立 V2 文件；正式计算结果目录不参与
替换，因此回滚不影响已有结果。

## 执行结果

已完成，最终状态如下：

- reusable-model 实现已成为唯一的
  `scripts/run_main/run_unit_cell_band.py`，独立 V2 文件已移除。
- 主程序按每组几何复用一个 COMSOL Model；缓存点不创建 Model，未缓存点复用同一
  runner；mesh 参数变化时只重建 mesh。
- 保留公开名称 `SimulationRun` 作为 reusable runner 的兼容别名，已有直接调用或测试
  替换点无需改用新的类名。
- 输出目录统一为
  `unitcell_band_<unit_cell_case_parameter_label>`，metadata 的执行引擎统一为
  `reused_model`。
- `unit_cell_eigenmode_count` 继续从 `scripts/parameter.json` 读取，当前值为 6。
- `pyproject.toml` 与同步后的本地环境只保留 `comsol-unit-cell-band`；README 已更新为
  唯一入口及 Model/mesh 复用说明。
- 独立测试已迁移为 `tests/test_run_unit_cell_band_model_reuse.py`，并补充唯一入口、标准
  输出名、Model 复用和 mesh 重建检查。

验证结果：

- `tests/test_run_unit_cell_band.py` 与
  `tests/test_run_unit_cell_band_model_reuse.py`：99 passed。
- 调整最终唯一入口断言后再次执行对应测试：1 passed。
- `py_compile` 通过。
- `uv sync` 通过，console scripts 审计结果仅为
  `['comsol-unit-cell-band']`。
- 活动代码、测试、README 和 `pyproject.toml` 中不存在旧 V1/V2 入口、V2 输出前缀或
  V2 execution-engine 标识。
- `git diff --check` 通过；仅输出工作区既有的 LF/CRLF 提示。

本次未重新运行昂贵 COMSOL 求解，也未移动、覆盖或删除任何正式计算结果。测试产生的
`_pytest_replacement*` 临时目录已清理，无遗留问题。
