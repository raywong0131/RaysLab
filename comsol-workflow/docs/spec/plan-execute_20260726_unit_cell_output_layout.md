# Unit-cell 正式输出布局：计划与执行记录

日期：2026-07-26

## 方案概述（实施前）

当前 `run_unit_cell_band.py` 仍直接在系列根目录生成 `cavity/`、`cladding/`、
`config.json` 和比较图。既有正式结果采用 `01_results/10_overview/99_config`
布局，但该布局此前依赖计算完成后的人工整理。用户已批准：

1. 当前完整 unit-cell band 计算完成后，将结果整理为正式布局；
2. 修改主程序，使后续运行从一开始就生成正式布局。

当前 COMSOL 进程继续使用旧式活动路径，实施期间不得移动其正在读写的目录。

## 正式目录与数据流

```text
unitcell_band_<parameter-label>/
  01_results/
    cavity/k_points/
    cladding/k_points/
  10_overview/
    cavity/       # geometry、connected bands、composition、Q fits、selected bands
    cladding/
    unit_cell_band_comparison.png
    p_bands_q_vs_k_comparison.png
  80_logs/        # 仅在后台运行确实保留日志时存在
  99_config/
    config.json
    cavity/config.json
    cladding/config.json
```

主程序在运行开始时创建需要使用的正式目录；k-point 求解和缓存始终位于
`01_results/<cell>/k_points`，因此中途停止与恢复不会依赖后续搬运。

## 实施范围

- 在 `run_unit_cell_band.py` 内集中定义 unit-cell 布局路径。
- 修改 `run_case()` 和 `run_band_workflow()` 的写入目标，不改变求解、模式选择、
  band tracking、Q 拟合或绘图算法。
- 在同一主程序中提供受控的旧布局迁移函数，仅用于当前已经启动的旧式结果。
- 更新现有测试，验证正式布局、缓存复用及旧布局迁移的文件完整性。
- 当前运行结束后再执行迁移，迁移前检查未知文件和目标冲突；迁移后核对文件数、
  总字节数及主要结果文件。

## 兼容与风险控制

- `run_case(case_name, params, output_root)` 的调用方式保持不变。
- `run_band_workflow()` 的返回 DataFrame 与现有调用方保持不变。
- 旧布局迁移不覆盖已存在的正式目标；出现未知文件、缺失主要产物或目标冲突时
  立即停止。
- 不移动或修改其他 unit-cell 结果目录。
- 不影响当前已加载旧代码的 COMSOL 进程；该进程完成后单独迁移。

## 验证与回滚

- 运行 unit-cell 聚焦测试、相关回归测试、`py_compile` 和
  `git diff --check`。
- 测试两次运行仍能复用完整 k points。
- 对迁移测试及当前真实结果分别比较迁移前后文件数和总字节数。
- 回滚仅涉及 `run_unit_cell_band.py`、对应测试和本记录，不触碰计算数据。

## 实际改动与结果

### 已完成的程序改动

- 在 `run_unit_cell_band.py` 中集中定义 `01_results`、`10_overview` 和
  `99_config` 路径，并由 `run_case()` 直接写入对应目录。
- k-point 原始结果与缓存写入 `01_results/<cell>/k_points/`；几何图、连接能带、
  模式成分、选中能带和 Q 拟合写入 `10_overview/<cell>/`。
- 两张 cavity/cladding 比较图写入 `10_overview/`；全局及分 cell 配置写入
  `99_config/`。
- 增加受控的 `migrate_legacy_unit_cell_layout()`：迁移前拒绝未知项、缺失产物和
  正式目录冲突，迁移后核对文件数与总字节数。
- 更新根 `README.md`，记录正式输出结构和各目录职责。

### 已完成的验证

- `python -m py_compile scripts/run_main/run_unit_cell_band.py`：通过。
- `tests/test_run_unit_cell_band.py` 与相关几何转换回归：85 项通过。
- 新测试覆盖主流程直接生成正式目录、缓存复用，以及旧布局迁移前后文件数和
  总字节数不变。
- `git diff --check`：通过；仅报告工作树既有的 LF/CRLF 转换提醒。

### 当前正式结果迁移

目标结果：
`scripts/.out/unit_cell_band/unitcell_band_cav(245-0.96-1.155)_clad(244-0.94-1)_mesh9/`。

- cavity 和 cladding 均完成 76/76 个 k 点，运行 stderr 为 0 字节，主进程与
  mphserver 正常退出。
- 迁移前检查确认两张比较图、两组 connected bands、selected bands、模式成分、
  Q 拟合和配置文件全部存在。
- 迁移前后文件数均为 2073，总字节数均为 160428064；迁移没有丢失或改写数据。
- 迁移后顶层仅为 `01_results/`、`10_overview/`、`80_logs/` 和 `99_config/`。
- `01_results/cavity/k_points/` 与 `01_results/cladding/k_points/` 均保留 76 个
  k-point 目录；旧式根目录 `cavity/`、`cladding/`、根配置和根比较图均已消失。
- 两张正式比较图现位于 `10_overview/`，每个 cell 的派生表格和几何图位于
  `10_overview/<cell>/`，配置位于 `99_config/`。

最终状态：两项批准内容均已完成；当前结果已整理，后续主程序将直接生成相同的
正式布局，不再依赖人工搬运。
