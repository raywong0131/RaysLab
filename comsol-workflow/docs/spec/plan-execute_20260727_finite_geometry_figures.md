# Finite / finite-quarter 几何图统一：计划与执行记录

日期：2026-07-27
状态：已完成

## 实施前方案概述

### 目标

统一 `run_finite.py` 与 `run_finite_quarter.py` 的模型级几何图输出。两种计算都只保留以下三张完整结构图，并与 MPH 模型及其进度日志一起写入每个 case 的 `00_model/`：

1. `full_finite_simulation.png`：保持现有绘图内容不变；
2. `full_lattice_modulation_vectors.png`：删除 strip-cut 参考线，并将 cladding shift 箭头显示倍率改为 `×10`；
3. `full_lattice_normal_fan_regions.png`：删除 strip-cut 参考线。

不再由 finite 或 finite-quarter 主流程生成 `full_lattice_modulated_holes*.png` 与 `full_finite_simulation_annotated.png`。

### 范围

- 修改 `scripts/run_main/run_finite.py` 的公共绘图函数和 full finite 接线；
- 修改 `scripts/run_main/run_finite_quarter.py`，复用同一组绘图函数生成完整结构图；
- 更新现有输出布局测试与用户文档；
- 更新最新 finite-quarter `shift0.060/00_model/` 的候选图，删除本次预览后确认不再保留的三张图；
- 不启动 COMSOL，不修改 MPH、场导出、far-field、Fourier 或模式评分结果。

### 数据流

```text
parameter.json + 当前 shift/profile
  -> build_full_lattice_holes()
  -> 完整 cell / hole records 与 finite 六边形边界
  -> 公共三图生成函数
  -> case/00_model/
       full_finite_simulation.png
       full_lattice_modulation_vectors.png
       full_lattice_normal_fan_regions.png
```

finite-quarter 的 eigensolve 仍只使用第一象限裁剪后的 records；完整 records 仅用于模型级几何图，不进入 COMSOL quarter geometry。

### 兼容策略

- 保留现有公共绘图函数，避免复制绘图实现；
- `full_finite_simulation.png` 的函数和参数保持不变，只改变目标目录；
- 不迁移或删除历史 case 中既有结果；只按用户明确要求整理本次最新结构中刚生成的候选图；
- shift-profile 诊断仍保存在 `10_overview/`，汇总 CSV、k-grid 和 BZ 图的归类不变；
- geometry-only quarter 调用也应生成三张完整结构图，但不创建或运行 COMSOL 模型。

### 验证

1. 测试 full finite 与 finite-quarter 都把三张图写入 `00_model/`；
2. 测试两个主流程不再请求生成其余三张候选图；
3. 测试 vector 图箭头倍率为 `10`，且 vector/normal-fan 图不再调用 strip-cut 线；
4. 运行相关 pytest、`py_compile` 与 `git diff --check`；
5. 对最新结构重绘并逐张打开检查坐标比例、边界、图例和文字。

### 回滚

回滚本次对两个 runner、测试和文档的增量修改即可恢复旧行为。已存在的 MPH 与数值结果不受影响；最新 case 中的三张 PNG 可从共享参数和 shift-profile 重新生成。

## 执行完成后记录

### 实际改动

1. `scripts/run_main/run_finite.py`
   - 将 cladding shift 箭头显示倍率从 `×4` 改为 `×10`；
   - 从 modulation-vector 与 normal-fan 图中移除 strip-cut 参考线；
   - 新增公共入口 `save_finite_model_geometry_figures()`，固定且只生成三张模型级完整结构图；
   - full finite 主流程改为调用该入口，并把三张图写入 `00_model/`；
   - 取消 full finite 主流程对 annotated simulation 与 modulated-hole 两组图的生成调用。
2. `scripts/run_main/run_finite_quarter.py`
   - `build_quarter_case_geometry()` 同时保留完整 records 与第一象限 records；
   - 第一象限 records 继续且仅用于 quarter COMSOL geometry；
   - 完整 records 只传给公共几何绘图入口，使 quarter 与 full finite 生成完全相同的三张完整结构图；
   - 输出目标统一为每个 symmetry case 的 `00_model/`。
3. 测试与文档
   - 更新 finite 输出与 quarter geometry-only 测试，验证固定文件集、自然路径和完整/quarter records 分流；
   - 更新根 README、scripts README 与 MEMORY，固化三图名称、目录、箭头倍率和无 strip-cut 规则；
   - 保留并兼容实施前工作区中尚未提交的 cladding shift-profile 改动。
4. 最新结果整理
   - 在最新 `finite_quarter_10-10_cav(245-0.96-1.155)_clad(240.2-0.967-1)_mesh9_shiftprof-tanhpow-l4-p2/shift0.060/00_model/` 重绘三张正式图；
   - 删除本次预览阶段生成、最终不保留的 `full_finite_simulation_annotated.png` 与两张 `full_lattice_modulated_holes*.png`；
   - MPH、progress log 与所有数值结果未修改。

### 测试与图像检查

- `py_compile scripts/run_main/run_finite.py scripts/run_main/run_finite_quarter.py`：通过；
- 聚焦测试：`50 passed`；
- 全量 pytest：`269 passed, 3 warnings`；三条 warning 均为既有 `test_strip_1d_fourier_quotient.py` 测试函数返回非 `None`；
- `git diff --check`：通过；
- 三张最新 PNG 均已逐张打开检查：物理纵横比正常、文字与图例未裁切；vector 图的图例和底部注释均显示 `×10`，vector 与 normal-fan 图中均无 strip-cut 线；`full_finite_simulation.png` 的绘图内容保持原样；
- pytest 临时目录已清理。

### 遗留问题

- 历史计算目录不会自动迁移或清理；这是为了遵守不修改既有计算结果的兼容原则。需要时可以从其保存配置单独重绘。
- `save_full_lattice_plot()` 与 annotated simulation 的底层函数仍保留，供显式诊断或其他工作流复用；finite 与 finite-quarter 主入口不再调用它们。

### 最终状态

实现、测试、文档和最新结果重绘均已完成。后续 full finite 与 finite-quarter 新计算会在 `00_model/` 生成一致且仅包含三张的模型级完整结构图集。
