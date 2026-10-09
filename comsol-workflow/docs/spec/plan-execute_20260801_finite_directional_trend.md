# Finite directional cladding-shift trend 兼容与执行记录

## 方案概述

扩展既有 `scripts/analysis/finite_trend.py`，使其在保留历史 scalar
`shiftX.XXX` 分析能力的同时，能够读取 finite / finite-quarter ellipse 系列使用的
`shiftxX.XXX_shiftyY.YYY` case。分析只读取已经保存的 `run_summary.json` 和
`10_overview/mode_scores.csv`，不读取 MPH、不启动 COMSOL，也不改写已有求解结果。

本次指定数据集为：

```text
finite_quarter_10-10_cav(245-0.96-1.156)_clad(240.1-0.967-1)_mesh5_ellipse-shift
```

## 范围

- 修改既有 `scripts/analysis/finite_trend.py`；
- 在既有 `tests/test_analyze_finite_trend.py` 中补充 directional 回归测试；
- 更新 `scripts/README.md` 的 trend 输入和横轴说明；
- 在指定 series 的既有结果树下生成 `finite_trend/` CSV 与两张 PNG；
- 不修改 finite/quarter 求解入口、共享参数、几何实现或 COMSOL 模型。

## 数据流与指标定义

```text
run_summary.json + series 根目录全部 shift case
  -> 按 (x factor, y factor, symmetry_id) 合并、去重
  -> 读取各 case 的 10_overview/mode_scores.csv
  -> 筛选 gamma_p_px_weight_fraction > 0.80
  -> 选择频率绝对距离 target frequency 最近的 mode
  -> finite_trend.csv
  -> py-mode weight/frequency 图 + Q 图
```

- ratio 扫描中 `shift_factor` 与 `cladding_x_shift_factor` 都表示基准 factor；
- directional CSV 额外保存 `cladding_x_shift_factor`、
  `cladding_y_shift_factor` 和可定义时的 `cladding_y_over_x_shift_ratio`；
- directional 图以 x factor 为横轴；同一 series 内存在多个 y/x ratio 时分线绘制，
  `(0, 0)` 原点作为各 ratio 曲线的共同起点；
- 同一 x、不同 y 是两个独立 case，不按重复值丢弃，也不把不同 ratio 强行连成一条线。

## 兼容策略

- 历史 scalar case 的字段、路径解析、CSV 必需列和画图入口保持不变；
- directional 元数据与目录名同时存在时逐轴核对，发生不一致则拒绝静默分析；
- `cladding_shift_base_factor` 存在时必须与 x factor 一致；显式 x/y pair 没有 base
  时以 x factor 作为 trend 基准；
- 仅 x 或 y 为 0 的合法 pair 保留，pair identity 仍由两个方向共同确定；
- series 根目录的补充 case 继续按既有规则并入分析，因此不会只依赖最后一次
  `run_summary.json` 而遗漏先前完成结果。

## 验证与回滚

- 聚焦 pytest 覆盖 legacy scalar、directional 同 x/不同 y、stale path、严格权重阈值、
  多 symmetry 和两张图输出；
- 对相关 Python 文件运行 `py_compile`；
- 对任务文件运行 `git diff --check`；
- 正式导出后检查 CSV 行数、所选 mode 和两张 PNG 的实际分辨率成图；
- 如需回滚，只需撤销本任务对 trend/测试/README 的修改并移除新生成的
  `finite_trend/` 分析目录；原始 case 结果不受影响。

## 实际改动与执行结果

- `finite_trend.py` 已同时支持 scalar 与 directional 路径/元数据；
- directional 去重键改为 `(x, y, symmetry_id)`，解决 `x=0.05` 下同时存在
  `y=0.08` 与 `y=0.10` 时的冲突；
- directional CSV 写入 x/y/ratio 审计列，图中 ratio=1.6 与 ratio=2 分线显示；
- 聚焦测试：`11 passed`；
- `py_compile`：通过；
- `git diff --check`：通过，仅报告 Windows 工作树 LF/CRLF 提示；
- 指定 series 共分析 8 个完成 case，全部选择 `mode_idx=1`，无缺失文件或无效模式告警；
- 输出：`finite_trend/finite_trend.csv`、`finite_trend.png`、
  `finite_q_trend.png`；两张 PNG 均已人工检查。

## 遗留问题与最终状态

本次没有启动 COMSOL，也没有创建新的求解 case。ratio=2 主扫描的频率随 x factor
单调上升、Q 单调下降；series 中历史显式点 `(0.05, 0.08)` 作为 ratio=1.6 分支单独
显示。所有请求范围内工作均已完成。

## 仅 ratio=2 的 trend 与 flatness 补充执行

用户进一步要求只分析 `y/x=2`。为避免手工删行导致不可复现，`finite_trend.py`
增加可选参数 `--y-over-x-ratio`：指定 ratio 时仅保留匹配的 directional pair，同时将
`(0,0)` 保留为扫描基线。使用 `--y-over-x-ratio 2` 原位重建本 series 的
`finite_trend.csv` 与两张 trend 图，最终清单为 7 个 case：

```text
(0.00, 0.00)  (0.01, 0.02)  (0.02, 0.04)  (0.03, 0.06)
(0.04, 0.08)  (0.05, 0.10)  (0.06, 0.12)
```

随后以这份 7 行 CSV 作为唯一模式清单运行 `finite_flatness_trend.py`，在
`finite_trend/flatness_trend/` 生成 summary/shell/cell 三张 CSV 与六张 PNG。
ratio=1.6 的 `(0.05,0.08)` 没有参与选模、trend 或 flatness。

补充验证结果：

- trend + flatness 聚焦测试：`19 passed`；
- 相关 Python 文件 `py_compile` 通过；
- `git diff --check` 通过，仅有 Windows LF/CRLF 提示；
- 7 个 case 的 flatness 输入全部存在；最大 Parseval 误差 `7.72e-16`、最大重构误差
  `8.36e-16`、最大 parquet 插值归一化 RMSE `0.00957`；
- 两张 trend 图和六张 flatness 图均已按实际导出分辨率检查。
- flatness 图会根据 `shift_kind` 自动使用 directional 横轴
  `CLADDING X SHIFT / A`，历史 scalar 输入仍保留原标签。

补充数据观察：

- cavity intensity CV 随 x 从 `45.06%` 单调降至 `8.24%`，但 Q 同时从
  `1458.58` 单调降至 `851.35`；
- `(0.05,0.10)` 是最接近目标频率的已算点，py-mode 权重也在此达到主扫描最大值
  `0.96576`；
- `(0.05,0.10)` 的第一层 cladding/cavity 最外层均值比为 `0.9993`，仍无向外回升；
  `(0.06,0.12)` 该比值升至 `1.0543`，且局部 cladding 热点达到 cavity 均值的
  `1.4446`，说明界面护栏已明显恶化；
- 因此当前 ratio=2 的 Pareto 转折点位于 `(0.05,0.10)` 附近；若继续优化目标频率，
  建议在 x `0.050–0.055`、y `0.100–0.110` 间加密，而不是直接继续增大到 `0.06`。

## Flatness canonical cell-map 提升计划

用户确认此前批准的中心 cell 归一化线性全结构图应成为标准输出。实施范围限定为
`finite_flatness_trend.py` 的输出契约、对应聚焦测试和文档，不改变 flatness 数值指标、
选模或 COMSOL 结果。

标准绘图合同：

- canonical 文件：`cavity_cladding_intensity_maps.png`；
- 兼容别名：`intensity_maps_normH.png`，与 canonical PNG 内容完全相同；
- 数据：每个 case 的 unit-cell mean `Hz` intensity，以该 case 的 cavity 中心 cell
  `I0` 归一化；
- panel：按 `shift_factor` 升序，使用紧凑自适应网格，x-y 几何等比例；
- colormap：`RdYlBu_r`；
- normalization：线性 `0–2`，`1` 为色标中点，所有 panel 共享；
- colorbar ticks：`0, 0.5, 1.0, 1.5, 2+`；
- 不标记或加粗中心 cell，不做分位数裁剪或逐 panel 归一化。

兼容与回滚：标准流程不再用旧 log-scale 全结构图写入 canonical 文件，但保留
`save_full_cell_maps()` 函数供显式诊断调用；如需回滚，恢复 `run_analysis()` 对该函数的
调用即可。实施后将用当前 ratio=2 的 7 行 selected-modes CSV 原位重跑 flatness，检查
canonical 与兼容别名逐字节一致，并人工检查全部输出 PNG。

### Canonical 提升执行结果

- `run_analysis()` 现由 `save_full_cell_maps_linear_center_normalized()` 直接写入
  `cavity_cladding_intensity_maps.png`，旧 log-scale `save_full_cell_maps()` 已退出标准
  调用链但函数仍保留；
- `intensity_maps_normH.png` 由 canonical 文件复制得到，保持兼容；
- colormap、normalization、标题、colorbar 标签/刻度/刻度文字和两个输出文件名均提升为
  共享常量，并新增实际 Matplotlib figure 合同测试；
- `scripts/README.md` 已加入 flatness 分析入口与 canonical 输出说明；
- trend + flatness 聚焦测试：`20 passed`，`py_compile` 与 `git diff --check` 通过；
- 当前 ratio=2 的 7 个 case 已原位重新运行 flatness；summary 仍为 7 行且唯一非零 ratio
  为 `2.0`；
- canonical 与兼容别名均为 `889264` bytes，SHA-256 同为
  `ADEBB7BD52482E77B70D745DB8F7ACF8D1A7FC7B266106DBAFCCCB45AA329F96`；
- 六张 flatness PNG 均已按实际分辨率检查，canonical 图满足共享线性 `0–2`、中点 `1`、
  顶端 `2+`、等比例与 7-panel 合同。
