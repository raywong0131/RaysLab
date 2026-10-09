# Strip 1D 输出结构重构实施记录

## 方案概述

将 `scripts/run_main/run_strip_1d.py` 的正式结果从 case 根目录扁平输出，重构为
已经用现有 10–10/mesh9 结果组验证过并继续简化的 layout v3：模型、跨 mode 汇总、per-mode
结果和配置分别进入 `00_model`、`10_overview`、`01_results/modeN` 与
`99_config`。底层 COMSOL 导出和 strip bulk Fourier 暂存到 case-local
`.staging`，全部成功并通过冲突预检后再整理到正式位置。

计算期间使用 Matplotlib `Agg` 非交互后端，不调用 `plt.show()`；默认不生成几何、
shift-profile 或 Fourier 诊断图。`Hz_Re`、`Hz_Im` 与 `Wem` 场图属于正式 per-mode
结果，继续写入 `11_simulation_exports`，但统一基于已读取的场数据由 Python
`save_field_plot` 生成；不再创建或运行 COMSOL `PlotGroup2D` 与 `Image Export`。

## 范围

- 修改 `scripts/run_main/run_strip_1d.py` 的路径规划、staging、最终化、完成判定和
  `run_summary.json` schema。
- 修改 `scripts/analysis/srip1d_trend.py`，优先读取 layout v3，同时保留旧扁平结果
  的只读兼容。
- 更新 strip 输出布局、扫描、trend 聚焦测试，以及 README、AGENTS 与 MEMORY 中的
  正式规则。
- `run_strip_1d.py` 不再承载 full-finite 正式运行；full finite 统一由
  `scripts/run_main/run_finite.py` 输出到 `finite_cavity/`。
- 将现有 10–10/mesh9 扫描组和本轮新完成的 20–20/mesh5 case 迁移为 layout v3；
  不迁移或改写其他工作流的正式结果组。

## 输出合同

```text
strip1d_<series>/
  cut_60/                              # 0° 时为 cut_00
    run_summary.json
    strip1d_trend/
    shiftX.XXX/
      00_model/
        strip_1d.mph
        comsol_progress.log
      01_results/                      # 仅 gamma_subspace_p > 0.9
        modeN/
          10_overview/
            eigenfrequency.csv
            mode_score.csv
          11_simulation_exports/
            Hz_center.parquet
            Hz_Re_2d.png
            Hz_Im_2d.png
            Wem_2d.png
            cutline_Wem.png            # x=0，沿 y 轴
          13_strip_bulk_fourier_Hz/
      10_overview/
        eigenfrequencies.csv           # 保留全部 mode 及筛选判断
        mode_scores.csv
        strip_bulk_fourier_summary.csv
      99_config/
        config.json
        objective.json
```

cut 文件夹使用两位角度；case 子目录不重复 cut 信息。case 根和 mode 根只放目录；
mode 使用不补零的自然编号。`gamma_subspace_p > 0.9` 是 strip 正式有效性判断，只有
严格超过阈值的 mode 进入 `01_results`；case 级汇总保留全部本征 mode 和筛选证据。
没有产生结果的分析目录和空 `80_logs` 不创建。

`cutline_Wem.png` 绘图合同：从 center-plane `ewfd.Wav` 数据在 `x=0` 上插值，横轴为
`y (µm)`，纵轴为单模最大值归一化后的 `Normalized Wem`，标题为
`Wem cutline at x = 0`，无图例、无交互窗口，使用克制的单条深色实线并输出 220 dpi。

## 数据流

```text
parameter.json + strip-only controls
  -> cut_<angle>/shiftX.XXX -> prepare_strip_case_output
  -> 00_model 中构建/保存 MPH 与 progress log
  -> .staging/simulation_exports
  -> .staging/strip_bulk_fourier_hz
  -> .staging/mode_scores.csv + objective.json
  -> 以 gamma_subspace_p > 0.9 标记有效 mode
  -> finalize_strip_case_output 预检 mode、未知文件和目标冲突，只保留有效 mode
  -> layout v3 正式目录并删除 .staging
  -> cut_XX/run_summary.json
  -> 可选 strip1d_trend/
```

## 兼容策略

- 新计算只写 layout v3。
- 历史 `cut60_shiftX.XXX` 可迁移为 `cut_60/shiftX.XXX`；迁移前检查目标冲突，保留
  case 级全部汇总，只删除/不复制不符合 gamma 阈值的 per-mode 派生结果。
- 完成判定识别 layout v3，也只读识别历史 layout v1/v2，避免重算既有结果。
- trend 优先读取 `10_overview/mode_scores.csv`，找不到时回退到 case 根
  `mode_scores.csv`。
- `run_summary.json` 使用 `output_layout_version: 3`，并位于对应 `cut_XX/` 中。
- 不覆盖已有 `01_results`、正式汇总或非空 `.staging`。

## 验证方案

- 精确测试目录名、mode 自然编号、文件分类、单模 CSV、缺失分析目录不创建。
- 测试 `cut_60`/`cut_00` 和无冗余 cut case 名；阈值严格使用 `> 0.9`，等于 0.9
  不通过；`01_results` 不包含无效 mode。
- 测试 `cutline_Wem.png` 的 x=0 插值、归一化、标题、坐标标签、文件名与 `Agg` 后端。
- 测试未知 staging 文件、未知 mode、重复目标和已有目标均在移动前失败。
- 测试成功后 `.staging` 清除，失败时保留诊断现场。
- 测试新版与旧版完成判定、trend 双布局读取、`Agg` 后端和关闭诊断绘图。
- 运行 strip 聚焦 pytest、相关 `py_compile` 和 `git diff --check`。
- 实际 COMSOL 启动前单独核对当前参数、输出目录、mesh、中心频率、shift、恢复与
  后处理开关。

## 回滚方案

代码回滚时只回退本任务涉及的源码、测试和文档，不触碰用户其他未提交修改。若新计算
在最终化前失败，保留 `00_model` checkpoint 与 `.staging`；若最终化成功，则正式结果
完整位于 layout v3，不依赖临时目录。

## 执行结果

- 2026-07-31：strip 的 `Hz_Re`、`Hz_Im` 与 `Wem` 图片已从 COMSOL
  `PlotGroup2D`/`Image Export` 改为 Python `save_field_plot`。文件名、目录、频率/Q
  标注、`Hz_center.parquet` 与 `cutline_Wem.png` 合同保持不变，避免每个结构完成后
  因 COMSOL 图片导出弹出绘图窗口；full finite 和共享 COMSOL 导出接口未修改。
- 2026-07-31：`RUN_GAMMA_TREND_AFTER_SCAN=None` 作为默认自动策略；cladding shift
  数量大于 3 时，在 `run_summary.json` 写入完成后自动生成 `strip1d_trend`。显式
  `True`/`False` 可覆盖默认策略，单点计算仍不会运行 trend。
- 2026-07-31：补算少量 shift 时，主流程会重新发现同一 cut 目录内全部完整
  `shiftX.XXX` case，并按 shift 去重、排序后写回 `run_summary.json`；自动 trend
  判断使用合并后的完整点数，避免新一轮少量参数覆盖此前扫描汇总。

- 已在 `run_strip_1d.py` 中实现 layout v3 路径规划、case-local `.staging`、未知文件/
  mode/目标冲突预检、mode-centric 最终化、摘要路径改写和 staging 清理。
- 新计算把 MPH/progress 直接写入 `00_model`，aggregate CSV 写入 case 级
  `10_overview`，配置/objective 写入 `99_config`，per-mode 场与 strip Fourier 分别
  写入 `11_simulation_exports` 和 `13_strip_bulk_fourier_Hz`。
- 完成判定支持 layout v3，并只读兼容历史 layout v1/v2；trend 优先读取 v3 并回退旧路径。
- `RUN_FINITE_FULL=True` 现在明确报错并引导使用 `run_finite.py`。
- Matplotlib 使用 `Agg`，正式默认 `RUN_DIAGNOSTIC_PLOTS=False`；shift profile、几何和
  Fourier 诊断绘图均不进入默认计算路径。正式 `Hz_Re`、`Hz_Im`、`Wem` 场图继续静默
  导出。
- 聚焦验证：99 tests passed；相关 `py_compile` 通过；`git diff --check` 通过。既有
  Fourier quotient 测试保留 3 条 `PytestReturnNotNoneWarning`，与本次修改无关。
- 2026-07-30 已按启动前预检完成新的 20–20、mesh5、cut60、uniform shift0.000 case：
  `scripts/.out/strip_1d/strip1d_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/`
  原始 `cut60_shift0.000/`。COMSOL mesh 共 223371 个单元，最小单元质量 0.1909；
  1588190 自由度的特征频率求解、场导出、Fourier 和评分全部成功，运行时 stderr 为空，
  诊断绘图关闭。四个模式的 `gamma_subspace_p` 分别为 0.0029333、0.9879232、
  0.0016116、0.9951392，因此 layout v3 仅保留 mode1 和 mode3。
- 新增 `cutline_Wem.png`：从已有 MPH 的 `ewfd.Wav` 求值生成，不重新执行 eigensolve，
  使用 `Agg` 后端且不弹窗。
- 迁移 10–10/mesh9 的前三个 case 时，发现 mode1 曾被旧 geometry-valid 规则提前跳过；
  已从现有 `Hz_center.parquet` 补做 Fourier 和 gamma 评分，不启动 COMSOL。三个
  `gamma_subspace_p` 分别为 0.9875324、0.9878282、0.9879621，均严格超过 0.9，
  因此十个 shift case 最终都保留 mode1 和 mode3。
- 共迁移 11 个正式 case，生成 22 张 `cutline_Wem.png`；10–10 的十点 trend 已在
  `cut_60/strip1d_trend/` 原位重建。两组 `run_summary.json` 均为 layout v3，旧
  `cut60_shiftX.XXX`、series 根 trend/summary 和迁移暂存目录均已清理。
- 最终聚焦验证为 99 tests passed；相关 `py_compile`、`git diff --check` 和当前
  20–20 case 的 `strip_case_is_complete(..., 0.0)` 均通过。测试仍有 3 条既有
  `PytestReturnNotNoneWarning`，与本任务无关。实际抽查两张截线 PNG，标题、坐标、
  归一化、边界和版面均正常。
- 最终状态：代码、测试和文档重构完成；正式 strip 结果按预检后的无冲突路径迁移到
  `cut_60/shiftX.XXX`，case 汇总保留全部模式，`01_results` 仅保留 gamma-valid 模式。
