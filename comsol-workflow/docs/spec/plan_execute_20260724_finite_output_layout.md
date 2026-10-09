# Finite 正式输出框架：计划与执行记录

日期：2026-07-24

## 实施前方案概述

### 背景与目标

当前 `run_finite.py` 与 `run_finite_quarter.py` 仍按旧工作流将 COMSOL 场导出、
far-field FFT、finite-lattice Fourier、评分结果和配置散写到 case 根目录及多个
处理类型目录中。现有结果已经人工整理出一套确认可用的正式框架，本次将该框架固化
到两个 main 入口，使后续新计算直接生成一致结果，不再依赖人工搬运。

正式 case 框架为：

```text
shiftX.XXX/
  00_model/
    finite_cavity.mph 或 finite_quarter.mph
    comsol_progress.log
  01_results/
    mode0/
      10_overview/
        eigenfrequency.csv
        mode_score.csv（若存在）
      11_simulation_exports/
        Hz_center.parquet
        Hz_Re_2d.png
        Hz_Im_2d.png
        Wem_2d.png
        E_air.parquet（有效模式且启用 far-field 时）
      12_farfield_FFT/（若产生结果）
        A_overview.png
        B_polarization.png
        C_cutlines.png
        D_gauss_fit.png
        farfield_summary.json
      13_lattice_fourier_Hz/（若产生结果）
        finite_lattice_fourier_hz.npz
        k_weight_first_bz.png
        p_subspace_mode_decomposition.csv
        top_k_peaks.csv
        k_weight_tops_norm.png
        k_weight_tops_phase.png
    mode1/
      ...
  10_overview/
    几何诊断图
    eigenfrequencies.csv
    mode_scores.csv
    farfield_summary.csv
    finite_lattice_fourier_summary.csv
    index_r_cavity.png
    index_k_1stBZ.png
  80_logs/
  99_config/
    config.json
    air_field_metadata.json（若存在）
    farfield_config.json（若存在）
    objective.json（若存在）
```

series 根目录继续直接保存 `run_summary.json`。mode 使用自然编号 `mode0`、
`mode1`、...、`mode10`，不补零。`11~13` 内文件直接平铺，不增加
`data/metrics/figures` 子目录。

### 实施范围

1. 在 `scripts/run_main/run_finite.py` 中增加可复用的 finite case 路径、临时输出
   和最终整理函数；不新增单独模块，以控制文件数量。
2. `run_finite.py` 使用该函数生成完整 finite cavity 的正式输出。
3. `run_finite_quarter.py` 复用同一函数。默认 `SYMMETRY_IDS=[1]` 时直接使用
   `shiftX.XXX/`；非默认或多 ID 时，在 shift 下增加 symmetry case 层以避免覆盖。
4. 更新现有测试，覆盖目录结构、mode 自然编号、文件归类、JSON 路径修正、公共图
   去重和临时目录清理。
5. 不启动 COMSOL，不重算数值结果，不再次迁移现有正式结果。

### 数据流

```text
COMSOL 与既有后处理模块
  -> case 内专用 .staging/
     simulation_exports/
     farfield_fft/
     finite_lattice_fourier_hz/
     mode_scores.csv / objective.json
  -> 全部计算成功
  -> 预检所有源文件和目标冲突
  -> 按 mode 整理到 01_results/modeN/10~13
  -> 公共汇总进入 10_overview
  -> 配置与元数据进入 99_config
  -> 修正结果 JSON/CSV 中的最终路径
  -> 删除空的 .staging
```

MPH checkpoint 和 `comsol_progress.log` 从计算开始即写入 `00_model/`，不经过
staging，保证长计算监控与 resume 路径稳定。几何诊断图直接写入 shift 的
`10_overview/`。

### 归类规则

- `10_overview`：跨 mode 汇总，以及只依赖几何、BZ、基底或采样网格的公共信息；
- mode `10_overview`：该 mode 的频率行和评分行；
- `11_simulation_exports`：COMSOL 或对称恢复直接产生的该 mode 场数据与场图；
- `12_farfield_FFT`：该 mode 的 far-field 派生结果；
- `13_lattice_fourier_Hz`：该 mode 的 lattice Fourier 派生结果；
- `99_config`：配置与元数据，且按已确认规则包含 `objective.json`；
- `80_logs`：除 `comsol_progress.log` 外的日志；当前流程没有其他持久日志时允许为空。

只创建确实产生结果的 `12`、`13` 目录。公共
`index_r_cavity.png` 和 `index_k_1stBZ.png` 只在 shift
`10_overview` 各保留一份，不复制进 mode。

### 兼容策略

- 底层 `farfield_fft.py` 和 `lattice_fourier_postprocess.py` 的输入输出接口保持不变，
  继续写 staging；由 main 层统一整理，避免影响独立分析脚本。
- `export_full_finite_results()` 与 quarter 对称恢复导出的原始命名保持不变，整理时
  去掉 `00_` 等 mode 前缀后进入对应 mode。
- `farfield_summary.json`、`farfield_summary.csv` 和
  `finite_lattice_fourier_summary.csv` 内的 staging 路径改写为最终路径。
- resume 使用 `00_model/*.mph` 与 `99_config/config.json`。若发现未清理的旧 staging
  或已有目标冲突，先报错，不静默覆盖正式结果。
- 现有分析模块仍可独立使用旧平铺接口；本次只改变两个 main 入口的最终归档行为。

### 验证计划

- 使用临时目录构造两位数 mode（如 `mode11`），验证自然编号和去前缀逻辑；
- 验证 mode 根目录只含目录，`11~13` 内无额外子目录；
- 验证公共 lattice 图只出现在 shift `10_overview`；
- 验证 `objective.json` 位于 `99_config`；
- 验证 far-field JSON 与两个 aggregate CSV 的路径都指向存在的最终文件；
- 验证没有结果的分析目录不创建；
- 验证整理完成后 `.staging` 被删除；
- 运行相关 pytest、`py_compile` 与 `git diff --check`，不启动 COMSOL。

### 回滚方案

代码回滚只涉及本次新增的路径/整理函数、两个 main 入口接线和测试。由于本次不迁移
现有正式结果，回滚不需要恢复计算数据。整理函数在执行移动前完成源文件白名单与目标
冲突预检；出现未知文件时保留 staging 并中止，不做部分覆盖。

## 执行完成后记录

### 实际改动

1. `scripts/run_main/run_finite.py`
   - 新增 `mode_centric_v6` 正式目录常量和 case 路径解析函数；
   - 新增 case 准备、staging 白名单预检、目标冲突检查和最终整理函数；
   - COMSOL 导出与既有 far-field/Fourier/评分模块继续写 `.staging/`；
   - 将 eigenfrequency、mode score、场导出、far-field、lattice Fourier、公共图、
     配置和 objective 按正式框架整理；
   - mode 目录改为不补零的自然编号；
   - 改写 per-mode far-field JSON、far-field aggregate CSV 和 lattice aggregate
     CSV 中的最终路径；
   - MPH、progress log、config 和几何概览从开始即写入最终稳定位置；
   - 正式结果存在、staging 非空或发现未知文件时拒绝静默覆盖；
   - 完整 finite 主流程和 `run_summary.json` 已接入正式布局。
2. `scripts/run_main/run_finite_quarter.py`
   - 复用 full finite 的路径与整理实现，没有复制新的整理逻辑；
   - 默认 `SYMMETRY_IDS=[1]` 直接输出到 `shiftX.XXX/`；
   - 多 ID 或非默认 ID 输出到独立 symmetry case；
   - quarter 仍只保存重建后的完整场结果，MPH 保存到 `00_model/`；
   - geometry-only 模式也建立正式固定目录，但不创建 staging。
3. 测试与文档
   - 在既有测试文件中增加 mode0/mode11、平铺目录、公共图去重、objective 位置、
     JSON/CSV 路径、未知 staging 文件原子拒绝、staging 清理和 quarter 默认目录测试；
   - 更新根 README、scripts README 和 MEMORY，固化正式规则；
   - 未新增额外 Python 模块，符合减少文件数量和复用现有入口的原则。

### 验证结果

- `py_compile`：`run_finite.py`、`run_finite_quarter.py` 通过；
- 聚焦输出、far-field、Fourier 和既有 finite 流程测试：`45 passed`；
- 全量测试：`218 passed, 5 failed, 3 warnings`。5 个失败均为本次改动前已存在的
  参数基准不一致：测试仍硬编码旧参数
  `cav(245-0.96-1.156)_clad(242-0.98-0.930)`，当前 `parameter.json` 为
  `cav(245-0.96-1.155)_clad(243.7-0.98-0.928)`；本次未越权修改这些物理参数测试；
- 测试使用专用 `--basetemp`，结束后已验证并清理过程目录；
- 未启动 COMSOL，未修改或重新迁移现有正式计算结果。

### 遗留问题

- 底层独立后处理 API 仍保留旧的平铺输入/输出接口，这是有意的兼容策略；只有两个
  finite main 入口负责最终整理。
- 旧参数硬编码测试应在单独任务中改为与共享参数配置同步，或明确固定成独立 reference
  fixture；不属于本次输出协议重构。
- 实际 COMSOL 长计算尚未执行；首次正式运行应核对输出目录为空、resume 设置和
  后处理开关，并在 mesh 完成后确认 `00_model/comsol_progress.log` 正常更新。

### 最终状态

代码实现、无 COMSOL 测试和文档记录已完成。正式输出协议为
`mode_centric_v6`，现有人工整理结果无需再次移动。

## 2026-07-29 `A_overview.png` 标准追加更新

`12_farfield_FFT/A_overview.png` 的四个面板固定为：

1. `Near field intensity (NA=XX)`；
2. `Far field Intensity (k-space, linear)`，`kx`、`ky` 均固定为 `[-0.2,+0.2]`；
3. `Far field intensity (1st BZ, log)`，沿用原第一 BZ 面板范围；
4. `Far field Intensity (k-space, log)`，沿用中心及六个相邻 BZ 面板范围。

图级标题固定为 `f=XXX THz, Q=XXX`，不再包含 mode 编号。两张 log 图共享
`[-5,0]` 色标，并将刻度显示为 `10^-5` 至 `10^0`；两图均叠加 NA 圆，NA 圆使用
低对比灰色虚线，BZ 边界保持白色实线。原始 `|E|² output` 面板不再输出。

共享绘图代码和 far-field 回归测试已同步更新。指定的 mode 1 样例只原位覆盖
`A_overview.png`，B、C、D 图未改写；15 项 far-field 聚焦测试通过，未启动 COMSOL，
随后包含 finite 整理流程在内的 26 项聚焦回归全部通过，相关入口 `py_compile` 与
`git diff --check` 通过；未创建新输出目录或额外结果文件。
