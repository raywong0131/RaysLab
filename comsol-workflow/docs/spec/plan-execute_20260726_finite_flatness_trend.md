# Finite cavity Hz 平整度趋势：计划与执行记录

日期：2026-07-26  
状态：已完成

## 目标

在不启动 COMSOL、不修改已有求解结果和既有 `finite_trend.csv` 选模逻辑的前提下，
读取 finite trend 已选出的逐 shift 目标模式，评价：

1. 全部 cavity unit cell（含最外层）的复数 `Hz_center` Gamma 平整度；
2. 整个 unit cell 上的平均强度 `mean(|Hz|^2)` 平整度；
3. cavity 外各 cladding shell 的归一化强度、边界热点和向外回升；
4. `0.000` 至 `0.090` 十个 shift 的趋势。

## 范围与数据源

- 模式选择唯一来源：现有 `finite_trend/finite_trend.csv`。
- cavity cell 复场来源：逐模式
  `13_lattice_fourier_Hz/finite_lattice_fourier_hz.npz`。
- cladding 场来源：逐模式
  `11_simulation_exports/Hz_center.parquet`。
- cell 中心、shell、周期与层数来源：逐 case `99_config/config.json`。
- NPZ 中的 `rho_points` 作为 cavity 和 cladding 共用的完整 unit-cell 采样点。
- PML、结构外区域和不完整边界 cell 不纳入评价。

## 指标

对 cavity cell 强度 `I_R = mean_rho(|Hz(R + rho)|^2)` 定义：

- `gamma_k_weight_fraction`：NPZ 中 `t=0` 的 Fourier 权重；
- `cavity_intensity_cv = std(I_R) / mean(I_R)`；
- `cavity_intensity_flatness = 1 / (1 + cavity_intensity_cv^2)`；
- `cavity_intensity_dmax = max(|I_R / mean(I_R) - 1|)`；
- `cavity_outer_mean_ratio`：cavity 最外 shell 均值除以全 cavity 均值。

对第 `s` 层 cladding 定义 `g_s = mean(I_R in shell s) / mean(I_R in cavity)`，
并汇总：

- `cladding_hotspot_max = max(I_R in cladding) / mean(I_R in cavity)`；
- `cladding_hotspot_p95`：同一归一化量的 95% 分位数；
- `cladding_upward_variation = sum(max(0, g_(s+1) - g_s))`；
- `cladding_max_shell_cv`：所有 cladding shell 内 CV 的最大值。

cladding 指标是独立护栏，不与 cavity 平整度混为一个区域平均值。

## 数据流与兼容策略

```text
finite_trend.csv（逐 shift 已选 mode_idx）
  -> 定位 case/mode 目录
  -> NPZ 反变换恢复 cavity cell 场并计算 cavity 指标
  -> config.json 读取 cladding cell 与 shell
  -> Hz_center.parquet 去重后线性插值到相同 rho_points
  -> 计算 cladding cell/shell 指标
  -> 输出 summary、shell profile、cell detail CSV 和静态趋势图
```

分析入口只消费 selected-modes CSV，不重新按频率或权重选择模式，因此与现有 finite
trend 保持严格一致。输出放在 `finite_trend/flatness_trend/`，不覆盖已有 CSV 和图片。

## 数值验证

1. NPZ `t=0` 权重必须与 trend CSV 的 `gamma_k_weight_fraction` 一致；
2. NPZ 反变换重构 cavity cell 场；
3. parquet 插值到 cavity cell 后，与 NPZ cavity 强度比较归一化 RMSE；
4. 十个 shift、每个 selected mode 的必要输入必须全部存在；
5. 所有汇总指标必须有限，shell 编号连续；
6. 运行聚焦 pytest、`py_compile`、`git diff --check`；
7. 打开并检查最终图片。

## 回滚

删除本任务新增的分析入口、测试、本文档和新增的
`finite_trend/flatness_trend/` 输出即可；原始求解数据和现有 trend 产物不受影响。

## 执行记录

### 实际改动

- 新增 `scripts/analysis/finite_flatness_trend.py`，直接消费现有
  `finite_trend.csv`，对其中十个 selected mode 计算 cavity 与 cladding 指标；
- 新增 `tests/test_finite_flatness_trend.py`，覆盖六角 shell、CV/flatness 恒等式和
  cladding 回升指标；
- 在当前 series 的 `finite_trend/flatness_trend/` 下生成 summary、shell、cell 三张
  CSV 与五张静态图；
- 未启动 COMSOL，未重选模式，未覆盖任何原始求解结果或既有 trend 产物。

### 验证结果

- `pytest -q -p no:cacheprovider --basetemp tests/.tmp_finite_flatness_pytest_run
  tests/test_analyze_finite_trend.py tests/test_finite_flatness_trend.py`：13 passed；
- `py_compile scripts/analysis/finite_flatness_trend.py`：通过；
- `git diff --check`：通过，仅报告既有用户文件的 LF/CRLF 提示；
- 十组 Parseval 最大误差 `8.06e-16`，重构最大误差 `8.36e-16`；
- parquet 插值回算 cavity 与 NPZ cavity 强度的归一化 RMSE 最大 `0.00270`；
- 五张输出图片均已打开检查，轴、标签、统一色标和十组 small multiples 可读。

首次合并 pytest 使用系统临时目录时，旧测试的 `tmp_path` fixture 因 Windows 临时目录
权限报错；改用仓库内受控 `--basetemp` 后全部通过，临时目录已清理。

### 当前数据观察

- `shift=0.050` 同时给出最高 Gamma 权重 `0.97115` 和最低 cavity intensity CV
  `11.60%`，但最大 cavity cell 偏差仍为 `32.90%`；当前十组没有达到此前提出的
  A 级平整度标准。
- `shift=0.040` 的 Gamma 权重 `0.96470`、CV `12.32%`，但 cladding 最大热点仅为
  cavity 均值的 `1.038`；`shift=0.050` 的该热点已增至 `1.435`。
- 所有 shift 的 cladding shell 均值都向外单调下降，因此异常不是径向 shell 均值回升；
  真正问题是第一层 cladding 的界面过冲和同一 shell 上的方位性热点。
- `shift>=0.060` 后左右边界热点快速增强；`shift=0.090` 的最大 cladding cell 达到
  cavity 均值的 `4.86` 倍，最外 cavity shell 均值也增至 cavity 均值的 `1.467`。
- 从 cavity 平整度、cladding 护栏和 Q 综合看，`0.040–0.050` 是当前 Pareto 转折区，
  下一轮适合在该区间细扫，而不是继续增加 shift。

### 遗留问题

- 当前 cladding 护栏基于 `|Hz|^2`；最终仍需导出可积分的电磁能量密度数据，复用同一
  cell/shell 指标验证界面热点；现有 `Wem_2d.png` 只能目视检查，不能做严格积分。

### 2026-07-27 可视化补充

- 增加以 cavity 中心 cell `(i,j)=(0,0)` 强度 `I0` 逐 case 归一化的全结构 cell map；
- 色值直接使用线性的 `I_R / I0`，十个 panel 共用从 0 到全局最大值的公共色标，不做
  对数变换或分位数裁剪；
- 中心 cell 仅作为数值归一化基准，不在图中使用边框或其他标记；
- 在 `clad(240.2-0.967-1)` 十组数据上生成测试图，公共范围为
  `0–2`，其中 `1.0` 位于色条正中，所有不小于 2 的值使用顶端颜色并将最高刻度标为
  `2+`；聚焦测试 5 项通过。
