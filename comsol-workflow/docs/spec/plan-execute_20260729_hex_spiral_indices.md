# Finite lattice 与 k-grid 零起始同心六边形编号：计划与执行记录

日期：2026-07-29

## 实施前方案概述

### 目标

将 finite lattice 与 finite reciprocal k-grid 的对外编号统一为零起始的同心
六边形顺序：

1. 中心点编号为 `0`；
2. 第一圈从第一象限方向的上方点开始，按逆时针方向编号 `1` 至 `6`；
3. 第二圈从同一方向开始，按逆时针方向编号 `7` 至 `18`；
4. 后续各圈从内向外连续编号，每圈均从第一象限方向的上方点开始并逆时针增加。

本次只改变对外显示和审计编号，不改变 finite-lattice Fourier 使用的循环群索引、
FFT 数组顺序、离散波矢生成、权重、模式分解或重建公式。

### 范围

- 在 `comsol_workflow/finite_lattice_fourier.py` 增加共享的规则六边形点集编号函数；
- 更新 lattice index 与 finite k-grid 诊断图中的数字标记和说明文字；
- 更新 finite `top_p_abs` / `top_p_phase` 面板显示的 k-grid 编号；
- 在 finite 的 CSV、NPZ 和 summary 中增加明确的外部编号字段，同时保留既有
  `cyclic_indices`、`t`、`xi_values` 等内部字段；
- 同步仍保留在 `scripts/run_main/run_finite.py` 中的旧兼容实现，避免手工调用时产生
  不同编号；
- 增加编号方向、逐圈连续性、输入顺序无关性和 Fourier 数值不变性测试。

不改动 strip 1D 的 row/m 编号，因为它不是二维 finite hexagonal k-grid。

### 数据流与索引边界

```text
lattice axial coordinates (i, j)
  -> cyclic_indices x(i,j)       [内部 FFT 重排，保持不变]
  -> lattice_indices             [外部零起始同心六边形编号]

xi_values[t]
  -> folded_k[t]                 [内部 t 顺序，保持不变]
  -> k_grid_indices[t]           [外部零起始同心六边形编号]
```

共享编号函数从规则六边形点集的中心和六个最近邻恢复二维晶格基矢，再将每个点映射为
整数轴坐标。圈层使用 hex distance，圈内使用相对起始方向的极角排序。函数必须验证：

- 点数符合 `3L(L+1)+1`；
- 中心点唯一；
- 每圈恰有 `6L` 个点；
- 点能够在容差内重建到恢复的规则三角晶格；
- 输出是 `0..N-1` 的完整排列。

### 兼容策略

- `hex_cyclic_indices()` 与 `hex_xi_values()` 不改签名、不改结果；
- 既有 NPZ 字段和 CSV 的 `t` 列继续表示内部 FFT 索引；
- 新增 `lattice_indices`、`k_grid_indices`、`k_index`、`top_k_index`，不复用或覆盖
  内部索引字段；
- 正式输出文件名和目录结构保持不变，避免影响 mode-centric 整理逻辑和已有消费者；
- `k_weight_Hz` 的权重位置、颜色、插图及数值均不改变。

### 验证与回滚方案

- 对一至多圈直接晶格验证中心、起点、逆时针顺序和圈起始编号；
- 对 finite k-grid 验证相同规则，并验证标签映射不改变 `folded_k` 和 Fourier 结果；
- 对打乱输入顺序的点集验证编号只由几何位置决定；
- 运行 lattice/Fourier、finite output layout 和 finite workflow 聚焦测试；
- 运行相关 `py_compile` 与 `git diff --check`；
- 如需回滚，只移除新增外部编号字段与显示映射，内部 FFT 数据无需迁移或恢复。

## 执行完成后记录

### 实际改动

1. `comsol_workflow/finite_lattice_fourier.py` 新增 `hex_spiral_indices()`：
   - 从规则六边形点集恢复中心、第一象限上方起始方向和三角晶格基矢；
   - 按 hex distance 从内向外分圈，每圈按极角逆时针编号；
   - 严格验证点数、中心、六个最近邻、整数晶格重建、每圈点数和最终排列。
2. `comsol_workflow/lattice_fourier_postprocess.py`：
   - lattice 诊断图改为显示 `lattice index`，中心为 0；
   - finite k-grid 诊断图使用相同的零起始同心六边形编号；
   - finite `top_p_abs` / `top_p_phase` 面板的 `k=` 改为外部 k-grid 编号；
   - `p_subspace_mode_decomposition.csv` 和 `top_k_peaks.csv` 新增 `k_index`；
   - mode NPZ 新增 `lattice_indices` 与 `k_grid_indices`；
   - mode summary 新增 `top_k_index`。
3. `scripts/run_main/run_finite.py` 中保留的旧兼容实现同步采用同一 helper 和字段，避免
   手工调用旧入口时出现不同编号。
4. 现有 `cyclic_indices`、`t`、`xi_values`、`folded_k`、FFT 调用和输出文件名均保持
   不变；绘图函数还会验证传入的 cyclic index 仍是完整 FFT 排列。
5. `tests/test_lattice_polygon_fourier_utils.py` 增加前两圈精确顺序、输入乱序不变性、
   finite k-grid 同规则以及输入数组不被改写的测试。

### 验证结果

- 聚焦回归：
  `test_lattice_polygon_fourier_utils.py`、`test_strip_1d_fourier_quotient.py`、
  `test_run_finite_farfield.py`、`test_finite_cavity_field_exports.py`，共 23 项通过；
  仅有 3 条既有的 pytest return-value warning。
- `py_compile`：共享 Fourier 模块、共享后处理、full finite 和 quarter finite 入口通过。
- `git diff --check` 通过；仅有工作区既有的 LF/CRLF 提示。
- 对当前 `20-20`、`shift0.000`、mode 1 样例只读核对：lattice 与 k-grid 均形成
  `0..1260` 完整排列；原主要内部 `t = 0, 1, 1260, 1200, 61, 1199` 分别映射为
  外部 `k_index = 0, 6, 3, 1, 4, 2`，权重和排序未改变。
- 在系统临时目录生成 L=3 lattice/k-grid 核对图并完成视觉检查，确认中心为 0、每圈
  从第一象限上方点开始逆时针增加；临时图已删除，未创建项目输出目录。
- 未启动 COMSOL，未改写任何正式计算结果。

### 遗留问题与最终状态

实现和验证完成。内部 cyclic index 仍只用于 FFT；对外使用零起始同心六边形编号。

### 2026-07-29 输出名追加更新

经用户确认，两个 shift 级公共索引图进一步固定为：

- 实空间 cavity lattice：`10_overview/index_r_cavity.png`；
- 第一布里渊区 finite k-grid：`10_overview/index_k_1stBZ.png`。

共享后处理、full finite 保留的旧兼容入口、mode-centric staging 整理映射、输出布局文档
和测试已同步更新。旧的 `lattice_bulk_cyclic_index.png` 与
`lattice_finite_k_grid_first_bz.png` 不再作为正式输出名。

重命名后运行 lattice、finite output layout 与 finite far-field 共 19 项聚焦测试，全部
通过；相关入口 `py_compile` 与 `git diff --check` 通过。指定的 `20-20`、
`shift0.000/10_overview` 样例已用新编号规则重绘为两个新文件，并在确认新文件完整后
删除两个旧名称文件；未创建新目录、未启动 COMSOL。

### 2026-07-29 索引图可读性追加更新

两个索引图进一步统一为同一套 `viridis` 连续编号色标：

- colorbar 通过 axes divider 与主坐标框等高，固定标出最小编号 `0` 和最大编号；
- colorbar 说明文字使用 `rotation=270`；
- k-grid 的所有点按 `k_grid_indices` 连续着色，不再使用单一蓝色；
- 数字字体由原来的 6.2/6.5 pt 缩小为 5.0 pt；
- 仅标记前四层，即中心层和第 1–3 圈，对应编号 `0..36`；更外层仍保留几何、点位和
  色彩编码；
- `run_finite.py` 中的旧兼容绘图函数改为直接调用共享绘图实现，避免再次发生样式分叉。

修改后运行相关 20 项聚焦测试，全部通过；`py_compile` 与 `git diff --check` 通过。
指定样例的两个现有 PNG 已原位重绘并完成视觉检查，未增加目录或其他结果文件。

### 2026-07-29 top-k 模式分布标准追加更新

经用户确认，`13_lattice_fourier_Hz` 中两张 top-k 模式分布图的正式输出名固定为：

- `k_weight_tops_norm.png`；
- `k_weight_tops_phase.png`。

本次追加更新覆盖上文“正式输出文件名保持不变”的历史兼容说明。共享后处理与
`run_finite.py` 保留入口均调用同一绘图实现；finite 模式重绘时会移除旧的
`top_p_abs.png`、`top_p_phase.png`。

幅值图对六个入选 top-k profile 使用同一个全局最大值归一化，色标固定为 `0..1`，
因此全图最大值严格为 1，同时保留不同 k 分量之间的相对幅值。相位图固定为
`[-π, +π]`，仅标记 `-π`、`0`、`+π`。两图的 colorbar 使用跨两行的独立 GridSpec
轴，上端与第一排黑色线框上缘平齐，下端与第二排黑色线框下缘平齐；说明文字统一使用
`rotation=270`，分别为 `Normalized |u_k(r)|` 与 `arg u_k(r)`。

指定的 `20-20`、`shift0.000`、mode 1 样例已在原有
`13_lattice_fourier_Hz` 目录内重绘并完成视觉检查；两个旧名 PNG 已在确认新图完整后
删除。lattice/Fourier、finite output layout 与 finite far-field 共 21 项聚焦测试通过，
相关入口 `py_compile` 与 `git diff --check` 通过；未创建新目录、未启动 COMSOL。

### 2026-07-29 top-k 图配色与标题追加更新

经用户进一步确认：

- `k_weight_tops_norm.png` 使用 `plasma`；
- `k_weight_tops_phase.png` 使用 `twilight_shifted`；
- 总标题分别固定为 `Top weight finite-cavity profiles: Intensity` 与
  `Top weight finite-cavity profiles: Phase`，不再包含 `mode XX:`；
- 各面板标题固定为 `k-index=XX, P(k)=XX`，替代原来的 `k=XX, frac=XX`。

指定样例已在原目录内原位重绘并完成视觉检查，原有归一化、相位范围以及跨两行
colorbar 对齐规则保持不变。新增标题、面板文字和 colormap 回归检查后，相关 22 项
聚焦测试全部通过；`py_compile` 与 `git diff --check` 通过。
