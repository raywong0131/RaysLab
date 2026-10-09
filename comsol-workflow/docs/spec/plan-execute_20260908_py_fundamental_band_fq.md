# p_y 基模 band–f-Q 组合图

## 目标

使用现有 unit-cell band 数据和已合并的 finite-cavity 本征频率表，在 finite 结果
`overview` 中生成一张与用户参考插图结构一致的三面板科研图，不启动 COMSOL。

输出文件：
`overview/p_y_fundamental_band_f_Q.png`

保留现有 `unit_cell_band_comparison_ylim185_210_THz_halfM_gamma_halfK.png` 和
`overview/f_Q.png`，不原位覆盖。

## 模式身份

- 用户可见的 `p_y` 对应项目内部 `px`。
- 必须从 finite 结果的 `analysis_mode_selection.csv` 读取唯一
  `analysis_selected=True` 且 `target_internal_mode=px` 的基模，不按最高 Q 猜测。
- 当前记录指向 `symmetry_1_xPEC_yPMC/mode19`：
  frequency = 197.9530009044521 THz，Q = 7676.733033102247，
  `fundamental_eligible=True`，内部 `mode_px=0.9989540631720737`。
- 图中只显示用户标签 `p_y`，不显化内部模式名、symmetry ID 或 PEC/PMC。

## 绘图合同

- 三面板顺序：`Cavity`、`Cladding`、`Quality factor`。
- 按用户参考图匹配画幅：两张 band 面板高/宽约 1.62，Q 面板约 1.73，
  三栏宽度约 1.065:1.055:1。
- 左两栏从 cavity/cladding `selected_bands.csv` 重绘，沿
  `M -> Gamma -> K` 截取 `|k| <= 0.3`，共享 190–210 THz 线性频率轴。
- band 横坐标三个标记依次为 `0.3MΓ`、`Γ`、`0.3ΓK`；不再显示第二排方向文字，
  并为 x/y 标签保留独立间距，禁止文字重叠。
- band 颜色沿用项目合同：用户 `p_y` 红、`p_x` 蓝，
  `d_xy` 浅灰、`d_(x^2+y^2)` 黑，并保留连续混合颜色。
- 右栏使用合并 `eigenfrequencies.csv` 的全部有限 frequency/Q 行；横轴 Q、
  纵轴 Frequency，均为线性坐标。普通模式统一蓝色散点。
- 在三栏相同频率处绘制红色水平虚线；右栏以红点和 `p_y` 标注突出基模，
  标注固定置于 Q 面板右上角。
- 红色基模参考线必须是跨过两个面板间隙的一条连续线，不得分成三段；线置于曲线和散点下方。
- 白底、刻度向内、无网格和冗余图例；只保留 band 成分图例。

## 验证与执行记录

- 校验基模选择唯一、内部/用户模式映射、与合并表 frequency/Q 一致。
- 校验 band 数据列、合并表中的全部有限 f-Q 点、线性坐标、输出 PNG 可解码。
- 目检标题、坐标、图例、虚线、标签和散点不裁切、不重叠。

### 执行结果

- 已生成 `overview/p_y_fundamental_band_f_Q.png`，尺寸 3600×1941；原有 band 图和
  `overview/f_Q.png` 均未覆盖。
- 基模从现有 selection 表唯一解析为用户 `p_y`：mode 19、197.9530009044521 THz、
  Q=7676.733033102247；与合并本征频率表逐项一致。
- 右栏绘制 216 个有限 frequency/Q 点，Q 范围 0–8000；三栏 frequency 范围均为
  190–210 THz，所有坐标均为 linear。
- PNG 已实际目检：标题、band 图例、Gamma/路径标签、红色基模标签、水平虚线和散点均未裁切或
  冲突；最终图不显示内部 symmetry/PEC/PMC。
- 按参考图复核 band/Q 面板高宽比和三栏宽度；横轴最终显示
  `0.3MΓ / Γ / 0.3ΓK`，没有第二排路径文字，x/y 文字无重叠。
- 三栏间距进一步压缩；基模标注移至 Q 面板右上角，参考虚线连续跨越三栏。
- 按用户复核要求再次整体放大字体：标题 16 pt，坐标轴/刻度/基模说明 14 pt，模式图例 12 pt；
  模式图例移至整图右上角并避开 `Quality factor` 标题和基模说明。
- Q 面板移除底部重复的 `Q` 标签；散点改为深色 58 pt² 底圈、43 pt² 实心蓝色主体和
  2 pt² 偏置高光，消除白心气泡感；
  红色基模点采用同样的三层结构。频率轴只标 `190/195/200/205/210`，Q 轴只标
  `0/4000/8000`；面板上缘下移，为右上角模式图例增加独立留白。
- `py_compile`、数据身份断言、PNG 解码和本任务 `git diff --check` 通过；未启动 COMSOL。
