# 科研绘图与复现

本文件是科研绘图样式及输出验证的规则所有者；长期协作偏好见 [Memory](../MEMORY.md)。当前用户的最新明确要求优先；
若新要求改变长期习惯，应同步更新 Memory。某张既有图的精确文件名、标题、数值范围、
色标端点和布局以当前代码、测试与本次请求为准，不把历史任务参数复制成新图默认值。

- 用户说“试一下”“以这个数据为例”或明确要求审批后再改代码时，视为样例预览阶段：
  只使用用户指定的现有数据，在原结果目录内重绘指定文件，不得新建
  `approval_preview`、`preview_v2` 等目录，不得提前修改共享生成器。
- 在上述预览阶段，用户明确说“批准修改”“确定为当前形式”后，把样例形式固化到共享绘图模块和
  保留的兼容入口，并同步聚焦测试与输出文档。不得让样例脚本和正式生成器形成两套
  实现。
- 用户直接要求修改生成器或已明确授权预览与实现两个阶段时，完成实现和验证，不额外设置样例审批。
- 修改多面板图前先列出最小绘图合同（已确定内容无需重问）：目标文件、面板顺序、精确标题与大小写、数据源、
  归一化、colormap、坐标范围与 aspect、colorbar 范围/刻度/标签/方向、叠加线和输出
  路径。字符串按用户给定形式处理，不擅自统一 `Intensity/intensity` 等大小写。
- 输出按格式分类：overview 目录只放 PNG；PDF 集中到同级编号 pdf 文件夹，任务 Markdown 报告放同级 `12_reports`；工程 Spec 放仓库 `docs/spec`。
  修改保存路径时同步活动清单和报告中的相对链接；历史结果仅在用户授权范围内迁移。
- 严格限制重绘范围。用户只要求 A 图时，不得顺带覆盖 B–D；只要求一个 mode 时，
  不得批量改写其他 mode。已有 CSV、JSON、NPZ、parquet 足够时不得重跑 COMSOL。
- 可比较的面板必须使用明确且科学合理的 scale。需要共享归一化时，应在绘图前对完整
  比较集合一次完成；不得把每个面板分别归一化后伪装成可比较。具体坐标范围、clim、
  归一化目标、刻度端点、marker/font 大小、精度和箭头长度每次结合数据与物理含义判断，
  不沿用上一张图的数值。当前请求明确给定数值时，才为该输出增加对应数值测试。
- 保持 x-y、kx-ky 图等比例，不拉伸数据。colorbar 与对应线框等高；跨两行的共享
  colorbar 必须从上排线框上缘准确延伸到下排线框下缘。当前反向竖排说明文字使用
  已批准的反向阅读方向，具体 rotation 参数按当前布局实现。log 色标优先按 `10^n`
  显示，不写 `lg`、`log10` 或裸指数；指数范围由当前数据决定。
- 审美默认是克制、清晰、可发表的 Nature-style 科研图：高值显眼、低值安静，相位使用
  循环色图；不添加无信息量的边框、虚线、逐点百分数、中心点加粗或重复图例。用户已
  指定 colormap 时必须使用指定值，不得以“更符合 Nature”为由自行替换。
- 频率、Q、NA、P(k) 等数值必须从当前结果 metadata 读取，不得硬编码样例值。显示
  编号与内部 FFT/cyclic 索引分离；改变编号标签不得改变数组顺序、Fourier 权重或重建。
- 导出后必须用实际 PNG 做视觉检查，确认标题不碰撞、线框和 colorbar 对齐、刻度完整、
  箭头/NA/BZ 不遮挡数据。修改共享生成器时，用聚焦合同测试覆盖本次涉及的标题、colormap、轴范围、clim、
  colorbar ticks/labels、overlay 与文件名；未变项目复用已有覆盖。随后运行相关 `py_compile`
  和 `git diff --check`。仅用现有函数重绘单张图时，核对数据与实际 PNG，不为临时预览新建测试套件。
- 交付时说明原位更新了哪些文件、删除或替代了哪些旧文件、是否创建目录、是否启动
  COMSOL，以及测试结果。仅创建本次需求所需目录，不覆盖无关正式结果。

## 通用新图样式与实现位置

- 新绘图脚本先在 `tests/<topic>/plot_<capability>.py` 打通；可复用绘图算法先在 tmp 试验，
  按 [晋升通道](../AGENTS.md#新算法与新脚本的两条晋升通道) 接入正式所有者。
- 字体统一为无衬线，优先 Arial，缺失用 DejaVu Sans；同一图/批次用同一字体族。
  新 PNG 默认不低于 300 dpi，PDF 保留矢量文字/线条；既有批准格式/DPI 不因文档更新重绘。
  同一图各格式同 stem；跨图物理面板尺寸、字号层级、线宽和符号含义在生成前确定并保持一致。
- 图保存使用已有 field_plotting 或 figure_output；位置由 output_paths 分类。
  预览/合成检查输出进试验区 .work，不写 archive；用户要求现有图原位预览时遵守该明确范围。
- 视觉样式可统一，phase、单位、数据范围、NA/BZ、归一化对象和颜色的科学含义须单独核验，
  不能为版面好看而改变数据、遮掉边界外强度或逐面板旋转相位。

## 其他既有显示约定

- “Nature style” means a clean scientific figure with restrained scaffolding,
  equal geometric aspect, deliberate whitespace, concise text, and no redundant
  decoration. Prefer direct, exact scientific labels over explanatory prose in
  the panel. Do not add legends, percentage labels, dashed guides, borders, or
  emphasis marks unless they carry requested information.

- All visible axis tick marks default to pointing inward. Apply this to both
  major and minor ticks on every axis; change the direction only when the user
  explicitly requests a figure-specific exception.

- For sequential intensity maps, the maximum should be immediately visible and
  the minimum visually quiet. Use a cyclic map for phase. An explicitly approved
  colormap always overrides a general palette preference; do not silently replace
  it with another “Nature-like” palette.

- Comparable panels share a defensible normalization when their values are
  meant to be compared. Never normalize each comparable panel independently
  without saying so. Choose the normalization from the current data, physical
  question, and explicit request rather than copying a previous figure's
  numerical convention.

- Keep colorbars aligned to the plotted frame. A single-panel colorbar must not
  extend past its axes; a two-row shared colorbar spans exactly from the upper
  row frame edge to the lower row frame edge. Show required endpoints explicitly
  when they are meaningful. Use the approved reversed reading direction for
  vertical colorbar labels without treating a past rotation value as a general
  numerical default.

- Display logarithmic intensity ticks as powers of ten (`10^n`), not as `lg`,
  `log10`, or raw exponents. Numerical transforms may use `log10` internally,
  but the reader-facing label follows the approved power notation and the
  exponent range comes from the current data.

- Populate frequency, Q, NA, fractions, and other numbers from the current
  result metadata. Never hard-code values copied from an example. Preserve exact
  requested capitalization and symbols, including `Intensity` versus
  `intensity`, `P(k)`, `P(Γ)`, and `H_z`.

- Do not promote past task-specific numbers into general plotting preferences.
  Axis ranges, color limits, normalization targets, tick endpoints, marker and
  font sizes, decimal precision, shell/layer counts, arrow lengths, and inset
  coordinates must be chosen afresh from the current data and physical meaning.
  Exact values for an existing approved output belong in its code, tests, and
  result metadata, not in durable memory.

- Useful non-numerical patterns from the approved finite figures may guide new
  plots: keep comparable panels visually aligned; use concise data-derived
  titles instead of redundant mode identifiers; keep main and inset palettes
  visually coherent without forcing them to be identical; place insets where
  they minimally obscure the data and avoid a conspicuous border; give equal
  semantic items equal marker emphasis; distinguish structural boundaries from
  secondary guides; keep reference overlays subdued; and place measurement
  arrows outside curves, pointing toward the measured feature without crossing
  the data.

- Label every plot axis with its quantity and unit. Use equal physical aspect
  for x-y and kx-ky images; never stretch field data to fit a panel.

- Keep side-by-side cavity/cladding comparisons aligned and use shared frequency
  or Q scales where comparison is intended.

- Center the shared frequency window of the unit-cell band comparison on the
  cavity `p2@Gamma` frequency so the two cells remain directly comparable.

- In every user-facing band plot, use red for the user's `py` and blue for the
  user's `px`, with a continuous red-white-blue transition for mixtures. This
  means project-internal `px` is rendered red and project-internal `py` blue.
  Display `d_{xy}` in light gray and `d_{x^2+y^2}` in black, with a continuous
  gray-black transition. Band-comparison curves use lines without point markers.

- Far-field k-space intensity plots use total `|E|^2`, preserve data outside
  drawn Brillouin-zone boundaries, and mark the relevant BZ and numerical
  aperture without distorting the axes. Gaussian-fit figures report the fitted
  k-space width and corresponding full-angle divergence without obscuring the
  data. Select linear/log display ranges from the actual dynamic range and the
  comparison being made; a range approved for an earlier figure is not a
  default for a new plot.

## Hz 场图标准

用户已指定 `results/finite_cavity_reference_results/shaolei_hex_side2_forward/singularity_original_r19_c20/selected_mode_14_plots/Hz_Im_2d.png` 的 colorbar 为后续默认标准：

- Hz 实部/虚部二维场图使用 `RdBu_r`：负值蓝、零值白、正值红；固定 clim 为 `[-1, 1]`，刻度为 `-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00`。
- 单独场图用当前分量有限样本的最大绝对值归一化；零场保持零。分量标为 `Re(Hz)` / `Im(Hz)`，色标明确写 `normalized`、无量纲（`1`），不再标 `A/m`。这种图用于比较形状，不能据此比较独立图的绝对幅值；需要幅值比较时另用完整比较集合的共同 scale。
- 共享 `field_plotting.py` 负责 PNG、SVG、hybrid SVG 与 PNG 嵌入 HTML 的一致性。控制器及旧版 spatial 的 Hz 二维图片导出也复用它；原始 TXT/parquet/复场及求解配置保持物理量不变。独立后处理包使用相同的 Matplotlib 色带与数值合同。
- Hz 强度、模长、相位和 Fourier 权重具有不同物理含义，不套用有符号场的 `[-1, 1]` 归一化。既有结果不会因修改默认生成器而自动重绘。
