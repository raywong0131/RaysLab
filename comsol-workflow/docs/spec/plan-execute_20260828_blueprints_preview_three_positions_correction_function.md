# Blueprints 三位置与修正函数预览执行计划

## 方案概述

重构 `preview.png` 的右侧信息区：删除 `Fabrication summary` 文本，以三个沿同一径向 cell strip 的代表三角形和一幅实际空间修正函数图替代。三个位置固定为中心、cavity 最外层和 cladding 最外层，并在下方径向带使用相同位置标记建立一一对应。

## 范围

- 修改 `blueprints/workflow/output/preview.py` 的代表 feature 选择、复合布局和修正曲线绘制。
- 修改 `blueprints/workflow/cli.py`，向预览上下文提供完整 `spatial_envelope` 配置。
- 修改严格 export schema 和 `main/parameter.json`，增加显式 `preview_layout` 版本，确保新预览使用非冲突 identity。
- 更新 `blueprints/tests/test_dxf.py`、配置/CLI identity 测试和 `blueprints/README.md`。
- 生成一个新 case；不覆盖或删除 `083368ca` 及更早 case。

## 数据流

```text
corrected LayoutPlan + nominal LayoutPlan
  -> 选择完整 shell 0..max_shell 的单 cell 宽径向带
  -> 同一 local_index 选取 shell 0 / cavity_layers / max_shell 三个三角形
  -> 三幅共享 local 坐标尺度的 nominal/corrected 对比

fabrication.spatial_envelope
  -> 沿同一径向带连续采样 radial_even_polynomial
  -> edge offset 转换为 nm
  -> 标记 shell 0 / cavity_layers / max_shell

三位置 shell
  -> 在下方径向带用同色位置标记
```

## 最小绘图合同

- 目标文件：新 case 内的 `preview.png`；仍只生成 DXF、PNG、report 三个文件。
- 面板顺序：左侧 `Corrected layout`；中列从上到下为
  `Center (shell 0) — nominal / corrected`、
  `Cavity boundary (shell N) — nominal / corrected`、
  `Cladding boundary (shell M) — nominal / corrected`；右列为
  `Applied correction along radial strip`；底部为现有 radial cell strip。
- 数据源：三个三角形严格来自下方所选径向 parent-cell ray，使用同一 `local_index`；当前正式参数中 N=20、M=40。
- 三角形对比：nominal 灰、corrected 蓝；转换为各自 parent cell 的 local 坐标；三个面板共享同一 x/y 范围和 equal aspect，不对每个面板独立缩放。
- 修正函数：按配置中的 `radial_even_polynomial` 在径向 parent-cell 中心线上连续采样；横轴为 parent shell，纵轴为 `Triangle side-length correction (nm)`。数值使用 `2√3 × (region-local edge offset + spatial-envelope edge offset)`，因此同时反映局部统一尺寸修正和空间包络修正；零线为克制的浅灰实线。
- 对应标记：中心、cavity boundary、cladding boundary 在曲线和下方径向带使用相同颜色与文本；保留 cavity/cladding 几何分界，不添加 colorbar。
- 不修改左侧整体图、下方孔形配色、参数标记、DXF 几何或 correction 数值。
- QA 样例：当前 20–20、圆角 3.0 µm、SEM radial correction 的正式参数。

## 参数与兼容策略

- `export.preview_layout` 使用严格枚举值 `radial_three_position_correction_v2`。
- 该字段参与配置 identity，使绘图布局修订生成新 case；DXF 几何参数和 fabrication revision 不作伪修改。
- `save_preview` 保留现有 positional API；预览上下文新增 `spatial_envelope` 字段。无法获得完整径向带或修正配置时，相应面板显示明确 unavailable 信息，不猜测数据。

## 验证

- 单元测试锁定三位置 shell 与同一 local index 选择。
- 锁定三个精确标题、共享坐标范围、底部三个位置 marker、修正曲线标题/坐标标签/采样数/标记数，并确认不存在 `Fabrication summary`。
- 运行聚焦预览测试、blueprints 完整 pytest、相关 `py_compile`、CLI preflight 和 `git diff --check`。
- 生成新 case，确认 DRC 0 error / 0 warning、严格三文件合同，并按实际 PNG 原分辨率检查面板排列、文字碰撞、三位置对应和修正曲线。
- 不启动 COMSOL。

## 回滚方案

恢复上一版 preview 生成器和 `export.preview_layout` 参数即可回到两位置加 summary 的布局；旧 `083368ca` case 保留，可作为未覆盖的视觉基准。

## 执行记录

### 紧凑布局修订

用户在首版 v2 PNG 基础上要求：放大中间三幅代表三角孔、压缩整体留白，并从可见
文字中移除 `nominal/corrected`。本轮保持三处 shell、共享 local 坐标尺度、灰/蓝
叠加、修正函数、下方径向带、DXF 几何和工艺数值不变；仅调整 Matplotlib 画布比例、
GridSpec 宽高比、constrained-layout padding、内层行距与标题。

精确标题修订为 `Center (shell 0)`、`Cavity boundary (shell 20)`、
`Cladding boundary (shell 40)` 和 `Radial cell strip`。为确保布局变更不覆盖
v2 case，`export.preview_layout` 升级为
`radial_three_position_correction_compact_v3`；配置校验继续接受历史 v2 值。

- 已重构 `blueprints/workflow/output/preview.py`：从同一条完整径向 parent-cell ray
  解析 shell 0、cavity 边界和 cladding 边界三处三角孔；三个 nominal/corrected
  面板使用共同 local 坐标范围，并在下方径向带加入同色 shell 标记。
- 已删除旧 `Fabrication summary` 信息区，改为连续采样真实
  `radial_even_polynomial` 的 `Edge offset (nm)` 曲线；shell 0/20/40 的曲线
  marker 与三角孔面板、下方径向带一一对应。
- `blueprints/workflow/cli.py` 已把完整 `fabrication.spatial_envelope` 传入预览；
  严格 export schema 和 `main/parameter.json` 已加入
  `preview_layout=radial_three_position_correction_v2`，且该字段参与 identity。
- 新 identity 为 `ba679d63`；新 case 为
  `dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_ba679d63`。
  旧 `083368ca` case 经生成前检查仍存在，未覆盖或删除。
- 聚焦测试：36 passed；`blueprints/tests` 全量测试：84 passed；相关三个源码文件
  通过 `py_compile`。
- 正式参数 preflight 和生成均得到 29,648 个 DXF entity，DRC 为 0 error /
  0 warning。新 case 严格只含 DXF、`layout_report.json` 和 `preview.png`
  三个文件。
- 已检查 3520×1650 正式 PNG 的整体图及原像素右侧/底部裁剪：三个标题不碰撞，
  local 尺度一致，修正曲线完整，shell 0/20/40 标签和 cavity/cladding 分界清晰。
- 未启动 COMSOL；未修改 DXF 几何、SEM 修正数值、dose 或 release 状态。

最终状态：本计划范围已完成，无代码遗留问题。当前配置仍是
`engineering_preview`，工艺标定限制维持不变。

### 紧凑布局修订结果

- 已将画布从 16.0×7.5 inch 调整为 14.8×8.2 inch，并缩小 constrained-layout
  padding、三角孔面板行距和底部径向带占比；中间列与修正函数列等宽，使三幅
  local 三角孔相对 v2 明显增大，正式 PNG 从 3520×1650 调整为 3256×1804。
- 中间三幅可见标题已改为 `Center (shell 0)`、
  `Cavity boundary (shell 20)`、`Cladding boundary (shell 40)`；
  底部标题改为 `Radial cell strip`。这些位置不再显示
  `nominal/corrected`，灰/蓝叠加本身保持不变。
- compact v3 identity 为 `8b2743f1`；正式 case 严格生成 DXF、
  `layout_report.json`、`preview.png` 三个文件。此前
  `083368ca` 和 `ba679d63` case 均保留。
- compact v3 聚焦测试 37 passed，`blueprints/tests` 全量测试 85 passed；
  preflight 与正式生成均为 29,648 features，DRC 0 error / 0 warning。
- 已检查正式 PNG 整体图：中间三图尺寸、标题间距、三处共享坐标尺度、修正曲线和
  底部 shell 0/20/40 对应均正确，无标题碰撞或裁切。未启动 COMSOL。

### 两侧面板加宽修订

用户确认 compact v3 后先要求右侧修正函数图更宽，随后补充要求左侧整体预览也放大。
本轮只调整横向布局：画布宽度由 14.8 inch 增至 16.5 inch，GridSpec 列权重由
`(1.85, 1.35, 1.35)` 改为 `(2.10, 1.35, 1.65)`。相对 compact v3，
左侧整体图实际宽度增加约 13%，中间三角孔面板近似保持，右侧曲线宽度增加约 22%；
高度、纵轴范围、曲线数据、标题、marker、底部径向带及 DXF 均不变。

为避免覆盖 compact v3，`export.preview_layout` 升级为
`radial_three_position_correction_compact_wide_v4`，配置校验保留 v2/v3
兼容值。完成后需运行聚焦/全量测试、preflight、严格三文件生成、DRC/DXF 读回，并
检查正式 PNG 的左右面板宽度与标签完整性。

wide v4 已完成：identity 为 `ad99cbcd`，正式 PNG 为 3630×1804。视觉检查确认
左侧整体版图相对 compact v3 放大、右侧修正函数横向空间增加、中间三图尺寸保持，
标题、三处 marker 和底部径向带均无碰撞或裁切。聚焦测试 37 passed，
`blueprints/tests` 全量测试 85 passed；preflight 与正式生成均为 29,648
features，DRC 0 error / 0 warning，严格三文件合同保持。compact v3 及更早 case
均保留，未启动 COMSOL。

### 压缩中列空白并继续加宽右图

用户检查后指出右图仍偏窄，且中间三角孔列两侧空白过多。原因是中间 GridSpec 列宽
明显大于由行高限制的方形 local 坐标轴；继续增大总画布不能有效解决该空白。

本轮将画布调整为 15.7×8.2 inch，列权重从
`(2.10, 1.35, 1.65)` 改为 `(2.10, 0.85, 1.90)`。按实际横向分配计算，
左列约 6.80 inch、与 v4 基本相同；中列容器由约 4.37 inch 收至约 2.75 inch，
但方形三角孔坐标轴仍由行高限制，因此实际图形尺寸近似不变；右列由约 5.34 inch
增至约 6.15 inch，增加约 15%。释放中列无效空白而不是单纯扩宽总画布。

本次按用户最新要求视为 wide v4 的预览原位细调，不再升级
`export.preview_layout`，也不创建新 case 或重写 DXF。修改共享 preview 生成器后，
在内存中用当前正式参数重建 nominal/corrected layout，只覆盖现有 `ad99cbcd`
case 的 `preview.png`，并同步更新 `layout_report.json` 中的 preview SHA；
DXF 文件和 DXF SHA 必须保持不变。完成后运行聚焦/全量测试、DRC/DXF 读回和正式 PNG
视觉检查，确认中列空白缩小、右图变宽且中间三图没有缩小。

原位细调已完成：正式 PNG 为 3454×1804，中列容器空白明显缩小，右侧修正函数
横向空间进一步增加，中间三幅方形图与左侧整体图的实际尺寸保持，标题和标记无碰撞
或裁切。只覆盖 `ad99cbcd/preview.png` 并同步
`layout_report.json.outputs.sha256["preview.png"]`；preview SHA 从
`1c9a2053…b83bd` 更新为 `6431c386…57755`。DXF 未重新生成，其 SHA
`5add8007…e41d` 前后完全一致，case 仍严格只有三个文件。聚焦测试 37 passed，
`blueprints/tests` 全量测试 85 passed，DRC 保持 0 error / 0 warning，未启动
COMSOL。

## 用户复核修订：右侧曲线纳入整体三角孔边长修正

用户检查 `93602506` 后指出，右侧曲线只显示空间包络 edge offset，遗漏了已实施的
三角孔边长统一减小 `10 nm`，导致中心看起来接近 `0`。共享 preview 现改为显示实际
三角形边长变化：先把 cavity/cladding region 的 local edge offset 与连续空间包络
edge offset 相加，再乘以等边三角形换算因子 `2√3`。纵轴同步改为
`Triangle side-length change (nm)`。

当前参数的三个标记值为：center `-10.341325 nm`、cavity boundary
`-12.420242 nm`、cladding boundary `-96.963687 nm`。center 与精确 `-10 nm` 的
`-0.341325 nm` 差异来自原有 radial polynomial 在中心的非零拟合常数，而不是整体
修正遗漏。

只原位重绘 `93602506/preview.png` 并更新报告中的 PNG SHA；DXF SHA 保持
`6dd9578ad7ca276d4d7b1790d002c98106c8d45abf008c287465fc0ea8d3dcc0`。新 preview SHA
为 `08f1b566316fb462148b5942169b0226847d5ac96abc04990ba96ec1047bb05c`，更新后 report
SHA 为 `498b8c75a00e5df55f59885028390058446155c0b1cb213a99a29c185c3786c9`。聚焦测试
`27 passed`，完整 blueprints 测试 `87 passed`；实际 PNG 视觉检查确认中心 marker
位于约 `-10 nm`、三处标记和曲线未裁切。未重写 DXF，未启动 COMSOL。

### 用户术语确认：显示边长修正量

用户进一步明确右图表达的是“边长修正”，不使用“总变化”表述。SEM 修正表中的
`recommended_side_compensation_nm` 已确认最外层包络边长修正为
`-86.963687 nm`；叠加独立的统一边长修正 `-10 nm` 后，右图最外层显示
`-96.963687 nm`，数值无需改变。标题改为
`Applied side-length correction along radial strip`，纵轴改为
`Triangle side-length correction (nm)`。

只再次原位更新 `93602506/preview.png` 及报告中的 PNG SHA。最新 preview SHA 为
`e98ce8406edd2c609ff9c896c84a1cd82ca7a7ffb67a8ad8748a7ec5d9131053`，report SHA 为
`4047d333400369b99be6280791bfecfc65d70ac54d088a23d0335cffa9003c4e`；DXF SHA 仍为
`6dd9578ad7ca276d4d7b1790d002c98106c8d45abf008c287465fc0ea8d3dcc0`。preview 聚焦
测试 `7 passed`，完整 blueprints 测试 `87 passed`，实际 PNG 已检查且未裁切。
