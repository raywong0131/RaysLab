# Blueprints 外围刻蚀开窗与参数标记执行计划

日期：2026-08-27
状态：工程实现与验证完成；尺寸待用户按预览审批

## 1. 目标

在现有 `blueprints/` 独立版图流水线中补充两类同属空气刻蚀的二维结构：

1. 围绕 finite 六边形结构六条边布置的 6 个大面积分段刻蚀开窗；
2. 位于完整结构下方、由 cavity/cladding 紧凑参数自动生成的老式计算器七段式参数标记。

最终输出合同不变：每个成功 case 仍严格只包含一个 DXF、一个综合 `preview.png` 和一个
`layout_report.json`。本任务不启动 COMSOL，不修改既有 `.out` case。

## 2. 工程默认值与显式参数

首版根据用户提供的 SEM 参考图采用以下可调工程默认值；这些数值只用于形态审批，不视为
已完成加工标定：

- 外围开窗引用 `finite_boundary`；内侧离结构边界 `2.5 um`；径向宽 `8.0 um`；
  每条边内端缩进 `3.0 um`；外端相对内端向两侧扩展 `1.5 um`；圆角 `1.0 um`；
  圆角每四分之一圆使用 8 段离散；layer 为 `ETCH`。
- 参数标记引用 `geometry.unit_cell.cavity/cladding`，不手写重复参数；格式为 b0 固定 1 位、
  eta 固定 2 位、zeta 最多 3 位并去除尾零，显示文本为
  `245.0 0.96 1.156 — 242.0 0.98 0.93`。
- 标记以外围开窗整体为竖直参考，水平居中，顶端离最下方开窗 `3.0 um`；首轮视觉检查
  后将数字调整为高 `3.0 um`、宽 `1.7 um`、段宽 `0.28 um`，段间隙 `0.10 um`，
  倒角 `0.09 um`，数字间距 `0.28 um`，小数点尺寸 `0.28 um`，参数间距
  `1.30 um`，两组与长横线间距 `1.60 um`；layer 为 `ETCH`。
- `fabrication.local_correction` 增加按 region 的显式配置：`cavity/cladding`、
  `outer_etch_window`、`parameter_marker` 可各自标定，不在源码内隐式复用加工 bias。

## 3. 数据流与模块边界

```text
parameter.json
  -> hex geometry source 构建 cavity/cladding nominal layout
  -> auxiliary geometry builder
       -> 6 segmented outer etch windows
       -> derived seven-segment parameter label
  -> 统一 dose/envelope + 按 region local correction
  -> layout transform + writer grid snap
  -> DRC
       -> 所有 feature 做合法性、grid、overlap、bridge、layer 检查
       -> finite boundary containment 只检查内部 cavity/cladding 孔
  -> DXF / preview.png / layout_report.json
```

辅助结构通过中性 `LayoutFeature` 接入。DXF writer 不增加六角或字体特例。外围开窗和参数
标记在 metadata 中明确声明为 finite-boundary 外部 feature，因此不参加原 finite footprint
containment；这不等于跳过 DRC，它们仍参加其余几何与加工规则。

## 4. 兼容与输出策略

- `outer_etch_windows.kind = none`、`markers.kind = none` 保留关闭能力；官方
  `parameter.json` 显式启用两项。
- 新字段进入 layout identity，旧 case 不覆盖；生成新的工程 case。
- 七段字体仅支持本合同需要的 `0-9`、小数点和长横线，字符均转为闭合 polygon，不在
  DXF 中写 `TEXT/MTEXT`。
- preview 的中心/边缘 nominal-corrected 对比只从允许比较的 cavity/cladding feature 中
  选择；全局面板显示全部空气开口。

## 5. 验证

- 配置严格性：缺字段、非法尺寸、未知 kind/字段拒绝；显示文本精确匹配目标字符串。
- 几何单测：6 个开窗、开窗相互不交叠、七段 glyph feature 合法、label 水平居中并位于
  开窗下方。
- DRC：外部 feature 不触发 finite-boundary containment，但仍能发现 overlap/grid/layer
  问题。
- CLI：小配置仍维持三文件输出；完整配置 preflight 通过并产生新的 feature/region 统计。
- DXF 独立回读：只包含闭合 `LWPOLYLINE`、微米单位和 `ETCH` layer。
- 实际检查新 `preview.png` 的六开窗形态、参数文字可读性、相对位置和无裁切。
- 运行 standalone 全套测试、COMSOL 几何 parity、相关 `py_compile` 和
  `git diff --check`；不启动 COMSOL。

## 6. 回滚

将 `outer_etch_windows` 与 `markers` 的 `kind` 改回 `none` 即可回到仅纳米孔版图；代码
模块和旧 case 均无需删除。若首版比例未获审批，只调整 `parameter.json` 中上述显式尺寸并
生成新 identity，不覆盖旧样例。

## 7. 执行记录

### 实际改动

- 新增 `blueprint_core/auxiliary_geometry.py`：从 finite 六边形 boundary 生成 6 个显式参数化
  梯形圆角开窗；从 cavity/cladding 紧凑参数自动格式化目标字符串，并将 0--9、小数点和
  长横线展开为闭合七段式 polygon。
- `parameter.json` 新增并启用 `outer_etch_windows` 与七段式 `markers`，所有位置、尺寸、
  圆角离散、格式精度、间距和 layer 均显式；`source_trace` 已同步当前
  `scripts/parameter.json` SHA-256。
- `fabrication.local_correction` 扩展为 `default + by_region` 严格合同；外围开窗和参数
  标记拥有独立零 bias 占位，后续可分别标定。
- CLI 在 hex nominal geometry 之后、统一修正之前组合辅助结构，并把 kind、显示文本和
  feature count 写入 layout report；辅助字段进入 artifact identity。
- DRC 新增通用 per-feature `boundary_check` 范围：外部开窗/标记不参加 finite footprint
  containment，但仍参加 polygon、grid、overlap、bridge、keep-out 和 layer 检查。
- preview 增加外围开窗/参数标记配色，并从 nominal/corrected 局部对比候选中排除大开窗
  和文字段，中心/边缘面板仍比较 photonic-crystal 小孔。
- README 已补充参数入口、关闭方式、按 region 修正和 DXF 无 `TEXT/MTEXT` 合同。

### 测试与样例

- 新增辅助几何测试，并扩充 config、DRC、CLI fixtures；standalone 全套加 COMSOL 纯
  Python geometry parity 共 `65 passed`。
- 完整 20/20 preflight：29,648 个闭合 feature，其中 cavity 7,566、cladding 21,960、
  outer windows 6、parameter-marker segments 116；DRC 0 error、0 warning。
- 最终 layout bounds 为 `[-41.63, -45.734, 41.63, 39.734] um`；最小孔宽
  `0.1928759266 um`、最小 polygon edge `0.1664331698 um`、最小桥宽
  `0.1155427295 um`、内部 finite-boundary clearance `0.471 um`。
- ezdxf 独立回读：R2010/AC1024、micrometers、29,648 个全部闭合 `LWPOLYLINE`、
  89,478 vertices、唯一 `ETCH` layer；DXF polygon 总面积 `2180.3678360002145 um^2`。
- 报告中的显示字符串精确为
  `245.0 0.96 1.156 — 242.0 0.98 0.93`，source trace 匹配当前共享参数文件。
- 实际 PNG 已在完整分辨率输出基础上完成视觉 QA：六开窗对称、与结构不接触、参数标记
  水平居中且无裁切；首轮偏细字体经检查后已加粗，偏细中间 case 已删除。
- 最终 case：
  `blueprints/.out/blueprint_hex_engineering_preview_cav20-clad20_shiftx0-y0_engineering-unassigned-r0_8076bc19/`，
  严格只含 `engineering_preview.dxf`、`preview.png` 和 `layout_report.json`。
- 相关 `py_compile`、standalone dependency scan 和 `git diff --check` 通过；未启动 COMSOL，
  未覆盖或删除任何既有正式结果。

### 遗留问题与最终状态

当前状态为 **engineering preview ready**。外围开窗和字体尺寸仍需要用户根据新预览确认；
production release 仍受 fabrication（包括三个 region 的独立 bias）、dose、polarity 与加工
DRC 标定门槛约束。若只调整形态，修改 `parameter.json` 会生成新 identity，不覆盖本次
样例。
