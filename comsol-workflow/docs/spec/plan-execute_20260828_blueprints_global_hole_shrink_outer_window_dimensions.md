# Blueprints 三角孔整体缩小与外窗尺寸固定执行计划

## 方案概述

在现有 `sem-air-hole-220nm-radial-v1` 空间包络修正基础上，为 cavity/cladding
三角孔叠加统一的边长缩小 `10 nm`；同时将六个外层梯形开窗改为显式、彼此独立的
`13.0 um` 高度、`3.8 um` 结构净距和 `5.0 um` 相邻平行斜边净距。生成新的稳定
identity 和非冲突工程预览 case，不覆盖既有 DXF，不启动 COMSOL。

## 范围

- 修改 `blueprints/main/parameter.json` 中的局部修正和外窗正式参数。
- 修改 `blueprints/workflow/config.py`，把外窗平行斜边净距纳入严格 schema。
- 修改 `blueprints/workflow/geometry/auxiliary.py`，独立构造给定净距的梯形侧边。
- 更新 `blueprints/tests/` 的配置、几何和 identity 合同。
- 更新 `blueprints/README.md` 的尺寸及修正语义。
- 运行独立测试、预检和实际生成，回读 DXF 并检查实际 PNG。

## 修正换算与数据流

等边三角孔由三条边做等距平行偏置。三角形边长变化与单边偏置的关系为
`delta_side = 2*sqrt(3)*edge_offset`，因此整体边长减小 `10 nm` 对应：

```text
edge_offset_nm = -10 / (2*sqrt(3))
               = -2.886751345948129 nm
```

该值写入 `fabrication.local_correction.default`，并与现有
`fabrication.spatial_envelope` 的位置相关 edge offset 相加。`outer_etch_window` 与
`parameter_marker` 继续使用各自的零偏置覆盖，因此不会被整体缩小。

```text
nominal cavity/cladding triangle
  -> local edge offset -2.886751345948129 nm
  +  radial SEM envelope edge offset
  -> grid snap -> DRC -> DXF / preview / report
```

## 外窗几何合同

正式参数为：

- `radial_width_um = 13.0`：沿对应六边形边外法线测得的梯形高；
- `inner_clearance_um = 3.8`：梯形近结构底边与 finite 边界的法向净距；
- `parallel_side_gap_um = 5.0`：相邻梯形两条互相平行的径向斜边之间的法向净距；
- `inner_corner_radius_um = 3.0`：只作用于近结构底边的两个角。

对于每个共享六边形顶点，先由 `inner_clearance_um` 得到相邻两条径向边的基础间距，
再沿两侧底边切向对称收短，使径向平行线的净距达到
`parallel_side_gap_um`。本组规则六个顶点等价，每个底边端部收短约 `0.693 um`。
用户已明确允许放宽零间距拼环能力，因此新字段要求为正，不保留
`parallel_side_gap_um=0` 的完整六边形环特殊合同。

## 预览绘图合同

- 目标：新 case 内唯一的 `preview.png`。
- 面板、标题、配色和现有 `radial_three_position_correction_compact_wide_v4` 布局不改。
- 左侧展示新完整布局；右侧中心/cavity 边界/cladding 边界三角孔与下方径向带保持对应；
  修正函数继续展示局部常量缩小与空间包络叠加后的结果。
- 使用同次生成的 nominal/corrected polygon；x-y 保持等比例。
- 以实际 PNG 原像素检查外窗、圆角、七段标记、三角孔和裁切。

## 验证与输出

- 配置：新字段必填、正数、拒绝未知字段；正式参数精确为 13/3.8/5 um。
- 几何：6 个合法且互不相交的开窗；相邻斜边平行；逐对法向距离为 5 um；高度为
  13 um；近结构底边净距为 3.8 um；只圆滑两个内角。
- 修正：中心和边界样本的 effective edge offset 等于统一负偏置与包络值之和；外围开窗
  和标记不受 10 nm 边长缩小影响。
- 运行聚焦测试、完整 `blueprints/tests`、相关 `py_compile`、CLI preflight、实际生成、
  DXF readback、`git diff --check` 和 PNG 视觉检查。
- 新 case 严格只包含 DXF、`preview.png`、`layout_report.json`，不删除旧 case。

## 回滚方案

将 `local_correction.default` 恢复为零偏置，并把外窗 13/3.8/5 um 参数恢复为上一版值，
即可恢复旧设计输入；旧 identity case 始终保留。代码回滚仅涉及本计划列出的
geometry/config/tests/README 文件和本计划文档。

## 执行记录

### 实际改动

- `blueprints/main/parameter.json` 将 cavity/cladding 默认 edge offset 改为
  `-2.886751345948129 nm`；外窗和参数标记的 region override 保持 `0 nm`。
- 外窗正式尺寸改为 `radial_width_um=13.0`、`inner_clearance_um=3.8`、
  `parallel_side_gap_um=5.0`，内角圆角继续使用 `3.0 um`。
- `blueprints/workflow/config.py` 将 `parallel_side_gap_um` 加入 v2 必填严格字段并要求
  正数，未知或零/负值均拒绝。
- `blueprints/workflow/geometry/auxiliary.py` 根据实际 finite 六边形逐顶点计算结构净距
  导致的基础间距，再沿相邻底边切向对称收短，使两条共享径向方向的平行斜边达到目标
  法向净距；metadata 记录高、结构净距、目标间距、基础间距及两端收短量。
- 按用户补充要求取消零间距拼环合同；当前只支持不小于 clearance-induced gap 的正
  `parallel_side_gap_um`，避免负向延长和相邻开窗接触。
- 更新 config/geometry/CLI identity 测试与 README；既有 case 未覆盖或删除。

### 数值与测试验证

- 当前实际 finite 边界的六组斜边间距范围：
  `4.999999999999998–5.0 um`。
- 外窗法向高度：`13.0 um`；近结构底边法向净距：`3.8 um`。
- 六边形等价顶点的单端切向收短范围：
  `0.692820323027551–0.6928203230275516 um`。
- `2*sqrt(3)*(-2.886751345948129 nm) = -10.0 nm`，局部统一缩小与 SEM 包络由既有
  correction pipeline 相加；辅助 region 的 effective local offset 仍为零。
- 聚焦 config/geometry/corrections：`43 passed`。
- 完整独立测试：`87 passed in 7.10s`。
- 相关 `py_compile`：通过。
- CLI preflight：29,648 features；DRC `passed=True`，0 error / 0 warning；未写文件。
- DXF 回读：29,648 个 `LWPOLYLINE`，全部闭合，唯一 layer 为 `ETCH`；版本
  `AC1024`/R2010，`INSUNITS=13`（micrometers）。
- 未启动 COMSOL。

### 新生成 case

identity：
`936025065d620c79dd17b6a1979e2ba8d4978f3521f138da7da3d7ee6a51e681`

目录：
`blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_93602506/`

严格包含三个文件：

- DXF：5,622,851 bytes；SHA-256
  `6dd9578ad7ca276d4d7b1790d002c98106c8d45abf008c287465fc0ea8d3dcc0`；
- `preview.png`：1,240,510 bytes；SHA-256
  `a5bd291fbf5a022aee08b6b8d34bb15ee958b141d0903d458ddf4d1552c77d0a`；
- `layout_report.json`：13,338 bytes；SHA-256
  `0a3c104940c179d2245a5ac9ff9bbd979e09d8dd15a117c60abf8a40f5e3665b`。

实际 PNG 已按完整分辨率数据生成，并通过内存缩放视图检查：六个开窗对称，5 um
分段间隔一致；靠结构两角圆滑、外角尖锐；七段参数标记完整；中心/cavity 边界/
cladding 边界三角孔、修正函数和下方径向 cell 带均完整且未裁切。视觉检查未写入任何
额外文件。

### 遗留问题与最终状态

本次参数、几何、测试、预检、正式生成和视觉检查均完成。输出仍为
`engineering_preview`：dose 未分配、polarity 未指定，生产 DRC 阈值和 SEM
`process_gain=1.0` 二次闭环尚未标定；因此该 DXF 不应描述为生产发布版图。未 commit
或 push。

### 用户复核修订：preview 右侧显示实际边长变化

用户指出首版 `93602506` preview 的右侧曲线只画空间包络 edge offset，没有把
local correction 中的统一边长减小 `10 nm` 纳入，因此中心看起来接近零。共享绘图已
改为 `2√3 × (local edge offset + envelope edge offset)`，纵轴改为
`Triangle side-length change (nm)`。当前三个标记值依次为 center
`-10.341325 nm`、cavity boundary `-12.420242 nm`、cladding boundary
`-96.963687 nm`；中心相对 `-10 nm` 的小残差来自既有包络拟合常数。

本次只原位更新 case 内 `preview.png` 和 `layout_report.json` 的 preview SHA。DXF SHA
仍为 `6dd9578ad7ca276d4d7b1790d002c98106c8d45abf008c287465fc0ea8d3dcc0`；新 preview
SHA 为 `08f1b566316fb462148b5942169b0226847d5ac96abc04990ba96ec1047bb05c`，新 report SHA
为 `498b8c75a00e5df55f59885028390058446155c0b1cb213a99a29c185c3786c9`。聚焦测试
`27 passed`、完整测试 `87 passed`，实际 PNG 已检查；未启动 COMSOL。

### 用户术语确认：右图为边长修正而非总变化

根据用户最新表述，右图数值定义不变，但统一使用“边长修正量”语义：局部统一边长
修正与 SEM `recommended_side_compensation_nm` 相加。最外层包络修正
`-86.963687 nm` 加统一修正 `-10 nm`，得到图中的 `-96.963687 nm`。标题和纵轴分别
改为 `Applied side-length correction along radial strip` 与
`Triangle side-length correction (nm)`。

最新 preview SHA 为
`e98ce8406edd2c609ff9c896c84a1cd82ca7a7ffb67a8ad8748a7ec5d9131053`，report SHA 为
`4047d333400369b99be6280791bfecfc65d70ac54d088a23d0335cffa9003c4e`；DXF SHA 保持
`6dd9578ad7ca276d4d7b1790d002c98106c8d45abf008c287465fc0ea8d3dcc0`。聚焦测试
`7 passed`、完整测试 `87 passed`，未启动 COMSOL。

## 参数接口扩展计划：整体边长修正与包络调制强度

### 目标与参数合同

将当前由 `local_correction.default.edge_offset=-2.886751... nm` 隐式表达的三角孔整体
边长修正，提升为 `main/parameter.json` 中直接以物理量书写的显式接口：

```json
"triangle_side_correction": {
  "kind": "equilateral_parallel_edge_v1",
  "global_delta_nm": -10.0,
  "envelope_modulation_strength": 1.0,
  "regions": ["cavity", "cladding"]
}
```

- `global_delta_nm`：每个三角孔单条边的统一修正量；负值缩小、正值放大。
- `envelope_modulation_strength`：无量纲包络强度；`0` 关闭包络、`1` 使用当前 SEM
  标定幅度、大于 `1` 墪强，必须非负且有限。本版固定为 `1.0`。
- `regions`：整体修正和包络强度共同作用的 region；正式值为 cavity/cladding，并与
  `spatial_envelope.regions` 严格一致，避免两个作用域产生歧义。
- `kind` 显式声明等边三角形三边平行偏置模型，为以后接入其他孔形保留分支。

`local_correction` 继续作为通用 region 级 edge offset/scale 接口，但正式 default
恢复为零，不再承担整体三角边长修正。外窗和参数标记的零修正覆盖保持不变。

### 数据流与数值兼容

```text
global_delta_nm
  -> 除以 2√3 -> global edge offset

SEM envelope edge offset × envelope_modulation_strength
  -> modulated envelope edge offset

global edge offset + modulated envelope + generic local/dose correction
  -> polygon correction -> grid snap -> DRC -> DXF
```

当 `global_delta_nm=-10.0` 且强度为 `1.0` 时，理论修正 polygon 应与当前
`93602506` 相同；由于参数 schema/identity 改变，仍生成新的非冲突 case，不覆盖旧
输出。preview 直接读取新接口并显示
`global_delta_nm + 2√3 × strength × envelope_edge_offset`，不再从低层 local edge
offset 反推用户参数。

### 实施范围与验证

- 更新 `workflow/config.py` 的严格 schema 和 region 一致性校验。
- 更新 `workflow/fabrication/corrections.py` 的显式参数解析、包络缩放、metadata 与
  correction pipeline 报告。
- 更新 `workflow/cli.py` 和 `workflow/output/preview.py` 的 preview context 与曲线数据。
- 更新正式参数、README、config/corrections/preview/identity 测试。
- 运行聚焦测试、完整 blueprints 测试、py_compile、preflight、正式生成、DXF 回读、
  新旧 polygon/DXF 几何对照、实际 PNG 检查和 `git diff --check`。
- 新 case 仍严格只含 DXF、preview、report；不启动 COMSOL，不删除旧 case。

### 回滚

恢复原 default edge-offset operator 并移除 `triangle_side_correction` 即可回到旧输入
表达；旧 `93602506` case 保留作为数值对照。实施完成后在本节后回填新 identity、
测试、哈希和兼容对照。

### 参数接口扩展执行结果

- `main/parameter.json` 已增加 `fabrication.triangle_side_correction`：当前
  `global_delta_nm=-10.0`、`envelope_modulation_strength=1.0`，作用域为
  cavity/cladding；通用 `local_correction.default` 已恢复为 `0 nm`，避免重复修正。
- `workflow/config.py` 已将该块纳入严格 schema：强度必须有限且非负，region 非空且唯一，
  并要求与 `spatial_envelope.regions` 完全一致。
- correction pipeline 先把单边长整体修正除以 `2√3` 换算为平行边偏移，再将原始 SEM
  包络乘以调制强度；feature metadata 同时记录原始/调制包络、整体偏移和实际单边长修正。
- preview 已改为直接读取新接口，右图绘制单个三角孔边长修正：
  `global_delta_nm + 2√3 × strength × envelope edge offset`，不再由通用 local correction
  反推。
- 新 identity 为
  `a02f4d079011efcfac477c745d4aa5bd244974210a3c8dc1d855c9b269fa8ca5`；新 case 为
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_a02f4d07/`，
  严格只含 DXF、`preview.png` 和 `layout_report.json`。
- 三处曲线值为：center `-10.341325000 nm`、cavity boundary `-12.420242338 nm`、
  cladding boundary `-96.963686911 nm`；中心相对 `-10 nm` 的差值来自 SEM 拟合常数项。
- 新旧 DXF 均含 29,648 个闭合 `LWPOLYLINE`、唯一 `ETCH` layer；逐实体、逐顶点回读
  几何完全一致，证明接口重构没有改变强度为 `1.0` 时的实际版图。
- preflight：29,648 features，DRC `0 error / 0 warning`；完整独立测试：`91 passed`；
  相关 `py_compile` 与 `git diff --check` 通过。实际 PNG 已检查，整体结构、六个圆角开窗、
  参数标记、三处三角孔、修正曲线和径向 cell strip 均完整且未裁切。
- 新 DXF SHA-256 为
  `ba985afdfdbc54608b450f24b52ff1fbf3eff32f583f9d448a1d7a1dd7923f0d`，preview 为
  `e98ce8406edd2c609ff9c896c84a1cd82ca7a7ffb67a8ad8748a7ec5d9131053`，report 为
  `556b8e2526417c5853a6241bf304ebb6a47599ebc46c850e3768f1e47daa19aa`。
- 未启动 COMSOL，旧 `93602506` case 未覆盖、未删除；当前结果仍为
  `engineering_preview`，dose、polarity、生产 DRC 阈值及二次 SEM 闭环仍未标定。

## 显式输出命名调整计划

用户要求最终 case 和 DXF 不再显示 identity 短哈希，改用
`_(global_delta_nm,envelope_modulation_strength)` 作为末尾后缀；当前值必须精确显示为
`_(-10,1)`。本次范围包括 `workflow/output/reporting.py`、命名测试、README 和
`blueprints/AGENTS.md`，并重新生成一个非覆盖的新名称 case。完整 64 位配置 identity
继续写入 `layout_report.json`，用于来源追溯；case-local 临时 staging 可继续使用短哈希，
最终化后必须删除。若不同完整配置映射到同一显式名称，沿用现有 no-overwrite gate，拒绝
覆盖而不是静默复用。验证包括完整 blueprints 测试、py_compile、preflight、正式生成、
DXF 回读、PNG 检查、严格三文件合同和 `git diff --check`。回滚时仅恢复 reporting 命名
函数、测试和文档；既有哈希目录与新显式目录均保留。

### 显式输出命名执行结果

- `workflow/output/reporting.py` 已停止把 `identity_sha256[:8]` 写入最终 case/DXF 名称，
  改为从 `fabrication.triangle_side_correction` 直接格式化
  `_(global_delta_nm,envelope_modulation_strength)`；当前精确后缀为 `_(-10,1)`。
- 正式 case 为
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_(-10,1)/`；
  DXF 为
  `dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_sem-air-hole-220nm-radial-v1-unassigned-r1_(-10,1).dxf`。
- 完整 identity 仍为
  `a02f4d079011efcfac477c745d4aa5bd244974210a3c8dc1d855c9b269fa8ca5`，只保留在
  `layout_report.json`；测试显式锁定其前 8 位不出现在最终 case 和 DXF 名中。
- preflight：29,648 features，DRC `0 error / 0 warning`；完整 blueprints 测试
  `91 passed`，相关 `py_compile` 通过。新旧 DXF 回读均为 29,648 个闭合
  `LWPOLYLINE`、唯一 `ETCH` layer，逐实体逐顶点几何完全一致。
- 新 DXF SHA-256 为
  `27082e74385382008294924aa1936ff74a6471eba2db9f81f28884f0cd221325`，preview 为
  `e98ce8406edd2c609ff9c896c84a1cd82ca7a7ffb67a8ad8748a7ec5d9131053`，report 为
  `343068282a87aa919c4d00d855a38b1008a856f7fccbc663198f25dd0ecb9b95`。
- 实际 PNG 已检查且与命名前版本内容一致；旧 `93602506`、`a02f4d07` case 均未删除或
  覆盖。未启动 COMSOL，release 状态仍为 `engineering_preview`。

### 整体边长修正参数变体生成

按用户要求，在 `envelope_modulation_strength=1.0` 不变的条件下，分别生成
`global_delta_nm=-15.0` 与 `-20.0` 两个非覆盖 case；生成后正式
`main/parameter.json` 已恢复为 `-10.0`。

- `(-15,1)`：目录
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_(-15,1)/`；
  identity `3423ea105dd9c71afdeb700803351c232b0dc97fa5172f4bbe11630a56395ffa`；
  DXF SHA-256 `7d89c90ffeb621985a5a2412d07fc00a64a3f15e24fc86429d9818227328a96f`。
- `(-20,1)`：目录
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_(-20,1)/`；
  identity `4ce9e3e246501c29b8693ce88bad36945f87b074e650568a023a4e92a0527f62`；
  DXF SHA-256 `73b092be6297da673d5e7d16c3b9362fbc5c948c9e08d84908e8f4f3aa198f19`。
- 两个 case 均严格包含 DXF、`preview.png`、`layout_report.json` 三个文件；各自 DXF
  均为 29,648 个闭合 `LWPOLYLINE`、唯一 `ETCH` layer，DRC 均为
  `0 error / 0 warning`。两张实际 preview 已检查，整体版图、三处三角孔、修正曲线和
  径向 cell strip 均完整且未裁切。未启动 COMSOL，既有 case 未覆盖或删除。

### 单个 DXF 工作流耗时基线

2026-08-28 在当前 Windows 工作站、29,648 features 的正式完整流程中：

- `(-15,1)`：内部各阶段累计 `31.241785 s`；
- `(-20,1)`：内部各阶段累计 `30.576670 s`；
- 两次平均：`30.909227 s`，因此当前规模可记录为单个完整 DXF case 约 `31 s`。

这里的完整流程包含 nominal/辅助几何构建、加工修正、DRC、DXF 写入与回读、preview
生成；不启动 COMSOL。若只统计 DXF 序列化与回读阶段，两次分别为 `7.215936 s` 与
`7.007932 s`，平均 `7.111934 s`。该数值是当前工作站与当前 29,648-feature 配置的实测
基线，不作为不同硬件、不同 feature 数或并发负载下的固定性能保证。

## 15-15 层、9 um 外窗参数系列执行计划

用户要求把正式 finite 参数改为 cavity/cladding 各 15 层，把外层梯形开窗法向高度改为
`9.0 um`，并在 `envelope_modulation_strength=1.0` 下分别生成
`global_delta_nm=-10.0/-15.0/-20.0` 三组结果。正式 `main/parameter.json` 最终保持
`15-15、9.0 um、-10.0 nm` 基准值；两个额外整体修正只用于生成非覆盖变体。

本次不改变现有 SEM 包络的物理空间尺度 `half_size_um=[32.8,32.8]`、拟合系数、开窗
内净距 `3.8 um`、相邻斜边间距 `5.0 um` 或内角圆角 `3.0 um`。同步更新严格参数测试、
正式命名测试和 README 当前尺寸说明。每组先 preflight，再生成严格三文件 case；验证
feature 数、DRC、DXF 闭合/layer/单位、报告内实际参数、PNG 视觉结果以及 no-overwrite。
完成后运行完整 blueprints 测试、相关 py_compile、`git diff --check` 并清理缓存，不启动
COMSOL。回滚时把 cavity/cladding 层数恢复为 20、开窗高度恢复为 13 um；已生成 case
保持不动。

### 15-15 层、9 um 外窗系列执行结果

- 正式 `main/parameter.json` 已改为 cavity/cladding 各 15 层、
  `outer_etch_windows.radial_width_um=9.0`，并保持基准
  `global_delta_nm=-10.0`、`envelope_modulation_strength=1.0`。SEM 包络
  `half_size_um=[32.8,32.8]` 及其他开窗尺寸未改变。
- 已生成三个非覆盖 case，目录分别以 `(-10,1)`、`(-15,1)`、`(-20,1)` 结尾；每个
  case 严格只含参数化 DXF、`preview.png`、`layout_report.json`。
- 三组各含 16,868 features；DXF 回读均为 16,868 个闭合 `LWPOLYLINE`、唯一
  `ETCH` layer、单位 `micrometers`；DRC 均为 `0 error / 0 warning`。
- 三组 identity 依次为
  `602064840cf497c1251fa36e0e6b65ed5304ad5f3c389bbe15a98afb1e83194e`、
  `9853e9e11fb67d43a063a8cc884e3a2c5477f5b7a9f0b3f25228a002cbccbcc5`、
  `4b17cc62ec72503143777cea8d4b78f12776bb175824e308bb579cdbe7240b2a`。
- DXF SHA-256 依次为
  `19ab55973a63e41eaf2c0899965f9109dd3d48847f12557e58286e4c59341707`、
  `e8bee9a386c9962c3abd087b4ee9aa4f2052ef41f9ad632c0decba5f1a04e1b5`、
  `d87fdd9ed7a3f8561971bf1eb1e5d3795ae480d0bc0da556a6a1e5edef4cdb06`。
- 完整阶段累计耗时分别为 `18.108954 s`、`18.069656 s`、`18.235433 s`，平均
  `18.138014 s`；纯 DXF 写入回读平均约 `3.674 s`。完整 blueprints 测试
  `91 passed`，相关 py_compile 通过。
- 三张实际 PNG 已联合检查：15/15 边界、9 um 六段外窗、圆角、参数标记、中心/cavity
  边界/cladding 边界三角孔、修正曲线及径向 strip 均完整且未裁切。未启动 COMSOL，
  既有 20-20 case 未覆盖或删除。

## Window height 线性层数关系执行计划

### 方案与数据流

根据 `20-20 -> 13 um` 与 `15-15 -> 9 um` 两个标定点，将当前 v2 外窗高度从固定
`radial_width_um` 改为显式线性规则：

```text
finite_total_layers = cavity_layers + cladding_layers
window_height_um = 0.4 * finite_total_layers - 3.0
```

因此 15+15 解析为 9 um，20+20 解析为 13 um；未来 cavity/cladding 不相等时仍按总层数
唯一确定高度。正式 `main/parameter.json` 写入规则类型、斜率和截距，不同时保留固定高度，
避免两套参数源。配置加载阶段联合校验 finite 层数与解析高度，几何构建阶段复用同一解析器，
并在 feature/layout metadata 和最终 report 中记录规则、总层数和实际高度。

### 范围、兼容与验证

- 更新严格 schema、外窗高度解析器、辅助几何汇总、正式参数、README 和相关测试。
- 当前 v2 正式入口采用线性规则；历史 v1 外窗及直接传入固定 `radial_width_um` 的低层调用
  保持兼容，但固定值与规则不得同时出现。
- 验证 15-15 为 9 um、20-20 为 13 um、20-15 为 11 um，并拒绝未知字段、非法规则、
  非正解析高度和缺少 geometry 的规则调用。
- 运行聚焦及完整 blueprints 测试、相关 `py_compile`、CLI preflight 和
  `git diff --check`；不启动 COMSOL、不覆盖或删除已有正式 case。本次需求不要求重新生成
  DXF，避免与已有显式命名目录冲突。

### 回滚方案

将正式参数恢复为固定 `radial_width_um=9.0`，并回滚本节涉及的 schema、解析器、测试和
文档即可；已有 15-15/20-20 输出保持不动。

### 执行结果

- `main/parameter.json` 已用 `linear_from_total_finite_layers_v1` 显式写入斜率
  `0.4 um/total-layer` 与截距 `-3.0 um`，不再重复写固定高度。
- `workflow/config.py` 提供唯一解析器并执行严格联合校验；规则与固定值互斥，规则字段、
  数值有限性、正斜率及当前层数下的正高度均受检验。
- `workflow/geometry/auxiliary.py` 使用同一解析器构建六个开窗；每个外窗和 layout 汇总
  metadata 均记录实际 `radial_width_um`、规则内容与 `finite_total_layers`，可进入最终
  report。
- 数值测试确认 15-15、20-20、20-15 分别解析为 `9.0`、`13.0`、`11.0 um`；当前
  15-15 外窗几何法向高度仍精确为 `9.0 um`。
- 聚焦测试 `43 passed`；完整 blueprints 测试 `99 passed in 6.78s`；相关
  `py_compile` 通过。CLI preflight 为 16,868 features，DRC
  `passed=True`、`0 error / 0 warning`，没有写入结果文件。
- `git diff --check` 通过；未启动 COMSOL，未生成、覆盖或删除任何已有 DXF case。
