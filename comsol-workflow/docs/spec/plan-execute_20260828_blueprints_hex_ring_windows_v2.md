# Blueprints 六边形环分段开窗 v2 执行计划

## 方案概述

将现有基于“内边端点缩进、外边切向延长、四角统一圆滑”的外层开窗 v1，重构为严格由 finite 六边形边界派生的六边形环分段开窗 v2。每个开窗的内侧参考边与对应 finite 六边形边等长且平行；相邻开窗在共享六边形顶点一侧的斜边沿同一径向方向，因而互相平行。只对靠近结构的两个内角做解析圆角，外侧两个角保持尖角。

## 范围

- 修改 `blueprints/workflow/geometry/auxiliary.py`，增加 v2 环段构造和选择性内角圆角。
- 修改 `blueprints/workflow/config.py`，增加严格的 v2 参数 schema，同时保留 v1 读取兼容。
- 将 `blueprints/main/parameter.json` 的正式工程配置切换为 v2，删除 v1 专属参数。
- 更新 `blueprints/tests/` 中的几何合同和配置合同。
- 更新 `blueprints/README.md` 的参数语义。
- 使用新 identity 生成一个新工程 case；不覆盖或删除旧 case。

## 几何合同与数据流

```text
finite 六边形边界
  -> 逐边计算外法线、边切向和中心到端点的径向
  -> 内边按 inner_clearance_um 沿外法线平移，长度保持不变
  -> 外端点沿两端径向推进，使法向宽度为 radial_width_um
  -> 仅在两个内端点应用 inner_corner_radius_um 解析圆角
  -> 6 个 outer_etch_window LayoutFeature
  -> correction / dose / grid / DRC / DXF / preview
```

当 `inner_clearance_um=0` 且 `inner_corner_radius_um=0` 时，相邻环段共享内外端点，六块的并集必须成为连续且无分段间隙的完整六边形环。正间距时，相邻斜边仍平行，并自然留下分段间隔。

## 参数与兼容策略

v2 使用以下显式参数：

- `kind = segmented_hex_ring_windows_v2`
- `reference = finite_boundary`
- `inner_clearance_um`
- `radial_width_um`
- `inner_corner_radius_um`
- `arc_segments_per_quadrant`
- `dxf_layer`

v2 不再接受 `inner_end_clearance_um`、`outer_tangential_extension_um` 和 `corner_radius_um`。v1 的严格 schema 与旧构造保留，确保历史参数仍可读取；仓库正式参数改用 v2。identity 由完整参数派生，因此新几何生成新目录，不覆盖旧输出。

## 预览绘图合同

- 目标文件：新 case 内的 `preview.png`。
- 数据源：同一次新参数流水线产生的 nominal/corrected polygon。
- 面板、配色、标题、归一化、坐标比例和径向放大带保持现有共享预览合同不变。
- 本轮只改变六个外层开窗 polygon；参数标记和纳米孔展示不改。
- 导出后按实际 PNG 检查六块对称性、内边位置、相邻侧边方向、内圆角和尖锐外角。

## 验证

- v2 schema：允许零内间距，要求正径向宽度、非负内圆角，拒绝未知/v1 专属字段。
- 几何合同：6 个合法 polygon；未圆角参考内边与 finite 边等长；相邻斜边平行。
- 零间距/零圆角：六块无面积重叠，并集等于完整六边形环。
- 选择性圆角：两个原内角被替换，两个外角精确保留，metadata 明确记录圆角范围。
- 当前正间距配置：六块互不相交。
- 运行 blueprints 完整 pytest、相关 `py_compile`、CLI preflight、实际生成、DRC、PNG 视觉检查和 `git diff --check`。
- 不启动 COMSOL。

## 回滚方案

将正式参数的 `outer_etch_windows` 恢复为 v1 字段即可复现历史开窗构造；旧输出目录保持不变。代码回滚仅涉及本计划列出的 geometry/config/tests/README 文件和本计划文档。

## 用户复核修订：增强内角圆角可见性

用户检查首版 `f4deb04d` 后确认 `1.0 µm` 内角圆角不够明显。保持 v2 几何算法、
内边/斜边合同和 `8.0 µm` 径向宽度不变，仅将正式参数
`inner_corner_radius_um` 从 `1.0` 调整为 `3.0`，并生成新 identity。修订后重新运行
聚焦测试、完整 blueprints 测试、preflight、实际生成、DXF 顶点核对和原像素 PNG 检查；
旧 case 保留不覆盖。

## 执行记录

### 实际改动

- `blueprints/workflow/geometry/auxiliary.py` 新增
  `segmented_hex_ring_windows_v2`：内边按外法线整体平移且保持原边长度，外端点
  沿六边形顶点径向推进；相邻环段侧边因此成对平行。
- 新增解析式选择性 fillet，只替换两个内侧角；两个外侧顶点原样保留。feature metadata
  记录参考边长、内角半径、圆角数量和外角半径。
- 保留 `segmented_hex_etch_windows_v1` 的历史构造与严格配置读取。
- `blueprints/workflow/config.py` 增加 v1/v2 分支 schema。v2 的
  `inner_clearance_um` 和 `inner_corner_radius_um` 允许非负值，径向宽度和圆弧离散
  仍要求为正；v2 严格拒绝 v1 专属字段。
- `blueprints/main/parameter.json` 切换到 v2，保留 2.5 um 内间距、8.0 um
  径向宽度和 1.0 um 内角圆角，删除端点缩进和切向延长参数。
- 更新 geometry/config/CLI identity 测试及 README。正式参数的新 identity 为
  `f4deb04d`。

### 验证结果

- 聚焦 geometry/config 测试：`29 passed`。
- blueprints 完整独立测试：`82 passed in 6.67s`。
- 相关 `py_compile`：通过。
- CLI preflight：29,648 features；DRC `passed=True`，0 error，0 warning；未写文件。
- `git diff --check`：通过。首次未限定测试目录的命令误收集外层仓库测试，
  因 blueprints 独立环境不含 pandas/scipy/jpype 而在收集期失败；按独立项目边界改为
  显式 `blueprints/tests` 后完整通过，该次环境错误不属于代码失败。
- 未启动 COMSOL。

### 新生成 case

`blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_f4deb04d/`

严格包含 3 个文件：

- DXF SHA-256：
  `47c574f122f6ee00ca67f7767535a775dcac03a1f343578373f9a32680bf2d1c`
- `layout_report.json` SHA-256：
  `e34e834ed4b331945baa65bfd423047de6bdb087ec39fe7f4291d1bb83e4e244`
- `preview.png` SHA-256：
  `4ed41822ae48d5a727dae853714511156d3ff3e98e09e14e5b11584409d20a5a`

实际 PNG 整图和原像素顶部开窗裁图均已检查：六块对称；内边位置与 finite 六边形
一致；相邻斜边平行；两个内角圆滑；两个外角尖锐；参数标记和径向孔带未改变。没有
创建视觉检查临时文件，也没有覆盖或删除旧 case。

### 遗留问题与最终状态

几何重构、参数切换、测试、预检、实际生成和视觉验收均已完成。当前 DXF 仍是
`engineering_preview`：dose 未分配、polarity 未指定，SEM `process_gain=1.0`
尚未完成二次闭环标定，因此不得描述为生产发布版图。v1 兼容路径保留，暂无本任务范围内
的代码遗留项。

### 用户复核修订执行结果

- 正式参数 `inner_corner_radius_um` 已从 `1.0` 改为 `3.0 µm`；
  v2 构造算法、内间距和径向宽度不变。
- 新 identity：`083368ca`。
- 聚焦 geometry/config/CLI 测试：`38 passed`。
- blueprints 完整独立测试：`82 passed in 6.48s`。
- preflight 与实际生成均为 29,648 features，DRC 0 error / 0 warning。
- 新 case 严格包含 DXF、`layout_report.json`、`preview.png` 三个文件：
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_083368ca/`。
- DXF SHA-256：
  `1e6f15b7ec7c18044d128d2ed5b2f9fb8e45506586cc92c9f11e830947b2d59a`。
- `layout_report.json` SHA-256：
  `d55c460647d6c38939d3f6ee6561348b4a1531863320b5a4d3dc8149cd3ca2ec`。
- `preview.png` SHA-256：
  `6b351f5d75eaefaaa48e2377330d6a8bc8477fc486fcd94f243e70804bbd92fb`。
- 直接回读 DXF 后，六个大开窗均为 16 顶点 LWPOLYLINE；圆角离散点跨度已扩大为
  约 `3.0 µm`。原像素顶部开窗裁图确认两个短边内角圆角清晰可见，外角仍为尖角，
  相邻斜边保持平行。
- 旧 `f4deb04d` 和更早 case 均保留；未启动 COMSOL，未 commit 或 push。
- 当前输出仍为 `engineering_preview`，生产标定遗留状态不变。
