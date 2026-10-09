# Blueprints DXF 逐孔几何覆盖 workflow 执行计划

日期：2026-09-01
状态：工程实现与验证完成；等待用户提供正式 DXF

## 1. 目标

在保持当前六角 finite-cavity 拓扑、region、外围开窗、加工修正、DRC 与三文件输出合同的
前提下，允许用户提供一个不区分 cavity/cladding layer 的二维 DXF。DXF 中每个闭合
polyline 表示一个尚未施加 Blueprints 整体边长修正和空间包络修正的等边三角孔；每个孔
可以相对参考结构独立改变位置、边长和旋转角。

程序以现有紧凑参数构建参考 `LayoutPlan`，将 DXF 三角孔与参考孔按质心做严格一对一匹配，
继承参考孔的稳定 `feature_id`、cavity/cladding region、parent cell、shell/layer 和 local
index，再用 DXF polygon 覆盖 nominal geometry。之后继续执行现有外围开窗、参数标记、
整体三角形边长修正、空间包络、dose、layout transform、writer-grid snapping、DRC、DXF
写入回读和综合预览。

本任务不读取 `.mph` 或 COMSOL Java，不启动 COMSOL，不修改或删除既有正式结果。

## 2. 首版输入合同

- DXF 坐标必须与参考六角结构处于同一坐标系，配置显式声明坐标单位为 `um`；不进行隐式
  平移、旋转、镜像或自动配准。
- 只读取 model space；其中所有绘图实体都必须是闭合、二维、零宽、无 bulge 的
  `LWPOLYLINE`，每个实体去除重复闭合点后必须恰有 3 个不同顶点。
- DXF 无需按 layer 区分 cavity/cladding；所有 layer 一视同仁。region 由参考孔匹配后继承。
- DXF 孔数必须与参考布局孔数完全相同；不允许新增、删除、重复匹配或未匹配孔。
- 每个三角形必须在显式相对容差内保持等边；每个导入质心与匹配参考质心的距离必须不超过
  显式阈值。程序还必须由当前参考几何计算最小孔中心间距，并要求该阈值严格小于间距的一半，
  以排除无 layer/ID 输入发生孔身份互换。
- 输入 DXF 的 SHA-256、用户可读 `model_id` 和 revision 必须写入配置。生成前校验实际
  SHA-256；不匹配时拒绝运行。
- DXF 表示待加工修正的 model nominal geometry，不能预先包含
  `triangle_side_correction` 或 `spatial_envelope`，避免双重修正。

## 3. 配置与兼容策略

- 保留现有 `geometry.schema_version=1` 及其全部行为和输出名称不变。
- 新增 `geometry.schema_version=2`；仍保留原有 `lattice/unit_cell/finite/cladding_shift`
  参考结构字段，并增加严格的 `model_geometry`：

```json
{
  "kind": "dxf_triangle_holes_v1",
  "model_id": "optimized-cavity-v1",
  "revision": "1",
  "path": "../data/models/optimized-cavity-v1.dxf",
  "sha256": "<64 hex characters>",
  "coordinate_unit": "um",
  "layer_policy": "all_modelspace",
  "matching": {
    "kind": "nearest_reference_centroid_v1",
    "max_centroid_distance_um": 0.1,
    "equilateral_relative_tolerance": 0.01
  }
}
```

- 配置路径相对当前参数文件所在目录解析；读入文件可以位于项目内其他目录，但不得由生成器
  修改。
- geometry/model 配置与 DXF 哈希进入完整 identity。case 和 DXF 名增加用户可读
  `_model(<model-id>-r<revision>)`，避免不同逐孔优化模型拥有相同基础紧凑参数时发生名称
  冲突；完整哈希仍只写入报告。
- 官方 `main/parameter.json` 在用户尚未提供正式 DXF 前保持 schema v1，不写入失效路径。

## 4. 数据流与模块边界

```text
parameter.json schema v2
  -> 严格配置校验
  -> 以 v1 参考字段构建 base hex LayoutPlan
  -> 读取并校验 DXF model-space triangles + source SHA-256
  -> STRtree 最近参考质心匹配
       -> 距离阈值
       -> 一对一双射
       -> 完整数目
  -> imported polygon 覆盖 base feature nominal/current polygon
       -> 保留 base identity/region/topology metadata
       -> 记录逐孔质心、边长和角度差异
  -> 添加现有 outer windows 与 parameter marker
  -> 现有 correction/dose/transform/grid pipeline
  -> 现有通用 DRC
  -> DXF readback + preview.png + layout_report.json
```

DXF reader 和 override builder 只负责 nominal geometry，不实现加工修正，不依赖外层
`scripts`、`comsol_workflow` 或 `mph`。后续若增加 `.mph/.java` 支持，应先独立导出满足本
合同的 DXF 或中性快照，再复用本 workflow。

## 5. 审核、报告和失败策略

- preflight 不写文件，但输出 model ID、实体数、匹配数、最大/平均质心位移、边长变化范围
  和最大等边误差。
- `layout_report.json` 记录输入路径的配置形式、实际 SHA-256、DXF header 单位/layer、匹配
  参数和聚合 QA；不把上万个逐孔记录复制到额外 sidecar。
- 输入包含其他 entity、开放 polyline、bulge、非平面坐标、非三角形、无效 polygon、哈希
  不符、计数不符、重复/缺失匹配、距离超限或非等边超限时立即失败。
- 参考 finite boundary 继续作为内部孔 containment 边界；外围开窗和标记沿用现有
  `boundary_check=false` 合同。
- 当前参数标记仍显示参考 compact parameters；报告明确标记实际孔形来自 model DXF。
- release gate 不变；首版按 engineering preview 验证，不把未标定工艺描述为生产结果。

## 6. 验证计划

- DXF reader：正常三角形、乱序实体、任意 layer；拒绝 hash mismatch、open/bulged/nonplanar、
  非三角形、退化或非等边 polygon。
- 匹配：平移/缩放/旋转后的每孔 polygon 正确继承 base feature identity 与 region；拒绝计数
  不同、重复 nearest reference、距离超限和非双射。
- correction：导入 polygon 成为 nominal，整体和包络修正在其实际质心位置继续施加；外围
  开窗正常加入且不受三角孔专用修正。
- CLI：schema v1 原有名称与行为不变；schema v2 名称包含 model ID/revision；preflight 不写
  文件；成功 case 仍严格只有参数化 DXF、`preview.png`、`layout_report.json`。
- 运行 Blueprints 独立完整 pytest、外层 geometry parity 聚焦测试、相关 `py_compile`、实际
  小型 DXF preflight、`git diff --check`；不启动 COMSOL。

## 7. 回滚方案

删除 schema v2 分派、DXF reader/override builder、model 命名后缀和对应测试即可；所有
schema v1 代码路径与既有 case 保持不变。官方参数文件默认不切换到 v2，因此回滚不需要
迁移正式配置或结果。

## 8. 执行记录

### 实际改动

- 新增 workflow/geometry/dxf_override.py：只读加载 model-space DXF，校验输入 SHA-256、
  单位、entity 类型、闭合/平面/零宽/无 bulge 合同、三顶点和显式等边容差。
- 新增 versioned geometry dispatcher。schema v1 继续调用原 HexGeometrySource；
  schema v2 先生成完全相同的 hex reference LayoutPlan，再按导入质心通过 STRtree 最近邻
  建立严格双射，用 DXF polygon 覆盖 nominal/current polygon，并保留原 feature ID、
  cavity/cladding、parent indices、shell/layer 和 local index。
- 匹配器从实际参考布局计算最小孔中心间距，配置阈值必须严格小于其一半；同时检查实体数
  完全相同、无重复 nearest reference、无遗漏和逐孔距离不超限。聚合记录最大/平均位移、
  边长变化范围、最大等边误差、DXF layer/header 和 source hash。
- CLI 改用统一 build_geometry()，相对输入路径按参数文件目录解析。preflight 增加 model
  ID/revision、匹配数和最大/平均质心距离输出；layout report 的 nominal/corrected summary
  均保留 model import QA。
- model ID/revision 进入 case 与 DXF 可读名称；完整 hash 仍只进入 identity/report。
  schema v1 原名称测试和官方 preflight 保持不变。
- 新增长路径文件系统适配，DXF 写入/回读、hash 和 staging 内容检查在 Windows 使用
  extended-length path；同时统一当前实际 staging 前缀 .s- 的创建与异常清理检查。
- README 和 Blueprints AGENTS 已记录正式输入、匹配、修正、命名和后续 MPH/Java adapter
  边界。官方 main/parameter.json 未切换 schema，也未写入虚构 DXF 路径。

### 验证结果

- 新增 8 项 DXF workflow 测试：乱序和任意 layer、逐孔平移/边长覆盖、source hash、
  缺孔、非双射、身份安全阈值、非三角形、schema 严格性，以及 schema v2 CLI preflight
  和三文件实际导出。
- schema v2 集成 fixture：114 个导入孔与 reference 全部匹配，加入 6 个外围开窗后共
  120 个 feature；现有整体边长/包络 correction、grid、DRC、DXF readback 均通过，成功
  case 严格只有参数化 DXF、preview.png 和 layout_report.json。
- 实际 preview PNG 为 3454 x 1804；已完成视觉 QA：孔阵列、6 个开窗、finite boundary、
  三个局部比较面板、径向修正函数和 radial strip 均完整，无裁切或标题碰撞。
- Blueprints 独立完整测试：107 passed；外层纯 Python geometry parity：2 passed；
  相关 py_compile 通过。
- 官方 schema v1 preflight：16868 features，DRC 0 error、0 warning，不写文件，名称保持
  dxf_15-15_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_(-10,1)。
- git diff --check 通过；未启动 COMSOL，未修改或删除既有正式结果。

### 遗留问题与最终状态

当前状态为 DXF holes-only engineering workflow ready。正式 main/parameter.json 仍保留
schema v1；获得用户真实 DXF 后，需要核对坐标单位/坐标系、孔数、实际 SHA-256、最大逐孔
位移和等边误差，再选择不超过 identity-safe 上限的匹配阈值并切换 schema v2。生产 release
仍受原 fabrication profile、dose、polarity 和 DRC 标定门槛约束。
