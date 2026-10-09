# finite 层数与 20-20 外窗基准修正

## 方案概述

当前 hex finite 几何把 `cavity_layers=N` 解释为晶格 shell `0..N`，实际包含
`N+1` 个以中心为起点的层。用户要求改成严格的 1-based 层数语义：中心 cell 是第 1 层，
因此 cavity 的最大 shell 必须为 `N-1`；紧邻 cavity 的第一层 cladding 位于 shell `N`，
其 cladding shift profile 的局部层号仍为 1。

外围开窗以 20-20 为最小参考结构。对总 finite 层数
`cavity_layers + cladding_layers <= 40` 的结构，窗口参考边界、位置和径向宽度全部固定为
20-20 基准；总层数大于 40 时，参考边界随真实结构增大，径向宽度继续使用现有线性公式。
选择总层数作为判据，是因为当前 hex 外边界和既有窗口线性规则都只依赖 cavity 与
cladding 的层数和，而不依赖二者如何分配。

## 范围

- 修正 `comsol_workflow` 的 hex finite plan 与 finite/strip 主入口中的直接 hex 构建路径；
- 同步修正独立 `blueprints` hex reference，使 schema v1 与 holes-only DXF override
  使用同一新参考拓扑；
- 为 v2 外窗线性规则增加显式 `minimum_total_layers`，避免在源码内隐藏 40 层常数；
- 更新 Blueprints 预览中的 cavity boundary shell 选择、测试、README 和参数示例；
- 不改变 square finite 的矩形裁剪语义，不启动 COMSOL，不覆盖既有结果。

## 数据流与几何合同

```text
cavity_layers=N
  -> cavity shell 0..N-1（对外第 1..N 层）
  -> cladding shell N..N+M-1
  -> cladding local layer = shell - N + 1

actual_total_layers = N + M
window_reference_total_layers = max(actual_total_layers, minimum_total_layers=40)
  -> <=40：20-20 reference boundary + 20-20 radial width
  -> >40：actual finite boundary + 原线性 radial width
```

窗口的 `inner_clearance_um`、`parallel_side_gap_um`、圆角和 DXF layer 保持不变。
`minimum_total_layers` 缺失的历史配置继续按旧式无下限线性规则解析；正式参数显式写入 40。

## 兼容与缓存策略

- schema v2 DXF 输入仍须与当前 reference holes 等数并严格一一匹配；层数语义改变后，旧
  holes-only DXF 若按旧拓扑导出会因 feature count 或质心不匹配而被明确拒绝，不静默迁移。
- `LayoutFeature.layer` 与 `parent_layer` 元数据统一表示从中心开始的绝对 1-based 层号；
  `parent_shell` 保持 0-based 晶格坐标。cladding shift profile 的局部层号独立记录为
  `cladding_shift_layer`，避免把局部调制层号与全结构绝对层号混为一谈。
- 不复用、移动或删除现有 finite/strip/blueprint 结果；本轮只做无 COMSOL 验证。

## 验证计划

- 单元测试锁定 `1-1`、`2-2` 的 cavity/cladding cell count、shell 范围和两边 parity；
- 外窗测试锁定 15-15、20-20 的六个 polygon 完全相同，以及 20-21 按规则增大；
- schema v2 fixture 随新 reference count 更新并验证严格匹配；
- 运行 Blueprints 完整测试、相关外层 finite/parity/strip 测试、`py_compile`、官方参数
  preflight 和 `git diff --check`；若生成 PNG，实际打开检查。

## 回滚方案

回滚本计划涉及的源码、测试、README 与正式 Blueprints 参数字段即可恢复旧行为。既有输出
未被修改，因此无需数据回滚。

## 执行记录

### 实际改动

- Blueprints 与 `comsol_workflow.finite_geometry` 的 hex cavity 改为 shell `0..N-1`，
  cladding 改为 shell `N..N+M-1`；finite/strip 主入口的 `BULK_RADIUS` 现在明确为
  `cavity_layers - 1`，直接 hex 构建、normal-fan、validity boundary、Fourier radius
  均沿用实际最大 shell。
- Blueprints `LayoutFeature.layer` 与几何 record 的 `parent_layer` 使用绝对 1-based
  层号；cladding profile 局部层号保存在 `cladding_shift_layer`，第一层仍为 1。
- 外窗线性规则新增显式 `minimum_total_layers`。正式参数设为 40；窗口构造器使用共享
  hex footprint 几何生成精确参考边界，而非近似比例放大。15-15 与 20-20 的六个窗口
  polygon、位置和面积相同；20-21 的参考总层数为 41、径向宽度为 13.4 um。
- Blueprints 径向预览的 cavity boundary 选择从 `cavity_layers` shell 修正为
  `cavity_layers - 1`。README、根/子项目 AGENTS 及测试合同已同步。
- schema v2 DXF override 继续从当前 reference 构建，因此 1-1 fixture 从旧 114 个孔更新
  为 42 个孔；旧层数语义导出的 holes-only DXF 会被等数/质心匹配门禁明确拒绝。
- finite-quarter 四模入口增加 cavity shell/边缘 cell 数量门禁，并在生成的 metadata 中显式
  保存 `cavity_layers`；以当前 `cavity_layers=20` 为例，强制 cavity shell 为 `0..19`，
  每条边为 20 个 cell。
- 四模入口增加受 `scripts/.out/finite_cavity` 约束的 `--output-dir` 覆盖参数，以便错误
  结果清理后在用户指定的新目录中重算，同时保留其他历史 finite 结果。

### 验证结果

- Blueprints 完整测试：`110 passed`。
- 根项目完整测试：`525 passed, 3 warnings`；三条 warning 均为既有
  `test_strip_1d_fourier_quotient.py` 测试函数返回非 None。
- 官方 `blueprints/main/parameter.json` preflight：15-15 新结构 `15788 features`，
  DRC `0 errors / 0 warnings`，未写输出。
- 相关模块 `py_compile`、两份参数 JSON 解析与 `git diff --check` 通过。
- 已检查测试实际生成的 schema v2 preview PNG：孔阵列、finite 边界、六个窗口、径向带
  和局部对比均未裁切；用于检查的仓库内临时 PNG 已删除。
- 未启动 COMSOL；未创建、覆盖、移动或删除任何正式计算/版图结果。
- 本轮针对中止的旧四模任务额外完成了只读 preflight 逻辑核对；旧结果中的 `bulk_radius=20`
  不被恢复，后续新任务必须通过 20-cell/边门禁。

### 遗留与最终状态

- 本次规则只改变当前 hex finite/strip；square finite 的矩形裁剪层数语义保持原样。
- 既有结果记录的是旧几何，保持只读历史状态。新任务必须重新构建，不能把旧 MPH/cache
  当作新 1-based 结构继续恢复。
- 本计划已完成。
