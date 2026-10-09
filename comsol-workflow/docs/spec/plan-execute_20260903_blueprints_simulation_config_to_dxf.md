# simulation config 直接生成 Blueprints DXF

## 方案概述

将 `results/finite_cavity_reference_results/finite_cavity_xy_fourier26_bulk5_cladding6_wide_bounds_mma_move0p01_eval0031/full_simulation_config.json`
作为展开后的 polygon 几何源接入 Blueprints。源文件 `layers[0].holes` 已包含 full cavity
对称展开后的 2382 个等边三角孔，不再先导出临时 holes-only DXF；适配器直接读取 JSON、
校验哈希/单位/材料/层/三角形，并与当前 6-6 reference 做严格质心双射，继承
cavity/cladding、parent cell、shell 和稳定 feature identity。

2382 个孔对应 legacy bulk radius 5（shell 0..5）和 6 层 cladding（shell 6..11）；按当前
1-based 语义配置为 `cavity_layers=6`、`cladding_layers=6`。源 JSON 的矩形 footprint 是
COMSOL 计算域，只记录 provenance，不作为版图外窗边界；外围开窗继续使用当前小于等于
20-20 时固定为 20-20 基准的规则。

## 范围

- 在 schema v2 `model_geometry` 下新增 `simulation_config_triangles_v1`；
- 抽取 DXF override 已有的通用严格三角孔匹配器，由 DXF 和 simulation config 共用；
- 新增独立、可迁移的 JSON geometry source，不导入外层 `scripts`、`comsol_workflow` 或
  `mph`；
- 新增任务级 Blueprints 参数快照，运行 preflight 后生成正式三文件工程预览；
- 更新严格配置校验、dispatch、公开接口、README、测试和执行记录。

## 数据流

```text
full_simulation_config.json
  -> SHA256 + JSON schema/source contract validation
  -> layers[0].holes -> 2382 equilateral polygons
  -> current 6-6 hex reference strict centroid bijection
  -> inherit region/topology; imported polygon becomes nominal geometry
  -> 20-20 baseline outer windows
  -> current global side correction + SEM spatial envelope
  -> grid snap + DRC
  -> parameterized DXF + preview.png + layout_report.json
```

参数标记对该 Fourier 逐孔结构并不能完整表达几何，因此任务快照中关闭 markers，避免把
reference compact parameters 误当成完整模型描述。case/DXF 名仍包含可读 model ID 与
revision，完整源路径、实际哈希、layer 选择、源 footprint 和匹配 QA 写入报告。

## 输入合同

- `kind="simulation_config_triangles_v1"`；
- `coordinate_unit="um"` 且源 JSON `length_unit="um"`；
- 显式 `layer_index=0`、`layer_material="silicon"`、`hole_material="air"`；
- `footprint_policy="reference_hex"`，明确不把 COMSOL 矩形计算域当作版图边界；
- SHA256 必须匹配当前文件字节；每个 hole 必须恰有三个有限且不同的二维顶点、有效且在
  配置容差内等边；
- 孔数必须等于 reference，质心匹配必须一一对应，匹配阈值严格小于 reference 最小孔距
  的一半。

当前源的只读预检数据：2382/2382 唯一匹配，最大质心距离约 0.10983 um，身份安全上限
约 0.1312 um，因此任务配置使用 0.12 um；最大非等边误差约 8.34e-15。

## 兼容与回滚

- 既有 `dxf_triangle_holes_v1` 字段、metadata 和行为保持不变；
- schema v1 紧凑参数路径保持不变；
- 删除新增 source/dispatch/config 分支、任务参数和测试即可回滚；不修改源结果 JSON，
  不覆盖既有 Blueprints case。

## 验证计划

- 合成 fixture 覆盖成功导入、哈希、单位、layer/material、退化/非等边和严格匹配错误；
- 对目标真实 JSON 运行 preflight，核对 2382 matched、最大距离、DRC 和输出身份；
- 生成三文件 case，回读 DXF，核对实体数、报告源哈希、外窗数、修正 pipeline；
- 实际查看 preview PNG；运行 Blueprints 完整测试、相关外层 parity、`py_compile`、JSON
  解析和 `git diff --check`；不启动 COMSOL。

## 执行记录

### 实际改动

- 新增 simulation_config_triangles_v1 严格输入规格与 geometry source，直接读取指定
  simulation config layer 的 holes；校验文件哈希、um 单位、layer/material、三顶点数值、
  有限性、有效性和等边误差；
- 与 dxf_triangle_holes_v1 共用现有 reference centroid 严格双射逻辑，保留 schema v1
  和 DXF schema v2 行为；新增 builder dispatch、公开导出和按 kind 分支的严格 config
  schema；
- 新增独立 fixture 测试，覆盖成功匹配、源合同错误、身份安全、CLI preflight 和三文件
  导出；README 与 Blueprints AGENTS 同步记录 JSON 输入合同；
- 新增 blueprints/main/parameter_eval0031.json，使用 6-6 reference、关闭 markers，
  保留当前 global_delta_nm=-10、envelope strength 1 和 minimum total layers 40 的
  20-20 基准外围开窗。

### 真实数据执行结果

- 当前源 SHA256：
  f947c149d2e37b578e68b626d799894236577f1d757f2b918a7d97093c87a033；
- preflight：2382/2382 matched；最大质心距离 0.10982991494065933 um，平均
  0.04004694506020869 um；加入 6 个外围开窗后共 2388 features；DRC 0 error /
  0 warning；
- 正式 case：
  blueprints/main/.out/dxf_6-6_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_model(finite-cavity-xy-fourier26-eval0031-r1)_(-10,1)；
- case 严格包含参数化 DXF、preview.png 和 layout_report.json。DXF 回读为 AC1024
  (R2010)、INSUNITS=13 (micrometers)、2388 entities；报告 region counts 为 cavity
  546、cladding 1836、outer_etch_window 6；
- 已实际检查 PNG：完整结构、6 个固定基准外窗、中心/cavity/cladding 三个孔面板、
  shell 0..11 径向条带和修正曲线均完整，无标题碰撞或截断。

### 验证

- 聚焦：17 passed；
- Blueprints 完整 pytest：119 passed in 16.67s；
- 修改模块与新增测试 py_compile 通过；
- 任务参数标准 JSON 解析通过；
- git diff --check 通过，仅显示仓库既有 Windows CRLF 转换提示；
- 未启动 COMSOL，未修改源 simulation config 或既有正式结果；视觉检查临时副本已删除。

最终状态：完成。
