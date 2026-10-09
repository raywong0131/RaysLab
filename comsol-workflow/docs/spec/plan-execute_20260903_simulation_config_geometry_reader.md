# simulation config 几何读取接口

## 方案概述

为当前 COMSOL workflow 增加一个可复用的 `*_simulation_config.json` 几何读取接口，
用于读取文件中的二维计算 footprint 和指定 layer 的显式 hole polygons。接口只消费
几何字段，不把源文件中的材料、物理场、边界条件、频率 shift 或 mode count 带入当前
求解流程；这些配置继续由 `scripts/parameter.json` 和现有 runner 决定。

本次完成底层读取与验证契约，并为后续 finite/quarter 几何构建提供稳定的数据结构。
没有启动 COMSOL，也没有修改共享 `parameter.json` 或既有正式结果。

## 范围

- 新增 `comsol_workflow.simulation_config_geometry` 模块。
- 支持读取当前 reference 目录中的 full/quarter 两类 JSON；默认读取
  `layers[0].holes`。
- 验证 JSON、长度单位、footprint、孔洞顶点、方向、退化情况和孔洞是否落在
  footprint 内。
- 返回只读语义明确的 NumPy 几何数组、footprint 和来源 SHA-256 摘要，供 COMSOL
  几何构建器使用。
- 增加针对真实 quarter 配置和错误输入的聚焦测试。

## 非范围与兼容策略

- 不读取或覆盖源 JSON 的 `materials`、`physics`、`boundary_conditions`、
  `shift_frequency`、`mode_count` 等求解配置。
- 不改变现有 compact geometry、layer 数、mesh、频率、输出目录和分析开关的参数源。
- 不把 JSON 中的矩形 footprint 自动解释为 fabrication outer window；其语义仅作为
  source geometry 的计算 footprint 返回。
- 源 JSON 不包含 region/cell 标签，因此读取器不猜测 cavity/cladding 归属，也不生成
  Fourier 所需的 cell 元数据；后续 runner 接入必须显式提供或验证对应拓扑，避免把
  `parameter.json` 的不匹配 cell 元数据静默套到外部几何上。

## 数据流

```text
simulation_config.json
  -> JSON/geometry contract validation
  -> SimulationConfigGeometry
  -> explicit footprint + holes for geometry construction
  -> parameter.json supplies physical/solver/postprocess configuration
```

## 实际改动

- 新增 `comsol_workflow/simulation_config_geometry.py`：提供
  `load_simulation_config_geometry(path, layer_index=0)`、
  `SimulationConfigGeometry` 和 `SimulationConfigGeometryError`。
- 新增 `tests/test_simulation_config_geometry.py`：覆盖真实 quarter 配置的 607 个三角
  孔洞、只读数组、源 hash、求解字段隔离和非法几何输入。
- 未接入 COMSOL 主入口；当前接口已经可以直接为 `LayerSpec.holes` 提供显式多边形，
  runner 接入时仍需单独处理 source geometry 与 parameter.json 拓扑元数据的一致性。

## 验证结果

- `uv run python -m pytest -q tests\\test_simulation_config_geometry.py`：9 passed。
- `uv run python -m pytest -q tests\\test_finite_cavity_field_exports.py tests\\test_run_finite_quarter_four_modes.py`：23 passed。
- 未启动 COMSOL，未生成或覆盖正式计算结果。
- 临时修复副本和测试占位文件已清理。

## 遗留问题与最终状态

读取接口已完成，可读取用户指定的 `quarter_simulation_config.json` 并返回 607 个有效
三角孔洞。若要实际启动该历史几何的 finite-quarter 求解，还需要将该接口接入具体
runner，并为后处理显式提供/验证 cavity-cladding cell 拓扑；本次不以不匹配的
`parameter.json` 几何元数据替代这一步。

最终状态：完成底层几何读取接口与测试，COMSOL 计算待后续 runner 接入后执行。
