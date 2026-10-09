# Unit-cell band pair/solo 入口重构计划与执行记录

## 方案概述

将现有同时计算 cavity 与 cladding 的 `run_unit_cell_band.py` 明确重命名为
`run_band_pair.py`，并新增 `run_band_solo.py`。solo 入口只从
`scripts/parameter.json` 解析 cavity compact cell、unit-cell eigensolver 数量和
mesh，不求解或比较 cladding。

## 范围

- 重命名双结构主入口及其活动测试、CLI 和 sweep 导入引用。
- 新增 cavity-only 主入口和聚焦测试。
- solo 结果仍写入 `scripts/.out/unit_cell_band/`，系列目录名为
  `unit_cell_(<cavity-b0>-<eta>-<zeta>)_mesh<mesh>`。
- 更新根 README、scripts README 和 pyproject CLI 说明。
- 不修改、移动、覆盖或删除既有计算结果；不启动 COMSOL。

## 数据流

```text
scripts/parameter.json
  -> parameter_config.ACTIVE_PARAMETERS
  -> run_band_pair: cavity + cladding -> pair comparison outputs
  -> run_band_solo: cavity only -> single-cell band/Q outputs
  -> scripts/.out/unit_cell_band/<derived-series-name>/
```

solo 复用 pair 模块的单 case 求解、band tracking、Q 拟合和基础绘图逻辑，避免复制
COMSOL 与数值实现。solo 的配置 metadata 只记录 cavity case；共享参数快照仍用于
可复现性，但 cladding 不参与求解、采样点构建、case 目录或图形输出。

## 兼容策略

- `run_band_pair.py` 保持原双结构输出布局和成功后用 cavity `p2@Gamma` 更新共享中心
  频率的行为。
- sweep 改为导入 `run_band_pair`，保留其现有算法与输出。
- 删除旧 Python 入口名，避免三个名称表达同一工作流；CLI 提供明确的 pair/solo
  命令。
- solo 使用新的非冲突系列名，不复用旧 pair 缓存目录。

## 绘图合同

- pair 图保持既有双面板比较合同不变。
- solo band 图只含 cavity 一个面板，数据来自 cavity `selected_bands.csv`，频率窗口
  仍以 cavity `p2@Gamma` 居中，颜色和线型沿用既有 unit-cell 科研绘图约定。
- solo Q 图只含 cavity 数据，横轴、颜色、拟合及注释规则沿用既有单 case 逻辑。
- 输出位于系列 `10_overview/`，不生成 cavity/cladding comparison 文件。

## 验证与回滚

- 聚焦测试覆盖入口唯一性、sweep/CLI 引用、solo 参数隔离、精确系列名、单 case
  编排、输出布局与中心频率更新。
- 运行 pair/solo 测试、相关 `py_compile` 和 `git diff --check`。
- 如需回滚，只回退本计划列出的入口、测试、CLI 与活动文档；计算结果不参与回滚。

## 执行记录

- 将 `scripts/run_main/run_unit_cell_band.py` 重命名为 `run_band_pair.py`，同步
  重命名原有 pair 测试，并更新五个 unit-cell sweep 的导入。
- pair 的 band/Q 绘图抽为可渲染一个或多个 case 的公共函数；原 pair 包装函数、
  双面板标题、共享 scale、颜色和输出名保持兼容。
- 新增 `run_band_solo.py`：只编排 cavity case，配置 metadata 不记录 cladding case，
  不创建 cladding 结果或比较图；系列名从 cavity compact 参数和 mesh 动态派生为
  `unit_cell_(<b0>-<eta>-<zeta>)_mesh<mesh>`。
- solo 输出单面板 `unit_cell_band.png` 与 `p_bands_q_vs_k.png`，完整成功后仍以 cavity
  `p2@Gamma` 回写中心频率，来源标记为 `run_band_solo`；回写一致性检查不包含未参与
  solo 计算的 cladding。
- CLI 更新为 `comsol-band-pair` 与 `comsol-band-solo`；活动 README、核心模块 README
  和长期 Memory 已同步。`parameter.json` 中既有 `center_frequency_source` 是历史运行
  记录，未伪造为新入口。
- 无 COMSOL 聚焦测试：160 passed。受影响入口与 sweep 的 `py_compile` 通过，
  `git diff --check` 通过。
- 使用测试数据生成了 70,592-byte solo PNG；图像查看器因 Windows sandbox helper
  连接超时无法显示，自动绘图合同测试通过，两个任务临时 PNG 均已删除。
- 未启动 COMSOL，未创建或修改 `scripts/.out/` 正式结果，未改动用户已有 finite
  cavity 工作树内容。

## 最终状态

实施完成。
