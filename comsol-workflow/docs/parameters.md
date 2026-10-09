# 主程序参数入口

所有主程序继续读取 [scripts/parameter.json](../scripts/parameter.json)，
由 [parameter_config.py](../scripts/run_main/parameter_config.py) 校验。
当前格式为 `schema_version: 2`，只存标准 JSON 数据，不放注释或说明字段。

## 分区与影响范围

| 路径 | 内容 | 使用入口 |
| --- | --- | --- |
| `mesh_size` | 1–9 的统一 mesh 等级 | 所有主程序；finite/quarter 使用对应手动 mesh 等级 |
| `cells.cavity`、`cells.cladding` | 各自的 `b0_nm`、`eta`、`zeta` | band pair、strip、full/quarter finite；band solo 只使用 cavity |
| `structure.cavity_layers`、`structure.cladding_layers` | 腔与包层层数 | strip、full/quarter finite |
| `structure.center_frequency_thz` | 求解和默认模式选择的目标频率，单位 THz | strip、full/quarter finite |
| `structure.center_frequency_source` | 中心频率的来源记录，不是另一份生效参数 | 供追溯；band 主程序成功后更新 |
| `structure.cladding_shift_factors` | 远端位移因子列表，即使单点也用列表 | strip 直接扫描；finite 默认作为 x 因子扫描 |
| `structure.cladding_shift_profile` | 逐层绝对位移分布：`uniform` 或 `tanh_power` | strip、full/quarter finite 共用 |
| `unit_cell_band.eigenmode_count` | band 求解模式数 | band pair/solo |
| `unit_cell_2d` | 独立单胞几何、模式对数、二维 k 网格、跟踪与分析、valley 加密参数 | 2D 主程序及对应扫描、分析 |
| `strip_1d` | `eigenmode_count`、`cut_angle_deg`（0 或 60）、`kx` | 仅 strip；kx 沿用现有 Floquet 波矢约定 |
| `finite_cavity.geometry` | 完整 footprint：`hex` 或 `square` | full/quarter finite；缺失时仍为 hex |
| `finite_cavity.eigenmode_count` | 每次求解的模式数 | full/quarter finite 和四扇区入口 |
| `finite_cavity.cladding_shift_geometry` | `groups` 或 `ellipse` 方向位移规则 | 仅 full/quarter finite；缺失时仍为 groups |
| `finite_cavity.cladding_y_over_x_shift_ratio` | 从共享列表派生 y 因子的比例 | full/quarter finite 的默认扫描输入 |
| `finite_cavity.quarter_symmetry_ids` | 非空、无重复的 1–4 整数列表 | 普通 quarter 入口；四扇区专用入口仍固定求解 1→2→3→4 |

`cells` 和 `structure` 的共享是现有工作流的明确约定：同一组腔/包层设计依次用于
band、strip 和 finite 验证。2D 单胞有自己的几何，修改 `cells.cavity` 不会改变它。
修改共享 mesh 会影响全部入口；模式数则分别归属 band、2D、strip、finite。

Band 和 2D 的求解 shift 仍固定为 `c_const/1.55[um]`，不读取结构目标频率。
`unit_cell_2d.eigenmode_pair_count=M` 表示每点请求 M 对、预期 2M 个本征值。
2D 网格、内部 band 名、偏振参数和加密字段见 [2D 合同](agent-unit-cell-2d.md)。

## 常见修改

- 改腔/包层设计：编辑 `cells` 中相应三元组；b0 使用 nm，eta/zeta 无量纲。
- 改有限结构大小或目标频率：编辑 `structure`。Hex 的 `cavity_layers=N` 包含中心
  cell，对应 shell 0…N−1；逐层 shift 为相对理想晶格的绝对位移，不逐层累加。
- 改某入口的模式数：编辑对应工作流的 `eigenmode_count`，2D 则编辑自己的模式对数。
- 改 strip 边界方向：编辑 `strip_1d.cut_angle_deg`；不会旋转 finite 或单胞几何。
- 改普通 quarter 扇区：编辑 `finite_cavity.quarter_symmetry_ids`；1=PEC/PMC、
  2=PMC/PEC、3=PEC/PEC、4=PMC/PMC，顺序为 x/y 轴边界。

默认 finite 位移扫描保持 `x=factor`、`y=factor*ratio`。
如需单个显式 x/y 位移，在 `finite_cavity` 中移除 `cladding_y_over_x_shift_ratio`，
同时加入以下两项（示例值）：

```json
"cladding_x_shift_factor": 0.02,
"cladding_y_shift_factor": 0.04
```

显式 x/y 和 ratio 互斥；显式 pair 只影响 finite，strip 仍扫描
`structure.cladding_shift_factors`。位移 profile 仍从 `structure` 读取。

## 兼容、回写与运行边界

加载器继续支持 schema v1 的历史任务快照，不自动迁移历史文件。
当前 v2 文件拒绝旧根级参数、未知分区字段和缺失的 2D 几何，避免拼错参数后静默使用默认值。
Python 中的 `SharedParameters` 属性、结果元数据字段、目录名和缓存身份沿用现有合同；
文件分区变化不会使已有计算缓存失效。

只有完整 band pair/solo 主流程成功后，才把 cavity `p2@Gamma` 更新至
`structure.center_frequency_thz` 并更新其来源。回写前继续检查参与 band 的几何、mesh
与模式数是否发生变化；v1 快照则仍写回 v1 路径。无关参数和同期结构层数修改予以保留。
扫描只生成任务级快照，不修改共享文件；用环境变量 `COMSOL_WORKFLOW_PARAMETER_PATH`
指定快照的方式不变。

本轮只调整现有参数归属并接入 strip/quarter 的上述入口参数。
Band 路径点数、频率显示窗口，以及各入口的执行、恢复、场导出和后处理开关仍位于脚本中；
JSON 不会因此触发计算。启动真实计算仍按 [仿真合同](agent-simulation.md) 核对实际配置。
