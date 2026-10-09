# `comsol_workflow/` 模块说明

本目录保存本机 COMSOL 控制器的可复用底层能力。`scripts/` 负责参数和任务编排，
这里负责几何、晶格、mesh、COMSOL 模型控制、模式识别与后处理。
模块应被导入使用，任务编排及命令入口位于 `scripts/`。

## 边界、dipolar 与统一输出

| 所有者 | 责任 |
| --- | --- |
| `boundary_integrals.py` / `boundary_plotting.py` | 原 S4 求积、Maxwell 项、奇偶性及绘图；保留共轭 SI 采样合同。 |
| `dipolar_radiation.py` | 复投影、能量 Gram、FFT/回溯、体辐射和周期模式匹配；保留原生场约定。 |
| `dipolar_decomposition.py` / `dipolar_lattice_analysis.py` | Gamma/REST 分解、共享节点投影、缝隙检查及晶格窗口/谱。 |
| `dipolar_inputs.py` | 显式来源清单和历史任务预设；不读取科学场数据或启动求解。 |
| `result_io.py` | 哈希、严格 JSON、原子 JSON 发布和复数 CSV 列。 |
| `output_paths.py` / `figure_output.py` | 按职责分类、旧布局读取及共享图保存；显示范围和归一化由图生成器决定。 |

旧 `s4_boundary`、`s4_plotting`、`s4_paths`、`s5_paths` 指向同一规范模块对象，
不复制业务实现。完整映射和验证见 [整合 Spec](../docs/spec/plan-execute_20261009_workflow_integration.md)。

## 三条主工作流的依赖

| 工作流 | 主要底层模块 |
| --- | --- |
| Unit-cell band | `geometry_utils.py`、`basis_utils.py`、`energy_recovery.py`、`band_connector.py`、`simulation_spatial/hexagon_unit_cell.py` |
| 1D strip | `geometry_utils.py`、`basis_utils.py`、`hex_lattice_utils.py`、`polygon_utils.py`、`simulation_utils.py`、`energy_recovery.py`、`finite_lattice_fourier.py`、`lattice_fourier_postprocess.py` |
| Finite cavity | `geometry_utils.py`、`finite_geometry.py`、`finite_geometry_comsol.py`、`hex_lattice_utils.py`、`polygon_utils.py`、`simulation_utils.py`、`mesh_constrcution.py`、`energy_recovery.py`、`finite_lattice_fourier.py`、`lattice_fourier_postprocess.py`、`farfield_fft.py` |
| Quarter finite cavity | 与 finite cavity 相同；`polygon_utils.py` 额外负责第一象限凸区域裁剪，`simulation_utils.py` 额外负责轴对称边界和四象限场重构。 |

```text
几何参数
  -> geometry/lattice 模块生成孔、cell 与有限边界
  -> simulation runner + mesh builder 构建并求解 COMSOL 模型
  -> 导出频率、Q、场、CSV/parquet/JSON
  -> band/Fourier/far-field 模块分析与绘图
```

## 几何、多边形与晶格

| 模块 | 功能与主要接口 | 定位 |
| --- | --- | --- |
| [`geometry_utils.py`](./geometry_utils.py) | 构造、旋转三角形，计算间距，生成包含六个三角孔的六边形 unit cell，并绘制几何预览。主要接口：`rotate_points()`、`create_inner_triangle()`、`create_single_sector()`、`create_hexagon_design()`、`visualize_hexagon_design()`。 | 三条主工作流共用的基础几何模块。 |
| [`polygon_utils.py`](./polygon_utils.py) | 多边形面积、最小宽度、偏移、点包含、线段相交、边界拼接、区域裁剪及裁剪结果验证。主要接口：`polygon_area()`、`polygon_min_width()`、`offset_polygon()`、`clip_polygon_to_outside_region()`、`clip_polygon_to_convex_region()`、`validate_clipped_hole()`。 | strip/finite 跨 cavity-cladding 边界孔裁剪，以及 quarter 第一象限裁剪的核心。 |
| [`hex_lattice_utils.py`](./hex_lattice_utils.py) | 六角晶格坐标、shell 编号、unit-cell 多边形、bulk/cladding 点、组合边界和连续层边界。主要接口：`LatticePoint`、`uv_to_xy()`、`lattice_points()`、`bulk_points()`、`boundary_from_lattice_cells()`、`retained_cladding_points()`。 | 当前 strip 与 finite cavity 的有限六角晶格生成器。 |
| [`finite_geometry.py`](./finite_geometry.py) | 不依赖 `mph` 的 finite geometry dispatcher；生成兼容接口的 hex full-cell plan 或 square half-cell plan，完成理想晶格归类、逐层 shift、footprint、validity region 和 analysis cells。 | Full finite 与 quarter finite 的 shape-aware 几何合同。 |
| [`finite_geometry_comsol.py`](./finite_geometry_comsol.py) | 把纯 `FiniteGeometryPlan` 编译为现有 `BoundarySpec`/`LayerSpec` 输入；full square 使用 `rectangle`，quarter square 使用 `quarter_rectangle`。 | 复用同一 COMSOL runner，不创建第二套 finite 求解器。 |
| [`finite_patch_lattice.py`](./finite_patch_lattice.py) | 生成 bulk、四侧边条带与四角块组成的有限 patch，支持正交或六角基。主要接口：`BulkSideCornerConfig`、`generate_lattice()`、`side_strips()`、`corner_blocks()`。 | 未被当前三入口调用，但作为今后新有限几何的可复用构建能力保留。 |

## 模式、能带与 Fourier 后处理

| 模块 | 功能与主要接口 | 定位 |
| --- | --- | --- |
| [`basis_utils.py`](./basis_utils.py) | 构造实正交 Fourier/旋转对称基，计算最小非零分量，并在 standard 与 Fourier basis 之间转换。 | 六孔参数化和模式成分的基础数学模块。 |
| [`energy_recovery.py`](./energy_recovery.py) | 对非结构化二维场插值和 Delaunay 积分，通过旋转自相关恢复 s/p/d/f 子空间能量，并把场映射为六维向量。 | unit-cell、strip 和 finite 模式识别的公共数值基础。 |
| [`band_connector.py`](./band_connector.py) | 读取各 k 点复数 Hz 场与模式成分，结合频率差和场重叠连接相邻 k 点的同一条 band。主要接口：`load_hz_center()`、`load_mode_composition()`、`assign_mode_candidates()`、`field_overlap()`、`BandConnector`。 | `run_band_pair.py` 的 band selection/tracking 模块。 |
| [`finite_lattice_fourier.py`](./finite_lattice_fourier.py) | 定义直接/倒易晶格、第一 Brillouin zone、有限六角 cyclic quotient，以及 square cavity 任意有限点集的直接 DFT。 | 有限晶格 Fourier 与 BZ 绘图的底层数学模块。 |
| [`lattice_fourier_postprocess.py`](./lattice_fourier_postprocess.py) | 从中心 Hz 场提取 strip rows 或 finite cells，执行晶格 Fourier 分解、绘制 k-space/模式图，并根据 Gamma 权重评分模式。 | strip 与 finite cavity 的主要模式后处理入口。 |
| [`farfield_fft.py`](./farfield_fft.py) | 读取 xy-air 复电场，在物理正方形窗口中补零并执行中心化 FFT/iFFT；计算 NA 收集、Gaussian FWHM 和发散角，输出 real/k-space/cut/fit 图与摘要。 | `run_finite.py` 的 valid-mode far-field 后处理。 |

## COMSOL 控制与 mesh

| 模块 | 功能与主要接口 | 定位 |
| --- | --- | --- |
| [`mesh_constrcution.py`](./mesh_constrcution.py) | 依据 named selections 构建 finite cavity 手动 mesh；把 `mesh_auto_size=1..9` 映射到最大增长率、曲率因子和狭窄区域分辨率，并配置材料尺寸与 PML sweep。主要接口：`resolve_manual_mesh_size_parameters()`、`construct_manual_mesh()`。 | 当前 finite cavity mesh builder。文件名保留既有拼写 `constrcution`。 |
| [`simulation_utils.py`](./simulation_utils.py) | 通用多层 COMSOL runner；以 `LayerSpec`、`BoundarySpec` 和 `SimulationConfig` 描述材料、边界和求解，支持 full/quarter hexagon 与 rectangle footprint，并负责 named selections、材料、top PML、周期/散射/PEC/PMC 边界、求解、日志和场导出。 | 当前 strip、完整 finite cavity 与 quarter finite cavity 的核心 COMSOL 接口。 |
| [`simulation_spatial/hexagon_unit_cell.py`](./simulation_spatial/hexagon_unit_cell.py) | 六角 unit cell 的 Floquet runner，使用位置 named selections 构建 slab、air、PML 和周期边界，求解 `(kx, ky)` 并导出频率、中心场和 polarization。 | 当前 unit-cell 主入口的首选 runner。 |
| [`simulation_spatial/square_unit_cell.py`](./simulation_spatial/square_unit_cell.py) | 针对正方形晶格的 Floquet unit-cell runner，具有与 hexagon runner 相近的构建、求解和导出接口。 | 当前三入口未调用；作为新正方形 unit-cell 计算的本机控制模板保留。 |
| [`simulation_spatial/rectangle_finite_size.py`](./simulation_spatial/rectangle_finite_size.py) | 构建 `Lx x Ly` 矩形有限 slab、任意多边形孔、air 和 PML，可只构建几何或完整求解并导出频率和场。 | 当前三入口未调用；作为新矩形有限结构计算的本机控制模板保留。 |

## 使用与维护约定

- 新需求优先复用这些模块；可复用数值逻辑放在本目录，具体参数和 orchestration
  放在 `scripts/` 对应分类目录。
- 本目录是标准 Python package。包内依赖使用相对导入（例如
  `from .energy_recovery import ...`），工作流和外部调用使用完整包名（例如
  `from comsol_workflow.geometry_utils import ...`）。新增模块应遵循同一规则，不要
  恢复 flat import 或把 `comsol_workflow/` 子目录直接加入 `sys.path`。
- 包根 `__init__.py` 和 `simulation_spatial/__init__.py` 保持轻量，不主动导入
  `mph`、runner 或其他重依赖；这使几何和分析模块可以独立复用，也便于继续增加
  新的几何、mesh、runner 与后处理模块。
- 创建 `SimulationRun`、调用 mesh builder 或实际求解会启动本机 COMSOL；普通
  单元测试应避免无意触发这些操作。
- 后处理优先读取已保存的 CSV、parquet、JSON 和场数据，不为重绘重复 eigensolve。
