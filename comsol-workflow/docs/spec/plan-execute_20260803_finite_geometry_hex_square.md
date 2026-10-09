# Finite hex/square 几何构建切换：实施前方案与执行记录

日期：2026-08-03
状态：核心功能已实施；full/quarter geometry+mesh 已验收，未运行 eigensolve

## 1. 背景与目标

当前 full finite 与 finite-quarter 只支持六角有限晶格：cavity 使用六角 shell，
cladding 继续沿六角 shell 向外扩展，最终计算 footprint 也是六角形。用户已经审核通过
`cav5_clad5_geometry_concept.png` 中的方形 fragment 结构，要求把它作为与现有六角结构
并列、互斥且可切换的第二种 finite 几何。

本次实施目标：

1. 在 `scripts/parameter.json` 增加 finite 几何选择，允许 `hex` 与 `square`；
2. `hex` 保持当前几何、shift、后处理、输出命名和 checkpoint 兼容；
3. `square` 严格复现已批准的 Group 1/Group 2 half-cell 装配合同；
4. square 同时兼容现有所有 cladding shift 输入、几何和逐层 profile；
5. full finite 与 finite-quarter 读取同一选择并复用同一个完整几何计划；
6. far-field、finite lattice Fourier、mode validity、resume 和几何诊断均正确理解实际
   footprint，不得把 square 结果继续按 hex 处理；
7. 不修改、移动、覆盖或删除现有正式结果，新的 square 输出使用非冲突目录。

本次只改变 finite/full-quarter。Unit-cell 不读取 finite 结构层数或 finite footprint；
strip 继续使用现有 1D cut 几何和 shift 扫描，不由该选项改变。

## 2. 已批准的 square 几何合同

### 2.1 晶格与组织单元

保持现有三角晶格周期 `A`，定义：

```text
a1 = (A, 0)
a2 = (A/2, sqrt(3) A/2)
h  = sqrt(3) A/2
R(m,n) = m a1 + n a2
       = (A(m+n/2), n h)
```

每个理想晶格点 `R=(x,y)` 对应一个仅用于归属判断的 brick 组织单元：

```text
[x-A/2, x+A/2] x [y-h/2, y+h/2]
```

这个矩形不是 COMSOL 物理域，不得进入最终 geometry。它只用于判断哪些 unit-cell
fragment 属于 cavity、cladding 或结构外部。

### 2.2 Group 1/Group 2

沿 `x=x_center` 把组织单元分成两个互不重叠的 half-cell fragment：

```text
Group 1 support = [x-A/2, x]     x [y-h/2, y+h/2]
Group 2 support = [x, x+A/2]     x [y-h/2, y+h/2]
```

以当前 `create_hexagon_design()` 返回的六个孔索引为准：

```text
Group 1 = hole indices {2, 3, 4}  # 左半边三孔
Group 2 = hole indices {0, 1, 5}  # 右半边三孔
```

索引分组必须由单元几何中心验证，而不是只依赖硬编码：测试同时断言三孔中心的 x
符号和精确索引，防止上游孔排序改变后静默生成错误结构。

### 2.3 cavN 与 cladM 的范围

令 `Nc=cavity_layers`、`Nd=cladding_layers`。

方形 cavity fragment support 的闭合范围为：

```text
|x| <= Nc A
|y| <= (Nc + 1/2) h
```

完整光子晶体结构的 fragment support 范围为：

```text
|x| <= (Nc + Nd) A
|y| <= (Nc + Nd + 1/2) h
```

fragment 归属规则按以下顺序执行：

1. support 完整落在 cavity 范围内：`region=cavity`；
2. support 完整落在结构范围内但不属于 cavity：`region=cladding`；
3. 其他 fragment：丢弃，不创建孔。

支持范围使用统一数值容差，边界判断采用包含端点的闭区间。每个 `(m,n,group_id)`
只生成一次，因此闭区间不会造成重复实体。

左右 cavity 边界经过 `x=+/-Nc A` 上的晶格点中心。在右边界，Group 1 可以属于
cavity、Group 2 属于 cladding；左边界则相反。上下边界位于相邻晶格行之间，保留完整
行，不需要再拆上下 Group。这与用户批准的 `cav5_clad5` 示意图一致。

### 2.4 cav5_clad5 固定验收样例

使用当前 `A=0.82 um`、当前 cavity/cladding 单元参数、shift=0 时，批准样例必须得到：

```text
total fragments             = 840
cavity fragments            = 220
cladding fragments          = 620
parent lattice points       = 431
full cavity parent cells    = 105
mixed cavity/clad parents   = 10
outer half-cell parents     = 22
hole polygons               = 2520
```

其中每个 fragment 含三个孔。`105` 个 parent 的两个 fragment 都属于 cavity，可作为完整
cavity unit cell 进入 unit-cell basis/Fourier 分析；10 个 mixed parent 参与 COMSOL 几何，
但不能伪装成完整 cavity unit cell。

### 2.5 外部物理 footprint

square 光子晶体结构范围外，在 Cartesian x/y 四侧各增加 `A/2` 的无孔均匀介质平板：

```text
physical half width  = (Nc + Nd) A + A/2
physical half height = (Nc + Nd + 1/2) h + A/2
```

外矩形四个侧面施加 SBC。cavity 边界和 cladding 外边界都只是逻辑/诊断边界，不得成为
COMSOL 物理边界。首版 square 保持批准示意图的 `lateral PML=off`；现有顶部空气层和
z 向 PML 保持不变。以后若增加横向 PML，必须另立方案处理侧向 PML domain partition
和 mesh，不能把它混入本次重构。

## 3. parameter.json 合同

新增顶层字段：

```json
"finite_geometry": "square"
```

只允许精确小写字符串：

```text
hex
square
```

含义：

| 值 | full finite | finite-quarter | strip | unit-cell |
| --- | --- | --- | --- | --- |
| `hex` | 当前六角完整结构 | 当前 quarter hex | 不影响 | 不影响 |
| `square` | 批准的 fragment 方形结构 | quarter rectangle | 不影响 | 不影响 |

兼容规则：

- 当前活动 `parameter.json` 在实施时显式写入 `"finite_geometry": "square"`；
- 旧参数文件缺失该字段时回退到 `hex`，保证历史调用和 fixture 可读取；
- 未知值、非字符串、大小写变体都明确报错；
- 该字段写入 `SharedParameters.to_metadata()`，但不加入 unit-cell identity，改变它不得
  阻止或触发 unit-cell `p2@Gamma` 中心频率回写；
- strip 的系列命名和运行行为不得因该字段变化。

`finite_geometry` 与现有 `cladding_shift_geometry` 是两个独立维度，禁止混淆：

```text
finite_geometry          = hex | square
cladding_shift_geometry = groups | ellipse
```

四种组合全部必须可运行：

```text
hex    + groups
hex    + ellipse
square + groups
square + ellipse
```

## 4. 目标代码架构

### 4.1 纯 Python 几何计划

新增 `comsol_workflow/finite_geometry.py`，只依赖 NumPy 和现有几何工具，不导入 `mph`
或 JPype。建议数据模型：

```text
FiniteGeometrySpec
  shape: hex | square
  period
  cavity_layers
  cladding_layers
  empty_buffer

CellFragmentTemplate
  cell_kind: cavity | cladding
  group_id: 1 | 2
  hole_polygons
  local_support

FragmentPlacement
  m, n
  parent ideal center
  group_id
  region
  ideal support
  square layer / hex shell layer
  outside axes
  shift vector
  final hole polygons

FiniteFootprint
  logical cavity boundary
  logical structured boundary
  physical boundary
  COMSOL BoundarySpec parameters

FiniteGeometryPlan
  spec
  placements
  full cavity analysis cells
  metadata
  footprint
```

`run_finite.py` 不再自行枚举 lattice point、孔 polygon 和外边界。它只负责：读取共享
参数、选择 geometry planner、调用 COMSOL runner、输出和后处理。

### 4.2 shape dispatcher

使用单一入口：

```text
build_finite_geometry_plan(...)
  if shape == hex:
      build_hex_geometry_plan(...)
  elif shape == square:
      build_square_geometry_plan(...)
```

两种 geometry plan 具有相同的下游接口，但 shape-specific 数学实现不互相调用。
选择是严格互斥的；一次 case 不得同时生成 hex 与 square fragment。

### 4.3 hex 适配原则

hex planner 第一阶段只是把当前稳定实现适配为 `FiniteGeometryPlan`：

- `bulk_points()`、`lattice_points()`、六角 shell 和 `inner_boundary` 计算不变；
- 当前 `normal_fan_feature()`、corner 补偿、ellipse 径向位移和逐层 profile 不变；
- 当前 hole polygon 数值和顺序保持；
- 当前 full/quarter 六角 footprint、far-field hex mask 和 hex cyclic Fourier 保持；
- 当前 hex 输出系列目录不增加新 suffix；
- 缺失 `finite_geometry` 的历史 config 按 hex 归一化后允许 resume；
- 用代表性 fixture 比较迁移前后 polygon 顶点、point、shift、metadata 和 boundary，要求
  数值逐项一致。

## 5. square cladding shift 合同

### 5.1 shift 执行顺序

square 必须按以下固定顺序处理：

1. 在理想晶格上生成 parent point；
2. 按理想 support 判断 Group 1/2 的 cavity/cladding/丢弃归属；
3. 为 cladding fragment 计算 square layer；
4. 根据 `cladding_shift_geometry`、x/y factor 和 profile 计算 shift；
5. 将 shift 施加到该 fragment 的三个孔；
6. 不因 shift 后的位置重新分类、裁剪或累计到下一层。

因此逐层 shift 始终是相对理想晶格的绝对位移。两个都属于 cladding 的同一 parent
Group 1/2 使用同一个 layer 和同一个 shift，保持完整 cladding cell 刚性；mixed parent
只移动 cladding fragment，cavity fragment 保持理想位置。

### 5.2 square layer

parent 理想中心为 `(x,y)`，定义：

```text
ux = |x| / A
uy = |y| / h

Lx = ceil(max(0, ux - Nc))
Ly = ceil(max(0, uy - Nc))
L  = max(1, Lx, Ly)
```

计算时在 `ceil` 前使用统一 epsilon，避免理论整数因浮点误差落入下一层。所有保留的
cladding fragment 必须满足 `1 <= L <= Nd`；超界立即报错，不静默 clamp。

该定义按 parent 分层，因此完整 cladding cell 的两个 fragment 不会因 half support
中心不同而获得不同 profile envelope。

### 5.3 profile 兼容

现有 `CladdingShiftProfile.envelope(L)` 原样复用：

- `uniform`；
- 当前支持的非 uniform profile（包括 `tanh_power`）；
- x/y 使用同一个 layer envelope；
- profile 只调制幅值，不改变 fragment 归属和方向；
- `scale_layers`、`power` 等参数合同继续由共享 profile 模块验证。

### 5.4 ellipse 兼容

square 的 ellipse 模式保留当前公式：

```text
theta = atan2(y, x)
f(theta) = fx cos(theta)^2 + fy sin(theta)^2
delta = -A g(L) f(theta) (x,y) / sqrt(x^2+y^2)
```

要求：

- 使用 parent 理想中心而不是 shift 后中心；
- 同一完整 cladding parent 的两个 fragment 共用 delta；
- mixed parent 只对 cladding fragment 应用 delta；
- `fx=0`、`fy=0` 或两者同时为 0 均合法；
- x/y 镜像后位移严格满足 C2v 对称，保证 quarter 可用；
- 不继承 hex groups 的 corner `2/sqrt(3)` 补偿。

### 5.5 groups 兼容

hex 继续使用当前 12 个 normal-fan side/corner 分组。square 新增独立的
`square_axis_components_v1`，不得调用 hex 的 `normal_fan_feature()`。

对 cladding fragment，利用其理想 support 判断 cavity 外侧轴：

```text
outside_x = support 在 x 方向超出 cavity support 范围
outside_y = support 在 y 方向超出 cavity support 范围
```

位移分量定义：

```text
delta_x = -sign(x_anchor) A fx g(L)   if outside_x else 0
delta_y = -sign(y_anchor) A fy g(L)   if outside_y else 0
delta   = (delta_x, delta_y)
```

其中 anchor 优先使用 parent center；若对应坐标理论上为 0，则使用 fragment support
中心确定符号。side fragment 只有一个非零分量；corner fragment 同时使用 x/y 分量。

square corner 不使用 hex 的 `2/sqrt(3)` 补偿。其角点位移大小自然为
`A g(L) sqrt(fx^2+fy^2)`。该差异必须通过 formula version 和 metadata 明确记录，禁止把
square groups checkpoint 当成 hex groups checkpoint 恢复。

### 5.6 shift 输入矩阵

以下输入在 square 下必须全部覆盖：

- ratio scan：`cladding_shift_factors + cladding_y_over_x_shift_ratio`；
- explicit pair：`cladding_x_shift_factor + cladding_y_shift_factor`；
- `groups` 与 `ellipse`；
- `uniform` 与所有现有非 uniform profile；
- `fx=fy`、`fx!=fy`、单轴为 0、双轴为 0；
- full finite 与 finite-quarter。

## 6. COMSOL geometry 编译

### 6.1 reference 与正式 compiler

纯 Python plan 同时提供：

1. flatten reference：把每个 placement 展开为最终 polygon，用于绘图、测试和与旧
   `LayerSpec.holes` backend 做等价比较；
2. template compiler input：保留 cavity/cladding、Group 1/2 template 和 translation，
   用于减少 COMSOL 几何树中的重复定义。

正式 square COMSOL builder 的目标流程：

```text
4 source templates
  cavity Group 1
  cavity Group 2
  cladding Group 1
  cladding Group 2
    -> placement Copy/Move 或可证明等价的 Array cohort
    -> hole Union, internal boundaries off
    -> slab footprint Difference holes
    -> air/top-PML extrusion
    -> final geometry union
```

由于 ellipse 和非 uniform profile 可能让每个 parent translation 不同，compiler 只在
translation 完全相同且晶格规则可表达时使用 Array；其他 placement 使用 template 的
Copy/Move，不得为了减少 feature 数近似 shift。

正式实现前先用 `cav1_clad1` geometry-only 模型验证 COMSOL 6.3 中 Array/Copy、Union、
`intbnd=off` 和 generated selection 的实际 Java API。验证代码不进入正式结果目录。

### 6.2 不创建逻辑矩形物理边界

cavity 和 structured rectangle 只存在于 Python metadata/diagnostic plot。COMSOL 中只
创建最终外部 physical footprint 和孔。不得把 cavity rectangle 与 cladding ring 分别
拉伸后再保留内部界面；最终 silicon slab 是统一 slab layer 减去 hole domains。

### 6.3 named selections

继续遵循位置/named selection 原则，至少提供：

```text
sel_layer_<i>_mat_dom
sel_layer_<i>_hole_doms
sel_physical_dom
sel_top_pml_dom
sel_pml_dom                 # top PML 的兼容 union
sel_bottom
sel_top_face
sel_side_x
sel_side_y
```

square side wall 用 physical rectangle 四侧的 Ball seed 加 `groupcontang=on`，再组合为
x/y 两组 SBC。禁止依赖 domain/boundary 固定 ID。

首版 square 只有 top PML，因此现有 `sel_pml_dom` 和 PML sweep 语义保持；命名中增加
`sel_top_pml_dom` 是为了以后扩展 lateral PML 时不混淆 mesh 方向。

### 6.4 通用 runner 接口

`SimulationRun` 保留现有 explicit polygon API供 unit-cell/strip 使用，并增加有限结构
geometry plan/backend 入口。材料、physics、study、求解和导出仍由同一 runner 负责，
避免再复制一个 finite-only 求解器。

`BoundarySpec` 扩展：

```text
hexagon_boundary
quarter_hexagon_boundary
rectangle
quarter_rectangle
```

full square 使用 rectangle；quarter square 使用 quarter_rectangle。

## 7. finite-quarter

### 7.1 共享 full plan

finite-quarter 必须先调用与 full finite 完全相同的 `build_finite_geometry_plan()`：

```text
parameter.json
  -> full hex/square geometry plan
  -> mirror-symmetry validation
  -> first-quadrant clipping
  -> quarter COMSOL geometry
```

禁止 quarter 单独重算 lattice、fragment region、layer 或 shift。

### 7.2 square quarter footprint

full physical rectangle 为：

```text
[-Lx/2,Lx/2] x [-Ly/2,Ly/2]
```

quarter rectangle 为：

```text
[0,Lx/2] x [0,Ly/2]
```

边界条件：

- `x=0` 与 `y=0`：继续按 symmetry case 使用 PEC/PMC；
- `x=Lx/2` 与 `y=Ly/2`：SBC；
- bottom/top：继续复用当前 PMC/PEC 与 top SBC/PML 合同。

跨越 x/y 对称轴的孔 polygon 继续使用精确 convex clipping；不得仅按中心点删除。现有
四象限 scalar/vector field 重建公式不依赖外 footprint，可继续复用。

### 7.3 对称性门槛

每种 square shift 组合在构建 quarter 前必须通过：

- hole polygon 关于 x/y 轴的镜像集合一致；
- fragment region、layer 和 shift vector 镜像一致；
- 单轴 factor 为 0 时仍保持对应镜像；
- 非对称输入或未来不对称 profile 明确拒绝 quarter，不静默运行。

## 8. 后处理迁移

### 8.1 footprint 与 far-field

新增通用 footprint metadata 和 `contains_xy()`：

```text
HexFootprint.contains_xy(x,y)
RectangleFootprint.contains_xy(x,y)
```

far-field 继续使用等间距 square FFT grid。采样窗口边长取 physical footprint 的
`max(Lx,Ly)`，窗口中 rectangle 外的点置 0；hex 保持当前 hex mask。输出 metadata 从
只写 `finite_hex_a_um` 改为同时写通用字段：

```text
finite_geometry
footprint_shape
footprint_Lx_um
footprint_Ly_um
sampling_window_um
domain_mask_kind
```

hex 兼容字段在 hex case 中继续保留，square 不伪造 `finite_hex_a`。

quarter 重建 air field 使用同一个 full footprint mask，保证 full/quarter far-field 输入
定义一致。

### 8.2 mode validity

mode validity 不再隐式调用 hex `bulk_boundary()`。geometry plan 显式提供 validity
polygon：

- hex：保持当前 bulk 加半 cladding 的六角边界；
- square：使用完整 cavity parent cell core，按已批准矩形逻辑向内排除 mixed half-cell
  和指定 cladding margin。

具体 margin 继续使用当前 `floor(cladding_layers/2)` 语义，但由 shape-specific planner
生成实际 polygon并写入 metadata。

### 8.3 finite lattice Fourier

hex 继续使用当前 finite cyclic quotient：

```text
N = 1 + 3L(L+1)
hex_spiral_indices
hex_cyclic_indices
```

square 不满足该群结构，禁止调用上述函数。新增 arbitrary finite lattice DFT：

```text
F(k) = sum_j c_j exp(-i k dot R_j)
```

其中 `R_j` 只取两个 fragment 都属于 cavity 的完整 parent cell；mixed boundary cell 参与
COMSOL 求解，但不进入完整 unit-cell basis decomposition。k-space 使用现有 reciprocal
basis 和第一 BZ定义，输出保持 `P(k)`、top-k peaks、p-subspace decomposition 和标准图，
但 index 审计明确写 `index_kind=arbitrary_square_cavity_v1`。

square 与 hex 的 Fourier summary 必须记录不同 transform kind，分析代码不得把两者合并
成同一趋势而不分组。

## 9. 输出、命名与 resume

### 9.1 系列目录

为了保留既有结果路径：

- `hex` 使用当前系列目录名，不增加 suffix；
- `square` 在 cell/mesh label 后增加稳定 `_square` suffix；
- profile 和 ellipse suffix 继续追加在 shape suffix 之后；
- shift case 名保持 `shiftxX.XXX_shiftyY.YYY`；
- quarter 继续使用 `finite_quarter_...` prefix。

示例：

```text
finite_5-5_cav(... )_clad(... )_mesh5_square_ellipse-shift/
  shiftx0.050_shifty0.080/

finite_quarter_5-5_cav(... )_clad(... )_mesh5_square_ellipse-shift/
  shiftx0.050_shifty0.080/
```

实际实现不在名称中保留示例空格。

### 9.2 config metadata

square config 至少记录：

```text
finite_geometry = square
finite_geometry_schema = square_fragments_v1
lattice_vectors_um
organizational_cell_size_um
group_hole_indices
cavity_fragment_bounds_um
structured_fragment_bounds_um
physical_footprint_um
empty_buffer_um
fragment_counts
parent_cell_counts
mixed_parent_count
square_layer_formula
cladding_shift_formula
analysis_cell_policy
```

hex config 写 `finite_geometry=hex` 和对应 hex schema，同时保留所有现有字段。

### 9.3 resume 合同

resume 比较项增加：

```text
finite_geometry
finite_geometry_schema
footprint shape/Lx/Ly
fragment/parent counts
group_hole_indices
square layer formula version
shift formula version
```

允许：

- 缺失 `finite_geometry` 的历史 checkpoint 归一化为 hex；
- 新 hex 配置恢复匹配的旧 hex checkpoint。

禁止：

- hex 与 square 互相恢复；
- square groups 与 square ellipse 互相恢复；
- square formula version 不同仍恢复；
- full 与 quarter checkpoint 互用。

## 10. 修改范围

预计新增：

- `comsol_workflow/finite_geometry.py`：纯 geometry plan、hex/square dispatcher；
- `comsol_workflow/finite_geometry_comsol.py`：COMSOL template compiler；
- `tests/test_finite_geometry.py`：shape 与 fragment 核心合同。

预计修改：

- `scripts/parameter.json`；
- `scripts/run_main/parameter_config.py`；
- `scripts/run_main/run_finite.py`；
- `scripts/run_main/run_finite_quarter.py`；
- `comsol_workflow/simulation_utils.py`；
- `comsol_workflow/mesh_constrcution.py`；
- `comsol_workflow/finite_lattice_fourier.py`；
- `comsol_workflow/lattice_fourier_postprocess.py`；
- 相关 far-field、resume、输出和 quarter 测试；
- 根 README、scripts README、comsol_workflow README、`AGENTS.md` 和 `MEMORY.md`。

明确不修改：

- unit-cell eigensolver 基准及 unit-cell 几何；
- strip cut 几何、strip shift 扫描和 strip 输出布局；
- 已有正式结果目录和文件；
- 当前参数中的 cell、mesh、频率与 shift 数值，除新增
  `"finite_geometry": "square"`。

## 11. 分阶段执行

### 阶段 A：参数与纯几何

1. 扩展 `SharedParameters`，验证/序列化 `finite_geometry`；
2. active parameter 写 `square`，旧 fixture 缺失时回退 hex；
3. 实现 geometry dataclass 和 dispatcher；
4. 把当前 hex 构建适配为 plan 并做逐点等价测试；
5. 实现 square fragment、parent、layer、footprint 和 flatten reference；
6. 用 cav5_clad5 固定计数和批准 PNG 锁定合同；
7. 此阶段不导入 mph、不启动 COMSOL。

### 阶段 B：所有 square shift

1. 把 shift 计算从 `run_finite.py` 移入 shape-aware geometry plan；
2. 保留 hex groups/ellipse bit-for-bit；
3. 实现 square groups component formula；
4. 复用 ellipse 与 profile，保证完整 cladding parent 刚性；
5. 测试 ratio/explicit、groups/ellipse、uniform/nonuniform、零轴；
6. 测试 full geometry x/y 镜像和 quarter 前置校验。

### 阶段 C：COMSOL geometry-only

1. 扩展 BoundarySpec 和 finite geometry backend；
2. 用 cav1_clad1 探针确认 Copy/Move/Array/Union API；
3. 用 cav2_clad2 构建 geometry-only full square；
4. 检查 cavity/clad 逻辑矩形未成为物理 boundary；
5. 检查 slab、hole、air、top PML 与四侧 selection；
6. 构建 mesh 但不 eigensolve，记录 feature/domain/boundary/mesh 数；
7. 导出 COMSOL 实际 geometry PNG人工检查。

任何 `SimulationRun` 或 geometry build 都可能启动 COMSOL。执行阶段 C 前重新核对本机
COMSOL 6.3、license、输出位置和用户授权，不把它当普通 smoke test。

### 阶段 D：full/quarter 工作流

1. full finite 消费统一 geometry plan；
2. quarter 从 full plan clip，增加 quarter_rectangle selections；
3. 保持材料、physics、study、checkpoint 和输出布局；
4. 更新 geometry diagnostic，让 square 图复现批准结构；
5. 更新命名和 resume 比较。

### 阶段 E：后处理

1. 引入通用 footprint mask；
2. 迁移 full/quarter air-field 和 far-field metadata；
3. 迁移 validity boundary；
4. 实现 square arbitrary finite lattice DFT；
5. 保持 hex cyclic Fourier 输出不变；
6. 更新 trend/summary 对 geometry shape 的分组和拒绝规则。

### 阶段 F：正式验证与文档

1. 运行全部纯 Python 聚焦测试；
2. 运行相关全仓回归；
3. 运行 `py_compile` 与 `git diff --check`；
4. 更新 README、AGENTS、Memory 和本文件执行记录；
5. 用户审核 geometry-only 后，才决定是否运行 small eigensolve；
6. full cav5_clad5 或正式 cav10_clad10 eigensolve 另行 preflight，不在本重构验证中自动
   启动。

## 12. 验证矩阵

### 12.1 配置

- `hex`、`square` 合法；
- 缺失字段回退 hex；
- 大写、未知值、非字符串拒绝；
- finite/quarter 读取一致；
- strip/unit-cell identity 和目录不变；
- square 系列目录不与 hex 冲突。

### 12.2 square 纯几何

- cav5_clad5 固定计数；
- Group 1/2 孔索引和 support；
- 无重复 `(m,n,group_id)`；
- 每个 fragment 恰好 cavity/cladding/丢弃之一；
- 边界 mixed parent 精确为批准结构；
- 所有孔位于对应 shifted fragment，shift=0 时等于模板平移；
- x/y 镜像集合一致；
- physical footprint 相对 structured bounds 四侧均为 `A/2`。

### 12.3 shift

- hex 原有测试全部通过；
- square layer 在 1..Nd 且 parent 内一致；
- groups side/corner 分量；
- ellipse 主轴、45 度和零轴；
- profile 不累计；
- mixed parent 只移动 cladding fragment；
- 完整 cladding parent 两个 fragment 同位移；
- quarter mirror symmetry。

### 12.4 COMSOL

- logical rectangle 不产生内部物理 face；
- named selections 非空、无错误重叠并覆盖目标实体；
- silicon selection 等于 slab extent 减 hole domains；
- top PML selection 不包含 physical domain；
- 四侧 SBC、顶部 SBC、底部 TE/PMC 正确；
- manual mesh 使用 material/physical/top-PML selections；
- geometry-only checkpoint 不含 mesh/solution；
- explicit flatten reference 与 template compiler 最终实体等价。

### 12.5 full/quarter 与后处理

- square full/quarter 频率在小模型对称 case 中一致；
- reconstructed field 的坐标和 parity 正确；
- far-field rectangle mask 正确置零；
- square 不出现 `finite_hex_a` 伪字段；
- square Fourier 只使用 full cavity parents；
- hex cyclic Fourier 数值回归不变；
- resume 跨 shape 明确拒绝。

## 13. 回滚方案

最低风险回滚只需把活动参数改为：

```json
"finite_geometry": "hex"
```

hex branch、旧输出路径和旧 checkpoint 合同保持。若 square 实现本身需要代码回滚，可删除
square planner/compiler、quarter_rectangle 和 arbitrary DFT 分支，不修改或迁移任何旧
结果。square 新目录保留为独立结果，不覆盖到 hex。

实现期间每个阶段保持独立可测试；若阶段 C 的 COMSOL template compiler 不稳定，可暂时
保留纯 Python flatten reference 作为 geometry-only 对照，但不得把它标记为最终完成，
也不得绕过逻辑矩形无物理边界的验收要求。

## 14. 风险与处理

1. **COMSOL feature 数仍过多**：template + cohort compiler；先用小模型测 feature 和
   build 时间，不直接运行 cav10_clad10。
2. **shift 后孔跨越逻辑边界**：归属只按 ideal support；metadata 明确，不重新裁剪。
3. **square groups corner 语义与 hex 不同**：独立 formula version，不共享 checkpoint。
4. **quarter clipping 产生过窄 polygon**：复用 polygon 宽度/面积验证，小于数值门槛明确
   报错或按正式既有过滤合同处理，不静默保留退化实体。
5. **far-field FFT window 改变采样**：统一 square sampling window并记录 mask/spacing；
   full/quarter 使用同一合同。
6. **旧分析误读新结构**：所有 summary/config 写 shape；聚合时 shape 必须分组。
7. **工作树已有修改**：只编辑本任务文件，不回滚或覆盖用户当前改动；实施前再次检查
   `git status`。

## 15. 执行记录

当前状态：核心功能已实施，纯 Python、full COMSOL geometry+mesh 与 quarter COMSOL
geometry+mesh 验收通过；正式大尺寸 eigensolve 未启动。

### 实际改动

- `scripts/parameter.json` 新增 `"finite_geometry": "square"`；共享参数严格验证
  `hex|square`，缺失回退 `hex`，且不进入 unit-cell identity。
- 新增 `comsol_workflow/finite_geometry.py`，实现 `FiniteGeometryPlan`、hex/square
  dispatcher、square Group 1/2 fragment、理想晶格归属、parent layer、groups/ellipse
  shift、rectangle footprint、validity boundary 和完整 cavity analysis cells。
- 新增 `comsol_workflow/finite_geometry_comsol.py`，把 plan 编译到现有
  `BoundarySpec`/`LayerSpec` named-selection backend；full square 使用 `rectangle`，
  quarter square 使用 `quarter_rectangle`，未创建 cavity/structured 逻辑矩形实体。
- `run_finite.py` 与 `run_finite_quarter.py` 已按 shape 选择边界、几何图、resume
  metadata、far-field footprint mask 和 Fourier transform；square 系列增加 `_square`，
  hex 系列名称保持不变。
- `simulation_utils.py` 新增 `quarter_rectangle`、square quarter 对称轴/外侧 wall seed、
  `sel_top_pml_dom` 与兼容 `sel_pml_dom` union。
- `finite_lattice_fourier.py` 与 `lattice_fourier_postprocess.py` 新增
  `arbitrary_square_cavity_v1` 直接 DFT；只消费两个 fragment 均为 cavity 的 parent，
  Gamma 明确为 index 0，hex 继续使用 cyclic quotient。
- full/quarter air-field 导出写通用 footprint metadata；square 使用 rectangle mask，
  不伪造 `finite_hex_a`。
- 更新根 README、scripts README、comsol_workflow README、`AGENTS.md` 与
  `MEMORY.md`。

### 测试结果

- `uv run python -m pytest -q`：`385 passed`，3 个既有
  `PytestReturnNotNoneWarning`，无失败。
- 新增/扩展测试覆盖：配置缺失回退与非法值、cav5_clad5 固定计数、fragment 去重、
  half-period buffer、groups/ellipse、uniform/tanh profile、完整 cladding parent 刚性、
  mixed parent、rectangle/quarter_rectangle compiler、square far-field mask、任意点 DFT、
  COMSOL boundary seed。
- 主入口及新增/修改模块 `py_compile` 通过。
- `scripts/parameter.json` 标准 JSON 读取检查通过；`git diff --check` 通过，仅有仓库
  现有 LF/CRLF 转换提示。

### COMSOL 运行与结果文件

实际启动 COMSOL 6.3 两次，仅执行小模型 geometry+manual mesh，未运行 eigensolve：

1. full `cav1_clad1 square+ellipse`：40 fragments、120 holes、rectangle
   `4.1 x 4.3707041555 um`；8 个关键 named selections 存在；mesh tags 为
   `size`、`size_layer0`、`ftet1`、`swe1`。
2. quarter 同一结构：32 个裁剪后 holes、`quarter_rectangle`；10 个关键 named
   selections 存在；PEC/PMC symmetry axes、outer rectangle SBC seeds 与 manual mesh
   均构建成功。

两次探针生成的临时脚本、MPH、progress log 和目录均已删除；未创建、移动、覆盖或删除
任何正式计算结果。批准样例仍为：

```text
scripts/.out/finite_cavity/cav5_clad5_geometry_concept.png
```

### 遗留问题与最终状态

功能合同已经完成：hex/square 可切换，四种 shape/shift 组合在纯几何层可运行，square
full/quarter COMSOL geometry+mesh 已通过，far-field、Fourier、命名和 resume 已按 shape
分离。

本次没有运行任何 square eigensolve，也没有对正式 `cav10_clad10` 做耗时/内存基准。
COMSOL compiler 当前使用已经过验证的 flattened polygon `LayerSpec` backend，而没有采用
原方案设想的 GeometryPart/Copy cohort 压缩；这不改变几何和 selection 合同，但
`cav10_clad10` 的 9840 个孔可能带来 geometry build 性能风险。正式大尺寸计算前应先运行
`cav2_clad2` 小本征对照，并单独完成 active case 的资源 preflight。

## 16. 2026-08-03 square overview cell rendering follow-up

### Plan

- Keep the square cavity classification and COMSOL geometry unchanged.
- Replace the rectangular organizational supports in
  `full_lattice_modulation_vectors.png` and
  `full_lattice_ellipse_shift.png` with the actual hexagonal parent-cell
  polygons from the triangular lattice.
- Merge two compatible half-cell placements into one complete hexagonal cell.
  Keep mixed cavity/cladding parents and outer half-cell parents split along
  the actual hexagon's left/right halves.
- Emit one modulation arrow for a rigidly shifted complete cladding parent and
  one arrow for each independently shifted cladding half-cell.
- Redraw only the two named PNG files in the existing square result directory;
  do not run COMSOL or overwrite other formal outputs.
- Verify polygon topology and merge counts with focused tests, inspect both
  exported PNG files, run `py_compile`, and finish with `git diff --check`.

### Rollback

The plotting helper and square-only branch are isolated from geometry
construction. Reverting this follow-up restores the old rectangular overview
without changing any COMSOL model, field data, or numerical result.

### Execution record

- Added `square_fragment_cell_polygon()` and `square_cell_plot_records()` to
  the pure finite-geometry module.
- Complete cavity/cladding parents now render as six-vertex hexagons. Mixed
  parents and outer half-cell parents render as four-vertex left/right halves
  of the same hexagon.
- The square overview generator consumes those records. A rigid complete
  cladding parent produces one modulation arrow instead of two duplicate
  half-cell arrows.
- Rebuilt the existing `cav20_clad20`, `square+ellipse`,
  `shiftx0.050_shifty0.080` plan from its saved `99_config/config.json` and
  redrew only the two requested PNG files in `00_model`.
- Visual inspection at the exported resolution confirmed the staggered
  triangular lattice, hexagonal cells, half-cell boundaries, equal x/y aspect,
  intact physical footprint, readable legend, and nonblank output.
- No COMSOL process, geometry build, mesh, or eigensolve was started.

### Verification

- `uv run python -m pytest -q tests/test_finite_geometry.py tests/test_finite_cavity_field_exports.py`:
  `22 passed`.
- `python -m py_compile` passed for `finite_geometry.py`, `run_finite.py`, and
  `test_finite_geometry.py`.
- `git diff --check` passed; only the repository's existing LF/CRLF warnings
  were reported.
