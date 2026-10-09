# COMSOL finite geometry 到微纳加工 DXF 版图：实施前方案与执行记录

日期：2026-08-04
状态：v1 工程版实现与验证完成；production 发布等待加工参数标定和外部 CAD 验收

## 1. 背景与目标

当前仓库已经能够从共享参数构建 photonic-crystal unit cell、strip 和 finite cavity，
并通过 COMSOL 完成本征频率、场、Fourier 与 far-field 分析。最后缺少的环节是把已经
选定的 finite cavity 设计转换为可交付给微纳加工流程的二维蚀刻版图。

本任务的目标不是把 COMSOL 输入多边形原样另存为 DXF，而是建立一条独立、可复现、
可审计的制造版图编译链：

1. 把所有独立几何参数显式写入 blueprint 配置，不依赖脚本常量或隐式默认；
2. 读取与仿真相同的 cavity/cladding、层数、finite geometry 和 cladding shift；
3. 使用与 full/quarter finite 完全相同的完整名义几何；
4. 在名义几何上应用有明确物理含义和来源的孔局部修正与全局空间包络修正；
5. 把 dose 建立为一等 process variable，并为后续 dose response/function 预留稳定接口；
6. 对补偿后的几何执行版图规则检查；
7. 导出单位、图层、极性和坐标约定明确的 `.dxf`；
8. 以精简输出为原则，只保存生产所需 DXF、一个综合快速预览 PNG 和一个必要的机器可读
   汇总文件，使加工文件仍能追溯到参数、工艺 profile、代码版本和输出哈希。

本任务不启动 COMSOL，不读取或修改 `.mph`，不移动、覆盖或删除现有
`scripts/.out/` 正式结果。

## 2. v1 范围与非目标

### 2.1 v1 范围

- v1 以当前六角晶格 finite cavity 为唯一主实现和正式验收对象；即使设计通过
  finite-quarter 求解，版图也必须从同一个完整结构生成。
- 从 v1 起就把 lattice、unit-cell pattern、finite footprint 和 correction model 分离为
  带 `kind/schema_version` 的接口，不能把“六角晶格 + 六个三角孔”写死在 DXF writer
  或 DRC 中。square、其他 Bravais lattice、其他孔形和任意 imported polygons 作为后续
  geometry source 接入，不要求在 v1 完成生产验收。
- 支持当前 `cladding_shift_geometry = groups | ellipse`、方向性 x/y shift 和逐层
  shift profile。
- 支持一个明确选择的 finite shift case；当共享设计参数包含扫描列表时，版图配置必须
  显式选中一个 case，不允许默认导出全部扫描点。
- 支持两级核心几何修正：对每个孔应用同一类 feature-local correction，以及按孔在整个
  写场中的位置应用 spatial envelope correction；同时支持可选角点策略、全局 x/y 比例
  补偿、版图坐标变换、设备网格量化、die/marker/keep-out 配置和几何 DRC。
- dose 在 v1 作为显式可调 process variable 进入 correction context 和输出汇总；v1 不要求
  建立完整连续 dose-response 模型，但接口不得阻止后续实现 `D(x,y)` 和
  `geometry_response(dose, position, feature)`。
- 生成工程预览与正式发布两种输出状态。只有正式发布状态且所有强制 DRC 通过时，才允许
  生成名为 `production.dxf` 的文件。

### 2.2 v1 非目标

- 不做三维刻蚀、侧壁角度、刻蚀深度或材料堆栈仿真。
- 不根据主观经验自动推断蚀刻 bias、圆角、proximity correction、dose function 或剂量补偿。
- 不把 COMSOL 的 SBC、PML、空气域或 simulation footprint 自动当作芯片/曝光边界。
- 不在 v1 自动排布整片 wafer、多 die 阵列或多器件 reticle。
- 不在 v1 实现 GDSII/OASIS；但几何计划与 writer 分离，为后续 writer 预留接口。
- 不把 mesh、中心频率、本征模式数等仿真专属参数写入制造版图 identity。

## 3. 当前代码基础与必须先解决的问题

### 3.1 可直接复用的基础

`comsol_workflow/finite_geometry.py` 中的纯 Python `FiniteGeometryPlan` 已经包含当前项目
可作为对照和迁移来源的：

- 完整 finite 结构的 hole polygons；
- cavity/cladding、parent lattice index、fragment group、layer 和 shift 信息；
- hex/square footprint 与完整结构边界；
- full/quarter 共用的几何语义；
- 不依赖 `mph` 的纯 NumPy 几何构建路径。

实施期应充分复用这些已经验证的公式和测试样例，但 standalone blueprint core 不得把
`FiniteGeometryPlan`、`scripts.run_main` 或 `mph` 类型作为公共 API。项目适配层负责把
当前 finite plan 转换为 blueprint 自己的中性 `NominalLayout`；迁出仓库后，只需替换或
带走 geometry source，不改 correction、dose、DRC、preview 或 DXF writer。

### 3.2 必须消除的隐藏/重复几何

当前晶格常数 `A=0.82 um`、参考尺寸和“紧凑参数 -> 六个三角孔”的转换仍分散在
`run_unit_cell_band.py`、`run_strip_1d.py` 和 `run_finite.py`。制造版图不能 import
这些主入口，也不能依赖仅存在于脚本常量中的物理尺寸。

实施时先把这些逻辑迁移/适配到独立 blueprint geometry source，并保留与当前项目的
parity adapter，例如：

```text
blueprints/src/blueprint_core/geometry/
  base.py                 # 通用 geometry source protocol
  hex_lattice.py          # 当前六角晶格坐标与 finite placement
  six_triangle_cell.py    # compact params -> six triangular holes
  comsol_parity.py        # 仓库内开发/一致性校验适配器，非 standalone runtime 必需
```

迁移时优先提取当前纯几何实现，不重新推导另一套公式。固定样例必须断言六个孔的顶点、
顺序、面积和中心以及完整 finite placements 逐项一致。blueprint 的生产运行路径只依赖
自身 package；`comsol_parity` 仅在仍位于当前仓库时用于交叉验证，未来可随测试 fixture
一起移除。

### 3.3 晶格常数参数化策略

所有独立几何参数必须显式进入 `blueprints/parameter.json`。其中包括但不限于：

- length unit、坐标原点和方向；
- lattice kind、period 或完整 basis vectors、lattice orientation；
- unit-cell pattern kind、参数化版本、cavity/cladding 的全部孔形输入；
- cavity/cladding 层数、finite footprint kind 和外部几何范围；
- cladding shift kind、x/y factor、逐层 profile 及其全部参数；
- die、marker、keep-out 和 correction field 的几何范围。

派生量不需要重复作为自由输入，但其公式必须版本化，且不能依赖未写明的 magic number。
例如由 period 推导的 `R0` 应记录公式版本，而不是再提供一个可能冲突的独立值。

当前集成期仍建议在共享设计参数中显式增加 `lattice_period_um`，当前活动值写为 `0.82`。
兼容策略为：

- 当前活动 `scripts/parameter.json` 显式保存该值；
- 历史仿真 fixture 缺失该字段时，现有仿真入口可暂时回退 `0.82` 并记录 legacy 来源；
- blueprint 的独立参数文件必须显式包含该值，正式发布模式拒绝任何隐式回退；
- 归一化参考常数可继续作为参数化实现细节，但 resolved geometry 必须记录最终实际单位和
  真实孔尺寸。

## 4. 建议目录与模块结构

```text
blueprints/
  pyproject.toml
  parameter.json
  run_blueprint.py
  README.md
  src/
    blueprint_core/
      __init__.py
      config.py
      models.py
      corrections.py
      dose.py
      drc.py
      dxf.py
      preview.py
      geometry/
        __init__.py
        base.py
        hex_lattice.py
        six_triangle_cell.py
        comsol_parity.py
  tests/
    test_config.py
    test_hex_geometry.py
    test_corrections.py
    test_dose_interface.py
    test_drc.py
    test_dxf.py
  .out/

tests/
  test_blueprint_comsol_geometry_parity.py
```

职责划分：

- `src/blueprint_core/geometry/base.py`：定义与 lattice/孔形无关的 geometry source protocol
  和中性 `NominalLayout`；DXF、correction 和 DRC 只消费该中性模型。
- `hex_lattice.py` 与 `six_triangle_cell.py`：v1 当前六角晶格和六三角孔实现，迁移并验证
  当前项目的纯几何公式，不导入 `mph` 或主工作流。
- `comsol_parity.py`：仅用于仍处于当前仓库期间，与当前 `FiniteGeometryPlan` 做一致性检查；
  不属于 standalone 生产路径。
- `corrections.py`：孔局部修正、空间包络修正、坐标变换和 grid snap。
- `dose.py`：dose context、dose field 和 process-response protocol；v1 可使用 scalar dose 和
  pass-through response。
- `drc.py`：对中性 layout 执行几何和加工规则检查。
- `dxf.py`：只负责把已经验证的 layout 写为 DXF，并负责回读校验。
- `preview.py`：生成单张综合快速预览，不承担几何计算。
- `run_blueprint.py`：唯一用户入口；解析配置、预检、选择 shift case、组织 staging、写
  单一汇总报告并最终化输出。
- `blueprints/parameter.json`：自包含几何、工艺/dose、输出合同和 DRC 合同。

standalone core 的生产依赖只能是标准库及明确声明的通用数值/CAD 依赖，不得反向依赖
`scripts`、`comsol_workflow`、`mph` 或当前仓库的输出目录。总项目文件夹固定命名为
`blueprints/`，内部采用独立 `pyproject.toml` 和 `src/blueprint_core` 标准 package 结构，
不把子目录直接加入 `sys.path`。功能稳定后，整个 `blueprints/` 目录可直接迁往独立仓库。
当前仓库内的命令为：

```powershell
uv run --project blueprints python blueprints\run_blueprint.py `
  --config blueprints\parameter.json
```

## 5. 参数与来源合同

### 5.1 自包含参数原则

根据“所有几何参数都要显化”和未来独立迁出的要求，blueprint 的生产入口以自身
`parameter.json` 为唯一运行输入，不在运行时读取当前项目的 `scripts/parameter.json`。
参数关系为：

```text
blueprints/parameter.json
  - 全部独立 nominal geometry 参数
  - 唯一 resolved shift pair
  - fabrication profile
  - dose variable/interface configuration
  - local correction + spatial envelope correction
  - layout coordinate contract
  - layer map、die、markers、keep-outs
  - DRC 与 release policy

optional source_trace
  - 原 scripts/parameter.json 的路径/哈希
  - 只用于集成期来源审计和 parity preflight
  - 缺失时不影响 standalone 运行
```

为降低手工转录风险，集成期可提供只读 importer，把当前 `scripts/parameter.json` 解析为
blueprint 配置草案；但 importer 生成后必须显式展开全部几何值，生产入口不得继续追随源
文件的后续变化。运行所需参数、DRC、摘要和哈希合并写入唯一 `layout_report.json`，不再
分别输出 resolved parameters、manifest 和 DRC JSON。

### 5.2 blueprint 参数块

| 参数块 | 必需内容 | 说明 |
| --- | --- | --- |
| `geometry` | units、lattice、unit cell、finite placement、resolved shift | 全部独立几何量显式，不允许生产模式 fallback |
| `parameterization` | geometry schema 与公式版本 | 区分六三角孔、未来其他孔形/晶格实现 |
| `source_trace` | 可选 source path/hash | 只作来源审计，不是运行依赖 |
| `fabrication` | `profile_id`、local correction、spatial envelope、grid、corner、scale | 所有长度带 `_nm` 或 `_um` 单位后缀 |
| `dose` | scalar/field kind、value、unit、response model | v1 允许 scalar + pass-through，接口必须稳定 |
| `layout` | origin、rotation、mirror、translation | 设计坐标到加工坐标的唯一变换 |
| `layers` | 逻辑 feature 到 DXF layer 的映射 | 可把 cavity/cladding 合并或保持独立 |
| `die` | outline、edge keep-out、是否导出 | 与 COMSOL footprint 独立 |
| `markers` | marker 类型、坐标、尺寸、layer | 未确认前默认不生成 |
| `drc` | min hole width、min edge、min bridge、boundary clearance | 违反强制规则时拒绝正式输出 |
| `export` | mode、DXF version、flatten、output policy | `engineering_preview` 或 `released` |

解析器要求严格字段、严格单位和有限数值，拒绝未知 enum、单边 x/y 配置、未 resolved 的
scan、隐藏默认、负的最小尺寸及 conflicting fields。`parameter.json` 必须保持为无注释
标准 JSON。

## 6. 名义几何数据流

```text
blueprints/parameter.json
  -> validate all explicit geometry/process/dose fields
  -> compact cavity/cladding params -> six local triangular holes
  -> selected GeometrySource(kind=hex_lattice_v1)
  -> standalone full NominalLayout
  -> tagged nominal features
  -> feature-local correction
  -> spatial envelope correction(position)
  -> dose/process response context
  -> DRC-clean LayoutPlan
  -> flattened DXF + one preview PNG + one layout_report.json
```

blueprint 生产入口不得 import `ACTIVE_PARAMETERS`、`run_finite.py`、`FiniteGeometryPlan`
或 `mph`。仅 parity test/adapter 可以读取当前项目对象，并把它转换为相同的中性模型。

每个孔建立稳定 ID，建议编码：

```text
<region>/<parent_m>/<parent_n>/<group_id>/<local_hole_index>
```

运行时中性 feature 至少携带 nominal/corrected polygon、region、parent index、layer、
cladding shift、local correction、envelope value、dose context、grid displacement、最终
DXF layer 和 DRC 状态。默认不把逐孔表写成独立 CSV；仅在显式 debug 模式下临时导出。

## 7. 制造修正合同

### 7.1 修正顺序

修正顺序必须固定并写入唯一 `layout_report.json`：

1. 从独立 hex geometry source 得到已经包含电磁设计 modulation 的 nominal hole polygons；
2. 在每个孔中心位置计算 spatial envelope value 和 dose context；
3. 由 base local parameters、envelope contribution 和 dose response 解析出该孔唯一的
   effective correction parameters；
4. 以孔自身 centroid/local frame 应用一次 feature-local geometry correction；
5. 应用角点策略；
6. 若存在经过标定的 placement distortion，再修正 feature center；该过程与孔尺寸修正分离；
7. 应用版图旋转、镜像、平移；
8. 把 die、marker 和其他加工结构变换到同一坐标合同；
9. 按 writer address grid 量化；
10. 删除量化造成的相邻重复顶点，但不得静默删除整个孔；
11. 执行最终 DRC；
12. 只有通过 release gate 后才写 production DXF。

### 7.2 第一级：每孔局部整体修正

每个孔都在自身 local frame 内调用同一个、版本化的 `FeatureCorrection`。v1 数据模型预留
以下 operator，具体启用集合和顺序必须显式：

- `edge_offset_nm`：沿孔边界法线做等距 offset；
- `isotropic_scale`：相对孔 centroid 做统一尺寸缩放；
- `anisotropic_scale_xy`：在明确的 local/global axis 中做 x/y 尺寸缩放；
- `corner_strategy`：sharp、rounded 或后续的 corner pre-compensation。

其中 edge bias 作用于最终孔边界，不通过修改 `b0_nm` 间接实现。正 bias 定义为沿孔外法线
扩大空气孔/蚀刻开口，负 bias 缩小空气孔。cavity 与 cladding 可以使用不同的 base
correction，但若没有区域差异标定，推荐引用同一个 profile。offset 必须显式声明 join style
和 miter limit；任何自交、退化或小于最小尺寸的结果均报错，不允许自动丢孔。

### 7.3 第二级：整个结构的空间包络修正

写场中心与边缘的加工差异用独立 `SpatialEnvelope` 表示。对第 i 个孔，以其名义中心
`(x_i, y_i)` 计算：

```text
e_i = E(x_i, y_i; write_field_geometry, coefficients)
d_i = D(x_i, y_i; dose_parameters)
c_i = c_base + c_envelope(e_i) + response(d_i, x_i, y_i, feature_i)
corrected_i = FeatureCorrection(nominal_i, c_i)
```

这里 `c_i` 是有效几何修正参数，不是直接把 polygon 坐标乘以一个来源不明的系数。所有
envelope 坐标范围、中心、归一化、系数、单位、作用对象和公式版本都必须显式。

v1 建议至少支持：

- `constant`：用于关闭空间变化或做对照；
- `radial_even_polynomial`：围绕写场中心的偶次径向包络；
- `separable_even_polynomial_xy`：x/y 非对称的中心到边缘修正；
- 后续预留 `sampled_map`：读取预实验/SEM 标定网格并插值。

默认建议 envelope 只调制孔的 size/bias，不移动 lattice center。若实际设备还存在 placement
distortion，应使用独立 vector field `delta_position(x,y)`，避免把 CD 修正和坐标畸变混为
一谈。此默认仍需用户确认。

### 7.4 第三级：dose 变量与响应接口

dose 在 schema 和运行模型中是一等变量，而不是写在备注里的数值。接口分成三层：

1. `DoseField D(x,y)`：v1 至少支持全版图 scalar dose；后续支持 spatial dose field；
2. `ProcessResponse`：把 dose、位置、feature 类型和 profile 转换为几何 correction contribution；
3. `DoseEncoder`：决定 dose 由外部曝光 recipe 处理，还是编码为 DXF layer/bucket。

v1 不虚构连续 dose-response 公式。允许的最小实现是：配置一个显式 scalar dose 和单位，
使用 `pass_through` response，几何修正仍由已知 base/envelope 参数给出；dose 值写入单一
report 并进入 layout identity。获得预实验数据后，再增加 table/interpolation response：

```text
dose -> measured CD/edge bias/corner response -> effective geometry correction
```

DXF 本身通常不完整表达曝光 dose。正式实现前必须确认下游设备是从外部 job recipe 读取
dose，还是使用 layer name/number 区分 dose bucket；未确认前不得擅自把 dose 编入图层。

### 7.5 角点与 proximity correction

v1 至少保留 `sharp` 和可配置 `rounded` 的数据结构。是否实现 rounded production geometry
由加工方实际角点模型决定。serif、局部预补偿、密度相关 shape bias 和剂量 PEC 只作为
版本化 rule 扩展点；没有 SEM/test coupon 标定时不得启用。

### 7.6 设计边界与加工边界分离

- `FiniteFootprint` 可用于 QA 和 bounds 检查，但不自动进入 etch layer。
- die outline、scribe lane、edge exclusion、write-field 和 alignment marker 只来自 blueprint
  配置或经过批准的模板。
- construction/reference/annotation 不得混入 production etch layer。

## 8. DXF 输出合同

在加工方未给出特定合同前，工程预览建议采用以下临时默认；正式发布前必须确认：

- 坐标数值以 `um` 表示，并正确设置 DXF `INSUNITS`；
- 工艺孔使用 closed `LWPOLYLINE`；
- production 文件 flatten 全部几何，不依赖 block/insert；
- etch layer 不使用 spline、hatch、开放 polyline 或文字；
- 每个 polygon 使用统一 winding，量化后不含连续重复点；
- cavity/cladding 在内部保持逻辑身份，最终物理 layer 由配置映射；
- DXF version 暂按 `R2010` 设计，若加工方要求 R12/R14 则以其导入合同为准；
- 生产文件只包含真实加工图层，辅助边界只出现在 QC PNG/JSON，不放入 production DXF。

建议使用 `ezdxf` 负责写入和回读；使用 `shapely>=2` 完成可靠 polygon offset、validity、
空间索引和邻近 DRC。具体依赖版本通过 `uv` 锁定，不在方案阶段安装。

## 9. 输出布局、命名与原子性

默认每个成功 case 只保留三个文件：

```text
blueprints/.out/
  blueprint_<shape>_<cavN-cladM>_<shift-label>_<process-id>_<hash8>/
    engineering_preview.dxf 或 production.dxf   # 两者只存在一个
    preview.png
    layout_report.json
```

`preview.png` 是一张综合图，至少组合全版图、中心/边缘代表孔、nominal/corrected 对比、
spatial envelope 摘要、dose/profile、关键尺寸和 DRC 结论；不为每个诊断维度生成单独 PNG。

`layout_report.json` 是唯一默认 JSON，合并以下必要信息：完整 resolved input、source trace、
geometry/process/dose schema version、counts/bbox、correction 范围、DRC 结论、DXF 回读摘要、
代码版本和文件 SHA-256。默认不生成 `resolved_parameter.json`、`manifest.json`、
`drc_report.json` 或逐孔 CSV。只有显式 `--debug-artifacts` 才允许在临时 debug 子目录输出
详细逐孔数据；debug 产物不得混入正式交付目录。

命名 identity 包含 period、cavity/cladding 参数、层数、finite geometry、shift pair 和
fabrication profile revision；不包含 mesh、eigenmode count 或 center frequency。

入口使用 case-local `.staging`：

1. 先检查目标目录和已知文件；
2. 在 staging 中只写当前模式的 DXF、综合 preview 和单一 report；
3. 回读 DXF 并完成最终 QA；
4. 确认无冲突后一次性最终化；
5. 默认拒绝覆盖已有 case；相同 resolved identity 可明确报告已存在。

`blueprints/.out/` 是独立的制造输出根，不在 `scripts/.out/` 下创建第四个一级目录。
实施时同步更新根 README、相关 README、AGENTS 输出边界说明和 `.gitignore` 验证。

## 10. DRC 与发布门槛

### 10.1 名义几何一致性

- 同一参数下 blueprint 六角晶格与当前 finite workflow 生成相同孔数、孔顶点、bbox 和
  region counts；
- groups/ellipse、zero/directional shift 的六角代表样例全部覆盖；
- 用至少一个非六角 synthetic geometry source 验证 correction/DRC/DXF 不依赖六角专属类型，
  但 v1 不宣称该结构已达到生产验收；
- quarter 来源仍导出完整、对称的 full layout。

### 10.2 修正后几何检查

- polygon 非退化、无自交、面积方向统一；
- 无重复孔、孔间重叠或跨越禁止区域；
- 最小孔宽、最短边、最小 silicon bridge 和 die clearance 满足配置；
- marker/die 与器件 keep-out 不冲突；
- 所有生产坐标位于 address grid；
- 记录 grid snapping 的最大、平均位移；
- 记录 nominal/corrected 面积、bbox 和孔尺寸变化分布；
- DRC 不通过时仍可生成报告和 PNG，但不得生成 `production.dxf`。

### 10.3 DXF 回读检查

用独立 read-back 路径重新打开 DXF，并断言：

- DXF version 与 units；
- 允许的 layer 集合；
- 生产 layer entity 全部为 closed polyline；
- entity 数量、顶点数量、bbox 与 writer 输入一致；
- 回读面积在 grid tolerance 内一致；
- 不存在未允许的 text、block、spline、hatch 或开放路径；
- 最终文件 SHA-256 写入 `layout_report.json`。

### 10.4 发布状态

- `engineering_preview`：允许使用尚未 release 的工艺 profile，文件名必须明确为 preview，
  不得同时生成 `production.dxf`。
- `released`：必须有非空 `profile_id/revision/provenance`、全部数值化工艺规则、明确 DXF
  合同、通过全部强制 DRC，并在 `layout_report.json` 标记 release 状态。

## 11. 测试与验收方案

### 11.1 自动测试

- parameter parser：合法配置、未知字段、单位错误、shift 歧义和 release gate。
- hex geometry parity：独立实现与当前项目的六孔顶点和完整 finite placement 严格一致。
- finite parity：小型/实际层数 hex、groups/ellipse、zero/directional shift。
- generic source test：非六角 synthetic polygons 可经过 correction/DRC/DXF 全链路。
- bias analytic tests：已知三角形 offset 后边到原边的法向距离等于目标 bias。
- local correction tests：centroid scale、operator 顺序、region mapping 和 correction provenance。
- envelope tests：constant、radial、separable x/y 的中心/边缘值、归一化和对称性。
- dose interface tests：scalar pass-through、单位校验、layout identity 和未支持 encoder 拒绝。
- grid tests：量化误差有界、重复点清理、退化孔拒绝。
- DRC fixtures：最小孔、最小 bridge、重叠、keep-out 和 die clearance 的正反样例。
- DXF round trip：units、layers、closed entities、counts、bbox、area、禁止实体。
- orchestration test：mock writer/geometry，不启动 COMSOL，不写 `scripts/.out/`。

### 11.2 聚焦命令

计划实施后至少运行：

```powershell
uv run --project blueprints python -m pytest -q `
  blueprints\tests\test_config.py `
  blueprints\tests\test_hex_geometry.py `
  blueprints\tests\test_corrections.py `
  blueprints\tests\test_dose_interface.py `
  blueprints\tests\test_drc.py `
  blueprints\tests\test_dxf.py

uv run python -m pytest -q tests\test_blueprint_comsol_geometry_parity.py

uv run --project blueprints python -m py_compile blueprints\run_blueprint.py
uv run --project blueprints python -m py_compile blueprints\src\blueprint_core\corrections.py
uv run --project blueprints python -m py_compile blueprints\src\blueprint_core\dose.py
uv run --project blueprints python -m py_compile blueprints\src\blueprint_core\dxf.py
git diff --check
```

随后用一个小层数 fixture 生成工程预览，实际检查 PNG 和 DXF 回读报告。正式目标结构只在
参数与加工合同确认后生成；本任务不需要 COMSOL license。

### 11.3 外部验收

自动 round-trip 不能替代加工方导入验收。正式文件至少需在加工方指定 CAD/mask 工具中
验证一次：单位、朝向、极性、layer、孔数量、边界尺寸和无隐式缩放。若加工方提供 import
截图或检查报告，将其路径/编号记录到 `layout_report.json` 或 release note。

## 12. 分阶段执行计划

### 阶段 0：锁定输入与加工合同

- 回答第 15 节仍未决的 P0 问题；
- 确定 v1 输出范围和一个目标 finite case；
- 锁定独立 geometry schema、local/envelope/dose 三层 correction contract；
- 至少获得可用于 engineering preview 的 provisional bias、grid、minimum feature/bridge 和
  DXF import 合同；
- 在本文件记录最终决定，才进入代码实施。

验收：不存在影响几何、单位、极性或输出内容的未决 P0 项。

### 阶段 1：建立可独立迁出的纯 Python 名义几何

- 在 `blueprints/src/blueprint_core/geometry/` 建立通用 source protocol；
- 迁移/复用当前六角晶格、六三角孔和 shift/profile 纯几何公式；
- 显式参数化所有独立几何量，不保留 `A/R0/B0` 类隐藏物理常数；
- 建立仅供仓库集成期使用的 COMSOL parity adapter；
- 保证 standalone 生产路径不 import 当前项目 package。

验收：代表参数下当前项目与独立实现的六个孔和完整 finite placement 逐点一致；将
整个 `blueprints/` 复制到隔离环境后仍可运行其纯几何测试。

### 阶段 2：建立 blueprint 配置与 nominal layout

- 新建独立 `blueprints/` 子项目、参数文件和入口；
- 建立自包含 geometry schema，并提供可选的 current-project importer/source trace；
- 显式展开所有几何参数并写入一个 resolved shift pair；
- 从 standalone hex source 建立 tagged nominal features；
- 生成 identity 和 preview-mode 输出目录。

验收：小型与实际层数结构的 nominal counts/bbox 与 finite workflow 一致，不启动 COMSOL。

### 阶段 3：孔局部修正、空间包络、dose 框架与 DRC

- 实现 feature-local correction operators；
- 实现 constant/radial/separable spatial envelope；
- 实现 scalar dose、pass-through response 和 future encoder protocol；
- 实现 corner strategy、scale/coordinate transform、grid snap；
- 实现 die/marker/keep-out 数据模型；
- 实现空间索引 DRC 和 release gate；
- 输出单一综合 preview 与内存/单一 report 中的必要摘要。

验收：中心/边缘样例得到精确可解释的 effective correction；改变 scalar dose 可以稳定进入
process context 和 identity；所有故障 fixture 在 DXF 写入前失败。

### 阶段 4：DXF writer 与原子输出

- 实现 flattened production entities、units 和 layer map；
- 实现 DXF 回读验证与 SHA-256；
- 实现 staging/finalize 和拒绝覆盖；
- 建立 engineering/released 两种文件命名。
- 默认严格只保留一个 DXF、一个 preview PNG 和一个 `layout_report.json`。

验收：round-trip tests 全部通过，production 文件不含禁止实体。

### 阶段 5：工程样例审批

- 只用确认的一个目标 design case 和 provisional/released process profile；
- 生成工程 DXF、综合 preview 和单一 report；
- 检查中心/边缘孔、cavity-cladding 过渡、shift、die/marker 位置和整体朝向；
- 用户批准后才把 profile 标记为 released。

验收：用户确认视觉与数值合同，加工方导入测试无单位/layer/极性问题。

### 阶段 6：正式输出与文档回填

- 生成唯一 `production.dxf`；
- 在唯一 `layout_report.json` 中保存 resolved parameters、correction/dose、DRC 和文件哈希；
- 更新 README 使用方式；
- 在本文件回填实际改动、测试、外部验收、遗留问题和最终状态。

验收：正式 DXF 可从 `layout_report.json` 完整复现，且没有覆盖任何既有正式结果。

## 13. 风险与控制

| 风险 | 后果 | 控制方式 |
| --- | --- | --- |
| 仿真与独立版图几何逐渐漂移 | 实际加工结构偏离仿真 | 迁移已验证公式 + 固定 fixture + 集成期 parity tests |
| 独立模块仍暗中 import 当前项目 | 将来无法迁出 | standalone dependency test，项目适配器与 core 分层 |
| 隐藏的 `A=0.82` 或其他 magic number | 整体尺寸错误 | 所有独立几何参数显式化，production 禁止 fallback |
| bias 正负号误解 | 孔整体放大/缩小方向错误 | 定义空气孔正 bias，解析样例与单一 report 明示 |
| envelope 同时混入尺寸与坐标畸变 | 修正无法解释或重复应用 | size/bias envelope 与 placement vector field 分开 |
| 没有数据却拟合 dose response | 产生虚假的加工精度 | v1 scalar/pass-through；只有标定数据才能启用 response |
| 把 simulation footprint 当 die | 多出错误加工边界 | die/marker 完全独立配置 |
| DXF unit 或 block 解释差异 | 导入缩放、漏图元 | INSUNITS + flatten + 外部工具验收 |
| 扫描参数未选定 | 一次输出多个或错误 case | 必须显式选择唯一 shift pair |
| grid snap 造成窄桥或退化孔 | 不可制造 | snap 后执行最终 DRC，不静默丢孔 |
| 30k 级孔数导致 O(N^2) DRC | 运行不可接受 | STRtree/空间索引，仅检查邻近候选 |
| 辅助图层被误加工 | 非预期结构 | production layer allowlist，辅助信息只进综合 PNG/report |
| 默认输出过多中间文件 | 交付混乱、难以识别正式文件 | 默认严格三文件；详细逐孔数据仅显式 debug 且不进交付目录 |
| 工作树已有大量用户修改 | 覆盖当前开发内容 | 逐文件小步修改，保留现有 diff，不做 reset/checkout |

## 14. 回滚方案

- 删除新增的 `blueprints/` 子项目、专用模块、测试和命令行入口；
- 若为 parity 临时调整当前项目 import，则恢复原调用，但不改变参数值或独立配置；
- `lattice_period_um` 为兼容新增字段，回滚时可保留而不影响旧逻辑；
- 所有 blueprint 输出位于独立 `.out`，不涉及 `scripts/.out` 迁移；
- 不执行 COMSOL、不修改 `.mph`，因此回滚不涉及计算结果恢复。

## 15. 需要确认与明确的问题

### 15.1 已确认的约束

| 编号 | 已确认事项 | 对方案的影响 |
| --- | --- | --- |
| C-1 | 所有几何参数必须显化 | blueprint 参数自包含；production 禁止隐藏常数和 fallback |
| C-2 | v1 围绕六角晶格，但要兼容未来其他结构 | 六角是首个 source；correction/DRC/DXF 只依赖通用 polygon/layout protocol |
| C-3 | 充分复用当前项目几何代码，但目标是独立模块，完善后直接迁出 | core 不依赖当前项目；只保留可删除的 parity adapter |
| C-4 | 修正分为每孔局部整体修正和全结构空间包络修正 | correction pipeline 明确拆成 local operator 与 spatial envelope |
| C-5 | dose 是显式可调变量，v1 不完整求解但必须预留函数框架 | 建立 DoseField、ProcessResponse、DoseEncoder 三层接口 |
| C-6 | 输出必须精简 | 默认严格为一个 DXF、一个综合 PNG、一个 `layout_report.json` |

### 15.2 P0：进入核心实施前仍需确认

| 编号 | 问题 | 当前建议/默认 | 为什么必须确认 |
| --- | --- | --- | --- |
| P0-1 | 首个验收结构的 cavity/cladding 层数、compact hole 参数和唯一 x/y shift pair 是什么？ | 把当前最终选定值完整复制为 blueprint 自包含参数，不保留 scan list | 决定 parity fixture 和首个版图 identity |
| P0-2 | 晶格常数是否精确为 `0.82 um`？六角晶格的 origin、basis orientation 和 finite footprint 采用哪一版本？ | 沿用当前项目公式，但把数值和公式版本全部显式化 | “六角晶格”仍可能有旋转和边界语义差异 |
| P0-3 | v1 每孔局部修正采用 edge offset、centroid scale，还是二者按固定顺序组合？ | 优先 edge offset；若同时启用，必须显式顺序且各有物理来源 | 两种操作对三角孔面积和尖角影响不同 |
| P0-4 | cavity 与 cladding 是否使用同一 local correction？初始 provisional 数值是多少？ | 没有分区标定时使用同一 profile | 避免人为改变 cavity-cladding 光学过渡 |
| P0-5 | spatial envelope 只调制孔 size/bias，还是也移动孔中心？ | v1 只调制 size/bias；位置畸变另设 vector field | 防止把 CD 差异与 writer 坐标畸变混合 |
| P0-6 | v1 包络选择径向偶次、多项式 x/y 分离，还是二者都支持？写场中心和有效范围是多少？ | 实现 constant + radial even + separable even；参数全部显式 | 决定 envelope schema、归一化和中心/边缘测试 |
| P0-7 | envelope 输出是 additive edge bias（nm）还是 multiplicative size scale？ | 优先 additive bias，单位 nm；需要 scale 时作为独立 operator | 必须保证系数的物理意义可解释 |
| P0-8 | v1 scalar dose 的数值、单位和适用曝光设备是什么？ | dose 明确写值和单位，response 暂用 `pass_through` | dose 不同设备可能使用完全不同单位/语义 |
| P0-9 | DXF 中只需要蚀刻孔，还是还要 die outline、alignment marker、scribe/write-field？ | 首版只输出确认过的蚀刻孔，其他逐项启用 | simulation boundary 不能代替加工边界 |
| P0-10 | 蚀刻/曝光极性是什么；闭合孔轮廓代表 clear/open 还是 dark/protected？ | 在 layer contract 和单一 report 中明确 | DXF 本身通常不携带完整 mask 极性语义 |
| P0-11 | writer/grid 最小地址单位是多少？ | 由设备给定，禁止把显示精度当 address grid | 决定坐标量化和最终尺寸误差 |
| P0-12 | 最小 hole width、最短边、最小 silicon bridge、die clearance 是多少？ | 全部作为强制 DRC 数值 | 没有阈值无法判断可制造性 |
| P0-13 | 加工方接受的 DXF version、单位约定和 layer 命名是什么？ | 暂定 R2010、um、flatten closed LWPOLYLINE | 必须匹配下游 CAD/mask importer |
| P0-14 | 设计坐标到加工坐标是否需要旋转、镜像、平移？坐标原点在哪里？ | 器件中心 `(0,0)`，不旋转、不镜像 | 错误镜像对方向性 shift 尤其危险 |

### 15.3 P1：生成 production DXF 前必须确认

| 编号 | 问题 | 当前建议/默认 |
| --- | --- | --- |
| P1-1 | 工艺 profile 的名称、revision、负责人和标定来源是什么？ | 使用不可变 `profile_id + revision`，记录 SEM/test coupon 或 foundry rule 来源 |
| P1-2 | x/y writer magnification 或 placement vector field 是否需要预补偿？ | 无标定时关闭，不能伪装为 released calibration |
| P1-3 | dose 由外部 job recipe 指定，还是需要按 DXF layer/dose bucket 编码？ | 优先外部 recipe；未确认前 DoseEncoder 不写 layer |
| P1-4 | cavity/cladding 是否最终合并到同一个物理 etch layer？ | 内部保留身份，生产层默认合并，除非加工方要求不同层 |
| P1-5 | 三角孔保持 sharp，还是输出 rounded/pre-compensated corner？若折线近似，最大 chord error 是多少？ | 没有标定时 sharp 只能作为 engineering profile |
| P1-6 | marker 模板、绝对坐标、尺寸、layer 和 keep-out 是什么？ | 未确认则不输出 marker |
| P1-7 | die 尺寸、edge exclusion 和器件相对 die 的位置是什么？ | 未确认则只验证孔阵列 bbox，不输出 die outline |
| P1-8 | 是否允许复用相同 identity 的已有输出，还是每次生成新 revision？ | 不覆盖；相同 identity 报告已存在，改变 profile revision 产生新目录 |
| P1-9 | 哪些 DRC 可降级为 warning，哪些必须 error？ | 几何非法、重叠、最小 bridge、单位/layer mismatch 必须 error |
| P1-10 | 谁批准 engineering profile 转为 released，批准信息如何记录？ | 在配置和单一 report 中记录 approver、日期、profile revision |
| P1-11 | 外部导入验收使用什么软件或加工方流程？ | 至少一次加工方指定工具导入并记录结果 |

### 15.4 P2：为后续版本预留，但不阻塞 v1

- 是否需要同一 die 上排布多个 cavity 参数或 shift 变体；
- 是否需要 unit-cell/strip/process-control test coupons；
- 是否需要 EBL write-field stitching、field boundary keep-out 或 spatial dose map；
- 是否需要基于局部图形密度的 shape PEC；
- 是否需要多次曝光、多 etch depth 或 overlay layer；
- 是否需要 GDSII/OASIS writer；
- 是否需要 wafer map、die ID、文字转轮廓和自动序列号；
- 是否需要从 SEM 测量自动反演新的 bias/corner profile；
- 是否需要把 `layout_report.json` 与实验数据库或样品编号关联。

## 16. 决策记录

| 日期 | 决策 | 状态 |
| --- | --- | --- |
| 2026-08-04 | 独立子项目总文件夹定名为 `blueprints/`，输出到其自身 `.out/` | 已实施 |
| 2026-08-04 | 内部 Python 包定名为 `blueprint_core/`，与外层 `blueprints/` 项目目录区分 | 已实施 |
| 2026-08-04 | 所有独立几何参数必须显式；blueprint 参数必须可自包含运行 | 已实施并通过严格配置校验 |
| 2026-08-04 | v1 以六角晶格为主，但 core API 必须允许未来其他结构 | 已实施中性 layout/GeometrySource 边界；v1 CLI 仅启用 hex source |
| 2026-08-04 | 复用/迁移当前几何实现，但 blueprint core 最终必须能直接迁出当前项目 | 已迁移纯 Python 公式并完成 parity 测试 |
| 2026-08-04 | 几何修正分为每孔局部整体修正和全结构空间包络修正 | 框架已实施；首个工程配置的标定值仍为零占位 |
| 2026-08-04 | dose 作为显式变量并预留 field/response/encoder 框架 | 已实施；v1 使用 unassigned/scalar + pass-through |
| 2026-08-04 | 默认输出精简为一个 DXF、一个综合 PNG、一个必要 JSON | 已实施并以 staging 三文件合同验证 |

## 17. 执行记录

### 实际改动

已新增可独立迁出的 `blueprints/` 项目，内部包名为 `blueprint_core/`。生产代码不导入
`scripts`、`comsol_workflow` 或 `mph`，也不启动 COMSOL。主要实现包括：

- 严格、自包含、拒绝未知字段的 `parameter.json` 配置合同，以及 engineering/released
  release gate；所有当前六角 finite geometry 独立参数、writer grid、修正、dose、layout、
  DRC 和 export 参数均已显式；
- `GeometrySource` 和中性 `LayoutPlan/LayoutFeature`，以及迁移后的六角晶格、六三角孔、
  groups/ellipse cladding shift 和 uniform/tanh-power profile；
- feature-local edge offset/isotropic/anisotropic 修正、constant/radial/separable spatial
  envelope、全局 placement/scale/layout transform 和 writer-grid snap；
- `DoseField`、`ProcessResponse`、`DoseEncoder` 稳定接口，v1 实现 unassigned/scalar、
  pass-through response 和 external-recipe encoder；
- polygon 合法性、重复/相交、最小孔宽/边/桥、boundary clearance、grid、keep-out 和 layer
  DRC；未赋值的加工阈值仍计算测量值，但 production release 必须全部赋值；
- flat DXF writer/readback（闭合 `LWPOLYLINE`、微米单位、单 `ETCH` layer）、一个综合预览
  PNG，以及包含 resolved config、source trace、源码树/Git 版本、修正、dose、DRC、DXF
  readback、输出哈希与 timing 的单一 `layout_report.json`；
- atomic staging、已有 case 拒绝覆盖、工程/生产 artifact identity 分离、完整三文件 manifest
  和 Windows 长路径规避；
- 独立 `pyproject.toml`、`uv.lock`、CLI wrapper、README、54 项 standalone tests，以及当前
  COMSOL geometry 对照的 2 项 parity tests。

当前自包含工程配置镜像 period `0.82 um`、cavity `(b0=245 nm, eta=0.96,
zeta=1.156)`、cladding `(b0=242 nm, eta=0.98, zeta=0.93)`、20/20 层、groups +
uniform shift、x/y factor `0`。由于实际 bias、dose、极性和加工 DRC 阈值尚未标定，配置
明确保持 `engineering_preview`，不会伪装成可投产的 `production.dxf`。

### 测试结果

- `blueprints/tests` 与 `tests/test_blueprint_comsol_geometry_parity.py`：`56 passed`；
- 当前完整 20/20 配置 preflight：29,526 个 feature，DRC `passed=True`，0 error，0 warning；
- 独立 DXF 回读：29,526 个闭合 `LWPOLYLINE`、88,578 vertices、唯一 `ETCH` layer、
  `$INSUNITS=13`（micrometers）、R2010/AC1024；
- 修正后孔 bbox：`[-33.133, -28.763, 33.133, 28.763] um`，DXF polygon 总面积
  `752.2203480002111 um^2`；测得最小孔宽 `0.1928759266 um`、最小桥宽
  `0.1155427295 um`、最小 boundary clearance `0.471 um`；
- `source_trace` 与当前 `scripts/parameter.json` SHA-256 匹配；production source 中未发现
  `comsol_workflow`、`scripts.run_main` 或 `mph` 引用；
- CLI/config 聚焦测试、`py_compile`、`git diff --check` 通过；实际 `preview.png` 已完成
  视觉检查；
- 全过程未启动 COMSOL，未修改、移动或覆盖 `scripts/.out/` 正式结果。

最终工程样例位于：

`blueprints/.out/blueprint_hex_engineering_preview_cav20-clad20_shiftx0-y0_engineering-unassigned-r0_2772286d/`

该 case 恰好包含 `engineering_preview.dxf`、`preview.png` 和 `layout_report.json`。

### 遗留问题与最终状态

v1 工程版链路已经完成，最终状态为 **engineering preview ready / production blocked by
unassigned process calibration**。转为生产发布前仍需确认并写入：实际 local/envelope
修正函数及系数、dose 数值/单位、clear/dark 极性、最小孔宽/边/桥和 boundary clearance
阈值、die/marker/keep-out（如需要），并完成目标 CAD/EBL 软件导入、单位、图层、闭合性、
极性、bbox 和关键尺寸验收。release gate 在这些字段未赋值时拒绝 `production.dxf`。

未来新增非六角结构时，下游 correction/dose/DRC/DXF/preview 无需改写；仍需为新结构增加
geometry source、对应严格配置 parser 和 case-label 分发。当前 synthetic polygon tests 已
验证下游模块不依赖六角专属类型。

执行中曾生成一个旧 identity 的工程验证 case 和若干测试缓存；自动清理命令被当前工具的
删除审批器拒绝，因此它们仍留在被 Git 忽略的本地目录中，不属于最终 artifact。确认删除
授权后应移除旧 case
`blueprint_hex_cav20-clad20_shiftx0-y0_engineering-unassigned-r0_199303db/` 及
`.test_tmp_*`/`.pytest_cache`/`__pycache__`，保留上述最终 case。
