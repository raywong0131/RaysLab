# 分层 tanh Cladding Shift：计划与执行记录

日期：2026-07-27
状态：已完成（代码与非 COMSOL 验证）

## 背景与目标

当前 finite、finite-quarter 和 strip 几何把同一个 `cladding shift` 施加到全部
cladding shell。cavity 最外层的位移为零，而第一层 cladding 直接达到目标位移，
因此位移场在 interface 处发生一步跳变。现有实现随后按 normal-fan 规则确定六边形
边区和角区的向内位移方向。

本任务把 cladding shift 改为由相对层数控制的绝对位移包络，在不改变现有
normal-fan 方向、cavity/cladding 单胞紧凑几何和六重对称性的前提下：

1. 以 cavity 最外层为第 `0` 层，cladding 从第 `1` 层开始编号；
2. interface 附近从零 shift 平滑增长；
3. 远端 cladding 回归当前扫描参数指定的目标 shift；
4. 减小相邻层位移跳变及其高空间频率成分，为降低边界散射提供可验证的几何方案；
5. 保留 uniform-shift 兼容模式，以便与历史结果进行严格对照。

本计划编写阶段只定义实现和验证方案；获得用户批准后已按此实施代码与共享配置
改造。实施过程未启动 COMSOL，也未修改已有正式计算结果。

## 已确认的当前行为

- `LatticePoint.shell` 使用六角轴坐标
  `max(abs(i), abs(j), abs(i + j))` 定义，适合作为保持 C6 对称的离散层编号。
- cavity 包含 `shell=0..BULK_RADIUS`；cladding 包含
  `shell=BULK_RADIUS+1..BULK_RADIUS+CLADDING_LAYERS`。
- 当前 `cladding_shift_factors` 是无量纲目标值 `f_target`，物理边区位移为
  `A * f_target`。
- 每个 cladding 单胞中的整组孔做相同的刚性平移；cavity 孔不平移，单胞内部孔间距
  不因 shift 改变。
- 边区沿对应六边形边的 inward normal 平移；角区沿相邻 inward normal 的角平分线
  平移，角区幅值保持为边区幅值的 `2/sqrt(3)`。
- strip 在完整晶格生成后连同 `modulation_shift` 向量一起刚性变换到 cut 坐标系；
  finite-quarter 复用 finite 的完整几何。

## 采纳的数学定义

### 相对层编号

对任一 cladding lattice point 定义：

```text
layer = point.shell - BULK_RADIUS
```

概念上的 cavity 最外层为 `layer=0`；实际 cladding 点的有效范围为
`layer=1..CLADDING_LAYERS`。不得按列表次序、欧氏半径或生成顺序推断层号。

### 目标值与逐层绝对位移

现有 `cladding_shift_factors` 继续表示远端目标值，不改成每层增量：

```text
f_target = 当前 shift case 的 cladding_shift_factor
f_layer  = f_target * envelope(layer)
```

逐层几何必须始终由理想晶格坐标直接构造：

```text
shifted_position(layer) = ideal_position(layer) + absolute_shift(layer)
```

禁止把本层 shift 累加到上一层已经移动后的坐标，否则远端位移会不断累计，无法恢复
为局部均匀的 cladding 晶格。

### tanh-power 包络

采纳的渐变函数为：

```text
envelope(layer) = tanh((layer / scale_layers) ** power)
```

首版活动参数为：

```text
kind         = "tanh_power"
scale_layers = 4.0
power        = 2.0
```

其性质为：

- `envelope(0)=0`；
- `power>1` 时连续模型在 interface 处的一阶导数为零；
- 单调趋近 `1`，但不通过有限层归一化强制等于 `1`；
- `scale_layers` 控制渐变长度，`power` 控制 interface 平坦程度与中部斜率集中程度。

不使用 `0.5 * (1 + tanh(...))`，因为它不能自然满足第 `0` 层严格为零；也不以
最外层值归一化整个 profile，因为归一化会掩盖“cladding 层数不足以到达远端目标”
的问题。

当前 `CLADDING_LAYERS=10`、`scale_layers=4`、`power=2` 时，预期包络约为：

```text
layer:     0      1      2      3      4      5      6      7      8      9     10
envelope: 0.000  0.062  0.245  0.510  0.762  0.916  0.978  0.996  0.999  1.000  1.000
```

### normal-fan 位移

保持现有方向与角区补偿，只把目标幅值乘以逐层包络：

```text
side_shift(layer)
  = A * f_target * envelope(layer) * side_inward_unit_normal

corner_shift(layer)
  = A * f_target * envelope(layer)
    * (2 / sqrt(3)) * corner_inward_bisector_unit_vector
```

同一 shell 上所有点使用相同 `envelope(layer)`；边区/角区只改变方向和既有幅值比例，
不得引入按方位角任意变化的标量包络。

## 配置方案

在 `scripts/parameter.json` 增加共享配置：

```json
"cladding_shift_profile": {
  "kind": "tanh_power",
  "scale_layers": 4.0,
  "power": 2.0
}
```

现有字段保持不变：

```json
"cladding_shift_factors": [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
```

读取和兼容规则：

1. `cladding_shift_profile` 缺失时使用 `uniform`，从而精确复现历史行为：
   `envelope(0)=0`，对所有实际 cladding 层有 `envelope(layer)=1`。
2. `kind` 首版只接受 `uniform` 和 `tanh_power`；未知值明确报错。
3. `tanh_power.scale_layers` 必须为有限正数，`power` 必须为有限且严格大于 `1`
   的数；活动默认值为 `4.0/2.0`。
4. `cladding_shift_factors` 仍只控制远端目标幅值；profile 参数不复制到每个
   shift case 中。
5. unit-cell 工作流继续忽略结构层数和 cladding shift profile；profile 变化不得
   改变 unit-cell 输出目录，也不得阻止 unit-cell 完成后的中心频率安全回写。
6. 不因新增可选字段升级现有 `schema_version=1`；旧配置可继续加载，但当前活动
   `parameter.json` 在实施时显式写入 `tanh_power` 配置。

配置加载后解析完整逐层包络，并计算以下诊断但不静默归一化：

- `outer_envelope = envelope(CLADDING_LAYERS)`；
- `bulk_start_envelope = envelope(max(1, CLADDING_LAYERS - 2))`；
- `max_first_difference = max(abs(f_layer - f_(layer-1)))`；
- `max_second_difference = max(abs(f_(layer+1) - 2*f_layer + f_(layer-1)))`。

当前默认配置的预期护栏为：最外层 `envelope>=0.999`，且最后约三层已经接近目标
shift。参数合法性与远端饱和度分开处理：前者不合法时拒绝加载；后者不足时在昂贵
计算 preflight 中明确报告，由用户决定是否调整渐变尺度，不得自动改变输入参数。

## 代码结构与数据流

新增一个不导入 `mph` 的轻量共享模块，例如：

```text
comsol_workflow/cladding_shift_profile.py
```

该模块集中负责：

- profile 数据结构、校验、序列化和稳定目录标签；
- `relative_cladding_layer(point, bulk_radius)`；
- `shift_envelope(layer, profile)`；
- `resolved_shift_profile(cladding_layers, target_factor, period, profile)`；
- uniform 与 tanh-power 的纯函数实现。

`parameter_config.py` 负责从 JSON 加载该 profile，并将它加入
`SharedParameters.to_metadata()`；finite 与 strip 必须调用同一共享实现，禁止各自复制
tanh 公式。

几何数据流为：

```text
parameter.json
  -> target shift list + cladding_shift_profile
  -> 当前 case 选择 f_target
  -> point.shell 解析 relative layer
  -> 共享 helper 计算 envelope(layer) 与 f_layer
  -> 现有 normal-fan 计算 side/corner 方向
  -> 形成该 lattice point 的 absolute modulation_shift
  -> 同一 cladding unit cell 的全部 hole polygon 刚性平移
  -> strip 坐标变换或 finite/quarter 几何构建
```

计划涉及：

- `comsol_workflow/cladding_shift_profile.py`：新增纯 profile 实现；
- `scripts/run_main/parameter_config.py`：加载、校验、元数据和目录标签；
- `scripts/parameter.json`：显式启用 `tanh_power/4.0/2.0`；
- `scripts/run_main/run_strip_1d.py`：逐 point 使用层包络；
- `scripts/run_main/run_finite.py`：逐 point 使用同一层包络；
- `scripts/run_main/run_finite_quarter.py`：不复制公式，只验证继承 finite 后的几何和
  metadata；
- 聚焦测试、README、`AGENTS.md` 和 `MEMORY.md`：同步参数语义与维护约束。

## 输出命名、元数据与历史结果隔离

profile 会改变实际几何，必须进入 finite/strip 的系列身份。否则新渐变几何可能落入
已有 uniform 结果的同一 `shift0.XXX` 目录，并被 completion/reuse 逻辑错误复用。

命名规则：

- `uniform` 保持现有系列目录名，保证历史兼容；
- 非 uniform profile 在系列名末尾增加稳定标签，例如：

```text
finite_10-10_<cell-and-mesh-label>_shiftprof-tanhpow-l4-p2
strip1d_10-10_<cell-and-mesh-label>_shiftprof-tanhpow-l4-p2
```

- case 名继续使用 `shift0.XXX` 或 `cut<angle>_shift0.XXX`；这里的 shift 始终是
  `f_target`，不是第一层或平均位移。
- unit-cell 系列名不增加 profile 标签。

每个 case 的 config、summary 和 resume metadata 至少记录：

```text
cladding_inward_shift_factor        # 远端无量纲目标 f_target
cladding_inward_shift               # 兼容字段，语义明确为 A*f_target
cladding_shift_profile              # kind/scale_layers/power/稳定标签
resolved_cladding_shift_by_layer    # layer/envelope/effective factor/side/corner 位移
outer_envelope
bulk_start_envelope
max_first_difference
max_second_difference
```

completion/reuse 与 finite resume 必须比较 profile identity 和解析后的逐层值；缺失或
不一致时拒绝复用，不得只检查原有标量 `cladding_inward_shift`。

已有 `scripts/.out/` 结果不迁移、不覆盖、不改写、不删除。新 profile 使用独立系列
目录，与 uniform 基线并存。

## 几何诊断产物

在现有 case overview/config 布局内增加或更新以下诊断，不创建新的一级输出目录：

- `cladding_shift_profile.csv`：逐层记录 envelope、有效 factor、side/corner 物理位移、
  一阶差分和二阶差分；
- `cladding_shift_profile.png`：同时显示逐层绝对位移及相邻层差分；
- 现有 normal-fan/shift-arrow 图继续绘制每个 point 的实际逐层向量；
- 图中文字改为 target shift、profile 标签和逐层范围，不能继续误写成所有 cladding
  共用一个 side/corner shift；
- config 中保留理想 `LatticePoint` 坐标和实际 `modulation_shift`，便于下游按 shell
  分析。

## 实施步骤

1. 在任何修改前复核工作树和现有正式结果边界，确认只改本计划涉及的源文件、测试
   和文档。
2. 新增共享 profile 数据结构及纯函数，实现 `uniform` 与 `tanh_power`，先以单元测试
   固定数学行为。
3. 扩展共享参数加载、序列化和 series label；保持 unit-cell identity 与输出目录不受
   profile 影响。
4. 将 finite 与 strip 的 cladding shift 幅值改为按 `point.shell` 解析；保持现有
   normal-fan 方向、corner 比例和 strip 整体坐标变换。
5. 让 finite-quarter 只继承 finite 的 profile 实现，并验证 full geometry 镜像对称。
6. 扩展 config、run summary、completion/reuse 和 resume 校验，阻止不同 profile
   之间交叉复用。
7. 增加逐层 CSV/PNG 与更新后的几何标注。
8. 更新共享参数说明、主入口说明和长期维护规则。
9. 完成全部非 COMSOL 验证；实施完成后在本文档回填实际改动、测试结果、已知问题和
   最终状态。

## 自动化验证计划

新增聚焦测试文件，例如 `tests/test_cladding_shift_profile.py`，并补充现有 finite、
strip、parameter-config 测试，至少覆盖：

1. `layer = point.shell - BULK_RADIUS`，cavity 外层为 `0`、第一层 cladding 为 `1`；
2. `tanh_power/4/2` 的已知逐层数值、边界、单调性和有限性；
3. 缺失 profile 时 `uniform` 精确复现当前所有 cladding point 的 polygon 和 shift；
4. target factor 为零时，所有 profile 的实际 shift 都严格为零；
5. 同一 shell 的 envelope 一致，C6 旋转和镜像对应点的向量按对称关系变换；
6. corner/side 幅值比始终为 `2/sqrt(3)`；
7. 同一 cladding 单胞的全部孔获得相同绝对位移，cavity 孔保持零位移；
8. finite 与 strip 对同一未旋转 lattice point 得到相同 shift；strip 旋转后 polygon、
   point 与 modulation vector 保持一致刚性变换；
9. finite-quarter 的 full geometry mirror-symmetry 验证继续通过；
10. 配置类型、未知 kind、非有限/非正 scale、`power<=1` 的拒绝路径；
11. uniform series 名保持不变，tanh profile 得到独立稳定 series 名；
12. completion/reuse 和 resume 在 profile 或逐层解析值不一致时拒绝复用；
13. config、summary、逐层 CSV 的字段、数值和排序；
14. 现有三位小数 `shift0.XXX` case 命名与 trend 横轴语义保持不变。

实施后的非 COMSOL 验证命令计划为：

```powershell
uv run python -m pytest -q `
  tests\test_cladding_shift_profile.py `
  tests\test_run_strip_shift_scan.py `
  tests\test_run_unit_cell_band.py

uv run python -m py_compile `
  comsol_workflow\cladding_shift_profile.py `
  scripts\run_main\parameter_config.py `
  scripts\run_main\run_strip_1d.py `
  scripts\run_main\run_finite.py `
  scripts\run_main\run_finite_quarter.py

git diff --check
```

不得把创建 `SimulationRun`、mesh、COMSOL model 或 eigensolve 当作自动化 smoke test。

## 后续物理验证矩阵（不在本轮自动执行）

代码和纯几何验证通过后，再由用户单独确认是否启动计算。第一轮应固定相同的
`f_target` 比较 profile，而不是同时大范围扫描幅值和函数形状：

```text
uniform
tanh_power(scale_layers=3, power=2)
tanh_power(scale_layers=4, power=2)  # 首选基线
tanh_power(scale_layers=5, power=2)
```

优先使用当前 Pareto 转折区的 `f_target=0.04` 和 `0.05`；若需要函数族对照，再加入
具有严格端点约束的 smoothstep，但不在首版实现中默认启用。

验证顺序：

1. 只生成并检查几何诊断，确认孔间最小距离、无重叠/异常 clipping、对称性和远端
   饱和度；
2. strip 1D 用于筛查 interface 模式、频率、局域长度和带隙是否保持，但不把 strip
   结果当作有限六边形角散射的完整结论；
3. finite-quarter/full finite 比较同一物理模式的 Q、Gamma 权重、频率、局域长度、
   cladding hotspot、Fourier/light-cone 和 far-field 辐射；
4. 跨 profile 比较时不能只固定 `mode_idx`，应结合频率、模式权重和场重叠跟踪同一
   模态；
5. 同时报告 Q 与模体积/局域性，避免把模式过度扩展造成的表观 Q 上升误判为纯散射
   改善。

## 风险与非目标

- 分层包络只平滑 shift 幅值，不消除现有 normal-fan 边区/角区之间的方位方向切换。
  若后续热点集中在六角形角区，再把连续 closest-point direction field 作为独立任务，
  不与本次径向 profile 同时修改。
- cavity 与 cladding 的 `b0/eta/zeta` 仍在 interface 处切换；本任务不渐变单胞紧凑
  参数，以免同时改变局域带结构和界面模形成机制。
- `scale_layers` 太小会把跳变移到中间少数层；太大会使最外层尚未恢复目标 shift。
  因此必须同时检查一阶差分、二阶差分和远端饱和度，不能只查看连续函数曲线。
- profile 可能改变孔间距、局域晶格常数、模式频率和模体积；Q 不保证随渐变长度单调
  提高。
- 不改现有 trend 的横轴定义；横轴继续表示远端目标 `f_target`。跨 profile 的结果
  必须按 profile 分系列比较，不能混合为一条单变量 shift 曲线。
- 本任务不启动参数扫描、不自动选择最佳 `scale_layers/power`、不修改既有正式结果。

## 回滚方案

若实现需要回滚：

1. 恢复 finite/strip 使用 uniform 目标幅值的原路径；
2. 删除新增共享 profile 模块及对应测试和文档改动；
3. 从共享配置移除可选 `cladding_shift_profile`；缺失字段的 uniform 兼容行为保证旧
   配置仍可读取；
4. 恢复非 uniform series label 扩展，但不得移动、覆盖或删除已经生成的独立 profile
   结果；
5. 已有 uniform 结果和 unit-cell 结果始终保持原位，不参与回滚写操作。

## 执行记录（实施后填写）

### 实际改动

- 新增 `comsol_workflow/cladding_shift_profile.py`，集中实现 `uniform` 与
  `tanh_power` 包络、相对层编号、逐层绝对位移解析、远端饱和预检及 CSV/PNG
  几何诊断。
- 在 `scripts/parameter.json` 启用已批准的默认值：
  `tanh_power(scale_layers=4.0, power=2.0)`；配置缺失该字段时仍回退到历史
  `uniform` 行为。
- `parameter_config.py` 负责统一解析 profile，并将非 uniform profile 的稳定标签加入
  finite/strip 系列目录；unit-cell identity 不包含 profile。
- `run_finite.py` 与 `run_strip_1d.py` 现在都按
  `layer=point.shell-BULK_RADIUS` 解析包络，并将它乘到原 normal-fan 位移幅值上；
  位移始终相对于理想晶格绝对施加，corner/side 比例保持 `2/sqrt(3)`。
- `run_finite_quarter.py` 继续复用 finite 完整几何，并同步保存 profile 诊断与 summary
  metadata。
- 扩展 config、summary、finite resume 与 strip completion/reuse 契约，使不同 profile
  或不同逐层解析值不能交叉复用；原 `shift0.XXX` case 名仍表示远端目标值。
- 更新根目录与 `scripts/` 使用说明、长期维护约束和聚焦测试。

### 验证结果

- 聚焦回归：
  `uv run python -m pytest -q tests/test_cladding_shift_profile.py tests/test_run_strip_shift_scan.py tests/test_run_unit_cell_band.py`
  通过，结果为 `129 passed`。
- 全套非 COMSOL 回归：`uv run python -m pytest -q` 通过，结果为
  `268 passed, 3 warnings`。3 条 warning 来自既有
  `tests/test_strip_1d_fourier_quotient.py` 测试函数返回非 `None`，与本改动无关。
- Python 编译检查通过：共享 profile 模块、配置模块及 finite/quarter/strip 三个入口
  均可 `py_compile`。
- 纯几何检查在目标 factor `0.09` 下验证 finite 与 strip 各生成 `7566` 条孔记录，
  相同 lattice point 的最大位移向量差为 `0`；finite-quarter 完整几何镜像对称检查
  通过。
- `git diff --check` 通过。

### 计算状态与遗留问题

代码与非 COMSOL 验证已经完成；未启动 COMSOL，未生成新的正式求解结果，也未迁移、
覆盖或删除 `scripts/.out/` 中的已有结果。profile 对 Q、局域长度、Fourier/light-cone
权重、far-field 辐射及边界散射的实际改善仍需按“后续物理验证矩阵”单独计算验证，
不能由几何平滑性直接推断。
