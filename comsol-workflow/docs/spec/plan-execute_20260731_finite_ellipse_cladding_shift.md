# Finite cladding 径向椭圆 shift 执行计划

## 方案概述

在现有 finite / finite-quarter 的 cladding x/y 分组平移之外，增加一个可选的
径向平滑模式。配置使用：

```json
"cladding_shift_geometry": {"kind": "groups"}
```

或：

```json
"cladding_shift_geometry": {"kind": "ellipse"}
```

`groups` 精确保留当前 12 个 normal-fan 区域的 x/y 分组、side/corner 方向和 corner
幅值补偿。`ellipse` 不再按扇区决定方向；每个 cladding cell 从理想晶格中心
`p=(x,y)` 沿指向 `(0,0)` 的径向作刚性平移。

设 `theta=atan2(y,x)`，x/y 目标参数分别为 `f_x`、`f_y`，层包络为 `g(L)`，
晶格周期为 `A`。`ellipse` 的平滑角度权重和位移定义为：

```text
f(theta) = f_x cos(theta)^2 + f_y sin(theta)^2
delta_p = -A g(L) f(theta) p / |p|
```

因此 0/180 度位置的位移大小严格为 `A*f_x`，90/270 度位置严格为 `A*f_y`，
主轴之间连续过渡。这里的兼容基准是用户确认的 `A*factor` 参数定义；ellipse 不继承
`groups` 在 x 轴 corner 上使用的 `2/sqrt(3)` 补偿，因此两种模式在 x 轴的实际物理
位移不相等。该定义允许 `f_x` 或 `f_y` 单独为零，但两个 factor 都必须为有限非负数。
两个 factor 同时为零时所有 cladding cell 保持理想位置。

## 范围

- 修改共享参数加载和 cladding-shift 纯函数模块。
- 修改 `scripts/run_main/run_finite.py`。
- 修改 `scripts/run_main/run_finite_quarter.py`，继续复用 full finite 几何。
- 更新 `scripts/parameter.json`，启用 `ellipse`，初值为 x=0.05、y=0.08，profile
  保持 `uniform`。
- 增加 finite / finite-quarter 的配置、几何、对称性、命名、metadata、resume 和
  诊断测试。
- 不修改 strip 的 shift 几何、扫描列表、case 命名或输出系列。
- 不启动 COMSOL，不改写或删除已有计算结果。

## 数据流

```text
scripts/parameter.json
  -> SharedParameters.cladding_shift_geometry (missing => groups)
  -> run_finite.CLADDING_SHIFT_GEOMETRY
  -> cladding_inward_shift(point)
       groups  -> normal_fan region/group shift（现有实现）
       ellipse -> ideal point angle weight + inward radial shift（新实现）
  -> full cladding hole records（cell 内所有孔共享同一 translation）
  -> finite full / finite quarter geometry
  -> config.json、profile diagnostics、geometry diagnostics
```

## 配置与命名

- `cladding_shift_geometry.kind` 仅允许 `groups` 和 `ellipse`。
- 旧参数文件缺少 `cladding_shift_geometry` 时回退到 `groups`。
- `groups` 系列名保持当前形式。
- `ellipse` 的 finite 和 finite-quarter 系列名末尾增加稳定后缀
  `_ellipse-shift`。
- case 名继续为 `shiftx0.050_shifty0.080`。
- strip 的系列名不增加该后缀，避免 finite 专用参数改变 strip 输出位置。

## Metadata 与诊断

- 共同记录 `cladding_shift_geometry`、x/y factor、profile 和公式版本。
- `groups` 保留 region-group、side/corner 补偿和逐层记录。
- `ellipse` 记录 x/y 主轴逐层目标、角度权重公式、原点和每个 cell 的理想中心、
  角度、有效 factor 与位移向量；不写入虚假的 corner 补偿信息。
- resume 契约包含 geometry kind 和对应的解析 metadata，禁止用 `groups` checkpoint
  恢复 `ellipse`，反之亦然。
- 几何图根据模式显示分组扇区或 ellipse 的连续径向位移。

## 验证

1. 参数解析：合法 kind、未知字段、缺省 `groups`、非负 factor、单轴零值。
2. 命名：只有 finite / finite-quarter 的 `ellipse` 系列增加 `_ellipse-shift`。
3. 几何：0/90 度主轴锚点、45 度权重、所有向量指向原点、零轴行为。
4. 对称性：关于 x/y 轴镜像后的位移满足 C2v；quarter 裁切检查继续通过。
5. 刚性：同一 cell 的所有孔使用同一位移向量。
6. 兼容：`groups` 的现有 12 区域和 corner 补偿测试继续通过；strip 测试不受影响。
7. 运行聚焦 pytest、`py_compile` 和 `git diff --check`；不使用几何构建或
   eigensolve 作为 smoke test。

## 回滚方案

把参数改为：

```json
"cladding_shift_geometry": {"kind": "groups"}
```

即可恢复当前 12 区域行为和原系列命名。代码回滚时可移除 ellipse 分支、专属
metadata/诊断和系列后缀；缺省 `groups` 兼容不依赖修改旧结果。

## 执行记录

### 实际改动

- 新增 `CladdingShiftGeometry`，只接受 `groups`/`ellipse`；字段缺失时回退到
  `groups`，未知字段和未知 kind 明确报错。
- x/y finite factor 改为有限非负数校验，允许任一或两者为 0。
- `ellipse_shift_weight` 与 `ellipse_inward_shift` 实现确认公式；0/180 度严格为
  `A*fx`，90/270 度严格为 `A*fy`，不继承 groups corner 补偿。
- `run_finite.py` 按 geometry kind 分派旧 normal-fan groups 或新连续径向位移；
  每个 cell 内所有孔继续共享同一刚性平移向量。
- finite-quarter 继续从 full finite records 裁切，因 ellipse 权重只含 x/y 偶次项而
  保持 x/y 镜像对称。
- groups 保留原逐层 group metadata/诊断；ellipse 新增主轴逐层 metadata、公式版本、
  每 cell 的理想坐标/角度/有效 factor/位移向量，以及不含 corner 字段的 CSV/PNG。
- resume 契约按 geometry kind 比较对应解析值；旧 metadata 缺失 geometry 时按
  groups 归一化。
- ellipse 的 full/quarter 系列名增加 `_ellipse-shift`；case 名仍为
  `shiftx0.050_shifty0.080`；strip 系列名不变。
- ellipse 模型级第三张几何图改为 `full_lattice_ellipse_shift.png`；groups 继续使用
  `full_lattice_normal_fan_regions.png`。
- `parameter.json` 已设为 x=0.05、y=0.08、geometry=ellipse、profile=uniform。
- 同步更新 `AGENTS.md`、`MEMORY.md`、根 README 和 scripts README 的长期配置、
  兼容、命名与诊断规则。

### 验证结果

- 聚焦测试：`50 passed`。
- unit-cell/strip/quarter 相关回归：`138 passed`。
- 全套测试：`349 passed, 3 warnings`；3 个 warning 均为既有
  `test_strip_1d_fourier_quotient.py` 测试函数返回非 `None`。
- `py_compile`：四个修改过的核心 Python 文件通过。
- `git diff --check`：通过，仅报告 Windows 工作树既有 LF/CRLF 提示。
- 未启动 COMSOL，未创建或改写 `scripts/.out/` 计算结果。

### 遗留问题与最终状态

- `ellipse` 是用户确认的平滑 `cos^2/sin^2` 径向权重，不是严格椭圆极径公式；
  这是为了自然支持单轴 factor 为 0。
- groups 的 x 轴 cell 属于 corner，实际位移含 `2/sqrt(3)`；ellipse 明确以
  `A*factor` 为兼容基准，因此两种模式在 x 轴的实际位移不同。
- 实施与非 COMSOL 验证完成，可以进入后续 geometry-only 预览或正式 COMSOL 计算，
  但两者均需用户另行明确授权。
