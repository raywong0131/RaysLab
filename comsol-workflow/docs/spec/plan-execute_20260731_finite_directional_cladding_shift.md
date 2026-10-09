# Finite 分区 Cladding Shift：计划与执行记录

日期：2026-07-31
状态：已完成（代码与非 COMSOL 验证）

## 目标

将 finite 与 finite-quarter 的 12 个 normal-fan 边界区域分成两组，并从
`scripts/parameter.json` 分别读取 shift factor：

- `cladding_y_shift_factor`：上下 6 区，包括 top/bottom 两个水平 side 和四个上下
  corner；
- `cladding_x_shift_factor`：左右 6 区，包括四个斜 side 和最左/最右两个 corner。

本任务只改变每个 normal-fan 区域选择的目标幅值。现有 inward direction、corner
幅值补偿 `2/sqrt(3)`、逐层 `cladding_shift_profile` 包络和“相对理想晶格的绝对位移”
语义保持不变。

## 精确分区

沿用 `ideal_bulk_hex_vertices()` 的索引：右、右上、左上、左、左下、右下为 corner
`0..5`；side `i` 连接 corner `i` 与 `i+1`。

```text
y group:
  side   = {1, 4}
  corner = {1, 2, 4, 5}

x group:
  side   = {0, 2, 3, 5}
  corner = {0, 3}
```

这里的 x/y 表示左右区域组和上下区域组，不表示把位移向量强制投影到笛卡尔 x/y
轴。每个区域仍沿原 normal-fan inward direction 移动。

## 参数与兼容

新增两个 finite 专用共享标量：

```json
"cladding_x_shift_factor": 0.05,
"cladding_y_shift_factor": 0.12
```

实施时采用用户明确指定的初值 `x=0.05`、`y=0.12`。`cladding_shift_factors`
继续供 strip 扫描使用，本任务不改变 strip 几何或输出。

旧参数文件缺少 x/y 字段时，两者都回退到 `cladding_shift_factors[0]`，从而复现原
finite uniform-factor 行为。显式字段必须是有限数值。

Unit-cell identity 不包含 x/y shift，因此这些字段变化不得影响 unit-cell 输出目录或
中心频率安全回写。

## 数据流与命名

```text
parameter.json
  -> SharedParameters.cladding_x/y_shift_factor
  -> run_finite normal_fan_feature(kind, idx)
  -> cladding_shift_group(kind, idx) = x | y
  -> target factor * layer profile envelope
  -> existing normal-fan direction and corner compensation
  -> finite-quarter reuses the full finite geometry
```

新的 finite/quarter case 名必须同时编码两个 factor：

```text
shiftx0.050_shifty0.120
```

不得落入历史 `shift0.XXX` case。Series 目录仍由 cell、层数、mesh 和 layer profile
决定；x/y factor 属于 case identity。

Config 和 run summary 至少记录：

```text
cladding_x_shift_factor
cladding_y_shift_factor
cladding_x_inward_shift
cladding_y_inward_shift
cladding_shift_region_groups
resolved_cladding_shift_by_group.x
resolved_cladding_shift_by_group.y
```

finite resume 必须比较 x/y factor、分组契约及两组逐层解析值。诊断 CSV/PNG 同时显示
x/y 两组，避免继续写成单一 target。

## 修改范围

- `scripts/parameter.json`
- `scripts/run_main/parameter_config.py`
- `comsol_workflow/cladding_shift_profile.py`
- `scripts/run_main/run_finite.py`
- `scripts/run_main/run_finite_quarter.py`
- 相关测试、README、`AGENTS.md` 和 `MEMORY.md`

不修改 `run_strip_1d.py` 的 shift 计算，不启动 COMSOL，不移动、覆盖或删除已有结果。
保留工作树中与 strip output layout 等其他任务有关的现有修改。

## 验证方案

1. 参数加载：显式 x/y、旧文件回退、非有限值拒绝、unit-cell identity 不变；
2. 分区映射：2 side + 4 corner 属于 y，4 side + 2 corner 属于 x，12 区无遗漏重叠；
3. 幅值：各组 side 使用 `A*f_group*envelope`，corner 保持 `2/sqrt(3)`；
4. 几何：同一 cell 内孔保持刚性位移，cavity 不移，x/y 轴镜像对称保持；
5. 命名：case 同时包含 x/y，且不等于旧 `shiftX.XXX`；
6. metadata、resume 和诊断产物同时包含两组；
7. finite-quarter 复用同一 full geometry，四种 PEC/PMC case 不改变几何；
8. 运行相关 pytest、全套非 COMSOL pytest、`py_compile` 与 `git diff --check`。

## 回滚

移除 x/y 参数字段和 directional metadata；finite/quarter 恢复遍历
`cladding_shift_factors` 并使用单一 `CLADDING_INWARD_SHIFT_FACTOR`。不得改写已经生成的
directional case 或历史 uniform case。

## 执行记录

### 实际改动

- `parameter.json` 新增 finite/quarter 专用
  `cladding_x_shift_factor=0.05`、`cladding_y_shift_factor=0.12`；原
  `cladding_shift_factors` 保留给 strip 扫描。
- `SharedParameters` 解析、验证并序列化 x/y factor；旧配置缺失字段时两者回退到
  `cladding_shift_factors[0]`，unit-cell identity 不包含它们。
- 共享 profile 模块新增精确 12 区分组、x/y 双 profile metadata 和双组 CSV/PNG
  诊断；原 strip 单 factor 接口保持不变。
- finite 根据 normal-fan 的 `(kind, idx)` 选择 x 或 y factor，再乘逐层 profile
  envelope；normal-fan 方向、corner `2/sqrt(3)` 补偿和绝对位移语义保持不变。
- finite main 从三点 scalar 扫描改为运行一个显式 x/y case，case 名为
  `shiftx0.050_shifty0.120`；metadata、图注、region 图和 resume 契约同步为双组。
- finite-quarter 继续复用 finite 完整几何，所有 symmetry case 使用相同 x/y factor，
  并沿用 directional case 名与双组 metadata/诊断。
- 更新 README、维护规则、长期记忆和相关测试；未修改 strip 几何实现。

### 验证结果

- 聚焦测试通过：
  `tests/test_cladding_shift_profile.py`、`tests/test_finite_cavity_field_exports.py`、
  `tests/test_run_unit_cell_band.py`、`tests/test_run_strip_shift_scan.py` 共
  `154 passed`。
- 纯几何检查确认 12 区完整且无重叠：x/y 各 6 区。当前 uniform profile 下实际幅值
  为 x-side `0.041 um`、x-corner `0.047342722074 um`、y-side `0.0984 um`、
  y-corner `0.113622532977 um`，符合 `A=0.82 um` 和 corner 补偿。
- Full finite 生成 `7566` 条孔记录，quarter 截取 `1912` 条；完整几何 x/y 轴镜像
  对称检查通过。
- 新 finite 和当前首个 quarter 输出路径均不存在；新命名不会落入历史 scalar case。
- 参数 JSON 校验、相关模块 `py_compile` 和任务文件 `git diff --check` 通过。
- 全套非 COMSOL pytest 结果为 `331 passed, 4 failed, 3 warnings`。4 个失败全部来自
  当前工作树中既有的 strip `STRIP_CUT_ANGLE=0` 与仍硬编码 `cut_60` 的旧测试：
  `test_analyze_strip_gamma_summary.py` 1 个、`test_run_strip_bulk_rotation.py` 1 个、
  `test_run_strip_output_layout.py` 2 个；与本任务的 finite/quarter 改动无关。本任务未
  擅自改写正在进行的 strip output-layout 工作。

### 计算状态与遗留问题

未启动 COMSOL，未生成新的正式结果，也未覆盖或删除已有结果。现有
`scripts/analysis/finite_trend.py` 仍按单变量 `shiftX.XXX` 组织趋势，不应直接把一个
directional x/y case 当作旧 scalar shift 趋势；若后续扫描 x 或 y，需要先明确固定轴与
横轴语义，再单独扩展趋势分析。
