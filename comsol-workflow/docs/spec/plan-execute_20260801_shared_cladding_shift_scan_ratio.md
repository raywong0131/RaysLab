# Finite 复用 cladding shift 扫描与 y/x 比例执行计划

## 目标

默认去除 `parameter.json` 中重复的 finite 专用 x/y shift 输入，让 strip、finite 和
finite-quarter 共同复用：

```json
"cladding_shift_factors": [0.00, 0.01, ..., 0.09]
```

finite/quarter 只额外读取一个比例：

```json
"cladding_y_over_x_shift_ratio": 1.6
```

对列表中的每个扫描值 `f`，解析为：

```text
cladding_x_shift_factor = f
cladding_y_shift_factor = f * cladding_y_over_x_shift_ratio
```

因此 `f=0.05` 时仍得到 `x=0.05`、`y=0.08`。同时保留非默认的显式 pair 输入：

```json
"cladding_x_shift_factor": 0.0,
"cladding_y_shift_factor": 0.08
```

显式 pair 模式只运行这一组 finite/quarter case，可表示 ratio 模式无法表达的
`x=0,y>0`。现有 `groups`/`ellipse` 几何、逐层 profile 和 `A*factor` 语义不改变。

## 范围与数据流

```text
parameter.json
  -> 默认：cladding_shift_factors + cladding_y_over_x_shift_ratio
     或非默认：cladding_x_shift_factor + cladding_y_shift_factor
  -> SharedParameters.finite_cladding_shift_factor_pairs
  -> run_finite / run_finite_quarter 逐 pair 扫描
  -> shiftxX.XXX_shiftyY.YYY case
```

- strip 继续逐项使用 `cladding_shift_factors`，行为与命名不变。
- finite 和 finite-quarter 改为逐项运行同一列表派生出的 x/y pair。
- `cladding_shift_geometry`、`cladding_shift_profile` 和 `_ellipse-shift` 系列后缀不变。
- 不启动 COMSOL，不写入或删除已有计算结果。

## 参数与兼容

- 新比例必须是有限非负数；0 合法，表示所有 finite y 主轴目标为 0。
- `cladding_shift_factors` 现在同时服务 strip 和 finite，因此每一项必须有限非负。
- 默认配置删除 `cladding_x_shift_factor` 和 `cladding_y_shift_factor`，避免三份目标输入。
- 若比例缺失但 x/y 字段同时存在，则选择 `xy` 输入模式，finite/quarter 只运行该
  pair；两个值都必须有限非负，且任一个均可为 0。
- 比例与任一显式 x/y 字段同时出现时直接报冲突，不允许静默选择优先级。
- 只出现一个显式 x/y 字段时报错，必须成对提供。
- 若比例和显式字段都缺失，回退 `scan_ratio` 模式和比例 1，从而保留历史各向同性扫描。

## 命名、metadata 与恢复

- 系列目录保持当前命名；ellipse 仍增加 `_ellipse-shift`。
- 每个 case 继续使用实际解析值命名，例如
  `shiftx0.050_shifty0.080`，不同扫描点不会互相覆盖。
- 共享 metadata 记录输入模式、原始 factor 列表、可选 y/x 比例和完整解析 pair。
- case metadata 继续记录实际 x/y factor；ratio 扫描时增加扫描基值和比例，显式 pair
  则记录输入模式并允许比例为 `null`。
- resume 比较共享扫描配置以及当前 case 的实际 x/y 值，禁止比例或列表变化后误用旧
  checkpoint。

## 验证

1. 配置加载：默认比例、显式比例、显式 x/y、互斥冲突、缺失配对、非法 factor。
2. 派生值：比例 1.6 时 0.05 映射到 0.08，0 映射到 0。
3. finite/quarter main：按完整 factor 列表调用，顺序稳定，汇总包含全部 case。
4. case 与系列命名：实际 x/y 值进入 case，ellipse 后缀不变，strip 不受影响。
5. geometry：派生后的每个 pair 继续走现有 groups/ellipse 实现。
6. 运行聚焦测试、全套 pytest、`py_compile` 和 `git diff --check`。

## 回滚

移除 ratio 扫描分支并让 finite/quarter main 只运行显式 pair；不需要移动或删除任何
既有结果目录。

## 执行记录

### 实际改动

- `SharedParameters` 新增互斥的 finite shift 输入解析：
  - `scan_ratio`：使用 `cladding_shift_factors` 与
    `cladding_y_over_x_shift_ratio` 派生完整 x/y pair 列表；
  - `xy`：ratio 缺失且 x/y 字段成对出现时，只解析一个显式 pair。
- ratio 与显式 x/y 同时出现、只出现一个 x/y、负数或非有限数都会在参数加载阶段
  报错；显式 pair 支持 `x=0,y>0`、`x>0,y=0` 和 `x=y=0`。
- `run_finite.main()` 逐项运行 `finite_cladding_shift_factor_pairs`；
  `run_finite_quarter.main()` 对每个 pair 展开所有指定 symmetry case。
- finite preset 同时支持结构化 input metadata、顶层 ratio 扫描输入、历史
  `cladding_inward_shift_factors` 以及显式 x/y pair。
- case metadata 和 summary 记录实际 x/y、输入 kind、可选扫描基值及可选 ratio；resume
  契约同步包含这些字段，并为旧显式 x/y metadata 填充兼容默认值。
- `parameter.json` 已删除重复 x/y，当前默认设置为
  `cladding_y_over_x_shift_ratio=1.6`，保留原 `cladding_shift_factors` 列表。
- strip 代码和扫描规则未改；groups/ellipse、uniform/profile、case 名和
  `_ellipse-shift` 系列后缀均保持原定义。
- 同步更新 `AGENTS.md`、`MEMORY.md`、根 README 与 scripts README。

### 当前解析结果

当前 10 个 finite/quarter pair 为：

```text
(0.000, 0.000)  (0.010, 0.016)  (0.020, 0.032)
(0.030, 0.048)  (0.040, 0.064)  (0.050, 0.080)
(0.060, 0.096)  (0.070, 0.112)  (0.080, 0.128)
(0.090, 0.144)
```

系列目录仍为
`finite_10-10_cav(245-0.96-1.156)_clad(240.1-0.967-1)_mesh5_ellipse-shift`；
各 case 使用实际 pair，例如 `shiftx0.050_shifty0.080`。

### 验证结果

- factor/ratio/显式 pair、finite/quarter main、几何和相关回归聚焦测试：
  `189 passed`（后续新增的 preset 显式单轴零值测试已包含在全套测试）。
- 全套测试：`358 passed, 3 warnings`；warning 均为既有
  `test_strip_1d_fourier_quotient.py` 测试返回非 `None`。
- 四个核心 Python 文件通过 `py_compile`。
- `git diff --check` 通过，仅有 Windows 工作树 LF/CRLF 提示。
- 未启动 COMSOL，未创建、覆盖或删除 `scripts/.out/` 结果。

### 最终状态与使用规则

- 默认扫描：保留 ratio，删除显式 x/y。
- 单 pair：删除 ratio，同时设置两个显式 x/y；strip 仍独立扫描 factor 列表。
- 由于输出层已有防覆盖检查，命中的既有 case 不会被静默覆盖；是否 resume 仍由
  `RESUME_FROM_EXISTING_MPH` 明确控制。
- 实施和非 COMSOL 验证完成。
