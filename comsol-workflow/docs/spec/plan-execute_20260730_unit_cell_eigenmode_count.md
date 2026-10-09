# unit-cell 本征模求解数量共享参数化

日期：2026-07-30

## 实施前方案

### 目标与范围

把 unit-cell 工作流请求 COMSOL 求解的本征模数量统一设为 6，并从
`scripts/parameter.json` 读取，使 v1、v2 及复用其 solver 配置的 unit-cell 扫描
不再各自硬编码数量。

本次新增顶层标准 JSON 字段：

```json
"unit_cell_eigenmode_count": 6
```

finite、finite-quarter 和 strip 的 `EIGENMODE_COUNT` 用途及计算规模不同，继续保留
各自主程序配置，不纳入本次统一。

### 数据流

```text
scripts/parameter.json
  -> parameter_config.load_shared_parameters()
  -> SharedParameters.unit_cell_eigenmode_count
  -> run_unit_cell_band v1/v2 的 EIGENMODE_COUNT
  -> 依赖 band.EIGENMODE_COUNT 的 unit-cell sweep SimulationConfig
```

`optimize_cladding_p2_d2_alignment.py` 当前单独硬编码 12，将改为直接读取同一个
`ACTIVE_PARAMETERS.unit_cell_eigenmode_count`。

### 校验与缓存兼容

- `unit_cell_eigenmode_count` 必须为正整数；布尔值、浮点数、0 和负数均拒绝。
- 将该字段加入共享 metadata 和 unit-cell 运行期间的参数身份检查；运行期间被改动
  时，不允许用旧参数结果回写中心频率。
- 新增 `unit_cell_case_parameter_label`，在原 cell/mesh 标签后增加稳定
  `_modes6` 后缀。v1/v2 使用该标签生成新目录，避免 6-mode 错误复用现有 4/12-mode
  正式缓存。
- finite/strip 的 `case_parameter_label` 与系列目录不增加 unit-cell mode 后缀，避免
  unit-cell 求解数量改变无关地切换 finite/strip 系列。
- unit-cell 扫描中自行维护缓存、且输出名未包含求解数量的入口，将增加 modes 后缀
  或改用 unit-cell 专属标签。

### 兼容与回滚

- v1/v2 的求解、模式分析和后处理代码不改变，只替换参数来源和输出身份。
- `parameter.json` 保持无注释标准 JSON。
- 回滚时移除共享字段与解析字段，并把各 unit-cell 常量恢复为局部值即可；已有正式
  结果不移动、不覆盖、不删除。

### 验证计划

- 验证共享文件解析值为 6，v1/v2 及 unit-cell 扫描配置均为 6。
- 验证缺失/非法字段被明确拒绝，metadata 与 identity 包含该值。
- 验证 v1/v2 独立输出标签包含 `_modes6`，finite/strip 系列标签保持原结构。
- 运行 parameter/unit-cell/sweep 聚焦测试、`py_compile` 和 `git diff --check`。
- 本次只改变参数传递，不启动新的 COMSOL eigensolve。

## 实施记录与结果

### 实际改动

- 在 `scripts/parameter.json` 新增
  `"unit_cell_eigenmode_count": 6`。
- `SharedParameters` 新增同名字段；loader 要求其为正整数，并将其写入 identity、
  metadata 和 unit-cell 运行期间的参数一致性检查。
- 新增 `unit_cell_case_parameter_label`，当前值为
  `cav(245-0.96-1.156)_clad(242-0.986-0.917)_mesh5_modes6`；原
  `case_parameter_label` 保持不变，finite/strip 系列名未受影响。
- v1、v2 的 `EIGENMODE_COUNT` 均改为
  `ACTIVE_PARAMETERS.unit_cell_eigenmode_count`，独立输出目录分别为
  `unitcell_band_<unit-cell-label>` 和 `unitcell_band_v2_<unit-cell-label>`。
- `scan_unit_cell_p2_gamma.py`、`cladding_bandgap_alignment.py`、
  `optimize_cavity_cladding_cells.py`、`optimize_cladding_doublet_alignment.py` 的
  mode-dependent 输出/源目录增加 `modes6` 身份。
- `optimize_cladding_p2_d2_alignment.py` 移除本地硬编码 12，改读共享值，并在输出
  目录加入 `modes6`。
- README 和参数/入口/扫描测试同步更新。

### 验证结果

- 运行时核对：shared、v1、v2、p2/d2 alignment 的值均为 6；其他 unit-cell
  扫描通过 `band.EIGENMODE_COUNT` 使用同一值。
- v1/v2/parameter/sweep/finite-quarter 相关聚焦测试：`112 passed`。
- 旧式缺省 shift-profile 参数文件兼容测试：`1 passed`。
- `py_compile`：所有本次修改的入口通过。
- `git diff --check`：通过。
- 另行运行整个 `test_cladding_shift_profile.py` 时有 4 个既有配置相关失败：测试
  假定 `tanh_power` 和非零 shift，而当前正式参数为 `uniform`、shift=0；与本次
  eigenmode count 参数化无关，未在本任务中改动该行为。
- 按计划未启动 COMSOL，也未写入、覆盖或删除正式计算结果。

### 最终状态

unit-cell 求解数已经统一由 `parameter.json` 控制，当前为 6。改变该字段后重新
启动 Python 入口即可读取新值；unit-cell 主流程输出标签随数值变化，避免跨求解数
复用缓存。finite/strip 的自身求解数量及系列目录保持原状。
## 2026-07-30 命名回滚补充

用户要求撤回输出目录中的 mode-count 后缀。实施范围如下：

- `unit_cell_eigenmode_count` 继续作为 `parameter.json` 中的可调求解参数，并继续写入
  metadata；不回滚参数化本身。
- unit-cell 主入口恢复使用 `case_parameter_label`，可见系列目录以 `_meshN` 结束，
  不再追加 `_modesN`。
- 同次参数化中为 unit-cell sweep、优化和缓存目录加入的 `modesN` 后缀一并移除。
- 为避免无后缀目录跨求解数量误用缓存，每个新求解 k 点写入 solver identity；缓存
  复用时核对 `eigenmode_count`。不符合或缺少 identity 的点视为需要重算。
- 已有 `_modesN` 正式结果目录保持原样，不移动、不覆盖、不删除。

验证与回滚：更新动态命名/缓存测试，运行 unit-cell 与相关 sweep 聚焦测试、
`py_compile` 和 `git diff --check`。如需回滚，只恢复可见命名表达式和测试；正式结果
不参与回滚操作。

### 命名回滚执行结果

- 主入口实际目录已恢复为
  `unitcell_band_cav(... )_clad(... )_mesh5` 形式，不含 `_modes6`。
- unit-cell scan、cell-pair optimization、doublet alignment 和 p2/d2 alignment 的
  `modesN` 后缀均已移除。
- 新求解点在 k-point 目录写入 `solver_config.json`；复用前核对求解数、shift、mesh、
  mode type 和 wavelength。各优化器的 candidate summary 和已接受 cavity 缓存也核对
  `eigenmode_count`，求解数变化不会误用同名目录中的不兼容缓存。
- 相关命名与缓存合同测试：13 passed；扩大聚焦范围结果为 108 passed、2 failed，两个
  failure 均来自工作区当前频率窗口为 175–220 THz、而既有测试仍断言 185–205 THz，
  与本次命名/缓存改动无关，未擅自回滚该窗口。
- 7 个相关 Python 入口 `py_compile` 通过，`git diff --check` 通过。
- 未启动 COMSOL；未移动、覆盖或删除已有 `_modesN` 正式结果目录。
