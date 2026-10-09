# Cladding 简并双重态对齐优化：计划与执行记录

日期：2026-07-26

## 方案概述（实施前）

本任务保持当前 cavity cell 不变，构造并优化满足 `zeta = 1`、`eta < 1`
的 cladding cell。用户侧目标定义为：

1. cladding 的简并 `p1/p2@Gamma` 对齐 cavity 的 `p1@Gamma`；
2. cladding 的简并 `d1/d2@Gamma` 对齐 cavity 的 `d2@Gamma`；
3. 仅优化 cladding 的 `b0` 与 `eta`，cavity 参数不变；
4. 最终使用完整 unit-cell band 计算验证候选参数。

实施前排查发现，三个主工作流的 Fourier 参数转换会把绝对值低于固定阈值的
系数强制截断。该行为会把 `zeta = 1` 对应的 `b_square_f3 = 0` 改成
`0.03`，破坏六个三角形等边和 Gamma 点的严格简并。用户已明确要求取消整个
程序的最小截断限制。因此，本次先全局移除该限制，再实施参数优化。

## 实施范围

- 修改 unit-cell、strip 1D 和 finite 主程序中的 `get_hole_params()`，使所有
  Fourier 系数按输入原值传递；finite-quarter 复用 finite，因此自动继承。
- 更新现有测试，验证零值和小系数不再被改写，并验证 `zeta = 1` 产生六个
  相同边长的三角形。
- 更新项目长期规则，删除“保留最小截断”的旧约束。
- 新增已获用户批准的 cladding 双重态对齐扫描脚本，复用现有 unit-cell
  Gamma 点求解、模式分解和完整能带验证接口。

## 数据流

```text
共享 cavity 参数 + 固定 cladding zeta=1
  -> 精确 Fourier 参数转换（不截断）
  -> Gamma 点候选求解与 p/d 模式选择
  -> b0/eta 自适应扫描
  -> 双目标误差排序
  -> 最优 cladding 参数
  -> 完整 unit-cell band 验证
```

## 兼容与风险控制

- 紧凑参数映射公式保持不变；只取消对 Fourier 系数的人为最小绝对值改写。
- 不修改已有计算结果，不启动 strip 或 finite 计算。
- 在启动昂贵 COMSOL 扫描前，重新核对 cavity 目标频率、扫描范围、mesh、输出
  目录和缓存行为。
- 若精确零系数暴露底层数值接口问题，应停止计算并报告，不恢复隐式截断。

## 验证与回滚

- 运行聚焦 pytest、`py_compile` 和 `git diff --check`。
- 数值验证 `zeta=1` 时六条边长一致，且显式的小 Fourier 系数原样生效。
- 回滚时可独立恢复本记录列出的源文件和测试；不触碰既有输出数据。

## 实际改动与结果（实施后填写）

### 已完成改动

- 删除 `run_unit_cell_band.py`、`run_strip_1d.py` 和 `run_finite.py` 中的
  `clip_min_abs()`，四类 Fourier 系数均按输入原值传递；quarter 工作流通过
  finite 继承该行为。
- 更新几何预览测试中的同源转换，避免测试继续模拟已废止的截断规则。
- 更新 `MEMORY.md`：精确零和小非零系数不得再被改写，`zeta=1` 必须保持
  `b_square_f3=0` 和六个等边三角形。
- 新增 `scripts/run_sweep/optimize_cladding_doublet_alignment.py`。脚本固定当前
  cavity 和 cladding `zeta=1`，在批准范围内以 `b0=0.1 nm`、`eta=0.001`
  的最终精度进行受限 Jacobian 自适应搜索和局部精修。
- 优化目标使用 cladding p/d 双重态的平均频率；同时记录两组数值劈裂，并仅在
  两个对齐误差和两个劈裂均不超过 `0.02 THz` 时接受候选。
- 扫描完成后自动执行完整 cavity/cladding unit-cell band 验证；成功后删除
  中间候选缓存，保留汇总、最优参数和完整验证结果。

### 验证结果

- `zeta=1` 的三条活动几何转换均得到六条 `243.7 nm` 的相同边长；旧截断
  造成的约 `10.888 nm` 边长差已消失。
- 小的 `r/theta/b_square/phi` Fourier 系数原样进入变换。
- 聚焦及相关回归测试：`97 passed`。
- `py_compile` 通过；`git diff --check` 无格式错误。

### 计算状态

- 2026-07-26 已启动 mesh9 自适应 Gamma 点扫描。
- 首个候选 `b0=243.7 nm, eta=0.980, zeta=1` 已完成：p 双重态平均频率
  `198.107881530 THz`，d 双重态平均频率 `200.623844409 THz`；相对目标误差
  分别为 `+1.462732351 THz` 和 `-2.662128002 THz`。p/d 数值劈裂分别为
  `0.014910231 THz`、`0.013609269 THz`，满足简并容差。
- 自适应及局部扫描共完成 38 个候选。按 `0.1 nm/0.001` 精度得到的最近候选为
  `b0=244.0 nm, eta=0.940, zeta=1`：p/d 双重态平均频率误差分别为
  `+0.016861776 THz` 和 `-0.020155586 THz`。用户接受该最近候选。
- 已将接受的 cladding 参数写入 `scripts/parameter.json`，cavity 与 mesh9 保持
  不变；其他共享参数未改动。
- 完整 unit-cell band 已启动到
  `unitcell_band_cav(245-0.96-1.155)_clad(244-0.94-1)_mesh9`。因 cavity 参数和
  mesh 与既有完整结果相同，复用其 76 个 k 点（1021 个文件、76557711 bytes）；
  cladding 的 76 个 k 点全部重新计算。
- 完整能带结果和最终状态待计算完成后回填。
# 2026-07-26 目标修正：d doublet 对齐 cavity d1@Gamma

## 修正说明（实施前）

用户修正 cladding C6 优化目标：固定 cavity 与原有 cladding 约束不变，cladding
`p1/p2@Gamma` 的均值继续对齐 cavity `p1@Gamma`；cladding `d1/d2@Gamma` 的
均值由原先对齐 cavity `d2@Gamma`，改为对齐 cavity `d1@Gamma`。

当前 mesh9 cavity 正式结果给出的目标为：

- cavity `p1@Gamma = 196.64514917832147 THz`；
- cavity `d1@Gamma = 200.5023632843541 THz`。

实施时保留 `zeta=1`、`eta<1`、`b0=235--255 nm`、`eta=0.900--0.999`、
`b0` 精度 0.1 nm、`eta` 精度 0.001，以及原有对齐和简并容差。新目标使用独立
输出目录，避免覆盖旧的 d2 目标结论；相同 `(b0, eta)` 的 Gamma 原始求解结果
与目标无关，可从旧优化缓存只读复用。修改后先运行聚焦测试和静态检查，再执行
自适应搜索；实际参数与结果在本节末尾补记。

## 修正执行结果

- 优化脚本现读取 cavity `p1@Gamma` 与 `d1@Gamma`，结果字段、缓存有效性检查、
  preflight 和输出元数据均使用新的 d1 目标语义。
- 新任务输出使用独立目录
  `cladding_doublet_alignment_to_cavity_p1_d1_cav(245-0.96-1.155)_zeta1_mesh9/`；
  旧 d2 目标结果未被覆盖。
- 相同几何可只读复用旧任务的 Gamma 原始结果；新目标摘要不会误用旧目标摘要。
- 聚焦测试 85 项通过，脚本 `py_compile` 与 `git diff --check` 通过。
- 实际扫描 33 个候选，最终最优 cladding 参数为
  `(b0, eta, zeta) = (240.2 nm, 0.967, 1.0)`。
- cavity `p1@Gamma = 196.64514917832147 THz`，最优 cladding p doublet 均值为
  `196.64631420591434 THz`，有符号误差 `+0.00116502759286 THz`。
- cavity `d1@Gamma = 200.5023632843541 THz`，最优 cladding d doublet 均值为
  `200.51231105483797 THz`，有符号误差 `+0.00994777048388 THz`。
- p/d 数值分裂分别为 `0.0166972062270 THz` 和 `0.0120563329692 THz`；
  两项对齐误差和两项分裂均小于 `0.02 THz` 门限，候选状态为 `accepted`。
- 扫描结束后停止了尚未确认的自动完整能带验证，并删除其不完整目录及候选缓存；
  仅保留 `scan_config.json`、`scan_summary.csv`、`best_parameters.json` 和
  `80_logs/`。
- 完整能带验证改为默认关闭，只有显式将 `RUN_FULL_BAND_VALIDATION=True` 时才执行；
  避免参数确认前自动启动昂贵计算。
