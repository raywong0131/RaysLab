# 理论能带与全波仿真参数拟合计划

## 当前状态

- 状态：理论模型、频率单位、波矢换算、目标 `C2v` cell 和 `C6v` eta 扫描范围
  已批准；eta 步长及近 Gamma 采样规模待继续讨论，尚未实施。
- 本计划只设计 theory-simulation fitting workflow；本轮不启动 COMSOL，不修改
  `scripts/parameter.json`，不创建或覆盖 `scripts/.out/` 正式结果。
- 第一版只覆盖当前 unit-cell workflow 已稳定识别的 four-band
  `p/d` 子空间；six-band 在 four-band 验证通过后再扩展。

### 已批准合同

- 数值 Hamiltonian 使用以 Gamma 四模平均频率为零点的普通频率 THz 表示；
  `t0/t1/mu` 的数值单位为 THz，`v` 的数值单位为 THz um。同一拟合内不混用
  angular-frequency 参数。
- COMSOL 保存的坐标是 `k/G`，其中 `G=4*pi/(sqrt(3)*a)`。代入 S17 前统一转换为
  `k[1/um]=(k/G)*4*pi/(sqrt(3)*a)`。
- `t0/t1` 不是先验常数。首版先对固定 `b0=245 nm`、`zeta=1` 的 `C6v` 结构扫描
  `eta in [0.94, 1.04]`，标定 `t0(eta)/t1(eta)/v(eta)`；再在目标 `eta=0.96`
  处冻结带不确定度的标定结果，拟合当前 `zeta=1.156` cavity 的
  `alpha/beta/gamma/mu`。多 zeta `C2v` 联合拟合在该单目标闭环验收后实施。

## 目标

建立一条可复现的两阶段工作流：

1. 用 COMSOL unit-cell 计算生成带模式标签和场重叠证据的目标能带数据；
2. 在不启动 COMSOL 的离线程序中，先用同一几何族的 `C6v` eta scan 确定
   `t0(eta)`、`t1(eta)`、`v(eta)`，再使用目标 eta 处的标定值拟合 `C2v` 结构的
   `alpha`、`beta`、`gamma`、`mu`；
3. 输出残差、参数不确定度、相关系数、Jacobian 奇异值、可辨识性和留出波矢验证，
   而不是只输出一组最优参数。

不把 COMSOL eigensolve 放进优化器的目标函数。优化迭代只调用纯 NumPy/SciPy
理论前向模型；仿真数据通过稳定的数据集合同一次生成、反复复用。

## 当前 workflow 核对结果

### 已可复用

- `scripts/run_main/run_band_pair.py` 已能为每个固定 unit-cell 几何复用一个 COMSOL
  Model，逐 `k` 点只重求解；`run_band_solo.py` 复用同一单 case 实现。
- COMSOL 中实际使用
  `G = 4*pi/(sqrt(3)*a)`，保存的 `kx`、`ky` 是 `G` 的无量纲倍数，当前
  `a = 0.82 um`。拟合程序必须显式保存这一换算，不能把 CSV 中的无量纲 `k`
  直接当作物理波数。
- 每个有效本征模已保存复数中心面
  `eigenmodes/<mode_idx>_Hz_center.parquet`；`band_connector.py` 已提供复场重叠和
  Hungarian 唯一匹配。
- `eigenfrequencies.csv` 和 `selected_bands.csv` 已包含频率、Q、有效性、
  `p/d` 权重、`px/py/dx/dy` 成分、模式编号、场重叠和匹配状态。
- 当前 cavity 样例
  `unitcell_band_cav(245-0.96-1.156)_clad(240.1-0.967-1)_mesh5` 中四条带在
  `Gamma-M`、`Gamma-K` 两条正向支路均完整连接；可作为读取、预处理和拟合程序的
  只读开发样例，无需重跑 COMSOL。
- 项目已经依赖 NumPy、Pandas 和 SciPy，具备 `least_squares`、SVD 和 Hungarian
  matching 所需依赖。

### 现有缺口

- `docs/SI_thoery.md` 的 S16/S17 已给出需要拟合的同一个 first-order `C2v`
  four-band Hamiltonian。S16 使用旋转基底，S17 使用 standing-wave 基底；两者由
  S10 的固定酉变换连接，不是两套独立拟合模型。
- 当前 `p1/p2/d1/d2` 是 Gamma 点各 `p/d` 子空间内按频率排序的运行标签，不等同于
  理论基底标签，不能直接作为拟合列名。根据当前 `basis_utils.py` 的
  `cos(theta)/sin(theta)/cos(2theta)/sin(2theta)` 定义和 SI 的 S7，映射应固定为：
  `internal px -> p_y`、`internal py -> p_x`、
  `internal dx -> d_(x2-y2)`、`internal dy -> d_xy`；其中可能出现的整体负号只是
  基底相位，不影响本征频率。该映射仍需用解析向量测试和 Gamma 场对称性锁定。
- 当前标准采样只有正向 `Gamma-K`（`+kx`）和 `Gamma-M`（`+ky`）支路，没有负向轴和
  对角留出点。对 S17 的频率拟合，负向点由时间反演给出重复信息，可作为符号/实现
  QA 而不必进入首轮正式拟合；对角点仍是检验完整二维耦合结构的必要留出集。
- 现有结果中没有固定 `b0=245 nm`、`zeta=1` 且覆盖目标 `eta=0.96` 与
  Dirac reference `eta=1` 的 `C6v` 标定扫描；已有 `zeta=1` cladding 的 `b0/eta`
  不同，不能直接并入本次 `t0(eta)/t1(eta)` 标定。

## 关键设计决定

### 1. 分离采集、数据整理与拟合

建议新增三个层次：

- `scripts/run_sweep/collect_theory_fit_bands.py`：唯一允许启动 COMSOL 的采集入口；
- `scripts/analysis/fit_theory_band_parameters.py`：只读取既有数据并执行离线拟合；
- `comsol_workflow/` 下的纯数值模块：承载理论前向模型、数据合同、标定、残差、
  优化和不确定度计算，供脚本与测试共同复用。

`collect` 不调用 `run_band_solo.main()`，从而不会在采集后隐式回写共享中心频率。
拟合和重绘也不得修改 `scripts/parameter.json`。

### 2. 第一版固定为 four-band

当前正式后处理稳定选择两条 `p` 带和两条 `d` 带，因此第一版数据形状固定为
`(n_k, 4)`。six-band 需要先定义额外两条带的物理基底、Gamma 标签和追踪规则，
不与第一版混做，避免把模型扩展问题混入参数拟合验收。

### 3. 配置不写入共享参数文件

优化专属配置使用独立标准 JSON，建议为
`scripts/theory_band_fit.json`。它只保存：

- 参考与目标结构选择；
- 近 Gamma 采样、拟合集/验证集划分；
- 参数边界、权重、robust loss、多起点随机种子；
- 理论模型版本和基底映射版本；
- 输出系列标签。

cell、mesh、unit-cell eigensolver 等共享值仍从 `scripts/parameter.json` 读取。
运行时把两份配置解析后的完整快照写入输出 manifest。

### 4. 输出位置

所有新结果继续位于允许的一级目录
`scripts/.out/unit_cell_band/`，建议使用两个非冲突系列：

```text
scripts/.out/unit_cell_band/
  theoryfit_dataset_<cell-and-mesh-label>/
    calibration_zeta1/eta<value>/
    target_zeta<value>/
    99_config/
  theoryfit_result_<cell-and-model-label>/
    10_calibration/
    20_fit/zeta<value>/
    30_validation/
    99_config/
```

数据集 manifest 记录源目录、参数快照、solver identity、采样点、文件校验信息和
完成状态。拟合结果引用原始数据，不复制全部场 parquet。已有正式结果只读复用，
不迁移、不覆盖。

## 数据合同

### 仿真数据集

整理后的长表每行表示一个 `(structure_id, k_point, theory_mode)`，至少包含：

- `structure_id`、`b0_nm`、`eta`、`zeta`、`mesh_size`；
- `kx_over_G`、`ky_over_G`、`k_abs_over_G`、`split`；
- `simulation_mode_idx`、`runtime_band_label`、`theory_mode_label`；
- `frequency_thz`、`frequency_imag_thz`、`q`、`is_valid`；
- `p_weight`、`d_weight`、`px_weight`、`py_weight`、`dx_weight`、`dy_weight`；
- `match_status`、`field_overlap`、原始结果相对路径。

`theory_mode_label` 只能由已批准的基底映射生成。任何 Gamma 模式纯度不足、支路
unmatched 或场重叠低于阈值的点保留在审计表中，但默认不进入拟合残差。

### 理论前向模型

新增 `comsol_workflow/c2v_band_model.py`。以 SI 的 S17 为唯一 canonical 数值实现，
standing-wave 基底顺序严格固定为

```text
(d_xy, d_(x2-y2), p_y, p_x)
```

并使用 S13/S14 的组合参数：

```text
mx  = gamma * (t0 - t1)
my  = ((4*beta - gamma)*t0 - (2*alpha + gamma)*t1) / 3
chi = (4*alpha - gamma) / 3
```

公开稳定接口：

```python
def calculate_bands(
    k_points_over_g,
    *,
    lattice_constant_um,
    t0_thz,
    t1_thz,
    v_thz_um,
    alpha,
    beta,
    gamma,
    mu_thz,
) -> BandSolution:
    ...
```

COMSOL 保存的是 `k/G`，其中 `G=4*pi/(sqrt(3)*a)`；进入 S17 前转换为
`k[1/um] = (k/G)*G`，使 `v*k` 与其他 Hamiltonian 元素同为 THz。

`BandSolution` 同时返回本征频率和本征矢，以便理论侧也从 Gamma 向外做确定性追踪，
而不是在每次优化迭代中仅按本征值排序。Hamiltonian 必须逐点验证 Hermitian；
参数为实数时优先使用 `numpy.linalg.eigh/eigvalsh`。

S16 不另建可独立漂移的生产实现。测试中按 S14 构造 S16，并严格验证
`H_sw = U @ H_spin @ U.conj().T` 与 S17 相等，以及二者本征值逐点一致。

数值拟合建议统一使用相对 Dirac/四模中心的普通频率 THz：即把 SI 中的 `omega`
整体除以 `2*pi` 后再定义 `t0/t1/v/mu`。这样 COMSOL 的 `frequency_thz` 可直接居中后
进入残差，且 `alpha/beta/gamma` 不受单位选择影响。manifest 必须记录该约定，禁止
在同一拟合中混用 rad/s 与 THz。

必须固定以下元数据：

- 理论来源：`SI_thoery.md` S10、S13、S14、S16、S17；
- standing-wave 基底顺序和内部 Fourier 标签映射版本；
- 输入 `k/G`、实际 `G`、转换后的 `k[1/um]` 和 `a[um]`；
- `t0/t1/mu` 使用 THz、`v` 使用 THz um；
- `zeta=1` 时严格回归 `alpha=beta=gamma=1, mu=0`。

## 执行阶段

### 阶段 A：锁定理论模型合同

1. 将 `SI_thoery.md` S17 及其 S13/S14 参数定义写入模块 docstring，S16/S10 作为
   等价表示和交叉验证来源。
2. 建立解析极限测试：
   - Gamma 点按 Appendix A 的 A2 严格返回
     `(mx-mu, my+mu, -my+mu, -mx-mu)`；
   - `C6v` 参数极限恢复预期双重简并；
   - `H(k) = H(k)^dagger`；
   - S16 经 S10 变换后逐元素等于 S17；
   - 满足时间反演下 `f(k)=f(-k)`；
   - 当前 Fourier basis 向量与 S7 标签映射逐向量相等到整体相位。
3. 在这些测试通过前不进入 COMSOL 采集。

### 阶段 B：抽取通用近 Gamma 采样器

1. 把当前仅支持 `gamma_m/gamma_k` 的采样构造扩展为由
   `(direction vector, magnitudes, split)` 定义的通用射线；现有 pair/solo 默认输出和
   文件名保持不变。
2. 第一轮建议采样：
   - 拟合集：Gamma、`+delta k` 与 `+2 delta k` 的 x/y 轴点；
   - 验证集：`(delta k, delta k)`、`(delta k, -delta k)`；
   - QA 集：少量负向轴点，只验证时间反演、COMSOL Bloch 符号和标签连续性，不作为
     独立统计样本重复加权；
   - 可选稳健性集：更密的近 Gamma 轴向点，用于改变 `k_max` 和采样密度检查。
3. `delta k` 先以 `k/G` 配置。最终数值应在查看理论有效尺度和现有 mesh 收敛后批准，
   不直接照搬当前 Q 拟合点。
4. 每条射线从 Gamma 向外独立使用现有复场重叠与 Hungarian matching；Gamma 重复点
   在汇总时按 solver identity 去重。

### 阶段 C：生成并冻结标定/目标数据

1. 标定结构族固定 `b0=245 nm`、`zeta=1`，在已批准的
   `eta in [0.94, 1.04]` 内扫描；网格必须精确包含目标 `eta=0.96` 和 Dirac
   reference `eta=1`，步长在正式运行前批准。
2. 目标结构第一版使用一个 `zeta != 1`；建议先以当前 cavity
   `(245 nm, 0.96, 1.156)` 为样例。
3. 昂贵运行前打印并保存 dry-run manifest，核对 geometry、mesh、中心频率、
   eigenmode count、全部 k 点、缓存命中路径和输出目录。
4. COMSOL 运行后至少检查 mesh 完成、首个完整点、Gamma 四模纯度和所有支路追踪
   状态，再冻结 dataset manifest。
5. 当前已有 `zeta=1.156` cavity 结果可先开发离线读取器；正式拟合是否复用该结果，
   由 manifest 的完整 solver identity 和新增采样覆盖共同决定。

### 阶段 D：C6v eta scan 与 t0/t1 标定

1. 对每个 `eta_j` 独立减去 Gamma 四模平均频率，并由 Gamma 双重态计算
   `m(eta_j) = (mean(f_d) - mean(f_p)) / 2 = t0(eta_j)-t1(eta_j)`。
2. 用两条轴向支路共同拟合 S12 的有质量色散
   `f_offset(k)=+/-sqrt(m**2+(v*k)**2)`：若 `m≈0`，它自然退化为线性色散；若
   `m!=0`，禁止用 Gamma 点零斜率代替完整拟合。两方向应共享同一个 `v(eta_j)`，
   方向差异作为 `C6v` 各向同性诊断。
3. 对每个扫描点计算
   `t1(eta_j)=2*v(eta_j)/a`、`t0(eta_j)=t1(eta_j)+m(eta_j)`，并保存
   `m/v/t0/t1` 的协方差与拟合窗口敏感性。这里的 `t0/t1` 是从全波数据确定的待定
   参数，不是输入常数。
4. 若 Gamma 处 `s/f` 模式能被可靠识别，则用 S8 做独立交叉检查：令
   `c_sf=(f_f-f_s)/2=2*t0+t1`，得到
   `t0=(c_sf+m)/3`、`t1=(c_sf-2*m)/3`。该 remote-band 估计只作验证，除非
   `s/f` 谱隔离和 six-site 截断误差通过检查。
5. 可选地用完整 C6v six-site S4 对近 Gamma 数据直接拟合 `t0/t1`，与
   `m+v` 提取结果比较；显著差异视为 nearest-neighbour/first-order 模型误差，
   不通过增加权重掩盖。
6. 先保留每个 eta 的独立估计，再建立带不确定度的平滑
   `t0(eta)/t1(eta)/v(eta)` 曲线；不预设线性或单调函数。`eta=1` 用于检查
   `t0=t1`，不在看到数据前硬性强制。
7. 输出 `eta_calibration.csv`、`calibration.json`、逐点残差、标定曲线、方向一致性、
   `s/f` 交叉检查和 `m≈0` 诊断。目标拟合读取 `eta=0.96` 的冻结标定分布，而不是
   无误差的单个常数。

### 阶段 E：目标数据预处理与解析初值

1. 用 Gamma 模式成分映射到理论四基底，并沿各射线保留现有场重叠追踪结果。
2. 对仿真和理论分别减去本结构 Gamma 四模平均频率，消除公共频率平移。
3. 按附带文本公式计算 `mx_target`、`my_target`、`mu0`、`gamma0`、`alpha0`、
   `beta0`。
4. 当 `abs(t0-t1)` 低于配置阈值时，不使用 `mx/m` 初始化 `gamma`，改用近 Gamma
   速度估计，并在结果中记录初始化分支。
5. 所有解析初值先投影到批准边界内；原始值与投影值都写入审计 JSON。

### 阶段 F：多起点非线性最小二乘

1. 使用 `scipy.optimize.least_squares`，默认 `loss="soft_l1"`，固定随机种子生成
   可复现的多起点。
2. 残差使用
   `sqrt(weight) * (predicted - observed) / frequency_scale`。若配置中的 `weight`
   表示平方损失权重，必须取平方根；不能直接乘 `weight` 后再次平方。
3. Gamma、`delta k`、`2 delta k` 和不同模式权重全部来自配置，并在输出中展开为逐点
   实际权重。
4. 基础硬边界只保证 `alpha,beta,gamma > 0`、`mu >= 0`。关于 `zeta` 的理论趋势先
   作为可开关的弱约束/软惩罚，避免用理论预期掩盖模型与全波结果的系统偏差。
5. 若不同起点收敛到显著不同的参数簇，再触发可选
   `differential_evolution -> least_squares` 后备流程；不默认支付全局搜索成本。

### 阶段 G：不确定度、可辨识性与验证

1. 对最优点残差 Jacobian 做 SVD，输出全部奇异值、数值秩、条件数和最弱可辨识
   参数组合。
2. 协方差使用 SVD 伪逆，不直接计算 `(J.T @ J)^-1`；报告自由度、残差方差和截断
   阈值。对 robust loss 明确标注这是局部线性近似。
3. 同时报告多起点解的离散程度；必要时增加基于 k 点重采样的 bootstrap，避免只依赖
   高斯线性近似。
4. 分开报告“给定 `t0/t1/v` 的条件不确定度”和“传播 eta 标定协方差后的总不确定度”。
   后者通过 calibration posterior/近似高斯样本重复拟合传播，避免把第一阶段待定参数
   当成精确常数。
5. 输出：
   - `fit_result.json`；
   - `parameters.csv`、`residuals.csv`；
   - `covariance.csv`、`correlation.csv`、`jacobian_singular_values.csv`；
   - 拟合/仿真能带对比、残差图和参数相关图；
   - 拟合集与对角留出集的独立 RMSE/最大误差。
6. 改变 `k_max`、采样密度、权重方案和模式纯度阈值，形成稳定性表；比较
   `gamma_mass` 与 `gamma_velocity`。

### 阶段 H：多 zeta 扩展

four-band 单目标验收后再实施：

1. 对每个 `zeta_j` 独立拟合并保留完整诊断；
2. 以独立解为初值进行联合平滑；
3. `alpha/beta` 增、`gamma` 减、`mu>=0` 和导数绝对值次序默认作为软惩罚；
4. 同时输出无趋势约束与有趋势约束结果，量化先验对结论的影响；
5. 用未参与联合拟合的 zeta 或 k 方向做外部验证。

## 验收标准

- 纯数值单元测试不导入或启动 `mph`。
- 已知合成参数能从无噪声 four-band 数据恢复到数值容差内；带噪声测试覆盖置信区间、
  参数相关和病态 Jacobian 告警。
- 模式 crossing/inversion 合成测试确认标签来自本征矢/子空间连续性，而非逐点频率排序。
- `C6v` 合成极限返回 `alpha=beta=gamma=1, mu=0`，并正确处理双重简并。
- 合成 eta scan 能分别恢复预设的 `t0(eta)`、`t1(eta)` 和 `v(eta)`；`eta=1` 检查
  `t0=t1`，`m+v` 提取与可选 S4/S8 交叉估计在设定容差内一致。
- 数据读取测试锁定当前 `selected_bands.csv` 列、Gamma 去重、标签映射、无效点过滤和
  train/validation split。
- 采集入口测试锁定任意射线、负向/对角点、缓存 identity、非冲突目录以及绝不回写
  `scripts/parameter.json`。
- 运行聚焦 pytest、相关 `py_compile` 和 `git diff --check`；绘图实施时按 MEMORY 的
  科研绘图规则检查实际 PNG。

## 风险与回滚

- 最大风险不是优化器，而是理论基底与仿真模式的错误映射；映射未批准时禁止给出
  物理参数结论。
- `m≈0` 时 `gamma` 的质量估计不可辨识；必须依赖速度信息并显式报告条件数。
- 只使用 Gamma `p/d` 频率只能得到 `t0-t1`，无法分别确定 `t0/t1`；标定必须包含
  近 Gamma 速度，且 `m+v`、S4/S8 估计不一致时应报告模型误差。
- 理论模型若只在极近 Gamma 有效，加入远端现有 31 点全 BZ 数据会产生有偏参数；
  默认只读取批准的 `k_max`，全 BZ 数据仅作模型失效诊断。
- robust loss 能隐藏少量坏点，但不能替代模式匹配审计；每个被降权或剔除的点必须
  可追溯。
- 回滚只删除本任务新增源码、测试、配置和新建的 theory-fit 输出；不触碰既有
  unit-cell、strip 或 finite 正式结果。

## 需要用户确认的决策

1. 确认已批准 `eta in [0.94, 1.04]` 范围内的扫描步长，以及是否采用粗扫后局部
   加密；网格必须包含 `eta=0.96` 与 `eta=1`。
2. 确认第一轮近 Gamma 的 `delta k`、`k_max` 和允许的 COMSOL 新计算规模；负向点仅
   作 QA，对角点作正式留出验证。

## 执行记录

- 2026-08-19：读取拟合设想和当前 pair/solo workflow；确认已有复场、模式成分、
  Hungarian tracking 与可复用 cavity 样例，确认同几何族 `zeta=1` 标定数据缺失。
  仅新增本计划，未实施代码，未启动 COMSOL，未改动正式结果。
- 2026-08-19：补读 `docs/SI_thoery.md` S2 和 Appendix A；将 S17 固定为拟合主实现、
  S16/S10 固定为等价性测试，并解析出当前 internal Fourier mode 到 SI standing-wave
  基底的精确映射。仍未实施代码或启动 COMSOL。
- 2026-08-19：用户批准以 Gamma 四模中心化后的 THz 作为数值频率单位，并批准将
  COMSOL 的 `k/G` 按 `G=4*pi/(sqrt(3)*a)` 转换为 `k[1/um]` 后代入 S17。
- 2026-08-19：用户初步批准单个 `eta=0.96`、`zeta=1` `C6v` 参考结构，加当前
  `zeta=1.156` cavity 的单目标 four-band 拟合闭环；多 zeta 后置。
- 2026-08-19：用户进一步明确 `t0/t1` 也是待定参数，必须通过 `zeta=1` 的 eta scan
  确认。计划据此把标定阶段扩展为 `t0(eta)/t1(eta)/v(eta)` 独立提取、交叉验证与
  不确定度传播；目标 `C2v` cell 保持不变。
- 2026-08-19：用户批准 `C6v` 标定扫描范围为 `eta in [0.94, 1.04]`；步长与是否局部
  加密仍待确认。
- 2026-08-19：用户批准开始实施，并采用粗扫步长 `0.02`，即
  `0.94, 0.96, 0.98, 1.00, 1.02, 1.04`；只有数据出现显著曲率时才追加 `0.01`
  局部加密。
- 2026-08-19：完成 `comsol_workflow/c2v_band_model.py`。S17 是 canonical 数值实现，
  S16 仅用于 S10 酉变换等价性验证；所有接口固定使用中心化普通频率 THz、
  `v[THz um]` 与显式 `k/G -> 1/um` 换算。
- 2026-08-19：完成 `comsol_workflow/theory_band_fitting.py`，包含 C6v 单 eta 的
  signed mass/velocity/t0/t1 标定、standing-wave 基底锚定的 Hungarian 本征矢追踪、
  deterministic multistart least-squares、SVD covariance/correlation/秩/条件数诊断，
  以及 `t0/t1/v` 标定协方差的 Monte Carlo 重拟合传播。
- 2026-08-19：新增独立 `scripts/theory_band_fit.json`、默认 dry-run 的
  `scripts/run_sweep/collect_theory_fit_bands.py` 和纯离线
  `scripts/analysis/fit_theory_band_parameters.py`。`run_band_pair.run_case` 只增加可选
  `branches` 参数；原 pair/solo 默认采样和输出不变。未修改 `scripts/parameter.json`。
- 2026-08-19：生成正式 dry-run manifest：6 个 C6v 标定结构加 1 个 C2v 目标结构，
  预计 95 个去重 eigensolve、7 次固定几何 Model/mesh 构建；快照为 mesh 5、10 模、
  `c_const/1.55[um]`、中心频率 `198.38683823113445 THz`。manifest 位于
  `scripts/.out/unit_cell_band/theoryfit_dataset_b245_eta0.94-1.04_target_eta0.96_zeta1.156_mesh5/99_config/dry_run_manifest.json`。
  本轮没有导入或启动 COMSOL，也没有写入任何 eigensolve 正式结果。
- 2026-08-19：聚焦验证通过：Hamiltonian、合成标定/拟合、病态 Jacobian、模式连续性、
  manifest/CSV 合同以及受影响 pair/solo 共 `117 passed`；相关文件 `py_compile`
  通过，`git diff --check` 通过（仅报告仓库既有文件的 Windows 行尾提示）。Ruff 未运行，
  因当前 uv 环境未安装该程序。

## 当前实施状态与剩余步骤

- 已完成：阶段 A、B 的通用采样接口、D、E、F、G 的纯数值实现，以及阶段 C 的
  dry-run manifest。
- 待启动：用户复核 95 次 eigensolve 的 manifest 后，显式执行 collector；完成 mesh、
  首个完整点、Gamma 四模纯度和 branch tracking 检查后冻结 dataset manifest。
- 待数据：COMSOL 数据完成后运行离线拟合，检查 eta 曲率决定是否在局部追加 `0.01`
  扫描，并输出留出对角点、负向 QA、残差、协方差和标定传播区间。
- 回滚：本轮新增源码、测试、配置和 theoryfit dry-run 目录可独立移除；既有 unit-cell、
  strip、finite 正式结果均未移动、覆盖或删除。

## Appendix C 二阶 direct-projection 扩展（2026-08-20）

- 理论合同更新为 `docs/SI_thoery_v2.md` Appendix C，新增 C2--C5 的 direct
  quadratic correction；S17 仍保留为一阶基线，二阶 canonical 实现使用 standing-wave
  基底 C5，并通过 C4 加 S10 变换做逐元素验证。
- 二阶项不新增拟合自由度：`Lambda_x/Lambda_y/Lambda_xy` 只读取已冻结的
  `a,t1,alpha,gamma`，单位继续为普通频率 THz；`beta,mu` 只通过原 Gamma/一阶项进入。
- 数值测试必须锁定 C2 系数、Gamma 退化为 S17、Hermitian、C4/C5 等价、时间反演谱，
  以及 `Lambda_xy` 在主轴消失而在 off-axis 非零。
- 使用现有完整 `(245,0.96,1.156)` run-band cavity CSV，保持 M--Gamma--K 路径和
  45 THz 频率窗，新建 `full_band_second_order_comparison.png`；同时绘制 COMSOL、
  first-order S17、second-order C5，并保留 `|k|/G<=0.02` 拟合窗标记。
- 分开报告一阶和二阶在 Gamma-M、Gamma-K 的 full-path RMSE。Appendix C 的 C6 明确
  排除 remote `s/f` Lowdin self-energy，因此 direct k^2 仍失败时不得继续把误差归因于
  原四个调制参数，也不得在没有新理论合同的情况下添加经验多项式。
- 2026-08-19 实际启动：首次显式执行在导入 COMSOL 前暴露 `_unique_case_points`
  仍硬编码 `gamma_m/gamma_k` 的兼容遗漏；未产生 eigensolve。已改为遍历任意 branch，
  新增回归测试并保留 `collector_*.log.failed_branch_key` 审计日志。
- 修复后后台 collector 与 `comsolmphserver` 已启动。eta=0.94 的 Gamma、
  `kx/G=0.002` 和 `kx/G=0.004` 已连续完成，stderr 为空；Gamma 目标双重态的 p/d
  纯度约为 `0.9992--0.9993`。manifest 当前为 `running`，任务继续在本机后台执行。
- 2026-08-19：COMSOL 数据集完成，manifest 为 `complete`，95/95 个去重 k 点、
  1197 个有效模式场全部落盘；7 个结构均为 0 unmatched、0 invalid，collector stderr
  为空。最小 matched overlap 为 eta=1 的 `0.5706`，目标 C2v 为 `0.9314`。
- 2026-08-19：完成离线标定和 S17 拟合。eta=0.96 标定得到
  `m=2.388981 THz`、`v=10.173604 THz um`、`t0=27.202650 THz`、
  `t1=24.813669 THz`；eta=1 得到 `m=-0.005436 +/- 0.008224 THz`，因此
  `t0=t1` 在误差内成立。目标拟合得到 `alpha=0.923957`、`beta=0.894576`、
  `gamma=0.848670`、`mu=1.137586 THz`；Jacobian rank 为 4，条件数 `112.58`，
  未判为病态。
- 两个对角非 Gamma 留出点 RMSE 分别为 `0.0471`、`0.0454 THz`，最大误差小于
  `0.062 THz`；两个负向 QA 点 RMSE 分别为 `0.0187`、`0.0347 THz`。结果目录补齐
  `parameters.csv`、`residuals.csv`、`validation_metrics.csv` 和
  `jacobian_singular_values.csv`。其中 conditional covariance 与 calibration-only
  propagation 分开报告，后者不得误读为包含拟合残差后的总置信区间。
- 拟合后验证：纯数值、数据合同和 workflow 聚焦测试 `11 passed`，相关文件
  `py_compile` 与 `git diff --check` 通过；本阶段未启动 COMSOL，未修改
  `scripts/parameter.json`。
- 2026-08-19：基于 `(245, 0.96, 1.156)` 已冻结拟合结果生成
  `band_fit_comparison.png`。面板顺序固定为 `Gamma-K (+kx)`、`Gamma-M (+ky)`，
  横轴为 `k/G in [0, 0.02]`，纵轴为恢复 Gamma 四模中心后的绝对频率 THz；COMSOL
  用半透明实线，S17 用虚线，模式颜色固定为 `p_y` 红、`p_x` 蓝、`d_xy` 浅灰、
  `d_x2-y2` 黑。图中非 Gamma RMSE 分别为 `0.057`、`0.081 THz`。
  实际 2181x992 PNG 已完成视觉检查并修正横轴刻度重叠；绘图相关聚焦测试
  `12 passed`，`py_compile` 与 `git diff --check` 通过，未启动 COMSOL，也未覆盖
  其他正式结果。
- 2026-08-19：进一步复用较新的完整 run-band cavity 数据
  `unitcell_band_cav(245-0.96-1.156)_clad(240.1-0.967-1)_mesh5`，生成
  `full_band_fit_comparison.png`。图使用与 `run_band` 相同的 `M-Gamma-K`
  归一化路径和以 `p2@Gamma` 为中心的 45 THz 频率窗，并用浅灰区域明确
  `|k|/G <= 0.02` 的实际拟合窗。完整路径没有 unmatched/invalid 数据。
  一阶 S17 强制外推后的全路径 RMSE 为 Gamma-M `10.94 THz`、Gamma-K
  `5.80 THz`，相对于近 Gamma 的 `0.06--0.08 THz` 显著增大；该差异记录为
  first-order/k.p 截断模型的适用域诊断，不通过重新调权掩盖。实际 1927x1337 PNG
  已视觉检查；聚焦测试 `12 passed`，`py_compile` 和 `git diff --check` 通过，未重跑
  COMSOL，未覆盖近 Gamma 比较图或其他正式结果。
- 2026-08-20：完成 Appendix C direct-projection 二阶扩展。数值实现新增 C2 的
  `Lambda_x/Lambda_y/Lambda_xy`、standing-wave 基底 C5，以及 C4 经 S10 变换到 C5
  的逐元素交叉验证；二阶项只使用冻结的 `a/t1/alpha/gamma`，没有新增拟合自由参数，
  一阶 S17 默认行为保持不变。Gamma 恢复、Hermitian、时间反演、主轴/off-axis 和
  C4/C5 等价测试均通过。
- 2026-08-20：使用同一完整 run-band cavity CSV 生成独立的
  `full_band_second_order_comparison.png`，未覆盖既有一阶图。图中 COMSOL 为半透明
  实线，一阶 S17 为点线，二阶 C5 为虚线，路径、模式颜色、45 THz 频率窗和
  `|k|/G <= 0.02` 拟合窗均与既有比较合同一致。全路径一阶/二阶 RMSE 分别为：
  Gamma-M `10.939/17.121 THz`，Gamma-K `5.798/31.229 THz`；最大绝对误差分别为
  Gamma-M `30.288/61.605 THz`，Gamma-K `16.666/127.022 THz`。因此 direct `k^2`
  项没有改善 full-BZ 外推，尤其在 Gamma-K 上显著恶化；该结果与 Appendix C6 未包含
  remote `s/f` 子空间 Löwdin self-energy 的限制一致，不能据此重新调整已冻结的四个
  调制参数或添加无理论合同的经验多项式。
- 新图实际为 `2110x1337` RGBA PNG，已检查图例、线型、频率窗、拟合窗标注和裁切；
  视窗外的二阶曲线保持裁切以诚实复用 COMSOL 的 45 THz 比较尺度。相关聚焦测试
  `15 passed`，`py_compile` 与 `git diff --check` 通过。本阶段未启动 COMSOL，未修改
  `scripts/parameter.json`，也未覆盖或删除其他正式结果。
