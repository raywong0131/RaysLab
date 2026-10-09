# S5：实际有限基模的全局 Γ 分量与中心远场——执行 Spec

- 日期：2026-10-06。
- run_id：`DS-20261006T073408Z-s5-gamma-origin-6c7bd2`。
- 正式位置：`I:/codeXproject/comsol-workflow/docs/spec/plan-execute_20261006_s5_gamma_origin.md`。
- 执行项目：`I:/codeXproject/comsol-workflow`，SSH alias `qd228-comsol`。
- 状态：G0/G1 与外域 mesh5 已完成；FULL 已启动，GAMMA/REST 和综合报告已接入自动顺序流程，尚未验收最终辐射结果。
- 科学方案源：本请求的 `source_plan.md`，对应 Vault 中的英文 `S5 Gamma Component Radiation Simulation Plan.md`。

## 1. 用户授权与问题

本次直接授权原文：

> 按照远端Project的要求，在指定位置构建Spec文档，推送给SSH端，开始进行执行计算

该指令授权本 Spec 范围内的脚本实现、必要聚焦测试、COMSOL 启动、保存模型载入、补充场导出、外域建模、三组频域求解和后处理。
英文源方案原先的“尚未授权执行”是编写计划时的历史状态，已由上述指令对本次范围取代。
不自动 commit/push，不修改 SI、共享 `scripts/parameter.json` 或既有结果。

本次只验证第一项机制：实际有限基模的全局 index-0 Γ 分量，是否解释法向中心的复振幅和中心亮斑。
三组输入为 FULL（完整实际场）、GAMMA（全局 index-0）和 REST（完整余项）。
每组均保留实际 cladding；本次不开展 cladding A/B/C 结构替换。
已确认的零点位置一致性作为前提，不重新要求网格、频率、源窗口、几何或尺寸扫描。
结果不得预设 GAMMA 必须成功解释中心；REST 的有限贡献必须如实报告。

## 2. 输入身份与固定参数

结果根目录以下简称 `R`：
`I:/codeXproject/comsol-workflow/results/S5_dipolarSingularity_analysis`。

| 输入 | 作用 |
| --- | --- |
| `R/00_model/02_gamma_cladding/A_solved.mph` | 实际结构的已保存实频驱动响应 |
| `R/01_results/02_gamma_cladding/A_cavity_EH.npz` | 1141 个完整 cavity cell 的共用 cell 内场采样 |
| `R/01_results/02_gamma_cladding/A_slab_EH.npz` | 原薄板截面场 |
| `R/01_results/02_gamma_cladding/A_air_EH.npz` | 原主空气平面的 E/H |
| `R/01_results/02_gamma_cladding/A_air_check_EH.npz` | 已有独立空气平面 |
| `R/99_config/02_gamma_cladding/A_simulation_config.json` | 已验证 footprint、holes 和完整几何来源 |
| `R/99_config/02_gamma_cladding/geometry_and_run_config.json` | 实频响应、源和资源配置来源 |
| `R/80_logs/02_gamma_cladding/source_and_projection.csv` | A 与目标模式对应关系 |

固定：user `py` / internal `px`，原 mode19；`f0=197953000904452.1 Hz`；mesh5；`a=820 nm`；薄板厚度 `200 nm`；折射率 `3.3`。
cavity：`b=245 nm, eta=0.96, zeta=1.156`，1141 cells，shell 0–19。
cladding：`b=242 nm, eta=0.98, zeta=0.93`，3540 cells，20 层，保留实际孔洞、位移和 profile。
保存对称性：quarter ID1，xPEC/yPMC，z=0 PMC。
以实际保存配置为准核对 pad、材料、空气域、PML、source support；不得从当前共享参数补齐缺失项。
实际几何经 `load_simulation_config_geometry()` 加载，不猜测 cell 拓扑。

A 是与目标本征模高度相符的实频响应，不是原复本征频率场本身。
本实验采用同一实频 A 来保证外域重放可比较；不得将复频本征场混入实频算子并宣称严格相等。
载入时查询实际 study/solver/dataset；预检为 `std1/freq, sol1/s1, dset1`，驱动解使用 `solnum=1`，不是 mode19。
旧元数据中的 `results/XX_dipolarSingularity_analysis` 只保留为历史来源；实际文件在当前 S5 根目录解析并哈希验证，不重建旧 alias。

## 3. 分解定义：全局 index-0，而非带包络的局部 Γ 投影

令 `N=1141`，`R_j` 为完整 cavity cell 的平移，`r,z` 为所有 cell 共用的局部坐标，`E_A,H_A` 为保存的 A 复场。
定义每个局部坐标上的复平均：

$$
\overline{\mathbf E}_\Gamma(\mathbf r,z)
=\frac{1}{N}\sum_{j=1}^{N}\mathbf E_A(\mathbf R_j+\mathbf r,z),
\qquad
\overline{\mathbf H}_\Gamma(\mathbf r,z)
=\frac{1}{N}\sum_{j=1}^{N}\mathbf H_A(\mathbf R_j+\mathbf r,z).
$$

把相同的平均 profile 重复到每个 cavity cell；完整余项为实际场减去该 profile。
不做逐 cell 相位对齐，不保留逐 cell 的包络乘子，不独立归一化两组。
复平均包含所有 cell 内的场分量，不限于单一周期性 py 模板。
不得将历史 19 个 k 点 ×4 个周期态的 76 个 basis 当作 1141 个 cell 的完整分解。

复用 `comsol_workflow/finite_lattice_fourier.py` 中 `hex_cyclic_indices`、`hex_xi_values`、`finite_lattice_fourier_components`。
验证正交归一化的 index-0 正变换为 `sqrt(N)` 乘复平均，逆变换贡献为平均 profile。
REST 包含所有非零 index 的完整余项，不只是选出的邻近 k 点。
检查 FULL=GAMMA+REST，且 REST 在共用局部坐标的 cell 平均为零。
该定义是有限 cell 的分解，并不把截断 Γ 场认作无限周期 Maxwell 本征态。

## 4. 闭合输入面与场导出

以实际完整 cavity cell union 的 footprint 构建棱柱，其闭合边界记为 `S`。
上下 cap 初值 `z=±300 nm`，必须位于实际空气域内、PML 外，并包住全部原外加电流。
侧面沿实际 cell union 边界；不得移除任何 cladding 孔或 cladding 物理域。
如 300 nm 不符合真实几何，选择一个有效值并记录原因；不开展窗口扫描。

补导出 S 全部侧面和两 cap 上的复 E/H、坐标、向外法线、面标签、积分权重、所属 cell 与材料侧信息。
同时在所有 cavity cell 的共用 cap/side 局部采样坐标求值，用于获得同一 Γ profile。
已有 `A_cavity_EH.npz` 的 3078 个体积分点不是闭合面采样；已有 `finite_surface_EH*.npz` 是单平面，不可替代 S。

初始采样：cap 用 61×61 共用局部 XY 网格，裁剪到 unit cell 且包含边界节点；每个 unit-cell edge 61 点；侧面 z 间距不大于 5 nm，并含薄板界面和 cap。
在材料界面保留合适的一侧值和分面插值，不跨越法向跳变做平滑。
批量插值、利用已验证镜面对称性避免重复取点，再恢复完整 1141 cells 后求平均。
复用 `NativeSampler`，六个 E/H 分量分别使用极向量和轴向量反射规则。
zPMC：E 的平面分量偶、Ez 奇；H 的平面分量奇、Hz 偶。
`axis='x'` helper 是关于 x 轴镜像，实际翻转 y；应使用显式坐标反射矩阵验证映射。

先固定共用的、边角及 unit-cell seam 兼容的切向场表示，再定义边界输入。
单独采用 cell index 的 tie-break 不足以保证 seam 上的切向场兼容。
所有相邻面共用一致的 edge/corner 数据；记录 projector/trace mapping，不为提高 Γ 贡献而修改定义。

令 `n` 为 S 的单位法线，`E_t=E-n(n·E)` 为切向电场。
三组边界输入定义为 `g_FULL=E_A,t|S`，`g_GAMMA=E_GAMMA,t|S`，`g_REST=g_FULL-g_GAMMA`。
REST 在同一边界表示上相减，不能独立插值破坏恒等式。

## 5. 同一外域算子下的三组求解

复制 A 到本次新模型，几何分区后只将 S 内部从 EWFD 求解域排除。
保留外部实际 cladding、pad、材料、空气、PML 和可用对称性；禁用原 impressed current。
新内边界采用 COMSOL `Electric Field` 条件，只规定切向 E，不同时规定 E 和 H。
H 由外域求解得到，用于独立重放检查。
各平面插值使用两个局部坐标，不使用退化的三维共面点云；侧面在实际材料界面处分 chart。
每组可用三个实部、三个虚部函数构造复场；COMSOL 中 REST 直接写 FULL 减 GAMMA。

这里是规定完整切向场的外域问题，不是平面波照射或独立 Γ 入射效率。
FULL trace 已含原结构反馈，因此删去内域仍可重放外域。
GAMMA/REST 是该固定边界分解下的贡献，不是另行本征求解后的模式。
三组仅换边界输入，几何、材料、mesh、频率、观察面、相位、振幅和远场算子完全相同。
不做等源功率或等储能重标定，不为每组另找共振。

官方接口依据：
- [Electric Field 条件](https://doc.comsol.com/6.3/doc/com.comsol.help.woptics/woptics_ug_optics.6.22.html)。
- [Far-field 理论](https://doc.comsol.com/6.3/doc/com.comsol.help.rf/rf_ug_theory.06.12.html)。

## 6. 阶段、程序职责与停止条件

| 阶段 | 操作及通过条件 |
| --- | --- |
| G0 audit | 不导入 mph、不启动 COMSOL；核对参数/源路径/哈希、1141 cells、输出非冲突、磁电分量、原源支撑域。缺失闭合场必须记 missing。 |
| G1 export/decompose | 启动 COMSOL 后只读取 A 已保存解，补导出闭合面，计算完整 Γ/REST；输入恒等式、REST 复平均和 DFT index-0 检查通过。 |
| G2 FULL pilot | 建一个外域 mesh，求 FULL；与 A 的同坐标空气 E/H、中心复振幅及中心图样比较。不独立调相、缩放；若不能重放，先修实现，不解释 Γ 来源。 |
| G3 GAMMA/REST | 同一外域和 mesh 顺序求两个余下 RHS；以同一 run/config 的 FULL 重放通过记录为入口。 |
| G4 report | 三组相同空气平面及远场算子后处理；核对复振幅线性相加，报告中心和中心邻域，绘图并回填本 Spec。 |

建议程序接口，须先实现再执行，普通 `audit` 不得隐式启动 COMSOL：

```powershell
uv run python -m scripts.analysis.prepare_s5_gamma_origin --config CONFIG --stage audit
uv run python -m scripts.run_main.run_s5_gamma_origin --config CONFIG --stage export --execute
uv run python -m scripts.analysis.prepare_s5_gamma_origin --config CONFIG --stage decompose
uv run python -m scripts.run_main.run_s5_gamma_origin --config CONFIG --stage full --execute
uv run python -m scripts.analysis.analyze_s5_gamma_origin --config CONFIG --stage replay-check
uv run python -m scripts.run_main.run_s5_gamma_origin --config CONFIG --stage components --execute
uv run python -m scripts.analysis.analyze_s5_gamma_origin --config CONFIG --stage report
```

`CONFIG` 为本次 `99_config/03_gamma_origin/<run_id>/replay_config.json` 的绝对路径。
程序组织可按真实复用关系减少文件，但保留无求解审计、导出、重放、门控和后处理职责。
复用 `run_s5_volume_export.py` 的 NativeSampler、已有保存几何加载及兼容 study/export helper。
复用 `analyze_s5_cladding_scattering.py` 的 `plane_farfield`。
不整体调用旧 cladding runner，其历史 preflight 存在旧路径依赖；不复用旧体电流模板作为新 Γ 定义。
记录 `exp(+iωt)` 时间约定和原 `exp(+iq·r)` 空间约定，观察 origin 沿用 A。
主空气平面初值为原 `z=1.1 μm`，81×81 的原 ±10° 观察窗口必须含精确 q=0。
异质 S 是输入面，不直接在 S 上套用均匀介质远场公式。

## 7. 机制判断与验收

令 `q` 为出射平面波矢，`F_s(q)` 为每组同一约定下的复远场振幅，`s` 为 FULL/GAMMA/REST；法向为 `q=0`。
先检查 `F_FULL(q)=F_GAMMA(q)+F_REST(q)`，再比较强度。
中心报告 x/y 振幅的实部、虚部、模和相位，并报告 `|F_REST,y(0)|/|F_FULL,y(0)|`。
不得把独立强度当作相加功率份额，必须保留 `2 Re(F_GAMMA·F_REST*)` 的相干交叉项。

FULL 重现 A、GAMMA 重现中心复振幅且 REST 在中心消去：支持本定义下 Γ 为中心来源。
REST 若存在实质有限中心振幅：不能写“完全来自 Γ”；保留数据、解释范围，不更改 projector 或归一化来获得预期结论。
仅精确中心一致而亮斑邻域不一致：只论证法向振幅，不扩大为整个亮斑。
FULL 重放失败或输出不满足线性关系：归类 implementation-unresolved，不作机制判定。
不以本次实验单独宣称“几乎完全来自 cladding 散射”；那属于下一项结构对照。
不新增零点精度、频率/网格/窗口稳健性研究作为本任务完成前提。

## 8. 资源、结果与安全约束

预算：一次已保存 A 载入/补导出、一个外域 mesh、三个顺序实频 RHS；复用分解可行时复用，否则同一 mesh 顺序求解。
初始沿用 A 的 16 cores、最多 768 GiB 与 24 小时单任务预算；真实可用资源/许可核对后可向下收敛，不增加并行 COMSOL。
不做新本征计算、1141 个逐输入响应或未经授权的结构扫描。
Windows 本机执行，PATH 加 `C:/Program Files/COMSOL/COMSOL63/Multiphysics/bin/win64`；不引入 Slurm、Linux launcher 或远程 scratch。
启动前遵守远端 AGENTS、仿真/绘图合同和相关 Memory，检查本机 Git dirty 状态并保留已有变更。
新 Python 跑聚焦测试、py_compile 与 `git diff --check`；不自动提交。

所有新批次文件位于 R 的既有格式分类下 `03_gamma_origin/<run_id>/`；不新建 S5 顶层运行长目录，不写 `.out` 外其他结果根。
主 PNG 直接平铺 `R/10_overview/`，使用 `s5_gamma_origin_<run_id>_` 唯一前缀，匹配 PDF 遵守既有格式分类。

| 分类 | 产物 |
| --- | --- |
| `99_config` | input_audit、gamma_definition、replay_config、authorization、源文件哈希及 frozen trace mapping |
| `01_results` | closed_surface_EH、folded_cell_profiles、boundary_inputs、三组空气 E/H 和远场复数组 |
| `00_model` | exterior base、FULL/GAMMA/REST 已求解模型 |
| `80_logs` | stdout/stderr、progress、source_identity、linearity_checks、normal_amplitudes.csv |
| `12_reports` | gamma_origin_report.md，区分 supported / not-supported / implementation-unresolved |
| `10_overview` | 边界分解、三组同标度远场、中心复振幅相加、两个角度 cut 及 cladding 响应分布 |

大型 MPH/NPZ 留在 SSH 端，仅明确打包且校验的轻量结果后续导回 Vault。
不得覆盖 `02_gamma_cladding` 或其他已存结果；resume 必须所有 model/mesh/source/frequency 身份一致。
后台启动后核对实际 PID、stderr、progress；记录模型载入/导出/mesh/solve 的真实阶段，不把派发或排队写作已求解。

## 9. 执行与验证回填

### 2026-10-07 14:27（北京时间）：外域网格完成，FULL 已启动

- 实现：scripts/analysis/prepare_s5_gamma_origin.py、scripts/run_main/run_s5_gamma_origin.py、scripts/run_main/run_s5_gamma_replay.py、scripts/analysis/analyze_s5_gamma_origin.py；聚焦测试 tests/test_s5_gamma_origin.py 共 5 项通过，相关入口 py_compile 通过。
- G0/G1 已完成：三组面采样文件均为 1141/1141。E/H 体场重构误差 0，REST 复平均误差约 9e-16；源模型 SHA256 与本批次 input_audit.json 一致。
- 闭合面独立核验见 80_logs/03_gamma_origin/<run_id>/boundary_input_checks.json：E/H 分别重构误差 0；cap–side 切向 E 差异 FULL=2.4274e-12、Γ=1.0075e-12；侧面共享角点 Ez 差异 0；材料界面单侧采样切向 E 差异 FULL=1.7572e-6、Γ=1.0598e-5。cap 上 Qh 对原始 Γ 的面积加权切向 E/H 改变量分别约 0.12868% 和 0.17898%。
- G2 外域网格已完成：1733 个内部域从 EWFD 排除，5333 个外部域保留；119 张内边界由 101 个插值分区恰好覆盖一次。已保存 exterior_partition.mph、exterior_base.mph，算子配置见 99_config/03_gamma_origin/<run_id>/exterior_operator.json。
- 边界选择修复：COMSOL 6.3 Finalize 不支持 intbnd 属性；Java 选择对象须使用 component.selection(tag)，整数参数显式 JInt。位置容差沿用保存模型 selection_tol=1 nm。为排除转角 Ball seed 同时选中相邻面，最终依据实际边界全部顶点的平面和高度归属动态创建 named selections；无硬编码边界编号。增加相邻面误选回归测试。
- 通过受控交互进程复用已载入检查点完成建模，之后直接调用同一 runner 的 solve_cases(model,cfg,'all',client,util)。当前 Python PID=21368、COMSOL PID=928，16 cores；FULL 于 14:27:30 开始编译频域方程。原参数、频率、mesh、cladding 不变。
- 常规恢复命令：.venv/Scripts/python.exe -B -u -m scripts.run_main.run_s5_gamma_replay --config CONFIG --stage all --execute；CONFIG 为第 6 节同一路径。已有有效 mesh/解可恢复，但不得与当前任务并行启动。
- 顺序流程已接入：FULL 求解与保存 → 同坐标 E/H、中心和远场复振幅门控 → GAMMA → REST → 共同远场与复振幅分析 → 五张 PNG/PDF、CSV/NPZ/JSON、中文 gamma_origin_report.md。FULL 门控不通过则停止归因，不调整相位或缩放来通过。
- FULL 重放与三组输出线性关系：等待真实求解结果，尚未验收。中心和亮斑邻域机制：未判定。
- 当前图：10_overview/s5_gamma_origin_<run_id>_01_boundary_inputs.png（已视觉检查）；每行内 FULL/Γ/REST 共用线性色标，cap 与 side 分行标注范围。
- 正式辐射图目标：00_replay_validation、02_common_farfields、03_center_phasors、04_cuts_and_exterior；Ex/Ey 每种偏振内三组共享线性色标，分别标注范围。
- 历史失败日志 replay_attempt1–5 保留；实时状态 execution.json 与 FULL_COMSOL_progress.log 为当前运行依据。前述失败均发生在 FULL 求解之前，不构成物理机制的反例。
