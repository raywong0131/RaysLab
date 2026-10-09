# `scripts/` 本机工作流

本目录负责参数与任务编排，底层几何、COMSOL 控制和数值处理由
[`comsol_workflow/`](../comsol_workflow/) 提供。

- [`run_main/`](run_main/)：band pair/solo、unit-cell 2D、strip、full/quarter finite 和四扇区批处理。
- [`run_sweep/`](run_sweep/)：COMSOL 参数扫描、优化和 refinement 求解。
- [`analysis/`](analysis/)：已有结果分析、重绘，以及显式请求的迁移/恢复工具。
- [`parameter.json`](parameter.json)：按共享数据与工作流分区的 schema v2 参数；由 [`parameter_config.py`](run_main/parameter_config.py) 校验，字段归属见 [参数指南](../docs/parameters.md)。

从仓库根目录运行下面的文件入口；安装命令以
[`pyproject.toml`](../pyproject.toml) 的 `[project.scripts]` 为准。
参数、模式命名、对称扇区和结果布局统一参见
[项目说明](../README.md#configuration) 与
[仿真合同](../docs/agent-simulation.md)，本页只维护脚本用途和特有选项。

计算结果固定在 `scripts/.out/` 的 `unit_cell_band/`、`unit_cell_2D/`、
`strip_1d/`、`finite_cavity/` 四个目录。分析优先使用已保存数据；
refinement、恢复或迁移工具的具体副作用须以入口及其合同为准。

## 活动 Python 程序

| 文件 | 功能与作用 | 使用方法 |
| --- | --- | --- |
| [`run_main/run_band_pair.py`](./run_main/run_band_pair.py) | 在一致的 Gamma-M/Gamma-K 网格上计算 cavity 与 cladding unit cell；同一组几何参数复用一个 COMSOL Model，后续 k 点只更新 Floquet 参数并重求，mesh 改变时只重构 mesh；选择并跟踪 2 条 p band 和 2 条 d band，输出模式成分、并排能带以及两方向 p-band Q-k 拟合。求解数量不足以覆盖某段目标带时保留 `unmatched` 证据、断开对应图线并继续其余输出。完整主流程成功后把 cavity `p2@Gamma` 回写为共享中心频率。 | 先核对 `parameter.json` 的 cell、mesh、`unit_cell_band.eigenmode_count`，以及文件顶部的频率窗口、band 和 Q-fit 取点，再运行 `uv run python scripts/run_main/run_band_pair.py` 或 `uv run comsol-band-pair`。其求解基准不读取共享中心频率；若 Gamma 点本身没有求出两条 p 和两条 d，程序会保留诊断性错误。 |
| [`run_main/run_band_solo.py`](./run_main/run_band_solo.py) | 只读取并计算共享配置中的 cavity unit cell，复用 pair 的 k 点求解、band tracking、Q 拟合和绘图实现；不创建 cladding case 或比较图。 | 运行 `uv run python scripts/run_main/run_band_solo.py` 或 `uv run comsol-band-solo`；当前参数输出到 `scripts/.out/unit_cell_band/unit_cell_(245-0.96-1.156)_mesh5/`。 |
| [`run_main/run_strip_1d.py`](./run_main/run_strip_1d.py) | 构建 0°/60° 物理边界方向的 1D cut，验证 cavity/cladding 组合并扫描远端目标 cladding shift；逐层幅值与 finite 共用 `comsol_workflow.cladding_shift_profile`，程序将完整晶格和位移向量刚性变换到 cut 局部坐标后裁剪，并把结果最终化为 `cut_XX/shiftX.XXX` 的 mode-centric layout v3；`01_results` 仅保留 `gamma_subspace_p > 0.9` 的模式，并为其导出 `cutline_Wem.png`。场图读取 COMSOL 数据后由 Python 生成，不调用 COMSOL `PlotGroup2D` 或 `Image Export`。 | 核对 `parameter.json` 的 cell/layers/mesh/中心频率/shift/profile，以及 `RUN_STRIP_CUT`、固定关闭的 `RUN_DIAGNOSTIC_PLOTS` 和 `strip_1d` 中的模式数、`cut_angle_deg`、kx，然后运行主脚本。`RUN_GAMMA_TREND_AFTER_SCAN=None` 时，shift 数量大于 3 会在扫描完成后自动运行 trend；可临时设为 `True/False` 覆盖。Full finite 请单独使用 `run_finite.py`。 |
| [`run_main/run_finite.py`](./run_main/run_finite.py) | 按 `finite_cavity.geometry=hex\|square` 构建并求解有限 cavity；按默认共享 factor+ratio 扫描或非默认显式 x/y pair 解析 case，再应用 shape-aware `groups` 或 `ellipse` 位移。Hex 保留 cyclic quotient Fourier；square 只对完整 cavity parent cell 执行 arbitrary finite DFT。默认在有效模式中选择与目标频率绝对差最小的唯一 mode，只让该 mode 进入 far-field、Fourier、标准 `intensity_maps_normH.png` 和评分。 | 核对 `parameter.json` 的 `finite_cavity.geometry`、cell/layers/mesh/`finite_cavity.eigenmode_count`/中心频率/shift 输入模式/geometry/profile、单模选择和各 `RUN_*` 开关，然后运行 `uv run python scripts/run_main/run_finite.py`。 |
| [`run_main/run_finite_quarter.py`](./run_main/run_finite_quarter.py) | 复用 finite 的完整 geometry plan 和 shift 分类，仅构建并求解 `x >= 0, y >= 0` 的第一象限；hex 使用 quarter hexagon，square 使用 quarter rectangle。通过 `SYMMETRY_IDS` 选择沿 x/y 轴的 PEC/PMC 组合，按场的奇偶性恢复四象限完整场，再执行相同的 shape-aware far-field 与 Fourier 后处理。 | 核对 `parameter.json` 的 `finite_cavity.geometry`、shift 输入模式、geometry、`quarter_symmetry_ids`、恢复和后处理开关，然后运行 `uv run python scripts/run_main/run_finite_quarter.py` 或 `uv run comsol-finite-quarter`。 |
| [`run_main/run_finite_quarter_all_symmetries.py`](./run_main/run_finite_quarter_all_symmetries.py) | 每个 shift pair 只构建一个 quarter COMSOL model 和一次 mesh，固定按 symmetry ID `1 -> 2 -> 3 -> 4` 切换内部 PEC/PMC feature 并求解；每轮复制、校验和保存 `solsymN/dsetsymN` 解副本后才进入下一轮。四轮全部 eigenpair 汇入同一个 `01_results`，mode 目录以 `modeN_<id>_xPEC|PMC_yPEC|PMC` 消歧；`10_overview/f_Q.png` 用统一点样式和线性双轴显示完整频率-Q 分布，不显化 symmetry 身份。 | 参数统一来自 `parameter.json` 的 finite geometry/layers/mesh/中心频率/`finite_cavity.eigenmode_count` 和后处理开关，不读取 unit-cell 目标。series 以 `finite_quarter_all_<cavity>-<cladding>...` 命名。先运行 `--preflight-only` 检查几何、四种边界顺序、预计 `4N` 模式、唯一 MPH 与输出冲突；确认后去掉该参数。可用 `--output-dir <series-root>` 指定非冲突输出。 |
| [`run_sweep/scan_unit_cell_p2_gamma.py`](./run_sweep/scan_unit_cell_p2_gamma.py) | 复用 unit-cell 主入口，按 coarse/fine/final 阶段扫描 cavity `zeta`，寻找目标 p band 的最大 Q 位于 Gamma 附近的候选；可复用已完成 k 点。 | 核对 `B0_NM`、`ETA`、`ZETA_*`、k 点、gap 与模式阈值，然后运行 `uv run python scripts/run_sweep/scan_unit_cell_p2_gamma.py`。 |
| [`run_sweep/cladding_bandgap_alignment.py`](./run_sweep/cladding_bandgap_alignment.py) | 复用 unit-cell 计算，按 `b0 -> zeta -> eta` 优先级优化 cladding；首要目标是 Gamma 目标 p 模式与 cavity 对齐，midgap 为软目标。 | 核对 cavity reference、`B0_*`、`ZETA_*`、`ETA_*`、容差和输出目录，然后运行 `uv run python scripts/run_sweep/cladding_bandgap_alignment.py`。 |
| [`run_sweep/optimize_cavity_cladding_cells.py`](./run_sweep/optimize_cavity_cladding_cells.py) | 使用 mesh9 联合优化 cell pair：先在 cavity 的二维 `eta/zeta` 网格中寻找红色 `py` 的 Q 最大值位于 Gamma 且频率漂移受限的候选，再按 `b0 -> zeta -> eta` 对齐 cladding 蓝色 `px@Gamma`，最后执行完整共同 k 路径验证；所有候选均可按 k 点恢复。 | 核对 cavity/cladding 范围、`CAVITY_REFERENCE_FREQUENCY_THZ`、频率与对齐容差、screen/final k 点及 `OUT_DIR`，然后运行 `uv run python scripts/run_sweep/optimize_cavity_cladding_cells.py`。 |
| [`run_sweep/optimize_cladding_p2_d2_alignment.py`](./run_sweep/optimize_cladding_p2_d2_alignment.py) | 固定 `cav(245-0.960-1.156)` 与 mesh5，从 `clad(242-0.980-0.930)` 出发，仅计算 Gamma 点；按 `zeta -> eta -> b0` 快速对齐 cladding `p2` 与 cavity `p2`，并以 cavity `d1`/cladding `d2` 频差作为软约束。 | 先运行 `uv run python scripts/run_sweep/optimize_cladding_p2_d2_alignment.py --dry-run` 核对范围、共享 `unit_cell_band.eigenmode_count`、固定求解基准和首批 zeta 候选；确认后去掉 `--dry-run`。共享参数文件不会被改写。 |
| [`run_sweep/optimize_finite_quarter_q.py`](./run_sweep/optimize_finite_quarter_q.py) | 在批准的 20–20、mesh9 finite-quarter 配置下固定 `SYMMETRY_ID=1` 和 cladding `eta`，按 `zeta -> b0` 方向搜索用户定义 py-mode 的 Q；每个候选只读取两个本征解的频率/Q/有效性，不生成场、几何、far-field 或 Fourier 图。 | 先运行 `uv run python scripts/run_sweep/optimize_finite_quarter_q.py --dry-run` 核对基线与首批候选；确认后运行不带 `--dry-run` 的同一命令。结果写入 `scripts/.out/finite_cavity/finite_quarter_q_optimization_.../10_overview/`，已完成参数自动复用。 |
| [`analysis/srip1d_trend.py`](./analysis/srip1d_trend.py) | 读取已有 strip 扫描结果；先筛选 `gamma_subspace_p > 0.9`，再选择最接近目标频率的模式，在 cut 目录的 `strip1d_trend/` 子目录中生成 selected-modes CSV 和 trend 图，不运行 COMSOL，也不生成 issues TXT。 | 运行 `uv run python scripts/analysis/srip1d_trend.py --out-root scripts/.out/strip_1d/<series-dir>/cut_60`；可用 `--target-frequency` 覆盖频率或用 `--selected-modes-csv` 直接重绘。 |
| [`analysis/migrate_strip_output_layout_v3.py`](./analysis/migrate_strip_output_layout_v3.py) | 将已完成的 layout v2 strip case 预检并迁移到 `cut_XX/shiftX.XXX`；补齐旧规则跳过的 Fourier/gamma 评分，从保存的 MPH 只读取 `ewfd.Wav` 生成 `cutline_Wem.png`，再按严格阈值收敛 `01_results`。不会重新 mesh 或 eigensolve。 | 先传入一个或多个 series 路径做默认 dry-run；确认后增加 `--execute`。目标 `cut_XX` 已存在或结果不完整时会拒绝执行。 |
| [`analysis/finite_trend.py`](./analysis/finite_trend.py) | 读取 full finite 或 finite-quarter 系列；合并 `run_summary.json` 与系列根目录的全部 scalar `shiftX.XXX` 或 directional `shiftxX.XXX_shiftyY.YYY` case，筛选 `gamma_p_px_weight_fraction > 0.80` 后选择最接近目标频率的 py-mode，生成 selected-modes CSV、py-mode 权重/频率双轴图和同一目标模式的 Q-shift 图。Directional CSV 同时记录 x/y factor，固定 ratio 扫描以 x factor 为横轴，不同 y/x ratio 分线显示；全程不运行 COMSOL。 | 运行 `uv run python scripts/analysis/finite_trend.py --out-root scripts/.out/finite_cavity/<series-dir>`；用 `--y-over-x-ratio 2` 可只保留指定 directional ratio（同时保留零 shift 基线），多 symmetry quarter 系列需传 `--symmetry-id`，也可用 `--selected-modes-csv` 同时重绘两张图。 |
| [`analysis/finite_flatness_trend.py`](./analysis/finite_flatness_trend.py) | 严格消费 `finite_trend.csv` 已选模式，从保存的 lattice-Fourier NPZ、`Hz_center.parquet` 和 case config 计算 cavity flatness 与 cladding guardrail。标准全结构 cell-map 为 `cavity_cladding_intensity_maps.png`：逐 case 以 cavity 中心 cell 强度归一化，使用共享线性 `0–2` 色标、`1` 为中点、顶端 `2+`；`intensity_maps_normH.png` 保留为逐字节相同的兼容别名。 | 运行 `uv run python scripts/analysis/finite_flatness_trend.py --selected-modes-csv scripts/.out/finite_cavity/<series-dir>/finite_trend/finite_trend.csv`；输出位于同目录的 `flatness_trend/`，不启动 COMSOL。 |

## 推荐顺序

新增统一入口由 [pyproject.toml](../pyproject.toml) 注册，保留原任务参数解析器：

| 入口 | 功能与边界 |
| --- | --- |
| `comsol-finite-quarter-all` | 四种对称边界完整谱；`--preflight-only` 只预检。 |
| `comsol-boundary-analysis` / `comsol-boundary-scan` | 边界贡献分析与七种扫描；求解阶段可启动 COMSOL。 |
| `comsol-dipolar-analyze <stage>` | 15 个保存数据分析阶段；`--describe-inputs` 只打印来源，`-- --help` 查看阶段参数。 |
| `comsol-dipolar-reference` / `comsol-dipolar-volume-export` / `comsol-dipolar-boundary-replay` | 周期参考、已有解导出和边界重放；`--describe-inputs` 不计算，实际运行按原阶段授权边界执行。 |

共享数值所有者和完整新旧名称见 [整合 Spec](../docs/spec/plan-execute_20261009_workflow_integration.md)。
来源 JSON、输出兼容规则与离线验证命令见 [根 README](../README.md#configuration)。
规范实现使用 `boundary` / `dipolar` 能力名称，旧 S4/S5 Python 入口继续可用。

```text
run_band_pair
  -> unit-cell 参数扫描（按需）
  -> run_strip_1d + Gamma trend
  -> run_finite 或 run_finite_quarter + finite 后处理
```

Quarter 对称扇区和四扇区批处理的区别见
[finite 工作流说明](../README.md#finite-cavity-and-symmetry-sectors)。

所有命令默认从仓库根目录执行。需要 COMSOL 的脚本必须能在当前 PowerShell
会话中找到本机 COMSOL 6.3；纯分析脚本优先复用 `scripts/.out/` 中已保存的
CSV、JSON、parquet 和场数据，避免重复 eigensolve。

启动昂贵计算前至少核对：模式命名对应、cavity/cladding 参数、结构尺寸、
cladding shift 远端目标及逐层 profile、mesh 的实际映射、中心频率或 k/扫描范围、缓存/恢复策略、输出
目录、后处理开关、内存和 license。
