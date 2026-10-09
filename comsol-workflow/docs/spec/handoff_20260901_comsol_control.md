# COMSOL Control 项目交接文档

> 更新时间：2026-09-01
> 用途：供新的 Codex 对话直接接手 `comsol-workflow` 的本机 COMSOL 控制、计算、分析、绘图与结果管理工作。
> 本文总结本对话中已经确认的长期规则、物理语义、程序结构、数据规范、历史决策、当前结果和未完成任务。执行任何新任务前，仍应读取项目根目录的 `AGENTS.md`、`MEMORY.md` 与相关源码，以代码和当前配置的实际状态为准。

## 1. 项目定位

`comsol-workflow` 最初为服务器/cluster 运行而构建，目前的目标已经明确调整为：

- 作为本机 COMSOL 控制器；
- 由 Codex 构建、启动、监控和终止 COMSOL 计算；
- 保留几何构建、网格构建、COMSOL 控制、求解、导出、后处理与绘图等核心能力；
- 不再考虑与服务器、Slurm 或 cluster 的通信与调度；
- 服务器相关内容和其他判定为无关的内容先移入归档目录，最终是否删除由用户决定；
- 归档目录中的程序，除非用户明确点名，否则不要修改。

项目的长期方向是“文件数量能少则少、入口清晰、参数统一、结果可复现、输出目录稳定”。

## 2. 用户协作习惯与执行规则

### 2.1 一般工作方式

- 用户常用“确认”“批准”“开始执行”授权已经讨论清楚的方案；获得批准后直接实施，不重复询问同一事项。
- 用户说“按照这个配置帮我设置参数”时，默认采用最小修改：只修改名称中明确给出的参数，其他运行开关与行为保持不变。
- 用户要求“先确认参数”时，应列出实际将使用的几何、层数、mesh、shift、中心频率、模式数和运行行为，等待确认后再开始昂贵计算。
- 对状态询问要检查进程、日志更新时间、COMSOL 进度与输出文件，不应只根据是否存在进程作判断。
- 用户明确说由其自行启动计算时，只修改配置和代码，不代为启动 COMSOL。
- 用户明确要求终止时，可以结束对应计算进程；应避免误杀无关 COMSOL 会话。
- 用户倾向直接实施，不需要为普通、小范围改动写计划。

### 2.2 昂贵计算前必须重新核对

开始 COMSOL 全带、finite、quarter、strip sweep 或其他昂贵任务前，至少核对：

1. cavity 几何参数；
2. cladding 几何参数；
3. cavity/cladding 层数；
4. mesh size；
5. cladding shift 列表；
6. 中心频率或 eigenfrequency shift；
7. 求解模式数；
8. strip cut 方向或 quarter symmetry ID；
9. 是否续算已有 MPH；
10. 输出目录是否会覆盖、合并或跳过已有参数；
11. 是否存在其他正在运行的 COMSOL 实例以及内存/端口冲突风险。

### 2.3 大改动记录规则

重构程序、创建新的 main 程序、统一配置体系、正式修改输出架构等大改动，必须：

1. 在 `docs/spec` 中创建带时间戳的中文 `plan_execute_YYYYMMDD_*.md`；
2. 实施前写方案概述、目标、范围、兼容性与验证方法；
3. 实施完成后在同一文档记录实际改动、验证结果、遗留问题；
4. 此规则已要求写入 `AGENTS.md`，新对话应检查并遵循。

### 2.4 临时文件规则

- 任务完成后不保留临时脚本、临时预览图、探针文件和无意义缓存。
- 即使过程文件必须暂时保留，也不要放在项目主目录或源码目录。
- 临时内容优先放系统临时目录或正式结果目录中明确的日志/配置子目录。
- 不要为了方便在根目录散落一次性 Python 文件。

### 2.5 Git 规则

- 用户会在需要时要求 commit 或 push；未明确要求时，不擅自提交或推送。
- 工作树中可能存在用户自己的未提交修改，禁止 reset、checkout 或覆盖无关改动。
- 当前对话中使用过远端 `https://github.com/jia-yli/comsol-workflow.git`，分支为 `raywong/finite-cavity`；执行 push 前必须重新检查当前 remote、branch、status。
- Git Credential Manager 曾提示无法使用 `wincredman` 持久保存凭据，但 push 本身可能已经成功并显示 `Everything up-to-date`。
- 当前 Git push 的具体操作方法已要求记录在 `AGENTS.md`；新对话应读取其中的最新命令，不要仅依赖本文的历史描述。
- `blueprints` 后来还被要求作为独立项目推送到 `https://github.com/raywong0131/dxf_blueprints`。处理主仓库与独立 blueprints 仓库时必须先确认工作目录和 remote，避免推错仓库。

## 3. 项目目录与程序分类

### 3.1 `scripts` 的正式分类

`scripts` 下的程序按以下职责组织：

- `scripts/run_main`：核心计算入口；
- `scripts/run_sweep`：调用 COMSOL 的参数扫描与优化；
- `scripts/analysis`：不负责 COMSOL 建模求解，只读取已有结果做分析和绘图；
- `Archive` 或项目实际采用的归档目录：不再主动维护的程序，除非用户明确提及不要操作。

曾经存在 `scripts/sweep`、零散脚本等旧布局；当前状态应以源码和 `scripts/README.md` 为准。

### 3.2 核心 main 程序

#### `run_unit_cell_band.py`

- 负责 cavity/cladding unit-cell 建模与完整能带计算；
- unit-cell 的求解基准频率始终保持 `c_const/1.55[um]`，不读取 finite 中心频率作为其自身基准；
- 完整 unit-cell 计算结束后，默认把 **cavity 的 p2@Γ** 写入统一参数中的中心频率，供 strip 和 finite 使用；
- 不是 cladding 的 p2@Γ。

#### `run_strip_1d.py`

- 负责 1D cut/strip 结构的验证计算；
- cavity 与 cladding 几何、层数、mesh、cladding shift、中心频率从统一配置读取；
- 已确认有限六边形的六个边界只需要 `0°` 和 `60°` 两类 cut：上下边为一类，另外四个侧边为同一类；
- 不能再通过简单旋转 cavity 几何来假装侧边 cut，`60°` cut 必须按真实边界方向构建；
- 已要求取消通过 COMSOL `PlotGroup2D` 和 `Image Export` 弹出/导出额外图像；
- cladding shift 数量大于 3 时，默认自动执行 strip trend 分析。

#### `run_finite.py`

- 核心完整有限尺寸 cavity 计算程序；
- cavity/cladding 几何、层数、mesh、cladding shift、中心频率、模式数由统一配置管理；
- `RESUME_FROM_EXISTING_MPH` 表示是否从已有 MPH 模型继续，不能仅因目录存在就盲目续算。

#### `run_finite_quarter.py`

- 使用第一象限 `x >= 0, y >= 0` 的 1/4 模型计算完整 finite cavity；
- 通过沿 x/y 轴的对称边界条件求解，之后在后处理中利用对称性恢复完整场；
- 保留 `finite_quarter.mph`；不保留不必要的 1/4 场导出，只保留恢复后的完整场结果；
- 最终分析类型与完整 `run_finite` 保持一致；
- 每个 symmetry ID 独立保存；
- 曾创建 `run_finite_quarter_four_modes.py` 等辅助入口，实际职责和是否仍保留需查源码。

### 3.3 Package 标准化

- 三个主程序移动到子目录后曾出现 VS Code“无法解析导入”；根因是脚本式路径注入与 package 布局不一致。
- 已决定标准化为 Python package，并保留可扩展性。
- `scripts/run_main/__init__.py`、`scripts/__init__.py` 的必要性应结合 package 入口和 `pyproject.toml` 判断，不应随意删除。
- 项目支持 `uv run` 和 `pyproject.toml` 中定义的命令入口，例如曾确认 `uv run comsol-finite-quarter` 可以启动相应 entry point。

## 4. 统一参数体系

### 4.1 文件职责

- 统一数据文件名：`scripts/parameter.json`；
- `parameter_config.py` 负责加载、校验、格式化或向各入口提供配置接口；
- 曾讨论把 JSON 和 Python 合并，但最终决定维持 `parameter.json + parameter_config.py` 的现状；
- 修改参数时优先只改 JSON，不要把运行参数重新散落回各 main 程序。

### 4.2 统一管理的参数

统一配置至少包含：

- cavity：`b0`、`eta/η`、`zeta/ζ`；
- cladding：`b0`、`eta/η`、`zeta/ζ`；
- cavity layer count；
- cladding layer count；
- mesh size；
- cladding shift factors；
- finite/strip 中心频率；
- `EIGENMODE_COUNT`；
- quarter 的 `SYMMETRY_IDS` 或相关运行选项。

### 4.3 各程序读取规则

- `run_unit_cell_band`：读取 cavity/cladding 几何与 mesh；无需读取层数和 cladding shift；其 eigenfrequency 基准固定为 `c_const/1.55[um]`。
- `run_strip_1d`：读取几何、两类层数、mesh、cladding shift、中心频率。
- `run_finite` / `run_finite_quarter`：读取几何、两类层数、mesh、cladding shift、中心频率和模式数。
- unit-cell 完整计算结束后，默认用 **cavity p2@Γ** 更新中心频率。

### 4.4 参数值与数据集不能混淆

本对话中计算过多组参数，至少包括：

- `cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5`：当前最新综合图使用的数据集；
- `cav(245-0.96-1.155)_clad(243.7-0.98-0.928)_mesh9`：finite-quarter 输出结构重构时使用的数据集；
- `cav(245-0.96-1.155)_clad(240.2-0.967-1)_mesh9`：后续 cladding 优化及 trend 分析使用的数据集；
- 更早的 `cav(235-0.960-1.156)_clad(231.8-0.960-0.830)` 等历史结果。

不要根据对话历史猜测 `parameter.json` 的当前值；新对话必须直接读取该文件确认。

## 5. 物理语义与模式定义

### 5.1 px、py 与 p1、p2

- 用户规定：所有关于 px/py 的表达与能带标注都按用户定义；
- **红色能带对应 py，蓝色能带对应 px**；
- p1 与 p2 应有严格的对应关系；
- px 与 py 的定义存在坐标/几何方向导致的互换问题，不能把 p1/p2 与 px/py 永久硬编码为同一映射；
- 报告或绘图时，先保持 p1/p2 的严格追踪，再按用户规定解释 px/py。

### 5.2 历史优化目标

本对话中曾完成或讨论多类 cell 优化：

1. cavity 红色 py 的 Q 最大值位于 Γ 点，同时尽量保持其频率，允许约 ±0.1 THz；
2. cladding 蓝色 px@Γ 与优化后的 cavity py@Γ 对齐；
3. 对 ζ=1、η<1 的 cladding：p1/p2 在 Γ 简并，d1/d2 在 Γ 简并；
4. 修正后的目标之一为 cavity p1@Γ 与 cladding p1/p2@Γ 一致；
5. 后续又修正 cladding d1/d2@Γ 需要与 cavity d1@Γ 一致；
6. 曾取消优化程序的“最小截断”限制。

这些是不同阶段的目标，不应混为一个当前目标。再次优化前必须让用户确认当前使用哪一套目标函数。

### 5.3 Strip Γ 点模式选择

strip trend 的模式选择不能只看频率差绝对值，必须：

1. `gamma_subspace_p > 0.99`；
2. 在满足该条件的候选中选择最接近基准频率的模式。

该规则源于历史上 shift=0.09 误选模式的问题。

## 6. Mesh 规则

### 6.1 Finite 手动网格映射

finite cavity 强制使用手动网格，但仍由 `MESH_AUTO_SIZE` 数字控制三项参数：

| Mesh size | 最大单元增长率 | 曲率因子 | 狭窄区域分辨率 |
|---:|---:|---:|---:|
| 1 | 1.30 | 0.20 | 1.00 |
| 2 | 1.35 | 0.30 | 0.85 |
| 3 | 1.40 | 0.40 | 0.70 |
| 4 | 1.45 | 0.50 | 0.60 |
| 5 | 1.50 | 0.60 | 0.50 |
| 6 | 1.60 | 0.70 | 0.40 |
| 7 | 1.70 | 0.80 | 0.30 |
| 8 | 1.85 | 0.90 | 0.20 |
| 9 | 2.00 | 1.00 | 0.10 |

历史上直接设置 `hgradactive` 会触发 COMSOL “未知属性”异常。已批准的修复方向是只设置当前 COMSOL API 实际支持的属性，不要重新引入 `hgradactive`。

### 6.2 内存与 COMSOL Java heap

- 本机内存约 6 TB，用户明确不考虑优化几何表达，倾向提高 `mphserver` Java heap 以充分利用本机内存；
- 但当前对话末尾的 `0xc0000142`/资源不足是 Windows 无法创建新进程，不能简单等同于 COMSOL 求解堆不足；
- 同时运行多个 COMSOL 可以实现，但必须使用独立端口、输出目录和资源配额，并确认内存、CPU、许可证与 mphserver 生命周期。

## 7. Quarter 对称边界条件

### 7.1 定义

`SYMMETRY_IDS` 中 x/y 边界指沿 x/y 轴的两条切割边界，编号为：

| ID | x 轴边界 | y 轴边界 |
|---:|---|---|
| 1 | PEC | PMC |
| 2 | PMC | PEC |
| 3 | PEC | PEC |
| 4 | PMC | PMC |

- 基模 py-mode 的默认情况已确认使用 `SYMMETRY_ID = 1`；
- 默认每次计算只考虑一种边界情况，但程序应保留遍历所有 ID 并与完整结构对应的扩展能力；
- `SYMMETRY_ID = 1` 在当前 COMSOL 设置中实际可能返回两个解：通常 mode0 的 Q 明显较低，mode1 才是目标 py-mode；不要因此擅自改变 ID。
- `EIGENMODE_COUNT` 已要求改为从 `parameter.json` 输入。

### 7.2 是否重新建模

修改 `SYMMETRY_IDS` 是否重新建模取决于程序的模型复用逻辑、输出目录存在性与 `RESUME_FROM_EXISTING_MPH`。新对话应直接检查 `run_finite_quarter.py`，不能仅凭配置名回答。

## 8. 标准输出根目录

`scripts/.out` 下长期只保留三个一级结果目录：

1. `unit_cell_band`；
2. `strip_1d`；
3. `finite_cavity`。

之后的正式计算和分析结果都必须归入这三类，不新增零散一级结果目录。

## 9. 结果文件夹命名

### 9.1 Unit-cell

标准形式：

```text
unitcell_band_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
```

### 9.2 Strip 1D

标准形式：

```text
strip1d_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
```

`20-20` 分别对应 cavity 与 cladding 层数。

### 9.3 Full finite

标准形式：

```text
finite_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
```

### 9.4 Finite quarter

标准形式：

```text
finite_quarter_10-10_cav(245-0.96-1.155)_clad(243.7-0.98-0.928)_mesh9
```

- shift 通常作为结构根目录下的子目录；
- 历史结果中同时存在 `shift0.000` 和 `shiftx0.000_shifty0.000` 两套路径，分析时必须按实际数据目录处理；
- 用户曾明确要求保留 `shift0.050` 这种名称，不要擅自改为 `shift00` 或其他缩写。

## 10. Finite / Quarter 正式输出布局

经过多轮整理，正式框架原则为“以 mode 为主、层级不过度嵌套、结果目录只放对应结果”。典型结构：

```text
finite_quarter_.../
├── run_summary.json
└── shift0.050/
    ├── 00_model/
    │   ├── finite_quarter.mph
    │   └── comsol_progress.log
    ├── 01_results/
    │   ├── mode0/
    │   │   ├── 11_simulation_exports/
    │   │   ├── 12_farfield_FFT/
    │   │   └── 13_lattice_fourier_Hz/
    │   ├── mode1/
    │   └── ...
    ├── 10_overview/
    ├── 80_logs/
    └── 99_config/
```

关键规则：

- 结构根目录不再额外创建 `99_log`；`run_summary.json` 可以直接放根目录；
- `00_model` 放 MPH 和 `comsol_progress.log`，二者始终同目录；
- `01_results` 名称保持不变，避免 `modes` 与 `model` 混淆；
- mode 命名采用 `mode0`、`mode1`……，不使用 `mode00`；若自然编号达到两位数则直接为 `mode11`；
- mode 文件夹下只允许结果子目录，不散放额外文件；
- `11`–`13` 内不再细分 `data/figures/metrics`，避免嵌套过深；
- `10_overview` 存放几何/shift/所有模式共享的汇总 CSV 和 PNG；
- `bulk_cyclic_index.png`、`finite_k_grid_first_hz.png` 等几何结构通用信息只在 `10_overview` 保留，不应在每个 mode 的 `13_*` 中重复；
- `80_logs` 存放日志；
- `99_config` 主要存 JSON 配置和元数据，`objective.json` 放在这里；
- `pari_validation.json` 也应生成并放在 `99_config`；
- `11_simulation_exports` 只放该 mode 的场导出；
- `12_farfield_FFT` 只放该 mode 的远场 FFT 结果；
- `13_lattice_fourier_Hz` 只放该 mode 的晶格 Fourier Hz 结果。

该正式布局曾记录在：

`docs/spec/plan_execute_20260724_finite_output_layout.md`

## 11. 几何结构图规则

finite 与 finite-quarter 均只保留以下三类结构图，并生成在对应 shift 的 `00_model`：

1. `full_finite_simulation`：保持原样；
2. `full_lattice_modulation_vectors`：删除 strip cut 线；cladding shift arrow 显示倍率改为 `×10`；
3. `full_lattice_normal_fan_regions`：删除 strip cut 线。

其余冗余几何图不再生成。Quarter 的图应使用恢复/对应的完整几何表达，而不是只展示 1/4 输出。

## 12. Strip trend 分析规范

- 分析入口曾由 `analyze_strip_gamma_trend` 改名为用户指定的 `srip1d_trend`，但结果目录/文件后来采用 `strip1d_trend`；新对话应检查 `pyproject.toml` 与实际脚本名，必要时统一拼写，不要猜测。
- 不再生成 `*_issues.txt`；
- 在对应 strip 数据根目录下创建 `strip1d_trend/`；
- 结果命名为：

```text
strip1d_trend/strip1d_trend.png
strip1d_trend/strip1d_trend.csv
```

- 当 cladding shift 点数大于 3 时，strip 计算默认自动触发 trend；
- 用户补充新数据后，trend 应重新读取全部现有点并覆盖更新结果。

## 13. Finite trend 与 Flatness trend

- 已对多组 finite-quarter 数据运行 finite trend 和 Flatness trend；
- 用户更新数据点后应重新扫描目录，而不是复用旧缓存；
- 用户曾更新 Flatness trend 的图像输出代码，若输出没有改变，应确认实际调用的是最新 package/entry point，而不是旧脚本或安装缓存；
- trend 分析通常只使用已有结果，不应启动 COMSOL。

## 14. Unit-cell 输出与绘图规范

### 14.1 Unit-cell 结果整理

- 历史 `unit_cell_band/cell_pair_opt_mesh9/best_pair_full_band` 已要求改成标准参数命名；
- unit-cell 结果也应按当前编号体系整理；
- 应生成 `10_overview`；
- `pari_validation.json` 放入 `99_config`，以后也直接生成到此处。

### 14.2 能带绘图约定

- 红色能带 = py；蓝色能带 = px；
- 图题使用 `Cavity` 和 `Cladding`；
- 默认不显示网格；
- 默认不显示 Γ 点竖直虚线；
- 所有可见主/次刻度朝内；
- 字体应适合论文图，标注清楚且略大；
- 图例不要遮挡能带，可放到图外；
- 横坐标不标 `wave-vector path`；
- 当前批准的路径表达为 `0.3MΓ–Γ–0.3ΓK`；
- `0.3` 表示从 Γ 沿 Γ–M、Γ–K 路径取到对应完整路径长度的 30%，不是波矢的绝对值；
- 历史参考图的纵坐标范围为 185–210 THz；
- fundamental py-mode 可用文本标识，但当前最新要求是不在注释中显示具体频率数字。

## 15. 当前最新综合图任务

### 15.1 正式图片

```text
scripts/.out/finite_cavity/
finite_quarter_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
shift0.000/10_overview/unit_cell_band_with_all_modes_Q.png
```

已生成的图由三个面板组成：

1. Cavity unit-cell bands；
2. Cladding unit-cell bands；
3. 全部 finite-quarter 模式的频率–Q 散点。

现有图的主要属性：

- 3334 × 1867，300 DPI；
- 包含全部 80 个 finite-quarter 模式；
- 每个 symmetry sector 原始数据有 20 个模式，但图中不区分 ID，普通点统一为蓝色；
- Q 为横坐标且采用线性轴；
- 三个面板共享纵坐标频率；
- 频率范围 185–210 THz；
- band path 范围为 ±0.3；
- 无网格、无 Γ 竖直虚线；
- 基模为全局最大 Q 点，约：
  - `f = 197.982801277 THz`
  - `Q = 8648.382282`
- 以红点和贯穿三个面板的红色水平虚线高亮。

### 15.2 数据来源

Cavity bands：

```text
scripts/.out/unit_cell_band/
unitcell_band_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
10_overview/cavity/selected_bands.csv
```

Cladding bands：

```text
scripts/.out/unit_cell_band/
unitcell_band_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
10_overview/cladding/selected_bands.csv
```

Finite modes：读取下列四个目录中的 `10_overview/eigenfrequencies.csv`：

```text
scripts/.out/finite_cavity/
finite_quarter_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
shiftx0.000_shifty0.000/symmetry_1_xPEC_yPMC/10_overview/eigenfrequencies.csv
```

以及同层级的 `symmetry_2_*`、`symmetry_3_*`、`symmetry_4_*`，共 80 个解。

参考能带图：

```text
scripts/.out/unit_cell_band/
unitcell_band_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5/
10_overview/unit_cell_band_comparison_ylim185_210_THz_halfM_gamma_halfK.png
```

### 15.3 当前唯一未完成的图像修改

当前注释为两行：

```text
fundamental py-mode
(197.98 THz)
```

用户最新要求：

- 删除 `(197.98 THz)`；
- 只保留 `fundamental py-mode`；
- 红点、红色水平虚线、80 个模式、线性 Q 轴以及所有其他版式不变；
- 覆盖原 PNG；
- 使用已有 CSV 重新绘图，不启动 COMSOL；
- 不使用生成式图像编辑改动科学图。

同一目录当前可能还有：

```text
unit_cell_band_with_all_modes_Q.txt
```

新对话应先检查该 TXT 是否包含绘图数据、说明或可复现脚本，再决定复用方式。

## 16. 当前系统阻塞

本对话末尾，本地 Windows 无法创建新的子进程：

- `cmd.exe` 返回 `0xc0000142`；
- PowerShell 返回 `0xc0000142`；
- Git Bash 无法正常执行命令；
- 尝试绕过 shell 直接启动项目 Python时，系统报告：

```text
Not enough memory resources are available to process this command.
```

因此这不是绘图脚本、PowerShell 或 cmd 的单独问题，而是 Codex/Windows 宿主当前无法分配创建进程所需资源。

建议恢复步骤：

1. VS Code 执行 `Developer: Reload Window`；
2. 若无效，完全退出全部 VS Code 窗口；
3. 结束残留的 `Code.exe`、Codex、Python、COMSOL、Java/mphserver 进程；
4. 重新打开项目并建立新 Codex 对话；
5. 若普通终端仍无法启动，重启 Windows。

仅新建聊天而不重启扩展宿主，可能无法解决 `0xc0000142`。

## 17. 新对话接手顺序

1. 读取根目录 `AGENTS.md` 与 `MEMORY.md`；
2. 读取本文件；
3. 用简单 PowerShell/Python 命令确认子进程已经恢复；
4. 检查 `git status`，保护所有用户已有改动；
5. 检查 `scripts/parameter.json`，记录当前参数但不要为当前绘图任务修改它；
6. 查看 `unit_cell_band_with_all_modes_Q.png` 和同目录 TXT；
7. 找到或重建绘图逻辑，只删除注释中的频率数字；
8. 使用四个 symmetry CSV，验证图中仍有 80 个解；
9. 验证 Q 轴为线性、ID 未区分、红点与水平参考线保留；
10. 覆盖正式 PNG；
11. 清理临时脚本、预览和缓存；
12. 向用户提供正式图片的绝对可点击路径；
13. 不 commit、不 push、不启动 COMSOL，除非用户后续明确要求。

## 18. 执行时的优先级

发生冲突时按以下优先级处理：

1. 用户当前对话中的最新明确要求；
2. 根目录 `AGENTS.md`；
3. 根目录 `MEMORY.md`；
4. 已批准的 `docs/spec/plan_execute_*.md`；
5. 本 handoff 文档；
6. 历史输出目录和旧脚本中的行为。

本文用于交接，不应取代对当前源码、配置和结果文件的实际检查。
