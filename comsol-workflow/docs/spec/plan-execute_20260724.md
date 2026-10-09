# 主工作流统一参数配置：计划与执行记录

日期：2026-07-24  
状态：已完成

## 实施前方案概述

### 目标

为 `scripts/run_main/` 下的四个主入口建立单一参数源，避免 cavity、cladding、
结构层数、mesh、中心频率和 cladding shift 在不同脚本中重复定义并产生漂移。
共享数据使用 `scripts/parameter.json` 保存；Python 代码负责加载、验证、派生内部
Fourier 参数和生成一致的结果目录标签。

涉及入口：

- `run_unit_cell_band.py`
- `run_strip_1d.py`
- `run_finite.py`
- `run_finite_quarter.py`

### 共享参数与读取规则

`parameter.json` 保存：

- cavity 紧凑几何参数 `(b0_nm, eta, zeta)`；
- cladding 紧凑几何参数 `(b0_nm, eta, zeta)`；
- `cavity_layers` 和 `cladding_layers`，均为正整数；
- `mesh_size`；
- `center_frequency_thz`；
- `cladding_shift_factors`，始终为非空列表。

读取矩阵：

| 参数 | Unit cell | Strip 1D | Finite | Finite quarter |
| --- | --- | --- | --- | --- |
| cavity 几何 | 读取 | 读取 | 读取 | 读取 |
| cladding 几何 | 读取 | 读取 | 读取 | 读取 |
| cavity/cladding 层数 | 不读取 | 读取 | 读取 | 通过 Finite 读取 |
| mesh size | 读取 | 读取 | 读取 | 读取 |
| center frequency | 不作为求解基准 | 读取 | 读取 | 读取 |
| cladding shift | 不读取 | 读取 | 读取 | 读取 |

Unit-cell eigensolver 的基准频率永久保持
`c_const/1.55[um]`，不读取共享中心频率。完整 unit-cell band 成功后，程序默认
选择 cavity 在 Gamma 点严格定义的 `p2`（两条已选 p band 中 Gamma 频率较高者），
并把该频率回写为 `parameter.json.center_frequency_thz`。可复用的底层 band 函数、
参数扫描、单点计算和仅重绘操作不得触发回写。

### 参数转译与命名

共享配置统一执行既有紧凑参数转译：

```text
r_f0        = eta
b_square_f0 = (b0_nm / 230.0)^2
b_square_f3 = (zeta^2 - 1.0) / 2.0
```

公共目录标签由实际解析值生成：

```text
cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
```

strip、finite 和 finite-quarter 的单个 case 继续显式包含 `shift0.XXX`。结构层数
写入系列目录名；strip cut angle、quarter symmetry IDs、unit-cell k 路径、模式数量
和后处理开关仍属于各工作流，不进入共享参数文件。

系列目录命名为：

```text
finite_<cavity_layers>-<cladding_layers>_<公共参数标签>
strip1d_<cavity_layers>-<cladding_layers>_<公共参数标签>
```

例如层数均为 20 时：

```text
finite_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
strip1d_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh9
```

系列目录下 finite case 使用 `shift0.XXX`；strip case 使用
`cut<angle>_shift0.XXX`。已有结果目录不迁移。

### 架构与兼容策略

1. 新增轻量共享配置模块，禁止导入 `mph` 或启动 COMSOL。
2. 加载阶段严格验证 schema、数值有限性、正值约束、层数正整数、mesh 取值范围
   和 shift 列表。
3. 四个入口保留现有常量别名，以最小化几何和后处理函数的改动；别名改由共享
   配置派生，避免一次性重写大量稳定逻辑。
4. 每次输出的 `config.json` 记录紧凑参数、内部派生参数、公共目录标签和中心频率。
5. Unit-cell 回写前确认完整计算成功、cavity `p2@Gamma` 存在且频率有效，并检查
   本次运行使用的几何与 mesh 仍和磁盘上的 `parameter.json` 一致，防止长任务结束
   后覆盖用户中途修改的新配置。
6. 回写采用同盘安全替换；过程文件只放在
   `scripts/.out/unit_cell_band/.tmp/`，成功或失败后清理。

### 验证计划

- 配置文件读取、验证、内部参数转译和标签格式单元测试；
- 四个主入口读取同一几何与 mesh 的契约测试，以及 strip/finite 读取共享层数测试；
- unit-cell 固定求解基准频率测试；
- cavity `p2@Gamma` 识别、只更新中心频率及配置冲突保护测试；
- strip/finite/quarter 的 shift、中心频率和目录命名回归测试；
- `py_compile`、相关 pytest、完整 pytest 和 `git diff --check`；
- 不启动昂贵 COMSOL eigensolve。

### 回滚方案

若迁移失败，可恢复各入口原有常量定义并移除共享加载；已有 `scripts/.out/` 结果
不移动、不覆盖、不删除。新目录名与旧结果目录并存，避免破坏已完成计算。

## 执行完成后记录

### 实际改动

1. 新增 `scripts/parameter.json`，作为四个主入口的唯一共享参数文件。初始
   cavity/cladding/mesh 为 `245-0.96-1.156`、`242-0.98-0.93` 和 mesh9；初始
   cavity/cladding 层数均为 10；初始
   中心频率采用已有完整 unit-cell 结果中的 cavity
   `p2@Gamma = 198.4493726243485 THz`；共享 shift 为 `[0.03]`。
2. 新增轻量 `scripts/run_main/parameter_config.py`：
   - 严格加载和验证 JSON schema、cell 参数、结构层数、mesh、中心频率和 shift 列表；
   - 统一生成内部 Fourier 参数、公共目录标签和带结构层数的系列目录标签；
   - 提供同盘安全替换函数，只更新 cavity `p2@Gamma` 中心频率；
   - 回写前检查本次 unit-cell 使用的 cell/mesh 与磁盘配置是否一致；
   - 临时文件位于 `scripts/.out/unit_cell_band/.tmp/` 并在操作结束后清理。
3. `run_unit_cell_band.py` 改为读取共享 cell/mesh，输出目录由共享标签生成；
   eigensolver shift 保持 `c_const/1.55[um]`。只有正式 `main()` 在完整 band、拟合
   和绘图全部成功后回写 cavity `p2@Gamma`；`run_band_workflow()` 本身不写共享
   参数，保证 sweep 和复用调用无副作用。
4. `run_strip_1d.py` 改为读取共享 cell、cavity/cladding 层数、mesh、中心频率和
   shift 列表。系列目录使用
   `strip1d_<layers>_<公共参数标签>`，其下 case 目录使用
   `cut<angle>_shift0.XXX`。迁移时原本地默认层数 20/20 改为共享文件中的 10/10。
5. `run_finite.py` 的活动配置改由共享参数派生；保留原 preset 解析函数以兼容
   已有测试和调用，但活动几何的 Fourier 转译只使用共享配置模块。中心频率、
   shift、mesh、cavity/cladding 层数和目录标签均来自 `parameter.json`。系列目录
   使用 `finite_<layers>_<公共参数标签>`，其下 case 目录使用 `shift0.XXX`。
6. `run_finite_quarter.py` 继续复用 `run_finite.py` 的稳定几何和后处理实现，因而
   自动使用同一个共享参数对象和结构层数；quarter 专属 `SYMMETRY_IDS` 保持本地
   配置。其系列目录同样显式记录层数和公共参数标签。
7. `pyproject.toml` 把 `parameter.json` 声明为 `scripts` package data；README、
   `AGENTS.md` 和 `MEMORY.md` 已同步共享参数的数据流与维护规则。
8. 原有 `scripts/.out/` 结果未移动、覆盖或删除；新的目录标签与旧目录并存。

### 验证结果

- `py_compile`：共享配置模块和四个主入口通过；
- 聚焦回归：106 项通过；
- 全仓库 pytest：219 项通过；
- pytest 仍报告 3 条既有 `PytestReturnNotNoneWarning`，与本次改动无关；
- 共享 JSON 加载、结构层数验证、参数转译、标签、中心频率原子回写、未知字段
  保留、临时文件清理和长任务配置冲突保护均有自动化测试；
- Unit-cell 不消费结构层数：层数在 unit-cell 长任务期间变化不会阻止中心频率
  回写，且新层数会原样保留；strip、finite 和 finite-quarter 均从共享配置读取层数；
- finite、strip 和 quarter 系列目录及其 shift/cut 子目录通过导入与命名回归测试；
- Unit-cell 固定基准频率及严格 cavity `p2@Gamma` 提取均有自动化测试；
- `parameter.json` 经 JSON 解析和主入口导入验证；
- `git diff --check` 通过；
- 本次配置重构未启动 COMSOL 或 eigensolve。

### 遗留问题与最终状态

实现已完成，可作为后续主工作流的统一参数入口。旧结果目录不会被自动识别为
新参数 case，也不会自动迁移。Unit-cell 自动回写已通过纯 Python 测试，首次正式
运行时仍应按昂贵计算 preflight 核对当前 `parameter.json`、输出目录和 cavity
`p2` 定义；运行成功后检查文件中的 `center_frequency_thz` 与本次
`cavity/selected_bands.csv` 一致。

## 三类输出根目录补充实施（2026-07-24）

### 实施前约束

为避免不同主程序与 sweep 在 `scripts/.out/` 下持续增加并列目录，输出一级目录
固定且只能保留三项：

```text
unit_cell_band/
strip_1d/
finite_cavity/
```

- unit-cell band、cavity/cladding 优化和全部 unit-cell sweep 进入
  `unit_cell_band/`；
- strip case、扫描汇总和 trend 分析进入 `strip_1d/`；
- full finite、quarter finite 和 finite 专属分支进入 `finite_cavity/`；
- 参数回写临时文件放在 `unit_cell_band/.tmp/` 并在结束后清理，不能在 `.out/`
  下临时创建第四个目录；
- 现有结果只做同盘目录迁移，不删除或覆盖。

### 实际改动

- 在 `parameter_config.py` 中集中定义三个输出根目录和唯一允许集合；
- 四个主入口、三个 unit-cell sweep、strip trend 分析及 finite 的旧 unit-cell
  验证分支均改用集中路径；
- finite-quarter 归入 `finite_cavity/`，不再使用独立 `finite_quarter/` 一级目录；
- 更新 README、AGENTS 和 MEMORY，使后续新增程序继续遵守三目录约束；
- 现有 `scripts/.out/` 已收敛为三项，正式计算结果保留。

### 验证结果

- 主入口和 sweep 的无 COMSOL 导入验证通过；
- 路径契约聚焦测试：115 项通过；
- 全仓库 pytest：220 项通过，仍有 3 条既有 `PytestReturnNotNoneWarning`；
- 测试临时文件在 `unit_cell_band/.tmp/` 内生成并在结束后清理；
- 磁盘审计确认 `scripts/.out/` 一级目录有且只有 `unit_cell_band/`、
  `strip_1d/`、`finite_cavity/`；
- `git diff --check` 通过；
- 未启动 COMSOL 或 eigensolve。

## Finite 结果按模式归档重构（2026-07-24）

### 实施前方案概述

#### 目标与范围

将 full finite 与 finite-quarter 的最终结果从“按后处理方法集中、同一模式分散”改为
“按 shift case 和 mode 集中”。先迁移用户明确指定的
`finite_quarter_10-10_cav(245-0.96-1.155)_clad(243.7-0.98-0.928)_mesh9`，
再修改两个 finite 主入口。Strip 与 unit-cell 输出不在本次范围内。

#### 最终目录契约

默认 quarter symmetry ID 1 和 full finite 使用如下布局：

```text
<series>/
  run_summary.json
  shift0.000/
    00_model/
      finite_quarter.mph 或 finite_cavity.mph
      comsol_progress.log
      几何/模型诊断图（若该工作流生成）
    01_results/
      mode0/
        10_overview/
        11_simulation_exports/
        12_farfield_FFT/
        13_lattice_fourier_Hz/
      mode1/
        ...
    10_overview/
      shift 级 CSV 与 PNG
    80_logs/
      其他 shift 级日志
    99_config/
      config.json
      objective.json
      其他 shift 级 JSON
```

- mode 目录使用不补零的自然编号：`mode0`、`mode1`、...、`mode9`、
  `mode10`、`mode11`；其根目录下只允许存在目录，不允许直接存在文件。
- `.mph` 只保存一份，位于 shift case 的 `00_model/`；`comsol_progress.log`
  始终与其对应的 MPH 位于同一目录。
- 所有 mode 目录统一位于 shift case 的 `01_results/`；保留 `01_results`
  这个名称，不改为容易与 `model` 混淆的 `modes`。
- mode 专属的摘要、场、far-field 和 finite-lattice Fourier 分别进入
  `10_overview/`、`11_simulation_exports/`、`12_farfield_FFT/` 和
  `13_lattice_fourier_Hz/`。
- `11~13` 内的文件直接平铺，不再增加 `data/`、`metrics/`、`figures/` 子目录，
  避免文件夹嵌套增加查看复杂度。与 mode 无关的几何、BZ、基底和采样网格信息
  只在 shift 的 `10_overview/` 保留一份。没有产生结果时不创建空分析目录。
- 结果摘要属于对应分析，例如 far-field 的摘要命名为
  `12_farfield_FFT/farfield_summary.json`，不归入 `99_config/`。
- shift 级 CSV/PNG 进入 `10_overview/`；其他日志进入 `80_logs/`；
  shift 级 JSON 配置与元数据进入 `99_config/`，`objective.json` 明确保留
  在此目录中。
- series 的 `run_summary.json` 直接位于 series 根目录，不额外创建 series `99_log/`。
- 两位数字前缀具有唯一且稳定的功能语义，不得把同一编号复用于不同类型目录。
- shift case 名称保持原有三位小数格式：`0 -> shift0.000`、
  `0.05 -> shift0.050`、`0.07 -> shift0.070`、`0.075 -> shift0.075`。
- quarter 默认仅计算 `SYMMETRY_IDS=[1]` 时使用上述直接布局；启用非默认或多个 ID 时，
  case 名增加 symmetry 标识，防止不同边界条件互相覆盖。

#### 实现与兼容策略

1. COMSOL 导出和现有后处理先写入 case 专用的任务级 staging 目录；全部成功后，
   由一个可测试的整理函数同盘移动到 mode-centric 最终布局并清空 staging。
2. 模型 checkpoint 和 progress log 从开始即写入 `00_model/`；case config 写入
   `99_config/`，以保证 resume 路径稳定。
3. 迁移正式旧结果时先核对源/目标绝对路径、拒绝覆盖已有目标，并在移动完成后验证
   mode 根目录无文件、源目录无未识别残留，再删除空旧目录。
4. 继续保留聚合 summary，便于跨模式比较；同时在各 mode 的 `10_overview/` 中保存该模式
   对应的 eigenfrequency/score 行。
5. 不重新启动 COMSOL，不重算现有结果；迁移只进行同盘目录整理和小型共享图复制。

#### 验证与回滚

- 为 shift 命名、mode 目录生成、旧布局迁移、full/quarter 路径和 mode 根目录约束增加测试。
- 运行 finite field export、far-field、finite-lattice Fourier 和 finite workflow 聚焦测试，
  再运行 `git diff --check`。
- 迁移采用“新目标完整生成并验证后再清理旧空目录”；若发现未知文件则停止并保留源，
  不做覆盖式删除。

### 执行完成后记录

当前已完成指定现有结果目录的人工整理，并确定 `11~13` 采用平铺结构；运行程序的
输出重构仍暂缓，待单独实施：

- 已保留 shift 名称 `shift0.000`、`shift0.050`、`shift0.070`；
- 已保留 `01_results/`，并将其中的 `mode00`~`mode03` 改为不补零的
  `mode0`~`mode3`；
- 已确认每个 mode 根目录只包含子目录；
- 已将 `objective.json` 保留在各 shift 的 `99_config/`；
- 已同步各 shift 的 `config.json` 和 series 根目录的 `run_summary.json`，其中
  mode 输出路径使用 `01_results/mode{mode_idx}`；
- 已逐一校验哈希，并从各 mode 的 `13_lattice_fourier_Hz/` 删除重复的
  `bulk_cyclic_index.png` 和 `finite_k_grid_first_bz.png` 共 24 份；每个 shift
  的 `10_overview/` 继续保留唯一的几何/采样网格概览副本；
- 曾按 `data/`、`metrics/`、`figures/` 试排 12 个 mode 的 `11~13` 内容；经查看
  后确认嵌套过多，已撤销该层分类并将 180 个文件平铺回对应分析目录；
  `summary.json` 仍统一使用更明确的名称 `farfield_summary.json`，其中原先指向
  staging 目录的源数据和绘图路径已更新到最终目录；无实际结果的空 `12/13`
  已删除；
- 各 shift 的 `config.json` 与 series 的 `run_summary.json` 已标记为
  `mode_centric_v6`，并完成所有 JSON、文件路径、目录类别和残留旧文件检查；
- 已完成目录、JSON 路径和旧 `mode0X` 残留检查。

尚未修改运行程序，也尚未调整 `11_simulation_exports`、`12_farfield_FFT`、
`13_lattice_fourier_Hz` 内部的进一步分类。
