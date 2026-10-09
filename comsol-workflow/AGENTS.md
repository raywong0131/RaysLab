# AGENTS.md

本仓库是 Windows 本机 COMSOL 6.3 控制器；默认用中文沟通。按当前任务读取专用规则，
不要求每次编辑都阅读全部文档、Memory 或整个代码库。

## 工作方式与授权

- 开始修改前检查 `git status --short --branch`，保留用户已有修改；读取相关代码及调用关系。
- “帮我改/设置/修复/实现”授权完成该范围内的实现、必要文件、聚焦验证及修复，不因新增
  `.py` 文件或已有审批而重复询问。普通实现选择自行判断；影响物理含义、计算规模或交付范围的
  未决选择才需澄清，先完成不依赖该答案的工作。
- 完成包含实现、适当验证、清理本轮临时产物和报告结果。仅写方案或完成第一版不算完成；
  遇到真实阻塞说明已完成部分及最小缺失条件。
- 跨模块重构、新算法/新计算脚本、目录迁移或配置/数据流大改动，在 `docs/spec/` 先记录中文
  `plan-execute_YYYYMMDD_<topic>.md`，完工后回填执行与验证。小修、措辞和参数修改无需独立计划。
- 技能服务当前任务，不能用简化偏好缩减明确请求。按需读取正文和支持材料。
  不因技能有多代理示例而自动委派，遵守当前会话设置和明确授权。
- 已批准具体方案后继续执行；确需暂停时引用造成暂停的文件、原文及适用原因。

## 入口与按需读取

| 当前任务 | 入口或参考 |
| --- | --- |
| Unit-cell band / Q-k | `scripts/run_main/run_band_pair.py`、`scripts/run_main/run_band_solo.py`；[仿真合同](docs/agent-simulation.md) 的 band 部分 |
| Unit-cell 2D、ζ 扫描、valley refinement、A–D 重绘 | [Unit-cell 2D](docs/agent-unit-cell-2d.md)；`scripts/run_main/run_unit_cell_2d.py` |
| Strip、finite、quarter、几何/mesh/输出、显式孔洞 | [仿真合同](docs/agent-simulation.md)；`scripts/run_main/run_strip_1d.py`、`scripts/run_main/run_finite.py`、`scripts/run_main/run_finite_quarter.py` |
| 模式命名、物理解释、compact 参数映射 | [Context 科学语义](CONTEXT.md#language) |
| 科研绘图 | [绘图合同](docs/agent-plotting.md)；长期偏好见 [Memory](MEMORY.md#plotting-output-and-naming-principles)；2D A–D 另读 2D 文档 |
| 用户要求 commit / push | [本机 Git 流程](docs/agent-git.md) |
| 独立版图或后处理子项目 | 对应兄弟项目 `../blueprints/AGENTS.md` 或 `../comsol-data-analysis/AGENTS.md`，仅进入该项目时适用 |

专用细节只在对应任务适用。准确命令以 `pyproject.toml` 和当前入口为准，参数以当前参数文件
及该次输出配置为准；历史计划不是当前运行配置。

## 项目共同约束

- 使用本机 COMSOL，不引入 Slurm、SSH、远程 scratch、服务器同步或 Linux 后台 launcher，除非用户改变定位。
- 成熟可复用几何/求解/分析位于 `comsol_workflow/`；成熟主入口在 `scripts/run_main/`，
  COMSOL 扫描在 `scripts/run_sweep/`，无求解分析在 `scripts/analysis/`。新增能力先按下述两条通道试验。
- 标准 package：包内相对导入，正式脚本用 `comsol_workflow.*` 或 `scripts.*`；
  试验入口按通道使用 `tmp.*` / `tests.*`，正式代码不反向依赖它们。
  不恢复 flat import，不向 `sys.path` 加包子目录。包初始化保持轻量，不导入 `mph` 或启动 COMSOL。
- 共享参数源 `scripts/parameter.json` 保持无注释的标准 JSON。Unit-cell 2D 用独立
  `unit_cell_2d` 对象及根级 mesh，可用 `COMSOL_WORKFLOW_PARAMETER_PATH` 快照，详见 2D 合同。
  不互借不同入口的几何/模式数；扫描、重绘和底层函数不得隐式改写共享参数。
- 保存的 `*_simulation_config.json` 用
  `comsol_workflow.simulation_config_geometry.load_simulation_config_geometry()` 读取；只取已验证的
  footprint/holes，不用历史求解字段替代当前参数，不猜测缺失 cell 拓扑。
- 活动配置测试从已加载的 `SharedParameters` 派生预期值，硬编码仅用于明确的合成/历史 fixture。
- 保持 `scripts/parameter.json` 的声明数据与 `scripts/run_main/parameter_config.py` 的校验/更新逻辑分离。

## 新算法与新脚本的两条晋升通道

统一使用 `tests/`，不另建 `test/`。先拆清算法和编排职责，再建立一个能力主题 `<topic>`：

| 类型 | 首次开发位置 | 允许依赖 | 验证后归属 |
| --- | --- | --- | --- |
| 可复用数值/几何/分析算法 | `tmp/<topic>/<algorithm>.py` | 现有 `comsol_workflow` 和已安装库 | 必要时并入现有 `comsol_workflow` 所有者，确有独立职责再建模块 |
| 调用算法完成计算/导出/分析/画图的入口 | `tests/<topic>/<verb>_<capability>.py` | 核心包及该主题已明确标注的试验算法 | `scripts/run_main`、`scripts/run_sweep` 或 `scripts/analysis` |
| 纯离线验证 | `tests/<topic>/checks/test_<capability>.py` | 被验证实现、合成数据或授权只读 fixture | 留在 `tests`，随晋升改为验证正式所有者 |

- 新需求先写 Spec，再在上述目录完成实际可运行实现与适当验证，不直接把新算法放入核心包，
  不直接把新计算脚本放入正式 scripts。既有错误修复、参数/路径适配和文档维护在原位置完成，
  新增数值方法不能作为“小修改”跳过试验。单纯复用已有算法的任务只新增候选入口，不创建空算法模块。
- 混合脚本把可复用算法抽到 `tmp/<topic>`，在 `tests/<topic>` 保留调用、参数、日志和输出编排；
  复用既有函数，不复制成熟算法。临时算法禁止从正式核心包或正式脚本反向导入。
- 从仓库根用 `uv run python -B -m tmp.<topic>.<algorithm>` 或
  `uv run python -B -m tests.<topic>.<verb>_<capability>` 调用；使用包导入，不向 sys.path 添加子目录。
  这些命令是路径形式，不授权实际 COMSOL 计算；导入不得启动求解。
- `tmp` 是算法试验源码，不是可全量清空的临时垃圾；`tests` 同时承载候选入口和纯验证，
  只有 `test_*.py` 作为 pytest 检查，`run_`/`export_`/`analyze_`/`plot_` 候选必须显式调用。
- 合成输出、验证缓存及 staging 分别进入 `tmp/.work/<task_id>/` 或 `tests/.work/<task_id>/`，
  不写源码同目录或 archive。真实结果仍进入已授权结果根，不因脚本处于试验区改变结果保护规则。
- 晋升门槛：本次需求已实现；来源/配置身份、单位、相位、模式映射和适用域明确；必要数值验证及
  输出检查通过；绘图查看实际 PNG；COMSOL 工作遵守实际计算边界，不以离线通过冒充科学实测。
- 达到门槛且需求适合长期维护时，在本次明确范围内自行晋升，无须为移动文件再问批准。
  算法优先合入已有所有者；脚本按求解/扫描/只读分析归类。同步消费者、验证、README/目录和必要 CLI，
  复用同一实现，试验区只保留有真实调用者的兼容入口或归档记录，禁止两套业务实现继续分叉。
- 仅用于本次研究的能力可留在试验区；Spec 记录位置、已验证范围、未决项和是否晋升。
  晋升不自动授权更多计算、历史结果迁移或发布。

`tmp/archive/` 与 `tests/archive/` 保存本轮之前的内容。归档保留原名称和目录，不默认递归执行，
不作为新输出目录、不当作待删除缓存。旧 tests 回归由 `scripts/check_offline.py` 显式读取；
可做必要路径维护，但不削弱断言。旧 tmp 准备清单含历史绝对路径时只作记录，不能静默改写哈希续跑。
以后仅归档明确完成且在授权范围内的主题；不得每次清空两个工作区。

## 新文件命名、位置和绘图约束

- 主题、模块和函数使用能力名称及 `lower_snake_case`，例如 `boundary_sampling`、`dipolar_radiation`；
  不用 `new`、`final`、`v2` 或 S4/S5 等任务编号作为新共享能力名称。科学批次/历史标识及旧兼容名保留。
- 脚本前缀表达真实职责：`prepare_` 准备，`run_` 求解/计算，`export_` 导出，`analyze_` 分析，
  `plot_` 绘图，`validate_` 校验。含求解的脚本不能仅标为 export。验证文件只用 `test_*.py`。
- 根目录仅维护约定的项目入口文档/配置；新工程 Spec 为 `docs/spec/plan-execute_YYYYMMDD_<topic>.md`，
  完工回填；领域指南为 `docs/<topic>.md`，长期决定按需为 `docs/adr/NNNN_<topic>.md`。
  任务计算报告为 `<result_root>/12_reports/[batch/]<topic>_report.md`，不把工程 Spec 写入结果 reports。
  `AGENTS.md`、`CONTEXT.md`、`MEMORY.md`、`README.md`、`__init__.py` 是保留名例外。
- 新数据/图片/日志/配置调用 `output_paths.artifact_path`、`result_io` 和已有图保存所有者；
  新绘图保存复用 `figure_output.save_figure_formats`，已有 field 图继续复用 `field_plotting`。
  不再各脚本自建一套分类/保存器。JSON 必须明确 data/log/config 职责。
- 通用新输出按职责放入 `00_model`、`01_results`、`10_overview`（PNG）、`11_pdf`、`12_reports`、
  `13_svg`、`80_logs`、`90_history`、`99_config`；根 HTML 索引是允许的导航文件。
  文件为 `<capability>[_<mode>][_<quantity>].<ext>`，同一图各格式共用 stem，用户模式名用 px/py。
  批次位于格式目录内部，主 PNG 平铺；当前 run_id 及参数派生 case 名依现有入口身份合同。
  配置/来源记录进 99_config，科学数值 CSV/NPZ/JSON 进 01_results，运行状态和检查进 80_logs。
- 已有精确输出合同优先（例如 finite v6 总览 CSV、per-mode 11/12/13 子目录、S4 compact 布局），
  不因通用规则移动旧产物或改写保存字段。变更这些合同另做显式兼容迁移。
- 绘图统一执行 [绘图合同](docs/agent-plotting.md)：克制 Nature 风格、无衬线字体、主次刻度朝内、
  等比例几何坐标、对齐的 colorbar、可靠共同 scale、数据来源及单位清楚，采用既定模式/Hz 色标。
  样式统一不能改变数值、原始场、显示 phase 或科学归一化；用户给定具体图合同优先。

## COMSOL 与结果边界

- 实际计算须在用户请求范围内。启动前从实际配置核对并报告用户/内部模式名、几何、mesh、
  频率/k 路径或扫描范围、求解规模、缓存恢复、输出目录和后处理；检查当前进程、PATH、license
  及目标 config identity。配置冲突时不覆盖 config 或混用 cache。
- `SimulationRun`、mesh builder、eigensolve 都可能启动 COMSOL，不作为普通 smoke test。
  已有完整数据足以重绘时使用分析入口。valley 新增点/修复点先准备计划，再按 2D 合同取得明确批准。
- 正式结果位于 `scripts/.out/` 的 `unit_cell_band/`、`unit_cell_2D/`、`strip_1d/`、
  `finite_cavity/` 四个一级目录。新计算从实际参数派生非冲突路径，不创建第五个一级目录。
- 未经明确授权，不移动、覆盖或删除既有计算结果。指定图原位重绘只覆盖该图，不扩到其他 mode 或面板。
  结果重组先确认解析后的源/目标在授权目录内并检查冲突，结束后比较文件数、总字节及主要产物。
- 长任务启动后检查 PID、stderr、progress，并报告日志及结果路径；通常监控到 mesh 完成或首个完整点后
  允许后台继续，除非用户要求持续监控。进程退出、错误解释清楚、summary 和产物满足合同后才报告完成。
- 用户要求停止时，只终止目标任务相关的 Python/COMSOL 进程并确认退出。后台 helper 用隐藏窗口。

从仓库根执行 `uv run <pyproject 中的命令>`；仅在环境缺失或依赖变动时 `uv sync`。
COMSOL PATH 包含 `C:\Program Files\COMSOL\COMSOL63\Multiphysics\bin\win64`；
找不到 COMSOL 时在同一 PowerShell 会话核对 `where.exe comsol` 和 PATH。

## 验证与收尾

- 验证与变更风险匹配：纯指令/Markdown 检查内容、链接和 diff；配置检查解析及对应合同；
  Python 行为修改运行聚焦测试和相关 `py_compile`。全仓测试用于广泛影响或明确提交检查需求。
  通过后仅在新修改、失败或未决风险出现时追加/重复测试。
- 已确认不启动 COMSOL、不读写正式结果的本地测试可直接执行，修复本次修改引起的失败。
  不为可逆措辞修改新增镜像实现的测试，不把 Windows 文件锁/权限失败误报成代码失败。
- 对本次修改运行 `git diff --check`。生成器合同测试和实际 PNG 检查见绘图文档。
- 只清理本轮临时目录、模型、日志和缓存；必要诊断/回滚资料保存在任务级临时目录并报告路径。
  不清空既有 cache，不向 ignore 加 `.test-tmp/` 或 `.slurm_out/` 隐藏应清理产物。
- 不访问或重整 `.archives/`，除非用户指定；仓库清理中的非临时删除候选按已授权范围移入归档，
  交由用户最终删除。偏好少量源文件和浅层输出目录，新 helper 应有明确职责或复用价值。
- 仅暂存当前任务文件。普通修改不自动 commit/push；明确要求“push 当前版本”时按 Git 流程完成必要提交与推送。
- 交付说明实际修改、验证及未验证部分；涉及输出时说明更新文件、目录变化及是否启动 COMSOL。

## Agent skills

Matt Pocock 的指定 11 个 skill 位于 workspace 根的 `../.agents/skills/`，由
`../skills-lock.json` 记录来源。调用时以本仓库为工作目录；第三方正文保持原样，
项目适配集中在下列文件。只读与当前任务相关的配置。

### Issue tracker

规格与票据使用 **GitHub Issues：`jia-yli/comsol-workflow`**；外部 PR 不作为 triage 请求入口。
`to-spec`、`to-tickets`、`triage`、`implement`、`code-review` 涉及票据时读取
[docs/agents/issue-tracker.md](docs/agents/issue-tracker.md)，其中规定认证、命令和依赖关系。

### Triage labels

使用 `needs-triage`、`needs-info`、`ready-for-agent`、`ready-for-human`、`wontfix` 五个默认标签。
映射和 COMSOL 任务准备要求见 [docs/agents/triage-labels.md](docs/agents/triage-labels.md)。
标签不能扩大计算授权或绕过本文件的结果保护规则。

### Domain docs

采用 **single-context**：根 [CONTEXT.md](CONTEXT.md) 维护项目结构与科学语义，`docs/adr/` 按需记录长期决策。
现有 Memory、领域合同与 `docs/spec/` 执行记录保留；按需读取规则见
[docs/agents/domain.md](docs/agents/domain.md)。`grill-with-docs`、`domain-modeling` 及其他工程
skill 应沿用这些科学定义；`tdd` / `diagnosing-bugs` 仍遵守无求解验证和实际计算边界。

## 已批准专题结果根

- 边界分析汇总保留 `results/S4_boundary_analysis/`，dipolar 分析保留
  `results/S5_dipolarSingularity_analysis/`；它们是通用 scripts/.out 四类工作流之外已批准的专题根。
- S5 顶层按格式分类，主 PNG 平铺；模型、数据、PDF、报告、日志、配置内部共享批次名
  （历史 `02_gamma_cladding` 等不改名）。整批运行兼容入口仅置于 90_history。
- 新程序不得重新创建旧 `results/XX_dipolarSingularity_analysis` 别名，不批量改写旧元数据来源。
  是否存在活动计算或 junction 必须当场检查，不从历史 PID/一次性清理记录推断状态。
- 运行中的模型、插值源及日志保持物理位置；移动前确认进程退出、绝对路径依赖与结果完整性。
