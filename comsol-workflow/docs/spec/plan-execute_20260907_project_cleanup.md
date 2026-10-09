# COMSOL workflow 清理与保功能简化

日期：2026-09-07。状态：已完成本次保功能清理与安全简化。

## 范围与保护边界

用户授权审查、删除不必要临时文件、优化结构和简化框架，要求保留所有功能与计算结果。
本次覆盖安装、入口编排、几何与数值模块、COMSOL 生命周期、结果兼容、分析、测试和文档。
沿用 2026-09-02 审计的职责划分，重新核对当前入口、符号引用、重复函数及文件用途。
不执行该旧审计中退役功能的建议：直接调用、诊断、历史结果读取和迁移能力继续保留。

当前四扇区 finite-quarter 计算仍在运行（Python PID 19856/40496，COMSOL server PID 33352）。
不改其求解、几何、配置、恢复或结果发布代码，不停止进程，不启动新的 COMSOL。
不调整依赖版本或整体同步正在使用的环境；仅修复已迁移 blueprint 开发依赖的定位。
用户已有工作树改动及上一任务迁出两个项目产生的状态保留，不提交或推送。

## 本次切割与验证

1. 删除已核实的空 `.codex*`/`.test-tmp` 目录、`.test_bp_final` 合成测试产物、
   `.pytest_cache`、源码树 `__pycache__` 和仅打印 ready 的 `.codex_process_probe.py`。
   不进入 `.git`、`.venv`、`.understand-anything`、`.archives`、正式输出或参考结果进行清理。
   保留安装元数据、编辑器配置、计算缓存、staging 和失败恢复日志。
2. 将开发依赖 `microfab-blueprints` 从 `blueprints` 改为 `../blueprints`，同步锁文件，
   保留跨项目几何一致性测试；验证锁文件与导入，再执行该测试。
3. strip 新旧布局分析入口复用现有模式筛选，并共享同一 CSV 行构建实现。
   保留函数签名、严格阈值、错误消息、结果列、历史路径回退和绘图行为；复用现有测试。
4. 文档以根 README 的工作流和结果布局为入口，修复 scripts README 的三/四目录冲突，
   更新迁出项目链接，保留历史决策记录。

基线：strip 分析两份测试共 23 passed。已使用 Windows 长路径盘点
55,024 个结果文件、146,180,091,419 字节。运行中的结果允许由原计算进程继续更新；
完工后检查非活动结果的路径、大小、mtime 及源配置哈希，单独列出活动目录变化。
备份、结果清单与验证日志存于系统临时目录 `comsol-simplify-20260907`。

## 保留项与盲点

- full/quarter、strip 诊断和历史调用接口：有不同边界条件或直接消费者，不能以重复为由退役。
- staging/finalize、模型复用、缓存身份与原子写入：承担数据保护，不简化其校验。
- 历史 layout 迁移和外部导入：外部数据/消费者无法穷尽，继续保留。
- 大型模块：本次检查入口、依赖和重复函数，不宣称对每行科学算法重新认证。
- 科学求解与实际渲染：通过已有离线/模拟测试验证；不以新 COMSOL 计算作为清理测试。

## 执行回执

### 实际变更

- 清理 15 个明确路径，共 185 个文件、7,197,922 字节（约 6.86 MiB）。
  其中 `.test_bp_final` 为 `small_cli_contract` 合成测试数据，非加工正式产物；
  `tests/__pycache__` 中的 `.pyc.<PID>` 是 9 月 2 日 pytest 原子写入残留，已一并清理。
  其余为 Python/pytest 缓存、空临时目录和只打印 `ready` 的诊断探针。
- `pyproject.toml` 与 `uv.lock` 改为 `../blueprints`；锁文件仅两处路径变化。
  离线重装同版本 editable blueprint 开发包，实际导入指向兄弟项目 `workflow/`。
  未升级运行依赖，未整体 `uv sync`，未修改兄弟项目源码。
- `srip1d_trend.py` 的两条历史/当前布局路径共享 `_trend_row()`，只维护一份导出行 schema；
  频率解析直接导入已有实现，保留原模块可调用名称。模式筛选、读取异常处理与行转换的
  异常边界均不变。源码增加 10 行，但消除了两份 CSV schema 和一个纯转发函数的维护义务。
- `scripts/README.md` 从 167 行收敛到 59 行，保留脚本目录和独有选项，通用合同指向
  根 README/仿真合同，消除“三个目录/禁止第四个”与实际四个工作流目录的矛盾。
  根 README 补充兄弟项目关系；AGENTS 仅修正迁出项目的两个路径。

### 审查覆盖与保留结论

| 范围 | 当前证据与处理 |
| --- | --- |
| 安装与入口 | 核对 8 个注册命令、直接脚本、开发依赖及包发现；修复失效的本地依赖，保留全部入口 |
| 参数与科学模型 | 核对参数所有者、源文件哈希及入口依赖；共享 JSON、Fourier 映射、几何、mesh、模式数未改 |
| band/2D/扫描 | 核对入口目录、导入、公共函数与旧审计；模型复用、valley 求解、扫描入口继续保留 |
| finite/quarter/strip | 追踪直接调用和重复函数；finite 与 strip 的相同函数仍依赖不同模块参数，未强行合并 |
| 生命周期与结果 | 模型、恢复数据、缓存身份、staging、原子最终化均为有效保护边界；未改代码或数据 |
| 离线分析 | 完整读取并合并 strip 行 schema；有限腔、Fourier、趋势及历史路径兼容由现有测试覆盖 |
| 临时与持久化文件 | 检查根目录和源码树；正式输出/参考结果只做长路径 metadata 盘点，其他用户工具和环境目录保留 |
| 文档与历史 | 收敛当前重复说明，保留历史设计理由、迁移工具和直接调用能力 |
| 不可穷尽部分 | 未盘点仓库外任意消费者、外部历史结果，也未重新执行科学求解；不能据此退役其兼容接口 |

### 验证

- 修改前 strip 分析基线：23 passed。
- 修改后现有离线测试：543 passed，3 个既有 `PytestReturnNotNoneWarning`，85.52 秒。
  测试运行显式拒绝 `mph.start()`，并排除三个需显式运行的 COMSOL 示例模块：
  `test_simulation_spatial.py`、`test_simulation_utils_spatial_comsol.py`、`test_simulation_plotting.py`。
  包含跨项目 blueprint 几何一致性、当前四扇区 mocked 流程及全部现有分析单元测试。
- `uv lock --check --offline`、源码内存编译、原公开函数名称保留检查、文档本地路径检查、
  全工作树 `git diff --check` 通过。
- 对照修改前源码哈希，科学源码与共享 JSON 未变；唯一变化的 Python 文件是
  `scripts/analysis/srip1d_trend.py`。用户先前源码改动未被覆盖。
- 结果文件数保持 55,024；无缺失、无新增。唯一变化是活动计算的
  `00_model/comsol_progress.log` 增长 531 字节；其他结果路径、大小、mtime 全部不变。
  该核查是文件 metadata 对比，不宣称对 146 GB 数据做了逐字节哈希。
- 最终确认 Python PID 19856/40496、COMSOL server PID 33352 仍在运行。
  本次没有新启动 COMSOL，也没有执行正式数据重绘、迁移或结果删除。

### 回退与证据

恢复目录：`C:/Users/Administrator/AppData/Local/Temp/comsol-simplify-20260907/`。
其中保留 `before/` 源码备份、`diff/` 本次差异、`cleanup-receipt.json`、结果前后清单、
源码初始哈希和验证回执；已清除本轮 pytest 与 Matplotlib 临时产物。
若需撤销源码修改，应逐项反向应用这里的本次 diff，保留任务开始前的用户改动；
不应使用整仓 `git reset/checkout`。临时缓存无需恢复，必要时由工具自行重建。
没有提交、推送或修改 Git 配置。
