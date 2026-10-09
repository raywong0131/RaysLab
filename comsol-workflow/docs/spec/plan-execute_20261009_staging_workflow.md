# 算法与脚本分层试验、晋升及目录归档

日期：2026-10-09。状态：已完成文档、归档及离线验证。

## 授权与目标

用户要求根据当前项目更新 AGENTS.md、CONTEXT.md、MEMORY.md、README.md，明确后续
命名、文件生成位置、绘图风格，并将 tmp 与 tests 当前内容分别放入各自 archive。
统一使用实际目录名 `tests/`，不另建 `test/` 或 `scipts/`。

两条通道：可复用新算法先在 `tmp/<topic>/` 验证，再按需晋升到 `comsol_workflow/`；
调用算法完成计算、分析、绘图的新脚本先在 `tests/<topic>/` 打通，再按需晋升到
`scripts/run_main/`、`scripts/run_sweep/` 或 `scripts/analysis/`。
普通既有实现的小修在原位置完成；新增算法/新脚本不得绕过试验阶段。

## 实施边界

- 不启动 COMSOL、不求解/重绘真实结果、不修改共享参数、核心科学计算、既有科学结论，
  不 commit/push；保留已有未提交修改。不操作兄弟项目和根 `.archives/`。
- 归档是保留文件的同盘移动，移动前检查绝对路径、目标冲突、运行进程、文件锁与链接。
  保留旧名称和子目录；不改 MPH 内绝对路径、准备清单、源哈希及历史元数据。
- 旧 tmp 任务归档后作为历史记录，含旧绝对路径的准备清单不自动续跑、不默默修复身份。
- tests 当前含正式回归测试及产物，必须保留其离线验证能力：必要时适配项目根/输出路径，
  仅改变位置关系、不改断言或算法。归档前保存原件及哈希，记录适配范围。
- 默认 pytest 不递归 archive，不自动执行候选计算脚本；离线检查显式读取归档回归，
  保留真实 COMSOL 示例排除和启动保护，避免“归档后测试数为零但看似通过”。

## 文档和目录职责

AGENTS 维护必须执行的流程、位置、命名、晋升门槛和保护规则；CONTEXT 维护架构及科学语义；
MEMORY 维护长期偏好并链接规则所有者；README 提供使用和开发操作步骤。
绘图细节由现有 docs/agent-plotting.md 维护，纠正工程 Spec 与产出报告的位置混用。

源码/主题使用能力名称和 lower_snake_case；脚本动词区分 prepare/run/export/analyze/plot/validate。
工程计划为 docs/spec/plan-execute_YYYYMMDD_<topic>.md；长期决定按需进 docs/adr/；
任务报告进入结果 12_reports。新数值数据/PNG/PDF/日志/配置使用共享输出路径职责分类，
保持既有精确合同及历史文件名。

tmp 与 tests 的试验源码不作为可随意清空的缓存；临时运行产物进入各任务 .work。
晋升前验证身份、数值/科学合同、输出和实际 PNG，晋升时复用已有所有者，更新入口和文档，
移除重复实现、保留必要兼容，完成同等检查并回填 Spec。晋升不附带新的计算授权。

## 验证与恢复

归档前后比较文件清单、大小、哈希、链接和总字节；检查归档外正式结果与共享参数未变。
检查全部受影响 Markdown 链接、现行指令一致性、pytest 收集范围、原离线回归和 Python 编译，
确认原 674 项离线检查继续可运行。
必要迁移代码修改运行聚焦检查；保留归档前文件清单和四份文档、配置/检查入口原件。
回退只还原本轮文件/目录，不使用整仓 reset/checkout。

## 执行回执

- 文档分工已落地：AGENTS 维护强制工作流、两条晋升通道、命名/位置和保护边界；
  CONTEXT 维护架构和科学定义；MEMORY 保留长期偏好并引用规则所有者；README 给出操作步骤。
  同步 docs/agent-plotting.md、docs/agents/domain.md，并新增 tmp/README.md、tests/README.md。
  科学定义从原 MEMORY 原样迁入 CONTEXT，旧 MEMORY 锚点保留。
- 归档全部原顶层条目并保留层级：tmp 的 33 个条目、173 个文件、838,760,865 字节进入
  tmp/archive；tests 的 68 个条目、140 个文件、5,139,700 字节进入 tests/archive，
  其中含 66 个测试模块。两个 archive 分别生成 archive_manifest_20261009.json。
  移动完成、路径适配之前，逐文件哈希、大小及修改时间与基线一致；没有目录链接。
- 为保留回归能力，17 个归档测试文件仅适配项目根、验证输出或移除冗余 sys.path 设置。
  适配前原件、前后 SHA256 和完整 diff 已保存；各文件所有 ast.Assert 比对不变。
  其余归档文件保持原哈希、大小和修改时间；旧 tmp 准备清单、绝对路径和源哈希未改。
- pytest 默认排除 archive、.work 和缓存，仅收集 test_*.py；候选计算脚本需显式调用。
  scripts/check_offline.py 显式读取旧回归和活动检查，保留核心/科学分进程、真实 COMSOL
  示例排除，并拒绝空测试组。启动保护统一至 tests/conftest.py；验证输出由 tests/_paths.py
  定向至独立 tests/.work 任务目录。tmp/tests 包初始化保持轻量，未新增依赖。
- 验证通过：默认 pytest 收集 1 个活动检查；聚焦检查 27 项；全量离线检查 675 项
  （核心 587、科学 88），即原 674 项加 1 个收集/归档规则检查。
  3 条 PytestReturnNotNoneWarning 来自既有 strip quotient 测试，不涉及本轮断言变化。
  受影响文档 65 个本地链接/锚点检查通过，231 个 Python 文件内存编译通过，
  git diff --check 通过。
- 核心算法和正式业务入口源码与本轮开始时一致；生产源码仅修改既有离线检查入口。
  scripts/parameter.json 和 uv.lock 的 SHA256 未变；63,892 个正式结果文件的路径清单、
  大小、修改时间及链接目标与基线一致。未启动 COMSOL、未重绘真实结果、未 commit/push。
- 本轮基线、原件、迁移 diff 和验证日志保存在
  I:\codex-task-temp\staging-workflow-20261009；聚焦验证产物移至该目录的 test-artifacts/focus，
  全量 runner 的任务目录已自行清理，试验区不保留本轮运行产物。

恢复时先确认无相关运行进程、源/目标路径均在本仓库并且无冲突，再按 baseline.json 清单
把两处 archive 中的原条目移回原位。17 个适配文件使用 before/tests 中的原件还原，
本轮修改的文档、pyproject.toml 和离线入口使用 before 中的对应原件还原；
docs/agents/domain.md 的原件位于 before/domain.md。本轮新增文件可移入上述任务备份目录保留。
只处理本轮清单内文件，不使用整仓 reset/checkout，不触及原有未提交工作和正式结果。
