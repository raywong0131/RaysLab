# GPT-6 Astra 持久指令优化

日期：2026-09-07。状态：已完成修改与静态验证。

## 依据与目标

已实际读取：

- Eric Provencher：[Rethinking skills and prompts for GPT-6 Astra](https://x.com/pvncher/status/2095991462416490862)。
- OpenAI：[Using GPT-6 Astra / Prompting best practices](https://developers.openai.com/api/docs/guides/latest-model)。
- 官方 [AGENTS.md 加载规则](https://developers.openai.com/codex/guides/agents-md) 与 [Skills](https://developers.openai.com/codex/skills)。
- 本机 `openai-docs/SKILL.md`、`skill-creator/SKILL.md` 及相关支持文档。

采用短而明确的触发描述、按需读取、与风险相称的验证和明确完成条件。
保留科学约束、用户授权范围和正式结果保护，不把官方示例原样堆进项目指令。

## 实际来源与范围

- 根 `AGENTS.md` 已出现在本轮上下文；磁盘版本已核对。
- 用户级 `C:/Users/Administrator/.codex/AGENTS.md` 为空；未发现 override。
- `MEMORY.md` 是项目引用的按需文档，不因文件名自动加载；用户级 memories 目录为空。
- 当前配置为 `gpt-6-astra`、`xhigh`、`priority`、HTTP Responses；保留这些设置。
- 两个个人技能 `nature-data`、`simplify-codebase` 可维护；本轮审计读取正文不等于执行其业务工作流。
- 已审计技能目录、可见描述和相关正文。系统/第三方缓存由分发系统管理，保留原件；通过个人持久约定明确任务边界。
- Understand 系列在 `.codex/skills` 与 `.agents/skills` 均为指向同一 checkout 的 junction；不视为两份可独立改写的源码。
- `handoff` 存在于磁盘但未出现在本轮可用技能列表；不声称它本轮已加载。
- `blueprints/AGENTS.md` 和未跟踪的 `comsol-data-analysis/AGENTS.md` 仅属于各自子项目；只审计，不改写。
- `.codex/rules/default.rules` 是历史命令审批规则；保留其权限语义，不把已批准命令视为当前任务授权。
- 系统/开发者注入指令和托管权限不是可编辑的用户文件，本任务不会修改它们。

## 方案、数据流与兼容策略

1. 根 AGENTS 保留项目范围、入口索引、共同约束、授权/完成条件与按需文档路由。
2. 将 Unit-cell 2D、仿真/输出、绘图和 Git push 细节迁入 `docs/agent-*.md`；迁移时保留有效科学合同。
3. Memory 保留长期偏好与物理命名；实现细节进入对应工作流文档，消除重复和已确认冲突。
4. 精简个人技能触发及必读流程，保留支持资料和默认自动发现；不新增依赖或修改插件缓存。
5. 用户级 AGENTS 保存简短的持续执行、按需技能与比例验证约定。

读取路径变为：根 AGENTS → 当前任务需要的参考文档/Memory 节 → 当前源码、参数与结果 metadata。
这只改变指令组织；不修改 Python、参数、求解流程、正式结果和 Git 提交。

## 已核实的问题

- 根 AGENTS 的 `run_unit_cell_band.py` 已不存在；当前 pyproject 注册 `run_band_pair.py`、`run_band_solo.py`。
- Memory 同时写 finite/quarter 默认扫描 ratio 派生列表与“只运行单组 x/y”；当前配置读取器支持前者及显式单组替代。
- Memory 同时保留 finite 的 `shiftX.XXX` 旧名和 directional 名；当前 runner 使用 `shiftxX.XXX_shiftyY.YYY`。
- Memory 对任何新 `.py` 文件逐项审批，超过一般已授权实现所需；改为按范围与后果判断。
- 绘图预览审批只适用于用户指定的预览阶段；明确直接修改/批准两阶段时继续完成。
- nature-data 对简单声明措辞强制加载 manifest 与三份 core，并跑全套清单；改为按请求深度选择。
- Ponytail 的 `ANY coding task`、`ACTIVE EVERY RESPONSE`、三行输出与“交付简版再询问”可能干扰当前任务；个人约定禁止其缩减已请求交付或跨任务强制生效。

## 验证与回滚

- 修改前保存任务级原件与哈希；仓库原有脏文件另记哈希，不写入其内容。
- 校验 UTF-8、Markdown 本地链接、技能 frontmatter、支持资源和剩余旧引用。
- 用代表性任务逐项检查授权、预览、COMSOL、技能路由、验证边界；不启动 COMSOL，不做无关全仓测试。
- 运行 `git diff --check`，核对最终 diff 与既有改动哈希。
- 仓库修改可按本次 diff 回退；用户级文件以任务备份逐文件恢复，不覆盖之后的新改动。

## 执行结果

共修改/新增 7 个仓库文档和 7 个用户级指令文件；未修改 Python、共享参数、模型配置、插件启用状态或命令审批规则。

| 文件 | 原状态 → 当前状态 | 实际变化 |
| --- | --- | --- |
| [根 AGENTS](../../AGENTS.md) | 475 → 83 行；27,746 → 7,354 字节 | 改为任务路由与必要边界，修正 band 入口，取消逐文件审批，校准完成和验证条件 |
| [Memory](../../MEMORY.md) | 625 → 92 行；41,229 → 6,292 字节 | 保留长期偏好与物理命名，移出专用实现和数值，清除冲突 |
| [Unit-cell 2D](../agent-unit-cell-2d.md) | 新增 | 保留 K/N、Quarter/full、identity、A–D、扫描、真实 refinement 和完成判据 |
| [仿真合同](../agent-simulation.md) | 新增 | 保留几何、mesh、参数、模式筛选、输出布局、显式孔洞读取合同 |
| [绘图合同](../agent-plotting.md) | 新增 | 限定预览审批，保留绘图复现和科学 scale，按改动选择合同测试 |
| [Git 流程](../agent-git.md) | 新增 | 保留已授权 push 的本机提交/验证/推送流程 |
| 用户级 `.codex/AGENTS.md` | 空 → 20 行 | 持续执行、按需技能、比例验证和第三方技能适用范围 |
| `nature-data/SKILL.md` | 62 → 31 行；description 831 → 197 字符 | 单句修改可直接完成，完整包才读工作流，保留事实与标识符真实性要求 |
| nature-data `manifest.yaml`、`static/core/stance.md`、`static/core/workflow.md` | 原位更新 | 取消 always_load，按需引用，按目标期刊适用政策，取消固定四段输出 |
| `simplify-codebase/SKILL.md` | description 449 → 192 字符 | 按范围读代码/文档，基线与证明深度按风险决定 |
| simplify-codebase `references/execution-and-recovery.md` | 原位更新 | 验证环节改为适用项，单个阻塞不阻止独立工作，小改动无需完整回执模板 |

simplify-codebase 正文与支持文档略有增长，用于澄清原有停止条件；未以行数下降替代语义正确性。
系统和第三方 SKILL 原件保留。没有卸载/禁用插件，也没有批量改写缓存或改变技能隐式调用策略。

### 验证结果

- 14 个涉及文件通过 UTF-8、Markdown 链接/锚点、代码围栏与实际入口路径静态检查。
- 两个个人技能通过未修改的官方 `quick_validate.py`；manifest 引用文件全部存在，隐式调用仍开启。
- Windows 的旧 Anaconda 默认 GBK 导致首轮校验失败；改用同一解释器/PyYAML，仅在校验进程内令
  `Path.read_text` 缺省 UTF-8 后通过。曾尝试把旧 PyYAML 加入新解释器，因 `collections.Hashable`
  兼容问题放弃；没有安装依赖或修改任一 Python 环境/校验器。
- 按原件与新文档对照关键技术字面量；已确认移除旧 `run_unit_cell_band.py`，其余有效合同已迁移或合并。
- 对原有 7 个已跟踪脏文件逐一比较 SHA256，内容未变化；未操作两个既有未跟踪目录和已有方案文档。
- `git diff --check` 通过。纯文档/指令变更未运行业务 pytest、py_compile、图像重绘或 COMSOL。

### 代表性场景审查（规则核对，非新模型实跑）

| 场景 | 当前规则支持的行为 |
| --- | --- |
| 简单修复需要新增 helper/test | 当前实现授权覆盖必要文件，运行相关检查，无逐文件审批 |
| “试一下这张图” | 仅指定数据/图原位预览；共享生成器等待该阶段批准 |
| “直接修改生成器并重绘 A” | 完成生成器及 A 的合同验证，不再次请求样例审批，不覆盖 B–D |
| valley 新增采样/不收敛修复 | 先形成具体计划，仅在明确批准后启动真实 COMSOL |
| 文档措辞或数据声明翻译 | 只做必要引用/事实检查，不运行求解或完整 FAIR/仓库清单 |
| 一个删除候选证据不足 | 暂停该候选，继续其他独立且已授权的改动 |
| 普通修复触及 Ponytail/业务分析关键词 | 按任务能力和范围选择技能，不自动缩减交付或引入外部发布 |
| 模式命名或 finite shift | 查询 Memory/当前配置，保留内部字段和 directional case identity |

### 回滚与生效

原件、逐文件前后哈希、指标和可审阅 diff 保存在：

`C:/Users/Administrator/AppData/Local/Temp/codex-astra-instructions-20260907-qwvwkb8s/`

- `manifest.json` 映射原路径与备份路径，包含修改前/后的 SHA256 与新增文档路径。
- `changes.patch` 包含 9 个既有文件与 4 个迁移文档的差异；本执行报告由仓库 diff 单独审阅。
- 恢复某文件前先比对 after_sha256，确认未发生后续编辑，再用对应备份逐文件恢复；新增文档仅在
  确认无后续改动后移除。不要 reset 工作树或回滚既有代码修改。

未生成结果目录，未启动 COMSOL，未 commit/push。临时脚本通过 stdin 执行，未留在源码目录；
仅保留上述任务级回滚资料。

AGENTS 在新 run/session 构建指令链；建议新建会话以加载精简后的内容。本轮已注入的旧文本不会因磁盘修改
追溯消失。技能通常自动发现变更，未刷新时重启 Codex。未启动额外模型会话，因此本次证明的是文件和规则一致性，
不宣称已量化后续模型时延、token 节省或行为改善。
