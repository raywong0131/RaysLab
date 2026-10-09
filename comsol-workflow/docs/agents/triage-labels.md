# Triage Labels

采用上游默认的五个角色名；标签描述任务准备状态，不提供额外执行权限。

| Skill 角色 | 本项目标签 / Status | 含义 |
| --- | --- | --- |
| `needs-triage` | `needs-triage` | 待核实范围、证据和影响 |
| `needs-info` | `needs-info` | 缺少必要复现信息、科学定义或计算授权 |
| `ready-for-agent` | `ready-for-agent` | 范围、验收与授权明确，可由 agent 实施 |
| `ready-for-human` | `ready-for-human` | 需要人的科学判断、外部访问或手工操作 |
| `wontfix` | `wontfix` | 明确不处理，记录理由 |

同一票据只保留一个当前 triage 角色。一般从 `needs-triage` 分流；补齐信息后回到
`needs-triage` 复核。完成状态由 tracker 的关闭/完成字段表达，不增加第六个 triage 角色。

仿真票据应写明用户/内部模式、参数来源、几何/mesh、求解规模、输出与恢复策略、验收证据。
`ready-for-agent` 不自动授权昂贵求解、扩大扫描、覆盖结果、commit 或 push。
验证不能以启动 `SimulationRun`、mesh builder 或 eigensolve 作为普通 smoke test。
已有数据能够回答问题时使用对应分析入口；细节见 [AGENTS.md](../../AGENTS.md)。
