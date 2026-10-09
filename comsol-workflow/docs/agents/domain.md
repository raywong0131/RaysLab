# Domain Docs

本仓库采用 **single-context**。领域是 Windows 本机 COMSOL 6.3 光子晶体仿真与结果分析；
workspace 下的独立兄弟项目不构成本仓库的多上下文 monorepo。

## 按任务读取

- 先读仓库根 [CONTEXT.md](../../CONTEXT.md) 中与任务相关的术语。
- 科学命名和不变量以 [CONTEXT.md](../../CONTEXT.md#language)
  为依据；不把用户模式名直接当作内部数据字段，不从历史计算推断新任务目标。
- 参数任务读 [参数合同](../parameters.md)；几何、求解及输出读
  [仿真合同](../agent-simulation.md)；Unit-cell 2D 读
  [2D 合同](../agent-unit-cell-2d.md)；绘图读 [绘图合同](../agent-plotting.md)。
  只读取当前任务相关材料，不要求每次遍历全部文档。
- `docs/adr/` 如存在，只读相关决策；缺失时直接继续。首次解决确有长期取舍的决策时，
  才创建 `docs/adr/NNNN_<topic>.md`。不为本次安装制造 ADR 或回填假定的历史决策。

## 维护边界

- `CONTEXT.md` 维护项目结构、算法/脚本职责及科学语义；执行规则由 AGENTS 维护；当前参数和求解设置以实际配置、入口及结果 metadata 为准。
- 沿用 `docs/spec/plan-execute_YYYYMMDD_<topic>.md` 记录项目规则要求的方案、执行和验证。
  ADR 记录长期决策，执行记录与票据互相链接；不搬迁或批量重写历史文档。
- 输出中的领域名沿用词汇表；遇到词汇表、Memory、合同或相关 ADR 冲突时明确指出，
  不擅自改变科学定义、模式映射或既有结果。
- 改变物理假设、计算规模或运行 COMSOL 时，继续遵守 [AGENTS.md](../../AGENTS.md)。
