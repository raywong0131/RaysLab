# Project Memory

Last reviewed: 2026-10-09

记录长期偏好和科学命名，按任务读取，不保存运行状态或重复手册。当前明确要求优先于旧偏好；
精确实现、输出名和数值从当前代码、测试与结果 metadata 核对。

## Collaboration style

- 默认中文，先说结果，准确简洁；技术细节按问题需要展开。
- “确认”批准紧邻的明确方案，不重复询问。明确要求修改/实现时完成必要代码和验证，新文件本身
  不另设审批；改变科学假设、计算规模或既有结果时按具体后果判断授权。
- “按照这个配置帮我设置参数”授权最小必要配套修改，使目标流程可运行；不扩建配置系统。
  没有目标 checkpoint 时，可在用户要求的全新运行中调整 strict resume，不清空已完成 cache。
- 问设置是否正确时检查实际代码、配置或产物；引号内文件名、目录名和标识符按原拼写处理。
- 用户要求先看样例/布局再改生成器时按该顺序；明确授权两个阶段时连续完成。
- 简单任务不额外写计划；大改动的方案和执行记录要求见 [AGENTS](AGENTS.md)。

## Scientific vocabulary and invariant mappings

科学语义和不变量统一维护在 [CONTEXT.md](CONTEXT.md#language)，包含 cavity/cladding、
用户/内部 px/py、p1/p2、compact/Fourier 映射和 cavity-only 参考适用域。
本节保留旧链接入口，不再维护第二份科学定义。

## Development workflow preferences

- 新需求先做可验证的试验：算法在 tmp，调用脚本在 tests；成熟后按必要性晋升，
  不为单次研究强行扩大核心包。细则与命名在 [AGENTS](AGENTS.md#新算法与新脚本的两条晋升通道)。
- tmp 中的试验源码具有保留价值，不能作为普通缓存全量删除；archive 是保留记录，
  不参与自动开发/计算，不作为新输出位置。真实结果保护始终适用。
- 保留功能、科学合同、来源和复现性优先于形式上的少文件/少行数；复用同一算法所有者。

## Plotting, output, and naming principles

- 默认克制、清晰的 Nature 科研风格；几何坐标等比例，主次刻度朝内，colorbar 与线框对齐，
  标签准确且注明单位，比较面板采用科学上合理的共同 scale。
- 尊重用户指定 palette、标题和布局；某次批准的数值范围、字号/尺寸、箭头长度不自动推广。
  样式统一不得改变原始场、phase、单位或科学归一化。
- Hz 有符号场、用户 px/py 能带颜色、预览授权及实际 PNG 检查统一见
  [绘图合同](docs/agent-plotting.md)，本文件不再复制参数表。
- 工程 Spec 放 docs/spec，计算报告放结果 12_reports；源文件和输出命名/位置见
  [AGENTS](AGENTS.md#新文件命名位置和绘图约束)。既有科学文件名与精确输出合同保持兼容。

## Source of truth and maintenance

- 共享配置：`scripts/parameter.json`；任务 identity：该次 `99_config/config.json`；
  行为/名称：源码、归档及活动验证、`pyproject.toml`；状态：进程、日志、summary、Git diff。
- 只在任务涉及时读取相应节，不全量加入每个任务。按主题合并/替换旧规则，不累积矛盾例外。
- 用户明确改变长期偏好时同步维护；不记录当前 mesh、频率、扫描值、PID、commit、分支或一次性布局。
  专用技术规则放入可路由的参考文档，保持 UTF-8 和有效链接。
