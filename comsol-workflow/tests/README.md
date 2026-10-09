# 候选计算脚本与离线验证

本目录是 `scripts/` 的前置区。新需求的计算、导出、分析、绘图入口先写入
`<topic>/<verb>_<capability>.py`，调用核心包或该主题明确的 tmp 试验算法。
脚本负责参数、来源、任务编排、日志和输出；可复用计算实现放 tmp，不内嵌第二套算法。
验证后按职责晋升到 scripts/run_main、run_sweep 或 analysis，验证转向正式所有者。
完整流程见 [AGENTS](../AGENTS.md#新算法与新脚本的两条晋升通道)。

候选文件前缀为 prepare_/run_/export_/analyze_/plot_/validate_，从仓库根按包显式执行。
`<topic>/checks/test_<capability>.py` 只做合成或已授权只读来源验证，不能启动 COMSOL；
它们留在 tests，不晋升为计算入口。默认 pytest 仅收集 test_*.py，不递归 archive 或 .work。
conftest 同时阻止真实 mph.start，并将验证输出放入任务 .work；mock 求解测试可继续运行。

```powershell
# 当前候选区的纯检查
uv run python -B -m pytest -q
# 归档回归与活动检查，分进程隔离核心/科学环境；排除真实 COMSOL 示例
uv run python -B scripts/check_offline.py
# 指定历史离线回归
uv run python -B -m pytest -q tests/archive/test_farfield_fft.py
```

`.work/<task_id>/` 存验证输出、日志和缓存；真实计算仍使用明确授权结果根。
本次原 66 个测试模块及 .out/缓存移入 archive；归档回归不会被当作完成的新需求实现。
archive_manifest_20261009.json 记录原文件身份及必要路径适配，现有断言和科学计算保持。
既有源码可做明确的回归维护，禁止新任务写入 archive 或修改历史结果/清单以放宽验证。

后续已按用户授权删除归档字节码及四张可重建测试 PNG，66 个原模块全部保留。
P1/P2/P3 分别为 bulk/cladding 边界比较、bulk/side/corner 分区预览和孔洞裁切预览；
现有算法已在核心，候选的价值在于通用编排与可视核验，尚未改名、移动或晋升。
Fourier 人工演示有转为纯检查的价值；三份实际 COMSOL 诊断仍排除自动执行。
详见 [候选审查](../docs/staging_cleanup.md) 与 [逐项分类](../docs/spec/staging_cleanup_inventory_20261009.json)；
原归档清单作为历史身份记录保留，不表示清理后的完整文件集。
