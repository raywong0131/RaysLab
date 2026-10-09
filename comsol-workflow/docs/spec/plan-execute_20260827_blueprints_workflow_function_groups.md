# Blueprints workflow 功能分组执行计划

日期：2026-08-27
状态：已完成

## 目标结构

在保持 `workflow/` 物理目录扁平映射为逻辑包 `blueprint_core` 的前提下，按职责建立：

```text
workflow/
  __init__.py
  cli.py               # 流水线编排与命令行入口
  config.py            # 严格参数配置
  models.py            # 共享中性 layout 数据合同
  geometry/            # 晶格、孔形、外围开窗、参数标记
  fabrication/         # 局部/包络修正与 dose/process 接口
  quality/             # DRC 与版图质量门禁
  output/              # DXF、PNG、报告与 source provenance
```

## 迁移映射

- `auxiliary_geometry.py` -> `geometry/auxiliary.py`；
- `corrections.py`、`dose.py` -> `fabrication/`；
- `drc.py` -> `quality/`；
- `dxf.py`、`preview.py`、`reporting.py` -> `output/`；
- `cli.py`、`config.py`、`models.py` 保持根层，分别承担编排、参数源和跨功能共享合同。

不保留旧模块路径的重复 wrapper；仓库内调用和测试统一切换到新的功能包路径。公开项目
根包名仍为 `blueprint_core`。

## 验证与回滚

- 更新所有相对导入、测试导入、setuptools 子包声明、README 与 AGENTS；
- `uv sync` 后验证各功能包导入；
- 运行独立完整 pytest、`compileall`、默认 preflight 和 `git diff --check`；
- preflight 不生成 `.out`，不启动 COMSOL；
- 回滚时反向移动模块并恢复导入及包声明。

## 执行记录

### 实际改动

- `geometry/` 接收 `auxiliary.py`，统一容纳晶格、孔形、外围开窗和参数标记构建；
- 新建 `fabrication/`，容纳 `corrections.py` 和 `dose.py`；
- 新建 `quality/`，容纳 `drc.py`；
- 新建 `output/`，容纳 `dxf.py`、`preview.py` 和 `reporting.py`；
- `cli.py`、`config.py`、`models.py` 保持根层，承担编排、参数和共享数据合同；
- 更新全部内部相对导入、测试导入、setuptools 子包声明、README 和 AGENTS；
- 未保留旧路径 wrapper 或重复模块。

### 验证结果

- `uv sync --project blueprints`：通过；
- 新增功能包 `geometry`、`fabrication`、`quality`、`output` 全部导入通过；
- 独立完整测试：64 passed；
- `compileall`：通过；
- 默认 preflight：29,648 features，DRC 0 error、0 warning，未写文件；
- COMSOL 未启动，DXF/PNG/JSON 未生成。

最终状态：workflow 已按功能分类，公开根包名仍为 `blueprint_core`。
