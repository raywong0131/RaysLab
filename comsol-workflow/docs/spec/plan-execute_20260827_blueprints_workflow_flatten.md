# Blueprints workflow 目录扁平化执行计划

日期：2026-08-27
状态：已完成

## 目标

删除 `blueprints/workflow/blueprint_core/` 这一层，使核心文件直接位于
`blueprints/workflow/`，几何子包位于 `blueprints/workflow/geometry/`。

## 兼容与数据流

- 保留公开 Python 包名 `blueprint_core` 和 `blueprint_core.*` 导入路径；
- 使用 setuptools 显式 package mapping，将逻辑包 `blueprint_core` 映射到物理目录
  `workflow/`，并将 `blueprint_core.geometry` 映射到 `workflow/geometry/`；
- CLI 中项目根路径从 `cli.py` 的新物理位置重新推导；
- source provenance、README、AGENTS 和测试期望全部切换到扁平路径；
- `main/parameter.json`、`main/run_blueprint.py` 与 `main/.out/` 合同不变。

## 验证与回滚

- 重新执行 `uv sync`，确认 `blueprint_core` 仍可导入；
- 运行独立完整 pytest、源码编译、默认 preflight 和 `git diff --check`；
- 确认 `workflow/blueprint_core/` 不存在且没有生成 `.out`；
- 回滚时将核心文件重新移入 `workflow/blueprint_core/`，并恢复 setuptools 与路径推导。

## 执行记录

### 实际改动

- 将 `workflow/blueprint_core/` 内全部 Python 文件及 `geometry/` 上移到 `workflow/`；
- 删除移动后确认已为空的 `workflow/blueprint_core/`；
- setuptools 改为显式映射：逻辑包 `blueprint_core` 对应物理目录 `workflow/`，公开导入名不变；
- CLI 项目根推导、source provenance、README、AGENTS 和路径测试已同步到扁平结构；
- `main/` 入口、参数和 `.out` 合同未改变。

### 验证结果

- `uv sync --project blueprints`：构建及安装通过；
- `import blueprint_core, blueprint_core.cli, blueprint_core.geometry`：通过，包文件来自
  `blueprints/workflow/__init__.py`；
- 独立完整测试：64 passed；
- `compileall`：通过；
- 默认 preflight：29,648 features，DRC 0 error、0 warning，未写文件；
- COMSOL 未启动，DXF/PNG/JSON 未生成。

最终状态：`workflow/blueprint_core/` 已不存在，核心程序直接位于 `workflow/`。
