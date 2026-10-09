# Blueprints 项目目录重构执行计划

日期：2026-08-27
状态：已完成

## 1. 目标

将 `blueprints/` 重构为可独立迁出的微纳加工版图项目，并形成清晰的源码、运行入口与
结果边界：

- 将原 `src/` 更名为 `workflow/`，核心可复用程序统一位于
  `workflow/`；
- 建立 `main/`，集中放置 `run_blueprint.py`、`parameter.json` 和运行产生的
  `.out/`；
- 删除现有版图结果以及项目内缓存、测试临时文件，不保留旧路径兼容副本；
- 新建 `blueprints/AGENTS.md`，固化独立项目的目录、参数、输出、加工安全与验证规则；
- 保持 `blueprint_core` Python 包名及 `blueprint-generate` 命令不变，使内部导入接口
  稳定。

## 2. 目标目录

```text
blueprints/
  AGENTS.md
  README.md
  pyproject.toml
  uv.lock
  workflow/
    ...核心几何、修正、DRC、DXF、预览与 CLI 模块
  main/
    run_blueprint.py
    parameter.json
    .out/                 # 运行时创建，不纳入版本控制
  tests/
    ...独立测试
```

根目录不再保留 `run_blueprint.py`、`parameter.json`、`.out/` 或 `src/`。测试仍单独
保留在 `tests/`，不与生产程序混放。

## 3. 路径与数据流

```text
main/parameter.json
  -> workflow/config.py 严格校验
  -> geometry + auxiliary geometry
  -> correction/dose/envelope
  -> DRC
  -> main/.out/<case-identity>/
       engineering_preview.dxf | production.dxf
       preview.png
       layout_report.json
```

`export.output_root` 继续采用相对参数文件的路径，因此配置值保持 `.out` 即自然落入
`main/.out/`。CLI 默认配置改为 `main/parameter.json`。源码 provenance 从
`workflow/` 收集，并包含根目录打包文件与 `main/run_blueprint.py`。

## 4. 兼容策略

- 保留 `blueprint_core.*` 导入路径，现有模块调用方只需重新安装项目，不改 Python
  import；
- `pyproject.toml` 的 setuptools package root 从 `src` 改为 `workflow`；
- 不保留旧入口和旧参数文件的双份副本，避免两套配置发生漂移；README 和测试全部切换
  到新路径；
- `source_trace.path` 相对参数文件重新计算到主项目 `scripts/parameter.json`；其摘要需与
  当前源文件一致。

## 5. 删除与清理范围

允许删除且不迁移：

- `blueprints/.out/` 下的全部既有 DXF、PNG、JSON case；
- `blueprints/.pytest_cache/`、所有 `__pycache__/`、`*.pyc`、`.test_tmp_*` 等缓存或
  临时文件；
- 迁移完成后为空的旧 `blueprints/src/` 路径。

不删除 `blueprints/.venv/`：它是独立项目的本机依赖环境，不是版图结果；验证过程仍需
复用它。若测试产生新的缓存，交付前再次清理。

## 6. 验证

- 检查根目录不再存在 `src/`、旧 `parameter.json`、旧 `run_blueprint.py` 和旧
  `.out/`；
- 检查 `main/.out/` 在清理后不存在，且 preflight 不创建结果目录；
- 运行 `blueprints` 独立完整测试；
- 对所有生产 Python 文件运行 `py_compile`；
- 运行 CLI `--preflight-only`，确认默认配置与新入口可用且不生成结果；
- 校验 `source_trace` SHA-256、包依赖独立性、`git diff --check`。

## 7. 回滚

目录迁移可以通过反向移动恢复：`workflow/` 回到 `src/`，`main/run_blueprint.py` 和
`main/parameter.json` 回到项目根，并恢复 setuptools、CLI、README 与测试路径。
已按用户要求删除的旧 `.out` 和缓存不纳入回滚，不能从本次重构恢复。

## 8. 执行记录

### 实际改动

- 将整个 `blueprints/src/` 原位迁移为 `blueprints/workflow/`，保留
  `blueprint_core` 包名和全部模块边界；
- 新建 `blueprints/main/`，将正式入口和唯一参数文件迁入；CLI 默认参数切换为
  `main/parameter.json`，配置中的 `.out` 因相对参数文件解析而落到 `main/.out/`；
- setuptools package root 切换为 `workflow`，独立 `uv sync` 已重新构建并安装成功；
- source provenance 改为收集 `workflow/**/*.py` 和
  `main/run_blueprint.py`；测试同步锁定新路径；
- `source_trace.path` 改为 `../../scripts/parameter.json`，摘要保持与当前外层参数源一致；
- 新建 `blueprints/AGENTS.md`，固化项目定位、目录、显式参数、几何扩展、工艺修正、
  dose、输出精简、release gate、安全与验证规则；
- README 已切换到新目录和默认命令；`.gitignore` 增加 editable-install egg-info 规则。

### 清理结果

- 已删除旧 `blueprints/.out/` 及其中全部既有 DXF、PNG、JSON case；这些结果按用户要求
  删除，不可由本次重构直接恢复；
- 已删除项目源码与测试中的 `.pytest_cache/`、`__pycache__/`、`*.pyc`、egg-info 和
  本轮 patch 临时文件；
- 保留 `blueprints/.venv/` 作为独立依赖环境，不把它视为版图或任务临时结果；最终验证后
  再次清理测试生成的缓存。

### 验证结果

- `uv sync --project blueprints`：通过，包从 `workflow` 成功构建安装；
- 从 `blueprints/` 运行完整独立 pytest：64 passed；
- `compileall`：通过；
- `main/run_blueprint.py --preflight-only`：通过，29,648 features，DRC 0 error、0 warning，
  明确未写文件；
- `source_trace` SHA-256：与当前 `scripts/parameter.json` 匹配；
- 生产源码外层依赖扫描：未发现 `scripts`、`comsol_workflow` 或 `mph` import；
- COMSOL 未启动，实际 DXF/PNG/JSON 未重新生成。

最终状态：目录重构已完成；当前仍为 `engineering_preview`，工艺 bias、dose、polarity 和
生产 DRC 阈值的标定状态不变。
