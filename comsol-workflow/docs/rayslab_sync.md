# RaysLab 项目同步记录

本次将 COMSOL Workflow 当前项目快照同步至
[`raywong0131/RaysLab/comsol-workflow`](https://github.com/raywong0131/RaysLab/tree/main/comsol-workflow)，
目标分支为 `main`。源项目基准为
[`1d591a49`](https://github.com/jia-yli/comsol-workflow/commit/1d591a494fe9a090c8de606ab3c3d0731899d6f4)。

## 同步范围

- 原项目 434 个版本文件：全部代码、文档、共享配置、依赖及锁文件、编辑器配置、
  离线测试、`tmp/archive/` 和 `tests/archive/` 中的现有项目记录。
- 一份原先被忽略、但离线几何回归必需的配置快照：
  `results/finite_cavity_reference_results/finite_cavity_xy_fourier26_bulk5_cladding6_wide_bounds_mma_move0p01_eval0031/quarter_simulation_config.json`。
  该文件为 179,388 字节的只读配置与 607 个孔洞多边形，不含求解模型或场结果；按原路径、原字节纳入。
- 共复制 435 个源文件，9,524,623 字节；源工作区检查时没有未提交或未跟踪项目文件。
- 新增本报告、同步 Spec 和文件身份清单。保留目标根 `.gitignore` 原内容，追加
  `!comsol-workflow/`，项目自身忽略规则继续约束输出位置。
- RaysLab 原先没有 `comsol-workflow/`；其他已有项目保持原样。源项目的 Git、远端及历史保持原样，
  源 `.git` 不进入目标目录。此处导入项目快照，不复制或重写源提交历史。

## 明确排除项

用户已确认历史求解模型与数据保留本地。扫描的历史结果根如下，必要配置快照是唯一例外：

| 结果根 | 扫描文件数 | 字节数 | MPH 模型数 |
| --- | ---: | ---: | ---: |
| `results/` | 4,660 | 29,664,239,604 | 118 |
| `scripts/.out/` | 59,232 | 192,867,070,946 | 212 |

上述历史目录合计约 222.5 GB，包含求解模型、场数据、NPZ/Parquet、结果图、报告、日志和运行配置。
目录链接未跟随，因此数量不包含 junction 指向的重复或外部目录。文件数量和大小为只读扫描记录，
没有为全部大文件计算内容哈希。

另外排除 `.git`、嵌套 Git 元数据、虚拟环境、Python/pytest 缓存、包安装元数据、
`.understand-anything/` 的生成状态、`tmp/.work/`、`tests/.work/`、运行临时文件、密钥及凭据。
`uv.lock` 是必要依赖锁文件；`tmp/archive/`、`tests/archive/` 是项目记录，均保留。
上传内容未检出私钥、真实令牌、凭据 URL 或敏感配置文件。

## 模型与数据用途、LFS

基础仿真由代码和 `scripts/parameter.json` 创建模型，项目启动不要求携带历史 MPH。
重绘、恢复和部分 boundary/dipolar 分析需要对应历史输入，此次没有把这些模型与数据上传。
例如 dipolar 预设的 finite 模型约 12 GB，Fourier 输入约 65.1 MB；具体路径、用途和本机存在状态
记录在[同步清单](spec/rayslab_sync_manifest_20261009.json)的 `historical_input_dependencies` 中。

历史预设的 `boundary_case_dir` 在源项目中已不存在；同步保留这一现有状态。运行相应历史流程前，
需提供有效的 `--source-manifest` / `COMSOL_WORKFLOW_DIPOLAR_INPUTS` 并准备对应模型、场数据和来源配置。
这份代码快照不代表历史数据已具备。

已检查 Git LFS 3.7.1 和目标仓库 LFS 端点。GitHub 普通 Git 单文件限制为 100 MiB，
[LFS 上限随方案为 2/4/5 GB](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage)。
本机有 108 个历史文件超过普通 Git 限制，其中 11 个超过 5 GiB，最大约 21.1 GB。
这类文件不能直接通过 GitHub LFS 保存完整单文件。本次没有拆分模型、重写历史或上传 LFS 对象；
纳入的文件最大 677,423 字节，无需 LFS。超限文件清单附于同步清单。

## 验证与使用条件

- 对全部源文件逐一校验 SHA-256；提交前同时核对 Git 暂存内容，确保复制未改变原始字节。
- 230 个 Python 文件编译检查、50 个源 JSON 配置解析和 1 个 TOML 解析通过。
- 离线回归使用副本中的代码及本机已有 Python 环境：587 项核心检查和 88 项科学检查通过，
  合计 675 项；不启动 COMSOL。
- 原回归包含三条已有 `PytestReturnNotNoneWarning`。首次检查发现缺失的只读几何配置及 Windows
  长临时路径限制，补入配置并使用较短的外部测试临时路径复验；未修改测试断言或项目代码。
- 新增同步文档和目标根忽略例外通过空白检查。原快照有 17 个文件的 126 处既有空白提示，
  本次按原字节保留，不将同步变成格式整理。
- `microfab-blueprints` 是源项目已有的兄弟开发依赖，不属于此次目录同步；完整开发环境需另行
  提供 `../blueprints`，生产依赖可使用 `uv sync --no-dev`。实际仿真仍需 COMSOL 6.3 和有效许可证。

执行边界及最终验证见[同步 Spec](spec/plan-execute_20261009_rayslab_sync.md)。
逐文件身份与排除范围见[同步清单](spec/rayslab_sync_manifest_20261009.json)。
提交与推送只在独立临时克隆中执行，不强制推送、不删除 RaysLab 现有内容。
