# COMSOL workflow 保功能整合与命名规范

日期：2026-10-09。状态：已按先写规格、后改代码的顺序完成整合及离线验证。

## 目标与版本范围

用户已授权按扫描方案执行，并要求先写 spec、保持当前功能完整。
基线为再上一次 push `39421be6`（2026-09-03 14:53），当前 HEAD 为
`d1902092`；范围包含其后的未提交及未跟踪代码。保留任务开始前全部修改，
不恢复已迁至兄弟项目的 blueprints，不自动提交或推送。

目标是让实际复用能力有唯一所有者、专题程序采用能力名称、任务来源显式可配置、
后续输出使用一致的分类与保存接口。保留 band pair/solo、unit-cell 2D/refinement、
strip、full/quarter/all-symmetry finite、边界分析与全部 dipolar 分析能力。

## 不可改变的合同

- 不启动 COMSOL、不运行 mesh/eigensolve；不移动、覆盖、删除正式结果、参考场或 MPH。
- 不改变 compact/Fourier 参数映射、px/py 内外命名对应、电场坐标、模式筛选、相位、
  求积权重、单位、镜像重建、辐射相干交叉项和科学验收状态。
- S4 共轭 SI 采样与 S5 原生采样是不同合同，不合并为隐式转换。
- 保留 schema v1/v2 读取、缓存/config 身份校验、源码与来源哈希、staging、
  原子发布、mesh 复用及四扇区解归档。旧准备清单不得因改名被静默改写或放宽校验。
- 旧模块/命令保留兼容入口并指向同一实现；monkeypatch 的模块状态也须共享。
  历史结果路径与字段保持可读取。旧记录中的程序名和绝对来源不批量改写。
- 已有确切输出合同优先于通用分类；对新格式提供统一路径入口，旧格式保留读取适配。
  S4/S5 已批准的专题结果根保留，主 PNG 平铺、格式在外、批次在内。

## 实施步骤

1. 保存源码基线、Git 状态、结果 metadata 与当前进程信息；检查离线测试的求解边界。
2. 将 S4/S5 实现改为能力名称，更新活动源码调用；旧路径用轻量模块别名及命令适配保留。
   保留本身已表达物理职责的 dipolar、finite 迁移/恢复程序，不为改名而改名。
3. 提取已被多个入口调用的结果 JSON/CSV/哈希操作及纯辐射/投影操作到包内。
   复用现有几何、Fourier、far-field 和 field_plotting；不引入插件框架或新依赖。
4. 为专题输入提供显式 run/source/output 参数或来源清单；历史默认值继续可复现。
   数据不足时明确失败，不猜测 cell 拓扑，不隐式使用共享求解字段替代历史几何。
5. 集中输出分类、figure 保存与产物记录，接入现有生成器；数值 CSV/JSON 按职责分类。
   通过格式版本及读取回退保留历史路径。仅修改生成器，不重绘正式结果。
6. 接入少量稳定公共 CLI，更新活动入口目录、功能对应表与输出合同。
7. 完成聚焦和必要的全离线测试、源码编译、旧入口与新入口等价检查、diff 检查；
   回填实际执行范围、保留项、未验证边界与恢复路径。

## 命名映射

以下是规范实现位置；完整实际映射在执行回执记录。旧位置保留兼容适配。

| 原位置/名称 | 规范位置/名称 |
| --- | --- |
| `comsol_workflow/s4_boundary.py` | `comsol_workflow/boundary_integrals.py` |
| `comsol_workflow/s4_plotting.py` | `comsol_workflow/boundary_plotting.py` |
| `comsol_workflow/s4_paths.py`、`s5_paths.py` | 公共输出路径模块及旧布局适配 |
| `run_main/run_s4_boundary.py` | `run_main/run_boundary_analysis.py` |
| `run_main/run_s4_mesh_convergence.py` | `run_sweep/run_boundary_mesh_convergence.py` |
| `run_main/run_s4_surface_control.py` | `run_main/run_boundary_surface_control.py` |
| `run_main/run_s4_surface_zeta.py` | `run_sweep/run_boundary_surface_zeta_scan.py` |
| `run_sweep/run_s4_*` | `run_sweep/run_boundary_*`（使用实际扫描职责命名） |
| `analysis/*s4*` | `analysis/*boundary*`（保留 compare/decompose/export/plot/reprocess/validate 动词） |
| `run_main/run_s5_cladding_scattering.py` | `run_main/run_dipolar_cladding_scattering.py` |
| `run_main/run_s5_periodic_reference.py` | `run_main/run_dipolar_periodic_reference.py` |
| `run_main/run_s5_gamma_origin.py` | `run_main/export_dipolar_gamma_boundary.py` |
| `run_main/run_s5_gamma_replay.py` | `run_main/run_dipolar_boundary_replay.py` |
| `run_main/run_s5_volume_export.py` | `run_main/run_dipolar_volume_export.py`（阶段明确是否求解） |
| `run_main/run_s5_mesh_radiation.py` | `run_main/run_dipolar_mesh_radiation.py`（包含周期本征求解） |
| `run_main/run_s5_cap_height.py` | `run_main/export_dipolar_cap_height.py` |
| `analysis/*s5*` | `analysis/*dipolar*`（描述实际分析职责） |
| `run_main/run_finite_quarter_four_modes.py` | `run_main/run_finite_quarter_all_symmetries.py` |

## 输出原则

统一接口区分模型、数值数据、PNG、PDF、报告、日志、配置与历史资料；默认数值数据
进入 `01_results`，PNG 进入 `10_overview`，PDF 进入 `11_pdf`，报告进入 `12_reports`，
日志进入 `80_logs`，参数/来源进入 `99_config`。不以扩展名猜测 JSON 的科学职责。
已保存的旧布局以及有明确字段/路径合同的旧接口继续支持，兼容分支标明布局版本。
共同风格仅负责保存和显示约定，不把某张图的归一化、clim 或科学标签推广为通用值。

## 验收与恢复

- 新旧模块可导入并共享同一状态，旧 `python -m`/必要的直接文件入口继续可用。
- 现有模式/求积/FFT/辐射测试保护数值与科学合同；补充跨模块复用、来源配置、
  输出分类和兼容入口的行为检查，不以复制实现生成预期值。
- 离线测试显式拒绝 `mph.start()`；排除已知真实 COMSOL 示例入口。mock 求解测试
  仍须验证四扇区次数、独立 dataset/solution、恢复身份、异常清理和结果完整性。
- 图保存修改使用合成数据生成 PNG/PDF 并检查实际 PNG；不重绘用户结果。
- 检查共享参数及正式结果未被本任务改写；活动进程产生的变化单独报告。
- 源码先备份任务开始时的实际工作树；回退只恢复本任务文件，不能用整仓 reset/checkout。
  本轮没有数据迁移，因此不需要结果恢复操作。

## 执行回执

2026-10-09 完成。本轮未启动 COMSOL、未重绘真实结果、未迁移结果、未 commit/push。
实现前的离线基线为 569 项核心测试和 88 项科学测试，共 657 项；最终为
586 项核心测试和 88 项科学测试，共 **674 项通过**。

### 实际范围与共享所有者

- 37 个实现采用能力名称，39 个旧模块保留兼容别名（两套旧路径模块共享一个新所有者）。
  兼容层通过 `sys.modules` 指向实际模块对象，导入、全局状态和 monkeypatch 共享；
  旧 `python -m` 和直接文件调用经 `runpy` 转到规范实现。
- 22 个数值函数抽入 `dipolar_radiation`、`dipolar_decomposition`、`dipolar_lattice_analysis`。
  函数签名、默认值及主体 AST 与任务前完全相同；原脚本继续导入并暴露原公共名称。
  数值消费者直接依赖包内所有者，保留原 phase、sign、单位、mesh 和求积合同。
- `result_io` 统一 SHA256、JSON 读取/原子写入及复数 CSV 列处理；边界准备清单增加
  实际共享 IO 所有者的源码哈希。旧准备清单因源码身份变化仍会被拒绝，须重新 prepare，
  不静默修改已保存清单。
- `output_paths` 统一按职责分类，`table_path` 保留批次并支持旧 CSV 读取；
  `table_links` 只改生成报告中的 CSV 链接。`figure_output` 已接入 8 个主要 dipolar 图保存入口，
  原图的坐标、clim、phase、布局和各格式 DPI 设置保留，Gamma PDF 的批次分组也保留。
- `dipolar_inputs` 集中历史来源预设并提供 schema v1 清单；原始场目录与 Fourier 总览目录
  分别配置。`workflow_cli` 只做白名单分发、参数转发及环境恢复，目标程序保留原解析器。
  禁止参数缩写，防止目标 `--source` 被错误消费。
- 更新根 README、脚本目录、包目录及活动边界文档；旧科学规格、历史元数据和结果名称未改。

任务前生产源码为 111 个文件、54,438 行；任务后为 158 个文件、55,232 行。
增加的文件主要是规范实现、兼容入口和显式来源/输出校验，总行数净增 794 行。
本轮简化的是共享能力的所有权和任务间依赖，不以删除能力或兼容性换取行数下降。
与实际工作树快照相比修改 57 个文件，增加 48 个文件；原文件均保留。

### 完整名称对应

旧结果中的 S4/S5 标记、物理模式名、模型 feature/tag、数据列和历史绝对来源保持原义。
`run_dipolar_mesh_radiation` 包含周期本征求解，因此使用 `run_`，不标为纯 export。

| 历史入口 | 当前唯一实现 |
| --- | --- |
| `comsol_workflow/s4_boundary.py` | `comsol_workflow/boundary_integrals.py` |
| `comsol_workflow/s4_paths.py` | `comsol_workflow/output_paths.py` |
| `comsol_workflow/s4_plotting.py` | `comsol_workflow/boundary_plotting.py` |
| `comsol_workflow/s5_paths.py` | `comsol_workflow/output_paths.py` |
| `scripts/analysis/analyze_s5_cladding_scattering.py` | `scripts/analysis/analyze_dipolar_cladding_scattering.py` |
| `scripts/analysis/analyze_s5_gamma_origin.py` | `scripts/analysis/analyze_dipolar_gamma_origin.py` |
| `scripts/analysis/compare_s4_ey_hz_boundary.py` | `scripts/analysis/compare_boundary_ey_hz.py` |
| `scripts/analysis/decompose_s4_c6v.py` | `scripts/analysis/decompose_boundary_c6v.py` |
| `scripts/analysis/export_s4_ey_area_scan.py` | `scripts/analysis/export_boundary_ey_area_scan.py` |
| `scripts/analysis/plot_s4_boundary.py` | `scripts/analysis/plot_boundary_analysis.py` |
| `scripts/analysis/plot_s4_cached_zeta.py` | `scripts/analysis/plot_boundary_cached_zeta.py` |
| `scripts/analysis/plot_s4_control.py` | `scripts/analysis/plot_boundary_control.py` |
| `scripts/analysis/plot_s4_dual_validation.py` | `scripts/analysis/plot_boundary_dual_validation.py` |
| `scripts/analysis/plot_s4_signed_scan.py` | `scripts/analysis/plot_boundary_signed_scan.py` |
| `scripts/analysis/plot_s5_gamma_mechanism.py` | `scripts/analysis/plot_dipolar_gamma_mechanism.py` |
| `scripts/analysis/plot_s5_volume_radiation.py` | `scripts/analysis/plot_dipolar_volume_radiation.py` |
| `scripts/analysis/prepare_s5_gamma_origin.py` | `scripts/analysis/prepare_dipolar_gamma_origin.py` |
| `scripts/analysis/reprocess_s4_boundary.py` | `scripts/analysis/reprocess_boundary_analysis.py` |
| `scripts/analysis/validate_s4_local_fields.py` | `scripts/analysis/validate_boundary_local_fields.py` |
| `scripts/analysis/validate_s5_bloch_radiation.py` | `scripts/analysis/validate_dipolar_bloch_radiation.py` |
| `scripts/analysis/visualize_s5_broadening.py` | `scripts/analysis/visualize_dipolar_broadening.py` |
| `scripts/analysis/visualize_s5_broadening_2d.py` | `scripts/analysis/visualize_dipolar_broadening_2d.py` |
| `scripts/run_main/run_finite_quarter_four_modes.py` | `scripts/run_main/run_finite_quarter_all_symmetries.py` |
| `scripts/run_main/run_s4_boundary.py` | `scripts/run_main/run_boundary_analysis.py` |
| `scripts/run_main/run_s4_mesh_convergence.py` | `scripts/run_sweep/run_boundary_mesh_convergence.py` |
| `scripts/run_main/run_s4_surface_control.py` | `scripts/run_main/run_boundary_surface_control.py` |
| `scripts/run_main/run_s4_surface_zeta.py` | `scripts/run_sweep/run_boundary_surface_zeta_scan.py` |
| `scripts/run_main/run_s5_cap_height.py` | `scripts/run_main/export_dipolar_cap_height.py` |
| `scripts/run_main/run_s5_cladding_scattering.py` | `scripts/run_main/run_dipolar_cladding_scattering.py` |
| `scripts/run_main/run_s5_gamma_origin.py` | `scripts/run_main/export_dipolar_gamma_boundary.py` |
| `scripts/run_main/run_s5_gamma_replay.py` | `scripts/run_main/run_dipolar_boundary_replay.py` |
| `scripts/run_main/run_s5_mesh_radiation.py` | `scripts/run_main/run_dipolar_mesh_radiation.py` |
| `scripts/run_main/run_s5_periodic_reference.py` | `scripts/run_main/run_dipolar_periodic_reference.py` |
| `scripts/run_main/run_s5_volume_export.py` | `scripts/run_main/run_dipolar_volume_export.py` |
| `scripts/run_sweep/run_s4_dual_validation.py` | `scripts/run_sweep/run_boundary_dual_validation.py` |
| `scripts/run_sweep/run_s4_full_refined_scan.py` | `scripts/run_sweep/run_boundary_full_refined_scan.py` |
| `scripts/run_sweep/run_s4_peak_resampling.py` | `scripts/run_sweep/run_boundary_peak_resampling.py` |
| `scripts/run_sweep/run_s4_px_mesh_test.py` | `scripts/run_sweep/run_boundary_px_mesh_test.py` |
| `scripts/run_sweep/run_s4_surface_s_scan.py` | `scripts/run_sweep/run_boundary_surface_signed_scan.py` |

### 公共命令与来源清单

原有 band、unit-cell 2D、strip、finite 和 quarter 命令保留，新增：

| 命令 | 用途 |
| --- | --- |
| `comsol-finite-quarter-all` | 四种边界对称性完整谱，原 four-modes 能力 |
| `comsol-boundary-analysis` | prepare / run / report / check-model |
| `comsol-boundary-scan <stage>` | surface-zeta、mesh-convergence、dual-validation、full-refined、peak-resampling、px-mesh、signed-surface |
| `comsol-dipolar-analyze <stage>` | 15 个保存数据分析阶段，`--help` 列出全部阶段 |
| `comsol-dipolar-reference` | 周期参考，实际 `--run` 可求解 |
| `comsol-dipolar-volume-export` | 保存解的体场导出，实际运行可打开 COMSOL |
| `comsol-dipolar-boundary-replay` | 按原 config / stage / execute 合同重放 |

命令通过 `pyproject.toml` 注册；已用 `uv pip install --no-deps -e .` 刷新当前可编辑安装。
未添加依赖，`uv.lock` 与任务前一致。

```powershell
uv run comsol-dipolar-analyze broadening-2d --describe-inputs
uv run comsol-dipolar-analyze broadening-2d -- --help
uv run comsol-dipolar-analyze broadening-2d --source-manifest .\sources.json --describe-inputs
```

`--describe-inputs` 仅读取清单并打印已解析路径，不读取科学数据、不创建目录。
也可在进程导入任务模块前设置 `COMSOL_WORKFLOW_DIPOLAR_INPUTS`。
以下是路径结构示例，不是本轮运行配置；路径相对于清单自身解析：

```json
{
  "schema_version": 1,
  "series_dir": "saved-series",
  "mode_dir": "saved-series/10_overview/symmetry_1_xPEC_yPMC_mode19",
  "field_mode_dir": "saved-series/01_results/saved-case/mode19",
  "source_config_dir": "saved-series/99_config/source_runs/saved-case",
  "finite_model": "saved-series/00_model/saved-case/finite_quarter.mph",
  "run_id": "review-run"
}
```

额外字段为 `fourier_source`、`reference_dir`、`boundary_case_dir`、`output_dir`、
`request_dir`、`source_batch`、`gamma_group`。改变 series 必须同时显式给出四个关联路径；
输出路径限制在已授权 S5 根内。位置可配置不代表任意几何/频率均通过现有科学校验，
原固定科学任务检查继续执行。

### 输出与历史兼容

- 公共 dipolar 命令设置 `COMSOL_WORKFLOW_OUTPUT_LAYOUT=unified`，数值 CSV/NPZ 进入
  `01_results`；PNG/PDF/报告/检查日志/配置分别进入 10/11/12/80/99 分类。
  原程序直接调用默认保持 v1 布局。交互式 broadening 导出继续使用保留的 HTML 模板，
  修正了改名后的模板定位，并将模板纳入 wheel 包资源；HTML 仍内嵌 PNG/PDF/CSV/JSON。
  JSON 在通用接口中必须显式指定职责；
  历史 `output_path` 中 summary/check JSON 仍归日志、preflight JSON 归配置。
- 数值 CSV 读者优先当前布局并回退另一布局，读取不创建目录。批次目录在各格式目录内部
  保留；主 PNG 平铺。Gamma 全局共相位及 PDF 运行分组均保留。
- 既有 finite v6 总览 CSV、S4 compact/legacy 数据位置及全部明确输出合同未改。
  S4 SI/共轭采样与 S5 native 采样未合并；finite/strip 中依赖不同全局状态的相似函数未合并。
- 本轮只修改生成器和读取代码，没有改写、迁移或删除任何保存结果。

### 验证与结果保护

1. `scripts/check_offline.py`：最终 **586 + 88 = 674 passed**。
   核心与科学测试在独立进程执行，阻止真实 `mph.start()`，仅排除三个已知真实 COMSOL 示例：
   `test_simulation_spatial.py`、`test_simulation_utils_spatial_comsol.py`、`test_simulation_plotting.py`。
   保留 mocked 求解器验证。出现的三个 `PytestReturnNotNoneWarning` 在基线中已存在。
2. 新增跨模块测试验证全部 39 个旧/新模块对象相同、可共享 monkeypatch 状态；
   验证默认来源、相对清单、错误清单不写入、清单到实际入口的进程级传播、
   CLI 参数转发/异常恢复/describe-only、输出角色/批次/历史读取及严格 JSON/复数 CSV。
   最后修复后的聚焦测试 40 项通过，包含使用合成来源的完整 HTML 导出，
   验证资源内嵌、来源哈希保持和安全数据替换。
3. 合成数据实际生成 PNG/PDF，并查看了共享保存器生成的 PNG；
   图像可打开、标题完整，数值 clim/坐标范围保持；未读取真实场重绘。
4. 224 个 Python 文件内存编译通过；22 个抽出数值函数与原始 AST 完全相同；
   11 项新/旧命令 help 或来源展示通过。`git diff --check` 通过；
   另存针对实际任务前快照的 diff，涵盖 Git 未跟踪文件。
   `uv build --wheel` 通过；检查 wheel 含原 HTML 模板的完整内容及七个新增命令。
5. 全部 63,892 个正式结果文件的大小、mtime_ns 与快照一致，
   总字节 **222,531,310,550**；共享参数、依赖锁、既有测试、科学映射/求解核心及用户指令
   与任务前哈希一致。没有对 222 GB 结果重新做全量内容哈希；核验边界为完整文件清单、大小和时间。

运行与科学实测边界：COMSOL 启动、mesh、eigensolve、真实导出、重放和正式结果重绘
均未执行；本轮完成的是源码/输入输出兼容及完整离线回归，不能据此宣称新的科学结果已验收。

### 恢复与收尾

恢复快照：`C:/Users/Administrator/AppData/Local/Temp/comsol-integration-20261009-ew3gnwjo`。
`before/` 保存任务前实际工作树（含原有未提交代码），`baseline.json` 保存源哈希/Git 状态/结果清单，
`naming-map.json` 保存名称对应，`verification.json` 保存核验结果，`task-source.diff` 保存任务实际改动，
`final-tests.log` 保存最终完整离线测试输出，`wheel/` 保存验证过的 wheel。

需要回退时，仅将 `verification.json.changed` 中对应文件从 `before/` 恢复，
新增文件按 `added` 清单移入该任务快照内的回退保留目录，再按恢复后的 pyproject 刷新可编辑安装。
不要执行整仓 reset/checkout，不删除用户已有未提交文件，不触碰正式结果或 `.archives`。
本轮临时测试目录已移出项目，验证文件保留在恢复快照的 `test-artifacts/`，
剩余 pytest 临时目录/链接以及本轮 wheel 构建目录保留于
`I:/codex-task-temp/workflow-integration-20261009/`。
工具策略拦截了递归删除，因此采用可恢复归档；跨盘移动遇到 pytest 符号链接后，
其余目录通过同盘移动完成。检查图另存为快照内 `synthetic_figure.png`。
现有用户缓存、历史临时目录和既有 egg-info 保留。
