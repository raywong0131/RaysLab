# COMSOL Workflow

本文件维护项目结构、概念边界和科学语义；必须执行的流程在 [AGENTS.md](AGENTS.md)，
使用方法在 [README.md](README.md)，长期协作偏好在 [MEMORY.md](MEMORY.md)。
当前配置来自 scripts/parameter.json 和任务快照，当前状态来自真实进程/日志/summary。

## Structure and ownership

| 位置 | 内容与边界 |
| --- | --- |
| `tmp/<topic>/` | 新的可复用算法试验；可使用核心包，核心包不得反向依赖试验区 |
| `tests/<topic>/` | 新计算/分析/绘图入口的前置区；调用算法，避免内嵌第二套算法实现 |
| `tests/<topic>/checks/` | 纯离线验证，`test_*.py`；不自动调用 COMSOL |
| `comsol_workflow/` | 已验证、适合长期维护的算法与 COMSOL/后处理能力 |
| `scripts/run_main/`、`run_sweep/`、`analysis/` | 已晋升的主计算、COMSOL 扫描和保存数据分析入口 |
| `tmp/.work/`、`tests/.work/` | 有任务边界的运行产物与缓存；不存唯一试验源码 |
| `tmp/archive/`、`tests/archive/` | 已归档源码/产物；旧测试可显式回归，旧清单不自动恢复运行 |
| `docs/spec/`、`docs/adr/` | 工程计划/执行回执，以及必要的长期设计决定 |
| `scripts/.out/`、已批准 `results/` 专题根 | 受保护的正式科学结果；位置不由试验源码所在目录决定 |

算法回答“如何计算”，脚本回答“本次计算什么、用何来源、怎样记录/保存”。
晋升是把已验证实现交给正式所有者并更新调用关系，不是复制一个新版本。
完整步骤、命名和生成位置见 [AGENTS](AGENTS.md#新算法与新脚本的两条晋升通道)。

## Shared capabilities

- 几何、求解、模式识别、Fourier/far-field 沿用核心包；边界贡献使用 boundary_integrals / boundary_plotting。
- dipolar 投影/辐射、Gamma 分解和晶格分析分别由 dipolar_radiation、dipolar_decomposition、
  dipolar_lattice_analysis 维护；来源由 dipolar_inputs 显式配置。
- 路径、序列化、图保存分别由 output_paths、result_io、figure_output 维护；field 图使用 field_plotting。
- 旧 S4/S5 模块仍为同一实现的兼容入口；新能力使用语义名称，历史科学批次标识保持原义。

## Language

以下科学定义从原 Memory 原样集中到本文件，后续只在此维护，不推断新的物理假设。

- For the two in-plane photonic-crystal regions, use `cavity` and `cladding`
  as fixed English technical terms: `cavity` is the central region and
  `cladding` is the surrounding peripheral photonic crystal. Keep these terms
  untranslated in Chinese conversation, reports, and figure labels.

- Cavity periodic reference fields may cover only cavity cells. Restrict both
  the internal-field projection and basis-radiation integrals to cavity; do
  not extend these fields into cladding. Keep the complete finite-structure
  solution as an independent reference, with unexpanded regions in the remainder.

- The compact `(b0, eta, zeta)` unit-cell convention is a compatibility input
  for the existing Fourier geometry. Its invariant input mapping is
  `r_f0 = eta`, `b_square_f0 = (b0_nm / 230.0)^2`, and
  `b_square_f3 = (zeta^2 - 1.0) / 2.0`, where `230.0 nm` is the Fourier
  reference `B_0`. The nominal left/right triangle side is `b0*zeta`; the other
  four nominal sides are
  `b0*sqrt((3-zeta^2)/2)` to preserve fill factor. Apply this contract through
  the conversion entry point used by the workflow and retain the existing
  parameterization and clearance checks. Fourier coefficients, including exact
  zeros and small nonzero values, must be passed through without a minimum
  absolute-value clip. In particular, `zeta = 1` must keep
  `b_square_f3 = 0`, produce six equal triangle side lengths, and preserve the
  intended C6-symmetric geometry. Treat the run's recorded resolved hole
  parameters, triangles, or equivalent realized-geometry metadata, rather than
  the nominal compact inputs alone, as the realized geometry.

- The user's external `px`/`py` convention is reversed relative to the project:
  user-facing `px` corresponds to project-internal `py`, and user-facing `py`
  corresponds to project-internal `px`.

- In the user's convention, the primary electric-field area integrals for
  `px` and `py` are respectively integral Ex dS and integral Ey dS. Internal
  mode-name conversion does not swap electric-field coordinates. Ex of a py
  mode is a transverse reference, not a separate px-mode result.
  User-facing mode names and figure filenames use px/py, never qx/qy.
  Cartesian boundary-integral fields keep their internal definitions; their
  x/y suffixes must not be used as mode labels.

- Preserve project-internal `px` and `py` names in code, decomposition results,
  CSV/JSON fields, and cached data. All user-facing conversation, documentation,
  plot labels, legends, and scientific interpretation use the user's convention.
  Mention the project-internal component only when it is needed to diagnose or
  implement the mapping; never silently rename stored internal fields.

- Before an expensive calculation targeting one p component, state the target
  in the user's convention and verify its project-internal component during the
  technical preflight. Do not permanently infer the research target from an
  earlier run.

- Within the two selected p-dominant bands at Gamma, `p1` is strictly the
  lower-frequency Gamma band and `p2` is strictly the higher-frequency Gamma
  band. Away from Gamma, band tracking preserves those identities. Neither
  `p1` nor `p2` is permanently equivalent to project `px` or `py`; inspect the
  actual Gamma composition.

- Keep internal d-component fields unchanged. In user-facing mathematical
  labels, project `dy_weight` is displayed as `d_{xy}`, and project `dx_weight`
  is displayed as `d_{x^2+y^2}`.

- In the strip workflow, “bulk cell” means the interior cavity-type
  `cavity_p_bic` cell, while “cladding cell” means the surrounding
  `bulk_gap_centered` cell.

- “cladding shift”指 cladding 位移目标，旧称“geometry shift”指同一量。逐层值是相对理想晶格的
  绝对位移，不逐 shell 累加；groups/ellipse、hex/square、ratio 扫描/显式 x/y 见 [仿真合同](docs/agent-simulation.md)。

- **quarter / C2v 展示重建**：对称性约化求解与由对称性生成的展示视图是不同概念。
  Unit-cell 2D 镜像展示点不是独立 COMSOL 求解点，插值/镜像不得计作新增求解证据。
- S4 共轭 SI 边界采样与 S5 native 场采样为不同合同，不能为统一接口而隐式互换。

## Contract routing

参数见 [docs/parameters.md](docs/parameters.md)；几何、mesh、求解与恢复见
[docs/agent-simulation.md](docs/agent-simulation.md)；2D 扫描/valley 见
[docs/agent-unit-cell-2d.md](docs/agent-unit-cell-2d.md)；绘图见
[docs/agent-plotting.md](docs/agent-plotting.md)。既有精确输出/科学合同优先于通用新文件分类。
不得用历史参数、mesh、频率或 PID 替代本次配置和预检。
