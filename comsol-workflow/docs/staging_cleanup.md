# tmp / tests 分类与晋升候选

2026-10-09。范围包含两个目录的 archive；[实施 Spec](spec/plan-execute_20261009_staging_cleanup.md) 记录边界与回执。
已直接删除 175 个冗余文件，合计 838,867,553 字节（约 839 MB / 800 MiB）；
清理 12 个空目录。675 项离线回归通过；四项候选保留原件，等待用户决定，尚未晋升。

## 分类依据

对全部 322 个文件建立大小、修改时间和 SHA256 清单，检查核心/正式入口的导入、动态加载、
路径及文档引用；对 Python 做静态解析并比较函数 AST。没有相关 Python/COMSOL 运行进程。
调用图、归档本身或测试文件名不能单独证明冗余；回归、真实求解诊断和来源快照有独立职责。

| 分类 | 数量 | 处理及依据 |
| --- | ---: | --- |
| Python/字体缓存 | 74 | 删除；保留源码和正常缓存机制 |
| 合成测试图 | 4 | 删除；test_geometry.py 的输出，不是测试输入 |
| 已发布数据的相同副本 | 62 | 删除；与保留的正式结果逐文件 SHA256 一致 |
| 相同的中间核验副本 | 5 | 删除 run05 副本；run04 的相同历史状态保留 |
| 旧任务启动日志/PID/活动指针 | 18 | 删除；任务已结束，科学审计与配置保留 |
| 旧排版/后处理绘图暂存 | 11 | 删除；不作计算输入，原数值表及当前正式图保留 |
| 旧固定案例包装脚本 | 1 | 删除 find_px_fundamental.py；核心筛选算法及 quarter 调用、回归保留 |
| 持续回归模块 | 59 | 保留；保护几何、求解编排、数值方法、配置、输出及兼容合同 |
| 工作流支撑文件 | 9 | 保留；目录说明、包身份、输出隔离、启动保护、当前收集检查及原归档清单 |
| 晋升候选 | 4 | 保留；见下表，交由用户判断 |
| 人工 COMSOL 诊断入口 | 3 | 保留；真实求解/多层/新旧实现/导出一致性诊断，不由离线 runner 执行 |
| 人工 Fourier 检查演示 | 1 | 保留；值得转成纯检查，未被默认 pytest 收集 |
| 独有历史证据 | 71 | 保留；参数/manifest、科学失败与收敛结论、旧数值版本及完整回滚快照 |

两份原 archive_manifest 是归档当时的身份记录，不再表示当前完整文件集。
[逐项清单](spec/staging_cleanup_inventory_20261009.json) 列出删除原因、替代文件和恢复方式；历史日志/文档的旧路径保持为历史信息，
不能按旧 PID、哈希或绝对路径自动续跑。

## 功能对比与晋升候选（待用户决定）

| 编号 | 当前文件与价值 | 已有正式所有者及重复部分 | 建议方向、必要工作与限制 |
| --- | --- | --- | --- |
| P1，优先 | [bulk/cladding 预览](../tests/archive/test_bulk_cladding_lattice_generation.py)：比较 envelope offset、miter offset、连续 lattice shell；不同周期与 gap 下展示 cladding 首层及保留单元 | hex_lattice_utils / polygon_utils 已有边界、shell、相交及保留点算法；脚本中的 cladding_selection 与核心职责重叠 | 可做 tests/lattice_geometry/plot_cladding_boundary_comparison.py，验证后晋升 scripts/analysis。只保留比较编排；先核验接触/凹边界判定，不能盲目替换选择规则。参数化固定样例，统一输出、≥300 dpi 与字体，并补纯检查、查看实际 PNG |
| P2，中 | [bulk/side/corner 预览](../tests/archive/test_bulk_side_corner_lattice_generation.py)：正交/hex 的三分区、间距/数量校验及 PNG/CSV | finite_patch_lattice 已有 generate_lattice、坐标转换和分区算法；preview 不属于新核心算法 | 可做 tests/lattice_geometry/plot_lattice_partitions.py → scripts/analysis。复用核心；整理 CLI 与 case identity，统一绘图/输出，检查正交和 hex 两种基底 |
| P3，中 | [孔洞裁切预览](../tests/archive/test_cavity_hole_polygon_clipping.py)：gap 扫描、保留/移除孔洞以及面积/外侧顶点的可视核验 | hex_lattice_utils / polygon_utils 已有算法；10 个定义与正式函数 AST 完全相同，含 Fourier 参数转换及单元几何包装 | 可做 tests/lattice_geometry/plot_clipped_holes.py → scripts/analysis。复用现有几何计划与裁切 API；AST 相同仍须核验全局常量和输入身份，不能盲目删除包装。固定科学参数和裁切语义不能随重命名改变。提取纯检查，整理 CLI 与图合同 |
| P4，条件性 | [CSV 单点评审](../tmp/archive/s4_20260916-gamma-v2/review.py)：跨积分分辨率变化、边界迹差异、air-height 的绝对/相对误差及有限场量级比较 | 正式 boundary 入口已写 identity、parity、trace 等指标，但此脚本另含跨分辨率/高度评审编排；不是新的 COMSOL 求解算法 | 如需通用评审，将可复用评估先放 tmp/boundary_validation，CSV 入口放 tests/boundary_validation/analyze_boundary_convergence.py，再考虑核心/analysis。当前固定 mode=1、旧行数、旧字段和“预期失败”断言，仅用于历史案例；须明确数据 schema、门槛与近零误差分母并用合成/授权只读 fixture 验证 |

本轮不移动、改名或运行以上候选；保留原件，避免以“已标记”冒充验证或晋升。
可优先把 P1–P3 做成同一几何诊断主题下的独立入口，共用正式算法与绘图保存，不新增一套框架。

## 保留但不直接晋升

- [Fourier 手工演示](../tests/archive/test_fourier_basis.py) 的正交、重建、Parseval 和最小非零基值检查有价值；
  它目前只在 __main__ 执行，多个失败仅打印文本。后续适合转入 tests/fourier_basis/checks，
  加确定性的断言；算法已在 basis_utils，无须再合入核心或 scripts。
- test_simulation_spatial.py、test_simulation_utils_spatial_comsol.py、test_simulation_plotting.py
  是实际 COMSOL 诊断，覆盖 unified/spatial、多层场重叠、频率/Q 与导出；mock 回归不能替代它们。
  本轮没有执行或宣称可运行；目前排除自动调用，后续若统一真实诊断入口再单独评估。
- test_strip_1d_fourier_quotient.py 中四个 helper 与正式函数相同，但所在模块还验证 strip 商空间的
  计数、正交和重建，不能因局部重复删除整份回归。59 个持续回归均保留，本轮不削弱断言。
- tmp 中保留的 S31 中间数值与 S4 旧版 backup 不是当前算法候选；它们记录数值不稳定、
  边界迹/厚度导数不足以及过去数据版本。新 final 并不与所有旧结果相同，不能按目录名一并删除。
- 参数快照即使字节相同也承载不同 manifest 的身份；保留完整快照关系。正式结果的历史
  metadata 可能记录已归档的 tmp 路径，本轮不重写来源或哈希。

## 验证及恢复

- tmp 删除 101 个文件，保留 75 个；tests 删除 74 个文件，保留 72 个。清理后原文件
  合计约 5.1 MB；随后仅补充目录说明。全部 66 个原测试/预览模块及活动检查仍保留，断言未改。
- 保留文件的原件 SHA256、大小和修改时间在清理后核验一致；随后编辑的三个 README 单独备份并记录。
  核心和正式脚本源码、共享参数、pyproject.toml、uv.lock 未改。
  63,892 个正式结果的路径、大小、修改时间及链接目标与基线一致；用于替代已删副本的文件再次核验 SHA256。
- `.\.venv\Scripts\python.exe -B scripts/check_offline.py` 检查通过：675 项
  （核心 587、科学 88），与清理前一致；3 条既有 PytestReturnNotNoneWarning 未改动。
  本轮未启动 COMSOL、未重绘科学结果、未自行晋升、未 commit/push。
- 44 个本地文档链接/锚点与逐项 JSON 清单检查通过，230 个剩余 Python 文件内存编译通过，
  git diff --check 通过；离线 runner 已清理本轮 .work 产物，测试完成后再次核验正式结果清单未变。
- 恢复资料位于 I:\codex-task-temp\staging-cleanup-20261009：baseline.json 为原件清单，
  classification.json 为执行前分类，deleted.jsonl 为实际删除回执，before_deletion 保存已删脚本、
  旧渲染快照和启动日志，before_docs 保存本轮修改的三个 README。字节相同的副本从逐项清单
  replacement 指向的保留文件恢复；缓存/合成测试 PNG 按原生成器重建，不必恢复旧字节码。
  恢复仅作用于清单内原路径，先检查冲突；不覆盖正式结果或还原整仓工作树。

覆盖限制：外部未知个人调用无法穷尽；删除的包装器只用于固定旧案例，能力仍由核心筛选函数
及 run_finite_quarter 消费。三份真实 COMSOL 诊断与晋升候选没有本轮运行验收，因而保留；
历史清单中的旧绝对路径不被当作当前可执行合同。
