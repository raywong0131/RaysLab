# unit_cell_2D ζ 参数扫描计划与执行记录

## 目标

在当前 scripts/parameter.json 的 unit_cell_2d 配置基础上，仅扫描
zeta = 1.151, 1.161。其余几何、mesh、二维 k 网格、模式追踪阈值与分析开关
保持不变，并在构建完成后直接启动本机 COMSOL 串行计算。

## 范围与数据流

- 新增 scripts/run_sweep/run_unit_cell_2d_zeta_scan.py 作为可恢复的扫描入口。
- 主参数文件只作为只读基准，不在扫描过程中改写。
- 每个 ζ 生成独立的完整参数快照，并通过
  COMSOL_WORKFLOW_PARAMETER_PATH 传给现有 run_unit_cell_2d。
- 两个 case 使用独立 Python 子进程串行执行，避免同一 Python 会话重复创建
  mph client；单个 case 内继续复用一个 COMSOL model 扫描全部 k 点。
- 每个 case 沿用既有 canonical series 命名、点级 cache identity、默认 Analysis
  和正式输出合同。
- 扫描级目录只保存 manifest、状态摘要和 stdout/stderr；数值结果仍位于各自
  unit_cell_2D_(...zeta...) series。

## 兼容与恢复

- 不改变 run_unit_cell_2d.py 的单 case 行为和 CLI。
- 已完成且 cache identity 一致的 k 点由原工作流自动跳过。
- 中断后重新运行同一扫描入口会复用各 case 参数快照和点级缓存。
- 若已有参数快照与本次基准配置冲突，则拒绝覆盖并停止。

## 验证

- 纯 Python 测试锁定：只修改 unit_cell_2d.zeta、case 目录互不冲突、
  主参数文件不变、子进程严格按输入顺序串行。
- 启动前运行聚焦 pytest、相关 py_compile 与 git diff --check。
- 启动后检查 PID、stdout/stderr、COMSOL 子进程和首个完整 k 点。

## 回滚

- 删除新增扫描入口、测试和 pyproject CLI 项即可回滚代码。
- 不通过删除正式计算结果回滚；已生成的 case 结果保留并可独立恢复。

## 执行记录

- 状态：实施中。
- 实际改动：
  - 新增可恢复的串行扫描入口和 pyproject 命令
    comsol-unit-cell-2d-zeta-scan。
  - 每个 ζ case 使用独立 Python 子进程、完整参数快照和既有点级缓存；
    主参数文件未修改。
  - 新增聚焦测试，并将过期的 0.1/19/361 网格断言同步为当前
    0.05/17/289 正式配置。
- 验证结果：unit-cell 2D 工作流、分析和扫描入口联合测试 47 passed；
  相关 py_compile 与 git diff --check 通过。
- 扫描目录：
  scripts/.out/unit_cell_2D/unit_cell_2D_zeta_scan_(245-0.96)_mesh5_kmax0.05_uniform17_band-p2_zeta1.151-1.161
- case 目录：
  - zeta=1.151：unit_cell_2D_(245-0.96-1.151)_mesh5_kmax0.05_uniform17_band-p2
  - zeta=1.161：unit_cell_2D_(245-0.96-1.161)_mesh5_kmax0.05_uniform17_band-p2
- 启动时间：2026-08-22 10:43；后台父进程 PID 32372。
- 启动检查：扫描级 stderr 为 0 字节；COMSOL mphserver 已启动；
  zeta=1.151 已生成 Gamma 和首批非 Gamma 点，监控时已有 3 个完整
  point_metadata.json。
- 当前状态：后台运行中，正在执行第 1/2 个 case。

## 2026-08-22 kmax0.15 / uniform25 三组扫描

- 用户指定保持其他参数不变，扫描 ζ=1.155、1.156、1.157。
- 主参数基准：`b0_nm=245`、`eta=0.96`、mesh5、目标 `p2`、
  `q_max_over_G=0.15`、`q_points_per_axis=25`。
- 每组为 25×25=625 个 k 点，三组共 1875 点；每点请求 2 对本征解
  （预期 4 个），eigensolver shift 固定为 `c_const/1.55[um]`。
- 每组使用一个固定几何 COMSOL model 串行复用全部 k 点；三组 case 之间继续
  使用独立 Python 子进程串行执行。点级缓存和默认 Analysis 均启用。
- 扫描目录预定为
  `scripts/.out/unit_cell_2D/unit_cell_2D_zeta_scan_(245-0.96)_mesh5_kmax0.15_uniform25_band-p2_zeta1.155-1.156-1.157`；
  启动前检查确认扫描目录和三个 case 目录均不存在。
- 启动前补齐既定 `00_model` 输出合同：新完成的 case 在释放复用模型前保存
  `geometry.png` 和 `unit_cell_2D.mph`。历史结果不在本次补写。
- 当前状态：启动前验证中；启动 PID、日志、首个完整点和最终结果将在下方继续回填。

### 启动记录

- 启动前聚焦测试共 `49 passed`；相关 `py_compile` 与
  `git diff --check` 通过。
- 2026-08-22 22:28:05 以隐藏后台进程启动；扫描父进程 PID 22892，
  当前 case Python PID 16212，COMSOL mphserver PID 55916。
- stdout：
  `80_logs/zeta_scan_20260822_222805_stdout.log`；
  stderr：
  `80_logs/zeta_scan_20260822_222805_stderr.log`。
- 启动检查时 stderr 为 0 字节，ζ=1.155 已完成 Gamma 与首批 7 个非 Gamma 点，
  共 8 个 `point_metadata.json`；全部为 `gamma/matched`，没有
  `missing/ambiguous`。
- 当前状态：后台正常运行，正在执行第 1/3 个 case；其余 case 按
  ζ=1.156、1.157 的顺序等待串行执行。
