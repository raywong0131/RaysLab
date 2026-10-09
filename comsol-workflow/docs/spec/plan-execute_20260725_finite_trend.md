# Finite Gamma 趋势分析：计划与执行记录

日期：2026-07-25  
修订：2026-07-26（按审核意见调整 py-mode 权重判据并增加 Q 趋势图）  
状态：已完成

## 实施前方案概述

### 目标

在 `scripts/analysis/` 下新增 `finite_trend.py`，对一组 full finite 或
finite-quarter 的 cladding inward shift 计算结果执行与
`srip1d_trend.py` 等价的 Gamma 趋势分析：

1. 从系列级 `run_summary.json` 定位本轮声明的 shift case；
2. 从每个 case 的 `10_overview/mode_scores.csv` 中筛选
   `gamma_p_px_weight_fraction > 0.80` 的候选模式。该列在本工作流中实际对应
   py-mode 的模式权重，是 finite trend 的核心筛选指标；
3. 在候选模式中选择频率最接近目标频率的模式，频率距离相同时按
   `mode_idx` 由小到大决定；
4. 汇总所选模式并生成两张图：
   - `finite_trend.png`：沿用 `strip1d_trend.png` 的双纵轴布局，绘制
     `gamma_p_px_weight_fraction`、frequency 与目标频率参考线；
   - `finite_q_trend.png`：绘制同一批目标模式的 Q 与 cladding shift；
5. 支持直接从已生成的趋势 CSV 同时重绘两张图，不启动 COMSOL、不重算 Fourier
   或 mode score。

本任务只消费已有结果，不修改、迁移、覆盖或删除正式计算结果。

### 已确认的输入差异

`srip1d_trend.py` 不能原样用于 finite，原因如下：

| 项目 | Strip 1D | Full finite | Finite quarter |
| --- | --- | --- | --- |
| 系列索引 | `run_summary.json` | `run_summary.json` | `run_summary.json` |
| case 目录 | cut/shift case | `shiftX.XXX` | `shiftX.XXX`，多 symmetry 时再嵌套 symmetry 目录 |
| score 文件 | case 根目录下 `mode_scores.csv` | `10_overview/mode_scores.csv` | `10_overview/mode_scores.csv` |
| summary 中的 shift | 显式字段 | 显式字段 | 当前每条 run 没有显式字段 |
| 特有维度 | cut angle | 无 | `symmetry_id` |

当前 finite-quarter 正式样例还表明：系列目录中可能保留未写入
`run_summary.json` 的其他完整 shift 结果。分析使用 summary 提供系列类型、目标频率和
symmetry 元数据，同时扫描同一系列根目录下所有严格匹配 `shiftX.XXX` 且具有正式输出
结构的 case；二者按 `(shift_factor, symmetry_id)` 去重合并。这样多次分批计算不会因
最后一次 summary 只覆盖部分 shift 而丢点，非 shift 目录不会被纳入。

### 范围与非目标

实施范围：

- 新增 `scripts/analysis/finite_trend.py`；
- 抽取或参数化 strip/finite 共用的模式选择、CSV 校验和绘图逻辑，避免复制数值实现；
- 新增 finite trend 聚焦测试，并确保 strip trend 现有输出契约不变；
- 补充 `scripts/README.md` 中的离线分析调用示例；
- 实施完成后在本文档回填实际改动、测试结果和最终状态。

本轮不包括：

- 启动 COMSOL、补算缺失 shift 或重新生成 `mode_scores.csv`；
- 修改 `scripts/parameter.json`；
- 自动把 trend 接入 `run_finite.py` 或 `run_finite_quarter.py` 的计算尾部；
- 实现跨 shift 的 overlap/band tracking；
- 迁移旧 finite 结果目录。

Q 图必须使用 Gamma py-mode 权重筛选后、按目标频率选中的同一个 mode，不允许为了
获得更高 Q 或补齐 Q 而进行第二次独立选模。`q` 是趋势 CSV 的必需输出列。

### 公共逻辑复用方案

将以下与数据源无关的逻辑整理为 `scripts/analysis/` 内的轻量公共 helper，或把
`srip1d_trend.py` 中现有函数参数化后由两个入口共同调用：

- 目标频率表达式解析；
- `mode_scores.csv` 必需列及数值有效性校验；
- 可配置“权重列 + 严格阈值”的过滤和最近频率选模：strip 继续使用原判据，finite
  使用 `gamma_p_px_weight_fraction > 0.80`；
- 趋势 CSV 的数值列验证、排序和重绘；
- 双纵轴趋势图绘制，标题、权重列和输出 stem 由调用方传入；
- finite 专用的 Q–shift 单纵轴绘图。

优先采用独立的私有公共 helper，避免让 finite 分析反向依赖名称和语义均属于
strip 的 `srip1d_trend.py`。`srip1d_trend.py` 的公开函数、CLI 参数、
`strip1d_trend.csv` 和 `strip1d_trend.png` 保持兼容。

### Finite 输入解析与兼容策略

#### 系列类型识别

默认根据 `run_summary.json` 自动识别：

- 顶层含 `symmetry_ids`，或 run 含 `symmetry_id`：`finite_quarter`；
- run 的 `case` 为 `finite_cavity`：full `finite`；
- 无法可靠判断时给出明确错误，不凭目录名称猜测物理类型。

CLI 显式传入的 series 路径始终优先。未传入时，默认使用当前
`run_finite.OUT_DIR`，以保持和 strip trend 类似的便捷入口；分析历史或 quarter
结果时要求使用 `--out-root <series-dir>`。

#### Shift 恢复

shift-factor 按以下顺序获取：

1. run 中的 `cladding_inward_shift_factor`；
2. 从 run 的 `out_dir` 路径段中严格解析 `shift<数字>`；
3. 两者都不存在或二者不一致时，将该 run 记为问题并跳过，不静默选择值。

这使当前 full finite summary 与缺少显式 shift 字段的 finite-quarter summary
都可使用，同时保留对未来 summary schema 改进的兼容性。

#### Case 与 score 路径恢复

case 路径解析顺序：

1. summary 中可用的绝对或相对 `out_dir`；
2. 绝对路径已失效时，在 `run_summary.json` 所在系列目录下按相对尾部恢复；
3. score 优先读取 run 中有效 `overview_dir` 下的 `mode_scores.csv`；
4. 回退到 `<case>/10_overview/mode_scores.csv`；
5. 最后兼容旧布局 `<case>/mode_scores.csv`。

所有回退都必须落在用户指定的 series 目录或 summary 明确声明的 case 内，且只读。
缺文件、空 CSV、列缺失、无合格 py-mode 或所选模式的 Q 无效时记录 warning；若所有 case 都不可用，整体
返回失败且不生成误导性空图。

选模使用严格不等式 `gamma_p_px_weight_fraction > 0.80`：恰好等于 `0.80` 的模式不
通过门槛。筛选阶段只使用 py-mode 权重和频率；选中模式后再校验该行的 `q` 是否为
有限数值。Q 缺失不会触发改选另一个 mode，而会将该 case 记为不可用，从而保证两张图
始终描述完全相同的目标模式集合。

#### Quarter symmetry

默认正式配置只有 `SYMMETRY_IDS=[1]`，此时输出和 full finite 使用相同的单序列图。
若 summary 同时包含多个 symmetry ID，不能把同一 shift 的不同边界条件当成一条
模式分支。首版接口应要求使用 `--symmetry-id` 明确选择一个 symmetry；未指定且检测到
多个 ID 时返回清晰错误。趋势 CSV 保留 `source_kind` 和 `symmetry_id` 字段以便审计。

### 目标频率规则

构建趋势时目标频率按以下优先级决定：

1. CLI `--target-frequency`；
2. `run_summary.json` 中
   `finite_cavity_preset_config.eigenfrequency_shift`；
3. 当前 `run_finite.EIGENFREQUENCY_SHIFT`。

summary 内嵌频率优先于当前共享参数，避免用新的 `parameter.json` 配置分析历史系列。
若 CLI 值与 summary 值不同，允许显式覆盖，并在终端输出实际使用值。

从已有趋势 CSV 重绘时，若未传 CLI 值，则使用 CSV 的 `target_frequency`；CSV 中目标
频率缺失或不唯一时拒绝推断。

### CLI 与输出契约

计划保持与 strip trend 对齐的主要参数：

```powershell
uv run python scripts\analysis\finite_trend.py `
  --out-root scripts\.out\finite_cavity\<finite-or-quarter-series>

uv run python scripts\analysis\finite_trend.py `
  --selected-modes-csv <path-to-finite_trend.csv>
```

参数：

- `--out-root`：finite/full-quarter 系列目录；
- `--output-dir`：可选输出目录；默认 `<out-root>/finite_trend/`；
- `--target-frequency`：可选目标频率覆盖；
- `--symmetry-id`：多 symmetry quarter 系列的显式选择；
- `--selected-modes-csv`：跳过原始 score 汇总，直接重绘。

默认产物：

```text
<series>/finite_trend/
  finite_trend.csv
  finite_trend.png
  finite_q_trend.png
```

CSV 至少包含：

- `shift_factor`
- `folder`
- `source_kind`
- `symmetry_id`
- `mode_idx`
- `frequency`
- `target_frequency`
- `frequency_distance_to_target`
- `gamma_p_px_weight_fraction`
- `q`
- `mode_scores_path`

`gamma_k_weight_fraction`、`gamma_subspace_p` 等其他 score 列可作为审计信息保留，
但不得参与 finite 的门槛判断。输出按 `shift_factor` 排序；重复
`(shift_factor, symmetry_id)` 视为 summary 歧义并报错。

`finite_trend.png` 维持 strip trend 的视觉结构：横轴为
`CLADDING_INWARD_SHIFT / A`，左轴为 `gamma_p_px_weight_fraction`，右轴为
`frequency (THz)`，并绘制目标频率虚线。标题改为 finite/finite-quarter 语义，不写死
Strip。`finite_q_trend.png` 使用相同横轴，纵轴为所选目标模式的 `Q`，按 shift 排序后
使用折线和数据点绘制；首版使用线性纵轴，不对 Q 做对数变换或归一化。

### 数据流

```text
series/run_summary.json
  -> 识别 full finite / finite-quarter
  -> 解析 run、shift、symmetry 与 case 路径
  -> case/10_overview/mode_scores.csv
  -> 校验列和数值
  -> gamma_p_px_weight_fraction > 0.80
  -> 选择最接近目标频率的 mode
  -> 校验该 mode 的 q
  -> finite_trend/finite_trend.csv
  -> finite_trend/finite_trend.png + finite_q_trend.png
```

### 测试计划

新增 `tests/test_analyze_finite_trend.py`，使用临时目录和小型 CSV/JSON fixture，至少覆盖：

1. full finite 从 `10_overview/mode_scores.csv` 读取并正确选模；
2. finite-quarter 从 `shiftX.XXX` 路径恢复缺失的 shift-factor；
3. summary 中相对路径、失效绝对路径和 `overview_dir` 的恢复；
4. 合并 summary runs 与系列根目录中的全部合法 `shiftX.XXX` case，不遗漏分批计算结果，
   同时忽略 `finite_trend/` 等非 shift 目录；
5. `gamma_p_px_weight_fraction > 0.80` 才能入选，并明确验证 `0.80` 边界值被拒绝；
6. 即使 `gamma_subspace_p` 或 `gamma_k_weight_fraction` 很高，py-mode 权重不超过 80%
   的模式仍被拒绝；
7. 缺列、缺文件、无合格 py-mode、所选模式 Q 缺失或非数值时的错误/告警；
8. 等频率距离时按 `mode_idx` 稳定打破平局；
9. 目标频率优先级：CLI、summary、当前配置；
10. 单 symmetry quarter 正常生成，多 symmetry 未指定 ID 时拒绝混合，指定后正确筛选；
11. 趋势 CSV 的列、排序、必需 Q 和来源审计字段；
12. Q 图与 Gamma/frequency 图使用完全相同的 `mode_idx` 行，不独立追逐最大 Q；
13. 从已有 CSV 同时重绘两张图，且不重新读取 raw `mode_scores.csv`；
14. 输出路径固定为 `finite_trend/finite_trend.csv`、`finite_trend.png` 和
    `finite_q_trend.png`；
15. `srip1d_trend.py` 现有测试全部通过，确认可配置公共 helper 没有改变 strip 原判据。

实施后的非 COMSOL 验证命令：

```powershell
uv run python -m pytest -q tests\test_analyze_finite_trend.py `
  tests\test_analyze_strip_gamma_trend.py `
  tests\test_analyze_strip_gamma_summary.py
uv run python -m py_compile scripts\analysis\finite_trend.py `
  scripts\analysis\srip1d_trend.py
git diff --check
```

在单元测试通过后，可对当前 finite-quarter 样例执行一次只读 smoke test，输出到新的
`finite_trend/` 子目录，并确认两张 PNG 使用同一组目标 mode；该步骤不读取 MPH、不启动 COMSOL。若正式结果目录写入需要单独
确认，则只在临时 fixture 上完成验证。

### 回滚方案

回滚仅涉及删除本任务新增的分析脚本、公共 helper、测试和文档变更，并恢复
`srip1d_trend.py` 的原有 helper 实现。由于分析只读已有 case，回滚不涉及任何
COMSOL 结果、MPH、parquet 或 summary 的恢复。

## 执行完成后记录

### 实际改动

1. 新增 `scripts/analysis/_gamma_trend_common.py`：
   - 提供共享的频率解析、可配置权重门槛选模、趋势 CSV 加载和双纵轴绘图；
   - 数值筛选同时拒绝 NaN 和无穷值；
   - 门槛使用严格大于关系，并以频率距离、`mode_idx` 的顺序稳定排序。
2. 新增 `scripts/analysis/finite_trend.py`：
   - 自动识别 full finite 与 finite-quarter summary；
   - 以 `gamma_p_px_weight_fraction > 0.80` 筛选 py-mode，再选择最接近目标频率的模式；
   - 支持 quarter 从 `shiftX.XXX` 路径恢复 shift，并支持失效绝对路径的 series 内恢复；
   - 多 symmetry summary 要求显式传入 `--symmetry-id`；
   - 目标频率优先级为 CLI、summary 内嵌配置、当前 finite 配置；
   - 输出 `finite_trend.csv`、`finite_trend.png` 和 `finite_q_trend.png`；
   - Q 无效时不会改选其他 mode，确保两张图始终对应同一目标模式集合；
   - 支持 `--selected-modes-csv` 同时重绘两张图。
3. `srip1d_trend.py` 改为复用公共 helper，但保留原有
   `gamma_subspace_p > 0.99` 判据、公开函数、CLI 和产物命名。
4. 新增 `tests/test_analyze_finite_trend.py`，覆盖 80% 严格边界、其他 Gamma 指标不参与
   finite 筛选、full/quarter 路径、summary 范围、多 symmetry、Q 一致性、目标频率优先级
   和双图输出。
5. 更新 `scripts/README.md`，加入 finite trend 的用途、命令和多 symmetry 使用说明。
6. 未修改用户已有的 `scripts/parameter.json`，未启动 COMSOL，也未读取或重写 MPH。
7. 2026-07-26 补充修复：发现正式 quarter 系列有 10 个 shift case，而最近一次
   `run_summary.json` 只记录其中 5 个。分析索引改为“summary runs + 系列根目录合法
   shift case”合并，summary 仍作为物理元数据来源，目录扫描只补齐遗漏 case。

### 测试结果

- 聚焦 trend 测试：`31 passed`。
- trend 加 finite 输出相关回归：`40 passed`。
- 全仓库测试：`233 passed, 3 warnings`；3 条均为已有
  `PytestReturnNotNoneWarning`，位于 `tests/test_strip_1d_fourier_quotient.py`。
- `py_compile`：`_gamma_trend_common.py`、`finite_trend.py` 和
  `srip1d_trend.py` 全部通过。
- `git diff --check` 通过；仅显示工作区既有 CRLF 提示，无空白错误。
- 使用现有
  `finite_quarter_10-10_cav(245-0.96-1.155)_clad(243.7-0.98-0.928)_mesh9`
  系列完成离线 smoke test：
  - 目标频率来自 summary，为 `198.424452450838 THz`；
  - 初次发现 summary 只声明 5 个 shift，因而只生成 5 点；修复索引策略后，从系列目录
    补齐 0.00–0.09 共 10 个 shift；
  - 0.00 和 0.05、0.07 选中 mode 3，其余 shift 选中 mode 1；
  - 所选 `gamma_p_px_weight_fraction` 均严格高于 0.80；
  - CSV、Gamma/frequency 图和 Q-shift 图均成功生成并完成视觉检查；
  - smoke test 未启动 COMSOL。

### 遗留问题与最终状态

实现与验证均已完成。首版按审核方案保持为独立离线分析入口，尚未自动接入
`run_finite.py` 或 `run_finite_quarter.py` 的计算尾部；若未来需要多 symmetry 同图比较、
Q 对数坐标或 overlap-based 跨 shift band tracking，应作为独立需求扩展，不能改变当前
py-mode 80% 门槛和同一 mode 双图契约。
