# Unit-cell 2D C2v 第一象限扫描与绘图对称映射计划

- 日期：2026-08-26
- 状态：**代码实施与纯 Python 验证完成；实机 C2v 对照验证待运行**
- 目标入口：`scripts/run_main/run_unit_cell_2d.py`
- 参数源：`scripts/parameter.json` 的 `unit_cell_2d`
- 输出根目录：`scripts/.out/unit_cell_2D/`
- 本文同时记录已批准的重构合同和实际执行结果；本轮未迁移既有结果，也未启动 COMSOL。

## 一、目标与边界

当前 unit-cell 结构及目标模式在研究窗口内满足 C2v 对称性。因此标准计算不再独立求解
四个象限，而是：

1. 只在第一象限 $k_x\ge 0,\ k_y\ge 0$ 求解真实 COMSOL 本征场和 $c(\mathbf k)$；
2. 只在绘图和闭合回路求值时使用 C2v 镜像映射，不生成、不保存其他象限的数据点；
3. 通过只读的对称场访问器生成现有 A--D 全部标准图、偏振椭圆和 winding number；
4. 保留 `isQuarter=0` 的完整 k 空间独立求解选项，作为兼容、数值对照和对称性验证入口；
5. 输入方式仍然只有 k 空间边界 `q_max_over_G` 和一个方向上的取点数
   `q_points_per_axis`，不恢复人工坐标列表。

本次重构不改变几何、mesh、四本征解合同、Gamma 模式识别、Hz-overlap band tracking、
air-plane 投影定义、B 图共同归一化、绘图风格和既有输出一级目录。

## 二、参数合同

在 `scripts/parameter.json` 的 `unit_cell_2d` 中增加：

```json
"isQuarter": 1
```

严格只接受整数 `0` 或 `1`，不接受字符串、浮点数或其他整数。

### 2.1 `isQuarter=1`：新的标准模式

设

```text
q_max_over_G = K
q_points_per_axis = N
```

则 COMSOL 真实求解轴为

$$
q_i/G=\frac{iK}{N-1},\qquad i=0,1,\ldots,N-1,
$$

即 `linspace(0, K, N)`。两方向各有 $N$ 点，共求解并保存 $N^2$ 个 k 点。绘图时使用的
完整坐标轴为

$$
[-K,-K+\Delta q,\ldots,0,\ldots,K-\Delta q,K],
\qquad \Delta q=\frac{K}{N-1},
$$

共 $2N-1$ 个坐标。这里的完整坐标轴只属于绘图视图，不据此创建
$(2N-1)^2$ 行结果表。

因此 `kmax0.15_uniform25_quarter` 的含义固定为：

- 第一象限每个方向25点；
- COMSOL求解 $25^2=625$ 点；
- 间隔为 $0.15/24=0.00625G$；
- 绘图分辨率对应 $49\times49$ 个坐标位置，但磁盘中仍只保存625个真实求解点。

`N` 必须是整数且 $N\ge2$；Gamma 和 $K$ 都由端点明确包含，因此不再要求奇数。

### 2.2 `isQuarter=0`：保留的完整扫描模式

为保持旧工作流和旧参数快照的语义，真实求解轴继续使用

```text
linspace(-K, K, N)
```

并真实求解 $N^2$ 点，不生成镜像行。此时 `N` 必须是大于等于3的奇数，保证Gamma唯一。
例如 `kmax0.15_uniform25_full` 仍是原来的25×25、625点完整扫描，间隔为0.0125G。

同一个 `N` 在 quarter/full 下表示不同的完整网格分辨率；这是有意保留的兼容合同。若要
以相同 k 间隔验证 `quarter25`，应使用 `full49`，而不是 `full25`。

### 2.3 缺失字段的兼容策略

- 当前 `scripts/parameter.json` 显式增加 `"isQuarter": 1`；
- 历史参数快照中若缺失 `isQuarter`，解析为 legacy full，即 `isQuarter=0`，不得把旧
  `uniform25` 静默解释成第一象限25点；
- `q_axis_over_G` 继续拒绝；
- `unit_cell_2d_metadata()` 必须写出输入模式、真实求解轴、绘图轴、真实求解点数、绘图
  坐标数和间隔，并明确后者不是 COMSOL 数据点数。

## 三、C2v 绘图映射的数学合同

定义镜像操作

$$
M_x:(k_x,k_y)\mapsto(-k_x,k_y),\qquad
D(M_x)=\begin{pmatrix}-1&0\\0&1\end{pmatrix},
$$

$$
M_y:(k_x,k_y)\mapsto(k_x,-k_y),\qquad
D(M_y)=\begin{pmatrix}1&0\\0&-1\end{pmatrix},
$$

以及 $C_2=M_xM_y$，其二维极向量表示为 $D(C_2)=-I$。

对同一条非简并、连续追踪的 Bloch band，Jones 向量满足

$$
\mathbf c(g\mathbf k)
=e^{i\phi_{g,n}(\mathbf k)}D(g)\mathbf c(\mathbf k),
\qquad g\in\{E,M_x,M_y,C_2\}.
$$

$e^{i\phi_{g,n}}$ 是 band 的公共 sewing phase；它同时乘在 $c_x,c_y$ 上，不影响强度、
Stokes 参数、偏振椭圆或 formal winding。绘图访问器采用明确的
`symmetry_gauge_v1`，令镜像求值时的公共 sewing phase 为1：

| 目标点 | `cx` | `cy` |
|---|---:|---:|
| $(+k_x,+k_y)$ | $c_x$ | $c_y$ |
| $(-k_x,+k_y)$ | $-c_x$ | $+c_y$ |
| $(+k_x,-k_y)$ | $+c_x$ | $-c_y$ |
| $(-k_x,-k_y)$ | $-c_x$ | $-c_y$ |

该规则只在内存中对 `cx_raw/cy_raw` 和 `cx_norm/cy_norm` 求值，并直接交给当前绘图或
winding 计算，不把返回值追加到 DataFrame、CSV、Parquet、point cache 或 COMSOL 目录。
频率、Q、`normalization_value` 等标量图若需要完整画幅，则在绘图函数中通过
$f(k_x,k_y)=f(|k_x|,|k_y|)$ 读取同一个真实源点，也不复制记录。

在当前 Stokes 定义

$$
S_0=|c_x|^2+|c_y|^2,\quad
S_1=|c_x|^2-|c_y|^2,\quad
S_2=2\operatorname{Re}(c_xc_y^*),\quad
S_3=-2\operatorname{Im}(c_xc_y^*)
$$

下，单次镜像保持 $S_0,S_1$，并使 $S_2,S_3$ 变号；$C_2$ 保持全部 Stokes 参数。
偏振椭圆的位置、长短轴和箭头手性必须由绘图时映射得到的 Jones 向量重新计算，不能
复制 PNG 或只复制角度。

### 3.1 坐标轴、Gamma 与只读访问器

对任一绘图目标坐标 $(k_x,k_y)$，访问器执行：

1. 将其规范化到真实源坐标 $(|k_x|,|k_y|)$；
2. 以现有6位小数坐标身份查找唯一的第一象限真实记录；
3. 根据目标坐标所在象限，仅对本次返回的 Jones 向量应用 `E/Mx/My/C2`；
4. 返回临时标量或矢量值，调用结束后不保留镜像记录。

Gamma 及两条轴仍只对应磁盘中的一个真实源记录。访问器必须检查源轴严格递增、Cartesian
product 完整、无重复、无非法负零；目标轴也必须唯一且严格递增。任何源点缺失都直接
阻止绘图，不能用插值或生成行补齐。

## 四、数据来源与可追溯性

新结果只保留一层权威点数据：`polarization_grid` 中每一行都必须对应真实 COMSOL 求解点
或同一 identity 下恢复的真实缓存点。quarter 模式下它严格只有 $N^2$ 行，且坐标均满足
$k_x\ge0,k_y\ge0$；full 模式下它有 $N^2$ 行完整空间真实数据。不得写出
`symmetry_expanded` 行，也不得另存完整镜像表。

表和元数据至少记录：

- `is_comsol_solved=true`；
- `sampling_domain=quarter/full`；
- `formal_band_label`：用户名称 `px/py`，内部 `p1/p2` 只保留在兼容元数据中；
- 真实的 `gamma/matched/repaired` 状态及 band-tracking 证据。

点目录、`point_metadata.json`、Hz parquet 和 MPH 也只为这些真实求解点创建。绘图访问器
版本、C2v 规则和 gauge 版本写入配置与绘图摘要，但不写逐镜像点 provenance；因为这些
点从未成为结果数据。

## 五、求解顺序、band tracking 与缓存

### 5.1 Quarter 模式

- Gamma仍为第一个点；
- 其余第一象限点继续按方形半径、欧氏半径和坐标确定性排序；
- Hz-overlap predecessor 只能引用已真实求解的第一象限状态；
- 绘图对称映射只在全部真实点完成模式匹配和偏振提取之后启用，不参与 COMSOL band
  tracking；
- 现有一次 repair pass 只修复真实第一象限点，之后直接更新权威真实点表。

### 5.2 Full 模式

保留现有完整方形径向调度、Gamma锚定、多前驱匹配和repair行为；绘图直接使用真实完整
表，不经过对称访问器。

### 5.3 Identity 与恢复

以下字段必须进入 cache identity：

- `isQuarter`；
- `solved_q_axis_over_G` 与 `analysis_q_axis_over_G`；
- 真实求解点数与绘图坐标数；
- 新的 scan definition version；
- `symmetry_group=C2v`；
- `symmetry_render_version` 与 `symmetry_gauge_version`。

quarter/full 必须使用不同 series 名称，禁止共享同一目录。历史 full 结果保持只读，不移动、
不覆盖，也不因部分坐标重合而跨 identity 复用缓存。

## 六、分析与绘图合同

标准 Analysis 统一读取只含真实求解点的 `polarization_grid`。quarter 模式创建一个只读
`C2vQuarterFieldView`；full 模式创建直接读取真实完整网格的 `FullFieldView`。两者向绘图
函数提供相同的按坐标求值接口，因此现有 A--D 图的画幅、标题、colorbar、归一化和输出
文件名保持不变：

- A：完整 k 空间的 $|c(\mathbf k)|$ linear/log；
- B：完整 k 空间的 $|c_x|,|c_y|$ linear/log；
- C：相位/偏振椭圆等现有输出；
- D：完整闭合回路的偏振椭圆与 BIC winding。

B 图在真实求解表上计算同一个

$$
C_{\max}=\max_{\mathbf k}\sqrt{|c_x^{\rm raw}|^2+|c_y^{\rm raw}|^2}
$$

作为共同分母。由于镜像不改变模长，quarter 真实表的最大值与完整画幅最大值严格相同。

`ordered_square_loop()` 和 shrinking-loop 扫描通过 `C2vQuarterFieldView` 在完整闭合路径
上逐点求值，不能只在第一象限拼接一个非闭合路径，也不能为该回路落盘镜像点表。正式
winding 继续以 `stokes-axis`/偏振椭圆为物理输出；
`real/imag` 属于所选 symmetry gauge 下的诊断结果，不作为跨 gauge 的物理结论。

绘图前增加真实源网格门禁：每个 band 必须恰有 $N^2$ 行，quarter 的源轴为 $[0,K]$，
full 的源轴为 $[-K,K]$，且 Cartesian product 完整、每个坐标唯一。quarter 绘图访问器还
必须能覆盖 $(2N-1)^2$ 个目标坐标；这只是访问器覆盖测试，不产生对应行。未通过时拒绝
生成 A--D，不能用插值填补缺失源点。

## 七、输出目录与命名

新结果仍与其他 unit-cell 2D 结果并列放在 `scripts/.out/unit_cell_2D/`。建议稳定命名：

```text
unit_cell_2D_(<b0>-<eta>-<zeta>)_mesh<mesh>_kmax<K>_uniform<N>_quarter_band-<formal>
unit_cell_2D_(<b0>-<eta>-<zeta>)_mesh<mesh>_kmax<K>_uniform<N>_full_band-<formal>
```

其中 `uniform<N>` 始终记录输入的 `q_points_per_axis`，`quarter/full` 明确其解释；用户可见
mode 使用 `px/py`，不得只以内部 `p1/p2` 命名新 series 或图。

目录内容：

```text
00_model/
  geometry.png
  unit_cell_2D.mph
01_results/
  k_points/                         # 仅真实求解点
  polarization_grid.parquet         # 仅真实求解点
10_overview/
  A--D 现有标准图
80_logs/
  polarization_grid.csv              # 仅真实求解点
  band_tracking_2d.csv               # 仅真实点的tracking证据
99_config/
  config.json
  run_summary.json
  winding_summary.json
  symmetry_render_summary.json
```

所有 CSV 继续放入 `80_logs`。`symmetry_render_summary.json` 只记录真实点数、绘图坐标
数、访问器版本、目标轴、完整性检查和源数据哈希，不包含逐镜像点记录。

## 八、zeta 扫描、独立 Analysis 与 valley refinement

### 8.1 zeta 扫描

`scripts/run_sweep/run_unit_cell_2d_zeta_scan.py` 必须把 `isQuarter` 写入每个 task-local
参数快照和系列摘要；quarter/full case 不得汇入同一系列目录。

### 8.2 独立重绘

`scripts/analysis/unit_cell_2D.py` 按参数解析 quarter/full series 名；若直接给定结果目录，
则以该目录 `99_config/config.json` 为权威，不从当前 `parameter.json` 猜测扫描类型。

### 8.3 COMSOL valley refinement

现有 refinement loader 假定输入表覆盖完整空间。quarter 新结果只包含第一象限真实点，
重构必须同步执行以下兼容修改：

- valley 几何检测通过同一个只读 C2v 场访问器观察完整画幅，不创建完整扩展表；
- 所有候选求解坐标先规范化为第一象限 $(|k_x|,|k_y|)$ 并按6位坐标去重；
- band tracking 只从 `polarization_grid` 及真实第一象限 Hz cache 恢复状态；
- 新增局部 COMSOL 点求解后追加为真实第一象限点，完整 B4 仍在绘图时按 C2v 映射；
- refinement summary 只报告 canonical solved 点和绘图覆盖范围，不报告生成点；
- 在这部分适配完成前，quarter source 必须被旧 refinement 入口明确拒绝，不能静默把
  目标象限坐标当成真实 Hz 文件路径。

既有709点 legacy full refinement 结果保持原状，不迁移到新 quarter series。

## 九、实现范围

计划批准后修改：

- `scripts/parameter.json`
- `scripts/run_main/parameter_config.py`
- `scripts/run_main/run_unit_cell_2d.py`
- `comsol_workflow/unit_cell_2d_analysis.py`
- `scripts/analysis/unit_cell_2D.py`
- `scripts/run_sweep/run_unit_cell_2d_zeta_scan.py`
- `comsol_workflow/unit_cell_2d_valley_refinement.py`
- `scripts/run_sweep/unit_cell_2d_valley_solver.py`
- 相关测试、README、长期 Memory 和本 spec 的执行回填

不修改 unit-cell band、strip、finite/finite-quarter 的扫描参数和求解语义。

## 十、验证计划

### 10.1 纯 Python 单元测试

1. `isQuarter` 只接受0/1；缺失字段回退 legacy full；
2. quarter `N=3` 只保存9个真实点，访问器覆盖5×5目标坐标且不改变表行数；
3. full `N=3` 保留3×3=9个真实点且不调用对称访问器；
4. quarter允许偶数N，full拒绝偶数N；
5. Gamma、两条轴和内部目标坐标分别正确映射到唯一真实源点；
6. `Mx^2=My^2=C2^2=E`，且 `MxMy=MyMx=C2`；
7. $S_0,S_1$ 镜像不变，$S_2,S_3$ 单镜像变号；
8. raw/norm Jones 向量、频率、Q和normalization的临时映射正确且不落盘；
9. quarter/full名称、identity、缓存互不冲突；
10. 第一象限真实表缺点、重复、冲突或来源字段缺失时 Analysis 拒绝；
11. winding loop 使用只读访问器完成完整闭合求值，且不生成镜像数据表；
12. valley refinement 只把 canonical 第一象限坐标交给 COMSOL。

### 10.2 小网格 COMSOL 对称性验证

在批准后先使用小 N：

1. 运行 quarter 小网格；
2. 使用 `isQuarter=0` 在对应代表点独立求解；
3. 比较频率、$1/Q$、$|c_x|,|c_y|$、归一化 Stokes 和 Jones 向量的最佳公共相位对齐误差；
4. 覆盖一个内部点、两条坐标轴和Gamma附近点；
5. 将实际最大偏差写入 `symmetry_validation.json`，通过后才允许正式 quarter25。

raw 本征矢的任意整体复尺度不能直接逐复数比较；验证必须使用 planar-L2 Jones 向量并
消除一个公共复相位，或使用 gauge-invariant Stokes/模长。

### 10.3 回归命令

```powershell
uv run python -m pytest -q tests\test_run_unit_cell_2d.py
uv run python -m pytest -q tests\test_unit_cell_2d_analysis.py
uv run python -m pytest -q tests\test_run_unit_cell_2d_zeta_scan.py
uv run python -m pytest -q tests\test_unit_cell_2d_valley_refinement.py
uv run python -m py_compile scripts\run_main\run_unit_cell_2d.py
uv run python -m py_compile scripts\analysis\unit_cell_2D.py
git diff --check
```

## 十一、完成门禁

只有同时满足以下条件才把 quarter run 标为 `complete`：

- $N^2$ 个第一象限真实点全部完成或由同 identity cache 恢复；
- 每个目标 band 在真实点上均通过模式匹配门禁；
- 第一象限真实表点数、Cartesian product、坐标和 provenance 完整；
- C2v 绘图访问器覆盖完整目标轴，没有修改源表或写出镜像点；
- A--D 分析只读取真实表并通过访问器全部成功；
- winding loop 闭合且摘要来自正式 `stokes-axis`；
- `run_summary.json` 分开记录真实求解点数和绘图目标坐标数，并注明后者不是数据点；
- 实际 PNG 完成视觉检查；
- 聚焦测试、语法检查和 `git diff --check` 通过。

任一条件失败时保持 `incomplete`，不得绕过访问器把第一象限表当作完整空间，也不得用
插值、复制 PNG 或镜像数据行补齐结果。

## 十二、回滚策略

- `isQuarter=0` 保留完整扫描实现，可作为即时计算回退；
- legacy full 参数快照继续按旧语义解析；
- 新 series 使用独立 `_quarter/_full` 后缀，回滚代码不会覆盖旧结果；
- 对称绘图访问发生在真实点求解完成后，失败时第一象限 COMSOL cache 仍可保留并恢复；
- 不移动、不改写任何既有 `uniform25`、历史 B4 或709点 refinement 结果。

## 十三、待审核的关键解释

本计划采用以下解释，需在实施前由用户确认：

1. `isQuarter=1, uniform25` 表示第一象限每方向25个真实点；绘图覆盖49×49个坐标，但
   始终只保存25×25个真实数据点；
2. `isQuarter=0, uniform25` 保留历史完整25×25语义，而不是完整49×49；
3. 新目录必须带 `_quarter` 或 `_full`，避免相同 `kmax/uniformN` 发生缓存碰撞；
4. C2v 只在绘图/闭合回路求值时对 complex Jones 向量采用显式
   `symmetry_gauge_v1`，不生成或保存镜像数据点；正式物理结论使用不依赖公共 sewing
   phase 的模长、Stokes、偏振椭圆与 winding；
5. valley refinement 同步改成第一象限 canonical 求解，其他象限只通过同一只读访问器
   显示。

## 十四、实际执行记录

### 14.1 已完成改动

- 在共享参数中加入严格的 `isQuarter=0/1`，当前默认值为 `1`；历史快照缺失该字段时仍按
  legacy full 解释；
- quarter 模式只构造、求解和保存第一象限的 $N\times N$ 个真实点，full 模式保留原有
  $[-K,K]$ 的 $N\times N$ 独立求解语义；
- 新增只读 `C2vQuarterFieldView`，在 A--D 绘图和闭合 winding 回路求值时临时完成
  `Mx/My/C2` 映射，不向 CSV、Parquet、point cache、MPH 目录或 DataFrame 写入镜像点；
- series、cache identity、配置快照和运行摘要均区分 quarter/full，并分开记录真实点数与
  绘图坐标数；用户可见模式名使用 `py/px`；
- zeta 扫描与 valley refinement 已同步 quarter 合同；refinement 的 COMSOL 候选坐标先
  规范化到第一象限，完整画幅只在检测与绘图访问时形成；
- minimum-curvature、topology reconstruction 和历史 B5 实验入口不会接受 quarter 数据；
- README、Memory 和聚焦测试已同步。
- quarter-series 的标准 B4 保持原始数据、共同归一化和共享平滑合同，仅将 logarithmic
  colorbar 固定为 $10^{-4}$ 到 $10^{-1}$；full-series 继续使用数据派生范围。
- quarter-series 的标准 B2 同样保持数据与绘图合同，仅将 logarithmic colorbar 固定为
  $10^{-4}$ 到 $10^0$；full-series B2 继续使用数据派生范围。
- quarter-series 的标准 B1 保持数据与绘图合同，仅将 linear colorbar 固定为 $0$ 到
  $1$；full-series B1 继续使用数据派生范围。
- quarter-series 的 A1/A2 改为绘制归一化
  $|c(k)|=\sqrt{|c_x(k)|^2+|c_y(k)|^2}$；为保持输出兼容，文件名仍保留既有
  `A1_ck_intensity_linear.png` 和 `A2_ck_intensity_log.png`。full-series A1/A2
  保留历史归一化 intensity 定义。
- quarter-series 的 D 图只保留未缩放的 winding-loop 偏振椭圆与对应 winding-number
  曲线，删除历史 `c_y\times10` 的两幅下排面板；full-series D 保留兼容的两行布局。
- quarter-series 的全部标准 kx--ky 面板（A1/A2、B1--B4、C1/C2、D 左图）明确标出
  坐标最小值和最大值，并以 $0.05G$ 为主刻度间隔；D 右图等非 k 空间坐标不修改，
  full-series 继续使用数据派生刻度。
- quarter-series 的正负坐标刻度统一去掉无意义的尾随零，例如显示 `0`、`0.1`、
  `-0.1`、`0.05` 和 `-0.05`，不显示 `0.00`、`0.10` 或 `-0.10`。

### 14.2 验证结果

- unit-cell 2D 聚焦测试：`91 passed`；
- 完整项目回归测试：`517 passed, 3 warnings`；三条 warning 均来自既有 strip 测试返回
  非 `None`，与本次改动无关；
- 相关 Python 文件 `py_compile` 通过；
- 合成 quarter 数据验证：3×3 真实源表保持 9 行，只读视图覆盖 5×5 坐标；Jones/Stokes
  奇偶性正确，闭合 complex-axis winding 为 `+1`；
- 已检查合成数据生成的 A2、C1、D PNG：完整四象限、等比例坐标、colorbar 和面板布局
  均正常，且渲染后源表行数未变化；
- `git diff --check` 通过，仅报告工作树既有的 LF/CRLF 转换提示。

### 14.3 未执行与首次正式运行门禁

本轮没有启动 COMSOL，也没有生成或覆盖正式计算结果。首次正式 quarter25 计算前，仍按
第 10.2 节执行小网格 quarter/full 实机对照，并将频率、$1/Q$、Jones/Stokes 的实际最大
偏差写入 `symmetry_validation.json`。因此当前状态表示“程序重构完成”，不表示新的
quarter COMSOL 数据已经求解完成。

### 14.4 C2 相位图坐标轴修正（2026-08-27）

- 以同参数 full-space `uniform25` 的 `C2_ck_phase.png` 和原始 Jones/Stokes 数据为
  对照，确认 quarter 与 full 的第一象限在重合采样点上一致；相位差来自对称绘图，
  不来自 COMSOL 求解、模式选择或 cyclic colormap。
- 象限内部继续遵守单次镜像使 polarization director angle 变号、两次镜像恢复原值的
  C2v 规则；`kx=0` 与 `ky=0` 上不再套用内点变号，而是沿对应正半轴复制同一 director
  branch。这样不生成镜像数据行，也不会让轴上 branch choice 经 doubled-angle 插值扩散
  成错误的十字相位拼接。
- ζ=1.156 的 625 个 full/quarter 重合坐标上，修正后 director-angle 的平均周期偏差为
  `1.35e-3 rad`；其中象限内平均偏差为 `1.21e-3 rad`，坐标轴平均偏差为
  `3.04e-3 rad`。偏差来自独立 full 求解的数值误差，而不是额外拟合。
- 修正严格限定于 quarter C2 的 `psi_stokes` 显示矩阵；Jones 向量、C1、D、winding、
  full-series 和源 CSV/Parquet 均不修改，也不启动 COMSOL。
- 已原位重绘 ζ=`1.154、1.155、1.156、1.157` 四组正式 quarter C2；实际 PNG 的标题、
  `twilight_shifted` cyclic colorbar、$[-\pi/2,\pi/2]$ 范围、等比例坐标和 $0.05G$
  刻度均通过视觉检查。
- 聚焦回归结果：`tests/test_unit_cell_2d_analysis.py` 为 `44 passed`；相关 `py_compile`
  与 `git diff --check` 通过。

### 14.5 C1 偏振椭圆密度与旋向配色（2026-08-27）

绘图合同：

- 目标文件仍为 `10_overview/C1_ck_vector.png`，标题仍为
  `Polarization ellipses (py mode)`；坐标范围、等比例画幅和 $0.05G$ 刻度不变。
- quarter25 的只读 49×49 对称视图仅在绘图时均匀抽取 25×25 个位置；每个位置画一个
  带箭头的真实 Jones 偏振椭圆。抽样不创建或保存镜像 CSV/Parquet 行。
- 使用当前定义 $S_3=-2\operatorname{Im}(c_xc_y^*)$ 判定旋向：$S_3>0$ 对应椭圆参数
  正向的逆时针旋转并使用红色，$S_3<0$ 对应顺时针旋转并使用蓝色；箭头继续作为非颜色
  的冗余方向编码。精确 $S_3=0$ 使用中性灰色。
- 修改只作用于 quarter C1；full C1、C2、D、winding、Jones 数据与 COMSOL 求解均保持
  不变。

执行结果：

- 已原位重绘 ζ=`1.154、1.155、1.156、1.157` 四组正式 Quarter
  `C1_ck_vector.png`，每张均绘制 25×25=625 个带箭头椭圆；实际 PNG 已完成视觉检查，
  图例位于图框上方且不遮挡数据。
- `tests/test_unit_cell_2d_analysis.py`：`45 passed`；相关 `py_compile` 与
  `git diff --check` 通过。未启动 COMSOL，也未修改任何源 CSV/Parquet。
