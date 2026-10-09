# unit_cell_2D low-value valley 的 COMSOL 局部加密重算 Spec

## 文档状态

- 状态：**实施中；ζ=1.156 采样计划已确认，正式 COMSOL 求解运行中**。
- 已完成参数解析、纯数据 valley 导航、首轮采样计划、可恢复 arbitrary-k 求解编排、
  真实点合并和 Delaunay 分片线性绘图实现。
- 已创建独立的新结果目录和采样预览；未启动 COMSOL，未覆盖任何既有 B4、CSV、JSON
  或 MPH。
- 用户已确认当前 189 点首轮采样计划，并明确允许与现有 COMSOL GUI 并发求解。

## 一、问题与结论

现有 `B4_cy_magnitude_log_valley_reconstruction.png` 通过 low-value 检测、拓扑
分类、支线拟合、distance suppression 和人为 line floor 改写显示场。即使所有步骤
可审计，它仍然是在稀疏数据上加入结构先验后的图像重建，不是对真实连续场或精确
nodal line 的独立求解。该路线不能作为新的物理结果继续使用。

新路线遵循以下边界：

1. 已有 uniform25 数据只用于估计 **应该在哪里增加 COMSOL 采样点**；
2. valley 的最终幅值只能来自已有 COMSOL 点或新求解的 COMSOL 点；
3. 低谷位置函数、切向、法向和局部加密网格不参与修改任何物理幅值；
4. 最终图只允许对真实采样值做分片线性插值和 Matplotlib 渲染；
5. 禁止 conic、固定 X/hyperbola 拓扑、五阶支线、minimum-curvature、Gaussian、
   人为 floor、中心恢复、支线压暗、幅值拟合或其他后验图像修饰；
6. 如果局部 COMSOL 采样仍不足以解析 valley，则结果标记为未收敛并继续加密或停止，
   不用数值修图掩盖采样不足。

“真实物理数值”在本文中指：在指定 COMSOL mesh、eigensolver、场平面和模式跟踪
合同下直接求得的离散数值。它不等价于数学上的无限精度连续解；因此仍要保留 mesh、
模式跟踪和局部采样收敛证据。

## 二、目标与非目标

### 目标

- 从已有 uniform25 的 py 模式数据自动提取 low-value valley 的分段位置函数；
- 围绕位置函数构造局部、可恢复、无重复的 arbitrary-k 加密点集合；
- 复用一个固定几何、固定 mesh 的 COMSOL model 逐点求解；
- 用 Hz 场重叠和频率连续性恢复同一 py 模式；
- 保存每个新点的本征频率、Q、中心面 Hz、`cx_norm/cy_norm` 和跟踪证据；
- 合并原 uniform25 与新 COMSOL 点，用分片线性插值生成新的 B4；
- 输出完整采样计划、运行清单、缓存、日志、收敛证据和可复现配置。

### 非目标

- 不从 PNG 反推数据；
- 不声称位置函数本身是精确 nodal line；
- 不用位置函数、拟合式或先验模型生成最终幅值；
- 不重新计算完整 uniform49 或更大的全平面均匀网格；
- 不覆盖标准 `B4_cy_magnitude_log.png`、历史参考图或被否决的历史输出；
- 不修改 unit-cell band、strip 或 finite-cavity 工作流。

## 三、入口与模块边界

用户要求的完整入口拟定为：

```text
scripts/analysis/unit_cell_2D_valley_refine.py
```

该脚本提供一个完整、可恢复的命令行工作流：

```text
读取源 series
  -> 验证 coarse 数据与 solver identity
  -> 提取 valley 位置函数
  -> 生成/审核局部加密点
  -> 调用 COMSOL 求解新点
  -> 合并真实采样
  -> 线性插值并绘制新 B4
```

虽然入口按要求位于 `scripts/analysis/`，但 COMSOL 底层实现不复制到该脚本。计划将
arbitrary-k 点求解复用为 `scripts/run_sweep/` 或 `comsol_workflow/` 中的公共函数，入口
只负责阶段编排。这样避免复制几何、mesh、模式跟踪、缓存和场导出实现。

拟新增 pyproject 命令：

```toml
comsol-unit-cell-2d-valley-refine = "scripts.analysis.unit_cell_2D_valley_refine:main"
```

默认命令执行完整的“prepare -> solve -> plot”流程；同时提供：

- `--prepare-only`：只生成采样计划，不启动 COMSOL；
- `--plot-only`：只在已完成求解后重新合并、绘图；
- `--resume`：默认启用，复用逐点缓存；
- `--allow-concurrent-comsol`：显式允许与另一个本项目 COMSOL 工作流并发；默认拒绝。

## 四、输入数据与物理标量

### 4.1 源 series

必选位置参数：

```powershell
uv run python scripts\analysis\unit_cell_2D_valley_refine.py `
  <source-series-dir> --prepare-only
```

源目录必须包含：

- `99_config/config.json`；
- `80_logs/polarization_grid.csv` 或 `01_results/polarization_grid.parquet`；
- 每个有效 coarse k 点的 `point_metadata.json`；
- 被跟踪 mode 的 `eigenmodes/<mode>_Hz_center.parquet`；
- 完整的几何、mesh、eigensolver 和 polarization definition identity。

实现前先验证：

- 数据是完整、无重复的 Cartesian uniform grid；
- 当前目标为内部 `p2`、用户正式名称 `py`；
- `match_status` 只能是 `gamma/matched/repaired`；
- 所有新旧点使用相同几何、mesh、air plane、4 个预期本征值和
  `c_const/1.55[um]` shift；
- 源文件哈希在整个流程前后不变。

### 4.2 使用 `cx_norm/cy_norm`，不跨批次拼接 raw 幅值

COMSOL 本征矢允许任意整体复数倍数。即使几何和 mesh 不变，不同 k 点、不同求解批次
的 `cx_raw/cy_raw` 也不能仅靠一个全局最大值消除本征矢幅度自由度。新流程使用现有
air-plane L2 归一化：

$$
\tilde c_x=\frac{c_x^{\rm raw}}
{\sqrt{\int_{\rm air}(|E_x|^2+|E_y|^2)\,dA}},\qquad
\tilde c_y=\frac{c_y^{\rm raw}}
{\sqrt{\int_{\rm air}(|E_x|^2+|E_y|^2)\,dA}}.
$$

用于定位和绘图的标量为 $|\tilde c_y|$。合并 coarse 与 refined 点后，只做一次共同
归一化：

$$
C_{\max}=\max_{\rm all\ valid\ points}
\sqrt{|\tilde c_x|^2+|\tilde c_y|^2},\qquad
M_y(\mathbf k)=\frac{|\tilde c_y(\mathbf k)|}{C_{\max}}.
$$

该决定会使新图成为独立命名的物理加密结果，而不是覆盖使用 raw 合同的历史标准 B4。
这是本 spec 最重要的审核项之一。

## 五、low-value valley 位置函数提取

位置提取只产生 COMSOL 扫描坐标，不产生最终幅值。

### 5.1 coarse 候选

在原始 uniform25 的 $|\tilde c_y|$ 网格上：

1. 逐个固定 $q_y$ 截面寻找严格离散局部最小值；
2. 逐个固定 $q_x$ 截面做同样操作，补充近水平或转向段；
3. 不使用 Otsu 或手写幅值阈值；只要一个点低于其同截面左右邻点，它就是采样候选；
4. 不预先规定 X 型、双曲线、相交、分裂方向或固定支线数；
5. 如果同一截面出现多个有效局部最低点，全部保留，增加 COMSOL 点而不主观删支线。

边界截面不做亚网格外推。内部三点可用简单抛物线插值细化 **位置坐标**，但拟合出的
幅值不得进入最终数据，也不得作为 line floor。

### 5.2 分段位置函数

- 相邻截面的候选只通过空间邻近和 mutual-nearest 连成边；
- 不确定配对同时保留，不强迫一对一连接；
- 每个连通段使用分片线性函数表示，局部可写成 $q_x=f_j(q_y)$ 或
  $q_y=g_j(q_x)$；
- 转向点由 row/column 两组分段共同覆盖；
- 局部切向由相邻线段平均得到，法向为其正交方向；
- 只允许对坐标做分片线性插值，不做 conic、全局多项式或拓扑模型拟合。

输出：

- `80_logs/valley_seed_candidates.csv`；
- `80_logs/valley_position_segments.csv`；
- `00_model/valley_refinement_sampling_plan.png`。

采样计划图叠加 coarse 点、位置分段、法向＋切向十字 stencil 和 junction 小网格，但不作为
最终物理图。`--prepare-only` 完成后，用户可先审核该图和点数。

## 六、局部 COMSOL 加密策略

### 6.1 拟写入 `scripts/parameter.json` 的参数组

在 `unit_cell_2d` 内增加：

```json
"valley_refinement": {
  "normal_half_width_over_G": 0.003125,
  "normal_points_per_round": 3,
  "anchor_stride": 3,
  "tangent_half_width_over_G": 0.0125,
  "tangent_points_per_round": 5,
  "longitudinal_max_gap_over_G": 0.00625,
  "refinement_rounds": 2,
  "refinement_ratio": 4.0,
  "junction_half_width_over_G": 0.0125,
  "junction_points_per_axis": 7,
  "max_edge_recenters": 2,
  "maximum_new_points": 1000,
  "display_points_per_axis": 577
}
```

这些数值针对当前 uniform25：原 coarse 间隔为
$\Delta q/G=0.0125$。

- 普通 valley 分支按连通段每 3 个 coarse candidate 选择一个采样锚点，端点保留；
- 第一轮法向范围为 $\pm0.003125G$，3 点间隔为 $0.003125G$；
- 第一轮切向范围为 $\pm0.0125G$，5 点间隔为 $0.00625G$；法向与切向中心去重，
  因而每个内部锚点最多仍形成 7 个十字 stencil 点，但采样权重转向 valley 纵向；
- 对每条分片线性 valley segment 检查其中点到已有新计划点的距离；若超过
  `0.00625G`，且该中点不与 uniform25 或其他计划点在 6 位坐标精度下重复，则加入一个
  仅用于纵向覆盖的真实 COMSOL gap-fill 点；
- 第二轮以第一轮实际最低 COMSOL 点为中心，半宽除以 4，5 点间隔为
  $0.00078125G$；
- junction 第一轮为 $[-0.0125,0.0125]^2$ 的 7×7 小网格；第二轮同样缩小 4 倍；
- 参数必须从源 coarse 间隔做一致性检查；改变源网格时不能静默沿用不合理绝对值。

### 6.2 普通 valley 段

对稀疏选择的位置函数锚点，沿估计法向和切向生成十字形奇数点 stencil。每一轮：

1. COMSOL 真实求解 stencil 上每个新 k 点；
2. 在实际求解值中直接选择 $|\tilde c_y|$ 最低的离散点；
3. 下一轮以该实际最低点为中心缩小扫描范围；
4. 不用抛物线或其他拟合生成幅值；
5. 如果最低值位于 stencil 端点，沿该方向平移扫描窗并真实补算；
6. 超过 `max_edge_recenters` 仍在端点时，将该段标记为 unresolved，不生成伪值。

### 6.3 交叉、近分裂和切向不稳定区

当候选图出现多分支汇聚、分支间距小于一个 coarse 间隔，或局部切向无法稳定定义时，
不使用任意法向。改为围绕该区域生成二维 Cartesian 小网格。该小网格中的每个点同样
由 COMSOL 求解。是否相交、是否存在奇点、低谷向哪个方向分裂，最终由这些真实二维
点自己显示。

### 6.4 点集合与边界

- 所有点必须位于源窗口 `[-q_max,+q_max]^2`；
- coarse 已有点直接复用，不重复求解；
- refined 点在 6 位 k 字符串精度下必须唯一；
- 生成超过 `maximum_new_points` 时 prepare 阶段失败并报告，不自动截断；
- `sampling_plan.json` 固定所有坐标及其来源，正式求解过程中不静默改变计划；
- 边界重心或额外 recenter 产生的新点通过 append-only 计划版本记录。

## 七、COMSOL 求解与模式跟踪

### 7.1 可复用内容

复用现有：

- `band.ReusableSimulationRun`：一个几何、一个 mesh、多 arbitrary-k 点；
- `persist_k_point_solution()`：保存本征频率和有效 mode 的中心 Hz；
- `compute_polarization_data()`：一次 air-plane batch 导出 raw 与 planar-L2 数据；
- `select_point_modes()`、Hz overlap、频率容差和 ambiguity gap；
- 点级 `solver_config.json`、`point_metadata.json` 和 atomic writes。

需要把当前私有 helper 整理为公共复用接口，禁止在新脚本中复制求解循环。

### 7.2 求解合同

- 几何：从 `scripts/parameter.json:unit_cell_2d` 读取，并与源 `config.json` 逐项一致；
- mesh：与源 series 相同，当前为 mesh5；
- 本征值：`eigenmode_pair_count=2`，即期望 4 个本征解；
- shift：固定 `c_const/1.55[um]`；
- k：arbitrary `(kx/G,ky/G)`，通过现有 `kx*G, ky*G` 参数写入；
- 模式：内部 `p2`，用户正式名称 `py`；
- 每个新点使用最近的有效 coarse/refined Hz 场作为 overlap predecessor；
- 先按位置函数顺序和离源 anchor 的距离排序，保证局部连续跟踪；
- ambiguous/missing 点做一次真实 re-solve repair；仍失败则标记无效并阻止正式 B4
  被宣称为 complete。

### 7.3 恢复、模型与并发

- 每个 k 点独立缓存，身份包含源 hash、geometry、mesh、solver、采样计划版本和坐标；
- 中断后重新执行默认跳过完整点，只补 pending/failed 点；
- 一个 case 内只创建一个 COMSOL model，geometry/mesh 各构建一次；
- 完成或安全中断前保存
  `00_model/unit_cell_2D_valley_refinement.mph`；
- 使用 case-local lock 和 PID 元数据；检测到本项目另一 COMSOL workflow 时默认拒绝
  启动，只有 `--allow-concurrent-comsol` 才允许并发；
- `--prepare-only` 和 `--plot-only` 永不创建 mph client。

## 八、基于真实新旧数据绘图

### 8.1 合并规则

- 合并 coarse 与 refined 的全部有效 `cx_norm/cy_norm`；
- 相同坐标必须数值一致，否则报告冲突，不选择平均值；
- 共同归一化只在完整合并集合上执行一次；
- 每一个存储节点均保留原始 COMSOL 值和来源标签 `coarse/refined/round`；
- 不把任何插值像素回写为数据点。

### 8.2 唯一允许的连续显示操作

对不规则真实点集合建立 Delaunay 三角剖分，并在每个三角形内部做重心分片线性插值：

$$
M(\mathbf q)=\lambda_1M_1+\lambda_2M_2+\lambda_3M_3,
\qquad \lambda_i\ge0,\quad \sum_i\lambda_i=1.
$$

因此插值值始终位于该三角形三个真实 COMSOL 节点的最小值和最大值之间，不产生 cubic
overshoot。禁止：

- Gaussian filter；
- bicubic/spline/RBF/Gaussian-process；
- minimum-curvature；
- 沿 valley 的幅值拟合；
- line floor、distance suppression、中心像素恢复；
- 人工连接、断开或压暗低谷。

源 uniform25 覆盖完整方形凸包，因此不需要边界外推。

### 8.3 新图合同

新图独立命名：

```text
10_overview/B4_cy_magnitude_log_comsol_refined.png
```

- 标题：`Normalized |c_y(k)| (py)`；
- colormap：`magma`；
- norm：`LogNorm`；
- 坐标：源窗口，0.05 主刻度，`kx-ky` equal aspect；
- colorbar：与图框等高；为与现有系列比较，拟继续显示
  `10^-5` 到 `10^-1`；这只是显示 clipping，不修改数值；
- 不叠加位置函数、采样点、拟合线、中心 marker 或拓扑文字；
- 图外 metadata 明确注明 `coarse + COMSOL-local-refinement + piecewise-linear`。

## 九、输出目录与命名

新计算不写入被否决的历史 series，也不创建第四个 `.out` 一级目录。拟使用：

```text
scripts/.out/unit_cell_2D/
  unit_cell_2D_valley_refined_(245-0.96-1.156)_mesh5_
  kmax0.15_source-uniform25_band-py/
```

目录合同：

```text
00_model/
  geometry.png
  valley_refinement_sampling_plan.png
  unit_cell_2D_valley_refinement.mph
01_results/
  refined_points/<point>/eigenfrequencies.csv
  refined_points/<point>/eigenmodes/*_Hz_center.parquet
  refined_points/<point>/point_metadata.json
  refined_polarization_points.parquet
  merged_polarization_points.parquet
10_overview/
  B4_cy_magnitude_log_comsol_refined.png
80_logs/
  valley_seed_candidates.csv
  valley_position_segments.csv
  refinement_points.csv
  band_tracking_refined.csv
  comsol_progress.log
  stdout/stderr logs
99_config/
  config.json
  sampling_plan.json
  convergence_summary.json
  run_summary.json
```

内部 band 仍写 `p2` 以兼容模式跟踪；目录、图名和标题使用用户正式名称 `py`。

## 十、当前被否决方案的处理

批准实施后：

1. 新流程不导入或调用 topology-adaptive、minimum-curvature 或 B5 幅值重建函数；
2. 当前 `--valley-reconstruction-b4-only` 不再作为受支持的正式路径；
3. 文档顶部把现有重建方法标记为 rejected/superseded；
4. 历史 PNG、CSV、JSON 保持原位作为过程记录，除非用户另行明确要求删除；
5. 标准 `B4_cy_magnitude_log.png` 保持不变；
6. 新 COMSOL-refined 图使用新文件名，不能无备份覆盖旧图。

是否在实施时彻底删除旧函数和 CLI，还是只将它们移入 archived 区域，作为本次审核的
一个明确决策点。推荐：移除公开 CLI 和正式调用，保留历史文档与结果，不删除已有数据。

## 十一、实现阶段

审核通过后按以下顺序实施：

1. 在 `parameter_config.py` 增加 `valley_refinement` dataclass、解析和验证；
2. 从现有 unit-cell 2D runner 提取 arbitrary-point 公共求解与模式跟踪接口；
3. 实现 pure-Python coarse valley candidate、位置分段和采样计划；
4. 实现 analysis 目录的完整编排脚本和 pyproject 命令；
5. 实现点级缓存、append-only recenter、lock、resume 和 summary；
6. 实现 coarse/refined 合并及 Delaunay 分片线性绘图；
7. 更新本方法文档和 README 命令；
8. 运行单元测试、fake-runner 集成测试、`py_compile` 和 `git diff --check`；
9. 只对 ζ=1.156 执行 `--prepare-only`，提交采样计划 PNG、点数和预计工作量供审核；
10. 用户再次确认采样计划后，才启动正式 COMSOL；
11. 首轮完成后检查 unresolved/ambiguous、最低点是否在 stencil 内部和第二轮变化；
12. 完整通过后生成新 B4，再决定是否批量处理其他 ζ。

## 十二、验证与验收标准

### 12.1 不启动 COMSOL 的测试

- coarse 数据完整性、哈希和 source identity；
- row/column 局部最低点和分段连接的确定性；
- 不预设拓扑、分支数或中心连接；
- 普通法向点、junction 小网格、边界裁切和坐标去重；
- round/recenter 只选择真实计算点，不拟合幅值；
- 最大点数保护、参数错误和计划 hash；
- cache identity、resume、冲突拒绝和 lock；
- coarse/refined 合并不平均冲突数据；
- 线性三角插值在节点精确、无 overshoot、无外推；
- 图名、标题、colormap、clim、ticks、equal aspect 和 colorbar 对齐；
- 标准 B4、历史参考图及源 CSV/JSON/MPH 哈希不变。

### 12.2 fake COMSOL 集成测试

- 一个 geometry、一个 mesh、多个 arbitrary-k solve；
- 每点 4 个预期本征解；
- py 模式使用最近 coarse/refined Hz predecessor；
- ambiguous repair 必须真实 re-solve；
- 中断后只补缺失点；
- `--prepare-only` 与 `--plot-only` 不调用 mph。

### 12.3 正式 COMSOL 验收

- sampling plan 与实际点集合逐点一致，额外 recenter 有 append-only 记录；
- 所有正式绘图点均有 COMSOL 来源或位于真实节点组成的线性三角形内部；
- 不存在 missing/ambiguous/unresolved 点；若存在则 status 不能为 complete；
- 每个普通 valley 截面的最低点不位于最终 stencil 边界；
- 第二轮相对第一轮的最低点位移不超过第一轮采样间隔；
- 关键区域可按需要增加第三轮做收敛复核，但不能用图像处理替代；
- ζ=1.156 的中心奇点和 valley 连接形态由 junction COMSOL 小网格自行给出；
- 最终 PNG 逐张视觉检查，确认无插值断层、colorbar/坐标/标题问题；
- `run_summary.json` 为 complete，stderr 为空或全部解释，MPH 与所有点级结果存在。

## 十三、风险与控制

### 模式误跟踪

局部低值区和 Γ 附近可能存在简并或频率接近。控制方式是保存全部 4 个候选频率和 Hz，
使用多个最近 predecessor，并把 ambiguity 作为硬失败；不能因为目标图需要低值而选择
幅值更低的另一个 mode。

### refined 点数量失控

位置函数的所有歧义都保留可能增大点数。prepare 阶段必须显示普通点、junction 点、
去重点和预计 solve 数；超过 1000 直接失败，由用户调整参数。

### 不规则三角剖分的长瘦三角形

完整 uniform25 coarse 网格保证全域基本覆盖；refined 点只在局部增加节点。验收时记录
最大三角边长，不允许跨越 coarse 网格空洞。若出现退化三角形，优先增加真实 COMSOL
点，不做平滑修补。

### 不同批次归一化

新流程强制使用 planar-L2 `cx_norm/cy_norm` 并验证 normalization definition version。
若源 series 缺失这些字段或 definition 不一致，拒绝运行，不能回退 raw 拼接。

### 与现有计算进程冲突

默认 lock 和进程预检拒绝并发 COMSOL。只有用户确认资源允许并显式传入
`--allow-concurrent-comsol` 才继续。

## 十四、回滚

- 在 COMSOL 启动前回滚：删除尚未批准的新增源码、测试、CLI 和空的新 series；
- 启动后不通过删除正式计算数据来回滚；新目录保持独立，可继续恢复或只读审计；
- 源 uniform25、标准 B4、历史参考与旧重建结果始终不被覆盖，因此无需恢复；
- 若新路线失败，只将新 series 标为 failed/incomplete，不回退到人为 floor 图作为物理结果。

## 十五、待审核确认项

请重点审核以下四项：

1. 新旧批次统一使用 planar-L2 `cx_norm/cy_norm`，新图不再沿用 raw B4 的跨点幅度合同；
2. 当前使用稀疏锚点的法向＋切向十字 stencil；法向间隔为
   `0.003125G -> 0.00078125G`，junction 使用两轮 7×7；
3. 新结果写入 `scripts/.out/unit_cell_2D/` 下的独立 series；源 uniform25 目录只读；
4. 旧公开 reconstruction CLI 从正式路径移除，但历史结果和文档保留，不删除数据。

审核通过后，先实现代码和测试，再只生成 ζ=1.156 的 `--prepare-only` 采样计划；该采样
计划还需一次人工确认，之后才实际启动 COMSOL。

## 十六、实施记录（2026-08-25）

### 已完成源码

- `parameter_config.py` 和 `scripts/parameter.json`：加入并验证
  `unit_cell_2d.valley_refinement`；旧参数快照缺少该组时使用同一组审核默认值。
- `comsol_workflow/unit_cell_2d_valley_refinement.py`：实现源数据/哈希验证、严格局部
  最小值、仅坐标抛物线细化、mutual-nearest 分段、数据驱动 junction、六位坐标去重、
  首轮采样计划、真实点冲突检查、共同 planar-L2 归一化及 Delaunay 线性绘图。
- `scripts/analysis/unit_cell_2D_valley_refine.py`：实现完整 CLI；`--prepare-only` 在任何
  COMSOL/mph 模块导入前退出。
- `scripts/run_sweep/unit_cell_2d_valley_solver.py`：实现点级缓存、最近 Hz predecessor、
  ambiguity 真实重算、真实最低点驱动第二轮、append-only edge recenter、MPH 保存、合并
  与最终绘图编排。
- `pyproject.toml`：加入 `comsol-unit-cell-2d-valley-refine`。
- `tests/test_unit_cell_2d_valley_refinement.py`：覆盖默认参数、源网格、严格最小值、
  junction、坐标去重、冲突拒绝和线性插值无 overshoot。

### ζ=1.156 prepare-only 实际结果

输入为旧 `scripts/.out/unit_cell_2D/` 下的 uniform25 完成结果及其 task-local
`scan_parameter.json`。输出位于 unit-cell 2D 专属的 `scripts/.out/unit_cell_2D/`：

```text
unit_cell_2D_valley_refined_(245-0.96-1.156)_mesh5_
kmax0.15_source-uniform25_band-py/
```

- coarse 间隔：`0.0125 G`；
- strict local-minimum candidates：66；
- mutual-nearest segments：46；
- 数据驱动二维 junction：1；
- 第一轮去重后新点：189，其中 20 个 valley 锚点产生 56 个法向点、72 个切向点，
  12 个 segment gap-fill 纵向点，中心 7×7 junction 产生 49 点；
- 按两轮和最多两次 recenter 的硬上限估计：1000 点；
- COMSOL 启动：否；
- 源结果与标准/历史 B4：未覆盖。

当前唯一待办是用户审核
`00_model/valley_refinement_sampling_plan.png`。确认后才允许正式求解，并在真实首轮结果
产生后继续核对 edge recenter、模式歧义和第二轮收敛。

### 输出位置修正

首次 prepare-only 错误写入了 `scripts/.out/unit_cell_band/`。经用户确认，Cartesian
unit-cell 2D 扫描与局部加密必须统一位于 `scripts/.out/unit_cell_2D/`。入口使用
`UNIT_CELL_2D_OUTPUT_ROOT`，首次生成的独立 refinement 目录整体移至正确位置；源
uniform25 结果未移动、未覆盖。

### 点数上限修正

根据用户审核意见，`maximum_new_points` 从 1200 收紧为 1000。prepare、第二轮动态
加点和 append-only edge recenter 均使用同一个硬上限；任何阶段若将唯一新增 k 点推到
1000 以上，流程立即失败并保留已有真实结果，不截断点集、不继续启动额外求解。

### 第一轮纵向采样修正

根据用户审核意见，原逐 candidate 法向 stencil 点数过多且沿 valley 切向覆盖不明确。
当前改为每 3 个连通 candidate 选择一个锚点，在锚点处同时放置法向 3 点和切向 5 点，
中心重复点去重；junction 从 9×9 调整为 7×7。ζ=1.156 第一轮因此从 313 点降为
177 点，并明确包含 valley 纵向采样。位置分段仍只用于选点，不生成或修改物理幅值。

随后按用户要求把每个锚点的权重进一步转向纵向：法向由 5 点减为 3 点，切向由 3 点
增为 5 点。两者的当前间隔均为 `0.00625G`；去重后的第一轮总数仍为 177，但实际构成
由“92 法向 + 36 切向”变为“56 法向 + 72 切向”，junction 保持 49 点。

用户随后圈出锚点衔接处仍存在的纵向空档。实现新增 segment midpoint coverage 检查，
仅当中点距离当前新计划点超过 `0.00625G` 时补点。本例自动识别并补入 12 个位置，
第一轮由 177 增至 189 点。最终 CSV 内部重复坐标为 0，与原 uniform25 坐标重合数也为
0。gap-fill 仍只定义 COMSOL 扫描坐标，不提供任何拟合幅值，也不作为第二轮最低值
搜索 stencil。

用户进一步指出普通法向端点过于接近原 uniform25 网格，而没有充分贴近估计 nodal
line。所有普通锚点的法向半宽因此从 `±0.00625G` 统一缩小为 `±0.003125G`，仍使用
3 点；纵向 5 点、segment gap-fill 和中心 junction 均保持不变。该调整只移动新 COMSOL
扫描坐标，不改变 valley 位置函数或任何幅值。

### 正式求解启动

用户确认最终采样预览后，正式流程以 `--allow-concurrent-comsol` 启动。启动合同为
ζ=1.156、mesh5、内部 `p2`/正式 `py`、每点 2 对即 4 个本征解、首轮 189 个唯一新点、
全部动态新增点硬上限 1000。后台日志与活动进程元数据写入 case 的 `80_logs/`。
启动检查确认 case lock 和 `status=solving` 已建立、stderr 为空，并已产生至少 3 个
完整 refined 点；抽查点的模式状态为 `matched`，使用 planar-L2 air-plane
`cx_norm/cy_norm`，因此流程已越过 Python 初始化并进入真实 COMSOL 求解。

### 正式求解完成与绘图修复（2026-08-26）

正式求解最终扩展为 709 个唯一新增 k 点，全部 709 点均完成，模式状态全部为
`matched`，`unresolved_point_count=0`；模型已保存为
`00_model/unit_cell_2D_valley_refinement.mph`。COMSOL 求解无需重跑。

首次自动绘图在 coarse/refined 合并表写入 Parquet 时失败。原因是 coarse CSV 将
`kx_str/ky_str` 推断为数值，而 refined Parquet 保留为字符串，合并后的 object 列含有
混合 Python 类型。修复方式是从权威数值坐标 `qx_over_G/qy_over_G` 以六位精度重新生成
两个字符串标签；这只统一存储 schema，不改变坐标、偏振系数、归一化或绘图数据。

同时修复 `--plot-only`：该分支现在直接由源 series 和参数快照解析已有 refinement
series，再读取已完成表绘图，不再调用 `prepare_valley_refinement()`。因此不会把包含
第二轮与 recenter 的 709 点动态计划覆盖回 189 点首轮计划，也不会启动 COMSOL。

修复后仅执行 plot-only，得到：

```text
10_overview/B4_cy_magnitude_log_comsol_refined.png
01_results/merged_polarization_points.parquet
```

最终合并节点为 1334 个（625 个 coarse + 709 个 refined），`run_summary.status` 已更新为
`complete`。绘图继续严格使用统一 planar-L2 归一化、Delaunay 分片线性插值、`magma`、
`LogNorm(10^-5, 10^-1)`、equal aspect 和等高 colorbar；未加入拟合幅值、拓扑约束、
平滑修补或 COMSOL 重算。

验证结果：

- `tests/test_unit_cell_2d_valley_refinement.py`：9 passed；
- mixed numeric/string k 标签可成功写入 Parquet；
- `--plot-only` 回归测试确认不会重建采样计划；
- 相关三个 Python 文件通过 `py_compile`；
- 实际 PNG 已按导出分辨率完成视觉检查。

## 十七、709 点结果的否决性复核与纠偏方案（2026-08-26）

首次绘图后的数据质量复核否决了上一节的 `complete` 结论。709 个点确实都完成了
COMSOL eigensolve 和模式匹配，但“求解成功”不等于“局部 valley 采样已经收敛”。当前
结果必须标记为 `incomplete`，已有 COMSOL 点全部保留并允许在纠偏后继续复用。

### 17.1 已确认的实现缺陷

1. `valley_cross_stencil` 的自适应中心在法向和切向点的并集上取最低值。首轮 20 个
   锚点中有 14 个、第二轮 20 个锚点中有 13 个因此被切向端点牵引，中心沿 valley
   方向漂移；切向点本来只用于纵向覆盖，不能驱动法向 recenter。
2. 第二轮三点法向 stencil 的中心通常已由上一轮求解。坐标去重正确地避免了重复
   COMSOL solve，但 stencil 判定没有把该中心作为只读复用节点纳入，因此 20/20 个
   第二轮首批法向 stencil 在表中都只剩两个端点，边缘判定必然失真。
3. recenter 仅检查最新新增 batch，没有联合中心节点和同一窗口的已有节点。去重后最新
   batch 往往只是一侧半窗，程序由此连续向外漂移，达到两次上限后仍没有真正的内部
   最低点证据。
4. 正式流程没有生成合同要求的 `convergence_summary.json`，也没有 edge-minimum、
   窗口完整性和三角剖分质量门禁；只要模式都匹配就把状态改成 `complete`。
5. refined B4 使用逐点 air-plane L2 后的 `cx_norm/cy_norm`，而标准 B1--B4 合同使用
   `cx_raw/cy_raw` 和完整集合的一个共同 `max|c(k)|`。本数据的 air-plane L2 因子从
   4697.5 到 1,313,491.2，跨度 279.6 倍；它改变了跨 k 的 B4 幅度关系并造成大面积
   饱和块。相距约 `1.74e-4 G` 的 coarse/refined 近邻显示 raw 系数连续，而逐点 L2
   标量在部分轴向近邻出现额外跳变。

### 17.2 代码纠偏合同

- cross stencil 只以法向节点决定最低点、edge 和 recenter；切向节点只参与最终空间
  覆盖和插值；junction 仍使用真实二维网格最低点。
- 坐标可以只求解一次，但必须允许作为多个 stencil 的只读成员复用。自适应判断使用
  “当前新点 + 已求解中心点”的完整逻辑窗口，不能把去重误解为缺失中心。
- 每次 recenter 检查完整逻辑窗口而不是最新半窗；达到上限仍在边缘时状态必须为
  `incomplete`，不得生成或宣称 final B4。
- refined B4 恢复标准 B-panel 合同：全部 coarse/refined `cx_raw/cy_raw` 使用一个共同
  `max sqrt(|cx_raw|^2+|cy_raw|^2)`；不修改任何节点值，仍只做 Delaunay 分片线性显示。
- 正式完成必须存在 `convergence_summary.json`，并同时满足：模式全部有效、所有最终
  法向/junction 窗口最低点位于内部、逻辑 stencil 完整、点数不超过1000、三角剖分
  无退化和跨 coarse 空洞长边。任一项失败则保持 `incomplete`。

### 17.3 当前数据的恢复策略

709 个已求解点不删除、不覆盖。纠偏程序先从首轮有效窗口重新建立 normal-only 自适应
链，所有已存在的六位坐标均作为 cache hit/只读节点复用，只把缺失坐标写入独立的
`convergence_repair_points.csv`。第一批纠偏后重新做收敛审计；若仍需加点，继续
append-only，但“原709点 + 纠偏点”的唯一新增总数仍受1000点硬上限保护。

纠偏采样预览必须再次人工审核。审核前只允许修改代码、生成计划和用已有点重绘诊断图，
不得启动 COMSOL。

### 17.4 已完成的纠偏实现与验证

- cross 自适应搜索已改为只读取 normal 节点；tangent 节点仅保留为空间覆盖样本；
- 去重中心可作为只读逻辑节点参与法向和 junction 窗口判断；
- 缺失或未通过 `convergence_summary.json` 时，输出只能标为 `incomplete`，图只能记录为
  `provisional_figure`；
- refined B4 已恢复 raw B-panel 共同归一化合同，并修复混合数值/字符串 `kx_str` 导致的
  parquet 写入失败；
- 首批纠偏计划包含 95 个唯一新点：48 个 `round2_junction_refine`、24 个
  `round2_normal_refine`、23 个 `round1_normal_recenter`；与已有709点零重复，预计总数
  804，小于1000点硬上限；
- 独立保存并哈希 `convergence_repair_targets.csv`，覆盖20个 cross 窗口和1个 junction；
  即使某个窗口所需坐标全部已经存在、未向95点CSV贡献新行，也不能逃逸收敛审计；
- 聚焦测试共 16 项通过，覆盖 tangent 不得移动法向中心、去重中心复用、raw 幅值派生、
  纠偏点唯一性、总点数上限，以及审批恢复入口不得重建动态计划。

当前没有重新启动 COMSOL。下一步必须先审核
`00_model/valley_convergence_repair_sampling_plan.png`，批准后只求解
`80_logs/convergence_repair_points.csv` 中缺失的95点；709个已有缓存不得重算或覆盖。

批准后的唯一恢复入口为：

```powershell
uv run python scripts\analysis\unit_cell_2D_valley_refine.py `
  "scripts\.out\unit_cell_2D\unit_cell_2D_(245-0.96-1.156)_mesh5_kmax0.15_uniform25_band-p2" `
  --parameter-file "scripts\.out\unit_cell_2D\unit_cell_2D_(245-0.96-1.156)_mesh5_kmax0.15_uniform25_band-p2\99_config\scan_parameter.json" `
  --run-approved-repair
```

若预检明确发现的是另一个已知且允许并行的 COMSOL 项目，才额外加入
`--allow-concurrent-comsol`。该入口先核对审批 CSV 的 SHA-256，再 append-only 写入动态
计划；逐点更新 parquet 断点，最后生成正式收敛审计。真实数据 dry-run 已确认：两个审批
哈希均匹配、709+95=804、21个审计目标完整保留；当前95点未求解时有11个窗口不完整、
17个窗口最低值仍位于边缘，因此门禁拒绝 `complete`，已有点三角剖分质量通过。
