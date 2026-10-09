# Unit-cell 2D B4 valley reconstruction：已退役方案记录

## 当前状态

截至 2026-09-02，topology-adaptive valley reconstruction、minimum-curvature B4 和
quadratic B5 fit 均已从生产代码、离线 CLI 和活动测试中删除。它们不是受支持的物理结果
生成器，也不能从当前 `scripts/analysis/unit_cell_2D.py` 重新生成。

当前权威 owner 是：

- 标准 A--D：`comsol_workflow/unit_cell_2d_analysis.py:run_analysis()`；
- 标准 B4：直接使用 COMSOL 节点上的 `cx_raw/cy_raw`、全网格共享归一化和标准显示插值；
- 真实 valley 改善：`scripts/analysis/unit_cell_2D_valley_refine.py`，设计合同见
  `docs/spec/plan-execute_20260825_unit_cell_2D_comsol_valley_refinement.md`。

已有历史 PNG、CSV 和 JSON 保持原位，只能作为比较证据读取；本次退役不移动、不覆盖、
不删除任何 `scripts/.out/` 文件。

## 被退役的三种方法

### Topology-adaptive valley reconstruction

该方法从既有 25×25 full Cartesian 数据中检测低值截面极小值，以中心不定二次型和
bootstrap 判定 X-like、left-right 或 upper-lower 拓扑，再用双曲线或奇次多项式分支、
分支距离抑制、人工 line floor 和 Gamma 恢复构造连续图。

它的主要问题不是数值实现错误，而是把“两条谷线、中心二次型、分支连接方式、floor 和
中心恢复”作为先验写入显示场。粗网格没有独立证明这些拓扑和振幅，因此结果不能作为
`c_y(k)` 的新物理解。

### Minimum-curvature B4

该方法先识别低值样本中的局部异常并调整 valley floor，再在加密网格上最小化有界薄板
弯曲能量，生成平滑的标量延拓。另一个较早版本保持全部 25×25 样本精确不变。

这类结果回答“满足选定约束的最平滑场是什么”，而不是“未采样位置的 Maxwell 解是什么”。
平滑目标、低值修正和 box constraint 都是分析先验，不能确定精确零线、BIC charge 或
Gamma Hessian。

### Quadratic B5 fit

该方法在 Gamma 附近 `|kx/G|,|ky/G|<=0.05` 的原始样本中选择低值点，拟合
`A(kx/G)^2+B(ky/G)^2=C`；历史 zeta=1.156 示例固定 `C=0`，并把拟合分支叠加在
minimum-curvature 背景上。

它可以作为局部模型比较，但固定中心、二次形式和拟合窗口会直接决定零线形状。拟合残差
也只是模型残差，不是新的场求解误差，所以 B5 不应进入正式物理结果链。

## 退役原因与替代方案

三种方法共同尝试解决粗 Cartesian 网格下低值 valley 显示不连续的问题，但都通过图像或
模型重建补充了 COMSOL 未求解的振幅。2026-08-25 的决策改为：已有极小值只用于规划新增
坐标；每个新增 `cx_raw/cy_raw` 必须由相同几何、mesh、eigensolver 和 mode tracking 的
COMSOL 求解得到；最终只允许对真实节点做 Delaunay 分片线性显示。

因此本次删除主动放弃以下能力：

- 从历史 full-grid CSV 重建三类比较 PNG；
- 重新生成对应 candidate/sample/relaxation CSV 和 diagnostics JSON；
- 从命令行选择上述拟合或平滑路线。

标准 A--D、legacy full-grid 的正常读取/重绘、Quarter C2v 只读视图、winding 和真实
COMSOL valley refinement 均不受影响。

## 保留的历史证据

- `docs/spec/plan-execute_20260821_unit_cell_2D.md` 保留当时的实验过程和输出记录；
- `docs/spec/plan-execute_20260825_unit_cell_2D_comsol_valley_refinement.md` 记录替代方案和
  为什么旧 reconstruction 不再受支持；
- Git 历史保留被删除实现和测试；
- 已有输出文件保持原位，其名称不再代表当前可调用接口。

## 重新引入条件

只有出现新的、明确授权的非物理比较需求时，才可从 Git 历史恢复为隔离研究工具，并必须
同时满足：名称显式标注 model/display reconstruction；不得进入标准 A--D；不得覆盖已有
结果；不得把拟合值写成 COMSOL 样本；必须重新建立独立 spec、消费者和验收测试。

若目标是改善正式 valley 结果，应扩展真实 COMSOL refinement，而不是恢复这些方法。