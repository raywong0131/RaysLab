# Blueprints SEM 三角孔径向修正数据执行计划

## 方案概述

基于 `blueprints/data/SEM_air_hole_size.csv` 中未进行空间包络修正的真实 SEM
三角孔边长数据，构建一套可审计、可复现的工程级相对版图补偿数据。当前任务不把结果
直接写入正式 `main/parameter.json`，不将其声明为生产标定，也不生成加工 DXF。

## 范围

- 新增独立的 SEM 标定解析、校验、径向拟合和 CSV/PNG 输出模块；
- 新增仓库内直接运行入口；
- 生成一份完整 40 壳层修正 CSV 和一张 QA 图；
- 增加聚焦测试与 README 使用说明；
- 不改变现有 nominal geometry、dose、DRC、DXF 和 release gate；
- 不启动 COMSOL，不覆盖 `main/.out/` 中任何正式 case。

## 数据合同与假设

1. `层数_外向内=1` 对应当前 20+20 finite layout 的最外层
   `parent_shell=40`，一般映射为 `parent_shell = total_shells + 1 - layer`；
2. `三角形N_nm` 和 `去极值平均_nm` 表示三角孔边长；程序会从原始样本复算并核对
   平均值和去极值平均值；
3. `layer >= 22` 的内区样本中位数定义为零补偿参考 CD；这个阈值作为显式 CLI 参数
   并写入输出表；
4. 未测但位于已测区间内的壳层使用一维线性插值；更靠中心的未测壳层使用参考 CD
   平坦外推；
5. 当前只掌握一次未修正版图的 SEM 数据，缺少多轮 design-to-SEM 响应斜率，因此使用
   显式 `process_gain=1.0` 构建相对工程补偿；生产使用前必须用二次实验重新标定；
6. 等边三角形平行边偏置与边长变化满足
   `edge_offset = side_length_change / (2*sqrt(3))`；
7. 推荐平滑包络使用现有 `radial_even_polynomial` 形式，并固定中心零项为 0：
   `offset_nm(rho)=sum(c[k]*rho^(2k))`。拟合只使用真实测量行，不用插值行。

## 数据流

```text
SEM_air_hole_size.csv
  -> 严格字段/数值/统计量校验
  -> 内区参考 CD + 40 壳层完整 profile
  -> 相对边长补偿与 polygon edge-offset 换算
  -> 固定 c0=0 的二阶径向偶多项式拟合
  -> SEM_air_hole_correction.csv + SEM_air_hole_correction.png
```

## 输出合同

- `blueprints/data/SEM_air_hole_correction.csv`：40 行完整壳层数据，包含原始/插值状态、
  参考 CD、表格补偿、推荐拟合补偿、预计修正后 CD、拟合系数和 source hash；
- `blueprints/data/SEM_air_hole_correction.png`：原始 SEM、补偿曲线和预计修正后 CD 三面板
  QA 图；
- 不生成额外 JSON、DXF 或调试文件。

## 验证

- 校验 29 个输入采样行、40 个输出壳层、层号唯一与映射正确；
- 校验原始平均值/去极值平均值可由样本复算；
- 校验三角形 edge-offset 换算、中心零补偿、拟合系数与误差；
- 运行 blueprints 完整独立 pytest、相关 `py_compile`、CLI 生成、`git diff --check`；
- 打开实际 PNG 检查标题、坐标方向、数据点、拟合线与布局。

## 回滚方案

删除新增的 SEM 标定模块、入口、测试、两项派生产物及 README 对应章节即可；现有正式
参数、版图生成器和 `main/.out/` case 均不受影响。

## 执行结果

已完成：

- 新增 `workflow/fabrication/sem_calibration.py`，实现严格 CSV 校验、内区参考、缺层
  插值、三角边长到 polygon edge offset 换算、固定中心零项的径向偶多项式拟合及
  CSV/PNG 输出；
- 新增 `main/build_sem_correction.py` 和 5 项聚焦测试，并更新 README；
- 生成 40 行 `data/SEM_air_hole_correction.csv` 与 220 DPI
  `data/SEM_air_hole_correction.png`；未生成额外 JSON 或 DXF；
- 29 行真实测量、6 行区间插值、5 行内区参考保持；参考 CD 为
  `220.341325 nm`；
- 推荐二阶 edge-offset 包络系数（nm）为
  `[0.0, 5.134538035, -30.140260019]`，最外圈推荐偏置
  `-25.005721983 nm`；拟合 RMSE 为 `0.566357622 nm` edge offset，最大绝对残差
  `1.334088211 nm`；
- 正式 PNG 已按实际分辨率视觉检查，面板、标题、图例、坐标方向和数据曲线无碰撞或
  截断；
- blueprints 完整独立测试 `70 passed`，相关 `py_compile`、任务文件
  `git diff --check` 均通过；主版图 preflight 为 29,648 features、DRC 0 error/
  0 warning，且没有写文件；
- 未修改正式 `main/parameter.json`，未覆盖或删除 `main/.out/` case，未启动 COMSOL。

最终状态：工程级相对修正数据已完成并可复现。由于 `process_gain=1.0` 尚未通过二次
SEM 闭环标定，该数据不得直接声明为生产工艺修正。
