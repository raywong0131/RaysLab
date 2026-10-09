# 11_simulation_exports 场图样式统一

## 方案概述

将完整有限腔与四分之一对称恢复的 `Hz_Re_2d.png`、`Hz_Im_2d.png`
和 `Wem_2d.png` 统一交给共享的 Matplotlib 绘图函数生成。这样透明背景、
COMSOL 参考色表、频率/Q/输出量标注、单位和 colorbar 几何关系只维护一份。

## 范围与数据流

- 输入仍来自现有 COMSOL 解：`get_2d_fields()` 提供坐标和场值，
  `get_eigenfrequencies()` 提供 frequency 与 `ewfd.Qfactor`。
- `Hz_center.parquet`、`eigenfrequencies.csv`、文件名和最终目录布局不变。
- 仅替换 `11_simulation_exports` 内三类 PNG 的绘制，不改求解、模式编号、
  有效性判定、远场或 Fourier 后处理。
- 对现有 `scripts/.out/finite_cavity` 结果提供一次性批量更新：Hz 从已保存的
  `Hz_center.parquet` 重绘，Wem 从对应 `finite_quarter.mph` 的已保存解读取，
  不重新求解。

## 兼容策略

- Hz 保持关于零点对称的线性色标；Wem 保持从零开始的线性色标。
- 颜色节点从用户提供的 COMSOL `Wave` 与 `HeatCamera` 参考图采样固化，
  不在运行时依赖参考图文件。
- PNG 保留原文件名与尺寸类别，并启用真正的 alpha 透明通道。

## 验证与回滚

- 单元测试检查调用方传入的 frequency、Q、输出量、单位与缩放模式。
- 绘图测试检查色表关键颜色、PNG alpha、顶栏文本，以及主图黑框和
  colorbar 的高度/垂直位置完全一致。
- 运行聚焦 pytest、`py_compile` 和 `git diff --check`。
- 回滚时移除共享绘图模块，并恢复两个导出函数原来的 COMSOL/局部绘图调用。
- 批量更新先在隐藏 staging 目录完整生成并检查全部图片；正式替换前复制旧图，
  任一替换失败时恢复全部旧图，成功后才清理 staging 与临时备份。

## 执行结果

- 新增 `comsol_workflow.field_plotting`，固化参考图采样得到的 `Wave` 与
  `HeatCamera` 色表，并统一生成带 alpha 通道的场图。
- 完整有限腔和四分之一恢复导出均传入对应 mode 的 frequency、
  `ewfd.Qfactor`、输出量名称与单位；图片内不再显示内部 mode 编号。
- colorbar 使用主图等比例调整后的最终轴位置，实测上下边与黑框完全对齐。
- 使用现有 `Hz_center.parquet` 生成临时 Hz/Wem 样图完成视觉检查；临时文件
  已清理。
- 新增 `scripts/analysis/update_finite_field_plots.py`，并对现有 finite-cavity
  结果执行事务式批量更新：56 个 mode、25 个已保存 MPH、168 张场图全部
  成功替换；没有重新 mesh 或求解。
- 批量更新后检查 168/168 张图片：均非空、尺寸均为 `1540×1320`、左上角
  alpha 均为 0；隐藏 staging 与临时备份已在成功后清理。
- `py_compile` 通过；新增/聚焦测试 13 项通过；完整测试 276 项通过，只有 3 个
  仓库既有的 `PytestReturnNotNoneWarning`。
- 现有正式结果和后续新运行生成的 `11_simulation_exports` 均已采用新样式。
- 根据后续版式反馈，顶栏调整为单行：frequency 左对齐、Q 位于主图 60% 宽度、
  输出量右对齐、单位位于 colorbar 上方；该行下移靠近黑框，colorbar 同时左移。
- 主图黑框恢复 `x/y (μm)` 几何尺度、外置数字和朝内刻度线。已仅用
  `mesh9_shiftprof-tanhpow-l4-p2/shift0.090/mode1` 的三张场图更新测试，
  其余既有 mode 尚未批量套用这一轮版式。
- frequency 显示固定保留两位小数，例如 `198.21 THz`。
- 为指定测试 mode 新增 `Hz_Re_2d.html`、`Hz_Im_2d.html`、`Wem_2d.html`：
  HTML 以 base64 内嵌当前透明 PNG，单文件、响应式、无网络依赖，因此现有
  顶栏、坐标尺度、色表和 colorbar 结构逐像素保持不变。
- 矢量格式测试：同一 mode 的全矢量 SVG 每张约 462.9 MiB；600 DPI hybrid
  SVG 仅将场分布栅格化，其他元素保持矢量，三张分别约 9.7、10.0、6.1 MiB。
- 已将 hybrid SVG 作为有限腔场图的常规矢量交付约定写入 `MEMORY.md`：场分布
  以 600 DPI 栅格化，黑框、几何刻度、文字与 colorbar 保持矢量，并使用
  `*_hybrid.svg` 后缀；超大的全矢量场图不作为常规输出。
- PNG、全矢量 SVG 与 hybrid SVG 的保存入口均改为按实际可见内容紧边界裁切，
  仅保留 `0.03 in` 防截断余量。指定测试 mode 的 PNG 由 `1540×1320` 缩小为
  Hz `1433×1088`、Wem `1413×1088`；hybrid SVG 的 viewBox 分别缩小为
  Hz `469.897×355.887`、Wem `463.193×355.887`，且每张仍只有一个栅格场图层。
- 同步重建三份 HTML，使其按图片自然尺寸显示并仅在视口较窄时缩放，避免紧边界
  PNG 被固定宽度再次放大。
- 紧边界改动完成后，`py_compile` 与 `git diff --check` 通过，完整测试为
  `280 passed`；仅保留 3 个既有的 `PytestReturnNotNoneWarning` 和 1 个本地
  pytest 缓存目录警告。
