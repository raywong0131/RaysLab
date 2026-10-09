# Blueprints 径向单胞放大带执行计划

## 目标

在现有 `preview.png` 下方增加一条单胞宽的径向放大带，从中心 parent cell 连续展示到
最外层 parent cell，使 SEM 空间包络造成的 nominal/corrected 孔形差异可直接观察。
DXF、参数、DRC 和 case 三文件合同保持不变，只原位重绘当前 case 的 `preview.png`。

## 最小绘图合同

- 输出：当前 `f4ea43b5` case 中的 `preview.png`；不新建预览目录。
- 面板顺序：上部保持整体版图、中心孔、边缘孔和 fabrication summary；下部新增跨页径向带。
- 数据源：同一次生成流程中的 nominal/corrected `LayoutPlan`，不读取或修改 DXF 几何。
- 径向路径：使用 feature metadata 中的 parent center/shell，选择从写场中心指向最外层
  parent cell 的一条完整晶格射线；只含 `cavity`、`cladding`。
- 范围：纵向严格限制为一个 cell 的横向宽度，横向从 shell 0 到最大 shell；保持几何等比例。
- 编码：nominal 为中性灰填充/深灰边，corrected 为蓝色半透明填充/深蓝边；不只依赖颜色。
- 标注：标题明确 nominal/corrected，横轴用 parent shell（中心到边缘），标出 cavity/cladding
  分界；不添加无信息量装饰。
- QA：使用当前 20+20、220 nm SEM 修正 case 的最终 PNG，在实际分辨率检查孔形、边缘灰色
  差分带、标题、刻度和布局。

## 实施

1. 在 `workflow/output/preview.py` 增加径向路径选择、坐标投影和绘制函数。
2. 将 figure grid 扩为底部跨栏面板，同时保持上部既有面板合同。
3. 增加聚焦测试，锁定 shell 选择、单胞范围、标题和输出 PNG。
4. 运行完整 blueprints 测试、`py_compile`、`git diff --check`。
5. 仅原位重绘当前 case 的 `preview.png`，同步更新 report 中该 PNG 的 SHA-256，保持 case
   仍严格三个文件且不重写 DXF。
6. 视觉检查完成后清理测试缓存，并回填本计划执行结果。

## 回滚

恢复旧的两行 preview grid 和相关测试，再使用现有 nominal/corrected 参数流程原位重绘
`preview.png` 即可；DXF 和 layout report 的非预览内容不受影响。

## 执行结果

- 已在 `workflow/output/preview.py` 增加 metadata 驱动的径向 parent-cell 路径选择、局部坐标
  投影和 nominal/corrected 叠加绘制。
- 当前 20+20 case 自动选择 shell 0–40，共 41 个 cell、246 个三角孔；实际周期为
  `0.82 µm`，单胞横向半宽为 `0.473427 µm`，cavity/cladding 分界为 shell 20/21。
- 最终图保持 x-y 等比例；nominal 使用灰色填充和深灰轮廓，corrected 使用蓝色半透明填充
  和深蓝轮廓，边缘 cell 可见随半径增加的灰色差分边带。
- 已原位更新当前 `f4ea43b5` case 的 `preview.png`，未修改 DXF；report 中 preview SHA-256
  更新为 `655764332365867a29160d0c4bb5c65de6269a2a7e6dc07ec65bba8a5d873e4d`，DXF SHA-256
  保持 `e53ee6f21ca5ff1a6831565076b6b2ea646eb9a0aa88bc8a010e5f9e341a5393`。
- 完整测试 `75 passed`；实际参数路径选择、整图和底部原始像素裁剪均已视觉检查。
- 未启动 COMSOL，case 仍严格只有 DXF、`preview.png`、`layout_report.json` 三个文件。

最终状态：已完成。
