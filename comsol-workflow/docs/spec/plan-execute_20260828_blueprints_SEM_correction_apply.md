# Blueprints 220 nm SEM 修正方案应用计划

## 目标

将 `SEM_air_hole_size.csv` 对应的未修正版图明确解释为设计三角边长 `220 nm`，构建并
保存可追踪的 SEM 径向修正方案，将方案写入当前 `main/parameter.json`，随后按当前
20+20 finite-cavity 参数生成新的 engineering-preview DXF case。

## 修正模型

- 实测代表值：每层 `去极值平均_nm`；
- 目标 CD：显式 `220 nm`；
- 内区参考：`layer >= 22` 实测值中位数；
- 三角边长到 polygon 平行边偏置：`offset = side_change/(2*sqrt(3))`；
- 过程增益：当前只有一轮数据，显式使用 engineering `process_gain=1.0`；
- 包络：`c0 + c1*rho^2 + c2*rho^4`；`c0` 由目标 CD 与内区参考直接确定，只拟合
  `c1/c2`，避免无中心测量时常数项漂移；
- 归一化半径：当前 40 壳层乘晶格周期，即 `32.8 um`；
- 作用域：仅 `cavity`、`cladding` region；外围开窗和参数标记不施加 SEM 包络。

## 程序与配置改动

1. 扩展 SEM calibration 模块，显式接受并输出 `target_cd_nm=220`；
2. 扩展 `fabrication.spatial_envelope.regions` 严格 schema 和修正流水线；
3. 更新测试，锁定目标 CD、系数、region scope 和非目标结构零包络；
4. 重新生成 `SEM_air_hole_correction.csv/png`；
5. 将 profile id、revision、provenance、radial coefficients、半径和作用 region 写入
   `main/parameter.json`；
6. 运行完整测试、py_compile、preflight、git diff check；
7. 运行正式 blueprint generator，检查新 case 中 DXF、preview 和 report，并用实际 PNG
   做视觉 QA。

## 发布边界

该方案仍使用未闭环标定的 `process_gain=1.0`，dose、polarity 和生产 DRC 未赋值，因此
配置与产物保持 `engineering_preview`，不得描述为生产发布版图。

## 输出与回滚

- 输出只增加一个参数 identity 对应的标准三文件 case；不覆盖旧 case；
- 回滚时恢复 `main/parameter.json` 的零包络及 engineering-unassigned profile，移除
  region-scope/schema/calibration 改动和对应测试即可；旧版图 case 不受影响。

## 执行结果

- 目标 CD 固定为 `220 nm`；内区参考中位数为 `220.341325 nm`。
- 保存的径向偶次多项式系数（polygon edge offset，nm）为
  `[-0.09853204031556974, 5.134538035287714, -30.140260018739795]`；
  拟合 RMSE 为 `0.566357622 nm`，最外圈修正为 `-25.104254024 nm`。
- 修正数据已原位更新为 `blueprints/data/SEM_air_hole_correction.csv` 与
  `blueprints/data/SEM_air_hole_correction.png`；QA 图已做实际视觉检查。
- 当前参数已登记 `profile_id=sem-air-hole-220nm-radial-v1`、`revision=1` 和数据
  provenance；空间包络仅作用于 `cavity`、`cladding`，开窗和参数标记保持零包络。
- 新 case：
  `blueprints/main/.out/dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_f4ea43b5`。
  case 严格只包含 DXF、`preview.png` 和 `layout_report.json` 三个文件，未覆盖旧 case。
- 版图共 `29,648` 个闭合 DXF entity；DXF 回读数量一致。DRC 结果为
  `0 error / 0 warning`，预览图已确认六角孔阵、6 个外围开窗和下方参数标记完整。
- 验证：blueprints 完整测试 `73 passed`；聚焦测试 `33 passed`；`py_compile`、
  preflight 和 `git diff --check` 均通过；测试与编译缓存已清理。
- 未启动 COMSOL。由于 `process_gain=1.0`、dose、polarity 和生产 DRC 阈值尚未闭环
  标定，最终状态保持 `engineering_preview`。

最终状态：已完成。
