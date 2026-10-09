# Hz 场图统一色标

以 `results/finite_cavity_reference_results/shaolei_hex_side2_forward/singularity_original_r19_c20/selected_mode_14_plots/Hz_Im_2d.png` 为标准。

- Hz 实部/虚部场图使用 Matplotlib `RdBu_r`，负蓝、零白、正红。
- 每张单独场图按当前分量的有限样本最大绝对值归一化，固定范围 `[-1, 1]`，九个刻度、两位小数。零场保持零；图上标识归一化，不再标 A/m。原始场数据保持物理单位。
- 覆盖共享 PNG/SVG/HTML 绘图、控制器及兼容入口的 Hz 二维图、独立后处理模块的 Hz 场图。频率/Q、几何范围和输出文件名沿用当前数据与合同。
- Hz 强度、模长、相位及 Fourier 权重不是有符号分量，不套用有符号归一化。后续若做幅值比较，应使用完整比较集合的共同 scale。
- 不重跑 COMSOL、不批量重绘既有结果。用合成数据、接口替身及实际 PNG 检查验证。

执行与验证：已完成。

- 共享绘图函数增加固定 Hz 色标，PNG / SVG / hybrid SVG / PNG 嵌入 HTML 复用；COMSOL 控制器及三个 spatial 兼容入口的二维 Hz 图片接入同一函数，原始 TXT 导出保留。
- 独立后处理包的 mode-fields 和 tetramer 场图均接入同一归一化合同。
- 主项目 38 项聚焦测试通过；排版修正后重新运行受影响的 12 项绘图测试通过。独立后处理 5 项测试通过。
- 相关源码与测试通过 py_compile；两个仓库的本轮 diff 空白检查通过。
- 使用合成场实际导出 PNG 并检查：色带、九个刻度、等高 colorbar、等比例坐标和无重叠的频率/Q/分量标签。归一化说明移至色标旁，避免与 Q 碰撞。
- 未启动 COMSOL，因此未做真实 COMSOL 运行验证；其图导出路由用接口替身验证。未覆盖、移动或删除既有结果；验证用临时 PNG 已清理。
