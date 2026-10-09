# finite quarter 两套结果合并与高 Q 模式归档

## 目标与边界

将以下两套已完成结果在数据层合并，并原位更新旧结果根目录的 `overview`：

- 旧结果（目标）：`fq4m_retry_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f197.980_e15_20260903`
- 新结果（追加）：`finite_quarter_all_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f203.600_e6_20260906_retry1`

本任务不启动 COMSOL、不改写或拼接 `.mph`，也不搬移原始模式目录。两套原始结果通过带来源字段的统一本征频率表组成一个逻辑数据集；满足 `Q > 1000` 的模式只将 Hz/Wem 图片汇总到目标 `overview`。

## 数据与命名合同

- 旧结果从四个 `symmetry_*` case 的 `10_overview/eigenfrequencies.csv` 读取，共 120 行。
- 新结果从统一 case 的 `10_overview/eigenfrequencies.csv` 读取，共 48 行。
- 目标 `overview/eigenfrequencies.csv` 应为 168 行，保留两表全部列，并新增：
  - `record_uid`：由来源运行名和来源 mode UID 组成的全局唯一身份；
  - `source_run`、`source_mode_uid`、`source_mode_dir`：记录来源和原始模式目录。
- 旧模式的来源 UID 补齐为 `modeN_<id>_xPEC|PMC_yPEC|PMC`；新模式沿用已有 `mode_uid`。
- 不按频率近似去重。两次本征求解返回的每一行均保留，避免误删相邻或重复本征模。

## f-Q 图合同

- 分析问题：两次求解合并后，全部模式的频率与 Q 如何分布。
- 图形：单系列散点图；每个有限频率/Q 记录对应一个点，预计 168 点。
- 横轴 `Frequency (THz)`，纵轴 `Q`；两轴均为线性坐标。
- 所有模式使用相同颜色、marker 和大小；不按来源、symmetry ID、PEC/PMC 或有效性分色，不显示这些内部身份的图例。
- 复用 `run_finite_quarter_four_modes.save_frequency_q_plot()`，覆盖目标 `overview/f_Q.png`。

## 高 Q 归档合同

- 条件严格为有限数值且 `Q > 1000`，不额外按 `is_valid` 过滤。
- 不复制完整 mode 目录，也不保留 `overview/q_gt_1000_modes/`。
- 每个筛选模式只汇总 `Hz_Im_2d.png` 和 `Wem_2d.png`，顶层名称遵循现有
  `symmetry_<id>_x<boundary>_y<boundary>_modeN_<field>.png` 规则。
- 重建 `overview/q_gt_1000_modes.csv`，包含来源路径、overview 前缀和两张目标图片路径。
- 当前预检结果为旧解 6 个、新解 1 个，共 7 个。

## 安全、验证与执行记录

- 写入前校验两套表、所有来源 mode 目录、168 个唯一 `record_uid`，并检查顶层图片冲突。
- 在仓库同盘的短 `tmp/fqmerge_*` 路径先生成临时输出并完成复制，以避开 Windows
  `MAX_PATH`；成功后再替换 CSV/PNG，最后原子发布新增的顶层图片。
- 验证合并行数、有限点数、Q 条件、每个清单图片路径、复制图片大小、PNG 可解码和线性坐标绘图摘要。
- 原有 `overview` 中单模式预览图保留不变。

### 执行结果

- 已生成 `overview/eigenfrequencies.csv`：旧解 120 行、新解 48 行，共 168 行；
  `record_uid` 168 个且全部唯一，所有来源模式目录存在。
- 已重绘 `overview/f_Q.png`：168 个有限点，x/y 均为 linear，统一散点样式，无来源、
  symmetry 或 PEC/PMC 图例；PNG 可解码，尺寸 1144×880，并已目检标签、边界和点分布。
- 严格按 `Q > 1000` 筛得 7 个模式；`overview` 现有 7 张 `Hz_Im` 和 7 张 `Wem`
  图片。旧解 12 张图保持不变，新解 `mode6_3_xPEC_yPEC` 的两张图按相同规则新增并通过大小校验。
- 早先误建的 `overview/q_gt_1000_modes/` 整模式目录已删除且未保留。
- 已写入 `q_gt_1000_modes.csv` 和 `merge_summary.json`；摘要中的 f-Q 路径指向正式文件，
  清单记录每个模式的 overview 前缀及两张图片路径，临时 staging 和 QA 副本均已清理。
- 首次按整目录理解复制时暴露 Windows `MAX_PATH`；最终实现改为同盘短 staging 和仅复制图片，
  原始两套模式目录和 MPH 均未改动，未启动 COMSOL。
- `py_compile`、结果断言、Pillow PNG 校验和本任务 `git diff --check` 均通过。

### 202.5 THz 增量更新

- 追加结果：`finite_quarter_all_20-20_c245-0.96-1.156_cl242-0.98-0.93_m5_f202.500_e6_20260908`。
- 合并器从目标 `overview/eigenfrequencies.csv` 继续追加，保留旧结果 120 行和 203.6 THz
  结果 48 行，再加入 202.5 THz 结果 48 行；最终共 216 行且 `record_uid` 全部唯一。
- 202.5 THz 新结果没有 `Q > 1000` 模式；高 Q 汇总仍为原有 7 个模式和 overview 外层
  14 张 Hz/Wem 图片，没有创建高 Q 子目录，也没有覆盖这些图片。
- `overview/f_Q.png` 已用全部 216 个有限点原位更新，横轴 Frequency、纵轴 Q，均为线性坐标。
