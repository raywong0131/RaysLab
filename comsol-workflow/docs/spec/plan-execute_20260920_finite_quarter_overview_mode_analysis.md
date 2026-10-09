# finite-quarter overview 高 Q 模式完整分析

## 目标

对整合结果
`finite_quarter_all_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5`
的 `10_overview/q_gt_1000_modes.csv` 中全部 7 个模式建立独立目录：

`10_overview/<overview_prefix>/`

每个目录遵循 finite cavity 的 mode-centric 结果合同，包含：

- `10_overview/eigenfrequency.csv`、`mode_score.csv`
- `11_simulation_exports/Hz_center.parquet`、`Hz_Re_2d.png`、`Hz_Im_2d.png`、
  `Wem_2d.png`、`E_air.parquet`
- `12_farfield_FFT/` 的 A–D 四图和 `farfield_summary.json`
- `13_lattice_fourier_Hz/` 的 NPZ、四张 Fourier 图和两张 CSV

模式目录使用现有 `overview_prefix`，即当前 overview 图片前缀，不显化新的外部命名体系。

## 模式与复用策略

完整结果直接复用并复制到 overview：

- `symmetry_1_xPEC_yPMC_mode19`
- `symmetry_3_xPEC_yPEC_mode6`

以下模式已有 `Hz_center.parquet`，但缺少完整后处理，需要从保存的 solution 补齐：

- symmetry 1：`symmetry_1_xPEC_yPMC_mode25`
- symmetry 2：`symmetry_2_xPMC_yPEC_mode28`
- symmetry 3：`symmetry_3_xPEC_yPEC_mode21`
- symmetry 4：`symmetry_4_xPMC_yPMC_mode21`、`symmetry_4_xPMC_yPMC_mode26`

同一 symmetry 的多个模式共用一次 `.mph` 加载。只 attach 已保存 solution 并读取场，不重新建模、
mesh 或 eigensolve。频率必须与来源 `eigenfrequencies.csv` 逐模式一致。

## 数值合同

- 场重建继续调用 `run_finite_quarter.export_reconstructed_selected_results()`。
- far-field 使用已完成结果的合同：801×801 空气面、2001×2001 FFT、NA=0.9、自动读取空气面高度。
- finite-lattice Fourier 使用 `a=0.82 um`、bulk radius 19、rho grid 41、top-k 6、
  `hex_cyclic_quotient_v1`；晶格点从各 symmetry 的保存 config 读取。
- Q 与频率从当前 eigenfrequency metadata 读取，不硬编码模式数值。

## 安全、发布与验证

1. dry-run 校验 7 个唯一模式、来源路径、config、模型、频率表和目标冲突。
2. 在 `finite_cavity/.fq_mode_analysis_staging` 生成全部结果；失败时不发布半成品。
3. 所有模式必须具备规定的 19 个文件，far-field 状态必须为 `ok`。
4. 校验 JSON/CSV/parquet/NPZ 可读，PNG 可解码，频率和 Q 与选择表一致。
5. 全部通过后逐目录同盘原子发布，并写
   `10_overview/finite_mode_analysis_manifest.json`；随后清理 staging。
6. 运行聚焦 dry-run、`py_compile`、实际 PNG 检查与 `git diff --check`。

## 执行记录

- dry-run 确认选择表包含 7 个唯一 `Q > 1000` 模式；其中 mode19（symmetry 1）和 mode6
  （symmetry 3）已有完整结果并复用，其余 5 个模式从四套保存 solution 补齐。
- 首次启动受沙箱对 `~/.comsol/v63/tomcat/logs` 的写权限限制，COMSOL 尚未加载模型即退出；
  未发布任何结果。核对并清理 320.5 MiB staging 后，以提升权限重新执行成功。
- 四组保存解依次完成，未重新 mesh 或 eigensolve。新增模式的组级 Gamma score 分别为：
  symmetry 2 mode28 `0.000665636`，symmetry 4 mode21/mode26 组 `0.01898`，
  symmetry 3 mode21 `0.000232545`，symmetry 1 mode25 `0.115057`。
- 已在根级 `10_overview` 发布 7 个以 `overview_prefix` 命名的模式目录；每个目录恰好
  19 个文件、11 张 PNG，并包含完整的 10/11/12/13 分段结果。
- 全部 far-field summary 状态为 `ok`；JSON、CSV、parquet、NPZ 和 PNG 均通过读取/解码，
  频率与 Q 和选择表一致。代表性 mode25 的 Hz、far-field A 面板及 k-weight 图已经目检。
- `finite_mode_analysis_manifest.json` 状态为 `complete`，staging 已清理，Python/COMSOL
  进程均退出；`py_compile` 与 `git diff --check` 通过。
