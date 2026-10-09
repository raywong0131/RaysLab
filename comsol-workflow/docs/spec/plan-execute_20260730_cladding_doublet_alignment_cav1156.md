# cavity p1/d1 对齐 cladding 双重态优化计划与执行记录

## 目标

以用户更正后的
`cav(245-0.96-1.156)_clad(240.2-0.967-1)_mesh5` 为出发点，固定 cavity
全部参数、cladding `zeta=1`、mesh5、材料与 unit-cell 求解设置，只优化 cladding
的 `b0_nm` 和 `eta`，使 Gamma 点：

- cladding `p1/p2` 双重态对齐 cavity `p1`；
- cladding `d1/d2` 双重态对齐 cavity `d1`。

上一轮 cavity `zeta=1.155` 源于输入笔误，其结果不作为本任务结论；对应独立结果目录
予以保留，不删除或覆盖。

## 复用方案与数据流

复用 `scripts/run_sweep/optimize_cladding_doublet_alignment.py` 的二维局部 Jacobian
搜索和 0.1 nm / 0.001 邻域细化，不新建第二套数值实现。

```text
固定 cavity(245, 0.96, 1.156) Gamma 求解或兼容缓存复用
  -> 识别 cavity p1@Gamma 与 d1@Gamma
  -> 从 cladding(240.2, 0.967, 1.0) 出发
  -> b0/eta Jacobian 更新
  -> 0.1 nm / 0.001 局部网格细化
  -> 按 p/d 双重态对齐误差和简并劈裂排序
```

## 固定项与可变量

- 固定 cavity：`b0=245 nm, eta=0.96, zeta=1.156`。
- 固定 cladding：`zeta=1.0`；仅 `b0_nm`、`eta` 可变。
- 起点：`b0=240.2 nm, eta=0.967`。
- 固定 mesh：COMSOL auto mesh size 5。
- unit-cell 求解数读取当前 `parameter.json`，本次启动前核对值；求解 shift 保持
  `c_const/1.55[um]`。
- 不运行 full-band validation，优化阶段只求 Gamma 点。

## 缓存与输出

- 新结果写入包含正确 cavity `1.156`、zeta1 和 mesh5 的非冲突目录，不增加
  mode-count 后缀。
- cavity reference 和运行期间的每个 cladding candidate 均使用 k-point solver
  identity 校验；求解数、mesh 或 shift 不一致时不复用。
- 仅复用几何和 solver identity 都兼容的 Γ 点结果；目标频率始终使用正确 cavity
  `1.156` 重新计算并重排。现有 sweep 在成功汇总后会清理体积较大的 candidate
  中间目录，保留 cavity reference、完整扫描表和最佳参数。
- 不移动、覆盖或删除上一轮 `1.155` 或其他正式结果。

## 验证

- 检查固定参数、起点、目标定义、搜索精度、mesh、求解数和输出目录。
- 运行 doublet-alignment 聚焦测试、`py_compile`、`git diff --check`。
- 启动 COMSOL 前检查现有进程和 license/内存占用；最终读取
  `best_parameters.json` 与 `scan_summary.csv`。

## 回滚

如需回滚，仅恢复该 sweep 的任务常量；新结果保留在正确 cavity 的独立输出目录，
不影响其他 unit-cell、strip 或 finite 工作流。

## 执行结果

- 已将固定 cavity 更正为 `245/0.960/1.156`；cladding 起点保持
  `240.2/0.967/1`，mesh5。本次从 `parameter.json` 读取到的 unit-cell 本征模式
  设定数为 10，shift 为 `c_const/1.55[um]`。
- 独立求解并保存正确 cavity Gamma reference，得到
  `p1=196.564396032325 THz`、`d1=200.461469058620 THz`。
- 自适应步骤从起点移动到 `240.1/0.967/1`；随后完成 25 点局部网格及 2 个
  网格外 Jacobian 探针，共 27 个可行候选。未启动 full-band validation。
- 唯一通过 `0.02 THz` 对齐及劈裂阈值的最优点为
  `b0=240.1 nm, eta=0.967, zeta=1`。
- 最优 cladding 结果：`p1=196.577459605006 THz`、
  `p2=196.583138261819 THz`，双重态中心误差 `+0.015902901087 THz`、
  劈裂 `0.005678656813 THz`；`d1=200.441638256705 THz`、
  `d2=200.455819077792 THz`，双重态中心误差 `-0.012740391372 THz`、
  劈裂 `0.014180821088 THz`。最大中心对齐误差为 `0.015902901087 THz`。
- 结果写入
  `scripts/.out/unit_cell_band/cladding_doublet_alignment_to_cavity_p1_d1_cav(245-0.96-1.156)_zeta1_mesh5/`；
  关键文件为 `best_parameters.json`、`scan_summary.csv`、`scan_config.json` 和
  `cavity_reference/config.json`。误写 cavity `1.155` 的旧目录予以保留，但不作为
  当前结论。
- 验证：4 个 doublet-alignment 聚焦测试通过；优化脚本 `py_compile` 通过；
  `git diff --check` 通过。
- 最终状态：更正后的优化完成，最终 cladding 相对起点仅将 `b0` 从 240.2 nm
  调整为 240.1 nm；未修改 `scripts/parameter.json`，未覆盖其他正式结果，原有
  COMSOL PID 12252 未被中止。
