# parameter.json 分区重构

## 目标与边界

统一参数入口按归属分区，保留当前物理数值、COMSOL 求解基准、结果目录及缓存身份。
以当前工作区为基线，保留已有修改；不启动 COMSOL，不修改正式结果。

## 实施计划

1. schema v2 保留根级 `mesh_size` 和独立 `unit_cell_2d`；新增 `cells`、
   `structure`、`unit_cell_band`、`strip_1d`、`finite_cavity` 分区。
2. 几何与结构共享关系保持现状；finite 的 footprint 和方向位移归 finite；
   strip 的模式数、切角、kx 和 quarter 对称性从原脚本默认值接入 JSON。
3. 在统一加载器解析分区，复用已有数值校验和 `SharedParameters`；兼容 v1 快照，
   v2 拒绝混用旧路径及拼错字段，2D 缺失时不得借用 cavity。
4. 中心频率回写保持原格式并仅更新目标字段；同步参数快照生成、直接读取者及文档。
5. 对比重构前后解析结果、目录与缓存身份；运行不启动 COMSOL 的配置/入口/扫描回归、
   相关 Python 编译检查及 `git diff --check`，清理本轮临时文件。

## 执行与验证

已完成。

- 迁移当前 parameter.json 至 v2 分区；保留当前 cavity/cladding、层数、mesh、
  中心频率、位移及全部 2D 数值。新增入口字段采用原脚本默认值：strip 模式数 2、
  切角 60°、kx=0；普通 quarter symmetry IDs=[1]。
- 加载器将 v2 分区映射为既有运行时属性；v1 快照仍受支持。v2 拒绝旧根级路径、
  未知字段、非法新入口参数及缺失 2D 对象。finite 显式 x/y 与 ratio 的互斥校验复用原实现。
- 更新中心频率回写、finite 优化快照与 theory-fit manifest 的参数读取；2D ζ 扫描原有路径
  不变，其快照与恢复测试通过。四扇区专用入口保留固定 1→2→3→4。
- 已同步 README、scripts/README、仿真与 2D 合同，新增 [参数指南](../parameters.md)。

验证：

- 本轮保存的重构前参数与元数据基线，与重构后 SharedParameters、to_metadata、
  2D cache identity、2D/strip/finite/quarter series 名逐项相等。
- 聚焦回归覆盖 test_parameter_config、run_band_pair、run_band_pair_model_reuse、
  run_band_solo、cladding_shift_profile、run_unit_cell_2d、run_unit_cell_2d_zeta_scan、
  optimize_finite_quarter_q、theory_fit_workflow、run_strip_output_layout、
  run_strip_bulk_rotation、finite_cavity_field_exports、run_finite_quarter_four_modes。
  共 284 个用例通过；首轮 283 通过，唯一失败为既有 alignment 测试临时 CSV 的
  WinError 5 原子替换错误，该用例独立重试通过（1 passed）。没有为文件锁修改业务代码。
- 测试通过 uv run --no-sync python 调用 pytest.main，临时将 mph.start/mph.Client
  替换为抛错函数，确保不启动 COMSOL。10 个相关 Python 文件通过 py_compile。
- 新文档 UTF-8 编码及本地链接检查通过；git diff --check 通过。
- 本轮 tmp/parameter-layout-20260908 基线、测试输出和编译产物已清理；正式结果目录
  未改动。保留用户原有工作区修改，没有 commit/push。

验证边界：未进行真实 COMSOL 求解；未跑与本次参数流无关的全仓数值/绘图测试。
