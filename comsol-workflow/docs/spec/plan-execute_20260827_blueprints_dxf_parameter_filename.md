# Blueprints DXF 参数化文件名执行计划

日期：2026-08-27
状态：已完成

## 目标

将 DXF 文件从通用的 `engineering_preview.dxf` / `production.dxf` 改为包含关键设计与
工艺参数的稳定文件名。当前参数的目标格式为：

```text
dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_engineering-unassigned-r0_8076bc19.dxf
```

## 命名字段

```text
dxf_<cavity layers>-<cladding layers>
_cav(<b0 nm>-<eta>-<zeta>)
_clad(<b0 nm>-<eta>-<zeta>)
_shift(<x factor>-<y factor>)
_<fabrication profile>-<dose kind>-r<revision>_<identity[:8]>.dxf
```

- 浮点数最多保留 6 位小数并去除无意义尾零；小数点按用户示例保留；
- 负号转换为 `m`，避免负数符号与字段分隔符混淆；
- 非数字工艺字段继续通过安全 label 规范化；profile 已以 `-<dose kind>` 结尾时不重复追加 dose；
- identity 短哈希保证未直接展示但会影响输出的参数仍不会发生文件名冲突；
- case 目录复用相同设计参数前缀并以 identity 短哈希结尾；staging 保持短名，完整实际路径控制在 Windows 限制内；
- `preview.png` 和 `layout_report.json` 保持精简固定命名。

release gate 仍由严格配置校验控制，不依赖文件名中是否出现 `production`。报告中的
`release_status` 继续明确区分工程预览与生产发布。

## 实施与验证

- 在 reporting 模块建立唯一 `dxf_artifact_name()`，CLI 与测试共同调用；
- 更新三文件输出合同、README、AGENTS 和相关测试；
- 测试精确锁定当前参数的示例文件名，并验证工程/发布配置得到不同名称；
- 运行完整独立 pytest、`compileall`、preflight 和 `git diff --check`；
- 删除当前旧命名 case 后，以同一参数重新生成并只保留新命名结果；
- 校验 case 仅含三个文件、DXF 回读与报告哈希一致，不启动 COMSOL。

## 回滚

恢复 CLI 的固定 DXF 名称及对应测试/文档，并重新生成 case。按本次请求删除的旧命名
结果不单独保留备份。

## 执行记录

### 实际改动

- 新增 `dxf_artifact_name()`，从 layers、cavity/cladding 紧凑参数、shift、fabrication、dose、revision 和 identity 构造 DXF 名；
- 当前参数精确生成
  `dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_engineering-unassigned-r0_8076bc19.dxf`；
- case 同步改为
  `dxf_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_shift(0-0)_8076bc19`；
- staging 保持短名；测试通过 extended-length path 检查超长 pytest 临时路径，实际输出完整路径为 222 字符；
- profile/dose 相邻重复已去除，发布状态继续由严格配置和报告控制；
- CLI、三文件合同、README、AGENTS 和测试已同步；旧固定 DXF 名不再使用。

### 验证结果

- 聚焦 CLI 测试：9 passed；独立完整测试：65 passed；
- `compileall`、preflight 和 `git diff --check`：通过；
- 新 case：29,648 features，DRC 0 error、0 warning；
- DXF：29,648 entities、89,478 vertices、micrometers、唯一 `ETCH` layer；
- DXF/PNG SHA-256 与 `layout_report.json` 一致；case 严格只有三个文件；
- 旧命名 case 在新 case 验证后删除，不保留重复结果；
- COMSOL 未启动。

最终状态：参数化 DXF 与 case 目录命名均已成为当前正式输出合同，当前产物仍为 engineering preview。
