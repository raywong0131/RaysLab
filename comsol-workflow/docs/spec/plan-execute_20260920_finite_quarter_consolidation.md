# finite-quarter 三套结果无损整合

## 目标与边界

将用户给出的四项输入按规范化绝对路径去重后，把三套唯一 finite-quarter 结果复制到：

`scripts/.out/finite_cavity/finite_quarter_all_20-20_cav(245-0.96-1.156)_clad(242-0.98-0.93)_mesh5`

其中 `f202.500_e6_20260908` 在原始输入中重复两次，只处理一次。其余两个唯一来源为
`f203.600_e6_20260906_retry1` 和 `f197.980_e15_20260903`。不启动 COMSOL，不移动、
删除或改写任一源目录，也不按频率相近程度删除物理模式。

“剔除重复数据”在本次整合中指删除重复输入路径；不同求解运行的结果即使模式频率相近，也保留其
独立来源和全部文件，避免误删简并或相邻本征模。

## 目标布局

目标沿用 finite 结果的两位数字目录约定，并用运行标签隔离同名文件：

```text
目标目录/
  00_model/<run_tag>/[<symmetry_case>/]
  01_results/<run_tag>/[<symmetry_case>/]
  10_overview/
    <旧 overview 的 19 个根级汇总文件>
    source_runs/<run_tag>/[<symmetry_case>/]
  80_logs/<run_tag>/[<symmetry_case>/]
  99_config/
    integration_manifest.json
    source_runs/<run_tag>/[<symmetry_case>/]
```

- 两套 all-mode 运行的 case 级标准目录分别映射到对应 `run_tag`。
- 旧 `fq4m_retry` 运行的四个 `symmetry_*` case 保持各自层级。
- 旧运行根级 `overview` 直接映射为新根级 `10_overview`，不额外复制一份。
- 旧运行根级 `four_mode_run_summary.json` 放入其来源标签对应的 `99_config/source_runs`。
- 空日志树不强制创建；没有文件的目录不属于需保留的数据。

## 安全与验证

1. 复制器拒绝已存在的正式目标和 staging，严格检查三套源的已知布局与标签唯一性。
2. 先在目标同级的短路径 staging 中复制，逐映射比较文件数与总字节数。
3. 写入 `99_config/integration_manifest.json`，记录四项原始输入、三项唯一来源、重复项、映射和源统计。
4. 所有映射合计必须等于三套唯一源的文件数与字节数；验证后同盘原子改名为正式目标。
5. 完成后复核源目录统计不变、目标文件总数、主要 CSV/JSON 可解析、PNG 可解码、无 staging 残留。
6. 运行复制器的 dry-run、`py_compile` 与本任务文件的 `git diff --check`。

## 预检统计

| 来源 | 文件数 | 字节数 |
| --- | ---: | ---: |
| `f202.500_e6_20260908` | 422 | 23,827,817,864 |
| `f203.600_e6_20260906_retry1` | 285 | 22,694,256,140 |
| `f197.980_e15_20260903` | 581 | 54,090,123,637 |
| 合计 | 1,288 | 100,612,197,641 |

正式目标预期为 1,289 个文件（源文件 1,288 个，加 1 个整合清单）；总字节数为源字节数加清单大小。

## 执行记录

- 首轮复制在 `f203.600_e6_20260906_retry1/80_logs` 遇到历史失败导出目录中的超长源路径，
  `shutil.copy2` 返回 `WinError 3`；正式目标尚未发布。普通 `Path.rglob()` 也漏计了该文件，解释了
  284/285 的预检差异。复制器随后为 Windows 的枚举、计数和复制统一增加 `\\?\` extended-path
  前缀；失败 staging 经路径边界核对后清理，再从头执行，避免保留半成品。
- 第二轮已通过超长路径复制，但因现代来源先创建 `10_overview/source_runs`，旧根 `overview` 随后不能
  再以整树方式创建其父目录。复制计划调整为先复制旧根 `overview`，再复制各来源 overview 子树；
  第二个未发布 staging 同样在边界核对后清理。
- 旧 `overview` 中 6 个历史更新文件的 ACL 拒绝沙箱身份读取；核验具体文件后，以提升权限运行同一
  受限复制器完成正式复制。27 个映射均通过逐树文件数/字节数校验，staging 最终原子改名为目标。
- 正式目标共 1,289 个文件、100,612,209,522 字节；三套源复核仍分别为 422/285/581 个文件，
  23,827,817,864 / 22,694,256,140 / 54,090,123,637 字节。
- 根级 `10_overview` 的 19 个直接文件与旧 `overview` 文件名和大小逐一一致；另外保留各来源的
  `source_runs` overview。目标内 77 个 JSON、334 个 CSV 和 659 个 PNG 均通过解析/解码检查，
  `integration_manifest.json` 状态为 `complete`，无 `.fq_integrate_staging_*` 残留。
- 未启动 COMSOL；长路径针对性复制、入口实际执行、Python 编译检查和 `git diff --check` 均通过。
