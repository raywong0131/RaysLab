# Issue tracker: GitHub

本项目的规格和实现票据使用 [jia-yli/comsol-workflow 的 GitHub Issues](https://github.com/jia-yli/comsol-workflow/issues)。
使用 `gh` CLI，命令在 `I:/codeXproject/comsol-workflow` 执行；显式指定
`--repo jia-yli/comsol-workflow`，避免 workspace 根目录或兄弟仓库导致目标错误。
不在 `.scratch/` 建立第二套正式 tracker。

## 访问与网络

- 首次使用先检查 `gh auth status --hostname github.com` 和
  `gh repo view jia-yli/comsol-workflow --json nameWithOwner,hasIssuesEnabled,viewerPermission`。
- 无认证时由用户运行 `gh auth login --hostname github.com --web --git-protocol https`，
  使用有目标仓库权限的账号；不要在对话或配置文件中粘贴 token。
- 本机已配置的系统代理为 `127.0.0.1:7897`。命令行未继承代理而连接失败时，
  在当前 PowerShell 会话设置 `$env:HTTPS_PROXY='http://127.0.0.1:7897'`；
  代理地址以当前主机实际设置为准，不改全局连接配置。
- 认证、权限或网络失败时明确报告阻塞；可以准备本地草稿，但不能称已发布，
  不擅自切换 tracker、公开私有仓库或改变仓库设置。

## 操作约定

以下 PowerShell 命令中的编号、标题、文件名替换为实际值：

```powershell
gh issue list --repo jia-yli/comsol-workflow --state open --json number,title,body,labels
gh issue view 42 --repo jia-yli/comsol-workflow --comments --json number,title,body,labels,comments,state
gh issue create --repo jia-yli/comsol-workflow --title '实际标题' --body-file '实际正文文件.md' --label ready-for-agent
gh issue comment 42 --repo jia-yli/comsol-workflow --body-file '实际评论文件.md'
gh issue edit 42 --repo jia-yli/comsol-workflow --remove-label needs-triage --add-label ready-for-agent
gh issue close 42 --repo jia-yli/comsol-workflow
```

- “publish to the issue tracker”表示创建 GitHub issue；“fetch ticket”表示读取 issue 正文、
  当前标签和评论。多行正文写 UTF-8 文件并用 `--body-file`，不拼接 shell 字符串。
- Spec 为一个 issue；实现票据一项一个 issue，链接原始 spec、相关 `docs/spec/` 记录和验收条件。
  按依赖顺序创建票据；`to-spec` / `to-tickets` 产出的已明确票据使用 `ready-for-agent`。
- 阻塞关系优先用 GitHub 原生 issue dependencies。添加时执行
  `gh api --method POST repos/jia-yli/comsol-workflow/issues/<child-number>/dependencies/blocked_by -F issue_id=<blocker-database-id>`；
  database id 由 `gh api repos/jia-yli/comsol-workflow/issues/<blocker-number> --jq .id` 获取，
  不是 issue number 或 node_id。接口不可用时正文写 `Blocked by: #...`，全部 blocker 关闭后才可实施。
- 现有 `docs/spec/plan-execute_*.md` 保留为仓库内方案和执行证据，不批量迁移为 issue，
  不因 setup 创建演示票据、评论或提交。

## Triage labels

标签名和语义见 [triage-labels.md](triage-labels.md)。首次认证后先用
`gh label list --repo jia-yli/comsol-workflow --limit 100 --json name,description,color`
核对，再用 `gh label create <name> --repo jia-yli/comsol-workflow --description <description>`
补齐缺少的五个默认标签；不使用 `--force` 覆盖已有标签。

## Pull requests as a triage surface

**PRs as a request surface: no.**

默认只 triage issues；仅在用户明确改变该约定后才纳入外部 PR。
`code-review` 仍可按用户指定范围审查 PR。发布规格/票据不意味着授权提交、push 或 merge。
