# 本机提交与 push

仅在用户要求提交/推送时读取；普通编辑的状态检查见 [AGENTS](../AGENTS.md)。

用户明确要求“git push 当前版本”时，在新对话中也按以下流程执行：

1. Git 仓库目录固定为 `I:\codeXproject\comsol-workflow`。先运行
   `git status -sb`、`git diff --stat` 和必要的具体 diff，确认当前分支、改动范围及
   是否混入无关文件。用户说“当前版本”并且工作树是一套连续、已确认的项目改动时，
   可以整体提交；存在来源不明或明显无关的修改时，必须先向用户确认，不得直接
   `git add -A`。
2. 提交前按改动风险运行测试、`git diff --check`，并清理本轮测试临时目录。大范围
   当前版本优先运行根项目完整 pytest；包含 `blueprints/` 时，同时使用其独立环境运行
   blueprint 测试。若一次测试因 Windows 文件锁或权限瞬态失败，应核对原因并重跑，
   不得把环境失败描述为代码失败。
3. 保持当前非默认分支，不因 push 自动新建分支。写入 `.git/index` 的 `git add`、
   `git commit`、`git commit --amend` 在 Codex sandbox 中可能报
   `index.lock: Permission denied`；遇到该环境时应直接使用
   `sandbox_permissions="require_escalated"` 执行，不要反复在只读 sandbox 中重试。
4. 使用简洁、覆盖完整 diff 的 commit message。提交后确认工作树干净且分支仅领先
   预期提交；若 `git diff --cached --check` 发现 whitespace 问题，先清理并 amend，
   再 push。
5. 普通 push 不依赖 GitHub CLI `gh`。即使本机没有安装 `gh`，也应继续使用原生 Git：

   ```powershell
   $branch = git branch --show-current
   git push -u origin $branch
   ```

   push 需要网络与 Git 元数据写权限时，使用
   `sandbox_permissions="require_escalated"`；可采用已批准的
   `prefix_rule=["git", "push"]`。只有用户明确要求创建或操作 PR 时才需要 `gh` 或
   GitHub connector；单纯 push 不创建 PR。
6. 推送完成后运行 `git status -sb` 与 `git log -1 --oneline --decorate`，确认
   `HEAD`、当前分支和 `origin/<branch>` 指向同一提交，并向用户报告分支、commit、
   测试结果和远端同步状态。
