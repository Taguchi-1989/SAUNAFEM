---
description: "Worktree-based task execution: create worktree → implement → test → commit → merge → cleanup in one flow"
user_invocable: true
---

# Worktree Task Skill

Execute a complete feature/fix cycle in an isolated git worktree.

## Usage
`/worktree-task <ticket-or-description>`

## Steps

1. **Create worktree** — Branch name from ticket/description (e.g., `feat/yaml-schema-validation`)
   ```bash
   git worktree add ../<branch-name> -b <branch-name>
   ```

2. **Implement** — Work in the worktree directory. Follow all CLAUDE.md rules.

3. **Test** — Run `pytest tests/ -x -q` in the worktree. All tests must pass.
   - If tests fail: fix and re-run. Do NOT proceed until green.

4. **Commit** — `git status` first, then stage and commit with descriptive message.
   - If pre-commit hook fails: fix issues, create NEW commit (never amend).

5. **Merge to main** — Switch to main, merge the branch:
   ```bash
   cd <original-repo>
   git merge <branch-name>
   ```
   - Run tests again on main after merge to confirm no regression.

6. **Cleanup** — Remove the worktree and branch:
   ```bash
   git worktree remove ../<branch-name>
   git branch -d <branch-name>
   ```

## Rules
- NEVER work directly on main — always use a worktree branch
- One worktree per task — don't mix unrelated changes
- If merge conflicts occur, resolve them carefully and re-test
- Report final status: tests passed count, files changed, merge result
