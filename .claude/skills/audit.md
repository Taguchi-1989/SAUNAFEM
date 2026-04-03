---
name: audit
description: "Full project quality audit. Runs lint, type check, security scan, dead code detection, and test coverage. Use when asked to review, audit, or check project health."
---

# Audit Skill

## When to Use
- User asks to "audit", "review project", "check quality", or "health check"
- Before major releases or milestones
- Periodically for maintenance

## Execution Steps

### Step 1: Automated Tool Checks (parallel)

Run these commands in parallel:

```bash
# Lint
ruff check src/ 2>&1; echo "EXIT:$?"
```

```bash
# Format check
ruff format --check src/ 2>&1; echo "EXIT:$?"
```

```bash
# Security scan
bandit -r src/ -c pyproject.toml --quiet 2>&1; echo "EXIT:$?"
```

```bash
# Dead code
vulture src/ --min-confidence 80 2>&1; echo "EXIT:$?"
```

```bash
# Tests with coverage
python -m pytest tests/ -x -q --cov=src --cov-report=term-missing 2>&1; echo "EXIT:$?"
```

### Step 2: Sub-Agent Deep Review (if issues found)

Launch sub-agents in parallel for:

#### Agent 1: Code Quality & Bug Detection
Scan all src/ files for logic errors, security issues, missing encoding="utf-8", dead code.

#### Agent 2: Uncommitted Changes Review
Run git diff, check for bugs/regressions in changed files.

### Step 3: Consolidate Report

```
## Audit Report: SaunaFlow

| Check          | Status | Details           |
|----------------|--------|-------------------|
| Ruff lint      | OK/NG  | N errors          |
| Ruff format    | OK/NG  | N files           |
| Bandit security| OK/NG  | H:N M:N L:N      |
| Vulture        | OK/NG  | N unused items    |
| Tests          | OK/NG  | N/N passed        |
| Coverage       | XX%    | threshold: 60%    |

### Critical Issues (must fix)
- ...

### Recommendations
1. ...
```

### Step 4: Fix Critical Issues
- Fix all BUG-severity issues immediately
- Run `pytest tests/ -x -q` after each fix
- Confirm all tests pass before declaring audit complete
