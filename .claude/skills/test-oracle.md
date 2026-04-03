---
name: test-oracle
description: "Run test suite as quality gate. Validates all tests pass before any commit or major change. The test oracle is the single source of truth for project health."
---

# Test Oracle Skill

## When to Use
- Before every git commit
- After implementing any new feature
- After modifying any existing code
- When validating that a change didn't regress anything

## Execution Steps

1. **Run full test suite**:
   ```bash
   pytest tests/ -x -q
   ```
   - `-x` stops on first failure (fail fast)
   - `-q` for concise output

2. **If tests fail**:
   - Read the failure output carefully
   - Fix the root cause (do NOT skip or weaken tests)
   - Re-run until all pass
   - Record what failed and why in CHANGELOG.md if it was a non-obvious issue

3. **If all pass**:
   - Report count: "N tests passed"
   - Proceed with commit or next task

4. **Coverage check** (optional, for major features):
   ```bash
   pytest tests/ --cov=src --cov-report=term-missing -q
   ```
   Coverage threshold: 60% minimum

## Success Criteria
- Exit code 0 from pytest
- No warnings that indicate real issues
- Test count should only increase, never decrease
