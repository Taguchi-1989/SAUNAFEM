---
name: implement-feature
description: "Structured feature implementation following test-driven development. Plan → Test → Implement → Validate → Document."
---

# Implement Feature Skill

## Workflow

### Phase 1: Plan (before writing any code)
1. Read CLAUDE.md for project rules and architecture
2. Read CHANGELOG.md for known limitations and failed approaches
3. Identify which files need to change
4. Identify what tests need to be written
5. Check: does this conflict with any known limitation?

### Phase 2: Write Tests First
1. Write failing tests that define the expected behavior
2. Run test oracle — confirm new tests fail (they should, feature doesn't exist yet)
3. Confirm existing tests still pass (new tests should only ADD failures)

### Phase 3: Implement
1. Make the minimal changes to pass the new tests
2. Follow project conventions:
   - Absolute imports from src/
   - UTF-8 encoding on all file opens
   - Type hints where practical
3. Run test oracle after each meaningful change

### Phase 4: Validate (Ralph Loop)
1. Run full test suite: `pytest tests/ -x -q`
2. If failures: fix and retry (up to 10 iterations)
3. If all pass: proceed to Phase 5

### Phase 5: Document
1. Update CHANGELOG.md with:
   - What was added
   - Any architectural decisions made
   - Any failed approaches (if iterations > 1)
2. Commit with descriptive message

## Subagent Delegation
For complex features, spawn subagents for independent subtasks:

| Subtask | Agent Type | When |
|---------|-----------|------|
| Research existing code | Explore | Before implementation |
| Write new module | general-purpose (worktree) | Independent module |
| Run tests | general-purpose | After implementation |

## Anti-Patterns to Avoid
- Writing code without tests first
- Committing without running test oracle
- Retrying the same failed approach
- Skipping CHANGELOG.md updates for non-trivial changes
- Over-engineering (make it work, then make it right)
- Hand-editing generated OpenFOAM case files (always go through case_builder)
