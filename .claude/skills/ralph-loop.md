---
name: ralph-loop
description: "Iterative completion loop. Keeps working on a task until success criteria are met, up to max iterations. Prevents 'agentic laziness' by re-checking claimed completion."
---

# Ralph Loop Skill

## Concept
Named after the pattern from Anthropic's long-running Claude research.
When a task has clear success criteria, iterate until those criteria are met.
Do not declare "done" until the oracle confirms it.

## When to Use
- Implementing a feature that must pass specific tests
- Debugging until a pipeline produces correct output
- Refactoring until all tests pass again
- Any task with measurable completion criteria

## Execution Pattern

```
for iteration in 1..max_iterations:
    1. Attempt the task (write code, fix bug, etc.)
    2. Run test oracle (pytest tests/ -x -q)
    3. If all success criteria met → DONE, break
    4. If not:
       a. Analyze what failed
       b. Record failed approach in CHANGELOG.md
       c. Adjust approach
       d. Continue to next iteration
```

## Parameters
- **goal**: What we're trying to achieve (e.g., "KPI calculation passes for all probe points")
- **max_iterations**: Maximum attempts before escalating to user (default: 10)
- **oracle_command**: How to check success (default: `pytest tests/ -x -q`)

## Rules
1. Never skip the oracle check
2. Never weaken tests to make them pass
3. Record each failed approach — don't retry the same thing
4. If hitting max_iterations, report what was tried and what remains
5. Each iteration should try a meaningfully different approach

## Example Usage
```
Goal: Implement YAML schema validation for case definitions
Max iterations: 10
Oracle: pytest tests/ -x -q (all tests pass + new schema tests pass)

Iteration 1: Used jsonschema for validation → test_missing_field FAILED (wrong error message format)
Iteration 2: Fixed error reporting → test_missing_field PASSED, test_range_check FAILED
Iteration 3: Added range validation → ALL TESTS PASSED → DONE
```
