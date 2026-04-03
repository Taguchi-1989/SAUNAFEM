---
name: deploy-check
description: "Pre-run comprehensive checklist. Runs all quality gates, validates case definitions, and confirms OpenFOAM environment readiness. Use before running a production simulation case."
---

# Deploy Check Skill

## When to Use
- Before running a full simulation case (especially M2+ mesh levels)
- Before experimental validation runs
- Before milestone deliverables
- User says "deploy check", "ready to run?", or "pre-flight check"

## Execution Steps

### Phase 1: Quality Gates (must all pass)

Run sequentially — stop on first failure:

```bash
echo "=== 1/5 Ruff Lint ==="
ruff check src/
if [ $? -ne 0 ]; then echo "BLOCKED: Fix lint errors"; exit 1; fi

echo "=== 2/5 Ruff Format ==="
ruff format --check src/
if [ $? -ne 0 ]; then echo "BLOCKED: Run ruff format src/"; exit 1; fi

echo "=== 3/5 Tests ==="
python -m pytest tests/ -x -q
if [ $? -ne 0 ]; then echo "BLOCKED: Fix failing tests"; exit 1; fi

echo "=== 4/5 Security ==="
bandit -r src/ -c pyproject.toml -ll
if [ $? -ne 0 ]; then echo "WARNING: Medium+ security issues found"; fi

echo "=== 5/5 Coverage ==="
python -m pytest tests/ -q --cov=src --cov-report=term-missing --cov-fail-under=60
if [ $? -ne 0 ]; then echo "BLOCKED: Coverage below 60%"; exit 1; fi

echo ""
echo "=== ALL GATES PASSED ==="
```

### Phase 2: Case Validation

```bash
# Validate the target case YAML
python -m harness.cli validate configs/cases/<target_case>.yaml

# Check OpenFOAM environment
which blockMesh && which buoyantSimpleFoam && echo "OpenFOAM OK" || echo "BLOCKED: OpenFOAM not found"
```

### Phase 3: Pre-Run Checklist

Report to user:

```
## Simulation Readiness Report

### Quality Gates
- [ ] Ruff lint: PASS (0 errors)
- [ ] Ruff format: PASS
- [ ] Tests: PASS (N/N)
- [ ] Security: PASS (0 medium+)
- [ ] Coverage: XX% (≥60%)

### Case Validation
- [ ] YAML schema valid
- [ ] Geometry dimensions reasonable
- [ ] Probe points within domain
- [ ] Mesh level appropriate for objective
- [ ] Solver settings match phase requirements

### Environment
- [ ] OpenFOAM installed and accessible
- [ ] Sufficient disk space for results
- [ ] Sufficient RAM for target mesh level
- [ ] Python dependencies up to date

### Manual Checks (user)
- [ ] Case YAML reviewed by human
- [ ] Expected run time is acceptable
- [ ] Results directory is clean or backed up
```

## Success Criteria
- All 5 quality gates pass
- Case YAML validates against schema
- OpenFOAM environment confirmed
- All checklist items addressed
