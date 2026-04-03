---
name: observe
description: "Multi-perspective observation agent. Dispatches specialized sub-agents that each audit from a distinct viewpoint: physical model fidelity, data pipeline integrity, case reproducibility, edge-case robustness, and cross-module contract violations."
---

# Observe Skill (Multi-Perspective Audit)

## Philosophy
Standard audits scan line-by-line for code defects. This skill operates at a higher level: it launches sub-agents that each adopt a **distinct perspective** and trace concerns *across module boundaries*.

## When to Use
- After standard audit has been run and fixes applied
- When "there must be more bugs" -- need fresh angles
- Before milestone validation or experimental comparison

## Perspectives (Sub-Agents)

### 1. Physical Model Fidelity Auditor
**Question:** "Does this implementation faithfully represent the physical phenomena?"
```
Focus areas:
- Boundary conditions: Are wall heat fluxes, heater output, steam source rates physically reasonable?
- Unit consistency: All temperatures in K or °C? All velocities in m/s? All times in seconds?
- Conservation: Does the mesh + solver setup conserve energy/mass within tolerance?
- Probe placement: Are probe coordinates inside the mesh domain?
- Mesh sensitivity: Are results at current mesh level trustworthy or under-resolved?
- Solver settings: Are relaxation factors, time steps, turbulence models appropriate?
Report deviations from physical expectations.
```

### 2. Data Pipeline Tracer
**Question:** "Can data flow through every path without corruption or loss?"
```
Trace these end-to-end pipelines:
- YAML case definition → schema validation → OpenFOAM case directory → solver execution → probe output → KPI
- Experimental CSV → timestamp alignment → comparison with CFD probes → validation report
- Multiple cases → batch comparison → differential KPI table → report
For each pipeline: feed a non-trivial case and trace what gets lost, corrupted,
or misinterpreted at each transformation boundary.
```

### 3. Edge Case & Boundary Tester
**Question:** "What happens at the extremes?"
```
Test these specific scenarios by reading the code (NOT running it):
- Case YAML with missing optional fields
- Case YAML with invalid geometry (negative dimensions, overlapping objects)
- Probe point outside domain
- Zero heat output
- Löyly with duration=0 or vapor_rate=0
- Aufguss with zero velocity
- Empty experimental CSV
- Experimental CSV with different timestamp resolution than CFD
- Mesh level M0 (very coarse) — do KPIs still compute?
- Solver that doesn't converge — how does harness handle it?
For each: trace the code path and report what actually happens.
```

### 4. Cross-Module Contract Verifier
**Question:** "Do modules agree on their shared interfaces?"
```
Check these interfaces for contract mismatches:
- YAML schema ↔ case_builder: Do all YAML fields get consumed? Are any ignored?
- case_builder ↔ OpenFOAM templates: Do template placeholders match builder output?
- solver_runner ↔ probe_parser: Does probe output format match parser expectations?
- probe_parser ↔ kpi.py: Do parsed data structures match KPI calculation inputs?
- kpi.py ↔ validation.py: Do KPI formats match validation comparison inputs?
- validation.py ↔ reporting.py: Do validation results match report generation expectations?
For each interface: list the "contract" vs "expectation" and flag mismatches.
```

### 5. Reproducibility & Determinism Verifier
**Question:** "Can every result be reproduced from the same inputs?"
```
Check:
- Are all random seeds fixed or controlled?
- Does case_builder produce identical output for identical YAML input?
- Are OpenFOAM solver settings fully specified (no reliance on defaults that may change)?
- Are probe output paths deterministic (no timestamp in directory names)?
- Can a result directory be traced back to its exact input YAML?
- Are Python package versions pinned in pyproject.toml?
- Is the OpenFOAM version recorded in result metadata?
Report any source of non-determinism.
```

## Execution

Launch all 5 sub-agents in parallel. Each should:
1. Read the relevant source files (not just scan -- understand the logic)
2. Report only **real issues** (not style, not theoretical)
3. Include file:line references
4. Classify as: BUG (wrong behavior), GAP (missing feature), MISMATCH (contract violation)
