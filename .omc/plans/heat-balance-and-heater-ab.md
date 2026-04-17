# Plan: Heat Balance Auto-Aggregation + Heater Model A/B Comparison

**Date:** 2026-04-16
**Scope:** 2 phases across ~12 files (new + modified)
**Estimated complexity:** MEDIUM
**Dependencies:** Phase 2 depends on Phase 1

---

## Context

OpenFOAM results currently show upper bench at 207C (should be ~80-100C). Before refining the mesh or tweaking wall h, we need:
1. Visibility into where energy goes (heat balance)
2. A controlled comparison of two heater models to isolate whether surface flux is the root cause

The project already has Jinja2 templates, a case builder, probe parsing, KPI calculation, batch runner, and reporting. This plan extends all of those to support heat balance diagnostics and A/B heater comparison.

---

## Phase 1: Heat Balance Auto-Aggregation

### Objective
Add OpenFOAM functionObjects that automatically output energy budget quantities every write step, parse those outputs in the harness, and include a heat balance summary in reporting.

### Step 1.1: Add functionObjects to controlDict.j2

**File:** `foam_templates/base_case/system/controlDict.j2`

Add three functionObjects inside the existing `functions { }` block, after `fieldAverage`:

```
    volAverageT
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        writeFields     false;
        log             true;
        operation       volAverage;
        fields          (T);
    }

    heaterHeatFlux
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        writeFields     false;
        log             true;
        operation       areaIntegrate;
        regionType      patch;
        name            heater_wall;
        fields          (wallHeatFlux);     // requires wallHeatFlux functionObject or phi-based
    }
```

**However**, `wallHeatFlux` is not a default field -- it requires either:
- A separate `wallHeatFlux` functionObject to compute it, OR
- Using `phi` (mass flux) and enthalpy to compute wall heat transfer

**Recommended approach:** Add a `wallHeatFlux` functionObject first, then integrate over patches:

```
    wallHeatFlux
    {
        type            wallHeatFlux;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        patches         (floor ceiling heater_wall heater_wall_surround
                         opposite_wall front back);
        // Writes wallHeatFlux field + per-patch integrated values to postProcessing/
    }

    volAverageT
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        writeFields     false;
        log             true;
        operation       volAverage;
        fields          (T);
    }
```

The `wallHeatFlux` functionObject (OpenFOAM v2312+) automatically:
- Computes `wallHeatFlux` field on all wall patches
- Writes per-patch integrated heat flux to `postProcessing/wallHeatFlux/0/surfaceFieldValue.dat`
- Positive = heat into domain, negative = heat out of domain

For ventilation cases, add conditional surfaceFieldValue for vent patches:
```
{% if ventilation | default(false) %}
    ventHeatLoss
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        writeFields     false;
        log             true;
        operation       sum;
        regionType      patch;
        name            exhaust_vent;
        fields          (phi);
    }
{% endif %}
```

**Acceptance criteria:**
- `controlDict.j2` renders with `wallHeatFlux` and `volAverageT` functionObjects
- Generated case directory contains these in `system/controlDict`
- OpenFOAM run produces `postProcessing/wallHeatFlux/` and `postProcessing/volAverageT/` directories

### Step 1.2: Add heat balance parser to harness

**New file:** `src/harness/heat_balance_parser.py`

Create a parser that reads the postProcessing output from the functionObjects above:

- `parse_wall_heat_flux(case_dir: Path) -> dict[str, list[tuple[float, float]]]`
  - Reads `postProcessing/wallHeatFlux/<startTime>/surfaceFieldValue.dat`
  - Returns dict mapping patch name to list of (time, integrated_flux_W) tuples
  - Handles OpenFOAM's wallHeatFlux output format (tab-separated, header with #)

- `parse_vol_average_t(case_dir: Path) -> list[tuple[float, float]]`
  - Reads `postProcessing/volAverageT/<startTime>/volFieldValue.dat`
  - Returns list of (time, volume_averaged_T) tuples

- `compute_heat_balance(wall_fluxes: dict, heater_input: float) -> HeatBalance`
  - Dataclass `HeatBalance` with fields:
    - `heater_input_W: float` (from heater_wall patch, should be positive)
    - `wall_loss_W: float` (sum of all non-heater wall patches, should be negative)
    - `vent_loss_W: float` (optional, from vent patches)
    - `imbalance_W: float` (heater + wall_loss + vent_loss; ~0 at steady state)
    - `imbalance_pct: float` (imbalance / heater_input * 100)

**Acceptance criteria:**
- Parser reads sample postProcessing output (test fixtures)
- `HeatBalance` dataclass computed correctly from mock data
- Handles missing ventilation gracefully

### Step 1.3: Integrate heat balance into reporting

**Modify:** `src/harness/reporting.py`

Add a `heat_balance_to_markdown(balance: HeatBalance) -> str` function that generates a table:

```markdown
## Heat Balance Summary

| Component | Value [W] | % of Input |
|-----------|-----------|------------|
| Heater input | +18000 | 100.0% |
| Wall losses | -17500 | 97.2% |
| Vent losses | -200 | 1.1% |
| **Imbalance** | **+300** | **1.7%** |

Volume-averaged T: 358.2 K (85.1 C)
```

**Acceptance criteria:**
- Markdown table renders correctly
- Percentages are relative to heater input
- Volume-averaged T shown in both K and C

### Step 1.4: Tests for Phase 1

**New file:** `tests/unit/test_heat_balance_parser.py`

Tests:
- Parse wallHeatFlux output with known values
- Parse volAverageT output with known values
- Compute heat balance from mock patch data
- Handle missing vent patches
- Handle empty / malformed postProcessing data

**Modify:** `tests/unit/test_reporting.py` -- add test for `heat_balance_to_markdown`

**Acceptance criteria:**
- All new tests pass
- `pytest tests/ -x -q` still green

---

## Phase 2: Heater Model A/B Comparison

### Objective
Support two heater models (surface flux vs volume source) selectable from YAML config, create variant case YAMLs, and run side-by-side comparison.

### Step 2.1: Extend YAML schema for heater model type

**Modify:** `configs/schemas/case_schema.json`

Add to `boundary_conditions.heater.properties`:
```json
"model": {
  "type": "string",
  "enum": ["surface_flux", "volume_source"],
  "default": "surface_flux",
  "description": "Heater model: surface_flux (externalWallHeatFluxTemperature on patch) or volume_source (fvOptions scalarSemiImplicitSource in cellZone)"
}
```

Also add optional `depth` property to heater (needed for volume source cellZone):
```json
"depth": {
  "type": "number",
  "exclusiveMinimum": 0,
  "default": 0.3,
  "description": "Heater depth in x-direction [m] (used for volume_source model cellZone)"
}
```

**Acceptance criteria:**
- Schema validates both `surface_flux` and `volume_source` heater models
- Existing YAML files still validate (default is `surface_flux`)

### Step 2.2: Add volume source support to templates

**Modify:** `foam_templates/base_case/constant/fvOptions.j2`

Add conditional block for volume heater source:
```
{% if heater_model == "volume_source" %}
heaterSource
{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   cellZone;
    cellZone        heaterZone;

    scalarSemiImplicitSourceCoeffs
    {
        selectionMode   cellZone;
        cellZone        heaterZone;
        volumeMode      absolute;    // total power, not per-unit-volume

        sources
        {
            h           ({{ heater_power_W }} 0);   // enthalpy source [W]
        }
    }
}
{% endif %}
```

Note: `buoyantPimpleFoam` uses enthalpy `h`, not temperature `T`, for the energy equation. The source must be in `h` (J/s = W).

**Modify:** `foam_templates/base_case/0/T.j2`

Make heater_wall BC conditional on heater model:
```
    heater_wall
    {
{% if heater_model == "volume_source" %}
        // Volume source model: heater wall is just a regular wall
        type            externalWallHeatFluxTemperature;
        mode            coefficient;
        h               uniform {{ wall_htc | default(1.0) }};
        Ta              uniform {{ T_walls }};
        kappaMethod     fluidThermo;
        value           uniform {{ T_wall_inner | default(T_walls) }};
{% else %}
        type            externalWallHeatFluxTemperature;
        mode            flux;
        q               uniform {{ heat_flux }};
        kappaMethod     fluidThermo;
        value           uniform {{ T_heater | default(473.15) }};
{% endif %}
    }
```

**New file:** `foam_templates/base_case/system/topoSetDict.j2`

Create a topoSetDict template for defining the heaterZone cellZone:
```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}

actions
(
{% if heater_model == "volume_source" %}
    {
        name    heaterZone;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (0 {{ heater_y0 }} {{ heater_z0 }}) ({{ heater_depth }} {{ heater_y1 }} {{ heater_z1 }});
    }
{% endif %}
);
```

**Acceptance criteria:**
- `fvOptions` file contains `heaterSource` block when `heater_model == "volume_source"`
- `T.j2` renders heater_wall as regular wall for volume source, as flux BC for surface flux
- `topoSetDict.j2` defines `heaterZone` box correctly matching heater geometry
- Existing cases (no `heater.model` key) default to surface flux behavior

### Step 2.3: Update case_builder.py

**Modify:** `src/harness/case_builder.py`

Changes to `build_case()`:
1. Read `heater.model` from YAML (default `"surface_flux"`)
2. Read `heater.depth` from YAML (default `0.3`)
3. Add to template context:
   - `heater_model`: `"surface_flux"` or `"volume_source"`
   - `heater_power_W`: `power_kw * 1000` (for fvOptions absolute mode)
   - `heater_depth`: depth in meters
4. Conditionally skip `topoSetDict.j2` rendering when model is `surface_flux`
5. For volume source model, do NOT skip `topoSetDict.j2`

Changes to `compute_heater_params()`:
- Add `heater_power_W` to returned dict

Update `scripts/run_openfoam_wsl.sh` (or create a new variant):
- For volume source cases, run `topoSet` before `blockMesh` (or after blockMesh, before solver)
- Execution order: `blockMesh` -> `topoSet` -> `buoyantPimpleFoam`

**Acceptance criteria:**
- `build_case()` produces correct output for both heater models
- Volume source case has `topoSetDict` with correct heaterZone box
- Surface flux case has no `topoSetDict` (or empty actions)
- Template context includes all new variables

### Step 2.4: Create A/B case YAML configs

**New file:** `configs/cases/dry_sauna_steady_surfflux.yaml`
- Copy of `dry_sauna_steady.yaml` with explicit `heater.model: surface_flux`
- Case name: `dry_sauna_steady_surfflux`

**New file:** `configs/cases/dry_sauna_steady_volsource.yaml`
- Copy with `heater.model: volume_source` and `heater.depth: 0.3`
- Case name: `dry_sauna_steady_volsource`

Both should have identical geometry, solver settings, and probes for fair comparison.

**Acceptance criteria:**
- Both YAML files pass schema validation
- `build_case()` succeeds for both
- Generated OpenFOAM directories differ only in heater-related files

### Step 2.5: Extend batch runner for OpenFOAM A/B comparison

**Modify:** `src/harness/batch.py`

The current batch runner only uses `simple_solver`. Add support for OpenFOAM result comparison:

- Add `BatchOpenFOAMResult` dataclass with:
  - `case_name: str`
  - `probe_values: dict[str, float]` (from probe_parser)
  - `heat_balance: HeatBalance | None` (from heat_balance_parser)
  - `kpis: list[KPIResult]`

- Add `compare_openfoam_results(case_dirs: list[Path], probe_names: list[str]) -> str`
  - Reads probe data and heat balance from each case directory
  - Returns Markdown comparison table showing:
    - Probe temperatures for each case
    - Heat balance for each case
    - KPI comparison

This function does NOT run OpenFOAM (that is manual via WSL). It only parses results from existing case directories.

**Acceptance criteria:**
- Can parse and compare results from two OpenFOAM case directories
- Comparison table clearly shows differences between A and B
- Works with missing heat balance data (graceful fallback)

### Step 2.6: Tests for Phase 2

**Modify:** `tests/unit/test_case_builder.py`
- Test `build_case` with `heater.model: volume_source`
- Test that `topoSetDict` is generated with correct box coordinates
- Test that `T` file uses coefficient BC for volume source
- Test that `fvOptions` contains `heaterSource` for volume source
- Test backward compatibility (no `heater.model` key defaults to surface flux)

**New file:** `tests/unit/test_heat_balance_parser.py` (if not already created in Phase 1)

**Modify:** `tests/unit/test_batch.py`
- Test `compare_openfoam_results` with mock directories

**Modify:** `tests/unit/test_schema.py`
- Test schema accepts `heater.model` and `heater.depth`
- Test schema rejects invalid heater model values

**Acceptance criteria:**
- All tests pass
- `pytest tests/ -x -q` returns 0

---

## Execution Order

```
Phase 1 (Heat Balance) — no external dependencies
  1.1  controlDict.j2 functionObjects
  1.2  heat_balance_parser.py (new)
  1.3  reporting.py integration
  1.4  Tests

Phase 2 (Heater A/B) — depends on Phase 1 for heat balance comparison
  2.1  Schema extension
  2.2  Template changes (fvOptions.j2, T.j2, topoSetDict.j2)
  2.3  case_builder.py updates
  2.4  A/B YAML configs
  2.5  Batch comparison extension
  2.6  Tests
```

## Files Summary

### New files
| File | Purpose |
|------|---------|
| `src/harness/heat_balance_parser.py` | Parse wallHeatFlux + volAverageT postProcessing |
| `foam_templates/base_case/system/topoSetDict.j2` | cellZone definition for volume source heater |
| `configs/cases/dry_sauna_steady_surfflux.yaml` | Case A: explicit surface flux |
| `configs/cases/dry_sauna_steady_volsource.yaml` | Case B: volume source |
| `tests/unit/test_heat_balance_parser.py` | Tests for heat balance parsing |

### Modified files
| File | Change |
|------|--------|
| `foam_templates/base_case/system/controlDict.j2` | Add wallHeatFlux + volAverageT functionObjects |
| `foam_templates/base_case/constant/fvOptions.j2` | Add conditional heaterSource block |
| `foam_templates/base_case/0/T.j2` | Conditional heater_wall BC by model type |
| `src/harness/case_builder.py` | heater_model context, topoSetDict skip logic |
| `src/harness/reporting.py` | heat_balance_to_markdown |
| `src/harness/batch.py` | OpenFOAM result comparison |
| `configs/schemas/case_schema.json` | heater.model + heater.depth properties |
| `tests/unit/test_case_builder.py` | Volume source build tests |
| `tests/unit/test_reporting.py` | Heat balance report test |
| `tests/unit/test_batch.py` | OpenFOAM comparison test |
| `tests/unit/test_schema.py` | Schema extension tests |
| `scripts/run_openfoam_wsl.sh` | Add topoSet step for volume source cases |

---

## Open Questions

See `.omc/plans/open-questions.md`

## Success Criteria

1. `pytest tests/ -x -q` passes with all new tests
2. `build_case()` succeeds for both surface flux and volume source YAML configs
3. Heat balance functionObjects render correctly in controlDict
4. Generated OpenFOAM case directories are structurally valid
5. Comparison report clearly shows energy budget and temperature differences between A/B cases
6. Existing tests remain green (backward compatibility)
