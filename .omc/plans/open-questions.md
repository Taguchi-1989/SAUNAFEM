## openfoam-template-alignment - 2026-04-05
- [ ] Initial turbulence field values (k=0.1, omega=5.0) are not documented in governing_equations.tex -- should they be added for reproducibility?
- [ ] delta_t default in case_builder.py is 0.1s but doc Section 6.1 says 0.05s -- should the code default be updated to match?
- [ ] Phase 2 multiComponent branch lacks `div(phi,Yi_h)` scheme in fvSchemes -- should this be added now (low risk) or deferred to Phase 2 work?
- [ ] Ventilation patches for Phase 2 will require significant blockMeshDict restructuring -- should a separate mesh template be created or should the existing one be extended with conditionals?

## heat-balance-and-heater-ab - 2026-04-16

- [ ] wallHeatFlux functionObject availability in OpenFOAM v2312 -- The `wallHeatFlux` type was introduced in ESI OpenFOAM v2006+. Need to confirm the exact WSL OpenFOAM version supports it. If not available, fall back to computing heat flux manually from `phi` and `h` fields using `surfaceFieldValue` with `areaIntegrate` on `phi*h`.
- [ ] fvOptions `h` vs `e` for enthalpy source -- `buoyantPimpleFoam` solves for enthalpy `h` in OpenFOAM v2312. The `scalarSemiImplicitSource` field name should be `h` (not `e` or `T`). Verify by checking the solver's energy equation field name in the specific installed version.
- [ ] topoSet execution in WSL pipeline -- Currently `run_openfoam_wsl.sh` runs only `buoyantPimpleFoam`. For volume source cases, `topoSet` must run after `blockMesh` and before the solver. Decision needed: modify the existing script with a conditional flag, or create a separate script for volume source cases?
- [ ] heaterZone box precision -- The `boxToCell` source in topoSetDict selects cells whose centers fall inside the box. If the heater depth (0.3m) doesn't align well with cell boundaries at M0 resolution (cell size ~0.125m in x), the selected zone may not match the intended volume exactly. This is acceptable for A/B comparison purposes but should be documented.
- [ ] Volume source power normalization -- `scalarSemiImplicitSource` with `volumeMode absolute` applies the total power to all cells in the zone. With `volumeMode specific`, it applies per unit volume. Using `absolute` is simpler and matches the 18 kW specification directly. Confirm this is the desired approach.
