## openfoam-template-alignment - 2026-04-05
- [ ] Initial turbulence field values (k=0.1, omega=5.0) are not documented in governing_equations.tex -- should they be added for reproducibility?
- [ ] delta_t default in case_builder.py is 0.1s but doc Section 6.1 says 0.05s -- should the code default be updated to match?
- [ ] Phase 2 multiComponent branch lacks `div(phi,Yi_h)` scheme in fvSchemes -- should this be added now (low risk) or deferred to Phase 2 work?
- [ ] Ventilation patches for Phase 2 will require significant blockMeshDict restructuring -- should a separate mesh template be created or should the existing one be extended with conditionals?
