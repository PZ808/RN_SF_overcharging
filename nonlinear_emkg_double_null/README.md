# NonlinearEMKGDoubleNull

Julia scaffold for a spherically symmetric Einstein-Maxwell-charged-scalar code.

The near-term target is to combine:

- Baake & Rinne, arXiv:1610.08352v2: fully coupled Einstein-Maxwell-Klein-Gordon equations.
- Gelles & Pretorius, arXiv:2503.04881v2: compactified double-null MRT coordinates.

This first scaffold implements the coordinate machinery, initial data, a fixed-background charged scalar electrodynamics update, and a first nonlinear MRT-style state/update for a charged complex scalar with dynamic `r`, `f`, `A_u`, `A_v`, and enclosed charge `Q`.

The nonlinear charged module is scaffolding, not yet a validated research solver. The main convention is fixed in `src/NonlinearEMKG.jl`:

```text
ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2
```

The next important work is checking every factor against the chosen action normalization, then adding constraint solves on the two initial null legs and convergence tests.

## Layout

- `src/Coordinates.jl`: RN metric functions and compactified MRT coordinates.
- `src/Fields.jl`: parameters, state containers, current definitions, constraints.
- `src/InitialData.jl`: ingoing compact-pulse data on two null legs.
- `src/Evolution.jl`: second-order double-null finite difference update.
- `src/NonlinearEMKG.jl`: nonlinear MRT-style charged scalar scaffolding.
- `src/Diagnostics.jl`: Maxwell residuals and simple conserved-quantity checks.
- `examples/run_smoke.jl`: small fixed-background evolution.
- `docs/formulation_notes.md`: implementation notes and equation provenance.

## Run

```bash
cd nonlinear_emkg_double_null
julia --project=. examples/run_smoke.jl
julia --project=. examples/run_nonlinear_smoke.jl
julia --project=. examples/check_mrt_initial_normalization.jl
julia --project=. examples/check_electrovac_mass_conservation.jl 0.02 531 10.0 0.02
julia --project=. examples/check_uncharged_mass_flux.jl 0.02 20.0 0.02 531 fixed
julia --project=. examples/check_uncharged_mass_flux.jl 0.02 20.0 0.02 531 split
julia --project=. examples/check_uncharged_mass_flux.jl 0.02 20.0 0.02 531 split 0.005
julia --project=. examples/check_uncharged_bondi_mass.jl 0.02 0.02 0.02 150.0 0.1
julia --project=. examples/check_uncharged_bondi_mass.jl 0.02 0.02 0.02 150.0 0.02
julia --project=. examples/check_uncharged_bondi_mass.jl 0.02 0.02 0.02 150.0 0.02 0.000025 0.05
julia --project=. examples/check_uncharged_bondi_mass.jl 0.02 0.02 0.02 150.0 0.02 0.000025 0.05 531
julia --project=. -e 'using Pkg; Pkg.test()'
```
