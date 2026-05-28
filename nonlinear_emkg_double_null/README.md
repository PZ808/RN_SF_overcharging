# NonlinearEMKGDoubleNull

Julia scaffold for a spherically symmetric Einstein-Maxwell-charged-scalar code.

The near-term target is to combine:

- Baake & Rinne, arXiv:1610.08352v2: fully coupled Einstein-Maxwell-Klein-Gordon equations.
- Gelles & Pretorius, arXiv:2503.04881v2: compactified double-null MRT coordinates.
- Gelles & Pretorius, arXiv:2602.11256v1: fully nonlinear charged scalar evolution in double-null MRT gauge.

This first scaffold implements the coordinate machinery, initial data, a fixed-background charged scalar electrodynamics update, and a first nonlinear MRT-style state/update for a charged complex scalar with dynamic `r`, `f`, `A_u`, `A_v`, and enclosed charge `Q`.

The nonlinear charged module is under validation. The main metric convention is fixed in `src/NonlinearEMKG.jl`:

```text
ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2
MRT branch:    Phi = sqrt(32*pi) * phi_GP
GP2026 branch: Psi = r*Phi = sqrt(32*pi) * r*phi_GP
```

Here `phi_GP` denotes the canonically normalized scalar used in
arXiv:2602.11256. The GP2026 driver evolves the paper's reduced field
`r*phi` with `reduced_scalar=true` and converts back to `Phi` only when
computing stress-energy and Maxwell current sources. Conservative charge and
mass-balance diagnostics are now available. A GP2026 paper-style driver also
marches complete fixed-`U` rows with `Delta U=C/f_GP(U,Vmax)` and the
published hyperbolic `Q` equation, with the older constraint march retained
as a comparison mode. A diagnostic `max-row` step limiter is available because
the outer-boundary rule can miss an interior large-`f` peak after apparent
horizon formation. The next important work is controlling the evolving-metric
long-run behavior before horizon-accumulation validation.

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
julia --project=. examples/check_nonlinear_charged_balance.jl
julia --project=. examples/check_gp2026_initial_data.jl
julia --project=. examples/check_gp2026_short_balance.jl
julia --project=. examples/check_gp2026_long_evolution.jl
julia --project=. examples/check_gp2026_u_refinement.jl
julia --project=. examples/diagnose_gp2026_nonfinite.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```
