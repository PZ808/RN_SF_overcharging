# Formulation Notes

## References

Baake/Rinne `arXiv:1610.08352v2` evolves the fully coupled Einstein-Maxwell-Klein-Gordon system in spherical symmetry on CMC hyperboloidal slices. The paper gives:

- conformal metric variables `Omega, Ntilde, X, Pi`;
- CMC/isotropic gauge equations;
- Maxwell variables `(ar, Er)` in temporal gauge;
- conformally rescaled scalar `phitilde` and `psitilde`;
- source terms for the Einstein constraints and slicing equation.

Gelles/Pretorius `arXiv:2503.04881v2` evolves charged scalar electrodynamics on a fixed RN background using compactified double-null MRT coordinates. The code starts from this coordinate system because the requested end state is "Baake/Rinne physics in Gelles/Pretorius coordinates."

## Current Implemented System

There are now two implementation tracks.

### Fixed-background SED track

The fixed-background fields are the Gelles/Pretorius SED variables:

- `xi = Re(r phi)`
- `pi = Im(r phi)`
- `Au, Av`
- enclosed charge `Q`

The quasi-Lorenz gauge condition is used:

```text
del_mu(r^-2 A^mu) = 0
```

In double-null coordinates this reduces to:

```text
Au,V + Av,U = 0
```

The fixed-background track follows the Gelles/Pretorius potential sign
convention used in their Appendix D.

The cell update follows the Burko-Ori corner scheme described in Appendix D of `arXiv:2503.04881v2`.

### Nonlinear MRT charged-scalar scaffold

`src/NonlinearEMKG.jl` starts the fully nonlinear charged problem in the Murata-Reall-Tanahashi convention:

```text
ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2
```

Its evolved variables are:

- `r`
- `logf`
- `phi_re, phi_im`
- `Au, Av`
- enclosed charge `Q`

The scaffold uses gauge-covariant scalar derivatives in the metric sources:

```text
D_a phi = partial_a phi - i e A_a phi
```

and evolves the enclosed charge from the null Maxwell constraints:

```text
Q_u =  4 pi r^2 J_u
Q_v = -4 pi r^2 J_v
```

`src/StressEnergy.jl` now isolates the matter sources in this convention:

```text
T_uu = 2 |D_u phi|^2
T_vv = 2 |D_v phi|^2
T_uv = alpha^2 / f
T_theta theta = 4 r^2 Re[(D_u phi)^* D_v phi] / f + 2 r^2 alpha^2 / f^2
J_a = 2 e Im[phi^* D_a phi]
alpha = Q f / r^2
F_UV = A_V,U - A_U,V = alpha
```

These are the canonical complex-scalar and Maxwell stress components for
`ds^2 = -f du dv + r^2 dOmega^2`, before choosing any extra overall
Einstein-equation normalization. The helper functions
`outgoing_constraint_source` and `ingoing_constraint_source` provide the
slots used by the two null Raychaudhuri constraints, normalized so the
single-real-scalar uncharged limit matches the MRT `1/4 r f^-1 phi_a^2`
constraint source.

## Missing Validation Work

The next physics/code steps are:

- verify every factor of 2 and 4 pi against the chosen action normalization;
- derive and solve the nonlinear constraints on the two initial null legs;
- validate that the uncharged limit reproduces Murata-Reall-Tanahashi;
- validate that the fixed-metric small-amplitude limit reproduces Gelles/Pretorius;
- define apparent horizon, Bondi mass, and horizon charge diagnostics in the double-null gauge.

## Current Diagnostics

`examples/residuals_mrt_equations.jl` evaluates the MRT equations on exact
extreme RN data in the Appendix A coordinates. The printed Eq. (4), Eq. (5),
Eq. (6), and Raychaudhuri signs are internally consistent: the current-sign
residuals decrease under refinement, while the old flipped `logf` sign leaves
an O(1) residual.

`examples/check_nonlinear_vacuum_rn.jl` checks that the nonlinear Burko-Ori
cell update preserves exact vacuum RN over a short run. It uses the optional
RN-background defect subtraction in `evolve_nonlinear!`, which subtracts the
known local truncation error of the exact MRT RN solution without changing the
continuum equations. Use enough fixed-point iterations for the nonlinear cell
solve; otherwise the remaining fixed-point error dominates this diagnostic.
This is the first metric regression to run after changing signs, factors, or
initial normalizations.

`examples/check_uncharged_decay.jl` is the intended MRT uncharged decay
target. After fixing the MRT signs in the metric equation, the current
uniform-`U` grid is not stable enough at the long times needed for this
diagnostic, so treat it as a target rather than a pass. MRT Appendix C uses
Burko-Ori adaptive mesh refinement in the `U` direction and notes that the
effective `U` resolution needed for late-time apparent-horizon work can be
very large.

`examples/check_charged_horizon_density.jl` is the charged-sector target
from Gelles/Pretorius. For extremal `eQ0=0.6`, the expected late-time
horizon charge-density exponent is `1 - 2s = 0`, i.e. a plateau. The
current scaffold does not pass this check yet; it is kept as a research
diagnostic to drive the next round of charged-sector corrections.

The Baake/Rinne equations cannot be pasted directly because their variables and gauge are CMC hyperboloidal, not compactified double-null. They are still the right source for matter stress tensor, charge conventions, and comparison diagnostics.
