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

`src/AdaptiveMRT.jl` is the first step toward the MRT/Burko-Ori adaptive
scheme: it introduces fixed-`V` nonlinear slices with independent `U` grids,
linear prolongation, cell-flag refinement helpers, and a first
`advance_adaptive_slice` primitive. That primitive prolongates the previous
slice onto the chosen next `U` grid and reuses the existing Burko-Ori cell
update to fill the new fixed-`V` slice from the ingoing boundary point.
`evolve_adaptive` loops over prescribed ingoing-boundary data, retaining the
`U` grid when no refinement policy is supplied so it can be tested directly
against the rectangular solver. Its Burko-Ori-style point-splitting policy
checks completed fixed-width `V` bands for relative variation in `r` or `f`,
absolute complex-scalar variation `|Delta Phi|`, and changes in a cell-slope
estimate of `|Delta partial_U Phi|`. The scalar thresholds are absolute to
avoid singular behavior at scalar zero crossings; diagnostics set them as
fractions of the injected amplitude. Inserted points use four-point
polynomial interpolation where four neighboring points are available,
following Burko-Ori's point-splitting description. An optional horizon-chop
policy evaluates `r,V` from two consecutive slices and retains one cell inside
the first detected `r,V <= 0` crossing, avoiding unnecessary evolution deeper
inside a trapped region. A horizon-local refinement policy can additionally
bisect cells near that crossing until a requested `Delta U` is reached.
Combining that policy with chopping requires a coordinate-width interior
buffer: retaining only one cell would shrink the physical buffer when that
cell is split and can let the horizon move beyond the stored domain before the
next refinement band. Point removal remains to be added.

For extreme RN validation, the analytic MRT helpers include the regular
future-horizon endpoint `U=0`, where `r=1` and `f=2`. This lets decay checks
distinguish exact-horizon sampling from a nearby finite-`U` extraction.

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

`examples/check_uncharged_decay.jl` targets the MRT Fig. 13 decay experiment:
an outgoing pulse on `Sigma_1`, sampled at the extreme-RN horizon endpoint
`U=0`. This differs from the Appendix B ingoing-wave data supplied by
`initialize_mrt2013_uncharged_ingoing!`, whose nonlinear endpoint is
non-extreme RN and therefore is not a `1/V` check. The script accepts optional
dimensionless thresholds for geometry, `Delta Phi / amplitude`, and
`Delta partial_U Phi / amplitude`, making threshold-refinement comparisons
repeatable. MRT Appendix C uses Burko-Ori adaptive mesh refinement in the `U`
direction and notes that the effective `U` resolution needed for late-time
horizon work can be very large.

The outgoing-wave helper now implements MRT Eq. (24) directly, instead of an
earlier convenient smooth substitute. The previous numerical tail-threshold
scan predates that correction and must not be used as validation evidence; it
should be repeated after the Bondi-mass convergence work below.

`examples/check_uncharged_bondi_mass.jl` targets MRT Fig. 7 for scalar charge
`e=0` on the electromagnetic background `Q=1`. It tunes `f0` by imposing the
degenerate initial apparent horizon condition of Eq. (28), evaluates the
renormalized Hawking mass of Eq. (13), and approximates `M_B(U)` by the
profile at `V approximately 150` as in the paper. The diagnostic subtracts
the exact extreme-RN stencil defect from the numerically differentiated mass;
without this subtraction, the finite-difference background error is comparable
to the physical `O(epsilon^2)` mass excess.

The initial-data tuner reproduces MRT's small-amplitude relation. The seeded
outgoing leg now evaluates the known analytic derivative of the Eq. (24)
wavepacket when integrating Eq. (25), instead of taking finite differences of
the sampled pulse:

| `epsilon` | `(f0-2)/epsilon^2` | `(Mi-1)/epsilon^2` | MRT targets |
| ---: | ---: | ---: | --- |
| 0.02 | 0.74045 | 0.78919 | 0.740, 0.789 |
| 0.03 | 0.74049 | 0.78908 | 0.740, 0.789 |

For `Delta V=0.1` and adaptive thresholds `0.02`, the scaled mass profiles
for `epsilon=0.02` and `0.03` nearly coincide in the wavepacket region:

| `U` sample | `epsilon=0.02` | `epsilon=0.03` |
| ---: | ---: | ---: |
| about -5.0 | 0.78916 | 0.78905 |
| about -4.0 | 0.71447 | 0.71443 |
| about -3.0 | 0.32931 | 0.32921 |
| about -2.0 | 0.28480 | 0.28492 |
| about -0.05 | 0.03665 | 0.03620 |

The last row is approaching the Fig. 7 final-mass target
`epsilon^-2(Mf-1) approximately 0.789 * 0.0375 = 0.0296`. Extracting the
value directly at the numerically located apparent horizon is not converged
yet: with `epsilon=0.02`, halving the adaptive thresholds at `Delta V=0.1`
changes the inferred retained fraction from `0.0391` to `0.0269`.

Before correcting the discrete Eq. (25) seeding, enabling horizon chopping
did not alter the `Delta V=0.1`, threshold `0.02` observables to shown
precision, while reducing the final stored `U` points from 2387 to 1897.
Keeping the adaptive thresholds fixed at `0.02`, that time-step scan gave:

| `Delta V` | final `U` points | `U_AH(V approximately 150)` | `epsilon^-2(M_B(-0.05)-1)` | horizon retained fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0.10 | 1897 | -0.00497 | 0.03665 | 0.03915 |
| 0.05 | 1891 | -0.00534 | 0.03616 | 0.05072 |
| 0.02 | 1889 | -0.00543 | 0.03602 | 0.05396 |

This separated the exterior mass profile from the apparent-horizon extraction
error. Still using the pre-correction seed and holding `Delta V=0.02` fixed,
horizon-local point splitting gave:

| requested horizon `Delta U` | interior `U` buffer | attained horizon `Delta U` | final `U` points | `U_AH(V approximately 150)` | horizon retained fraction |
| ---: | ---: | ---: | ---: | ---: | ---: |
| disabled | 0.00 | 4.7484e-5 | 1889 | -0.00542952 | 0.05395993 |
| disabled | 0.05 | 4.7484e-5 | 2165 | -0.00542952 | 0.05395993 |
| 2.5e-5 | 0.05 | 2.3742e-5 | 2191 | -0.00542976 | 0.05395833 |
| 1.25e-5 | 0.05 | 1.1871e-5 | 2233 | -0.00542979 | 0.05395797 |

The buffer-only row verifies that retaining enough interior domain does not
change that observable. Halving the actual horizon-cell width twice produced
only an `O(2e-6)` change in the retained fraction.

`examples/check_mrt_initial_normalization.jl` now evaluates Eq. (28) directly
on the seeded `V=0` leg. MRT's continuum degenerate-horizon data have
`r_V(U=0,V=0)=0`, but the original 110-point initial grid does not even
contain `U=0`, and its sampled minimum is not close to zero:

| initial `U` points | initial `Delta U` | location of minimum `r_V` | `min(r_V)/epsilon^2` |
| ---: | ---: | ---: | ---: |
| 110 | 4.8624e-2 | 0.00550 | 0.53935 |
| 107 | 5.0000e-2 | 0.00000 | 0.51821 |
| 531 | 1.0000e-2 | 0.00000 | 0.020738 |
| 1061 | 5.0000e-3 | 0.00000 | 0.005185 |
| 5301 | 1.0000e-3 | 0.00000 | 0.000207 |

Thus the initial degenerate horizon is represented only to second-order
accuracy and must be included in convergence studies. Repeating the late-time
run after the Eq. (25) correction, with `Delta V=0.02`, interior buffer
`0.05`, and requested horizon `Delta U=2.5e-5`, gives:

| initial `U` points | attained horizon `Delta U` | `epsilon^-2(M_B(-0.05)-1)` | horizon retained fraction |
| ---: | ---: | ---: | ---: |
| 110 | 2.3742e-5 | 0.0356756 | 0.0534719 |
| 531 | 1.9531e-5 | 0.0338520 | 0.0538673 |
| 1061 | 1.9531e-5 | 0.0337517 | 0.0536832 |

The initial resolution materially changes the exterior mass profile, but it
does not bring the horizon-extracted retained fraction toward MRT's `0.0375`.
Moreover, the refined-initial-grid runs develop a small `M_B < 1` region and
a non-monotone late profile, inconsistent with the BPS bound and Eq. (14).

The evolved-mass checks now separate this failure further. In
`examples/check_electrovac_mass_conservation.jl`, scalar-free data with the
same `f0 > 2` must evolve as a nearby non-extreme RN spacetime with constant
mass `M_i`. At `V approximately 10`, `Delta V=0.02`, and 531 fixed `U`
points, the maximum mass error is:

| evolution defect correction | mass diagnostic correction | `max abs(varpi-M_i)` | relative to `M_i-1` |
| :---: | :---: | ---: | ---: |
| off | off | 2.75e-5 | 8.70e-2 |
| off | on | 1.94e-5 | 6.14e-2 |
| on | off | 3.82e-5 | 1.21e-1 |
| on | on | 8.45e-7 | 2.68e-3 |

Thus the paired extreme-RN stencil corrections improve, rather than corrupt,
the nearby electrovac mass in this short control.

`examples/check_uncharged_mass_flux.jl` evaluates MRT Eq. (14),
`varpi_U = -r^2 r_V |phi_U|^2/(2f)`, directly from the evolved scalar fields
and compares it with the derivative of the reconstructed mass. For
`epsilon=0.02`, `V approximately 20`, `Delta V=0.02`, and the corrected mass
diagnostic:

| initial `U` points | adaptive mode | final `U` points | max exterior `abs(flux residual)` |
| ---: | :--- | ---: | ---: |
| 531 | fixed | 531 | 1.55e-5 |
| 531 | chop only | 531 | 1.55e-5 |
| 531 | point splitting only | 948 | 1.58e-3 |
| 531 | point splitting plus chop | 821 | 1.58e-3 |
| 1061 | fixed | 1061 | 3.96e-6 |
| 1061 | point splitting only | 1404 | 1.45e-3 |

The fixed-grid residual reduces by approximately four when the initial
resolution doubles, as expected for the second-order stencil. In this
`V approximately 20` control, chopping alone leaves the grid unchanged and
cannot be the source of the displayed failure. Point insertion is necessary
for the large residual, remains two orders of magnitude worse, and does not
improve when starting from the finer grid. The next implementation target is
therefore constraint-preserving prolongation for refined fixed-`V` slices,
most directly by reconstructing `log(f)` from the `U` constraint after
interpolating `r` and `phi`.

Constraint-preserving prolongation has now been added. When a parent cell is
split, `r`, the scalar, potentials, and charge are midpoint-interpolated as
before, but the inserted `log(f)` is reconstructed from the `U` Raychaudhuri
constraint using both parent endpoints. Existing evolved grid values are not
altered; projecting all existing values from one boundary was tested and
discarded because it accumulated a larger slice-wide integration error.

At the original point-splitting threshold `0.02`, this local constraint
projection changes the corrected 531-point flux residual only from the
pre-projection value `1.58e-3` to `1.54e-3`. Threshold refinement gives:

| initial `U` points | point-splitting threshold | final `U` points | max exterior `abs(flux residual)` |
| ---: | ---: | ---: | ---: |
| 531 | fixed grid | 531 | 1.55e-5 |
| 531 | 0.020 | 948 | 1.54e-3 |
| 531 | 0.010 | 1479 | 5.23e-4 |
| 531 | 0.005 | 2876 | 1.45e-4 |

The projected prolongation is therefore retained as the consistent insertion
rule, but it is not the dominant repair: at the thresholds affordable so far,
adaptive evolutions still violate the flux diagnostic more strongly than a
fixed grid. Before returning to the `V=150` Fig. 7 comparison, the next
controlled step is to improve or further tighten the splitting criterion
around the observed residual peak (`U approximately -0.23` in this run), and
verify convergence of the flux residual.

The conservative diagnostic provides a more direct resolution of the Fig. 7
comparison. `uncharged_flux_integrated_mass_profile` anchors `varpi` in the
outer electrovac region and integrates MRT Eq. (14) inward:

`varpi(U) = varpi(U0) + integral[-r^2 r_V |phi_U|^2/(2f) dU]`.

For the same `epsilon=0.02`, `Delta V=0.02`, 531-point aligned initial grid,
point-splitting threshold `0.02`, horizon `Delta U <= 2.5e-5`, and
`V approximately 150` evolution, the two extractions give:

| extraction | `epsilon^-2(M_f-1)` | `(M_f-1)/(M_i-1)` | maximum exterior upward mass step |
| :--- | ---: | ---: | ---: |
| geometric Eq. (13) from differentiated metric data | 0.0428763 | 0.0543293 | 2.98e-5 |
| flux-integrated Eq. (14) | 0.0295846 | 0.0374871 | 0 |
| MRT Fig. 7 target | about 0.0296 | about 0.0375 | monotone |

Thus the evolved scalar/metric solution carries the correct integrated mass
loss for this validation run; the prior Fig. 7 failure came from evaluating
the derivative-sensitive algebraic mass on split data. The immediate
conservative follow-up is to treat flux-integrated `varpi` as the primary
uncharged Bondi-mass diagnostic and derive the analogous charged scalar
mass/charge balance laws before relying on horizon charge-accumulation
measurements.

`examples/check_charged_horizon_density.jl` is the charged-sector target
from Gelles/Pretorius. For extremal `eQ0=0.6`, the expected late-time
horizon charge-density exponent is `1 - 2s = 0`, i.e. a plateau. The
current scaffold does not pass this check yet; it is kept as a research
diagnostic to drive the next round of charged-sector corrections.

The Baake/Rinne equations cannot be pasted directly because their variables and gauge are CMC hyperboloidal, not compactified double-null. They are still the right source for matter stress tensor, charge conventions, and comparison diagnostics.
