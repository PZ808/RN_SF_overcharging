# Formulation Notes

## References

Baake/Rinne `arXiv:1610.08352v2` evolves the fully coupled Einstein-Maxwell-Klein-Gordon system in spherical symmetry on CMC hyperboloidal slices. The paper gives:

- conformal metric variables `Omega, Ntilde, X, Pi`;
- CMC/isotropic gauge equations;
- Maxwell variables `(ar, Er)` in temporal gauge;
- conformally rescaled scalar `phitilde` and `psitilde`;
- source terms for the Einstein constraints and slicing equation.

Gelles/Pretorius `arXiv:2503.04881v2` evolves charged scalar electrodynamics on a fixed RN background using compactified double-null MRT coordinates. The code starts from this coordinate system because the requested end state is "Baake/Rinne physics in Gelles/Pretorius coordinates."

Gelles/Pretorius `arXiv:2602.11256v1`, submitted on 11 February 2026, now
provides that end state directly: fully nonlinear Einstein-Maxwell-charged
scalar evolution in double-null MRT gauge. The nonlinear module is being
aligned against its Eqs. (3)-(10), while the MRT uncharged tests remain the
baseline regression.

The current equation mapping is:

| Gelles/Pretorius equation | Stored-variable conversion | Code status |
| :--- | :--- | :--- |
| `ds^2=-2 f_GP dU dV+r^2 dOmega^2` | `f_code=2 f_GP` | implemented |
| Eq. (3), `r_UV` | Coulomb coefficient becomes `-f_code(1-Q^2/r^2)/(4r)` | implemented and RN tested |
| Eq. (4), `f_UV` | evolve `log(f_code)`; sources convert `Psi=sqrt(32 pi) r phi_GP` to `Phi=Psi/r` | implemented |
| Eqs. (7)-(8), null metric constraints | `r T_UU/(8 f_code)` and `r T_VV/(8 f_code)` with `T_aa=2|D_a Phi|^2` | implemented |
| Scalar evolution in arXiv:2503.04881 Appendix D | GP2026 path evolves `Psi=sqrt(32 pi) r phi_GP` with the published `bar xi, bar Pi` equations | implemented |
| Eq. (10), renormalized mass | `M=r[1+4r_Ur_V/f_code+Q^2/r^2]/2` | implemented and conservatively checked |
| Maxwell constraints from arXiv:2503.04881 Appendix D | `Q_U=r^2J_U/8`, `Q_V=-r^2J_V/8`, `F_UV=-Qf_code/(2r^2)` | implemented and conservatively checked |

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

The older MRT regression branch retains the normalized scalar:

```text
Phi = sqrt(32 pi) phi_GP
D_a Phi = partial_a Phi - i e A_a Phi
```

The GP2026 production branch instead stores and evolves the paper's reduced
scalar when called with `reduced_scalar=true`:

```text
Psi = r Phi = sqrt(32 pi) r phi_GP
Psi_re,UV = (r_UV/r) Psi_re + e^2 A_U A_V Psi_re
            - e(A_V Psi_im,U + A_U Psi_im,V)
Psi_im,UV = (r_UV/r) Psi_im + e^2 A_U A_V Psi_im
            + e(A_V Psi_re,U + A_U Psi_re,V)
```

Stress-energy and Maxwell current evaluation converts `Psi` back to `Phi`
and its derivatives. Since the paper uses `ds^2=-2 f_GP dU dV`, while this
module stores `f=2 f_GP`, its Faraday component and Maxwell constraints are:

```text
F_UV = -Q f / (2 r^2)
J_a = 2 e Im[Phi^* D_a Phi] = 32 pi J_GP,a
Q_U =  r^2 J_U / 8
Q_V = -r^2 J_V / 8
```

`src/StressEnergy.jl` isolates the reduced scalar sources in this convention:

```text
T_UU = 2 |D_U Phi|^2
T_VV = 2 |D_V Phi|^2
T_uv = alpha^2 / f
T_theta theta = 4 r^2 Re[(D_U Phi)^* D_V Phi] / f + 2 r^2 alpha^2 / f^2
alpha = -Q f / (2 r^2)
F_UV = A_V,U - A_U,V = alpha
```

The Maxwell entries omit the Gaussian-unit `1/(4 pi)` coefficient because the
metric update already carries its Coulomb terms explicitly through `Q`; use
`maxwell_weight=1/(4*pi)` when requesting physical Maxwell stress output. The
helper functions
`outgoing_constraint_source` and `ingoing_constraint_source` provide the
slots used by the two null Raychaudhuri constraints, normalized so the
single-real-scalar uncharged limit matches the MRT `1/4 r f^-1 Phi_a^2`
constraint source.

## Missing Validation Work

The next physics/code steps are:

- converge the charged conservative charge and mass-balance residuals;
- complete charged initial-data construction for the intended MRT experiments;
- validate that the uncharged limit reproduces Murata-Reall-Tanahashi;
- validate that the fixed-metric small-amplitude limit reproduces Gelles/Pretorius;
- validate horizon charge accumulation against the fixed and nonlinear Gelles/Pretorius results.

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

The charged nonlinear diagnostic now follows `arXiv:2602.11256v1`. In the
stored scalar convention, the conservative fixed-`V` laws implemented in
`src/Diagnostics.jl` are:

```text
Q_U = r^2 J_U / 8
varpi_U = -r^2 r_V T_UU / (4 f) + Q Q_U / r
```

The initializer now integrates the `Q_V=-r^2 J_V/8` Maxwell constraint on its
ingoing matter leg before constructing the associated quasi-Lorenz potential.
For direct continuation of the validated MRT experiment,
`initialize_mrt2013_charged_outgoing_wave!` seeds MRT's initially real
outgoing pulse. In quasi-Lorenz gauge, it sets `A_U=0` on the outgoing
initial leg and `A_V=0` on the ingoing initial leg, while integrating the
transverse potential components required by the Coulomb field. Because the
outgoing scalar pulse is initially real and `A_U=0` along its support, its
initial `Q_U` source vanishes, so the uncharged degenerate-horizon `f0`
tuning remains applicable at the initial slice even when `e` is nonzero.
`examples/check_nonlinear_charged_balance.jl` is a short-run regression that
reports both algebraic-versus-flux-integrated `Q` and
algebraic-versus-flux-integrated `varpi`; these residuals should be converged
before using the nonlinear charged code for horizon scaling claims.

This diagnostic exposed an evolution inconsistency immediately: the original
nonlinear cell update treated `A_U,V` and `A_V,U`, which are first-order
Maxwell constraints, as if they were mixed second derivatives in a diamond
update. The current repair imposes those two potential equations with
centered corner derivatives and restores the quasi-Lorenz initial Coulomb
potential. For charge, the production data prescribe `Q(V)` on `U=U0`, so
the interior march advances

```text
Q_11 = Q_01 + Delta U (Q_U)_corner.
```

The complementary constraint `Q_V=-r^2 J_V/8` is now available as
`charged_charge_flux_v_profile` and is not enforced by this update.
For the charged MRT balance regression with `epsilon=0.02`, `e=0.6`, and
`V approximately 4`, the production-oriented update gives:

| initial `U` points | `Delta V` | `max abs(Q_U-source)` | `max abs(varpi_U-source)` | `max abs(Q-Q_flux)` | `max abs(varpi-varpi_flux)` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 531 | 0.020 | 1.76e-7 | 1.89e-7 | 1.41e-7 | 7.93e-8 |
| 1061 | 0.020 | 1.77e-7 | 1.75e-7 | 1.41e-7 | 7.66e-8 |
| 531 | 0.010 | 8.83e-8 | 1.04e-7 | 7.05e-8 | 4.18e-8 |

These residuals remain small compared with the physical fluxes
(`max abs(Q_U) approximately 5.35e-5`,
`max abs(varpi_U) approximately 2.66e-4`) and reduce with `Delta V`.
The previously recorded parts-in-`10^9` charge numbers came from a symmetric
charge corner trial that is inconsistent with an ingoing boundary already
carrying a nonconstant `Q(V)` profile; that trial is not used.

`initialize_gp2026_single_pulse!` now implements the 2026 production family:
ingoing data on `N_B`, `A0=0.01`, `omega_tilde=1`, `(U0,V0)=(-1,0)`,
`M0=1`, and a default near-threshold background charge
`Q0=1.0033218`. The code uses the paper's metric convention through
`f_code=2 f_GP` and stores the boundary pulse directly as
`Psi=sqrt(32*pi) r phi_GP`. For the initial radius gauge, Appendix A/B contain mutually
inconsistent signs; following the main text and the quoted numerical domain,
the implementation uses

```text
r(U,V0) = M0 - U/2,  r_U = -1/2,
r_star(r(U0,V)) = r_star(M0-U0/2) + (V-V0)/2.
```

This yields `r(1.6,0)=0.2`, as stated in the paper. A boundary diagnostic
then found and fixed a separate sign error in our transcription of
`d[A(V) cos(omega V)]/dV`. At the paper parameters and `Delta V=0.08`,
the corrected initializer gives `Delta Q=0.01477994` and
`max abs(Q_V-source)=1.64e-6`; the latter decreases by approximately four
under each halving of `Delta V`. The short production evolution in
`examples/check_gp2026_short_balance.jl` gives:

| `Delta V` | `max abs(Q_U-source)` | `max abs(Q_V-source)` | `max abs(Q-Q_flux)` |
| ---: | ---: | ---: | ---: |
| 0.080 | 8.86e-7 | 1.68e-6 | 1.49e-7 |
| 0.040 | 4.41e-7 | 4.47e-7 | 7.51e-8 |

`examples/convergence_gp2026_charge_residuals.jl` now performs a cleaner
two-dimensional convergence study by halving both `Delta U` and `Delta V` on
the fixed short domain `U in [-1,-0.8]`, `V in [0,20]`. With the GP2026
production hyperbolic charge update, sampled at `U=-0.9` and `V approximately
10`, the charge residuals converge at second order:

| `Delta U` | `Delta V` | `max abs(Q_U-source)` | rate | `max abs(Q_V-source)` | rate | `max abs(Q-Q_flux)` | rate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01000 | 0.080 | 1.089e-6 | - | 1.651e-6 | - | 1.830e-7 | - |
| 0.00500 | 0.040 | 2.757e-7 | 1.98 | 4.128e-7 | 2.00 | 4.723e-8 | 1.95 |
| 0.00250 | 0.020 | 6.936e-8 | 1.99 | 1.032e-7 | 2.00 | 1.199e-8 | 1.98 |
| 0.00125 | 0.010 | 1.740e-8 | 2.00 | 2.580e-8 | 2.00 | 3.021e-9 | 1.99 |

The row-marched production path shows the same behavior on the same short
domain when the paper-AMR constant `C` is halved with `Delta V`:

| `Delta V` | `C` | U rows | `max abs(Q_U-source)` | rate | `max abs(Q_V-source)` | rate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.080 | 0.100 | 8 | 2.420e-6 | - | 1.655e-6 | - |
| 0.040 | 0.050 | 15 | 6.210e-7 | 1.96 | 4.138e-7 | 2.00 |
| 0.020 | 0.025 | 29 | 1.572e-7 | 1.98 | 1.035e-7 | 2.00 |

These tests check the Maxwell sector in a region that reaches the same
physical endpoint at every resolution. They do not by themselves validate the
near-horizon critical scaling, where the Eq. (9) grid accumulation changes the
meaning of a fixed `Umax` comparison.

Re-running these checks after converting the GP2026 evolution path from
`Phi` to the published reduced field `Psi=r*Phi` leaves the quoted short
balance residuals and the long-run transition below unchanged at the reported
precision. The representation discrepancy is therefore repaired but is not
the cause of the threshold disagreement.

The first longer run is not yet a reproduction of the paper's critical
solution. `examples/check_gp2026_long_evolution.jl` reports horizon
appearance and loss of finiteness rather than silently returning late-time
diagnostics. Using the implemented internally consistent initial-radius gauge,
`Delta U=0.01`, and `Delta V=0.08` gives:

| `Q0` | requested `Vmax` | first trapped `V` | first nonfinite `V` | outcome |
| ---: | ---: | ---: | ---: | :--- |
| 1.0033218 | 100 | 6.56 | 50.64 | apparent horizon forms during the pulse; evolution fails later |
| 1.0200000 | 400 | none | none | finite, dispersive control run |

The reduced-field diagnostic `gp2026_horizon_rphi_series` now directly
samples the paper's observable `|r phi_GP|` at detected apparent horizons.
For the `Q0=1.0033218` run above it decreases from `1.0485e-2` at `V=6.56`
to `1.3148e-4` at `V=50.56` over 551 finite horizon samples.

For the finite `Q0=1.02` control at `V approximately 400`,
`r_V` remains positive in the range `0.49366` to `0.49468` and
`max abs(Q_U-source)=2.92e-11`. The alleged threshold value from the paper
therefore does not agree with this implementation yet. This discrepancy is
already present during the transient and is not explained by late-time
horizon resolution.

The next controlled implementation step follows the row-oriented production
evolution direction in arXiv:2602.11256. `evolve_gp2026_u_adaptive` stores
rows at fixed `U` and fills each row across the prescribed `V` mesh. The
literal paper step rule is

```text
Delta U = C/f_GP(U,Vmax) = 2C/f_code(U,Vmax).
```

The solver keeps that rule available as `step_control=:outer`, but the default
uses `step_control=:local`, the smaller of `Delta U=2C/max_row(f_code)` and a
local areal-radius limiter

```text
|r_U| Delta U <= C |Delta r|.
```

This follows the stability logic of Gundlach/Baumgarte/Hilditch
arXiv:1908.05971, whose double-null method notes that there is no causal
Courant condition in these coordinates but still imposes a geometric step
bound. It also provides the Appendix-D hyperbolic evolution equation for `Q`,
converted to `Psi=sqrt(32*pi) r phi_GP` and `f_code=2 f_GP`. The prior
`Q_U` constraint march remains available as a comparison diagnostic, not as
the production default. `examples/check_gp2026_u_refinement.jl` runs either
branch and can promote arithmetic to `BigFloat`.

At `Q0=1.0033218`, `A0=0.01`, `e Q0=0.6`, and `Vmax=100`, the hyperbolic
charge branch with the literal `outer` step control gives:

| `Delta V` | `C` | first trapped `V` | `max abs(Q_U-source)` | `max abs(Q_V-source)` |
| ---: | ---: | ---: | ---: | ---: |
| 0.080 | 0.600 | 51.36 | 8.64e-5 | 4.74e-6 |
| 0.040 | 0.300 | 45.88 | 1.13e-5 | 1.08e-6 |
| 0.020 | 0.150 | 24.80 | 6.76e-6 | 2.16e-7 |

Thus the charge residuals decrease under this joint refinement, but the
apparent-horizon time has not converged and every listed run eventually
encounters a nonfinite row. The paper-style marching implementation is useful
for the next debugging comparison; it is not yet a validated reproduction of
the nonlinear threshold solution.

`examples/diagnose_gp2026_nonfinite.jl` localizes the first coarse
`Delta V=0.08`, `C=0.6` nonfinite row. The immediate failure is not a missing
Maxwell source or a horizon regularity singularity: after the apparent-horizon
region forms, the largest `f_code` value moves away from `Vmax`, while
`f_code(U,Vmax)` collapses from about `1.7e4` to `12`. The paper step rule
then jumps from `Delta U=6.9e-5` to `Delta U=0.10`; at `V approximately 50.4`
the fixed-point cell update drives `r_11` negative and then overflows on the
next cell. A diagnostic `step_control=:max_row` option instead uses
`Delta U=2C/max_row(f_code)` and reaches `U=1.6` without a nonfinite row for
that same run, confirming that the blow-up is caused by row step control.
The stricter `local` controller also avoids the invalid row and gives
`max abs(Q_U-source)=8.01e-6` at the 1000-row cap for the same coarse run,
but it does not yet reach `U=1.6`; it stops cleanly near `U=-0.202`. This is
the current default fix. It is not yet a physics result by itself: the
apparent-horizon location and charge residuals still require a convergence
study with local refinement and point removal/insertion near the stiff layer.

This failure mode is sensitive to the outer boundary. The paper's production
domain uses `Vmax=400`, not the shorter `Vmax=100` stress test above. With
`Vmax=400`, `Delta V=0.08`, `C=0.6`, and the literal Eq. (9) rule, the outer
boundary metric grows monotonically in the expected way and the coordinate
spacing collapses instead of producing an invalid row over the first 120 rows:

| row | `U` | `Delta U` | `f_code(U,Vmax)` |
| ---: | ---: | ---: | ---: |
| 2 | -0.869658 | 1.30e-1 | 9.21e0 |
| 41 | -0.206532 | 2.14e-4 | 5.60e3 |
| 81 | -0.205396 | 1.49e-7 | 8.07e6 |
| 111 | -0.205395 | 3.59e-10 | 3.34e9 |

The correct paper-reproduction interpretation is therefore not "does the
coordinate reach `U=1.6`?" in a horizon-forming run. Eq. (9) is intended to
accumulate rows near the event-horizon/throat location. The diagnostic
`examples/check_gp2026_paper_amr.jl` records this limiting-`U` behavior
directly and should be used alongside the shorter stress tests.

For the critical-phenomena comparison in Sec. IIIA, the direct observables are
not the limiting `U` itself but the trapped-surface formation time and late
horizon geometry:

```text
V_trap proportional to (Q_* - Q0)^(-1/2),
1 - Q/M proportional to (Q_* - Q0),
1 - r/M proportional to (Q_* - Q0)^(1/2).
```

`examples/scan_gp2026_threshold.jl` now scans the paper's one-parameter
family and prints both the direct Sec. IIIA columns (`slice_trap_*`,
`final_AH_*`) and the limiting-surface proxy columns (`last_U`, `next_du`,
`outer_f`, `min_rv`, `max_rho`). A first 120-row `Vmax=400` scan gives a
smooth proxy trend on the BH side,

| `Q0` | `last_U` | `min r_V` on final row | `max rho` | direct `V_trap` |
| ---: | ---: | ---: | ---: | :--- |
| 1.0000000 | -0.250919 | 8.87e-4 | 2.08 | missing |
| 1.0020000 | -0.224237 | 9.04e-4 | 2.21 | missing |
| 1.0030000 | -0.210081 | 9.13e-4 | 2.29 | missing |
| 1.0033218 | -0.205395 | 9.16e-4 | 2.31 | missing |

This is suggestive of the expected limiting surface as `Q0` approaches the
quoted `Q_*`, but it is not yet a measurement of the Sec. IIIA power laws.
The current horizon/trapped-surface extraction does not see trapped cells in
these paper-AMR rows before Eq. (9) has accumulated at the marginal surface.
The next debugging target is therefore the trapped-region detector and/or the
ability to continue through the accumulated event-horizon layer, not another
claim of critical scaling.

To make the near-horizon throat explicit, the row diagnostics now include

```text
y = r - |Q|,        rho = -log(y/|Q|).
```

`rho` is a stretched extremal-throat coordinate: large `rho` means the row is
close to the would-be AdS2/JT region. The optional `step_control=:throat`
limits row-to-row changes in this coordinate, while the default `:local`
controller takes the minimum of the largest-`f`, geometric-`r`, and
throat-`rho` limits. The same diagnostic reports a matching candidate and
the full band where `rho >= rho_match`; for the coarse `Vmax=100`,
`Delta V=0.08`, `C=0.6` run at the 1000-row cap, `rho_match=2` gives a band
from `V=0` to `V=4.88`, with `max rho=2.33`. This is the data needed for a
future matched full-system plus near-AdS2/JT patch: choose a `rho=rho_match`
interface, pass `r`, `Q`, `Psi`, and fluxes across it, and evolve the deeper
throat with effective near-horizon variables.

A first check of the lapse in `rho` coordinates compares `log f` with

```text
log f_(U rho) = log f_(U V) - log |rho_V|.
```

This is only a coordinate diagnostic on the existing `V` grid. For the
default local-controller run at the 1000-row cap, the throat band
`rho >= 2` already has a small raw `log f` range, about `0.42`; the
`U-rho` coefficient has range about `0.59`. For the failing literal
outer-step run, the last valid row has a throat-band raw `log f` range of
about `8.54`, reduced only modestly to `7.62` by the `rho` transform, while
the full-row transformed range is worse. The important signal is instead
`max |Delta rho| approximately 11.5` on that last valid outer-step row:
the violent behavior is an unresolved throat-coordinate gradient. A simple
post-processing coordinate change does not tame it; a real `rho`-adapted mesh
or matched near-horizon patch is needed.

`examples/check_charged_horizon_density.jl` is the charged-sector target
from Gelles/Pretorius. For extremal `eQ0=0.6`, the expected late-time
horizon charge-density exponent is `1 - 2s = 0`, i.e. a plateau. The
current scaffold does not pass this check yet; it is kept as a research
diagnostic to drive the next round of charged-sector corrections.

The Baake/Rinne equations cannot be pasted directly because their variables and gauge are CMC hyperboloidal, not compactified double-null. They are still the right source for matter stress tensor, charge conventions, and comparison diagnostics.
