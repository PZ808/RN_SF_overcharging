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

### GP2026 Equation Audit

The June 2026 audit of arXiv:2602.11256 found one displayed-equation
inconsistency that should not be copied into the code. With the paper's metric
`ds^2=-2 f_GP dU dV+r^2 dOmega^2`, Eq. (4) is printed with
`+f_GP^2/(2r^2)-2 f_GP^2 Q^2/r^4` in the `f_GP,UV` equation. Converting this
literal equation to the stored variable `f_code=2 f_GP` would give

```text
(log f_code)_UV = f_code/(4r^2) + 2 r_U r_V/r^2
                  - f_code Q^2/r^4 - scalar source.
```

That formula fails the exact RN electrovac residual test. The implemented
formula,

```text
(log f_code)_UV = f_code/(2r^2) + 2 r_U r_V/r^2
                  - f_code Q^2/r^4 - scalar source,
```

is the one that converges on exact RN. A direct residual comparison on extreme
RN gave, at resolutions `(60,180)`, `(120,360)`, `(240,720)`,

| formula | residuals |
| :--- | :--- |
| implemented `+f/(2r^2)-fQ^2/r^4` | `2.99e-1`, `1.55e-1`, `5.07e-2` |
| literal GP Eq. (4) after `f_code=2f_GP` | `1.36`, `1.93`, `2.21` |
| extra Coulomb factor `-2fQ^2/r^4` in `f_code` | `2.88`, `4.33`, `5.01` |

Thus Eq. (4)'s displayed `f^2/(2r^2)` coefficient is best treated as a typo or
normalization slip. The mass identity in Eq. (10) is internally consistent with
the paper metric and remains our reference for diagnostics.

Appendix C also appears to drop a minus sign in its displayed `r_U` formula.
Differentiating the stated MRT-gauge relation
`r_*(r_+ + V/2)+r_*(r_+ - U/2)=r_*(r)` gives a negative `r_U`; the main text
also states that the GP single-pulse gauge has `r_U=-1/2` on `N_A`. The code
uses the negative sign.

Appendix A is a convention trap rather than a confirmed typo. For
super-extremal data, GP switch to "extremal gauge": the initial-leg radius is
the extremal reference radius, while the corner lapse is normalized with the
super-extremal `Q0` through Eq. (A4). Their displayed Eq. (A3), taken
literally with `(U0,V0)=(-1,0)`, is singular at the stated corner. A shift of
the extremal tortoise-coordinate origin, or the equivalent subtraction of
its corner value, is therefore implicit. With that corner normalization and
the main-text conditions `r_U=-1/2` and `r(1.6,0)=0.2`, the pulse-leg gauge
reduces to the production `:areal_affine` data:

```text
r(U,V0) = 1 - U/2,
r(U0,V) = 3/2 + V/2.
```

The initializer sets `f_code=2 f_GP` and obtains `f_code(U0,V0)` from the
Hawking-mass normalization. The statement that pure ingoing data allow the
metric on `N_A` to be exact RN should therefore be read as a gauge-fixed
characteristic construction, not as permission to replace the `N_A` data
with subextremal MRT-gauge RN formulae.

`gp2026_initial_constraint_residuals` and
`examples/check_gp2026_initial_data.jl` directly audit the data written by
the production initializer, using the same midpoint and trapezoidal
quadratures. At `Q0=1.0033218` and `Delta V=0.08`, the corner mass is one
exactly, the lapse residual is below `4e-16`, the charge-increment residual
is below `1.2e-16`, and the potential residuals are below `2e-15`.

`examples/diagnose_gp2026_initial_data.jl` additionally separates the two
possible readings of `r_V` on `N_A`. With
`Q0=1.0033218`, `A0=0.01`, and `eQ0=0.6`, using the mass-compatible RN value
`r_V=f F_Q(r)/2` gives `max |M-M0|=2.22e-16`. Using the tempting but wrong
extremal-reference value `r_V=F_ext(r)/2` away from the corner gives
`max |M-M0|=2.60e-2`. The same script checks the pulse leg `N_B`; at
`Delta V=0.16, 0.08, 0.04`, the residuals
`max |Q_V+r^2J_V/8|` are `4.94e-7`, `1.23e-7`, `3.09e-8`, and the log-lapse
constraint residuals are `7.81e-6`, `1.95e-6`, `5.81e-7`. Thus the initial
data are internally consistent to the expected finite-difference order. The
remaining threshold shift is unlikely to come from a simple constraint
integration error, though the Appendix-A coordinate-origin ambiguity remains.

There is a separate ambiguity in the radius gauge on the pulse-carrying
initial leg `N_B`. The displayed MRT relation in GP Appendix A cannot pass
through the stated corner `(U0,V0)=(-1,0)` while also satisfying the main-text
conditions `r_U=-1/2` and `r(1.6,0)=0.2`. The production initializer currently
uses the corner-compatible areal-affine reading

```text
r(U0,V) = r0 + V/2.
```

The earlier EF-affine reading remains available with
`pulse_leg_gauge=:ef_affine`:

```text
r_star(r(U0,V)) = r_star(r0) + V/2.
```

`examples/compare_gp2026_initial_gauges.jl` compares that choice with the
EF-affine reading. At `Q0=1`, the post-pulse characteristic data have
`M-Q=5.3509e-3` (areal affine) and `1.5968e-2` (EF affine). At the quoted
`Q0=1.0033218`, the corresponding gaps are `2.0831e-3` and `1.2779e-2`.
The areal-affine result is much closer to the scale of the reported threshold
offset `Qstar-1=3.3218e-3`, which is why it is now the production default.

### Direct-lapse Newton update

GP solve all fields at the unknown cell corner simultaneously with
Newton-Raphson and discretize `f` directly. The production row path now
matches that design with `cell_solver=:newton_direct`; the previous fixed-count
Picard update of `log(f)` remains available as `:picard_log`. Newton uses a
scaled finite-difference Jacobian, positivity-preserving damping for `r` and
`f`, and defaults to `(rtol,atol)=(1e-13,1e-15)`.

On `U in [-1,-0.8]`, `V in [0,20]`, the maximum unscaled residual over all
seven discrete cell equations is `1.78e-14`, `5.31e-15`, and `1.15e-12` at
joint resolutions `(Delta U,Delta V)=(0.01,0.08)`, `(0.005,0.04)`, and
`(0.0025,0.02)`. The independent charge constraints converge at second order:
`Q_U` rates are `1.98,1.99`, and `Q_V` rates are `2.00,2.00`.

Before root-step rejection was added, the corrected horizon controls did not
produce a completed negative `r_V` row. With the literal Eq. (9) rule:

| `Q0` | limiting `U` | closest positive `r_V` | status |
| ---: | ---: | ---: | :--- |
| 1.0 | -0.156800 | 5.51e-3 | Float64 coordinate stall |
| 1.0033218 | -0.082402 | 4.76e-3 | Float64 coordinate stall |

The local controller gets much closer for `Q0=1`: after 5,500 rows it reaches
`U=-0.12313525` and `min r_V=2.25e-5` at `V=144.24`, with a smallest realized
`Delta U=1.89e-14`. The Newton cells remain finite, while the derivative-based
`Q_U` diagnostic loses precision; the transverse `Q_V` residual remains
`1.73e-6`.

A separate uniform-`U` experiment forces the grid through the Eq. (9)
limiting coordinate. It confirms that this surface is not a physical endpoint,
but Newton fails in the interior before a complete trapped row is obtained.
At `Delta U=0.002`, the last complete row is `U=-0.114` with
`min r_V=2.22e-3`; at `Delta U=0.0005`, it is `U=-0.111` with
`min r_V=5.39e-4`. The near factor-four reduction tracks the step refinement.
This identifies a resolution-dependent failure in the pulse evolution, but
does not by itself show that the MRT variables or bulk update cannot cross a
horizon.

### Analytic electrovacuum horizon crossing

The exact extremal-RN solution provides a stronger separation test. In the
GP2026 normalization of the uncompactified MRT chart,

```text
rstar(r) = rstar(M-U/2) + (V-V0)/2,
f_code = F(r)/F(M-U/2),
ds^2 = -f_code dU dV.
```

The future event horizon is the regular grid line `U=0`, where `r=M` and
`f_code=1`. `initialize_gp2026_exact_extremal_rn!` seeds this exact solution
on both initial characteristic legs, and
`examples/convergence_gp2026_electrovacuum.jl` evolves from `U=-0.4` through
the horizon to `U=0.2`.

| `Delta U` | `Delta V` | max `|delta r|` | rate | max `|delta log(f)|` | rate | horizon `|delta r|` | horizon `|delta f|` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.020 | 0.100 | 3.445e-5 | - | 1.334e-3 | - | 3.413e-5 | 4.061e-4 |
| 0.010 | 0.050 | 8.613e-6 | 2.000 | 3.331e-4 | 2.001 | 8.532e-6 | 1.015e-4 |
| 0.005 | 0.025 | 2.153e-6 | 2.000 | 8.326e-5 | 2.000 | 2.133e-6 | 2.538e-5 |

The direct cell residual is at roundoff and `Q` remains exactly constant.
Therefore the direct-Newton bulk equations do penetrate the horizon with
second-order convergence. The missing trapped surface in the scalar runs must
instead involve the pulse-leg initial data, adaptive approach to the
near-marginal layer, or their interaction; it is not a generic coordinate or
electrovacuum evolution obstruction.

The zero-pulse initial-leg comparison also distinguishes the two available
gauge readings. On `U0=-1`, `V in [0,10]`, `:ef_affine` agrees with the exact
analytic chart to about `5e-13`, while `:areal_affine` differs by `3.82` in
`r` and `2.20` in `log(f)`. The areal-affine data may describe RN in a
different null reparameterization, but they are not the analytic
GP-normalized MRT chart above. Reproduction runs should therefore explicitly
compare the two choices; the production default is not changed by this test.

A controlled charged comparison at `Q0=1`, `eQ0=0.6`, `A0=0.01`,
`Vmax=400`, `Delta V=0.08`, and a 6,000-row local-controller budget changes
only this pulse-leg gauge. Neither branch produces a trapped row:

| pulse-leg gauge | rows | last `U` | status | closest `V` | closest positive `r_V` | `r-Q` there |
| :--- | ---: | ---: | :--- | ---: | ---: | ---: |
| `:areal_affine` | 6000 | -0.123135 | max rows | 164.64 | 2.03e-5 | 0.1397 |
| `:ef_affine` | 5295 | -0.249278 | precision stalled | 172.72 | 2.11e-5 | 0.2180 |

Thus selecting the analytically matched electrovacuum gauge does not recover
the GP trapped surface and actually reaches the Float64 limiting layer sooner
in this test. This is not caused by a low-order charged update:
`examples/convergence_gp2026_charge_residuals.jl` with `:ef_affine` gives
rates `1.98,1.99,2.00` for the `Q_U` residual and `2.00,2.00,2.00` for
`Q_V`; the integrated charge mismatch approaches second order as well. The
infinite late-run `Q_U` diagnostic is derivative roundoff once adjacent
`U` rows differ at machine precision, while the transverse `Q_V` residual
remains finite at `1.70e-6`.

`cell_equation_residual_summary` and
`examples/check_gp2026_cell_residuals.jl` audit the evolved rectangular-cell
equations directly. On the GP2026 production branch
(`Q0=1.0033218`, `A0=0.01`, `eQ0=0.6`, `U in [-1,-0.8]`, `V in [0,20]`,
hyperbolic charge update), the implemented centered equations close at the
solver/roundoff level: at `Delta U=0.01`, `Delta V=0.08`,
`max |r_UV-RHS|=1.66e-12`, `max |logf_UV-RHS|=9.27e-13`,
`max |Psi_UV-RHS|=2.83e-14`, `max |Q_UV-RHS|=4.16e-13`, and the Lorenz/Faraday
potential residuals are below `3e-15`. The unevolved charge constraints
remain finite-difference diagnostics and converge at second order:
`max |Q_U-source|` goes `2.52e-6 -> 6.34e-7 -> 1.59e-7 -> 3.98e-8`, while
`max |Q_V-source|` goes `1.68e-6 -> 4.21e-7 -> 1.05e-7 -> 2.63e-8`.
The negative-control metric variants remain large: the literal printed GP
Eq. (4) conversion gives `O(1.9e-1)`, and an extra Coulomb factor gives
`O(2.6e-1)`. This pushes the missing trapped-surface crossing away from a
local cell-equation sign/factor error and toward global/gauge/discretization
issues in the stiff throat evolution.

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
uses `step_control=:local`, the minimum of four caps:

- `Delta U=2C/max_row(f_code)`;
- a local areal-radius limiter

```text
|r_U| Delta U <= C |Delta r|.
```

- throat-coordinate limiters on row-to-row changes in
  `rho=-log((r-|Q|)/|Q|)` and `eta=1-|Q|/r`.

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
family and prints both the direct Sec. IIIA columns (`direct_vtrap_*`,
`slice_trap_*`, `final_AH_*`) and the limiting-surface proxy columns
(`last_U`, `next_du`, `outer_f`, `min_rv`, `max_rho`,
`vtrap_proxy_*`). The helper

```julia
vtrap_diagnostic(rows; missing_status)
```

defines `direct_vtrap_*` as the earliest row-local apparent-horizon crossing
in `V`, where the centered row diagnostic `r_V` reaches zero. If no crossing
is present, `vtrap_status` records why the run stopped and `vtrap_proxy_*`
reports the row point with the smallest still-positive `r_V`.

A first 120-row `Vmax=400` scan gives a smooth proxy trend on the BH side,

| `Q0` | `last_U` | `vtrap_status` | proxy `V` | proxy `min r_V` | `max rho` | direct `V_trap` |
| ---: | ---: | :--- | ---: | ---: | ---: | :--- |
| 1.0000000 | -0.250919 | `max_rows` | 104.08 | 8.87e-4 | 2.08 | missing |
| 1.0020000 | -0.224237 | `max_rows` | 103.68 | 9.04e-4 | 2.21 | missing |
| 1.0030000 | -0.210081 | `max_rows` | 103.28 | 9.13e-4 | 2.29 | missing |
| 1.0033218 | -0.205395 | `max_rows` | 103.12 | 9.16e-4 | 2.31 | missing |

This is suggestive of the expected limiting surface as `Q0` approaches the
quoted `Q_*`, but it is not yet a measurement of the Sec. IIIA power laws.
The direct trapped-surface extraction does not see trapped cells in these
paper-AMR rows before Eq. (9) has accumulated at the marginal surface. This is
not just a Float64 detector miss: a 160-bit `Q0=1.0`, 220-row check still has
positive proxy `min r_V=7.63e-4` at `V=240.64`. The next debugging target is
therefore why the row evolution approaches a positive-expansion limiting
surface instead of producing the Sec. IIIA trapped-surface curve, not another
claim of critical scaling.

`examples/check_gp2026_stiffness.jl` tests whether this is a local stepping
problem. It evolves to a chosen base row, then advances to the same target `U`
with one row step and with repeated substeps. At `Q0=1.0033218`, `Vmax=400`,
`Delta V=0.08`, `C=0.6`, and 40 paper-AMR rows, the outer-boundary rule selects
`Delta U=2.14e-4`, while the local geometric limiter would select
`8.45e-6`. Using the full outer step, the one-step row is finite and Picard
iterations converge, but it differs from a 16-substep reference by
`max |Delta r|=9.21e-4` and `max |Delta log f|=4.27e-4`. Repeating the same
base-row probe with `Delta U` reduced by a factor `0.04` gives
`max |Delta r|=3.95e-8` and `max |Delta log f|=1.94e-8`. This points to a
local stiffness/step-size problem in the row march, not primarily to a lack of
fixed-point iterations inside the cell solve.

`examples/analyze_gp2026_trapped_surfaces.jl` adds the corresponding
trapped-surface analysis. It follows the row point where `r_V` is smallest and
prints both null expansions, the renormalized Hawking mass, and the horizon
function

```text
H = 1 - 2M/r + Q^2/r^2 = -4 r_U r_V/f.
```

For the 120-row paper-AMR run above, the closest point is on the last row at
`V=103.12`: `r_V=9.16e-4`, `r_U=-3.09e5`, `H=1.83e-3`, and
`r-r_+=6.27e-3`. The identity for `H` closes to roundoff, so the detector and
mass diagnostic agree that this point is still outside the apparent horizon.
A 1000-row local-controller run reaches `U=-0.202078` and gives its closest
point at `V=18.32`, with `r_V=9.93e-4`, `H=1.95e-3`, and `r-r_+=8.45e-3`.
Thus the missing `Vtrap` is not currently a sign convention error in the
trapped-surface detector: both expansion and Hawking-mass diagnostics remain
positive before the run stops or exhausts the row budget.

The row marcher now supports adaptive substepping inside a selected macro row
step. For example, `step_control=:outer, substep_control=:local` keeps the
paper-style Eq. (9) macro target but advances toward it with smaller stored
rows satisfying the local geometric, `rho`, and `eta` throat caps. This
substantially changes the near-marginal behavior but does not yet reproduce
the GP trapped-surface curve.
For the quoted-critical run with `max_rows=6000`, the closest point is
`V=144.40`, `r_V=6.13e-6`, `H=1.22e-5`, and `r-r_+=5.33e-5`; for a BH-side
`Q0=1.0` run it reaches `V=139.04`, `r_V=5.62e-6`, `H=1.12e-5`, and
`r-r_+=4.53e-5`. A middle BH-side sample at `Q0=1.002` gives the same pattern:
`V=144.40`, `r_V=5.90e-6`, `H=1.18e-5`, and `r-r_+=4.97e-5`. These runs
remain future-untrapped (`r_U<0`, `r_V>0`) and
eventually stall in the limiting layer. Substepping therefore confirms the
outer row step was too crude, but it also exposes a remaining physics or
coordinate/equation discrepancy: the solution approaches a marginal surface
from outside instead of crossing to the Sec. IIIA apparent horizon.

Tail fits in the trapped-surface analyzer make this more quantitative. On the
last 300 proxy samples, a floor model `y(V)=y_inf+a/V` gives positive limits
and no finite crossing. For the quoted-critical run, the fitted floors are
`r_V -> 6.02e-6`, `H -> 1.20e-5`, and `r-r_+ -> 5.23e-5`; for the BH-side
`Q0=1.0` run they are `r_V -> 5.57e-6`, `H -> 1.11e-5`, and
`r-r_+ -> 4.49e-5`; and for `Q0=1.002` they are `r_V -> 5.80e-6`,
`H -> 1.16e-5`, and `r-r_+ -> 4.88e-5`. A zero-asymptote power-law fit returns exponents near
zero, another sign of plateauing. This does not support a simple "GP
extrapolated a finite trapped surface from our kind of tail" explanation:
with the current equations/data, both the direct signs and the fitted tails
remain outside the future apparent horizon. Negative values of `H` that appear
near the small-`V` initial corner are classified separately by expansion signs
and are not future trapped (`r_U` and `r_V` are both positive there).

To make the near-horizon throat explicit, the row diagnostics now include

```text
y = r - |Q|,
rho = -log(y/|Q|),
eta = 1 - |Q|/r,
zeta = |Q|/y.
```

`rho` is a stretched extremal-throat coordinate: large `rho` means the row is
close to the would-be AdS2/JT region. `eta` is a bounded rational throat
distance, while `zeta` is the direct reciprocal distance. The optional
`step_control=:throat` limits row-to-row changes in `rho`, and
`step_control=:eta` limits row-to-row changes in the compact rational
coordinate `eta`, capped by the largest-`f` row rule so it remains a
refinement criterion. The default `:local` controller takes the minimum of the
largest-`f`, geometric-`r`, throat-`rho`, and throat-`eta` limits. The same
diagnostic reports a matching candidate and the full band where
`rho >= rho_match`; for the coarse `Vmax=100`,
`Delta V=0.08`, `C=0.6` run at the 1000-row cap, `rho_match=2` gives a band
from `V=0` to `V=4.88`, with `max rho=2.33`. This is the data needed for a
future matched full-system plus near-AdS2/JT patch: choose a `rho=rho_match`
interface, pass `r`, `Q`, `Psi`, and fluxes across it, and evolve the deeper
throat with effective near-horizon variables.

The first fixed-`rho` boundary-layer extractor is now implemented. The
diagnostic functions

```julia
throat_boundary_sample(row; rho_match=2.0)
throat_boundary_series(rows; rho_match=2.0)
```

interpolate the largest-`V` crossing of a requested `rho_match` on each
constant-`U` GP2026 row. `examples/extract_gp2026_throat_boundary.jl` evolves
the paper-AMR row system and prints a TSV time series of the matching surface:
`U`, `V(rho_match)`, `r`, `Q`, `r-|Q|`, `log f`, the stored scalar
`Psi=sqrt(32*pi) r phi_GP`, the GP reduced amplitude `r phi_GP`, gauge
potentials, `V` derivatives, and the local `J_V`, `T_VV`, and `Q_V` residual.
The helper

```julia
throat_boundary_observables(sample, ep)
```

collects the gauge-invariant throat quantities used by that extractor:
`|Psi|`, `|r phi_GP|`, the covariant phase gradient
`theta_V - e A_V = Im(Psi^* D_V Psi)/|Psi|^2`, `|D_V(r phi_GP)|`,
`J_V`, `T_VV`, the constraint source terms, `Q/r`, `1-|Q|/r`, and the
coordinate-diagnostic
`log f_(U rho) = log f_(U V) - log |rho_V|`. The raw scalar phase is still
printed for mode tracking, but it is gauge dependent. Conversely,
`theta_V-e A_V` is gauge invariant but divides by `|Psi|^2`, so mode fits
should apply an amplitude threshold before using it in the wave tail.
For example, at `Q0=1.0033218`, `Vmax=400`, `Delta V=0.08`, `C=0.6`,
`rho_match=2`, and a 40-row cap, the extractor finds 25 fixed-`rho` samples,
with the outer matching surface moving from `V approximately 0.20` to
`V approximately 4.40` over rows 16 through 36. This is intended as the first
interface data stream for comparing the full nonlinear dynamics to a near-AdS2
or resonant-QNM reduced model.

`examples/convergence_gp2026_throat_boundary.jl` checks this extractor under
joint paper-AMR refinement. At `rho_match=2`, `Vmax=400`, and
`(Delta V,C,max_rows)=(0.08,0.6,80),(0.04,0.3,160),(0.02,0.15,320)`, the
fixed-`rho` local `Q_V` residual decreases as

| `Delta V` | `C` | samples | `max abs(Q_V residual)` |
| ---: | ---: | ---: | ---: |
| 0.080 | 0.600 | 65 | 1.884e-6 |
| 0.040 | 0.300 | 129 | 5.531e-7 |
| 0.020 | 0.150 | 258 | 1.434e-7 |

The residual rates are approximately `1.77` and `1.95`. Comparing the
extracted boundary fields as functions of `U`, with the finer run
interpolated onto the coarser samples, gives rates about `1.86` for
`V(rho_match)`, `1.55` for `r` and `Q`, and `2.01` for `|r phi_GP|` between
the last two refinement pairs. The derivative columns in the extractor now use
centered row derivatives interpolated to the fixed-`rho` crossing; raw
single-cell secants made the local residual look only first-order.

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

`examples/check_gp2026_throat_stiffness_scaling.jl` is the second-generation
stiffness diagnostic. It evolves the row system, samples the minimum-`r_V`
point, and reports `rho`, `H`, `logf_rho`, row-to-row changes, candidate step
limits, and local one-step versus two-half-step probes. In the default
`Q0=1.0033218`, `Vmax=400`, `Delta V=0.08`, `max_rows=1200` run with
`step_control=:outer` and `substep_control=:local`, the local stiffness grows
rapidly as the throat layer is approached. The max-norm row derivative of the
throat coordinate increases from `rho_U_inf=1.22` near row 2 to `1.04e3` by
row 1200, and `logf_rho_U_inf` increases from `1.48` to `9.88e2`. The GP
outer step stays about `25` times larger than the geometric local cap in the
tail: at row 1200,

```text
Delta U_outer = 1.37e-4,
Delta U_geometric = 5.41e-6,
Delta U_throat = 2.41e-4.
```

A local probe that advances the final row by one GP outer step and compares it
to two half steps remains finite, but it produces
`max |Delta r|=5.36e-4`, `max |Delta logf|=3.04e-4`,
`max |Delta rho|=6.68e-5`, `max |Delta logf_rho|=3.00e-4`, and
`max |Delta H|=4.82e-4`. Thus the stiff layer is not primarily a Picard
nonconvergence: it is a variable/mesh conditioning problem where invariant or
near-throat variables change more gently than raw metric fields, but the GP
outer step is still far outside the local geometric resolution scale.

The same diagnostic now compares `rho` with two rational throat variables:

```text
eta = 1 - |Q|/r = (r-|Q|)/r,
zeta = |Q|/(r-|Q|).
```

Here `eta` is the useful compact rational coordinate. In the default run at
row 1200 it gives `eta_U_inf=2.59e2`, compared with `rho_U_inf=1.04e3`,
while the direct reciprocal gives `zeta_U_inf=1.73e3`. The one-step versus
two-half-step probe shows the same ordering:
`max |Delta eta|=1.65e-5`, `max |Delta rho|=6.68e-5`, and
`max |Delta zeta|=6.50e-5`. The direct inverse is therefore not a good
marching coordinate by itself; the bounded rational distance `eta`, or a
closely related compactified near-horizon coordinate, is the more promising
candidate for a throat-adapted mesh.

The first implementation of that idea is an `eta` row-step limiter. It
estimates `max |partial_U eta|` from the two most recent rows and sets
`Delta U_eta=max_delta_eta/max |partial_U eta|`, with
`max_delta_eta=0.025` by default. `step_control=:eta` uses this as a refinement
cap on top of the largest-`f` row rule, selecting
`min(Delta U_maxrow, Delta U_eta)`. The conservative `step_control=:local`
uses the minimum of the largest-`f`, geometric, `rho`, and `eta` caps. A
short `Vmax=80`, `step_control=:eta` smoke run reaches the requested
`U=1.6` without a nonfinite row and detects a row-local trapped-surface
crossing, but the late row is already in a floor-dominated throat/interior
region where `Delta U_eta` becomes huge and the max-row cap takes over. This
is useful controller plumbing, not yet a replacement for the conservative
geometric-plus-throat local controller.

The next semi-implicit step is row rejection with backtracking. Passing
`backtrack=true` to `evolve_gp2026_u_adaptive` changes a proposed row advance
from "take the step and inspect it later" to "try the step, measure the
realized future-row change, and reject it if needed." The accept/reject
diagnostic is `realized_row_change_summary(rows, candidate)`, which reports
the actual row-to-row maxima of `r`, `logf`, `rho`, `eta`, and
`H=-4r_Ur_V/f`. The default backtracking caps are the existing
`max_delta_rho` and `max_delta_eta`; `r`, `logf`, and `H` caps are available
but default to `Inf`. A rejected row is retried with
`Delta U <- backtrack_factor * Delta U`, with `backtrack_factor=0.5` and
`max_backtracks=20` by default. This is not a full Newton solve for the
entire row, but it is semi-implicit in the practical sense that the controller
now depends on the candidate future row instead of only on past-row speed
estimates. A short `Vmax=40`, `step_control=:outer`, `backtrack=true` smoke
run with realized caps `Delta rho <= 0.05` and `Delta eta <= 0.02` reaches
`U=-0.5` with finite rows and 23 accepted rows.

The first Berger-Oliger-inspired AMR piece is now implemented as a row local
truncation-error probe. `berger_oliger_row_lte(previous, target_u, ep)`
computes the next row once with a full `Delta U` step and once through the
midpoint with two half steps, then estimates the second-order Richardson LTE
from their difference. `row_lte_error` forms a normalized max norm over
`r`, `logf`, `|Psi|`, `Q`, and `eta` by default, and
`buffered_flag_intervals` clusters the flagged `V` points into candidate
refinement patches. The default clustering is deliberately Hamade-Stewart-like
and returns one buffered patch spanning all flagged points; `cluster=:components`
keeps separate connected components.

The July 2026 BO tests confirm that the full-step/two-half-step discrepancy
has the expected cubic local-defect scaling for a second-order scheme. On the
unit-test grid, halving `Delta U` gives a physical-sector error ratio `8.55`,
or observed power `3.10`; individual `r`, `logf`, `|Psi|`, `Q`, and `eta`
rates are all close to this value. Including the gauge-dependent `A_U` instead
gives power about `1.84` and can flag the entire late-time row, which is why
the default norm now excludes `A_U,A_V`. They remain available explicitly
through the `fields` keyword.

The older single-row prototype remains available for isolated tests.
`berger_oliger_refine_patch` refines one selected parent `V` interval by a
factor of four by default, interpolates the lower row onto that child grid,
evolves the child through the two half-`U` steps, and injects coincident child
values into both the parent midpoint and target rows. It then re-integrates
the coarse `V` suffix after the patch endpoint, following Hamade-Stewart's
downstream synchronization step. Unit tests cover child-grid dimensions,
finite evolution, exact coincident-point injection, and suffix re-integration.
A factor-four smoke over parent points `100:180` produces 321 child points, a
finite corrected parent row, and maximum raw physical correction `1.93e-7`.

### Persistent Hamade-Stewart hierarchy

`src/StewartAMR.jl` now implements the characteristic hierarchy algorithm of
Hamadé and Stewart, section 3 of `gr-qc/9506044`, extended to multiple
buffered sibling patches per parent. It supports:

1. persistent child grids between revisions;
2. factor-four refinement in both `U` and `V` by default;
3. revision after a configurable fixed number of level steps;
4. dynamic child creation, rebuilding, and destruction;
5. disjoint connected-component sibling patches with configurable gap merging
   and a patch-count cap;
6. recursively many child generations;
7. parent-to-child boundary interpolation during `U` subcycling;
8. left-to-right finest-to-coarsest injection at synchronized `U`;
9. outward `V` reintegration after every injected patch endpoint;
10. optional rejection when the finest-level LTE still exceeds tolerance.

`evolve_gp2026_u_adaptive(...; bo_amr=true)` now uses this persistent
hierarchy. The disposable `advance_u_row_berger_oliger` path is retained only
as a local regression utility. Unit tests exercise two-level persistence,
multiple evolved siblings, recursive subcycling, patch rebuilding and
destruction, coincident-point injection, causal suffix reintegration,
finest-LTE rejection, topology-only rollback, and driver integration.

`examples/convergence_gp2026_stewart_amr.jl` evolves exact extremal RN through
`U=0`. With one factor-four child level it gives:

| root `Delta U` | root `Delta V` | max `|delta r|` | rate | max `|delta log(f)|` | rate |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.040 | 0.200 | 8.610e-6 | - | 3.331e-4 | - |
| 0.020 | 0.100 | 2.153e-6 | 2.000 | 8.326e-5 | 2.000 |
| 0.010 | 0.050 | 5.383e-7 | 2.000 | 2.081e-5 | 2.000 |

The electrovacuum charge error remains below `7e-14`.
`examples/convergence_gp2026_stewart_charge.jl` forces an active child on the
charged pulse and finds second-order independent Maxwell residuals:

| root `Delta U` | root `Delta V` | max `Q_U` residual | rate | max `Q_V` residual | rate |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.020 | 0.160 | 1.594e-6 | - | 6.573e-6 | - |
| 0.010 | 0.080 | 3.997e-7 | 1.996 | 1.645e-6 | 1.999 |
| 0.005 | 0.040 | 1.002e-7 | 1.995 | 4.113e-7 | 2.000 |

The Stewart LTE norm uses `atol + rtol*max(|y_coarse|,|y_fine|)`, without
the previous artificial unit floor that made small throat variables such as
`eta=1-|Q|/r` invisible. `examples/check_gp2026_stewart_amr.jl` reports
per-level LTE, patch ranges, revisions, injections, and reintegrations.

On quoted-critical data at `Vmax=400`, `Delta V=0.08`, and a 1200-row budget,
the active setting `(atol,rtol)=(1e-10,1e-8)` performs 420 injections,
24 downstream reintegrations, and 105 dynamic child create/destroy cycles,
reaching three levels and `rho=9.48`. It remains finite but does not yet form
a trapped row; the closest positive expansion is `r_V=3.76e-3` at `V=17.12`.
A looser 6000-row `(1e-8,1e-5)` run does not request late-time refinement and
also remains untrapped. The hierarchy implementation therefore passed its
algorithmic and convergence checks, but the hierarchy alone did not resolve
the GP physics discrepancy.

The missing numerical operation was transactional rejection of a failed
coarse root step. `advance_stewart_hierarchy!` now evolves a trial hierarchy,
halves `Delta U` after a Newton failure, and commits only a completely
synchronized recursive step. This prevents a failed coarse solve from
terminating the evolution before the active child can resolve the throat.
Using the literal outer-boundary Eq. (9) controller,
`(atol,rtol)=(1e-10,1e-8)`, three levels, and `Vmax=40`, the complete
evolution to `Umax=1.6` now gives:

| `Q0` | `Delta V` | `Vtrap` | `U` at minimum crossing |
| ---: | ---: | ---: | ---: |
| 1.001 | 0.080 | 7.740115 | 0.037091 |
| 1.002 | 0.080 | 10.292215 | 0.052061 |
| 1.003 | 0.080 | 12.316536 | 0.176653 |
| 1.0032 | 0.080 | 12.662902 | 0.182329 |

The monotonic increase of `Vtrap` as `Q0` approaches the quoted threshold
from the black-hole side is the expected Section IIIA trend. For `Q0=1.001`,
the joint `Delta V=0.08,0.04,0.02` sequence is
`7.7401149, 7.7393168, 7.7391419`; the two error differences have ratio
`4.56`, corresponding to observed order `2.19`. These are dynamically
formed horizons because the initial corner is super-extremal. The reported
`Vtrap` is the minimum crossing over all evolved rows, as in GP; stopping at
the first trapped row gives a different and substantially larger value.

This recovers horizon formation but does not yet reproduce the paper's
quoted threshold quantitatively. At `Q0=1.0033218`, the minimum crossing is
`12.88002` at `Delta V=0.08` and `12.87952` at `Delta V=0.04`, so the finite
crossing is not explained by root-grid truncation error. At `Delta V=0.08`,
`Q0=1.005` traps at `V=16.55368`, whereas `Q0=1.006` is untrapped through
`Umax=1.6` on `Vmax=40`; the finite-domain effective threshold is therefore
shifted into that interval.
A fixed `p=1/2` fit to the current low-`Vtrap` points gives approximately
`Qstar=1.0049`, but these crossings lie in the pulse transient and should not
be used as a final critical fit.

Hierarchy-depth saturation is visible near the quoted threshold. Stopping at
the first trapped row gives `V=16.9981` with three levels and `V=17.0750`
with four factor-four levels. The crossing survives the extra generation,
but the fourth-level normalized LTE is still about `15.6`; its precise
location is not tolerance-resolved.

The cell Newton path now owns one reusable workspace per row advance. Residual
evaluation, finite-difference Jacobian construction, damping trials, and the
pivoted `7x7` linear solve all operate in place. A 500-cell direct-Newton row
allocates about 177 bytes per cell including construction of the output row.
Rejected Stewart root steps no longer `deepcopy` all field vectors: a
topology-only snapshot shares the immutable completed `NLRow` objects and
copies only hierarchy nodes and statistics.

The quoted-critical three-level run, including 195 rejected root steps,
completes in `27.9 s`, allocates `4.30 GB` across all recursively generated
rows, and reproduces `Vtrap=12.8800186240`. The four-level
stop-at-first-crossing run completes in `136.7 s` and reproduces
`V=17.07502566`; this was roughly twice as fast as the previous observed
runtime. Its `20.9 GB` allocation total is now dominated by the millions of
intermediate full/half rows and LTE arrays rather than cell-local Newton
scratch space.

### Sibling patches and acceptance controls

LTE flags are now split into buffered connected components. Overlapping or
nearby components are merged and at most eight sibling patches are retained
by default. Siblings are advanced and injected from small to large `V`; after
each injection the downstream coarse suffix is re-integrated before the next
sibling obtains its parent boundary data. This preserves the causal ordering
of the characteristic equations.

Localization is conditional on the error actually being local. On the first
GP step with the strict all-field norm, `r`, `log(f)`, and `eta` are above
tolerance at 500 of 501 points, so a full-domain child is required. At later
times the `Q0=1.001`, `Delta V=0.08` run splits its finest level into two
patches and uses 5,998 points instead of 8,001. At `Delta V=0.04`, the
geometry error remains global and the code correctly retains a full patch.
A charge-only negative control localizes to one 289-point child and runs in
8.8 seconds, but moves the first crossing from `U≈0.023` to `U≈-0.056`; it is
therefore not a valid production refinement norm.

The optional `reject_on_finest_lte` setting rejects the complete root step
when a revision on the maximum level remains above tolerance. The hierarchy
is rolled back and the root step is halved transactionally. The current
`Q0=1.001` control results are:

| control | `Vtrap` | conclusion |
| :--- | ---: | :--- |
| all fields, `Delta V=0.08,0.04,0.02`, strict norm without LTE rejection | `7.739510, 7.739158, 7.739127` | strong root-`V` convergence |
| two levels, `Delta V=0.08` | `7.107288` | rejected as under-resolved |
| three levels, tolerance relaxed by ten | `7.747665` | `1.1e-3` fractional tolerance effect |
| accepted LTE, `C=0.6`, `Delta V=0.08,0.04` | `7.72110, 7.71993` | `1.6e-4` root-`V` effect |
| accepted LTE, `Vmax=40,80` | `7.72110, 7.72214` | `1.3e-4` outer-boundary effect |
| accepted LTE, `C=0.6,0.4` | `7.72110, 7.76051` | `5.1e-3` controller effect; dominant error |

The accepted runs have maximum finest-level normalized LTE below one. The
strict `(atol,rtol)=(1e-10,1e-8)` all-field acceptance run is prohibitively
expensive at three levels; the control table uses `(1e-9,1e-7)`. The remaining
`C` dependence means that horizon formation is resolved qualitatively but
`Vtrap` is not yet a sub-percent precision observable.

`refined_vtrap_sample` fits the nonuniform three-row apparent-horizon crossing
curve near its minimum. `trapped_surface_invariants` then uses the
marginal-surface identity `M=(r^2+Q^2)/(2r)` to report gauge-invariant
quantities without coordinate derivatives. For the accepted `Q0=1.001`
runs, the nucleating surface has

```text
M = 1.00290--1.00294,
1-Q/M = 3.2e-11--5.3e-11,
|1-r/M| = 0.8e-5--1.0e-5.
```

This confirms formation through an almost degenerate marginal pair. The
Reissner-Nordstrom surface-gravity proxy is `0.8e-5--1.0e-5`, but its
remaining controller sensitivity prevents using it for a critical scaling
fit. Conservative refluxing is still absent because the characteristic cell
equations are not finite-volume flux balances.

`examples/check_charged_horizon_density.jl` is the charged-sector target
from Gelles/Pretorius. For extremal `eQ0=0.6`, the expected late-time
horizon charge-density exponent is `1 - 2s = 0`, i.e. a plateau. The
current scaffold does not pass this check yet; it is kept as a research
diagnostic to drive the next round of charged-sector corrections.

The Baake/Rinne equations cannot be pasted directly because their variables and gauge are CMC hyperboloidal, not compactified double-null. They are still the right source for matter stress tensor, charge conventions, and comparison diagnostics.
