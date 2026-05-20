using NonlinearEMKGDoubleNull

function quantile_like(x, q)
    xs = sort(collect(x))
    k = clamp(Int(round(1 + q * (length(xs) - 1))), 1, length(xs))
    return xs[k]
end

# Nontrivial physics regression target:
# in the e=0, real-scalar, extremal-RN limit, Fig. 8 and the associated
# discussion of arXiv:1307.6800 imply late-time horizon decay phi ~ V^-1
# for the scalar itself.

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = 1.0e-4,
    omega = 0.0,
    center = 15.0,
    width = 3.0,
)

grid = mrt2013_grid(; nu=180, nv=900, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=120.0)
state = NLState(grid)
initialize_mrt2013_uncharged_ingoing!(state, grid, ep)

evolve_nonlinear!(state, grid, ep; iterations=5)

vef = grid.v
phi_abs = abs.(state.phi_re[end, :])
finite = isfinite.(vef) .& isfinite.(phi_abs) .& (vef .> 0)
vef = vef[finite]
phi_abs = phi_abs[finite]

vmin = quantile_like(vef, 0.65)
slope, intercept, nfit = fit_power_law(vef, phi_abs; xmin=vmin)

println("samples = ", length(vef))
println("fit samples = ", nfit)
println("fit window V_EF >= ", vmin)
println("late-time horizon-adjacent |phi| slope = ", slope)
println("target slope for phi on extremal horizon: -1")
