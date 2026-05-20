using NonlinearEMKGDoubleNull

function quantile_like(x, q)
    xs = sort(collect(x))
    k = clamp(Int(round(1 + q * (length(xs) - 1))), 1, length(xs))
    return xs[k]
end

# Charged-sector research diagnostic from Gelles/Pretorius arXiv:2503.04881:
# on extremal RN, the horizon charge density scales as
#   rho_Q ~ V_EF^(1 - 2s)
# and for |e Q0| >= 1/2, s=1/2 so rho_Q approaches a nonzero constant.
#
# This is intentionally not a unit test yet. At the current scaffold stage it
# is expected to expose whether the charged horizon evolution is faithful.

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.6,
    amplitude = 1.0e-4,
    omega = 0.0,
    center = 20.0,
    width = 4.0,
)

v0 = compact_v_from_ef_v(0.0, ep.rn)
v1 = compact_v_from_ef_v(900.0, ep.rn)
grid = compact_mrt_grid(ep.rn; nu=90, nv=360, u0=-1.0, v0=v0, u1=-1.0e-4, v1=v1)
state = initialize_state(grid, ep)
evolve!(state, grid, ep)
ru, rv = maxwell_residuals(state, grid, ep)

vef, rho = horizon_charge_density_series(state, grid, ep)
finite = isfinite.(vef) .& isfinite.(rho) .& (vef .> 0) .& (abs.(rho) .> 0)
vef = vef[finite]
rho_abs = abs.(rho[finite])

vmin = max(40.0, quantile_like(vef, 0.65))
slope, intercept, nfit = fit_power_law(vef, rho_abs; xmin=vmin)
s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)
target = 1 - 2s
late = rho_abs[vef .>= vmin]
plateau_ratio = maximum(late) / minimum(late)

println("samples = ", length(vef))
println("fit samples = ", nfit)
println("fit window V_EF >= ", vmin)
println("eQ0 = ", ep.scalar_charge * ep.rn.Q0)
println("s = ", s)
println("late-time |rho_Q| slope = ", slope)
println("linear target slope 1 - 2s = ", target)
println("late-window plateau max/min = ", plateau_ratio)
println("max |Maxwell-U residual| = ", maximum(abs, ru))
println("max |Maxwell-V residual| = ", maximum(abs, rv))
