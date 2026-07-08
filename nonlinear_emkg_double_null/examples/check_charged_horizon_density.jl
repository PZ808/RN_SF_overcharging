using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))
string_argument(index, default) = argument(index, default)

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

q0 = real_argument(1, 1.0)
eQ0 = real_argument(2, 0.6)
amplitude = real_argument(3, 1.0e-4)
vmax = real_argument(4, 900.0)
nu = integer_argument(5, 90)
nv = integer_argument(6, 360)
fit_quantile = real_argument(7, 0.65)
fit_vmin_floor = real_argument(8, 40.0)
u1 = real_argument(9, -1.0e-4)
fit_xmax_quantile = real_argument(10, 0.95)
grid_mode = string_argument(11, "compact")

ep = EvolutionParams(
    rn = RNParams(1.0, q0),
    scalar_charge = eQ0 / q0,
    amplitude = amplitude,
    omega = 0.0,
    center = 20.0,
    width = 4.0,
)

grid = if grid_mode == "compact"
    v0 = compact_v_from_ef_v(0.0, ep.rn)
    v1 = compact_v_from_ef_v(vmax, ep.rn)
    compact_mrt_grid(ep.rn; nu, nv, u0=-1.0, v0=v0, u1, v1)
elseif grid_mode == "ef-uniform"
    ef_v_mrt_grid(ep.rn; nu, nv, u0=-1.0, u1, Vef0=0.0, Vef1=vmax)
else
    throw(ArgumentError("grid_mode must be compact or ef-uniform"))
end
state = initialize_state(grid, ep)
evolve!(state, grid, ep)
ru, rv = maxwell_residuals(state, grid, ep)

vef, rho = horizon_charge_density_series(state, grid, ep)
finite = isfinite.(vef) .& isfinite.(rho) .& (vef .> 0) .& (abs.(rho) .> 0)
vef = vef[finite]
rho = rho[finite]
rho_abs = abs.(rho)

vmin = max(fit_vmin_floor, quantile_like(vef, fit_quantile))
vmax_fit = quantile_like(vef, fit_xmax_quantile)
slope, intercept, nfit = fit_power_law(vef, rho_abs; xmin=vmin, xmax=vmax_fit)
s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)
target = 1 - 2s
late_mask = (vef .>= vmin) .& (vef .<= vmax_fit)
late = rho_abs[late_mask]
plateau_ratio = maximum(late) / minimum(late)
late_signed = rho[late_mask]
dominant_sign = count(>=(0), late_signed) >= count(<(0), late_signed) ? 1.0 : -1.0
sign_mask = late_mask .& (dominant_sign .* rho .> 0)
signed_slope = count(sign_mask) >= 2 ?
               fit_power_law(vef[sign_mask], dominant_sign .* rho[sign_mask])[1] :
               NaN
signed_late = dominant_sign .* rho[sign_mask]
signed_ratio = isempty(signed_late) ? NaN : maximum(signed_late) / minimum(signed_late)

println("# q0 = ", q0, ", eQ0 = ", eQ0, ", amplitude = ", amplitude,
        ", Vmax = ", vmax, ", nu = ", nu, ", nv = ", nv,
        ", fit_quantile = ", fit_quantile, ", fit_vmin_floor = ",
        fit_vmin_floor, ", u1 = ", u1, ", fit_xmax_quantile = ",
        fit_xmax_quantile, ", grid_mode = ", grid_mode)
println("samples = ", length(vef))
println("fit samples = ", nfit)
println("fit window V_EF = ", (vmin, vmax_fit))
println("eQ0 = ", ep.scalar_charge * ep.rn.Q0)
println("s = ", s)
println("late-time |rho_Q| slope = ", slope)
println("late-time sign-coherent rho_Q slope = ", signed_slope)
println("linear target slope 1 - 2s = ", target)
println("late-window plateau max/min = ", plateau_ratio)
println("sign-coherent late-window max/min = ", signed_ratio)
println("max |Maxwell-U residual| = ", maximum(abs, ru))
println("max |Maxwell-V residual| = ", maximum(abs, rv))
