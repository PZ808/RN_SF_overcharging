using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))
string_argument(index, default) = argument(index, default)

function safe_fit_power_law(x, y; xmin=nothing, xmax=nothing)
    try
        return fit_power_law(x, y; xmin, xmax)
    catch err
        err isa ArgumentError || rethrow()
        return NaN, NaN, 0
    end
end

q0 = real_argument(1, 1.0)
eQ0 = real_argument(2, 0.4)
amplitude = real_argument(3, 1.0e-4)
vmax = real_argument(4, 600.0)
nu = integer_argument(5, 260)
nv = integer_argument(6, 700)
u1 = real_argument(7, -1.0e-4)
extractor_mode = string_argument(8, "row")
envelope_mode = string_argument(9, "gp-bump")
pulse_width = real_argument(10, 4.0)
fit_vmin = real_argument(11, 40.0)
fit_vmax = real_argument(12, 150.0)

extractor_mode in ("row", "extrapolated") ||
    throw(ArgumentError("extractor_mode must be row or extrapolated"))

ep = EvolutionParams(
    rn=RNParams(1.0, q0),
    scalar_charge=eQ0 / q0,
    amplitude=amplitude,
    omega=0.0,
    center=20.0,
    width=pulse_width,
)

envelope = if envelope_mode == "gaussian"
    gaussian_envelope
elseif envelope_mode == "gp-bump"
    gp2025_bump_envelope
else
    throw(ArgumentError("envelope_mode must be gaussian or gp-bump"))
end

grid = compact_mrt_grid(
    ep.rn;
    nu,
    nv,
    u0=-1.0,
    u1,
    v0=compact_v_from_ef_v(0.0, ep.rn),
    v1=compact_v_from_ef_v(vmax, ep.rn),
)
state = initialize_state(grid, ep; envelope)
evolve_passive_scalar!(state, grid, ep)
reconstruct_passive_charge!(state, grid, ep)

vef, rho = if extractor_mode == "row"
    horizon_charge_density_series(state, grid, ep)
else
    horizon_charge_density_extrapolated_series(state, grid, ep; rows=4)
end
finite = isfinite.(vef) .& isfinite.(rho) .& (vef .> 0) .& (abs.(rho) .> 0)
vef = vef[finite]
rho = rho[finite]
rho_abs = abs.(rho)
s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)
slope, _, nfit = safe_fit_power_law(vef, rho_abs; xmin=fit_vmin, xmax=fit_vmax)
mask = (vef .>= fit_vmin) .& (vef .<= fit_vmax)

println("# passive charge-density audit")
println("# q0 = ", q0, ", eQ0 = ", eQ0, ", amplitude = ", amplitude,
        ", Vmax = ", vmax, ", nu = ", nu, ", nv = ", nv,
        ", u1 = ", u1, ", extractor_mode = ", extractor_mode,
        ", envelope_mode = ", envelope_mode, ", pulse_width = ", pulse_width)
println("s = ", s)
println("target rho_Q slope 1 - 2s = ", 1 - 2s)
println("fit window V_EF = ", (fit_vmin, fit_vmax))
println("fit samples = ", nfit)
println("passive |rho_Q| slope = ", slope)
println("late-window max/min = ",
        count(mask) == 0 ? NaN : maximum(rho_abs[mask]) / minimum(rho_abs[mask]))
println("max |Q-Q0| = ", maximum(abs.(state.Q .- ep.rn.Q0)))
for target_v in (40.0, 80.0, 120.0, 150.0, 200.0, 300.0, 600.0)
    target_v <= maximum(vef) || continue
    _, j = findmin(abs.(vef .- target_v))
    println("  V=", vef[j], " rho_Q=", rho[j])
end

