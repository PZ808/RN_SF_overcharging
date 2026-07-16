using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))
string_argument(index, default) = argument(index, default)

q0 = real_argument(1, 1.0)
eQ0 = real_argument(2, 0.4)
amplitude = real_argument(3, 1.0e-4)
vmax = real_argument(4, 600.0)
nu = integer_argument(5, 260)
nv = integer_argument(6, 700)
u1 = real_argument(7, -1.0e-4)
mode = string_argument(8, "passive")
envelope_mode = string_argument(9, "gp-bump")
pulse_width = real_argument(10, 4.0)

mode in ("passive", "coupled") ||
    throw(ArgumentError("mode must be passive or coupled"))

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
if mode == "passive"
    evolve_passive_scalar!(state, grid, ep)
else
    evolve!(state, grid, ep)
end

vef = [NonlinearEMKGDoubleNull.ef_v_from_mrt(v, ep.rn) for v in grid.v]
P = sqrt.(state.xi[end, :] .^ 2 .+ state.pi[end, :] .^ 2)
s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)

println("# fixed-background scalar-tail audit")
println("# q0 = ", q0, ", eQ0 = ", eQ0, ", amplitude = ", amplitude,
        ", Vmax = ", vmax, ", nu = ", nu, ", nv = ", nv,
        ", u1 = ", u1, ", mode = ", mode,
        ", envelope_mode = ", envelope_mode, ", pulse_width = ", pulse_width)
println("s = ", s)
println("target horizon |r phi| slope = ", -s)
println("max |Q-Q0| = ", maximum(abs.(state.Q .- ep.rn.Q0)))
for (lo, hi) in ((40.0, 150.0), (80.0, 300.0), (150.0, 300.0),
                 (150.0, 600.0))
    mask = isfinite.(vef) .& isfinite.(P) .& (vef .>= lo) .&
           (vef .<= hi) .& (P .> 0)
    count(mask) < 4 && continue
    slope, _, nfit = fit_power_law(vef[mask], P[mask])
    println("fit V_EF = ", (lo, hi),
            ": samples = ", nfit,
            ", slope = ", slope,
            ", max/min = ", maximum(P[mask]) / minimum(P[mask]))
end

