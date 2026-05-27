using NonlinearEMKGDoubleNull

function real_argument(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

function first_trapped_slice(state)
    for j in 2:length(state.slices)
        lower, upper = state.slices[j - 1], state.slices[j]
        fields_finite = all(isfinite, upper.r) && all(isfinite, upper.logf) &&
                        all(isfinite, upper.Q)
        fields_finite || return nothing
        rv = adaptive_outgoing_expansion(lower, upper)
        crossing = findfirst(value -> value <= 0, rv)
        if !isnothing(crossing)
            return (upper.v, upper.u[crossing])
        end
    end
    return nothing
end

function last_finite_index(state)
    return findlast(slice -> all(isfinite, slice.r) && all(isfinite, slice.logf) &&
                            all(isfinite, slice.Q), state.slices)
end

q0 = real_argument(1, 1.0033218)
vmax = real_argument(2, 100.0)
dv = real_argument(3, 0.08)
amplitude = real_argument(4, 0.01)
du = real_argument(5, 0.01)
use_chopping = length(ARGS) < 6 || ARGS[6] == "chop"

U0 = -1.0
U1 = 1.6
ep = EvolutionParams(
    rn = RNParams(1.0, q0),
    scalar_charge = 0.6 / q0,
    amplitude = amplitude,
    omega = 1.0,
)
nu = Int(round((U1 - U0) / du)) + 1
nv = Int(round(vmax / dv)) + 1
grid = gp2026_grid(; nu, nv, U0, V0=0.0, U1, V1=vmax)
seed = NLState(grid)
initialize_gp2026_single_pulse!(seed, grid, ep)

state = if use_chopping
    chopping = HorizonChoppingConfig(
        start_v = dv,
        band_width = dv,
        interior_buffer_cells = 1,
        interior_buffer_width = du,
    )
    evolve_adaptive(
        slice_from_rectangular(seed, grid, 1),
        west_boundary_from_rectangular(seed, grid),
        ep;
        iterations=10,
        horizon_chopping=chopping,
        reduced_scalar=true,
    )
else
    evolve_nonlinear!(seed, grid, ep; iterations=10, reduced_scalar=true)
    adaptive_state_from_rectangular(seed, grid)
end

finite_index = last_finite_index(state)
isnothing(finite_index) && error("initial slice is nonfinite")
first_bad_v = finite_index == length(state.slices) ? nothing :
              state.slices[finite_index + 1].v
finite_state = AdaptiveNLState(state.slices[1:finite_index])
trap = first_trapped_slice(finite_state)
last_v = last(finite_state.slices).v

println("Gelles-Pretorius 2026 longer nonlinear evolution")
println("Q0 = ", q0)
println("e Q0 = ", ep.scalar_charge * q0)
println("A0 = ", amplitude)
println("requested Vmax = ", vmax)
println("Delta U = ", du)
println("Delta V = ", dv)
println("horizon chopping = ", use_chopping)
println("first trapped (V,U-cell) = ", trap)
println("first nonfinite V = ", first_bad_v)
println("last finite V = ", last_v)
println("last finite U points = ", length(last(finite_state.slices).u))
println("last retained U = ", last(last(finite_state.slices).u))

if length(finite_state.slices) >= 2
    target_v = last_v - dv / 2
    target_u = min(-0.5, minimum(last(slice.u) for slice in finite_state.slices))
    u, sampled_v, rv = outgoing_expansion_profile(finite_state; target_v)
    _, _, mass = renormalized_hawking_mass_profile(finite_state.slices[end - 1],
                                                   finite_state.slices[end])
    _, _, _, _, q_u_residual =
        charged_charge_flux_u_profile(finite_state, ep; target_v, reduced_scalar=true)
    _, _, _, _, q_v_residual =
        charged_charge_flux_v_profile(finite_state, ep; target_u, reduced_scalar=true)
    _, final_rphi = gp2026_rphi_profile(last(finite_state.slices))
    horizon_v, _, horizon_rphi = gp2026_horizon_rphi_series(finite_state)
    println("diagnostic V = ", sampled_v)
    println("diagnostic U for Q_V = ", target_u)
    println("r_V extrema = ", extrema(rv))
    println("apparent-horizon U = ", apparent_horizon_location(u, rv))
    println("mass extrema = ", extrema(mass))
    println("last-slice max |r phi_GP| = ", maximum(final_rphi))
    if !isempty(horizon_rphi)
        println("horizon |r phi_GP| first/last = ",
                (first(horizon_rphi), last(horizon_rphi)))
        println("horizon amplitude samples and V range = ",
                (length(horizon_rphi), first(horizon_v), last(horizon_v)))
    end
    println("max |Q_U - flux source| = ", maximum(abs, q_u_residual))
    println("max |Q_V - flux source| = ", maximum(abs, q_v_residual))
end
