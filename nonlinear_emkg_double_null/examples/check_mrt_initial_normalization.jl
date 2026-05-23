using NonlinearEMKGDoubleNull

epsilon = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.02
resolutions = length(ARGS) >= 2 ? parse.(Int, ARGS[2:end]) :
              [110, 107, 531, 1061, 2651, 5301]

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = epsilon,
)
f0 = mrt2013_degenerate_horizon_f0(ep)
initial_mass = mrt2013_initial_bondi_mass(f0)

println("epsilon = ", epsilon)
println("(f0 - 2) / epsilon^2 = ", (f0 - 2) / epsilon^2)
println("(Mi - 1) / epsilon^2 = ", (initial_mass - 1) / epsilon^2)
println("MRT small-epsilon targets: 0.740 and 0.789")
println("discrete Eq. (28) reconstruction on the seeded V=0 leg:")

for nu in resolutions
    grid = mrt2013_grid(; nu, nv=2, U0=-5.1, V0=0.0, U1=0.2, V1=0.02)
    state = NLState(grid)
    initialize_mrt2013_outgoing_wave!(state, grid, ep; f0)
    rv = mrt2013_initial_rv_profile(state, grid)
    k = argmin(rv)
    zero_index = findfirst(iszero, grid.u)
    at_zero = isnothing(zero_index) ? "not sampled" :
              string(rv[zero_index] / epsilon^2)
    println("  nu = ", nu,
            ", Delta U = ", grid.u[2] - grid.u[1],
            ", U(min r_V) = ", grid.u[k],
            ", min(r_V)/epsilon^2 = ", rv[k] / epsilon^2,
            ", r_V(0)/epsilon^2 = ", at_zero)
end
