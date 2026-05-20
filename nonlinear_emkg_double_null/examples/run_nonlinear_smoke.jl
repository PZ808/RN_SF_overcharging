using NonlinearEMKGDoubleNull

ep = EvolutionParams(
    rn = RNParams(1.0, 0.999),
    scalar_charge = 0.6,
    amplitude = 1.0e-5,
    omega = 0.6,
    center = 20.0,
    width = 4.0,
)

grid = compact_mrt_grid(ep.rn; nu=32, nv=32, u0=-1.0, v0=0.05, u1=-1.0e-2, v1=1.25)
state = initialize_nonlinear_state(grid, ep)
evolve_nonlinear!(state, grid, ep)

println("grid = ", size(grid))
println("r range = ", extrema(state.r))
println("f range = ", extrema(exp.(state.logf)))
println("Q range = ", extrema(state.Q))
println("max |phi| = ", maximum(sqrt.(state.phi_re.^2 .+ state.phi_im.^2)))

