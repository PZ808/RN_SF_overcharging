using NonlinearEMKGDoubleNull

ep = EvolutionParams(
    rn = RNParams(1.0, 0.999),
    scalar_charge = 0.4,
    amplitude = 1.0e-4,
    omega = 0.4,
    center = 20.0,
    width = 4.0,
)

grid = compact_mrt_grid(ep.rn; nu=48, nv=48, u0=-1.0, v0=0.05, u1=-1.0e-3, v1=pi / 2 - 1.0e-2)
state = initialize_state(grid, ep)
evolve!(state, grid, ep)
ru, rv = maxwell_residuals(state, grid, ep)

println("grid = ", size(grid))
println("Q range = ", extrema(state.Q))
println("max |Maxwell-U residual| = ", maximum(abs, ru))
println("max |Maxwell-V residual| = ", maximum(abs, rv))

