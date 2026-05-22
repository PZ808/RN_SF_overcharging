using NonlinearEMKGDoubleNull

function run_vacuum_rn_case(nu, nv; V1=20.0)
    ep = EvolutionParams(
        rn = RNParams(1.0, 1.0),
        scalar_charge = 0.0,
        amplitude = 0.0,
        omega = 0.0,
        center = 15.0,
        width = 3.0,
    )
    grid = mrt2013_grid(; nu=nu, nv=nv, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=V1)
    state = NLState(grid)
    initialize_mrt2013_uncharged_ingoing!(state, grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations=20, subtract_rn_background=true)

    exact_r = [mrt2013_areal_radius(U, V, ep.rn) for U in grid.u, V in grid.v]
    exact_logf = [log(mrt2013_metric_f(U, V, ep.rn)) for U in grid.u, V in grid.v]
    return maximum(abs.(state.r .- exact_r)), maximum(abs.(state.logf .- exact_logf))
end

resolutions = ((60, 180), (120, 360), (240, 720))
errors = [run_vacuum_rn_case(nu, nv) for (nu, nv) in resolutions]

println("resolutions = ", resolutions)
for (res, err) in zip(resolutions, errors)
    println("resolution ", res, ": max |r-r_RN| = ", err[1],
            ", max |logf-logf_RN| = ", err[2])
end

if all(err -> err[1] < 1.0e-9 && err[2] < 1.0e-9, errors)
    println("errors are at roundoff; refinement ratios are not meaningful")
else
    println("r error ratios    = ", errors[1][1] / errors[2][1], ", ",
            errors[2][1] / errors[3][1])
    println("logf error ratios = ", errors[1][2] / errors[2][2], ", ",
            errors[2][2] / errors[3][2])
end
