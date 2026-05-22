using NonlinearEMKGDoubleNull

function restrict_to(coarse_grid, fine_grid, fine_values)
    out = similar(coarse_grid.u, length(coarse_grid.u), length(coarse_grid.v))
    for i in eachindex(coarse_grid.u)
        ui = searchsortedlast(fine_grid.u, coarse_grid.u[i])
        ui = clamp(ui, firstindex(fine_grid.u), lastindex(fine_grid.u) - 1)
        u0, u1 = fine_grid.u[ui], fine_grid.u[ui + 1]
        tu = iszero(u1 - u0) ? zero(u0) : (coarse_grid.u[i] - u0) / (u1 - u0)
        for j in eachindex(coarse_grid.v)
            vj = searchsortedlast(fine_grid.v, coarse_grid.v[j])
            vj = clamp(vj, firstindex(fine_grid.v), lastindex(fine_grid.v) - 1)
            v0, v1 = fine_grid.v[vj], fine_grid.v[vj + 1]
            tv = iszero(v1 - v0) ? zero(v0) : (coarse_grid.v[j] - v0) / (v1 - v0)
            f00 = fine_values[ui, vj]
            f10 = fine_values[ui + 1, vj]
            f01 = fine_values[ui, vj + 1]
            f11 = fine_values[ui + 1, vj + 1]
            out[i, j] = (1 - tu) * (1 - tv) * f00 + tu * (1 - tv) * f10 +
                        (1 - tu) * tv * f01 + tu * tv * f11
        end
    end
    return out
end

function run_case(nu, nv)
    ep = EvolutionParams(
        rn = RNParams(1.0, 1.0),
        scalar_charge = 0.0,
        amplitude = 1.0e-4,
        omega = 0.0,
        center = 15.0,
        width = 3.0,
    )
    grid = mrt2013_grid(; nu=nu, nv=nv, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=20.0)
    state = NLState(grid)
    initialize_mrt2013_uncharged_ingoing!(state, grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations=12, subtract_rn_background=true)
    return grid, state
end

function maxdiff_ratio(a1, a2_on_1, a3_on_1)
    err12 = maximum(abs.(a1 .- a2_on_1))
    err23 = maximum(abs.(a2_on_1 .- a3_on_1))
    return err12, err23, err12 / err23
end

resolutions = [(60, 180), (120, 360), (240, 720)]
results = [run_case(nu, nv) for (nu, nv) in resolutions]

g1, st1 = results[1]
g2, st2 = results[2]
g3, st3 = results[3]

r2 = restrict_to(g1, g2, st2.r)
r3 = restrict_to(g1, g3, st3.r)
lf2 = restrict_to(g1, g2, st2.logf)
lf3 = restrict_to(g1, g3, st3.logf)
phi2 = restrict_to(g1, g2, st2.phi_re)
phi3 = restrict_to(g1, g3, st3.phi_re)

r_err = maxdiff_ratio(st1.r, r2, r3)
lf_err = maxdiff_ratio(st1.logf, lf2, lf3)
phi_err = maxdiff_ratio(st1.phi_re, phi2, phi3)

interior = (2:length(g1.u)-1, 2:length(g1.v)-1)
r_int = maxdiff_ratio(st1.r[interior...], r2[interior...], r3[interior...])
lf_int = maxdiff_ratio(st1.logf[interior...], lf2[interior...], lf3[interior...])
phi_int = maxdiff_ratio(st1.phi_re[interior...], phi2[interior...], phi3[interior...])

println("resolutions = ", resolutions)
println("r max diff coarse-medium    = ", r_err[1])
println("r max diff medium-fine      = ", r_err[2])
println("r convergence ratio         = ", r_err[3])
println("logf max diff coarse-medium = ", lf_err[1])
println("logf max diff medium-fine   = ", lf_err[2])
println("logf convergence ratio      = ", lf_err[3])
println("phi max diff coarse-medium  = ", phi_err[1])
println("phi max diff medium-fine    = ", phi_err[2])
println("phi convergence ratio       = ", phi_err[3])
println("r interior ratio            = ", r_int[3])
println("logf interior ratio         = ", lf_int[3])
println("phi interior ratio          = ", phi_int[3])
