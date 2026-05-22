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
        scalar_charge = 0.6,
        amplitude = 1.0e-4,
        omega = 0.0,
        center = 20.0,
        width = 4.0,
    )
    v0 = compact_v_from_ef_v(0.0, ep.rn)
    v1 = compact_v_from_ef_v(80.0, ep.rn)
    grid = compact_mrt_grid(ep.rn; nu=nu, nv=nv, u0=-1.0, v0=v0, u1=-1.0e-4, v1=v1)
    state = initialize_state(grid, ep)
    evolve!(state, grid, ep)
    ru, rv = maxwell_residuals(state, grid, ep)
    amp = sqrt.(state.xi.^2 .+ state.pi.^2)
    return grid, state, amp, maximum(abs, ru), maximum(abs, rv)
end

resolutions = [(45, 80), (90, 160), (180, 320)]
results = [run_case(nu, nv) for (nu, nv) in resolutions]

g1, st1, amp1, ru1, rv1 = results[1]
g2, st2, amp2, ru2, rv2 = results[2]
g3, st3, amp3, ru3, rv3 = results[3]

q2_on_1 = restrict_to(g1, g2, st2.Q)
q3_on_1 = restrict_to(g1, g3, st3.Q)
xi2_on_1 = restrict_to(g1, g2, st2.xi)
xi3_on_1 = restrict_to(g1, g3, st3.xi)
pi2_on_1 = restrict_to(g1, g2, st2.pi)
pi3_on_1 = restrict_to(g1, g3, st3.pi)
amp2_on_1 = restrict_to(g1, g2, amp2)
amp3_on_1 = restrict_to(g1, g3, amp3)

q_err_12 = maximum(abs.(st1.Q .- q2_on_1))
q_err_23 = maximum(abs.(q2_on_1 .- q3_on_1))
xi_err_12 = maximum(abs.(st1.xi .- xi2_on_1))
xi_err_23 = maximum(abs.(xi2_on_1 .- xi3_on_1))
pi_err_12 = maximum(abs.(st1.pi .- pi2_on_1))
pi_err_23 = maximum(abs.(pi2_on_1 .- pi3_on_1))
amp_err_12 = maximum(abs.(amp1 .- amp2_on_1))
amp_err_23 = maximum(abs.(amp2_on_1 .- amp3_on_1))

interior = (2:length(g1.u)-1, 2:length(g1.v)-1)
q_int_12 = maximum(abs.(st1.Q[interior...] .- q2_on_1[interior...]))
q_int_23 = maximum(abs.(q2_on_1[interior...] .- q3_on_1[interior...]))
xi_int_12 = maximum(abs.(st1.xi[interior...] .- xi2_on_1[interior...]))
xi_int_23 = maximum(abs.(xi2_on_1[interior...] .- xi3_on_1[interior...]))
pi_int_12 = maximum(abs.(st1.pi[interior...] .- pi2_on_1[interior...]))
pi_int_23 = maximum(abs.(pi2_on_1[interior...] .- pi3_on_1[interior...]))
amp_int_12 = maximum(abs.(amp1[interior...] .- amp2_on_1[interior...]))
amp_int_23 = maximum(abs.(amp2_on_1[interior...] .- amp3_on_1[interior...]))

println("resolutions = ", resolutions)
println("Q max diff coarse-medium = ", q_err_12)
println("Q max diff medium-fine   = ", q_err_23)
println("Q convergence ratio      = ", q_err_12 / q_err_23)
println("xi max diff coarse-medium = ", xi_err_12)
println("xi max diff medium-fine   = ", xi_err_23)
println("xi convergence ratio      = ", xi_err_12 / xi_err_23)
println("Pi max diff coarse-medium = ", pi_err_12)
println("Pi max diff medium-fine   = ", pi_err_23)
println("Pi convergence ratio      = ", pi_err_12 / pi_err_23)
println("P max diff coarse-medium = ", amp_err_12)
println("P max diff medium-fine   = ", amp_err_23)
println("P convergence ratio      = ", amp_err_12 / amp_err_23)
println("Q interior ratio         = ", q_int_12 / q_int_23)
println("xi interior ratio        = ", xi_int_12 / xi_int_23)
println("Pi interior ratio        = ", pi_int_12 / pi_int_23)
println("P interior ratio         = ", amp_int_12 / amp_int_23)
println("Maxwell residuals U      = ", (ru1, ru2, ru3))
println("Maxwell residuals V      = ", (rv1, rv2, rv3))
