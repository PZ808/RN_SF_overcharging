using Test
using NonlinearEMKGDoubleNull

@testset "RN geometry" begin
    p = RNParams(1.0, 0.8)
    @test metric_F(10.0, p) > 0
    @test isfinite(rstar(3.0, p))
    @test areal_radius(-0.8, 0.2, p) > first(NonlinearEMKGDoubleNull.horizons(p))
    @test isfinite(metric_ftilde(-0.8, 0.2, p))
end

@testset "small evolution" begin
    ep = EvolutionParams(rn=RNParams(1.0, 0.9), scalar_charge=0.1, amplitude=1.0e-5, omega=0.2)
    grid = compact_mrt_grid(ep.rn; nu=8, nv=8, u0=-1.0, v0=0.05, u1=-0.1, v1=0.8)
    state = initialize_state(grid, ep)
    evolve!(state, grid, ep)
    @test all(isfinite, state.xi)
    @test all(isfinite, state.pi)
    @test all(isfinite, state.Q)
end

@testset "small nonlinear scaffold evolution" begin
    ep = EvolutionParams(rn=RNParams(1.0, 0.9), scalar_charge=0.2, amplitude=1.0e-6, omega=0.2)
    grid = compact_mrt_grid(ep.rn; nu=6, nv=6, u0=-1.0, v0=0.05, u1=-0.2, v1=0.7)
    state = initialize_nonlinear_state(grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations=2)
    @test all(isfinite, state.r)
    @test all(isfinite, state.logf)
    @test all(isfinite, state.phi_re)
    @test all(isfinite, state.phi_im)
    @test all(isfinite, state.Q)
end

@testset "exact GP2026 extremal RN horizon crossing" begin
    p = RNParams(1.0, 1.0)
    ep = EvolutionParams(rn=p, scalar_charge=0.0, amplitude=0.0, omega=0.0)

    @test gp2026_exact_extremal_rn_radius(0.0, 7.0, p) == 1.0
    @test gp2026_exact_extremal_rn_fcode(0.0, 7.0, p) == 1.0
    @test gp2026_exact_extremal_rn_radius(-0.2, 0.0, p) ≈ 1.1
    @test gp2026_exact_extremal_rn_radius(0.2, 0.0, p) ≈ 0.9
    @test gp2026_exact_extremal_rn_fcode(-0.2, 0.0, p) ≈ 1.0
    @test gp2026_exact_extremal_rn_fcode(0.2, 0.0, p) ≈ 1.0

    errors = Float64[]
    logf_errors = Float64[]
    horizon_errors = Float64[]
    horizon_f_errors = Float64[]
    for (nu, nv) in ((31, 101), (61, 201))
        grid = gp2026_grid(; nu, nv, U0=-0.4, V0=0.0, U1=0.2, V1=10.0)
        state = NLState(grid)
        initialize_gp2026_exact_extremal_rn!(state, grid, ep)
        evolve_nonlinear!(
            state, grid, ep;
            iterations=15,
            reduced_scalar=true,
            hyperbolic_charge=true,
            cell_solver=:newton_direct,
        )

        exact_r = [
            gp2026_exact_extremal_rn_radius(U, V, p)
            for U in grid.u, V in grid.v
        ]
        exact_logf = [
            log(gp2026_exact_extremal_rn_fcode(U, V, p))
            for U in grid.u, V in grid.v
        ]
        push!(errors, maximum(abs.(state.r .- exact_r)))
        push!(logf_errors, maximum(abs.(state.logf .- exact_logf)))
        horizon_index = findfirst(iszero, grid.u)
        @test !isnothing(horizon_index)
        push!(
            horizon_errors,
            maximum(abs.(state.r[horizon_index, :] .- p.M)),
        )
        push!(
            horizon_f_errors,
            maximum(abs.(exp.(state.logf[horizon_index, :]) .- 1)),
        )
        @test all(isfinite, state.r)
        @test all(isfinite, state.logf)
        @test all(state.Q .== p.Q0)
    end

    @test 3.8 < errors[1] / errors[2] < 4.2
    @test 3.8 < logf_errors[1] / logf_errors[2] < 4.2
    @test 3.8 < horizon_errors[1] / horizon_errors[2] < 4.2
    @test 3.8 < horizon_f_errors[1] / horizon_f_errors[2] < 4.2
end

@testset "charged nonlinear initial Maxwell constraint" begin
    ep = EvolutionParams(rn=RNParams(1.0, 0.999), scalar_charge=0.6,
                         amplitude=1.0e-2, omega=0.6, center=4.0, width=1.0)
    grid = compact_mrt_grid(ep.rn; nu=8, nv=24, u0=-1.0,
                            v0=compact_v_from_ef_v(0.0, ep.rn), u1=-0.1,
                            v1=compact_v_from_ef_v(8.0, ep.rn))
    state = initialize_nonlinear_state(grid, ep)
    i = firstindex(grid.u)
    @test maximum(abs.(state.Q[i, :] .- ep.rn.Q0)) > 0
    checked_cells = 0
    for j in firstindex(grid.v)+1:lastindex(grid.v)
        dv = grid.v[j] - grid.v[j - 1]
        r = (state.r[i, j - 1] + state.r[i, j]) / 2
        f = exp((state.logf[i, j - 1] + state.logf[i, j]) / 2)
        q = state.Q[i, j - 1]
        phi_re = (state.phi_re[i, j - 1] + state.phi_re[i, j]) / 2
        phi_im = (state.phi_im[i, j - 1] + state.phi_im[i, j]) / 2
        phiv_re = (state.phi_re[i, j] - state.phi_re[i, j - 1]) / dv
        phiv_im = (state.phi_im[i, j] - state.phi_im[i, j - 1]) / dv
        Av = (state.Av[i, j - 1] + state.Av[i, j]) / 2
        source = stress_energy(r, f, q, phi_re, phi_im, zero(phiv_re), phiv_re,
                               zero(phiv_im), phiv_im, 0.0, Av, ep.scalar_charge)
        observed = (state.Q[i, j] - state.Q[i, j - 1]) / dv
        expected = -r^2 * source.Jv / 8
        if abs(dv * expected) > 16 * eps(ep.rn.Q0)
            @test observed ≈ expected atol=4 * eps(ep.rn.Q0) / dv
            checked_cells += 1
        end
    end
    @test checked_cells > 0
end

@testset "centered nonlinear Maxwell potential update" begin
    ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.6, amplitude=0.02)
    neutral = EvolutionParams(rn=ep.rn, scalar_charge=0.0, amplitude=ep.amplitude)
    grid = mrt2013_grid(; nu=6, nv=4, U0=-5.1, V0=0.0, U1=-4.9, V1=0.06)
    state = NLState(grid)
    initialize_mrt2013_charged_outgoing_wave!(state, grid, ep;
                                              f0=mrt2013_degenerate_horizon_f0(neutral))
    evolve_nonlinear!(state, grid, ep; iterations=12, subtract_rn_background=true)

    i, j = 1, 1
    du = grid.u[i + 1] - grid.u[i]
    dv = grid.v[j + 1] - grid.v[j]
    r = NonlinearEMKGDoubleNull.corner_average(state.r[i, j], state.r[i + 1, j],
                                               state.r[i, j + 1], state.r[i + 1, j + 1])
    f = exp(NonlinearEMKGDoubleNull.corner_average(state.logf[i, j],
                                                   state.logf[i + 1, j],
                                                   state.logf[i, j + 1],
                                                   state.logf[i + 1, j + 1]))
    q = NonlinearEMKGDoubleNull.corner_average(state.Q[i, j], state.Q[i + 1, j],
                                               state.Q[i, j + 1], state.Q[i + 1, j + 1])
    phi_re = NonlinearEMKGDoubleNull.corner_average(state.phi_re[i, j],
                                                    state.phi_re[i + 1, j],
                                                    state.phi_re[i, j + 1],
                                                    state.phi_re[i + 1, j + 1])
    phi_im = NonlinearEMKGDoubleNull.corner_average(state.phi_im[i, j],
                                                    state.phi_im[i + 1, j],
                                                    state.phi_im[i, j + 1],
                                                    state.phi_im[i + 1, j + 1])
    phiu_re = NonlinearEMKGDoubleNull.corner_du(state.phi_re[i, j],
                                                state.phi_re[i + 1, j],
                                                state.phi_re[i, j + 1],
                                                state.phi_re[i + 1, j + 1], du)
    phiv_re = NonlinearEMKGDoubleNull.corner_dv(state.phi_re[i, j],
                                                state.phi_re[i + 1, j],
                                                state.phi_re[i, j + 1],
                                                state.phi_re[i + 1, j + 1], dv)
    phiu_im = NonlinearEMKGDoubleNull.corner_du(state.phi_im[i, j],
                                                state.phi_im[i + 1, j],
                                                state.phi_im[i, j + 1],
                                                state.phi_im[i + 1, j + 1], du)
    phiv_im = NonlinearEMKGDoubleNull.corner_dv(state.phi_im[i, j],
                                                state.phi_im[i + 1, j],
                                                state.phi_im[i, j + 1],
                                                state.phi_im[i + 1, j + 1], dv)
    Au = NonlinearEMKGDoubleNull.corner_average(state.Au[i, j], state.Au[i + 1, j],
                                                state.Au[i, j + 1], state.Au[i + 1, j + 1])
    Av = NonlinearEMKGDoubleNull.corner_average(state.Av[i, j], state.Av[i + 1, j],
                                                state.Av[i, j + 1], state.Av[i + 1, j + 1])
    source = stress_energy(r, f, q, phi_re, phi_im, phiu_re, phiv_re, phiu_im,
                           phiv_im, Au, Av, ep.scalar_charge)
    auv, avu, _, _, _ = maxwell_rhs(r, f, q, source)
    numerical_auv = NonlinearEMKGDoubleNull.corner_dv(state.Au[i, j], state.Au[i + 1, j],
                                                       state.Au[i, j + 1],
                                                       state.Au[i + 1, j + 1], dv)
    numerical_avu = NonlinearEMKGDoubleNull.corner_du(state.Av[i, j], state.Av[i + 1, j],
                                                       state.Av[i, j + 1],
                                                       state.Av[i + 1, j + 1], du)
    @test numerical_auv ≈ auv atol=1.0e-12
    @test numerical_avu ≈ avu atol=1.0e-12
end

@testset "GP2026 single-pulse super-extremal initial data" begin
    q0 = 1.0033218
    ep = EvolutionParams(rn=RNParams(1.0, q0), scalar_charge=0.6 / q0,
                         amplitude=0.01, omega=1.0)
    grid = gp2026_grid(; nu=27, nv=201, U0=-1.0, V0=0.0, U1=1.6, V1=20.0)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)

    @test state.r[:, 1] ≈ 1 .- grid.u ./ 2
    @test state.r[end, 1] ≈ 0.2
    @test state.r[1, :] ≈ 1.5 .+ grid.v ./ 2
    @test all(iszero, state.phi_re[:, 1])
    @test all(iszero, state.phi_im[:, 1])
    @test all(state.Q[:, 1] .== q0)

    peak = findfirst(==(10.0), grid.v)
    physical_p = hypot(state.phi_re[1, peak], state.phi_im[1, peak]) / sqrt(32 * pi)
    @test physical_p ≈ ep.amplitude
    @test state.Q[1, end] > q0
    @test all(iszero, state.Au[:, 1])
    @test all(iszero, state.Av[1, :])
    @test maximum(abs, state.Au[1, :]) > 0
    @test maximum(abs, state.Av[:, 1]) > 0

    r0 = state.r[1, 1]
    ru0 = gp2026_extremal_gauge_ru(-1.0)
    rv0 = gp2026_extremal_gauge_rv(-1.0, 0.0)
    m0 = renormalized_hawking_mass(r0, exp(state.logf[1, 1]), ru0, rv0, q0)
    @test m0 ≈ 1.0
    @test maximum(abs.(state.logf[:, 1] .- state.logf[1, 1])) < 1.0e-14

    base_f_at_end = exp(state.logf[1, 1]) *
                    gp2026_extremal_gauge_rv(-1.0, last(grid.v)) / rv0
    @test exp(state.logf[1, end]) > base_f_at_end
    _, _, _, _, q_v_residual =
        charged_charge_flux_v_profile(adaptive_state_from_rectangular(state, grid), ep;
                                      target_u=-1.0, reduced_scalar=true)
    @test maximum(abs, q_v_residual) < 3.0e-6

    ef_state = NLState(grid)
    initialize_gp2026_single_pulse!(
        ef_state, grid, ep; pulse_leg_gauge=:ef_affine,
    )
    @test ef_state.r[1, end] < state.r[1, end]
    ef_rv0 = gp2026_extremal_gauge_rv(
        -1.0, 0.0; pulse_leg_gauge=:ef_affine,
    )
    ef_mass = renormalized_hawking_mass(
        ef_state.r[1, 1], exp(ef_state.logf[1, 1]), ru0, ef_rv0, q0,
    )
    @test ef_mass ≈ 1.0
    @test_throws ArgumentError initialize_gp2026_single_pulse!(
        NLState(grid), grid, ep; pulse_leg_gauge=:invalid,
    )
end

@testset "GP2026 short charged constraint evolution" begin
    q0 = 1.0033218
    ep = EvolutionParams(rn=RNParams(1.0, q0), scalar_charge=0.6 / q0,
                         amplitude=0.01, omega=1.0)
    grid = gp2026_grid(; nu=11, nv=101, U0=-1.0, V0=0.0, U1=-0.9, V1=20.0)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations=10, reduced_scalar=true)
    adaptive = adaptive_state_from_rectangular(state, grid)
    _, _, _, _, q_u_residual =
        charged_charge_flux_u_profile(adaptive, ep; target_v=10.0, reduced_scalar=true)
    _, _, _, _, q_v_residual =
        charged_charge_flux_v_profile(adaptive, ep; target_u=-0.95, reduced_scalar=true)
    @test maximum(abs, q_u_residual) < 1.2e-5
    @test maximum(abs, q_v_residual) < 1.2e-5
end

@testset "GP2026 centered cell residuals" begin
    q0 = 1.0033218
    ep = EvolutionParams(rn=RNParams(1.0, q0), scalar_charge=0.6 / q0,
                         amplitude=0.01, omega=1.0)
    grid = gp2026_grid(; nu=11, nv=101, U0=-1.0, V0=0.0,
                       U1=-0.9, V1=20.0)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations=10, reduced_scalar=true,
                      hyperbolic_charge=true)
    summary = cell_equation_residual_summary(
        state, grid, ep; reduced_scalar=true, hyperbolic_charge=true,
    )
    @test summary.cells == (length(grid.u) - 1) * (length(grid.v) - 1)
    @test summary.max_abs_r_uv < 1.0e-10
    @test summary.max_abs_logf_uv < 1.0e-10
    @test summary.max_abs_psi_re_uv < 1.0e-11
    @test summary.max_abs_psi_im_uv < 1.0e-11
    @test summary.max_abs_q_uv < 1.0e-10
    @test summary.max_abs_au_v < 1.0e-12
    @test summary.max_abs_av_u < 1.0e-12
    @test summary.max_abs_quasilorenz < 1.0e-12
    @test summary.max_abs_faraday < 1.0e-12
    @test summary.max_abs_q_u_constraint < 2.0e-5
    @test summary.max_abs_q_v_constraint < 2.0e-5
    @test summary.max_abs_logf_gp_literal > 1.0e-2
    @test summary.max_abs_logf_coulomb2 > 1.0e-2

    newton_state = NLState(grid)
    initialize_gp2026_single_pulse!(newton_state, grid, ep)
    evolve_nonlinear!(
        newton_state, grid, ep;
        iterations=12, reduced_scalar=true, hyperbolic_charge=true,
        cell_solver=:newton_direct,
    )
    newton_summary = cell_equation_residual_summary(
        newton_state, grid, ep;
        reduced_scalar=true, hyperbolic_charge=true,
        cell_solver=:newton_direct,
    )
    @test newton_summary.max_abs_cell_residual < 1.0e-10
    @test newton_summary.max_abs_r_uv < 1.0e-7
    @test newton_summary.max_abs_f_uv < 1.0e-6
    @test newton_summary.max_abs_psi_re_uv < 1.0e-7
    @test newton_summary.max_abs_psi_im_uv < 1.0e-7
    @test newton_summary.max_abs_q_uv < 1.0e-7
    @test newton_summary.max_abs_au_v < 1.0e-8
    @test newton_summary.max_abs_av_u < 1.0e-8
end

@testset "GP2026 U-step evolution" begin
    q0 = 1.0033218
    ep = EvolutionParams(rn=RNParams(1.0, q0), scalar_charge=0.6 / q0,
                         amplitude=0.01, omega=1.0)
    grid = gp2026_grid(; nu=2, nv=101, U0=-1.0, V0=0.0, U1=-0.99, V1=20.0)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)
    initial = row_from_rectangular(state, grid, 1)

    boundary = gp2026_na_boundary_point(-0.9, ep)
    @test boundary.r ≈ 1.45
    @test boundary.logf ≈ log(gp2026_fcorner_code(ep))
    @test iszero(boundary.phi_re)
    @test iszero(boundary.phi_im)
    @test iszero(boundary.Au)
    @test boundary.Q == q0

    C = 0.1
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C, iterations=10)
    expected_first_du = min(0.1, 2C / exp(maximum(initial.logf)))
    @test evolved.rows[2].u - evolved.rows[1].u ≈ expected_first_du
    @test last(evolved.rows).u ≈ -0.9
    @test all(row -> all(isfinite, row.r) && all(isfinite, row.logf) &&
                     all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
                     all(isfinite, row.Au) && all(isfinite, row.Av) &&
                     all(isfinite, row.Q) && all(>(0), row.r), evolved.rows)
    transposed = adaptive_state_from_u_rows(evolved)
    @test length(transposed.slices) == length(grid.v)
    @test last(transposed.slices).u == [row.u for row in evolved.rows]

    paper_step = evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C,
                                          iterations=10, step_control=:outer)
    @test last(paper_step.rows).u ≈ -0.9
    substepped_paper =
        evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C, iterations=10,
                                 step_control=:outer, substep_control=:local)
    @test last(substepped_paper.rows).u ≈ -0.9
    @test all(row -> all(isfinite, row.r) && all(isfinite, row.logf) &&
                     all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
                     all(isfinite, row.Au) && all(isfinite, row.Av) &&
                     all(isfinite, row.Q) && all(>(0), row.r),
              substepped_paper.rows)

    geometric_step = evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C,
                                              iterations=10, step_control=:geometric)
    @test last(geometric_step.rows).u ≈ -0.9
    eta_step = evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C,
                                        iterations=10, step_control=:eta,
                                        max_delta_eta=0.01)
    @test last(eta_step.rows).u ≈ -0.9
    eta_du_info = gp2026_row_step_du(eta_step.rows, C, :eta; max_delta_eta=0.01)
    @test isfinite(eta_du_info.eta_du)
    @test eta_du_info.eta_du > 0
    @test eta_du_info.selected == min(eta_du_info.max_row_du, eta_du_info.eta_du)
    @test buffered_flag_intervals([false, true, false, false, true, false];
                                  buffer_points=1) == [1:6]
    @test buffered_flag_intervals([false, true, false, false, true, false];
                                  buffer_points=1, cluster=:components) == [1:3, 4:6]
    patches = row_lte_patches(initial, [2:4])
    @test only(patches).first_index == 2
    @test only(patches).last_index == 4
    @test only(patches).first_v == grid.v[2]
    lte = berger_oliger_row_lte(initial, -0.95, ep;
                                iterations=10, atol=1.0e-8, rtol=1.0e-5,
                                buffer_points=1)
    @test lte.full.u ≈ -0.95
    @test lte.refined.u ≈ -0.95
    @test length(lte.error) == length(grid.v)
    @test lte.max_error == maximum(lte.error)
    @test all(isfinite, lte.error)
    @test all(>=(0), lte.error)
    @test row_lte_error(lte.refined, lte.refined) == zeros(length(grid.v))
    @test all(interval -> first(interval) >= firstindex(grid.v) &&
                           last(interval) <= lastindex(grid.v),
              lte.intervals)
    lte_half = berger_oliger_row_lte(initial, -0.975, ep;
                                     iterations=10, atol=1.0, rtol=0.0)
    raw_lte = maximum(row_lte_error(lte.full, lte.refined;
                                    atol=1.0, rtol=0.0))
    raw_lte_half = maximum(lte_half.error)
    lte_ratio = raw_lte / raw_lte_half
    @test 6.0 < lte_ratio < 10.0
    child = berger_oliger_refine_patch(
        initial, lte, ep;
        interval=2:21, refinement_factor=2, iterations=10,
    )
    @test child.interval == 2:21
    @test length(child.fine_v) == (length(child.interval) - 1) * 2 + 1
    @test child.reintegrated
    @test child.parent_midpoint.v == grid.v
    @test child.parent_target.v == grid.v
    @test all(isfinite, child.child_midpoint.r)
    @test all(isfinite, child.child_target.r)
    @test all(isfinite, child.parent_target.r)
    @test child.max_correction >= 0
    for (offset, parent_index) in enumerate(child.interval)
        child_index = 1 + 2 * (offset - 1)
        @test child.parent_midpoint.r[parent_index] ≈
              child.child_midpoint.r[child_index]
        @test child.parent_target.r[parent_index] ≈
              child.child_target.r[child_index]
        @test child.parent_target.Q[parent_index] ≈
              child.child_target.Q[child_index]
    end
    bo_evolved = evolve_gp2026_u_adaptive(
        initial, ep;
        Umax=-0.98, C, iterations=10, max_rows=20,
        step_control=:outer, bo_amr=true,
        bo_rtol=1.0e-8, bo_refinement_factor=2,
    )
    @test last(bo_evolved.rows).u ≈ -0.98
    @test all(row -> all(isfinite, row.r) && all(isfinite, row.logf),
              bo_evolved.rows)
    @test_throws ArgumentError evolve_gp2026_u_adaptive(
        initial, ep;
        Umax=-0.98, C, bo_amr=true, backtrack=true,
    )
    backtracked_step =
        evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C,
                                 iterations=10, step_control=:outer,
                                 backtrack=true,
                                 max_realized_delta_rho=0.05,
                                 max_realized_delta_eta=0.02,
                                 max_rows=1_000)
    @test last(backtracked_step.rows).u ≈ -0.9
    for k in 2:length(backtracked_step.rows)
        change = realized_row_change_summary(backtracked_step.rows[1:k - 1],
                                             backtracked_step.rows[k])
        @test change.finite
        @test change.max_abs_rho <= 0.05 + 100eps()
        @test change.max_abs_eta <= 0.02 + 100eps()
    end

    throat = throat_row_diagnostics(last(evolved.rows))
    @test length(throat.y) == length(grid.v)
    @test all(>(0), throat.y)
    @test all(isfinite, throat.rho)
    @test length(throat.eta) == length(grid.v)
    @test length(throat.zeta) == length(grid.v)
    @test all(isfinite, throat.eta)
    @test all(isfinite, throat.zeta)
    @test all(>=(0), throat.eta)
    @test all(>(0), throat.zeta)
    @test throat.max_rho == maximum(throat.rho)
    @test throat.max_abs_delta_rho >= 0
    @test throat.max_abs_delta_eta >= 0
    @test throat.max_abs_delta_zeta >= 0
    match = throat_matching_candidate(last(evolved.rows); rho_min=0.0)
    @test !isnothing(match)
    @test match.index == firstindex(last(evolved.rows).v)
    @test match.rho == throat.rho[match.index]
    band = throat_matching_band(last(evolved.rows); rho_min=0.0)
    @test !isnothing(band)
    rho_indices = findall(>=(0.0), throat.rho)
    @test band.first_index == first(rho_indices)
    @test band.last_index == last(rho_indices)
    @test band.count == length(rho_indices)
    @test band.component_count >= 1
    lapse = rho_lapse_diagnostics(last(evolved.rows); rho_min=0.0)
    @test length(lapse.rho) == length(last(evolved.rows).v)
    @test length(lapse.rho_v) == length(last(evolved.rows).v)
    @test length(lapse.logf_rho) == length(last(evolved.rows).v)
    @test all(isfinite, lapse.logf_range)
    @test all(isfinite, lapse.logf_rho_range)
    @test lapse.throat_count == length(rho_indices)
    @test range_width(lapse.logf_range) >= 0
    @test range_width(lapse.logf_rho_range) >= 0

    synthetic_v = [0.0, 1.0, 2.0]
    synthetic_q = ones(3)
    synthetic_r = 1 .+ exp.(-[3.0, 2.0, 1.0])
    synthetic_row = NLRow(
        0.0,
        synthetic_v,
        synthetic_r,
        zeros(3),
        [1.0, 2.0, 3.0],
        [-1.0, -2.0, -3.0],
        zeros(3),
        zeros(3),
        synthetic_q,
    )
    boundary_sample = throat_boundary_sample(synthetic_row; rho_match=2.0)
    @test !isnothing(boundary_sample)
    @test boundary_sample.v ≈ 1.0
    @test boundary_sample.rho ≈ 2.0
    @test boundary_sample.y ≈ exp(-2.0)
    @test boundary_sample.phi_re ≈ 2.0
    @test boundary_sample.phi_im ≈ -2.0
    @test boundary_sample.rho_v ≈ -1.0
    @test isnothing(throat_boundary_sample(synthetic_row; rho_match=4.0))
    boundary_series = throat_boundary_series(UAdaptiveNLState([synthetic_row]);
                                             rho_match=2.0)
    @test length(boundary_series) == 1
    @test boundary_series[1].row_index == 1
    neutral_observables =
        throat_boundary_observables(boundary_sample,
                                    EvolutionParams(rn=RNParams(1.0, 1.0),
                                                    scalar_charge=0.0))
    @test neutral_observables.psi_abs ≈ sqrt(8.0)
    @test neutral_observables.rphi_gp_abs ≈ sqrt(8.0) / sqrt(32pi)
    @test neutral_observables.raw_phase ≈ -pi / 4
    @test neutral_observables.rphi_gp_abs_v ≈ sqrt(2.0) / sqrt(32pi)
    @test neutral_observables.covariant_phase_v ≈ 0.0 atol=1.0e-14
    @test neutral_observables.covariant_dv_rphi_abs ≈ sqrt(2.0) / sqrt(32pi)
    @test neutral_observables.q_v_source ≈ 0.0
    @test neutral_observables.q_v_residual ≈ 0.0
    @test neutral_observables.one_minus_absq_over_r ≈ 1 - 1 / boundary_sample.r

    throat_step = evolve_gp2026_u_adaptive(initial, ep; Umax=-0.9, C,
                                           iterations=10, step_control=:throat)
    @test last(throat_step.rows).u ≈ -0.9
end

@testset "persistent Hamade-Stewart hierarchy" begin
    ep = EvolutionParams(
        rn=RNParams(1.0, 1.0),
        scalar_charge=0.6,
        amplitude=0.01,
        omega=1.0,
    )
    grid = gp2026_grid(
        ; nu=2, nv=41, U0=-1.0, V0=0.0, U1=-0.99, V1=20.0,
    )
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)
    initial = row_from_rectangular(state, grid, 1)

    target_boundary = gp2026_na_boundary_point(-0.98, ep)
    upper = advance_u_row(
        initial, target_boundary, ep;
        iterations=12,
        reduced_scalar=true,
        hyperbolic_charge=true,
        cell_solver=:newton_direct,
    )
    midpoint_boundary =
        interpolate_parent_boundary(initial, upper, grid.v[8], -0.99)
    lower_point = row_point(initial, 8)
    upper_point = row_point(upper, 8)
    @test midpoint_boundary.u == -0.99
    @test midpoint_boundary.v == grid.v[8]
    @test midpoint_boundary.r ≈ (lower_point.r + upper_point.r) / 2
    @test midpoint_boundary.Q ≈ (lower_point.Q + upper_point.Q) / 2

    persistent_config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=100,
        max_levels=2,
        atol=1.0,
        rtol=0.0,
        buffer_points=2,
    )
    persistent = initialize_stewart_hierarchy(
        initial;
        config=persistent_config,
    )
    patch = 8:25
    NonlinearEMKGDoubleNull.rebuild_stewart_child!(
        persistent.root,
        patch,
        persistent.config,
        persistent.stats,
    )
    persistent.root.steps_since_revision = 0
    child_identity = persistent.root.child
    first_result = advance_stewart_hierarchy!(
        persistent, -0.98, ep;
        iterations=12,
        cell_solver=:newton_direct,
    )
    second_result = advance_stewart_hierarchy!(
        persistent, -0.96, ep;
        iterations=12,
        cell_solver=:newton_direct,
    )
    @test persistent.root.child === child_identity
    @test first_result.depth == 2
    @test second_result.depth == 2
    @test second_result.row.u == -0.96
    @test persistent.root.child.current.u == -0.96
    @test persistent.stats.level_steps == [2, 8]
    @test persistent.stats.injections == 2
    @test persistent.stats.suffix_reintegrations == 2
    @test validate_stewart_hierarchy(persistent.root)
    for (offset, parent_index) in enumerate(patch)
        child_index = 1 + 4 * (offset - 1)
        @test persistent.root.current.r[parent_index] ≈
              persistent.root.child.current.r[child_index]
        @test persistent.root.current.Q[parent_index] ≈
              persistent.root.child.current.Q[child_index]
    end

    destruction_config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=1,
        max_levels=2,
        atol=1.0,
        rtol=0.0,
    )
    destruction = initialize_stewart_hierarchy(
        initial;
        config=destruction_config,
    )
    NonlinearEMKGDoubleNull.rebuild_stewart_child!(
        destruction.root,
        patch,
        destruction.config,
        destruction.stats,
    )
    destruction.root.steps_since_revision = 1
    destruction_result = advance_stewart_hierarchy!(
        destruction, -0.98, ep;
        iterations=12,
        cell_solver=:newton_direct,
    )
    @test destruction_result.depth == 1
    @test isnothing(destruction.root.child)
    @test destruction.stats.child_destructions == 1

    containment_config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=4,
        max_levels=3,
    )
    containment = initialize_stewart_hierarchy(
        initial;
        config=containment_config,
    )
    NonlinearEMKGDoubleNull.rebuild_stewart_child!(
        containment.root,
        5:30,
        containment.config,
        containment.stats,
    )
    NonlinearEMKGDoubleNull.rebuild_stewart_child!(
        containment.root.child,
        30:60,
        containment.config,
        containment.stats,
    )
    old_grandchild_range = extrema(containment.root.child.child.current.v)
    NonlinearEMKGDoubleNull.rebuild_stewart_child!(
        containment.root,
        20:25,
        containment.config,
        containment.stats,
    )
    @test stewart_hierarchy_depth(containment.root) == 3
    @test first(containment.root.child.current.v) <=
          first(old_grandchild_range)
    @test last(containment.root.child.current.v) >=
          last(old_grandchild_range)
    @test extrema(containment.root.child.child.current.v) ==
          old_grandchild_range
    @test validate_stewart_hierarchy(containment.root)

    recursive_config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=1,
        max_levels=3,
        atol=1.0e-14,
        rtol=1.0e-12,
        buffer_points=2,
    )
    recursive = initialize_stewart_hierarchy(
        initial;
        config=recursive_config,
    )
    recursive_result = advance_stewart_hierarchy!(
        recursive, -0.98, ep;
        iterations=12,
        cell_solver=:newton_direct,
    )
    @test recursive_result.depth == 3
    @test recursive.stats.level_steps == [1, 4, 16]
    @test recursive.stats.max_level_reached == 2
    @test recursive.stats.injections >= 5
    @test recursive.root.current.u ==
          recursive.root.child.current.u ==
          recursive.root.child.child.current.u
    @test validate_stewart_hierarchy(recursive.root)

    vacuum = EvolutionParams(
        rn=RNParams(1.0, 1.0),
        scalar_charge=0.0,
        amplitude=0.0,
        omega=0.0,
    )
    vacuum_errors = Float64[]
    vacuum_logf_errors = Float64[]
    for (du, dv) in ((0.04, 0.2), (0.02, 0.1))
        U0, U1, V1 = -0.2, 0.12, 4.0
        root_steps = Int(round((U1 - U0) / du))
        vacuum_grid = gp2026_grid(
            ; nu=2,
            nv=Int(round(V1 / dv)) + 1,
            U0,
            V0=0.0,
            U1=U0 + du,
            V1,
        )
        vacuum_state = NLState(vacuum_grid)
        initialize_gp2026_exact_extremal_rn!(
            vacuum_state, vacuum_grid, vacuum,
        )
        vacuum_initial = row_from_rectangular(vacuum_state, vacuum_grid, 1)
        vacuum_config = StewartAMRConfig(
            refinement_factor=4,
            revision_interval=4,
            max_levels=2,
            atol=1.0e-14,
            rtol=1.0e-10,
        )
        vacuum_hierarchy = initialize_stewart_hierarchy(
            vacuum_initial;
            config=vacuum_config,
        )
        vacuum_rows = NLRow[vacuum_initial]
        for step in 1:root_steps
            target_u = step == root_steps ? U1 : U0 + step * du
            result = advance_stewart_hierarchy!(
                vacuum_hierarchy, target_u, vacuum;
                U0,
                pulse_leg_gauge=:ef_affine,
                iterations=15,
                cell_solver=:newton_direct,
            )
            push!(vacuum_rows, result.row)
        end
        push!(
            vacuum_errors,
            maximum(
                abs(
                    row.r[j] -
                    gp2026_exact_extremal_rn_radius(
                        row.u, row.v[j], vacuum.rn,
                    ),
                )
                for row in vacuum_rows for j in eachindex(row.v)
            ),
        )
        push!(
            vacuum_logf_errors,
            maximum(
                abs(
                    row.logf[j] -
                    log(
                        gp2026_exact_extremal_rn_fcode(
                            row.u, row.v[j], vacuum.rn,
                        ),
                    ),
                )
                for row in vacuum_rows for j in eachindex(row.v)
            ),
        )
        @test maximum(
            abs(row.Q[j] - 1.0)
            for row in vacuum_rows for j in eachindex(row.Q)
        ) < 1.0e-12
    end
    @test 3.7 < vacuum_errors[1] / vacuum_errors[2] < 4.3
    @test 3.7 < vacuum_logf_errors[1] / vacuum_logf_errors[2] < 4.3
end

@testset "MRT ingoing initial leg normalization" begin
    ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=0.0)
    grid = mrt2013_grid(; nu=16, nv=48, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=20.0)
    state = NLState(grid)
    initialize_mrt2013_uncharged_ingoing!(state, grid, ep)

    exact_logf = [log(mrt2013_metric_f(grid.u[1], V, ep.rn)) for V in grid.v]
    @test maximum(abs.(state.logf[1, :] .- exact_logf)) < 1.0e-10
    @test mrt2013_areal_radius(0.0, 120.0, ep.rn) == 1.0
    @test mrt2013_metric_f(0.0, 120.0, ep.rn) == 2.0
    @test mrt2013_areal_radius(0.2, 0.0, ep.rn) ≈ 0.8
    @test mrt2013_metric_f(0.2, 0.0, ep.rn) ≈ 2.0
    @test 0.8 < mrt2013_areal_radius(0.2, 1.0, ep.rn) < 1.0
end

@testset "MRT outgoing-wave initial data" begin
    vacuum = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=0.0)
    grid = mrt2013_grid(; nu=32, nv=16, U0=-5.1, V0=0.0, U1=0.0, V1=1.0)
    state = NLState(grid)
    initialize_mrt2013_outgoing_wave!(state, grid, vacuum)
    exact_logf = [log(mrt2013_metric_f(U, 0.0, vacuum.rn)) for U in grid.u]
    @test maximum(abs.(state.logf[:, 1] .- exact_logf)) < 5.0e-12
    @test all(iszero, state.phi_re)

    perturbed = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=1.0e-3)
    initialize_mrt2013_outgoing_wave!(state, grid, perturbed)
    @test maximum(abs, state.phi_re[:, 1]) > 0
    @test all(iszero, state.phi_re[1, :])
    @test state.logf[end, 1] < log(2.0)

    f0 = mrt2013_degenerate_horizon_f0(perturbed)
    mass = mrt2013_initial_bondi_mass(f0)
    @test (f0 - 2) / perturbed.amplitude^2 ≈ 0.740 atol=2.0e-3
    @test (mass - 1) / perturbed.amplitude^2 ≈ 0.789 atol=2.0e-3

    aligned_coarse = mrt2013_grid(; nu=531, nv=2, U0=-5.1, V0=0.0, U1=0.2, V1=0.02)
    aligned_fine = mrt2013_grid(; nu=1061, nv=2, U0=-5.1, V0=0.0, U1=0.2, V1=0.02)
    coarse_state = NLState(aligned_coarse)
    fine_state = NLState(aligned_fine)
    initialize_mrt2013_outgoing_wave!(coarse_state, aligned_coarse, perturbed; f0)
    initialize_mrt2013_outgoing_wave!(fine_state, aligned_fine, perturbed; f0)
    coarse_zero = findfirst(iszero, aligned_coarse.u)
    fine_zero = findfirst(iszero, aligned_fine.u)
    coarse_rv = mrt2013_initial_rv_profile(coarse_state, aligned_coarse)[coarse_zero]
    fine_rv = mrt2013_initial_rv_profile(fine_state, aligned_fine)[fine_zero]
    @test abs(coarse_rv / fine_rv) ≈ 4 rtol=0.02

    charged = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.6,
                              amplitude=1.0e-3)
    charged_state = NLState(aligned_coarse)
    initialize_mrt2013_charged_outgoing_wave!(charged_state, aligned_coarse, charged;
                                              f0=mrt2013_degenerate_horizon_f0(perturbed))
    @test all(charged_state.Q .== charged.rn.Q0)
    @test all(iszero, charged_state.Au[:, 1])
    @test all(iszero, charged_state.Av[1, :])
    @test maximum(abs, charged_state.Au[1, :]) > 0
    @test maximum(abs, charged_state.Av[:, 1]) > 0
end

@testset "adaptive MRT slice utilities" begin
    u = [-2.0, -1.0, 1.0]
    slice = NLSlice(
        3.0,
        u,
        1 .+ u,
        2 .- u,
        3 .* u,
        -2 .* u,
        4 .+ 2 .* u,
        5 .- 3 .* u,
        fill(0.7, length(u)),
    )

    refined_u = refine_u_grid(u, [true, false])
    @test refined_u == [-2.0, -1.5, -1.0, 1.0]

    refined = interpolate_slice(slice, refined_u)
    @test refined.r ≈ 1 .+ refined_u
    @test refined.logf ≈ 2 .- refined_u
    @test refined.phi_re ≈ 3 .* refined_u
    @test refined.Q ≈ fill(0.7, length(refined_u))

    @test spacing_refinement_flags(u; max_du=1.5) == [false, true]
    @test variation_refinement_flags(slice; max_dr=1.5) == [false, true]
    splitting = PointSplittingConfig(; band_width=0.1, max_relative_r=0.4,
                                     max_relative_f=Inf)
    @test point_splitting_flags(slice, splitting) == [true, true]

    scalar_u = [-2.0, -1.0, 0.0, 1.0]
    scalar_fields = fill(1.0, length(scalar_u))
    scalar_slice = NLSlice(0.0, scalar_u, scalar_fields, zero.(scalar_u),
                           [0.0, 0.1, 0.1, 0.1], [0.0, 0.0, 0.2, 0.2],
                           zero.(scalar_u), zero.(scalar_u), scalar_fields)
    scalar_splitting = PointSplittingConfig(; max_dphi=0.15)
    @test point_splitting_flags(scalar_slice, scalar_splitting) ==
          [false, true, false]

    gradient_slice = NLSlice(0.0, scalar_u, scalar_fields, zero.(scalar_u),
                             [0.0, 0.0, 0.0, 1.0], zero.(scalar_u),
                             zero.(scalar_u), zero.(scalar_u), scalar_fields)
    gradient_splitting = PointSplittingConfig(; max_dphiu=0.5)
    @test point_splitting_flags(gradient_slice, gradient_splitting) ==
          [false, true, true]

    chop_u = collect(0.0:5.0)
    previous = NLSlice(0.0, chop_u, fill(2.0, 6), zeros(6), zeros(6), zeros(6),
                       zeros(6), zeros(6), ones(6))
    current = NLSlice(0.1, chop_u, 2 .+ [0.2, 0.1, 0.05, -0.1, -0.2, -0.3],
                      zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), ones(6))
    chopping = HorizonChoppingConfig(; band_width=0.1, interior_buffer_cells=1)
    @test adaptive_outgoing_expansion(previous, current) ≈
          [1.5, 0.75, -0.25, -1.5, -2.5]
    chopped = chop_inside_apparent_horizon(previous, current, chopping)
    @test length(chopped.u) == 5
    horizon_v, horizon_u, horizon_rphi =
        gp2026_horizon_rphi_series(AdaptiveNLState([previous, chopped]))
    @test length(horizon_v) == length(horizon_u) == length(horizon_rphi) == 1
    @test only(horizon_rphi) == 0
    width_chopping = HorizonChoppingConfig(; band_width=0.1, interior_buffer_cells=0,
                                           interior_buffer_width=0.5)
    @test length(chop_inside_apparent_horizon(previous, current, width_chopping).u) == 5
    horizon_refinement = HorizonRefinementConfig(; band_width=0.1, max_du=0.6,
                                                  exterior_cells=1, interior_cells=1)
    @test horizon_refinement_flags(previous, current, horizon_refinement) ==
          [false, true, true, true, false]
    @test length(refine_near_apparent_horizon(previous, current, horizon_refinement).u) == 9

    row_v = [0.0, 0.5, 1.0]
    row_r = [2.0, 2.25, 2.0]
    row = NLRow(-0.5, row_v, row_r, zeros(3), zeros(3), zeros(3),
                zeros(3), zeros(3), ones(3))
    @test row_outgoing_expansion(row) ≈ [0.5, 0.0, -0.5]
    minimum_sample = row_expansion_minimum(row; row_index=7)
    @test minimum_sample.row_index == 7
    @test minimum_sample.v ≈ 1.0
    @test minimum_sample.rv ≈ -0.5
    crossing = row_apparent_horizon_crossing(row; row_index=7)
    @test crossing.row_index == 7
    @test crossing.u ≈ -0.5
    @test crossing.v ≈ 0.5
    @test crossing.rv ≈ 0.0
    diagnostic = vtrap_diagnostic([row]; missing_status=:max_rows)
    @test diagnostic.status == :trapped
    @test diagnostic.trapped
    @test diagnostic.trap.v ≈ 0.5

    untrapped_row = NLRow(-0.4, row_v, [2.0, 2.2, 2.5], zeros(3), zeros(3),
                          zeros(3), zeros(3), zeros(3), ones(3))
    untrapped = vtrap_diagnostic([untrapped_row]; missing_status=:precision_stalled)
    @test untrapped.status == :precision_stalled
    @test !untrapped.trapped
    @test isnothing(untrapped.trap)
    @test untrapped.closest.rv > 0

    cubic_u = [-2.0, -1.0, 0.0, 1.0]
    cubic = cubic_u .^ 3
    cubic_slice = NLSlice(0.0, cubic_u, cubic, cubic, cubic, cubic, cubic, cubic,
                          fill(1.0, length(cubic_u)))
    split_cubic = refine_slice(cubic_slice, [true, false, true])
    @test split_cubic.r ≈ split_cubic.u .^ 3

    constrained_ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=0.0)
    constrained_u = [0.0, 0.5, 1.0]
    slope = 0.2
    constrained_r = 2 .- constrained_u
    constrained_phi = slope .* constrained_u
    constrained_logf = log(2.0) .-
                       slope^2 / 4 .* (2 .* constrained_u .- constrained_u .^ 2 ./ 2)
    constrained_slice = NLSlice(0.0, constrained_u, constrained_r, constrained_logf,
                                constrained_phi, zeros(3), zeros(3), zeros(3), ones(3))
    split_constrained = refine_slice_constrained(constrained_slice, [true, true],
                                                 constrained_ep)
    expected_logf = log(2.0) .-
                    slope^2 / 4 .* (2 .* split_constrained.u .-
                                     split_constrained.u .^ 2 ./ 2)
    @test split_constrained.logf ≈ expected_logf
    @test first(split_constrained.logf) == first(constrained_slice.logf)

    grid = Grid(u, [0.0, 1.0])
    state = NLState(grid)
    state.r[:, 2] .= slice.r
    state.logf[:, 2] .= slice.logf
    state.phi_re[:, 2] .= slice.phi_re
    state.phi_im[:, 2] .= slice.phi_im
    state.Au[:, 2] .= slice.Au
    state.Av[:, 2] .= slice.Av
    state.Q[:, 2] .= slice.Q
    from_rect = slice_from_rectangular(state, grid, 2)
    @test from_rect.v == 1.0
    @test from_rect.u == u
    @test from_rect.phi_im == slice.phi_im
    @test slice_point(from_rect, 2).r == from_rect.r[2]
end

@testset "MRT renormalized Hawking mass" begin
    r = 2.0
    f = 2.0
    ru = -1.0
    rv = 0.5 * (1 - 1 / r)^2
    @test renormalized_hawking_mass(r, f, ru, rv, 1.0) ≈ 1.0

    lower = NLSlice(0.0, [-2.0, -1.0, 0.0], [3.0, 2.0, 1.0],
                    fill(log(2.0), 3), zeros(3), zeros(3), zeros(3), zeros(3),
                    ones(3))
    upper = NLSlice(0.1, [-2.0, -1.0, 0.0], [3.01, 2.01, 1.01],
                    fill(log(2.0), 3), zeros(3), zeros(3), zeros(3), zeros(3),
                    ones(3))
    refined_upper = refine_slice(upper, [true, true])
    _, _, original_mass = renormalized_hawking_mass_profile(lower, upper)
    _, _, refined_mass = renormalized_hawking_mass_profile(lower, refined_upper)
    @test refined_mass ≈ original_mass

    p = RNParams(1.0, 1.0)
    u = [-5.1, -2.0, -0.2]
    exact_lower = NLSlice(0.0, u,
                          [mrt2013_areal_radius(U, 0.0, p) for U in u],
                          [log(mrt2013_metric_f(U, 0.0, p)) for U in u],
                          zeros(3), zeros(3), zeros(3), zeros(3), ones(3))
    exact_upper = NLSlice(0.2, u,
                          [mrt2013_areal_radius(U, 0.2, p) for U in u],
                          [log(mrt2013_metric_f(U, 0.2, p)) for U in u],
                          zeros(3), zeros(3), zeros(3), zeros(3), ones(3))
    _, _, corrected_mass =
        renormalized_hawking_mass_profile(exact_lower, exact_upper; rn_background=p)
    @test corrected_mass ≈ ones(2)
    _, _, mass_u, expected_mass_u, residual =
        uncharged_mass_flux_u_profile(exact_lower, exact_upper; rn_background=p)
    @test mass_u ≈ zeros(1)
    @test expected_mass_u == zeros(1)
    @test residual ≈ zeros(1)
    _, _, geometric_mass, flux_mass, balance_error =
        uncharged_flux_integrated_mass_profile(exact_lower, exact_upper; rn_background=p)
    @test geometric_mass ≈ ones(2)
    @test flux_mass ≈ ones(2)
    @test balance_error ≈ zeros(2)
    _, _, q_u, expected_q_u, q_residual =
        charged_charge_flux_u_profile(exact_lower, exact_upper, 0.6)
    @test q_u ≈ zeros(1)
    @test expected_q_u ≈ zeros(1)
    @test q_residual ≈ zeros(1)
    _, _, geometric_q, flux_q, q_balance_error =
        charged_flux_integrated_charge_profile(exact_lower, exact_upper, 0.6)
    @test geometric_q ≈ ones(2)
    @test flux_q ≈ ones(2)
    @test q_balance_error ≈ zeros(2)
    _, _, q_v, expected_q_v, q_v_residual =
        charged_charge_flux_v_profile(AdaptiveNLState([exact_lower, exact_upper]), 0.6;
                                      target_u=-1.0)
    @test q_v ≈ zeros(1)
    @test expected_q_v ≈ zeros(1)
    @test q_v_residual ≈ zeros(1)
    _, _, charged_geometric_mass, charged_flux_mass, charged_balance_error =
        charged_flux_integrated_mass_profile(exact_lower, exact_upper, 0.6; rn_background=p)
    @test charged_geometric_mass ≈ ones(2)
    @test charged_flux_mass ≈ ones(2)
    @test charged_balance_error ≈ zeros(2)

    @test apparent_horizon_location([-1.0, 0.0, 1.0], [0.5, 0.25, -0.25]) ≈ 0.5
    @test isnothing(apparent_horizon_location([-1.0, 0.0], [0.5, 0.25]))
end

@testset "adaptive MRT one-step advance" begin
    ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=0.0)
    u = collect(range(-5.1, -1.0e-3; length=8))
    v0 = 0.0
    v1 = 0.2

    exact_r0 = [mrt2013_areal_radius(U, v0, ep.rn) for U in u]
    exact_logf0 = [log(mrt2013_metric_f(U, v0, ep.rn)) for U in u]
    previous = NLSlice(v0, u, exact_r0, exact_logf0, zero.(u), zero.(u),
                       zero.(u), zero.(u), fill(ep.rn.Q0, length(u)))

    northwest = NLPoint(
        first(u),
        v1,
        mrt2013_areal_radius(first(u), v1, ep.rn),
        log(mrt2013_metric_f(first(u), v1, ep.rn)),
        0.0,
        0.0,
        0.0,
        0.0,
        ep.rn.Q0,
    )
    next = advance_adaptive_slice(previous, northwest, ep;
                                  iterations=20, subtract_rn_background=true)

    exact_r1 = [mrt2013_areal_radius(U, v1, ep.rn) for U in u]
    exact_logf1 = [log(mrt2013_metric_f(U, v1, ep.rn)) for U in u]
    @test maximum(abs.(next.r .- exact_r1)) < 1.0e-10
    @test maximum(abs.(next.logf .- exact_logf1)) < 1.0e-10
    @test all(next.Q .≈ ep.rn.Q0)
end

@testset "adaptive MRT driver matches rectangular evolution" begin
    ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=1.0e-5)
    grid = mrt2013_grid(; nu=10, nv=16, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=1.0)
    rectangular = NLState(grid)
    initialize_mrt2013_uncharged_ingoing!(rectangular, grid, ep)

    initial = slice_from_rectangular(rectangular, grid, 1)
    west_boundary = west_boundary_from_rectangular(rectangular, grid)
    adaptive = evolve_adaptive(initial, west_boundary, ep;
                               iterations=5, subtract_rn_background=true)

    evolve_nonlinear!(rectangular, grid, ep; iterations=5, subtract_rn_background=true)
    @test length(adaptive.slices) == length(grid.v)
    for j in eachindex(grid.v)
        slice = adaptive.slices[j]
        @test slice.u == grid.u
        @test slice.r ≈ rectangular.r[:, j]
        @test slice.logf ≈ rectangular.logf[:, j]
        @test slice.phi_re ≈ rectangular.phi_re[:, j]
        @test slice.Q ≈ rectangular.Q[:, j]
    end
end

@testset "adaptive MRT point splitting policy" begin
    ep = EvolutionParams(rn=RNParams(1.0, 1.0), scalar_charge=0.0, amplitude=0.0)
    grid = mrt2013_grid(; nu=6, nv=4, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=0.3)
    rectangular = NLState(grid)
    initialize_mrt2013_uncharged_ingoing!(rectangular, grid, ep)

    splitting = PointSplittingConfig(; band_width=0.1, max_relative_r=1.0e-12,
                                     max_relative_f=Inf)
    adaptive = evolve_adaptive(slice_from_rectangular(rectangular, grid, 1),
                               west_boundary_from_rectangular(rectangular, grid), ep;
                               iterations=5, subtract_rn_background=true,
                               point_splitting=splitting)

    @test length(adaptive.slices[1].u) == 6
    @test length(adaptive.slices[2].u) == 11
    @test length(adaptive.slices[3].u) == 21
    @test length(adaptive.slices[4].u) == 21
    @test all(isfinite, adaptive.slices[end].r)
    @test all(isfinite, adaptive.slices[end].logf)
end

@testset "stress energy sources" begin
    src = stress_energy(2.0, 3.0, 0.4, 0.1, -0.2, 0.3, -0.1, 0.05, 0.2, 0.7, -0.4, 0.6)
    @test isfinite(src.Tuu)
    @test isfinite(src.Tvv)
    @test isfinite(src.Tuv)
    @test isfinite(src.Tthth)
    @test isfinite(src.Ju)
    @test isfinite(src.Jv)
    @test isfinite(src.scalar_logf_source)
    @test src.Tuu >= 0
    @test src.Tvv >= 0

    neutral = stress_energy(2.0, 3.0, 0.4, 0.1, -0.2, 0.3, -0.1, 0.05, 0.2, 0.0, 0.0, 0.0)
    @test neutral.Ju == 0
    @test neutral.Jv == 0
    @test neutral.Tuv ≈ neutral.alpha^2 / 3.0
    @test neutral.alpha ≈ -0.4 * 3.0 / (2 * 2.0^2)
    @test neutral.scalar_logf_source ≈ (0.3 * -0.1 + 0.05 * 0.2) / 2

    auv, avu, _, _, _ = maxwell_rhs(2.0, 3.0, 0.4, neutral)
    @test avu - auv ≈ -0.4 * 3.0 / (2 * 2.0^2)
    @test auv + avu ≈ 0.0

    _, _, _, qu, qv = maxwell_rhs(2.0, 3.0, 0.4, src)
    @test qu ≈ 2.0^2 * src.Ju / 8
    @test qv ≈ -2.0^2 * src.Jv / 8

    r, ru, rv = 2.0, -0.3, 0.2
    psi_re, psi_im = r * 0.1, r * -0.2
    psi_re_u = r * 0.3 + ru * 0.1
    psi_re_v = r * -0.1 + rv * 0.1
    psi_im_u = r * 0.05 + ru * -0.2
    psi_im_v = r * 0.2 + rv * -0.2
    from_reduced = stress_energy_reduced_scalar(
        r, 3.0, 0.4, ru, rv, psi_re, psi_im, psi_re_u, psi_re_v,
        psi_im_u, psi_im_v, 0.7, -0.4, 0.6,
    )
    @test from_reduced.Tuu ≈ src.Tuu
    @test from_reduced.Tvv ≈ src.Tvv
    @test from_reduced.Ju ≈ src.Ju
    @test from_reduced.Jv ≈ src.Jv
    @test from_reduced.scalar_logf_source ≈ src.scalar_logf_source
end
