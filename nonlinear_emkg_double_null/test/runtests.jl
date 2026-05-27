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
    @test maximum(abs, q_u_residual) < 3.0e-6
    @test maximum(abs, q_v_residual) < 1.2e-5
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
    @test src.Tuu >= 0
    @test src.Tvv >= 0

    neutral = stress_energy(2.0, 3.0, 0.4, 0.1, -0.2, 0.3, -0.1, 0.05, 0.2, 0.0, 0.0, 0.0)
    @test neutral.Ju == 0
    @test neutral.Jv == 0
    @test neutral.Tuv ≈ neutral.alpha^2 / 3.0
    @test neutral.alpha ≈ -0.4 * 3.0 / (2 * 2.0^2)

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
end
