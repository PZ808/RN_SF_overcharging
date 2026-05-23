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
    @test length(chop_inside_apparent_horizon(previous, current, chopping).u) == 5

    cubic_u = [-2.0, -1.0, 0.0, 1.0]
    cubic = cubic_u .^ 3
    cubic_slice = NLSlice(0.0, cubic_u, cubic, cubic, cubic, cubic, cubic, cubic,
                          fill(1.0, length(cubic_u)))
    split_cubic = refine_slice(cubic_slice, [true, false, true])
    @test split_cubic.r ≈ split_cubic.u .^ 3

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

    auv, avu, _, _, _ = maxwell_rhs(2.0, 3.0, 0.4, neutral)
    @test avu - auv ≈ 0.4 * 3.0 / 2.0^2
    @test auv + avu ≈ 0.0
end
