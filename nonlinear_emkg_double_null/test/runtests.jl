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
