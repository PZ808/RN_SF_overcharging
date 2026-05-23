using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

# MRT Fig. 13 uses an outgoing pulse on Sigma_1 and samples the event horizon.
# The small amplitude keeps this first check close to the fixed-extreme-RN
# scalar problem while the nonlinear/adaptive machinery is exercised.
ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = 1.0e-8,
    omega = 0.0,
    center = 15.0,
    width = 3.0,
)

grid = mrt2013_grid(; nu=80, nv=1501, U0=-5.1, V0=0.0, U1=0.0, V1=300.0)
seed = NLState(grid)
initialize_mrt2013_outgoing_wave!(seed, grid, ep)

geometry_threshold = argument_or_default(1, 0.03)
phi_threshold_over_amplitude = argument_or_default(2, 0.03)
phiu_threshold_over_amplitude = argument_or_default(3, 0.03)
splitting = PointSplittingConfig(
    band_width = 1.0,
    max_relative_r = geometry_threshold,
    max_relative_f = geometry_threshold,
    max_dphi = phi_threshold_over_amplitude * ep.amplitude,
    max_dphiu = phiu_threshold_over_amplitude * ep.amplitude,
)
state = evolve_adaptive(
    slice_from_rectangular(seed, grid, 1),
    west_boundary_from_rectangular(seed, grid),
    ep;
    iterations=20,
    subtract_rn_background=true,
    point_splitting=splitting,
)

vef = [slice.v for slice in state.slices]
phi_horizon = [slice.phi_re[end] for slice in state.slices]
slope, intercept, nfit = fit_power_law(vef, abs.(phi_horizon); xmin=150.0, xmax=300.0)

last_slice = state.slices[end]
dphi_dr = (last_slice.phi_re[end] - last_slice.phi_re[end - 1]) /
          (last_slice.r[end] - last_slice.r[end - 1])
aretakis_h = last_slice.phi_re[end] + last_slice.r[end] * dphi_dr

println("relative geometry threshold = ", geometry_threshold)
println("Delta phi / amplitude threshold = ", phi_threshold_over_amplitude)
println("Delta phi_U / amplitude threshold = ", phiu_threshold_over_amplitude)
println("final U points = ", length(last_slice.u))
println("fit samples = ", nfit)
println("fit window V_EF in [150, 300]")
println("late-time horizon |phi| slope = ", slope)
println("target slope for phi on an extremal horizon: -1")
println("late approximate H/A = ", aretakis_h / ep.amplitude)
println("late r on sampled horizon = ", last_slice.r[end])
println("status: adaptive prototype; threshold convergence remains under study")
