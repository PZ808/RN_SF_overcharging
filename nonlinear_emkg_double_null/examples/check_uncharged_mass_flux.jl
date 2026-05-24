using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

epsilon = argument_or_default(1, 0.02)
target_v = argument_or_default(2, 20.0)
dv = argument_or_default(3, 0.02)
initial_nu = round(Int, argument_or_default(4, 531.0))
adaptive_mode = length(ARGS) >= 5 ? ARGS[5] : "both"
adaptive_mode in ("fixed", "split", "chop", "both") ||
    throw(ArgumentError("adaptive mode must be fixed, split, chop, or both"))
split_threshold = argument_or_default(6, 0.02)

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = epsilon,
)
f0 = mrt2013_degenerate_horizon_f0(ep)
nv = round(Int, target_v / dv) + 2
grid = mrt2013_grid(; nu=initial_nu, nv, U0=-5.1, V0=0.0, U1=0.2,
                    V1=dv * (nv - 1))
seed = NLState(grid)
initialize_mrt2013_outgoing_wave!(seed, grid, ep; f0)

splitting = PointSplittingConfig(
    band_width = 1.0,
    max_relative_r = split_threshold,
    max_relative_f = split_threshold,
    max_dphi = split_threshold * epsilon,
    max_dphiu = split_threshold * epsilon,
)
chopping = HorizonChoppingConfig(
    start_v = 1.0,
    band_width = 1.0,
    interior_buffer_cells = 1,
    interior_buffer_width = 0.05,
)
if adaptive_mode != "fixed"
    split_policy = adaptive_mode in ("split", "both") ? splitting : nothing
    chop_policy = adaptive_mode in ("chop", "both") ? chopping : nothing
    state = evolve_adaptive(
        slice_from_rectangular(seed, grid, 1),
        west_boundary_from_rectangular(seed, grid),
        ep;
        iterations=20,
        subtract_rn_background=true,
        point_splitting=split_policy,
        horizon_chopping=chop_policy,
    )
else
    rectangular = seed
    evolve_nonlinear!(rectangular, grid, ep; iterations=20, subtract_rn_background=true)
    state = adaptive_state_from_rectangular(rectangular, grid)
end

println("epsilon = ", epsilon, ", initial U points = ", initial_nu, ", Delta V = ", dv,
        ", adaptive mode = ", adaptive_mode, ", split threshold = ", split_threshold)
println("final U points = ", length(last(state.slices).u))
for diagnostic_correction in (false, true)
    background = diagnostic_correction ? ep.rn : nothing
    u, sampled_v, mass = bondi_mass_profile(state; target_v, rn_background=background)
    flux_u, _, mass_u, expected_mass_u, residual =
        uncharged_mass_flux_u_profile(state; target_v, rn_background=background)
    _, _, geometric_mass, flux_mass, balance_error =
        uncharged_flux_integrated_mass_profile(state; target_v, rn_background=background)
    expansion_u, _, rv = outgoing_expansion_profile(state; target_v)
    crossing = findfirst(value -> value <= 0, rv)
    exterior_limit = isnothing(crossing) ? last(u) : expansion_u[max(crossing - 1, 1)]
    mass_exterior = u .<= exterior_limit
    exterior = flux_u .<= exterior_limit
    exterior_indices = findall(exterior)
    max_residual_index = exterior_indices[argmax(abs.(residual[exterior]))]
    println("  diagnostic correction = ", diagnostic_correction,
            ", sampled V = ", sampled_v,
            ", minimum exterior varpi = ", minimum(mass[mass_exterior]),
            ", max exterior varpi_U = ", maximum(mass_u[exterior]),
            ", minimum exterior expected varpi_U = ", minimum(expected_mass_u[exterior]),
            ", max exterior expected varpi_U = ", maximum(expected_mass_u[exterior]),
            ", max exterior |flux residual| = ", maximum(abs, residual[exterior]),
            ", at U = ", flux_u[max_residual_index],
            ", min exterior flux-integrated varpi = ", minimum(flux_mass[mass_exterior]),
            ", max exterior |geometric - flux varpi| = ",
            maximum(abs, balance_error[mass_exterior]))
end
