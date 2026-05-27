using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

epsilon = argument_or_default(1, 0.02)
scalar_charge = argument_or_default(2, 0.6)
target_v = argument_or_default(3, 4.0)
dv = argument_or_default(4, 0.02)
initial_nu = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 531

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = scalar_charge,
    amplitude = epsilon,
)
neutral_ep = EvolutionParams(
    rn = ep.rn,
    scalar_charge = 0.0,
    amplitude = epsilon,
)
f0 = mrt2013_degenerate_horizon_f0(neutral_ep)
nv = Int(round(target_v / dv)) + 2
grid = mrt2013_grid(; nu=initial_nu, nv, U0=-5.1, V0=0.0, U1=0.2,
                    V1=(nv - 1) * dv)
state = NLState(grid)
initialize_mrt2013_charged_outgoing_wave!(state, grid, ep; f0)
evolve_nonlinear!(state, grid, ep; iterations=8, subtract_rn_background=true)
adaptive = adaptive_state_from_rectangular(state, grid)

_, sampled_v, _, expected_q_u, q_residual =
    charged_charge_flux_u_profile(adaptive, ep; target_v)
_, _, geometric_q, flux_q, charge_balance_error =
    charged_flux_integrated_charge_profile(adaptive, ep; target_v)
_, _, _, expected_mass_u, mass_residual =
    charged_mass_flux_u_profile(adaptive, ep; target_v, rn_background=ep.rn)
_, _, geometric_mass, flux_mass, mass_balance_error =
    charged_flux_integrated_mass_profile(adaptive, ep; target_v, rn_background=ep.rn)

println("stored scalar convention: Phi = sqrt(32*pi) * phi_GP")
println("epsilon = ", epsilon)
println("scalar charge e = ", scalar_charge)
println("f0 = ", f0)
println("sampled V = ", sampled_v)
println("initial U points = ", initial_nu)
println("Q range at sampled slice pair = ", extrema(geometric_q))
println("max |Q_U - flux source| = ", maximum(abs, q_residual))
println("max |Q - flux-integrated Q| = ", maximum(abs, charge_balance_error))
println("max |varpi_U - charged flux source| = ", maximum(abs, mass_residual))
println("max |varpi - flux-integrated varpi| = ", maximum(abs, mass_balance_error))
println("flux-integrated Q range = ", extrema(flux_q))
println("flux-integrated varpi range = ", extrema(flux_mass))
println("max |expected Q_U| = ", maximum(abs, expected_q_u))
println("max |expected varpi_U| = ", maximum(abs, expected_mass_u))
println("geometric varpi range = ", extrema(geometric_mass))
