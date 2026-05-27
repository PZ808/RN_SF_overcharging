using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

q0 = argument_or_default(1, 1.0033218)
amplitude = argument_or_default(2, 0.01)
dv = argument_or_default(3, 0.08)
target_u = argument_or_default(4, -0.9)
target_v = argument_or_default(5, 10.0)
U1 = -0.8
du = 0.01

ep = EvolutionParams(
    rn = RNParams(1.0, q0),
    scalar_charge = 0.6 / q0,
    amplitude = amplitude,
    omega = 1.0,
)

nu = Int(round((U1 + 1.0) / du)) + 1
nv = Int(round(20.0 / dv)) + 1
grid = gp2026_grid(; nu, nv, U0=-1.0, V0=0.0, U1, V1=20.0)
state = NLState(grid)
initialize_gp2026_single_pulse!(state, grid, ep)
evolve_nonlinear!(state, grid, ep; iterations=10, reduced_scalar=true)
adaptive = adaptive_state_from_rectangular(state, grid)

_, sampled_v, _, expected_q_u, q_u_residual =
    charged_charge_flux_u_profile(adaptive, ep; target_v, reduced_scalar=true)
_, _, _, _, integrated_q_error =
    charged_flux_integrated_charge_profile(adaptive, ep; target_v, reduced_scalar=true)
_, sampled_u, _, expected_q_v, q_v_residual =
    charged_charge_flux_v_profile(adaptive, ep; target_u, reduced_scalar=true)

println("Gelles-Pretorius 2026 short nonlinear charge-balance check")
println("Q0 = ", q0)
println("e Q0 = ", ep.scalar_charge * q0)
println("A0 = ", amplitude)
println("Delta U = ", du)
println("Delta V = ", dv)
println("sampled U for Q_V = ", sampled_u)
println("sampled V for Q_U = ", sampled_v)
println("max |expected Q_U| = ", maximum(abs, expected_q_u))
println("max |Q_U - flux source| = ", maximum(abs, q_u_residual))
println("max |Q - flux-integrated Q| = ", maximum(abs, integrated_q_error))
println("max |expected Q_V| = ", maximum(abs, expected_q_v))
println("max |Q_V - flux source| = ", maximum(abs, q_v_residual))
