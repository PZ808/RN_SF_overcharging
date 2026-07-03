using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

q0 = argument_or_default(1, 1.0033218)
amplitude = argument_or_default(2, 0.01)
omega_tilde = argument_or_default(3, 1.0)
dv = argument_or_default(4, 0.08)

ep = EvolutionParams(
    rn = RNParams(1.0, q0),
    scalar_charge = 0.6 / q0,
    amplitude = amplitude,
    omega = omega_tilde,
)

nu = 261
nv = Int(round(20.0 / dv)) + 1
grid = gp2026_grid(; nu, nv, U0=-1.0, V0=0.0, U1=1.6, V1=20.0)
state = NLState(grid)
initialize_gp2026_single_pulse!(state, grid, ep)
constraints = gp2026_initial_constraint_residuals(state, grid, ep)

r0 = state.r[1, 1]
ru0 = gp2026_extremal_gauge_ru(-1.0)
rv0 = gp2026_extremal_gauge_rv(-1.0, 0.0)
corner_mass = renormalized_hawking_mass(r0, exp(state.logf[1, 1]), ru0, rv0, q0)
peak_index = argmax([gp2026_single_pulse_envelope(V; amplitude) for V in grid.v])
peak_p = hypot(state.phi_re[1, peak_index], state.phi_im[1, peak_index]) /
         sqrt(32 * pi)
adaptive = adaptive_state_from_rectangular(state, grid)
_, _, _, _, q_v_residual =
    charged_charge_flux_v_profile(adaptive, ep; target_u=-1.0, reduced_scalar=true)

println("Gelles-Pretorius 2026 single-pulse initial data")
println("Q0 = ", q0)
println("e Q0 = ", ep.scalar_charge * q0)
println("A0 = ", amplitude)
println("omega_tilde = ", omega_tilde)
println("Delta V = ", dv)
println("corner mass M(U0,V0) = ", corner_mass)
println("r(U=1.6,V0) = ", state.r[end, 1])
println("r_U on N_A = ", (state.r[end, 1] - state.r[1, 1]) /
                          (grid.u[end] - grid.u[1]))
println("peak V sample = ", grid.v[peak_index])
println("peak |r phi_GP| = ", peak_p)
println("charge deposited on N_B = ", state.Q[1, end] - state.Q[1, 1])
println("max initial |Q_V - flux source| = ", maximum(abs, q_v_residual))
println("f_code corner = ", exp(state.logf[1, 1]))
println("f_GP corner = ", exp(state.logf[1, 1]) / 2)
println("corner mass error = ", constraints.corner_mass_error)
println("max N_A radius residual = ", constraints.na_radius)
println("max N_A lapse-constraint residual = ",
        constraints.na_lapse_constraint)
println("max N_A charge residual = ", constraints.na_charge)
println("max N_A scalar residual = ", constraints.na_scalar)
println("max N_A A_U residual = ", constraints.na_au)
println("max N_A A_V-constraint residual = ",
        constraints.na_av_constraint)
println("max N_B radius residual = ", constraints.nb_radius)
println("max N_B scalar residual = ", constraints.nb_scalar)
println("max N_B charge-constraint residual = ",
        constraints.nb_charge_constraint)
println("max N_B lapse-constraint residual = ",
        constraints.nb_lapse_constraint)
println("max N_B A_U-constraint residual = ",
        constraints.nb_au_constraint)
println("max N_B A_V residual = ", constraints.nb_av)
println("corner Faraday residual = ", constraints.faraday_corner)
