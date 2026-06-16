using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function real_argument(index, default)
    return parse(Float64, argument(index, string(default)))
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function step_control_argument(index, default)
    value = argument(index, default)
    value in ("outer", "max-row", "geometric", "throat", "eta", "local") ||
        throw(ArgumentError("step control must be outer, max-row, geometric, throat, eta, or local"))
    return value == "outer" ? :outer :
           value == "max-row" ? :max_row :
           value == "geometric" ? :geometric :
           value == "throat" ? :throat :
           value == "eta" ? :eta :
           :local
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function format_value(value)
    value isa Integer && return string(value)
    value isa Real && return string(value)
    return string(value)
end

function print_table(samples, ep; stride=1)
    headers = (
        "row", "U", "V", "rho", "r", "Q", "r_minus_absQ", "logf", "f_code",
        "logf_rho", "rho_V", "r_V", "Q_V",
        "Psi_re", "Psi_im", "abs_Psi",
        "rphi_GP_re", "rphi_GP_im", "abs_rphi_GP", "raw_phase_rphi_GP",
        "abs_rphi_GP_V", "gauge_phase_V", "covariant_DV_abs_rphi",
        "A_U", "A_V", "J_V", "T_VV", "outgoing_constraint_source",
        "Q_V_source", "Q_V_residual", "Q_over_r", "one_minus_absQ_over_r",
    )
    println(join(headers, '\t'))
    for sample in samples[1:stride:end]
        observables = throat_boundary_observables(sample, ep)
        f = exp(sample.logf)
        values = (
            sample.row_index,
            sample.u,
            sample.v,
            sample.rho,
            sample.r,
            sample.q,
            sample.y,
            sample.logf,
            f,
            observables.logf_rho,
            sample.rho_v,
            sample.r_v,
            sample.q_v,
            sample.phi_re,
            sample.phi_im,
            observables.psi_abs,
            observables.rphi_gp_re,
            observables.rphi_gp_im,
            observables.rphi_gp_abs,
            observables.raw_phase,
            observables.rphi_gp_abs_v,
            observables.covariant_phase_v,
            observables.covariant_dv_rphi_abs,
            sample.Au,
            sample.Av,
            observables.Jv,
            observables.Tvv,
            observables.outgoing_constraint_source,
            observables.q_v_source,
            observables.q_v_residual,
            observables.q_over_r,
            observables.one_minus_absq_over_r,
        )
        println(join((format_value(value) for value in values), '\t'))
    end
end

function main()
    q0 = real_argument(1, 1.0033218)
    rho_match = real_argument(2, 2.0)
    vmax = real_argument(3, 400.0)
    dv = real_argument(4, 0.08)
    C = real_argument(5, 0.6)
    max_rows = integer_argument(6, 120)
    Umax = real_argument(7, 1.6)
    amplitude = real_argument(8, 0.01)
    step_control = step_control_argument(9, "outer")
    stride = integer_argument(10, 1)
    max_delta_rho = real_argument(11, 0.25)
    max_delta_eta = real_argument(12, 0.025)
    stride >= 1 || throw(ArgumentError("stride must be at least 1"))

    ep = EvolutionParams(
        rn = RNParams(1.0, q0),
        scalar_charge = 0.6 / q0,
        amplitude = amplitude,
        omega = 1.0,
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0=-1.0, V0=0.0, U1=-0.99, V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    evolved = evolve_gp2026_u_adaptive(
        row_from_rectangular(seed, grid, 1), ep;
        Umax, C, iterations=10, max_rows, hyperbolic_charge=true,
        step_control, max_delta_rho, max_delta_eta,
    )
    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    samples = throat_boundary_series(rows; rho_match, boundary=:outer)

    println("# GP2026 fixed-rho throat boundary extractor")
    println("# Q0 = ", q0, ", eQ0 = ", ep.scalar_charge * q0,
            ", A0 = ", amplitude)
    println("# rho_match = ", rho_match, ", Vmax = ", vmax,
            ", Delta V = ", dv, ", C = ", C,
            ", step_control = ", step_control,
            ", max_delta_rho = ", max_delta_rho,
            ", max_delta_eta = ", max_delta_eta)
    println("# valid rows = ", length(rows), ", samples = ", length(samples),
            ", last valid U = ", last(rows).u)
    println("# The stored scalar Psi equals sqrt(32*pi) * r * phi_GP.")
    print_table(samples, ep; stride)
end

main()
