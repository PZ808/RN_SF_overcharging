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

function refinement_rate(previous, current)
    (isnothing(previous) || previous <= 0 || current <= 0) && return nothing
    return log(previous / current) / log(2)
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Integer && return string(value)
    value isa AbstractString && return value
    value isa Bool && return string(value)
    value isa Real && return string(value)
    return string(value)
end

function run_level(; q0, amplitude, U1, V1, du, dv, target_u, target_v,
                   iterations, hyperbolic_charge, cell_solver,
                   pulse_leg_gauge)
    ep = EvolutionParams(
        rn = RNParams(1.0, q0),
        scalar_charge = 0.6 / q0,
        amplitude = amplitude,
        omega = 1.0,
    )
    nu = Int(round((U1 + 1.0) / du)) + 1
    nv = Int(round(V1 / dv)) + 1
    grid = gp2026_grid(; nu, nv, U0=-1.0, V0=0.0, U1, V1)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(
        state, grid, ep; pulse_leg_gauge,
    )
    evolve_nonlinear!(state, grid, ep; iterations, reduced_scalar=true,
                      hyperbolic_charge, cell_solver)
    adaptive = adaptive_state_from_rectangular(state, grid)

    _, sampled_v, _, expected_q_u, q_u_residual =
        charged_charge_flux_u_profile(adaptive, ep; target_v, reduced_scalar=true)
    _, sampled_u, _, expected_q_v, q_v_residual =
        charged_charge_flux_v_profile(adaptive, ep; target_u, reduced_scalar=true)
    _, _, _, _, integrated_q_error =
        charged_flux_integrated_charge_profile(adaptive, ep; target_v, reduced_scalar=true)

    return (
        du=du,
        dv=dv,
        nu=nu,
        nv=nv,
        sampled_u=sampled_u,
        sampled_v=sampled_v,
        max_expected_q_u=maximum(abs, expected_q_u),
        max_q_u_residual=maximum(abs, q_u_residual),
        max_expected_q_v=maximum(abs, expected_q_v),
        max_q_v_residual=maximum(abs, q_v_residual),
        max_integrated_q_error=maximum(abs, integrated_q_error),
    )
end

function print_table(rows)
    previous_q_u = nothing
    previous_q_v = nothing
    previous_integrated = nothing
    headers = (
        "level", "Delta_U", "Delta_V", "nu", "nv", "sampled_U", "sampled_V",
        "max_Q_U_residual", "rate_Q_U",
        "max_Q_V_residual", "rate_Q_V",
        "max_integrated_Q_error", "rate_integrated_Q",
        "max_expected_Q_U", "max_expected_Q_V",
    )
    println(join(headers, '\t'))
    for (level, row) in enumerate(rows)
        rate_q_u = refinement_rate(previous_q_u, row.max_q_u_residual)
        rate_q_v = refinement_rate(previous_q_v, row.max_q_v_residual)
        rate_integrated = refinement_rate(previous_integrated, row.max_integrated_q_error)
        values = (
            level - 1, row.du, row.dv, row.nu, row.nv, row.sampled_u, row.sampled_v,
            row.max_q_u_residual, rate_q_u,
            row.max_q_v_residual, rate_q_v,
            row.max_integrated_q_error, rate_integrated,
            row.max_expected_q_u, row.max_expected_q_v,
        )
        println(join((format_value(value) for value in values), '\t'))
        previous_q_u = row.max_q_u_residual
        previous_q_v = row.max_q_v_residual
        previous_integrated = row.max_integrated_q_error
    end
end

function main()
    q0 = real_argument(1, 1.0033218)
    amplitude = real_argument(2, 0.01)
    base_du = real_argument(3, 0.01)
    base_dv = real_argument(4, 0.08)
    levels = integer_argument(5, 4)
    U1 = real_argument(6, -0.8)
    V1 = real_argument(7, 20.0)
    target_u = real_argument(8, -0.9)
    target_v = real_argument(9, 10.0)
    iterations = integer_argument(10, 12)
    charge_mode = argument(11, "hyperbolic")
    charge_mode in ("hyperbolic", "constraint") ||
        throw(ArgumentError("charge mode must be hyperbolic or constraint"))
    hyperbolic_charge = charge_mode == "hyperbolic"
    cell_solver_argument = argument(12, "newton-direct")
    cell_solver_argument in ("newton-direct", "picard-log") ||
        throw(ArgumentError("cell solver must be newton-direct or picard-log"))
    cell_solver = cell_solver_argument == "newton-direct" ?
                  :newton_direct : :picard_log
    pulse_leg_gauge_argument = argument(13, "areal-affine")
    pulse_leg_gauge_argument in ("areal-affine", "ef-affine") ||
        throw(ArgumentError("pulse-leg gauge must be areal-affine or ef-affine"))
    pulse_leg_gauge = pulse_leg_gauge_argument == "areal-affine" ?
                      :areal_affine : :ef_affine

    println("# GP2026 charge residual convergence")
    println("# Q0 = ", q0, ", eQ0 = 0.6, A0 = ", amplitude,
            ", charge_mode = ", charge_mode,
            ", cell_solver = ", cell_solver,
            ", pulse_leg_gauge = ", pulse_leg_gauge)
    println("# U range = [-1, ", U1, "], V range = [0, ", V1, "]")
    println("# target_u = ", target_u, ", target_v = ", target_v,
            ", iterations = ", iterations)
    rows = [
        run_level(; q0, amplitude, U1, V1,
                  du=base_du / 2.0^level,
                  dv=base_dv / 2.0^level,
                  target_u, target_v, iterations, hyperbolic_charge,
                  cell_solver, pulse_leg_gauge)
        for level in 0:levels-1
    ]
    print_table(rows)
end

main()
