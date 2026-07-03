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
    value isa Real && return isfinite(value) ? string(value) : string(value)
    return string(value)
end

function run_level(; q0, amplitude, U1, V1, du, dv, iterations,
                   cell_solver, pulse_leg_gauge)
    ep = EvolutionParams(
        rn=RNParams(1.0, q0),
        scalar_charge=0.6 / q0,
        amplitude=amplitude,
        omega=1.0,
    )
    nu = Int(round((U1 + 1.0) / du)) + 1
    nv = Int(round(V1 / dv)) + 1
    grid = gp2026_grid(; nu, nv, U0=-1.0, V0=0.0, U1, V1)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(
        state, grid, ep; pulse_leg_gauge,
    )
    evolve_nonlinear!(state, grid, ep; iterations, reduced_scalar=true,
                      hyperbolic_charge=true, cell_solver)
    summary = cell_equation_residual_summary(
        state, grid, ep; reduced_scalar=true, hyperbolic_charge=true,
        cell_solver,
    )
    return merge((du=du, dv=dv, nu=nu, nv=nv), summary)
end

function print_table(rows)
    headers = (
        "level", "Delta_U", "Delta_V", "nu", "nv", "cells",
        "cell_residual",
        "r_UV", "rate_r",
        "f_UV", "rate_f",
        "logf_UV", "rate_logf",
        "psi_re_UV", "rate_psi_re",
        "psi_im_UV", "rate_psi_im",
        "Q_UV", "rate_Q",
        "A_U,V", "A_V,U", "Lorenz", "Faraday",
        "Q_U_constraint", "Q_V_constraint",
        "logf_GP_literal", "logf_extra_Coulomb",
    )
    println(join(headers, '\t'))
    previous = nothing
    for (level, row) in enumerate(rows)
        rate_r = isnothing(previous) ? nothing :
                 refinement_rate(previous.max_abs_r_uv, row.max_abs_r_uv)
        rate_logf = isnothing(previous) ? nothing :
                    refinement_rate(previous.max_abs_logf_uv,
                                    row.max_abs_logf_uv)
        rate_f = isnothing(previous) ? nothing :
                 refinement_rate(previous.max_abs_f_uv, row.max_abs_f_uv)
        rate_psi_re = isnothing(previous) ? nothing :
                      refinement_rate(previous.max_abs_psi_re_uv,
                                      row.max_abs_psi_re_uv)
        rate_psi_im = isnothing(previous) ? nothing :
                      refinement_rate(previous.max_abs_psi_im_uv,
                                      row.max_abs_psi_im_uv)
        rate_q = isnothing(previous) ? nothing :
                 refinement_rate(previous.max_abs_q_uv, row.max_abs_q_uv)
        values = (
            level - 1, row.du, row.dv, row.nu, row.nv, row.cells,
            row.max_abs_cell_residual,
            row.max_abs_r_uv, rate_r,
            row.max_abs_f_uv, rate_f,
            row.max_abs_logf_uv, rate_logf,
            row.max_abs_psi_re_uv, rate_psi_re,
            row.max_abs_psi_im_uv, rate_psi_im,
            row.max_abs_q_uv, rate_q,
            row.max_abs_au_v, row.max_abs_av_u,
            row.max_abs_quasilorenz, row.max_abs_faraday,
            row.max_abs_q_u_constraint, row.max_abs_q_v_constraint,
            row.max_abs_logf_gp_literal, row.max_abs_logf_coulomb2,
        )
        println(join((format_value(value) for value in values), '\t'))
        previous = row
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
    iterations = integer_argument(8, 12)
    cell_solver_argument = argument(9, "newton-direct")
    cell_solver_argument in ("newton-direct", "picard-log") ||
        throw(ArgumentError("cell solver must be newton-direct or picard-log"))
    cell_solver = cell_solver_argument == "newton-direct" ?
                  :newton_direct : :picard_log
    pulse_leg_gauge_argument = argument(10, "areal-affine")
    pulse_leg_gauge_argument in ("areal-affine", "ef-affine") ||
        throw(ArgumentError("pulse-leg gauge must be areal-affine or ef-affine"))
    pulse_leg_gauge = pulse_leg_gauge_argument == "areal-affine" ?
                      :areal_affine : :ef_affine

    println("# GP2026 centered-cell residual audit")
    println("# Q0 = ", q0, ", eQ0 = 0.6, A0 = ", amplitude)
    println("# U range = [-1, ", U1, "], V range = [0, ", V1, "]")
    println("# iterations = ", iterations,
            ", cell_solver = ", cell_solver,
            ", pulse_leg_gauge = ", pulse_leg_gauge)
    rows = [
        run_level(; q0, amplitude, U1, V1,
                  du=base_du / 2.0^level,
                  dv=base_dv / 2.0^level,
                  iterations, cell_solver, pulse_leg_gauge)
        for level in 0:levels-1
    ]
    print_table(rows)
end

main()
