using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))

function refinement_rate(previous, current)
    previous > 0 && current > 0 || return nothing
    return log(previous / current) / log(2)
end

function run_level(
    level;
    q0,
    amplitude,
    base_du,
    base_dv,
    U1,
    V1,
    target_u,
    target_v,
    iterations,
)
    U0 = -1.0
    du = base_du / 2.0^level
    dv = base_dv / 2.0^level
    root_steps = Int(round((U1 - U0) / du))
    ep = EvolutionParams(
        rn=RNParams(1.0, q0),
        scalar_charge=0.6 / q0,
        amplitude=amplitude,
        omega=1.0,
    )
    grid = gp2026_grid(
        ; nu=2,
        nv=Int(round(V1 / dv)) + 1,
        U0,
        V0=0.0,
        U1=U0 + du,
        V1,
    )
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep)
    initial = row_from_rectangular(state, grid, 1)
    config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=4,
        max_levels=2,
        atol=1.0e-14,
        rtol=1.0e-10,
        buffer_points=4,
    )
    hierarchy = initialize_stewart_hierarchy(initial; config)
    rows = NLRow[initial]
    for step in 1:root_steps
        next_u = step == root_steps ? U1 : U0 + step * du
        result = advance_stewart_hierarchy!(
            hierarchy,
            next_u,
            ep;
            U0,
            pulse_leg_gauge=:areal_affine,
            iterations,
            cell_solver=:newton_direct,
        )
        push!(rows, result.row)
    end

    adaptive = adaptive_state_from_u_rows(UAdaptiveNLState(rows))
    _, sampled_v, _, _, q_u_residual = charged_charge_flux_u_profile(
        adaptive,
        ep;
        target_v,
        reduced_scalar=true,
    )
    _, sampled_u, _, _, q_v_residual = charged_charge_flux_v_profile(
        adaptive,
        ep;
        target_u,
        reduced_scalar=true,
    )
    return (
        level=level,
        du=du,
        dv=dv,
        sampled_u=sampled_u,
        sampled_v=sampled_v,
        max_q_u=maximum(abs, q_u_residual),
        max_q_v=maximum(abs, q_v_residual),
        depth=stewart_hierarchy_depth(hierarchy.root),
        level_steps=join(hierarchy.stats.level_steps, ","),
        injections=hierarchy.stats.injections,
        reintegrations=hierarchy.stats.suffix_reintegrations,
    )
end

format_value(value) = isnothing(value) ? "missing" : string(value)

function main()
    q0 = real_argument(1, 1.0)
    amplitude = real_argument(2, 0.01)
    base_du = real_argument(3, 0.02)
    base_dv = real_argument(4, 0.16)
    levels = integer_argument(5, 3)
    U1 = real_argument(6, -0.8)
    V1 = real_argument(7, 20.0)
    target_u = real_argument(8, -0.9)
    target_v = real_argument(9, 10.0)
    iterations = integer_argument(10, 12)

    rows = [
        run_level(
            level;
            q0,
            amplitude,
            base_du,
            base_dv,
            U1,
            V1,
            target_u,
            target_v,
            iterations,
        )
        for level in 0:levels-1
    ]
    println("# Charged GP2026 constraint convergence with active Stewart AMR")
    println(
        "# Q0=", q0,
        ", eQ0=0.6, A0=", amplitude,
        ", U=[-1,", U1, "], V=[0,", V1, "]",
    )
    println(
        join(
            (
                "level", "Delta_U", "Delta_V", "sampled_U", "sampled_V",
                "max_Q_U_residual", "rate_Q_U",
                "max_Q_V_residual", "rate_Q_V",
                "depth", "level_steps", "injections",
                "suffix_reintegrations",
            ),
            '\t',
        ),
    )
    previous = nothing
    for row in rows
        rate_q_u = isnothing(previous) ? nothing :
                   refinement_rate(previous.max_q_u, row.max_q_u)
        rate_q_v = isnothing(previous) ? nothing :
                   refinement_rate(previous.max_q_v, row.max_q_v)
        values = (
            row.level,
            row.du,
            row.dv,
            row.sampled_u,
            row.sampled_v,
            row.max_q_u,
            rate_q_u,
            row.max_q_v,
            rate_q_v,
            row.depth,
            row.level_steps,
            row.injections,
            row.reintegrations,
        )
        println(join(format_value.(values), '\t'))
        previous = row
    end
end

main()
