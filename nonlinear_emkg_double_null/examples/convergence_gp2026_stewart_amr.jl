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
    base_du,
    base_dv,
    U0,
    U1,
    V1,
    refinement_factor,
    revision_interval,
    max_levels,
    iterations,
)
    du = base_du / 2.0^level
    dv = base_dv / 2.0^level
    root_steps = Int(round((U1 - U0) / du))
    nv = Int(round(V1 / dv)) + 1
    p = RNParams(1.0, 1.0)
    ep = EvolutionParams(
        rn=p,
        scalar_charge=0.0,
        amplitude=0.0,
        omega=0.0,
    )
    grid = gp2026_grid(
        ; nu=2, nv, U0, V0=0.0, U1=U0 + du, V1,
    )
    state = NLState(grid)
    initialize_gp2026_exact_extremal_rn!(state, grid, ep)
    initial = row_from_rectangular(state, grid, 1)
    config = StewartAMRConfig(
        refinement_factor=refinement_factor,
        revision_interval=revision_interval,
        max_levels=max_levels,
        atol=1.0e-14,
        rtol=1.0e-10,
        buffer_points=4,
    )
    hierarchy = initialize_stewart_hierarchy(initial; config)
    rows = NLRow[initial]
    for step in 1:root_steps
        target_u = step == root_steps ? U1 : U0 + step * du
        result = advance_stewart_hierarchy!(
            hierarchy,
            target_u,
            ep;
            U0,
            V0=0.0,
            M0=1.0,
            pulse_leg_gauge=:ef_affine,
            iterations,
            cell_solver=:newton_direct,
        )
        push!(rows, result.row)
    end

    max_r = maximum(
        abs(
            row.r[j] -
            gp2026_exact_extremal_rn_radius(row.u, row.v[j], p),
        )
        for row in rows for j in eachindex(row.v)
    )
    max_logf = maximum(
        abs(
            row.logf[j] -
            log(gp2026_exact_extremal_rn_fcode(row.u, row.v[j], p)),
        )
        for row in rows for j in eachindex(row.v)
    )
    max_q = maximum(
        abs(row.Q[j] - p.Q0)
        for row in rows for j in eachindex(row.v)
    )
    return (
        level=level,
        du=du,
        dv=dv,
        root_steps=root_steps,
        max_r=max_r,
        max_logf=max_logf,
        max_q=max_q,
        depth=stewart_hierarchy_depth(hierarchy.root),
        level_steps=join(hierarchy.stats.level_steps, ","),
        revisions=join(hierarchy.stats.revisions, ","),
        creations=hierarchy.stats.child_creations,
        destructions=hierarchy.stats.child_destructions,
        injections=hierarchy.stats.injections,
        reintegrations=hierarchy.stats.suffix_reintegrations,
    )
end

format_value(value) = isnothing(value) ? "missing" : string(value)

function main()
    base_du = real_argument(1, 0.04)
    base_dv = real_argument(2, 0.2)
    levels = integer_argument(3, 3)
    U0 = real_argument(4, -0.4)
    U1 = real_argument(5, 0.2)
    V1 = real_argument(6, 10.0)
    refinement_factor = integer_argument(7, 4)
    revision_interval = integer_argument(8, 4)
    max_levels = integer_argument(9, 2)
    iterations = integer_argument(10, 15)

    rows = [
        run_level(
            level;
            base_du,
            base_dv,
            U0,
            U1,
            V1,
            refinement_factor,
            revision_interval,
            max_levels,
            iterations,
        )
        for level in 0:levels-1
    ]
    println("# Exact extremal-RN convergence with persistent Stewart AMR")
    println(
        "# refinement_factor=", refinement_factor,
        ", revision_interval=", revision_interval,
        ", max_levels=", max_levels,
        ", U=[", U0, ",", U1, "], V=[0,", V1, "]",
    )
    println(
        join(
            (
                "level", "Delta_U", "Delta_V", "root_steps",
                "max_r_error", "rate_r", "max_logf_error", "rate_logf",
                "max_Q_error", "depth", "level_steps", "revisions",
                "creations", "destructions", "injections",
                "suffix_reintegrations",
            ),
            '\t',
        ),
    )
    previous = nothing
    for row in rows
        rate_r = isnothing(previous) ? nothing :
                 refinement_rate(previous.max_r, row.max_r)
        rate_logf = isnothing(previous) ? nothing :
                    refinement_rate(previous.max_logf, row.max_logf)
        values = (
            row.level,
            row.du,
            row.dv,
            row.root_steps,
            row.max_r,
            rate_r,
            row.max_logf,
            rate_logf,
            row.max_q,
            row.depth,
            row.level_steps,
            row.revisions,
            row.creations,
            row.destructions,
            row.injections,
            row.reintegrations,
        )
        println(join(format_value.(values), '\t'))
        previous = row
    end
end

main()
