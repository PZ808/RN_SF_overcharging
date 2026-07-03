using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

function real_argument(index, default)
    return parse(Float64, argument(index, default))
end

function integer_argument(index, default)
    return parse(Int, argument(index, default))
end

function refinement_rate(previous, current)
    previous > 0 && current > 0 || return nothing
    return log(previous / current) / log(2)
end

function max_errors(state, grid, p)
    max_r = 0.0
    max_logf = 0.0
    max_q = 0.0
    max_r_exterior = 0.0
    max_r_interior = 0.0
    for i in eachindex(grid.u), j in eachindex(grid.v)
        exact_r = gp2026_exact_extremal_rn_radius(grid.u[i], grid.v[j], p)
        exact_logf =
            log(gp2026_exact_extremal_rn_fcode(grid.u[i], grid.v[j], p))
        r_error = abs(state.r[i, j] - exact_r)
        max_r = max(max_r, r_error)
        max_logf = max(max_logf, abs(state.logf[i, j] - exact_logf))
        max_q = max(max_q, abs(state.Q[i, j] - p.Q0))
        if grid.u[i] < 0
            max_r_exterior = max(max_r_exterior, r_error)
        elseif grid.u[i] > 0
            max_r_interior = max(max_r_interior, r_error)
        end
    end

    horizon_index = findfirst(iszero, grid.u)
    isnothing(horizon_index) &&
        throw(ArgumentError("the electrovacuum grid must contain U=0"))
    horizon_r_error = maximum(abs.(state.r[horizon_index, :] .- p.M))
    horizon_f_error = maximum(abs.(exp.(state.logf[horizon_index, :]) .- 1))
    residuals = cell_equation_residual_summary(
        state, grid,
        EvolutionParams(rn=p, scalar_charge=0.0);
        reduced_scalar=true,
        hyperbolic_charge=true,
        cell_solver=:newton_direct,
    )
    return (
        max_r=max_r,
        max_logf=max_logf,
        max_q=max_q,
        max_r_exterior=max_r_exterior,
        max_r_interior=max_r_interior,
        horizon_r_error=horizon_r_error,
        horizon_f_error=horizon_f_error,
        max_cell_residual=residuals.max_abs_cell_residual,
    )
end

function run_level(level; U0, U1, V1, base_du, base_dv, iterations)
    du = base_du / 2.0^level
    dv = base_dv / 2.0^level
    nu = Int(round((U1 - U0) / du)) + 1
    nv = Int(round(V1 / dv)) + 1
    grid = gp2026_grid(; nu, nv, U0, V0=0.0, U1, V1)
    any(iszero, grid.u) ||
        throw(ArgumentError("choose U bounds and Delta U so the grid contains U=0"))

    ep = EvolutionParams(
        rn=RNParams(1.0, 1.0),
        scalar_charge=0.0,
        amplitude=0.0,
        omega=0.0,
    )
    state = NLState(grid)
    initialize_gp2026_exact_extremal_rn!(state, grid, ep)
    evolve_nonlinear!(
        state, grid, ep;
        iterations,
        reduced_scalar=true,
        hyperbolic_charge=true,
        cell_solver=:newton_direct,
    )
    return merge(
        (level=level, du=du, dv=dv, nu=nu, nv=nv),
        max_errors(state, grid, ep.rn),
    )
end

function pulse_leg_errors(U0, V1, dv, gauge)
    nv = Int(round(V1 / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=0.0, U1=U0 + 0.01, V1)
    ep = EvolutionParams(
        rn=RNParams(1.0, 1.0),
        scalar_charge=0.0,
        amplitude=0.0,
        omega=0.0,
    )
    state = NLState(grid)
    initialize_gp2026_single_pulse!(
        state, grid, ep;
        U0,
        pulse_leg_gauge=gauge,
    )
    max_r = maximum(
        abs(
            state.r[1, j] -
            gp2026_exact_extremal_rn_radius(U0, grid.v[j], ep.rn),
        )
        for j in eachindex(grid.v)
    )
    max_logf = maximum(
        abs(
            state.logf[1, j] -
            log(gp2026_exact_extremal_rn_fcode(U0, grid.v[j], ep.rn)),
        )
        for j in eachindex(grid.v)
    )
    return (max_r=max_r, max_logf=max_logf)
end

function format_value(value)
    isnothing(value) && return "missing"
    return string(value)
end

function main()
    base_du = real_argument(1, "0.02")
    base_dv = real_argument(2, "0.1")
    levels = integer_argument(3, "3")
    U0 = real_argument(4, "-0.4")
    U1 = real_argument(5, "0.2")
    V1 = real_argument(6, "10.0")
    iterations = integer_argument(7, "15")

    rows = [
        run_level(level; U0, U1, V1, base_du, base_dv, iterations)
        for level in 0:levels-1
    ]
    println("# Exact extremal-RN evolution in GP2026 MRT coordinates")
    println("# U range = [", U0, ", ", U1, "], V range = [0, ", V1, "]")
    println("# Direct-f Newton cells; horizon U=0 lies inside every grid")
    for gauge in (:ef_affine, :areal_affine)
        leg = pulse_leg_errors(-1.0, V1, base_dv, gauge)
        println(
            "# zero-pulse ", gauge, " U0=-1 initial-leg errors: max_r=",
            leg.max_r, ", max_logf=", leg.max_logf,
        )
    end
    headers = (
        "level", "Delta_U", "Delta_V", "nu", "nv",
        "max_r_error", "rate_r",
        "max_logf_error", "rate_logf",
        "horizon_r_error", "rate_horizon_r",
        "horizon_f_error", "rate_horizon_f",
        "exterior_r_error", "interior_r_error",
        "max_Q_error", "max_cell_residual",
    )
    println(join(headers, '\t'))
    previous = nothing
    for row in rows
        rate_r = isnothing(previous) ? nothing :
                 refinement_rate(previous.max_r, row.max_r)
        rate_logf = isnothing(previous) ? nothing :
                    refinement_rate(previous.max_logf, row.max_logf)
        rate_horizon_r = isnothing(previous) ? nothing :
                         refinement_rate(previous.horizon_r_error,
                                         row.horizon_r_error)
        rate_horizon_f = isnothing(previous) ? nothing :
                         refinement_rate(previous.horizon_f_error,
                                         row.horizon_f_error)
        values = (
            row.level, row.du, row.dv, row.nu, row.nv,
            row.max_r, rate_r,
            row.max_logf, rate_logf,
            row.horizon_r_error, rate_horizon_r,
            row.horizon_f_error, rate_horizon_f,
            row.max_r_exterior, row.max_r_interior,
            row.max_q, row.max_cell_residual,
        )
        println(join(format_value.(values), '\t'))
        previous = row
    end
end

main()
