using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function boolean_argument(index, default)
    value = lowercase(argument(index, string(default)))
    value in ("true", "t", "1", "yes", "y", "on") && return true
    value in ("false", "f", "0", "no", "n", "off") && return false
    throw(ArgumentError("boolean argument must be true or false"))
end

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, argument(index, string(default)))
end

function step_control_argument(index, default)
    value = argument(index, default)
    value in ("none", "outer", "max-row", "geometric", "throat", "eta", "local") ||
        throw(ArgumentError("step control must be none, outer, max-row, geometric, throat, eta, or local"))
    value == "none" && return :none
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

function row_u_derivative(rows, index)
    length(rows) >= 2 || throw(ArgumentError("at least two rows are required"))
    row = rows[index]
    derivative = similar(row.r)
    if index == firstindex(rows)
        next = rows[index + 1]
        derivative .= (next.r .- row.r) ./ (next.u - row.u)
    elseif index == lastindex(rows)
        previous = rows[index - 1]
        derivative .= (row.r .- previous.r) ./ (row.u - previous.u)
    else
        previous = rows[index - 1]
        next = rows[index + 1]
        derivative .= (next.r .- previous.r) ./ (next.u - previous.u)
    end
    return derivative
end

function logf_rho_array(row::NLRow)
    rho = throat_row_diagnostics(row).rho
    rho_v = N.coordinate_derivative(rho, row.v)
    return [isfinite(rho_v[j]) && rho_v[j] != 0 ?
            row.logf[j] - log(abs(rho_v[j])) :
            typeof(row.logf[j])(Inf)
            for j in eachindex(row.v)]
end

function rational_throat_arrays(row::NLRow)
    throat = throat_row_diagnostics(row)
    return (; throat.eta, throat.zeta)
end

function logf_coordinate_array(row::NLRow, coordinate)
    coordinate_v = N.coordinate_derivative(coordinate, row.v)
    return [isfinite(coordinate_v[j]) && coordinate_v[j] != 0 ?
            row.logf[j] - log(abs(coordinate_v[j])) :
            typeof(row.logf[j])(Inf)
            for j in eachindex(row.v)]
end

function horizon_function_array(row::NLRow, ru, rv)
    return [begin
        mass = renormalized_hawking_mass(row.r[j], exp(row.logf[j]),
                                         ru[j], rv[j], row.Q[j])
        one(row.r[j]) - 2mass / row.r[j] + row.Q[j]^2 / row.r[j]^2
    end for j in eachindex(row.v)]
end

function min_rv_point(row::NLRow, ru)
    rv = row_outgoing_expansion(row)
    horizon = horizon_function_array(row, ru, rv)
    j = last(findmin(rv))
    throat = throat_row_diagnostics(row)
    logf_rho = logf_rho_array(row)
    rational = rational_throat_arrays(row)
    logf_eta = logf_coordinate_array(row, rational.eta)
    logf_zeta = logf_coordinate_array(row, rational.zeta)
    return (; index=j, V=row.v[j], r=row.r[j], Q=row.Q[j],
            y=throat.y[j], rho=throat.rho[j],
            eta=rational.eta[j], zeta=rational.zeta[j],
            r_U=ru[j], r_V=rv[j],
            H=horizon[j], logf=row.logf[j],
            logf_rho=logf_rho[j],
            logf_eta=logf_eta[j],
            logf_zeta=logf_zeta[j],
            max_rho=throat.max_rho)
end

function max_abs_difference(a, b)
    length(a) == length(b) || throw(ArgumentError("arrays must have equal length"))
    return maximum(abs(a[k] - b[k]) for k in eachindex(a)
                   if isfinite(a[k]) && isfinite(b[k]))
end

function row_conditioning(rows, index, C, step_control, max_delta_rho, max_delta_eta)
    row = rows[index]
    ru = row_u_derivative(rows, index)
    point = min_rv_point(row, ru)
    du_info = gp2026_row_step_du(rows[1:index], C, step_control;
                                 max_delta_rho, max_delta_eta)
    local_cap = minimum((du_info.max_row_du, du_info.geometric_du,
                         du_info.throat_du, du_info.eta_du))
    throat = throat_row_diagnostics(row)
    logf_rho = logf_rho_array(row)
    rational = rational_throat_arrays(row)
    logf_eta = logf_coordinate_array(row, rational.eta)
    logf_zeta = logf_coordinate_array(row, rational.zeta)
    rv = row_outgoing_expansion(row)
    horizon = horizon_function_array(row, ru, rv)

    if index > firstindex(rows)
        previous = rows[index - 1]
        du_prev = row.u - previous.u
        previous_throat = throat_row_diagnostics(previous)
        previous_logf_rho = logf_rho_array(previous)
        previous_rational = rational_throat_arrays(previous)
        previous_logf_eta = logf_coordinate_array(previous, previous_rational.eta)
        previous_logf_zeta = logf_coordinate_array(previous, previous_rational.zeta)
        previous_ru = row_u_derivative(rows, index - 1)
        previous_rv = row_outgoing_expansion(previous)
        previous_horizon = horizon_function_array(previous, previous_ru, previous_rv)
        delta_rho = max_abs_difference(throat.rho, previous_throat.rho)
        delta_eta = max_abs_difference(rational.eta, previous_rational.eta)
        delta_zeta = max_abs_difference(rational.zeta, previous_rational.zeta)
        delta_logf = max_abs_difference(row.logf, previous.logf)
        delta_logf_rho = max_abs_difference(logf_rho, previous_logf_rho)
        delta_logf_eta = max_abs_difference(logf_eta, previous_logf_eta)
        delta_logf_zeta = max_abs_difference(logf_zeta, previous_logf_zeta)
        delta_H = max_abs_difference(horizon, previous_horizon)
    else
        du_prev = zero(row.u)
        delta_rho = zero(row.u)
        delta_eta = zero(row.u)
        delta_zeta = zero(row.u)
        delta_logf = zero(row.u)
        delta_logf_rho = zero(row.u)
        delta_logf_eta = zero(row.u)
        delta_logf_zeta = zero(row.u)
        delta_H = zero(row.u)
    end

    return (
        row=index,
        U=row.u,
        Delta_U_prev=du_prev,
        V_min_rV=point.V,
        r=point.r,
        Q=point.Q,
        y=point.y,
        rho=point.rho,
        eta=point.eta,
        zeta=point.zeta,
        r_U=point.r_U,
        r_V=point.r_V,
        H=point.H,
        logf=point.logf,
        logf_rho=point.logf_rho,
        logf_eta=point.logf_eta,
        logf_zeta=point.logf_zeta,
        max_rho=point.max_rho,
        Delta_rho_inf=delta_rho,
        Delta_eta_inf=delta_eta,
        Delta_zeta_inf=delta_zeta,
        Delta_logf_inf=delta_logf,
        Delta_logf_rho_inf=delta_logf_rho,
        Delta_logf_eta_inf=delta_logf_eta,
        Delta_logf_zeta_inf=delta_logf_zeta,
        Delta_H_inf=delta_H,
        rho_U_inf=du_prev > 0 ? delta_rho / du_prev : zero(du_prev),
        eta_U_inf=du_prev > 0 ? delta_eta / du_prev : zero(du_prev),
        zeta_U_inf=du_prev > 0 ? delta_zeta / du_prev : zero(du_prev),
        logf_rho_U_inf=du_prev > 0 ? delta_logf_rho / du_prev : zero(du_prev),
        logf_eta_U_inf=du_prev > 0 ? delta_logf_eta / du_prev : zero(du_prev),
        logf_zeta_U_inf=du_prev > 0 ? delta_logf_zeta / du_prev : zero(du_prev),
        next_outer_Delta_U=du_info.outer_du,
        next_maxrow_Delta_U=du_info.max_row_du,
        next_geometric_Delta_U=du_info.geometric_du,
        next_throat_Delta_U=du_info.throat_du,
        next_eta_Delta_U=du_info.eta_du,
        outer_over_local_cap=isfinite(local_cap) && local_cap > 0 ?
                             du_info.outer_du / local_cap : typeof(row.u)(Inf),
    )
end

function advance_to_target(base::NLRow, target_u, ep::EvolutionParams;
                           pieces::Int, iterations::Int, hyperbolic_charge::Bool,
                           U0, V0, M0)
    pieces >= 1 || throw(ArgumentError("pieces must be positive"))
    target_u > base.u || throw(ArgumentError("target U must exceed base U"))
    current = base
    for k in 1:pieces
        next_u = base.u + (target_u - base.u) * k / pieces
        south = gp2026_na_boundary_point(next_u, ep; U0, V0, M0)
        current = advance_u_row(current, south, ep;
                                iterations, reduced_scalar=true,
                                hyperbolic_charge)
        finite_row(current) || return current
    end
    return current
end

function probe_difference(base::NLRow, one::NLRow, two::NLRow)
    rho_one = throat_row_diagnostics(one).rho
    rho_two = throat_row_diagnostics(two).rho
    rational_one = rational_throat_arrays(one)
    rational_two = rational_throat_arrays(two)
    logf_rho_one = logf_rho_array(one)
    logf_rho_two = logf_rho_array(two)
    logf_eta_one = logf_coordinate_array(one, rational_one.eta)
    logf_eta_two = logf_coordinate_array(two, rational_two.eta)
    logf_zeta_one = logf_coordinate_array(one, rational_one.zeta)
    logf_zeta_two = logf_coordinate_array(two, rational_two.zeta)
    du_one = one.u - base.u
    ru_one = (one.r .- base.r) ./ du_one
    ru_two = (two.r .- base.r) ./ (two.u - base.u)
    rv_one = row_outgoing_expansion(one)
    rv_two = row_outgoing_expansion(two)
    H_one = horizon_function_array(one, ru_one, rv_one)
    H_two = horizon_function_array(two, ru_two, rv_two)
    point_one = min_rv_point(one, ru_one)
    point_two = min_rv_point(two, ru_two)
    return (
        max_abs_r=max_abs_difference(one.r, two.r),
        max_abs_logf=max_abs_difference(one.logf, two.logf),
        max_abs_rho=max_abs_difference(rho_one, rho_two),
        max_abs_eta=max_abs_difference(rational_one.eta, rational_two.eta),
        max_abs_zeta=max_abs_difference(rational_one.zeta, rational_two.zeta),
        max_abs_logf_rho=max_abs_difference(logf_rho_one, logf_rho_two),
        max_abs_logf_eta=max_abs_difference(logf_eta_one, logf_eta_two),
        max_abs_logf_zeta=max_abs_difference(logf_zeta_one, logf_zeta_two),
        max_abs_H=max_abs_difference(H_one, H_two),
        one_min_rV=point_one.r_V,
        two_min_rV=point_two.r_V,
        one_H_min_rV=point_one.H,
        two_H_min_rV=point_two.H,
    )
end

function probe_row(rows, index, ep, C, max_delta_rho, max_delta_eta;
                   iterations, U0, V0, M0)
    base = rows[index]
    du_info = gp2026_row_step_du(rows[1:index], C, :outer;
                                 max_delta_rho, max_delta_eta)
    target_du = du_info.outer_du
    target_u = base.u + target_du
    one = advance_to_target(base, target_u, ep;
                            pieces=1, iterations, hyperbolic_charge=true,
                            U0, V0, M0)
    two = advance_to_target(base, target_u, ep;
                            pieces=2, iterations, hyperbolic_charge=true,
                            U0, V0, M0)
    diff = finite_row(one) && finite_row(two) ?
           probe_difference(base, one, two) :
           (; max_abs_r=Inf, max_abs_logf=Inf, max_abs_rho=Inf,
              max_abs_eta=Inf, max_abs_zeta=Inf,
              max_abs_logf_rho=Inf,
              max_abs_logf_eta=Inf, max_abs_logf_zeta=Inf,
              max_abs_H=Inf,
              one_min_rV=Inf, two_min_rV=Inf,
              one_H_min_rV=Inf, two_H_min_rV=Inf)
    local_cap = minimum((du_info.max_row_du, du_info.geometric_du,
                         du_info.throat_du, du_info.eta_du))
    return (
        row=index,
        U=base.u,
        target_Delta_U=target_du,
        outer_Delta_U=du_info.outer_du,
        maxrow_Delta_U=du_info.max_row_du,
        geometric_Delta_U=du_info.geometric_du,
        throat_Delta_U=du_info.throat_du,
        eta_Delta_U=du_info.eta_du,
        outer_over_local_cap=isfinite(local_cap) && local_cap > 0 ?
                             du_info.outer_du / local_cap : typeof(base.u)(Inf),
        one_finite=finite_row(one),
        two_half_finite=finite_row(two),
        max_abs_r=diff.max_abs_r,
        max_abs_logf=diff.max_abs_logf,
        max_abs_rho=diff.max_abs_rho,
        max_abs_eta=diff.max_abs_eta,
        max_abs_zeta=diff.max_abs_zeta,
        max_abs_logf_rho=diff.max_abs_logf_rho,
        max_abs_logf_eta=diff.max_abs_logf_eta,
        max_abs_logf_zeta=diff.max_abs_logf_zeta,
        max_abs_H=diff.max_abs_H,
        one_min_rV=diff.one_min_rV,
        two_half_min_rV=diff.two_min_rV,
        one_H_min_rV=diff.one_H_min_rV,
        two_half_H_min_rV=diff.two_H_min_rV,
    )
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Bool && return string(value)
    value isa Symbol && return string(value)
    value isa Integer && return string(value)
    value isa Real && return isfinite(value) ? string(Float64(value)) : string(value)
    return string(value)
end

function print_table(rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    println(join(string.(names), '\t'))
    for row in rows
        println(join((format_value(getproperty(row, name)) for name in names), '\t'))
    end
end

function tail_probe_indices(row_count, probe_tail_count)
    probe_tail_count <= 0 && return Int[]
    first_index = max(2, row_count - probe_tail_count + 1)
    return collect(first_index:row_count)
end

function run(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "400.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    Umax = real_argument(6, "1.6", T)
    max_rows = integer_argument(7, 1200)
    step_control = step_control_argument(8, "outer")
    substep_control = step_control_argument(9, "local")
    substep_C = real_argument(10, string(C), T)
    max_delta_rho = real_argument(11, "0.25", T)
    stride = integer_argument(12, 120)
    probe_tail_count = integer_argument(13, 4)
    max_substeps_per_row = integer_argument(14, 10_000)
    iterations = integer_argument(15, 10)
    precision_bits = integer_argument(16, 0)
    max_delta_eta = real_argument(17, "0.025", T)
    backtrack = boolean_argument(18, false)
    backtrack_factor = real_argument(19, "0.5", T)
    max_backtracks = integer_argument(20, 20)
    max_realized_delta_rho = real_argument(21, string(max_delta_rho), T)
    max_realized_delta_eta = real_argument(22, string(max_delta_eta), T)

    U0 = parse(T, "-1.0")
    V0 = zero(T)
    ep = EvolutionParams(
        rn=RNParams(one(T), q0),
        scalar_charge=parse(T, "0.6") / q0,
        amplitude=amplitude,
        omega=one(T),
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0, U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    evolved = evolve_gp2026_u_adaptive(
        initial, ep; Umax, C, iterations, max_rows,
        hyperbolic_charge=true, step_control, max_delta_rho,
        max_delta_eta, substep_control, substep_C, max_substeps_per_row,
        backtrack, backtrack_factor, max_backtracks,
        max_realized_delta_rho, max_realized_delta_eta,
    )
    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    length(rows) >= 2 || error("diagnostic needs at least two valid rows")

    sample_indices = unique(vcat(collect(2:stride:length(rows)), length(rows)))
    row_records = [row_conditioning(rows, index, C, step_control,
                                    max_delta_rho, max_delta_eta)
                   for index in sample_indices]
    probe_records = [probe_row(rows, index, ep, C, max_delta_rho, max_delta_eta;
                               iterations, U0, V0, M0=ep.rn.M)
                     for index in tail_probe_indices(length(rows), probe_tail_count)]

    println("# GP2026 throat-variable stiffness scaling")
    println("# numeric type = ", T, ", precision bits = ", precision(T),
            ", requested precision arg = ", precision_bits)
    println("# Q0 = ", q0, ", eQ0 = ", ep.scalar_charge * q0,
            ", A0 = ", amplitude)
    println("# Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C,
            ", Umax = ", Umax)
    println("# max_rows = ", max_rows, ", valid_rows = ", length(rows),
            ", last U = ", last(rows).u)
    println("# step_control = ", step_control,
            ", substep_control = ", substep_control,
            ", substep_C = ", substep_C,
            ", max_delta_rho = ", max_delta_rho,
            ", max_delta_eta = ", max_delta_eta,
            ", backtrack = ", backtrack,
            ", backtrack_factor = ", backtrack_factor,
            ", max_backtracks = ", max_backtracks,
            ", max_realized_delta_rho = ", max_realized_delta_rho,
            ", max_realized_delta_eta = ", max_realized_delta_eta)
    println("# Row throat conditioning series")
    print_table(row_records)
    println("# Local one-step versus two-half-step probes using the GP outer step")
    print_table(probe_records)
end

precision_bits = integer_argument(16, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run(BigFloat)
    end
else
    run(Float64)
end
