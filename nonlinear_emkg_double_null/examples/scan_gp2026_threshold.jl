using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, argument(index, string(default)))
end

function q_values_argument(index, default, ::Type{T}) where {T<:Real}
    return [parse(T, value) for value in split(argument(index, default), ",")]
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function linear_interpolate(x, y, xq)
    length(x) == length(y) || throw(ArgumentError("x and y lengths differ"))
    xq <= first(x) && return first(y)
    xq >= last(x) && return last(y)
    i = searchsortedlast(x, xq)
    i = min(max(i, firstindex(x)), lastindex(x) - 1)
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - t) * y[i] + t * y[i + 1]
end

function row_rv_profile(row::NLRow)
    n = length(row.v)
    rv = similar(row.r)
    rv[begin] = (row.r[begin + 1] - row.r[begin]) /
                (row.v[begin + 1] - row.v[begin])
    rv[end] = (row.r[end] - row.r[end - 1]) /
              (row.v[end] - row.v[end - 1])
    for j in firstindex(row.v)+1:lastindex(row.v)-1
        rv[j] = (row.r[j + 1] - row.r[j - 1]) /
                (row.v[j + 1] - row.v[j - 1])
    end
    return rv
end

function final_row_expansion_minimum(row::NLRow)
    rv = row_rv_profile(row)
    value, index = findmin(rv)
    return (value=value, V=row.v[index], r=row.r[index], Q=row.Q[index])
end

function first_row_trapped(rows)
    for row in rows
        rv = row_rv_profile(row)
        index = findfirst(<=(zero(eltype(rv))), rv)
        if !isnothing(index)
            return (U=row.u, V=row.v[index], r=row.r[index],
                    Q=row.Q[index], rv=rv[index])
        end
    end
    return nothing
end

function first_slice_trapped(state::AdaptiveNLState)
    for j in 2:length(state.slices)
        rv = adaptive_outgoing_expansion(state.slices[j - 1], state.slices[j])
        crossing = findfirst(<=(zero(eltype(rv))), rv)
        if !isnothing(crossing)
            u = (state.slices[j].u[crossing] + state.slices[j].u[crossing + 1]) / 2
            return (V=state.slices[j].v, U=u, rv=rv[crossing])
        end
    end
    return nothing
end

function final_horizon_properties(state::AdaptiveNLState)
    length(state.slices) >= 2 || return nothing
    lower, upper = state.slices[end - 1], state.slices[end]
    expansion_u, _, rv = outgoing_expansion_profile(lower, upper)
    uh = apparent_horizon_location(expansion_u, rv)
    isnothing(uh) && return nothing

    mass_u, _, mass = renormalized_hawking_mass_profile(lower, upper)
    r = linear_interpolate(upper.u, upper.r, uh)
    q = linear_interpolate(upper.u, upper.Q, uh)
    m = linear_interpolate(mass_u, mass, uh)
    return (U=uh, V=upper.v, r=r, Q=q, M=m,
            one_minus_Q_over_M=one(eltype(upper.r)) - q / m,
            one_minus_r_over_M=one(eltype(upper.r)) - r / m)
end

function run_case(q0, qstar, vmax, dv, amplitude, C, Umax, max_rows)
    T = promote_type(typeof(q0), typeof(qstar), typeof(vmax), typeof(dv),
                     typeof(amplitude), typeof(C), typeof(Umax))
    U0 = parse(T, "-1.0")
    ep = EvolutionParams(
        rn = RNParams(one(T), q0),
        scalar_charge = parse(T, "0.6") / q0,
        amplitude = amplitude,
        omega = one(T),
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T), U1=U0 + parse(T, "0.01"),
                       V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax, C, iterations=10,
                                       max_rows, hyperbolic_charge=true,
                                       step_control=:outer)

    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid for Q0=$q0")
    rows = evolved.rows[1:last_valid]
    final = last(rows)
    du = diff([row.u for row in rows])
    outer_f = exp(last(final.logf))
    next_du = 2C / outer_f
    final_min_rv = final_row_expansion_minimum(final)
    throat = throat_row_diagnostics(final)
    row_trap = first_row_trapped(rows)

    state = adaptive_state_from_u_rows(UAdaptiveNLState(rows))
    slice_trap = first_slice_trapped(state)
    horizon = final_horizon_properties(state)

    dq = qstar - q0
    vtrap = isnothing(slice_trap) ? nothing : slice_trap.V
    scaled_vtrap = !isnothing(vtrap) && dq > 0 ? vtrap * sqrt(dq) : nothing

    return (
        Q0=q0,
        dQ=dq,
        side=dq > 0 ? "BH-side" : dq < 0 ? "dispersive-side" : "quoted-critical",
        rows=length(evolved.rows),
        valid_rows=length(rows),
        invalid_row=last_valid < length(evolved.rows),
        last_U=final.u,
        reached_Umax=final.u == Umax,
        stalled=final.u + next_du == final.u,
        first_du=isempty(du) ? nothing : first(du),
        last_du=isempty(du) ? nothing : last(du),
        next_du=next_du,
        outer_f=outer_f,
        min_rv=final_min_rv.value,
        min_rv_V=final_min_rv.V,
        min_rv_r=final_min_rv.r,
        min_rv_Q=final_min_rv.Q,
        max_rho=throat.max_rho,
        min_throat_y=throat.min_y,
        row_trap_V=isnothing(row_trap) ? nothing : row_trap.V,
        row_trap_U=isnothing(row_trap) ? nothing : row_trap.U,
        slice_trap_V=isnothing(slice_trap) ? nothing : slice_trap.V,
        slice_trap_U=isnothing(slice_trap) ? nothing : slice_trap.U,
        scaled_vtrap=scaled_vtrap,
        final_AH_U=isnothing(horizon) ? nothing : horizon.U,
        final_AH_M=isnothing(horizon) ? nothing : horizon.M,
        final_AH_Q=isnothing(horizon) ? nothing : horizon.Q,
        final_AH_r=isnothing(horizon) ? nothing : horizon.r,
        one_minus_Q_over_M=isnothing(horizon) ? nothing : horizon.one_minus_Q_over_M,
        one_minus_r_over_M=isnothing(horizon) ? nothing : horizon.one_minus_r_over_M,
    )
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Bool && return string(value)
    value isa AbstractString && return value
    value isa Integer && return string(value)
    if value isa Real
        return isfinite(value) ? string(Float64(value)) : string(value)
    end
    return string(value)
end

function print_table(rows)
    names = propertynames(first(rows))
    println(join(string.(names), '\t'))
    for row in rows
        println(join((format_value(getproperty(row, name)) for name in names), '\t'))
    end
end

function run_scan(::Type{T}) where {T<:Real}
    qvalues = q_values_argument(1, "1.0,1.001,1.002,1.003,1.0032,1.0033218,1.004,1.02", T)
    qstar = real_argument(2, "1.0033218", T)
    vmax = real_argument(3, "400.0", T)
    dv = real_argument(4, "0.08", T)
    amplitude = real_argument(5, "0.01", T)
    C = real_argument(6, "0.6", T)
    Umax = real_argument(7, "1.6", T)
    max_rows = integer_argument(8, 180)

    println("# GP2026 Section IIIA threshold scan")
    println("# qstar = ", qstar, ", Vmax = ", vmax, ", Delta V = ", dv,
            ", C = ", C, ", max_rows = ", max_rows)
    println("# Columns with `slice_trap_*` and final_AH_* are the direct Section IIIA diagnostics.")
    println("# Missing trapped-surface columns mean this run only reached the limiting-surface proxy.")
    rows = [run_case(q0, qstar, vmax, dv, amplitude, C, Umax, max_rows)
            for q0 in qvalues]
    print_table(rows)
end

precision_bits = integer_argument(9, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_scan(BigFloat)
    end
else
    run_scan(Float64)
end
