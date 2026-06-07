using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, argument(index, string(default)))
end

function step_control_argument(index, default="outer")
    value = argument(index, default)
    value in ("none", "outer", "max-row", "geometric", "throat", "local") ||
        throw(ArgumentError("step control must be none, outer, max-row, geometric, throat, or local"))
    value == "none" && return :none
    return value == "outer" ? :outer :
           value == "max-row" ? :max_row :
           value == "geometric" ? :geometric :
           value == "throat" ? :throat :
           :local
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function row_u_derivative(rows, index)
    1 <= index <= length(rows) || throw(BoundsError(rows, index))
    length(rows) >= 2 || throw(ArgumentError("at least two rows are required"))
    row = rows[index]
    derivative = similar(row.r)
    if index == firstindex(rows)
        next = rows[index + 1]
        du = next.u - row.u
        derivative .= (next.r .- row.r) ./ du
    elseif index == lastindex(rows)
        previous = rows[index - 1]
        du = row.u - previous.u
        derivative .= (row.r .- previous.r) ./ du
    else
        previous = rows[index - 1]
        next = rows[index + 1]
        du = next.u - previous.u
        derivative .= (next.r .- previous.r) ./ du
    end
    return derivative
end

function first_crossing_coordinate(x, values; target=0)
    for j in firstindex(values):lastindex(values)-1
        left = values[j] - target
        right = values[j + 1] - target
        if left == 0
            return x[j]
        elseif right == 0
            return x[j + 1]
        elseif signbit(left) != signbit(right)
            t = -left / (right - left)
            return (one(t) - t) * x[j] + t * x[j + 1]
        end
    end
    return nothing
end

function point_record(rows, index; point=:min_rv)
    row = rows[index]
    rv = row_outgoing_expansion(row)
    ru = row_u_derivative(rows, index)
    mass = [renormalized_hawking_mass(row.r[j], exp(row.logf[j]),
                                      ru[j], rv[j], row.Q[j])
            for j in eachindex(row.v)]
    horizon_function = [one(row.r[j]) - 2mass[j] / row.r[j] +
                        row.Q[j]^2 / row.r[j]^2
                        for j in eachindex(row.v)]

    j = if point === :min_mass_gap
        last(findmin(horizon_function))
    elseif point === :min_r_minus_absq
        last(findmin(row.r .- abs.(row.Q)))
    else
        last(findmin(rv))
    end

    r = row.r[j]
    q = row.Q[j]
    m = mass[j]
    qabs = max(abs(q), sqrt(eps(float(one(r)))))
    y = r - qabs
    y_for_rho = max(y, sqrt(eps(float(one(r)))) * qabs)
    rho = -log(y_for_rho / qabs)
    theta_u = 2ru[j] / r
    theta_v = 2rv[j] / r
    expansion_class = ru[j] < 0 && rv[j] < 0 ? "future_trapped" :
                      ru[j] > 0 && rv[j] > 0 ? "past_trapped" :
                      ru[j] < 0 && rv[j] > 0 ? "exterior" :
                      "interior_white_hole_side"
    horizon_gap = horizon_function[j]
    gradient_gap = -4ru[j] * rv[j] / exp(row.logf[j])
    discriminant = m^2 - q^2
    rplus = discriminant >= 0 ? m + sqrt(discriminant) : typeof(m)(NaN)
    crossing = row_apparent_horizon_crossing(row; row_index=index)
    mass_gap_crossing_v = first_crossing_coordinate(row.v, horizon_function)

    return (
        row=index,
        U=row.u,
        point=String(point),
        V=row.v[j],
        r=r,
        Q=q,
        M=m,
        r_U=ru[j],
        r_V=rv[j],
        theta_U=theta_u,
        theta_V=theta_v,
        expansion_class=expansion_class,
        future_trapped=ru[j] < 0 && rv[j] < 0,
        horizon_function=horizon_gap,
        minus4_ru_rv_over_f=gradient_gap,
        horizon_identity_error=horizon_gap - gradient_gap,
        one_minus_absQ_over_r=one(r) - abs(q) / r,
        rho=rho,
        Q_over_M=q / m,
        one_minus_Q_over_M=one(m) - q / m,
        r_minus_rplus=isfinite(rplus) ? r - rplus : typeof(r)(NaN),
        direct_trap_V=isnothing(crossing) ? nothing : crossing.v,
        mass_gap_crossing_V=mass_gap_crossing_v,
        min_row_rV=minimum(rv),
        min_row_horizon_function=minimum(horizon_function),
    )
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Bool && return string(value)
    value isa AbstractString && return value
    value isa Symbol && return string(value)
    value isa Integer && return string(value)
    value isa Real && return isfinite(value) ? string(Float64(value)) : string(value)
    return string(value)
end

function print_table(rows)
    names = propertynames(first(rows))
    println(join(string.(names), '\t'))
    for row in rows
        println(join((format_value(getproperty(row, name)) for name in names), '\t'))
    end
end

function linear_fit(x, y)
    length(x) == length(y) || throw(ArgumentError("fit arrays must have the same length"))
    length(x) >= 2 || throw(ArgumentError("need at least two samples"))
    xbar = sum(x) / length(x)
    ybar = sum(y) / length(y)
    denom = sum((x[k] - xbar)^2 for k in eachindex(x))
    denom > 0 || throw(ArgumentError("fit abscissae must not all match"))
    slope = sum((x[k] - xbar) * (y[k] - ybar) for k in eachindex(x)) / denom
    intercept = ybar - slope * xbar
    residual = sqrt(sum((y[k] - (intercept + slope * x[k]))^2
                        for k in eachindex(x)) / length(x))
    return (; intercept, slope, residual)
end

function positive_tail_fit(x, y)
    positive = [(x[k], y[k]) for k in eachindex(x)
                if isfinite(x[k]) && isfinite(y[k]) && x[k] > 0 && y[k] > 0]
    length(positive) >= 2 || return nothing
    logx = log.([pair[1] for pair in positive])
    logy = log.([pair[2] for pair in positive])
    fit = linear_fit(logx, logy)
    return (; exponent=-fit.slope, intercept=fit.intercept,
            residual=fit.residual, samples=length(positive))
end

function floor_tail_fit(V, y)
    samples = [(V[k], y[k]) for k in eachindex(V)
               if isfinite(V[k]) && isfinite(y[k]) && V[k] > 0]
    length(samples) >= 2 || return nothing
    invV = [inv(pair[1]) for pair in samples]
    values = [pair[2] for pair in samples]
    fit = linear_fit(invV, values)
    crossing_V = fit.intercept * fit.slope < 0 ? -fit.slope / fit.intercept : nothing
    return (; y_inf=fit.intercept, invV_slope=fit.slope,
            residual=fit.residual, crossing_V, samples=length(samples))
end

function print_tail_fit(label, records)
    V = [record.V for record in records]
    values = [getproperty(record, label) for record in records]
    power = positive_tail_fit(V, values)
    floor = floor_tail_fit(V, values)
    println("# tail fit for ", label)
    println("#   positive power-law-to-zero fit = ", power)
    println("#   floor fit y(V)=y_inf+a/V = ", floor)
end

function run_analysis(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "400.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    Umax = real_argument(6, "1.6", T)
    max_rows = integer_argument(7, 120)
    step_control = step_control_argument(8)
    max_delta_rho = real_argument(9, "0.25", T)
    stride = integer_argument(10, 10)
    substep_control = step_control_argument(11, "none")
    substep_C = real_argument(12, string(C), T)
    max_substeps_per_row = integer_argument(13, 10_000)
    tail_fit_count = integer_argument(14, 200)
    precision_bits = integer_argument(15, 0)

    U0 = parse(T, "-1.0")
    ep = EvolutionParams(
        rn=RNParams(one(T), q0),
        scalar_charge=parse(T, "0.6") / q0,
        amplitude=amplitude,
        omega=one(T),
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T),
                       U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax, C,
                                       iterations=10, max_rows,
                                       hyperbolic_charge=true,
                                       step_control, max_delta_rho,
                                       substep_control, substep_C,
                                       max_substeps_per_row)
    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    length(rows) >= 2 || error("analysis needs at least two valid rows")

    stalled = last(rows).u + 2C / exp(last(last(rows).logf)) == last(rows).u
    missing_status = last_valid < length(evolved.rows) ? :invalid_row :
                     last(rows).u == Umax ? :reached_umax :
                     stalled ? :precision_stalled :
                     length(evolved.rows) >= max_rows ? :max_rows :
                     :stopped
    vtrap = vtrap_diagnostic(rows; missing_status)

    records = NamedTuple[]
    sample_indices = unique(vcat(collect(1:stride:length(rows)), length(rows)))
    for index in sample_indices
        push!(records, point_record(rows, index; point=:min_rv))
    end

    tail_count = min(tail_fit_count, length(rows))
    tail_indices = collect(length(rows)-tail_count+1:length(rows))
    tail_records = [point_record(rows, index; point=:min_rv)
                    for index in tail_indices]

    global_min = point_record(rows, vtrap.closest.row_index; point=:min_rv)
    final_min_mass_gap = point_record(rows, length(rows); point=:min_mass_gap)

    println("# GP2026 trapped-surface analysis")
    println("# numeric type = ", T, ", precision bits = ", precision(T),
            ", requested precision arg = ", precision_bits)
    println("# Q0 = ", q0, ", eQ0 = ", ep.scalar_charge * q0,
            ", A0 = ", amplitude)
    println("# Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C,
            ", step_control = ", step_control,
            ", substep_control = ", substep_control,
            ", substep_C = ", substep_C,
            ", max rows = ", max_rows)
    println("# valid rows = ", length(rows), ", last U = ", last(rows).u,
            ", Vtrap status = ", vtrap.status, ", direct trap = ", vtrap.trap)
    println("# global closest min-rV point = ", global_min)
    println("# final-row min horizon-function point = ", final_min_mass_gap)
    println("# tail fit samples = ", length(tail_records))
    print_tail_fit(:r_V, tail_records)
    print_tail_fit(:horizon_function, tail_records)
    print_tail_fit(:r_minus_rplus, tail_records)
    println("# Time series at each row's minimum r_V")
    print_table(records)
end

precision_bits = integer_argument(15, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_analysis(BigFloat)
    end
else
    run_analysis(Float64)
end
