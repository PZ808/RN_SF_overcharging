using NonlinearEMKGDoubleNull
using Printf

const N = NonlinearEMKGDoubleNull

const COLORS = (
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
    "#17becf", "#8c564b", "#7f7f7f",
)

struct PlotSeries{T<:Real}
    x::Vector{T}
    y::Vector{T}
    label::String
    color::String
end

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

function svg_escape(text)
    return replace(string(text), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
end

function tick_label(value)
    value == 0 && return "0"
    absvalue = abs(value)
    if absvalue < 1.0e-2 || absvalue >= 1.0e4
        return @sprintf("%.1e", value)
    elseif absvalue < 1
        return @sprintf("%.3f", value)
    end
    return @sprintf("%.3g", value)
end

function transformed(value, scale::Symbol)
    scale === :log10 && return log10(value)
    scale === :linear && return value
    throw(ArgumentError("scale must be :linear or :log10"))
end

function inverse_transformed(value, scale::Symbol)
    scale === :log10 && return 10.0^value
    scale === :linear && return value
    throw(ArgumentError("scale must be :linear or :log10"))
end

function linear_ticks(lo, hi; count=5)
    hi > lo || return [lo]
    return [lo + (hi - lo) * (k - 1) / (count - 1) for k in 1:count]
end

function log_ticks(lo, hi)
    lo > 0 || throw(ArgumentError("log ticks require positive lower bound"))
    first_power = ceil(Int, log10(lo))
    last_power = floor(Int, log10(hi))
    ticks = [10.0^power for power in first_power:last_power]
    if isempty(ticks)
        return [lo, sqrt(lo * hi), hi]
    end
    return ticks
end

function finite_domain(series::Vector{<:PlotSeries}, scale::Symbol, selector)
    values = Float64[]
    for item in series
        for value in selector(item)
            if isfinite(value) && (scale === :linear || value > 0)
                push!(values, Float64(value))
            end
        end
    end
    isempty(values) && throw(ArgumentError("no finite values available for plot"))
    lo, hi = extrema(values)
    if lo == hi
        delta = scale === :log10 ? lo / 10 : max(abs(lo), 1.0) * 0.05
        lo -= delta
        hi += delta
    end
    if scale === :log10
        lo /= 1.2
        hi *= 1.2
    else
        pad = 0.05 * (hi - lo)
        lo -= pad
        hi += pad
    end
    return lo, hi
end

function has_finite_pairs(series::Vector{<:PlotSeries}, xscale::Symbol, yscale::Symbol)
    for item in series
        for k in eachindex(item.x)
            x, y = item.x[k], item.y[k]
            isfinite(x) && isfinite(y) || continue
            (xscale === :linear || x > 0) && (yscale === :linear || y > 0) || continue
            return true
        end
    end
    return false
end

function write_empty_svg(path; title, message, width=920, height=560)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="18" font-weight="600">$(svg_escape(title))</text>""")
        println(io, """<text x="$(width / 2)" y="$(height / 2)" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#555">$(svg_escape(message))</text>""")
        println(io, "</svg>")
    end
    return path
end

function write_svg_line_plot(path, series::Vector{<:PlotSeries};
                             title, xlabel, ylabel,
                             xscale::Symbol=:linear, yscale::Symbol=:linear,
                             width=920, height=560)
    mkpath(dirname(path))
    has_finite_pairs(series, xscale, yscale) ||
        return write_empty_svg(path; title,
                               message="No finite positive data for this scale",
                               width, height)
    left, right, top, bottom = 82.0, 215.0, 60.0, 74.0
    plot_w = width - left - right
    plot_h = height - top - bottom
    xmin, xmax = finite_domain(series, xscale, item -> item.x)
    ymin, ymax = finite_domain(series, yscale, item -> item.y)
    txmin, txmax = transformed(xmin, xscale), transformed(xmax, xscale)
    tymin, tymax = transformed(ymin, yscale), transformed(ymax, yscale)

    sx(x) = left + (transformed(Float64(x), xscale) - txmin) / (txmax - txmin) * plot_w
    sy(y) = top + (tymax - transformed(Float64(y), yscale)) / (tymax - tymin) * plot_h

    xticks = xscale === :log10 ? log_ticks(xmin, xmax) : linear_ticks(xmin, xmax)
    yticks = yscale === :log10 ? log_ticks(ymin, ymax) : linear_ticks(ymin, ymax)

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(width / 2)" y="30" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="18" font-weight="600">$(svg_escape(title))</text>""")
        println(io, """<rect x="$left" y="$top" width="$plot_w" height="$plot_h" fill="#fafafa" stroke="#222" stroke-width="1"/>""")

        for tick in xticks
            x = sx(tick)
            println(io, """<line x1="$x" y1="$top" x2="$x" y2="$(top + plot_h)" stroke="#ddd" stroke-width="1"/>""")
            println(io, """<text x="$x" y="$(top + plot_h + 24)" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="12">$(tick_label(tick))</text>""")
        end
        for tick in yticks
            y = sy(tick)
            println(io, """<line x1="$left" y1="$y" x2="$(left + plot_w)" y2="$y" stroke="#ddd" stroke-width="1"/>""")
            println(io, """<text x="$(left - 10)" y="$(y + 4)" text-anchor="end" font-family="Helvetica, Arial, sans-serif" font-size="12">$(tick_label(tick))</text>""")
        end

        for item in series
            points = String[]
            for k in eachindex(item.x)
                x, y = item.x[k], item.y[k]
                isfinite(x) && isfinite(y) || continue
                (xscale === :linear || x > 0) && (yscale === :linear || y > 0) || continue
                push!(points, @sprintf("%.3f,%.3f", sx(x), sy(y)))
            end
            length(points) >= 2 || continue
            println(io, """<polyline points="$(join(points, " "))" fill="none" stroke="$(item.color)" stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/>""")
        end

        println(io, """<text x="$(left + plot_w / 2)" y="$(height - 20)" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="14">$(svg_escape(xlabel))</text>""")
        println(io, """<text x="20" y="$(top + plot_h / 2)" text-anchor="middle" transform="rotate(-90 20 $(top + plot_h / 2))" font-family="Helvetica, Arial, sans-serif" font-size="14">$(svg_escape(ylabel))</text>""")

        legend_x = left + plot_w + 24
        legend_y = top + 20
        for (index, item) in pairs(series)
            y = legend_y + 24 * (index - 1)
            println(io, """<line x1="$legend_x" y1="$y" x2="$(legend_x + 24)" y2="$y" stroke="$(item.color)" stroke-width="3"/>""")
            println(io, """<text x="$(legend_x + 32)" y="$(y + 4)" font-family="Helvetica, Arial, sans-serif" font-size="13">$(svg_escape(item.label))</text>""")
        end
        println(io, "</svg>")
    end
    return path
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

function trapped_tail_records(q0, ::Type{T}; vmax, dv, amplitude, C, Umax,
                              max_rows, stride, substep_C) where {T<:Real}
    ep = EvolutionParams(
        rn=RNParams(one(T), q0),
        scalar_charge=parse(T, "0.6") / q0,
        amplitude=amplitude,
        omega=one(T),
    )
    U0 = parse(T, "-1.0")
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T),
                       U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax, C,
                                       iterations=10, max_rows,
                                       hyperbolic_charge=true,
                                       step_control=:outer,
                                       substep_control=:local,
                                       substep_C)
    last_valid = findlast(finite_row, evolved.rows)
    rows = evolved.rows[1:last_valid]
    sample_indices = unique(vcat(collect(2:stride:length(rows)), length(rows)))
    records = NamedTuple[]
    for index in sample_indices
        row = rows[index]
        rv = row_outgoing_expansion(row)
        ru = row_u_derivative(rows, index)
        _, j = findmin(rv)
        r = row.r[j]
        q = row.Q[j]
        f = exp(row.logf[j])
        mass = renormalized_hawking_mass(r, f, ru[j], rv[j], q)
        horizon_function = one(r) - 2mass / r + q^2 / r^2
        discriminant = mass^2 - q^2
        rplus = discriminant >= 0 ? mass + sqrt(discriminant) : T(NaN)
        push!(records, (
            row=index,
            U=row.u,
            V=row.v[j],
            r_V=rv[j],
            horizon_function=horizon_function,
            r_minus_rplus=isfinite(rplus) ? r - rplus : T(NaN),
            max_rho=throat_row_diagnostics(row).max_rho,
        ))
    end
    return records
end

function finite_difference(values, coordinates)
    derivative = similar(values)
    derivative[begin] = (values[begin + 1] - values[begin]) /
                        (coordinates[begin + 1] - coordinates[begin])
    derivative[end] = (values[end] - values[end - 1]) /
                      (coordinates[end] - coordinates[end - 1])
    for j in firstindex(values)+1:lastindex(values)-1
        derivative[j] = (values[j + 1] - values[j - 1]) /
                        (coordinates[j + 1] - coordinates[j - 1])
    end
    return derivative
end

function extremal_F_prime(r, M)
    return 2M / r^2 - 2M^2 / r^3
end

function initial_pulse_sources(V, ep, M0, U0, V0, width)
    r = gp2026_extremal_gauge_initial_radius(U0, V; U0, V0, M0)
    rv = gp2026_extremal_gauge_rv(U0, V; U0, V0, M0)
    envelope = gp2026_single_pulse_envelope(V; amplitude=ep.amplitude, width)
    envelope_v = N.gp2026_single_pulse_envelope_derivative(V;
                                                           amplitude=ep.amplitude,
                                                           width)
    phase = ep.omega * V
    psi_re = sqrt(32pi) * envelope * cos(phase)
    psi_im = -sqrt(32pi) * envelope * sin(phase)
    psi_v_re = sqrt(32pi) * (envelope_v * cos(phase) -
                             ep.omega * envelope * sin(phase))
    psi_v_im = sqrt(32pi) * (-envelope_v * sin(phase) -
                             ep.omega * envelope * cos(phase))
    phi_re = psi_re / r
    phi_im = psi_im / r
    phi_v_re = (psi_v_re - rv * phi_re) / r
    phi_v_im = (psi_v_im - rv * phi_im) / r
    _, Jv = current_components(phi_re, phi_im,
                               zero(phi_v_re), zero(phi_v_im),
                               phi_v_re, phi_v_im, ep.scalar_charge)
    rv_v = extremal_F_prime(r, M0) * rv / 2
    logf_v = rv_v / rv + r * (phi_v_re^2 + phi_v_im^2) / (4rv)
    q_v = -r^2 * Jv / 8
    return (; q_v, logf_v)
end

function initial_data_residuals(::Type{T}; q0, vmax, dvs, amplitude, omega,
                                width) where {T<:Real}
    ep = EvolutionParams(rn=RNParams(one(T), q0),
                         scalar_charge=parse(T, "0.6") / q0,
                         amplitude=amplitude,
                         omega=omega)
    q_residual = T[]
    logf_residual = T[]
    U0 = parse(T, "-1.0")
    V0 = zero(T)
    for dv in dvs
        nv = Int(round(vmax / dv)) + 1
        grid = gp2026_grid(; nu=2, nv, U0, V0, U1=U0 + parse(T, "0.01"),
                           V1=vmax)
        state = NLState(grid)
        initialize_gp2026_single_pulse!(state, grid, ep; U0, V0, width,
                                        M0=one(T))
        q_v = finite_difference(collect(state.Q[1, :]), grid.v)
        logf_v = finite_difference(collect(state.logf[1, :]), grid.v)
        expected_q_v = similar(grid.v)
        expected_logf_v = similar(grid.v)
        for (j, V) in pairs(grid.v)
            sources = initial_pulse_sources(V, ep, one(T), U0, V0, width)
            expected_q_v[j] = sources.q_v
            expected_logf_v[j] = sources.logf_v
        end
        push!(q_residual, maximum(abs, q_v .- expected_q_v))
        push!(logf_residual, maximum(abs, logf_v .- expected_logf_v))
    end
    return q_residual, logf_residual
end

function plot_all(::Type{T}) where {T<:Real}
    output_dir = argument(1, "plots")
    qvalues = q_values_argument(2, "1.0,1.002,1.0033218", T)
    vmax = real_argument(3, "400.0", T)
    dv = real_argument(4, "0.08", T)
    max_rows = integer_argument(5, 6000)
    stride = integer_argument(6, 120)
    amplitude = parse(T, "0.01")
    C = parse(T, "0.6")
    Umax = parse(T, "1.6")
    substep_C = C

    mkpath(output_dir)
    println("# Writing plots to ", abspath(output_dir))

    records_by_q = Dict{T,Vector{NamedTuple}}()
    for q0 in qvalues
        println("# evolving Q0=", q0, " for trapped-surface proxy plots")
        records_by_q[q0] = trapped_tail_records(q0, T; vmax, dv, amplitude, C,
                                                Umax, max_rows, stride,
                                                substep_C)
    end

    function series_for(field)
        return [PlotSeries(
                    [record.V for record in records_by_q[q0]],
                    [getproperty(record, field) for record in records_by_q[q0]],
                    "Q0=$(q0)",
                    COLORS[mod1(index, length(COLORS))],
                )
                for (index, q0) in pairs(qvalues)]
    end

    rv_path = joinpath(output_dir, "gp2026_trapped_proxy_rv.svg")
    h_path = joinpath(output_dir, "gp2026_trapped_proxy_horizon_function.svg")
    rplus_path = joinpath(output_dir, "gp2026_trapped_proxy_rminus_rplus.svg")
    rho_path = joinpath(output_dir, "gp2026_trapped_proxy_maxrho.svg")
    write_svg_line_plot(rv_path, series_for(:r_V);
                        title="GP2026 Row Proxy: Minimum Outgoing Expansion",
                        xlabel="V at row minimum r_V", ylabel="min r_V",
                        xscale=:linear, yscale=:log10)
    write_svg_line_plot(h_path, series_for(:horizon_function);
                        title="GP2026 Row Proxy: Horizon Function",
                        xlabel="V at row minimum r_V", ylabel="H",
                        xscale=:linear, yscale=:log10)
    write_svg_line_plot(rplus_path, series_for(:r_minus_rplus);
                        title="GP2026 Row Proxy: Distance to r_+",
                        xlabel="V at row minimum r_V", ylabel="r - r_+",
                        xscale=:linear, yscale=:log10)
    write_svg_line_plot(rho_path, series_for(:max_rho);
                        title="GP2026 Row Proxy: Throat Coordinate",
                        xlabel="V at row minimum r_V", ylabel="max rho",
                        xscale=:linear, yscale=:linear)

    dvs = T[parse(T, "0.16"), parse(T, "0.08"), parse(T, "0.04")]
    q_residual, logf_residual =
        initial_data_residuals(T; q0=last(qvalues), vmax=parse(T, "80.0"),
                               dvs, amplitude, omega=one(T),
                               width=parse(T, "20.0"))
    residual_series = [
        PlotSeries(collect(dvs), q_residual, "Q_V constraint", COLORS[1]),
        PlotSeries(collect(dvs), logf_residual, "log f constraint", COLORS[2]),
    ]
    residual_path = joinpath(output_dir, "gp2026_initial_residual_convergence.svg")
    write_svg_line_plot(residual_path, residual_series;
                        title="GP2026 Initial Pulse Constraint Residuals",
                        xlabel="Delta V", ylabel="max residual",
                        xscale=:log10, yscale=:log10)

    println("# Plot files:")
    for path in (rv_path, h_path, rplus_path, rho_path, residual_path)
        println(abspath(path))
    end
end

plot_all(Float64)
