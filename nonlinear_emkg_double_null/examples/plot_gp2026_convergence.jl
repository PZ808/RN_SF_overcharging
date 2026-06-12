using NonlinearEMKGDoubleNull
using Printf

const SQRT32PI = sqrt(32 * pi)
const COLORS = (
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf",
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

function real_argument(index, default)
    return parse(Float64, argument(index, string(default)))
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

function log_ticks(lo, hi)
    first_power = ceil(Int, log10(lo))
    last_power = floor(Int, log10(hi))
    ticks = [10.0^power for power in first_power:last_power]
    isempty(ticks) ? [lo, sqrt(lo * hi), hi] : ticks
end

function linear_ticks(lo, hi; count=5)
    [lo + (hi - lo) * (k - 1) / (count - 1) for k in 1:count]
end

function write_svg_line_plot(path, series::Vector{<:PlotSeries};
                             title, xlabel, ylabel,
                             xscale::Symbol=:log10, yscale::Symbol=:log10,
                             width=920, height=560)
    mkpath(dirname(path))
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
                isfinite(x) && isfinite(y) && x > 0 && y > 0 || continue
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

function reference_series(x, y, order, label, color)
    finite = [(x[k], y[k]) for k in eachindex(x) if isfinite(y[k]) && x[k] > 0 && y[k] > 0]
    isempty(finite) && return PlotSeries(Float64[], Float64[], label, color)
    x0, y0 = first(finite)
    return PlotSeries(collect(x), [y0 * (xi / x0)^order for xi in x], label, color)
end

function refinement_rate(previous, current)
    (previous > 0 && current > 0) ? log(previous / current) / log(2) : NaN
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function run_charge_level(; q0, amplitude, U1, V1, du, dv, target_u,
                          target_v, iterations)
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
    initialize_gp2026_single_pulse!(state, grid, ep)
    evolve_nonlinear!(state, grid, ep; iterations, reduced_scalar=true,
                      hyperbolic_charge=true)
    adaptive = adaptive_state_from_rectangular(state, grid)
    _, _, _, _, q_u_residual =
        charged_charge_flux_u_profile(adaptive, ep; target_v, reduced_scalar=true)
    _, _, _, _, q_v_residual =
        charged_charge_flux_v_profile(adaptive, ep; target_u, reduced_scalar=true)
    _, _, _, _, integrated_q_error =
        charged_flux_integrated_charge_profile(adaptive, ep; target_v,
                                               reduced_scalar=true)
    return (
        h=dv,
        du=du,
        dv=dv,
        max_q_u_residual=maximum(abs, q_u_residual),
        max_q_v_residual=maximum(abs, q_v_residual),
        max_integrated_q_error=maximum(abs, integrated_q_error),
    )
end

function interpolate_series(x, y, xq)
    xq < first(x) || xq > last(x) && return nothing
    i = searchsortedlast(x, xq)
    i == length(x) && return y[end]
    i < 1 && return nothing
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - t) * y[i] + t * y[i + 1]
end

function build_throat_series(; q0, rho_match, vmax, dv, C, max_rows, Umax,
                             amplitude)
    ep = EvolutionParams(
        rn=RNParams(1.0, q0),
        scalar_charge=0.6 / q0,
        amplitude=amplitude,
        omega=1.0,
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0=-1.0, V0=0.0, U1=-0.99, V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    evolved = evolve_gp2026_u_adaptive(
        row_from_rectangular(seed, grid, 1), ep;
        Umax, C, iterations=10, max_rows, hyperbolic_charge=true,
        step_control=:outer,
    )
    last_valid = findlast(finite_row, evolved.rows)
    rows = evolved.rows[1:last_valid]
    samples = throat_boundary_series(rows; rho_match, boundary=:outer)
    u = [sample.u for sample in samples]
    return (
        h=dv,
        u=u,
        v=[sample.v for sample in samples],
        r=[sample.r for sample in samples],
        q=[sample.q for sample in samples],
        amplitude=[hypot(sample.phi_re, sample.phi_im) / SQRT32PI
                   for sample in samples],
        q_v_residual=[throat_boundary_observables(sample, ep).q_v_residual
                      for sample in samples],
        rows=length(rows),
        samples=length(samples),
    )
end

function max_difference(coarse, fine, field)
    common_u = [u for u in coarse.u if first(fine.u) <= u <= last(fine.u)]
    isempty(common_u) && return NaN
    coarse_values = getproperty(coarse, field)
    fine_values = getproperty(fine, field)
    differences = Float64[]
    for u in common_u
        coarse_value = interpolate_series(coarse.u, coarse_values, u)
        fine_value = interpolate_series(fine.u, fine_values, u)
        if !isnothing(coarse_value) && !isnothing(fine_value)
            push!(differences, abs(coarse_value - fine_value))
        end
    end
    isempty(differences) ? NaN : maximum(differences)
end

function plot_convergence()
    output_dir = argument(1, "plots")
    q0 = real_argument(2, 1.0033218)
    amplitude = real_argument(3, 0.01)
    charge_levels = integer_argument(4, 4)
    throat_levels = integer_argument(5, 3)
    mkpath(output_dir)

    println("# Charge residual convergence")
    charge_rows = [
        run_charge_level(; q0, amplitude, U1=-0.8, V1=20.0,
                         du=0.01 / 2.0^level,
                         dv=0.08 / 2.0^level,
                         target_u=-0.9, target_v=10.0,
                         iterations=10)
        for level in 0:charge_levels-1
    ]
    for row in charge_rows
        println(row)
    end
    h_charge = [row.h for row in charge_rows]
    charge_series = [
        PlotSeries(h_charge, [row.max_q_u_residual for row in charge_rows],
                   "Q_U residual", COLORS[1]),
        PlotSeries(h_charge, [row.max_q_v_residual for row in charge_rows],
                   "Q_V residual", COLORS[2]),
        PlotSeries(h_charge, [row.max_integrated_q_error for row in charge_rows],
                   "integrated Q", COLORS[3]),
        reference_series(h_charge, [row.max_q_v_residual for row in charge_rows],
                         2, "O(h^2)", COLORS[6]),
    ]
    charge_path = joinpath(output_dir, "gp2026_charge_residual_convergence.svg")
    write_svg_line_plot(charge_path, charge_series;
                        title="GP2026 Evolved Charge Constraint Convergence",
                        xlabel="Delta V", ylabel="max residual")

    println("# Fixed-rho throat-boundary convergence")
    throat_levels_data = [
        build_throat_series(; q0, rho_match=2.0, vmax=400.0,
                            dv=0.08 / 2.0^level,
                            C=0.6 / 2.0^level,
                            max_rows=80 * 2^level,
                            Umax=1.6, amplitude)
        for level in 0:throat_levels-1
    ]
    for row in throat_levels_data
        println((h=row.h, rows=row.rows, samples=row.samples,
                 max_qv_residual=maximum(abs, row.q_v_residual)))
    end
    pair_h = [throat_levels_data[i].h for i in 1:length(throat_levels_data)-1]
    d_v = [max_difference(throat_levels_data[i], throat_levels_data[i + 1], :v)
           for i in 1:length(throat_levels_data)-1]
    d_r = [max_difference(throat_levels_data[i], throat_levels_data[i + 1], :r)
           for i in 1:length(throat_levels_data)-1]
    d_q = [max_difference(throat_levels_data[i], throat_levels_data[i + 1], :q)
           for i in 1:length(throat_levels_data)-1]
    d_amp = [max_difference(throat_levels_data[i], throat_levels_data[i + 1],
                            :amplitude)
             for i in 1:length(throat_levels_data)-1]
    throat_series = [
        PlotSeries(pair_h, d_v, "V(rho=2)", COLORS[1]),
        PlotSeries(pair_h, d_r, "r(rho=2)", COLORS[2]),
        PlotSeries(pair_h, d_q, "Q(rho=2)", COLORS[3]),
        PlotSeries(pair_h, d_amp, "|r phi|(rho=2)", COLORS[4]),
        reference_series(pair_h, d_r, 2, "O(h^2)", COLORS[6]),
    ]
    throat_path = joinpath(output_dir, "gp2026_throat_boundary_convergence.svg")
    write_svg_line_plot(throat_path, throat_series;
                        title="GP2026 Fixed-rho Throat Boundary Convergence",
                        xlabel="coarse Delta V", ylabel="coarse-fine max difference")

    println("# Plot files:")
    println(abspath(charge_path))
    println(abspath(throat_path))
end

plot_convergence()
