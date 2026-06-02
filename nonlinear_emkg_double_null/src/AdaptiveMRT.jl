const NL_FIELD_NAMES = (:r, :logf, :phi_re, :phi_im, :Au, :Av, :Q)

"""
One nonlinear MRT slice at fixed V.

Unlike `NLState`, this stores its own U grid. This is the storage unit needed
for Burko-Ori/MRT refinement, where the number of U points can change from one
V slice to the next while the V step remains fixed.
"""
struct NLSlice{T<:Real}
    v::T
    u::Vector{T}
    r::Vector{T}
    logf::Vector{T}
    phi_re::Vector{T}
    phi_im::Vector{T}
    Au::Vector{T}
    Av::Vector{T}
    Q::Vector{T}
end

struct AdaptiveNLState{T<:Real}
    slices::Vector{NLSlice{T}}
end

"""
One nonlinear GP2026 row at fixed `U`.

This storage orientation supports the production rule
`Delta U = C/f_GP(U,Vmax)`, which selects the next outgoing hypersurface
after the current row has been evolved to the outer `V` boundary.
"""
struct NLRow{T<:Real}
    u::T
    v::Vector{T}
    r::Vector{T}
    logf::Vector{T}
    phi_re::Vector{T}
    phi_im::Vector{T}
    Au::Vector{T}
    Av::Vector{T}
    Q::Vector{T}
end

struct UAdaptiveNLState{T<:Real}
    rows::Vector{NLRow{T}}
end

struct ThroatRowDiagnostics{T<:Real}
    y::Vector{T}
    rho::Vector{T}
    qabs::Vector{T}
    min_y::T
    max_rho::T
    max_abs_delta_rho::T
end

struct ThroatMatchCandidate{T<:Real}
    index::Int
    u::T
    v::T
    r::T
    q::T
    y::T
    rho::T
end

struct ThroatMatchBand{T<:Real}
    first_index::Int
    last_index::Int
    count::Int
    component_count::Int
    u::T
    v_first::T
    v_last::T
    rho_minimum_in_band::T
    rho_maximum_in_band::T
end

struct RhoLapseDiagnostics{T<:Real}
    rho::Vector{T}
    rho_v::Vector{T}
    logf_rho::Vector{T}
    logf_range::Tuple{T,T}
    logf_rho_range::Tuple{T,T}
    throat_logf_range::Tuple{T,T}
    throat_logf_rho_range::Tuple{T,T}
    throat_count::Int
end

struct NLPoint{T<:Real}
    u::T
    v::T
    r::T
    logf::T
    phi_re::T
    phi_im::T
    Au::T
    Av::T
    Q::T
end

struct PointSplittingConfig{T<:Real}
    band_width::T
    max_relative_r::T
    max_relative_f::T
    max_dphi::T
    max_dphiu::T
end

struct HorizonChoppingConfig{T<:Real}
    start_v::T
    band_width::T
    interior_buffer_cells::Int
    interior_buffer_width::T
end

struct HorizonRefinementConfig{T<:Real}
    start_v::T
    band_width::T
    max_du::T
    exterior_cells::Int
    interior_cells::Int
end

function PointSplittingConfig(; band_width=1.0, max_relative_r=Inf, max_relative_f=Inf,
                              max_dphi=Inf, max_dphiu=Inf)
    band_width > 0 || throw(ArgumentError("band_width must be positive"))
    max_relative_r > 0 || throw(ArgumentError("max_relative_r must be positive"))
    max_relative_f > 0 || throw(ArgumentError("max_relative_f must be positive"))
    max_dphi > 0 || throw(ArgumentError("max_dphi must be positive"))
    max_dphiu > 0 || throw(ArgumentError("max_dphiu must be positive"))
    T = promote_type(typeof(band_width), typeof(max_relative_r), typeof(max_relative_f),
                     typeof(max_dphi), typeof(max_dphiu))
    return PointSplittingConfig{T}(convert(T, band_width), convert(T, max_relative_r),
                                   convert(T, max_relative_f), convert(T, max_dphi),
                                   convert(T, max_dphiu))
end

function HorizonChoppingConfig(; start_v=0.0, band_width=1.0, interior_buffer_cells::Int=1,
                               interior_buffer_width=0.0)
    start_v >= 0 || throw(ArgumentError("start_v must be nonnegative"))
    band_width > 0 || throw(ArgumentError("band_width must be positive"))
    interior_buffer_cells >= 0 ||
        throw(ArgumentError("interior_buffer_cells must be nonnegative"))
    interior_buffer_width >= 0 ||
        throw(ArgumentError("interior_buffer_width must be nonnegative"))
    T = promote_type(typeof(start_v), typeof(band_width), typeof(interior_buffer_width))
    return HorizonChoppingConfig{T}(convert(T, start_v), convert(T, band_width),
                                    interior_buffer_cells, convert(T, interior_buffer_width))
end

function HorizonRefinementConfig(; start_v=0.0, band_width=1.0, max_du,
                                 exterior_cells::Int=2, interior_cells::Int=1)
    start_v >= 0 || throw(ArgumentError("start_v must be nonnegative"))
    band_width > 0 || throw(ArgumentError("band_width must be positive"))
    max_du > 0 || throw(ArgumentError("max_du must be positive"))
    exterior_cells >= 0 || throw(ArgumentError("exterior_cells must be nonnegative"))
    interior_cells >= 0 || throw(ArgumentError("interior_cells must be nonnegative"))
    T = promote_type(typeof(start_v), typeof(band_width), typeof(max_du))
    return HorizonRefinementConfig{T}(convert(T, start_v), convert(T, band_width),
                                      convert(T, max_du), exterior_cells, interior_cells)
end

function NLPoint(u::Real, v::Real, r::Real, logf::Real, phi_re::Real, phi_im::Real,
                 Au::Real, Av::Real, Q::Real)
    T = promote_type(typeof(u), typeof(v), typeof(r), typeof(logf), typeof(phi_re),
                     typeof(phi_im), typeof(Au), typeof(Av), typeof(Q))
    return NLPoint{T}(convert(T, u), convert(T, v), convert(T, r),
                      convert(T, logf), convert(T, phi_re), convert(T, phi_im),
                      convert(T, Au), convert(T, Av), convert(T, Q))
end

function require_strictly_increasing(x::AbstractVector, name::AbstractString)
    length(x) >= 2 || throw(ArgumentError("$name must have at least two points"))
    for i in firstindex(x):lastindex(x)-1
        x[i] < x[i + 1] || throw(ArgumentError("$name must be strictly increasing"))
    end
    return x
end

function require_same_length(name::AbstractString, n::Int, values::AbstractVector...)
    for value in values
        length(value) == n || throw(ArgumentError("$name fields must all have length $n"))
    end
    return n
end

function NLSlice(v::Real, u::AbstractVector{<:Real}, r::AbstractVector{<:Real},
                 logf::AbstractVector{<:Real}, phi_re::AbstractVector{<:Real},
                 phi_im::AbstractVector{<:Real}, Au::AbstractVector{<:Real},
                 Av::AbstractVector{<:Real}, Q::AbstractVector{<:Real})
    n = length(u)
    require_same_length("NLSlice", n, r, logf, phi_re, phi_im, Au, Av, Q)
    require_strictly_increasing(u, "u")
    T = promote_type(typeof(v), eltype(u), eltype(r), eltype(logf), eltype(phi_re),
                     eltype(phi_im), eltype(Au), eltype(Av), eltype(Q))
    return NLSlice{T}(convert(T, v), collect(T, u), collect(T, r), collect(T, logf),
                      collect(T, phi_re), collect(T, phi_im), collect(T, Au),
                      collect(T, Av), collect(T, Q))
end

function NLSlice(v::Real, u::AbstractVector{<:Real})
    require_strictly_increasing(u, "u")
    T = promote_type(typeof(v), eltype(u))
    uu = collect(T, u)
    z = zeros(T, length(uu))
    return NLSlice{T}(convert(T, v), uu, copy(z), copy(z), copy(z), copy(z),
                      copy(z), copy(z), copy(z))
end

function AdaptiveNLState(slices::AbstractVector{<:NLSlice})
    isempty(slices) && throw(ArgumentError("AdaptiveNLState needs at least one slice"))
    T = promote_type(map(s -> typeof(s.v), slices)...)
    converted = NLSlice{T}[]
    for slice in slices
        push!(converted, NLSlice(slice.v, slice.u, slice.r, slice.logf, slice.phi_re,
                                slice.phi_im, slice.Au, slice.Av, slice.Q))
    end
    return AdaptiveNLState{T}(converted)
end

function NLRow(u::Real, v::AbstractVector{<:Real}, r::AbstractVector{<:Real},
               logf::AbstractVector{<:Real}, phi_re::AbstractVector{<:Real},
               phi_im::AbstractVector{<:Real}, Au::AbstractVector{<:Real},
               Av::AbstractVector{<:Real}, Q::AbstractVector{<:Real})
    n = length(v)
    require_same_length("NLRow", n, r, logf, phi_re, phi_im, Au, Av, Q)
    require_strictly_increasing(v, "v")
    T = promote_type(typeof(u), eltype(v), eltype(r), eltype(logf), eltype(phi_re),
                     eltype(phi_im), eltype(Au), eltype(Av), eltype(Q))
    return NLRow{T}(convert(T, u), collect(T, v), collect(T, r), collect(T, logf),
                    collect(T, phi_re), collect(T, phi_im), collect(T, Au),
                    collect(T, Av), collect(T, Q))
end

function UAdaptiveNLState(rows::AbstractVector{<:NLRow})
    isempty(rows) && throw(ArgumentError("UAdaptiveNLState needs at least one row"))
    T = promote_type(map(row -> typeof(row.u), rows)...)
    converted = NLRow{T}[]
    for row in rows
        push!(converted, NLRow(row.u, row.v, row.r, row.logf, row.phi_re,
                               row.phi_im, row.Au, row.Av, row.Q))
    end
    return UAdaptiveNLState{T}(converted)
end

function slice_point(slice::NLSlice, i::Int)
    1 <= i <= length(slice.u) || throw(BoundsError(slice.u, i))
    return NLPoint(slice.u[i], slice.v, slice.r[i], slice.logf[i], slice.phi_re[i],
                   slice.phi_im[i], slice.Au[i], slice.Av[i], slice.Q[i])
end

function slice_field(slice::NLSlice, field::Symbol)
    field === :r && return slice.r
    field === :logf && return slice.logf
    field === :phi_re && return slice.phi_re
    field === :phi_im && return slice.phi_im
    field === :Au && return slice.Au
    field === :Av && return slice.Av
    field === :Q && return slice.Q
    throw(ArgumentError("unknown nonlinear field: $field"))
end

function slice_from_rectangular(st::NLState, g::Grid, j::Int)
    1 <= j <= length(g.v) || throw(BoundsError(g.v, j))
    return NLSlice(g.v[j], g.u, st.r[:, j], st.logf[:, j], st.phi_re[:, j],
                   st.phi_im[:, j], st.Au[:, j], st.Av[:, j], st.Q[:, j])
end

function row_from_rectangular(st::NLState, g::Grid, i::Int)
    1 <= i <= length(g.u) || throw(BoundsError(g.u, i))
    return NLRow(g.u[i], g.v, st.r[i, :], st.logf[i, :], st.phi_re[i, :],
                 st.phi_im[i, :], st.Au[i, :], st.Av[i, :], st.Q[i, :])
end

function adaptive_state_from_rectangular(st::NLState, g::Grid)
    return AdaptiveNLState([slice_from_rectangular(st, g, j) for j in eachindex(g.v)])
end

function adaptive_state_from_u_rows(st::UAdaptiveNLState)
    rows = st.rows
    reference_v = first(rows).v
    for row in rows
        row.v == reference_v ||
            throw(ArgumentError("all U-adaptive rows must share the same V grid"))
    end
    u = [row.u for row in rows]
    return AdaptiveNLState([
        NLSlice(reference_v[j], u, [row.r[j] for row in rows],
                [row.logf[j] for row in rows], [row.phi_re[j] for row in rows],
                [row.phi_im[j] for row in rows], [row.Au[j] for row in rows],
                [row.Av[j] for row in rows], [row.Q[j] for row in rows])
        for j in eachindex(reference_v)
    ])
end

function throat_row_diagnostics(row::NLRow; charge_floor=nothing,
                                throat_floor=nothing)
    T = promote_type(eltype(row.r), eltype(row.Q))
    qfloor = isnothing(charge_floor) ? sqrt(eps(float(one(T)))) :
             convert(T, charge_floor)
    yfloor = isnothing(throat_floor) ? sqrt(eps(float(one(T)))) :
             convert(T, throat_floor)
    qabs = [max(abs(q), qfloor) for q in row.Q]
    y = [max(row.r[j] - qabs[j], yfloor * qabs[j]) for j in eachindex(row.r)]
    rho = [-log(y[j] / qabs[j]) for j in eachindex(y)]
    max_abs_delta_rho = length(rho) <= 1 ? zero(eltype(rho)) : maximum(abs, diff(rho))
    return ThroatRowDiagnostics(y, rho, qabs, minimum(y), maximum(rho),
                                max_abs_delta_rho)
end

function throat_matching_candidate(row::NLRow; rho_min=2.0, charge_floor=nothing,
                                   throat_floor=nothing)
    diagnostics = throat_row_diagnostics(row; charge_floor, throat_floor)
    index = findfirst(rho -> rho >= rho_min, diagnostics.rho)
    isnothing(index) && return nothing
    return ThroatMatchCandidate(index, row.u, row.v[index], row.r[index],
                                row.Q[index], diagnostics.y[index],
                                diagnostics.rho[index])
end

function throat_matching_band(row::NLRow; rho_min=2.0, charge_floor=nothing,
                              throat_floor=nothing)
    diagnostics = throat_row_diagnostics(row; charge_floor, throat_floor)
    indices = findall(rho -> rho >= rho_min, diagnostics.rho)
    isempty(indices) && return nothing
    component_count = one(Int)
    for k in firstindex(indices)+1:lastindex(indices)
        component_count += indices[k] == indices[k - 1] + 1 ? 0 : 1
    end
    first_index = first(indices)
    last_index = last(indices)
    band_rho = diagnostics.rho[indices]
    return ThroatMatchBand(first_index, last_index, length(indices),
                           component_count, row.u,
                           row.v[first_index], row.v[last_index],
                           minimum(band_rho), maximum(band_rho))
end

function finite_extrema(values::AbstractVector{T}) where {T<:Real}
    finite = [value for value in values if isfinite(value)]
    isempty(finite) && return (T(NaN), T(NaN))
    return extrema(finite)
end

range_width(range::Tuple{<:Real,<:Real}) = range[2] - range[1]

function coordinate_derivative(values::AbstractVector{T},
                               coordinates::AbstractVector{T}) where {T<:Real}
    length(values) == length(coordinates) ||
        throw(ArgumentError("derivative arrays must have equal length"))
    length(values) >= 2 || throw(ArgumentError("at least two points are required"))
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

function rho_lapse_diagnostics(row::NLRow; rho_min=2.0, charge_floor=nothing,
                               throat_floor=nothing)
    throat = throat_row_diagnostics(row; charge_floor, throat_floor)
    rho_v = coordinate_derivative(throat.rho, row.v)
    logf_rho = [isfinite(rho_v[j]) && rho_v[j] != 0 ?
                row.logf[j] - log(abs(rho_v[j])) :
                typeof(row.logf[j])(Inf)
                for j in eachindex(row.logf)]
    throat_indices = findall(rho -> rho >= rho_min, throat.rho)
    throat_logf_range = isempty(throat_indices) ?
                         finite_extrema(eltype(row.logf)[]) :
                         finite_extrema(row.logf[throat_indices])
    throat_logf_rho_range = isempty(throat_indices) ?
                             finite_extrema(eltype(logf_rho)[]) :
                             finite_extrema(logf_rho[throat_indices])
    return RhoLapseDiagnostics(
        throat.rho,
        rho_v,
        logf_rho,
        finite_extrema(row.logf),
        finite_extrema(logf_rho),
        throat_logf_range,
        throat_logf_rho_range,
        length(throat_indices),
    )
end

function throat_row_du(rows::AbstractVector{<:NLRow}, C::Real; max_delta_rho=0.25)
    T = promote_type(typeof(C), typeof(first(rows).u), typeof(max_delta_rho))
    length(rows) >= 2 || return T(Inf)
    current = rows[end]
    previous = rows[end - 1]
    current.v == previous.v ||
        throw(ArgumentError("throat step requires matching V grids"))
    du_previous = current.u - previous.u
    du_previous > 0 || throw(ArgumentError("row U coordinates must increase"))
    current_rho = throat_row_diagnostics(current).rho
    previous_rho = throat_row_diagnostics(previous).rho
    rho_u = [(current_rho[j] - previous_rho[j]) / du_previous
             for j in eachindex(current_rho)]
    max_speed = maximum(abs, rho_u)
    isfinite(max_speed) && max_speed > 0 || return T(Inf)
    return convert(T, max_delta_rho) / max_speed
end

function west_boundary_from_rectangular(st::NLState, g::Grid)
    i = firstindex(g.u)
    return [NLPoint(g.u[i], g.v[j], st.r[i, j], st.logf[i, j], st.phi_re[i, j],
                    st.phi_im[i, j], st.Au[i, j], st.Av[i, j], st.Q[i, j])
            for j in eachindex(g.v)]
end

function interpolate_linear(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))
    require_strictly_increasing(x, "x")
    scale = max(abs(first(x)), abs(last(x)), abs(xq), one(float(xq)))
    tol = 32 * eps(float(scale)) * scale
    xq < first(x) - tol && throw(ArgumentError("query point is below interpolation range"))
    xq > last(x) + tol && throw(ArgumentError("query point is above interpolation range"))
    xq <= first(x) && return y[firstindex(y)]
    xq >= last(x) && return y[lastindex(y)]

    i = searchsortedlast(x, xq)
    x[i] == xq && return y[i]
    i == lastindex(x) && return y[lastindex(y)]
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - t) * y[i] + t * y[i + 1]
end

function interpolate_values(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                            xq::AbstractVector{<:Real})
    return [interpolate_linear(x, y, xqi) for xqi in xq]
end

function interpolate_slice(slice::NLSlice, new_u::AbstractVector{<:Real})
    require_strictly_increasing(new_u, "new_u")
    return NLSlice(slice.v, new_u,
                   interpolate_values(slice.u, slice.r, new_u),
                   interpolate_values(slice.u, slice.logf, new_u),
                   interpolate_values(slice.u, slice.phi_re, new_u),
                   interpolate_values(slice.u, slice.phi_im, new_u),
                   interpolate_values(slice.u, slice.Au, new_u),
                   interpolate_values(slice.u, slice.Av, new_u),
                   interpolate_values(slice.u, slice.Q, new_u))
end

function refine_u_grid(u::AbstractVector{<:Real}, refine_cells::AbstractVector{Bool})
    require_strictly_increasing(u, "u")
    length(refine_cells) == length(u) - 1 ||
        throw(ArgumentError("refine_cells must have length length(u)-1"))
    T = eltype(u)
    refined = T[]
    sizehint!(refined, length(u) + count(refine_cells))
    push!(refined, first(u))
    for i in firstindex(refine_cells):lastindex(refine_cells)
        if refine_cells[i]
            push!(refined, (u[i] + u[i + 1]) / 2)
        end
        push!(refined, u[i + 1])
    end
    return refined
end

function local_polynomial_value(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                                xq::Real, cell::Int)
    n = length(x)
    length(y) == n || throw(ArgumentError("x and y must have the same length"))
    width = min(4, n)
    first_stencil = clamp(cell - 1, 1, n - width + 1)
    value = zero(promote_type(eltype(y), typeof(xq)))
    for a in first_stencil:first_stencil+width-1
        basis = one(value)
        for b in first_stencil:first_stencil+width-1
            a == b && continue
            basis *= (xq - x[b]) / (x[a] - x[b])
        end
        value += basis * y[a]
    end
    return value
end

function local_polynomial_derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                                     xq::Real, cell::Int)
    n = length(x)
    length(y) == n || throw(ArgumentError("x and y must have the same length"))
    width = min(4, n)
    first_stencil = clamp(cell - 1, 1, n - width + 1)
    derivative = zero(promote_type(eltype(y), typeof(xq)))
    for a in first_stencil:first_stencil+width-1
        basis_derivative = zero(derivative)
        for c in first_stencil:first_stencil+width-1
            c == a && continue
            term = inv(x[a] - x[c])
            for b in first_stencil:first_stencil+width-1
                (b == a || b == c) && continue
                term *= (xq - x[b]) / (x[a] - x[b])
            end
            basis_derivative += term
        end
        derivative += basis_derivative * y[a]
    end
    return derivative
end

function split_interpolate_values(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                                  refine_cells::AbstractVector{Bool})
    length(refine_cells) == length(x) - 1 ||
        throw(ArgumentError("refine_cells must have length length(x)-1"))
    T = promote_type(eltype(x), eltype(y))
    values = T[]
    sizehint!(values, length(y) + count(refine_cells))
    push!(values, y[1])
    for i in eachindex(refine_cells)
        if refine_cells[i]
            midpoint = (x[i] + x[i + 1]) / 2
            push!(values, local_polynomial_value(x, y, midpoint, i))
        end
        push!(values, y[i + 1])
    end
    return values
end

function spacing_refinement_flags(u::AbstractVector{<:Real}; max_du::Real)
    require_strictly_increasing(u, "u")
    return [u[i + 1] - u[i] > max_du for i in firstindex(u):lastindex(u)-1]
end

function variation_refinement_flags(slice::NLSlice; max_dr=Inf, max_dlogf=Inf,
                                    max_dphi=Inf, max_dA=Inf, max_dQ=Inf)
    ncell = length(slice.u) - 1
    flags = falses(ncell)
    thresholds = (
        (:r, max_dr),
        (:logf, max_dlogf),
        (:phi_re, max_dphi),
        (:phi_im, max_dphi),
        (:Au, max_dA),
        (:Av, max_dA),
        (:Q, max_dQ),
    )
    for (field, threshold) in thresholds
        isfinite(threshold) || continue
        values = slice_field(slice, field)
        for i in 1:ncell
            flags[i] |= abs(values[i + 1] - values[i]) > threshold
        end
    end
    return flags
end

function point_splitting_flags(slice::NLSlice, config::PointSplittingConfig)
    ncell = length(slice.u) - 1
    flags = falses(ncell)
    for i in 1:ncell
        relative_r = abs(slice.r[i + 1] - slice.r[i]) /
                     max(abs(slice.r[i]), eps(float(slice.r[i])))
        # This is |(f[i+1]-f[i])/f[i]| without exponentiating f itself.
        relative_f = abs(expm1(slice.logf[i + 1] - slice.logf[i]))
        dphi = hypot(slice.phi_re[i + 1] - slice.phi_re[i],
                     slice.phi_im[i + 1] - slice.phi_im[i])
        flags[i] = relative_r > config.max_relative_r ||
                   relative_f > config.max_relative_f ||
                   dphi > config.max_dphi
    end

    if isfinite(config.max_dphiu) && ncell > 1
        phiu_re = [(slice.phi_re[i + 1] - slice.phi_re[i]) /
                   (slice.u[i + 1] - slice.u[i]) for i in 1:ncell]
        phiu_im = [(slice.phi_im[i + 1] - slice.phi_im[i]) /
                   (slice.u[i + 1] - slice.u[i]) for i in 1:ncell]
        for i in 1:ncell-1
            dphiu = hypot(phiu_re[i + 1] - phiu_re[i],
                          phiu_im[i + 1] - phiu_im[i])
            if dphiu > config.max_dphiu
                flags[i] = true
                flags[i + 1] = true
            end
        end
    end
    return flags
end

function refine_slice(slice::NLSlice, refine_cells::AbstractVector{Bool})
    refined_u = refine_u_grid(slice.u, refine_cells)
    return NLSlice(slice.v, refined_u,
                   split_interpolate_values(slice.u, slice.r, refine_cells),
                   split_interpolate_values(slice.u, slice.logf, refine_cells),
                   split_interpolate_values(slice.u, slice.phi_re, refine_cells),
                   split_interpolate_values(slice.u, slice.phi_im, refine_cells),
                   split_interpolate_values(slice.u, slice.Au, refine_cells),
                   split_interpolate_values(slice.u, slice.Av, refine_cells),
                   split_interpolate_values(slice.u, slice.Q, refine_cells))
end

"""
Evaluate the `U` Raychaudhuri constraint integrand inside one parent cell.

For `y = r_U / f`, the convention used by the nonlinear evolution gives

    y_U = -r T_UU / (8f),

and equivalently

    (log f)_U = (log |r_U|)_U + r T_UU / (8r_U).
"""
function u_constraint_data(slice::NLSlice, ep::EvolutionParams, u::Real, cell::Int;
                           reduced_scalar::Bool=false)
    r = local_polynomial_value(slice.u, slice.r, u, cell)
    logf = local_polynomial_value(slice.u, slice.logf, u, cell)
    phi_re = local_polynomial_value(slice.u, slice.phi_re, u, cell)
    phi_im = local_polynomial_value(slice.u, slice.phi_im, u, cell)
    Au = local_polynomial_value(slice.u, slice.Au, u, cell)
    Av = local_polynomial_value(slice.u, slice.Av, u, cell)
    q = local_polynomial_value(slice.u, slice.Q, u, cell)
    ru = local_polynomial_derivative(slice.u, slice.r, u, cell)
    phiu_re = local_polynomial_derivative(slice.u, slice.phi_re, u, cell)
    phiu_im = local_polynomial_derivative(slice.u, slice.phi_im, u, cell)
    scale = max(abs(ru), one(float(ru)))
    threshold = 64 * eps(float(scale)) * scale
    abs(ru) > threshold ||
        throw(ArgumentError("constraint-preserving insertion cannot cross r_U = 0"))
    source = if reduced_scalar
        stress_energy_reduced_scalar(r, exp(logf), q, ru, zero(ru),
                                     phi_re, phi_im, phiu_re, zero(phiu_re),
                                     phiu_im, zero(phiu_im), Au, Av,
                                     ep.scalar_charge)
    else
        stress_energy(r, exp(logf), q, phi_re, phi_im,
                      phiu_re, zero(phiu_re), phiu_im, zero(phiu_im),
                      Au, Av, ep.scalar_charge)
    end
    return ru, r * source.Tuu / (8 * ru)
end

"""
Refine a fixed-`V` slice while preserving its evolved coarse-grid data.

All fields except `logf` use four-point midpoint interpolation. For each
inserted midpoint, `logf` is instead reconstructed from the two parent
endpoints using the `U` Raychaudhuri constraint. Existing coarse-grid points
are left unchanged, avoiding a projection error accumulated across the slice.
"""
function refine_slice_constrained(slice::NLSlice, refine_cells::AbstractVector{Bool},
                                  ep::EvolutionParams; reduced_scalar::Bool=false)
    any(refine_cells) || return slice
    refined = refine_slice(slice, refine_cells)
    logf = copy(refined.logf)
    refined_index = 1
    for i in eachindex(refine_cells)
        if refine_cells[i]
            left_u = slice.u[i]
            right_u = slice.u[i + 1]
            mid_u = (left_u + right_u) / 2
            half_du = mid_u - left_u
            ru_left, source_left = u_constraint_data(slice, ep, left_u, i; reduced_scalar)
            ru_mid, source_mid = u_constraint_data(slice, ep, mid_u, i; reduced_scalar)
            ru_right, source_right = u_constraint_data(slice, ep, right_u, i; reduced_scalar)
            sign(ru_left) == sign(ru_mid) == sign(ru_right) ||
                throw(ArgumentError("constraint-preserving insertion requires one sign of r_U"))
            from_left = slice.logf[i] + log(abs(ru_mid / ru_left)) +
                        half_du * (source_left + source_mid) / 2
            from_right = slice.logf[i + 1] + log(abs(ru_mid / ru_right)) -
                         half_du * (source_mid + source_right) / 2
            logf[refined_index + 1] = (from_left + from_right) / 2
            refined_index += 1
        end
        refined_index += 1
    end
    return NLSlice(refined.v, refined.u, refined.r, logf, refined.phi_re, refined.phi_im,
                   refined.Au, refined.Av, refined.Q)
end

function truncate_slice(slice::NLSlice, last_point::Int)
    2 <= last_point <= length(slice.u) ||
        throw(ArgumentError("last_point must retain between two and all slice points"))
    return NLSlice(slice.v, slice.u[1:last_point], slice.r[1:last_point],
                   slice.logf[1:last_point], slice.phi_re[1:last_point],
                   slice.phi_im[1:last_point], slice.Au[1:last_point],
                   slice.Av[1:last_point], slice.Q[1:last_point])
end

function adaptive_outgoing_expansion(previous::NLSlice, current::NLSlice)
    current.v > previous.v ||
        throw(ArgumentError("current slice must have larger V than previous"))
    previous_on_current = previous.u == current.u ? previous :
                          interpolate_slice(previous, current.u)
    dv = current.v - previous.v
    return [((current.r[i] - previous_on_current.r[i]) +
             (current.r[i + 1] - previous_on_current.r[i + 1])) / (2dv)
            for i in 1:length(current.u)-1]
end

function chop_inside_apparent_horizon(previous::NLSlice, current::NLSlice,
                                      config::HorizonChoppingConfig)
    rv = adaptive_outgoing_expansion(previous, current)
    crossing = findfirst(value -> value <= 0, rv)
    isnothing(crossing) && return current
    last_point = min(length(current.u), crossing + 1 + config.interior_buffer_cells)
    if config.interior_buffer_width > 0
        target_u = current.u[crossing + 1] + config.interior_buffer_width
        last_point = max(last_point, min(length(current.u),
                                         searchsortedfirst(current.u, target_u)))
    end
    return last_point == length(current.u) ? current : truncate_slice(current, last_point)
end

function horizon_refinement_flags(previous::NLSlice, current::NLSlice,
                                  config::HorizonRefinementConfig)
    rv = adaptive_outgoing_expansion(previous, current)
    flags = falses(length(current.u) - 1)
    crossing = findfirst(value -> value <= 0, rv)
    isnothing(crossing) && return flags

    lo = max(firstindex(flags), crossing - config.exterior_cells)
    hi = min(lastindex(flags), crossing + config.interior_cells)
    for i in lo:hi
        flags[i] = current.u[i + 1] - current.u[i] > config.max_du
    end
    return flags
end

function refine_near_apparent_horizon(previous::NLSlice, current::NLSlice,
                                      config::HorizonRefinementConfig)
    flags = horizon_refinement_flags(previous, current, config)
    return any(flags) ? refine_slice(current, flags) : current
end

function refine_near_apparent_horizon(previous::NLSlice, current::NLSlice,
                                      config::HorizonRefinementConfig,
                                      ep::EvolutionParams; reduced_scalar::Bool=false)
    flags = horizon_refinement_flags(previous, current, config)
    return any(flags) ? refine_slice_constrained(current, flags, ep; reduced_scalar) : current
end

function require_same_coordinate(a::Real, b::Real, name::AbstractString)
    T = promote_type(typeof(a), typeof(b))
    aa = convert(T, a)
    bb = convert(T, b)
    scale = max(abs(aa), abs(bb), one(float(aa)))
    abs(aa - bb) <= 64 * eps(float(scale)) * scale ||
        throw(ArgumentError("$name mismatch: $aa != $bb"))
    return aa
end

function fill_first_column!(st::NLState, slice::NLSlice)
    st.r[:, 1] .= slice.r
    st.logf[:, 1] .= slice.logf
    st.phi_re[:, 1] .= slice.phi_re
    st.phi_im[:, 1] .= slice.phi_im
    st.Au[:, 1] .= slice.Au
    st.Av[:, 1] .= slice.Av
    st.Q[:, 1] .= slice.Q
    return st
end

function fill_northwest_boundary!(st::NLState, point::NLPoint)
    st.r[1, 2] = point.r
    st.logf[1, 2] = point.logf
    st.phi_re[1, 2] = point.phi_re
    st.phi_im[1, 2] = point.phi_im
    st.Au[1, 2] = point.Au
    st.Av[1, 2] = point.Av
    st.Q[1, 2] = point.Q
    return st
end

function gp2026_na_boundary_point(U::Real, ep::EvolutionParams; U0=-1.0, V0=0.0,
                                  M0=ep.rn.M)
    U >= U0 || throw(ArgumentError("GP2026 N_A boundary requires U >= U0"))
    r0 = gp2026_extremal_gauge_initial_radius(U0, V0; U0, V0, M0)
    r = gp2026_extremal_gauge_initial_radius(U, V0; U0, V0, M0)
    fcorner = gp2026_fcorner_code(ep; U0, V0, M0)
    # Along N_A the scalar vanishes and f_code is constant. Integrating
    # A_V,U=-Q*f_code/(4r^2), with r_U=-1/2, gives this exact boundary value.
    Av = -ep.rn.Q0 * fcorner / 2 * (inv(r) - inv(r0))
    return NLPoint(U, V0, r, log(fcorner), zero(r), zero(r), zero(r), Av, ep.rn.Q0)
end

function advance_u_row(previous::NLRow, south::NLPoint, ep::EvolutionParams;
                       iterations::Int=5, subtract_rn_background::Bool=false,
                       reduced_scalar::Bool=false, hyperbolic_charge::Bool=false)
    south.u > previous.u ||
        throw(ArgumentError("south.u must be larger than the previous row U"))
    require_same_coordinate(first(previous.v), south.v, "south V")
    grid = Grid([previous.u, south.u], previous.v)
    st = NLState(grid)
    st.r[1, :] .= previous.r
    st.logf[1, :] .= previous.logf
    st.phi_re[1, :] .= previous.phi_re
    st.phi_im[1, :] .= previous.phi_im
    st.Au[1, :] .= previous.Au
    st.Av[1, :] .= previous.Av
    st.Q[1, :] .= previous.Q
    st.r[2, 1] = south.r
    st.logf[2, 1] = south.logf
    st.phi_re[2, 1] = south.phi_re
    st.phi_im[2, 1] = south.phi_im
    st.Au[2, 1] = south.Au
    st.Av[2, 1] = south.Av
    st.Q[2, 1] = south.Q
    for j in 1:length(previous.v)-1
        step_nonlinear_cell!(st, grid, ep, 1, j;
                             iterations, subtract_rn_background, reduced_scalar,
                             hyperbolic_charge)
    end
    return row_from_rectangular(st, grid, 2)
end

function geometric_row_du(rows::AbstractVector{<:NLRow}, C::Real)
    T = promote_type(typeof(C), typeof(first(rows).u))
    length(rows) >= 2 || return T(Inf)
    current = rows[end]
    previous = rows[end - 1]
    current.v == previous.v ||
        throw(ArgumentError("geometric step requires matching V grids"))
    du_previous = current.u - previous.u
    du_previous > 0 || throw(ArgumentError("row U coordinates must increase"))
    ru = [(current.r[j] - previous.r[j]) / du_previous for j in eachindex(current.r)]
    candidates = promote_type(eltype(current.r), typeof(C), typeof(du_previous))[]
    for j in 1:length(current.r)-1
        dR = abs(current.r[j + 1] - current.r[j])
        speed = (abs(ru[j]) + abs(ru[j + 1])) / 2
        if isfinite(dR) && isfinite(speed) && dR > 0 && speed > 0
            push!(candidates, C * dR / speed)
        end
    end
    isempty(candidates) && return T(Inf)
    return minimum(candidates)
end

"""
Evolve GP2026 initial data with a row-wise `U` step criterion.

The paper uses `ds^2=-2*f_GP*dU*dV`, while the solver stores
`f_code=2*f_GP`. Therefore its `Delta U=C/f_GP(U,Vmax)` rule is
`Delta U=2C/f_code(U,Vmax)` here. In practice the stiff `f_code` peak can
move away from `Vmax` after apparent-horizon formation, so the default
`step_control=:local` takes the smaller of the largest-`f_code` limiter and
the local geometric condition `|r_U| Delta U <= C |Delta r|`, following the
spirit of Gundlach/Baumgarte/Hilditch arXiv:1908.05971. The literal paper
rule remains available as `step_control=:outer`; `step_control=:max_row`
uses only the largest-`f_code` limiter.
"""
function evolve_gp2026_u_adaptive(initial::NLRow, ep::EvolutionParams; Umax=1.6, C=0.6,
                                  U0=initial.u, V0=first(initial.v), M0=ep.rn.M,
                                  iterations::Int=10, max_rows::Int=100_000,
                                  hyperbolic_charge::Bool=true,
                                  step_control::Symbol=:local,
                                  max_delta_rho=0.25)
    C > 0 || throw(ArgumentError("C must be positive"))
    Umax > initial.u || throw(ArgumentError("Umax must exceed initial U"))
    max_rows >= 2 || throw(ArgumentError("max_rows must be at least 2"))
    step_control in (:outer, :max_row, :geometric, :throat, :local) ||
        throw(ArgumentError("step_control must be :outer, :max_row, :geometric, :throat, or :local"))
    rows = [initial]
    while last(rows).u < Umax && length(rows) < max_rows
        previous = last(rows)
        outer_du = 2C / exp(last(previous.logf))
        max_row_du = 2C / exp(maximum(previous.logf))
        geometric_du = geometric_row_du(rows, C)
        throat_du = throat_row_du(rows, C; max_delta_rho)
        du = if step_control === :outer
            outer_du
        elseif step_control === :max_row
            max_row_du
        elseif step_control === :geometric
            isfinite(geometric_du) ? geometric_du : max_row_du
        elseif step_control === :throat
            isfinite(throat_du) ? throat_du : max_row_du
        else
            minimum((max_row_du, geometric_du, throat_du))
        end
        isfinite(du) && du > 0 || break
        next_u = min(Umax, previous.u + du)
        next_u > previous.u || break
        south = gp2026_na_boundary_point(next_u, ep; U0, V0, M0)
        next = advance_u_row(previous, south, ep;
                             iterations, reduced_scalar=true, hyperbolic_charge)
        push!(rows, next)
        all(isfinite, next.r) && all(isfinite, next.logf) &&
            all(isfinite, next.phi_re) && all(isfinite, next.phi_im) &&
            all(isfinite, next.Au) && all(isfinite, next.Av) &&
            all(isfinite, next.Q) && all(>(zero(eltype(next.r))), next.r) ||
            break
    end
    return UAdaptiveNLState(rows)
end

"""
Advance one adaptive MRT slice by one fixed V step.

`previous` is the completed slice at `V_j`. `northwest` supplies the boundary
point at `(first(new_u), V_{j+1})` on the ingoing initial/boundary leg. The
previous slice is first prolongated onto `new_u`, then the existing nonlinear
Burko-Ori cell update fills the new slice from left to right.
"""
function advance_adaptive_slice(previous::NLSlice, northwest::NLPoint, ep::EvolutionParams;
                                new_u::AbstractVector{<:Real}=previous.u,
                                iterations::Int=5,
                                subtract_rn_background::Bool=false,
                                reduced_scalar::Bool=false)
    require_strictly_increasing(new_u, "new_u")
    require_same_coordinate(first(new_u), northwest.u, "first new_u")
    northwest.v > previous.v ||
        throw(ArgumentError("northwest.v must be larger than the previous slice V"))

    T = promote_type(eltype(new_u), typeof(previous.v), typeof(northwest.v))
    uu = collect(T, new_u)
    prev_on_new = interpolate_slice(previous, uu)
    grid = Grid(uu, T[previous.v, northwest.v])
    st = NLState(grid)

    fill_first_column!(st, prev_on_new)
    fill_northwest_boundary!(st, northwest)
    for i in 1:length(uu)-1
        step_nonlinear_cell!(st, grid, ep, i, 1;
                             iterations, subtract_rn_background, reduced_scalar)
    end
    return slice_from_rectangular(st, grid, 2)
end

"""
Evolve through a prescribed western boundary using slice storage.

Without point splitting or horizon chopping this retains the same U grid on
every slice, providing a regression check against rectangular evolution. With
policies it can refine around, then discard points inside, a detected apparent
horizon and apply Burko-Ori-style U point splitting at completed V bands before
advancing the next slice.
"""
function evolve_adaptive(initial::NLSlice, west_boundary::AbstractVector{<:NLPoint},
                         ep::EvolutionParams; iterations::Int=5,
                         subtract_rn_background::Bool=false,
                         point_splitting::Union{Nothing,PointSplittingConfig}=nothing,
                         horizon_chopping::Union{Nothing,HorizonChoppingConfig}=nothing,
                         horizon_refinement::Union{Nothing,HorizonRefinementConfig}=nothing,
                         reduced_scalar::Bool=false)
    isempty(west_boundary) &&
        throw(ArgumentError("west_boundary needs at least its initial point"))
    require_same_coordinate(first(initial.u), first(west_boundary).u, "initial west U")
    require_same_coordinate(initial.v, first(west_boundary).v, "initial west V")

    slices = [initial]
    next_split_v = isnothing(point_splitting) ? Inf : initial.v + point_splitting.band_width
    next_chop_v = isnothing(horizon_chopping) ? Inf :
                  max(initial.v + horizon_chopping.band_width, horizon_chopping.start_v)
    next_horizon_refine_v = isnothing(horizon_refinement) ? Inf :
                            max(initial.v + horizon_refinement.band_width,
                                horizon_refinement.start_v)
    for j in 2:length(west_boundary)
        west_boundary[j].v > west_boundary[j - 1].v ||
            throw(ArgumentError("west_boundary V values must be strictly increasing"))
        if !isnothing(horizon_refinement) && length(slices) >= 2
            scale_v = max(abs(last(slices).v), abs(next_horizon_refine_v),
                          one(float(last(slices).v)))
            refine_due = last(slices).v >= next_horizon_refine_v -
                         64 * eps(float(scale_v)) * scale_v
            if refine_due
                slices[end] = refine_near_apparent_horizon(slices[end - 1], last(slices),
                                                            horizon_refinement, ep;
                                                            reduced_scalar)
                next_horizon_refine_v += horizon_refinement.band_width
            end
        end
        if !isnothing(horizon_chopping) && length(slices) >= 2
            scale_v = max(abs(last(slices).v), abs(next_chop_v),
                          one(float(last(slices).v)))
            chop_due = last(slices).v >= next_chop_v -
                       64 * eps(float(scale_v)) * scale_v
            if chop_due
                slices[end] = chop_inside_apparent_horizon(slices[end - 1], last(slices),
                                                            horizon_chopping)
                next_chop_v += horizon_chopping.band_width
            end
        end
        if !isnothing(point_splitting)
            scale_v = max(abs(last(slices).v), abs(next_split_v), one(float(last(slices).v)))
            split_due = last(slices).v >= next_split_v - 64 * eps(float(scale_v)) * scale_v
            if split_due
                flags = point_splitting_flags(last(slices), point_splitting)
                if any(flags)
                    slices[end] = refine_slice_constrained(last(slices), flags, ep;
                                                           reduced_scalar)
                end
                next_split_v += point_splitting.band_width
            end
        end
        next = advance_adaptive_slice(last(slices), west_boundary[j], ep;
                                      iterations, subtract_rn_background, reduced_scalar)
        push!(slices, next)
    end
    return AdaptiveNLState(slices)
end
