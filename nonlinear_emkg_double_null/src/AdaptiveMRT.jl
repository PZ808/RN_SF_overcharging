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

function HorizonChoppingConfig(; start_v=0.0, band_width=1.0, interior_buffer_cells::Int=1)
    start_v >= 0 || throw(ArgumentError("start_v must be nonnegative"))
    band_width > 0 || throw(ArgumentError("band_width must be positive"))
    interior_buffer_cells >= 0 ||
        throw(ArgumentError("interior_buffer_cells must be nonnegative"))
    T = promote_type(typeof(start_v), typeof(band_width))
    return HorizonChoppingConfig{T}(convert(T, start_v), convert(T, band_width),
                                    interior_buffer_cells)
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

function adaptive_state_from_rectangular(st::NLState, g::Grid)
    return AdaptiveNLState([slice_from_rectangular(st, g, j) for j in eachindex(g.v)])
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
    return last_point == length(current.u) ? current : truncate_slice(current, last_point)
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
                                subtract_rn_background::Bool=false)
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
        step_nonlinear_cell!(st, grid, ep, i, 1; iterations, subtract_rn_background)
    end
    return slice_from_rectangular(st, grid, 2)
end

"""
Evolve through a prescribed western boundary using slice storage.

Without point splitting or horizon chopping this retains the same U grid on
every slice, providing a regression check against rectangular evolution. With
policies it can discard points inside a detected apparent horizon and apply
Burko-Ori-style U point splitting at completed V bands before advancing the
next slice.
"""
function evolve_adaptive(initial::NLSlice, west_boundary::AbstractVector{<:NLPoint},
                         ep::EvolutionParams; iterations::Int=5,
                         subtract_rn_background::Bool=false,
                         point_splitting::Union{Nothing,PointSplittingConfig}=nothing,
                         horizon_chopping::Union{Nothing,HorizonChoppingConfig}=nothing)
    isempty(west_boundary) &&
        throw(ArgumentError("west_boundary needs at least its initial point"))
    require_same_coordinate(first(initial.u), first(west_boundary).u, "initial west U")
    require_same_coordinate(initial.v, first(west_boundary).v, "initial west V")

    slices = [initial]
    next_split_v = isnothing(point_splitting) ? Inf : initial.v + point_splitting.band_width
    next_chop_v = isnothing(horizon_chopping) ? Inf :
                  max(initial.v + horizon_chopping.band_width, horizon_chopping.start_v)
    for j in 2:length(west_boundary)
        west_boundary[j].v > west_boundary[j - 1].v ||
            throw(ArgumentError("west_boundary V values must be strictly increasing"))
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
                    slices[end] = refine_slice(last(slices), flags)
                end
                next_split_v += point_splitting.band_width
            end
        end
        next = advance_adaptive_slice(last(slices), west_boundary[j], ep;
                                      iterations, subtract_rn_background)
        push!(slices, next)
    end
    return AdaptiveNLState(slices)
end
