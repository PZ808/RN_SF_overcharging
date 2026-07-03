"""
Configuration for the simplified Hamade-Stewart Berger-Oliger hierarchy.

`max_levels` counts the root grid. Each child refines both null directions by
`refinement_factor`; Stewart used a factor of four and one buffered error
cluster per level.
"""
struct StewartAMRConfig{T<:Real}
    refinement_factor::Int
    revision_interval::Int
    max_levels::Int
    fields::Tuple
    atol::T
    rtol::T
    order::Int
    buffer_points::Int
    reintegrate::Bool
end

function StewartAMRConfig(;
    refinement_factor::Int=4,
    revision_interval::Int=4,
    max_levels::Int=4,
    fields=DEFAULT_ROW_LTE_FIELDS,
    atol=1.0e-8,
    rtol=1.0e-5,
    order::Int=2,
    buffer_points::Int=4,
    reintegrate::Bool=true,
)
    refinement_factor >= 2 ||
        throw(ArgumentError("refinement_factor must be at least two"))
    revision_interval >= 1 ||
        throw(ArgumentError("revision_interval must be positive"))
    max_levels >= 1 || throw(ArgumentError("max_levels must be positive"))
    order >= 1 || throw(ArgumentError("order must be positive"))
    buffer_points >= 0 ||
        throw(ArgumentError("buffer_points must be nonnegative"))
    atol > 0 || throw(ArgumentError("atol must be positive"))
    rtol >= 0 || throw(ArgumentError("rtol must be nonnegative"))
    T = promote_type(typeof(atol), typeof(rtol))
    return StewartAMRConfig{T}(
        refinement_factor,
        revision_interval,
        max_levels,
        Tuple(fields),
        convert(T, atol),
        convert(T, rtol),
        order,
        buffer_points,
        reintegrate,
    )
end

mutable struct StewartAMRStats
    level_steps::Vector{Int}
    revisions::Vector{Int}
    last_lte::Vector{Float64}
    max_lte::Vector{Float64}
    child_creations::Int
    child_destructions::Int
    injections::Int
    suffix_reintegrations::Int
    precision_fallbacks::Int
    max_level_reached::Int
end

StewartAMRStats() =
    StewartAMRStats(Int[], Int[], Float64[], Float64[], 0, 0, 0, 0, 0, 0)

mutable struct StewartAMRLevel{T<:Real}
    level::Int
    parent_interval::UnitRange{Int}
    current::NLRow{T}
    child::Union{Nothing,StewartAMRLevel{T}}
    steps_since_revision::Int
end

mutable struct StewartAMRHierarchy{T<:Real,C<:Real}
    root::StewartAMRLevel{T}
    config::StewartAMRConfig{C}
    stats::StewartAMRStats
    initial_u::T
    initial_v::T
end

function initialize_stewart_hierarchy(
    initial::NLRow{T};
    config::StewartAMRConfig=StewartAMRConfig(),
) where {T<:Real}
    root = StewartAMRLevel{T}(
        0,
        firstindex(initial.v):lastindex(initial.v),
        initial,
        nothing,
        config.revision_interval,
    )
    return StewartAMRHierarchy(
        root,
        config,
        StewartAMRStats(),
        initial.u,
        first(initial.v),
    )
end

function ensure_stat_level!(values::Vector{Int}, level::Int)
    while length(values) <= level
        push!(values, 0)
    end
    return values
end

function record_level_step!(stats::StewartAMRStats, level::Int)
    ensure_stat_level!(stats.level_steps, level)
    stats.level_steps[level + 1] += 1
    stats.max_level_reached = max(stats.max_level_reached, level)
    return stats
end

function record_revision!(
    stats::StewartAMRStats,
    level::Int,
    max_error,
)
    ensure_stat_level!(stats.revisions, level)
    stats.revisions[level + 1] += 1
    while length(stats.last_lte) <= level
        push!(stats.last_lte, NaN)
        push!(stats.max_lte, 0.0)
    end
    value = isnothing(max_error) ? NaN : Float64(max_error)
    stats.last_lte[level + 1] = value
    if isfinite(value)
        stats.max_lte[level + 1] =
            max(stats.max_lte[level + 1], value)
    end
    return stats
end

function row_point_at_v(row::NLRow, v::Real)
    first(row.v) <= v <= last(row.v) ||
        throw(ArgumentError("requested V lies outside the row"))
    index = searchsortedfirst(row.v, v)
    if index <= length(row.v) && row.v[index] == v
        return row_point(row, index)
    end
    cell = clamp(
        searchsortedlast(row.v, v),
        firstindex(row.v),
        lastindex(row.v) - 1,
    )
    values = (
        local_polynomial_value(row.v, row.r, v, cell),
        local_polynomial_value(row.v, row.logf, v, cell),
        local_polynomial_value(row.v, row.phi_re, v, cell),
        local_polynomial_value(row.v, row.phi_im, v, cell),
        local_polynomial_value(row.v, row.Au, v, cell),
        local_polynomial_value(row.v, row.Av, v, cell),
        local_polynomial_value(row.v, row.Q, v, cell),
    )
    return NLPoint(row.u, v, values...)
end

function interpolate_parent_boundary(
    lower::NLRow,
    upper::NLRow,
    v::Real,
    target_u::Real,
)
    lower.v == upper.v ||
        throw(ArgumentError("parent boundary rows must share their V grid"))
    lower.u <= target_u <= upper.u ||
        throw(ArgumentError("boundary interpolation target lies outside parent step"))
    lower_point = row_point_at_v(lower, v)
    target_u == lower.u && return lower_point
    upper_point = row_point_at_v(upper, v)
    target_u == upper.u && return upper_point
    weight = (target_u - lower.u) / (upper.u - lower.u)
    blend(a, b) = (one(weight) - weight) * a + weight * b
    return NLPoint(
        target_u,
        v,
        blend(lower_point.r, upper_point.r),
        blend(lower_point.logf, upper_point.logf),
        blend(lower_point.phi_re, upper_point.phi_re),
        blend(lower_point.phi_im, upper_point.phi_im),
        blend(lower_point.Au, upper_point.Au),
        blend(lower_point.Av, upper_point.Av),
        blend(lower_point.Q, upper_point.Q),
    )
end

function stewart_row_lte_error(
    full_step::NLRow,
    half_step::NLRow,
    config::StewartAMRConfig,
)
    full_step.v == half_step.v ||
        throw(ArgumentError("LTE comparison rows must share the same V grid"))
    richardson = 2^config.order - 1
    T = promote_type(
        eltype(full_step.r),
        eltype(half_step.r),
        typeof(config.atol),
        typeof(config.rtol),
    )
    error = zeros(T, length(full_step.v))
    for field in config.fields
        full_values = row_field(full_step, field)
        half_values = row_field(half_step, field)
        for j in eachindex(error)
            difference =
                abs(full_values[j] - half_values[j]) / richardson
            scale = convert(T, config.atol) +
                    convert(T, config.rtol) *
                    max(abs(full_values[j]), abs(half_values[j]))
            value = isfinite(difference) && isfinite(scale) && scale > 0 ?
                    difference / scale : T(Inf)
            error[j] = max(error[j], value)
        end
    end
    return error
end

function stewart_row_lte(
    previous::NLRow,
    target_u::Real,
    boundary_point,
    ep::EvolutionParams,
    config::StewartAMRConfig;
    iterations::Int=10,
    hyperbolic_charge::Bool=true,
    cell_solver::Symbol=:newton_direct,
    newton_rtol=1.0e-13,
    newton_atol=1.0e-15,
)
    target_u > previous.u ||
        throw(ArgumentError("target U must exceed the previous row"))
    midpoint_u = (previous.u + target_u) / 2
    previous.u < midpoint_u < target_u || return nothing
    midpoint_boundary = boundary_point(midpoint_u)
    target_boundary = boundary_point(target_u)
    full = advance_u_row(
        previous, target_boundary, ep;
        iterations,
        reduced_scalar=true,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
    midpoint = advance_u_row(
        previous, midpoint_boundary, ep;
        iterations,
        reduced_scalar=true,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
    half = advance_u_row(
        midpoint, target_boundary, ep;
        iterations,
        reduced_scalar=true,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
    error = stewart_row_lte_error(full, half, config)
    flags = [value > one(value) for value in error]
    intervals = buffered_flag_intervals(
        flags;
        buffer_points=config.buffer_points,
        cluster=:all,
    )
    return (
        full=full,
        half=half,
        midpoint=midpoint,
        error=error,
        max_error=maximum(error),
        flags=flags,
        intervals=intervals,
    )
end

function destroy_stewart_child!(
    level::StewartAMRLevel,
    stats::StewartAMRStats,
)
    isnothing(level.child) && return level
    level.child = nothing
    stats.child_destructions += 1
    return level
end

function stewart_descendant_v_ranges(level::StewartAMRLevel{T}) where {T<:Real}
    ranges = Tuple{T,T}[]
    current = level.child
    while !isnothing(current)
        push!(ranges, (first(current.current.v), last(current.current.v)))
        current = current.child
    end
    return ranges
end

function parent_interval_covering_v_range(
    v::AbstractVector{<:Real},
    first_v::Real,
    last_v::Real,
)
    first(v) <= first_v <= last_v <= last(v) ||
        throw(ArgumentError("descendant V range lies outside rebuilt parent"))
    lo = clamp(searchsortedlast(v, first_v), firstindex(v), lastindex(v) - 1)
    hi = clamp(searchsortedfirst(v, last_v), lo + 1, lastindex(v))
    return lo:hi
end

function rebuild_stewart_child!(
    level::StewartAMRLevel{T},
    interval::UnitRange{Int},
    config::StewartAMRConfig,
    stats::StewartAMRStats,
) where {T<:Real}
    level.level + 1 < config.max_levels || return level
    descendant_ranges = isnothing(level.child) ?
                        Tuple{T,T}[] :
                        stewart_descendant_v_ranges(level.child)
    patch = normalized_patch_interval(interval, length(level.current.v))
    if !isempty(descendant_ranges)
        descendant_patch = parent_interval_covering_v_range(
            level.current.v,
            first(first(descendant_ranges)),
            last(first(descendant_ranges)),
        )
        patch = (
            min(first(patch), first(descendant_patch)):
            max(last(patch), last(descendant_patch))
        )
    end
    fine_v = refined_v_patch_grid(
        level.current.v, patch;
        refinement_factor=config.refinement_factor,
    )
    child_row = interpolate_row(level.current, fine_v)
    if !isnothing(level.child)
        stats.child_destructions += 1
    end
    level.child = StewartAMRLevel{T}(
        level.level + 1,
        patch,
        child_row,
        nothing,
        config.revision_interval,
    )
    stats.child_creations += 1
    stats.max_level_reached = max(stats.max_level_reached, level.level + 1)

    parent = level.child
    for (first_v, last_v) in descendant_ranges
        parent.level + 1 < config.max_levels || break
        descendant_patch = parent_interval_covering_v_range(
            parent.current.v, first_v, last_v,
        )
        descendant_v = refined_v_patch_grid(
            parent.current.v,
            descendant_patch;
            refinement_factor=config.refinement_factor,
        )
        descendant_row = interpolate_row(parent.current, descendant_v)
        parent.child = StewartAMRLevel{T}(
            parent.level + 1,
            descendant_patch,
            descendant_row,
            nothing,
            config.revision_interval,
        )
        stats.child_creations += 1
        stats.max_level_reached =
            max(stats.max_level_reached, parent.level + 1)
        parent = parent.child
    end
    return level
end

function stewart_hierarchy_depth(level::StewartAMRLevel)
    return isnothing(level.child) ? level.level + 1 :
           stewart_hierarchy_depth(level.child)
end

function stewart_hierarchy_intervals(level::StewartAMRLevel)
    intervals = UnitRange{Int}[level.parent_interval]
    current = level
    while !isnothing(current.child)
        current = current.child
        push!(intervals, current.parent_interval)
    end
    return intervals
end

function validate_stewart_hierarchy(level::StewartAMRLevel)
    require_strictly_increasing(level.current.v, "Stewart level V")
    child = level.child
    isnothing(child) && return true
    patch = normalized_patch_interval(
        child.parent_interval, length(level.current.v),
    )
    expected = refined_v_patch_grid(
        level.current.v, patch;
        refinement_factor=(
            length(child.current.v) - 1
        ) ÷ (length(patch) - 1),
    )
    child.current.v == expected ||
        throw(ArgumentError("child grid is not nested in its parent patch"))
    child.current.u == level.current.u ||
        throw(ArgumentError("Stewart levels are not synchronized in U"))
    return validate_stewart_hierarchy(child)
end

function stewart_subtree_has_error(
    level::StewartAMRLevel,
    parent_lower::NLRow,
    parent_upper::NLRow,
    ep::EvolutionParams,
    config::StewartAMRConfig;
    iterations::Int=10,
    hyperbolic_charge::Bool=true,
    cell_solver::Symbol=:newton_direct,
    newton_rtol=1.0e-13,
    newton_atol=1.0e-15,
)
    refinement_factor = config.refinement_factor
    target_u = level.current.u +
               (parent_upper.u - parent_lower.u) / refinement_factor
    target_u > level.current.u || return true
    child_v0 = first(level.current.v)
    boundary_point = u -> interpolate_parent_boundary(
        parent_lower, parent_upper, child_v0, u,
    )
    estimate = stewart_row_lte(
        level.current,
        target_u,
        boundary_point,
        ep,
        config;
        iterations,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
    isnothing(estimate) && return true
    !isempty(estimate.intervals) && return true
    isnothing(level.child) && return false
    return stewart_subtree_has_error(
        level.child,
        level.current,
        estimate.full,
        ep,
        config;
        iterations,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
end

function advance_stewart_level!(
    level::StewartAMRLevel,
    target_u::Real,
    boundary_point,
    ep::EvolutionParams,
    config::StewartAMRConfig,
    stats::StewartAMRStats;
    iterations::Int=10,
    hyperbolic_charge::Bool=true,
    cell_solver::Symbol=:newton_direct,
    newton_rtol=1.0e-13,
    newton_atol=1.0e-15,
)
    lower = level.current
    target_u > lower.u ||
        throw(ArgumentError("Stewart target U must exceed current level U"))
    revise = level.steps_since_revision >= config.revision_interval
    estimate = revise ? stewart_row_lte(
        lower,
        target_u,
        boundary_point,
        ep,
        config;
        iterations,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    ) : nothing

    upper = if isnothing(estimate)
        if revise
            stats.precision_fallbacks += 1
        end
        advance_u_row(
            lower, boundary_point(target_u), ep;
            iterations,
            reduced_scalar=true,
            hyperbolic_charge,
            cell_solver,
            newton_rtol,
            newton_atol,
        )
    else
        estimate.full
    end

    if revise
        record_revision!(
            stats,
            level.level,
            isnothing(estimate) ? nothing : estimate.max_error,
        )
        child_has_error = !isnothing(level.child) &&
                          !isnothing(estimate) &&
                          stewart_subtree_has_error(
            level.child,
            lower,
            upper,
            ep,
            config;
            iterations,
            hyperbolic_charge,
            cell_solver,
            newton_rtol,
            newton_atol,
        )
        if level.level + 1 >= config.max_levels ||
           isnothing(estimate) ||
           (isempty(estimate.intervals) && !child_has_error)
            destroy_stewart_child!(level, stats)
        elseif !isempty(estimate.intervals)
            rebuild_stewart_child!(
                level,
                only(estimate.intervals),
                config,
                stats,
            )
        end
        level.steps_since_revision = 0
    end

    child = level.child
    if !isnothing(child)
        refinement_factor = config.refinement_factor
        child_start_u = child.current.u
        child_start_u == lower.u ||
            throw(ArgumentError("child and parent are not synchronized before evolution"))
        child_v0 = first(child.current.v)
        parent_boundary = u -> interpolate_parent_boundary(
            lower, upper, child_v0, u,
        )
        parent_du = target_u - child_start_u
        for substep in 1:refinement_factor
            child_target_u = substep == refinement_factor ?
                             target_u :
                             child_start_u +
                             parent_du * (substep / refinement_factor)
            advance_stewart_level!(
                child,
                child_target_u,
                parent_boundary,
                ep,
                config,
                stats;
                iterations,
                hyperbolic_charge,
                cell_solver,
                newton_rtol,
                newton_atol,
            )
        end
        upper = inject_row_patch(
            upper,
            child.current,
            child.parent_interval;
            refinement_factor,
        )
        stats.injections += 1
        if config.reintegrate &&
           last(child.parent_interval) < lastindex(upper.v)
            upper = reintegrate_row_suffix(
                lower,
                upper,
                last(child.parent_interval),
                ep;
                iterations,
                hyperbolic_charge,
                cell_solver,
                newton_rtol,
                newton_atol,
            )
            stats.suffix_reintegrations += 1
        end
    end

    level.current = upper
    level.steps_since_revision += 1
    record_level_step!(stats, level.level)
    return (
        row=upper,
        estimate=estimate,
        revised=revise,
        has_child=!isnothing(level.child),
        level=level.level,
    )
end

function advance_stewart_hierarchy!(
    hierarchy::StewartAMRHierarchy,
    target_u::Real,
    ep::EvolutionParams;
    U0=hierarchy.initial_u,
    V0=hierarchy.initial_v,
    M0=ep.rn.M,
    pulse_leg_gauge::Symbol=:areal_affine,
    iterations::Int=10,
    hyperbolic_charge::Bool=true,
    cell_solver::Symbol=:newton_direct,
    newton_rtol=1.0e-13,
    newton_atol=1.0e-15,
)
    boundary_point = u -> gp2026_na_boundary_point(
        u, ep;
        U0,
        V0,
        M0,
        pulse_leg_gauge,
    )
    result = advance_stewart_level!(
        hierarchy.root,
        target_u,
        boundary_point,
        ep,
        hierarchy.config,
        hierarchy.stats;
        iterations,
        hyperbolic_charge,
        cell_solver,
        newton_rtol,
        newton_atol,
    )
    validate_stewart_hierarchy(hierarchy.root)
    return merge(
        result,
        (
            depth=stewart_hierarchy_depth(hierarchy.root),
            stats=hierarchy.stats,
        ),
    )
end
