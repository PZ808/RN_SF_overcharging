struct FixedBackgroundVRefinementConfig{T<:Real}
    max_passes::Int
    max_points::Int
    vmin::T
    vmax::T
    charge_relative_threshold::T
    energy_relative_threshold::T
    component_relative_threshold::T
    scalar_relative_threshold::T
    residual_threshold::T
    floor::T
    signal_floor_fraction::T
end

function FixedBackgroundVRefinementConfig(;
    max_passes::Int=3,
    max_points::Int=1200,
    vmin=0.0,
    vmax=Inf,
    charge_relative_threshold=0.5,
    energy_relative_threshold=0.5,
    component_relative_threshold=0.75,
    scalar_relative_threshold=0.5,
    residual_threshold=Inf,
    floor=1.0e-30,
    signal_floor_fraction=1.0e-4,
)
    T = promote_type(typeof(vmin), typeof(vmax), typeof(charge_relative_threshold),
                     typeof(energy_relative_threshold),
                     typeof(component_relative_threshold),
                     typeof(scalar_relative_threshold), typeof(residual_threshold),
                     typeof(floor), typeof(signal_floor_fraction))
    return FixedBackgroundVRefinementConfig{T}(
        max_passes,
        max_points,
        convert(T, vmin),
        convert(T, vmax),
        convert(T, charge_relative_threshold),
        convert(T, energy_relative_threshold),
        convert(T, component_relative_threshold),
        convert(T, scalar_relative_threshold),
        convert(T, residual_threshold),
        convert(T, floor),
        convert(T, signal_floor_fraction),
    )
end

struct FixedBackgroundVRefinementSummary{T<:Real}
    pass::Int
    nv::Int
    flagged_intervals::Int
    max_indicator::T
    max_charge_indicator::T
    max_energy_indicator::T
    max_component_indicator::T
    max_scalar_indicator::T
    max_residual_indicator::T
end

struct FixedBackgroundPatchResult{T<:Real}
    state::State{T}
    grid::Grid{T}
    coarse_start::Int
    coarse_stop::Int
    v_refinement_factor::Int
    u_refinement_factor::Int
end

function relative_indicator(a, b, floor)
    if !(isfinite(a) && isfinite(b))
        return zero(promote_type(typeof(a), typeof(b), typeof(floor)))
    end
    return abs(b - a) / max(abs(a), abs(b), floor)
end

function positive_relative_indicator(a, b, floor)
    if !(isfinite(a) && isfinite(b)) || a < 0 || b < 0
        return zero(promote_type(typeof(a), typeof(b), typeof(floor)))
    end
    return relative_indicator(a, b, floor)
end

function window_signal_floor(values, vef, config::FixedBackgroundVRefinementConfig)
    samples = [
        abs(values[j]) for j in eachindex(values)
        if isfinite(values[j]) && config.vmin <= vef[j] <= config.vmax
    ]
    return isempty(samples) ? config.floor :
           max(config.floor, config.signal_floor_fraction * maximum(samples))
end

function fixed_background_v_refinement_indicators(
    st::State,
    g::Grid,
    ep::EvolutionParams,
    config::FixedBackgroundVRefinementConfig=FixedBackgroundVRefinementConfig(),
)
    nv = length(g.v)
    scores = zeros(eltype(g.v), nv - 1)
    charge_scores = zeros(eltype(g.v), nv - 1)
    energy_scores = zeros(eltype(g.v), nv - 1)
    component_scores = zeros(eltype(g.v), nv - 1)
    scalar_scores = zeros(eltype(g.v), nv - 1)
    residual_scores = zeros(eltype(g.v), nv - 1)

    vef, rho_q = horizon_charge_density_series(st, g, ep)
    _, rho_e = horizon_energy_density_series(st, g, ep)
    _, e_qr, e_pr, e_qv, e_pv = horizon_energy_density_divided_components(st, g, ep)
    P = sqrt.(st.xi[end, :] .^ 2 .+ st.pi[end, :] .^ 2)
    maxwell_u, maxwell_v = maxwell_residuals(st, g, ep)

    charge_floor = window_signal_floor(rho_q, vef, config)
    energy_floor = window_signal_floor(rho_e, vef, config)
    scalar_floor = window_signal_floor(P, vef, config)
    qr_floor = window_signal_floor(e_qr, vef, config)
    pr_floor = window_signal_floor(e_pr, vef, config)
    qv_floor = window_signal_floor(e_qv, vef, config)
    pv_floor = window_signal_floor(e_pv, vef, config)

    for j in 1:nv-1
        vmid = (vef[j] + vef[j + 1]) / 2
        if !(config.vmin <= vmid <= config.vmax)
            continue
        end

        charge_scores[j] =
            relative_indicator(abs(rho_q[j]), abs(rho_q[j + 1]), charge_floor)
        energy_scores[j] =
            positive_relative_indicator(rho_e[j], rho_e[j + 1], energy_floor)
        scalar_scores[j] =
            positive_relative_indicator(P[j], P[j + 1], scalar_floor)

        component_scores[j] = maximum((
            positive_relative_indicator(e_qr[j], e_qr[j + 1], qr_floor),
            positive_relative_indicator(e_pr[j], e_pr[j + 1], pr_floor),
            positive_relative_indicator(e_qv[j], e_qv[j + 1], qv_floor),
            positive_relative_indicator(e_pv[j], e_pv[j + 1], pv_floor),
        ))

        residual_scores[j] = max(maximum(abs, maxwell_u[:, j]),
                                 maximum(abs, maxwell_v[:, j]))
        residual_indicator =
            isfinite(config.residual_threshold) && config.residual_threshold > 0 ?
            residual_scores[j] / config.residual_threshold :
            zero(residual_scores[j])

        scores[j] = maximum((
            charge_scores[j] / config.charge_relative_threshold,
            energy_scores[j] / config.energy_relative_threshold,
            component_scores[j] / config.component_relative_threshold,
            scalar_scores[j] / config.scalar_relative_threshold,
            residual_indicator,
        ))
    end

    return scores, charge_scores, energy_scores, component_scores, scalar_scores,
           residual_scores
end

function finite_max(x)
    values = [v for v in x if isfinite(v)]
    return isempty(values) ? NaN : maximum(values)
end

function fixed_background_v_refinement_flags(
    st::State,
    g::Grid,
    ep::EvolutionParams,
    config::FixedBackgroundVRefinementConfig=FixedBackgroundVRefinementConfig(),
)
    scores, charge_scores, energy_scores, component_scores, scalar_scores,
        residual_scores = fixed_background_v_refinement_indicators(st, g, ep, config)
    flags = scores .> 1
    return flags, scores, charge_scores, energy_scores, component_scores,
           scalar_scores, residual_scores
end

function refine_fixed_background_v_grid(
    g::Grid,
    ep::EvolutionParams,
    flags::AbstractVector{Bool},
    scores::AbstractVector,
    max_points::Int,
)
    vef = [ef_v_from_mrt(v, ep.rn) for v in g.v]
    selected = copy(flags)
    available = max_points - length(g.v)
    if available <= 0
        selected .= false
    elseif count(selected) > available
        flagged = findall(selected)
        order = sortperm(scores[flagged]; rev=true)
        keep = Set(flagged[order[1:available]])
        selected .= [j in keep for j in eachindex(selected)]
    end

    refined_vef = eltype(vef)[]
    for j in 1:length(vef)-1
        push!(refined_vef, vef[j])
        selected[j] && push!(refined_vef, (vef[j] + vef[j + 1]) / 2)
    end
    push!(refined_vef, last(vef))
    refined_v = [compact_v_from_ef_v(V, ep.rn) for V in refined_vef]
    return Grid(copy(g.u), refined_v), selected
end

function evolve_fixed_background_v_adaptive(
    ep::EvolutionParams;
    nu::Int=120,
    nv::Int=360,
    u0=-1.0,
    u1=-1.0e-5,
    Vef0=0.0,
    Vef1=600.0,
    envelope=gaussian_envelope,
    config::FixedBackgroundVRefinementConfig=FixedBackgroundVRefinementConfig(),
)
    u = collect(range(u0, u1; length=nu))
    vef = collect(range(Vef0, Vef1; length=nv))
    grid = Grid(u, [compact_v_from_ef_v(V, ep.rn) for V in vef])
    summaries = FixedBackgroundVRefinementSummary{eltype(grid.v)}[]
    state = initialize_state(grid, ep; envelope)
    evolve!(state, grid, ep)

    for pass in 0:config.max_passes
        flags, scores, charge_scores, energy_scores, component_scores,
            scalar_scores, residual_scores =
            fixed_background_v_refinement_flags(state, grid, ep, config)
        push!(summaries, FixedBackgroundVRefinementSummary(
            pass,
            length(grid.v),
            count(flags),
            finite_max(scores),
            finite_max(charge_scores),
            finite_max(energy_scores),
            finite_max(component_scores),
            finite_max(scalar_scores),
            finite_max(residual_scores),
        ))

        pass == config.max_passes && break
        count(flags) == 0 && break
        length(grid.v) >= config.max_points && break

        grid, selected = refine_fixed_background_v_grid(
            grid,
            ep,
            flags,
            scores,
            config.max_points,
        )
        count(selected) == 0 && break
        state = initialize_state(grid, ep; envelope)
        evolve!(state, grid, ep)
    end

    return state, grid, summaries
end

function linear_interpolate_series(x::AbstractVector, y::AbstractVector, xq)
    length(x) == length(y) || throw(ArgumentError("x and y must have matching lengths"))
    lo, hi = x[firstindex(x)], x[lastindex(x)]
    scale = max(abs(float(lo)), abs(float(hi)), abs(float(xq)), 1.0)
    tol = 100 * eps(scale)
    if xq < lo
        xq >= lo - tol ||
            throw(ArgumentError("query point outside interpolation range"))
        return y[firstindex(y)]
    elseif xq > hi
        xq <= hi + tol ||
            throw(ArgumentError("query point outside interpolation range"))
        return y[lastindex(y)]
    end
    idx = searchsortedlast(x, xq)
    if idx == lastindex(x)
        return y[idx]
    end
    idx = max(idx, firstindex(x))
    x0, x1 = x[idx], x[idx + 1]
    x1 == x0 && return y[idx]
    weight = (xq - x0) / (x1 - x0)
    return (1 - weight) * y[idx] + weight * y[idx + 1]
end

function refined_patch_vef(vef::AbstractVector, start_index::Int, stop_index::Int,
                           refinement_factor::Int)
    1 <= start_index < stop_index <= length(vef) ||
        throw(ArgumentError("patch indices must satisfy 1 <= start < stop <= nv"))
    refinement_factor >= 1 || throw(ArgumentError("refinement_factor must be positive"))
    refined = eltype(vef)[]
    for j in start_index:stop_index-1
        left, right = vef[j], vef[j + 1]
        for sub in 0:refinement_factor-1
            push!(refined, left + (right - left) * sub / refinement_factor)
        end
    end
    push!(refined, vef[stop_index])
    return refined
end

function refined_axis(values::AbstractVector, refinement_factor::Int)
    refinement_factor >= 1 || throw(ArgumentError("refinement_factor must be positive"))
    refined = eltype(values)[]
    for i in firstindex(values):lastindex(values)-1
        left, right = values[i], values[i + 1]
        for sub in 0:refinement_factor-1
            push!(refined, left + (right - left) * sub / refinement_factor)
        end
    end
    push!(refined, last(values))
    return refined
end

function fixed_background_patch_indices(g::Grid, ep::EvolutionParams, vmin, vmax)
    vef = [ef_v_from_mrt(v, ep.rn) for v in g.v]
    vmin <= vmax || throw(ArgumentError("vmin must be <= vmax"))
    start_index = searchsortedlast(vef, vmin)
    start_index = clamp(start_index, firstindex(vef), lastindex(vef) - 1)
    stop_index = searchsortedfirst(vef, vmax)
    stop_index = clamp(stop_index, start_index + 1, lastindex(vef))
    return start_index, stop_index
end

function evolve_fixed_background_v_patch(
    coarse_state::State,
    coarse_grid::Grid,
    ep::EvolutionParams;
    vmin,
    vmax,
    refinement_factor::Int=4,
    u_refinement_factor::Int=1,
    envelope=gaussian_envelope,
)
    start_index, stop_index = fixed_background_patch_indices(coarse_grid, ep, vmin, vmax)
    coarse_vef = [ef_v_from_mrt(v, ep.rn) for v in coarse_grid.v]
    patch_vef = refined_patch_vef(coarse_vef, start_index, stop_index,
                                  refinement_factor)
    patch_u = refined_axis(coarse_grid.u, u_refinement_factor)
    patch_grid = Grid(patch_u,
                      [compact_v_from_ef_v(V, ep.rn) for V in patch_vef])
    boundary_vef = sort(unique(vcat(coarse_vef[begin:start_index], patch_vef)))
    boundary_grid = Grid(copy(coarse_grid.u),
                         [compact_v_from_ef_v(V, ep.rn) for V in boundary_vef])
    boundary_state = initialize_state(boundary_grid, ep; envelope)
    patch_state = State(promote_type(eltype(coarse_grid.u), eltype(coarse_state.xi)),
                        length(patch_grid.u), length(patch_grid.v))

    fields = (:xi, :pi, :Au, :Av, :Q)
    for field in fields
        coarse_values = getfield(coarse_state, field)
        patch_values = getfield(patch_state, field)
        boundary_values = getfield(boundary_state, field)
        for i in eachindex(patch_grid.u)
            patch_values[i, firstindex(patch_grid.v)] =
                linear_interpolate_series(coarse_grid.u, coarse_values[:, start_index],
                                          patch_grid.u[i])
        end
        for j in eachindex(patch_grid.v)
            patch_values[firstindex(patch_grid.u), j] =
                linear_interpolate_series(boundary_vef,
                                          boundary_values[firstindex(boundary_grid.u), :],
                                          patch_vef[j])
        end
    end

    evolve!(patch_state, patch_grid, ep)
    return FixedBackgroundPatchResult(
        patch_state,
        patch_grid,
        start_index,
        stop_index,
        refinement_factor,
        u_refinement_factor,
    )
end

function fixed_background_patch_parent_errors(
    patch::FixedBackgroundPatchResult,
    coarse_state::State,
)
    coarse_rows = axes(coarse_state.xi, 1)
    coarse_cols = patch.coarse_start:patch.coarse_stop
    patch_rows = firstindex(patch.grid.u):patch.u_refinement_factor:lastindex(patch.grid.u)
    patch_cols = firstindex(patch.grid.v):patch.v_refinement_factor:lastindex(patch.grid.v)

    length(patch_rows) == length(coarse_rows) ||
        throw(ArgumentError("patch parent U nodes do not match coarse U nodes"))
    length(patch_cols) == length(coarse_cols) ||
        throw(ArgumentError("patch parent V nodes do not match coarse V nodes"))

    return (;
        xi=maximum(abs.(patch.state.xi[patch_rows, patch_cols] .-
                        coarse_state.xi[coarse_rows, coarse_cols])),
        pi=maximum(abs.(patch.state.pi[patch_rows, patch_cols] .-
                        coarse_state.pi[coarse_rows, coarse_cols])),
        Au=maximum(abs.(patch.state.Au[patch_rows, patch_cols] .-
                        coarse_state.Au[coarse_rows, coarse_cols])),
        Av=maximum(abs.(patch.state.Av[patch_rows, patch_cols] .-
                        coarse_state.Av[coarse_rows, coarse_cols])),
        Q=maximum(abs.(patch.state.Q[patch_rows, patch_cols] .-
                       coarse_state.Q[coarse_rows, coarse_cols])),
    )
end
