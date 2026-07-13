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
