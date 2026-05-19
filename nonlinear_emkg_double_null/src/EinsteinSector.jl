struct MetricState{T<:Real}
    r::Matrix{T}
    f::Matrix{T}
end

function MetricState(g::Grid{T}, p::RNParams{T}) where {T<:Real}
    nu, nv = size(g)
    r = zeros(T, nu, nv)
    f = zeros(T, nu, nv)
    for i in 1:nu, j in 1:nv
        r[i, j] = areal_radius(g.u[i], g.v[j], p)
        f[i, j] = metric_f(g.u[i], g.v[j], p)
    end
    return MetricState(r, f)
end

function einstein_backreaction_rhs!(duv_metric::MetricState, metric::MetricState,
                                    matter::State, g::Grid, ep::EvolutionParams)
    throw(ErrorException(
        "The fully nonlinear Einstein sector is not implemented yet. " *
        "Next step: derive double-null equations for r(u,v) and f(u,v) " *
        "from the EMKG stress tensor and replace the fixed RN metric calls."
    ))
end

