"""
Legacy fixed-background metric container from the first scaffold.

The active fully nonlinear double-null code stores `r`, `logf`, scalar fields,
potentials, and charge together in `NLState`, and advances them through
`evolve_nonlinear!` / `evolve_gp2026_u_adaptive`.
"""
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
        "`einstein_backreaction_rhs!` is a legacy scaffold and is not the " *
        "active nonlinear evolution path. Use `NLState` with " *
        "`evolve_nonlinear!` or the GP2026 row marcher " *
        "`evolve_gp2026_u_adaptive`; the nonlinear metric equations are " *
        "implemented there via `metric_rhs` and `step_nonlinear_cell!`."
    ))
end
