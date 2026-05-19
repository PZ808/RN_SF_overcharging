function maxwell_residuals(st::State, g::Grid, ep::EvolutionParams)
    nu, nv = size(g)
    ru = zeros(eltype(g.u), nu - 1, nv - 1)
    rv = zeros(eltype(g.u), nu - 1, nv - 1)
    e = ep.scalar_charge
    for i in 1:nu-1, j in 1:nv-1
        du = g.u[i + 1] - g.u[i]
        dv = g.v[j + 1] - g.v[j]
        uc = (g.u[i] + g.u[i + 1]) / 2
        vc = (g.v[j] + g.v[j + 1]) / 2
        r = areal_radius(uc, vc, ep.rn)

        xi = (st.xi[i, j] + st.xi[i + 1, j] + st.xi[i, j + 1] + st.xi[i + 1, j + 1]) / 4
        pii = (st.pi[i, j] + st.pi[i + 1, j] + st.pi[i, j + 1] + st.pi[i + 1, j + 1]) / 4
        Au = (st.Au[i, j] + st.Au[i + 1, j] + st.Au[i, j + 1] + st.Au[i + 1, j + 1]) / 4
        Av = (st.Av[i, j] + st.Av[i + 1, j] + st.Av[i, j + 1] + st.Av[i + 1, j + 1]) / 4

        xi_u = ((st.xi[i + 1, j] - st.xi[i, j]) + (st.xi[i + 1, j + 1] - st.xi[i, j + 1])) / (2du)
        pi_u = ((st.pi[i + 1, j] - st.pi[i, j]) + (st.pi[i + 1, j + 1] - st.pi[i, j + 1])) / (2du)
        xi_v = ((st.xi[i, j + 1] - st.xi[i, j]) + (st.xi[i + 1, j + 1] - st.xi[i + 1, j])) / (2dv)
        pi_v = ((st.pi[i, j + 1] - st.pi[i, j]) + (st.pi[i + 1, j + 1] - st.pi[i + 1, j])) / (2dv)
        q_u = ((st.Q[i + 1, j] - st.Q[i, j]) + (st.Q[i + 1, j + 1] - st.Q[i, j + 1])) / (2du)
        q_v = ((st.Q[i, j + 1] - st.Q[i, j]) + (st.Q[i + 1, j + 1] - st.Q[i + 1, j])) / (2dv)

        Ju = current_u(xi, pii, xi_u, pi_u, Au, e) / r^2
        Jv = current_v(xi, pii, xi_v, pi_v, Av, e) / r^2
        ru[i, j] = q_u - 4pi * r^2 * Ju
        rv[i, j] = q_v + 4pi * r^2 * Jv
    end
    return ru, rv
end

function horizon_phi_series(st::NLState, g::Grid, ep::EvolutionParams; i::Int=lastindex(g.u))
    vef = [ef_v_from_mrt(v, ep.rn) for v in g.v]
    phi_abs = sqrt.(st.phi_re[i, :].^2 .+ st.phi_im[i, :].^2)
    return vef, phi_abs
end

function fit_power_law(x, y; xmin=nothing, xmax=nothing)
    idx = trues(length(x))
    if xmin !== nothing
        idx .&= x .>= xmin
    end
    if xmax !== nothing
        idx .&= x .<= xmax
    end
    idx .&= y .> 0
    xs = log.(x[idx])
    ys = log.(y[idx])
    length(xs) >= 2 || throw(ArgumentError("Need at least two positive samples for a power-law fit"))

    xbar = sum(xs) / length(xs)
    ybar = sum(ys) / length(ys)
    denom = sum((xx - xbar)^2 for xx in xs)
    denom > 0 || throw(ArgumentError("Degenerate x samples"))
    slope = sum((xs[k] - xbar) * (ys[k] - ybar) for k in eachindex(xs)) / denom
    intercept = ybar - slope * xbar
    return slope, intercept, count(idx)
end
