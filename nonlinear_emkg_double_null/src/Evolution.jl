function cell_metric_source(u, v, ep::EvolutionParams)
    p = ep.rn
    r = areal_radius(u, v, p)
    fp = 2p.M / r^2 - 2p.Q0^2 / r^3
    return -metric_f(u, v, p) * fp / (2r)
end

function evolve!(st::State, g::Grid, ep::EvolutionParams)
    nu, nv = size(g)
    for i in 1:nu-1
        for j in 1:nv-1
            step_cell!(st, g, ep, i, j)
        end
    end
    return st
end

function background_Au(u::Real, v::Real, g::Grid, p::RNParams)
    rarg = mrt_arg_u(u, p)
    return p.Q0 * sec(u)^2 / (2 * metric_F(rarg, p)) *
           (inv(areal_radius(u, first(g.v), p)) - inv(areal_radius(u, v, p)))
end

function background_Av(u::Real, v::Real, g::Grid, p::RNParams)
    rarg = mrt_arg_v(v, p)
    return p.Q0 * sec(v)^2 / (2 * metric_F(rarg, p)) *
           (inv(areal_radius(first(g.u), v, p)) - inv(areal_radius(u, v, p)))
end

function fill_background_potentials!(st::State, g::Grid, ep::EvolutionParams)
    for i in eachindex(g.u), j in eachindex(g.v)
        st.Au[i, j] = background_Au(g.u[i], g.v[j], g, ep.rn)
        st.Av[i, j] = background_Av(g.u[i], g.v[j], g, ep.rn)
        st.Q[i, j] = ep.rn.Q0
    end
    return st
end

function evolve_passive_scalar!(st::State, g::Grid, ep::EvolutionParams)
    fill_background_potentials!(st, g, ep)
    nu, nv = size(g)
    for i in 1:nu-1
        for j in 1:nv-1
            step_passive_scalar_cell!(st, g, ep, i, j)
        end
    end
    return st
end

function passive_charge_source(st::State, g::Grid, ep::EvolutionParams, i::Int, j::Int)
    du = g.u[i + 1] - g.u[i]
    dv = g.v[j + 1] - g.v[j]
    uc = (g.u[i] + g.u[i + 1]) / 2
    vc = (g.v[j] + g.v[j + 1]) / 2
    e = ep.scalar_charge
    f = metric_f(uc, vc, ep.rn)
    r = areal_radius(uc, vc, ep.rn)

    xi00, xi10, xi01, xi11 =
        st.xi[i, j], st.xi[i + 1, j], st.xi[i, j + 1], st.xi[i + 1, j + 1]
    pi00, pi10, pi01, pi11 =
        st.pi[i, j], st.pi[i + 1, j], st.pi[i, j + 1], st.pi[i + 1, j + 1]
    au00, au10, au01, au11 =
        st.Au[i, j], st.Au[i + 1, j], st.Au[i, j + 1], st.Au[i + 1, j + 1]
    av00, av10, av01, av11 =
        st.Av[i, j], st.Av[i + 1, j], st.Av[i, j + 1], st.Av[i + 1, j + 1]

    xic = (xi00 + xi10 + xi01 + xi11) / 4
    pic = (pi00 + pi10 + pi01 + pi11) / 4
    auc = (au00 + au10 + au01 + au11) / 4
    avc = (av00 + av10 + av01 + av11) / 4
    qc = ep.rn.Q0

    xiu = ((xi10 - xi00) + (xi11 - xi01)) / (2du)
    xiv = ((xi01 - xi00) + (xi11 - xi10)) / (2dv)
    piu = ((pi10 - pi00) + (pi11 - pi01)) / (2du)
    piv = ((pi01 - pi00) + (pi11 - pi10)) / (2dv)

    return 8pi * e * (-e * f * qc * (xic^2 + pic^2) / (2r^2) -
                      e * auc * (xic * xiv + pic * piv) +
                      e * avc * (xic * xiu + pic * piu) -
                      xiu * piv + xiv * piu)
end

function reconstruct_passive_charge!(st::State, g::Grid, ep::EvolutionParams)
    st.Q[firstindex(g.u), :] .= ep.rn.Q0
    st.Q[:, firstindex(g.v)] .= ep.rn.Q0
    nu, nv = size(g)
    for i in 1:nu-1
        for j in 1:nv-1
            du = g.u[i + 1] - g.u[i]
            dv = g.v[j + 1] - g.v[j]
            st.Q[i + 1, j + 1] =
                st.Q[i + 1, j] + st.Q[i, j + 1] - st.Q[i, j] +
                du * dv * passive_charge_source(st, g, ep, i, j)
        end
    end
    return st
end

function step_passive_scalar_cell!(st::State, g::Grid, ep::EvolutionParams,
                                   i::Int, j::Int)
    du = g.u[i + 1] - g.u[i]
    dv = g.v[j + 1] - g.v[j]
    uc = (g.u[i] + g.u[i + 1]) / 2
    vc = (g.v[j] + g.v[j + 1]) / 2
    e = ep.scalar_charge
    gg = cell_metric_source(uc, vc, ep)

    xi00, xi10, xi01 = st.xi[i, j], st.xi[i + 1, j], st.xi[i, j + 1]
    pi00, pi10, pi01 = st.pi[i, j], st.pi[i + 1, j], st.pi[i, j + 1]
    au00, au10, au01 = st.Au[i, j], st.Au[i + 1, j], st.Au[i, j + 1]
    av00, av10, av01 = st.Av[i, j], st.Av[i + 1, j], st.Av[i, j + 1]
    au11, av11 = st.Au[i + 1, j + 1], st.Av[i + 1, j + 1]

    T = promote_type(eltype(g.u), eltype(st.xi))
    residual_function! = function (residual, values)
        xi11, pi11 = values
        xic = (xi00 + xi10 + xi01 + xi11) / 4
        pic = (pi00 + pi10 + pi01 + pi11) / 4
        auc = (au00 + au10 + au01 + au11) / 4
        avc = (av00 + av10 + av01 + av11) / 4

        xiu = ((xi10 - xi00) + (xi11 - xi01)) / (2du)
        xiv = ((xi01 - xi00) + (xi11 - xi10)) / (2dv)
        piu = ((pi10 - pi00) + (pi11 - pi01)) / (2du)
        piv = ((pi01 - pi00) + (pi11 - pi10)) / (2dv)

        xiuv = gg * xic + e^2 * xic * auc * avc - e * (avc * piu + auc * piv)
        piuv = gg * pic + e^2 * pic * auc * avc + e * (avc * xiu + auc * xiv)

        residual[1] = xi11 - xi10 - xi01 + xi00 - du * dv * xiuv
        residual[2] = pi11 - pi10 - pi01 + pi00 - du * dv * piuv
        return residual
    end

    guess = T[xi10 + xi01 - xi00, pi10 + pi01 - pi00]
    base = zeros(T, 2)
    col1 = zeros(T, 2)
    col2 = zeros(T, 2)
    residual_function!(base, guess)
    trial = copy(guess)
    trial[1] += one(T)
    residual_function!(col1, trial)
    trial .= guess
    trial[2] += one(T)
    residual_function!(col2, trial)
    jacobian = hcat(col1 .- base, col2 .- base)
    correction = jacobian \ base
    st.xi[i + 1, j + 1], st.pi[i + 1, j + 1] = guess .- correction
    return st
end

function step_cell!(st::State, g::Grid, ep::EvolutionParams, i::Int, j::Int)
    du = g.u[i + 1] - g.u[i]
    dv = g.v[j + 1] - g.v[j]
    uc = (g.u[i] + g.u[i + 1]) / 2
    vc = (g.v[j] + g.v[j + 1]) / 2
    e = ep.scalar_charge
    f = metric_f(uc, vc, ep.rn)
    r = areal_radius(uc, vc, ep.rn)
    gg = cell_metric_source(uc, vc, ep)

    xi00, xi10, xi01 = st.xi[i, j], st.xi[i + 1, j], st.xi[i, j + 1]
    pi00, pi10, pi01 = st.pi[i, j], st.pi[i + 1, j], st.pi[i, j + 1]
    au00, au10, au01 = st.Au[i, j], st.Au[i + 1, j], st.Au[i, j + 1]
    av00, av10, av01 = st.Av[i, j], st.Av[i + 1, j], st.Av[i, j + 1]
    q00, q10, q01 = st.Q[i, j], st.Q[i + 1, j], st.Q[i, j + 1]

    xi11 = xi10 + xi01 - xi00
    pi11 = pi10 + pi01 - pi00
    au11 = au10 + (g.v[j + 1] - g.v[j]) * q00 * metric_f(uc, vc, ep.rn) / (2r^2)
    av11 = av01 - (g.u[i + 1] - g.u[i]) * q00 * metric_f(uc, vc, ep.rn) / (2r^2)
    q11 = q10 + q01 - q00

    workspace = NewtonCellWorkspace(promote_type(eltype(g.u), eltype(st.xi)), 5)
    workspace.values .= (xi11, pi11, au11, av11, q11)
    residual_function! = function (residual, values)
        xi11, pi11, au11, av11, q11 = values
        xic = (xi00 + xi10 + xi01 + xi11) / 4
        pic = (pi00 + pi10 + pi01 + pi11) / 4
        auc = (au00 + au10 + au01 + au11) / 4
        avc = (av00 + av10 + av01 + av11) / 4
        qc = (q00 + q10 + q01 + q11) / 4

        xiu = ((xi10 - xi00) + (xi11 - xi01)) / (2du)
        xiv = ((xi01 - xi00) + (xi11 - xi10)) / (2dv)
        piu = ((pi10 - pi00) + (pi11 - pi01)) / (2du)
        piv = ((pi01 - pi00) + (pi11 - pi10)) / (2dv)

        xiuv = gg * xic + e^2 * xic * auc * avc - e * (avc * piu + auc * piv)
        piuv = gg * pic + e^2 * pic * auc * avc + e * (avc * xiu + auc * xiv)

        faraday = qc * f / (2r^2)
        quv = 8pi * e * (-e * f * qc * (xic^2 + pic^2) / (2r^2) -
                          e * auc * (xic * xiv + pic * piv) +
                          e * avc * (xic * xiu + pic * piu) -
                          xiu * piv + xiv * piu)

        residual[1] = xi11 - xi10 - xi01 + xi00 - du * dv * xiuv
        residual[2] = pi11 - pi10 - pi01 + pi00 - du * dv * piuv
        residual[3] = au11 - au10 - dv * faraday
        residual[4] = av11 - av01 + du * faraday
        residual[5] = q11 - q10 - q01 + q00 - du * dv * quv
        return residual
    end
    solved = damped_newton_cell!(
        residual_function!,
        workspace;
        max_iterations=12,
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    if solved.converged
        xi11, pi11, au11, av11, q11 = workspace.values
    else
        xi11, pi11, au11, av11, q11 = workspace.values
    end

    st.xi[i + 1, j + 1] = xi11
    st.pi[i + 1, j + 1] = pi11
    st.Au[i + 1, j + 1] = au11
    st.Av[i + 1, j + 1] = av11
    st.Q[i + 1, j + 1] = q11
    return st
end
