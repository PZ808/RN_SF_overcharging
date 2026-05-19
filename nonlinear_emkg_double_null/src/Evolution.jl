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
    au11 = au10 + au01 - au00
    av11 = av10 + av01 - av00
    q11 = q10 + q01 - q00

    for _ in 1:4
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

        auv = qc * f / (2r^2)
        avu = -qc * f / (2r^2)
        quv = 8pi * e * (-e * f * qc * (xic^2 + pic^2) / (2r^2) -
                          e * auc * (xic * xiv + pic * piv) +
                          e * avc * (xic * xiu + pic * piu) -
                          xiu * piv + xiv * piu)

        xi11 = xi10 + xi01 - xi00 + du * dv * xiuv
        pi11 = pi10 + pi01 - pi00 + du * dv * piuv
        au11 = au10 + au01 - au00 + du * dv * auv
        av11 = av10 + av01 - av00 + du * dv * avu
        q11 = q10 + q01 - q00 + du * dv * quv
    end

    st.xi[i + 1, j + 1] = xi11
    st.pi[i + 1, j + 1] = pi11
    st.Au[i + 1, j + 1] = au11
    st.Av[i + 1, j + 1] = av11
    st.Q[i + 1, j + 1] = q11
    return st
end

