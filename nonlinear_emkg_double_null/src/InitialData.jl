function gaussian_envelope(Vef::Real, ep::EvolutionParams)
    return ep.amplitude * exp(-((Vef - ep.center) / ep.width)^2)
end

function ef_v_from_mrt(v::Real, p::RNParams)
    return 2 * rstar(mrt_arg_v(v, p), p)
end

function initialize_state(g::Grid, ep::EvolutionParams)
    st = State(g)
    p = ep.rn
    st.Q[:, :] .= p.Q0

    for j in eachindex(g.v)
        Vef = ef_v_from_mrt(g.v[j], p)
        amp = gaussian_envelope(Vef, ep)
        st.xi[1, j] = amp * cos(ep.omega * Vef)
        st.pi[1, j] = -amp * sin(ep.omega * Vef)
    end

    for j in 2:length(g.v)
        dv = g.v[j] - g.v[j - 1]
        qsrc = 8pi * ep.scalar_charge * ep.omega * st.xi[1, j - 1]^2
        st.Q[1, j] = st.Q[1, j - 1] + dv * qsrc
        u = g.u[1]
        v = g.v[j - 1]
        r = areal_radius(u, v, p)
        st.Au[1, j] = st.Au[1, j - 1] + dv * metric_f(u, v, p) * st.Q[1, j - 1] / (2r^2)
    end

    for i in eachindex(g.u)
        st.Q[i, 1] = p.Q0
    end

    return st
end

