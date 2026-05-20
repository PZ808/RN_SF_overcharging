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

    phase = zero(eltype(g.v))
    prev_vef = ef_v_from_mrt(g.v[firstindex(g.v)], p)

    for j in eachindex(g.v)
        Vef = ef_v_from_mrt(g.v[j], p)
        if j > firstindex(g.v)
            dVef = Vef - prev_vef
            rprev = areal_radius(g.u[1], g.v[j - 1], p)
            # Gauge rotation from Lorenz-style monochromatic data to the
            # quasi-Lorenz boundary condition Av|NB=0.  In EF coordinates
            # A_V^Lorenz = -Q/(2r), so alpha_V = Q/(2r).
            phase += ep.scalar_charge * p.Q0 * dVef / (2rprev)
            prev_vef = Vef
        end
        amp = gaussian_envelope(Vef, ep)
        theta = phase - ep.omega * Vef
        st.xi[1, j] = amp * cos(theta)
        st.pi[1, j] = amp * sin(theta)
    end

    for j in 2:length(g.v)
        dVef = ef_v_from_mrt(g.v[j], p) - ef_v_from_mrt(g.v[j - 1], p)
        amp2 = st.xi[1, j - 1]^2 + st.pi[1, j - 1]^2
        qsrc = 8pi * ep.scalar_charge * ep.omega * amp2
        st.Q[1, j] = st.Q[1, j - 1] + dVef * qsrc
        u = g.u[1]
        v = g.v[j - 1]
        r = areal_radius(u, v, p)
        dv = g.v[j] - g.v[j - 1]
        st.Au[1, j] = st.Au[1, j - 1] + dv * metric_f(u, v, p) * st.Q[1, j - 1] / (2r^2)
    end

    for i in eachindex(g.u)
        st.Q[i, 1] = p.Q0
    end

    for i in 2:length(g.u)
        du = g.u[i] - g.u[i - 1]
        u = g.u[i - 1]
        v = g.v[1]
        r = areal_radius(u, v, p)
        st.Av[i, 1] = st.Av[i - 1, 1] - du * metric_f(u, v, p) * st.Q[i - 1, 1] / (2r^2)
    end

    return st
end
