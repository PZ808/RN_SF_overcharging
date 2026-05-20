"""
State for the fully nonlinear spherical Einstein-Maxwell-charged-scalar system.

Convention used in this file:

    ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2

This is the Murata-Reall-Tanahashi 2013 convention. It differs by a factor of
two from the fixed-background SED convention used in `Evolution.jl`.
"""
struct NLState{T<:Real}
    r::Matrix{T}
    logf::Matrix{T}
    phi_re::Matrix{T}
    phi_im::Matrix{T}
    Au::Matrix{T}
    Av::Matrix{T}
    Q::Matrix{T}
end

function NLState(::Type{T}, nu::Integer, nv::Integer) where {T<:Real}
    z = () -> zeros(T, nu, nv)
    return NLState(z(), z(), z(), z(), z(), z(), z())
end

NLState(g::Grid{T}) where {T<:Real} = NLState(T, length(g.u), length(g.v))

"RN background metric coefficient in the MRT 2013 convention ds^2=-f du dv+r^2 dOmega^2."
mrt2013_background_f(u::Real, v::Real, p::RNParams) = 2 * metric_f(u, v, p)

function mrt2013_grid(; nu::Int=300, nv::Int=1200, U0=-5.1, V0=0.0, U1=0.95, V1=200.0)
    return Grid(collect(range(U0, U1; length=nu)), collect(range(V0, V1; length=nv)))
end

function mrt2013_bump(x, xmin, xmax; alpha=4.0, amplitude=1.0e-2)
    (xmin < x < xmax) || return zero(x)
    mid = (xmin + xmax) / 2
    halfwidth = (xmax - xmin) / 2
    z = (x - mid) / halfwidth
    return amplitude * exp(-alpha * z^2 / (1 - z^2))
end

function initialize_mrt2013_uncharged_ingoing!(st::NLState, g::Grid, ep::EvolutionParams;
                                               Vini=0.0, Vfin=5.9, alpha=4.0)
    fill!(st.r, zero(eltype(st.r)))
    fill!(st.logf, zero(eltype(st.logf)))
    fill!(st.phi_re, zero(eltype(st.phi_re)))
    fill!(st.phi_im, zero(eltype(st.phi_im)))
    fill!(st.Au, zero(eltype(st.Au)))
    fill!(st.Av, zero(eltype(st.Av)))
    fill!(st.Q, one(eltype(st.Q)))

    U0 = g.u[firstindex(g.u)]
    V0 = g.v[firstindex(g.v)]
    r0 = 1 - U0

    # Sigma_1: V=0, Eq. (21) and Appendix B: phi=0, f=2.
    j0 = firstindex(g.v)
    for i in eachindex(g.u)
        st.r[i, j0] = 1 - g.u[i]
        st.logf[i, j0] = log(2)
    end

    # Sigma_2: U=U0, use Eq. (22) to integrate r_V.
    i0 = firstindex(g.u)
    st.r[i0, j0] = r0
    for j in j0+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        rprev = st.r[i0, j - 1]
        rv = 0.5 * (1 - inv(rprev))^2
        st.r[i0, j] = rprev + dv * rv
    end

    for j in eachindex(g.v)
        st.phi_re[i0, j] = mrt2013_bump(g.v[j], Vini, Vfin; alpha, amplitude=ep.amplitude)
    end

    # Appendix B: f=2 on Sigma_1; determine f on Sigma_2 by constraint C1.
    st.logf[i0, j0] = log(2)
    y = (0.5 * (1 - inv(st.r[i0, j0]))^2) / 2
    for j in j0+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        rprev = st.r[i0, j - 1]
        phiv = (st.phi_re[i0, j] - st.phi_re[i0, j - 1]) / dv
        # C1 => d_V(r_V/f) = r phi_V^2 / (4 f).
        # With y=r_V/f and f=r_V/y: y_V = r phi_V^2 y / (4 r_V).
        rv = 0.5 * (1 - inv(rprev))^2
        y += dv * rprev * phiv^2 * y / (4 * max(abs(rv), eps(typeof(rv))))
        st.logf[i0, j] = log(abs(rv / y))
    end

    return st
end

function initialize_nonlinear_state(g::Grid, ep::EvolutionParams)
    st = NLState(g)
    p = ep.rn
    for i in eachindex(g.u), j in eachindex(g.v)
        st.r[i, j] = areal_radius(g.u[i], g.v[j], p)
        st.logf[i, j] = log(mrt2013_background_f(g.u[i], g.v[j], p))
        st.Q[i, j] = p.Q0
    end

    seed_ingoing_scalar!(st, g, ep)
    seed_quasilorenz_potential!(st, g, ep)
    solve_initial_logf_constraints!(st, g, ep)
    return st
end

function seed_ingoing_scalar!(st::NLState, g::Grid, ep::EvolutionParams)
    p = ep.rn
    i = firstindex(g.u)
    for j in eachindex(g.v)
        Vef = ef_v_from_mrt(g.v[j], p)
        amp = gaussian_envelope(Vef, ep)
        st.phi_re[i, j] = amp * cos(ep.omega * Vef) / max(st.r[i, j], eps(eltype(g.v)))
        st.phi_im[i, j] = -amp * sin(ep.omega * Vef) / max(st.r[i, j], eps(eltype(g.v)))
    end
    return st
end

function seed_quasilorenz_potential!(st::NLState, g::Grid, ep::EvolutionParams)
    i0 = firstindex(g.u)
    j0 = firstindex(g.v)
    for j in j0+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        f = exp(st.logf[i0, j - 1])
        st.Au[i0, j] = st.Au[i0, j - 1] + dv * st.Q[i0, j - 1] * f / (2 * st.r[i0, j - 1]^2)
    end
    for i in i0+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        f = exp(st.logf[i - 1, j0])
        st.Av[i, j0] = st.Av[i - 1, j0] - du * st.Q[i - 1, j0] * f / (2 * st.r[i - 1, j0]^2)
    end
    return st
end

function solve_initial_logf_constraints!(st::NLState, g::Grid, ep::EvolutionParams)
    solve_outgoing_leg_logf_constraint!(st, g, ep)
    solve_ingoing_leg_logf_constraint!(st, g, ep)
    return st
end

function solve_outgoing_leg_logf_constraint!(st::NLState, g::Grid, ep::EvolutionParams)
    i = firstindex(g.u)
    length(g.v) < 2 && return st

    dv1 = g.v[2] - g.v[1]
    rv0 = (st.r[i, 2] - st.r[i, 1]) / dv1
    y = rv0 / exp(st.logf[i, 1])

    for j in firstindex(g.v)+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        rv = (st.r[i, j] - st.r[i, j - 1]) / dv
        phiv_re = (st.phi_re[i, j] - st.phi_re[i, j - 1]) / dv
        phiv_im = (st.phi_im[i, j] - st.phi_im[i, j - 1]) / dv
        fprev = exp(st.logf[i, j - 1])
        source = stress_energy(st.r[i, j - 1], fprev, st.Q[i, j - 1],
                               st.phi_re[i, j - 1], st.phi_im[i, j - 1],
                               zero(rv), phiv_re, zero(rv), phiv_im,
                               st.Au[i, j - 1], st.Av[i, j - 1], ep.scalar_charge)
        y += dv * source.Tvv * y^2 / (8 * max(rv^2, eps(typeof(rv))))
        st.logf[i, j] = log(abs(rv / y))
    end
    return st
end

function solve_ingoing_leg_logf_constraint!(st::NLState, g::Grid, ep::EvolutionParams)
    j = firstindex(g.v)
    length(g.u) < 2 && return st

    du1 = g.u[2] - g.u[1]
    ru0 = (st.r[2, j] - st.r[1, j]) / du1
    y = ru0 / exp(st.logf[1, j])

    for i in firstindex(g.u)+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        ru = (st.r[i, j] - st.r[i - 1, j]) / du
        phiu_re = (st.phi_re[i, j] - st.phi_re[i - 1, j]) / du
        phiu_im = (st.phi_im[i, j] - st.phi_im[i - 1, j]) / du
        fprev = exp(st.logf[i - 1, j])
        source = stress_energy(st.r[i - 1, j], fprev, st.Q[i - 1, j],
                               st.phi_re[i - 1, j], st.phi_im[i - 1, j],
                               phiu_re, zero(ru), phiu_im, zero(ru),
                               st.Au[i - 1, j], st.Av[i - 1, j], ep.scalar_charge)
        y += du * source.Tuu * y^2 / (8 * max(ru^2, eps(typeof(ru))))
        st.logf[i, j] = log(abs(ru / y))
    end
    return st
end

function evolve_nonlinear!(st::NLState, g::Grid, ep::EvolutionParams; iterations::Int=5)
    nu, nv = size(g)
    for i in 1:nu-1
        for j in 1:nv-1
            step_nonlinear_cell!(st, g, ep, i, j; iterations)
        end
    end
    return st
end

function corner_average(a00, a10, a01, a11)
    return (a00 + a10 + a01 + a11) / 4
end

function corner_du(a00, a10, a01, a11, du)
    return ((a10 - a00) + (a11 - a01)) / (2du)
end

function corner_dv(a00, a10, a01, a11, dv)
    return ((a01 - a00) + (a11 - a10)) / (2dv)
end

function metric_rhs(r, f, ru, rv, q, source::StressEnergyComponents)
    # These are the MRT uncharged-real-scalar equations with charged matter
    # inserted in the same slots. The coefficients multiplying T_ab still need
    # a final normalization pass against the chosen action.
    scalar_uv_source = source.Tthth * f / (4r^2)
    ruv = (-ru * rv - f * (1 - q^2 / r^2) / 4) / r
    logfuv = -2 * ru * rv / r^2 - f / (2r^2) + q^2 * f / r^4 - scalar_uv_source
    return ruv, logfuv
end

function charged_scalar_rhs(r, ru, rv, phi_re_u, phi_re_v, phi_im_u, phi_im_v,
                            phi_re, phi_im, Au, Av, e)
    # Real/imaginary split of D_a D^a phi + (2/r) r_,a D^a phi = 0.
    # This assumes quasi-Lorenz Au,v + Av,u = 0 at the cell level.
    re_uv = -(ru * phi_re_v + rv * phi_re_u) / r -
            e * (Av * phi_im_u + Au * phi_im_v) -
            e * (ru * Av + rv * Au) * phi_im / r +
            e^2 * Au * Av * phi_re
    im_uv = -(ru * phi_im_v + rv * phi_im_u) / r +
            e * (Av * phi_re_u + Au * phi_re_v) +
            e * (ru * Av + rv * Au) * phi_re / r +
            e^2 * Au * Av * phi_im
    return re_uv, im_uv
end

function maxwell_rhs(r, f, q, source::StressEnergyComponents)
    auv = q * f / (2r^2)
    avu = -q * f / (2r^2)

    # Constraint-derived charge evolution. Coefficients should be rechecked
    # against the chosen 4pi/action normalization before production runs.
    q_uv = zero(r)
    q_v_constraint = -4pi * r^2 * source.Jv
    q_u_constraint = 4pi * r^2 * source.Ju
    return auv, avu, q_uv, q_u_constraint, q_v_constraint
end

function step_nonlinear_cell!(st::NLState, g::Grid, ep::EvolutionParams, i::Int, j::Int; iterations::Int=5)
    du = g.u[i + 1] - g.u[i]
    dv = g.v[j + 1] - g.v[j]
    e = ep.scalar_charge

    r00, r10, r01 = st.r[i, j], st.r[i + 1, j], st.r[i, j + 1]
    lf00, lf10, lf01 = st.logf[i, j], st.logf[i + 1, j], st.logf[i, j + 1]
    pr00, pr10, pr01 = st.phi_re[i, j], st.phi_re[i + 1, j], st.phi_re[i, j + 1]
    pi00, pi10, pi01 = st.phi_im[i, j], st.phi_im[i + 1, j], st.phi_im[i, j + 1]
    au00, au10, au01 = st.Au[i, j], st.Au[i + 1, j], st.Au[i, j + 1]
    av00, av10, av01 = st.Av[i, j], st.Av[i + 1, j], st.Av[i, j + 1]
    q00, q10, q01 = st.Q[i, j], st.Q[i + 1, j], st.Q[i, j + 1]

    r11 = r10 + r01 - r00
    lf11 = lf10 + lf01 - lf00
    pr11 = pr10 + pr01 - pr00
    pi11 = pi10 + pi01 - pi00
    au11 = au10 + au01 - au00
    av11 = av10 + av01 - av00
    q11 = q10 + q01 - q00

    for _ in 1:iterations
        r = corner_average(r00, r10, r01, r11)
        lf = corner_average(lf00, lf10, lf01, lf11)
        f = exp(lf)
        pr = corner_average(pr00, pr10, pr01, pr11)
        pii = corner_average(pi00, pi10, pi01, pi11)
        au = corner_average(au00, au10, au01, au11)
        av = corner_average(av00, av10, av01, av11)
        q = corner_average(q00, q10, q01, q11)

        ru = corner_du(r00, r10, r01, r11, du)
        rv = corner_dv(r00, r10, r01, r11, dv)
        pru = corner_du(pr00, pr10, pr01, pr11, du)
        prv = corner_dv(pr00, pr10, pr01, pr11, dv)
        piu = corner_du(pi00, pi10, pi01, pi11, du)
        piv = corner_dv(pi00, pi10, pi01, pi11, dv)

        source = stress_energy(r, f, q, pr, pii, pru, prv, piu, piv, au, av, e)

        ruv, lfuv = metric_rhs(r, f, ru, rv, q, source)
        pruv, piuv = charged_scalar_rhs(r, ru, rv, pru, prv, piu, piv, pr, pii, au, av, e)
        auv, avu, _, quc, qvc = maxwell_rhs(r, f, q, source)

        r11 = r10 + r01 - r00 + du * dv * ruv
        lf11 = lf10 + lf01 - lf00 + du * dv * lfuv
        pr11 = pr10 + pr01 - pr00 + du * dv * pruv
        pi11 = pi10 + pi01 - pi00 + du * dv * piuv
        au11 = au10 + au01 - au00 + du * dv * auv
        av11 = av10 + av01 - av00 + du * dv * avu
        q11 = (q10 + dv * qvc + q01 + du * quc) / 2
    end

    st.r[i + 1, j + 1] = r11
    st.logf[i + 1, j + 1] = lf11
    st.phi_re[i + 1, j + 1] = pr11
    st.phi_im[i + 1, j + 1] = pi11
    st.Au[i + 1, j + 1] = au11
    st.Av[i + 1, j + 1] = av11
    st.Q[i + 1, j + 1] = q11
    return st
end
