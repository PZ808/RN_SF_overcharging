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

function mrt2013_grid(; nu::Int=300, nv::Int=1200, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=200.0)
    return Grid(collect(range(U0, U1; length=nu)), collect(range(V0, V1; length=nv)))
end

function require_extreme_horizon_endpoint(p::RNParams)
    scale = max(abs(p.M), abs(p.Q0), one(float(p.M)))
    abs(abs(p.Q0) - p.M) <= 64 * eps(float(scale)) * scale ||
        throw(ArgumentError("the analytic U=0 MRT endpoint is currently implemented only for extreme RN"))
    return p
end

function mrt2013_interior_radius_from_rstar(target::Real, p::RNParams; tol=1.0e-12,
                                            maxiter=100)
    require_extreme_horizon_endpoint(p)
    rplus, _ = horizons(p)
    lo = eps(float(rplus)) * rplus
    hi = rplus * (1 - eps(float(rplus)))
    target >= rstar(lo, p) ||
        throw(ArgumentError("target lies outside the regular interior MRT branch"))
    for _ in 1:maxiter
        mid = (lo + hi) / 2
        if rstar(mid, p) < target
            lo = mid
        else
            hi = mid
        end
        abs(hi - lo) <= tol * max(one(mid), abs(mid)) && return (lo + hi) / 2
    end
    return (lo + hi) / 2
end

function mrt2013_areal_radius(U::Real, V::Real, p::RNParams; tol=1.0e-12, maxiter=100)
    rplus, _ = horizons(p)
    if iszero(U)
        require_extreme_horizon_endpoint(p)
        return rplus
    end
    target = rstar(rplus - U, p) + V / 2
    if U > 0
        return mrt2013_interior_radius_from_rstar(target, p; tol, maxiter)
    end
    return radius_from_rstar(target, p; tol, maxiter)
end

function mrt2013_metric_f(U::Real, V::Real, p::RNParams)
    if iszero(U)
        require_extreme_horizon_endpoint(p)
        return 2 * one(promote_type(typeof(U), typeof(V), typeof(p.M)))
    end
    rplus, _ = horizons(p)
    r = mrt2013_areal_radius(U, V, p)
    return 2 * metric_F(r, p) / metric_F(rplus - U, p)
end

function mrt2013_bump(x, xmin, xmax; alpha=4.0, amplitude=1.0e-2)
    (xmin < x < xmax) || return zero(x)
    width = xmax - xmin
    exponent = alpha * (inv(x - xmax) - inv(x - xmin) + 4 / width)
    return amplitude * exp(exponent)
end

function mrt2013_bump_derivative(x, xmin, xmax; alpha=4.0, amplitude=1.0e-2)
    (xmin < x < xmax) || return zero(x)
    phi = mrt2013_bump(x, xmin, xmax; alpha, amplitude)
    return alpha * phi * (-inv(x - xmax)^2 + inv(x - xmin)^2)
end

function trapezoidal_integral(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(ArgumentError("integration arrays must have equal length"))
    integral = zero(promote_type(eltype(x), eltype(y)))
    for i in firstindex(x):lastindex(x)-1
        integral += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2
    end
    return integral
end

function mrt2013_initial_bondi_mass(f0::Real; U0=-5.1)
    r0 = 1 - U0
    return r0 / 2 * (1 + inv(r0)^2 - 2 / f0 * (1 - inv(r0))^2)
end

"""
Reconstruct `r_V` on MRT's outgoing initial leg from Eq. (28).

This exposes how accurately a discrete initial `U` grid represents the
degenerate apparent-horizon condition `r_V(U=0,V=0)=0`.
"""
function mrt2013_initial_rv_profile(st::NLState, g::Grid)
    j0 = firstindex(g.v)
    n = length(g.u)
    rv = zeros(promote_type(eltype(g.u), eltype(st.r)), n)
    r0 = st.r[firstindex(g.u), j0]
    q0 = st.Q[firstindex(g.u), j0]
    rr_v = r0 / 2 * (1 - q0 / r0)^2
    rv[firstindex(g.u)] = rr_v / r0
    for i in firstindex(g.u)+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        r_left, r_right = st.r[i - 1, j0], st.r[i, j0]
        f_left, f_right = exp(st.logf[i - 1, j0]), exp(st.logf[i, j0])
        q_left, q_right = st.Q[i - 1, j0], st.Q[i, j0]
        source_left = f_left * (1 - q_left^2 / r_left^2)
        source_right = f_right * (1 - q_right^2 / r_right^2)
        rr_v -= du * (source_left + source_right) / 8
        rv[i] = rr_v / r_right
    end
    return rv
end

"""
Choose `f0` for MRT's degenerate-initial-apparent-horizon family.

This is the zero scalar-charge, `Q=M=1` outgoing-wave family of Sec. 3.1:
Eq. (28) is integrated until `U=0` and `f0` is selected so that
`r_V(0,0)=0`. The auxiliary quadrature is intentionally much finer than the
evolution grid because `f0 - 2` is of order `amplitude^2`.
"""
function mrt2013_degenerate_horizon_f0(ep::EvolutionParams; U0=-5.1, Uah=0.0,
                                       Uout=-5.0, Uin=0.9, alpha=4.0,
                                       quadrature_points::Int=20001)
    ep.scalar_charge == 0 ||
        throw(ArgumentError("the MRT Fig. 7 initial-data tuner is for scalar charge e=0"))
    ep.rn.M == 1 && ep.rn.Q0 == 1 ||
        throw(ArgumentError("the MRT Fig. 7 initial-data tuner assumes Q=M=1"))
    Uah == 0 ||
        throw(ArgumentError("the degenerate apparent horizon is fixed at U=0"))
    quadrature_points >= 3 ||
        throw(ArgumentError("quadrature_points must be at least three"))

    u = collect(range(U0, Uah; length=quadrature_points))
    r = 1 .- u
    phiu = [mrt2013_bump_derivative(U, Uout, Uin; alpha, amplitude=ep.amplitude)
            for U in u]

    scalar_integral = zeros(promote_type(eltype(u), typeof(ep.amplitude)), length(u))
    for i in 2:length(u)
        du = u[i] - u[i - 1]
        scalar_integral[i] = scalar_integral[i - 1] +
                             du * (r[i - 1] * phiu[i - 1]^2 + r[i] * phiu[i]^2) / 2
    end
    fshape = exp.(-scalar_integral ./ 4)
    vacuum_integrand = 1 .- inv.(r).^2
    integrand = vacuum_integrand .* fshape
    denominator = trapezoidal_integral(u, integrand)
    vacuum_denominator = trapezoidal_integral(u, vacuum_integrand)
    # Normalize by the same quadrature's vacuum value so the exact f0=2
    # background is preserved before extracting the O(amplitude^2) shift.
    return 2 * vacuum_denominator / denominator
end

function initialize_mrt2013_uncharged_ingoing!(st::NLState, g::Grid, ep::EvolutionParams;
                                               Vini=0.0, Vfin=5.9, alpha=4.0)
    fill!(st.r, zero(eltype(st.r)))
    fill!(st.logf, zero(eltype(st.logf)))
    fill!(st.phi_re, zero(eltype(st.phi_re)))
    fill!(st.phi_im, zero(eltype(st.phi_im)))
    fill!(st.Au, zero(eltype(st.Au)))
    fill!(st.Av, zero(eltype(st.Av)))
    fill!(st.Q, convert(eltype(st.Q), ep.rn.Q0))

    U0 = g.u[firstindex(g.u)]
    V0 = g.v[firstindex(g.v)]

    # Sigma_1: V=0, Eq. (21) and Appendix B: phi=0, f=2.
    j0 = firstindex(g.v)
    for i in eachindex(g.u)
        st.r[i, j0] = mrt2013_areal_radius(g.u[i], V0, ep.rn)
        st.logf[i, j0] = log(mrt2013_metric_f(g.u[i], V0, ep.rn))
    end

    # Sigma_2: U=U0, use the exact extremal RN radius from Appendix A.
    i0 = firstindex(g.u)
    for j in eachindex(g.v)
        st.r[i0, j] = mrt2013_areal_radius(U0, g.v[j], ep.rn)
    end

    for j in eachindex(g.v)
        st.phi_re[i0, j] = mrt2013_bump(g.v[j], Vini, Vfin; alpha, amplitude=ep.amplitude)
    end

    # Appendix B: f=2 on Sigma_1; determine f on Sigma_2 by constraint C1.
    st.logf[i0, j0] = log(2)
    y = (0.5 * metric_F(st.r[i0, j0], ep.rn)) / exp(st.logf[i0, j0])
    for j in j0+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        rprev = st.r[i0, j - 1]
        phiv = (st.phi_re[i0, j] - st.phi_re[i0, j - 1]) / dv
        # MRT Eq. (7): d_V(r_V/f) = -r phi_V^2 / (4f).
        # With y=r_V/f and f=r_V/y: y_V = -r phi_V^2 y / (4r_V).
        rv = 0.5 * metric_F(rprev, ep.rn)
        y -= dv * rprev * phiv^2 * y / (4 * max(abs(rv), eps(typeof(rv))))
        rv_current = 0.5 * metric_F(st.r[i0, j], ep.rn)
        st.logf[i0, j] = log(abs(rv_current / y))
    end

    return st
end

function initialize_mrt2013_outgoing_wave!(st::NLState, g::Grid, ep::EvolutionParams;
                                           Uout=-5.0, Uin=0.9, alpha=4.0, f0=2.0)
    fill!(st.r, zero(eltype(st.r)))
    fill!(st.logf, zero(eltype(st.logf)))
    fill!(st.phi_re, zero(eltype(st.phi_re)))
    fill!(st.phi_im, zero(eltype(st.phi_im)))
    fill!(st.Au, zero(eltype(st.Au)))
    fill!(st.Av, zero(eltype(st.Av)))
    fill!(st.Q, convert(eltype(st.Q), ep.rn.Q0))

    i0 = firstindex(g.u)
    j0 = firstindex(g.v)
    U0 = g.u[i0]
    V0 = g.v[j0]

    # MRT Eqs. (20), (23), and (24): outgoing wavepacket on Sigma_1.
    for i in eachindex(g.u)
        st.r[i, j0] = mrt2013_areal_radius(g.u[i], V0, ep.rn)
        st.phi_re[i, j0] = mrt2013_bump(g.u[i], Uout, Uin;
                                        alpha, amplitude=ep.amplitude)
    end

    # MRT Eqs. (20) and (23): trivial scalar data on Sigma_2.
    for j in eachindex(g.v)
        st.r[i0, j] = mrt2013_areal_radius(U0, g.v[j], ep.rn)
    end

    # MRT Eq. (25): solve C2 on Sigma_1.
    st.logf[i0, j0] = log(f0)
    integral = zero(eltype(st.logf))
    for i in i0+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        left_source = st.r[i - 1, j0] *
                      mrt2013_bump_derivative(g.u[i - 1], Uout, Uin;
                                              alpha, amplitude=ep.amplitude)^2
        right_source = st.r[i, j0] *
                       mrt2013_bump_derivative(g.u[i], Uout, Uin;
                                               alpha, amplitude=ep.amplitude)^2
        integral += du * (left_source + right_source) / 2
        st.logf[i, j0] = log(f0) - integral / 4
    end

    # MRT Eq. (26): solve C1 on Sigma_2, where phi=0.
    normalization = f0 / mrt2013_metric_f(U0, V0, ep.rn)
    for j in eachindex(g.v)
        st.logf[i0, j] = log(normalization * mrt2013_metric_f(U0, g.v[j], ep.rn))
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
        st.Au[i0, j] = st.Au[i0, j - 1] - dv * st.Q[i0, j - 1] * f / (2 * st.r[i0, j - 1]^2)
    end
    for i in i0+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        f = exp(st.logf[i - 1, j0])
        st.Av[i, j0] = st.Av[i - 1, j0] + du * st.Q[i - 1, j0] * f / (2 * st.r[i - 1, j0]^2)
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
        y -= dv * st.r[i, j - 1] * source.Tvv * y / (8 * rv)
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
        y -= du * st.r[i - 1, j] * source.Tuu * y / (8 * ru)
        st.logf[i, j] = log(abs(ru / y))
    end
    return st
end

function evolve_nonlinear!(st::NLState, g::Grid, ep::EvolutionParams; iterations::Int=5,
                           subtract_rn_background::Bool=false)
    nu, nv = size(g)
    for i in 1:nu-1
        for j in 1:nv-1
            step_nonlinear_cell!(st, g, ep, i, j; iterations, subtract_rn_background)
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
    # MRT Eq. (4): (log f)_UV = f/(2r^2) + 2 r_U r_V/r^2
    #                         - Q^2 f/r^4 - (1/2) phi_U phi_V.
    # `source.Tthth` includes Maxwell angular stress, while the MRT equation
    # already carries the Q^2 term explicitly, so subtract that part here.
    scalar_tthth = source.Tthth - 2 * q^2 / r^2
    scalar_uv_source = scalar_tthth * f / (8r^2)
    ruv = (-ru * rv - f * (1 - q^2 / r^2) / 4) / r
    logfuv = f / (2r^2) + 2 * ru * rv / r^2 - q^2 * f / r^4 - scalar_uv_source
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
    # With F=dA and F_UV = A_V,U - A_U,V = Q f/r^2, quasi-Lorenz
    # A_U,V + A_V,U = 0 gives the following split.
    auv = -q * f / (2r^2)
    avu = q * f / (2r^2)

    # Constraint-derived charge evolution. Coefficients should be rechecked
    # against the chosen 4pi/action normalization before production runs.
    q_uv = zero(r)
    q_v_constraint = -4pi * r^2 * source.Jv
    q_u_constraint = 4pi * r^2 * source.Ju
    return auv, avu, q_uv, q_u_constraint, q_v_constraint
end

function rn_background_update_defect(g::Grid, ep::EvolutionParams, i::Int, j::Int, du, dv)
    u0, u1 = g.u[i], g.u[i + 1]
    v0, v1 = g.v[j], g.v[j + 1]
    p = ep.rn

    rb00 = mrt2013_areal_radius(u0, v0, p)
    rb10 = mrt2013_areal_radius(u1, v0, p)
    rb01 = mrt2013_areal_radius(u0, v1, p)
    rb11 = mrt2013_areal_radius(u1, v1, p)
    lfb00 = log(mrt2013_metric_f(u0, v0, p))
    lfb10 = log(mrt2013_metric_f(u1, v0, p))
    lfb01 = log(mrt2013_metric_f(u0, v1, p))
    lfb11 = log(mrt2013_metric_f(u1, v1, p))

    rb = corner_average(rb00, rb10, rb01, rb11)
    lfb = corner_average(lfb00, lfb10, lfb01, lfb11)
    fb = exp(lfb)
    rub = corner_du(rb00, rb10, rb01, rb11, du)
    rvb = corner_dv(rb00, rb10, rb01, rb11, dv)
    source = stress_energy(rb, fb, p.Q0, zero(rb), zero(rb), zero(rb), zero(rb),
                           zero(rb), zero(rb), zero(rb), zero(rb), zero(rb))
    ruvb, lfuvb = metric_rhs(rb, fb, rub, rvb, p.Q0, source)

    r_defect = rb10 + rb01 - rb00 + du * dv * ruvb - rb11
    lf_defect = lfb10 + lfb01 - lfb00 + du * dv * lfuvb - lfb11
    return r_defect, lf_defect
end

function step_nonlinear_cell!(st::NLState, g::Grid, ep::EvolutionParams, i::Int, j::Int;
                              iterations::Int=5, subtract_rn_background::Bool=false)
    du = g.u[i + 1] - g.u[i]
    dv = g.v[j + 1] - g.v[j]
    e = ep.scalar_charge
    r_defect, lf_defect = subtract_rn_background ?
                           rn_background_update_defect(g, ep, i, j, du, dv) :
                           (zero(du), zero(du))

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

        r11 = r10 + r01 - r00 + du * dv * ruv - r_defect
        lf11 = lf10 + lf01 - lf00 + du * dv * lfuv - lf_defect
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
