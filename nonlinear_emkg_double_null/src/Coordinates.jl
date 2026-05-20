struct RNParams{T<:Real}
    M::T
    Q0::T
end

function horizons(p::RNParams)
    disc = p.M^2 - p.Q0^2
    disc < 0 && throw(ArgumentError("RN naked singularity: |Q0| must be <= M"))
    root = sqrt(max(zero(disc), disc))
    return p.M + root, p.M - root
end

metric_F(r::Real, p::RNParams) = 1 - 2p.M / r + p.Q0^2 / r^2

function rstar(r::Real, p::RNParams)
    rplus, rminus = horizons(p)
    if abs(p.Q0) ≈ p.M
        return r + 2p.M * log(abs(r / p.M - 1)) - p.M^2 / (r - p.M)
    end
    return r +
           rplus^2 / (rplus - rminus) * log(abs(r / rplus - 1)) -
           rminus^2 / (rplus - rminus) * log(abs(r / rminus - 1))
end

function radius_from_rstar(target::Real, p::RNParams; tol=1.0e-12, maxiter=100)
    rplus, _ = horizons(p)
    lo = rplus * (1 + 1.0e-12)
    hi = max(rplus + one(float(rplus)), 2p.M + abs(target) + one(float(p.M)))
    while rstar(hi, p) < target
        hi *= 2
    end
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

function compact_v_from_ef_v(Vef::Real, p::RNParams)
    rplus, _ = horizons(p)
    rarg = radius_from_rstar(Vef / 2, p)
    return atan(2 * (rarg - rplus))
end

function compact_mrt_grid(p::RNParams; nu::Int=80, nv::Int=80, u0=-1.2, v0=0.0, u1=-1.0e-3, v1=pi / 2 - 1.0e-3)
    u = collect(range(u0, u1; length=nu))
    v = collect(range(v0, v1; length=nv))
    return Grid(u, v)
end

struct Grid{T<:Real}
    u::Vector{T}
    v::Vector{T}
end

Base.size(g::Grid) = (length(g.u), length(g.v))

mrt_arg_u(u::Real, p::RNParams) = horizons(p)[1] - tan(u) / 2
mrt_arg_v(v::Real, p::RNParams) = horizons(p)[1] + tan(v) / 2

function areal_radius(u::Real, v::Real, p::RNParams; tol=1.0e-12, maxiter=80)
    target = rstar(mrt_arg_u(u, p), p) + rstar(mrt_arg_v(v, p), p)
    rplus, _ = horizons(p)
    lo = max(rplus * (1 + 1.0e-12), 1.0e-12)
    hi = max(2rplus + tan(v) - tan(u) + 10p.M, lo * 2)
    while rstar(hi, p) < target
        hi *= 2
    end
    for _ in 1:maxiter
        mid = (lo + hi) / 2
        if rstar(mid, p) < target
            lo = mid
        else
            hi = mid
        end
        abs(hi - lo) <= tol * max(one(hi), abs(mid)) && return (lo + hi) / 2
    end
    return (lo + hi) / 2
end

function metric_f(u::Real, v::Real, p::RNParams)
    r = areal_radius(u, v, p)
    den = 2 * metric_F(mrt_arg_u(u, p), p) * metric_F(mrt_arg_v(v, p), p)
    return metric_F(r, p) * sec(u)^2 * sec(v)^2 / den
end

metric_ftilde(u::Real, v::Real, p::RNParams) = metric_f(u, v, p) * cos(u)^2 * cos(v)^2

renormalized_R(u::Real, v::Real, p::RNParams) = 2 * areal_radius(u, v, p) / (tan(v) - tan(u))
