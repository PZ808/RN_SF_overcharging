struct EvolutionParams{T<:Real}
    rn::RNParams{T}
    scalar_charge::T
    amplitude::T
    omega::T
    center::T
    width::T
end

function EvolutionParams(; rn=RNParams(1.0, 0.8), scalar_charge=0.2, amplitude=1.0e-3,
                         omega=0.1, center=12.0, width=2.0)
    T = promote_type(typeof(rn.M), typeof(rn.Q0), typeof(scalar_charge),
                     typeof(amplitude), typeof(omega), typeof(center), typeof(width))
    return EvolutionParams{T}(
        RNParams(convert(T, rn.M), convert(T, rn.Q0)),
        convert(T, scalar_charge),
        convert(T, amplitude),
        convert(T, omega),
        convert(T, center),
        convert(T, width),
    )
end

struct State{T<:Real}
    xi::Matrix{T}
    pi::Matrix{T}
    Au::Matrix{T}
    Av::Matrix{T}
    Q::Matrix{T}
end

function State(::Type{T}, nu::Integer, nv::Integer) where {T<:Real}
    return State(zeros(T, nu, nv), zeros(T, nu, nv), zeros(T, nu, nv), zeros(T, nu, nv), zeros(T, nu, nv))
end

State(g::Grid{T}) where {T<:Real} = State(T, length(g.u), length(g.v))

function current_u(xi, pii, xi_u, pi_u, Au, e)
    return -2e * (xi_u * pii - xi * pi_u + e * Au * (xi^2 + pii^2))
end

function current_v(xi, pii, xi_v, pi_v, Av, e)
    return -2e * (xi_v * pii - xi * pi_v + e * Av * (xi^2 + pii^2))
end

function copy_cell_values!(dst::State, src::State, i::Int, j::Int)
    dst.xi[i, j] = src.xi[i, j]
    dst.pi[i, j] = src.pi[i, j]
    dst.Au[i, j] = src.Au[i, j]
    dst.Av[i, j] = src.Av[i, j]
    dst.Q[i, j] = src.Q[i, j]
    return dst
end
