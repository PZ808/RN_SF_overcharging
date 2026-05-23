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

function renormalized_hawking_mass(r, f, ru, rv, q)
    return r / 2 * (1 + 4 * ru * rv / f + q^2 / r^2)
end

function renormalized_hawking_mass_profile(lower::NLSlice, upper::NLSlice;
                                           rn_background::Union{Nothing,RNParams}=nothing)
    upper.v > lower.v || throw(ArgumentError("upper slice must have larger V"))
    upper_on_lower = lower.u == upper.u ? upper : interpolate_slice(upper, lower.u)

    ncell = length(lower.u) - 1
    dv = upper.v - lower.v
    u = [(lower.u[i] + lower.u[i + 1]) / 2 for i in 1:ncell]
    mass = similar(u)
    for i in 1:ncell
        du = lower.u[i + 1] - lower.u[i]
        r = corner_average(lower.r[i], lower.r[i + 1],
                           upper_on_lower.r[i], upper_on_lower.r[i + 1])
        logf = corner_average(lower.logf[i], lower.logf[i + 1],
                              upper_on_lower.logf[i], upper_on_lower.logf[i + 1])
        q = corner_average(lower.Q[i], lower.Q[i + 1],
                           upper_on_lower.Q[i], upper_on_lower.Q[i + 1])
        ru = corner_du(lower.r[i], lower.r[i + 1],
                       upper_on_lower.r[i], upper_on_lower.r[i + 1], du)
        rv = corner_dv(lower.r[i], lower.r[i + 1],
                       upper_on_lower.r[i], upper_on_lower.r[i + 1], dv)
        mass[i] = renormalized_hawking_mass(r, exp(logf), ru, rv, q)
        if !isnothing(rn_background)
            rb00 = mrt2013_areal_radius(lower.u[i], lower.v, rn_background)
            rb10 = mrt2013_areal_radius(lower.u[i + 1], lower.v, rn_background)
            rb01 = mrt2013_areal_radius(lower.u[i], upper.v, rn_background)
            rb11 = mrt2013_areal_radius(lower.u[i + 1], upper.v, rn_background)
            lfb00 = log(mrt2013_metric_f(lower.u[i], lower.v, rn_background))
            lfb10 = log(mrt2013_metric_f(lower.u[i + 1], lower.v, rn_background))
            lfb01 = log(mrt2013_metric_f(lower.u[i], upper.v, rn_background))
            lfb11 = log(mrt2013_metric_f(lower.u[i + 1], upper.v, rn_background))
            rb = corner_average(rb00, rb10, rb01, rb11)
            rub = corner_du(rb00, rb10, rb01, rb11, du)
            rvb = corner_dv(rb00, rb10, rb01, rb11, dv)
            lfb = corner_average(lfb00, lfb10, lfb01, lfb11)
            discrete_background_mass =
                renormalized_hawking_mass(rb, exp(lfb), rub, rvb, rn_background.Q0)
            mass[i] -= discrete_background_mass - rn_background.M
        end
    end
    return u, (lower.v + upper.v) / 2, mass
end

function bondi_mass_profile(st::AdaptiveNLState; target_v=150.0,
                            rn_background::Union{Nothing,RNParams}=nothing)
    length(st.slices) >= 2 ||
        throw(ArgumentError("Bondi-mass approximation needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return renormalized_hawking_mass_profile(st.slices[j], st.slices[j + 1];
                                             rn_background)
end

function outgoing_expansion_profile(lower::NLSlice, upper::NLSlice)
    upper.v > lower.v || throw(ArgumentError("upper slice must have larger V"))
    upper_on_lower = lower.u == upper.u ? upper : interpolate_slice(upper, lower.u)
    ncell = length(lower.u) - 1
    dv = upper.v - lower.v
    u = [(lower.u[i] + lower.u[i + 1]) / 2 for i in 1:ncell]
    rv = similar(u)
    for i in 1:ncell
        rv[i] = corner_dv(lower.r[i], lower.r[i + 1],
                          upper_on_lower.r[i], upper_on_lower.r[i + 1], dv)
    end
    return u, (lower.v + upper.v) / 2, rv
end

function outgoing_expansion_profile(st::AdaptiveNLState; target_v=150.0)
    length(st.slices) >= 2 ||
        throw(ArgumentError("expansion profile needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return outgoing_expansion_profile(st.slices[j], st.slices[j + 1])
end

function apparent_horizon_location(u::AbstractVector, rv::AbstractVector)
    length(u) == length(rv) ||
        throw(ArgumentError("U and r_V profiles must have the same length"))
    crossing = findfirst(value -> value <= 0, rv)
    isnothing(crossing) && return nothing
    crossing == firstindex(rv) && return u[crossing]
    i = crossing
    fraction = rv[i - 1] / (rv[i - 1] - rv[i])
    return u[i - 1] + fraction * (u[i] - u[i - 1])
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

function conformal_weight_s(eQ0)
    return real(0.5 + sqrt(complex(0.25 - eQ0^2)))
end

function horizon_charge_density_series(st::State, g::Grid, ep::EvolutionParams; i::Int=lastindex(g.u))
    p = ep.rn
    rplus, _ = horizons(p)
    vef = [ef_v_from_mrt(v, p) for v in g.v]
    rho = similar(vef)

    for j in eachindex(g.v)
        if i == firstindex(g.u)
            q_u_compact = (st.Q[i + 1, j] - st.Q[i, j]) / (g.u[i + 1] - g.u[i])
        elseif i == lastindex(g.u)
            q_u_compact = (st.Q[i, j] - st.Q[i - 1, j]) / (g.u[i] - g.u[i - 1])
        else
            q_u_compact = (st.Q[i + 1, j] - st.Q[i - 1, j]) / (g.u[i + 1] - g.u[i - 1])
        end

        if j == firstindex(g.v)
            q_vef = (st.Q[i, j + 1] - st.Q[i, j]) / (vef[j + 1] - vef[j])
        elseif j == lastindex(g.v)
            q_vef = (st.Q[i, j] - st.Q[i, j - 1]) / (vef[j] - vef[j - 1])
        else
            q_vef = (st.Q[i, j + 1] - st.Q[i, j - 1]) / (vef[j + 1] - vef[j - 1])
        end

        U = tan(g.u[i])
        du_dU = cos(g.u[i])^2
        if abs(p.Q0) ≈ p.M
            q_r = 0.5 * du_dU * q_u_compact
        else
            kappa = (horizons(p)[1] - horizons(p)[2]) / (2 * rplus^2)
            q_r = 0.5 * exp(-kappa * vef[j]) * du_dU * q_u_compact
        end
        rho[j] = (q_r + 0.5 * q_vef) / (4pi * rplus^2)
    end

    return vef, rho
end
