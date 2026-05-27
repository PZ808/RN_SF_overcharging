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

"""
Evaluate the MRT Eq. (14) `U`-flux law for an uncharged scalar:

    varpi_U = -r^2 r_V |phi_U|^2 / (2f).

The returned `mass_u` is differentiated from the reconstructed Hawking-mass
profile; `expected_mass_u` is evaluated directly from the evolved fields.
Passing `rn_background` applies the same diagnostic stencil-defect correction
as `renormalized_hawking_mass_profile`.
"""
function uncharged_mass_flux_u_profile(lower::NLSlice, upper::NLSlice;
                                       rn_background::Union{Nothing,RNParams}=nothing)
    u, v, mass = renormalized_hawking_mass_profile(lower, upper; rn_background)
    length(mass) >= 2 ||
        throw(ArgumentError("mass-flux check needs at least three U grid points"))
    upper_on_lower = lower.u == upper.u ? upper : interpolate_slice(upper, lower.u)
    expected_cells = similar(mass)
    for i in eachindex(mass)
        du = lower.u[i + 1] - lower.u[i]
        dv = upper.v - lower.v
        r = corner_average(lower.r[i], lower.r[i + 1],
                           upper_on_lower.r[i], upper_on_lower.r[i + 1])
        f = exp(corner_average(lower.logf[i], lower.logf[i + 1],
                               upper_on_lower.logf[i], upper_on_lower.logf[i + 1]))
        rv = corner_dv(lower.r[i], lower.r[i + 1],
                       upper_on_lower.r[i], upper_on_lower.r[i + 1], dv)
        phiu_re = corner_du(lower.phi_re[i], lower.phi_re[i + 1],
                            upper_on_lower.phi_re[i], upper_on_lower.phi_re[i + 1], du)
        phiu_im = corner_du(lower.phi_im[i], lower.phi_im[i + 1],
                            upper_on_lower.phi_im[i], upper_on_lower.phi_im[i + 1], du)
        expected_cells[i] = -r^2 * rv * (phiu_re^2 + phiu_im^2) / (2f)
    end
    centered_u = [(u[i] + u[i + 1]) / 2 for i in firstindex(u):lastindex(u)-1]
    mass_u = diff(mass) ./ diff(u)
    expected_mass_u = [(expected_cells[i] + expected_cells[i + 1]) / 2
                       for i in firstindex(expected_cells):lastindex(expected_cells)-1]
    return centered_u, v, mass_u, expected_mass_u, mass_u .- expected_mass_u
end

function uncharged_mass_flux_u_profile(st::AdaptiveNLState; target_v=150.0,
                                       rn_background::Union{Nothing,RNParams}=nothing)
    length(st.slices) >= 2 ||
        throw(ArgumentError("mass-flux check needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return uncharged_mass_flux_u_profile(st.slices[j], st.slices[j + 1]; rn_background)
end

"""
Construct an uncharged scalar mass profile by integrating MRT Eq. (14) in `U`.

The integrated profile is anchored to the geometrically reconstructed mass in
the first outer cell. Away from that anchor, it is a conservative alternative
to differentiating `r_U r_V/f`; a disagreement with `geometric_mass` measures
the accumulated failure of the geometric mass to satisfy the flux law.
"""
function uncharged_flux_integrated_mass_profile(lower::NLSlice, upper::NLSlice;
                                                rn_background::Union{Nothing,RNParams}=nothing)
    u, v, geometric_mass = renormalized_hawking_mass_profile(lower, upper; rn_background)
    _, _, _, expected_mass_u, _ =
        uncharged_mass_flux_u_profile(lower, upper; rn_background)
    flux_mass = similar(geometric_mass)
    flux_mass[1] = geometric_mass[1]
    for i in 2:length(flux_mass)
        flux_mass[i] = flux_mass[i - 1] + (u[i] - u[i - 1]) * expected_mass_u[i - 1]
    end
    return u, v, geometric_mass, flux_mass, geometric_mass .- flux_mass
end

function uncharged_flux_integrated_mass_profile(st::AdaptiveNLState; target_v=150.0,
                                                rn_background::Union{Nothing,RNParams}=nothing)
    length(st.slices) >= 2 ||
        throw(ArgumentError("flux-integrated mass needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return uncharged_flux_integrated_mass_profile(st.slices[j], st.slices[j + 1];
                                                  rn_background)
end

function charged_u_flux_cells(lower::NLSlice, upper::NLSlice, scalar_charge;
                              reduced_scalar::Bool=false)
    upper.v > lower.v || throw(ArgumentError("upper slice must have larger V"))
    upper_on_lower = lower.u == upper.u ? upper : interpolate_slice(upper, lower.u)

    ncell = length(lower.u) - 1
    dv = upper.v - lower.v
    u = [(lower.u[i] + lower.u[i + 1]) / 2 for i in 1:ncell]
    q = similar(u)
    expected_q_u = similar(u)
    expected_mass_u = similar(u)
    for i in eachindex(u)
        du = lower.u[i + 1] - lower.u[i]
        r = corner_average(lower.r[i], lower.r[i + 1],
                           upper_on_lower.r[i], upper_on_lower.r[i + 1])
        f = exp(corner_average(lower.logf[i], lower.logf[i + 1],
                               upper_on_lower.logf[i], upper_on_lower.logf[i + 1]))
        q[i] = corner_average(lower.Q[i], lower.Q[i + 1],
                              upper_on_lower.Q[i], upper_on_lower.Q[i + 1])
        ru = corner_du(lower.r[i], lower.r[i + 1],
                       upper_on_lower.r[i], upper_on_lower.r[i + 1], du)
        rv = corner_dv(lower.r[i], lower.r[i + 1],
                       upper_on_lower.r[i], upper_on_lower.r[i + 1], dv)
        phi_re = corner_average(lower.phi_re[i], lower.phi_re[i + 1],
                                upper_on_lower.phi_re[i], upper_on_lower.phi_re[i + 1])
        phi_im = corner_average(lower.phi_im[i], lower.phi_im[i + 1],
                                upper_on_lower.phi_im[i], upper_on_lower.phi_im[i + 1])
        phiu_re = corner_du(lower.phi_re[i], lower.phi_re[i + 1],
                            upper_on_lower.phi_re[i], upper_on_lower.phi_re[i + 1], du)
        phiu_im = corner_du(lower.phi_im[i], lower.phi_im[i + 1],
                            upper_on_lower.phi_im[i], upper_on_lower.phi_im[i + 1], du)
        Au = corner_average(lower.Au[i], lower.Au[i + 1],
                            upper_on_lower.Au[i], upper_on_lower.Au[i + 1])
        Av = corner_average(lower.Av[i], lower.Av[i + 1],
                            upper_on_lower.Av[i], upper_on_lower.Av[i + 1])
        source = if reduced_scalar
            stress_energy_reduced_scalar(r, f, q[i], ru, rv, phi_re, phi_im,
                                         phiu_re, zero(phiu_re), phiu_im,
                                         zero(phiu_im), Au, Av, scalar_charge)
        else
            stress_energy(r, f, q[i], phi_re, phi_im, phiu_re, zero(phiu_re),
                          phiu_im, zero(phiu_im), Au, Av, scalar_charge)
        end
        expected_q_u[i] = r^2 * source.Ju / 8
        expected_mass_u[i] = -r^2 * rv * source.Tuu / (4f) +
                             q[i] * expected_q_u[i] / r
    end
    return u, (lower.v + upper.v) / 2, q, expected_q_u, expected_mass_u
end

"""
Evaluate the nonlinear charged Maxwell constraint along a fixed-`V` slice:

    Q_U = r^2 J_U / 8.

For the GP2026 production path set `reduced_scalar=true`, because that state
stores `Psi=r*Phi`. In either representation this corresponds to the
Gelles-Pretorius law `Q_U=4*pi*r^2*J_GP,U`.
"""
function charged_charge_flux_u_profile(lower::NLSlice, upper::NLSlice, scalar_charge;
                                       reduced_scalar::Bool=false)
    u, v, q, expected_cells, _ =
        charged_u_flux_cells(lower, upper, scalar_charge; reduced_scalar)
    length(q) >= 2 ||
        throw(ArgumentError("charge-flux check needs at least three U grid points"))
    centered_u = [(u[i] + u[i + 1]) / 2 for i in firstindex(u):lastindex(u)-1]
    q_u = diff(q) ./ diff(u)
    expected_q_u = [(expected_cells[i] + expected_cells[i + 1]) / 2
                    for i in firstindex(expected_cells):lastindex(expected_cells)-1]
    return centered_u, v, q_u, expected_q_u, q_u .- expected_q_u
end

charged_charge_flux_u_profile(lower::NLSlice, upper::NLSlice, ep::EvolutionParams;
                              reduced_scalar::Bool=false) =
    charged_charge_flux_u_profile(lower, upper, ep.scalar_charge; reduced_scalar)

function charged_charge_flux_u_profile(st::AdaptiveNLState, scalar_charge; target_v=150.0,
                                       reduced_scalar::Bool=false)
    length(st.slices) >= 2 ||
        throw(ArgumentError("charge-flux check needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return charged_charge_flux_u_profile(st.slices[j], st.slices[j + 1], scalar_charge;
                                         reduced_scalar)
end

charged_charge_flux_u_profile(st::AdaptiveNLState, ep::EvolutionParams; target_v=150.0,
                              reduced_scalar::Bool=false) =
    charged_charge_flux_u_profile(st, ep.scalar_charge; target_v, reduced_scalar)

"""
Evaluate the complementary nonlinear Maxwell constraint at fixed `U`:

    Q_V = -r^2 J_V / 8.

The GP2026 production march prescribes `Q` on `U=U0` and advances it with
the `U` constraint. This profile therefore tests the unevolved Maxwell
constraint independently.
"""
function charged_charge_flux_v_profile(st::AdaptiveNLState, scalar_charge; target_u,
                                       reduced_scalar::Bool=false)
    length(st.slices) >= 2 ||
        throw(ArgumentError("charge-flux check needs at least two V slices"))

    ncell = length(st.slices) - 1
    v = similar([slice.v for slice in st.slices], ncell)
    q_v = similar(v)
    expected_q_v = similar(v)
    for j in 1:ncell
        lower, upper = st.slices[j], st.slices[j + 1]
        dv = upper.v - lower.v
        v[j] = (lower.v + upper.v) / 2
        lr = interpolate_linear(lower.u, lower.r, target_u)
        ur = interpolate_linear(upper.u, upper.r, target_u)
        lf = interpolate_linear(lower.u, lower.logf, target_u)
        uf = interpolate_linear(upper.u, upper.logf, target_u)
        lq = interpolate_linear(lower.u, lower.Q, target_u)
        uq = interpolate_linear(upper.u, upper.Q, target_u)
        lphi_re = interpolate_linear(lower.u, lower.phi_re, target_u)
        uphi_re = interpolate_linear(upper.u, upper.phi_re, target_u)
        lphi_im = interpolate_linear(lower.u, lower.phi_im, target_u)
        uphi_im = interpolate_linear(upper.u, upper.phi_im, target_u)
        r = (lr + ur) / 2
        rv = (ur - lr) / dv
        f = exp((lf + uf) / 2)
        q = (lq + uq) / 2
        phi_re = (lphi_re + uphi_re) / 2
        phi_im = (lphi_im + uphi_im) / 2
        phi_re_v = (uphi_re - lphi_re) / dv
        phi_im_v = (uphi_im - lphi_im) / dv
        Au = (interpolate_linear(lower.u, lower.Au, target_u) +
              interpolate_linear(upper.u, upper.Au, target_u)) / 2
        Av = (interpolate_linear(lower.u, lower.Av, target_u) +
              interpolate_linear(upper.u, upper.Av, target_u)) / 2
        source = if reduced_scalar
            stress_energy_reduced_scalar(r, f, q, zero(r), rv, phi_re, phi_im,
                                         zero(r), phi_re_v, zero(r), phi_im_v,
                                         Au, Av, scalar_charge)
        else
            stress_energy(r, f, q, phi_re, phi_im, zero(r), phi_re_v,
                          zero(r), phi_im_v, Au, Av, scalar_charge)
        end
        q_v[j] = (uq - lq) / dv
        expected_q_v[j] = -r^2 * source.Jv / 8
    end
    return v, target_u, q_v, expected_q_v, q_v .- expected_q_v
end

charged_charge_flux_v_profile(st::AdaptiveNLState, ep::EvolutionParams; target_u,
                              reduced_scalar::Bool=false) =
    charged_charge_flux_v_profile(st, ep.scalar_charge; target_u, reduced_scalar)

"""
Integrate the charged Maxwell constraint in `U`, anchored to the outer
geometrically sampled charge.
"""
function charged_flux_integrated_charge_profile(lower::NLSlice, upper::NLSlice, scalar_charge;
                                                reduced_scalar::Bool=false)
    u, v, geometric_q, _, _ =
        charged_u_flux_cells(lower, upper, scalar_charge; reduced_scalar)
    _, _, _, expected_q_u, _ =
        charged_charge_flux_u_profile(lower, upper, scalar_charge; reduced_scalar)
    flux_q = similar(geometric_q)
    flux_q[1] = geometric_q[1]
    for i in 2:length(flux_q)
        flux_q[i] = flux_q[i - 1] + (u[i] - u[i - 1]) * expected_q_u[i - 1]
    end
    return u, v, geometric_q, flux_q, geometric_q .- flux_q
end

charged_flux_integrated_charge_profile(lower::NLSlice, upper::NLSlice, ep::EvolutionParams;
                                       reduced_scalar::Bool=false) =
    charged_flux_integrated_charge_profile(lower, upper, ep.scalar_charge; reduced_scalar)

function charged_flux_integrated_charge_profile(st::AdaptiveNLState, scalar_charge;
                                                target_v=150.0,
                                                reduced_scalar::Bool=false)
    length(st.slices) >= 2 ||
        throw(ArgumentError("flux-integrated charge needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return charged_flux_integrated_charge_profile(st.slices[j], st.slices[j + 1],
                                                  scalar_charge; reduced_scalar)
end

charged_flux_integrated_charge_profile(st::AdaptiveNLState, ep::EvolutionParams;
                                       target_v=150.0, reduced_scalar::Bool=false) =
    charged_flux_integrated_charge_profile(st, ep.scalar_charge; target_v, reduced_scalar)

"""
Evaluate the charged renormalized Hawking-mass balance law:

    varpi_U = -r^2 r_V T_UU/(4f) + Q Q_U/r.

For the source scalar `T_UU=2|D_U Phi|^2`; on reduced-field runs `Phi` is
reconstructed from `Psi=r*Phi`. The first term reduces to MRT Eq. (14) when
the scalar is real and uncharged.
"""
function charged_mass_flux_u_profile(lower::NLSlice, upper::NLSlice, scalar_charge;
                                     rn_background::Union{Nothing,RNParams}=nothing,
                                     reduced_scalar::Bool=false)
    u, v, mass = renormalized_hawking_mass_profile(lower, upper; rn_background)
    _, _, _, _, expected_cells =
        charged_u_flux_cells(lower, upper, scalar_charge; reduced_scalar)
    length(mass) >= 2 ||
        throw(ArgumentError("mass-flux check needs at least three U grid points"))
    centered_u = [(u[i] + u[i + 1]) / 2 for i in firstindex(u):lastindex(u)-1]
    mass_u = diff(mass) ./ diff(u)
    expected_mass_u = [(expected_cells[i] + expected_cells[i + 1]) / 2
                       for i in firstindex(expected_cells):lastindex(expected_cells)-1]
    return centered_u, v, mass_u, expected_mass_u, mass_u .- expected_mass_u
end

charged_mass_flux_u_profile(lower::NLSlice, upper::NLSlice, ep::EvolutionParams;
                            rn_background::Union{Nothing,RNParams}=nothing,
                            reduced_scalar::Bool=false) =
    charged_mass_flux_u_profile(lower, upper, ep.scalar_charge; rn_background, reduced_scalar)

function charged_mass_flux_u_profile(st::AdaptiveNLState, scalar_charge; target_v=150.0,
                                     rn_background::Union{Nothing,RNParams}=nothing,
                                     reduced_scalar::Bool=false)
    length(st.slices) >= 2 ||
        throw(ArgumentError("mass-flux check needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return charged_mass_flux_u_profile(st.slices[j], st.slices[j + 1], scalar_charge;
                                       rn_background, reduced_scalar)
end

charged_mass_flux_u_profile(st::AdaptiveNLState, ep::EvolutionParams; target_v=150.0,
                            rn_background::Union{Nothing,RNParams}=nothing,
                            reduced_scalar::Bool=false) =
    charged_mass_flux_u_profile(st, ep.scalar_charge; target_v, rn_background, reduced_scalar)

"""
Integrate the charged mass-balance law in `U`, anchored to the outer
geometrically reconstructed mass.
"""
function charged_flux_integrated_mass_profile(lower::NLSlice, upper::NLSlice, scalar_charge;
                                              rn_background::Union{Nothing,RNParams}=nothing,
                                              reduced_scalar::Bool=false)
    u, v, geometric_mass = renormalized_hawking_mass_profile(lower, upper; rn_background)
    _, _, _, expected_mass_u, _ =
        charged_mass_flux_u_profile(lower, upper, scalar_charge; rn_background, reduced_scalar)
    flux_mass = similar(geometric_mass)
    flux_mass[1] = geometric_mass[1]
    for i in 2:length(flux_mass)
        flux_mass[i] = flux_mass[i - 1] + (u[i] - u[i - 1]) * expected_mass_u[i - 1]
    end
    return u, v, geometric_mass, flux_mass, geometric_mass .- flux_mass
end

charged_flux_integrated_mass_profile(lower::NLSlice, upper::NLSlice, ep::EvolutionParams;
                                     rn_background::Union{Nothing,RNParams}=nothing,
                                     reduced_scalar::Bool=false) =
    charged_flux_integrated_mass_profile(lower, upper, ep.scalar_charge;
                                         rn_background, reduced_scalar)

function charged_flux_integrated_mass_profile(st::AdaptiveNLState, scalar_charge; target_v=150.0,
                                              rn_background::Union{Nothing,RNParams}=nothing,
                                              reduced_scalar::Bool=false)
    length(st.slices) >= 2 ||
        throw(ArgumentError("flux-integrated mass needs at least two V slices"))
    vmid = [(st.slices[j].v + st.slices[j + 1].v) / 2
            for j in 1:length(st.slices)-1]
    _, j = findmin(abs.(vmid .- target_v))
    return charged_flux_integrated_mass_profile(st.slices[j], st.slices[j + 1],
                                                scalar_charge; rn_background, reduced_scalar)
end

charged_flux_integrated_mass_profile(st::AdaptiveNLState, ep::EvolutionParams; target_v=150.0,
                                     rn_background::Union{Nothing,RNParams}=nothing,
                                     reduced_scalar::Bool=false) =
    charged_flux_integrated_mass_profile(st, ep.scalar_charge;
                                         target_v, rn_background, reduced_scalar)

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

"""
Return the GP2026 reduced scalar amplitude `|r*phi_GP|` on one stored slice.

This assumes the slice was evolved with `reduced_scalar=true`, so its scalar
arrays contain `Psi=sqrt(32*pi)*r*phi_GP`.
"""
function gp2026_rphi_profile(slice::NLSlice)
    return slice.u, hypot.(slice.phi_re, slice.phi_im) ./ sqrt(32 * pi)
end

"""
Sample `|r*phi_GP|` at the apparent horizon on all GP2026 slices where a
zero of `r_V` is bracketed.
"""
function gp2026_horizon_rphi_series(st::AdaptiveNLState)
    T = eltype(first(st.slices).r)
    horizon_v = T[]
    horizon_u = T[]
    amplitude = T[]
    for j in 2:length(st.slices)
        lower, upper = st.slices[j - 1], st.slices[j]
        u = [(upper.u[i] + upper.u[i + 1]) / 2 for i in 1:length(upper.u)-1]
        rv = adaptive_outgoing_expansion(lower, upper)
        uh = apparent_horizon_location(u, rv)
        isnothing(uh) && continue
        psi_re = interpolate_linear(upper.u, upper.phi_re, uh)
        psi_im = interpolate_linear(upper.u, upper.phi_im, uh)
        push!(horizon_v, upper.v)
        push!(horizon_u, uh)
        push!(amplitude, hypot(psi_re, psi_im) / sqrt(32 * pi))
    end
    return horizon_v, horizon_u, amplitude
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
