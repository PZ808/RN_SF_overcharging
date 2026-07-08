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

struct TrappedSurfaceSample{T<:Real}
    row_index::Int
    u::T
    v::T
    r::T
    q::T
    rv::T
end

struct VTrapDiagnostic{T<:Real}
    status::Symbol
    trapped::Bool
    trap::Union{Nothing,TrappedSurfaceSample{T}}
    closest::TrappedSurfaceSample{T}
end

struct HorizonChargeDensitySample{T<:Real}
    row_index::Int
    u::T
    v::T
    r::T
    q::T
    rv::T
    q_v::T
    surface_density::T
    flux_density_v::T
    q_u::T
    r_u::T
    radial_density_proxy::T
end

"""
Gauge-invariant quasi-local quantities at a marginally trapped sphere.

At `r_V=0`, the renormalized Hawking identity reduces to
`M=(r^2+Q^2)/(2r)` without requiring coordinate derivatives. The
surface-gravity value is the Reissner-Nordstrom proxy constructed from this
instantaneous `(M,Q)` pair; it is a diagnostic, not an assumption that the
surface is stationary.
"""
function trapped_surface_invariants(sample::TrappedSurfaceSample)
    r = sample.r
    q = sample.q
    r > zero(r) || throw(ArgumentError("trapped-surface radius must be positive"))
    mass = (r^2 + q^2) / (2r)
    discriminant = max(mass^2 - q^2, zero(mass))
    root = sqrt(discriminant)
    outer_radius = mass + root
    kappa_proxy = root / outer_radius^2
    return (
        mass=mass,
        q_over_m=q / mass,
        r_over_m=r / mass,
        one_minus_q_over_m=one(mass) - q / mass,
        one_minus_r_over_m=one(mass) - r / mass,
        rn_surface_gravity_proxy=kappa_proxy,
    )
end

"""
Centered outgoing expansion proxy `r_V` on a single GP2026 `U` row.

The apparent-horizon condition used by Gelles-Pretorius is a sign flip of the
outgoing expansion, `theta_+ proportional to r_V`. This helper keeps the
row-oriented diagnostic consistent with the fixed-`V` slice detector.
"""
row_outgoing_expansion(row::NLRow) = coordinate_derivative(row.r, row.v)

function row_expansion_minimum(row::NLRow; row_index::Int=0)
    rv = row_outgoing_expansion(row)
    value, index = findmin(rv)
    return TrappedSurfaceSample(row_index, row.u, row.v[index], row.r[index],
                                row.Q[index], value)
end

function row_apparent_horizon_crossing(row::NLRow; row_index::Int=0)
    rv = row_outgoing_expansion(row)
    crossing = findfirst(value -> value <= zero(value), rv)
    isnothing(crossing) && return nothing
    if crossing == firstindex(rv)
        return TrappedSurfaceSample(row_index, row.u, row.v[crossing],
                                    row.r[crossing], row.Q[crossing],
                                    rv[crossing])
    end

    i = crossing
    fraction = rv[i - 1] / (rv[i - 1] - rv[i])
    return TrappedSurfaceSample(
        row_index,
        row.u,
        segment_value(row.v, i - 1, fraction),
        segment_value(row.r, i - 1, fraction),
        segment_value(row.Q, i - 1, fraction),
        zero(rv[i]),
    )
end

function row_horizon_charge_density_sample(row::NLRow; row_index::Int=0)
    rv = row_outgoing_expansion(row)
    crossing = findfirst(value -> value <= zero(value), rv)
    isnothing(crossing) && return nothing
    qv = coordinate_derivative(row.Q, row.v)

    if crossing == firstindex(rv)
        r = row.r[crossing]
        q = row.Q[crossing]
        return HorizonChargeDensitySample(
            row_index,
            row.u,
            row.v[crossing],
            r,
            q,
            rv[crossing],
            qv[crossing],
            q / (4pi * r^2),
            qv[crossing] / (4pi * r^2),
            oftype(r, NaN),
            oftype(r, NaN),
            oftype(r, NaN),
        )
    end

    i = crossing
    fraction = rv[i - 1] / (rv[i - 1] - rv[i])
    r = segment_value(row.r, i - 1, fraction)
    q = segment_value(row.Q, i - 1, fraction)
    q_v = segment_value(qv, i - 1, fraction)
    return HorizonChargeDensitySample(
        row_index,
        row.u,
        segment_value(row.v, i - 1, fraction),
        r,
        q,
        zero(rv[i]),
        q_v,
        q / (4pi * r^2),
        q_v / (4pi * r^2),
        oftype(r, NaN),
        oftype(r, NaN),
        oftype(r, NaN),
    )
end

function with_radial_charge_density_proxy(
    sample::HorizonChargeDensitySample,
    rows::AbstractVector{<:NLRow},
)
    i = sample.row_index
    if i <= firstindex(rows) && i < lastindex(rows)
        lower, upper = rows[i], rows[i + 1]
    elseif i >= lastindex(rows) && i > firstindex(rows)
        lower, upper = rows[i - 1], rows[i]
    elseif firstindex(rows) < i < lastindex(rows)
        lower, upper = rows[i - 1], rows[i + 1]
    else
        return sample
    end
    first(lower.v) <= sample.v <= last(lower.v) &&
        first(upper.v) <= sample.v <= last(upper.v) || return sample
    lower_point = row_point_at_v(lower, sample.v)
    upper_point = row_point_at_v(upper, sample.v)
    du = upper.u - lower.u
    du > 0 || return sample
    q_u = (upper_point.Q - lower_point.Q) / du
    r_u = (upper_point.r - lower_point.r) / du
    radial_q = r_u != 0 ? q_u / r_u : oftype(sample.r, NaN)
    proxy = isfinite(radial_q) ?
            (radial_q + sample.q_v / 2) / (4pi * sample.r^2) :
            oftype(sample.r, NaN)
    return HorizonChargeDensitySample(
        sample.row_index,
        sample.u,
        sample.v,
        sample.r,
        sample.q,
        sample.rv,
        sample.q_v,
        sample.surface_density,
        sample.flux_density_v,
        q_u,
        r_u,
        proxy,
    )
end

function apparent_horizon_charge_density_series(rows::AbstractVector{<:NLRow})
    samples = HorizonChargeDensitySample[]
    for (row_index, row) in pairs(rows)
        sample = row_horizon_charge_density_sample(row; row_index)
        isnothing(sample) ||
            push!(samples, with_radial_charge_density_proxy(sample, rows))
    end
    return samples
end

function quadratic_sample_value(left, center, right, x)
    x1 = left.u - center.u
    x3 = right.u - center.u
    xx = x - center.u
    weight1 = xx * (xx - x3) / (x1 * (x1 - x3))
    weight2 = (xx - x1) * (xx - x3) / (x1 * x3)
    weight3 = xx * (xx - x1) / (x3 * (x3 - x1))
    return weight1, weight2, weight3
end

function quadratic_trapped_surface_sample(
    left::TrappedSurfaceSample,
    center::TrappedSurfaceSample,
    right::TrappedSurfaceSample,
)
    left.u < center.u < right.u || return center
    x1 = left.u - center.u
    x3 = right.u - center.u
    y1 = left.v - center.v
    y3 = right.v - center.v
    determinant = x1^2 * x3 - x3^2 * x1
    determinant != zero(determinant) || return center
    curvature = (y1 * x3 - y3 * x1) / determinant
    slope = (x1^2 * y3 - x3^2 * y1) / determinant
    curvature > zero(curvature) || return center
    offset = -slope / (2curvature)
    x1 < offset < x3 || return center
    u = center.u + offset
    weights = quadratic_sample_value(left, center, right, u)
    blend(field) =
        weights[1] * getproperty(left, field) +
        weights[2] * getproperty(center, field) +
        weights[3] * getproperty(right, field)
    return TrappedSurfaceSample(
        center.row_index,
        u,
        center.v + slope * offset + curvature * offset^2,
        blend(:r),
        blend(:q),
        zero(center.rv),
    )
end

"""
Refine the discrete minimum of the apparent-horizon crossing curve in `U`.

The direct GP diagnostic minimizes crossing `V` over completed rows. When
crossings on both neighboring rows are available, this helper fits the local
nonuniform three-point parabola and returns its interior vertex. It falls back
to the discrete sample unless the fit is convex and bracketed.
"""
function refined_vtrap_sample(rows::AbstractVector{<:NLRow})
    crossings = TrappedSurfaceSample[]
    for (row_index, row) in pairs(rows)
        crossing = row_apparent_horizon_crossing(row; row_index)
        isnothing(crossing) || push!(crossings, crossing)
    end
    isempty(crossings) && return nothing
    index = argmin(sample.v for sample in crossings)
    1 < index < length(crossings) || return crossings[index]
    return quadratic_trapped_surface_sample(
        crossings[index - 1],
        crossings[index],
        crossings[index + 1],
    )
end

"""
Return the direct `Vtrap` diagnostic for GP2026 row data.

If a row contains `r_V <= 0`, `trap` is the earliest apparent-horizon point in
`V` among all rows. If no such point exists, `trap === nothing`, `status` is the
caller-provided missing status, and `closest` records the row point with the
smallest positive outgoing expansion. This makes a missing `Vtrap` distinguish
"not trapped yet" from "the run stalled before trapping".
"""
function vtrap_diagnostic(rows::AbstractVector{<:NLRow}; missing_status::Symbol=:untrapped)
    isempty(rows) && throw(ArgumentError("Vtrap diagnostic needs at least one row"))

    closest = row_expansion_minimum(first(rows); row_index=firstindex(rows))
    trap = nothing
    for (row_index, row) in pairs(rows)
        minimum_sample = row_expansion_minimum(row; row_index)
        if minimum_sample.rv < closest.rv
            closest = minimum_sample
        end

        crossing = row_apparent_horizon_crossing(row; row_index)
        if !isnothing(crossing) &&
           (isnothing(trap) || crossing.v < trap.v)
            trap = crossing
        end
    end

    status = isnothing(trap) ? missing_status : :trapped
    return VTrapDiagnostic(status, !isnothing(trap), trap, closest)
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

struct ThroatBoundaryObservables{T<:Real}
    psi_abs::T
    rphi_gp_re::T
    rphi_gp_im::T
    rphi_gp_abs::T
    raw_phase::T
    rphi_gp_abs_v::T
    covariant_phase_v::T
    covariant_dv_rphi_abs::T
    Jv::T
    Tvv::T
    q_v_source::T
    q_v_residual::T
    outgoing_constraint_source::T
    q_over_r::T
    one_minus_absq_over_r::T
    logf_rho::T
end

"""
Gauge-invariant and near-throat observables at a fixed-`rho` boundary sample.

The stored GP2026 scalar is `Psi=sqrt(32*pi)*r*phi_GP`. The returned
`rphi_gp_*` fields divide out this normalization. `raw_phase` is included for
mode-tracking convenience, while `covariant_phase_v` is the gauge-invariant
combination `theta_V - e A_V` computed from `Im(Psi^* D_V Psi)/|Psi|^2`.
"""
function throat_boundary_observables(sample::ThroatBoundarySample,
                                     ep::EvolutionParams)
    scale = sqrt(32 * pi)
    f = exp(sample.logf)
    source = stress_energy_reduced_scalar(
        sample.r, f, sample.q, zero(sample.r), sample.r_v,
        sample.phi_re, sample.phi_im,
        zero(sample.phi_re_v), sample.phi_re_v,
        zero(sample.phi_im_v), sample.phi_im_v,
        sample.Au, sample.Av, ep.scalar_charge,
    )

    psi2 = sample.phi_re^2 + sample.phi_im^2
    psi_abs = sqrt(psi2)
    rphi_gp_re = sample.phi_re / scale
    rphi_gp_im = sample.phi_im / scale
    rphi_gp_abs = psi_abs / scale
    raw_phase = atan(sample.phi_im, sample.phi_re)

    Dv_re = sample.phi_re_v + ep.scalar_charge * sample.Av * sample.phi_im
    Dv_im = sample.phi_im_v - ep.scalar_charge * sample.Av * sample.phi_re
    rphi_gp_abs_v = psi_abs > 0 ?
                    (sample.phi_re * sample.phi_re_v +
                     sample.phi_im * sample.phi_im_v) / (psi_abs * scale) :
                    zero(psi_abs)
    covariant_phase_v = psi2 > 0 ?
                        (sample.phi_re * Dv_im - sample.phi_im * Dv_re) / psi2 :
                        zero(psi2)
    covariant_dv_rphi_abs = hypot(Dv_re, Dv_im) / scale
    q_v_source = -sample.r^2 * source.Jv / 8
    logf_rho = isfinite(sample.rho_v) && sample.rho_v != 0 ?
               sample.logf - log(abs(sample.rho_v)) :
               typeof(sample.logf)(Inf)

    return ThroatBoundaryObservables(
        psi_abs,
        rphi_gp_re,
        rphi_gp_im,
        rphi_gp_abs,
        raw_phase,
        rphi_gp_abs_v,
        covariant_phase_v,
        covariant_dv_rphi_abs,
        source.Jv,
        source.Tvv,
        q_v_source,
        sample.q_v - q_v_source,
        outgoing_constraint_source(sample.r, f, source),
        sample.q / sample.r,
        one(sample.r) - abs(sample.q) / sample.r,
        logf_rho,
    )
end

function checked_grid_size(st::NLState, g::Grid)
    size(st.r) == size(st.logf) == size(st.phi_re) == size(st.phi_im) ==
        size(st.Au) == size(st.Av) == size(st.Q) == size(g) ||
        throw(ArgumentError("state field sizes must match the grid"))
    return size(g)
end

finite_abs(value) = isfinite(value) ? abs(value) : typeof(value)(Inf)

"""
Summarize centered-cell residuals for the fully nonlinear EMKG equations.

The residuals are evaluated on an already evolved rectangular `NLState`, using
the same corner averages and centered derivatives as `step_nonlinear_cell!`.
Set `reduced_scalar=true` for the GP2026 production branch, where
`phi_re/phi_im` store `Psi=sqrt(32*pi)*r*phi_GP`. The two `logf` variants are
negative controls for the most plausible paper-copying mistakes:

* `max_abs_logf_gp_literal` uses the literal printed GP Eq. (4) coefficient
  after converting to `f_code=2f_GP`.
* `max_abs_logf_coulomb2` uses an extra factor of two in the Coulomb term.

The `Q_UV` residual is meaningful for the GP2026 production mode
`hyperbolic_charge=true`; in constraint-marched runs it is just measured
against zero. For `cell_solver=:newton_direct`, `max_abs_cell_residual`
measures the unscaled seven-equation algebraic residual and `max_abs_f_uv`
measures the direct lapse equation. This avoids interpreting roundoff divided
by `Delta U*Delta V` as a truncation error.
"""
function cell_equation_residual_summary(st::NLState, g::Grid,
                                        ep::EvolutionParams;
                                        reduced_scalar::Bool=false,
                                        hyperbolic_charge::Bool=false,
                                        subtract_rn_background::Bool=false,
                                        cell_solver::Symbol=:picard_log,
                                        interior_only::Bool=false)
    hyperbolic_charge && !reduced_scalar &&
        throw(ArgumentError("hyperbolic charge residual requires reduced_scalar=true"))
    cell_solver in (:picard_log, :newton_direct) ||
        throw(ArgumentError("cell_solver must be :picard_log or :newton_direct"))
    cell_solver === :newton_direct && subtract_rn_background &&
        throw(ArgumentError(
            "RN background subtraction is not implemented for :newton_direct",
        ))
    nu, nv = checked_grid_size(st, g)
    nu >= 2 && nv >= 2 ||
        throw(ArgumentError("cell residuals need at least one cell"))

    i_range = if interior_only && nu > 3
        2:(nu - 2)
    else
        1:(nu - 1)
    end
    j_range = if interior_only && nv > 3
        2:(nv - 2)
    else
        1:(nv - 1)
    end
    T = promote_type(eltype(g.u), eltype(g.v), eltype(st.r),
                     typeof(ep.scalar_charge))
    maxima = Dict{Symbol,T}(
        :r_uv => zero(T),
        :f_uv => zero(T),
        :logf_uv => zero(T),
        :cell_residual => zero(T),
        :logf_gp_literal => zero(T),
        :logf_coulomb2 => zero(T),
        :psi_re_uv => zero(T),
        :psi_im_uv => zero(T),
        :q_uv => zero(T),
        :q_u_constraint => zero(T),
        :q_v_constraint => zero(T),
        :au_v => zero(T),
        :av_u => zero(T),
        :quasilorenz => zero(T),
        :faraday => zero(T),
    )

    function update!(key, value)
        maxima[key] = max(maxima[key], convert(T, finite_abs(value)))
    end

    cells = 0
    e = ep.scalar_charge
    for i in i_range, j in j_range
        du = g.u[i + 1] - g.u[i]
        dv = g.v[j + 1] - g.v[j]
        r00, r10, r01, r11 = st.r[i, j], st.r[i + 1, j],
                              st.r[i, j + 1], st.r[i + 1, j + 1]
        lf00, lf10, lf01, lf11 = st.logf[i, j], st.logf[i + 1, j],
                                  st.logf[i, j + 1], st.logf[i + 1, j + 1]
        pr00, pr10, pr01, pr11 = st.phi_re[i, j], st.phi_re[i + 1, j],
                                  st.phi_re[i, j + 1], st.phi_re[i + 1, j + 1]
        pi00, pi10, pi01, pi11 = st.phi_im[i, j], st.phi_im[i + 1, j],
                                  st.phi_im[i, j + 1], st.phi_im[i + 1, j + 1]
        au00, au10, au01, au11 = st.Au[i, j], st.Au[i + 1, j],
                                  st.Au[i, j + 1], st.Au[i + 1, j + 1]
        av00, av10, av01, av11 = st.Av[i, j], st.Av[i + 1, j],
                                  st.Av[i, j + 1], st.Av[i + 1, j + 1]
        q00, q10, q01, q11 = st.Q[i, j], st.Q[i + 1, j],
                              st.Q[i, j + 1], st.Q[i + 1, j + 1]

        r = corner_average(r00, r10, r01, r11)
        f00, f10, f01, f11 = exp(lf00), exp(lf10), exp(lf01), exp(lf11)
        lf = corner_average(lf00, lf10, lf01, lf11)
        f = cell_solver === :newton_direct ?
            corner_average(f00, f10, f01, f11) : exp(lf)
        pr = corner_average(pr00, pr10, pr01, pr11)
        pii = corner_average(pi00, pi10, pi01, pi11)
        au = corner_average(au00, au10, au01, au11)
        av = corner_average(av00, av10, av01, av11)
        q = corner_average(q00, q10, q01, q11)

        ru = corner_du(r00, r10, r01, r11, du)
        rv = corner_dv(r00, r10, r01, r11, dv)
        fu = corner_du(f00, f10, f01, f11, du)
        fv = corner_dv(f00, f10, f01, f11, dv)
        pru = corner_du(pr00, pr10, pr01, pr11, du)
        prv = corner_dv(pr00, pr10, pr01, pr11, dv)
        piu = corner_du(pi00, pi10, pi01, pi11, du)
        piv = corner_dv(pi00, pi10, pi01, pi11, dv)
        au_v = corner_dv(au00, au10, au01, au11, dv)
        av_u = corner_du(av00, av10, av01, av11, du)
        q_u = corner_du(q00, q10, q01, q11, du)
        q_v = corner_dv(q00, q10, q01, q11, dv)

        source = if reduced_scalar
            stress_energy_reduced_scalar(r, f, q, ru, rv, pr, pii, pru, prv,
                                         piu, piv, au, av, e)
        else
            stress_energy(r, f, q, pr, pii, pru, prv, piu, piv, au, av, e)
        end
        ruv, lfuv = metric_rhs(r, f, ru, rv, q, source)
        _, fuv = direct_lapse_rhs(r, f, fu, fv, ru, rv, q, source)
        pruv, piuv = if reduced_scalar
            charged_reduced_scalar_rhs(r, ruv, pru, prv, piu, piv,
                                       pr, pii, au, av, e)
        else
            charged_scalar_rhs(r, ru, rv, pru, prv, piu, piv,
                               pr, pii, au, av, e)
        end
        auv, avu, _, q_u_source, q_v_source = maxwell_rhs(r, f, q, source)
        quv = reduced_scalar && hyperbolic_charge ?
              charged_reduced_charge_rhs(r, f, q, pr, pii, pru, prv,
                                         piu, piv, au, av, e) :
              zero(q)

        r_defect, lf_defect = subtract_rn_background ?
                               rn_background_update_defect(g, ep, i, j, du, dv) :
                               (zero(du), zero(du))
        ruv_observed = (r11 - r10 - r01 + r00) / (du * dv)
        fuv_observed = (f11 - f10 - f01 + f00) / (du * dv)
        lfuv_observed = (lf11 - lf10 - lf01 + lf00) / (du * dv)
        pruv_observed = (pr11 - pr10 - pr01 + pr00) / (du * dv)
        piuv_observed = (pi11 - pi10 - pi01 + pi00) / (du * dv)
        quv_observed = (q11 - q10 - q01 + q00) / (du * dv)

        lfuv_gp_literal = f / (4r^2) + 2 * ru * rv / r^2 -
                           q^2 * f / r^4 - source.scalar_logf_source
        lfuv_coulomb2 = f / (2r^2) + 2 * ru * rv / r^2 -
                        2 * q^2 * f / r^4 - source.scalar_logf_source

        update!(:r_uv, ruv_observed - ruv + r_defect / (du * dv))
        update!(:f_uv, fuv_observed - fuv)
        update!(:logf_uv, lfuv_observed - lfuv + lf_defect / (du * dv))
        update!(:logf_gp_literal,
                lfuv_observed - lfuv_gp_literal + lf_defect / (du * dv))
        update!(:logf_coulomb2,
                lfuv_observed - lfuv_coulomb2 + lf_defect / (du * dv))
        update!(:psi_re_uv, pruv_observed - pruv)
        update!(:psi_im_uv, piuv_observed - piuv)
        update!(:q_uv, quv_observed - quv)
        update!(:q_u_constraint, q_u - q_u_source)
        update!(:q_v_constraint, q_v - q_v_source)
        update!(:au_v, au_v - auv)
        update!(:av_u, av_u - avu)
        update!(:quasilorenz, au_v + av_u)
        update!(:faraday, av_u - au_v + q * f / (2r^2))
        update!(
            :cell_residual,
            (r11 - r10 - r01 + r00) - du * dv * ruv + r_defect,
        )
        lapse_cell_residual = if cell_solver === :newton_direct
            (f11 - f10 - f01 + f00) - du * dv * fuv
        else
            (lf11 - lf10 - lf01 + lf00) - du * dv * lfuv + lf_defect
        end
        update!(:cell_residual, lapse_cell_residual)
        update!(
            :cell_residual,
            (pr11 - pr10 - pr01 + pr00) - du * dv * pruv,
        )
        update!(
            :cell_residual,
            (pi11 - pi10 - pi01 + pi00) - du * dv * piuv,
        )
        update!(:cell_residual, 2dv * (au_v - auv))
        update!(:cell_residual, 2du * (av_u - avu))
        charge_cell_residual = hyperbolic_charge ?
                               (q11 - q10 - q01 + q00) - du * dv * quv :
                               (q11 - q01) - du * q_u_source
        update!(:cell_residual, charge_cell_residual)
        cells += 1
    end

    return (
        cells=cells,
        max_abs_cell_residual=maxima[:cell_residual],
        max_abs_r_uv=maxima[:r_uv],
        max_abs_f_uv=maxima[:f_uv],
        max_abs_logf_uv=maxima[:logf_uv],
        max_abs_logf_gp_literal=maxima[:logf_gp_literal],
        max_abs_logf_coulomb2=maxima[:logf_coulomb2],
        max_abs_psi_re_uv=maxima[:psi_re_uv],
        max_abs_psi_im_uv=maxima[:psi_im_uv],
        max_abs_q_uv=maxima[:q_uv],
        max_abs_q_u_constraint=maxima[:q_u_constraint],
        max_abs_q_v_constraint=maxima[:q_v_constraint],
        max_abs_au_v=maxima[:au_v],
        max_abs_av_u=maxima[:av_u],
        max_abs_quasilorenz=maxima[:quasilorenz],
        max_abs_faraday=maxima[:faraday],
    )
end

"""
Audit the GP2026 characteristic initial data against the discrete constraints.

Residuals use the same midpoint and trapezoidal rules as
`initialize_gp2026_single_pulse!`. The stored metric coefficient is
`f_code=2f_GP` and the stored scalar is
`Psi=sqrt(32*pi)*r*phi_GP`.
"""
function gp2026_initial_constraint_residuals(
    st::NLState,
    g::Grid,
    ep::EvolutionParams;
    U0=first(g.u),
    V0=first(g.v),
    width=20.0,
    M0=ep.rn.M,
    pulse_leg_gauge::Symbol=:areal_affine,
)
    first(g.u) == U0 ||
        throw(ArgumentError("GP initial audit requires first U point U0"))
    first(g.v) == V0 ||
        throw(ArgumentError("GP initial audit requires first V point V0"))
    length(g.u) >= 2 && length(g.v) >= 2 ||
        throw(ArgumentError("GP initial audit requires at least two points per leg"))

    i0 = firstindex(g.u)
    j0 = firstindex(g.v)
    fcorner = gp2026_fcorner_code(
        ep;
        U0,
        V0,
        M0,
        pulse_leg_gauge,
    )
    ru0 = gp2026_extremal_gauge_ru(U0; M0)
    rv0 = gp2026_extremal_gauge_rv(
        U0, V0;
        U0,
        V0,
        M0,
        pulse_leg_gauge,
    )
    corner_mass = renormalized_hawking_mass(
        st.r[i0, j0],
        exp(st.logf[i0, j0]),
        ru0,
        rv0,
        st.Q[i0, j0],
    )

    T = promote_type(eltype(st.r), eltype(g.u), eltype(g.v))
    na_radius = zero(T)
    na_lapse = zero(T)
    na_charge = zero(T)
    na_scalar = zero(T)
    na_au = zero(T)
    na_av = zero(T)
    reference_ru_over_f = ru0 / exp(st.logf[i0, j0])
    for i in eachindex(g.u)
        expected_r = gp2026_extremal_gauge_initial_radius(
            g.u[i], V0;
            U0,
            V0,
            M0,
            pulse_leg_gauge,
        )
        na_radius = max(na_radius, abs(st.r[i, j0] - expected_r))
        ru_over_f = gp2026_extremal_gauge_ru(g.u[i]; M0) /
                    exp(st.logf[i, j0])
        na_lapse = max(na_lapse, abs(ru_over_f - reference_ru_over_f))
        na_charge = max(na_charge, abs(st.Q[i, j0] - ep.rn.Q0))
        na_scalar = max(
            na_scalar,
            hypot(st.phi_re[i, j0], st.phi_im[i, j0]),
        )
        na_au = max(na_au, abs(st.Au[i, j0]))
    end
    for i in i0+1:lastindex(g.u)
        du = g.u[i] - g.u[i - 1]
        left = -st.Q[i - 1, j0] * exp(st.logf[i - 1, j0]) /
               (4 * st.r[i - 1, j0]^2)
        right = -st.Q[i, j0] * exp(st.logf[i, j0]) /
                (4 * st.r[i, j0]^2)
        residual = st.Av[i, j0] - st.Av[i - 1, j0] -
                   du * (left + right) / 2
        na_av = max(na_av, abs(residual))
    end

    nb_radius = zero(T)
    nb_scalar = zero(T)
    nb_charge = zero(T)
    nb_lapse = zero(T)
    nb_au = zero(T)
    nb_av = zero(T)
    for j in eachindex(g.v)
        scalar = gp2026_initial_scalar_data(
            g.v[j], ep;
            U0,
            V0,
            width,
            M0,
            pulse_leg_gauge,
        )
        nb_radius = max(nb_radius, abs(st.r[i0, j] - scalar.r))
        nb_scalar = max(
            nb_scalar,
            abs(st.phi_re[i0, j] - scalar.psi_re),
            abs(st.phi_im[i0, j] - scalar.psi_im),
        )
        nb_av = max(nb_av, abs(st.Av[i0, j]))
    end
    for j in j0+1:lastindex(g.v)
        dv = g.v[j] - g.v[j - 1]
        midpoint = (g.v[j] + g.v[j - 1]) / 2
        scalar = gp2026_initial_scalar_data(
            midpoint, ep;
            U0,
            V0,
            width,
            M0,
            pulse_leg_gauge,
        )
        _, Jv = current_components(
            scalar.phi_re,
            scalar.phi_im,
            zero(scalar.phi_v_re),
            zero(scalar.phi_v_im),
            scalar.phi_v_re,
            scalar.phi_v_im,
            ep.scalar_charge,
        )
        expected_delta_q = -dv * scalar.r^2 * Jv / 8
        nb_charge = max(
            nb_charge,
            abs(st.Q[i0, j] - st.Q[i0, j - 1] - expected_delta_q),
        )

        rv_left = gp2026_extremal_gauge_rv(
            U0, g.v[j - 1];
            U0,
            V0,
            M0,
            pulse_leg_gauge,
        )
        rv_right = gp2026_extremal_gauge_rv(
            U0, g.v[j];
            U0,
            V0,
            M0,
            pulse_leg_gauge,
        )
        scalar_source = scalar.r *
                        (scalar.phi_v_re^2 + scalar.phi_v_im^2) /
                        (4 * scalar.rv)
        expected_delta_logf =
            log(rv_right / rv_left) + dv * scalar_source
        nb_lapse = max(
            nb_lapse,
            abs(
                st.logf[i0, j] - st.logf[i0, j - 1] -
                expected_delta_logf
            ),
        )

        left = st.Q[i0, j - 1] * exp(st.logf[i0, j - 1]) /
               (4 * st.r[i0, j - 1]^2)
        right = st.Q[i0, j] * exp(st.logf[i0, j]) /
                (4 * st.r[i0, j]^2)
        potential_residual = st.Au[i0, j] - st.Au[i0, j - 1] -
                             dv * (left + right) / 2
        nb_au = max(nb_au, abs(potential_residual))
    end

    du0 = g.u[i0 + 1] - g.u[i0]
    dv0 = g.v[j0 + 1] - g.v[j0]
    av_u = (st.Av[i0 + 1, j0] - st.Av[i0, j0]) / du0
    au_v = (st.Au[i0, j0 + 1] - st.Au[i0, j0]) / dv0
    expected_av_u = (
        -st.Q[i0, j0] * exp(st.logf[i0, j0]) /
        (4 * st.r[i0, j0]^2) -
        st.Q[i0 + 1, j0] * exp(st.logf[i0 + 1, j0]) /
        (4 * st.r[i0 + 1, j0]^2)
    ) / 2
    expected_au_v = (
        st.Q[i0, j0] * exp(st.logf[i0, j0]) /
        (4 * st.r[i0, j0]^2) +
        st.Q[i0, j0 + 1] * exp(st.logf[i0, j0 + 1]) /
        (4 * st.r[i0, j0 + 1]^2)
    ) / 2
    expected_faraday = expected_av_u - expected_au_v
    faraday_corner = av_u - au_v - expected_faraday

    return (
        corner_mass=corner_mass,
        corner_mass_error=corner_mass - M0,
        corner_fcode=exp(st.logf[i0, j0]),
        expected_corner_fcode=fcorner,
        na_radius=na_radius,
        na_lapse_constraint=na_lapse,
        na_charge=na_charge,
        na_scalar=na_scalar,
        na_au=na_au,
        na_av_constraint=na_av,
        nb_radius=nb_radius,
        nb_scalar=nb_scalar,
        nb_charge_constraint=nb_charge,
        nb_lapse_constraint=nb_lapse,
        nb_au_constraint=nb_au,
        nb_av=nb_av,
        faraday_corner=faraday_corner,
    )
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
            r_u_compact = (
                areal_radius(g.u[i + 1], g.v[j], p) -
                areal_radius(g.u[i], g.v[j], p)
            ) / (g.u[i + 1] - g.u[i])
        elseif i == lastindex(g.u)
            q_u_compact = (st.Q[i, j] - st.Q[i - 1, j]) / (g.u[i] - g.u[i - 1])
            r_u_compact = (
                areal_radius(g.u[i], g.v[j], p) -
                areal_radius(g.u[i - 1], g.v[j], p)
            ) / (g.u[i] - g.u[i - 1])
        else
            q_u_compact = (st.Q[i + 1, j] - st.Q[i - 1, j]) / (g.u[i + 1] - g.u[i - 1])
            r_u_compact = (
                areal_radius(g.u[i + 1], g.v[j], p) -
                areal_radius(g.u[i - 1], g.v[j], p)
            ) / (g.u[i + 1] - g.u[i - 1])
        end

        if j == firstindex(g.v)
            q_vef = (st.Q[i, j + 1] - st.Q[i, j]) / (vef[j + 1] - vef[j])
        elseif j == lastindex(g.v)
            q_vef = (st.Q[i, j] - st.Q[i, j - 1]) / (vef[j] - vef[j - 1])
        else
            q_vef = (st.Q[i, j + 1] - st.Q[i, j - 1]) / (vef[j + 1] - vef[j - 1])
        end

        q_r = q_u_compact / r_u_compact
        rho[j] = (q_r + 0.5 * q_vef) / (4pi * rplus^2)
    end

    return vef, rho
end
