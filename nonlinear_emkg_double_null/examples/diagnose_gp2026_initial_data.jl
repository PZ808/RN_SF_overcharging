using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, argument(index, string(default)))
end

function finite_difference(values, coordinates)
    length(values) == length(coordinates) ||
        throw(ArgumentError("values and coordinates must have equal length"))
    length(values) >= 2 || throw(ArgumentError("at least two points are required"))
    derivative = similar(values)
    derivative[begin] = (values[begin + 1] - values[begin]) /
                        (coordinates[begin + 1] - coordinates[begin])
    derivative[end] = (values[end] - values[end - 1]) /
                      (coordinates[end] - coordinates[end - 1])
    for j in firstindex(values)+1:lastindex(values)-1
        derivative[j] = (values[j + 1] - values[j - 1]) /
                        (coordinates[j + 1] - coordinates[j - 1])
    end
    return derivative
end

function rn_F(r, M, Q)
    return one(r) - 2M / r + Q^2 / r^2
end

function extremal_F_prime(r, M)
    return 2M / r^2 - 2M^2 / r^3
end

function initial_reduced_scalar_and_derivative(V, ep, M0, U0, V0, width)
    r = gp2026_extremal_gauge_initial_radius(U0, V; U0, V0, M0)
    rv = gp2026_extremal_gauge_rv(U0, V; U0, V0, M0)
    amplitude = gp2026_single_pulse_envelope(V; amplitude=ep.amplitude, width)
    derivative = N.gp2026_single_pulse_envelope_derivative(V;
                                                           amplitude=ep.amplitude,
                                                           width)
    phase = ep.omega * V
    psi_re = sqrt(32pi) * amplitude * cos(phase)
    psi_im = -sqrt(32pi) * amplitude * sin(phase)
    psi_v_re = sqrt(32pi) * (derivative * cos(phase) -
                             ep.omega * amplitude * sin(phase))
    psi_v_im = sqrt(32pi) * (-derivative * sin(phase) -
                             ep.omega * amplitude * cos(phase))
    phi_re = psi_re / r
    phi_im = psi_im / r
    phi_v_re = (psi_v_re - rv * phi_re) / r
    phi_v_im = (psi_v_im - rv * phi_im) / r
    return (; r, rv, psi_re, psi_im, psi_v_re, psi_v_im,
            phi_re, phi_im, phi_v_re, phi_v_im)
end

function diagnose_na(ep; U0, V0, U1, M0, samples)
    q0 = ep.rn.Q0
    fcorner = gp2026_fcorner_code(ep; U0, V0, M0)
    u = collect(range(U0, U1; length=samples))
    mass_with_rn_rv = similar(u)
    mass_with_extremal_rv = similar(u)
    rv_rn = similar(u)
    rv_ext = similar(u)
    av_constraint_error = similar(u)
    for (k, U) in pairs(u)
        point = gp2026_na_boundary_point(U, ep; U0, V0, M0)
        r = point.r
        f = exp(point.logf)
        ru = gp2026_extremal_gauge_ru(U; M0)
        rv_rn[k] = f * rn_F(r, M0, q0) / 2
        rv_ext[k] = gp2026_extremal_gauge_rv(U, V0; U0, V0, M0)
        mass_with_rn_rv[k] = renormalized_hawking_mass(r, f, ru, rv_rn[k], q0)
        mass_with_extremal_rv[k] = renormalized_hawking_mass(r, f, ru, rv_ext[k], q0)
        av_u = q0 * f * ru / (2r^2)
        av_constraint_error[k] = av_u + q0 * f / (4r^2)
    end
    return (;
        fcorner,
        max_abs_mass_error_rn_rv=maximum(abs.(mass_with_rn_rv .- M0)),
        max_abs_mass_error_extremal_rv=maximum(abs.(mass_with_extremal_rv .- M0)),
        max_abs_rv_difference=maximum(abs.(rv_ext .- rv_rn)),
        max_abs_av_constraint_error=maximum(abs, av_constraint_error),
        first_sample=(U=first(u), rv_rn=first(rv_rn), rv_ext=first(rv_ext),
                      mass_rn=first(mass_with_rn_rv),
                      mass_ext=first(mass_with_extremal_rv)),
        last_sample=(U=last(u), rv_rn=last(rv_rn), rv_ext=last(rv_ext),
                     mass_rn=last(mass_with_rn_rv),
                     mass_ext=last(mass_with_extremal_rv)),
    )
end

function diagnose_nb(::Type{T}, ep; U0, V0, Vmax, dv, M0, width) where {T<:Real}
    nv = Int(round(Vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0, U1=U0 + parse(T, "0.01"), V1=Vmax)
    state = NLState(grid)
    initialize_gp2026_single_pulse!(state, grid, ep; U0, V0, width, M0)

    v = grid.v
    logf = collect(state.logf[1, :])
    q = collect(state.Q[1, :])
    q_v = finite_difference(q, v)
    logf_v = finite_difference(logf, v)

    expected_q_v = similar(v)
    expected_logf_v = similar(v)
    for (j, V) in pairs(v)
        data = initial_reduced_scalar_and_derivative(V, ep, M0, U0, V0, width)
        Ju, Jv = current_components(data.phi_re, data.phi_im,
                                    zero(data.phi_v_re), zero(data.phi_v_im),
                                    data.phi_v_re, data.phi_v_im,
                                    ep.scalar_charge)
        expected_q_v[j] = -data.r^2 * Jv / 8
        rv_v = extremal_F_prime(data.r, M0) * data.rv / 2
        phiv2 = data.phi_v_re^2 + data.phi_v_im^2
        expected_logf_v[j] = rv_v / data.rv + data.r * phiv2 / (4data.rv)
    end

    active = findall(V -> zero(T) < V < width, v)
    all_indices = eachindex(v)
    return (;
        max_abs_q_v_residual=maximum(abs, q_v .- expected_q_v),
        max_abs_logf_v_residual=maximum(abs, logf_v .- expected_logf_v),
        active_max_abs_q_v_residual=isempty(active) ? zero(T) :
                                    maximum(abs, (q_v .- expected_q_v)[active]),
        active_max_abs_logf_v_residual=isempty(active) ? zero(T) :
                                       maximum(abs, (logf_v .- expected_logf_v)[active]),
        final_Q=last(q),
        final_logf=last(logf),
        min_Q=minimum(q[all_indices]),
        max_Q=maximum(q[all_indices]),
    )
end

function run(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "400.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    omega = real_argument(5, "1.0", T)
    eQ0 = real_argument(6, "0.6", T)
    U0 = real_argument(7, "-1.0", T)
    V0 = real_argument(8, "0.0", T)
    U1 = real_argument(9, "1.6", T)
    width = real_argument(10, "20.0", T)
    na_samples = integer_argument(11, 200)

    ep = EvolutionParams(rn=RNParams(one(T), q0),
                         scalar_charge=eQ0 / q0,
                         amplitude=amplitude,
                         omega=omega)
    na = diagnose_na(ep; U0, V0, U1, M0=one(T), samples=na_samples)
    nb = diagnose_nb(T, ep; U0, V0, Vmax=vmax, dv, M0=one(T), width)

    println("# GP2026 Appendix-A initial-data diagnostic")
    println("# Q0 = ", q0, ", eQ0 = ", eQ0, ", A0 = ", amplitude,
            ", omega = ", omega)
    println("# Vmax = ", vmax, ", Delta V = ", dv,
            ", U0 = ", U0, ", U1 = ", U1, ", width = ", width)
    println("# N_A")
    println("fcorner_code = ", na.fcorner)
    println("max |M[rV=f F_Q/2] - M0| = ", na.max_abs_mass_error_rn_rv)
    println("max |M[rV=F_ext/2] - M0| = ", na.max_abs_mass_error_extremal_rv)
    println("max |rV_ext - rV_mass_compatible| = ", na.max_abs_rv_difference)
    println("max |A_V,U + Q f/(4 r^2)| = ", na.max_abs_av_constraint_error)
    println("first N_A sample = ", na.first_sample)
    println("last N_A sample = ", na.last_sample)
    println("# N_B")
    println("max |Q_V + r^2 J_V/8| = ", nb.max_abs_q_v_residual)
    println("max active-pulse |Q_V + r^2 J_V/8| = ",
            nb.active_max_abs_q_v_residual)
    println("max |(log f)_V - expected| = ", nb.max_abs_logf_v_residual)
    println("max active-pulse |(log f)_V - expected| = ",
            nb.active_max_abs_logf_v_residual)
    println("final Q = ", nb.final_Q)
    println("final logf = ", nb.final_logf)
    println("Q range = ", (nb.min_Q, nb.max_Q))
end

run(Float64)
