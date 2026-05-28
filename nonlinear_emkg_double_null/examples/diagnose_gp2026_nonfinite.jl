using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, length(ARGS) >= index ? ARGS[index] : default)
end

function integer_argument(index, default)
    return length(ARGS) >= index ? parse(Int, ARGS[index]) : default
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q)
end

function first_bad_index(name, values)
    index = findfirst(!isfinite, values)
    isnothing(index) || println("first nonfinite ", name, " index = ", index)
    return index
end

function row_summary(row::NLRow)
    return (
        u = row.u,
        r_extrema = extrema(row.r),
        logf_extrema = extrema(row.logf),
        psi_abs_max = maximum(hypot.(row.phi_re, row.phi_im)),
        Au_extrema = extrema(row.Au),
        Av_extrema = extrema(row.Av),
        Q_extrema = extrema(row.Q),
    )
end

function averaged_rv_crossing(previous::NLRow, current::NLRow)
    dv = diff(current.v)
    rv = [((current.r[j + 1] - current.r[j]) +
           (previous.r[j + 1] - previous.r[j])) / (2dv[j])
          for j in 1:length(dv)]
    index = findfirst(value -> value <= 0, rv)
    return isnothing(index) ? nothing :
           (index = index, V = current.v[index + 1], r = current.r[index + 1],
            rv = rv[index])
end

function trace_cell_step!(st::NLState, g::Grid, ep::EvolutionParams, j::Int;
                          iterations::Int=10, hyperbolic_charge::Bool=true)
    du = g.u[2] - g.u[1]
    dv = g.v[j + 1] - g.v[j]
    e = ep.scalar_charge

    r00, r10, r01 = st.r[1, j], st.r[2, j], st.r[1, j + 1]
    lf00, lf10, lf01 = st.logf[1, j], st.logf[2, j], st.logf[1, j + 1]
    pr00, pr10, pr01 = st.phi_re[1, j], st.phi_re[2, j], st.phi_re[1, j + 1]
    pi00, pi10, pi01 = st.phi_im[1, j], st.phi_im[2, j], st.phi_im[1, j + 1]
    au00, au10, au01 = st.Au[1, j], st.Au[2, j], st.Au[1, j + 1]
    av00, av10, av01 = st.Av[1, j], st.Av[2, j], st.Av[1, j + 1]
    q00, q10, q01 = st.Q[1, j], st.Q[2, j], st.Q[1, j + 1]

    r11 = r10 + r01 - r00
    lf11 = lf10 + lf01 - lf00
    pr11 = pr10 + pr01 - pr00
    pi11 = pi10 + pi01 - pi00
    au11 = au10 - au01 + au00
    av11 = av01 - av10 + av00
    q11 = q10 + q01 - q00

    last_report = nothing
    for iteration in 1:iterations
        r = N.corner_average(r00, r10, r01, r11)
        lf = N.corner_average(lf00, lf10, lf01, lf11)
        f = exp(lf)
        pr = N.corner_average(pr00, pr10, pr01, pr11)
        pii = N.corner_average(pi00, pi10, pi01, pi11)
        au = N.corner_average(au00, au10, au01, au11)
        av = N.corner_average(av00, av10, av01, av11)
        q = N.corner_average(q00, q10, q01, q11)

        ru = N.corner_du(r00, r10, r01, r11, du)
        rv = N.corner_dv(r00, r10, r01, r11, dv)
        pru = N.corner_du(pr00, pr10, pr01, pr11, du)
        prv = N.corner_dv(pr00, pr10, pr01, pr11, dv)
        piu = N.corner_du(pi00, pi10, pi01, pi11, du)
        piv = N.corner_dv(pi00, pi10, pi01, pi11, dv)

        source = stress_energy_reduced_scalar(r, f, q, ru, rv, pr, pii, pru, prv,
                                              piu, piv, au, av, e)
        ruv, lfuv = metric_rhs(r, f, ru, rv, q, source)
        pruv, piuv = charged_reduced_scalar_rhs(r, ruv, pru, prv, piu, piv,
                                                pr, pii, au, av, e)
        auv, avu, _, quc, qvc = maxwell_rhs(r, f, q, source)
        quv = hyperbolic_charge ?
              charged_reduced_charge_rhs(r, f, q, pr, pii, pru, prv, piu, piv,
                                         au, av, e) :
              zero(q)

        next_r11 = r10 + r01 - r00 + du * dv * ruv
        next_lf11 = lf10 + lf01 - lf00 + du * dv * lfuv
        next_pr11 = pr10 + pr01 - pr00 + du * dv * pruv
        next_pi11 = pi10 + pi01 - pi00 + du * dv * piuv
        next_au11 = au10 - au01 + au00 + 2dv * auv
        next_av11 = av01 - av10 + av00 + 2du * avu
        next_q11 = hyperbolic_charge ? q10 + q01 - q00 + du * dv * quv :
                   q01 + du * quc

        last_report = (
            iteration = iteration,
            center = (r = r, logf = lf, f = f, q = q, ru = ru, rv = rv,
                      psi = hypot(pr, pii), Au = au, Av = av),
            source = (Tuu = source.Tuu, Tvv = source.Tvv, Tuv = source.Tuv,
                      Tthth = source.Tthth, Ju = source.Ju, Jv = source.Jv),
            rhs = (ruv = ruv, lfuv = lfuv, psi_re_uv = pruv, psi_im_uv = piuv,
                   Au_v = auv, Av_u = avu, Q_u_constraint = quc,
                   Q_v_constraint = qvc, Q_uv = quv),
            update = (r11 = next_r11, logf11 = next_lf11, phi_re11 = next_pr11,
                      phi_im11 = next_pi11, Au11 = next_au11, Av11 = next_av11,
                      Q11 = next_q11),
        )

        r11, lf11 = next_r11, next_lf11
        pr11, pi11 = next_pr11, next_pi11
        au11, av11 = next_au11, next_av11
        q11 = next_q11

        values = (r, lf, f, q, ru, rv, source.Tuu, source.Tvv, source.Ju,
                  source.Jv, ruv, lfuv, pruv, piuv, auv, avu, quc, qvc, quv,
                  r11, lf11, pr11, pi11, au11, av11, q11)
        if any(!isfinite, values)
            println("nonfinite during cell iteration ", iteration)
            println(last_report)
            break
        end
        if r11 <= 0
            println("nonpositive areal-radius update during cell iteration ", iteration)
            println(last_report)
            break
        end
    end

    st.r[2, j + 1] = r11
    st.logf[2, j + 1] = lf11
    st.phi_re[2, j + 1] = pr11
    st.phi_im[2, j + 1] = pi11
    st.Au[2, j + 1] = au11
    st.Av[2, j + 1] = av11
    st.Q[2, j + 1] = q11
    return last_report
end

function run_trace(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "100.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    Umax = real_argument(6, "1.6", T)
    max_rows = integer_argument(7, 100_000)
    charge_mode = length(ARGS) >= 9 ? ARGS[9] : "hyperbolic"
    hyperbolic_charge = charge_mode == "hyperbolic"
    charge_mode in ("hyperbolic", "constraint") ||
        throw(ArgumentError("charge mode must be hyperbolic or constraint"))
    step_control = length(ARGS) >= 10 ? ARGS[10] : "outer"
    step_control in ("outer", "max-row") ||
        throw(ArgumentError("step control must be outer or max-row"))

    U0 = parse(T, "-1.0")
    ep = EvolutionParams(
        rn = RNParams(one(T), q0),
        scalar_charge = parse(T, "0.6") / q0,
        amplitude = amplitude,
        omega = one(T),
    )
    nv = Int(round(vmax / dv)) + 1
    seed_grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T),
                            U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(seed_grid)
    initialize_gp2026_single_pulse!(seed, seed_grid, ep)
    rows = [row_from_rectangular(seed, seed_grid, 1)]

    println("Tracing GP2026 row march")
    println("numeric type = ", T)
    println("precision bits = ", precision(T))
    println("Q0 = ", q0, ", eQ0 = ", ep.scalar_charge * q0)
    println("Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C)
    println("charge evolution = ", charge_mode)
    println("step control = ", step_control)

    first_horizon = nothing
    while last(rows).u < Umax && length(rows) < max_rows
        previous = last(rows)
        fcode_outer = exp(last(previous.logf))
        fcode_step = step_control == "outer" ? fcode_outer : exp(maximum(previous.logf))
        isfinite(fcode_step) && fcode_step > 0 || error("bad f before row")
        du = 2C / fcode_step
        next_u = min(Umax, previous.u + du)
        next_u > previous.u || error("coordinate stalled before row")
        south = gp2026_na_boundary_point(next_u, ep; U0, V0=zero(T), M0=one(T))
        grid = Grid([previous.u, next_u], previous.v)
        st = NLState(grid)
        st.r[1, :] .= previous.r
        st.logf[1, :] .= previous.logf
        st.phi_re[1, :] .= previous.phi_re
        st.phi_im[1, :] .= previous.phi_im
        st.Au[1, :] .= previous.Au
        st.Av[1, :] .= previous.Av
        st.Q[1, :] .= previous.Q
        st.r[2, 1] = south.r
        st.logf[2, 1] = south.logf
        st.phi_re[2, 1] = south.phi_re
        st.phi_im[2, 1] = south.phi_im
        st.Au[2, 1] = south.Au
        st.Av[2, 1] = south.Av
        st.Q[2, 1] = south.Q

        bad_j = nothing
        bad_report = nothing
        for j in 1:length(previous.v)-1
            bad_report = trace_cell_step!(st, grid, ep, j; iterations=10,
                                          hyperbolic_charge)
            if !isfinite(st.r[2, j + 1]) || !isfinite(st.logf[2, j + 1]) ||
               !isfinite(st.phi_re[2, j + 1]) || !isfinite(st.phi_im[2, j + 1]) ||
               !isfinite(st.Au[2, j + 1]) || !isfinite(st.Av[2, j + 1]) ||
               !isfinite(st.Q[2, j + 1]) || st.r[2, j + 1] <= 0
                bad_j = j
                break
            end
        end

        next = row_from_rectangular(st, grid, 2)
        horizon = averaged_rv_crossing(previous, next)
        if isnothing(first_horizon) && !isnothing(horizon)
            first_horizon = (row = length(rows) + 1, U = next.u,
                             V = horizon.V, r = horizon.r, rv = horizon.rv)
            println("first apparent-horizon crossing = ", first_horizon)
        end

        if length(rows) == 1 || length(rows) % 10 == 0 || !finite_row(next) ||
           !isnothing(horizon)
            println("row ", length(rows) + 1, " summary: du=", du,
                    ", next_u=", next.u, ", f_outer=", fcode_outer,
                    ", f_step=", fcode_step,
                    ", ", row_summary(next))
        end

        if !isnothing(bad_j) || !finite_row(next)
            println("stopping at first invalid/nonfinite row")
            println("row index = ", length(rows) + 1)
            println("cell j = ", bad_j, ", V cell = ",
                    isnothing(bad_j) ? nothing : (grid.v[bad_j], grid.v[bad_j + 1]))
            first_bad_index("r", next.r)
            first_bad_index("logf", next.logf)
            first_bad_index("phi_re", next.phi_re)
            first_bad_index("phi_im", next.phi_im)
            first_bad_index("Au", next.Au)
            first_bad_index("Av", next.Av)
            first_bad_index("Q", next.Q)
            println("previous row summary: ", row_summary(previous))
            println("bad row summary: ", row_summary(next))
            println("last cell report: ", bad_report)
            return
        end
        push!(rows, next)
    end
    println("completed without a nonfinite row; rows = ", length(rows))
end

precision_bits = integer_argument(8, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_trace(BigFloat)
    end
else
    run_trace(Float64)
end
