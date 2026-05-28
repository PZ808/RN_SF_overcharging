using NonlinearEMKGDoubleNull

function real_argument(index, default, ::Type{T}) where {T<:Real}
    return parse(T, length(ARGS) >= index ? ARGS[index] : default)
end

function integer_argument(index, default)
    return length(ARGS) >= index ? parse(Int, ARGS[index]) : default
end

function first_trapped_slice(state)
    for j in 2:length(state.slices)
        rv = adaptive_outgoing_expansion(state.slices[j - 1], state.slices[j])
        crossing = findfirst(value -> value <= 0, rv)
        isnothing(crossing) || return (state.slices[j].v, state.slices[j].u[crossing])
    end
    return nothing
end

function run_comparison(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "20.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    Umax = real_argument(6, "1.6", T)
    max_rows = integer_argument(7, 100_000)
    charge_mode = length(ARGS) >= 9 ? ARGS[9] : "hyperbolic"
    hyperbolic_charge = charge_mode == "hyperbolic"
    charge_mode in ("hyperbolic", "constraint") ||
        throw(ArgumentError("charge mode must be hyperbolic or constraint"))
    step_control_argument = length(ARGS) >= 10 ? ARGS[10] : "outer"
    step_control_argument in ("outer", "max-row") ||
        throw(ArgumentError("step control must be outer or max-row"))
    step_control = step_control_argument == "outer" ? :outer : :max_row

    U0 = parse(T, "-1.0")
    ep = EvolutionParams(
        rn = RNParams(one(T), q0),
        scalar_charge = parse(T, "0.6") / q0,
        amplitude = amplitude,
        omega = one(T),
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T),
                       U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    raw = evolve_gp2026_u_adaptive(initial, ep; Umax, C, iterations=10, max_rows,
                                   hyperbolic_charge, step_control)

    last_finite = findlast(row -> all(isfinite, row.r) && all(isfinite, row.logf) &&
                                   all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
                                   all(isfinite, row.Au) && all(isfinite, row.Av) &&
                                   all(isfinite, row.Q), raw.rows)
    isnothing(last_finite) && error("initial GP row is nonfinite")
    finite_rows = UAdaptiveNLState(raw.rows[1:last_finite])
    state = adaptive_state_from_u_rows(finite_rows)
    du = diff([row.u for row in finite_rows.rows])
    trap = length(finite_rows.rows) >= 2 ? first_trapped_slice(state) : nothing
    last_fcode = exp(last(finite_rows.rows).logf[end])
    next_paper_du = 2C / last_fcode
    hit_nonfinite_row = last_finite < length(raw.rows)
    coordinate_stalled = last(finite_rows.rows).u + next_paper_du ==
                         last(finite_rows.rows).u

    println("Gelles-Pretorius 2026 paper U-step evolution")
    println("numeric type = ", T)
    println("precision bits = ", precision(T))
    println("Q0 = ", q0)
    println("e Q0 = ", ep.scalar_charge * q0)
    println("A0 = ", amplitude)
    println("requested Vmax = ", vmax)
    println("Delta V = ", dv)
    println("C = ", C)
    println("charge evolution = ", charge_mode)
    println("step control = ", step_control_argument)
    println("requested Umax = ", Umax)
    println("stored U rows = ", length(raw.rows))
    println("finite U rows = ", length(finite_rows.rows))
    println("last finite U = ", last(finite_rows.rows).u)
    println("reached Umax = ", last(finite_rows.rows).u == Umax)
    println("encountered nonfinite row = ", hit_nonfinite_row)
    println("coordinate stalled at precision = ", coordinate_stalled)
    println("Delta U extrema = ", isempty(du) ? nothing : extrema(du))
    println("last outer f_code = ", last_fcode)
    println("next paper Delta U = ", next_paper_du)
    println("first trapped (V,U-cell) = ", trap)

    if length(finite_rows.rows) >= 3
        target_v = vmax - dv / 2
        target_u = min(parse(T, "-0.5"), last(finite_rows.rows).u)
        _, _, _, _, q_u_residual =
            charged_charge_flux_u_profile(state, ep; target_v, reduced_scalar=true)
        _, _, _, _, q_v_residual =
            charged_charge_flux_v_profile(state, ep; target_u, reduced_scalar=true)
        horizon_v, _, horizon_rphi = gp2026_horizon_rphi_series(state)
        println("diagnostic U for Q_V = ", target_u)
        println("max |Q_U - flux source| = ", maximum(abs, q_u_residual))
        println("max |Q_V - flux source| = ", maximum(abs, q_v_residual))
        if !isempty(horizon_rphi)
            println("horizon |r phi_GP| first/last = ",
                    (first(horizon_rphi), last(horizon_rphi)))
            println("horizon amplitude samples and V range = ",
                    (length(horizon_rphi), first(horizon_v), last(horizon_v)))
        end
    end
end

precision_bits = integer_argument(8, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_comparison(BigFloat)
    end
else
    run_comparison(Float64)
end
