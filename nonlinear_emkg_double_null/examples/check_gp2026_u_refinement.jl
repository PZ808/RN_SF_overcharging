using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

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
    step_control_argument = length(ARGS) >= 10 ? ARGS[10] : "local"
    step_control_argument in ("outer", "max-row", "geometric", "throat", "local") ||
        throw(ArgumentError("step control must be outer, max-row, geometric, throat, or local"))
    step_control = if step_control_argument == "outer"
        :outer
    elseif step_control_argument == "max-row"
        :max_row
    elseif step_control_argument == "geometric"
        :geometric
    elseif step_control_argument == "throat"
        :throat
    else
        :local
    end
    max_delta_rho = real_argument(11, "0.25", T)
    match_rho = real_argument(12, "2.0", T)

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
                                   hyperbolic_charge, step_control, max_delta_rho)

    last_valid = findlast(row -> all(isfinite, row.r) && all(isfinite, row.logf) &&
                                  all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
                                  all(isfinite, row.Au) && all(isfinite, row.Av) &&
                                  all(isfinite, row.Q) &&
                                  all(>(zero(eltype(row.r))), row.r), raw.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    valid_rows = UAdaptiveNLState(raw.rows[1:last_valid])
    final_row = last(valid_rows.rows)
    state = adaptive_state_from_u_rows(valid_rows)
    du = diff([row.u for row in valid_rows.rows])
    trap = length(valid_rows.rows) >= 2 ? first_trapped_slice(state) : nothing
    last_fcode = exp(last(final_row.logf))
    max_fcode = exp(maximum(final_row.logf))
    next_paper_du = 2C / last_fcode
    next_max_row_du = 2C / max_fcode
    next_geometric_du = N.geometric_row_du(valid_rows.rows, C)
    next_throat_du = throat_row_du(valid_rows.rows, C; max_delta_rho)
    next_controlled_du = if step_control === :outer
        next_paper_du
    elseif step_control === :max_row
        next_max_row_du
    elseif step_control === :geometric
        next_geometric_du
    elseif step_control === :throat
        next_throat_du
    else
        minimum((next_max_row_du, next_geometric_du, next_throat_du))
    end
    hit_invalid_row = last_valid < length(raw.rows)
    coordinate_stalled = final_row.u + next_controlled_du == final_row.u

    println("Gelles-Pretorius 2026 U-step evolution")
    println("numeric type = ", T)
    println("precision bits = ", precision(T))
    println("Q0 = ", q0)
    println("e Q0 = ", ep.scalar_charge * q0)
    println("A0 = ", amplitude)
    println("requested Vmax = ", vmax)
    println("Delta V = ", dv)
    println("C = ", C)
    println("max Delta rho = ", max_delta_rho)
    println("matching rho = ", match_rho)
    println("charge evolution = ", charge_mode)
    println("step control = ", step_control_argument)
    println("requested Umax = ", Umax)
    println("stored U rows = ", length(raw.rows))
    println("valid U rows = ", length(valid_rows.rows))
    println("last valid U = ", final_row.u)
    println("reached Umax = ", final_row.u == Umax)
    println("encountered invalid row = ", hit_invalid_row)
    println("coordinate stalled at precision = ", coordinate_stalled)
    println("Delta U extrema = ", isempty(du) ? nothing : extrema(du))
    println("last outer f_code = ", last_fcode)
    println("last max-row f_code = ", max_fcode)
    println("next paper Delta U = ", next_paper_du)
    println("next max-row Delta U = ", next_max_row_du)
    println("next geometric Delta U = ", next_geometric_du)
    println("next throat Delta U = ", next_throat_du)
    println("next controlled Delta U = ", next_controlled_du)
    println("first trapped (V,U-cell) = ", trap)
    throat = throat_row_diagnostics(final_row)
    println("throat min(r-|Q|) = ", throat.min_y)
    println("throat max rho = ", throat.max_rho)
    println("throat max |Delta rho| = ", throat.max_abs_delta_rho)
    match = throat_matching_candidate(final_row; rho_min=match_rho)
    println("throat matching candidate = ", match)
    println("throat matching band = ", throat_matching_band(final_row; rho_min=match_rho))
    rho_lapse = rho_lapse_diagnostics(final_row; rho_min=match_rho)
    println("logf range = ", rho_lapse.logf_range)
    println("logf_rho range = ", rho_lapse.logf_rho_range)
    println("logf range width = ", range_width(rho_lapse.logf_range))
    println("logf_rho range width = ", range_width(rho_lapse.logf_rho_range))
    println("throat logf range = ", rho_lapse.throat_logf_range)
    println("throat logf_rho range = ", rho_lapse.throat_logf_rho_range)
    println("throat logf range width = ", range_width(rho_lapse.throat_logf_range))
    println("throat logf_rho range width = ", range_width(rho_lapse.throat_logf_rho_range))

    if length(valid_rows.rows) >= 3
        target_v = vmax - dv / 2
        target_u = min(parse(T, "-0.5"), final_row.u)
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
