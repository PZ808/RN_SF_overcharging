using NonlinearEMKGDoubleNull

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
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function run_paper_amr(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "400.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    Umax = real_argument(6, "1.6", T)
    max_rows = integer_argument(7, 120)

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
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax, C, iterations=10,
                                       max_rows, hyperbolic_charge=true,
                                       step_control=:outer)

    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    final = last(rows)
    u_values = [row.u for row in rows]
    du = diff(u_values)
    outer_f = [exp(last(row.logf)) for row in rows]
    next_du = 2C / last(outer_f)
    stalled = final.u + next_du == final.u
    hit_invalid_row = last_valid < length(evolved.rows)

    println("Gelles-Pretorius 2026 paper-style Eq. 9 AMR")
    println("numeric type = ", T)
    println("precision bits = ", precision(T))
    println("Q0 = ", q0, ", e Q0 = ", ep.scalar_charge * q0)
    println("A0 = ", amplitude)
    println("Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C)
    println("requested Umax = ", Umax)
    println("stored rows = ", length(evolved.rows))
    println("valid rows = ", length(rows))
    println("last valid U = ", final.u)
    println("reached Umax = ", final.u == Umax)
    println("encountered invalid row = ", hit_invalid_row)
    println("coordinate stalled at precision = ", stalled)
    println("Delta U first/last = ", isempty(du) ? nothing : (first(du), last(du)))
    println("Delta U extrema = ", isempty(du) ? nothing : extrema(du))
    println("outer f_code first/last = ", (first(outer_f), last(outer_f)))
    println("next Eq. 9 Delta U = ", next_du)

    throat = throat_row_diagnostics(final)
    println("throat min(r-|Q|) = ", throat.min_y)
    println("throat max rho = ", throat.max_rho)
    println("throat max |Delta rho| = ", throat.max_abs_delta_rho)
    println("interpretation = ",
            final.u == Umax ? "reached requested Umax" :
            hit_invalid_row ? "invalid row before clean Eq. 9 accumulation" :
            "Eq. 9 accumulated near a limiting U before reaching Umax")
end

precision_bits = integer_argument(8, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_paper_amr(BigFloat)
    end
else
    run_paper_amr(Float64)
end
