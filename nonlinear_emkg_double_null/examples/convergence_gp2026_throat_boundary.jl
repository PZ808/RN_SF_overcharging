using NonlinearEMKGDoubleNull

const SQRT32PI = sqrt(32 * pi)

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : default
end

function real_argument(index, default)
    return parse(Float64, argument(index, string(default)))
end

function integer_argument(index, default)
    return parse(Int, argument(index, string(default)))
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function refinement_rate(previous, current)
    (previous > 0 && current > 0) ? log(previous / current) / log(2) : NaN
end

function interpolate_series(x, y, xq)
    if xq < first(x) || xq > last(x)
        return nothing
    end
    i = searchsortedlast(x, xq)
    i == length(x) && return y[end]
    i < 1 && return nothing
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - t) * y[i] + t * y[i + 1]
end

function build_series(; q0, rho_match, vmax, dv, C, max_rows, Umax, amplitude)
    ep = EvolutionParams(
        rn = RNParams(1.0, q0),
        scalar_charge = 0.6 / q0,
        amplitude = amplitude,
        omega = 1.0,
    )
    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0=-1.0, V0=0.0, U1=-0.99, V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    evolved = evolve_gp2026_u_adaptive(
        row_from_rectangular(seed, grid, 1), ep;
        Umax, C, iterations=10, max_rows, hyperbolic_charge=true,
        step_control=:outer,
    )
    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    samples = throat_boundary_series(rows; rho_match, boundary=:outer)
    return ep, samples, rows
end

function sample_arrays(samples, ep)
    u = [sample.u for sample in samples]
    v = [sample.v for sample in samples]
    r = [sample.r for sample in samples]
    q = [sample.q for sample in samples]
    amplitude = [hypot(sample.phi_re, sample.phi_im) / SQRT32PI
                 for sample in samples]
    q_v_residual = [throat_boundary_observables(sample, ep).q_v_residual
                    for sample in samples]
    return (; u, v, r, q, amplitude, q_v_residual)
end

function compare_series(coarse, fine)
    common_u = [u for u in coarse.u if first(fine.u) <= u <= last(fine.u)]
    isempty(common_u) && return nothing

    function max_difference(field)
        coarse_values = getproperty(coarse, field)
        fine_values = getproperty(fine, field)
        differences = Float64[]
        for u in common_u
            coarse_value = interpolate_series(coarse.u, coarse_values, u)
            fine_value = interpolate_series(fine.u, fine_values, u)
            if !isnothing(coarse_value) && !isnothing(fine_value)
                push!(differences, abs(coarse_value - fine_value))
            end
        end
        return maximum(differences)
    end

    return (
        samples = length(common_u),
        u_min = first(common_u),
        u_max = last(common_u),
        max_delta_v = max_difference(:v),
        max_delta_r = max_difference(:r),
        max_delta_q = max_difference(:q),
        max_delta_amplitude = max_difference(:amplitude),
    )
end

function print_level_table(levels, arrays, rows, samples)
    println("level\tDelta_V\tC\tvalid_rows\tsamples\tU_first\tU_last\tmax_abs_QV_residual")
    for i in eachindex(levels)
        arr = arrays[i]
        println(join((
            i - 1,
            levels[i].dv,
            levels[i].C,
            length(rows[i]),
            length(samples[i]),
            first(arr.u),
            last(arr.u),
            maximum(abs, arr.q_v_residual),
        ), '\t'))
    end
end

function print_pair_table(arrays)
    println()
    println("pair\tcommon_samples\tUmin\tUmax\tmax_dV\trate_dV\tmax_dr\trate_dr\tmax_dQ\trate_dQ\tmax_d_abs_rphi\trate_d_abs_rphi")
    previous = nothing
    for i in 1:length(arrays)-1
        comparison = compare_series(arrays[i], arrays[i + 1])
        if isnothing(comparison)
            println("$(i - 1)-$i\tno-overlap")
            continue
        end
        rate_v = isnothing(previous) ? NaN :
                 refinement_rate(previous.max_delta_v, comparison.max_delta_v)
        rate_r = isnothing(previous) ? NaN :
                 refinement_rate(previous.max_delta_r, comparison.max_delta_r)
        rate_q = isnothing(previous) ? NaN :
                 refinement_rate(previous.max_delta_q, comparison.max_delta_q)
        rate_amplitude = isnothing(previous) ? NaN :
                         refinement_rate(previous.max_delta_amplitude,
                                         comparison.max_delta_amplitude)
        println(join((
            "$(i - 1)-$i",
            comparison.samples,
            comparison.u_min,
            comparison.u_max,
            comparison.max_delta_v,
            rate_v,
            comparison.max_delta_r,
            rate_r,
            comparison.max_delta_q,
            rate_q,
            comparison.max_delta_amplitude,
            rate_amplitude,
        ), '\t'))
        previous = comparison
    end
end

function main()
    q0 = real_argument(1, 1.0033218)
    rho_match = real_argument(2, 2.0)
    vmax = real_argument(3, 400.0)
    base_dv = real_argument(4, 0.08)
    base_C = real_argument(5, 0.6)
    base_rows = integer_argument(6, 80)
    levels_count = integer_argument(7, 3)
    Umax = real_argument(8, 1.6)
    amplitude = real_argument(9, 0.01)

    levels = [
        (dv = base_dv / 2.0^level,
         C = base_C / 2.0^level,
         max_rows = base_rows * 2^level)
        for level in 0:levels_count-1
    ]

    println("# GP2026 fixed-rho throat-boundary convergence")
    println("# Q0 = ", q0, ", eQ0 = 0.6, A0 = ", amplitude,
            ", rho_match = ", rho_match)
    println("# Vmax = ", vmax, ", base Delta V = ", base_dv,
            ", base C = ", base_C, ", base rows = ", base_rows)

    all_samples = []
    all_rows = []
    arrays = []
    for level in levels
        ep, samples, rows = build_series(;
            q0, rho_match, vmax, dv=level.dv, C=level.C,
            max_rows=level.max_rows, Umax, amplitude,
        )
        push!(all_samples, samples)
        push!(all_rows, rows)
        push!(arrays, sample_arrays(samples, ep))
    end

    print_level_table(levels, arrays, all_rows, all_samples)
    print_pair_table(arrays)
end

main()
