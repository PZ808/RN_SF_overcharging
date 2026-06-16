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

function integer_list_argument(index, default)
    return [parse(Int, value) for value in split(argument(index, default), ",")]
end

function step_control_argument(index)
    value = argument(index, "outer")
    value in ("outer", "max-row", "geometric", "throat", "eta", "local") ||
        throw(ArgumentError("step control must be outer, max-row, geometric, throat, eta, or local"))
    return value == "outer" ? :outer :
           value == "max-row" ? :max_row :
           value == "geometric" ? :geometric :
           value == "throat" ? :throat :
           value == "eta" ? :eta :
           :local
end

function finite_row(row::NLRow)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(zero(eltype(row.r))), row.r)
end

function max_abs(values)
    isempty(values) && return NaN
    return maximum(abs, values)
end

function row_qv_residual(row::NLRow, ep::EvolutionParams)
    rv = N.coordinate_derivative(row.r, row.v)
    qv = N.coordinate_derivative(row.Q, row.v)
    psi_re_v = N.coordinate_derivative(row.phi_re, row.v)
    psi_im_v = N.coordinate_derivative(row.phi_im, row.v)
    residual = similar(row.r)
    for j in eachindex(row.v)
        source = stress_energy_reduced_scalar(
            row.r[j], exp(row.logf[j]), row.Q[j], zero(row.r[j]), rv[j],
            row.phi_re[j], row.phi_im[j],
            zero(row.phi_re[j]), psi_re_v[j],
            zero(row.phi_im[j]), psi_im_v[j],
            row.Au[j], row.Av[j], ep.scalar_charge,
        )
        residual[j] = qv[j] + row.r[j]^2 * source.Jv / 8
    end
    return max_abs(residual)
end

function row_difference(row::NLRow, reference::NLRow)
    row.v == reference.v || throw(ArgumentError("rows must share the same V grid"))
    psi_diff = [hypot(row.phi_re[j] - reference.phi_re[j],
                      row.phi_im[j] - reference.phi_im[j])
                for j in eachindex(row.v)]
    return (
        max_abs_r=max_abs(row.r .- reference.r),
        max_abs_logf=max_abs(row.logf .- reference.logf),
        max_abs_psi=max_abs(psi_diff),
        max_abs_Au=max_abs(row.Au .- reference.Au),
        max_abs_Av=max_abs(row.Av .- reference.Av),
        max_abs_Q=max_abs(row.Q .- reference.Q),
    )
end

function row_summary(row::NLRow, ep::EvolutionParams)
    minimum_expansion = row_expansion_minimum(row)
    return (
        finite=finite_row(row),
        min_r=minimum(row.r),
        max_logf=maximum(row.logf),
        min_rv=minimum_expansion.rv,
        min_rv_V=minimum_expansion.v,
        max_qv_residual=row_qv_residual(row, ep),
    )
end

function advance_to_target(base::NLRow, target_u, ep::EvolutionParams;
                           pieces::Int, iterations::Int, hyperbolic_charge::Bool,
                           U0, V0, M0)
    pieces >= 1 || throw(ArgumentError("pieces must be positive"))
    target_u > base.u || throw(ArgumentError("target U must exceed base U"))
    current = base
    for k in 1:pieces
        next_u = base.u + (target_u - base.u) * k / pieces
        next_u > current.u || return current
        south = gp2026_na_boundary_point(next_u, ep; U0, V0, M0)
        current = advance_u_row(current, south, ep;
                                iterations, reduced_scalar=true, hyperbolic_charge)
        finite_row(current) || return current
    end
    return current
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Bool && return string(value)
    value isa Symbol && return string(value)
    value isa Integer && return string(value)
    value isa Real && return isfinite(value) ? string(Float64(value)) : string(value)
    return string(value)
end

function print_table(rows)
    names = propertynames(first(rows))
    println(join(string.(names), '\t'))
    for row in rows
        println(join((format_value(getproperty(row, name)) for name in names), '\t'))
    end
end

function run_stiffness(::Type{T}) where {T<:Real}
    q0 = real_argument(1, "1.0033218", T)
    vmax = real_argument(2, "400.0", T)
    dv = real_argument(3, "0.08", T)
    amplitude = real_argument(4, "0.01", T)
    C = real_argument(5, "0.6", T)
    base_rows = integer_argument(6, 40)
    du_factor = real_argument(7, "1.0", T)
    max_delta_rho = real_argument(8, "0.25", T)
    step_control = step_control_argument(9)
    pieces_list = integer_list_argument(10, "1,2,4,8,16")
    iterations_list = integer_list_argument(11, "1,2,4,8,16")
    reference_iterations = integer_argument(12, 20)
    precision_bits = integer_argument(13, 0)
    max_delta_eta = real_argument(14, "0.025", T)

    U0 = parse(T, "-1.0")
    Umax = parse(T, "1.6")
    ep = EvolutionParams(
        rn=RNParams(one(T), q0),
        scalar_charge=parse(T, "0.6") / q0,
        amplitude=amplitude,
        omega=one(T),
    )

    nv = Int(round(vmax / dv)) + 1
    grid = gp2026_grid(; nu=2, nv, U0, V0=zero(T),
                       U1=U0 + parse(T, "0.01"), V1=vmax)
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(seed, grid, ep)
    initial = row_from_rectangular(seed, grid, 1)
    evolved = evolve_gp2026_u_adaptive(initial, ep; Umax, C,
                                       iterations=reference_iterations,
                                       max_rows=base_rows,
                                       hyperbolic_charge=true,
                                       step_control,
                                       max_delta_rho,
                                       max_delta_eta)
    last_valid = findlast(finite_row, evolved.rows)
    isnothing(last_valid) && error("initial GP row is invalid")
    rows = evolved.rows[1:last_valid]
    base = last(rows)
    du_info = gp2026_row_step_du(rows, C, step_control;
                                 max_delta_rho, max_delta_eta)
    target_du = du_factor * du_info.selected
    isfinite(target_du) && target_du > 0 || error("selected target Delta U is not positive")
    target_u = min(Umax, base.u + target_du)
    target_u > base.u || error("target U did not advance")

    reference_pieces = maximum(pieces_list)
    reference = advance_to_target(base, target_u, ep;
                                  pieces=reference_pieces,
                                  iterations=reference_iterations,
                                  hyperbolic_charge=true, U0, V0=zero(T), M0=ep.rn.M)
    reference_summary = row_summary(reference, ep)
    vtrap = vtrap_diagnostic([base, reference]; missing_status=:stiffness_probe)

    println("# GP2026 row-step stiffness check")
    println("# numeric type = ", T, ", precision bits = ", precision(T),
            ", requested precision arg = ", precision_bits)
    println("# Q0 = ", q0, ", eQ0 = ", ep.scalar_charge * q0,
            ", A0 = ", amplitude)
    println("# Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C,
            ", step_control = ", step_control)
    println("# base rows requested/valid = ", (base_rows, length(rows)),
            ", base U = ", base.u)
    println("# max_delta_rho/max_delta_eta = ",
            (max_delta_rho, max_delta_eta))
    println("# Delta U selected/outer/max-row/geometric/throat/eta = ",
            (du_info.selected, du_info.outer_du, du_info.max_row_du,
             du_info.geometric_du, du_info.throat_du, du_info.eta_du))
    println("# target Delta U = ", target_u - base.u,
            ", target U = ", target_u,
            ", reference pieces = ", reference_pieces,
            ", reference iterations = ", reference_iterations)
    println("# reference finite/min_r/max_logf/min_rv/max_QV_residual = ",
            (reference_summary.finite, reference_summary.min_r,
             reference_summary.max_logf, reference_summary.min_rv,
             reference_summary.max_qv_residual))
    println("# direct Vtrap over base/reference rows = ", vtrap.trap,
            ", status = ", vtrap.status)

    substep_rows = NamedTuple[]
    for pieces in pieces_list
        row = advance_to_target(base, target_u, ep;
                                pieces, iterations=reference_iterations,
                                hyperbolic_charge=true, U0, V0=zero(T), M0=ep.rn.M)
        summary = row_summary(row, ep)
        diff = finite_row(row) && finite_row(reference) ?
               row_difference(row, reference) :
               (; max_abs_r=Inf, max_abs_logf=Inf, max_abs_psi=Inf,
                  max_abs_Au=Inf, max_abs_Av=Inf, max_abs_Q=Inf)
        push!(substep_rows, (
            mode="substeps",
            pieces=pieces,
            iterations=reference_iterations,
            finite=summary.finite,
            min_r=summary.min_r,
            max_logf=summary.max_logf,
            min_rv=summary.min_rv,
            min_rv_V=summary.min_rv_V,
            max_QV_residual=summary.max_qv_residual,
            max_abs_r_vs_ref=diff.max_abs_r,
            max_abs_logf_vs_ref=diff.max_abs_logf,
            max_abs_psi_vs_ref=diff.max_abs_psi,
            max_abs_Q_vs_ref=diff.max_abs_Q,
        ))
    end
    println("# Substep convergence at fixed target U")
    print_table(substep_rows)

    iteration_rows = NamedTuple[]
    for iterations in iterations_list
        row = advance_to_target(base, target_u, ep;
                                pieces=1, iterations,
                                hyperbolic_charge=true, U0, V0=zero(T), M0=ep.rn.M)
        summary = row_summary(row, ep)
        diff = finite_row(row) && finite_row(reference) ?
               row_difference(row, reference) :
               (; max_abs_r=Inf, max_abs_logf=Inf, max_abs_psi=Inf,
                  max_abs_Au=Inf, max_abs_Av=Inf, max_abs_Q=Inf)
        push!(iteration_rows, (
            mode="one-step-iterations",
            pieces=1,
            iterations=iterations,
            finite=summary.finite,
            min_r=summary.min_r,
            max_logf=summary.max_logf,
            min_rv=summary.min_rv,
            min_rv_V=summary.min_rv_V,
            max_QV_residual=summary.max_qv_residual,
            max_abs_r_vs_ref=diff.max_abs_r,
            max_abs_logf_vs_ref=diff.max_abs_logf,
            max_abs_psi_vs_ref=diff.max_abs_psi,
            max_abs_Q_vs_ref=diff.max_abs_Q,
        ))
    end
    println("# Picard-iteration sensitivity for a single target step")
    print_table(iteration_rows)
end

precision_bits = integer_argument(13, 0)
if precision_bits > 0
    setprecision(precision_bits) do
        run_stiffness(BigFloat)
    end
else
    run_stiffness(Float64)
end
