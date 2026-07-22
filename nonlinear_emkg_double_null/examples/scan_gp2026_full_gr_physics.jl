using NonlinearEMKGDoubleNull

# Usage:
#   julia --project=nonlinear_emkg_double_null \
#       nonlinear_emkg_double_null/examples/scan_gp2026_full_gr_physics.jl \
#       q_values amplitudes vmax dv C Umax max_rows bo_rtol max_levels \
#       revision_interval gauge step_control bo_atol stop_on_trap fields \
#       max_sibling_patches merge_gap_points reject_on_finest_lte root_stepper
#
# Example:
#   julia --project=nonlinear_emkg_double_null \
#       nonlinear_emkg_double_null/examples/scan_gp2026_full_gr_physics.jl \
#       1.001,1.0033218,1.006 0.01 160 0.12 0.6 1.2 100
#
# The script prints one TSV row per case with trapped-surface, throat,
# horizon charge-density, and Stewart-AMR controller diagnostics.

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))

function boolean_argument(index, default)
    value = lowercase(argument(index, default))
    value in ("true", "t", "1", "yes", "y", "on") && return true
    value in ("false", "f", "0", "no", "n", "off") && return false
    throw(ArgumentError("argument $index must be true or false"))
end

function real_list_argument(index, default)
    return [parse(Float64, value) for value in split(argument(index, default), ",")]
end

function symbol_argument(index, default, choices)
    value = argument(index, default)
    value in keys(choices) ||
        throw(ArgumentError("argument $index must be one of $(collect(keys(choices)))"))
    return choices[value]
end

function finite_row(row)
    return all(isfinite, row.r) && all(isfinite, row.logf) &&
           all(isfinite, row.phi_re) && all(isfinite, row.phi_im) &&
           all(isfinite, row.Au) && all(isfinite, row.Av) &&
           all(isfinite, row.Q) && all(>(0), row.r)
end

function quantile_like(x, q)
    xs = sort(collect(x))
    isempty(xs) && return nothing
    k = clamp(Int(round(1 + q * (length(xs) - 1))), 1, length(xs))
    return xs[k]
end

function relative_std(x)
    isempty(x) && return NaN
    mean = sum(x) / length(x)
    mean == 0 && return Inf
    variance = sum((value - mean)^2 for value in x) / length(x)
    return sqrt(variance) / abs(mean)
end

function finite_extrema_or_nothing(x)
    values = [value for value in x if isfinite(value)]
    isempty(values) && return nothing
    return extrema(values)
end

function hierarchy_point_counts(level, counts=Int[])
    while length(counts) <= level.level
        push!(counts, 0)
    end
    counts[level.level + 1] += length(level.current.v)
    for child in level.children
        hierarchy_point_counts(child, counts)
    end
    return counts
end

function hierarchy_patch_counts(level, counts=Int[])
    while length(counts) <= level.level
        push!(counts, 0)
    end
    counts[level.level + 1] += 1
    for child in level.children
        hierarchy_patch_counts(child, counts)
    end
    return counts
end

function charge_density_summary(rows)
    samples = apparent_horizon_charge_density_series(rows)
    pairs = sort(
        [
            (
                v=sample.v,
                density=abs(sample.surface_density),
                signed_density=sample.surface_density,
                flux_density=sample.flux_density_v,
                radial_density=abs(sample.radial_density_proxy),
            )
            for sample in samples
            if isfinite(sample.v) &&
               isfinite(sample.surface_density) &&
               abs(sample.surface_density) > 0
        ];
        by=pair -> pair.v,
    )
    v = [pair.v for pair in pairs]
    density = [pair.density for pair in pairs]
    late_vmin = quantile_like(v, 0.65)
    late_pairs = isnothing(late_vmin) ? typeof(pairs)() :
                 [pair for pair in pairs if pair.v >= late_vmin]
    late_v = [pair.v for pair in late_pairs]
    late_density = [pair.density for pair in late_pairs]
    late_slope = if length(late_density) >= 2
        fit_power_law(late_v, late_density)[1]
    else
        NaN
    end
    radial_pairs = [
        pair for pair in late_pairs
        if isfinite(pair.radial_density) && pair.radial_density > 0
    ]
    radial_slope = if length(radial_pairs) >= 2
        fit_power_law([pair.v for pair in radial_pairs],
                      [pair.radial_density for pair in radial_pairs])[1]
    else
        NaN
    end

    return (;
        samples=length(samples),
        finite_samples=length(pairs),
        v_min=isempty(v) ? nothing : minimum(v),
        v_max=isempty(v) ? nothing : maximum(v),
        late_vmin,
        late_slope,
        late_maxmin=isempty(late_density) ? nothing :
                    maximum(late_density) / minimum(late_density),
        late_relstd=relative_std(late_density),
        first_density=isempty(density) ? nothing : first(density),
        last_density=isempty(density) ? nothing : last(density),
        first_flux=isempty(pairs) ? nothing : first(pairs).flux_density,
        last_flux=isempty(pairs) ? nothing : last(pairs).flux_density,
        radial_late_slope=radial_slope,
    )
end

function run_case(;
    q0,
    amplitude,
    vmax,
    dv,
    C,
    Umax,
    max_rows,
    bo_rtol,
    bo_atol,
    max_levels,
    revision_interval,
    pulse_leg_gauge,
    step_control,
    stop_on_trap,
    refinement_fields,
    max_sibling_patches,
    merge_gap_points,
    reject_on_finest_lte,
    root_stepper,
    controller_target_lte,
    controller_safety,
    controller_min_factor,
    controller_max_factor,
)
    U0 = -1.0
    ep = EvolutionParams(
        rn=RNParams(1.0, q0),
        scalar_charge=0.6 / q0,
        amplitude=amplitude,
        omega=1.0,
    )
    grid = gp2026_grid(
        ; nu=2,
        nv=Int(round(vmax / dv)) + 1,
        U0,
        V0=0.0,
        U1=U0 + 0.01,
        V1=vmax,
    )
    seed = NLState(grid)
    initialize_gp2026_single_pulse!(
        seed, grid, ep;
        pulse_leg_gauge=pulse_leg_gauge,
    )
    initial = row_from_rectangular(seed, grid, 1)
    config = StewartAMRConfig(
        refinement_factor=4,
        revision_interval=revision_interval,
        max_levels=max_levels,
        atol=bo_atol,
        rtol=bo_rtol,
        order=2,
        buffer_points=4,
        merge_gap_points=merge_gap_points,
        max_sibling_patches=max_sibling_patches,
        fields=refinement_fields,
        reject_on_finest_lte=reject_on_finest_lte,
        reintegrate=true,
    )
    hierarchy = initialize_stewart_hierarchy(initial; config)
    rows = NLRow[initial]
    controller = StewartRootStepController(
        target_lte=controller_target_lte,
        safety=controller_safety,
        min_factor=controller_min_factor,
        max_factor=controller_max_factor,
        min_du=0.0,
        max_du=Inf,
    )
    next_du = gp2026_row_step_du(rows, C, step_control).selected
    accepted_lte = Float64[]
    accepted_du = Float64[]
    while last(rows).u < Umax && length(rows) < max_rows
        du = root_stepper === :lte ? next_du :
             gp2026_row_step_du(rows, C, step_control).selected
        isfinite(du) && du > 0 || break
        start_u = last(rows).u
        target_u = min(Umax, start_u + du)
        target_u > start_u || break

        result = advance_stewart_hierarchy!(
            hierarchy,
            target_u,
            ep;
            U0,
            V0=0.0,
            M0=1.0,
            pulse_leg_gauge=pulse_leg_gauge,
            iterations=12,
            hyperbolic_charge=true,
            cell_solver=:newton_direct,
        )
        push!(rows, result.row)
        finite_row(result.row) || break

        step_du = result.accepted_target_u - start_u
        step_lte = deepest_recent_lte(hierarchy.stats)
        push!(accepted_du, Float64(step_du))
        push!(accepted_lte, Float64(step_lte))
        if root_stepper === :lte
            next_du = next_stewart_root_du(
                controller,
                step_du,
                step_lte;
                order=config.order,
            )
        end
        if stop_on_trap &&
           !isnothing(row_apparent_horizon_crossing(
               result.row;
               row_index=length(rows),
           ))
            break
        end
    end

    last_valid = findlast(finite_row, rows)
    valid_rows = isnothing(last_valid) ? NLRow[] : rows[1:last_valid]
    isempty(valid_rows) && error("initial row invalid")
    final = last(valid_rows)
    next_step = gp2026_row_step_du(valid_rows, C, step_control)
    coordinate_stalled = final.u + next_step.selected == final.u
    termination = last_valid < length(rows) ? :invalid_row :
                  final.u == Umax ? :reached_umax :
                  coordinate_stalled ? :precision_stalled :
                  length(rows) >= max_rows ? :max_rows :
                  :stopped
    vtrap = vtrap_diagnostic(valid_rows; missing_status=termination)
    refined_vtrap = refined_vtrap_sample(valid_rows)
    trap = isnothing(refined_vtrap) ? vtrap.trap : refined_vtrap
    invariants = isnothing(trap) ? nothing : trapped_surface_invariants(trap)
    closest_invariants = trapped_surface_invariants(vtrap.closest)
    throat = throat_row_diagnostics(final)
    charge = charge_density_summary(valid_rows)
    du_extrema = finite_extrema_or_nothing(accepted_du)
    lte_extrema = finite_extrema_or_nothing(accepted_lte)

    return (;
        q0,
        eQ0=ep.scalar_charge * q0,
        amplitude,
        vmax,
        dv,
        rows=length(rows),
        valid_rows=length(valid_rows),
        last_U=final.u,
        termination,
        vtrap_status=vtrap.status,
        vtrap_V=isnothing(trap) ? nothing : trap.v,
        vtrap_U=isnothing(trap) ? nothing : trap.u,
        vtrap_r=isnothing(trap) ? nothing : trap.r,
        vtrap_Q=isnothing(trap) ? nothing : trap.q,
        vtrap_M=isnothing(invariants) ? nothing : invariants.mass,
        vtrap_q_over_m=isnothing(invariants) ? nothing : invariants.q_over_m,
        vtrap_one_minus_q_over_m=isnothing(invariants) ? nothing :
                                 invariants.one_minus_q_over_m,
        vtrap_r_over_m=isnothing(invariants) ? nothing : invariants.r_over_m,
        closest_V=vtrap.closest.v,
        closest_U=vtrap.closest.u,
        closest_r=vtrap.closest.r,
        closest_Q=vtrap.closest.q,
        closest_rv=vtrap.closest.rv,
        closest_q_over_m=closest_invariants.q_over_m,
        throat_min_y=throat.min_y,
        throat_max_rho=throat.max_rho,
        charge_samples=charge.samples,
        charge_finite_samples=charge.finite_samples,
        charge_V_min=charge.v_min,
        charge_V_max=charge.v_max,
        charge_late_Vmin=charge.late_vmin,
        charge_late_slope=charge.late_slope,
        charge_late_maxmin=charge.late_maxmin,
        charge_late_relstd=charge.late_relstd,
        charge_first=charge.first_density,
        charge_last=charge.last_density,
        charge_flux_first=charge.first_flux,
        charge_flux_last=charge.last_flux,
        charge_radial_late_slope=charge.radial_late_slope,
        hierarchy_depth=stewart_hierarchy_depth(hierarchy.root),
        hierarchy_point_counts=join(hierarchy_point_counts(hierarchy.root), ","),
        hierarchy_patch_counts=join(hierarchy_patch_counts(hierarchy.root), ","),
        level_steps=join(hierarchy.stats.level_steps, ","),
        level_revisions=join(hierarchy.stats.revisions, ","),
        max_lte=join(hierarchy.stats.max_lte, ","),
        last_lte=join(hierarchy.stats.last_lte, ","),
        accepted_du_min=isnothing(du_extrema) ? nothing : first(du_extrema),
        accepted_du_max=isnothing(du_extrema) ? nothing : last(du_extrema),
        accepted_lte_min=isnothing(lte_extrema) ? nothing : first(lte_extrema),
        accepted_lte_max=isnothing(lte_extrema) ? nothing : last(lte_extrema),
        child_creations=hierarchy.stats.child_creations,
        child_destructions=hierarchy.stats.child_destructions,
        injections=hierarchy.stats.injections,
        suffix_reintegrations=hierarchy.stats.suffix_reintegrations,
        rejected_root_steps=hierarchy.stats.rejected_root_steps,
        rejected_finest_lte_steps=hierarchy.stats.rejected_finest_lte_steps,
    )
end

function format_value(value)
    isnothing(value) && return "missing"
    value isa Symbol && return string(value)
    value isa Bool && return string(value)
    value isa AbstractString && return value
    value isa Integer && return string(value)
    value isa Real && return isfinite(value) ? string(Float64(value)) : string(value)
    return replace(string(value), '\t' => ' ')
end

function print_row(columns, row)
    println(join([format_value(getproperty(row, column)) for column in columns], '\t'))
end

function main()
    q_values = real_list_argument(1, "1.001,1.0033218,1.006")
    amplitudes = real_list_argument(2, "0.01")
    vmax = real_argument(3, 160.0)
    dv = real_argument(4, 0.12)
    C = real_argument(5, 0.6)
    Umax = real_argument(6, 1.2)
    max_rows = integer_argument(7, 100)
    bo_rtol = real_argument(8, 1.0e-5)
    max_levels = integer_argument(9, 3)
    revision_interval = integer_argument(10, 4)
    pulse_leg_gauge = symbol_argument(
        11,
        "areal-affine",
        Dict("areal-affine" => :areal_affine, "ef-affine" => :ef_affine),
    )
    step_control = symbol_argument(
        12,
        "local",
        Dict(
            "outer" => :outer,
            "max-row" => :max_row,
            "geometric" => :geometric,
            "throat" => :throat,
            "eta" => :eta,
            "local" => :local,
        ),
    )
    bo_atol = real_argument(13, 1.0e-8)
    stop_on_trap = boolean_argument(14, false)
    refinement_fields = symbol_argument(
        15,
        "all",
        Dict(
            "all" => (:r, :logf, :psi_abs, :Q, :eta),
            "matter" => (:psi_abs, :Q),
            "charge" => (:Q,),
        ),
    )
    max_sibling_patches = integer_argument(16, 8)
    merge_gap_points = integer_argument(17, 2)
    reject_on_finest_lte = boolean_argument(18, true)
    root_stepper = symbol_argument(
        19,
        "lte",
        Dict("paper" => :paper, "lte" => :lte),
    )
    controller_target_lte = real_argument(20, 0.7)
    controller_safety = real_argument(21, 0.9)
    controller_min_factor = real_argument(22, 0.5)
    controller_max_factor = real_argument(23, 1.5)

    columns = (
        :q0, :eQ0, :amplitude, :vmax, :dv, :rows, :valid_rows, :last_U,
        :termination, :vtrap_status, :vtrap_V, :vtrap_U, :vtrap_M,
        :vtrap_q_over_m, :vtrap_one_minus_q_over_m, :vtrap_r_over_m,
        :closest_V, :closest_rv, :closest_q_over_m, :throat_min_y,
        :throat_max_rho, :charge_samples, :charge_finite_samples,
        :charge_V_min, :charge_V_max, :charge_late_Vmin,
        :charge_late_slope, :charge_late_maxmin, :charge_late_relstd,
        :charge_first, :charge_last, :charge_flux_first, :charge_flux_last,
        :charge_radial_late_slope, :hierarchy_depth, :hierarchy_point_counts,
        :hierarchy_patch_counts, :level_steps, :level_revisions, :max_lte,
        :last_lte, :accepted_du_min, :accepted_du_max, :accepted_lte_min,
        :accepted_lte_max, :child_creations, :child_destructions, :injections,
        :suffix_reintegrations, :rejected_root_steps,
        :rejected_finest_lte_steps,
    )
    println("# GP2026 full-GR Stewart AMR physics scan")
    println("# q_values=", join(q_values, ","),
            ", amplitudes=", join(amplitudes, ","),
            ", Vmax=", vmax,
            ", DeltaV=", dv,
            ", C=", C,
            ", Umax=", Umax,
            ", max_rows=", max_rows,
            ", bo_rtol=", bo_rtol,
            ", max_levels=", max_levels,
            ", revision_interval=", revision_interval,
            ", gauge=", pulse_leg_gauge,
            ", step_control=", step_control,
            ", bo_atol=", bo_atol,
            ", stop_on_trap=", stop_on_trap,
            ", fields=", refinement_fields,
            ", max_siblings=", max_sibling_patches,
            ", merge_gap_points=", merge_gap_points,
            ", reject_finest_lte=", reject_on_finest_lte,
            ", root_stepper=", root_stepper)
    println(join(string.(columns), '\t'))
    for amplitude in amplitudes
        for q0 in q_values
            try
                row = run_case(;
                    q0,
                    amplitude,
                    vmax,
                    dv,
                    C,
                    Umax,
                    max_rows,
                    bo_rtol,
                    bo_atol,
                    max_levels,
                    revision_interval,
                    pulse_leg_gauge,
                    step_control,
                    stop_on_trap,
                    refinement_fields,
                    max_sibling_patches,
                    merge_gap_points,
                    reject_on_finest_lte,
                    root_stepper,
                    controller_target_lte,
                    controller_safety,
                    controller_min_factor,
                    controller_max_factor,
                )
                print_row(columns, row)
            catch err
                @warn "scan case failed" q0 amplitude err
            end
        end
    end
end

main()
