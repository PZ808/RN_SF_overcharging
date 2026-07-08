using NonlinearEMKGDoubleNull

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

function main()
    q0 = real_argument(1, 1.0)
    vmax = real_argument(2, 400.0)
    dv = real_argument(3, 0.08)
    amplitude = real_argument(4, 0.01)
    C = real_argument(5, 0.6)
    Umax = real_argument(6, 1.6)
    max_rows = integer_argument(7, 6000)
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
        "paper",
        Dict("paper" => :paper, "lte" => :lte),
    )
    controller_target_lte = real_argument(20, 0.7)
    controller_safety = real_argument(21, 0.9)
    controller_min_factor = real_argument(22, 0.5)
    controller_max_factor = real_argument(23, 1.5)

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
        du = if root_stepper === :lte
            next_du
        else
            gp2026_row_step_du(rows, C, step_control).selected
        end
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
    valid_rows = rows[1:last_valid]
    final = last(valid_rows)
    next_step = gp2026_row_step_du(valid_rows, C, step_control)
    coordinate_stalled = final.u + next_step.selected == final.u
    missing_status = last_valid < length(rows) ? :invalid_row :
                     final.u == Umax ? :reached_umax :
                     coordinate_stalled ? :precision_stalled :
                     length(rows) >= max_rows ? :max_rows :
                     :stopped
    vtrap = vtrap_diagnostic(valid_rows; missing_status)
    refined_vtrap = refined_vtrap_sample(valid_rows)
    throat = throat_row_diagnostics(final)
    horizon_charge = apparent_horizon_charge_density_series(valid_rows)
    charge_pairs = sort(
        [
            (v=sample.v, density=abs(sample.surface_density),
             flux_density=sample.flux_density_v,
             radial_density=abs(sample.radial_density_proxy))
            for sample in horizon_charge
            if isfinite(sample.v) &&
               isfinite(sample.surface_density) &&
               abs(sample.surface_density) > 0
        ];
        by=pair -> pair.v,
    )
    charge_v = [pair.v for pair in charge_pairs]
    charge_density_abs = [pair.density for pair in charge_pairs]
    charge_late_vmin = quantile_like(charge_v, 0.65)
    late_charge_density = isnothing(charge_late_vmin) ?
                          Float64[] :
                          [pair.density for pair in charge_pairs
                           if pair.v >= charge_late_vmin]
    charge_slope = if !isnothing(charge_late_vmin) &&
                      length(late_charge_density) >= 2
        finite_v = [pair.v for pair in charge_pairs
                    if pair.v >= charge_late_vmin]
        fit_power_law(finite_v, late_charge_density)[1]
    else
        NaN
    end
    late_radial_density = isnothing(charge_late_vmin) ?
                          Float64[] :
                          [pair.radial_density for pair in charge_pairs
                           if pair.v >= charge_late_vmin &&
                              isfinite(pair.radial_density) &&
                              pair.radial_density > 0]
    radial_slope = if !isnothing(charge_late_vmin) &&
                      length(late_radial_density) >= 2
        finite_v = [pair.v for pair in charge_pairs
                    if pair.v >= charge_late_vmin &&
                       isfinite(pair.radial_density) &&
                       pair.radial_density > 0]
        fit_power_law(finite_v, late_radial_density)[1]
    else
        NaN
    end

    println("# GP2026 persistent Hamade-Stewart AMR")
    println(
        "# Q0=", q0,
        ", eQ0=", ep.scalar_charge * q0,
        ", A0=", amplitude,
        ", Vmax=", vmax,
        ", DeltaV=", dv,
        ", C=", C,
    )
    println(
        "# gauge=", pulse_leg_gauge,
        ", step_control=", step_control,
        ", bo_atol=", bo_atol,
        ", bo_rtol=", bo_rtol,
        ", max_levels=", max_levels,
        ", revision_interval=", revision_interval,
        ", stop_on_trap=", stop_on_trap,
        ", fields=", refinement_fields,
        ", max_siblings=", max_sibling_patches,
        ", merge_gap_points=", merge_gap_points,
        ", reject_finest_lte=", reject_on_finest_lte,
        ", root_stepper=", root_stepper,
    )
    println("stored rows = ", length(rows))
    println("valid rows = ", length(valid_rows))
    println("last U = ", final.u)
    println("termination = ", missing_status)
    println("Vtrap status = ", vtrap.status)
    println("direct Vtrap = ", vtrap.trap)
    println("quadratic Vtrap = ", refined_vtrap)
    if !isnothing(vtrap.trap)
        println("Vtrap invariants = ",
                trapped_surface_invariants(vtrap.trap))
    end
    if !isnothing(refined_vtrap)
        println("quadratic Vtrap invariants = ",
                trapped_surface_invariants(refined_vtrap))
    end
    println("closest Vtrap proxy = ", vtrap.closest)
    println("horizon charge-density samples = ", length(horizon_charge))
    println("horizon charge-density V range = ",
            isempty(charge_v) ? nothing : extrema(charge_v))
    println("horizon charge-density late Vmin = ", charge_late_vmin)
    println("horizon |Q|/(4pi r^2) sorted first/last = ",
            isempty(charge_density_abs) ? nothing :
            (first(charge_density_abs), last(charge_density_abs)))
    println("horizon Q_V/(4pi r^2) sorted first/last = ",
            isempty(charge_pairs) ? nothing :
            (first(charge_pairs).flux_density, last(charge_pairs).flux_density))
    println("horizon |Q|/(4pi r^2) late slope = ", charge_slope)
    println("horizon |Q|/(4pi r^2) late max/min = ",
            isempty(late_charge_density) ? nothing :
            maximum(late_charge_density) / minimum(late_charge_density))
    println("horizon |Q|/(4pi r^2) late relative std = ",
            relative_std(late_charge_density))
    println("horizon derivative-density proxy late slope = ", radial_slope)
    println("horizon derivative-density proxy late max/min = ",
            isempty(late_radial_density) ? nothing :
            maximum(late_radial_density) / minimum(late_radial_density))
    println("horizon derivative-density proxy late relative std = ",
            relative_std(late_radial_density))
    println("throat min(r-|Q|) = ", throat.min_y)
    println("throat max rho = ", throat.max_rho)
    println("hierarchy depth = ", stewart_hierarchy_depth(hierarchy.root))
    println("hierarchy point counts = ", hierarchy_point_counts(hierarchy.root))
    println("hierarchy patch counts = ", hierarchy_patch_counts(hierarchy.root))
    println("hierarchy patches = ",
            stewart_hierarchy_intervals(hierarchy.root))
    println("level steps = ", hierarchy.stats.level_steps)
    println("level revisions = ", hierarchy.stats.revisions)
    println("last normalized LTE by level = ", hierarchy.stats.last_lte)
    println("maximum normalized LTE by level = ", hierarchy.stats.max_lte)
    println("child creations = ", hierarchy.stats.child_creations)
    println("child destructions = ", hierarchy.stats.child_destructions)
    println("injections = ", hierarchy.stats.injections)
    println("suffix reintegrations = ",
            hierarchy.stats.suffix_reintegrations)
    println("precision fallbacks = ", hierarchy.stats.precision_fallbacks)
    println("rejected root steps = ", hierarchy.stats.rejected_root_steps)
    println("rejected finest-LTE steps = ",
            hierarchy.stats.rejected_finest_lte_steps)
    println("accepted root Delta U extrema = ",
            isempty(accepted_du) ? nothing : extrema(accepted_du))
    println("accepted root LTE extrema = ",
            isempty(accepted_lte) ? nothing : extrema(accepted_lte))
    println("next LTE-controller Delta U = ", next_du)
    println("next paper/local Delta U = ", next_step.selected)
end

main()
