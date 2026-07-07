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
    while last(rows).u < Umax && length(rows) < max_rows
        step = gp2026_row_step_du(rows, C, step_control)
        du = step.selected
        isfinite(du) && du > 0 || break
        target_u = min(Umax, last(rows).u + du)
        target_u > last(rows).u || break
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
    println("next controlled Delta U = ", next_step.selected)
end

main()
