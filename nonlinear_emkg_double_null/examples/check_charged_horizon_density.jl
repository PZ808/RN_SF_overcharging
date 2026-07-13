using NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

real_argument(index, default) = parse(Float64, argument(index, default))
integer_argument(index, default) = parse(Int, argument(index, default))
string_argument(index, default) = argument(index, default)

function quantile_like(x, q)
    xs = sort(collect(x))
    k = clamp(Int(round(1 + q * (length(xs) - 1))), 1, length(xs))
    return xs[k]
end

function safe_fit_power_law(x, y; xmin=nothing, xmax=nothing)
    try
        return fit_power_law(x, y; xmin, xmax)
    catch err
        err isa ArgumentError || rethrow()
        return NaN, NaN, 0
    end
end

# Charged-sector research diagnostic from Gelles/Pretorius arXiv:2503.04881:
# on extremal RN, the horizon charge density scales as
#   rho_Q ~ V_EF^(1 - 2s)
# and for |e Q0| >= 1/2, s=1/2 so rho_Q approaches a nonzero constant.
#
# This is intentionally not a unit test yet. At the current scaffold stage it
# is expected to expose whether the charged horizon evolution is faithful.

q0 = real_argument(1, 1.0)
eQ0 = real_argument(2, 0.6)
amplitude = real_argument(3, 1.0e-4)
vmax = real_argument(4, 900.0)
nu = integer_argument(5, 90)
nv = integer_argument(6, 360)
fit_quantile = real_argument(7, 0.65)
fit_vmin_floor = real_argument(8, 40.0)
u1 = real_argument(9, -1.0e-4)
fit_xmax_quantile = real_argument(10, 0.95)
grid_mode = string_argument(11, "compact")
fit_vmax_override = length(ARGS) >= 12 ? real_argument(12, NaN) : NaN
uef1 = length(ARGS) >= 13 ? real_argument(13, 1000.0) : 1000.0
envelope_mode = length(ARGS) >= 14 ? string_argument(14, "gaussian") : "gaussian"
pulse_width = length(ARGS) >= 15 ? real_argument(15, 4.0) : 4.0
adaptive_passes = length(ARGS) >= 16 ? integer_argument(16, 2) : 2
adaptive_max_points = length(ARGS) >= 17 ? integer_argument(17, max(nv, 1200)) : max(nv, 1200)
extractor_mode = length(ARGS) >= 18 ? string_argument(18, "row") : "row"
extractor_rows = length(ARGS) >= 19 ? integer_argument(19, 4) : 4
patch_refinement_factor = length(ARGS) >= 20 ? integer_argument(20, 4) : 4
patch_vmin_override = length(ARGS) >= 21 ? real_argument(21, NaN) : NaN
patch_vmax_override = length(ARGS) >= 22 ? real_argument(22, NaN) : NaN
patch_u_refinement_factor = length(ARGS) >= 23 ? integer_argument(23, 1) : 1

ep = EvolutionParams(
    rn = RNParams(1.0, q0),
    scalar_charge = eQ0 / q0,
    amplitude = amplitude,
    omega = 0.0,
    center = 20.0,
    width = pulse_width,
)

envelope_function = if envelope_mode == "gaussian"
    gaussian_envelope
elseif envelope_mode == "gp-bump"
    gp2025_bump_envelope
else
    throw(ArgumentError("envelope_mode must be gaussian or gp-bump"))
end

adaptive_summaries = FixedBackgroundVRefinementSummary{Float64}[]
patch_result = nothing

if grid_mode == "adaptive"
    refinement_vmax = isfinite(fit_vmax_override) ? fit_vmax_override : vmax
    refinement_config = FixedBackgroundVRefinementConfig(
        max_passes=adaptive_passes,
        max_points=adaptive_max_points,
        vmin=fit_vmin_floor,
        vmax=refinement_vmax,
        charge_relative_threshold=0.35,
        energy_relative_threshold=0.35,
        component_relative_threshold=0.35,
        scalar_relative_threshold=0.35,
        residual_threshold=Inf,
    )
    state, grid, adaptive_summaries =
        evolve_fixed_background_v_adaptive(
            ep;
            nu,
            nv,
            u0=-1.0,
            u1,
            Vef0=0.0,
            Vef1=vmax,
            envelope=envelope_function,
            config=refinement_config,
        )
elseif grid_mode == "patch"
    refinement_vmax = isfinite(fit_vmax_override) ? fit_vmax_override : vmax
    patch_vmin = isfinite(patch_vmin_override) ? patch_vmin_override : fit_vmin_floor
    patch_vmax = isfinite(patch_vmax_override) ? patch_vmax_override : refinement_vmax
    v0 = compact_v_from_ef_v(0.0, ep.rn)
    v1 = compact_v_from_ef_v(vmax, ep.rn)
    coarse_grid = compact_mrt_grid(ep.rn; nu, nv, u0=-1.0, v0=v0, u1, v1)
    coarse_state = initialize_state(coarse_grid, ep; envelope=envelope_function)
    evolve!(coarse_state, coarse_grid, ep)
    patch_result = evolve_fixed_background_v_patch(
        coarse_state,
        coarse_grid,
        ep;
        vmin=patch_vmin,
        vmax=patch_vmax,
        refinement_factor=patch_refinement_factor,
        u_refinement_factor=patch_u_refinement_factor,
        envelope=envelope_function,
    )
    state = patch_result.state
    grid = patch_result.grid
else
    grid = if grid_mode == "compact"
    v0 = compact_v_from_ef_v(0.0, ep.rn)
    v1 = compact_v_from_ef_v(vmax, ep.rn)
    compact_mrt_grid(ep.rn; nu, nv, u0=-1.0, v0=v0, u1, v1)
elseif grid_mode == "ef-uniform"
    ef_v_mrt_grid(ep.rn; nu, nv, u0=-1.0, u1, Vef0=0.0, Vef1=vmax)
elseif grid_mode == "ef-uv"
    ef_uv_mrt_grid(ep.rn; nu, nv, Uef0=0.0, Uef1=uef1, Vef0=0.0, Vef1=vmax)
else
    throw(ArgumentError("grid_mode must be compact, ef-uniform, ef-uv, adaptive, or patch"))
end
    state = initialize_state(grid, ep; envelope=envelope_function)
    evolve!(state, grid, ep)
end
ru, rv = maxwell_residuals(state, grid, ep)

vef, rho = if extractor_mode == "row"
    horizon_charge_density_series(state, grid, ep)
elseif extractor_mode == "extrapolated"
    horizon_charge_density_extrapolated_series(state, grid, ep; rows=extractor_rows)
else
    throw(ArgumentError("extractor_mode must be row or extrapolated"))
end
vef_energy, rho_energy = if extractor_mode == "row"
    horizon_energy_density_series(state, grid, ep)
else
    horizon_energy_density_extrapolated_series(state, grid, ep; rows=extractor_rows)
end
vef_energy_direct, rho_energy_direct =
    horizon_energy_density_direct_series(state, grid, ep)
vef_components, e_qr, e_pr, e_qv, e_pv = if extractor_mode == "row"
    horizon_energy_density_divided_components(state, grid, ep)
else
    horizon_energy_density_extrapolated_components(state, grid, ep; rows=extractor_rows)
end
finite_charge = isfinite.(vef) .& isfinite.(rho) .& (vef .> 0) .& (abs.(rho) .> 0)
finite_energy = isfinite.(vef_energy) .& isfinite.(rho_energy) .&
                (vef_energy .> 0) .& (rho_energy .> 0)
finite_energy_direct = isfinite.(vef_energy_direct) .& isfinite.(rho_energy_direct) .&
                       (vef_energy_direct .> 0) .& (rho_energy_direct .> 0)
vef_charge = vef[finite_charge]
rho_charge = rho[finite_charge]
rho_abs = abs.(rho_charge)
vef_energy = vef_energy[finite_energy]
rho_energy = rho_energy[finite_energy]
vef_energy_direct = vef_energy_direct[finite_energy_direct]
rho_energy_direct = rho_energy_direct[finite_energy_direct]

vmin = max(fit_vmin_floor, quantile_like(vef_charge, fit_quantile))
vmax_fit = isfinite(fit_vmax_override) ? fit_vmax_override :
           quantile_like(vef_charge, fit_xmax_quantile)
slope, intercept, nfit = safe_fit_power_law(vef_charge, rho_abs; xmin=vmin, xmax=vmax_fit)
s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)
target = 1 - 2s
energy_target = isapprox(abs(ep.rn.Q0), ep.rn.M) ? 2 - 2s : -4s
energy_slope, energy_intercept, energy_nfit =
    safe_fit_power_law(vef_energy, rho_energy; xmin=vmin, xmax=vmax_fit)
energy_direct_slope, energy_direct_intercept, energy_direct_nfit =
    safe_fit_power_law(vef_energy_direct, rho_energy_direct; xmin=vmin, xmax=vmax_fit)
late_mask = (vef_charge .>= vmin) .& (vef_charge .<= vmax_fit)
late = rho_abs[late_mask]
plateau_ratio = maximum(late) / minimum(late)
late_signed = rho_charge[late_mask]
dominant_sign = count(>=(0), late_signed) >= count(<(0), late_signed) ? 1.0 : -1.0
sign_mask = late_mask .& (dominant_sign .* rho_charge .> 0)
signed_slope = count(sign_mask) >= 2 ?
               safe_fit_power_law(vef_charge[sign_mask],
                                  dominant_sign .* rho_charge[sign_mask])[1] :
               NaN
signed_late = dominant_sign .* rho_charge[sign_mask]
signed_ratio = isempty(signed_late) ? NaN : maximum(signed_late) / minimum(signed_late)

println("# q0 = ", q0, ", eQ0 = ", eQ0, ", amplitude = ", amplitude,
        ", Vmax = ", vmax, ", nu = ", nu, ", nv = ", nv,
        ", fit_quantile = ", fit_quantile, ", fit_vmin_floor = ",
        fit_vmin_floor, ", u1 = ", u1, ", fit_xmax_quantile = ",
        fit_xmax_quantile, ", grid_mode = ", grid_mode,
        ", fit_vmax_override = ", fit_vmax_override, ", Uef1 = ", uef1,
        ", envelope_mode = ", envelope_mode, ", pulse_width = ", pulse_width,
        ", adaptive_passes = ", adaptive_passes,
        ", adaptive_max_points = ", adaptive_max_points,
        ", extractor_mode = ", extractor_mode,
        ", extractor_rows = ", extractor_rows,
        ", patch_refinement_factor = ", patch_refinement_factor,
        ", patch_vmin_override = ", patch_vmin_override,
        ", patch_vmax_override = ", patch_vmax_override,
        ", patch_u_refinement_factor = ", patch_u_refinement_factor)
if !isempty(adaptive_summaries)
    println("adaptive summaries:")
    for summary in adaptive_summaries
        println("  pass=", summary.pass,
                " nv=", summary.nv,
                " flagged=", summary.flagged_intervals,
                " max_indicator=", summary.max_indicator,
                " max_charge=", summary.max_charge_indicator,
                " max_energy=", summary.max_energy_indicator,
                " max_component=", summary.max_component_indicator,
                " max_scalar=", summary.max_scalar_indicator,
                " max_residual=", summary.max_residual_indicator)
    end
end
if patch_result !== nothing
    println("patch summary:")
    println("  coarse_start = ", patch_result.coarse_start,
            ", coarse_stop = ", patch_result.coarse_stop,
            ", patch_nv = ", length(patch_result.grid.v),
            ", patch_nu = ", length(patch_result.grid.u),
            ", v_refinement_factor = ", patch_result.v_refinement_factor,
            ", u_refinement_factor = ", patch_result.u_refinement_factor)
end
println("charge samples = ", length(vef_charge))
println("charge fit samples = ", nfit)
println("energy samples = ", length(vef_energy))
println("energy fit samples = ", energy_nfit)
println("direct-energy samples = ", length(vef_energy_direct))
println("direct-energy fit samples = ", energy_direct_nfit)
println("fit window V_EF = ", (vmin, vmax_fit))
println("eQ0 = ", ep.scalar_charge * ep.rn.Q0)
println("s = ", s)
println("late-time |rho_Q| slope = ", slope)
println("late-time sign-coherent rho_Q slope = ", signed_slope)
println("linear target slope 1 - 2s = ", target)
println("late-time rho_E slope = ", energy_slope)
println("late-time rho_E direct-check slope = ", energy_direct_slope)
println("linear target energy slope = ", energy_target)
println("late-window plateau max/min = ", plateau_ratio)
println("sign-coherent late-window max/min = ", signed_ratio)
println("max |Maxwell-U residual| = ", maximum(abs, ru))
println("max |Maxwell-V residual| = ", maximum(abs, rv))
println("late samples:")
for target_v in (80.0, 120.0, 150.0, 200.0, 300.0, 400.0, 600.0, 900.0, 1200.0)
    target_v <= maximum(vef_charge) || continue
    _, jq = findmin(abs.(vef_charge .- target_v))
    _, je = findmin(abs.(vef_energy .- target_v))
    _, jd = findmin(abs.(vef_energy_direct .- target_v))
    _, jc = findmin(abs.(vef_components .- target_v))
    println("  V=", vef_charge[jq],
            " rho_Q=", rho_charge[jq],
            " rho_E=", rho_energy[je],
            " rho_E_direct_check=", rho_energy_direct[jd],
            " components=(Qr=", e_qr[jc],
            ", Pr=", e_pr[jc],
            ", Qv=", e_qv[jc],
            ", Pv=", e_pv[jc], ")")
end
