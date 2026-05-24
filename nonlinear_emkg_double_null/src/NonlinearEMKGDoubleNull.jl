module NonlinearEMKGDoubleNull

include("Coordinates.jl")
include("Fields.jl")
include("EinsteinSector.jl")
include("StressEnergy.jl")
include("NonlinearEMKG.jl")
include("AdaptiveMRT.jl")
include("InitialData.jl")
include("Evolution.jl")
include("Diagnostics.jl")

export RNParams, EvolutionParams, Grid, State
export compact_mrt_grid, initialize_state, evolve!, maxwell_residuals
export MetricState, einstein_backreaction_rhs!
export NLState, initialize_nonlinear_state, evolve_nonlinear!, mrt2013_background_f
export mrt2013_grid, mrt2013_areal_radius, mrt2013_metric_f, initialize_mrt2013_uncharged_ingoing!
export initialize_mrt2013_outgoing_wave!, mrt2013_initial_bondi_mass, mrt2013_degenerate_horizon_f0
export mrt2013_initial_rv_profile
export charged_scalar_rhs, metric_rhs, maxwell_rhs
export NLSlice, NLPoint, AdaptiveNLState, PointSplittingConfig, HorizonChoppingConfig
export HorizonRefinementConfig, slice_point
export slice_from_rectangular, adaptive_state_from_rectangular, west_boundary_from_rectangular
export interpolate_slice, refine_u_grid, refine_slice, truncate_slice
export spacing_refinement_flags, variation_refinement_flags, point_splitting_flags
export adaptive_outgoing_expansion, chop_inside_apparent_horizon
export horizon_refinement_flags, refine_near_apparent_horizon
export advance_adaptive_slice, evolve_adaptive
export StressEnergyComponents, stress_energy, covariant_scalar_derivatives, current_components
export outgoing_constraint_source, ingoing_constraint_source
export horizon_phi_series, fit_power_law
export renormalized_hawking_mass, renormalized_hawking_mass_profile, bondi_mass_profile
export uncharged_mass_flux_u_profile
export outgoing_expansion_profile, apparent_horizon_location
export horizon_charge_density_series, conformal_weight_s
export rstar, areal_radius, metric_F, metric_f, metric_ftilde, radius_from_rstar, compact_v_from_ef_v

end
