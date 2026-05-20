module NonlinearEMKGDoubleNull

include("Coordinates.jl")
include("Fields.jl")
include("EinsteinSector.jl")
include("StressEnergy.jl")
include("NonlinearEMKG.jl")
include("InitialData.jl")
include("Evolution.jl")
include("Diagnostics.jl")

export RNParams, EvolutionParams, Grid, State
export compact_mrt_grid, initialize_state, evolve!, maxwell_residuals
export MetricState, einstein_backreaction_rhs!
export NLState, initialize_nonlinear_state, evolve_nonlinear!, mrt2013_background_f
export mrt2013_grid, initialize_mrt2013_uncharged_ingoing!
export charged_scalar_rhs, metric_rhs, maxwell_rhs
export StressEnergyComponents, stress_energy, covariant_scalar_derivatives, current_components
export outgoing_constraint_source, ingoing_constraint_source
export horizon_phi_series, fit_power_law
export horizon_charge_density_series, conformal_weight_s
export rstar, areal_radius, metric_F, metric_f, metric_ftilde, radius_from_rstar, compact_v_from_ef_v

end
