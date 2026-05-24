using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

epsilon = argument_or_default(1, 0.02)
initial_nu = round(Int, argument_or_default(2, 531.0))
target_v = argument_or_default(3, 10.0)
dv = argument_or_default(4, 0.02)

reference = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = epsilon,
)
f0 = mrt2013_degenerate_horizon_f0(reference)
initial_mass = mrt2013_initial_bondi_mass(f0)

# The scalar is turned off while retaining f0 > 2. MRT's electrovac argument
# then identifies the entire evolved solution as RN with mass initial_mass.
electrovac = EvolutionParams(
    rn = reference.rn,
    scalar_charge = 0.0,
    amplitude = 0.0,
)
nv = round(Int, target_v / dv) + 2
grid = mrt2013_grid(; nu=initial_nu, nv, U0=-5.1, V0=0.0, U1=-0.01,
                    V1=dv * (nv - 1))

println("epsilon used to choose f0 = ", epsilon)
println("target electrovac mass Mi = ", initial_mass)
println("Delta V = ", dv, ", initial U points = ", initial_nu)
println("sampled target V = ", target_v)

for evolution_correction in (false, true)
    state = NLState(grid)
    initialize_mrt2013_outgoing_wave!(state, grid, electrovac; f0)
    evolve_nonlinear!(state, grid, electrovac; iterations=20,
                      subtract_rn_background=evolution_correction)
    adaptive = adaptive_state_from_rectangular(state, grid)
    for diagnostic_correction in (false, true)
        background = diagnostic_correction ? electrovac.rn : nothing
        u, sampled_v, mass = bondi_mass_profile(adaptive; target_v, rn_background=background)
        mass_error = mass .- initial_mass
        derivative = diff(mass) ./ diff(u)
        println("  evolution correction = ", evolution_correction,
                ", diagnostic correction = ", diagnostic_correction,
                ", sampled V = ", sampled_v,
                ", max |varpi - Mi| = ", maximum(abs, mass_error),
                ", relative to Mi - 1 = ", maximum(abs, mass_error) / (initial_mass - 1),
                ", min(varpi - Mi) = ", minimum(mass_error),
                ", max varpi_U = ", maximum(derivative))
    end
end
