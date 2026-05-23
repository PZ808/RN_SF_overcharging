using NonlinearEMKGDoubleNull

function argument_or_default(index, default)
    return length(ARGS) >= index ? parse(Float64, ARGS[index]) : default
end

function nearest_sample(u, values, target)
    _, i = findmin(abs.(u .- target))
    return u[i], values[i]
end

# MRT Fig. 7 uses the outgoing, degenerate-initial-apparent-horizon family.
# Scalar charge e=0, while the background electromagnetic charge remains Q=1.
epsilon = argument_or_default(1, 0.02)
geometry_threshold = argument_or_default(2, 0.02)
scalar_threshold = argument_or_default(3, 0.02)
target_v = argument_or_default(4, 150.0)
dv = argument_or_default(5, 0.2)

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.0,
    amplitude = epsilon,
    omega = 0.0,
    center = 15.0,
    width = 3.0,
)

f0 = mrt2013_degenerate_horizon_f0(ep)
initial_mass = mrt2013_initial_bondi_mass(f0)

vmax = target_v + dv
nv = round(Int, vmax / dv) + 1
grid = mrt2013_grid(; nu=110, nv, U0=-5.1, V0=0.0, U1=0.2, V1=dv * (nv - 1))
seed = NLState(grid)
initialize_mrt2013_outgoing_wave!(seed, grid, ep; f0)

splitting = PointSplittingConfig(
    band_width = 1.0,
    max_relative_r = geometry_threshold,
    max_relative_f = geometry_threshold,
    max_dphi = scalar_threshold * epsilon,
    max_dphiu = scalar_threshold * epsilon,
)
chopping = HorizonChoppingConfig(
    start_v = 1.0,
    band_width = 1.0,
    interior_buffer_cells = 1,
)
state = evolve_adaptive(
    slice_from_rectangular(seed, grid, 1),
    west_boundary_from_rectangular(seed, grid),
    ep;
    iterations=20,
    subtract_rn_background=true,
    point_splitting=splitting,
    horizon_chopping=chopping,
)

u, sampled_v, mass = bondi_mass_profile(state; target_v, rn_background=ep.rn)
expansion_u, _, rv = outgoing_expansion_profile(state; target_v)
scaled_mass = (mass .- 1) ./ epsilon^2
scaled_initial_mass = (initial_mass - 1) / epsilon^2
uah = apparent_horizon_location(expansion_u, rv)

if isnothing(uah)
    final_mass = last(mass)
    exterior_mass = mass
else
    exterior_indices = findall(value -> value <= uah, u)
    last_exterior = last(exterior_indices)
    next_index = min(last_exterior + 1, lastindex(u))
    if next_index == last_exterior
        final_mass = mass[last_exterior]
    else
        fraction = (uah - u[last_exterior]) / (u[next_index] - u[last_exterior])
        final_mass = mass[last_exterior] +
                     fraction * (mass[next_index] - mass[last_exterior])
    end
    exterior_mass = mass[firstindex(mass):last_exterior]
end

retained_fraction = (final_mass - 1) / (initial_mass - 1)
max_upward_step = maximum(diff(exterior_mass))
final_mass_label = isnothing(uah) ? "right-boundary" : "horizon-interpolated"

println("epsilon = ", epsilon)
println("f0 = ", f0)
println("(Mi - 1) / epsilon^2 = ", scaled_initial_mass)
println("MRT small-epsilon target for initial coefficient: 0.789")
println("Bondi profile sampled at V = ", sampled_v)
println("Delta V = ", dv)
println("mass extraction uses exact RN stencil-defect subtraction")
println("final U points = ", length(last(state.slices).u))
println("apparent-horizon U at sampled V = ", uah)
println("minimum and rightmost r_V = ", minimum(rv), ", ", last(rv))
println("minimum exterior sampled mass = ", minimum(exterior_mass))
println("maximum exterior upward mass step in U = ", max_upward_step)
println(final_mass_label, " (MB - 1) / epsilon^2 = ", (final_mass - 1) / epsilon^2)
println("(Mf - 1) / (Mi - 1) approximation = ", retained_fraction)
println("MRT Fig. 7 target for retained fraction: 0.0375")

println("sampled epsilon^-2 (MB - 1):")
for target in (-5.0, -4.0, -3.0, -2.0, -1.0, -0.05)
    ui, value = nearest_sample(u, scaled_mass, target)
    println("  U = ", ui, ": ", value)
end
