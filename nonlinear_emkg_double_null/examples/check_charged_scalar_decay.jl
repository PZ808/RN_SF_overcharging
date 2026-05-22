using NonlinearEMKGDoubleNull

# Gelles/Pretorius arXiv:2503.04881, Eq. (72):
# on extremal RN, the gauge-invariant horizon amplitude P=|r phi|
# decays as V_EF^-s.  For eQ0=0.6, s=1/2, so P ~ V_EF^-1/2.

ep = EvolutionParams(
    rn = RNParams(1.0, 1.0),
    scalar_charge = 0.6,
    amplitude = 1.0e-4,
    omega = 0.0,
    center = 20.0,
    width = 4.0,
)

v0 = compact_v_from_ef_v(0.0, ep.rn)
v1 = compact_v_from_ef_v(900.0, ep.rn)
grid = compact_mrt_grid(ep.rn; nu=90, nv=360, u0=-1.0, v0=v0, u1=-1.0e-4, v1=v1)
state = initialize_state(grid, ep)
evolve!(state, grid, ep)

vef = [NonlinearEMKGDoubleNull.ef_v_from_mrt(v, ep.rn) for v in grid.v]
P = sqrt.(state.xi[end, :].^2 .+ state.pi[end, :].^2)
finite = isfinite.(vef) .& isfinite.(P) .& (vef .> 0) .& (P .> 0)
vef = vef[finite]
P = P[finite]

s = conformal_weight_s(ep.scalar_charge * ep.rn.Q0)
target = -s

println("eQ0 = ", ep.scalar_charge * ep.rn.Q0)
println("s = ", s)
println("target P slope = ", target)
for cut in (20.0, 40.0, 80.0, 120.0)
    idx = vef .>= cut
    count(idx) > 5 || continue
    slope, _, nfit = fit_power_law(vef[idx], P[idx])
    println("fit V_EF >= ", cut, ": samples = ", nfit,
            ", slope = ", slope,
            ", max/min = ", maximum(P[idx]) / minimum(P[idx]))
end

