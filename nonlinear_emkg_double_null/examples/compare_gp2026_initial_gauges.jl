using NonlinearEMKGDoubleNull

const N = NonlinearEMKGDoubleNull

function argument(index, default)
    return length(ARGS) >= index ? ARGS[index] : string(default)
end

function pulse_leg_geometry(mode::Symbol, v, r0, reference)
    if mode === :ef_affine
        r = radius_from_rstar(rstar(r0, reference) + v / 2, reference)
        return r, metric_F(r, reference) / 2
    elseif mode === :areal_affine
        return r0 + v / 2, one(v) / 2
    end
    throw(ArgumentError("unknown pulse-leg gauge: $mode"))
end

function initial_row(mode::Symbol, ep::EvolutionParams, v)
    q0 = ep.rn.Q0
    r0 = ep.rn.M + one(eltype(v)) / 2
    reference = RNParams(ep.rn.M, ep.rn.M)
    scale = sqrt(32pi)

    r = similar(v)
    rv = similar(v)
    psi_re = similar(v)
    psi_im = similar(v)
    q = fill(convert(eltype(v), q0), length(v))
    logf = similar(v)
    Au = zeros(eltype(v), length(v))
    Av = zeros(eltype(v), length(v))

    for j in eachindex(v)
        r[j], rv[j] = pulse_leg_geometry(mode, v[j], r0, reference)
        amplitude = gp2026_single_pulse_envelope(v[j]; amplitude=ep.amplitude)
        phase = ep.omega * v[j]
        psi_re[j] = scale * amplitude * cos(phase)
        psi_im[j] = -scale * amplitude * sin(phase)
    end

    ru0 = -one(eltype(v)) / 2
    denominator = q0^2 + r0 * (r0 - 2ep.rn.M)
    fcorner = -4r0^2 * ru0 * rv[1] / denominator
    logf[1] = log(fcorner)
    logf_integral = zero(eltype(v))

    for j in firstindex(v)+1:lastindex(v)
        dv = v[j] - v[j - 1]
        vmid = (v[j] + v[j - 1]) / 2
        rmid, rvmid = pulse_leg_geometry(mode, vmid, r0, reference)
        amplitude = gp2026_single_pulse_envelope(vmid; amplitude=ep.amplitude)
        derivative = N.gp2026_single_pulse_envelope_derivative(
            vmid; amplitude=ep.amplitude,
        )
        phase = ep.omega * vmid
        z_re = amplitude * cos(phase)
        z_im = -amplitude * sin(phase)
        dz_re = derivative * cos(phase) - ep.omega * amplitude * sin(phase)
        dz_im = -derivative * sin(phase) - ep.omega * amplitude * cos(phase)
        phi_re = scale * z_re / rmid
        phi_im = scale * z_im / rmid
        phiv_re = scale * (dz_re / rmid - z_re * rvmid / rmid^2)
        phiv_im = scale * (dz_im / rmid - z_im * rvmid / rmid^2)
        _, Jv = current_components(
            phi_re, phi_im, zero(phi_re), zero(phi_im),
            phiv_re, phiv_im, ep.scalar_charge,
        )
        q[j] = q[j - 1] - dv * rmid^2 * Jv / 8
        logf_integral += dv * rmid * (phiv_re^2 + phiv_im^2) / (4rvmid)
        logf[j] = log(fcorner * rv[j] / rv[1]) + logf_integral
    end

    for j in firstindex(v)+1:lastindex(v)
        dv = v[j] - v[j - 1]
        left = q[j - 1] * exp(logf[j - 1]) / (4r[j - 1]^2)
        right = q[j] * exp(logf[j]) / (4r[j]^2)
        Au[j] = Au[j - 1] + dv * (left + right) / 2
    end

    return NLRow(-one(eltype(v)), v, r, logf, psi_re, psi_im, Au, Av, q),
           fcorner, rv
end

function initial_mass_profile(row::NLRow, rv)
    ru = similar(row.v)
    mass = similar(row.v)
    ru[1] = -one(eltype(row.v)) / 2
    mass[1] = renormalized_hawking_mass(
        row.r[1], exp(row.logf[1]), ru[1], rv[1], row.Q[1],
    )
    for j in firstindex(row.v)+1:lastindex(row.v)
        dv = row.v[j] - row.v[j - 1]
        rmid = (row.r[j] + row.r[j - 1]) / 2
        rvmid = (row.r[j] - row.r[j - 1]) / dv
        fmid = exp((row.logf[j] + row.logf[j - 1]) / 2)
        qmid = (row.Q[j] + row.Q[j - 1]) / 2
        a = rvmid / rmid
        b = fmid * (1 - qmid^2 / rmid^2) / (4rmid)
        ru[j] = ((1 - dv * a / 2) * ru[j - 1] - dv * b) /
                (1 + dv * a / 2)
        mass[j] = renormalized_hawking_mass(
            row.r[j], exp(row.logf[j]), ru[j], rv[j], row.Q[j],
        )
    end
    return mass
end

function boundary_point(u, ep, fcorner)
    r0 = ep.rn.M + one(u) / 2
    r = ep.rn.M - u / 2
    Av = -ep.rn.Q0 * fcorner / 2 * (inv(r) - inv(r0))
    return NLPoint(
        u, zero(u), r, log(fcorner), zero(u), zero(u), zero(u), Av, ep.rn.Q0,
    )
end

function run_case(mode, ep, v; C, Umax, max_rows)
    initial, fcorner, rv = initial_row(mode, ep, v)
    initial_mass = initial_mass_profile(initial, rv)
    rows = [initial]
    while last(rows).u < Umax && length(rows) < max_rows
        previous = last(rows)
        du = 2C / exp(last(previous.logf))
        target_u = min(Umax, previous.u + du)
        target_u > previous.u || break
        south = boundary_point(target_u, ep, fcorner)
        next = advance_u_row(
            previous, south, ep;
            iterations=10, reduced_scalar=true, hyperbolic_charge=true,
        )
        push!(rows, next)
        all(isfinite, next.r) && all(>(zero(eltype(next.r))), next.r) || break
    end

    diagnostic = vtrap_diagnostic(
        rows;
        missing_status=last(rows).u + 2C / exp(last(last(rows).logf)) ==
                       last(rows).u ? :precision_stalled : :max_rows,
    )
    return (
        mode=mode,
        initial_M=last(initial_mass),
        initial_Q=last(initial.Q),
        initial_M_minus_Q=last(initial_mass) - last(initial.Q),
        fcorner=fcorner,
        rows=length(rows),
        last_U=last(rows).u,
        status=diagnostic.status,
        trap=diagnostic.trap,
        closest=diagnostic.closest,
    )
end

q0 = parse(Float64, argument(1, 1.0))
vmax = parse(Float64, argument(2, 400.0))
dv = parse(Float64, argument(3, 0.08))
max_rows = parse(Int, argument(4, 300))
C = parse(Float64, argument(5, 0.6))
ep = EvolutionParams(
    rn=RNParams(1.0, q0),
    scalar_charge=0.6 / q0,
    amplitude=0.01,
    omega=1.0,
)
v = collect(range(0.0, vmax; length=Int(round(vmax / dv)) + 1))

println("# GP2026 corner-compatible pulse-leg gauge comparison")
println("# Q0 = ", q0, ", Vmax = ", vmax, ", Delta V = ", dv, ", C = ", C)
for mode in (:ef_affine, :areal_affine)
    println(run_case(mode, ep, v; C, Umax=1.6, max_rows))
end
