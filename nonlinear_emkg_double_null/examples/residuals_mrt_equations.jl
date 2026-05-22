using NonlinearEMKGDoubleNull

function cell_residuals(nu, nv; V1=20.0)
    p = RNParams(1.0, 1.0)
    g = mrt2013_grid(; nu=nu, nv=nv, U0=-5.1, V0=0.0, U1=-1.0e-3, V1=V1)
    r = [mrt2013_areal_radius(U, V, p) for U in g.u, V in g.v]
    lf = [log(mrt2013_metric_f(U, V, p)) for U in g.u, V in g.v]

    r_res = Float64[]
    lf_res = Float64[]
    lf_old_res = Float64[]

    for i in 1:length(g.u)-1, j in 1:length(g.v)-1
        du = g.u[i + 1] - g.u[i]
        dv = g.v[j + 1] - g.v[j]
        rc = (r[i, j] + r[i + 1, j] + r[i, j + 1] + r[i + 1, j + 1]) / 4
        lfc = (lf[i, j] + lf[i + 1, j] + lf[i, j + 1] + lf[i + 1, j + 1]) / 4
        fc = exp(lfc)

        ru = ((r[i + 1, j] - r[i, j]) + (r[i + 1, j + 1] - r[i, j + 1])) / (2du)
        rv = ((r[i, j + 1] - r[i, j]) + (r[i + 1, j + 1] - r[i + 1, j])) / (2dv)
        ruv = ((r[i + 1, j + 1] - r[i, j + 1]) - (r[i + 1, j] - r[i, j])) / (du * dv)
        lfuv = ((lf[i + 1, j + 1] - lf[i, j + 1]) - (lf[i + 1, j] - lf[i, j])) / (du * dv)

        push!(r_res, rc * ruv + ru * rv + fc * (1 - p.Q0^2 / rc^2) / 4)
        push!(lf_res, lfuv - (fc / (2rc^2) + 2 * ru * rv / rc^2 - p.Q0^2 * fc / rc^4))
        push!(lf_old_res, lfuv - (-fc / (2rc^2) - 2 * ru * rv / rc^2 + p.Q0^2 * fc / rc^4))
    end

    c1 = Float64[]
    for i in eachindex(g.u), j in 2:length(g.v)-1
        dvp = g.v[j + 1] - g.v[j]
        dvm = g.v[j] - g.v[j - 1]
        rv_p = (r[i, j + 1] - r[i, j]) / dvp
        rv_m = (r[i, j] - r[i, j - 1]) / dvm
        y_p = rv_p / exp((lf[i, j + 1] + lf[i, j]) / 2)
        y_m = rv_m / exp((lf[i, j] + lf[i, j - 1]) / 2)
        push!(c1, (y_p - y_m) / ((dvp + dvm) / 2))
    end

    c2 = Float64[]
    for i in 2:length(g.u)-1, j in eachindex(g.v)
        dup = g.u[i + 1] - g.u[i]
        dum = g.u[i] - g.u[i - 1]
        ru_p = (r[i + 1, j] - r[i, j]) / dup
        ru_m = (r[i, j] - r[i - 1, j]) / dum
        y_p = ru_p / exp((lf[i + 1, j] + lf[i, j]) / 2)
        y_m = ru_m / exp((lf[i, j] + lf[i - 1, j]) / 2)
        push!(c2, (y_p - y_m) / ((dup + dum) / 2))
    end

    return maximum(abs, r_res), maximum(abs, lf_res), maximum(abs, lf_old_res),
           maximum(abs, c1), maximum(abs, c2)
end

for res in ((60, 180), (120, 360), (240, 720))
    vals = cell_residuals(res...)
    println("resolution = ", res)
    println("  Eq5 r residual          = ", vals[1])
    println("  Eq4 logf residual       = ", vals[2])
    println("  old-sign logf residual  = ", vals[3])
    println("  C1 constraint residual  = ", vals[4])
    println("  C2 constraint residual  = ", vals[5])
end

