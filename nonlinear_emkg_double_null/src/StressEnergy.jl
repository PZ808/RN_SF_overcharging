"""
Matter sources for a spherically symmetric charged complex scalar.

Metric convention:

    ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2

so `g_uv = g_vu = -f/2` and `g^uv = g^vu = -2/f`.

Gauge convention:

    D_a phi = partial_a phi - i e A_a phi

For `phi = phi_re + i phi_im`, this gives

    Re(D_a phi) = partial_a phi_re + e A_a phi_im
    Im(D_a phi) = partial_a phi_im - e A_a phi_re

The default scalar stress tensor is the canonical complex-scalar one:

    T_ab^phi = (D_a phi)^* D_b phi + (D_b phi)^* D_a phi
               - g_ab (D_c phi)^* D^c phi

This means `T_uu = 2 |D_u phi|^2` for the raw complex field. If later we
choose `phi = (phi1 + i phi2)/sqrt(2)`, the scalar contribution should be
multiplied by `scalar_weight = 1/2`.

The Maxwell stress tensor here is written for L_EM = -F_ab F^ab / 4:

    T_ab^EM = F_a{}^c F_bc - (1/4) g_ab F_cd F^cd

No overall `1/(4*pi)` is included. Put that normalization into the Einstein
equations, or pass `maxwell_weight`, once the action convention is frozen.
"""

struct StressEnergyComponents{T<:Real}
    Tuu::T
    Tvv::T
    Tuv::T
    Tthth::T
    Ju::T
    Jv::T
    alpha::T
end

function covariant_scalar_derivatives(phi_re, phi_im, phi_re_u, phi_re_v,
                                      phi_im_u, phi_im_v, Au, Av, e)
    du_re = phi_re_u + e * Au * phi_im
    du_im = phi_im_u - e * Au * phi_re
    dv_re = phi_re_v + e * Av * phi_im
    dv_im = phi_im_v - e * Av * phi_re
    return du_re, du_im, dv_re, dv_im
end

function current_components(phi_re, phi_im, du_re, du_im, dv_re, dv_im, e)
    Ju = 2e * (phi_re * du_im - phi_im * du_re)
    Jv = 2e * (phi_re * dv_im - phi_im * dv_re)
    return Ju, Jv
end

function stress_energy(r, f, q, phi_re, phi_im, phi_re_u, phi_re_v,
                       phi_im_u, phi_im_v, Au, Av, e;
                       scalar_weight=one(r), maxwell_weight=one(r))
    du_re, du_im, dv_re, dv_im =
        covariant_scalar_derivatives(phi_re, phi_im, phi_re_u, phi_re_v,
                                     phi_im_u, phi_im_v, Au, Av, e)

    Ju, Jv = current_components(phi_re, phi_im, du_re, du_im, dv_re, dv_im, e)

    du2 = du_re^2 + du_im^2
    dv2 = dv_re^2 + dv_im^2
    du_dot_dv = du_re * dv_re + du_im * dv_im

    # Scalar source. The massless scalar has T_uv = 0 in these coordinates.
    scalar_Tuu = 2 * du2
    scalar_Tvv = 2 * dv2
    scalar_Tuv = zero(r)
    scalar_Tthth = 4 * r^2 * du_dot_dv / f

    # F = alpha du wedge dv with alpha = Q f / r^2 in this convention.
    alpha = q * f / r^2
    maxwell_Tuu = zero(r)
    maxwell_Tvv = zero(r)
    maxwell_Tuv = alpha^2 / f
    maxwell_Tthth = 2 * r^2 * alpha^2 / f^2

    return StressEnergyComponents(
        scalar_weight * scalar_Tuu + maxwell_weight * maxwell_Tuu,
        scalar_weight * scalar_Tvv + maxwell_weight * maxwell_Tvv,
        scalar_weight * scalar_Tuv + maxwell_weight * maxwell_Tuv,
        scalar_weight * scalar_Tthth + maxwell_weight * maxwell_Tthth,
        Ju,
        Jv,
        alpha,
    )
end

outgoing_constraint_source(r, f, source::StressEnergyComponents) = r * source.Tvv / (8f)
ingoing_constraint_source(r, f, source::StressEnergyComponents) = r * source.Tuu / (8f)
