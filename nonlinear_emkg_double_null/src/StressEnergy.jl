"""
Reduced matter sources for a spherically symmetric charged complex scalar.

Metric convention:

    ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2

so `g_uv = g_vu = -f/2` and `g^uv = g^vu = -2/f`.

The original MRT branch of the nonlinear solver stores the MRT-normalized scalar

    Phi = sqrt(32*pi) * phi_GP

where `phi_GP` is the canonically normalized field of Gelles-Pretorius
arXiv:2602.11256. This leaves the MRT uncharged equations in their original
normalization while allowing direct use of the charged equations after the
corresponding current rescaling.

For the GP2026 production branch, the state instead stores their evolved
reduced field

    Psi = r * Phi = sqrt(32*pi) * r * phi_GP.

`stress_energy_reduced_scalar` converts `Psi` and its coordinate derivatives
to `Phi` before evaluating the same source terms.

Gauge convention:

    D_a Phi = partial_a Phi - i e A_a Phi

For `Phi = phi_re + i phi_im`, this gives

    Re(D_a Phi) = partial_a phi_re + e A_a phi_im
    Im(D_a Phi) = partial_a phi_im - e A_a phi_re

The scalar components below are the canonical complex-scalar expression in
terms of `Phi`:

    T_ab^Phi = (D_a Phi)^* D_b Phi + (D_b Phi)^* D_a Phi
               - g_ab (D_c Phi)^* D^c Phi

This means `T_uu = 2 |D_u Phi|^2`. The nonlinear Raychaudhuri and metric
updates consume this reduced source with their MRT normalization factors.

The Faraday component is represented in the Gelles-Pretorius charge
normalization, converted to the solver metric coefficient:

    F_uv = -Q f/(2 r^2).

The Maxwell stress entries are provided algebraically for diagnostics. The
metric evolution carries the Coulomb terms explicitly through `Q`, rather
than feeding these entries back as an independent source.

The Maxwell expression below omits its Gaussian-unit overall `1/(4*pi)`:

    T_ab^EM = F_a{}^c F_bc - (1/4) g_ab F_cd F^cd

Supply `maxwell_weight=1/(4*pi)` when a physical Gaussian-unit Maxwell stress
is needed as output.
"""

struct StressEnergyComponents{T<:Real}
    Tuu::T
    Tvv::T
    Tuv::T
    Tthth::T
    Ju::T
    Jv::T
    alpha::T
    scalar_logf_source::T
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

    # F = alpha du wedge dv with the Gaussian-unit enclosed charge Q.
    alpha = -q * f / (2r^2)
    maxwell_Tuu = zero(r)
    maxwell_Tvv = zero(r)
    maxwell_Tuv = q^2 * f / (4r^4)
    maxwell_Tthth = q^2 / (2r^2)

    return StressEnergyComponents(
        scalar_weight * scalar_Tuu + maxwell_weight * maxwell_Tuu,
        scalar_weight * scalar_Tvv + maxwell_weight * maxwell_Tvv,
        scalar_weight * scalar_Tuv + maxwell_weight * maxwell_Tuv,
        scalar_weight * scalar_Tthth + maxwell_weight * maxwell_Tthth,
        Ju,
        Jv,
        alpha,
        scalar_weight * du_dot_dv / 2,
    )
end

function stress_energy_reduced_scalar(r, f, q, ru, rv, psi_re, psi_im,
                                      psi_re_u, psi_re_v, psi_im_u, psi_im_v,
                                      Au, Av, e;
                                      scalar_weight=one(r), maxwell_weight=one(r))
    phi_re = psi_re / r
    phi_im = psi_im / r
    phi_re_u = (psi_re_u - ru * phi_re) / r
    phi_re_v = (psi_re_v - rv * phi_re) / r
    phi_im_u = (psi_im_u - ru * phi_im) / r
    phi_im_v = (psi_im_v - rv * phi_im) / r
    return stress_energy(r, f, q, phi_re, phi_im, phi_re_u, phi_re_v,
                         phi_im_u, phi_im_v, Au, Av, e;
                         scalar_weight, maxwell_weight)
end

outgoing_constraint_source(r, f, source::StressEnergyComponents) = r * source.Tvv / (8f)
ingoing_constraint_source(r, f, source::StressEnergyComponents) = r * source.Tuu / (8f)
