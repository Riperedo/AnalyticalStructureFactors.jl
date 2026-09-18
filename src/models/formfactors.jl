# src/models/formfactors.jl

# ==============================================================================
# PARTICLE FORM FACTORS P(q) AND DECOUPLING FACTOR β(q)
# ==============================================================================
# This module implements analytical and semi-analytical form factors for
# various particle geometries (homogeneous spheres, core-shell, multi-shell,
# ellipsoids of revolution, and Y-shaped antibody macromolecules), as well
# as the Kotlarchyk-Chen decoupling factor β(q).
#
# References:
# - Kotlarchyk, M., & Chen, S. H. (1983). Analysis of small-angle neutron scattering
#   spectra from polydisperse interacting colloids. J. Chem. Phys., 79(5), 2461-2469.
# - Pedersen, J. S. (1997). Analysis of small-angle scattering data from colloids
#   and polymer solutions: modeling and least-squares fitting. Adv. Colloid Interface Sci., 70, 171-210.
# - Yearley, E.J. et al. (2013). Small-Angle Neutron Scattering Characterization of Monoclonal Antibody 
#   Conformations and Interactions at High Concentrations
# ==============================================================================

"""
    _f_sphere_amplitude(x::T) where {T<:AbstractFloat}

Calculates the normalized single-sphere scattering amplitude:
``f(x) = 3 \\frac{\\sin(x) - x \\cos(x)}{x^3}``
with a Taylor expansion for ``|x| < 10^{-4}`` to ensure numerical stability at ``x \\to 0``.
"""
function _f_sphere_amplitude(x::T) where {T<:AbstractFloat}
    x_abs = abs(x)
    if x_abs < T(1e-4)
        x2 = x * x
        # 3 * j1(x)/x = 1 - x^2/10 + x^4/280 - x^6/15120
        return one(T) - x2 / T(10.0) + (x2^2) / T(280.0) - (x2^3) / T(15120.0)
    else
        return T(3.0) * (sin(x) - x * cos(x)) / (x^3)
    end
end

# ------------------------------------------------------------------------------
# 1. HOMOGENEOUS SPHERE FORM FACTOR
# ------------------------------------------------------------------------------

@doc """
    P_sphere(q::T, R::T) where {T<:AbstractFloat}

Calculates the normalized form factor ``P(q)`` for a monodisperse homogeneous sphere of radius `R`.

``P(q) = \\left[ 3 \\frac{\\sin(q R) - q R \\cos(q R)}{(q R)^3} \\right]^2``

Normalized such that ``P(0) = 1``.

# Arguments
- `q::T`: Magnitude of the scattering wavevector.
- `R::T`: Sphere radius.

# Returns
- `T`: Normalized form factor ``P(q) \\in [0, 1]``.

# References
- Guinier, A., & Fournet, G. (1955). Small-Angle Scattering of X-Rays.
"""
function P_sphere(q::T, R::T) where {T<:AbstractFloat}
    if R < zero(T)
        throw(DomainError(R, "Sphere radius R must be non-negative."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end
    amp = _f_sphere_amplitude(q * R)
    return amp^2
end

"""
    P_sphere(R::T) where {T<:AbstractFloat}

Factory function returning a closure `q -> P_sphere(q, R)`.
"""
function P_sphere(R::T) where {T<:AbstractFloat}
    return q -> P_sphere(q, R)
end

"""
    P_sphere(q_values::AbstractVector{T}, R::T) where {T<:AbstractFloat}

Vectorized evaluation of `P_sphere` over an array of wavevectors `q_values`.
"""
function P_sphere(q_values::AbstractVector{T}, R::T) where {T<:AbstractFloat}
    return P_sphere.(q_values, R)
end

# ------------------------------------------------------------------------------
# 2. CORE-SHELL SPHERE FORM FACTOR
# ------------------------------------------------------------------------------

@doc """
    P_coreshell(q::T, R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}

Calculates the normalized form factor ``P(q)`` of a concentric core-shell sphere.

# Arguments
- `q::T`: Magnitude of the scattering wavevector.
- `R_core::T`: Radius of the core (``R_c``).
- `R_shell::T`: Overall radius of the particle including the outer shell (``R_s > R_c``).
- `sld_core::T`: Scattering length density (SLD) of the core (``\\rho_c``).
- `sld_shell::T`: Scattering length density of the shell (``\\rho_s``).
- `sld_solv::T`: Scattering length density of the surrounding solvent (``\\rho_{solv}``).

# Returns
- `T`: Normalized form factor ``P(q)`` with ``P(0) = 1``.

# References
- Pedersen, J. S. (1997). Adv. Colloid Interface Sci., 70, 171-210.
"""
function P_coreshell(q::T, R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}
    if R_core < zero(T) || R_shell < R_core
        throw(DomainError((R_core, R_shell), "Radii must satisfy 0 <= R_core <= R_shell."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end

    Vc = (T(4.0) * T(π) / T(3.0)) * (R_core^3)
    Vs = (T(4.0) * T(π) / T(3.0)) * (R_shell^3)

    contrast_c = sld_core - sld_shell
    contrast_s = sld_shell - sld_solv

    # Total scattering amplitude at q=0
    f0 = Vc * contrast_c + Vs * contrast_s
    if abs(f0) < eps(T)
        # Contrast matched system (f(0) = 0)
        fc = Vc * contrast_c * _f_sphere_amplitude(q * R_core)
        fs = Vs * contrast_s * _f_sphere_amplitude(q * R_shell)
        return (fc + fs)^2 / (Vs^2)
    end

    fq = Vc * contrast_c * _f_sphere_amplitude(q * R_core) + Vs * contrast_s * _f_sphere_amplitude(q * R_shell)
    return (fq / f0)^2
end

"""
    P_coreshell(R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}

Factory function returning a closure `q -> P_coreshell(q, R_core, R_shell, sld_core, sld_shell, sld_solv)`.
"""
function P_coreshell(R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}
    return q -> P_coreshell(q, R_core, R_shell, sld_core, sld_shell, sld_solv)
end

"""
    P_coreshell(q_values::AbstractVector{T}, R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}

Vectorized evaluation of `P_coreshell`.
"""
function P_coreshell(q_values::AbstractVector{T}, R_core::T, R_shell::T, sld_core::T, sld_shell::T, sld_solv::T) where {T<:AbstractFloat}
    return [P_coreshell(q, R_core, R_shell, sld_core, sld_shell, sld_solv) for q in q_values]
end

# ------------------------------------------------------------------------------
# 3. MULTI-SHELL CONCENTRIC SPHERES FORM FACTOR (N LAYERS)
# ------------------------------------------------------------------------------

@doc """
    P_multishell(q::T, radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}

Calculates the normalized form factor ``P(q)`` of a multi-layer concentric spherical particle
having ``N`` layers with outer radii ``R_1 < R_2 < \\dots < R_N`` and scattering length densities ``\\rho_1, \\rho_2, \\dots, \\rho_N``.

# Arguments
- `q::T`: Magnitude of the scattering wavevector.
- `radii::AbstractVector{T}`: Monotonically increasing outer radii of each shell ``[R_1, R_2, \\dots, R_N]``.
- `slds::AbstractVector{T}`: Scattering length densities of each shell ``[\\rho_1, \\rho_2, \\dots, \\rho_N]``.
- `sld_solv::T`: Solvent scattering length density ``\\rho_{solv}``.

# Returns
- `T`: Normalized form factor ``P(q)``.
"""
function P_multishell(q::T, radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}
    n_layers = length(radii)
    if length(slds) != n_layers || n_layers == 0
        throw(DimensionMismatch("radii and slds must have the same non-zero length."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end

    f0 = zero(T)
    fq = zero(T)

    for i in 1:n_layers
        R_i = radii[i]
        rho_i = slds[i]
        rho_next = (i == n_layers) ? sld_solv : slds[i+1]
        delta_rho = rho_i - rho_next
        V_i = (T(4.0) * T(π) / T(3.0)) * (R_i^3)

        f0 += V_i * delta_rho
        fq += V_i * delta_rho * _f_sphere_amplitude(q * R_i)
    end

    if abs(f0) < eps(T)
        V_outer = (T(4.0) * T(π) / T(3.0)) * (radii[end]^3)
        return (fq / V_outer)^2
    end

    return (fq / f0)^2
end

"""
    P_multishell(radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}

Factory function returning a closure `q -> P_multishell(q, radii, slds, sld_solv)`.
"""
function P_multishell(radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}
    return q -> P_multishell(q, radii, slds, sld_solv)
end

"""
    P_multishell(q_values::AbstractVector{T}, radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}

Vectorized evaluation of `P_multishell`.
"""
function P_multishell(q_values::AbstractVector{T}, radii::AbstractVector{T}, slds::AbstractVector{T}, sld_solv::T) where {T<:AbstractFloat}
    return [P_multishell(q, radii, slds, sld_solv) for q in q_values]
end

# ------------------------------------------------------------------------------
# 4. ELLIPSOID OF REVOLUTION FORM FACTOR AND KOTLARCHYK-CHEN DECOUPLING β(q)
# ------------------------------------------------------------------------------

@doc """
    P_ellipsoid(q::T, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Calculates the orientationally averaged form factor ``P(q)`` for an ellipsoid of revolution
with rotational semi-axis `Ra` and equatorial semi-axis `Rb` (aspect ratio ``\\nu = Ra/Rb``).

``P(q) = \\int_0^1 \\left[ 3 \\frac{\\sin(u) - u \\cos(u)}{u^3} \\right]^2 d\\mu``
where ``u = q R_b \\sqrt{1 + \\mu^2 ((Ra/Rb)^2 - 1)}``.

# Arguments
- `q::T`: Magnitude of the wavevector.
- `Ra::T`: Semi-axis of rotation (prolate if `Ra > Rb`, oblate if `Ra < Rb`).
- `Rb::T`: Equatorial semi-axis.
- `n_quad::Int=100`: Number of quadrature integration points over ``\\mu \\in [0, 1]``.

# Returns
- `T`: Normalized form factor ``P(q) \\in [0, 1]``.

# References
- Kotlarchyk, M., & Chen, S. H. (1983). J. Chem. Phys. 79, 2461.
"""
function P_ellipsoid(q::T, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    if Ra <= zero(T) || Rb <= zero(T)
        throw(DomainError((Ra, Rb), "Ellipsoid semi-axes Ra and Rb must be strictly positive."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end
    if q ≈ zero(T) || Ra ≈ Rb
        return P_sphere(q, Ra)
    end

    nu_sq_minus_one = (Ra / Rb)^2 - one(T)
    dmu = one(T) / T(n_quad)
    
    # Midpoint numerical quadrature
    integral = zero(T)
    for i in 1:n_quad
        mu = (T(i) - T(0.5)) * dmu
        u = q * Rb * sqrt(one(T) + (mu^2) * nu_sq_minus_one)
        amp = _f_sphere_amplitude(u)
        integral += amp^2
    end

    return integral * dmu
end

"""
    P_ellipsoid(Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Factory function returning a closure `q -> P_ellipsoid(q, Ra, Rb; n_quad=n_quad)`.
"""
function P_ellipsoid(Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    return q -> P_ellipsoid(q, Ra, Rb; n_quad=n_quad)
end

"""
    P_ellipsoid(q_values::AbstractVector{T}, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Vectorized evaluation of `P_ellipsoid`.
"""
function P_ellipsoid(q_values::AbstractVector{T}, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    return [P_ellipsoid(q, Ra, Rb; n_quad=n_quad) for q in q_values]
end

@doc """
    beta_ellipsoid(q::T, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Calculates the Kotlarchyk-Chen decoupling factor ``\\beta(q)`` for an ellipsoid of revolution:
``\\beta(q) = \\frac{\\langle F(q) \\rangle^2}{\\langle |F(q)|^2 \\rangle} = \\frac{\\left( \\int_0^1 f(u) d\\mu \\right)^2}{\\int_0^1 f(u)^2 d\\mu}``
where ``f(u) = 3 j_1(u)/u`` and ``u = q R_b \\sqrt{1 + \\mu^2 ((Ra/Rb)^2 - 1)}``.

Used in the decoupling approximation for non-spherical particles:
``I(q) = N P(q) \\left[ 1 + \\beta(q) (S_{eff}(q) - 1) \\right]``
where ``S_{eff}(q)`` is typically evaluated using the effective volume-equivalent hard-sphere radius
``R_{eff} = (Ra \\cdot Rb^2)^{1/3}``.

# Arguments
- `q::T`: Magnitude of the wavevector.
- `Ra::T`: Semi-axis of rotation.
- `Rb::T`: Equatorial semi-axis.
- `n_quad::Int=100`: Number of quadrature steps.

# Returns
- `T`: Decoupling factor ``\\beta(q) \\in [0, 1]``.

# References
- Kotlarchyk, M., & Chen, S. H. (1983). Analysis of small-angle neutron scattering
  spectra from polydisperse interacting colloids. J. Chem. Phys., 79(5), 2461-2469.
"""
function beta_ellipsoid(q::T, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    if Ra <= zero(T) || Rb <= zero(T)
        throw(DomainError((Ra, Rb), "Ellipsoid semi-axes Ra and Rb must be strictly positive."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end
    if q ≈ zero(T) || Ra ≈ Rb
        return one(T)
    end

    nu_sq_minus_one = (Ra / Rb)^2 - one(T)
    dmu = one(T) / T(n_quad)

    sum_F = zero(T)
    sum_F2 = zero(T)

    for i in 1:n_quad
        mu = (T(i) - T(0.5)) * dmu
        u = q * Rb * sqrt(one(T) + (mu^2) * nu_sq_minus_one)
        amp = _f_sphere_amplitude(u)
        sum_F += amp
        sum_F2 += amp^2
    end

    mean_F = sum_F * dmu
    mean_F2 = sum_F2 * dmu

    if mean_F2 < eps(T)
        return one(T)
    end

    return clamp((mean_F^2) / mean_F2, zero(T), one(T))
end

"""
    beta_ellipsoid(Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Factory function returning a closure `q -> beta_ellipsoid(q, Ra, Rb; n_quad=n_quad)`.
"""
function beta_ellipsoid(Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    return q -> beta_ellipsoid(q, Ra, Rb; n_quad=n_quad)
end

"""
    beta_ellipsoid(q_values::AbstractVector{T}, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}

Vectorized evaluation of `beta_ellipsoid`.
"""
function beta_ellipsoid(q_values::AbstractVector{T}, Ra::T, Rb::T; n_quad::Int=100) where {T<:AbstractFloat}
    return [beta_ellipsoid(q, Ra, Rb; n_quad=n_quad) for q in q_values]
end

"""
    effective_radius_ellipsoid(Ra::T, Rb::T) where {T<:AbstractFloat}

Calculates the equivalent spherical radius ``R_{eff} = (Ra \\cdot Rb^2)^{1/3}`` having the same volume as the ellipsoid.
"""
function effective_radius_ellipsoid(Ra::T, Rb::T) where {T<:AbstractFloat}
    return cbrt(Ra * (Rb^2))
end

# ------------------------------------------------------------------------------
# 5. DECOUPLED TOTAL SCATTERING INTENSITY
# ------------------------------------------------------------------------------

@doc """
    intensity_decoupled(P_val::T, beta_val::T, S_val::T) where {T<:AbstractFloat}

Calculates the effective structure factor contribution to scattering intensity under the decoupling approximation:
``I_{norm}(q) = P(q) \\left[ 1 + \\beta(q) (S(q) - 1) \\right]``

# Arguments
- `P_val::T`: Form factor value ``P(q)``.
- `beta_val::T`: Decoupling factor ``\\beta(q)``.
- `S_val::T`: Interparticle structure factor ``S(q)``.

# Returns
- `T`: Decoupled scattering intensity.
"""
function intensity_decoupled(P_val::T, beta_val::T, S_val::T) where {T<:AbstractFloat}
    return P_val * (one(T) + beta_val * (S_val - one(T)))
end
