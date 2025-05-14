# src/monodisperse.jl (o un archivo de utilidades si son muy generales)

# --- Helper Functions for Percus-Yevick (Hard Spheres) ---
# These are typically not exported and are used internally.
# Prefixing with an underscore `_` is a common convention for internal helper functions.

"""
    _alpha1_py(ϕ::T) where {T<:AbstractFloat}

Coefficient α₁ for the Percus-Yevick direct correlation function for hard spheres.
Depends on volume fraction `ϕ`.
"""
_alpha1_py(ϕ::T) where {T<:AbstractFloat} = -(one(T) + T(2.0) * ϕ)^2 / (one(T) - ϕ)^4

"""
    _alpha2_py(ϕ::T) where {T<:AbstractFloat}

Coefficient α₂ for the Percus-Yevick direct correlation function for hard spheres.
Depends on volume fraction `ϕ`.
"""
_alpha2_py(ϕ::T) where {T<:AbstractFloat} = T(6.0) * ϕ * (one(T) + T(0.5) * ϕ)^2 / (one(T) - ϕ)^4

"""
    _alpha3_py(ϕ::T) where {T<:AbstractFloat}

Coefficient α₃ for the Percus-Yevick direct correlation function for hard spheres.
Depends on volume fraction `ϕ` and `_alpha1_py`.
"""
_alpha3_py(ϕ::T) where {T<:AbstractFloat} = T(0.5) * ϕ * _alpha1_py(ϕ)


# --- Useful Integrals for Fourier Transform of Direct Correlation Function ---
# These are also typically internal helper functions.

"""
    _integral1_py(k::T) where {T<:AbstractFloat}

Integral ∫(0,1) dx x² j₀(kx)dx, where j₀ is the spherical Bessel function of order 0.
This is equivalent to ∫(0,1) dx x² (sin(kx)/kx) dx.
The diameter of the sphere is scaled out and re-introduced via `k` (i.e. k = q*sigma).
"""
function _integral1_py(k::T) where {T<:AbstractFloat}
    k_abs = abs(k) # Use abs(k) if k can be negative, though wavevector magnitude is usually positive
    if k_abs < eps(T) # Check if k is very close to zero
        return one(T) / T(3.0)
    else
        # sin(k_abs)/k_abs^3 - cos(k_abs)/k_abs^2 based on common forms,
        # or (sin(k_abs) - k_abs*cos(k_abs))/k_abs^3
        return (sin(k_abs) - k_abs * cos(k_abs)) / k_abs^3
    end
end

"""
    _integral2_py(k::T) where {T<:AbstractFloat}

Integral ∫(0,1) dx x³ j₀(kx)dx. (Original comment said x² xj₀(kx), which is x³)
"""
function _integral2_py(k::T) where {T<:AbstractFloat}
    k_abs = abs(k)
    if k_abs < eps(T)
        return one(T) / T(4.0)
    else
        # Original: (-(k^2 - 2.0)*cos(k) + 2.0*k*sin(k) - 2.0)/k^4
        # This is for ∫ x³ j₀(x k) dx from 0 to 1.
        # Let's verify or use a known form if available.
        # WolframAlpha for integral x^3 sin(k x)/(k x) dx from 0 to 1:
        # ( (2-k^2)cos(k) + 2k sin(k) - 2 ) / k^4
        # This matches the user's original form, just with k_abs
        return ( (T(2.0) - k_abs^2) * cos(k_abs) + T(2.0) * k_abs * sin(k_abs) - T(2.0) ) / k_abs^4
    end
end

"""
    _integral3_py(k::T) where {T<:AbstractFloat}

Integral ∫(0,1) dx x⁵ j₀(kx)dx. (Original comment said x² x³j₀(kx), which is x⁵)
"""
function _integral3_py(k::T) where {T<:AbstractFloat}
    k_abs = abs(k)
    if k_abs < eps(T)
        return one(T) / T(6.0)
    else
        # Original: (4.0*k*(k^2 - 6.0)*sin(k) - (k^4 - 12.0*k^2 + 24.0)*cos(k) + 24.0)/k^6
        # This is for ∫ x⁵ j₀(x k) dx from 0 to 1.
        # WolframAlpha for integral x^5 sin(k x)/(k x) dx from 0 to 1:
        # ( (24 - 12k^2 + k^4)cos(k) + (6k^2 - 24) 2k sin(k) - 24 + k^4 sin(k) ) / k^6
        # ( (k^4 - 12k^2 + 24)cos(k) + (4k^3 - 24k)sin(k) - 24 ) / k^6
        # This matches the user's original form (with a sign flip for cos term due to rearrangement)
        return ( T(4.0)*k_abs*(k_abs^2 - T(6.0))*sin(k_abs) - (k_abs^4 - T(12.0)*k_abs^2 + T(24.0))*cos(k_abs) + T(24.0) ) / k_abs^6
    end
end


# --- Public API Functions ---
# These functions would typically be exported by the module.

@doc """
    C_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}

Computes the dimensionless direct correlation function c(k*sigma) for hard spheres using the Percus-Yevick (PY) approximation.
The result is c(q*sigma), where q is the wavevector magnitude and sigma is the particle diameter.
The input `k` is expected to be the dimensionless wavevector k = q*sigma.

# Arguments
- `ϕ::T`: The volume fraction of spheres. Should be in the range `[0, 1)`.
- `k::T`: The dimensionless wavevector, k = q*sigma, where q is the wavevector and sigma is the sphere diameter.

# Returns
- Value of the direct correlation function c(k).

# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
function C_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}
    if !(zero(T) <= ϕ < one(T))
        @warn "Volume fraction ϕ = $ϕ is outside the typical physical range [0, 1). Results may be unphysical."
    end
    # It's assumed k = qσ (dimensionless wavevector, where σ is the hard-sphere diameter)
    # The integrals I1, I2, I3 are defined for k = qσ.
    # The final c(k) is also dimensionless.
    return _alpha1_py(ϕ) * _integral1_py(k) +
           _alpha2_py(ϕ) * _integral2_py(k) +
           _alpha3_py(ϕ) * _integral3_py(k)
end

@doc """
    IS_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}

Computes the inverse of the static structure factor 1/S(k) for hard spheres
using the direct correlation function c(k) from the Percus-Yevick (PY) approximation.
The inverse structure factor is given by 1/S(k) = 1 - n * c_tilde(q), where n is number density and c_tilde(q) is the FT of c(r).
For PY hard spheres, this is often written as 1/S(k) = 1 - 24*ϕ*c_PY(k),
where c_PY(k) is the dimensionless form computed by `C_HS_PY`.

# Arguments
- `ϕ::T`: The volume fraction of spheres. Should be in the range `[0, 1)`.
- `k::T`: The dimensionless wavevector, k = q*sigma.

# Returns
- Value of the inverse static structure factor 1/S(k).

# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
function IS_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}
    # The factor -24ϕ comes from n * integral of c_PY(r), where n = 6ϕ/(πσ³),
    # and the integral involves σ.
    # If C_HS_PY returns c(kσ) directly as used in S(k) = 1 / (1 - ρ ĉ(k))
    # where ĉ(k) is the FT of c(r).
    # The common form is 1/S(k) = 1 - ρ_norm * C_HS_PY(ϕ,k)
    # where ρ_norm = -24ϕ for the Wertheim/Thiele solution.
    # So the formula 1.0 - 24*ϕ*C_HS_PY(ϕ, k) is correct if C_HS_PY is defined as it is.
    return one(T) - T(24.0) * ϕ * C_HS_PY(ϕ, k)
end

@doc """
    S_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}

Returns the Static Structure Factor S(k) for hard spheres under the Percus-Yevick (PY) Closure.
The input `k` is expected to be the dimensionless wavevector k = q*sigma.

The formula used is:
S(k) = 1 / (1 - 24*ϕ*c(k))

# Arguments
- `ϕ::T`: The volume fraction. Should be in the range `[0, 1)`.
- `k::T`: The dimensionless wavevector, k = q*sigma.

# Returns
- Value of the static structure factor S(k).

# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
function S_HS_PY(ϕ::T, k::T) where {T<:AbstractFloat}
    # Input validation for ϕ is handled by C_HS_PY, which is called by IS_HS_PY.
    # Consider if S(k) can be zero or negative, which would be unphysical.
    # IS_HS_PY can be zero or negative.
    # If IS_HS_PY is zero, S(k) would be Inf. If negative, S(k) would be negative.
    # This can happen for unphysical ϕ or k values in the PY approx.
    inv_S = IS_HS_PY(ϕ, k)
    if inv_S <= zero(T) && !(abs(k) < eps(T) && abs(ϕ) < eps(T)) # Avoid warning for S(0) in ideal gas limit
        @warn "Inverse structure factor 1/S(k) = $inv_S is non-positive for ϕ=$ϕ, k=$k. Resulting S(k) will be unphysical (Inf or negative)."
    end
    return one(T) / inv_S
end

