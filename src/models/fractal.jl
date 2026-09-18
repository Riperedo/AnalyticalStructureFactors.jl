# src/models/fractal.jl

# ==============================================================================
# FRACTAL AGGREGATE STRUCTURE FACTORS S(q)
# ==============================================================================
# This module implements structure factors for fractal colloidal aggregates
# and mass-fractal polymer/protein networks based on the Teixeira-Chen model
# and hard-sphere modulated fractal clusters.
#
# References:
# - Chen, S.-H., & Teixeira, J. (1986). Structure and Fractal Dimension of
#   Protein-Detergent Complexes. Phys. Rev. Lett., 57(20), 2583-2586.
# - Teixeira, J. (1988). Small-angle scattering by fractal systems.
#   J. Appl. Crystallogr., 21(6), 781-785.
# - Freltoft, T., Kjems, J. K., & Sinha, S. K. (1986). Power-law correlations
#   and finite-size effects in silica particle aggregates studied by small-angle neutron scattering.
#   Phys. Rev. B, 33(1), 269-275.
# ==============================================================================

using SpecialFunctions: gamma

@doc """
    S_Fractal_Chen(q::T, R::T, D::T, ξ::T) where {T<:AbstractFloat}

Calculates the analytical static structure factor ``S(q)`` for mass-fractal aggregates
of spherical primary particles according to the Chen-Teixeira model.

``S(q) = 1 + \\frac{D \\cdot \\Gamma(D - 1) \\sin\\left((D - 1) \\arctan(q \\xi)\\right)}{(q R)^D \\left[1 + (q \\xi)^{-2}\\right]^{(D - 1)/2}}``

# Arguments
- `q::T`: Scattering wavevector magnitude (in units of ``\\text{length}^{-1}``).
- `R::T`: Primary particle radius (``R > 0``).
- `D::T`: Fractal dimension (typically ``1 < D < 3``, e.g., ``D \\approx 1.8`` for DLCA or ``D \\approx 2.1`` for RLCA).
- `ξ::T`: Characteristic aggregate correlation / cutoff length (``\\xi > R``).

# Returns
- `T`: Static structure factor ``S(q)``.

# Limiting Behaviors
- **``q \\to 0`` (Clustering Limit)**:
  ``\\lim_{q \\to 0} S(q) = 1 + D \\Gamma(D) \\left(\\frac{\\xi}{R}\\right)^D`` (related to the average aggregate mass).
- **``1/\\xi \\ll q \\ll 1/R`` (Fractal Power Law)**:
  ``S(q) \\propto (q R)^{-D}``.
- **``q \\gg 1/R`` (High-``q`` Asymptote)**:
  ``S(q) \\to 1``.

# References
- Chen, S.-H., & Teixeira, J. (1986). Phys. Rev. Lett., 57, 2583.
- Teixeira, J. (1988). J. Appl. Crystallogr., 21, 781.
"""
function S_Fractal_Chen(q::T, R::T, D::T, ξ::T) where {T<:AbstractFloat}
    if R <= zero(T)
        throw(DomainError(R, "Primary particle radius R must be strictly positive."))
    end
    if D <= one(T) || D >= T(3.0)
        throw(DomainError(D, "Fractal dimension D should typically be in the range (1, 3)."))
    end
    if ξ <= zero(T)
        throw(DomainError(ξ, "Correlation cutoff length ξ must be strictly positive."))
    end
    if q < zero(T)
        throw(DomainError(q, "Wavevector magnitude q must be non-negative."))
    end

    # Handle q -> 0 smoothly using Taylor limit
    if q * ξ < T(1e-4)
        # Using sin(alpha * atan(x)) ~ alpha * x for small x
        # S(0) = 1 + D * Gamma(D-1) * (D-1) * (xi/R)^D = 1 + D * Gamma(D) * (xi/R)^D
        return one(T) + D * gamma(D) * ((ξ / R)^D)
    end

    qR = q * R
    qξ = q * ξ
    alpha = D - one(T)

    # Numerically stable formulation:
    # [1 + (q*xi)^(-2)]^(-alpha/2) = (q*xi)^alpha / [1 + (q*xi)^2]^(alpha/2)
    # Term = D * Gamma(alpha) * (xi/R)^D * sin(alpha * atan(q*xi)) / [ q*xi * (1 + (q*xi)^2)^(alpha/2) ]
    numerator = D * gamma(alpha) * ((ξ / R)^D) * sin(alpha * atan(qξ))
    denominator = qξ * ((one(T) + qξ^2)^(alpha / T(2.0)))

    return one(T) + numerator / denominator
end

"""
    S_Fractal_Chen(R::T, D::T, ξ::T) where {T<:AbstractFloat}

Factory function returning a closure `q -> S_Fractal_Chen(q, R, D, ξ)`.
"""
function S_Fractal_Chen(R::T, D::T, ξ::T) where {T<:AbstractFloat}
    return q -> S_Fractal_Chen(q, R, D, ξ)
end

"""
    S_Fractal_Chen(q_values::AbstractVector{T}, R::T, D::T, ξ::T) where {T<:AbstractFloat}

Vectorized evaluation of `S_Fractal_Chen`.
"""
function S_Fractal_Chen(q_values::AbstractVector{T}, R::T, D::T, ξ::T) where {T<:AbstractFloat}
    return [S_Fractal_Chen(q, R, D, ξ) for q in q_values]
end
