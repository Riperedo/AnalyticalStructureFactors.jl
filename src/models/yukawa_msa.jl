# src/models/yukawa_msa.jl

# ==============================================================================
# SINGLE-YUKAWA FLUID - EXACT ANALYTICAL MEAN SPHERICAL APPROXIMATION (MSA)
# ==============================================================================
# This module implements the exact analytical solution of the Ornstein-Zernike (OZ)
# equation with Mean Spherical Approximation (MSA) closure for a monodisperse fluid
# of hard spheres with an attractive or repulsive Yukawa tail:
#
#   V(r) / (k_B T) = -K * exp(-z * (r - 1)) / r   (for r > 1, with σ = 1)
#
# The solution uses the Baxter Wiener-Hopf factorization method developed by
# Waisman (1973) and Høye & Blum (1977). The Baxter scaling parameter `d` is found
# as the physical root of a 4th-degree polynomial.
#
# References:
# - Waisman, E. (1973). The radial distribution function for a fluid of hard spheres at high densities:
#   mean spherical integral equation approach. Molecular Physics, 25(1), 45-48.
# - Høye, J. S., & Blum, L. (1977). Solution of the Yukawa closure of the Ornstein-Zernike equation.
#   J. Stat. Phys., 16(5), 399-413.
# - Cummings, P. T., & Smith, W. R. (1979). Chem. Phys. Lett., 62(1), 187-190.
# ==============================================================================

using LinearAlgebra: eigvals

"""
    _poly_roots(coeffs::AbstractVector{T}) where {T<:AbstractFloat}

Calculates all complex roots of a polynomial ``c_n x^n + c_{n-1} x^{n-1} + \\dots + c_0 = 0``
using the eigenvalues of its companion matrix.
"""
function _poly_roots(coeffs::AbstractVector{T}) where {T<:AbstractFloat}
    max_c = maximum(abs, coeffs)
    if max_c == zero(T)
        return Complex{T}[]
    end
    lead_idx = findfirst(x -> abs(x) > 1e-15 * max_c, coeffs)
    if isnothing(lead_idx)
        return Complex{T}[]
    end
    c = coeffs[lead_idx:end]
    deg = length(c) - 1
    if deg == 0
        return Complex{T}[]
    elseif deg == 1
        return [-Complex{T}(c[2] / c[1])]
    end
    A = zeros(Complex{T}, deg, deg)
    for i in 1:deg-1
        A[i+1, i] = one(T)
    end
    for i in 1:deg
        val = -c[deg - i + 2] / c[1]
        if !isfinite(val)
            c_big = BigFloat.(coeffs)
            return Complex{T}.(_poly_roots(c_big))
        end
        A[i, deg] = val
    end
    try
        return eigvals(A)
    catch
        if T == BigFloat
            return Complex{T}[]
        else
            c_big = BigFloat.(coeffs)
            return Complex{T}.(_poly_roots(c_big))
        end
    end
end

"""
    yukawa_msa_roots(ϕ::T, K::T, z::T) where {T<:AbstractFloat}

Calculates and returns all real roots for the Baxter scaling parameter `d` from the 4th-degree
polynomial equation of the Ornstein-Zernike equation with Mean Spherical Approximation (MSA).
"""
function yukawa_msa_roots(ϕ::T, K::T, z::T) where {T<:AbstractFloat}
    if abs(K) < eps(T)
        return T[zero(T)]
    end

    Ez = exp(z)
    E2z = exp(T(2.0) * z)
    E3z = exp(T(3.0) * z)
    E4z = exp(T(4.0) * z)

    # 4th-degree polynomial in d: para1*d^4 + para2*d^3 + para3*d^2 + para4*d + para5 = 0
    para1 = -T(36.0) * (z^4) * (ϕ^2) + T(72.0) * (z^4) * (ϕ^3) - T(36.0) * (z^4) * (ϕ^4)
    para2 = T(6.0) * (z^5) * ϕ + T(72.0) * (z^2) * (ϕ^2) - T(72.0) * Ez * (z^2) * (ϕ^2) +
            T(72.0) * Ez * (z^3) * (ϕ^2) - T(36.0) * (z^4) * (ϕ^2) - T(12.0) * (z^5) * (ϕ^2) +
            T(144.0) * (z^2) * (ϕ^3) - T(144.0) * Ez * (z^2) * (ϕ^3) + T(108.0) * (z^3) * (ϕ^3) +
            T(36.0) * Ez * (z^3) * (ϕ^3) + T(36.0) * (z^4) * (ϕ^3) + T(6.0) * (z^5) * (ϕ^3)
    para3 = -T(12.0) * E2z * K * (z^4) * ϕ + T(24.0) * E2z * K * (z^4) * (ϕ^2) - T(12.0) * E2z * K * (z^4) * (ϕ^3)
    para4 = E3z * K * (z^5) + T(12.0) * E2z * K * (z^2) * ϕ - T(12.0) * E3z * K * (z^2) * ϕ +
            T(12.0) * E2z * K * (z^3) * ϕ + T(6.0) * E3z * K * (z^4) * ϕ - T(2.0) * E3z * K * (z^5) * ϕ +
            T(24.0) * E2z * K * (z^2) * (ϕ^2) - T(24.0) * E3z * K * (z^2) * (ϕ^2) +
            6.0 * E2z * K * (z^3) * (ϕ^2) + 18.0 * E3z * K * (z^3) * (ϕ^2) -
            6.0 * E3z * K * (z^4) * (ϕ^2) + E3z * K * (z^5) * (ϕ^2)
    para5 = -E4z * (K^2) * (z^4) + T(2.0) * E4z * (K^2) * (z^4) * ϕ - E4z * (K^2) * (z^4) * (ϕ^2)

    all_roots = _poly_roots([para1, para2, para3, para4, para5])
    real_roots = T[]
    for r in all_roots
        if abs(imag(r)) < T(1e-4)
            push!(real_roots, real(r))
        end
    end
    return real_roots
end
yukawa_msa_roots(ϕ::Real, K::Real, z::Real) = yukawa_msa_roots(promote(Float64(ϕ), Float64(K), Float64(z))...)

"""
    yukawa_msa_beta_roots(ϕ::T, K::T, z::T) where {T<:AbstractFloat}

Calculates and returns all real roots for Cummings & Smith (1979) Baxter parameter
`beta = z * d / (2 * exp(z))` as a function of packing fraction `ϕ` and isotherm `K`.
"""
function yukawa_msa_beta_roots(ϕ::T, K::T, z::T) where {T<:AbstractFloat}
    d_roots = yukawa_msa_roots(ϕ, K, z)
    Ez = exp(z)
    return sort([z * d / (T(2.0) * Ez) for d in d_roots])
end
yukawa_msa_beta_roots(ϕ::Real, K::Real, z::Real) = yukawa_msa_beta_roots(promote(Float64(ϕ), Float64(K), Float64(z))...)

"""
    yukawa_msa_coefficients(ϕ::T, K::T, z::T) where {T<:AbstractFloat}

Solves the quartic algebraic polynomial for the Baxter parameter `d` and returns
the coefficients `(a, b, c, d)` of the analytical physical MSA solution.
"""
function yukawa_msa_coefficients(ϕ::T, K::T, z::T) where {T<:AbstractFloat}
    Ez = exp(z)
    E2z = exp(T(2.0) * z)
    E3z = exp(T(3.0) * z)
    E4z = exp(T(4.0) * z)

    if abs(K) < eps(T)
        a_val = (one(T) + T(2.0) * ϕ) / ((one(T) - ϕ)^2)
        b_val = -T(3.0) * ϕ / (T(2.0) * ((one(T) - ϕ)^2))
        return (a=a_val, b=b_val, c=zero(T), d=zero(T))
    end

    real_roots = yukawa_msa_roots(ϕ, K, z)

    if isempty(real_roots)
        d = zero(T)
    else
        # Select the physical root connecting continuously to hard-spheres (d -> 0 as K -> 0)
        d = real_roots[argmin(abs.(real_roots))]
    end

    # Denominator for a, b, c
    den = T(6.0) * d * ϕ * (-T(2.0) * E2z * K * (z^2) + T(2.0) * E3z * K * (z^2) - E3z * K * (z^3) +
          T(12.0) * E2z * K * ϕ - T(12.0) * E3z * K * ϕ + T(12.0) * E3z * K * z * ϕ -
          T(12.0) * (d^2) * (z^2) * ϕ + T(12.0) * (d^2) * Ez * (z^2) * ϕ - E2z * K * (z^2) * ϕ -
          T(5.0) * E3z * K * (z^2) * ϕ - T(6.0) * (d^2) * (z^3) * ϕ + E3z * K * (z^3) * ϕ +
          T(72.0) * (d^2) * (ϕ^2) - T(72.0) * (d^2) * Ez * (ϕ^2) + T(72.0) * (d^2) * z * (ϕ^2) +
          T(30.0) * (d^2) * (z^2) * (ϕ^2) + T(6.0) * (d^2) * Ez * (z^2) * (ϕ^2) + T(6.0) * (d^2) * (z^3) * (ϕ^2))

    num_a = -E4z * (K^2) * (z^3) + d * E3z * K * (z^4) + T(12.0) * d * E2z * K * z * ϕ -
            T(12.0) * d * E3z * K * z * ϕ + T(12.0) * d * E3z * K * (z^2) * ϕ +
            T(6.0) * E4z * (K^2) * (z^2) * ϕ - T(12.0) * (d^2) * E2z * K * (z^3) * ϕ -
            T(6.0) * d * E3z * K * (z^3) * ϕ - T(2.0) * E4z * (K^2) * (z^3) * ϕ +
            T(6.0) * (d^3) * (z^4) * ϕ + T(2.0) * d * E3z * K * (z^4) * ϕ +
            T(72.0) * (d^3) * z * (ϕ^2) - T(72.0) * (d^3) * Ez * z * (ϕ^2) +
            T(72.0) * (d^3) * (z^2) * (ϕ^2) + T(36.0) * (d^3) * (z^3) * (ϕ^2) -
            T(36.0) * (d^4) * (z^3) * (ϕ^2) + T(48.0) * (d^2) * E2z * K * (z^3) * (ϕ^2) +
            T(12.0) * (d^3) * (z^4) * (ϕ^2) - T(216.0) * (d^4) * (z^2) * (ϕ^3) -
            T(72.0) * (d^4) * (z^3) * (ϕ^3)

    num_b = T(2.0) * E4z * (K^2) * (z^2) - T(2.0) * d * E3z * K * (z^3) - T(24.0) * d * E2z * K * ϕ +
            T(24.0) * d * E3z * K * ϕ - T(24.0) * d * E2z * K * z * ϕ +
            T(24.0) * (d^2) * E2z * K * (z^2) * ϕ - T(12.0) * d * E3z * K * (z^2) * ϕ -
            T(8.0) * E4z * (K^2) * (z^2) * ϕ - T(12.0) * (d^3) * (z^3) * ϕ +
            T(8.0) * d * E3z * K * (z^3) * ϕ + T(3.0) * E4z * (K^2) * (z^3) * ϕ -
            T(3.0) * d * E3z * K * (z^4) * ϕ - T(144.0) * (d^3) * (ϕ^2) +
            T(144.0) * (d^3) * Ez * (ϕ^2) - T(288.0) * (d^3) * z * (ϕ^2) +
            T(144.0) * (d^3) * Ez * z * (ϕ^2) + T(144.0) * (d^2) * E2z * K * z * (ϕ^2) -
            T(216.0) * (d^3) * (z^2) * (ϕ^2) + T(72.0) * (d^4) * (z^2) * (ϕ^2) -
            T(96.0) * (d^2) * E2z * K * (z^2) * (ϕ^2) - T(96.0) * (d^3) * (z^3) * (ϕ^2) -
            T(36.0) * (d^2) * E2z * K * (z^3) * (ϕ^2) - T(18.0) * (d^3) * (z^4) * (ϕ^2) +
            T(864.0) * (d^4) * z * (ϕ^3) + T(576.0) * (d^4) * (z^2) * (ϕ^3) +
            T(108.0) * (d^4) * (z^3) * (ϕ^3)

    num_c = -d * Ez * (T(2.0) * E2z * K * (z^2) - T(2.0) * d * (z^3) - E2z * K * (z^3) -
            T(12.0) * E2z * K * ϕ + T(12.0) * E2z * K * z * ϕ + T(12.0) * (d^2) * (z^2) * ϕ -
            T(5.0) * E2z * K * (z^2) * ϕ - d * (z^3) * ϕ + E2z * K * (z^3) * ϕ -
            T(72.0) * (d^2) * (ϕ^2) + T(6.0) * (d^2) * (z^2) * (ϕ^2))

    den_c = -T(2.0) * E2z * K * (z^2) + T(2.0) * E3z * K * (z^2) - E3z * K * (z^3) +
            T(12.0) * E2z * K * ϕ - T(12.0) * E3z * K * ϕ + T(12.0) * E3z * K * z * ϕ -
            T(12.0) * (d^2) * (z^2) * ϕ + T(12.0) * (d^2) * Ez * (z^2) * ϕ - E2z * K * (z^2) * ϕ -
            T(5.0) * E3z * K * (z^2) * ϕ - T(6.0) * (d^2) * (z^3) * ϕ + E3z * K * (z^3) * ϕ +
            T(72.0) * (d^2) * (ϕ^2) - T(72.0) * (d^2) * Ez * (ϕ^2) + T(72.0) * (d^2) * z * (ϕ^2) +
            T(30.0) * (d^2) * (z^2) * (ϕ^2) + T(6.0) * (d^2) * Ez * (z^2) * (ϕ^2) + T(6.0) * (d^2) * (z^3) * (ϕ^2)

    if abs(den) < eps(T) || abs(den_c) < eps(T)
        # Fallback to Hard Sphere Percus-Yevick limits
        a_val = (one(T) + T(2.0) * ϕ) / ((one(T) - ϕ)^2)
        b_val = -T(3.0) * ϕ / (T(2.0) * ((one(T) - ϕ)^2))
        c_val = zero(T)
        d = zero(T)
    else
        a_val = num_a / den
        b_val = num_b / (T(2.0) * den)
        c_val = num_c / den_c
    end

    return (a=a_val, b=b_val, c=c_val, d=d)
end
const _yukawa_msa_coefficients = yukawa_msa_coefficients

"""
    _P_factor_yukawa_msa(k::T, ϕ::T, z::T, a::T, b::T, c::T, d::T) where {T<:AbstractFloat}

Calculates Baxter's complex factor ``P(k) = 1 - 2\\pi \\rho \\hat{Q}(k)`` for the Yukawa MSA model.
"""
function _P_factor_yukawa_msa(k::T, ϕ::T, z::T, a::T, b::T, c::T, d::T) where {T<:AbstractFloat}
    W = im
    # Threshold for evaluating direct formula vs quadratic continuation to prevent float roundoff at k -> 0
    k_eval = max(k, T(1e-3))
    Ez_inv = exp(-z)
    E_ik = exp(W * k_eval)

    term1 = one(Complex{T})
    term2 = (T(12.0) * d * ϕ) / (W * k_eval - z)
    term3 = (T(6.0) * a * ϕ + T(12.0) * b * ϕ - T(12.0) * c * ϕ + T(12.0) * c * Ez_inv * ϕ) / (-W * k_eval + z)
    term4 = (T(12.0) * a * ϕ - T(12.0) * a * E_ik * ϕ) / ((k_eval^2) * (-W * k_eval + z))
    term5 = (-T(12.0) * W * b * ϕ + T(12.0) * W * a * E_ik * ϕ + T(12.0) * W * b * E_ik * ϕ) / (k_eval * (-W * k_eval + z))
    term6 = (T(12.0) * W * a * z * ϕ - T(12.0) * W * a * E_ik * z * ϕ) / ((k_eval^3) * (-W * k_eval + z))
    term7 = (T(12.0) * b * z * ϕ - T(12.0) * a * E_ik * z * ϕ - T(12.0) * b * E_ik * z * ϕ) / ((k_eval^2) * (-W * k_eval + z))
    term8 = (T(6.0) * W * a * z * ϕ + T(12.0) * W * b * z * ϕ - T(12.0) * W * c * exp(W * k_eval - z) * z * ϕ + T(12.0) * W * c * Ez_inv * z * ϕ) / (k_eval * (-W * k_eval + z))

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8
end

# ------------------------------------------------------------------------------
# PUBLIC API: S_Yukawa_MSA, IS_Yukawa_MSA, C_Yukawa_MSA
# ------------------------------------------------------------------------------

@doc """
    IS_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}

Calculates the inverse static structure factor ``1 / S(k)`` for a single Yukawa fluid
within the exact analytical Mean Spherical Approximation (MSA).

``\\frac{1}{S(k)} = |P(k)|^2``

# Arguments
- `ϕ::T`: Volume fraction (``0 \\le \\phi < 1``).
- `K::T`: Dimensionless contact interaction energy ``K = \\epsilon / (k_B T)``. Positive for attraction, negative for repulsion.
- `z::T`: Dimensionless inverse screening length ``z = \\kappa \\sigma`` (``z > 0``).
- `k::T`: Dimensionless wavevector ``k = q \\sigma``.

# Returns
- `T`: Inverse static structure factor ``IS(k)``.

# References
- Waisman, E. (1973). Mol. Phys., 25(1), 45-48.
- Høye, J. S., & Blum, L. (1977). J. Stat. Phys., 16(5), 399-413.
"""
function IS_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}
    if ϕ < zero(T) || ϕ >= one(T)
        throw(DomainError(ϕ, "Volume fraction ϕ must satisfy 0 <= ϕ < 1."))
    end
    if z <= zero(T)
        throw(DomainError(z, "Screening parameter z must be strictly positive."))
    end
    if k < zero(T)
        throw(DomainError(k, "Wavevector k must be non-negative."))
    end

    if ϕ ≈ zero(T)
        return one(T)
    end

    coeffs = _yukawa_msa_coefficients(ϕ, K, z)
    Pk = _P_factor_yukawa_msa(k, ϕ, z, coeffs.a, coeffs.b, coeffs.c, coeffs.d)
    return abs2(Pk)
end

@doc """
    S_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}

Calculates the exact analytical static structure factor ``S(k)`` for a monodisperse Yukawa fluid
under the Mean Spherical Approximation (MSA).

``S(k) = \\frac{1}{|P(k)|^2}``

# Arguments
- `ϕ::T`: Volume fraction.
- `K::T`: Contact energy parameter (dimensionless).
- `z::T`: Inverse screening length (dimensionless).
- `k::T`: Dimensionless wavevector ``k = q\\sigma``.

# Returns
- `T`: Static structure factor ``S(k)``.
"""
function S_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}
    inv_S = IS_Yukawa_MSA(ϕ, K, z, k)
    return one(T) / inv_S
end

@doc """
    C_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}

Calculates the Fourier transform of the direct correlation function ``C(k)`` for the Yukawa MSA model:
``C(k) = \\frac{1 - 1/S(k)}{24 \\phi}``
"""
function C_Yukawa_MSA(ϕ::T, K::T, z::T, k::T) where {T<:AbstractFloat}
    if ϕ ≈ zero(T)
        return C_HS_PY(zero(T), k)
    end
    inv_S = IS_Yukawa_MSA(ϕ, K, z, k)
    return (one(T) - inv_S) / (T(24.0) * ϕ)
end

# ------------------------------------------------------------------------------
# CLOSURE FACTORY AND VECTORIZED METHODS
# ------------------------------------------------------------------------------

"""
    S_Yukawa_MSA(ϕ::T, K::T, z::T) where {T<:AbstractFloat}

Factory function returning a closure `k -> S_Yukawa_MSA(ϕ, K, z, k)`.
Precomputes the polynomial root coefficients `(a, b, c, d)` for maximum performance.
"""
function S_Yukawa_MSA(ϕ::T, K::T, z::T) where {T<:AbstractFloat}
    coeffs = _yukawa_msa_coefficients(ϕ, K, z)
    return k -> begin
        Pk = _P_factor_yukawa_msa(k, ϕ, z, coeffs.a, coeffs.b, coeffs.c, coeffs.d)
        return one(T) / abs2(Pk)
    end
end

"""
    S_Yukawa_MSA(ϕ::T, K::T, z::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Vectorized evaluation of `S_Yukawa_MSA`.
"""
function S_Yukawa_MSA(ϕ::T, K::T, z::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    f = S_Yukawa_MSA(ϕ, K, z)
    return f.(k_values)
end

"""
    IS_Yukawa_MSA(ϕ::T, K::T, z::T) where {T<:AbstractFloat}

Factory function returning a closure `k -> IS_Yukawa_MSA(ϕ, K, z, k)`.
"""
function IS_Yukawa_MSA(ϕ::T, K::T, z::T) where {T<:AbstractFloat}
    coeffs = _yukawa_msa_coefficients(ϕ, K, z)
    return k -> begin
        Pk = _P_factor_yukawa_msa(k, ϕ, z, coeffs.a, coeffs.b, coeffs.c, coeffs.d)
        return abs2(Pk)
    end
end

"""
    IS_Yukawa_MSA(ϕ::T, K::T, z::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Vectorized evaluation of `IS_Yukawa_MSA`.
"""
function IS_Yukawa_MSA(ϕ::T, K::T, z::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    f = IS_Yukawa_MSA(ϕ, K, z)
    return f.(k_values)
end
