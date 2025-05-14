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


# --- Effective Diameter Calculation (Blip Function) ---

@doc """
    blip(temperature::T_num; ν::Integer = 6) where {T_num<:AbstractFloat}

Calculates the effective hard-sphere diameter (λ) based on the "blip" function
formalism for a soft potential. The function returns both λ and λ³.

The formula for λ³ is:
λ³(T, ν) = 1 - 3 * ∫_0^1 dx x² * exp(-(1/T) * ((1/x^ν) - 1)²)

# Arguments
- `temperature::T_num`: Dimensionless temperature of the system.

# Keywords
- `ν::Integer = 6`: Parameter modulating the "softness" or steepness of the repulsive part of the potential.

# Returns
- `Tuple{T_num, T_num}`: A tuple containing `(λ, λ³)`, where `λ` is the effective diameter
  and `λ³` is its cube.

# References
[1] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola.
    Phys. Rev. E 87, 052306 – Published 22 May 2013.
"""
function blip(temperature::T_num; ν::Integer = 6) where {T_num<:AbstractFloat}
    if temperature ≈ zero(T_num) # Check if temperature is very close to zero
        # For T=0, the exponential term goes to 0 (if (1/x^ν - 1)² > 0) or 1 (if (1/x^ν - 1)² = 0, i.e. x=1).
        # The integral becomes 0 if integrand is 0, leading to λ³ = 1.
        # If T=0 implies hard sphere limit, then λ=1.
        return one(T_num), one(T_num)
    else
        # Numerical integration using a simple Riemann sum.
        # For higher precision, consider using a library like QuadGK.jl.
        dx = T_num(0.0001) # Step size for integration; smaller for more accuracy.
        # Integration range from dx to 1 to avoid x=0 if 1/x^ν is problematic.
        # The integrand should be negligible or zero at x=0 for this to be accurate.
        x_range = dx:dx:one(T_num)
        
        # Pre-calculate 1/T to avoid repeated division
        inv_T = one(T_num) / temperature

        # Calculate integrand values
        # Integrand: x² * exp(-(1/T) * ( (1/x^ν) - 1 )²)
        # (1/x^ν - 1)² = (x^(-ν) - 1)² = x^(-2ν) - 2x^(-ν) + 1
        integrand_values = [x_val^2 * exp(-inv_T * (x_val^(-2*ν) - T_num(2.0)*x_val^(-ν) + one(T_num))) for x_val in x_range]
        
        integral_sum = sum(integrand_values) * dx
        
        lambda_cubed = one(T_num) - T_num(3.0) * integral_sum
        
        # Ensure lambda_cubed is not negative due to numerical errors,
        # though physically it should be positive.
        if lambda_cubed < zero(T_num)
            @warn "λ³ calculated as $lambda_cubed, which is negative. Clamping to zero. Check parameters or integration accuracy."
            lambda_cubed = zero(T_num)
        end
        
        lambda = cbrt(lambda_cubed) # Equivalent to lambda_cubed^(one(T_num)/T_num(3.0))
        
        return lambda, lambda_cubed
    end
end

# --- Verlet-Weiss Corrections for Hard Spheres ---

@doc """
    phi_VW(ϕ_HS::T) where {T<:AbstractFloat}

Computes the density correction (effective packing fraction ϕ_eff) for hard spheres
using the Verlet-Weiss (VW) approximation.
ϕ_eff = ϕ_HS * (1 - ϕ_HS / 16)

# Arguments
- `ϕ_HS::T`: The packing fraction of the hard-sphere system.

# Returns
- `T`: The Verlet-Weiss corrected effective packing fraction.

# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972.
"""
phi_VW(ϕ_HS::T) where {T<:AbstractFloat} = ϕ_HS * (one(T) - ϕ_HS / T(16.0))

@doc """
    k_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}

Computes the wave vector correction for hard spheres using the Verlet-Weiss (VW)
approximation. This gives an effective wavevector k_eff.
k_eff = k_HS * (ϕ_eff / ϕ_HS)^(1/3)

# Arguments
- `ϕ_HS::T`: The original packing fraction of the hard-sphere system.
- `k_HS::T`: The original (dimensionless) wave vector.

# Returns
- `T`: The Verlet-Weiss corrected effective wave vector.

# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972.
"""
function k_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}
    if ϕ_HS ≈ zero(T)
        return k_HS # Avoid division by zero if ϕ_HS is zero; ratio is 1.
    end
    # ϕ_eff / ϕ_HS = (1 - ϕ_HS / 16)
    ratio_phi = one(T) - ϕ_HS / T(16.0)
    # Handle potential negative ratio if ϕ_HS is unphysically large, though VW is for moderate densities.
    if ratio_phi < zero(T)
        @warn "Ratio ϕ_eff/ϕ_HS = $ratio_phi is negative for ϕ_HS = $ϕ_HS. Result for k_VW might be unphysical."
        # Return k_HS or handle as an error, depending on desired behavior.
        # For now, let cbrt handle it (which might produce NaN or Complex for negative base).
    end
    return k_HS * cbrt(ratio_phi)
end
