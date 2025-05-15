# src/models/mixtures/vw_mixture.jl

# Este archivo contiene la lógica para la corrección de Verlet-Weiss (VW)
# aplicada a mezclas de esferas duras, utilizando la solución de Baxter como base.

# Se asume que S_HS_Baxter_mixture y/o IS_HS_Baxter_mixture están disponibles
# desde hardsphere_mixture.jl (incluido en el módulo principal).
# También se asume que phi_to_rho_mixture y rho_to_phi_mixture están disponibles
# desde utils.jl o un archivo similar.

@doc """
    VW_correction_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the Verlet-Weiss effective number densities and wavevector scaling factor
for a multi-component hard-sphere mixture.

This function implements the logic from the user's `VW_Baxter.m` script.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters for each component.
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions for each component.

# Returns
- `Tuple{Vector{T}, T}`: A tuple `(ρ_eff_vector, λk_eff_factor)` where:
    - `ρ_eff_vector`: Vector of effective number densities for each component.
    - `λk_eff_factor`: Effective wavevector scaling factor.
"""
function VW_correction_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}) where {T<:AbstractFloat}
    if length(ϕ_vector) != length(σ_vector)
        throw(DimensionMismatch("ϕ_vector and σ_vector must have the same length"))
    end

    # Calculate original number densities and mole fractions
    # Assuming phi_to_rho_mixture is available from utils.jl
    ρ_orig_vector = phi_to_rho_mixture(ϕ_vector, σ_vector)
    
    sum_ρ_orig = sum(ρ_orig_vector)
    if sum_ρ_orig ≈ zero(T) # Avoid division by zero if all densities are zero
        # If no particles, no correction needed, return original densities (all zero) and lambda_k=1
        return ρ_orig_vector, one(T)
    end
    ξ_vector = ρ_orig_vector ./ sum_ρ_orig # Mole fractions x_i

    # Total packing fraction
    ϕ_total = sum(ϕ_vector)

    # Verlet-Weiss correction for total packing fraction
    ϕ_total_eff = ϕ_total * (one(T) - ϕ_total / T(16.0))
    if ϕ_total_eff < zero(T)
        @warn "Effective total packing fraction ϕ_total_eff = $ϕ_total_eff is negative. Check input parameters."
        # Decide how to handle: clamp, error, or allow? For now, allow.
    end
    
    # Wavevector scaling factor
    # Avoid issues if ϕ_total is zero (implies ϕ_total_eff is also zero)
    λk_eff_factor = (ϕ_total ≈ zero(T)) ? one(T) : cbrt(ϕ_total_eff / ϕ_total)
    if ϕ_total ≈ zero(T) && !(ϕ_total_eff ≈ zero(T)) # Should not happen if ϕ_total_eff formula is correct
         @warn "ϕ_total is zero but ϕ_total_eff is not. λk_eff_factor might be ill-defined."
    end


    # Effective total number density
    # <σ³> = Σᵢ xᵢ σᵢ³
    mean_sigma_cubed = sum(ξ_vector .* (σ_vector.^3))
    if mean_sigma_cubed ≈ zero(T) && !(ϕ_total_eff ≈ zero(T))
        @warn "<σ³> is zero but ϕ_total_eff is not. Effective total density will be Inf/NaN."
        # This case implies all sigma are zero, but there's a packing fraction.
        # Or all xi are zero for non-zero sigma, which is impossible if sum_ρ_orig > 0.
        # If all sigma are zero, then ϕ_total should be zero too.
        ρ_total_eff = T(Inf) # Or handle as an error
    elseif mean_sigma_cubed ≈ zero(T) && ϕ_total_eff ≈ zero(T)
        ρ_total_eff = zero(T)
    else
        ρ_total_eff = (T(6.0) / T(π)) * ϕ_total_eff / mean_sigma_cubed
    end

    # Effective number densities for each component
    ρ_eff_vector = ξ_vector .* ρ_total_eff

    return ρ_eff_vector, λk_eff_factor
end


@doc """
    IS_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates the inverse of the partial structure factor matrix, S⁻¹(k), for a
multi-component hard-sphere mixture using the Percus-Yevick approximation
with Verlet-Weiss corrections.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of original hard-sphere diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of original packing fractions.
- `k_wavevector::T`: Magnitude of the wavevector.

# Returns
- `Matrix{Complex{T}}`: The n x n inverse structure factor matrix S⁻¹_ij(k) with VW corrections.
"""
function IS_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}
    ρ_eff_vector, λk_eff_factor = VW_correction_mixture(σ_vector, ϕ_vector)
    
    # Call the Baxter solution with effective parameters
    # IS_HS_Baxter_mixture should be available from hardsphere_mixture.jl
    return IS_HS_Baxter_mixture(σ_vector, ρ_eff_vector, k_wavevector * λk_eff_factor)
end

@doc """
    S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates the partial structure factor matrix S(k) for a multi-component
hard-sphere mixture using the Percus-Yevick approximation with Verlet-Weiss (VW) corrections.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of original hard-sphere diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of original packing fractions.
- `k_wavevector::T`: Magnitude of the wavevector.

# Returns
- `Matrix{T}`: The n x n real-valued partial structure factor matrix S_ij(k) with VW corrections.
"""
function S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}
    # S_HS_Baxter_mixture should be available from hardsphere_mixture.jl
    ρ_eff_vector, λk_eff_factor = VW_correction_mixture(σ_vector, ϕ_vector)
    
    # Call the Baxter solution with effective parameters
    return S_HS_Baxter_mixture(σ_vector, ρ_eff_vector, k_wavevector * λk_eff_factor)
    # Alternatively, using the IS_HS_VW_mixture:
    # inv_S_matrix_vw = IS_HS_VW_mixture(σ_vector, ϕ_vector, k_wavevector)
    # S_matrix_complex_vw = inv(inv_S_matrix_vw)
    # if any(abs.(imag.(S_matrix_complex_vw)) .> sqrt(eps(T)))
    #     @warn "Significant imaginary part found in S_HS_VW_mixture for k=$k_wavevector."
    # end
    # return real.(S_matrix_complex_vw)
end

# --- Versión "Fábrica" (Opcional) ---
@doc """
    S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}) where {T<:AbstractFloat}

Returns a function `f(k)` that calculates the VW-corrected partial structure factor matrix S(k)
for a multi-component hard-sphere mixture. System parameters are fixed.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of original hard-sphere diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of original packing fractions.

# Returns
- `Function`: A function that takes a wavevector `k` and returns the S_ij(k) matrix.
"""
function S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}) where {T<:AbstractFloat}
    # Pre-calculate effective parameters to capture in the closure
    ρ_eff_vector, λk_eff_factor = VW_correction_mixture(σ_vector, ϕ_vector)
    
    function structure_factor_k_matrix_vw(k_val::T)
        # S_HS_Baxter_mixture should be available
        return S_HS_Baxter_mixture(σ_vector, ρ_eff_vector, k_val * λk_eff_factor)
    end
    return structure_factor_k_matrix_vw
end

# --- Versión Vectorizada (Opcional) ---
@doc """
    S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the VW-corrected partial structure factor matrix S_ij(k) for a vector of k values.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of original diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of original packing fractions.
- `k_values::AbstractVector{T}`: A vector of wavevector magnitudes.

# Returns
- `Vector{Matrix{T}}`: A vector of S_ij(k) matrices, one for each k in `k_values`.
"""
function S_HS_VW_mixture(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    sk_function_vw = S_HS_VW_mixture(σ_vector, ϕ_vector) # Get the factory function
    return [sk_function_vw(k_val) for k_val in k_values]
end

