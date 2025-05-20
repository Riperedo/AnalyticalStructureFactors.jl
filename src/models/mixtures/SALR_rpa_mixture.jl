# src/models/mixtures/SALR_rpa_mixture.jl

# Este archivo contiene la lógica para la Random Phase Approximation (RPA)
# aplicada a sistemas multicomponente con interacciones SALR (Short-Range Attraction, Long-Range Repulsion).

# Se asume que IS_HS_VW_mixture está disponible desde vw_mixture.jl.
# Se asume que betaU_Yukawa está disponible desde potentials.jl.
# Se asume que phi_to_rho_mixture está disponible desde utils.jl.

@doc """
    S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                       k_wavevector::T, 
                       A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                       A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}

Calculates the partial structure factor matrix S_ij(k) for a multi-component mixture
using the Random Phase Approximation (RPA). The reference system is a hard-sphere
mixture with Verlet-Weiss corrections. The perturbation potential for each pair (i,j)
is assumed to be of SALR form, modeled as u_SALR = u_Yukawa_attractive - u_Yukawa_repulsive.

The RPA formula used is:
S⁻¹ = (S⁰)⁻¹ - βŨ_scaled
where βŨ_scaled_ij(k) = √(ρᵢρⱼ) * βũ_SALR_ij(k).

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters for each component (reference system).
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions for each component (reference system).
- `k_wavevector::T`: Magnitude of the wavevector.
- `A_attr_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aᵃᵢⱼ for the attractive Yukawa component.
- `Z_attr_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zᵃᵢⱼ for the attractive Yukawa component.
- `A_rep_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aʳᵢⱼ for the repulsive Yukawa component.
- `Z_rep_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zʳᵢⱼ for the repulsive Yukawa component.

# Returns
- `Matrix{T}`: The n x n real-valued partial structure factor matrix S_ij(k).
"""
function S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                            k_wavevector::T, 
                            A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                            A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}
    
    n_components = length(ϕ_vector)
    if length(σ_vector) != n_components || 
       size(A_attr_matrix) != (n_components, n_components) || 
       size(Z_attr_matrix) != (n_components, n_components) ||
       size(A_rep_matrix) != (n_components, n_components) ||
       size(Z_rep_matrix) != (n_components, n_components)
        throw(DimensionMismatch("Input vector/matrix dimensions are inconsistent."))
    end

    # 1. Get the inverse structure factor matrix of the reference system (HS-VW)
    inv_S0_matrix = IS_HS_VW_mixture(σ_vector, ϕ_vector, k_wavevector)

    # 2. Calculate number densities for scaling the perturbation
    ρ_vector = phi_to_rho_mixture(ϕ_vector, σ_vector)

    # 3. Construct the matrix of scaled Fourier-transformed perturbation potentials
    # βŨ_scaled_ij(k) = √(ρᵢρⱼ) * ( βũ_Yukawa_attr_ij(k) - βũ_Yukawa_rep_ij(k) )
    
    beta_U_pert_scaled_matrix = zeros(Complex{T}, n_components, n_components)
    for i in 1:n_components
        for j in 1:n_components
            # Parameters for the attractive Yukawa part
            A_attr_ij = A_attr_matrix[i,j]
            Z_attr_ij = Z_attr_matrix[i,j]
            
            # Parameters for the repulsive Yukawa part
            A_rep_ij = A_rep_matrix[i,j]
            Z_rep_ij = Z_rep_matrix[i,j]
            
            # Calculate Fourier transforms for attractive and repulsive parts
            beta_u_tilde_attr_ij = betaU_Yukawa(A_attr_ij, Z_attr_ij, k_wavevector)
            beta_u_tilde_rep_ij  = betaU_Yukawa(A_rep_ij, Z_rep_ij, k_wavevector)
            
            # Combine for SALR perturbation
            beta_u_tilde_salr_ij = beta_u_tilde_attr_ij - beta_u_tilde_rep_ij
            
            # Apply density scaling
            beta_U_pert_scaled_matrix[i,j] = sqrt(ρ_vector[i] * ρ_vector[j]) * beta_u_tilde_salr_ij
        end
    end

    # 4. Apply the RPA formula: S⁻¹ = (S⁰)⁻¹ + βŨ_scaled
    inv_S_matrix_RPA = inv_S0_matrix + beta_U_pert_scaled_matrix
    
    # 5. Invert to get S(k)
    S_matrix_RPA_complex = inv(inv_S_matrix_RPA)

    # The result S_ij(k) should be real. Small imaginary parts can arise from numerical precision.
    if any(abs.(imag.(S_matrix_RPA_complex)) .> sqrt(eps(T)))
        @warn "Significant imaginary part found in S_RPA_mixture_SALR for k=$k_wavevector. Max imag: $(maximum(abs.(imag.(S_matrix_RPA_complex))))."
    end
    
    return real.(S_matrix_RPA_complex)
end


# --- Versión "Fábrica" (Opcional) ---
@doc """
    S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                       A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                       A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}

Returns a function `f(k)` that calculates the RPA partial structure factor matrix S_ij(k)
for a multi-component SALR mixture. System and potential parameters are fixed.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions.
- `A_attr_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aᵃᵢⱼ for the attractive Yukawa component.
- `Z_attr_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zᵃᵢⱼ for the attractive Yukawa component.
- `A_rep_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aʳᵢⱼ for the repulsive Yukawa component.
- `Z_rep_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zʳᵢⱼ for the repulsive Yukawa component.

# Returns
- `Function`: A function that takes a wavevector `k` and returns the S_ij(k) matrix.
"""
function S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                            A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                            A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}
    
    function structure_factor_k_matrix_rpa_salr(k_val::T)
        return S_RPA_mixture_SALR(σ_vector, ϕ_vector, k_val, 
                                  A_attr_matrix, Z_attr_matrix,
                                  A_rep_matrix, Z_rep_matrix)
    end
    return structure_factor_k_matrix_rpa_salr
end

# --- Versión Vectorizada (Opcional) ---
@doc """
    S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                       A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                       A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}, 
                       k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the RPA partial structure factor matrix S_ij(k) for a vector of k values
for a multi-component SALR mixture.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions.
- `A_attr_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aᵃᵢⱼ for the attractive Yukawa component.
- `Z_attr_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zᵃᵢⱼ for the attractive Yukawa component.
- `A_rep_matrix::AbstractMatrix{T}`: Matrix of amplitudes Aʳᵢⱼ for the repulsive Yukawa component.
- `Z_rep_matrix::AbstractMatrix{T}`: Matrix of inverse screening lengths zʳᵢⱼ for the repulsive Yukawa component.
- `k_values::AbstractVector{T}`: A vector of wavevector magnitudes.

# Returns
- `Vector{Matrix{T}}`: A vector of S_ij(k) matrices, one for each k in `k_values`.
"""
function S_RPA_mixture_SALR(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                            A_attr_matrix::AbstractMatrix{T}, Z_attr_matrix::AbstractMatrix{T},
                            A_rep_matrix::AbstractMatrix{T}, Z_rep_matrix::AbstractMatrix{T}, 
                            k_values::AbstractVector{T}) where {T<:AbstractFloat}
    
    sk_function_rpa_salr = S_RPA_mixture_SALR(σ_vector, ϕ_vector, 
                                             A_attr_matrix, Z_attr_matrix,
                                             A_rep_matrix, Z_rep_matrix)
    return [sk_function_rpa_salr(k_val) for k_val in k_values]
end

