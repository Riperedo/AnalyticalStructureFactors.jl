# src/models/mixtures/squarewell_rpa_mixture.jl

# Este archivo contiene la lógica para la Random Phase Approximation (RPA)
# aplicada a sistemas multicomponente con interacciones de Pozo Cuadrado (Square-Well).

# Se asume que IS_HS_VW_mixture está disponible desde vw_mixture.jl.
# Se asume que betaU_SW está disponible desde potentials.jl.
# Se asume que phi_to_rho_mixture está disponible desde utils.jl.

@doc """
    S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                             k_wavevector::T, 
                             temperature_matrix::AbstractMatrix{T}, 
                             lambda_range_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}

Calculates the partial structure factor matrix S_ij(k) for a multi-component mixture
using the Random Phase Approximation (RPA). The reference system is a hard-sphere
mixture with Verlet-Weiss corrections. The perturbation potential for each pair (i,j)
is assumed to be of Square-Well form.

The RPA formula used is:
S⁻¹ = (S⁰)⁻¹ - βŨ_scaled
where βŨ_scaled_ij(k) = √(ρᵢρⱼ) * βũ_SquareWell_ij(k).

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters for each component (reference system).
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions for each component (reference system).
- `k_wavevector::T`: Magnitude of the wavevector.
- `temperature_matrix::AbstractMatrix{T}`: Matrix of dimensionless temperatures Tᵢⱼ for each pair (i,j).
                                          Tᵢⱼ is typically related to the well depth εᵢⱼ (e.g., k_B * T_abs / εᵢⱼ).
- `lambda_range_matrix::AbstractMatrix{T}`: Matrix of Square-Well ranges λᵢⱼ for each pair (i,j)
                                            (in units of the hard-sphere diameter σ of a reference component, or σᵢⱼ).

# Returns
- `Matrix{T}`: The n x n real-valued partial structure factor matrix S_ij(k).
"""
function S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                                  k_wavevector::T, 
                                  temperature_matrix::AbstractMatrix{T}, 
                                  lambda_range_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}
    
    n_components = length(ϕ_vector)
    if length(σ_vector) != n_components || 
       size(temperature_matrix) != (n_components, n_components) || 
       size(lambda_range_matrix) != (n_components, n_components)
        throw(DimensionMismatch("Input vector/matrix dimensions are inconsistent."))
    end

    # 1. Get the inverse structure factor matrix of the reference system (HS-VW)
    inv_S0_matrix = IS_HS_VW_mixture(σ_vector, ϕ_vector, k_wavevector)

    # 2. Calculate number densities for scaling the perturbation
    ρ_vector = phi_to_rho_mixture(ϕ_vector, σ_vector)

    # 3. Construct the matrix of scaled Fourier-transformed perturbation potentials
    # βŨ_scaled_ij(k) = √(ρᵢρⱼ) * βũ_SquareWell_ij(k)
    
    beta_U_pert_scaled_matrix = zeros(Complex{T}, n_components, n_components)
    for i in 1:n_components
        for j in 1:n_components
            # Get parameters for the ij pair
            temp_ij = temperature_matrix[i,j]
            lambda_ij = lambda_range_matrix[i,j]
            
            # Calculate βũ_SquareWell_ij(k)
            # Asumimos que betaU_SW toma la temperatura y el rango para el par específico.
            # El 'k' aquí es el k_wavevector global.
            beta_u_tilde_ij = betaU_SW(temp_ij, lambda_ij, k_wavevector)
            
            # Apply density scaling
            beta_U_pert_scaled_matrix[i,j] = sqrt(ρ_vector[i] * ρ_vector[j]) * beta_u_tilde_ij
        end
    end

    # 4. Apply the RPA formula: S⁻¹ = (S⁰)⁻¹ - βŨ_scaled
    inv_S_matrix_RPA = inv_S0_matrix - beta_U_pert_scaled_matrix
    
    # 5. Invert to get S(k)
    S_matrix_RPA_complex = inv(inv_S_matrix_RPA)

    # The result S_ij(k) should be real. Small imaginary parts can arise from numerical precision.
    if any(abs.(imag.(S_matrix_RPA_complex)) .> sqrt(eps(T)))
        @warn "Significant imaginary part found in S_RPA_mixture_SquareWell for k=$k_wavevector. Max imag: $(maximum(abs.(imag.(S_matrix_RPA_complex))))."
    end
    
    return real.(S_matrix_RPA_complex)
end


# --- Versión "Fábrica" (Opcional) ---
@doc """
    S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                             temperature_matrix::AbstractMatrix{T}, 
                             lambda_range_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}

Returns a function `f(k)` that calculates the RPA partial structure factor matrix S_ij(k)
for a multi-component Square-Well mixture. System and potential parameters are fixed.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions.
- `temperature_matrix::AbstractMatrix{T}`: Matrix of dimensionless temperatures Tᵢⱼ.
- `lambda_range_matrix::AbstractMatrix{T}`: Matrix of Square-Well ranges λᵢⱼ.

# Returns
- `Function`: A function that takes a wavevector `k` and returns the S_ij(k) matrix.
"""
function S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                                  temperature_matrix::AbstractMatrix{T}, 
                                  lambda_range_matrix::AbstractMatrix{T}) where {T<:AbstractFloat}
    
    function structure_factor_k_matrix_rpa_sw(k_val::T)
        return S_RPA_mixture_SquareWell(σ_vector, ϕ_vector, k_val, temperature_matrix, lambda_range_matrix)
    end
    return structure_factor_k_matrix_rpa_sw
end

# --- Versión Vectorizada (Opcional) ---
@doc """
    S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                             temperature_matrix::AbstractMatrix{T}, 
                             lambda_range_matrix::AbstractMatrix{T}, 
                             k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the RPA partial structure factor matrix S_ij(k) for a vector of k values
for a multi-component Square-Well mixture.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of diameters.
- `ϕ_vector::AbstractVector{T}`: Vector of packing fractions.
- `temperature_matrix::AbstractMatrix{T}`: Matrix of dimensionless temperatures Tᵢⱼ.
- `lambda_range_matrix::AbstractMatrix{T}`: Matrix of Square-Well ranges λᵢⱼ.
- `k_values::AbstractVector{T}`: A vector of wavevector magnitudes.

# Returns
- `Vector{Matrix{T}}`: A vector of S_ij(k) matrices, one for each k in `k_values`.
"""
function S_RPA_mixture_SquareWell(σ_vector::AbstractVector{T}, ϕ_vector::AbstractVector{T}, 
                                  temperature_matrix::AbstractMatrix{T}, 
                                  lambda_range_matrix::AbstractMatrix{T}, 
                                  k_values::AbstractVector{T}) where {T<:AbstractFloat}
    
    sk_function_rpa_sw = S_RPA_mixture_SquareWell(σ_vector, ϕ_vector, temperature_matrix, lambda_range_matrix)
    return [sk_function_rpa_sw(k_val) for k_val in k_values]
end

