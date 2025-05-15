# src/models/mixtures/hardsphere_mixture.jl

# Contiene la lógica para la solución de Baxter de la
# aproximación de Percus-Yevick para mezclas de esferas duras.

@doc """
    Qk_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates Baxter's Q-matrix for a multi-component hard-sphere mixture
within the Percus-Yevick approximation.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters for each component.
- `ρ_vector::AbstractVector{T}`: Vector of number densities for each component.
- `k_wavevector::T`: Magnitude of the wavevector.

# Returns
- `Matrix{Complex{T}}`: The n x n Baxter Q-matrix.
"""
function Qk_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}
    n_components = length(ρ_vector)
    if length(σ_vector) != n_components
        throw(DimensionMismatch("σ_vector and ρ_vector must have the same length"))
    end

    ζ2 = sum(ρ_vector .* (σ_vector.^2)) * (T(π)/T(6.0))
    ζ3 = sum(ρ_vector .* (σ_vector.^3)) * (T(π)/T(6.0)) # Total packing fraction ϕ_total

    Δ = one(T) - ζ3
    if Δ <= zero(T)
        @warn "Denominator (1 - ζ3) = $Δ is non-positive. Results may be unphysical."
    end
    Δ² = Δ^2

    Q_matrix = zeros(Complex{T}, n_components, n_components)
    k = k_wavevector # Alias for clarity

    for i in 1:n_components
        σ_i = σ_vector[i]
        ρ_i = ρ_vector[i]
        
        a_i = (one(T) - ζ3 + T(3.0)*σ_i*ζ2)/(Δ²)
        b_i = -(T(3.0)/T(2.0))*(σ_i^2)*ζ2/(Δ²)

        for j in 1:n_components
            σ_j = σ_vector[j]
            ρ_j = ρ_vector[j]

            delta_ij = (i == j ? one(T) : zero(T))

            R_limit = (σ_i + σ_j) / T(2.0)
            Sr_limit = (σ_i - σ_j) / T(2.0)

            if k ≈ zero(T)
                # Límite k=0 para la integral de Q_ij(r).
                # Q_ij(r) = a_i(r-R_limit) + 0.5*b_i(r-R_limit)^2, integrada de Sr_limit a R_limit.
                # Integral = -a_i * (σ_j^2) / 8 + b_i * (σ_j^3) / 48
                integral_Qij_r_k0 = -a_i * (σ_j^2) / T(8.0) + b_i * (σ_j^3) / T(48.0)
                Q_matrix[i,j] = delta_ij - T(2.0)*T(π)*sqrt(ρ_i*ρ_j) * integral_Qij_r_k0
            else
                ik = im * k
                ik2 = ik^2
                ik3 = ik^3

                # pR (límite superior R_limit)
                exp_ikR = exp(ik * R_limit)
                term_a1_R = (a_i / T(2.0)) * exp_ikR * (T(2.0) - T(2.0)*ik*R_limit - (k*R_limit)^2) / ik3
                term_a2_R = -(a_i / T(2.0)) * (R_limit^2) * exp_ikR / ik
                term_b1_R = b_i * exp_ikR * (ik*R_limit - one(T)) / ik2
                term_b2_R = -b_i * R_limit * exp_ikR / ik 
                pR = term_a1_R + term_a2_R + term_b1_R + term_b2_R

                # pS (límite superior Sr_limit)
                exp_ikSr = exp(ik * Sr_limit)
                term_a1_Sr = (a_i / T(2.0)) * exp_ikSr * (T(2.0) - T(2.0)*ik*Sr_limit - (k*Sr_limit)^2) / ik3
                term_a2_Sr = -(a_i / T(2.0)) * (Sr_limit^2) * exp_ikSr / ik # Este término usa Sr_limit^2
                term_b1_Sr = b_i * exp_ikSr * (ik*Sr_limit - one(T)) / ik2
                # CORRECCIÓN: El siguiente término usa R_limit, como en el código Matlab para pS
                term_b2_Sr = -b_i * R_limit * exp_ikSr / ik 
                pS = term_a1_Sr + term_a2_Sr + term_b1_Sr + term_b2_Sr
                
                Q_matrix[i,j] = delta_ij - T(2.0)*T(π)*sqrt(ρ_i*ρ_j)*(pR-pS)
            end
        end
    end
    return Q_matrix
end

@doc """
    IS_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates the inverse of the partial structure factor matrix, S⁻¹(k), for a
multi-component hard-sphere mixture using Baxter's Q-matrix formalism.
S⁻¹(k) = Qᴴ(k) * Q(k)
"""
function IS_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}
    Q = Qk_mixture(σ_vector, ρ_vector, k_wavevector)
    return Q'*Q 
end

@doc """
    S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates the partial structure factor matrix S(k) for a multi-component
hard-sphere mixture using Baxter's Q-matrix formalism (Percus-Yevick approximation).
S(k) = inv(Qᴴ(k) * Q(k))
"""
function S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}
    inv_S_matrix = IS_HS_Baxter_mixture(σ_vector, ρ_vector, k_wavevector)
    S_matrix_complex = inv(inv_S_matrix)
    
    max_imag_part = maximum(abs.(imag.(S_matrix_complex)))
    # Ajustar la tolerancia para la advertencia si es necesario, o hacerla relativa al valor real.
    # Por ejemplo, relativo al máximo elemento real de la diagonal.
    # CORRECCIÓN: Usar LinearAlgebra.diag
    diag_elements_real = real.(LinearAlgebra.diag(S_matrix_complex))
    diag_max_real = isempty(diag_elements_real) ? one(T) : max(one(T), maximum(abs.(diag_elements_real))) # Evitar error si diag_elements_real es vacío
    
    if max_imag_part > sqrt(eps(T)) * diag_max_real 
        @warn "Significant imaginary part found in S_HS_Baxter_mixture for k=$k_wavevector. Max imag part: $max_imag_part."
    end
    return real.(S_matrix_complex)
end

# --- Versión "Fábrica" (Opcional) ---
@doc """
    S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}) where {T<:AbstractFloat}

Returns a function `f(k)` that calculates the partial structure factor matrix S(k)
for a multi-component hard-sphere mixture.
"""
function S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}) where {T<:AbstractFloat}
    function structure_factor_k_matrix(k_val::T)
        return S_HS_Baxter_mixture(σ_vector, ρ_vector, k_val)
    end
    return structure_factor_k_matrix
end

# --- Versión Vectorizada (Opcional) ---
@doc """
    S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the partial structure factor matrix S_ij(k) for a vector of k values.
"""
function S_HS_Baxter_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    sk_function = S_HS_Baxter_mixture(σ_vector, ρ_vector)
    return [sk_function(k_val) for k_val in k_values]
end

