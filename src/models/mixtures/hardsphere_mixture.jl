# src/models/mixtures/hardsphere_mixture.jl

# Contiene la lógica para la solución de Baxter de la
# aproximación de Percus-Yevick para mezclas de esferas duras.

@doc """
    Qk_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k_wavevector::T) where {T<:AbstractFloat}

Calculates Baxter's Q-matrix for a multi-component hard-sphere mixture
within the Percus-Yevick approximation.

The formulas for `pR` and `pS` are based on the user's provided  implementation.

# Arguments
- `σ_vector::AbstractVector{T}`: Vector of hard-sphere diameters for each component.
- `ρ_vector::AbstractVector{T}`: Vector of number densities for each component.
- `k_wavevector::T`: Magnitude of the wavevector.

# Returns
- `Matrix{Complex{T}}`: The n x n Baxter Q-matrix.
"""
function Qk_mixture(σ_vector::AbstractVector{T}, ρ_vector::AbstractVector{T}, k::T) where {T<:AbstractFloat}
    n_components = length(ρ_vector)
    if length(σ_vector) != n_components
        throw(DimensionMismatch("σ_vector and ρ_vector must have the same length"))
    end

    # Corregido: Asegurar multiplicación por π/6
    ζ2 = sum(ρ_vector .* (σ_vector.^2)) * (T(π)/T(6.0))
    ζ3 = sum(ρ_vector .* (σ_vector.^3)) * (T(π)/T(6.0)) # Total packing fraction ϕ_total

    Δ = one(T) - ζ3
    if Δ <= zero(T) # Podría ser problemático si Δ es exactamente cero también.
        @warn "Denominator (1 - ζ3) = $Δ is non-positive or zero. Results may be unphysical."
    end
    Δ² = Δ^2 # Evitar Δ = 0 aquí. Si Δ es cero, Δ² es cero.

    # Generalizado a n componentes
    Q_matrix = zeros(Complex{T}, n_components, n_components)

    for i in 1:n_components
        σ_i = σ_vector[i]
        ρ_i = ρ_vector[i]
        
        # Coeficientes a_i y b_i según la forma del código 
        # Estos son parámetros para la función Q_ij(r) de Baxter.
        # Manejar el caso Δ² ≈ 0 para evitar división por cero si Δ es muy pequeño.
        if Δ² ≈ zero(T)
            # Asignar valores que eviten NaN/Inf o manejar como error.
            # Esto indica una densidad muy alta, probablemente no física para PY.
            a_i = T(Inf) # O algún otro valor de señalización
            b_i = T(Inf)
            @warn "Δ² is close to zero in Qk_mixture. Coefficients a_i, b_i might be Inf."
        else
            a_i = (one(T) - ζ3 + T(3.0)*σ_i*ζ2)/(Δ²)
            b_i = -(T(3.0)/T(2.0))*(σ_i^2)*ζ2/(Δ²)
        end


        for j in 1:n_components
            σ_j = σ_vector[j]
            ρ_j = ρ_vector[j]

            delta_ij = (i == j ? one(T) : zero(T))

            # Límites usados en las fórmulas de  para pR y pS
            Sr_limit =(σ_i - σ_j)/T(2.0)
            R_limit = (σ_i + σ_j)/T(2.0)
            
            if k ≈ zero(T)
                # Límite k=0 para la integral de Q_ij(r).
                # Asumiendo Q_ij(r) = a_i(r-R_limit) + 0.5*b_i(r-R_limit)^2, integrada de Sr_limit a R_limit.
                # Esta forma de Q_ij(r) es una interpretación común que podría llevar a las expresiones pR, pS.
                # Integral = -a_i * (σ_j^2) / 8 + b_i * (σ_j^3) / 48
                # (Esta fórmula para la integral debe ser consistente con la forma de Q_ij(r)
                # cuya transformada de Fourier da las expresiones pR y pS)
                integral_Qij_r_k0 = -a_i * (σ_j^2) / T(8.0) + b_i * (σ_j^3) / T(48.0)
                
                # El factor 2π está en la fórmula final de Q_matrix[i,j]
                Q_matrix[i,j] = delta_ij - T(2.0)*T(π)*sqrt(ρ_i*ρ_j) * integral_Qij_r_k0
            else
                # Usando las fórmulas exactas del código  proporcionado por el usuario para pR y pS
                ik = im * k
                ik2 = ik^2
                ik3 = ik^3

                # pR (límite superior R_limit)
                exp_ikR = exp(ik * R_limit)
                # El 'R' en bi*R*... se interpreta como R_limit aquí
                pR = (a_i/T(2.0))*exp_ikR*((T(2.0)-T(2.0)*ik*R_limit-(k^2)*(R_limit^2))/(ik3)) -
                     (a_i/T(2.0))*(R_limit^2)*(exp_ikR/(ik)) +
                     b_i*exp_ikR*((ik*R_limit - one(T))/(ik2)) -
                     b_i*R_limit*(exp_ikR/(ik))
            
                # pS (límite superior Sr_limit)
                exp_ikSr = exp(ik * Sr_limit)
                # El 'R' en bi*R*... se interpreta como R_limit aquí, como en el código 
                pS = (a_i/T(2.0))*exp_ikSr*((T(2.0)-T(2.0)*ik*Sr_limit-(k^2)*(Sr_limit^2))/(ik3)) -
                     (a_i/T(2.0))*(R_limit^2)*(exp_ikSr/(ik)) + # Error en : Debería ser Sr_limit^2 aquí? Usando R_limit^2 como en .
                     b_i*exp_ikSr*((ik*Sr_limit - one(T))/(ik2)) -
                     b_i*R_limit*(exp_ikSr/(ik)) # Usando R_limit como en el código 

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

