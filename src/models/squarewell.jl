# src/models/squarewell.jl

@doc """
    S_SW_RPA(ϕ::T, temperature::T, lambda_range::T, k::T) where {T<:AbstractFloat}

Calcula el factor de estructura estático para un fluido de pozo cuadrado
utilizando la Random Phase Approximation (RPA).
El sistema de referencia es de esferas duras con corrección de Verlet-Weiss.

La fórmula general de RPA es: S(k) = 1 / ( 1/S₀(k) - ρ * β * u_pert_tilde(k) )
donde S₀(k) es el factor de estructura del sistema de referencia,
ρ es la densidad numérica, y u_pert_tilde(k) es la transformada de Fourier
de la parte perturbativa del potencial.

# Arguments
- `ϕ::T`: Fracción de volumen de las esferas.
- `temperature::T`: Temperatura adimensional (ej. k_B * T_abs / ε).
- `lambda_range::T`: Rango del pozo cuadrado (λ en unidades de σ).
- `k::T`: Vector de onda adimensional (k = qσ).

# Returns
- `T`: Valor del factor de estructura S(k).
"""
function S_SW_RPA(ϕ::T, temperature::T, lambda_range::T, k::T) where {T<:AbstractFloat}
    # Factor de estructura del sistema de referencia (HS-VW)
    s0_k = S_HS_VW(ϕ, k) # S_HS_VW ya usa phi_VW y k_VW internamente si es necesario por su propia definición.
                         # O podrías usar IS_HS_VW(ϕ, k) para obtener 1/S₀(k) directamente.
    
    inv_s0_k = IS_HS_VW(ϕ, k) # Más directo: 1 / S₀(k)

    beta_u_pert_k = betaU_SW(temperature, lambda_range, k)

    
    density_factor = (T(6.0) * ϕ) / T(π)
    
    inv_s_k = inv_s0_k - density_factor * beta_u_pert_k
    
    if inv_s_k <= zero(T)
        @warn "Inverse structure factor 1/S(k) = $inv_s_k for SW-RPA is non-positive. Resulting S(k) will be unphysical."
    end

    return one(T) / inv_s_k
end
