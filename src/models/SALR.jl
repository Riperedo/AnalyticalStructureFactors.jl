# src/models/SALR.jl

# Este archivo contendrá la lógica para factores de estructura
# basados en el potencial SALR (Short-Range Attraction, Long-Range Repulsion).

# --- Random Phase Approximation (RPA) para Potencial SALR ---
@doc """
    S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T, k::T) where {T<:AbstractFloat}

Calcula el factor de estructura estático para un fluido con potencial SALR
utilizando la Random Phase Approximation (RPA).
El sistema de referencia es de esferas duras con corrección de Verlet-Weiss.
El potencial SALR se modela como la diferencia de dos potenciales de Yukawa:
u_SALR(r) = u_Yukawa_attr(r) - u_Yukawa_rep(r).

# Arguments
- `ϕ::T`: Fracción de volumen de las esferas.
- `amplitude_attr::T`: Amplitud efectiva del componente atractivo de Yukawa (K_A).
- `inv_screening_length_attr::T`: Longitud de apantallamiento inversa del componente atractivo (z_A).
- `amplitude_rep::T`: Amplitud efectiva del componente repulsivo de Yukawa (K_R).
- `inv_screening_length_rep::T`: Longitud de apantallamiento inversa del componente repulsivo (z_R).
- `k::T`: Vector de onda adimensional (k = qσ).

# Returns
- `T`: Valor del factor de estructura S(k).
"""
function S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T, k::T) where {T<:AbstractFloat}
    inv_s0_k = IS_HS_VW(ϕ, k) 

    beta_u_attr_k = betaU_Yukawa(amplitude_attr, inv_screening_length_attr, k)
    beta_u_rep_k  = betaU_Yukawa(amplitude_rep, inv_screening_length_rep, k)
    beta_u_pert_k = beta_u_attr_k - beta_u_rep_k 

    density_factor = (T(6.0) * ϕ) / T(π)
    
    inv_s_k = inv_s0_k - density_factor * beta_u_pert_k
    
    # La siguiente advertencia ha sido eliminada según la solicitud del usuario.
    # if inv_s_k <= zero(T)
    #     @warn "Inverse structure factor 1/S(k) = $inv_s_k for SALR-RPA is non-positive. Resulting S(k) will be unphysical."
    # end

    return one(T) / inv_s_k
end

# --- Versión "Fábrica" (Opcional, pero sigue el patrón de tu proyecto anterior) ---
@doc """
    S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T) where {T<:AbstractFloat}

Retorna una función `f(k)` que calcula el factor de estructura estático para un fluido
con potencial SALR utilizando la RPA. Los parámetros del sistema y del potencial SALR se fijan.

# Arguments
- `ϕ::T`: Fracción de volumen de las esferas.
- `amplitude_attr::T`: Amplitud del componente atractivo de Yukawa.
- `inv_screening_length_attr::T`: Longitud de apantallamiento inversa del componente atractivo.
- `amplitude_rep::T`: Amplitud del componente repulsivo de Yukawa.
- `inv_screening_length_rep::T`: Longitud de apantallamiento inversa del componente repulsivo.

# Returns
- `Function`: Una función que toma un vector de onda `k` y retorna `S(k)`.
"""
function S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T) where {T<:AbstractFloat}
    function structure_factor_k(k_val::T)
        return S_SALR_RPA(ϕ, amplitude_attr, inv_screening_length_attr, amplitude_rep, inv_screening_length_rep, k_val)
    end
    return structure_factor_k 
end

# --- Versión Vectorizada (Opcional, para conveniencia) ---
@doc """
    S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calcula el factor de estructura S(k) para un vector de valores de k para el potencial SALR.

# Arguments
- `ϕ::T`: Fracción de volumen.
- `amplitude_attr::T`: Amplitud del componente atractivo de Yukawa.
- `inv_screening_length_attr::T`: Longitud de apantallamiento inversa del componente atractivo.
- `amplitude_rep::T`: Amplitud del componente repulsivo de Yukawa.
- `inv_screening_length_rep::T`: Longitud de apantallamiento inversa del componente repulsivo.
- `k_values::AbstractVector{T}`: Un vector de valores de k adimensionales.

# Returns
- `Vector{T}`: Un vector de valores de S(k) correspondientes a cada k en `k_values`.
"""
function S_SALR_RPA(ϕ::T, amplitude_attr::T, inv_screening_length_attr::T, amplitude_rep::T, inv_screening_length_rep::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    sk_function = S_SALR_RPA(ϕ, amplitude_attr, inv_screening_length_attr, amplitude_rep, inv_screening_length_rep)
    return sk_function.(k_values)
end

