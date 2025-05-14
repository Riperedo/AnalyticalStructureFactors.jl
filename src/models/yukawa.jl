# src/models/yukawa.jl

# Este archivo contendrá la lógica para factores de estructura
# basados en el potencial de Yukawa.

# --- Random Phase Approximation (RPA) para Potencial de Yukawa ---

@doc """
    S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T, k::T) where {T<:AbstractFloat}

Calcula el factor de estructura estático para un fluido con potencial de Yukawa
utilizando la Random Phase Approximation (RPA).
El sistema de referencia es de esferas duras con corrección de Verlet-Weiss.

La fórmula general de RPA es: S(k) = 1 / ( 1/S₀(k) - ρ_factor * β * u_pert_tilde(k) )
donde S₀(k) es el factor de estructura del sistema de referencia,
ρ_factor es un factor relacionado con la densidad (usualmente 6ϕ/π para σ=1),
y u_pert_tilde(k) es la transformada de Fourier de la parte perturbativa del potencial.

# Arguments
- `ϕ::T`: Fracción de volumen de las esferas.
- `amplitude::T`: Amplitud efectiva del potencial de Yukawa (a menudo incluye β, ej., βε).
- `inv_screening_length::T`: Longitud de apantallamiento inversa adimensional (z).
- `k::T`: Vector de onda adimensional (k = qσ).

# Returns
- `T`: Valor del factor de estructura S(k).
"""
function S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T, k::T) where {T<:AbstractFloat}
    inv_s0_k = IS_HS_VW(ϕ, k) 
    beta_u_pert_k = betaU_Yukawa(amplitude, inv_screening_length, k)
    density_factor = (T(6.0) * ϕ) / T(π)
    
    inv_s_k = inv_s0_k - density_factor * beta_u_pert_k
    
    # La siguiente advertencia ha sido eliminada según la solicitud del usuario.
    # if inv_s_k <= zero(T)
    #     @warn "Inverse structure factor 1/S(k) = $inv_s_k for Yukawa-RPA is non-positive. Resulting S(k) will be unphysical."
    # end

    return one(T) / inv_s_k
end

# --- Versión "Fábrica" (Opcional, pero sigue el patrón de tu proyecto anterior) ---
@doc """
    S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T) where {T<:AbstractFloat}

Retorna una función `f(k)` que calcula el factor de estructura estático para un fluido
con potencial de Yukawa utilizando la RPA. Los parámetros del sistema y del potencial
(ϕ, amplitud, longitud de apantallamiento inversa) se fijan.

# Arguments
- `ϕ::T`: Fracción de volumen de las esferas.
- `amplitude::T`: Amplitud efectiva del potencial de Yukawa.
- `inv_screening_length::T`: Longitud de apantallamiento inversa adimensional (z).

# Returns
- `Function`: Una función que toma un vector de onda `k` y retorna `S(k)`.
"""
function S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T) where {T<:AbstractFloat}
    function structure_factor_k(k_val::T)
        return S_Yukawa_RPA(ϕ, amplitude, inv_screening_length, k_val)
    end
    return structure_factor_k
end

# --- Versión Vectorizada (Opcional, para conveniencia) ---
@doc """
    S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}

Calcula el factor de estructura S(k) para un vector de valores de k.

# Arguments
- `ϕ::T`: Fracción de volumen.
- `amplitude::T`: Amplitud del potencial de Yukawa.
- `inv_screening_length::T`: Longitud de apantallamiento inversa.
- `k_values::AbstractVector{T}`: Un vector de valores de k adimensionales.

# Returns
- `Vector{T}`: Un vector de valores de S(k) correspondientes a cada k en `k_values`.
"""
function S_Yukawa_RPA(ϕ::T, amplitude::T, inv_screening_length::T, k_values::AbstractVector{T}) where {T<:AbstractFloat}
    sk_function = S_Yukawa_RPA(ϕ, amplitude, inv_screening_length)
    return sk_function.(k_values)
end

