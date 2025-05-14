# src/models/hardsphere.jl

@doc """
    C_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}

Direct correlation function for Hard Spheres using Percus-Yevick
with Verlet-Weiss corrections.
Utiliza `phi_VW` y `k_VW` de `utils.jl`.
"""
function C_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}
    # phi_VW y k_VW deberían estar disponibles aquí porque utils.jl se incluye
    # en el módulo principal AnalyticalStructureFactors.
    # C_HS_PY también debería estar disponible.
    return C_HS_PY(phi_VW(ϕ_HS), k_VW(ϕ_HS, k_HS))
end

@doc """
    IS_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}

Inverse static structure factor for Hard Spheres using Percus-Yevick
with Verlet-Weiss corrections.
"""
function IS_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}
    return IS_HS_PY(phi_VW(ϕ_HS), k_VW(ϕ_HS, k_HS))
end

@doc """
    S_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}

Static structure factor for Hard Spheres using Percus-Yevick
with Verlet-Weiss corrections.
"""
function S_HS_VW(ϕ_HS::T, k_HS::T) where {T<:AbstractFloat}
    return S_HS_PY(phi_VW(ϕ_HS), k_VW(ϕ_HS, k_HS)) # O 1.0 / IS_HS_VW(ϕ_HS, k_HS)
end


# --- Aproximación WCA (Weeks-Chandler-Andersen) vía Blip Function ---
# Estas funciones usan las funciones HS-VW y la función blip de utils.jl

@doc """
    C_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}

Direct correlation function using the WCA approximation with a blip function
to map to an effective hard-sphere system with Verlet-Weiss corrections.
"""
function C_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}
    # blip debería estar disponible aquí desde utils.jl
    eff_lambda, eff_lambda_cubed = blip(temperature, ν=ν)
    # C_HS_VW debería estar disponible
    return C_HS_VW(eff_lambda_cubed * ϕ, eff_lambda * k)
end

@doc """
    IS_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}

Inverse static structure factor using the WCA approximation with a blip function.
"""
function IS_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}
    eff_lambda, eff_lambda_cubed = blip(temperature, ν=ν)
    return IS_HS_VW(eff_lambda_cubed * ϕ, eff_lambda * k)
end

@doc """
    S_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}

Static structure factor using the WCA approximation with a blip function.
"""
function S_WCA_blip(ϕ::T, temperature::T, k::T; ν::Integer = 6) where {T<:AbstractFloat}
    eff_lambda, eff_lambda_cubed = blip(temperature, ν=ν)
    return S_HS_VW(eff_lambda_cubed * ϕ, eff_lambda * k) # O 1.0 / IS_WCA_blip(...)
end

