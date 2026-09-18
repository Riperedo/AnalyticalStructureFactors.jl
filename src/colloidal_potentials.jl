# src/colloidal_potentials.jl

# ==============================================================================
# COLLOIDAL ELECTROSTATICS AND DLVO / GOCM POTENTIAL PARAMETERS
# ==============================================================================
# This module provides utilities to convert physical laboratory quantities
# (electrolyte concentration, salt valence, surface charge, macroion size, temperature,
# dielectric constant) into the dimensionless parameters (K, z) used by analytical
# Yukawa and DLVO static structure factors.
#
# References:
# - Derjaguin, B., & Landau, L. (1941). Acta Physicochim. URSS, 14, 633.
# - Verwey, E. J. W., & Overbeek, J. T. G. (1948). Theory of the Stability
#   of Lyophobic Colloids. Elsevier.
# - Belloni, L. (1998). Colloidal interactions. J. Phys.: Condens. Matter, 10(48), R549.
# - Hayter, J. B., & Penfold, J. (1981). Mol. Phys., 42(1), 109-118.
# ==============================================================================

# Physical constants in SI units
const KB_SI = 1.380649e-23        # Boltzmann constant (J/K)
const EPSILON_0_SI = 8.8541878128e-12 # Vacuum permittivity (F/m)
const ELEMENTARY_CHARGE_SI = 1.602176634e-19 # Elementary charge (C)
const AVOGADRO_SI = 6.02214076e23  # Avogadro constant (mol^-1)

@doc """
    ionic_strength(concentrations::AbstractVector{T}, valences::AbstractVector{T}) where {T<:AbstractFloat}

Calculates the molar ionic strength ``I`` of an electrolyte solution (in ``\\text{mol}/\\text{L}`` or ``\\text{M}``):
``I = \\frac{1}{2} \\sum_i c_i z_i^2``

# Arguments
- `concentrations::AbstractVector{T}`: Molar concentrations of each ionic species (in M = mol/L).
- `valences::AbstractVector{T}`: Valence (charge number) of each ionic species.

# Returns
- `T`: Molar ionic strength ``I`` (mol/L).
"""
function ionic_strength(concentrations::AbstractVector{T}, valences::AbstractVector{T}) where {T<:AbstractFloat}
    if length(concentrations) != length(valences)
        throw(DimensionMismatch("concentrations and valences must have the same length."))
    end
    return T(0.5) * sum(concentrations .* (valences .^ 2))
end

@doc """
    bjerrum_length(temperature_K::T, epsilon_r::T) where {T<:AbstractFloat}

Calculates the Bjerrum length ``\\lambda_B`` (in meters), which is the separation
at which the electrostatic interaction between two elementary charges equals the thermal energy ``k_B T``:
``\\lambda_B = \\frac{e^2}{4\\pi \\varepsilon_0 \\varepsilon_r k_B T}``

# Arguments
- `temperature_K::T`: Absolute temperature in Kelvin (``T > 0``).
- `epsilon_r::T`: Relative permittivity (dielectric constant) of the solvent (e.g., ``\\approx 78.4`` for water at 298 K).

# Returns
- `T`: Bjerrum length ``\\lambda_B`` in meters.
"""
function bjerrum_length(temperature_K::T, epsilon_r::T) where {T<:AbstractFloat}
    if temperature_K <= zero(T)
        throw(DomainError(temperature_K, "Temperature must be strictly positive."))
    end
    if epsilon_r <= zero(T)
        throw(DomainError(epsilon_r, "Relative permittivity must be strictly positive."))
    end
    num = (T(ELEMENTARY_CHARGE_SI)^2)
    den = T(4.0) * T(π) * T(EPSILON_0_SI) * epsilon_r * T(KB_SI) * temperature_K
    return num / den
end

@doc """
    debye_screening_length(temperature_K::T, epsilon_r::T, I_molar::T) where {T<:AbstractFloat}

Calculates the Debye screening length ``\\kappa^{-1}`` (in meters):
``\\kappa^{-1} = \\sqrt{\\frac{\\varepsilon_0 \\varepsilon_r k_B T}{2000 N_A e^2 I}} = \\frac{1}{\\sqrt{8\\pi \\lambda_B N_A \\cdot 1000 I}}``

# Arguments
- `temperature_K::T`: Temperature in Kelvin.
- `epsilon_r::T`: Relative permittivity of the solvent.
- `I_molar::T`: Molar ionic strength (in mol/L).

# Returns
- `T`: Debye length ``\\kappa^{-1}`` in meters.
"""
function debye_screening_length(temperature_K::T, epsilon_r::T, I_molar::T) where {T<:AbstractFloat}
    if I_molar <= zero(T)
        return T(Inf)
    end
    l_B = bjerrum_length(temperature_K, epsilon_r)
    # Number density in particles/m^3: n_i = 1000 * N_A * I
    kappa_sq = T(8.0) * T(π) * l_B * (T(AVOGADRO_SI) * T(1000.0) * I_molar)
    return one(T) / sqrt(kappa_sq)
end

@doc """
    dlvo_yukawa_parameters(Z_effective::T, diameter_m::T, temperature_K::T, epsilon_r::T, I_molar::T) where {T<:AbstractFloat}

Calculates the dimensionless DLVO Yukawa parameters ``(K, z)`` used in ``S_Yukawa_MSA`` and ``S_Yukawa_RPA``:
- Dimensionless screening parameter: ``z = \\kappa \\sigma``
- Dimensionless contact energy amplitude: ``K = \\frac{Z_{eff}^2 \\lambda_B}{\\sigma (1 + z/2)^2}``

# Arguments
- `Z_effective::T`: Effective surface charge of the macroion (in elementary charge units ``e``).
- `diameter_m::T`: Hard-core particle diameter ``\\sigma`` in meters.
- `temperature_K::T`: Temperature in Kelvin.
- `epsilon_r::T`: Relative permittivity of the medium.
- `I_molar::T`: Molar ionic strength in mol/L.

# Returns
- `NamedTuple{(:K, :z), Tuple{T, T}}`: Dimensionless Yukawa contact amplitude `K` and inverse screening parameter `z`.

# References
- Verwey, E. J. W., & Overbeek, J. T. G. (1948). Theory of the Stability of Lyophobic Colloids.
- Hayter, J. B., & Penfold, J. (1981). Mol. Phys., 42, 109.
"""
function dlvo_yukawa_parameters(Z_effective::T, diameter_m::T, temperature_K::T, epsilon_r::T, I_molar::T) where {T<:AbstractFloat}
    l_B = bjerrum_length(temperature_K, epsilon_r)
    inv_kappa = debye_screening_length(temperature_K, epsilon_r, I_molar)
    kappa = one(T) / inv_kappa
    z = kappa * diameter_m

    # DLVO Yukawa contact amplitude (Derjaguin-Landau-Verwey-Overbeek)
    K = (Z_effective^2 * l_B) / (diameter_m * ((one(T) + z / T(2.0))^2))

    return (K=K, z=z)
end
