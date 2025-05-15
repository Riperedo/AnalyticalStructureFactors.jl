# src/AnalyticalStructureFactors.jl

module AnalyticalStructureFactors

# --- Dependencias del Paquete (si las hay) ---
# using SpecialFunctions
# using QuadGK # Si decides usarlo para la función blip, por ejemplo
using LinearAlgebra # Necesario para operaciones matriciales en mezclas (inv, ')

# --- Inclusión de Archivos de Código Fuente ---

# Utilidades generales y correcciones
include("utils.jl") # Contiene blip, phi_VW (monodisp.), k_VW (monodisp.), phi_to_rho_mixture, rho_to_phi_mixture

# Transformadas de Fourier de potenciales de par
include("potentials.jl") # Contiene betaU_SW, betaU_Yukawa (para monodispersos o pares)

# --- Modelos Monodispersos ---
include(joinpath("models", "hardsphere.jl"))    # Lógica para HS-PY, HS-VW, WCA-blip
include(joinpath("models", "squarewell.jl"))  # Lógica para S(k) con potencial de pozo cuadrado (monodisperso)
include(joinpath("models", "yukawa.jl"))      # Lógica para S(k) con potencial Yukawa (monodisperso)
include(joinpath("models", "SALR.jl"))        # Lógica para S(k) con potencial SALR (monodisperso)

# --- Modelos de Mezclas ---
# Sistema de referencia para mezclas
include(joinpath("models", "mixtures", "hardsphere_mixture.jl")) # S_HS_Baxter_mixture, IS_HS_Baxter_mixture
include(joinpath("models", "mixtures", "vw_mixture.jl"))         # VW_correction_mixture, S_HS_VW_mixture, IS_HS_VW_mixture

# RPA para mezclas
include(joinpath("models", "mixtures", "yukawa_rpa_mixture.jl"))       # S_RPA_mixture_Yukawa
include(joinpath("models", "mixtures", "squarewell_rpa_mixture.jl")) # S_RPA_mixture_SquareWell
include(joinpath("models", "mixtures", "SALR_rpa_mixture.jl"))       # S_RPA_mixture_SALR


# --- Exportaciones de la API Pública ---

# Desde utils.jl
export blip
export phi_VW, k_VW # Para sistemas monodispersos
export phi_to_rho_mixture, rho_to_phi_mixture # Para mezclas

# Desde potentials.jl (exportar si los usuarios las necesitan directamente)
export betaU_SW, betaU_Yukawa

# Desde models/hardsphere.jl (Monodisperso)
export S_HS_PY, C_HS_PY, IS_HS_PY
export S_HS_VW, C_HS_VW, IS_HS_VW
export S_WCA_blip, C_WCA_blip, IS_WCA_blip

# Desde models/squarewell.jl (Monodisperso)
export S_SW_RPA # Asumiendo que este es el nombre de tu función principal aquí

# Desde models/yukawa.jl (Monodisperso)
export S_Yukawa_RPA # Asumiendo que este es el nombre de tu función principal aquí

# Desde models/SALR.jl (Monodisperso)
export S_SALR_RPA # Asumiendo que este es el nombre de tu función principal aquí

# --- Exportaciones para Mezclas ---

# Desde models/mixtures/hardsphere_mixture.jl
export Qk_mixture, S_HS_Baxter_mixture, IS_HS_Baxter_mixture

# Desde models/mixtures/vw_mixture.jl
export VW_correction_mixture, S_HS_VW_mixture, IS_HS_VW_mixture

# Desde models/mixtures/yukawa_rpa_mixture.jl
export S_RPA_mixture_Yukawa

# Desde models/mixtures/squarewell_rpa_mixture.jl
export S_RPA_mixture_SquareWell

# Desde models/mixtures/SALR_rpa_mixture.jl
export S_RPA_mixture_SALR


# Considera también las funciones "fábrica" que devuelven f(k) si mantienes ese patrón
# y quieres exportarlas explícitamente. Si una función tiene múltiples métodos (como las
# versiones fábrica y las que toman k), exportar el nombre de la función exporta todos sus métodos.

end # module
