# src/AnalyticalStructureFactors.jl

module AnalyticalStructureFactors

# --- Dependencias del Paquete (si las hay) ---
# using SpecialFunctions
# using QuadGK # Si decides usarlo para la función blip, por ejemplo

# --- Inclusión de Archivos de Código Fuente ---

# Utilidades generales y correcciones
include("utils.jl") # Contiene blip, phi_VW, k_VW

# Transformadas de Fourier de potenciales de par
include("potentials.jl") # Contiene betaU_SW, betaU_Yukawa

# Modelos específicos de factores de estructura
# Nota: Usamos joinpath para construir rutas de archivo de forma robusta
include(joinpath("models", "hardsphere.jl"))    # Lógica para HS-PY, HS-VW, WCA-blip
include(joinpath("models", "squarewell.jl"))  # Lógica para S(k) con potencial de pozo cuadrado
include(joinpath("models", "yukawa.jl"))      # (Futuro) Para modelos basados en Yukawa
include(joinpath("models", "SALR.jl"))
# include("mixtures.jl")                       # (Futuro) Para sistemas de mezclas

# --- Exportaciones de la API Pública ---

# Desde utils.jl
export blip, phi_VW, k_VW

# Desde potentials.jl (exportar si los usuarios las necesitan directamente)
export betaU_SW, betaU_Yukawa

# Desde models/hardsphere.jl
export S_HS_PY, C_HS_PY, IS_HS_PY  # Funciones básicas de Percus-Yevick
export S_HS_VW, C_HS_VW, IS_HS_VW  # Funciones con corrección Verlet-Weiss
export S_WCA_blip, C_WCA_blip, IS_WCA_blip # Factor de estructura WCA usando blip

# Desde models/
export S_SW_RPA # Ejemplo de nombre para S(k) de Pozo Cuadrado con RPA
export S_Yukawa_RPA # Esto exportará todos los métodos de S_Yukawa_RPA
export S_SALR_RPA


# Aquí podrías añadir más exportaciones a medida que desarrollas el paquete.
# Considera también las funciones "fábrica" que devuelven f(k) si mantienes ese patrón.

end # module
