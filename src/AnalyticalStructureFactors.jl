module AnalyticalStructureFactors

# Importar dependencias si las tienes (ej. SpecialFunctions, Parameters)
# using SpecialFunctions
# using Parameters

# Incluir los otros archivos de código fuente
include("utils.jl")
include("monodisperse.jl")
include("mixtures.jl")

# Exportar las funciones que quieres que sean públicas y accesibles para los usuarios
# Ejemplo de exportaciones para funciones de monodisperse.jl:
export S_HS_PY, C_HS_PY, IS_HS_PY

# Ejemplo de exportaciones para funciones de mixtures.jl:
#export S_ij_binary_hardsphere_py, S_total_binary_mixture # Nombres de ejemplo
# Exportar tipos de datos personalizados si los tienes (ver utils.jl)
# export HardSphereParams # Ejemplo

end # module