# AnalyticalStructureFactors.jl

```@meta
CurrentModule = AnalyticalStructureFactors
Bienvenido a la documentación de AnalyticalStructureFactors.jl.Este paquete proporciona herramientas para calcular factores de estructura estáticos S(q) utilizando soluciones analíticas para diversos modelos físicos, tanto para sistemas monodispersos como para mezclas. Está diseñado para ser una herramienta eficiente y fácil de usar para investigadores en física de la materia condensada, física química y ciencia de materiales.EmpezandoPara comenzar a usar el paquete, primero asegúrate de tenerlo instalado. Si aún no lo has hecho, puedes instalarlo desde el REPL de Julia:import Pkg
# Pkg.add("AnalyticalStructureFactors") # Descomenta cuando esté registrado
Pkg.add(url="[https://github.com/Riperedo/AnalyticalStructureFactors.jl](https://github.com/Riperedo/AnalyticalStructureFactors.jl)") # Para instalar desde GitHub
Luego, puedes empezar a usar las funciones importando el paquete:using AnalyticalStructureFactors

# Ejemplo (reemplaza con un ejemplo real y simple de tu paquete):
# ϕ_example = 0.3 # Fracción de volumen
# k_example = 1.0 # Vector de onda adimensional q*sigma
# sq_val_example = S_HS_PY(ϕ_example, k_example)

# Usando interpolación de cadenas de forma segura:
# println("S(k) para esferas duras (PY) con ϕ=$(ϕ_example), k=$(k_example) es: $(sq_val_example)")
Contenido de la DocumentaciónPages = [
    "index.md",
    # "man/monodisperse_theory.md", # Descomenta y ajusta si tienes estas páginas
    # "man/mixtures_theory.md",
    "lib/public_api.md"
]
Referencia Rápida de la APIUna lista completa de las funciones y tipos exportados se puede encontrar en la sección de Referencia de la API.A continuación, se muestra una selección de las funcionalidades principales (esto se generará automáticamente si configuras @autodocs correctamente aquí o en la página de la API):Modules = [AnalyticalStructureFactors]
Order   = [:module, :constant, :type, :macro, :function]
Private = false
Contribuciones y Reporte de ErroresLas contribuciones al paquete son bienvenidas. Si deseas contribuir o reportar un error, por favor visita el repositorio en GitHub y abre un "issue".