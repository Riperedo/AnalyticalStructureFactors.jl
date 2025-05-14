# docs/make.jl
using Documenter
using AnalyticalStructureFactors # Importa tu paquete

# Asegúrate de que el directorio actual sea docs/
cd(@__DIR__)

makedocs(
    sitename = "AnalyticalStructureFactors.jl",
    authors = "Ricardo Peredo Ortiz",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = String[], # Puedes añadir CSS personalizado aquí
    ),
    modules = [AnalyticalStructureFactors],
    pages = [
        "Home" => "index.md",
        # "Manual" => [
        #     "Guide Page 1" => "man/guide1.md",
        #     "Guide Page 2" => "man/guide2.md",
        # ],
        "API Reference" => "lib/public_api.md", # Asegúrate que este archivo exista en docs/src/lib/
    ],
    # --- AJUSTE IMPORTANTE AQUÍ ---
    # `strict` por defecto es true, lo que significa que errores como :missing_docs y :doctest detendrán la compilación.
    # Usamos `warnonly` para especificar qué errores deben ser solo advertencias y no detener la compilación.
    # strict = true, # Puedes omitir esta línea ya que true es el valor por defecto.
    warnonly = [:missing_docs] # Esto hará que los "docstrings no incluidos" solo muestren una advertencia.
                               # Otros errores, como los de doctest, seguirán deteniendo la compilación.
)

deploydocs(
    repo = "github.com/Riperedo/AnalyticalStructureFactors.jl.git", # URL corregida (sin formato Markdown)
    devbranch = "main",
    # push_preview = true,
)
