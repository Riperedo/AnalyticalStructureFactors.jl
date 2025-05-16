# examples/plot_monodisperse_SW_RPA.jl
#
# Este script calcula el factor de estructura estático S(k) para
# un sistema monocomponente con potencial de Pozo Cuadrado (Square-Well)
# usando la aproximación de Fase Aleatoria (RPA) sobre un sistema de
# referencia de esferas duras (S_SW_RPA).
# Los resultados se guardan en un archivo .dat.
#
# Para ejecutar desde la raíz del proyecto AnalyticalStructureFactors.jl:
# julia --project examples/plot_monodisperse_SW_RPA.jl

using AnalyticalStructureFactors
# using Gnuplot # Descomenta si quieres usar la sección de graficación

println("--- Iniciando script de ejemplo: Monodisperse Square-Well (RPA) ---")

# --- Definición de Parámetros ---

# Rango de k (vector de onda adimensional, k = qσ)
k_min = 0.01
k_max = 20.0 # Ampliamos un poco para ver más detalles del pozo cuadrado
k_step = 0.05
k_values = collect(k_min:k_step:k_max)
println("Rango de k: de $k_min a $k_max con paso $k_step. Total puntos: $(length(k_values))")

# Fracción de empaquetamiento (ϕ) fija para este ejemplo
phi_fixed = 0.2
println("Fracción de empaquetamiento fija (ϕ): $phi_fixed")

# Parámetros del potencial de Pozo Cuadrado a explorar
# Cada elemento de `sw_params` será una tupla: (temperatura_efectiva, rango_lambda)
# temperatura_efectiva: T* (relacionada con la profundidad del pozo, ej. kBT/ε)
# rango_lambda: λ (ancho del pozo en unidades de σ)
sw_parameter_sets = [
    (temperature=1.0, lambda_range=1.5, label="T=1.0_L=1.5"), # Caso 1
    (temperature=1.5, lambda_range=1.5, label="T=1.5_L=1.5"), # Caso 2
    (temperature=1.0, lambda_range=1.8, label="T=1.0_L=1.8"), # Caso 3
    (temperature=0.8, lambda_range=1.2, label="T=0.8_L=1.2")  # Caso 4 (pozo más profundo y corto)
]
println("Conjuntos de parámetros para Pozo Cuadrado (T*, λ) a explorar:")
for params in sw_parameter_sets
    println("  Temperatura efectiva: $(params.temperature), Rango λ: $(params.lambda_range)")
end

# --- Cálculo de S(k) ---

# Almacenaremos los resultados en una matriz.
# La primera columna será k, las siguientes S(k) para cada conjunto de parámetros del pozo.
s_k_results_sw = zeros(Float64, length(k_values), length(sw_parameter_sets))

println("Calculando S(k) para cada conjunto de parámetros del Pozo Cuadrado...")
for (i, params) in enumerate(sw_parameter_sets)
    temp_eff = params.temperature
    lambda_r = params.lambda_range
    println("  Calculando para T* = $temp_eff, λ = $lambda_r...")
    
    # Asumimos que S_SW_RPA está exportada y toma (phi, temp_eff, lambda_range, k)
    s_k_for_params = S_SW_RPA.(Ref(phi_fixed), Ref(temp_eff), Ref(lambda_r), k_values)
    s_k_results_sw[:, i] = s_k_for_params
end
println("Cálculos de S(k) para Pozo Cuadrado completados.")

# --- Preparar Datos para Guardar ---
data_to_save_sw = hcat(k_values, s_k_results_sw)

# --- Definir Cabecera para el Archivo .dat ---
header_info_sw = ["Resultados del Factor de Estructura para Pozo Cuadrado Monocomponente (RPA)"]
push!(header_info_sw, "Generado por: examples/plot_monodisperse_SW_RPA.jl")
push!(header_info_sw, "Fracción de empaquetamiento (phi) = $phi_fixed")
push!(header_info_sw, "Column 1: k (vector de onda adimensional qσ)")
for (i, params) in enumerate(sw_parameter_sets)
    push!(header_info_sw, "Column $(i+1): S(k) para T*=$(params.temperature), lambda=$(params.lambda_range) (Label: $(params.label))")
end

# --- Guardar Datos ---
output_directory = "examples_output" # Mismo directorio que el ejemplo anterior
if !isdir(output_directory)
    mkdir(output_directory)
    println("Directorio creado: $output_directory")
end
output_filepath_sw = joinpath(output_directory, "monodisperse_SW_RPA_Sk.dat")

try
    save_data(output_filepath_sw, data_to_save_sw, header_lines=header_info_sw, delimiter='\t')
catch e
    println(stderr, "Error al intentar guardar datos. Asegúrate de que la función 'save_data' esté disponible y exportada por el módulo AnalyticalStructureFactors.")
    showerror(stderr, e)
    println(stderr)
end

println("--- Script de ejemplo para Pozo Cuadrado finalizado ---")

# --- Opcional: Graficar con Gnuplot.jl ---
# println("Intentando graficar S(k) vs k para Pozo Cuadrado...")
# try
#     using Gnuplot
#     gp_plot_sw = Gnuplot.GPPlot()
#     for (i, params) in enumerate(sw_parameter_sets)
#         plot_legend = "S(k) T*=$(params.temperature), λ=$(params.lambda_range)"
#         if i == 1
#             Gnuplot.plot(gp_plot_sw, k_values, s_k_results_sw[:,i], legend=plot_legend, 韩式="l")
#         else
#             Gnuplot.plot!(gp_plot_sw, k_values, s_k_results_sw[:,i], legend=plot_legend, 韩式="l")
#         end
#     end
#     Gnuplot.xlabel(gp_plot_sw, "k (qσ)")
#     Gnuplot.ylabel(gp_plot_sw, "S(k)")
#     Gnuplot.title(gp_plot_sw, "Factor de Estructura Pozo Cuadrado (RPA), ϕ=$phi_fixed")
#     Gnuplot.grid(gp_plot_sw, "on")
#     Gnuplot.figure(gp_plot_sw) # Muestra la gráfica
#     # save(term="pngcairo", output=joinpath(output_directory, "monodisperse_SW_RPA_Sk.png"))
#     println("Gráfica generada (si Gnuplot está configurado).")
# catch e
#     println(stderr, "No se pudo graficar con Gnuplot.jl. Asegúrate de que esté instalado y configurado.")
#     println(stderr, "Error: $e")
# end

