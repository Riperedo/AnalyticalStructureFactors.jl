# examples/plot_monodisperse_HS.jl
#
# Este script calcula el factor de estructura estático S(k) para
# esferas duras usando la aproximación de Percus-Yevick (S_HS_PY)
# para diferentes fracciones de empaquetamiento y un rango de k.
# Los resultados se guardan en un archivo .dat.
#
# Para ejecutar desde la raíz del proyecto AnalyticalStructureFactors.jl:
# julia --project examples/plot_monodisperse_HS.jl

# Importar el paquete principal.
# Si el paquete está en desarrollo y no en el entorno global,
# asegúrate de que el entorno del proyecto esté activo.
# (Esto se hace con `julia --project ...` desde la raíz del paquete)
using AnalyticalStructureFactors
# Gnuplot.jl no es una dependencia de AnalyticalStructureFactors,
# pero si quisieras graficar directamente desde este script, lo añadirías aquí
# y al Project.toml de la carpeta examples si tuviera uno, o al global.
# using Gnuplot 

println("--- Iniciando script de ejemplo: Monodisperse Hard Spheres (PY) ---")

# --- Definición de Parámetros ---

# Rango de k (vector de onda adimensional, k = qσ)
k_min = 0.01
k_max = 15.0
k_step = 0.05
k_values = collect(k_min:k_step:k_max)
println("Rango de k: de $k_min a $k_max con paso $k_step. Total puntos: $(length(k_values))")

# Fracciones de empaquetamiento (ϕ) a explorar
phi_values = [0.1, 0.2, 0.3, 0.4]
println("Fracciones de empaquetamiento a explorar (ϕ): $phi_values")

# --- Cálculo de S(k) ---

# Almacenaremos los resultados en una matriz.
# La primera columna será k, las siguientes S(k) para cada phi.
# El número de columnas de datos de S(k) es length(phi_values).
# El número de filas es length(k_values).
s_k_results = zeros(Float64, length(k_values), length(phi_values))

println("Calculando S(k) para cada ϕ...")
for (i, phi) in enumerate(phi_values)
    println("  Calculando para ϕ = $phi...")
    # Usamos broadcasting para aplicar S_HS_PY a cada elemento de k_values
    # S_HS_PY(phi, k_val) o S_HS_PY.(Ref(phi), k_values)
    # Asumimos que S_HS_PY está exportada y toma (phi, k)
    s_k_for_phi = S_HS_PY.(Ref(phi), k_values)
    s_k_results[:, i] = s_k_for_phi
end
println("Cálculos de S(k) completados.")

# --- Preparar Datos para Guardar ---
# Concatenar el vector k_values como la primera columna
data_to_save = hcat(k_values, s_k_results)

# --- Definir Cabecera para el Archivo .dat ---
header_info = ["Resultados del Factor de Estructura para Esferas Duras (Percus-Yevick)"]
push!(header_info, "Generado por: examples/plot_monodisperse_HS.jl")
push!(header_info, "Column 1: k (vector de onda adimensional qσ)")
for (i, phi) in enumerate(phi_values)
    push!(header_info, "Column $(i+1): S(k) para phi = $phi")
end

# --- Guardar Datos ---
# Asegúrate de que el directorio 'examples_output' exista o ajusta la ruta.
output_directory = "examples_output"
if !isdir(output_directory)
    mkdir(output_directory)
    println("Directorio creado: $output_directory")
end
output_filepath = joinpath(output_directory, "monodisperse_HS_PY_Sk.dat")

# Usamos la función save_data de AnalyticalStructureFactors.utils
# (asumiendo que está exportada o la llamas con AnalyticalStructureFactors.save_data)
try
    save_data(output_filepath, data_to_save, header_lines=header_info, delimiter='\t')
catch e
    println(stderr, "Error al intentar guardar datos. Asegúrate de que la función 'save_data' esté disponible y exportada por el módulo AnalyticalStructureFactors.")
    showerror(stderr, e)
    println(stderr)
end

println("--- Script de ejemplo finalizado ---")

# --- Opcional: Graficar con Gnuplot.jl (si está instalado y configurado) ---
# Descomenta las siguientes líneas si quieres una gráfica rápida.
# Necesitarás tener Gnuplot instalado en tu sistema y Gnuplot.jl en tu entorno.
#
# println("Intentando graficar S(k) vs k...")
# try
#     using Gnuplot # Asegúrate de que Gnuplot.jl esté en el Project.toml de este script/entorno
#     gp_plot = Gnuplot.GPPlot() # Crea un objeto de gráfica
#     for (i, phi) in enumerate(phi_values)
#         # @gp gp_plot k_values s_k_results[:,i] "with lines title 'S(k) ϕ=$phi'" :-
#         # La macro @gp puede ser un poco quisquillosa con el scope en scripts.
#         # Una forma más programática:
#         if i == 1
#             Gnuplot.plot(gp_plot, k_values, s_k_results[:,i], legend="S(k) ϕ=$phi", 韩式="l") # 'l' for lines
#         else
#             Gnuplot.plot!(gp_plot, k_values, s_k_results[:,i], legend="S(k) ϕ=$phi", 韩式="l")
#         end
#     end
#     Gnuplot.xlabel(gp_plot, "k (qσ)")
#     Gnuplot.ylabel(gp_plot, "S(k)")
#     Gnuplot.title(gp_plot, "Factor de Estructura de Esferas Duras (PY)")
#     Gnuplot.grid(gp_plot, "on")
#     Gnuplot.figure(gp_plot) # Muestra la gráfica
#     # save(term="pngcairo", output=joinpath(output_directory, "monodisperse_HS_PY_Sk.png"))
#     println("Gráfica generada (si Gnuplot está configurado).")
# catch e
#     println(stderr, "No se pudo graficar con Gnuplot.jl. Asegúrate de que esté instalado y configurado.")
#     println(stderr, "Error: $e")
# end

