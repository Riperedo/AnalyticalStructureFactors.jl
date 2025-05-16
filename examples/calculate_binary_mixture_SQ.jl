# examples/calculate_binary_mixture_SQ.jl
#
# Este script calcula los factores de estructura parciales S_ij(k) y el
# factor de estructura total S_T(k) para una mezcla binaria utilizando la
# Random Phase Approximation (RPA) con un potencial de perturbación de Yukawa.
# El sistema de referencia es una mezcla de esferas duras con corrección de Verlet-Weiss.
# Los resultados se guardan en un archivo .dat.
#
# Para ejecutar desde la raíz del proyecto AnalyticalStructureFactors.jl:
# julia --project examples/calculate_binary_mixture_SQ.jl

using AnalyticalStructureFactors
# using Gnuplot # Descomenta si quieres usar la sección de graficación

println("--- Iniciando script de ejemplo: Mezcla Binaria con Potencial Yukawa (RPA) ---")

# --- Definición de Parámetros ---

# Rango de k (vector de onda adimensional, k = qσ_ref)
k_min = 0.01
k_max = 15.0
k_step = 0.05
k_values = collect(k_min:k_step:k_max)
println("Rango de k: de $k_min a $k_max con paso $k_step. Total puntos: $(length(k_values))")

# Parámetros de la Mezcla Binaria (Sistema de Referencia de Esferas Duras)
# Componente 1 y Componente 2
sigma_vector = [1.0, 0.8]  # Diámetros σ₁ , σ₂
phi_vector   = [0.15, 0.1] # Fracciones de empaquetamiento ϕ₁, ϕ₂
println("Parámetros de la mezcla:")
println("  σ = $sigma_vector")
println("  ϕ = $phi_vector (fracciones de empaquetamiento de componentes)")

# Parámetros del Potencial de Perturbación de Yukawa entre pares (i,j)
# Matriz de amplitudes A_ij (pueden ser positivas para repulsión, negativas para atracción)
# A_ij = A_ji
A_matrix = [-1.0  -0.5;  # A₁₁, A₁₂
            -0.5  -0.8]  # A₂₁, A₂₂
# Matriz de inversas de longitudes de apantallamiento Z_ij
# Z_ij = Z_ji
Z_matrix = [1.8  2.0;   # Z₁₁, Z₁₂
            2.0  2.2]   # Z₂₁, Z₂₂
println("Parámetros del potencial de Yukawa (perturbación):")
println("  Matriz de Amplitudes A_ij = \n$A_matrix")
println("  Matriz de Inversas de Long. de Apantallamiento Z_ij = \n$Z_matrix")

# --- Cálculo de Factores de Estructura Parciales S_ij(k) ---

println("Calculando S_ij(k) para la mezcla con RPA-Yukawa...")
# S_RPA_mixture_Yukawa devuelve un vector de matrices S_ij(k)
# (una matriz para cada valor de k en k_values)
s_ij_matrices = S_RPA_mixture_Yukawa(sigma_vector, phi_vector, 
                                     A_matrix, Z_matrix, 
                                     k_values)
println("Cálculos de S_ij(k) completados.")

# Extraer los componentes S11, S12 (S21), S22 para cada k
s11_values = [s_matrix[1,1] for s_matrix in s_ij_matrices]
s12_values = [s_matrix[1,2] for s_matrix in s_ij_matrices] # S12 = S21
s22_values = [s_matrix[2,2] for s_matrix in s_ij_matrices]

# --- Cálculo del Factor de Estructura Total S_T(k) ---
# S_T(k) = x₁S₁₁(k) + x₂S₂₂(k) + 2√(x₁x₂)S₁₂(k)
# donde x_i son las fracciones molares.

# Calcular densidades numéricas y fracciones molares
rho_vector = phi_to_rho_mixture(phi_vector, sigma_vector)
rho_total = sum(rho_vector)
x_vector = rho_vector ./ rho_total
x1 = x_vector[1]
x2 = x_vector[2]
println("Fracciones molares: x₁ = $x1, x₂ = $x2")

s_total_values = (x1 .* s11_values) .+ 
                 (x2 .* s22_values) .+ 
                 (2 * sqrt(x1 * x2) .* s12_values)
println("Cálculo de S_Total(k) completado.")

# --- Preparar Datos para Guardar ---
# Columnas: k, S11(k), S12(k), S22(k), S_Total(k)
data_to_save_mixture = hcat(k_values, s11_values, s12_values, s22_values, s_total_values)

# --- Definir Cabecera para el Archivo .dat ---
header_info_mixture = ["Resultados del Factor de Estructura para Mezcla Binaria (RPA con Potencial Yukawa)"]
push!(header_info_mixture, "Generado por: examples/calculate_binary_mixture_SQ.jl")
push!(header_info_mixture, "Parámetros de Esferas Duras (Referencia):")
push!(header_info_mixture, "  σ = $sigma_vector, ϕ_componentes = $phi_vector (ϕ_total = $(sum(phi_vector)))")
push!(header_info_mixture, "Parámetros del Potencial de Yukawa (Perturbación):")
push!(header_info_mixture, "  A_matrix = $A_matrix")
push!(header_info_mixture, "  Z_matrix = $Z_matrix")
push!(header_info_mixture, "Fracciones Molares: x1=$x1, x2=$x2")
push!(header_info_mixture, "Column 1: k (vector de onda adimensional qσ_ref)")
push!(header_info_mixture, "Column 2: S_11(k)")
push!(header_info_mixture, "Column 3: S_12(k)")
push!(header_info_mixture, "Column 4: S_22(k)")
push!(header_info_mixture, "Column 5: S_Total(k)")

# --- Guardar Datos ---
output_directory = "examples_output"
if !isdir(output_directory)
    mkdir(output_directory)
    println("Directorio creado: $output_directory")
end
output_filepath_mixture = joinpath(output_directory, "binary_mixture_Yukawa_RPA_Sk.dat")

try
    save_data(output_filepath_mixture, data_to_save_mixture, header_lines=header_info_mixture, delimiter='\t')
catch e
    println(stderr, "Error al intentar guardar datos. Asegúrate de que la función 'save_data' esté disponible y exportada.")
    showerror(stderr, e)
    println(stderr)
end

println("--- Script de ejemplo para Mezcla Binaria (Yukawa RPA) finalizado ---")

# --- Opcional: Graficar con Gnuplot.jl ---
# println("Intentando graficar S_ij(k) y S_T(k) vs k...")
# try
#     using Gnuplot
#     gp_plot_mix = Gnuplot.GPPlot()
#     Gnuplot.plot(gp_plot_mix, k_values, s11_values, legend="S_11(k)", 韩式="l")
#     Gnuplot.plot!(gp_plot_mix, k_values, s22_values, legend="S_22(k)", 韩式="l")
#     Gnuplot.plot!(gp_plot_mix, k_values, s12_values, legend="S_12(k)", 韩式="l")
#     Gnuplot.plot!(gp_plot_mix, k_values, s_total_values, legend="S_Total(k)", 韩式="l dashtype 2 linewidth 2") # Línea más gruesa y discontinua
    
#     Gnuplot.xlabel(gp_plot_mix, "k (qσ_ref)")
#     Gnuplot.ylabel(gp_plot_mix, "S_ij(k) / S_T(k)")
#     Gnuplot.title(gp_plot_mix, "Factores de Estructura Mezcla Binaria Yukawa (RPA)")
#     Gnuplot.grid(gp_plot_mix, "on")
#     Gnuplot.figure(gp_plot_mix) # Muestra la gráfica
#     # save(term="pngcairo", output=joinpath(output_directory, "binary_mixture_Yukawa_RPA_Sk.png"))
#     println("Gráfica generada (si Gnuplot está configurado).")
# catch e
#     println(stderr, "No se pudo graficar con Gnuplot.jl. Asegúrate de que esté instalado y configurado.")
#     println(stderr, "Error: $e")
# end

