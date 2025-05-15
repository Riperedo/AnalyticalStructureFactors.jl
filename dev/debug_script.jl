# dev/debug_script.jl
# Para ejecutar desde la raíz del proyecto AnalyticalStructureFactors.jl:
# julia --project dev/debug_script.jl

using Revise # Carga Revise primero
using AnalyticalStructureFactors # Carga tu paquete
# DelimitedFiles es necesario para save_data, pero save_data ya lo usa internamente.
# Si llamaras a writedlm/readdlm directamente aquí, necesitarías `using DelimitedFiles`.

println("--- Iniciando Script de Depuración y Comparación ---")

# --- Parámetros ---
# Monodisperso
phi_mono = 0.6
sigma_mono = 1.0 # Asumimos diámetro 1 para S_HS_PY

# Mezcla Binaria Equivalente
# Componentes idénticos, la suma de sus fracciones de volumen debe ser phi_mono
# y sus diámetros deben ser sigma_mono.
sigma_mixture = [sigma_mono, sigma_mono]   # Ambos componentes tienen diámetro 1.0
phi_mixture = [phi_mono / 2.0, phi_mono / 2.0] # Ej: [0.05, 0.05] para que sumen 0.1

# Rango de k común
k_range = collect(0.01:0.05:15.0) # Evitar k=0 exacto si causa problemas en alguna función,
                                 # o asegurar que la función lo maneje.

# --- Cálculo para Sistema Monodisperso ---
println("Calculando S(k) para sistema monodisperso (PY)...")
S_mono_values = S_HS_PY.(Ref(phi_mono), k_range)
println("Cálculo monodisperso completado. Primeros valores de S_mono: ", S_mono_values[1:min(5, end)])

# --- Cálculo para Mezcla Binaria Equivalente ---
println("\nCalculando S_ij(k) y S_Total(k) para mezcla binaria equivalente (Baxter)...")
# Convertir fracciones de volumen de la mezcla a densidades numéricas
rho_mixture = phi_to_rho_mixture(phi_mixture, sigma_mixture)

# Calcular fracciones molares x_i = ρ_i / ρ_total
rho_total_mixture = sum(rho_mixture)
x_mixture = rho_mixture ./ rho_total_mixture
x1 = x_mixture[1]
x2 = x_mixture[2]

S_mixture_matrices = S_HS_Baxter_mixture(sigma_mixture, rho_mixture, k_range)

# Extraer los componentes S11, S12, S22
S11_mixture_values = [S_mat[1,1] for S_mat in S_mixture_matrices]
S12_mixture_values = [S_mat[1,2] for S_mat in S_mixture_matrices] # S12 = S21
S22_mixture_values = [S_mat[2,2] for S_mat in S_mixture_matrices]

# Calcular el Factor de Estructura Total para la mezcla
# S_T(k) = x₁S₁₁(k) + x₂S₂₂(k) + 2√(x₁x₂)S₁₂(k)
S_total_mixture_values = (x1 .* S11_mixture_values) .+ 
                         (x2 .* S22_mixture_values) .+ 
                         (2 * sqrt(x1 * x2) .* S12_mixture_values)

println("Cálculo de mezcla completado.")
println("Primeros valores de S11_mix: ", S11_mixture_values[1:min(5, end)])
println("Primeros valores de S12_mix: ", S12_mixture_values[1:min(5, end)])
println("Primeros valores de S22_mix: ", S22_mixture_values[1:min(5, end)])
println("Primeros valores de S_Total_mix: ", S_total_mixture_values[1:min(5, end)])


# --- Guardar Datos para Comparación ---
# Columnas: k, S_mono(k), S11_mix(k), S12_mix(k), S22_mix(k), S_Total_mix(k)
data_to_save = hcat(k_range, S_mono_values, S11_mixture_values, S12_mixture_values, S22_mixture_values, S_total_mixture_values)

output_filepath = joinpath("dev", "comparison_mono_vs_mixture.dat")
header_lines = [
    "Comparison of Monodisperse PY S(k) vs Baxter Mixture S_ij(k) and S_Total(k) for equivalent system.",
    "Monodisperse: phi = $phi_mono, sigma = $sigma_mono",
    "Mixture: sigma = $sigma_mixture, phi_components = $phi_mixture (total_phi = $(sum(phi_mixture)))",
    "Mole fractions: x1 = $x1, x2 = $x2",
    "Column 1: k (wavevector)",
    "Column 2: S_monodisperse_PY(k)",
    "Column 3: S_11_mixture(k)",
    "Column 4: S_12_mixture(k)",
    "Column 5: S_22_mixture(k)",
    "Column 6: S_Total_mixture(k) = x1*S11 + x2*S22 + 2*sqrt(x1*x2)*S12"
]

save_data(output_filepath, data_to_save, header_lines=header_lines, delimiter='\t')

println("\n--- Script de Depuración Finalizado ---")
