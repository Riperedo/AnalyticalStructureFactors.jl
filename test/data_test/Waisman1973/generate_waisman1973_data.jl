# test/data_test/Waisman1973/generate_waisman1973_data.jl
#
# Generates the radial distribution function g(r) for high-density hard spheres (eta = 0.49)
# using AnalyticalStructureFactors.jl to validate against Table 1 of Waisman (1973)
# (Mol. Phys., 25(1), 45-48).

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__

# Physical parameters for hard spheres at high density
eta = 0.49
rho = eta / (π / 6.0)

# Waisman (1973) Yukawa MSA parameters matching Carnahan-Starling EOS
K = 1.20
z = 14.5

# High-resolution wavevector grid for accurate Fourier inversion
q_grid = collect(range(1e-4, 200.0, length=10000))

println("Evaluating static structure factors S(q) using AnalyticalStructureFactors.jl...")
Sq_waisman = [S_Yukawa_MSA(eta, K, z, q) for q in q_grid]
Sq_py = [S_HS_PY(eta, q) for q in q_grid]

# 1. Discrete evaluation at Table 1 distance points
r_table = [1.00, 1.06, 1.12, 1.18, 1.24, 1.30, 1.36, 1.42, 1.48, 1.54, 1.60]

# Continuous radial grid for smooth plotting
r_dense = collect(range(1.00, 3.00, length=500))

println("Performing 3D isotropic Fourier inversion sk_to_gr...")
gr_dense_waisman = sk_to_gr(q_grid, Sq_waisman, r_dense, rho)
gr_dense_py = sk_to_gr(q_grid, Sq_py, r_dense, rho)

# Analytical contact values g(1^+)
g_contact_CS = (1.0 - 0.5 * eta) / ((1.0 - eta)^3) # Carnahan-Starling -> 5.6916
g_contact_PY = (1.0 + 0.5 * eta) / ((1.0 - eta)^2) # Percus-Yevick -> 4.7866

# Set contact point in dense grid (since Fourier transform of discontinuous function yields Gibbs midpoint)
gr_dense_waisman[1] = g_contact_CS
gr_dense_py[1] = g_contact_PY

# Discrete evaluated points
gr_tab_waisman = sk_to_gr(q_grid, Sq_waisman, r_table, rho)
gr_tab_py = sk_to_gr(q_grid, Sq_py, r_table, rho)
gr_tab_waisman[1] = g_contact_CS
gr_tab_py[1] = g_contact_PY

# Write dense continuous curves
dense_file = joinpath(out_dir, "waisman_continuous.dat")
open(dense_file, "w") do io
    println(io, "# r/d  g_Waisman_MSA(r)  g_PercusYevick(r)")
    for i in 1:length(r_dense)
        @printf(io, "%.5f  %.5f  %.5f\n", r_dense[i], gr_dense_waisman[i], gr_dense_py[i])
    end
end
println("Generated continuous curve data -> $dense_file")

# Write tabular comparison
table_file = joinpath(out_dir, "waisman_table_evaluated.dat")
open(table_file, "w") do io
    println(io, "# x(r/d)  g_MC_Ref9  g_PY_tab  g_PY_calc  g_Waisman_tab  g_Waisman_calc")
    # Reference values from Waisman (1973) Table 1
    mc_vals = [5.75, 3.40, 2.18, 1.48, 1.08, 0.84, 0.66, 0.58, 0.57, 0.61, 0.67]
    py_tab = [4.79, 3.45, 2.43, 1.67, 1.15, 0.81, 0.62, 0.54, 0.53, 0.58, 0.65]
    tw_tab = [5.69, 3.39, 2.28, 1.58, 1.11, 0.81, 0.64, 0.56, 0.56, 0.60, 0.67]
    for i in 1:length(r_table)
        @printf(io, "%.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n",
                r_table[i], mc_vals[i], py_tab[i], gr_tab_py[i], tw_tab[i], gr_tab_waisman[i])
    end
end
println("Generated discrete table comparison -> $table_file")
