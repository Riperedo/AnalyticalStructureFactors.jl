# test/data_test/Menon1991/generate_sticky_limit_comparison.jl
#
# Generates comparison data between:
# 1. Menon et al. (1991) Physical SHS Mapping
# 2. Conventional Sticky Limit (B2-matching / Regnaut & Ravey / Barboy)
# 3. Pure Hard Sphere Reference Limit (PY)
# 4. Monte Carlo simulation data (Huang et al. 1984)

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__
QR_grid = range(0.001, 6.0, length=1000)

# ------------------------------------------------------------------------------
# 1. Figure 3: Homogeneous Fluid Phase
# ------------------------------------------------------------------------------
phi3 = 0.07
a_sigma3 = 1.1
u0_3 = 0.92

# Menon mapping
eps3 = 1.0 - 1.0 / a_sigma3
eta_menon3 = phi3 / ((1.0 - eps3)^3)
tau_menon3 = exp(-u0_3) / (12.0 * eps3)

# Conventional B2-matching sticky limit (Barboy 1974 / Regnaut & Ravey 1989 / Menon Eq. 6)
tau_conv3 = 1.0 / (4.0 * (exp(u0_3) - 1.0) * (a_sigma3^3 - 1.0))
eta_conv3 = phi3

fig3_comp_file = joinpath(out_dir, "menon_fig3_sticky_limit_comparison.dat")
open(fig3_comp_file, "w") do io
    println(io, "# QR  kappa  S_Menon_SHS  S_Conv_Sticky_B2  S_Pure_HS_phi  S_Pure_HS_eta")
    for QR in QR_grid
        kappa = 2.0 * QR
        s_menon = S_Baxter_SHS(eta_menon3, tau_menon3, kappa)
        s_conv  = S_Baxter_SHS(eta_conv3, tau_conv3, kappa)
        s_hs_phi = S_HS_PY(phi3, kappa)
        s_hs_eta = S_HS_PY(eta_menon3, kappa)
        @printf(io, "%.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n", QR, kappa, s_menon, s_conv, s_hs_phi, s_hs_eta)
    end
end
println("Saved Fig. 3 sticky limit comparison -> $fig3_comp_file")

# ------------------------------------------------------------------------------
# 2. Figure 4: Two-Phase Coexistence (Gas Phase)
# ------------------------------------------------------------------------------
phi4 = 0.07
a_sigma4 = 1.02
u0_4 = 3.83

# Menon mapping (Gas phase equilibrium density eta_g ≈ 0.0581)
eps4 = 1.0 - 1.0 / a_sigma4
eta_menon4_gas = 0.0581
tau_menon4 = exp(-u0_4) / (12.0 * eps4)

# Conventional B2-matching sticky limit
tau_conv4 = 1.0 / (4.0 * (exp(u0_4) - 1.0) * (a_sigma4^3 - 1.0))
eta_conv4_gas = 0.0581

fig4_comp_file = joinpath(out_dir, "menon_fig4_sticky_limit_comparison.dat")
open(fig4_comp_file, "w") do io
    println(io, "# QR  kappa  S_Menon_Gas_SHS  S_Conv_Gas_Sticky_B2  S_Pure_HS_eta_g  S_Pure_HS_phi")
    for QR in QR_grid
        kappa = 2.0 * QR
        s_menon = S_Baxter_SHS(eta_menon4_gas, tau_menon4, kappa)
        s_conv  = S_Baxter_SHS(eta_conv4_gas, tau_conv4, kappa)
        s_hs_gas = S_HS_PY(eta_menon4_gas, kappa)
        s_hs_phi = S_HS_PY(phi4, kappa)
        @printf(io, "%.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n", QR, kappa, s_menon, s_conv, s_hs_gas, s_hs_phi)
    end
end
println("Saved Fig. 4 sticky limit comparison -> $fig4_comp_file")
