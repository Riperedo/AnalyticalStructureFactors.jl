# src/AnalyticalStructureFactors.jl

module AnalyticalStructureFactors

# --- Dependencies ---
using LinearAlgebra      # Matrix operations and companion matrix eigenvalues
using SpecialFunctions   # Gamma and special functions

# --- Source Code Inclusions ---

# General utilities, corrections, and conversions
include("utils.jl") # blip, phi_VW, k_VW, phi_to_rho_mixture, rho_to_phi_mixture, save_data

# Fourier transforms of pair potentials (SW, Yukawa)
include("potentials.jl")

# Colloidal electrostatics and DLVO / GOCM parameters
include("colloidal_potentials.jl")

# 3D Spherical Fourier transforms and S(q) <-> g(r) inversion tools
include("fourier_tools.jl")

# --- Single-Component (Monodisperse) Structure Factors ---
include(joinpath("models", "hardsphere.jl")) # S_HS_PY, S_HS_VW, S_WCA_blip
include(joinpath("models", "sticky_hardsphere.jl")) # S_Baxter_SHS, IS_Baxter_SHS, S_Menon_SHS (Baxter 1968, Menon et al. 1991)
include(joinpath("models", "squarewell.jl")) # S_SW_RPA
include(joinpath("models", "yukawa.jl"))     # S_Yukawa_RPA
include(joinpath("models", "yukawa_msa.jl")) # S_Yukawa_MSA, IS_Yukawa_MSA, C_Yukawa_MSA (Exact Baxter-Høye-Blum solution)
include(joinpath("models", "SALR.jl"))       # S_SALR_RPA
include(joinpath("models", "SALR_msa.jl"))   # S_SALR_MSA, IS_SALR_MSA, C_SALR_MSA (Exact Two-Yukawa MSA)
include(joinpath("models", "fractal.jl"))    # S_Fractal_Chen (Chen-Teixeira mass-fractal aggregate model)

# --- Particle Form Factors P(q) and Decoupling Factor β(q) ---
include(joinpath("models", "formfactors.jl")) # P_sphere, P_coreshell, P_multishell, P_ellipsoid, beta_ellipsoid, intensity_decoupled

# --- Multicomponent Mixture Structure Factors ---
include(joinpath("models", "asakura_oosawa.jl")) # S_AO_cc, S_AO_cp, S_AO_pp, S_AO_matrix, S_AO_eff_cc, S_AO_eff_pp (Colloid-Polymer Mixtures)
include(joinpath("models", "mixtures", "hardsphere_mixture.jl")) # S_HS_Baxter_mixture, IS_HS_Baxter_mixture
include(joinpath("models", "mixtures", "vw_mixture.jl"))         # S_HS_VW_mixture, IS_HS_VW_mixture
include(joinpath("models", "mixtures", "yukawa_rpa_mixture.jl")) # S_RPA_mixture_Yukawa
include(joinpath("models", "mixtures", "squarewell_rpa_mixture.jl")) # S_RPA_mixture_SquareWell
include(joinpath("models", "mixtures", "SALR_rpa_mixture.jl"))       # S_RPA_mixture_SALR


# ==============================================================================
# PUBLIC API EXPORTS
# ==============================================================================

# From utils.jl
export blip
export phi_VW, k_VW
export phi_to_rho_mixture, rho_to_phi_mixture
export save_data

# From potentials.jl
export betaU_SW, betaU_Yukawa

# From colloidal_potentials.jl
export ionic_strength, bjerrum_length, debye_screening_length, dlvo_yukawa_parameters

# From fourier_tools.jl
export sk_to_gr, gr_to_sk

# From models/hardsphere.jl
export S_HS_PY, C_HS_PY, IS_HS_PY
export S_HS_VW, C_HS_VW, IS_HS_VW
export S_WCA_blip, C_WCA_blip, IS_WCA_blip

# From models/sticky_hardsphere.jl (Baxter 1968 & Menon et al. 1991)
export S_Baxter_SHS, IS_Baxter_SHS, S_Menon_SHS
export baxter_shs_lambda, baxter_shs_critical_tau

# From models/squarewell.jl
export S_SW_RPA

# From models/yukawa.jl
export S_Yukawa_RPA

# From models/yukawa_msa.jl (Exact MSA)
export S_Yukawa_MSA, IS_Yukawa_MSA, C_Yukawa_MSA
export yukawa_msa_coefficients, yukawa_msa_roots, yukawa_msa_beta_roots

# From models/SALR.jl & models/SALR_msa.jl
export S_SALR_RPA
export S_SALR_MSA, IS_SALR_MSA, C_SALR_MSA

# From models/fractal.jl
export S_Fractal_Chen

# From models/formfactors.jl
export P_sphere, P_coreshell, P_multishell, P_ellipsoid, beta_ellipsoid, effective_radius_ellipsoid, intensity_decoupled

# From models/asakura_oosawa.jl (Colloid-Polymer Mixtures)
export S_AO_cc, S_AO_cp, S_AO_pp, S_AO_matrix, S_AO_eff_cc, S_AO_eff_pp
export C_AO_cc, C_AO_cp, C_AO_star, C_AO_eff_cc, C_AO_HS
export free_volume_fraction_AO, spinodal_AO, polymer_reservoir_to_system_density_AO, polymer_system_to_reservoir_density_AO

# Multicomponent Mixtures
export Qk_mixture, S_HS_Baxter_mixture, IS_HS_Baxter_mixture
export VW_correction_mixture, S_HS_VW_mixture, IS_HS_VW_mixture
export S_RPA_mixture_Yukawa
export S_RPA_mixture_SquareWell
export S_RPA_mixture_SALR

end # module
