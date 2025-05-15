# AnalyticalStructureFactors.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Riperedo.github.io/AnalyticalStructureFactors.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Riperedo.github.io/AnalyticalStructureFactors.jl/dev/)
[![Build Status](https://github.com/Riperedo/AnalyticalStructureFactors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Riperedo/AnalyticalStructureFactors.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Riperedo/AnalyticalStructureFactors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Riperedo/AnalyticalStructureFactors.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

`AnalyticalStructureFactors.jl`: A Julia package for computing static structure factors $S(q)$ using analytical solutions and well-established approximations for various physical models. Ideal for rapid and precise calculations in condensed matter physics, chemical physics, and materials science.

## Features

* Calculation of $S(q)$ for **Monodisperse Systems**:
    * Hard Spheres (HS) using Percus-Yevick (PY) approximation (`S_HS_PY`).
    * Verlet-Weiss (VW) correction for HS (`S_HS_VW`).
    * Weeks-Chandler-Andersen (WCA) approximation via Blip Function (`S_WCA_blip`).
    * Random Phase Approximation (RPA) for:
        * Square-Well potential (`S_SW_RPA`).
        * Yukawa potential (`S_Yukawa_RPA`).
        * SALR (Short-Range Attraction, Long-Range Repulsion) potential (`S_SALR_RPA`).
* Calculation of partial structure factors $S_{ij}(k)$ for **Mixture Systems**:
    * Hard Sphere mixtures using Baxter's solution for Percus-Yevick (`S_HS_Baxter_mixture`).
    * Verlet-Weiss correction for HS mixtures (`S_HS_VW_mixture`).
    * RPA for mixtures with:
        * Yukawa pair potentials (`S_RPA_mixture_Yukawa`).
        * Square-Well pair potentials (`S_RPA_mixture_SquareWell`).
        * SALR pair potentials (`S_RPA_mixture_SALR`).
* Helper functions for:
    * Effective diameter calculation via Blip function (`blip`).
    * Verlet-Weiss corrections for monodisperse systems (`phi_VW`, `k_VW`).
    * Conversions between packing fraction and number density for mixtures (`phi_to_rho_mixture`, `rho_to_phi_mixture`).
    * Fourier transforms of pair potentials (`betaU_SW`, `betaU_Yukawa`).

## Validation Highlight: Monodisperse Limit of Mixtures

A key validation of the mixture formalism has been successfully performed:
* The partial structure factors $S_{ij}(k)$ for a binary mixture composed of **two identical components** were calculated using the Baxter solution (`S_HS_Baxter_mixture`).
* The total structure factor $S_T(k) = x_1 S_{11}(k) + x_2 S_{22}(k) + 2 \sqrt{x_1 x_2} S_{12}(k)$ was computed from these partials.
* This $S_T(k)$ was found to be **identical** to the $S(k)$ calculated for an equivalent monodisperse system using the Percus-Yevick approximation (`S_HS_PY`) across various total packing fractions.
* This confirms the consistency and correct reduction of the mixture code to the well-established monodisperse limit.

## Installation

The package is registered in the Julia General Registry. You can install it from the Julia REPL:
```julia
import Pkg
Pkg.add("AnalyticalStructureFactors")
If you want to install the latest development version directly from GitHub:import Pkg
Pkg.add(url="[https://github.com/Riperedo/AnalyticalStructureFactors.jl](https://github.com/Riperedo/AnalyticalStructureFactors.jl)")
Quick StartHere's a basic example of how to calculate the static structure factor for Hard Spheres using the Percus-Yevick approximation:using AnalyticalStructureFactors

# Define parameters
ϕ = 0.3  # Volume fraction
k = 2.0  # Dimensionless wavevector (qσ)

# Calculate S(k) for monodisperse system
sq_value = S_HS_PY(ϕ, k)
println("S(k) for HS-PY with ϕ=$ϕ, k=$k is: $sq_value")

# For a range of k values
k_array = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
sq_array = S_HS_PY.(Ref(ϕ), k_array) # Using broadcasting

# Example for WCA (monodisperse)
T_wca = 1.0 # Temperature
sq_wca = S_WCA_blip(ϕ, T_wca, k)
println("S(k) for WCA with ϕ=$ϕ, T=$T_wca, k=$k is: $sq_wca")

# Example for a binary mixture (HS Baxter)
σ_mix = [1.0, 1.0] # Diameters
ϕ_mix = [0.05, 0.05] # Packing fractions of components
ρ_mix = phi_to_rho_mixture(ϕ_mix, σ_mix)
S_ij_matrix = S_HS_Baxter_mixture(σ_mix, ρ_mix, k)
println("S_ij(k) matrix for HS mixture at k=$k: \\n", S_ij_matrix)
```
## Documentation

For more detailed information on the implemented models, API reference, and examples, please refer to the documentation.How to CiteIf you use AnalyticalStructureFactors.jl in your research or software, please cite it.You can use the information in the CITATION.cff file (if available) or cite it as:Ricardo Peredo-Ortiz. (2025). AnalyticalStructureFactors.jl (Version 0.0.1) [Software]. Available at: https://github.com/Riperedo/AnalyticalStructureFactors.jl(Please update the version number and year as appropriate for new releases).ContributingContributions to AnalyticalStructureFactors.jl are welcome! If you find any issues, have suggestions for improvements, or want to add new models, please feel free to:Open an issue on the GitHub repository.