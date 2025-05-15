# AnalyticalStructureFactors.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Riperedo.github.io/AnalyticalStructureFactors.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Riperedo.github.io/AnalyticalStructureFactors.jl/dev/)
[![Build Status](https://github.com/Riperedo/AnalyticalStructureFactors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Riperedo/AnalyticalStructureFactors.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Riperedo/AnalyticalStructureFactors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Riperedo/AnalyticalStructureFactors.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

`AnalyticalStructureFactors.jl`: A Julia package for computing static structure factors $S(q)$ using analytical solutions and well-established approximations for various physical models. Ideal for rapid and precise calculations in condensed matter physics, chemical physics, and materials science.

## Features

* Calculation of $S(q)$ for **Hard Spheres (HS)** using:
    * Percus-Yevick (PY) approximation (`S_HS_PY`).
    * Verlet-Weiss (VW) correction (`S_HS_VW`).
* **Weeks-Chandler-Andersen (WCA)** approximation via Blip Function for soft potentials, mapping to an effective HS-VW system (`S_WCA_blip`).
* **Random Phase Approximation (RPA)** for various interaction potentials:
    * Square-Well potential (`S_SW_RPA`).
    * Yukawa potential (`S_Yukawa_RPA`).
    * SALR (Short-Range Attraction, Long-Range Repulsion) potential (`S_SALR_RPA`).
* Helper functions for:
    * Effective diameter calculation via Blip function (`blip`).
    * Verlet-Weiss corrections (`phi_VW`, `k_VW`).
    * Fourier transforms of pair potentials (`betaU_SW`, `betaU_Yukawa`).

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

# Calculate S(k)
sq_value = S_HS_PY(ϕ, k)
println("S(k) for HS-PY with ϕ=$ϕ, k=$k is: $sq_value")

# For a range of k values
k_array = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
sq_array = S_HS_PY.(Ref(ϕ), k_array) # Using broadcasting

# Example using the "function factory" pattern for S_HS_PY
s_hs_py_func = S_HS_PY(ϕ) # Returns a function f(k)
sq_value_factory = s_hs_py_func(k)
println("S(k) from factory with ϕ=$ϕ, k=$k is: $sq_value_factory")
sq_array_factory = s_hs_py_func.(k_array)

# Example for WCA
T_wca = 1.0 # Temperature
sq_wca = S_WCA_blip(ϕ, T_wca, k)
println("S(k) for WCA with ϕ=$ϕ, T=$T_wca, k=$k is: $sq_wca")

# Example for Yukawa
beta_yuk = 1.0 # Inverse temperature
z_yuk = 2.0 # inv screening length
sq_yuk = S_Yukawa_RPA(ϕ, beta_yuk, z, k)
println("S(k) for Yukaya RPA with ϕ=$ϕ, β =$beta_yuk, z = $z_yuk k=$k is: $sq_yuk")
```

## Documentation

For more detailed information on the implemented models, API reference, and examples, please refer to the documentation.How to CiteIf you use AnalyticalStructureFactors.jl in your research or software, please cite it.You can use the information in the CITATION.cff file (if available) or cite it as:Ricardo Peredo-Ortiz. (2025). AnalyticalStructureFactors.jl (Version 0.0.1) [Software]. Available at: https://github.com/Riperedo/AnalyticalStructureFactors.jl.ContributingContributions to AnalyticalStructureFactors.jl are welcome! If you find any issues, have suggestions for improvements, or want to add new models, please feel free to:Open an issue on the GitHub repository.Submit a pull request with your changes.Please follow the contribution guidelines if available.LicenseAnalyticalStructureFactors.jl is licensed under the [MIT License](https://opensource.org/license/mit).
