# Welcome to AnalyticalStructureFactors.jl Wiki

**AnalyticalStructureFactors.jl** is a high-performance, purely analytical Julia package designed for calculating static structure factors $S(q)$, direct correlation functions $C(q)$, particle form factors $P(q)$, decoupling factors $\beta(q)$, and small-angle scattering intensities $I(q)$ for monodisperse colloids, globular proteins, fractal aggregates, and polydisperse multicomponent mixtures.

---

## Key Highlights

- **Purely Analytical Solutions**: Eliminates iterative numerical integral equation solvers (such as hypernetted-chain HNC or Percus-Yevick Picard iterations) for rapid, deterministic microstructural calculations ($< 0.1$ ms per evaluation).
- **Exact Mean Spherical Approximation (MSA)**:
  - **Single-Yukawa MSA**: Exact Baxter-Høye-Blum factorization and analytic quartic root solver ($S(q), c(r), S(0)$).
  - **Two-Yukawa / SALR MSA**: Exact algebraic reduction to a 1D scalar root-finder on the Baxter parameter $d_2$ with automated physical core-penalty branch selection.
- **Multicomponent Mixtures**: Full support for $p$-component Baxter hard-sphere mixtures, Verlet-Weis semi-empirical corrections, and Random Phase Approximation (RPA) for Yukawa, Square-Well, and SALR mixtures.
- **Form Factors & Decoupling Approximation**: Spheres, Core-Shell, Multi-Shell, and Prolate/Oblate Ellipsoids with orientational decoupling factor $\beta(q)$ for polydisperse/non-spherical scattering $I(q) = I_0 P(q) [1 + \beta(q)(S(q)-1)]$.
- **Fractal Aggregates & Gels**: Sow-Hsin Chen / Teixeira mass-fractal aggregate model $S_{\text{fractal}}(q)$ with cutoff length $\xi$, fractal dimension $D_f$, and monomer radius $R_0$.
- **Colloidal Electrostatics & DLVO**: Automated computation of ionic strength, Debye screening length $\kappa^{-1}$, Bjerrum length $\ell_B$, and mapping from surface charge/potential to Yukawa parameters.
- **3D Isotropic Fourier Inversion**: Robust sine-FFT routines converting $S(q) \leftrightarrow g(r)$ with strict asymptotic damping.
- **Validated Against Literature**: Tested and benchmarked against Waisman (1973), Cummings & Smith (1979), and Liu, Chen, & Chen (2005).

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Riperedo/AnalyticalStructureFactors.jl.git")
```

Or via the Julia package manager REPL (press `]`):
```julia
pkg> add https://github.com/Riperedo/AnalyticalStructureFactors.jl.git
```

---

## Quick Start Example

```julia
using AnalyticalStructureFactors

# 1. Monodisperse Hard Sphere Structure Factor (Percus-Yevick)
ϕ = 0.30      # Volume fraction
k = 2.5       # Dimensionless wavevector (q * σ)
S_hs = S_HS_PY(ϕ, k)
println("Hard Sphere S(k = 2.5): ", S_hs)

# 2. Exact SALR Two-Yukawa MSA Structure Factor
# Parameters: ϕ, K1 (attraction), Z1, K2 (repulsion), Z2, k
K1, Z1 = 6.0, 10.0   # Short-range attraction
K2, Z2 = -1.0, 0.5   # Long-range repulsion
S_salr = S_SALR_MSA(ϕ, K1, Z1, K2, Z2, k)
println("SALR MSA S(k = 2.5): ", S_salr)

# 3. Particle Form Factor and Decoupled Scattering Intensity
R = 15.0      # Particle radius (Å)
q = 0.05      # Scattering wavevector (1/Å)
P_q = P_sphere(q, R)
I_q = intensity_decoupled(q, S_salr, P_q; I0 = 1.0)
println("Scattering Intensity I(q = 0.05): ", I_q)
```

---

## Wiki Table of Contents

| Section | Description |
| :--- | :--- |
| **[Monodisperse Models](Monodisperse-Models)** | Hard Spheres (PY, Verlet-Weis, WCA), Square-Well RPA, Yukawa RPA, SALR RPA |
| **[Mean Spherical Approximation](Mean-Spherical-Approximation)** | Exact Single-Yukawa & Two-Yukawa (SALR) MSA analytical solvers |
| **[Multicomponent Mixtures](Multicomponent-Mixtures)** | $p$-component Baxter HS, Verlet-Weis corrections, Yukawa/SW/SALR RPA mixtures |
| **[Form Factors & Scattering](Form-Factors-and-Scattering)** | Spheres, core-shell, multi-shell, ellipsoids, orientational decoupling approximation |
| **[Fractal Structures](Fractal-Structures)** | Chen-Teixeira mass-fractal model for aggregates, gels, and cluster networks |
| **[Colloidal Electrostatics & DLVO](Colloidal-Electrostatics-and-DLVO)** | Debye length, Bjerrum length, charge renormalization, DLVO-to-Yukawa mapping |
| **[Fourier Transform Tools](Fourier-Transform-Tools)** | 3D isotropic Fourier inversion routines $S(q) \leftrightarrow g(r)$ |
| **[Benchmark Validation](Benchmark-Validation)** | Validation suites against Waisman (1973), Cummings & Smith (1979), and Liu et al. (2005) |
| **[API Reference](API-Reference)** | Comprehensive reference of all exported functions, signatures, and types |

---

## Citation

If you use `AnalyticalStructureFactors.jl` in academic research, please cite:

```bibtex
@software{PeredoOrtiz_AnalyticalStructureFactors_2025,
  author = {Peredo-Ortiz, Ricardo and Development Team},
  title = {{AnalyticalStructureFactors.jl: Analytical Static Structure Factors and Scattering Functions in Julia}},
  url = {https://github.com/Riperedo/AnalyticalStructureFactors.jl},
  year = {2025}
}
```
