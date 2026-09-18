# Colloidal Electrostatics and DLVO

`AnalyticalStructureFactors.jl` provides dedicated utilities to compute electrostatic properties of electrolyte solutions and convert physical colloidal/protein parameters (surface charge, ionic strength, temperature, dielectric constant) into dimensionless Yukawa parameters for structural solvers.

---

## 1. Fundamental Electrostatic Scales

### Bjerrum Length (`bjerrum_length`)
The Bjerrum length $\ell_B$ is the separation distance at which the electrostatic interaction between two elementary charges equals the thermal energy $k_B T$:
$$\ell_B = \frac{e^2}{4\pi \varepsilon_0 \varepsilon_r k_B T}$$
In water at $T = 298.15\text{ K}$ ($\varepsilon_r \approx 78.4$), $\ell_B \approx 7.14\text{ \AA}$ ($0.714\text{ nm}$).

```julia
using AnalyticalStructureFactors

T = 298.15      # Temperature (K)
eps_r = 78.4    # Relative permittivity of water
l_B = bjerrum_length(T; epsilon_r = eps_r)
println("Bjerrum length = ", l_B, " m")
```

---

### Ionic Strength (`ionic_strength`)
For a solution of electrolytes with molar concentrations $c_i$ and valencies $z_i$:
$$I = \frac{1}{2} \sum_{i} c_i z_i^2$$

```julia
# Example: 10 mM NaCl (1:1 salt)
c_salts = [0.010, 0.010]  # mol / L
valencies = [1, 1]
I_ionic = ionic_strength(c_salts, valencies)
```

---

### Debye Screening Length (`debye_screening_length`)
The characteristic screening distance $\kappa^{-1}$ of the electrostatic double layer:
$$\kappa = \sqrt{8\pi \ell_B N_A I}$$
$$\lambda_D = \kappa^{-1}$$

```julia
# Compute Debye length in meters
kappa_inv = debye_screening_length(I_ionic, T; epsilon_r = 78.4)
println("Debye screening length = ", kappa_inv * 1e9, " nm")
```

---

## 2. DLVO to Yukawa Mapping (`dlvo_yukawa_parameters`)

According to the Derjaguin-Landau-Verwey-Overbeek (DLVO) theory for screened electrostatic repulsion between charged spheres of diameter $\sigma = 2a$ and effective charge $Z_{\text{eff}}$:
$$V_{\text{DLVO}}(r) = \frac{Z_{\text{eff}}^2 e^2}{4\pi \varepsilon_0 \varepsilon_r} \left(\frac{e^{\kappa a}}{1 + \kappa a}\right)^2 \frac{e^{-\kappa r}}{r}, \quad r > \sigma$$

Expressed in the dimensionless Yukawa form:
$$\frac{V(r)}{k_B T} = -K_{\text{DLVO}} \frac{e^{-z(r/\sigma - 1)}}{r/\sigma}$$
where:
$$z = \kappa \sigma$$
$$K_{\text{DLVO}} = -\frac{Z_{\text{eff}}^2 \ell_B}{\sigma (1 + \kappa a)^2} = -\frac{Z_{\text{eff}}^2 \ell_B}{\sigma (1 + z/2)^2}$$

### Practical Example: Protein Electrostatic Repulsion
```julia
# Lysozyme protein in 10 mM NaCl
T = 298.15
sigma = 3.0e-9         # Diameter 3 nm (a = 1.5 nm)
Z_eff = +8.0           # Net charge +8 e
I = 0.010              # 10 mM NaCl

K_dlvo, z_dlvo = dlvo_yukawa_parameters(Z_eff, sigma, I, T; epsilon_r = 78.4)
println("Dimensionless K = ", K_dlvo, ", z = ", z_dlvo)

# Direct input into exact Single-Yukawa MSA solver
phi = 0.08
s_protein = S_Yukawa_MSA(phi, K_dlvo, z_dlvo, 2.0)
```
