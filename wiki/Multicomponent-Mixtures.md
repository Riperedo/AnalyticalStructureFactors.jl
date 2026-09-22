# Multicomponent Mixtures

`AnalyticalStructureFactors.jl` provides analytical solutions for $p$-component polydisperse mixtures of spherical particles.

---

## 1. Conventions and Matrix Definitions

For a mixture of $p$ species with hard-core diameters $\boldsymbol{\sigma} = (\sigma_1, \sigma_2, \dots, \sigma_p)$, number densities $\boldsymbol{\rho} = (\rho_1, \rho_2, \dots, \rho_p)$, and volume fractions $\boldsymbol{\phi} = (\phi_1, \phi_2, \dots, \phi_p)$ where $\phi_i = \frac{\pi}{6}\rho_i \sigma_i^3$:

- The **Ashcroft-Langreth partial structure factor matrix** $\mathbf{S}(k) \in \mathbb{R}^{p \times p}$ is defined as:
  $$S_{ij}(k) = \frac{1}{\sqrt{\rho_i \rho_j}} \left[ \left(\mathbf{I} - \mathbf{C}(k)\right)^{-1} \right]_{ij}$$
  where $\mathbf{C}(k)_{ij} = \sqrt{\rho_i \rho_j} \hat{c}_{ij}(k)$.

- Density conversion helpers:
  ```julia
  # Convert volume fractions to number densities
  rho_vec = phi_to_rho_mixture(phi_vec, sigma_vec)

  # Convert number densities to volume fractions
  phi_vec = rho_to_phi_mixture(rho_vec, sigma_vec)
  ```

---

## 2. Baxter Hard Sphere Mixture (`S_HS_Baxter_mixture`)

The exact analytical solution of the multicomponent Ornstein-Zernike equation with Percus-Yevick closure was derived by Baxter using matrix Wiener-Hopf factorization:
$$\hat{q}_{ij}(k) = \delta_{ij} - 2\pi \sqrt{\rho_i \rho_j} \int_0^{\infty} e^{i k r} Q_{ij}(r) \, dr$$

### Usage Example:
```julia
using AnalyticalStructureFactors

# Binary mixture: species 1 (diameter 1.0, phi 0.2), species 2 (diameter 1.5, phi 0.1)
sigma_vec = [1.0, 1.5]
phi_vec = [0.2, 0.1]
k = 2.0

# Returns the p x p partial structure factor matrix S_ij(k)
S_matrix = S_HS_Baxter_mixture(phi_vec, sigma_vec, k)
println("S_11 = ", S_matrix[1,1], ", S_12 = ", S_matrix[1,2], ", S_22 = ", S_matrix[2,2])

# Inverse structure factor matrix
inv_S_matrix = IS_HS_Baxter_mixture(phi_vec, sigma_vec, k)
```

---

## 3. Verlet-Weis Mixture Correction (`S_HS_VW_mixture`)

The semi-empirical Verlet-Weis correction generalizes to mixtures by scaling each component's volume fraction $\phi_i \to \phi_{w,i}$ and diameter $\sigma_i \to \sigma_{w,i}$, providing improved structure factor predictions near the contact peaks:

```julia
S_vw_matrix = S_HS_VW_mixture(phi_vec, sigma_vec, k)
```

---

## 4. RPA Potential Mixtures

Under the Random Phase Approximation (RPA), attractive and repulsive continuous potential tails $\hat{U}_{ij}(k)$ are added linearly to the direct correlation function matrix:
$$\mathbf{C}_{\text{RPA}}(k) = \mathbf{C}_{\text{HS}}(k) - \beta \hat{\mathbf{U}}(k)$$

### Yukawa Mixture (`S_RPA_mixture_Yukawa`)
Cross-interaction contact energies $\mathbf{K} \in \mathbb{R}^{p \times p}$ and screening parameters $\mathbf{Z} \in \mathbb{R}^{p \times p}$:
```julia
K_mat = [1.0 0.5; 0.5 1.2]   # Contact energies K_ij = \beta \epsilon_ij
Z_mat = [1.8 1.8; 1.8 1.8]   # Inverse screening lengths z_ij

S_yukawa_mix = S_RPA_mixture_Yukawa(phi_vec, sigma_vec, K_mat, Z_mat, k)
```

### Square-Well Mixture (`S_RPA_mixture_SquareWell`)
Cross-interaction well depths $\boldsymbol{\epsilon} \in \mathbb{R}^{p \times p}$ and ranges $\boldsymbol{\lambda} \in \mathbb{R}^{p \times p}$:
```julia
eps_mat = [1.0 0.8; 0.8 1.5]   # \epsilon_ij / (k_B * T)
lambda_mat = [1.2 1.2; 1.2 1.2] # \lambda_ij (well width factor)

S_sw_mix = S_RPA_mixture_SquareWell(phi_vec, sigma_vec, lambda_mat, eps_mat, k)
```

### SALR Mixture (`S_RPA_mixture_SALR`)
Competing two-Yukawa interactions across multiple components:
```julia
K1_mat = [6.0 4.0; 4.0 6.0]    # Short-range attraction
Z1_mat = [10.0 10.0; 10.0 10.0]
K2_mat = [-1.0 -0.5; -0.5 -1.0] # Long-range repulsion
Z2_mat = [0.5 0.5; 0.5 0.5]

S_salr_mix = S_RPA_mixture_SALR(phi_vec, sigma_vec, K1_mat, Z1_mat, K2_mat, Z2_mat, k)
```

---

## 5. Asakura-Oosawa (AO) Colloid-Polymer Model

The **Asakura-Oosawa (AO) model** describes mixtures of hard-sphere colloids (species $c$, packing fraction $\eta_c$) and non-adsorbing ideal polymer coils (species $p$, packing fraction $\eta_p$) with size ratio $q = R_p / R_c$.

Based on the density functional theory and Percus-Yevick analytical solution by **Schmidt et al. (2002)** (*J. Phys.: Condens. Matter* **14**, 9353):
- **Free Volume Fraction**: $\alpha(\eta_c, q) = (1-\eta_c)\exp(-[A\gamma + B\gamma^2 + C\gamma^3])$ where $\gamma = \eta_c / (1 - \eta_c)$.
- **Colloid-Colloid Direct Correlation**: $C_{cc}(k) = C_{\text{HS}}(\eta_c, k) + \eta_p C^*(\eta_c, q, k)$.
- **Cross Direct Correlation**: $C_{cp}(k)$ across the core and overlap shells.
- **Partial Structure Factors**: $S_{cc}(k), S_{cp}(k), S_{pp}(k)$.
- **Effective One-Component Structure Factor**: $S_{cc}^{\text{eff}}(k) = \frac{1}{1 - \rho_c C_{cc}^{\text{eff}}(k)}$.

### Usage Example:
```julia
using AnalyticalStructureFactors

eta_c = 0.20   # Colloid volume fraction
eta_p = 0.10   # Polymer volume fraction
q = 0.6        # Size ratio Rp / Rc
k = 2.5        # Dimensionless wavevector q * sigma_c

# Partial structure factors
s_cc = S_AO_cc(eta_c, eta_p, q, k)
s_cp = S_AO_cp(eta_c, eta_p, q, k)
s_pp = S_AO_pp(eta_c, eta_p, q, k)

# Full 2x2 Ashcroft-Langreth structure factor matrix
S_matrix = S_AO_matrix(eta_c, eta_p, q, k; VW=false)

# Effective one-component colloid structure factor
s_eff_cc = S_AO_eff_cc(eta_c, eta_p, q, k)

# Spinodal curve polymer packing fraction
eta_p_spinodal = spinodal_AO(eta_c, q)

# Free volume fraction alpha(eta_c, q) and reservoir conversion
alpha = free_volume_fraction_AO(eta_c, q)
eta_p_res = polymer_system_to_reservoir_density_AO(eta_p, eta_c, q)
```

---

## 6. Total Coherent Scattering from Mixtures

For polydisperse mixtures with scattering length densities $b_i$, the total coherent scattering intensity $I(q)$ is:
$$I(q) = \sum_{i=1}^p \sum_{j=1}^p \sqrt{\rho_i \rho_j} \, b_i(q) b_j(q) S_{ij}(q)$$

