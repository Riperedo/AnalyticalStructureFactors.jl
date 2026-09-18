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

## 5. Total Coherent Scattering from Mixtures

For polydisperse mixtures with scattering length densities $b_i$, the total coherent scattering intensity $I(q)$ is:
$$I(q) = \sum_{i=1}^p \sum_{j=1}^p \sqrt{\rho_i \rho_j} \, b_i(q) b_j(q) S_{ij}(q)$$
