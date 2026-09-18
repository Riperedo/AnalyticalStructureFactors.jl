# API Reference

Complete reference of all functions exported by `AnalyticalStructureFactors.jl`.

---

## 1. Monodisperse Hard Spheres & Potentials

### `S_HS_PY(ϕ, k)`
- **Description**: Static structure factor for hard spheres using the Percus-Yevick analytical solution.
- **Arguments**:
  - `ϕ::Real`: Volume fraction ($\phi \in (0, 1)$).
  - `k::Real`: Dimensionless wavevector ($k = q\sigma$).
- **Returns**: `Float64` ($S(k)$).

### `C_HS_PY(ϕ, k)`
- **Description**: 3D Fourier transform of direct correlation function $C(k) = \rho \hat{c}(k)$ for Percus-Yevick hard spheres.

### `IS_HS_PY(ϕ, k)`
- **Description**: Inverse structure factor $1 / S(k)$ for Percus-Yevick hard spheres.

### `S_HS_VW(ϕ, k)`
- **Description**: Hard sphere structure factor with Verlet-Weis semi-empirical correction.

### `S_WCA_blip(ϕ, k, T_star)`
- **Description**: Structure factor for Weeks-Chandler-Andersen (WCA) repulsive core with temperature-dependent blip function.

---

## 2. Monodisperse RPA & Exact MSA Solvers

### `S_SW_RPA(ϕ, λ, ϵ_kT, k)`
- **Description**: Square-Well structure factor under Random Phase Approximation.
- **Arguments**:
  - `ϕ::Real`: Hard core volume fraction.
  - `λ::Real`: Well width relative to diameter ($\lambda > 1$).
  - `ϵ_kT::Real`: Dimensionless well depth $\epsilon / (k_B T)$.
  - `k::Real`: Dimensionless wavevector ($q\sigma$).

### `S_Yukawa_RPA(ϕ, K, z, k)`
- **Description**: Single-Yukawa structure factor under Random Phase Approximation.
- **Arguments**:
  - `ϕ::Real`: Volume fraction.
  - `K::Real`: Contact energy $\beta\epsilon$.
  - `z::Real`: Dimensionless inverse screening length $\kappa\sigma$.
  - `k::Real`: Dimensionless wavevector ($q\sigma$).

### `S_Yukawa_MSA(ϕ, K, z, k)`
- **Description**: Exact analytical Mean Spherical Approximation structure factor for Single-Yukawa fluid via Baxter $Q$-factorization.

### `S_SALR_RPA(ϕ, K1, z1, K2, z2, k)`
- **Description**: Two-Yukawa / SALR structure factor under Random Phase Approximation.

### `S_SALR_MSA(ϕ, K1, z1, K2, z2, k)`
- **Description**: Exact analytical Mean Spherical Approximation structure factor for Two-Yukawa / SALR fluid via 1D scalar Baxter root-finder.
- **Arguments**:
  - `ϕ::Real`: Volume fraction.
  - `K1::Real`: Attraction contact amplitude ($K_1 > 0$).
  - `z1::Real`: Attraction screening parameter.
  - `K2::Real`: Repulsion contact amplitude ($K_2 < 0$).
  - `z2::Real`: Repulsion screening parameter.
  - `k::Real`: Dimensionless wavevector ($q\sigma$).

---

## 3. Fractal Aggregates

### `S_Fractal_Chen(q, R0, Df, xi)`
- **Description**: Sow-Hsin Chen / Teixeira mass-fractal static structure factor.
- **Arguments**:
  - `q::Real`: Scattering wavevector ($1/\text{\AA}$ or $1/\text{nm}$).
  - `R0::Real`: Primary monomer particle radius.
  - `Df::Real`: Fractal dimension ($1 < D_f < 3$).
  - `xi::Real`: Aggregate correlation cutoff length (cluster size).

---

## 4. Form Factors & Scattering Intensity

### `P_sphere(q, R)`
- **Description**: Normalized spherical form factor $P(q) = |F(q, R)|^2$.

### `P_coreshell(q, Rc, t, rho_c, rho_s, rho_solv)`
- **Description**: Form factor for a concentric core-shell sphere.

### `P_multishell(q, radii, slds, rho_solv)`
- **Description**: Form factor for multi-shell concentric spheres.

### `P_ellipsoid(q, a, b)`
- **Description**: Orientationally averaged form factor for a rotational ellipsoid ($a$: rotational axis, $b$: equatorial radius).

### `beta_ellipsoid(q, a, b)`
- **Description**: Orientational decoupling parameter $\beta(q) = |\langle F(q) \rangle|^2 / \langle |F(q)|^2 \rangle \in [0, 1]$.

### `intensity_decoupled(q, S_q, P_q; beta=1.0, I0=1.0, background=0.0)`
- **Description**: Decoupled coherent scattering intensity $I(q) = I_0 P(q) [1 + \beta(q)(S(q) - 1)] + \text{background}$.

---

## 5. Multicomponent Mixtures

### `S_HS_Baxter_mixture(phi_vec, sigma_vec, k)`
- **Description**: $p \times p$ partial structure factor matrix for hard-sphere mixture using Baxter Percus-Yevick solution.

### `S_HS_VW_mixture(phi_vec, sigma_vec, k)`
- **Description**: $p \times p$ partial structure factor matrix for hard-sphere mixture with Verlet-Weis corrections.

### `S_RPA_mixture_Yukawa(phi_vec, sigma_vec, K_mat, Z_mat, k)`
- **Description**: $p \times p$ partial structure factor matrix for multicomponent Yukawa mixture under RPA.

### `S_RPA_mixture_SquareWell(phi_vec, sigma_vec, lambda_mat, eps_mat, k)`
- **Description**: $p \times p$ partial structure factor matrix for multicomponent Square-Well mixture under RPA.

### `S_RPA_mixture_SALR(phi_vec, sigma_vec, K1_mat, Z1_mat, K2_mat, Z2_mat, k)`
- **Description**: $p \times p$ partial structure factor matrix for multicomponent SALR mixture under RPA.

---

## 6. Colloidal Electrostatics & DLVO

### `bjerrum_length(T; epsilon_r = 78.4)`
- **Description**: Calculates Bjerrum length $\ell_B$ (meters) in a solvent of dielectric constant $\varepsilon_r$ at temperature $T$ (Kelvin).

### `ionic_strength(c_salts, valencies)`
- **Description**: Computes molar ionic strength $I = \frac{1}{2} \sum c_i z_i^2$ (mol/L).

### `debye_screening_length(I, T; epsilon_r = 78.4)`
- **Description**: Calculates Debye screening length $\kappa^{-1}$ (meters).

### `dlvo_yukawa_parameters(Z_eff, sigma, I, T; epsilon_r = 78.4)`
- **Description**: Converts physical colloid charge and buffer salt concentration into dimensionless contact energy $K_{\text{DLVO}}$ and screening parameter $z = \kappa\sigma$.

---

## 7. Fourier Inversion Tools

### `sk_to_gr(q_vec, Sq_vec, rho; r_max = 10.0, nr = 1000, damping = 0.0)`
- **Description**: 3D isotropic Fourier inversion from $S(q)$ to radial distribution function $g(r)$.

### `gr_to_sk(r_vec, gr_vec, rho; q_max = 25.0, nq = 500)`
- **Description**: 3D isotropic forward Fourier transform from $g(r)$ to static structure factor $S(q)$.
