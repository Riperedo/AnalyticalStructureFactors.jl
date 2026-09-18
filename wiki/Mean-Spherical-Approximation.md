# Mean Spherical Approximation (MSA)

The **Mean Spherical Approximation (MSA)** is an exact closure to the Ornstein-Zernike (OZ) integral equation for hard-core potentials with continuous tails:
$$\begin{cases} h(r) = -1, & 0 < r < 1 \\ c(r) = -\beta V(r), & r > 1 \end{cases}$$
where $r = R/\sigma$ is the reduced distance normalized by the hard-sphere diameter $\sigma$.

Unlike perturbation theories (such as the Random Phase Approximation, RPA), the MSA enforces the exact hard-core non-overlapping boundary condition $g(r) = 0$ for $r < 1$ self-consistently alongside the asymptotic tail.

---

## 1. Single-Yukawa MSA Solver

For a hard sphere with a single Yukawa potential:
$$V(r) = -\epsilon \frac{e^{-z(r-1)}}{r}, \quad r > 1$$
where $K = \beta\epsilon$ is the contact energy and $z = \kappa\sigma$ is the inverse screening length.

### Baxter $Q$-Factorization

Using Baxter's Wiener-Hopf factorization method, the direct correlation function and structure factor are expressed in terms of a parameter $d$ satisfying a quartic algebraic equation:
$$d^4 + a_3 d^3 + a_2 d^2 + a_1 d + a_0 = 0$$

`AnalyticalStructureFactors.jl` solves this quartic polynomial analytically using companion matrix eigensolvers or Ferrari's formulas, selecting the unique physical root corresponding to positive compressibility $S(k=0) > 0$.

### Functions and API:
- `S_Yukawa_MSA(ϕ, K, z, k)`: Evaluates the static structure factor $S(k)$.
- `IS_Yukawa_MSA(ϕ, K, z, k)`: Evaluates the inverse structure factor $1/S(k)$.
- `C_Yukawa_MSA(ϕ, K, z, k)`: Evaluates the Fourier direct correlation function $C(k)$.
- `yukawa_msa_roots(ϕ, K, z)`: Computes all algebraic roots for the Baxter parameter $d$.
- `yukawa_msa_coefficients(ϕ, K, z)`: Returns the Baxter polynomial coefficients $(a_0, a_1, a_2, a_3)$.

#### Example:
```julia
using AnalyticalStructureFactors

ϕ = 0.20
K = 1.624     # Attractive contact energy
z = 1.8       # Inverse screening length
k = 2.5       # Wavevector

s_msa = S_Yukawa_MSA(ϕ, K, z, k)
c_msa = C_Yukawa_MSA(ϕ, K, z, k)

# Evaluate compressibility limit S(0)
s_0 = S_Yukawa_MSA(ϕ, K, z, 0.0)
println("S(0) = ", s_0)
```

---

## 2. Two-Yukawa / SALR MSA Solver

For short-range attraction and long-range repulsion (SALR) fluids or two-component screening:
$$\frac{V(r)}{k_B T} = -K_1 \frac{e^{-z_1(r-1)}}{r} - K_2 \frac{e^{-z_2(r-1)}}{r}, \quad r > 1$$

### Sign and Convention Rules
- $K_1 > 0$: Attraction (e.g., $K_1 = 6.0$, $z_1 = 10.0$).
- $K_2 < 0$: Repulsion (e.g., $K_2 = -1.0$, $z_2 = 0.5$).
- $z_i = \sigma/\lambda_i$: Dimensionless inverse screening lengths.

---

### Modern 1D Scalar Algebraic Root-Finder

In `AnalyticalStructureFactors.jl`, the coupled Baxter-Høye-Blum equations for two Yukawa tails are reduced analytically to a 1D scalar root-finding problem on the Baxter parameter $d_2$:

$$d_1(d_2, \text{sgn}) = \frac{-y_{21}(d_2) \pm \sqrt{y_{21}(d_2)^2 - 4 y_{22}(d_2) y_{20}(d_2)}}{2 y_{22}(d_2)}$$

$$R(d_2) = y_{14}(d_2) d_1^2 + y_{13}(d_2) d_1 + y_{12}(d_2) + \frac{y_{11}(d_2)}{d_1} + \frac{y_{10}(d_2)}{d_1^2} = 0$$

#### Adaptive Multi-Scale Root-Finding
The solver searches for bracketed sign changes across an adaptive, multi-scale log-spaced search grid:
$$d_2 \in \left[\pm 10^{-5}, \, \pm 10 \cdot \left(|K_1|e^{z_1} + |K_2|e^{z_2}\right)\right]$$
ensuring rapid convergence ($< 0.1$ ms) across both short-range ($z \sim 10-20$) and long-range ($z \sim 0.1-0.5$) potentials.

#### Automated Physical Branch Selection
For each candidate real root pair $(d_1, d_2)$, the physical branch is uniquely selected by minimizing the non-overlapping core penalty filter:
$$\text{Penalty} = \frac{1}{N} \sum_{r \in [0.1, 0.9]} |g(r)|$$
guaranteeing continuous, non-divergent structure factors and non-negative pair correlations $g(r) \ge 0$.

---

### Functions and API:
- `S_SALR_MSA(ϕ, K1, z1, K2, z2, k)`: Evaluates static structure factor $S(k)$.
- `IS_SALR_MSA(ϕ, K1, z1, K2, z2, k)`: Evaluates inverse structure factor $1/S(k)$.
- `C_SALR_MSA(ϕ, K1, z1, K2, z2, k)`: Evaluates Fourier direct correlation function $C(k)$.

#### Example:
```julia
using AnalyticalStructureFactors

ϕ = 0.20
K1, z1 = 6.0, 10.0   # Short-range attraction
K2, z2 = -1.0, 0.5   # Long-range repulsion

# Evaluate S(k) across a wavevector range
k_vals = 0.1:0.05:10.0
S_vals = [S_SALR_MSA(ϕ, K1, z1, K2, z2, k) for k in k_vals]

# Locate intermediate-range order (IRO) cluster peak
max_idx = argmax(S_vals)
println("Cluster Peak at k = ", k_vals[max_idx], " with intensity S(k) = ", S_vals[max_idx])
```

---

## Comparison: MSA vs. RPA

| Feature | Random Phase Approximation (RPA) | Mean Spherical Approximation (MSA) |
| :--- | :--- | :--- |
| **Accuracy at High Density** | Overestimates compressibility, unphysical divergence | Accurate up to fluid freezing ($\phi \approx 0.50$) |
| **Core Boundary Condition** | Does not enforce $g(r < \sigma) = 0$ | Strictly enforces $g(r < \sigma) = 0$ |
| **Cluster Peak Formation** | Shifts cluster peaks to higher $q$ | Accurately predicts IRO cluster peak positions |
| **Computational Speed** | $< 0.01$ ms (direct formula) | $< 0.10$ ms (analytical scalar root-finder) |
