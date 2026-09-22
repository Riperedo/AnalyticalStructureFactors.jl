# Monodisperse Models

Monodisperse models describe single-component systems of spherical particles with diameter $\sigma$, number density $\rho$, and volume fraction $\phi = \frac{\pi}{6}\rho\sigma^3$. In dimensionless units, the wavevector is normalized as $k = q\sigma$.

---

## 1. Hard Spheres

The hard-sphere interaction potential is defined as:
$$V_{\text{HS}}(r) = \begin{cases} \infty, & r < 1 \\ 0, & r \ge 1 \end{cases}$$
where $r = R/\sigma$.

### Percus-Yevick Solution (`S_HS_PY`)

The exact analytical solution of the Ornstein-Zernike equation with Percus-Yevick closure gives the direct correlation function $C(k)$:
$$C_{\text{HS}}(k) = -\frac{4\pi}{k^3} \left[ \alpha (\sin k - k\cos k) + \beta \left(\frac{2k\sin k + (2 - k^2)\cos k - 2}{k}\right) + \gamma \left(\frac{(4k^3 - 24k)\sin k - (k^4 - 12k^2 + 24)\cos k + 24}{k^3}\right) \right]$$
with density parameters:
$$\alpha = \frac{(1 + 2\phi)^2}{(1 - \phi)^4}, \quad \beta = \frac{-6\phi(1 + \phi/2)^2}{(1 - \phi)^4}, \quad \gamma = \frac{\phi(1 + 2\phi)^2}{2(1 - \phi)^4}$$

The static structure factor is:
$$S_{\text{HS}}(k) = \frac{1}{1 - \rho C_{\text{HS}}(k)}$$

#### Usage Example:
```julia
using AnalyticalStructureFactors

ϕ = 0.35
k = 3.2
s = S_HS_PY(ϕ, k)
c = C_HS_PY(ϕ, k)
inv_s = IS_HS_PY(ϕ, k)  # 1 / S(k)
```

---

### Verlet-Weis Semi-Empirical Correction (`S_HS_VW`)

The Verlet-Weis correction improves the Percus-Yevick approximation near the principal diffraction peak by rescaling the volume fraction and wavevector:
$$\phi_w = \phi \left(1 - \frac{\phi}{16}\right), \quad k_w = k \left(\frac{\phi_w}{\phi}\right)^{1/3}$$

The structure factor is evaluated as:
$$S_{\text{VW}}(k) = S_{\text{PY}}(k_w, \phi_w) + \Delta S_{\text{VW}}(k, \phi)$$

#### Usage Example:
```julia
s_vw = S_HS_VW(ϕ, k)
```

---

### Weeks-Chandler-Andersen Soft Sphere Blip Function (`S_WCA_blip`)

For continuous soft repulsive potentials (e.g., Lennard-Jones repulsive core), the WCA blip function maps the soft potential to an effective hard-sphere diameter $d_{\text{eff}} = \sigma (1 - B(T))$:
$$B(T) = \int_0^\infty \left(1 - e^{-\beta V_{\text{rep}}(r)}\right) dr - 1$$

#### Usage Example:
```julia
T_star = 1.2
B = blip(T_star)
s_wca = S_WCA_blip(ϕ, k, T_star)
```

---

## 2. Square-Well Fluid (RPA)

The square-well potential combines hard-sphere repulsion with an attractive well of depth $\epsilon$ and range $\lambda\sigma$:
$$V_{\text{SW}}(r) = \begin{cases} \infty, & r < 1 \\ -\epsilon, & 1 \le r \le \lambda \\ 0, & r > \lambda \end{cases}$$

Under the Random Phase Approximation (RPA):
$$C_{\text{SW}}(k) = C_{\text{HS}}(k) - \beta \hat{U}_{\text{SW}}(k)$$
where $\beta = 1/(k_B T)$ and $\hat{U}_{\text{SW}}(k)$ is the 3D Fourier transform of the attractive well:
$$\hat{U}_{\text{SW}}(k) = -\frac{4\pi \epsilon \sigma^3}{k^3} \left[ \sin(\lambda k) - \lambda k \cos(\lambda k) - \sin k + k\cos k \right]$$

### Usage Example:
```julia
ϕ = 0.20
λ = 1.5       # Well width (relative to σ)
ϵ_kT = 1.2    # Attraction depth ε / (k_B * T)
k = 2.0

s_sw = S_SW_RPA(ϕ, λ, ϵ_kT, k)
```

---

## 3. Screened Coulomb / Yukawa Fluid (RPA)

The single-Yukawa potential represents screened electrostatic or dispersion interactions:
$$V_{\text{Yukawa}}(r) = \begin{cases} \infty, & r < 1 \\ -\epsilon \dfrac{e^{-z(r-1)}}{r}, & r \ge 1 \end{cases}$$
where $z = \kappa\sigma$ is the dimensionless inverse screening length and $K = \beta\epsilon$ is the contact energy.

The Fourier transform of the attractive tail is:
$$\hat{U}_{\text{Yukawa}}(k) = \frac{4\pi \epsilon \sigma^3}{k} \left[ \frac{z \sin k + k \cos k}{z^2 + k^2} \right]$$

### Usage Example:
```julia
ϕ = 0.15
K = 2.0       # Contact energy βε
z = 1.8       # Inverse screening length κσ
k = 3.0

s_yukawa_rpa = S_Yukawa_RPA(ϕ, K, z, k)
```

---

## 4. SALR Fluid (RPA)

Short-range attraction and long-range repulsion (SALR) fluids feature competing interactions modeled by two Yukawa tails:
$$V_{\text{SALR}}(r) = \begin{cases} \infty, & r < 1 \\ -K_1 \dfrac{e^{-z_1(r-1)}}{r} - K_2 \dfrac{e^{-z_2(r-1)}}{r}, & r \ge 1 \end{cases}$$
where $K_1 > 0$ denotes short-range attraction ($z_1 \sim 5-10$) and $K_2 < 0$ denotes long-range repulsion ($z_2 \sim 0.5-2$).

Under the RPA:
$$S_{\text{SALR}}(k) = \frac{1}{1 - \rho \left[ C_{\text{HS}}(k) + \beta\hat{U}_{\text{Yukawa}}(k; K_1, z_1) + \beta\hat{U}_{\text{Yukawa}}(k; K_2, z_2) \right]}$$

### Usage Example:
```julia
ϕ = 0.20
K1, z1 = 6.0, 10.0   # Short-range attraction
K2, z2 = -1.0, 0.5   # Long-range repulsion
k = 1.5

s_salr_rpa = S_SALR_RPA(ϕ, K1, z1, K2, z2, k)
```

> **Tip**: For dense liquids and strong coupling regimes, use the exact **Mean Spherical Approximation (MSA)** instead of RPA, as described in [Mean Spherical Approximation](Mean-Spherical-Approximation).

---

## 5. Sticky Hard Spheres (Baxter & Menon Models)

The **Sticky Hard Sphere (SHS)** model describes hard spheres with infinitely narrow, infinitely deep attractive surface adhesion under the Percus-Yevick approximation.

### 5.1 Baxter Sticky Limit (`S_Baxter_SHS`)
Baxter (1968) solved the Percus-Yevick equation using Wiener-Hopf factorization:
$$S(k)^{-1} = A^2(\kappa) + B^2(\kappa)$$
where $\kappa = k a$, and $\lambda$ is the physical root of $A\lambda^2 + B\lambda + C = 0$:
$$A = \frac{\eta}{12}, \quad B = -\left(\tau + \frac{\eta^2}{1-\eta}\right), \quad C = \frac{1 + \eta/2}{(1-\eta)^2}$$

```julia
using AnalyticalStructureFactors

eta = 0.25   # Volume fraction
tau = 0.8    # Baxter stickiness parameter (smaller = more sticky, tau -> infty is pure HS)
k = 3.5      # Dimensionless wavevector

s_baxter = S_Baxter_SHS(eta, tau, k)
inv_s = IS_Baxter_SHS(eta, tau, k)

# Critical stickiness parameter (spinodal boundary)
tau_c = baxter_shs_critical_tau(eta)
```

### 5.2 Menon Physical Square-Well Mapping (`S_Menon_SHS`)
Menon et al. (1991) derived an exact physical mapping relating square-well potential parameters (hard-core diameter $\sigma$, well width $\Delta$, well depth $u_0/(k_B T)$, and physical volume fraction $\phi$) to the Baxter parameters:
- $a = \sigma + \Delta$
- $\epsilon = \Delta / a$
- $\eta = \phi / (1 - \epsilon)^3$
- $\tau = \frac{1}{12\epsilon} \exp\left(\frac{u_0}{k_B T}\right)$

```julia
phi = 0.20      # Physical volume fraction
u0_kT = -1.2    # Attractive well depth u0 / (kB * T)
delta = 0.05    # Well width Δ
sigma = 1.0     # Hard core diameter σ
k = 3.5         # Wavevector

s_menon = S_Menon_SHS(phi, u0_kT, delta, sigma, k)
```

### References
- R. J. Baxter, "Percus–Yevick Equation for Hard Spheres with Surface Adhesion", *J. Chem. Phys.* **49**(6), 2770–2774 (1968). DOI: [10.1063/1.1670482](https://doi.org/10.1063/1.1670482).
- S. V. G. Menon, C. Manohar, and K. S. Rao, "A new interpretation of the sticky hard sphere model", *J. Chem. Phys.* **95**(12), 9186–9190 (1991). DOI: [10.1063/1.461199](https://doi.org/10.1063/1.461199).

