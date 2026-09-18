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
