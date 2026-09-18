# Form Factors and Scattering

The total coherent scattering intensity measured in Small-Angle X-ray Scattering (SAXS) and Small-Angle Neutron Scattering (SANS) experiments depends on both the intra-particle structure (form factor $P(q)$) and inter-particle spatial correlations (structure factor $S(q)$).

---

## 1. Homogeneous Sphere (`P_sphere`)

For a homogeneous sphere of radius $R$:
$$F(q, R) = \frac{3\left(\sin(qR) - qR\cos(qR)\right)}{(qR)^3}$$
$$P(q, R) = |F(q, R)|^2$$

### Usage:
```julia
using AnalyticalStructureFactors

q = 0.05    # Wavevector (1/Å)
R = 15.0    # Sphere radius (Å)

P = P_sphere(q, R)
```

---

## 2. Core-Shell and Multi-Shell Spheres

### Core-Shell Sphere (`P_coreshell`)
For a particle with core radius $R_c$, total radius $R_s = R_c + t$ (thickness $t$), core scattering length density (SLD) $\rho_c$, shell SLD $\rho_s$, and solvent SLD $\rho_{\text{solv}}$:
$$V_c = \frac{4\pi}{3} R_c^3, \quad V_s = \frac{4\pi}{3} R_s^3$$
$$A(q) = (\rho_c - \rho_s) V_c F(q, R_c) + (\rho_s - \rho_{\text{solv}}) V_s F(q, R_s)$$
$$P_{\text{coreshell}}(q) = \frac{|A(q)|^2}{\left[(\rho_c - \rho_s)V_c + (\rho_s - \rho_{\text{solv}})V_s\right]^2}$$

```julia
Rc = 20.0       # Core radius
t = 5.0         # Shell thickness
rho_c = 1.0e-6  # Core SLD
rho_s = 3.0e-6  # Shell SLD
rho_solv = 0.0  # Solvent SLD

P_cs = P_coreshell(q, Rc, t, rho_c, rho_s, rho_solv)
```

### Multi-Shell Sphere (`P_multishell`)
Generalization to $N$ concentric spherical shells with radii $\mathbf{R} = (R_1, R_2, \dots, R_N)$ and SLDs $\boldsymbol{\rho} = (\rho_1, \rho_2, \dots, \rho_N)$:
```julia
radii = [15.0, 25.0, 30.0]
slds = [1.5e-6, 4.0e-6, 2.0e-6]
rho_solv = 0.0

P_multi = P_multishell(q, radii, slds, rho_solv)
```

---

## 3. Ellipsoids and Orientational Averaging

For a triaxial or rotational ellipsoid with semi-axes $a$ (rotational axis) and $b = c$ (equatorial radius), the scattering amplitude at angle $\theta$ ($\mu = \cos\theta$) is:
$$r(\mu) = \sqrt{a^2 \mu^2 + b^2 (1 - \mu^2)}$$
$$F(q, \mu) = \frac{3\left(\sin(q\,r(\mu)) - q\,r(\mu)\cos(q\,r(\mu))\right)}{(q\,r(\mu))^3}$$

### Orientational Averaging
- **Form Factor**:
  $$P_{\text{ellipsoid}}(q) = \langle |F(q, \mu)|^2 \rangle = \int_0^1 |F(q, \mu)|^2 \, d\mu$$
- **Decoupling Parameter $\beta(q)$**:
  $$\beta(q) = \frac{|\langle F(q, \mu) \rangle|^2}{\langle |F(q, \mu)|^2 \rangle} = \frac{\left[\int_0^1 F(q, \mu) \, d\mu\right]^2}{\int_0^1 |F(q, \mu)|^2 \, d\mu}$$
- **Effective Sphere Radius**:
  $$R_{\text{eff}} = \left(a b^2\right)^{1/3}$$

### Usage:
```julia
a = 25.0   # Semi-major axis (Å)
b = 12.0   # Semi-minor axis (Å)

# Form factor averaged over random orientations
P_ell = P_ellipsoid(q, a, b)

# Orientational decoupling parameter beta(q) in [0, 1]
beta = beta_ellipsoid(q, a, b)

# Effective hard sphere radius
R_eff = effective_radius_ellipsoid(a, b)
```

---

## 4. Decoupling Approximation for Scattering Intensity

For non-spherical particles, polydisperse systems, or interacting colloids, Kotlarchyk and Chen's decoupling approximation expresses the total coherent scattering intensity as:
$$I(q) = I_0 P(q) \left[ 1 + \beta(q) \left(S_{\text{eff}}(q) - 1\right) \right] + \text{background}$$

Where:
- For isotropic monodisperse spheres: $\beta(q) \equiv 1$, reducing to $I(q) = I_0 P(q) S(q)$.
- For non-spherical particles (ellipsoids, cylinders): $\beta(q) < 1$, dampening inter-particle interference at high $q$.

### API Function: `intensity_decoupled`
```julia
# Compute decoupled scattering profile
I_q = intensity_decoupled(q, S_val, P_val; beta = beta_val, I0 = 1.0, background = 0.0)
```
