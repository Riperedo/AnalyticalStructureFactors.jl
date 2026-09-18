# Fractal Structures

Disordered colloidal networks, colloidal gels, and cluster aggregates often exhibit self-similar fractal microstructures over intermediate spatial scales. `AnalyticalStructureFactors.jl` implements the widely used **Sow-Hsin Chen and José Teixeira (1986, 1988)** mass-fractal static structure factor.

---

## 1. Chen-Teixeira Mass-Fractal Model (`S_Fractal_Chen`)

For an aggregate of spherical primary particles (monomers) of radius $R_0$, mass-fractal dimension $D_f$, and macroscopic correlation cutoff length $\xi$ (characterizing the aggregate cluster size), the inter-particle correlation function decays exponentially:
$$g(r) - 1 \propto r^{D_f - 3} \exp(-r/\xi)$$

The Fourier transform yields the analytical static structure factor:
$$S_{\text{fractal}}(q) = 1 + \frac{D_f \, \Gamma(D_f - 1)}{(q R_0)^{D_f} \left[1 + (q \xi)^{-2}\right]^{(D_f - 1)/2}} \sin\left[(D_f - 1) \arctan(q \xi)\right]$$
where $\Gamma(x)$ is the Gamma function.

---

## 2. Physical Regimes and Asymptotic Limits

```
  S(q) ^
       |  Guinier regime: S(0) ~ (ξ / R0)^Df
       |  \
       |   \   Fractal Power-Law: S(q) ~ q^(-Df)
       |    \
       |     \_________________ Asymptotic limit: S(q) -> 1
       +----------------------------------------------------> q
            1/ξ             1/R0
```

1. **Guinier Regime ($q \ll 1/\xi$)**:
   At very low wavevectors, the scattering probes the entire aggregate:
   $$S(0) = 1 + \Gamma(D_f + 1) \left(\frac{\xi}{R_0}\right)^{D_f} \approx N_{\text{agg}}$$
   where $N_{\text{agg}}$ is the average number of monomers per cluster.

2. **Fractal Power-Law Scaling ($1/\xi \ll q \ll 1/R_0$)**:
   The intensity follows pure fractal self-similarity:
   $$S(q) \propto q^{-D_f}$$
   - Diffusion-Limited Cluster Aggregation (DLCA): $D_f \approx 1.75 - 1.85$ (open, branched fractal structures).
   - Reaction-Limited Cluster Aggregation (RLCA): $D_f \approx 2.05 - 2.15$ (denser, compact aggregates).
   - Percolation / Gel Networks: $D_f \approx 2.5$.

3. **Monomer Scale ($q \gg 1/R_0$)**:
   The intra-monomer structure dominates: $S(q) \to 1$.

---

## 3. Usage Example

```julia
using AnalyticalStructureFactors

q = 0.02          # Scattering vector (1/Å)
R0 = 15.0         # Monomer radius (Å)
Df = 2.1          # Fractal dimension (RLCA regime)
xi = 350.0        # Aggregate cluster size (Å)

# Evaluate fractal structure factor
s_frac = S_Fractal_Chen(q, R0, Df, xi)
println("S_fractal(q = 0.02) = ", s_frac)

# Full aggregate scattering profile I(q) = I0 * P_sphere(q, R0) * S_fractal(q)
P_q = P_sphere(q, R0)
I_q = intensity_decoupled(q, s_frac, P_q; I0 = 1.0)
```
