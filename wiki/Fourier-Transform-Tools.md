# Fourier Transform Tools

`AnalyticalStructureFactors.jl` provides fast, robust 3D isotropic Fourier transformation routines for converting back and forth between reciprocal-space static structure factors $S(q)$ and real-space pair correlation functions $g(r)$.

---

## 1. Direct Inversion: $S(q) \to g(r)$ (`sk_to_gr`)

In an isotropic 3D fluid with number density $\rho$, the radial distribution function $g(r)$ is related to the structure factor $S(q)$ by the inverse 3D Fourier transform:
$$g(r) = 1 + \frac{1}{2\pi^2 \rho r} \int_0^\infty q \left[S(q) - 1\right] \sin(qr) \, dq$$

To eliminate finite-range truncation oscillations (Gibbs ringing), `sk_to_gr` incorporates optional Gaussian or Lorch damping:
$$S_{\text{damped}}(q) = \left[S(q) - 1\right] \exp\left(-\alpha q^2\right)$$

### Function Signature:
```julia
r_grid, gr_vals = sk_to_gr(q_grid, Sq_vals, rho; r_max = 10.0, nr = 1000, damping = 0.0)
```

### Usage Example:
```julia
using AnalyticalStructureFactors

ϕ = 0.35
σ = 1.0
rho = 6 * ϕ / (π * σ^3)

# Generate fine q grid
q_grid = range(0.001, 50.0, length = 2000)
Sq_vals = [S_HS_PY(ϕ, q * σ) for q in q_grid]

# Invert to g(r)
r_grid, gr_vals = sk_to_gr(q_grid, Sq_vals, rho; r_max = 8.0, nr = 800)

println("Contact value g(r = 1.0+): ", gr_vals[findfirst(>=(1.0), r_grid)])
```

---

## 2. Forward Transform: $g(r) \to S(q)$ (`gr_to_sk`)

The forward 3D Fourier transform converts a radial distribution function $g(r)$ (e.g., obtained from Molecular Dynamics or Monte Carlo simulations) into the experimental static structure factor $S(q)$:
$$S(q) = 1 + \frac{4\pi \rho}{q} \int_0^\infty r \left[g(r) - 1\right] \sin(qr) \, dr$$

At $q = 0$, the compressibility limit is evaluated analytically:
$$S(0) = 1 + 4\pi \rho \int_0^\infty r^2 \left[g(r) - 1\right] dr$$

### Function Signature:
```julia
q_grid, Sq_vals = gr_to_sk(r_grid, gr_vals, rho; q_max = 25.0, nq = 500)
```

### Usage Example:
```julia
# Re-transform g(r) back to S(q)
q_rec, Sq_rec = gr_to_sk(r_grid, gr_vals, rho; q_max = 25.0, nq = 500)
```
