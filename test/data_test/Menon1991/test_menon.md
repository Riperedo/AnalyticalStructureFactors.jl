# Benchmark: Adhesive Hard Spheres (Menon et al., 1991)

## Source
- S. V. G. Menon, C. Manohar, and K. S. Rao, *"A new interpretation of the sticky hard sphere model"*, *The Journal of Chemical Physics* **95**(12), 9186–9190 (1991). DOI: [10.1063/1.461199](https://doi.org/10.1063/1.461199).
- R. J. Baxter, *"Percus–Yevick Equation for Hard Spheres with Surface Adhesion"*, *The Journal of Chemical Physics* **49**(6), 2770–2774 (1968). DOI: [10.1063/1.1670482](https://doi.org/10.1063/1.1670482).

## Overview
Menon et al. (1991) provided a physical mapping between physical square-well fluids of finite width $\Delta$ and core diameter $\sigma$, and Baxter's Adhesive/Sticky Hard Sphere (SHS) model. 

### Parameter Mappings
1. **Interaction Range & Perturbation Parameter**:
   $$a = \sigma + \Delta, \quad \epsilon = \frac{\Delta}{a} = 1 - \frac{\sigma}{a}$$
2. **Effective Volume Fraction & Stickiness Parameter**:
   $$\eta = \frac{\phi}{(1 - \epsilon)^3}, \quad \tau = \frac{1}{12\epsilon} \exp\left(\frac{u_0}{k_B T}\right)$$
3. **Dimensionless Wavevector**:
   $$\kappa = k a = 2 (Q R), \quad \text{where } R = \frac{a}{2}$$

## Reproduced Figures
- **Figure 2**: Phase diagram ($\tau$ vs. $\eta$) showing the analytical spinodal line $\tau_s(\eta) = \frac{\sqrt{\eta(1+\eta/2)/3} - \eta}{1-\eta}$, the critical point $(\eta_c \approx 0.1213, \tau_c \approx 0.0976)$, the digitized binodal points (`Fig2.dat`), and the thermodynamic state points for Figs. 3 and 4.
- **Figure 3**: Structure factor $S(Q)$ in the homogeneous fluid phase ($\phi = 0.07, a/\sigma = 1.1, -u_0/k_B T = 0.92 \implies \eta = 0.0931, \tau = 0.3653$) compared against digitized experimental/simulation data (`Fig3_circles.dat`, `Fig3_triangles.dat`).
- **Figure 4**: Structure factor $S(Q)$ in the two-phase coexistence region ($\phi = 0.07, a/\sigma = 1.02, -u_0/k_B T = 3.83 \implies \tau = 0.0923 < \tau_c$), evaluated at the equilibrium gas-phase density $\eta_g \approx 0.0581$ compared against scanned data (`Fig4_circles.dat`, `Fig4_triangles.dat`).

## Reproduction Command
To execute the automated reproduction pipeline:
```bash
./reproduce_figure.sh
```
