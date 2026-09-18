# Benchmark Validation

`AnalyticalStructureFactors.jl` contains three independent benchmark suites located under [`test/data_test/`](https://github.com/Riperedo/AnalyticalStructureFactors.jl/tree/main/test/data_test) that systematically reproduce landmark publications from the statistical physics and soft matter literature.

Each benchmark suite contains:
- Raw scanned digitized experimental/theoretical data points.
- Julia reproduction scripts using `AnalyticalStructureFactors.jl`.
- Gnuplot visual validation scripts.
- Automated bash pipeline (`reproduce_figure.sh`).
- Full publication-ready LaTeX validation reports (`.pdf`).

---

## 1. Waisman (1973): Single-Yukawa MSA Pair Correlation $g(r)$

*Reference: E. Waisman, Mol. Phys. 25, 45–48 (1973).*

- **System**: Hard core with single attractive Yukawa tail at $\eta = \phi = \pi/6 \approx 0.5236$, $z = 1.8$, and $K = 1.0827, 1.2180, 1.6240$.
- **Validation**: Verifies that the analytical Fourier inversion `sk_to_gr` applied to `S_Yukawa_MSA` reproduces the exact discrete $g(r)$ table values reported by Waisman.
- **Location**: [`test/data_test/Waisman1973/`](https://github.com/Riperedo/AnalyticalStructureFactors.jl/tree/main/test/data_test/Waisman1973)
- **Reproduction**:
  ```bash
  cd test/data_test/Waisman1973
  bash reproduce_figure.sh
  ```

---

## 2. Cummings & Smith (1979): Analytic Quartic Roots & $S(Q)$

*Reference: P. T. Cummings and E. R. Smith, Mol. Phys. 38, 997–1011 (1979).*

- **System**: Single-Yukawa MSA solver across four coupling states ($K = 0.8120, 1.0827, 1.2180, 1.6240$ at $\phi = 0.5236, z = 1.8$).
- **Validation**: Compares the continuous $S(Q)$ curves computed by `S_Yukawa_MSA` with scanned curves from Figure 1 of Cummings & Smith (1979).
- **Location**: [`test/data_test/CummingsSmith1979/`](https://github.com/Riperedo/AnalyticalStructureFactors.jl/tree/main/test/data_test/CummingsSmith1979)
- **Reproduction**:
  ```bash
  cd test/data_test/CummingsSmith1979
  bash reproduce_figure.sh
  ```

---

## 3. Liu, Chen, & Chen (2005): SALR Cluster Formation & Protein SANS

*Reference: Y. Liu, W.-R. Chen, and S.-H. Chen, J. Chem. Phys. 122, 044507 (2005).*

- **System**: Two-Yukawa / SALR fluid with competing short-range attraction and long-range repulsion under Mean Spherical Approximation.
- **Validation Scope**:
  - **Figures 1–5**: Evolution of the intermediate-range order (IRO) cluster peak in $S(Q)$ and peak intensity $I_{\text{cluster}}$ across variations in $K_1$, $1/Z_1$, $|K_2|$, $1/Z_2$, and volume fraction $\phi$.
  - **Figure 7**: Structure factor evolution $S(Q)$ at $\phi=0.15$ ($K_1 = 3.6, 4.5, 5.5, 6.72$).
  - **Figure 9**: Long-range attraction and short-range repulsion on logarithmic wavevector scale.
  - **Figures 11 & 12**: Decoupled scattering modeling $I(Q) = I_0 P(Q) S(Q)$ and interaction potential decomposition $V(r) = V_{\text{att}}(r) + V_{\text{rep}}(r)$ for aqueous Cytochrome C protein solutions (10.18 wt\% and 20.40 wt\%).
- **Location**: [`test/data_test/Liu2005/`](https://github.com/Riperedo/AnalyticalStructureFactors.jl/tree/main/test/data_test/Liu2005)
- **Reproduction**:
  ```bash
  cd test/data_test/Liu2005
  bash reproduce_figure.sh
  ```

### Quantitative Error Summary

| Benchmark Dataset | State / Curve Parameter | Points ($N$) | RMS Deviation |
| :--- | :--- | :---: | :---: |
| **Liu (2005) Fig. 1(a)** | $K_1 = 0.0, 3.0, 6.0, 10.0$ | 113 | 0.0392 |
| **Liu (2005) Fig. 2(a)** | $Z_1 = 14.0, 8.0, 4.0$ | 85 | 0.0298 |
| **Liu (2005) Fig. 3(a)** | $K_2 = -1.0, -0.1, -0.01$ | 90 | 0.0273 |
| **Liu (2005) Fig. 4(a)** | $Z_2 = 0.1, 2.0, 4.0, 8.0$ | 114 | 0.0226 |
| **Liu (2005) Fig. 5(a)** | $\phi = 0.05, 0.20, 0.40, 0.55$ | 110 | 0.0418 |
| **Liu (2005) Fig. 7** | $K_1 = 3.6, 4.5, 5.5, 6.72$ | 100 | 0.0222 |
| **Liu (2005) Fig. 9** | $K_1 = 0.0, 0.2, 0.4, 0.45$ | 95 | 0.0343 |
| **Liu (2005) Fig. 11** | Cytochrome C 10.18 wt\% SANS | 120 | 0.0245 |
| **Liu (2005) Fig. 12** | Cytochrome C 20.40 wt\% SANS | 120 | 0.0289 |
