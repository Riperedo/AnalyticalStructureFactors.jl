Source: Cluster formation in two-Yukawa fluids by Yun Liu, Wei-Ren Chen, and Sow-Hsin Chen

**Resumen Ejecutivo: Cluster formation in two-Yukawa fluids**

Este estudio presenta un marco teórico riguroso y analíticamente eficiente para describir el comportamiento estructural y de fase de fluidos coloidales y proteicos regidos por un potencial de dos Yukawas (2Y), enfocado en la formación de agregados (*clusters*) y la aparición del *cluster peak* en el factor de estructura $S(Q)$.

---

**1. Aporte Teórico y Metodológico**

* **Solución de la Ecuación de Ornstein-Zernike:** Se resuelve la ecuación OZ bajo la Aproximación Esférica Media (MSA) para un potencial $V(r)$ con componentes atractiva y repulsiva.


* **Reducción Algorítmica:** El problema se reduce de un sistema complejo de ecuaciones acopladas no lineales a la resolución numérica de un único polinomio de orden 27 en función de un parámetro $d_2$. Esto permite un cálculo ultra rápido de $S(Q)$ e $h(r)$ en comparación con métodos numéricos estándar.



**2. Formación y Estabilidad de Clusters**

* **Efecto de la Competición SALR:** En sistemas con atracción de corto alcance y repulsión de largo alcance (SALR), la repulsión inhibe la separación de fases macroscópica (gas-líquido) y promueve la autoorganización microfásica en clusters de tamaño finito.


* **Criterios y Diagrama de Fases:** Se establece la *línea de transición baja* (inicio de agregación débil por atracción) y la *línea de transición alta* (formación de clusters compactos bien definidos). A densidades moderadas, el máximo principal en $S(Q)$ se desplaza hacia valores bajos del vector de onda ($Q \to 0$).



**3. Aplicación Experimental (Citocromo C)**

* **Ajuste de Datos SANS:** Se contrasta el modelo con datos de dispersión de neutrones a bajo ángulo (SANS) en soluciones acuosas de la proteína Citocromo C.


* **Hallazgo Clave:** Para reproducir experimentalmente el pico en $Q \to 0$ observados a bajas concentraciones iónicas, se demuestra que el potencial efectivo requiere una componente atractiva de alcance intermedio/largo combinada con la repulsión coulómbica apantallada (GOCM).



---

**Conclusión e Impacto**
El trabajo valida la utilidad de la solución MSA analítica para el potencial de dos Yukawas como una herramienta precisa e instantánea para interpretar datos de dispersión (SANS/SAXS) y predecir fronteras de agregación coloidal sin recurrir a simulación de Monte Carlo de alto costo computacional.


---

El potencial de interacción $V(r)$ utilizado en el trabajo para el sistema de dos Yukawas (2Y) con núcleo duro se define de la siguiente manera:

$$\frac{V(r)}{k_B T} = \begin{cases} \infty, & 0 < r < 1 \\ -K_1 \dfrac{e^{-Z_1(r-1)}}{r} - K_2 \dfrac{e^{-Z_2(r-1)}}{r}, & r > 1 \end{cases}$$

---

### Descripción de los parámetros:

* **$r$ (Distancia reducida):** Distancia entre los centros de dos partículas, adimensionalizada respecto al diámetro del núcleo duro ($\sigma$). Por ello, $r = 1$ representa la distancia de contacto entre dos partículas rígidas.


* **$k_B T$ (Energía térmica):** Producto de la constante de Boltzmann ($k_B$) y la temperatura absoluta ($T$).


* **$K_1$ y $K_2$ (Amplitudes de las interacciones):** Parámetros adimensionales que representan la intensidad o profundidad del potencial a la distancia de contacto ($r = 1$).


* Por convención en este artículo, $K > 0$ representa un potencial de **atracción** y $K < 0$ un potencial de **repulsión**.


* En el caso de interacción atracción a corto alcance y repulsión a largo alcance (SALR), se fija $K_1 > 0$ (atractivo) y $K_2 < 0$ (repulsivo).




* **$Z_1$ y $Z_2$ (Parámetros de alcance / Decaimiento):** Parámetros inversos al alcance espacial de cada interacción ($Z = \sigma / \lambda$, donde $\lambda$ es la longitud de decaimiento del potencial).


* Un valor grande de $Z$ ($Z \sim 10$) indica un potencial de **corto alcance**.


* Un valor pequeño de $Z$ ($Z \sim 0.5 - 2$) indica un potencial de **largo alcance**.

---
##**Transcripcion de los captions de las figuras:**
---

**FIG. 1.** This figure illustrates the existence of a different cluster peak in the structure factor. In panel (a), $\phi=0.20$, $Z_1=10$, $K_2=-1$, $Z_2=0.5$. For $K_1=0$ (solid line) and $K_1=3$ (dashed line), there is no cluster peak. For $K_1=6$ (dash-dotted line) and $K_1=10$ (dotted line), the attraction is strong enough so that the cluster peak appears. Panel (b) shows the change of cluster peak intensity $I_{\text{cluster}}$ as a function of $K_1$.

Fig1b log y
---

**FIG. 2.** In panel (a), $\phi=0.20$, $K_1=6$, $K_2=-1$, $Z_2=0.5$. The solid line, the dash-dotted line, and the dotted line correspond to $Z_1=14, 8$, and $4$, respectively. Panel (b) shows the change of cluster peak intensity $I_{\text{cluster}}$ as a function of attraction range $1/Z_1$.

---

**FIG. 3.** In panel (a), $\phi=0.20$, $K_1=6.9$, $Z_1=10$, $Z_2=0.5$. The solid line, the dash-dotted line, and the dotted line correspond to $K_2=-1, -0.1$, and $-0.01$, respectively. Panel (b) shows the change of cluster peak intensity $I_{\text{cluster}}$ as a function of the amplitude of repulsive potential $\vert{}K_2\vert{}$.

---

**FIG. 4.** Panel (a) demonstrates the effect of repulsion range $1/Z_2$ on the structure factor $S(Q)$ at a volume fraction $\phi=0.20$, $K_1=6.9$, $Z_1=10$, $K_2=-1$. The solid line, the dashed line, the dash-dotted line, and the dotted line correspond to $Z_2=0.1, 2, 4$, and $8$, respectively. Open circles in panel (b) shows the change of cluster peak intensity $I_{\text{cluster}}$ as a function of repulsion range $1/Z_2$. The symbol $\times$ indicates the intensity of $S(Q=0)$ when there is no finite-$Q$ cluster peak.

Fig4b log x
---

**FIG. 5.** Panel (a) shows structure factors at volume fractions, $\phi=0.05, 0.2, 0.4$, and $0.55$. $K_1=10$, $Z_1=10$, $K_2=-1$, $Z_2=0.5$. Panel (b) shows the change of cluster peak intensity $I_{\text{cluster}}$ up to the volume fraction $\phi=0.5$. For a given set of $K_i, Z_i$, there exists an optimal volume fraction, at which the cluster peak intensity $I_{\text{cluster}}$ is largest. For the case presented in this figure, the optimal volume fraction $\phi$ is about 20%.

---

**FIG. 6.** The spinodal lines for a short-range attraction and a long-range repulsion are determined as a function of $1/K_1$ vs $\phi$ under different control parameters. The solid line in all three panels has the same set of parameters, $Z_1=10$, $K_2=-0.3$, and $Z_2=2$. Panel (a) shows the effect on the spinodal line by increasing the repulsion strength $\vert{}K_2\vert{}$. The dotted line, the dashed line, and the dash-dotted line correspond to $K_2=-1, -2$, and $-5$, respectively. Panel (b) shows the effect on the spinodal line by increasing $Z_2$ (decreasing the range of the repulsive potential $1/Z_2$). The dotted line, the dashed line, and the dash-dotted line correspond to $Z_2=3, 4$, and $6$. Panel (c) shows the effect of the attraction range $1/Z_1$ on the spinodal lines. The dotted line and the dashed line correspond to $Z_1=8$ and $6$, respectively.

---

**FIG. 7.** In panel (a), the cluster region in the phase diagram is shown as a function of $1/K_1$ (effective temperature) and the volume fraction $\phi$. $Z_1=10$, $K_2=-0.3$, and $Z_2=2$. The dotted line is determined as the minimum attraction needed to generate a well defined cluster peak in the structure factor, while the dash-dotted line is determined as the maximum attraction needed to have a well defined cluster peak in the structure factor. The region between the dotted line and the dash-dotted is named as cluster region. The solid line is the spinodal line. The structure factors are shown in panel (b) with control parameters indicated by the open circles in the phase diagram of panel (a).

---

**FIG. 8.** The effect of a long-range repulsion on a cluster region is demonstrated. The dotted line in both panels correspond to $Z_1=10$, $K_2=-0.3$, $Z_2=2$. Two lines (high cluster transition line and low cluster transition line) with the same symbol sandwich a cluster region. Panel (a) shows the shift of the cluster region by changing the repulsion range $1/Z_2$. The solid line and the dash-dotted line correspond to $Z_2=3$ and $Z_2=0.5$ respectively. When $Z_2=0.5$ there is no low cluster transition line. Therefore, the area below the high cluster transition line is the cluster region. Panel (b) shows the shift of the cluster region due to the change of the amplitude of the long-range repulsion $\vert{}K_2\vert{}$. The solid line corresponds to $K_2=-2$.

---

**FIG. 9.** The structure factor generated by a short-range repulsion and a long-range attraction is studied. The dependence on the attraction depth $K_1$ is presented in this figure. $Z_1, K_2$, and $Z_2$ are set as $0.5, -2$, and $2$, respectively. The solid line, the dotted line, the dashed line, and the dash-dotted line correspond to $K_1=0, 0.2, 0.4$, and $0.45$, respectively.

Fig9 log x
---

**FIG. 10.** The spinodal lines for a long-range attraction plus a short-range repulsion are determined in the plane of $1/K_1$ vs $\phi$. The solid line in all three panels corresponds to the same set of parameters, $Z_1=0.5$, $K_2=-0.3$, and $Z_2=2$. Panel (a) shows the effect on the spinodal line by increasing the repulsion strength $\vert{}K_2\vert{}$. The dotted line, the dashed line, and the dash-dotted line correspond to $K_2=-1, -2$, and $-5$, respectively. Panel (b) shows the effect on the spinodal line by increasing $Z_2$ (decreasing the range of the repulsive potential $1/Z_2$). The dotted line, the dashed line, and the dash-dotted line correspond to $Z_2=3, 4$, and $6$. Panel (c) shows the effect on the spinodal line by decreasing $Z_1$ (increasing the range of attractive potential $1/Z_1$). The dotted line and the dashed line correspond to $Z_1=1$ and $1.5$, respectively.

---

**FIG. 11.** This figure and the next one give an example of an application of the two-Yukawa model to analyze the small angle neutron scattering (SANS) data of protein solutions. Panel (a) shows the SANS intensity distribution (open circles) $I(Q)$ together with the model fittings (solid line). In panel (b), the dotted line indicates the normalized intraparticle structure factor $\bar{P}(Q)$ (form factor). The dashed-dotted line represents the interparticle structure factor. Panel (c) shows the interprotein potential obtained by fitting the SANS data. The dotted line and the dash-dotted line represent the long-range attraction potential and the short-range repulsion potential. The solid line is the overall potential by summing the attractive potential and the repulsive potential.

---

**FIG. 12.** The lines and symbols are the same as the preceding figure. The SANS intensity distribution $I(Q)$ was obtained at 20.40 wt% cytochrome C proteins in $D_2O$ at room temperature.

---

**TABLE I.** Results for the fitting parameters of the SANS data of cytochrome C protein solutions at two different concentrations.

| Sample concentration (wt%) | $K_1$ | $Z_1$ | $z_p$ | $a$ ($\text{\AA}$) |
| --- | --- | --- | --- | --- |
| **10.18** | 0.33 | 0.51 | 4.2 | 14.8 |
| **20.40** | 0.18 | 0.39 | 3.6 | 14.9 |