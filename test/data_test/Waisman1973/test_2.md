Source: The radial distribution function for a fluid of hard spheres at high densities Mean Spherical integral equation approach by Eduardo Waisman

resumen: En este trabajo, Waisman resuelve analíticamente la ecuación integral de la aproximación esférica media (MSA) para un fluido de esferas rígidas con una cola de potencial de Yukawa fuera del núcleo ($r>d$). El autor aprovecha esta solución exacta con dos propósitos clave: por un lado, desarrollar un enfoque semiempírico para corregir las deficiencias de la aproximación de Percus-Yevick (PY) en el fluido de esferas rígidas a altas densidades (ajustando los parámetros del potencial de Yukawa a la ecuación de estado de Carnahan-Starling y obteniendo excelente acuerdo en la función de distribución radial $g(r)$ frente a experimentos de simulación Monte Carlo); y por otro lado, reportar la solución analítica rigurosa del modelo esférico medio para un fluido cuyas moléculas interactúan mediante un potencial completo de esferas rígidas más una interacción de Yukawa atractiva o repulsiva de la forma $\omega(r) = -A r^{-1} \exp[-z(r-d)]$ para $r>d$.


caption table: $\xi = 0.49$

Los datos de la tabal se encuentran en `tabla.dat`.

Las variables representadas en la tabla corresponden a la estructura del fluido de esferas rígidas a alta densidad ($\eta = \xi = 0.49$):

* **$x$:** Es la distancia interpartícula reducida, dada por $x = r / d$, donde $r$ es la distancia entre los centros de las partículas y $d$ es el diámetro de las esferas rígidas. Por lo tanto, $x = 1$ corresponde al punto de contacto entre dos esferas.


* **$\xi$ (o $\eta$):** Es la fracción de volumen o densidad reducida del sistema, definida como $\xi = \pi \rho d^3 / 6$, donde $\rho$ es la densidad numérica de partículas. En el encabezado de esta tabla se fija en $\xi = 0.49$, que representa un estado de muy alta densidad para el fluido.


* **Valores tabulados ($g(x)$):** La tabla tabula los valores de la **función de distribución radial**, $g(x)$ (o $g(r)$), a la distancia $x$ especificada. Esta función mide la probabilidad relativa de encontrar el centro de una partícula a una distancia $x$ de otra partícula de referencia, comparada con una distribución completamente aleatoria (gas ideal).



Las tres columnas finales comparan los valores de $g(x)$ obtenidos por tres métodos distintos: los resultados de simulación por computadora (Monte Carlo), la solución teórica de Percus-Yevick (P.Y.) y el método propuesto por Waisman en este artículo (This work).

