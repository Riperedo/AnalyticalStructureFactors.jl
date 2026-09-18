# Gnuplot script for Figure 9 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 6.0in,4.0in
set output "fig9.pdf"

set title "{/:Bold FIG. 9: S(Q) for Short-Range Repulsion & Long-Range Attraction (\\phi=0.20, Z_1=0.5, K_2=-2, Z_2=2)}" font "Helvetica,11"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set logscale x
set xrange [0.05:10]
set yrange [0:3.5]
set grid
set key top right font "Helvetica,8"

plot "Fig9/theory_Fig9_K1_0.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory K_1=0", \
     "Fig9/Fig9_K1_0.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned K_1=0", \
     "Fig9/theory_Fig9_K1_0.2.dat" u 1:2 w l dt 3 lc rgb "#ff7f0e" lw 2 t "Theory K_1=0.2", \
     "Fig9/Fig9_K1_0.2" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned K_1=0.2", \
     "Fig9/theory_Fig9_K1_0.4.dat" u 1:2 w l dt 2 lc rgb "#2ca02c" lw 2 t "Theory K_1=0.4", \
     "Fig9/Fig8_K1_0.4.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned K_1=0.4", \
     "Fig9/theory_Fig9_K1_0.45.dat" u 1:2 w l dt 4 lc rgb "#d62728" lw 2 t "Theory K_1=0.45", \
     "Fig9/Fig9_K1_0.45.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned K_1=0.45"
