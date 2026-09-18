# Gnuplot script for Figure 1 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 7.5in,3.5in
set output "fig1.pdf"

set multiplot layout 1,2 title "{/:Bold FIG. 1: Cluster Peak in Structure Factor S(Q) and Peak Intensity I_{cluster}}" font "Helvetica,12"

# Panel (a)
set title "(a) S(Q) for \\phi=0.20, Z_1=10, K_2=-1, Z_2=0.5" font "Helvetica,10"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:3.5]
set grid
set key top right font "Helvetica,8"

plot "Fig1/theory_Fig1a_K1_0.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory K_1=0", \
     "Fig1/Fig1a_K1_0.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned K_1=0", \
     "Fig1/theory_Fig1a_K1_3.dat" u 1:2 w l dt 2 lc rgb "#ff7f0e" lw 2 t "Theory K_1=3", \
     "Fig1/Fig1a_K1_3.dat" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned K_1=3", \
     "Fig1/theory_Fig1a_K1_6.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory K_1=6", \
     "Fig1/Fig1a_K1_6.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned K_1=6", \
     "Fig1/theory_Fig1a_K1_10.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory K_1=10", \
     "Fig1/Fig1a_K1_10.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned K_1=10"

# Panel (b)
set title "(b) Cluster peak intensity I_{cluster} vs K_1" font "Helvetica,10"
set xlabel "K_1 (Attraction strength)" font "Helvetica,10"
set ylabel "I_{cluster}" font "Helvetica,10"
set logscale y
set xrange [5:21]
set yrange [0.5:25]
set grid
set key top left font "Helvetica,8"

plot "Fig1/theory_Fig1b.dat" u 1:2 w l lc rgb "#9467bd" lw 2 t "Theory (SALR MSA)", \
     "Fig1/Fig1b.dat" u 1:2 w p pt 7 ps 0.8 lc rgb "#d62728" t "Scanned data"

unset logscale y
unset multiplot
