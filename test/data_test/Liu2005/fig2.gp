# Gnuplot script for Figure 2 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 7.5in,3.5in
set output "fig2.pdf"

set multiplot layout 1,2 title "{/:Bold FIG. 2: Effect of Attraction Range 1/Z_1 on S(Q) and I_{cluster}}" font "Helvetica,12"

# Panel (a)
set title "(a) S(Q) for \\phi=0.20, K_1=6, K_2=-1, Z_2=0.5" font "Helvetica,10"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:2.5]
set grid
set key top right font "Helvetica,8"

plot "Fig2/theory_Fig2a_Z1_14.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory Z_1=14", \
     "Fig2/Fig2a_Z1_14.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned Z_1=14", \
     "Fig2/theory_Fig2a_Z1_8.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory Z_1=8", \
     "Fig2/Fig2a_Z1_8.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned Z_1=8", \
     "Fig2/theory_Fig2a_Z1_4.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory Z_1=4", \
     "Fig2/Fig2a_Z1_4.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned Z_1=4"

# Panel (b)
set title "(b) Cluster peak intensity I_{cluster} vs 1/Z_1" font "Helvetica,10"
set xlabel "1/Z_1 (Attraction range)" font "Helvetica,10"
set ylabel "I_{cluster}" font "Helvetica,10"
set xrange [0.05:0.26]
set yrange [0.5:4.5]
set grid
set key top left font "Helvetica,8"

plot "Fig2/theory_Fig2b.dat" u 1:2 w l lc rgb "#9467bd" lw 2 t "Theory (SALR MSA)", \
     "Fig2/Fig2b.dat" u 1:2 w p pt 7 ps 0.8 lc rgb "#d62728" t "Scanned data"

unset multiplot
