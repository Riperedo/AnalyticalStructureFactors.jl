# Gnuplot script for Figure 4 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 7.5in,3.5in
set output "fig4.pdf"

set multiplot layout 1,2 title "{/:Bold FIG. 4: Effect of Repulsion Range 1/Z_2 on S(Q) and I_{cluster}}" font "Helvetica,12"

# Panel (a)
set title "(a) S(Q) for \\phi=0.20, K_1=6.9, Z_1=10, K_2=-1" font "Helvetica,10"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:2.5]
set grid
set key top right font "Helvetica,8"

plot "Fig4/theory_Fig4a_Z2_0.1.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory Z_2=0.1", \
     "Fig4/Fig4a_Z2_0.1.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned Z_2=0.1", \
     "Fig4/theory_Fig4a_Z2_2.dat" u 1:2 w l dt 2 lc rgb "#ff7f0e" lw 2 t "Theory Z_2=2.0", \
     "Fig4/Fig4a_Z2_2.dat" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned Z_2=2.0", \
     "Fig4/theory_Fig4a_Z2_4.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory Z_2=4.0", \
     "Fig4/Fig4a_Z2_4.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned Z_2=4.0", \
     "Fig4/theory_Fig4a_Z2_8.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory Z_2=8.0", \
     "Fig4/Fig4a_Z2_8.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned Z_2=8.0"

# Panel (b)
set title "(b) Cluster peak intensity I_{cluster} vs 1/Z_2" font "Helvetica,10"
set xlabel "1/Z_2 (Repulsion range)" font "Helvetica,10"
set ylabel "I_{cluster}" font "Helvetica,10"
set logscale x
set xrange [0.1:10.0]
set yrange [0.5:2.5]
set grid
set key top right font "Helvetica,8"

plot "Fig4/theory_Fig4b_circles.dat" u 1:2 w l lc rgb "#9467bd" lw 2 t "Theory (Cluster peak)", \
     "Fig4/Fig4b_circles.dat" u 1:2 w p pt 6 ps 0.8 lc rgb "#1f77b4" t "Scanned I_{cluster}", \
     "Fig4/theory_Fig4b_stars.dat" u 1:2 w l dt 2 lc rgb "#d62728" lw 2 t "Theory S(Q=0)", \
     "Fig4/Fig4b_stars.dat" u 1:2 w p pt 2 ps 0.8 lc rgb "#d62728" t "Scanned S(Q=0)"

unset logscale x
unset multiplot
