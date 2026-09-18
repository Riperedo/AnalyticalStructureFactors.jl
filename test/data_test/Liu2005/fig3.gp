# Gnuplot script for Figure 3 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 7.5in,3.5in
set output "fig3.pdf"

set multiplot layout 1,2 title "{/:Bold FIG. 3: Effect of Repulsion Amplitude |K_2| on S(Q) and I_{cluster}}" font "Helvetica,12"

# Panel (a)
set title "(a) S(Q) for \\phi=0.20, K_1=6.9, Z_1=10, Z_2=0.5" font "Helvetica,10"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:3.5]
set grid
set key top right font "Helvetica,8"

plot "Fig3/theory_Fig3a_K2_1.0.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory K_2=-1.0", \
     "Fig3/Fig3a_K2_1.0.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned K_2=-1.0", \
     "Fig3/theory_Fig3a_K2_0.1.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory K_2=-0.1", \
     "Fig3/Fig3a_K2_0.1.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned K_2=-0.1", \
     "Fig3/theory_Fig3a_K2_0.01.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory K_2=-0.01", \
     "Fig3/Fig3a_K2_0.01.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned K_2=-0.01"

# Panel (b)
set title "(b) Cluster peak intensity I_{cluster} vs |K_2|" font "Helvetica,10"
set xlabel "|K_2| (Repulsion strength)" font "Helvetica,10"
set ylabel "I_{cluster}" font "Helvetica,10"
set xrange [0:2.0]
set yrange [0.5:9.0]
set grid
set key top right font "Helvetica,8"

plot "Fig3/theory_Fig3b.dat" u 1:2 w l lc rgb "#9467bd" lw 2 t "Theory (SALR MSA)", \
     "Fig3/Fig3b.dat" u 1:2 w p pt 7 ps 0.8 lc rgb "#d62728" t "Scanned data"

unset multiplot
