# Gnuplot script for Figure 5 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 7.5in,3.5in
set output "fig5.pdf"

set multiplot layout 1,2 title "{/:Bold FIG. 5: Volume Fraction Dependence of S(Q) and I_{cluster}}" font "Helvetica,12"

# Panel (a)
set title "(a) S(Q) for K_1=10, Z_1=10, K_2=-1, Z_2=0.5" font "Helvetica,10"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:3.5]
set grid
set key top right font "Helvetica,8"

plot "Fig5/theory_Fig5a_phi_0.05.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory {/Symbol f}=0.05", \
     "Fig5/Fig5a_phi_0.05.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned {/Symbol f}=0.05", \
     "Fig5/theory_Fig5a_phi_0.2.dat" u 1:2 w l dt 2 lc rgb "#ff7f0e" lw 2 t "Theory {/Symbol f}=0.20", \
     "Fig5/Fig5a_phi_0.2.dat" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned {/Symbol f}=0.20", \
     "Fig5/theory_Fig5a_phi_0.4.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory {/Symbol f}=0.40", \
     "Fig5/Fig5a_phi_0.4.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned {/Symbol f}=0.40", \
     "Fig5/theory_Fig5a_phi_0.55.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory {/Symbol f}=0.55", \
     "Fig5/Fig5a_phi_0.55.dat" u 1:2 w p pt 12 ps 0.6 lc rgb "#d62728" t "Scanned {/Symbol f}=0.55"

# Panel (b)
set title "(b) Cluster peak intensity I_{cluster} vs {/Symbol f}" font "Helvetica,10"
set xlabel "Volume fraction {/Symbol f}" font "Helvetica,10"
set ylabel "I_{cluster}" font "Helvetica,10"
set xrange [0:0.55]
set yrange [0.5:2.5]
set grid
set key top right font "Helvetica,8"

plot "Fig5/theory_Fig5b.dat" u 1:2 w l lc rgb "#9467bd" lw 2 t "Theory (Optimal {/Symbol f} {/Symbol \\approx} 0.20)", \
     "Fig5/Fig5b.dat" u 1:2 w p pt 7 ps 0.8 lc rgb "#d62728" t "Scanned data"

unset multiplot
