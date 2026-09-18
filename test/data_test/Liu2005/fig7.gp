# Gnuplot script for Figure 7 of Liu et al. (2005) - Structure Factors
set terminal pdfcairo enhanced color font "Helvetica,10" size 5.5in,3.8in
set output "fig7.pdf"

set title "{/:Bold FIG. 7: Static Structure Factor S(Q) at {/Symbol f}=0.15 for Indicated Attraction K_1}" font "Helvetica,11"
set xlabel "Q{/Symbol s}" font "Helvetica,10"
set ylabel "S(Q)" font "Helvetica,10"
set xrange [0:10]
set yrange [0:4.0]
set grid
set key top right font "Helvetica,9"

plot "Fig7/theory_Fig7b_K1_3.6.dat" u 1:2 w l dt 1 lc rgb "#1f77b4" lw 2 t "Theory K_1=3.6", \
     "Fig7/Fig7b_K1_3.6.dat" u 1:2 w p pt 6 ps 0.7 lc rgb "#1f77b4" t "Scanned K_1=3.6", \
     "Fig7/theory_Fig7b_K1_4.5.dat" u 1:2 w l dt 2 lc rgb "#ff7f0e" lw 2 t "Theory K_1=4.5", \
     "Fig7/Fig7b_K1_4.5.dat" u 1:2 w p pt 4 ps 0.7 lc rgb "#ff7f0e" t "Scanned K_1=4.5", \
     "Fig7/theory_Fig7b_K1_5.5.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory K_1=5.5", \
     "Fig7/Fig7b_K1_5.5.dat" u 1:2 w p pt 8 ps 0.7 lc rgb "#2ca02c" t "Scanned K_1=5.5", \
     "Fig7/theory_Fig7b_K1_6.72.dat" u 1:2 w l dt 3 lc rgb "#d62728" lw 2 t "Theory K_1=6.72", \
     "Fig7/Fig7b_K1_6.72.dat" u 1:2 w p pt 12 ps 0.7 lc rgb "#d62728" t "Scanned K_1=6.72"
