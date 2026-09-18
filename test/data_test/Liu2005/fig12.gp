# Gnuplot script for Figure 12 of Liu et al. (2005)
set terminal pdfcairo enhanced color font "Helvetica,10" size 10.0in,3.5in
set output "fig12.pdf"

set multiplot layout 1,3 title "{/:Bold FIG. 12: Application to SANS of Cytochrome C (20.40 wt%)}" font "Helvetica,12"

# Panel (a) SANS Intensity
set title "(a) SANS Intensity I(Q)" font "Helvetica,10"
set xlabel "Q (\\305^{-1})" font "Helvetica,10"
set ylabel "I(Q) (cm^{-1})" font "Helvetica,10"
set xrange [0:0.25]
set yrange [0:0.6]
set grid
set key top right font "Helvetica,8"

plot "Fig12/theory_Fig12a_Iq.dat" u 1:2 w l lc rgb "#1f77b4" lw 2 t "Model Fit", \
     "Fig12/Fig12a_open-circles.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#d62728" t "SANS Data"

# Panel (b) Form Factor & Structure Factor
set title "(b) P(Q) and S(Q)" font "Helvetica,10"
set xlabel "Q (\\305^{-1})" font "Helvetica,10"
set ylabel "P(Q), S(Q)" font "Helvetica,10"
set xrange [0:0.25]
set yrange [0:1.2]
set grid
set key top right font "Helvetica,8"

plot "Fig12/theory_Fig12b_Pq.dat" u 1:2 w l dt 3 lc rgb "#ff7f0e" lw 2 t "Theory P(Q)", \
     "Fig12/Fig12b_dotted.dat" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned P(Q)", \
     "Fig12/theory_Fig12b_Sq.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory S(Q)", \
     "Fig12/Fig12b_dashed-dotted.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned S(Q)"

# Panel (c) Interprotein Potential
set title "(c) Effective Interprotein Potential" font "Helvetica,10"
set xlabel "r (\\305)" font "Helvetica,10"
set ylabel "V(r) / k_B T" font "Helvetica,10"
set xrange [25:90]
set yrange [-0.5:1.6]
set grid
set key top right font "Helvetica,8"

plot "Fig12/theory_Fig12c_vtot.dat" u 1:2 w l lc rgb "#1f77b4" lw 2 t "Theory V_{tot}(r)", \
     "Fig12/Fig12c_solid.dat" u 1:2 w p pt 6 ps 0.6 lc rgb "#1f77b4" t "Scanned V_{tot}(r)", \
     "Fig12/theory_Fig12c_vrep.dat" u 1:2 w l dt 4 lc rgb "#2ca02c" lw 2 t "Theory V_{rep}(r)", \
     "Fig12/Fig12c_dashed-dotted.dat" u 1:2 w p pt 8 ps 0.6 lc rgb "#2ca02c" t "Scanned V_{rep}(r)", \
     "Fig12/theory_Fig12c_vatt.dat" u 1:2 w l dt 3 lc rgb "#ff7f0e" lw 2 t "Theory V_{att}(r)", \
     "Fig12/Fig12c_doted.dat" u 1:2 w p pt 4 ps 0.6 lc rgb "#ff7f0e" t "Scanned V_{att}(r)"

unset multiplot
