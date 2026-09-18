# test/data_test/Waisman1973/plot_waisman1973.gp
#
# Generates high-quality PDF and PNG plots comparing radial distribution function g(r)
# from AnalyticalStructureFactors.jl against Waisman (1973) Table 1 data.

# --- PDF Terminal Output ---
set terminal pdfcairo enhanced color font "Helvetica,12" size 7in,5in
set output "waisman1973_comparison.pdf"

set title "Hard Sphere Fluid at High Density ({/Symbol h} = 0.49): Radial Distribution Function g(r)\nValidation against Waisman (1973) Mol. Phys. 25, 45-48" font "Helvetica-Bold,13"
set xlabel "Reduced Distance x = r / {/Symbol s}" font "Helvetica-Bold,12"
set ylabel "Radial Distribution Function g(r)" font "Helvetica-Bold,12"

set grid lc rgb "#E0E0E0" dt 3 lw 1
set key top right box opaque spacing 1.3 font "Helvetica,11"

set xrange [0.98:2.0]
set yrange [0.3:6.2]

# Colors and styles
set style line 1 lc rgb "#0072BD" lw 2.2 dt 2 pt 6 ps 1.2 # PY (Blue Dashed)
set style line 2 lc rgb "#D9531E" pt 7 ps 1.4              # Monte Carlo (Red/Orange Circles)
set style line 3 lc rgb "#009E73" lw 2.5 pt 9 ps 1.2        # Waisman MSA (Green Solid)
set style line 4 lc rgb "#7E2F8E" pt 4 ps 1.2              # Table This Work (Purple Squares)

plot "waisman_continuous.dat" using 1:3 with lines ls 1 title "Percus-Yevick (PY) [AnalyticalStructureFactors.jl]", \
     "waisman_continuous.dat" using 1:2 with lines ls 3 title "Waisman MSA [AnalyticalStructureFactors.jl]", \
     "tabla.dat" using 1:2 with points ls 2 title "Monte Carlo Computer Exp. [Barker & Henderson 1971]", \
     "tabla.dat" using 1:3 with points pt 6 ps 1.2 lc rgb "#0072BD" title "Percus-Yevick (Table 1)", \
     "tabla.dat" using 1:4 with points pt 4 ps 1.2 lc rgb "#009E73" title "Waisman (Table 1 'This work')"

# --- PNG Terminal Output ---
set terminal pngcairo enhanced font "Helvetica,12" size 1000,750
set output "waisman1973_comparison.png"
replot
