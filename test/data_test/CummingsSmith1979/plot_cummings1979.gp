# Gnuplot script to compare analytical MSA beta roots with Cummings & Smith (1979) Fig. 1 data

set terminal pdfcairo enhanced color font 'Helvetica,12' size 6.0in, 4.5in
set output 'cummings1979_fig1_comparison.pdf'

set title "Cummings & Smith (1979) - Fig. 1: Parameter {/Symbol b} vs Density {/Symbol h} ({/Symbol x} = 2)" font "Helvetica-Bold,13"
set xlabel "Packing Fraction {/Symbol h}" font "Helvetica-Bold,12"
set ylabel "Baxter Parameter {/Symbol b}" font "Helvetica-Bold,12"

set xrange [0.0:0.7]
set yrange [0.0:2.5]
set grid linetype 1 linecolor rgb "#E0E0E0" linewidth 1.0

# Styles
# Colors: Blue (K=0.8120), Green (K=1.0827), Orange (K=1.2180), Red (K=1.6240)
set style line 1 lc rgb "#1f77b4" lw 2 dt 1  # K=0.8120 line
set style line 2 lc rgb "#1f77b4" pt 6 ps 0.8 lw 1.5 # K=0.8120 points

set style line 3 lc rgb "#2ca02c" lw 2 dt 1  # K=1.0827 line
set style line 4 lc rgb "#2ca02c" pt 4 ps 0.8 lw 1.5 # K=1.0827 points

set style line 5 lc rgb "#ff7f0e" lw 2 dt 1  # K=1.2180 line
set style line 6 lc rgb "#ff7f0e" pt 8 ps 0.8 lw 1.5 # K=1.2180 points

set style line 7 lc rgb "#d62728" lw 2 dt 1  # K=1.6240 line
set style line 8 lc rgb "#d62728" pt 12 ps 0.8 lw 1.5 # K=1.6240 points

set key top right box opaque spacing 1.2 font "Helvetica,10"

plot \
    'analytical_K_0.8120.dat' using 1:2 with lines ls 1 title "K = 0.8120 (Analytical)", \
    'analytical_K_0.8120.dat' using 1:3 with lines ls 1 notitle, \
    'Fig1_K_0.812_1.dat' using 1:2 with points ls 2 title "K = 0.8120 (Scanned)", \
    'Fig1_K_0.812_2.dat' using 1:2 with points ls 2 notitle, \
    \
    'analytical_K_1.0827.dat' using 1:2 with lines ls 3 title "K = 1.0827 (Analytical)", \
    'analytical_K_1.0827.dat' using 1:3 with lines ls 3 notitle, \
    'Fig1_K_1.0827_1 copy.dat' using 1:2 with points ls 4 title "K = 1.0827 (Scanned)", \
    'Fig1_K_1.0827_2.dat' using 1:2 with points ls 4 notitle, \
    \
    'analytical_K_1.2180.dat' using 1:2 with lines ls 5 title "K = 1.2180 (Analytical)", \
    'analytical_K_1.2180.dat' using 1:3 with lines ls 5 notitle, \
    'Fig1_K_1.2180_1.dat' using 1:2 with points ls 6 title "K = 1.2180 (Scanned)", \
    'Fig1_K_1.2180_2.dat' using 1:2 with points ls 6 notitle, \
    \
    'analytical_K_1.6240.dat' using 1:2 with lines ls 7 title "K = 1.6240 (Analytical)", \
    'analytical_K_1.6240.dat' using 1:3 with lines ls 7 notitle, \
    'Fig1_K_1.6240_1.dat' using 1:2 with points ls 8 title "K = 1.6240 (Scanned)", \
    'Fig1_K_1.6240_2.dat' using 1:2 with points ls 8 notitle

# Output PNG version
set terminal pngcairo enhanced color font 'Helvetica,12' size 900, 650
set output 'cummings1979_fig1_comparison.png'
replot
