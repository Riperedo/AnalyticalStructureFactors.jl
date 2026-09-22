# Gnuplot script to reproduce Figures 2, 3, and 4 from Menon et al. (1991)
# J. Chem. Phys. 95, 9186–9190.

set encoding utf8

# ==============================================================================
# Figure 2: Phase Diagram (tau vs eta)
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,11' size 5.8in, 4.4in
set output 'menon1991_fig2_phase_diagram.pdf'

set title "Menon et al. (1991) - Fig. 2: Adhesive Hard Sphere Phase Diagram" font "Helvetica-Bold,12"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,11"
set ylabel "Baxter Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,11"

set xrange [0.0:0.60]
set yrange [0.0:0.40]
set grid linetype 1 linecolor rgb "#E0E0E0" linewidth 0.8

set style line 1 lc rgb "#1f77b4" lw 2.2 dt 1           # Binodal curve
set style line 2 lc rgb "#000000" pt 7 ps 0.7 lw 1.2    # Scanned binodal points
set style line 3 lc rgb "#ff7f0e" lw 1.8 dt 2           # Spinodal boundary
set style line 4 lc rgb "#d62728" pt 4 ps 1.1 lw 1.5    # Fig 3 state point
set style line 5 lc rgb "#2ca02c" pt 6 ps 1.1 lw 1.5    # Fig 4 gas point
set style line 6 lc rgb "#2ca02c" pt 8 ps 1.0 lw 1.2    # Fig 4 nominal point
set style line 7 lc rgb "#9467bd" pt 9 ps 1.4 lw 1.8    # Critical point

set key top right box opaque spacing 1.2 font "Helvetica,9.5"

# Annotations
set label 1 at 0.13, 0.11 "Critical Point ({/Symbol h}_c=0.1213, {/Symbol t}_c=0.0976)" font "Helvetica-Bold,8.5" tc rgb "#9467bd"
set label 2 at 0.11, 0.365 "Fig. 3 (Homogeneous, {/Symbol t}=0.3653)" font "Helvetica,8.5" tc rgb "#d62728"
set label 3 at 0.08, 0.075 "Fig. 4 ({/Symbol t}=0.0923)" font "Helvetica,8.5" tc rgb "#2ca02c"

plot \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 1 title "Binodal / Coexistence Curve (Menon 1991)", \
    'Fig2.dat' using 1:2 with points ls 2 title "Binodal Scanned Data (Fig2.dat)", \
    'menon_fig2_spinodal.dat' using 1:2 with lines ls 3 title "Spinodal Boundary {/Symbol t}_s({/Symbol h}) (Analytical)", \
    '< grep Fig3_Homogeneous menon_fig2_statepoints.dat' using 2:3 with points ls 4 title "State Fig. 3 ({/Symbol h}=0.0931, {/Symbol t}=0.3653)", \
    '< grep Fig4_GasPhase menon_fig2_statepoints.dat' using 2:3 with points ls 5 title "State Fig. 4 (Gas {/Symbol h}=0.0581, {/Symbol t}=0.0923)", \
    '< grep Fig4_Nominal menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "State Fig. 4 (Nominal {/Symbol h}=0.0744)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 7 title "Critical Point"

unset label 1; unset label 2; unset label 3

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,11' size 850, 620
set output 'menon1991_fig2_phase_diagram.png'
replot

# ==============================================================================
# Figure 3: Structure Factor S(Q) in Homogeneous Fluid
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,11' size 5.5in, 4.2in
set output 'menon1991_fig3_structure_factor.pdf'

set title "Menon et al. (1991) - Fig. 3: Structure Factor S(Q) (Homogeneous Fluid Phase)" font "Helvetica-Bold,12"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,11"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,11"

set xrange [0.0:6.0]
set yrange [0.5:2.5]
set grid linetype 1 linecolor rgb "#E0E0E0" linewidth 0.8

set style line 1 lc rgb "#1f77b4" lw 2 dt 1          # Analytical line
set style line 2 lc rgb "#d62728" pt 6 ps 0.8 lw 1.2  # Scanned circles
set style line 3 lc rgb "#2ca02c" pt 8 ps 0.8 lw 1.2  # Scanned triangles

set key top right box opaque spacing 1.2 font "Helvetica,9.5"

plot \
    'menon_fig3_analytical.dat' using 1:3 with lines ls 1 title "Baxter PY Analytical ({/Symbol h}=0.0931, {/Symbol t}=0.3653)", \
    'Fig3_circles.dat' using 1:2 with points ls 2 title "Menon et al. (1991) Fig. 3 (Circles)", \
    'Fig3_triangles.dat' using 1:2 with points ls 3 title "Menon et al. (1991) Fig. 3 (Triangles)"

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,11' size 800, 600
set output 'menon1991_fig3_structure_factor.png'
replot

# ==============================================================================
# Figure 4: Structure Factor S(Q) in Coexistence Region
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,11' size 5.5in, 4.2in
set output 'menon1991_fig4_structure_factor.pdf'

set title "Menon et al. (1991) - Fig. 4: Structure Factor S(Q) (Phase Coexistence Region)" font "Helvetica-Bold,12"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,11"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,11"

set xrange [0.0:6.0]
set yrange [0.5:4.0]
set grid linetype 1 linecolor rgb "#E0E0E0" linewidth 0.8

plot \
    'menon_fig4_analytical.dat' using 1:3 with lines ls 1 title "Baxter PY Gas Phase ({/Symbol h}=0.0581, {/Symbol t}=0.0923)", \
    'Fig4_circles.dat' using 1:2 with points ls 2 title "Menon et al. (1991) Fig. 4 (Circles)", \
    'Fig4_triangles.dat' using 1:2 with points ls 3 title "Menon et al. (1991) Fig. 4 (Triangles)"

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,11' size 800, 600
set output 'menon1991_fig4_structure_factor.png'
replot

# ==============================================================================
# Combined 3-Panel Figure
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,10' size 7.5in, 7.5in
set output 'menon1991_combined_comparison.pdf'

set multiplot layout 2,2 rowsfirst title "Validation and Reproduction of Menon et al. (1991)\nAdhesive Hard Sphere Analytical Model" font "Helvetica-Bold,13"

# Panel 1: Fig 2 Phase Diagram
set title "(a) Phase Diagram ({/Symbol t} vs {/Symbol h})" font "Helvetica-Bold,11"
set xlabel "{/Symbol h}" font "Helvetica-Bold,10"
set ylabel "{/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.60]
set yrange [0.0:0.40]
set key top right box opaque spacing 1.1 font "Helvetica,7.8"
plot \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 1 title "Binodal (Menon 1991)", \
    'Fig2.dat' using 1:2 with points ls 2 title "Binodal Scanned", \
    'menon_fig2_spinodal.dat' using 1:2 with lines ls 3 title "Spinodal {/Symbol t}_s({/Symbol h})", \
    '< grep Fig3_Homogeneous menon_fig2_statepoints.dat' using 2:3 with points ls 4 title "State Fig. 3", \
    '< grep Fig4_GasPhase menon_fig2_statepoints.dat' using 2:3 with points ls 5 title "State Fig. 4 (Gas)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 7 title "Critical Pt."

# Panel 2: Fig 3 Homogeneous S(Q)
set title "(b) Homogeneous S(Q) ({/Symbol h}=0.0931, {/Symbol t}=0.3653)" font "Helvetica-Bold,11"
set xlabel "QR" font "Helvetica-Bold,10"
set ylabel "S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.5:2.5]
set key top right box opaque spacing 1.1 font "Helvetica,8"
plot \
    'menon_fig3_analytical.dat' using 1:3 with lines ls 1 title "Analytical", \
    'Fig3_circles.dat' using 1:2 with points ls 2 title "Scanned (Circles)", \
    'Fig3_triangles.dat' using 1:2 with points ls 3 title "Scanned (Triangles)"

# Panel 3: Fig 4 Coexistence S(Q)
set title "(c) Coexistence Gas S(Q) ({/Symbol h}=0.0581, {/Symbol t}=0.0923)" font "Helvetica-Bold,11"
set xlabel "QR" font "Helvetica-Bold,10"
set ylabel "S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.5:4.0]
set key top right box opaque spacing 1.1 font "Helvetica,8"
plot \
    'menon_fig4_analytical.dat' using 1:3 with lines ls 1 title "Analytical (Gas Phase)", \
    'Fig4_circles.dat' using 1:2 with points ls 2 title "Scanned (Circles)", \
    'Fig4_triangles.dat' using 1:2 with points ls 3 title "Scanned (Triangles)"

unset multiplot

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,10' size 1000, 1000
set output 'menon1991_combined_comparison.png'
set multiplot layout 2,2 rowsfirst title "Validation and Reproduction of Menon et al. (1991)\nAdhesive Hard Sphere Analytical Model" font "Helvetica-Bold,13"

set title "(a) Phase Diagram ({/Symbol t} vs {/Symbol h})" font "Helvetica-Bold,11"
set xlabel "{/Symbol h}" font "Helvetica-Bold,10"
set ylabel "{/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.60]
set yrange [0.0:0.40]
set key top right box opaque spacing 1.1 font "Helvetica,7.8"
plot \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 1 title "Binodal (Menon 1991)", \
    'Fig2.dat' using 1:2 with points ls 2 title "Binodal Scanned", \
    'menon_fig2_spinodal.dat' using 1:2 with lines ls 3 title "Spinodal {/Symbol t}_s({/Symbol h})", \
    '< grep Fig3_Homogeneous menon_fig2_statepoints.dat' using 2:3 with points ls 4 title "State Fig. 3", \
    '< grep Fig4_GasPhase menon_fig2_statepoints.dat' using 2:3 with points ls 5 title "State Fig. 4 (Gas)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 7 title "Critical Pt."

set title "(b) Homogeneous S(Q) ({/Symbol h}=0.0931, {/Symbol t}=0.3653)" font "Helvetica-Bold,11"
set xlabel "QR" font "Helvetica-Bold,10"
set ylabel "S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.5:2.5]
set key top right box opaque spacing 1.1 font "Helvetica,8"
plot \
    'menon_fig3_analytical.dat' using 1:3 with lines ls 1 title "Analytical", \
    'Fig3_circles.dat' using 1:2 with points ls 2 title "Scanned (Circles)", \
    'Fig3_triangles.dat' using 1:2 with points ls 3 title "Scanned (Triangles)"

set title "(c) Coexistence Gas S(Q) ({/Symbol h}=0.0581, {/Symbol t}=0.0923)" font "Helvetica-Bold,11"
set xlabel "QR" font "Helvetica-Bold,10"
set ylabel "S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.5:4.0]
set key top right box opaque spacing 1.1 font "Helvetica,8"
plot \
    'menon_fig4_analytical.dat' using 1:3 with lines ls 1 title "Analytical (Gas Phase)", \
    'Fig4_circles.dat' using 1:2 with points ls 2 title "Scanned (Circles)", \
    'Fig4_triangles.dat' using 1:2 with points ls 3 title "Scanned (Triangles)"

unset multiplot
