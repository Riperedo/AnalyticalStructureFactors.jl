# Gnuplot script to compare Menon et al. (1991) SHS model with the Sticky Limits
# (Conventional B2-matching, Pure Hard Spheres, and Monte Carlo data)

set encoding utf8

set terminal pdfcairo enhanced color font 'Helvetica,10' size 7.5in, 7.5in
set output 'menon1991_sticky_limit_comparison.pdf'

set multiplot layout 2,1 rowsfirst title "Comparison: Menon SHS Mapping vs Conventional Sticky Limit and Hard Spheres" font "Helvetica-Bold,13"

# ------------------------------------------------------------------------------
# Panel (a): Fig 3 Homogeneous Fluid Phase
# ------------------------------------------------------------------------------
set title "(a) Homogeneous Fluid Phase ({/Symbol f}=0.07, a/{/Symbol s}=1.10, -u_0/k_B T=0.92)" font "Helvetica-Bold,11"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,10"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,10"

set xrange [0.0:6.0]
set yrange [0.4:2.0]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

set style line 1 lc rgb "#1f77b4" lw 2.2 dt 1           # Menon SHS
set style line 2 lc rgb "#ff7f0e" lw 2.0 dt 2           # Conventional B2 Sticky
set style line 3 lc rgb "#7f7f7f" lw 1.5 dt 3           # Pure HS (phi=0.07)
set style line 4 lc rgb "#d62728" pt 6 ps 0.8 lw 1.2    # MC Circles
set style line 5 lc rgb "#2ca02c" pt 8 ps 0.8 lw 1.2    # MC Triangles

set key top right box opaque spacing 1.15 font "Helvetica,8.5"

plot \
    'menon_fig3_sticky_limit_comparison.dat' using 1:3 with lines ls 1 title "Menon SHS ({/Symbol h}=0.0931, {/Symbol t}=0.3653)", \
    'menon_fig3_sticky_limit_comparison.dat' using 1:4 with lines ls 2 title "Conv. Sticky Limit B_2 ({/Symbol h}=0.07, {/Symbol t}=0.5004)", \
    'menon_fig3_sticky_limit_comparison.dat' using 1:5 with lines ls 3 title "Pure Hard Sphere PY ({/Symbol h}=0.07)", \
    'Fig3_circles.dat' using 1:2 with points ls 4 title "MC Simulation (Circles)", \
    'Fig3_triangles.dat' using 1:2 with points ls 5 title "MC Simulation (Triangles)"

# ------------------------------------------------------------------------------
# Panel (b): Fig 4 Two-Phase Coexistence Gas Phase
# ------------------------------------------------------------------------------
set title "(b) Phase Coexistence Gas Phase ({/Symbol f}=0.07, a/{/Symbol s}=1.02, -u_0/k_B T=3.83)" font "Helvetica-Bold,11"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,10"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,10"

set xrange [0.0:6.0]
set yrange [0.4:4.2]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'menon_fig4_sticky_limit_comparison.dat' using 1:3 with lines ls 1 title "Menon SHS Gas ({/Symbol h}_g=0.0581, {/Symbol t}=0.0923)", \
    'menon_fig4_sticky_limit_comparison.dat' using 1:4 with lines ls 2 title "Conv. Sticky Gas ({/Symbol h}_g=0.0581, {/Symbol t}=0.0906)", \
    'menon_fig4_sticky_limit_comparison.dat' using 1:5 with lines ls 3 title "Pure Hard Sphere PY ({/Symbol h}_g=0.0581)", \
    'Fig4_circles.dat' using 1:2 with points ls 4 title "MC Simulation (Circles)", \
    'Fig4_triangles.dat' using 1:2 with points ls 5 title "MC Simulation (Triangles)"

unset multiplot

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,10' size 1000, 1000
set output 'menon1991_sticky_limit_comparison.png'

set multiplot layout 2,1 rowsfirst title "Comparison: Menon SHS Mapping vs Conventional Sticky Limit and Hard Spheres" font "Helvetica-Bold,13"

set title "(a) Homogeneous Fluid Phase ({/Symbol f}=0.07, a/{/Symbol s}=1.10, -u_0/k_B T=0.92)" font "Helvetica-Bold,11"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,10"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.4:2.0]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8
set key top right box opaque spacing 1.15 font "Helvetica,8.5"

plot \
    'menon_fig3_sticky_limit_comparison.dat' using 1:3 with lines ls 1 title "Menon SHS ({/Symbol h}=0.0931, {/Symbol t}=0.3653)", \
    'menon_fig3_sticky_limit_comparison.dat' using 1:4 with lines ls 2 title "Conv. Sticky Limit B_2 ({/Symbol h}=0.07, {/Symbol t}=0.5004)", \
    'menon_fig3_sticky_limit_comparison.dat' using 1:5 with lines ls 3 title "Pure Hard Sphere PY ({/Symbol h}=0.07)", \
    'Fig3_circles.dat' using 1:2 with points ls 4 title "MC Simulation (Circles)", \
    'Fig3_triangles.dat' using 1:2 with points ls 5 title "MC Simulation (Triangles)"

set title "(b) Phase Coexistence Gas Phase ({/Symbol f}=0.07, a/{/Symbol s}=1.02, -u_0/k_B T=3.83)" font "Helvetica-Bold,11"
set xlabel "Dimensionless Wavevector QR" font "Helvetica-Bold,10"
set ylabel "Structure Factor S(Q)" font "Helvetica-Bold,10"
set xrange [0.0:6.0]
set yrange [0.4:4.2]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'menon_fig4_sticky_limit_comparison.dat' using 1:3 with lines ls 1 title "Menon SHS Gas ({/Symbol h}_g=0.0581, {/Symbol t}=0.0923)", \
    'menon_fig4_sticky_limit_comparison.dat' using 1:4 with lines ls 2 title "Conv. Sticky Gas ({/Symbol h}_g=0.0581, {/Symbol t}=0.0906)", \
    'menon_fig4_sticky_limit_comparison.dat' using 1:5 with lines ls 3 title "Pure Hard Sphere PY ({/Symbol h}_g=0.0581)", \
    'Fig4_circles.dat' using 1:2 with points ls 4 title "MC Simulation (Circles)", \
    'Fig4_triangles.dat' using 1:2 with points ls 5 title "MC Simulation (Triangles)"

unset multiplot
