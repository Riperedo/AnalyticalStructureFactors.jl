# Gnuplot script for Bisection Spinodal & Unphysical State Analysis
# Sticky Hard Spheres (Menon et al. 1991 / Baxter 1968)

set encoding utf8

# ==============================================================================
# 1. Bisection Phase Diagram with Sampled States
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,11' size 6.2in, 4.8in
set output 'menon1991_bisection_phase_diagram.pdf'

set title "Successive Bisection Analysis: Spinodal and Unphysical Negative A(k) States" font "Helvetica-Bold,12"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,11"
set ylabel "Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,11"

set xrange [0.0:0.42]
set yrange [0.0:0.40]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

set style line 1 lc rgb "#2ca02c" pt 7 ps 0.45       # Physical samples (green)
set style line 2 lc rgb "#d62728" pt 4 ps 0.65 lw 1.1 # Negative A samples (red)
set style line 3 lc rgb "#7f7f7f" pt 2 ps 0.50       # Complex lambda samples (grey)
set style line 4 lc rgb "#1f77b4" lw 2.2 dt 1        # Analytical Spinodal (D=0)
set style line 5 lc rgb "#ff7f0e" lw 2.0 dt 2        # Analytical S(0) Divergence (A=0)
set style line 6 lc rgb "#9467bd" pt 9 ps 1.5 lw 2.0 # Critical Point
set style line 7 lc rgb "#333333" lw 1.5 dt 3        # Binodal Coexistence Curve

set key top right box opaque spacing 1.15 font "Helvetica,8.8"

# Annotations
set label 1 at 0.13, 0.11 "Critical Point ({/Symbol h}_c=0.1213, {/Symbol t}_c=0.0976)" font "Helvetica-Bold,8.5" tc rgb "#9467bd"
set label 2 at 0.22, 0.05 "Forbidden Complex {/Symbol l} (D < 0)" font "Helvetica,8" tc rgb "#555555"
set label 3 at 0.25, 0.14 "Negative A(k) / S(k) < 0 (A < 0)" font "Helvetica,8" tc rgb "#d62728"
set label 4 at 0.08, 0.25 "Physical Fluid (A > 0, D > 0)" font "Helvetica,8" tc rgb "#2ca02c"

plot \
    'bisection_samples_complex_lambda.dat' using 1:2 with points ls 3 title "Bisection Probes: Complex {/Symbol l} (D < 0)", \
    'bisection_samples_negative_A.dat' using 1:2 with points ls 2 title "Bisection Probes: Unphysical A(k) < 0", \
    'bisection_samples_physical.dat' using 1:2 with points ls 1 title "Bisection Probes: Physical Fluid (A > 0)", \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 7 title "Binodal Coexistence (Menon 1991)", \
    'analytical_spinodal_curve.dat' using 1:2 with lines ls 4 title "Spinodal Boundary {/Symbol t}_s({/Symbol h}) (D = 0)", \
    'analytical_compressibility_div_curve.dat' using 1:2 with lines ls 5 title "Compressibility Divergence (A(0) = 0)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "Critical Point"

unset label 1; unset label 2; unset label 3; unset label 4

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,11' size 900, 700
set output 'menon1991_bisection_phase_diagram.png'
replot

# ==============================================================================
# 2. Combined Overview Plot
# ==============================================================================
set terminal pdfcairo enhanced color font 'Helvetica,10' size 7.5in, 7.5in
set output 'menon1991_bisection_combined.pdf'

set multiplot layout 2,1 rowsfirst title "Bisection Spinodal and Unphysical State Analysis (Menon / Baxter SHS)" font "Helvetica-Bold,13"

# Top: Phase Diagram
set title "(a) Thermodynamic Phase Diagram with Bisection State Probes" font "Helvetica-Bold,11"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,10"
set ylabel "Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.42]
set yrange [0.0:0.35]
set key top right box opaque spacing 1.1 font "Helvetica,7.8"

plot \
    'bisection_samples_complex_lambda.dat' using 1:2 with points ls 3 title "Complex {/Symbol l} (D < 0)", \
    'bisection_samples_negative_A.dat' using 1:2 with points ls 2 title "Unphysical A(k) < 0", \
    'bisection_samples_physical.dat' using 1:2 with points ls 1 title "Physical Fluid (A > 0)", \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 7 title "Binodal (Menon)", \
    'analytical_spinodal_curve.dat' using 1:2 with lines ls 4 title "Spinodal (D = 0)", \
    'analytical_compressibility_div_curve.dat' using 1:2 with lines ls 5 title "A(0) = 0 Locus", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "Critical Pt."

# Bottom: Bisection Boundaries vs Analytical Predictions
set title "(b) Converged Bisection Boundaries vs Exact Analytical Limits" font "Helvetica-Bold,11"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,10"
set ylabel "Transition Parameter {/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.42]
set yrange [0.0:0.12]
set key top right box opaque spacing 1.1 font "Helvetica,8"

plot \
    'analytical_spinodal_curve.dat' using 1:2 with lines lc rgb "#1f77b4" lw 2.2 title "Exact Spinodal {/Symbol t}_s({/Symbol h}) (D=0)", \
    'analytical_compressibility_div_curve.dat' using 1:2 with lines lc rgb "#ff7f0e" lw 2.0 dt 2 title "Exact A(0)=0 Locus", \
    'bisection_found_boundaries.dat' using 1:2 with points pt 7 ps 0.8 lc rgb "#1f77b4" title "Bisection Converged (D=0)", \
    'bisection_found_boundaries.dat' using 1:3 with points pt 6 ps 0.8 lc rgb "#ff7f0e" title "Bisection Converged (A=0)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "Critical Point"

unset multiplot

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,10' size 1000, 1000
set output 'menon1991_bisection_combined.png'
set multiplot layout 2,1 rowsfirst title "Bisection Spinodal and Unphysical State Analysis (Menon / Baxter SHS)" font "Helvetica-Bold,13"

set title "(a) Thermodynamic Phase Diagram with Bisection State Probes" font "Helvetica-Bold,11"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,10"
set ylabel "Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.42]
set yrange [0.0:0.35]
set key top right box opaque spacing 1.1 font "Helvetica,7.8"

plot \
    'bisection_samples_complex_lambda.dat' using 1:2 with points ls 3 title "Complex {/Symbol l} (D < 0)", \
    'bisection_samples_negative_A.dat' using 1:2 with points ls 2 title "Unphysical A(k) < 0", \
    'bisection_samples_physical.dat' using 1:2 with points ls 1 title "Physical Fluid (A > 0)", \
    'menon_fig2_binodal_smooth.dat' using 1:2 with lines ls 7 title "Binodal (Menon)", \
    'analytical_spinodal_curve.dat' using 1:2 with lines ls 4 title "Spinodal (D = 0)", \
    'analytical_compressibility_div_curve.dat' using 1:2 with lines ls 5 title "A(0) = 0 Locus", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "Critical Pt."

set title "(b) Converged Bisection Boundaries vs Exact Analytical Limits" font "Helvetica-Bold,11"
set xlabel "Effective Volume Fraction {/Symbol h}" font "Helvetica-Bold,10"
set ylabel "Transition Parameter {/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.42]
set yrange [0.0:0.12]
set key top right box opaque spacing 1.1 font "Helvetica,8"

plot \
    'analytical_spinodal_curve.dat' using 1:2 with lines lc rgb "#1f77b4" lw 2.2 title "Exact Spinodal {/Symbol t}_s({/Symbol h}) (D=0)", \
    'analytical_compressibility_div_curve.dat' using 1:2 with lines lc rgb "#ff7f0e" lw 2.0 dt 2 title "Exact A(0)=0 Locus", \
    'bisection_found_boundaries.dat' using 1:2 with points pt 7 ps 0.8 lc rgb "#1f77b4" title "Bisection Converged (D=0)", \
    'bisection_found_boundaries.dat' using 1:3 with points pt 6 ps 0.8 lc rgb "#ff7f0e" title "Bisection Converged (A=0)", \
    '< grep Critical_Point menon_fig2_statepoints.dat' using 2:3 with points ls 6 title "Critical Point"

unset multiplot
