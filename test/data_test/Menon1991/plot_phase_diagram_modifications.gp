# Gnuplot script to plot Phase Diagram modifications in the Sticky Limit
# Sticky Hard Spheres (Menon et al. 1991 / Baxter 1968)

set encoding utf8

set terminal pdfcairo enhanced color font 'Helvetica,10' size 7.5in, 9.5in
set output 'menon1991_phase_diagram_modifications.pdf'

set multiplot layout 3,1 rowsfirst title "Modifications to the Phase Diagram in the Sticky Limit ({/Symbol e} {/Symbol \256} 0)" font "Helvetica-Bold,13"

# ==============================================================================
# Panel (a): (phi, tau) Plane as epsilon -> 0
# ==============================================================================
set title "(a) Coexistence Boundary Evolution in ({/Symbol f}, {/Symbol t}) Plane" font "Helvetica-Bold,11"
set xlabel "Physical Hard Core Volume Fraction {/Symbol f} = {/Symbol h}(1 - {/Symbol e})^3" font "Helvetica-Bold,10"
set ylabel "Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,10"

set xrange [0.0:0.60]
set yrange [0.0:0.12]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

set style line 1 lc rgb "#1f77b4" lw 2.2 dt 1  # Sticky limit (eps = 0)
set style line 2 lc rgb "#2ca02c" lw 2.0 dt 2  # eps = 0.02
set style line 3 lc rgb "#ff7f0e" lw 2.0 dt 4  # eps = 0.05
set style line 4 lc rgb "#d62728" lw 2.0 dt 5  # eps = 0.10

set key top right box opaque spacing 1.15 font "Helvetica,8.5"

plot \
    'phase_diagram_tau_vs_phi.dat' using 6:2 with lines ls 1 title "Sticky Limit {/Symbol e} = 0.00 ({/Symbol f}_c = 0.1213)", \
    'phase_diagram_tau_vs_phi.dat' using 5:2 with lines ls 2 title "Thin Well {/Symbol e} = 0.02 ({/Symbol f}_c = 0.1143)", \
    'phase_diagram_tau_vs_phi.dat' using 4:2 with lines ls 3 title "Intermediate {/Symbol e} = 0.05 ({/Symbol f}_c = 0.1040)", \
    'phase_diagram_tau_vs_phi.dat' using 3:2 with lines ls 4 title "Moderate Well {/Symbol e} = 0.10 ({/Symbol f}_c = 0.0884)"

# ==============================================================================
# Panel (b): Physical Temperature (phi, T*) Plane
# ==============================================================================
set title "(b) Physical Phase Coexistence in ({/Symbol f}, T^* = k_B T / |u_0|) Plane" font "Helvetica-Bold,11"
set xlabel "Physical Hard Core Volume Fraction {/Symbol f}" font "Helvetica-Bold,10"
set ylabel "Reduced Temperature T^* = k_B T / |u_0|" font "Helvetica-Bold,10"

set xrange [0.0:0.50]
set yrange [0.15:0.55]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'phase_diagram_T_vs_phi.dat' using 3:4 with lines lc rgb "#d62728" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.10, T^*_c=0.466)", \
    'phase_diagram_T_vs_phi.dat' using 3:5 with lines lc rgb "#d62728" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.10, T^*_c=0.484)", \
    'phase_diagram_T_vs_phi.dat' using 6:7 with lines lc rgb "#ff7f0e" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.05, T^*_c=0.352)", \
    'phase_diagram_T_vs_phi.dat' using 6:8 with lines lc rgb "#ff7f0e" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.05, T^*_c=0.358)", \
    'phase_diagram_T_vs_phi.dat' using 9:10 with lines lc rgb "#2ca02c" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.02, T^*_c=0.265)", \
    'phase_diagram_T_vs_phi.dat' using 9:11 with lines lc rgb "#2ca02c" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.02, T^*_c=0.266)"

# ==============================================================================
# Panel (c): Critical Point Trajectory vs epsilon
# ==============================================================================
set title "(c) Critical Point Trajectory: Critical Temperature T^*_c and Packing Fraction {/Symbol f}_c vs Well Width {/Symbol e}" font "Helvetica-Bold,11"
set xlabel "Dimensionless Well Width {/Symbol e} = {/Symbol D}/a = 1 - {/Symbol s}/a" font "Helvetica-Bold,10"
set ylabel "Critical Parameters" font "Helvetica-Bold,10"

set xrange [0.0:0.15]
set yrange [0.0:0.60]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'phase_diagram_critical_loci.dat' using 1:3 with lines lc rgb "#1f77b4" lw 2.2 title "Menon Critical Temp. T^*_c({/Symbol e})", \
    'phase_diagram_critical_loci.dat' using 1:4 with lines lc rgb "#ff7f0e" lw 1.8 dt 2 title "Conv. B_2 Critical Temp. T^*_c({/Symbol e})", \
    'phase_diagram_critical_loci.dat' using 1:2 with lines lc rgb "#2ca02c" lw 2.2 title "Critical Volume Fraction {/Symbol f}_c({/Symbol e})"

unset multiplot

# Output PNG
set terminal pngcairo enhanced color font 'Helvetica,10' size 1000, 1200
set output 'menon1991_phase_diagram_modifications.png'

set multiplot layout 3,1 rowsfirst title "Modifications to the Phase Diagram in the Sticky Limit ({/Symbol e} {/Symbol \256} 0)" font "Helvetica-Bold,13"

set title "(a) Coexistence Boundary Evolution in ({/Symbol f}, {/Symbol t}) Plane" font "Helvetica-Bold,11"
set xlabel "Physical Hard Core Volume Fraction {/Symbol f} = {/Symbol h}(1 - {/Symbol e})^3" font "Helvetica-Bold,10"
set ylabel "Stickiness Parameter {/Symbol t}" font "Helvetica-Bold,10"
set xrange [0.0:0.60]
set yrange [0.0:0.12]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8
set key top right box opaque spacing 1.15 font "Helvetica,8.5"

plot \
    'phase_diagram_tau_vs_phi.dat' using 6:2 with lines ls 1 title "Sticky Limit {/Symbol e} = 0.00 ({/Symbol f}_c = 0.1213)", \
    'phase_diagram_tau_vs_phi.dat' using 5:2 with lines ls 2 title "Thin Well {/Symbol e} = 0.02 ({/Symbol f}_c = 0.1143)", \
    'phase_diagram_tau_vs_phi.dat' using 4:2 with lines ls 3 title "Intermediate {/Symbol e} = 0.05 ({/Symbol f}_c = 0.1040)", \
    'phase_diagram_tau_vs_phi.dat' using 3:2 with lines ls 4 title "Moderate Well {/Symbol e} = 0.10 ({/Symbol f}_c = 0.0884)"

set title "(b) Physical Phase Coexistence in ({/Symbol f}, T^* = k_B T / |u_0|) Plane" font "Helvetica-Bold,11"
set xlabel "Physical Hard Core Volume Fraction {/Symbol f}" font "Helvetica-Bold,10"
set ylabel "Reduced Temperature T^* = k_B T / |u_0|" font "Helvetica-Bold,10"
set xrange [0.0:0.50]
set yrange [0.15:0.55]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'phase_diagram_T_vs_phi.dat' using 3:4 with lines lc rgb "#d62728" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.10, T^*_c=0.466)", \
    'phase_diagram_T_vs_phi.dat' using 3:5 with lines lc rgb "#d62728" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.10, T^*_c=0.484)", \
    'phase_diagram_T_vs_phi.dat' using 6:7 with lines lc rgb "#ff7f0e" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.05, T^*_c=0.352)", \
    'phase_diagram_T_vs_phi.dat' using 6:8 with lines lc rgb "#ff7f0e" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.05, T^*_c=0.358)", \
    'phase_diagram_T_vs_phi.dat' using 9:10 with lines lc rgb "#2ca02c" lw 2.2 dt 1 title "Menon Mapping ({/Symbol e}=0.02, T^*_c=0.265)", \
    'phase_diagram_T_vs_phi.dat' using 9:11 with lines lc rgb "#2ca02c" lw 1.5 dt 2 title "Conv. B_2 Matching ({/Symbol e}=0.02, T^*_c=0.266)"

set title "(c) Critical Point Trajectory: Critical Temperature T^*_c and Packing Fraction {/Symbol f}_c vs Well Width {/Symbol e}" font "Helvetica-Bold,11"
set xlabel "Dimensionless Well Width {/Symbol e} = {/Symbol D}/a = 1 - {/Symbol s}/a" font "Helvetica-Bold,10"
set ylabel "Critical Parameters" font "Helvetica-Bold,10"
set xrange [0.0:0.15]
set yrange [0.0:0.60]
set grid linetype 1 linecolor rgb "#EAEAEA" linewidth 0.8

plot \
    'phase_diagram_critical_loci.dat' using 1:3 with lines lc rgb "#1f77b4" lw 2.2 title "Menon Critical Temp. T^*_c({/Symbol e})", \
    'phase_diagram_critical_loci.dat' using 1:4 with lines lc rgb "#ff7f0e" lw 1.8 dt 2 title "Conv. B_2 Critical Temp. T^*_c({/Symbol e})", \
    'phase_diagram_critical_loci.dat' using 1:2 with lines lc rgb "#2ca02c" lw 2.2 title "Critical Volume Fraction {/Symbol f}_c({/Symbol e})"

unset multiplot
