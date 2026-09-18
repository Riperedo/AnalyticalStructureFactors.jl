#!/usr/bin/env bash
# ==============================================================================
# Step-by-step reproduction script for Cummings & Smith (1979) Fig. 1 benchmark
# ==============================================================================
set -e

# Change directory to the script's directory
cd "$(dirname "$0")"

echo "[Step 1/3] Generating analytical MSA roots using Julia..."
julia generate_cummings1979_data.jl

echo "[Step 2/3] Generating comparative PDF and PNG plots with Gnuplot..."
gnuplot plot_cummings1979.gp

echo "[Step 3/3] Compiling LaTeX benchmark report..."
pdflatex -interaction=nonstopmode report_cummings1979.tex > /dev/null
pdflatex -interaction=nonstopmode report_cummings1979.tex > /dev/null

echo "========================================================================"
echo "Reproduction completed successfully!"
echo "Outputs generated:"
echo "  - Data files:     analytical_K_*.dat"
echo "  - Plot (PNG):     cummings1979_fig1_comparison.png"
echo "  - Plot (PDF):     cummings1979_fig1_comparison.pdf"
echo "  - Report (PDF):   report_cummings1979.pdf"
echo "========================================================================"
