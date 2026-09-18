#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

echo "========================================================================"
echo "Reproducing Liu, Chen, and Chen (2005) Structure Factor Benchmark Report"
echo "========================================================================"

echo "Step 1: Generating theoretical structure factors with AnalyticalStructureFactors.jl..."
julia --project="$ROOT_DIR" "$SCRIPT_DIR/generate_liu2005_data.jl"

echo "Step 2: Plotting Structure Factor Figures using Gnuplot..."
cd "$SCRIPT_DIR"
for gp in fig1.gp fig2.gp fig3.gp fig4.gp fig5.gp fig7.gp fig9.gp fig11.gp fig12.gp; do
    echo "  Generating ${gp%.gp}.pdf from $gp..."
    gnuplot "$gp"
done

echo "Step 3: Compiling LaTeX report..."
pdflatex -interaction=nonstopmode report_liu2005.tex
pdflatex -interaction=nonstopmode report_liu2005.tex

echo "========================================================================"
echo "Reproduction successfully completed!"
echo "Report generated at: $SCRIPT_DIR/report_liu2005.pdf"
echo "========================================================================"
