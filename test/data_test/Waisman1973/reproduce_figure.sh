#!/usr/bin/env bash
# ==============================================================================
# reproduce_figure.sh - Waisman (1973) Table 1 & g(r) Benchmark Reproduction
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================================================"
echo "Starting Benchmark Reproduction for Waisman (1973) Table 1"
echo "========================================================================"

# Step 1: Generate analytical and inverted g(r) data
echo "[Step 1/3] Generating g(r) from AnalyticalStructureFactors.jl using Julia..."
julia generate_waisman1973_data.jl

# Step 2: Plot comparative figure with Gnuplot
echo "[Step 2/3] Generating comparative PDF and PNG plots with Gnuplot..."
gnuplot plot_waisman1973.gp

# Step 3: Compile LaTeX report
echo "[Step 3/3] Compiling LaTeX benchmark report..."
pdflatex -interaction=nonstopmode report_waisman1973.tex > /dev/null 2>&1
pdflatex -interaction=nonstopmode report_waisman1973.tex > /dev/null 2>&1

echo "========================================================================"
echo "Reproduction completed successfully!"
echo "Outputs generated:"
echo "  - Continuous data: waisman_continuous.dat"
echo "  - Table data:      waisman_table_evaluated.dat"
echo "  - Plot (PNG):      waisman1973_comparison.png"
echo "  - Plot (PDF):      waisman1973_comparison.pdf"
echo "  - Report (PDF):    report_waisman1973.pdf"
echo "========================================================================"
