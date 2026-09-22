#!/usr/bin/env bash
# test/data_test/Menon1991/reproduce_figure.sh
#
# Automated script to reproduce all analytical benchmark figures, sticky limit comparisons,
# phase diagram modifications, successive bisection analysis, and compiled LaTeX report.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

echo "========================================================================"
echo "Reproducing Menon et al. (1991) Benchmarks, Sticky Limits & Phase Diagrams"
echo "========================================================================"

echo "[1/7] Running Julia script to generate analytical curves..."
cd "${ROOT_DIR}"
julia --project=. "${SCRIPT_DIR}/generate_menon1991_data.jl"

echo "[2/7] Running Julia script for sticky limit comparisons..."
julia --project=. "${SCRIPT_DIR}/generate_sticky_limit_comparison.jl"

echo "[3/7] Running Julia script for phase diagram modifications..."
julia --project=. "${SCRIPT_DIR}/generate_phase_diagram_modifications.jl"

echo "[4/7] Running Julia script for successive bisection spinodal analysis..."
julia --project=. "${SCRIPT_DIR}/bisection_spinodal_analysis.jl"

echo "[5/7] Generating benchmark, sticky limit, phase diagram, and bisection plots with Gnuplot..."
cd "${SCRIPT_DIR}"
gnuplot plot_menon1991.gp
gnuplot plot_sticky_limit_comparison.gp
gnuplot plot_phase_diagram_modifications.gp
gnuplot plot_bisection_analysis.gp

echo "[6/7] Compiling LaTeX benchmark report (pass 1)..."
pdflatex -interaction=nonstopmode report_menon1991.tex > /dev/null 2>&1 || true

echo "[7/7] Compiling LaTeX benchmark report (pass 2)..."
pdflatex -interaction=nonstopmode report_menon1991.tex > /dev/null 2>&1

echo "========================================================================"
echo "Menon et al. (1991) benchmark reproduction completed successfully!"
echo "Generated files in ${SCRIPT_DIR}:"
echo "  - menon1991_fig2_phase_diagram.pdf / .png"
echo "  - menon1991_fig3_structure_factor.pdf / .png"
echo "  - menon1991_fig4_structure_factor.pdf / .png"
echo "  - menon1991_combined_comparison.pdf / .png"
echo "  - menon1991_sticky_limit_comparison.pdf / .png"
echo "  - menon1991_phase_diagram_modifications.pdf / .png"
echo "  - menon1991_bisection_phase_diagram.pdf / .png"
echo "  - menon1991_bisection_combined.pdf / .png"
echo "  - report_menon1991.pdf"
echo "========================================================================"
