# test/data_test/Menon1991/generate_menon1991_data.jl
#
# Reproduces the analytical benchmark data for:
# - Fig. 2: Phase diagram / coexistence curve (binodal) and spinodal boundary τ vs η
# - Fig. 3: Static structure factor S(Q) in the homogeneous fluid phase
# - Fig. 4: Static structure factor S(Q) in the phase coexistence (gas phase) regime
# as presented in Menon, Manohar, and Rao (1991), J. Chem. Phys. 95, 9186.
#
# Uses AnalyticalStructureFactors.jl to perform all analytical evaluations.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__

# ------------------------------------------------------------------------------
# 1. Figure 2: Phase Diagram (Binodal Coexistence & Spinodal Boundary)
# ------------------------------------------------------------------------------
println("Generating Figure 2 data (Phase diagram / Coexistence curve)...")

# A. Analytical Spinodal Boundary (where discriminant D = 0)
# tau_s(eta) = (sqrt(eta*(1 + eta/2)/3) - eta) / (1 - eta)
# Defined for eta in (0, 0.4142) where tau_s >= 0
eta_spinodal_grid = range(0.001, 0.414, length=1000)
fig2_spinodal_file = joinpath(out_dir, "menon_fig2_spinodal.dat")

open(fig2_spinodal_file, "w") do io
    println(io, "# eta  tau_spinodal")
    for eta in eta_spinodal_grid
        tau_s = baxter_shs_critical_tau(eta)
        @printf(io, "%.6f  %.6f\n", eta, tau_s)
    end
end
println("  Saved Fig. 2 spinodal curve -> $fig2_spinodal_file")

# B. Smooth Binodal / Coexistence Curve from Scanned Fig2.dat points
fig2_raw_file = joinpath(out_dir, "Fig2.dat")
pts = Tuple{Float64, Float64}[]
open(fig2_raw_file, "r") do f
    for line in eachline(f)
        s = strip(line)
        (startswith(s, "#") || isempty(s)) && continue
        parts = split(s)
        push!(pts, (parse(Float64, parts[1]), parse(Float64, parts[2])))
    end
end
sort!(pts, by = x -> x[1])

x_pts = [p[1] for p in pts]
y_pts = [p[2] for p in pts]
n = length(x_pts)
h = diff(x_pts)

# Natural cubic spline construction
alpha_vec = zeros(n)
for i in 2:n-1
    alpha_vec[i] = (3.0/h[i])*(y_pts[i+1] - y_pts[i]) - (3.0/h[i-1])*(y_pts[i] - y_pts[i-1])
end
l = zeros(n); mu = zeros(n); z = zeros(n)
l[1] = 1.0; mu[1] = 0.0; z[1] = 0.0
for i in 2:n-1
    l[i] = 2.0*(x_pts[i+1] - x_pts[i-1]) - h[i-1]*mu[i-1]
    mu[i] = h[i]/l[i]
    z[i] = (alpha_vec[i] - h[i-1]*z[i-1])/l[i]
end
l[n] = 1.0; z[n] = 0.0
c_coef = zeros(n); b_coef = zeros(n); d_coef = zeros(n)
for j in (n-1):-1:1
    c_coef[j] = z[j] - mu[j]*c_coef[j+1]
    b_coef[j] = (y_pts[j+1] - y_pts[j])/h[j] - h[j]*(c_coef[j+1] + 2.0*c_coef[j])/3.0
    d_coef[j] = (c_coef[j+1] - c_coef[j])/(3.0*h[j])
end

function eval_binodal_spline(x)
    if x <= x_pts[1]; return y_pts[1]; end
    if x >= x_pts[end]; return y_pts[end]; end
    idx = searchsortedlast(x_pts, x)
    idx = clamp(idx, 1, n-1)
    dx = x - x_pts[idx]
    return y_pts[idx] + b_coef[idx]*dx + c_coef[idx]*dx^2 + d_coef[idx]*dx^3
end

fig2_binodal_file = joinpath(out_dir, "menon_fig2_binodal_smooth.dat")
open(fig2_binodal_file, "w") do io
    println(io, "# eta  tau_binodal_coexistence")
    fine_eta = range(x_pts[1], x_pts[end], length=1000)
    for eta in fine_eta
        tau_val = eval_binodal_spline(eta)
        @printf(io, "%.6f  %.6f\n", eta, tau_val)
    end
end
println("  Saved Fig. 2 continuous binodal curve -> $fig2_binodal_file")

# C. Critical point & thermodynamic state points for Figs. 3 and 4
eta_c = 3.0 / sqrt(2.0) - 2.0
tau_c_exact = (2.0 - sqrt(2.0)) / 6.0

# Fig 3: phi = 0.07, a/sigma = 1.1 => epsilon = 1 - 1/1.1 = 0.090909 => eta = phi / (1-eps)^3 = 0.093077
#        -u0/kT = 0.92 => tau = exp(-0.92) / (12*eps) = 0.365313
eps_fig3 = 1.0 - 1.0 / 1.1
eta_fig3 = 0.07 / ((1.0 - eps_fig3)^3)
tau_fig3 = exp(-0.92) / (12.0 * eps_fig3)

# Fig 4: phi = 0.07, a/sigma = 1.02 => epsilon = 1 - 1/1.02 = 0.019608 => eta_nominal = 0.074351
#        -u0/kT = 3.83 => tau = exp(-3.83) / (12*eps) = 0.092305
# Gas phase equilibrium density: eta_g ≈ 0.0581
eps_fig4 = 1.0 - 1.0 / 1.02
eta_fig4_nominal = 0.07 / ((1.0 - eps_fig4)^3)
tau_fig4 = exp(-3.83) / (12.0 * eps_fig4)
eta_fig4_gas = 0.0581

statepoints_file = joinpath(out_dir, "menon_fig2_statepoints.dat")
open(statepoints_file, "w") do io
    println(io, "# label  eta  tau")
    @printf(io, "Critical_Point    %.6f  %.6f\n", eta_c, tau_c_exact)
    @printf(io, "Fig3_Homogeneous  %.6f  %.6f\n", eta_fig3, tau_fig3)
    @printf(io, "Fig4_Nominal      %.6f  %.6f\n", eta_fig4_nominal, tau_fig4)
    @printf(io, "Fig4_GasPhase     %.6f  %.6f\n", eta_fig4_gas, tau_fig4)
end
println("  Saved Fig. 2 state points -> $statepoints_file")

# ------------------------------------------------------------------------------
# 2. Figure 3: Structure Factor S(Q) in Homogeneous Fluid Phase
# ------------------------------------------------------------------------------
println("Generating Figure 3 data (Homogeneous S(Q))...")
# QR in [0.0, 6.0] => kappa = 2*QR in [0.0, 12.0]
QR_grid = range(0.001, 6.0, length=1000)
fig3_file = joinpath(out_dir, "menon_fig3_analytical.dat")

open(fig3_file, "w") do io
    println(io, "# QR  kappa  S(Q)")
    for QR in QR_grid
        kappa = 2.0 * QR
        s_val = S_Baxter_SHS(eta_fig3, tau_fig3, kappa)
        @printf(io, "%.6f  %.6f  %.6f\n", QR, kappa, s_val)
    end
end
println("  Saved Fig. 3 analytical curve -> $fig3_file")

# ------------------------------------------------------------------------------
# 3. Figure 4: Structure Factor S(Q) in Phase Coexistence (Gas Phase)
# ------------------------------------------------------------------------------
println("Generating Figure 4 data (Phase coexistence gas S(Q))...")
fig4_file = joinpath(out_dir, "menon_fig4_analytical.dat")

open(fig4_file, "w") do io
    println(io, "# QR  kappa  S(Q)")
    for QR in QR_grid
        kappa = 2.0 * QR
        s_val = S_Baxter_SHS(eta_fig4_gas, tau_fig4, kappa)
        @printf(io, "%.6f  %.6f  %.6f\n", QR, kappa, s_val)
    end
end
println("  Saved Fig. 4 analytical curve -> $fig4_file")

println("\nAll Menon et al. (1991) analytical benchmark data successfully generated!")
