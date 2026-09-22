# test/data_test/Menon1991/generate_phase_diagram_modifications.jl
#
# Investigates how the Sticky Hard Sphere phase diagram is modified:
# 1. In the (phi, tau) plane as well width epsilon -> 0 (the ideal sticky sphere limit).
# 2. In the physical temperature-density plane (phi, T* = kT / |u0|).
# 3. Comparing Menon's physical square-well mapping with conventional B2-matching.
# 4. Critical point trajectory (phi_c, T*_c) as a function of epsilon.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__

# ------------------------------------------------------------------------------
# 1. Load Universal Binodal in (eta, tau) from Scanned Fig2.dat
# ------------------------------------------------------------------------------
pts = Tuple{Float64, Float64}[]
open(joinpath(out_dir, "Fig2.dat")) do f
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

# Natural cubic spline interpolation for universal binodal
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

function eval_binodal(x)
    if x <= x_pts[1]; return y_pts[1]; end
    if x >= x_pts[end]; return y_pts[end]; end
    idx = searchsortedlast(x_pts, x)
    idx = clamp(idx, 1, n-1)
    dx = x - x_pts[idx]
    return y_pts[idx] + b_coef[idx]*dx + c_coef[idx]*dx^2 + d_coef[idx]*dx^3
end

eta_fine = range(x_pts[1], x_pts[end], length=1000)

# ------------------------------------------------------------------------------
# 2. Phase Diagram in (phi, tau) Plane for various epsilon -> 0
# ------------------------------------------------------------------------------
eps_list = [0.10, 0.05, 0.02, 0.00] # 0.00 is ideal sticky sphere limit

file_tau_phi = joinpath(out_dir, "phase_diagram_tau_vs_phi.dat")
open(file_tau_phi, "w") do io
    println(io, "# eta  tau  phi_eps0.10  phi_eps0.05  phi_eps0.02  phi_sticky_limit")
    for eta in eta_fine
        tau = eval_binodal(eta)
        phi_10 = eta * ((1.0 - 0.10)^3)
        phi_05 = eta * ((1.0 - 0.05)^3)
        phi_02 = eta * ((1.0 - 0.02)^3)
        phi_00 = eta # sticky limit
        @printf(io, "%.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n", eta, tau, phi_10, phi_05, phi_02, phi_00)
    end
end
println("Saved Phase Diagram (phi, tau) -> $file_tau_phi")

# ------------------------------------------------------------------------------
# 3. Phase Diagram in Physical Plane (phi, T* = kT / |u0|)
# ------------------------------------------------------------------------------
file_T_phi = joinpath(out_dir, "phase_diagram_T_vs_phi.dat")
open(file_T_phi, "w") do io
    println(io, "# eta  tau  phi_10  T_menon_10  T_conv_10  phi_05  T_menon_05  T_conv_05  phi_02  T_menon_02  T_conv_02")
    for eta in eta_fine
        tau = eval_binodal(eta)
        
        # eps = 0.10
        eps_10 = 0.10
        phi_10 = eta * ((1.0 - eps_10)^3)
        T_m_10 = 1.0 / log(1.0 / (12.0 * eps_10 * tau))
        term_c_10 = 4.0 * tau * (((1.0 / (1.0 - eps_10))^3) - 1.0)
        T_c_10 = 1.0 / log(1.0 + 1.0 / term_c_10)
        
        # eps = 0.05
        eps_05 = 0.05
        phi_05 = eta * ((1.0 - eps_05)^3)
        T_m_05 = 1.0 / log(1.0 / (12.0 * eps_05 * tau))
        term_c_05 = 4.0 * tau * (((1.0 / (1.0 - eps_05))^3) - 1.0)
        T_c_05 = 1.0 / log(1.0 + 1.0 / term_c_05)
        
        # eps = 0.02
        eps_02 = 0.02
        phi_02 = eta * ((1.0 - eps_02)^3)
        T_m_02 = 1.0 / log(1.0 / (12.0 * eps_02 * tau))
        term_c_02 = 4.0 * tau * (((1.0 / (1.0 - eps_02))^3) - 1.0)
        T_c_02 = 1.0 / log(1.0 + 1.0 / term_c_02)
        
        @printf(io, "%.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n",
                eta, tau, phi_10, T_m_10, T_c_10, phi_05, T_m_05, T_c_05, phi_02, T_m_02, T_c_02)
    end
end
println("Saved Phase Diagram (phi, T*) -> $file_T_phi")

# ------------------------------------------------------------------------------
# 4. Critical Point Trajectory (phi_c, T*_c) vs Well Width epsilon
# ------------------------------------------------------------------------------
eta_c = 3.0 / sqrt(2.0) - 2.0
tau_c = (2.0 - sqrt(2.0)) / 6.0

file_crit = joinpath(out_dir, "phase_diagram_critical_loci.dat")
open(file_crit, "w") do io
    println(io, "# epsilon  phi_c  T_c_menon  T_c_conv  Delta_a_over_sigma")
    for eps in range(0.001, 0.15, length=500)
        phi_c = eta_c * ((1.0 - eps)^3)
        T_c_m = 1.0 / log(1.0 / (12.0 * eps * tau_c))
        
        term_c = 4.0 * tau_c * (((1.0 / (1.0 - eps))^3) - 1.0)
        T_c_c = 1.0 / log(1.0 + 1.0 / term_c)
        
        a_over_sigma = 1.0 / (1.0 - eps)
        delta_sigma = a_over_sigma - 1.0
        
        @printf(io, "%.6f  %.6f  %.6f  %.6f  %.6f\n", eps, phi_c, T_c_m, T_c_c, delta_sigma)
    end
end
println("Saved Critical Points vs epsilon -> $file_crit")
