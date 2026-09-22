# test/data_test/Menon1991/bisection_spinodal_analysis.jl
#
# Numerical Spinodal and Unphysical State Analysis for Sticky Hard Spheres
# via Successive Bisection on S(k=0.01) and Baxter's Discriminant.
#
# Investigates:
# 1. Divergence of S(k=0.01) near the critical point and spinodal.
# 2. Emergence of unphysical negative A(k) (and S(k) < 0 in linearized closures).
# 3. Complex λ branch-cut boundary (D < 0).
# 4. Successive bisection tracking and phase diagram state mapping.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__

println("========================================================================")
println("Running Successive Bisection Spinodal & Unphysical State Analysis")
println("========================================================================")

k_probe = 0.01
eta_grid = range(0.01, 0.40, length=40)
n_bisection_steps = 30

# Data structures to store sampled points during bisection
physical_samples = Tuple{Float64, Float64}[]
negative_A_samples = Tuple{Float64, Float64, Float64}[] # (eta, tau, A_val)
complex_lambda_samples = Tuple{Float64, Float64}[]
bisection_boundaries = Tuple{Float64, Float64, Float64}[] # (eta, tau_found_D, tau_found_A0)

# Track convergence trace for 3 representative volume fractions
selected_etas = [0.06, 0.12132034, 0.25]
convergence_traces = Dict{Float64, Vector{Tuple{Int, Float64, Float64, Float64, Symbol}}}()

for eta in eta_grid
    # --- Bisection 1: Boundary where D <= 0 (Spinodal) ---
    tau_hi_D = 0.40
    tau_lo_D = 0.001
    for step in 1:n_bisection_steps
        tau_mid = 0.5 * (tau_hi_D + tau_lo_D)
        
        A_q = eta / 12.0
        B_q = -(tau_mid + eta / (1.0 - eta))
        C_q = (1.0 + eta / 2.0) / ((1.0 - eta)^2)
        disc = B_q^2 - 4.0 * A_q * C_q
        
        if disc < 0.0
            push!(complex_lambda_samples, (eta, tau_mid))
            tau_lo_D = tau_mid
        else
            lambda_val = (-B_q - sqrt(disc)) / (2.0 * A_q)
            A_k, B_k = AnalyticalStructureFactors._shs_AB_functions(eta, lambda_val, k_probe)
            
            if A_k < 0.0
                push!(negative_A_samples, (eta, tau_mid, A_k))
            else
                push!(physical_samples, (eta, tau_mid))
            end
            tau_hi_D = tau_mid
        end
    end
    tau_spinodal_found = 0.5 * (tau_hi_D + tau_lo_D)
    
    # --- Bisection 2: Boundary where A(k=0.01) <= 0 or D < 0 (Compressibility Divergence) ---
    tau_hi_A = 0.40
    tau_lo_A = 0.001
    
    trace_vec = Tuple{Int, Float64, Float64, Float64, Symbol}[]
    is_tracked = any(abs(eta - target) < 1e-4 for target in selected_etas)
    
    for step in 1:n_bisection_steps
        tau_mid = 0.5 * (tau_hi_A + tau_lo_A)
        
        A_q = eta / 12.0
        B_q = -(tau_mid + eta / (1.0 - eta))
        C_q = (1.0 + eta / 2.0) / ((1.0 - eta)^2)
        disc = B_q^2 - 4.0 * A_q * C_q
        
        status = :unknown
        A_k_val = NaN
        S_val = NaN
        
        if disc < 0.0
            status = :complex_lambda
            tau_lo_A = tau_mid
        else
            lambda_val = (-B_q - sqrt(disc)) / (2.0 * A_q)
            A_k, B_k = AnalyticalStructureFactors._shs_AB_functions(eta, lambda_val, k_probe)
            A_k_val = A_k
            S_val = 1.0 / (A_k^2 + B_k^2)
            
            if A_k <= 0.0
                status = :negative_A
                tau_lo_A = tau_mid
            else
                status = :physical
                tau_hi_A = tau_mid
            end
        end
        
        if is_tracked
            push!(trace_vec, (step, tau_mid, A_k_val, S_val, status))
        end
    end
    tau_A0_found = 0.5 * (tau_hi_A + tau_lo_A)
    
    push!(bisection_boundaries, (eta, tau_spinodal_found, tau_A0_found))
    
    if is_tracked
        convergence_traces[eta] = trace_vec
    end
end

# ------------------------------------------------------------------------------
# Save sampled points and boundaries to files
# ------------------------------------------------------------------------------

# 1. Physical sampled points
open(joinpath(out_dir, "bisection_samples_physical.dat"), "w") do io
    println(io, "# eta  tau")
    for pt in physical_samples
        @printf(io, "%.6f  %.6f\n", pt[1], pt[2])
    end
end

# 2. Negative A(k) / unphysical sampled points
open(joinpath(out_dir, "bisection_samples_negative_A.dat"), "w") do io
    println(io, "# eta  tau  A_val")
    for pt in negative_A_samples
        @printf(io, "%.6f  %.6f  %.6e\n", pt[1], pt[2], pt[3])
    end
end

# 3. Complex lambda sampled points
open(joinpath(out_dir, "bisection_samples_complex_lambda.dat"), "w") do io
    println(io, "# eta  tau")
    for pt in complex_lambda_samples
        @printf(io, "%.6f  %.6f\n", pt[1], pt[2])
    end
end

# 4. Bisection found boundaries
open(joinpath(out_dir, "bisection_found_boundaries.dat"), "w") do io
    println(io, "# eta  tau_spinodal_found(D=0)  tau_A0_found(A=0)")
    for pt in bisection_boundaries
        @printf(io, "%.6f  %.6f  %.6f\n", pt[1], pt[2], pt[3])
    end
end

# 5. Exact Analytical Curves: Spinodal (D=0) & Compressibility Divergence (A(0)=0)
open(joinpath(out_dir, "analytical_spinodal_curve.dat"), "w") do io
    println(io, "# eta  tau_spinodal_exact")
    for eta in range(0.005, 0.414, length=500)
        tau_s = baxter_shs_critical_tau(eta)
        @printf(io, "%.6f  %.6f\n", eta, tau_s)
    end
end

open(joinpath(out_dir, "analytical_compressibility_div_curve.dat"), "w") do io
    println(io, "# eta  tau_A0_divergence_exact")
    for eta in range(0.005, 0.50, length=500)
        lam0 = (1.0 + 2.0 * eta) / (eta * (1.0 - eta))
        tau_div = (eta / 12.0) * lam0 - eta / (1.0 - eta) + (1.0 + eta / 2.0) / (((1.0 - eta)^2) * lam0)
        @printf(io, "%.6f  %.6f\n", eta, tau_div)
    end
end

# 6. Convergence Trace
open(joinpath(out_dir, "bisection_convergence_trace.dat"), "w") do io
    println(io, "# eta  step  tau_mid  A_val  S_val  status")
    for (eta_val, trace) in convergence_traces
        for t in trace
            @printf(io, "%.6f  %2d  %.8f  %12.5e  %12.5e  %s\n", eta_val, t[1], t[2], t[3], t[4], string(t[5]))
        end
    end
end

println("Bisection analysis data successfully generated:")
println("  - bisection_samples_physical.dat ($(length(physical_samples)) points)")
println("  - bisection_samples_negative_A.dat ($(length(negative_A_samples)) points)")
println("  - bisection_samples_complex_lambda.dat ($(length(complex_lambda_samples)) points)")
println("  - bisection_found_boundaries.dat")
println("  - analytical_spinodal_curve.dat")
println("  - analytical_compressibility_div_curve.dat")
println("  - bisection_convergence_trace.dat")
