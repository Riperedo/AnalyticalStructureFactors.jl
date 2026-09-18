# test/data_test/CummingsSmith1979/generate_cummings1979_data.jl
#
# Reproduces the analytical curves for Baxter's parameter beta as a function of density eta
# for a single Yukawa fluid (xi = 2.0) along isotherms K = 0.812, 1.0827, 1.2180, and 1.6240,
# as presented in Cummings & Smith (1979) Fig. 1.
#
# This script directly utilizes AnalyticalStructureFactors.jl to evaluate the exact
# analytical MSA roots.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AnalyticalStructureFactors
using Printf

out_dir = @__DIR__
z = 2.0
isotherms = [0.812, 1.0827, 1.2180, 1.6240]
eta_grid = range(0.0001, 0.70, length=2000)

for K in isotherms
    K_str = @sprintf("%.4f", K)
    # Collect real roots directly using AnalyticalStructureFactors.yukawa_msa_beta_roots
    data_all = Tuple{Float64, Vector{Float64}}[]
    for eta in eta_grid
        r = yukawa_msa_beta_roots(eta, K, z)
        if !isempty(r)
            push!(data_all, (eta, r))
        end
    end

    # Write branch data files for clean plotting
    filename = joinpath(out_dir, "analytical_K_$(K_str).dat")
    open(filename, "w") do io
        println(io, "# eta  beta_root1  beta_root2")
        for (eta, r) in data_all
            if length(r) >= 2
                @printf(io, "%.6f  %.6f  %.6f\n", eta, r[1], r[2])
            elseif length(r) == 1
                @printf(io, "%.6f  %.6f  NaN\n", eta, r[1])
            end
        end
    end
    println("Generated analytical data for K = $K_str using AnalyticalStructureFactors.jl -> $filename")
end
