# ==============================================================================
# Benchmark Data Generator for Liu, Chen & Chen (2005) - Pure Julia
# AnalyticalStructureFactors.jl
# Focus: Static Structure Factors S(Q), Form Factors P(Q), and SANS Intensities I(Q)
# ==============================================================================

using AnalyticalStructureFactors
using DelimitedFiles
using Printf

const BASE_DIR = @__DIR__

"""
    write_dat(filepath, x, y; header="# x y")
Helper function to write two columns of data to a text file.
"""
function write_dat(filepath, x, y; header="# x y")
    open(filepath, "w") do io
        println(io, header)
        for (xi, yi) in zip(x, y)
            @printf(io, "%.8e\t%.8e\n", xi, yi)
        end
    end
end

"""
    get_Icluster(ϕ, K1, z1, K2, z2)
Calculates the cluster peak intensity max_{q in (0, 3.5)} S(q).
"""
function get_Icluster(ϕ, K1, z1, K2, z2)
    sq_func = S_SALR_MSA(ϕ, K1, z1, K2, z2)
    q_grid = range(0.01, 3.5, length=200)
    return maximum(sq_func.(q_grid))
end

println("Starting generation of theoretical structure factors for Liu (2005)...")

# ------------------------------------------------------------------------------
# FIGURE 1: Effect of Attraction Strength K1
# ------------------------------------------------------------------------------
println("Generating Fig 1 theory...")
q_grid_fig1 = range(0.02, 20.0, length=300)
for (K1, tag) in [(0.0, "K1_0"), (3.0, "K1_3"), (6.0, "K1_6"), (10.0, "K1_10")]
    sq_func = S_SALR_MSA(0.20, K1, 10.0, -1.0, 0.5)
    sq = sq_func.(q_grid_fig1)
    write_dat(joinpath(BASE_DIR, "Fig1", "theory_Fig1a_$(tag).dat"), q_grid_fig1, sq, header="# Q S(Q)")
end

k1_grid_fig1b = range(5.0, 21.0, length=70)
icl_fig1b = [get_Icluster(0.20, K1, 10.0, -1.0, 0.5) for K1 in k1_grid_fig1b]
write_dat(joinpath(BASE_DIR, "Fig1", "theory_Fig1b.dat"), k1_grid_fig1b, icl_fig1b, header="# K1 I_cluster")

# ------------------------------------------------------------------------------
# FIGURE 2: Effect of Attraction Range 1/Z1
# ------------------------------------------------------------------------------
println("Generating Fig 2 theory...")
q_grid_fig2 = range(0.02, 20.0, length=300)
for (Z1, tag) in [(14.0, "Z1_14"), (8.0, "Z1_8"), (4.0, "Z1_4")]
    sq_func = S_SALR_MSA(0.20, 6.0, Z1, -1.0, 0.5)
    sq = sq_func.(q_grid_fig2)
    write_dat(joinpath(BASE_DIR, "Fig2", "theory_Fig2a_$(tag).dat"), q_grid_fig2, sq, header="# Q S(Q)")
end

inv_z1_grid_fig2b = range(0.06, 0.26, length=60)
icl_fig2b = [get_Icluster(0.20, 6.0, 1.0 / inv_z1, -1.0, 0.5) for inv_z1 in inv_z1_grid_fig2b]
write_dat(joinpath(BASE_DIR, "Fig2", "theory_Fig2b.dat"), inv_z1_grid_fig2b, icl_fig2b, header="# 1/Z1 I_cluster")

# ------------------------------------------------------------------------------
# FIGURE 3: Effect of Repulsion Strength |K2|
# ------------------------------------------------------------------------------
println("Generating Fig 3 theory...")
q_grid_fig3 = range(0.02, 20.0, length=300)
for (K2, tag) in [(-1.0, "K2_-1"), (-0.1, "K2_-0.1"), (-0.01, "K2_-0.01")]
    sq_func = S_SALR_MSA(0.20, 6.9, 10.0, K2, 0.5)
    sq = sq_func.(q_grid_fig3)
    write_dat(joinpath(BASE_DIR, "Fig3", "theory_Fig3a_$(tag).dat"), q_grid_fig3, sq, header="# Q S(Q)")
end

k2_grid_fig3b = range(0.01, 2.0, length=60)
icl_fig3b = [get_Icluster(0.20, 6.9, 10.0, -k2, 0.5) for k2 in k2_grid_fig3b]
write_dat(joinpath(BASE_DIR, "Fig3", "theory_Fig3b.dat"), k2_grid_fig3b, icl_fig3b, header="# |K2| I_cluster")

# ------------------------------------------------------------------------------
# FIGURE 4: Effect of Repulsion Range 1/Z2
# ------------------------------------------------------------------------------
println("Generating Fig 4 theory...")
q_grid_fig4 = range(0.02, 20.0, length=300)
for (Z2, tag) in [(0.1, "Z2_0.1"), (2.0, "Z2_2"), (4.0, "Z2_4"), (8.0, "Z2_8")]
    sq_func = S_SALR_MSA(0.20, 6.9, 10.0, -1.0, Z2)
    sq = sq_func.(q_grid_fig4)
    write_dat(joinpath(BASE_DIR, "Fig4", "theory_Fig4a_$(tag).dat"), q_grid_fig4, sq, header="# Q S(Q)")
end

inv_z2_grid_fig4b = 10.0 .^ range(-1.0, 1.1, length=60)
icl_fig4b_circles = Float64[]
x_circles = Float64[]
sq0_fig4b_stars = Float64[]
x_stars = Float64[]

for inv_z2 in inv_z2_grid_fig4b
    z2 = 1.0 / inv_z2
    sq_func = S_SALR_MSA(0.20, 6.9, 10.0, -1.0, z2)
    q_test = range(0.01, 3.5, length=200)
    sq_test = sq_func.(q_test)
    max_idx = argmax(sq_test)
    if max_idx > 1 && q_test[max_idx] > 0.1
        push!(x_circles, inv_z2)
        push!(icl_fig4b_circles, sq_test[max_idx])
    else
        push!(x_stars, inv_z2)
        push!(sq0_fig4b_stars, sq_func(0.001))
    end
end
write_dat(joinpath(BASE_DIR, "Fig4", "theory_Fig4b_circles.dat"), x_circles, icl_fig4b_circles, header="# 1/Z2 I_cluster")
write_dat(joinpath(BASE_DIR, "Fig4", "theory_Fig4b_stars.dat"), x_stars, sq0_fig4b_stars, header="# 1/Z2 S(0)")

# ------------------------------------------------------------------------------
# FIGURE 5: Effect of Volume Fraction phi
# ------------------------------------------------------------------------------
println("Generating Fig 5 theory...")
q_grid_fig5 = range(0.02, 20.0, length=300)
for (phi, tag) in [(0.05, "phi_0.05"), (0.20, "phi_0.2"), (0.40, "phi_0.4"), (0.55, "phi_0.55")]
    sq_func = S_SALR_MSA(phi, 10.0, 10.0, -1.0, 0.5)
    sq = sq_func.(q_grid_fig5)
    write_dat(joinpath(BASE_DIR, "Fig5", "theory_Fig5a_$(tag).dat"), q_grid_fig5, sq, header="# Q S(Q)")
end

phi_grid_fig5b = range(0.02, 0.52, length=60)
icl_fig5b = [get_Icluster(phi, 10.0, 10.0, -1.0, 0.5) for phi in phi_grid_fig5b]
write_dat(joinpath(BASE_DIR, "Fig5", "theory_Fig5b.dat"), phi_grid_fig5b, icl_fig5b, header="# phi I_cluster")

# ------------------------------------------------------------------------------
# FIGURE 7: Static Structure Factors at Selected State Points (phi=0.15)
# ------------------------------------------------------------------------------
println("Generating Fig 7 theory...")
q_grid_fig7b = range(0.02, 15.0, length=300)
for (K1, tag) in [(3.6, "K1_3.6"), (4.5, "K1_4.5"), (5.5, "K1_5.5"), (6.72, "K1_6.72")]
    sq_func = S_SALR_MSA(0.15, K1, 10.0, -0.3, 2.0)
    sq = sq_func.(q_grid_fig7b)
    write_dat(joinpath(BASE_DIR, "Fig7", "theory_Fig7b_$(tag).dat"), q_grid_fig7b, sq, header="# Q S(Q)")
end

# ------------------------------------------------------------------------------
# FIGURE 9: Short-Range Repulsion and Long-Range Attraction
# ------------------------------------------------------------------------------
println("Generating Fig 9 theory...")
q_grid_fig9 = 10.0 .^ range(-1.3, 1.0, length=250)
for (K1, tag) in [(0.0, "K1_0"), (0.2, "K1_0.2"), (0.4, "K1_0.4"), (0.45, "K1_0.45")]
    sq_func = S_SALR_MSA(0.20, K1, 0.5, -2.0, 2.0)
    sq = sq_func.(q_grid_fig9)
    write_dat(joinpath(BASE_DIR, "Fig9", "theory_Fig9_$(tag).dat"), q_grid_fig9, sq, header="# Q S(Q)")
end

# ------------------------------------------------------------------------------
# FIGURES 11 & 12: Cytochrome C SANS Modeling (10.18 wt% and 20.40 wt%)
# ------------------------------------------------------------------------------
println("Generating Fig 11 & 12 Cytochrome C theory...")
function P_sphere(Q, a)
    qa = max(Q * a, 1e-6)
    return (3.0 * (sin(qa) - qa * cos(qa)) / (qa^3))^2
end

# Fig 11: 10.18 wt%, a = 14.8 A, sigma = 29.6 A, K1 = 0.33, Z1 = 0.51, K2 = -2.20, Z2 = 1.85, phi = 0.075
a11 = 14.8
sigma11 = 2.0 * a11
K1_11 = 0.33
Z1_11 = 0.51
K2_11 = -2.20
Z2_11 = 1.85
phi_11 = 0.075

q_grid_sans = range(0.002, 0.30, length=200)
p11 = [P_sphere(q, a11) for q in q_grid_sans]
sq_func_11 = S_SALR_MSA(phi_11, K1_11, Z1_11, K2_11, Z2_11)
s11 = [sq_func_11(q * sigma11) for q in q_grid_sans]
i11 = 1.064 .* p11 .* s11

write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11a_Iq.dat"), q_grid_sans, i11, header="# Q I(Q)")
write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11b_Pq.dat"), q_grid_sans, p11, header="# Q P(Q)")
write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11b_Sq.dat"), q_grid_sans, s11, header="# Q S(Q)")

r_grid_11 = range(sigma11, 100.0, length=150)
v_att_11 = [-K1_11 * exp(-Z1_11 * (r/sigma11 - 1.0)) / (r/sigma11) for r in r_grid_11]
v_rep_11 = [-K2_11 * exp(-Z2_11 * (r/sigma11 - 1.0)) / (r/sigma11) for r in r_grid_11]
v_tot_11 = v_att_11 .+ v_rep_11

write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11c_vatt.dat"), r_grid_11, v_att_11, header="# Distance(A) V_att/kBT")
write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11c_vrep.dat"), r_grid_11, v_rep_11, header="# Distance(A) V_rep/kBT")
write_dat(joinpath(BASE_DIR, "Fig11", "theory_Fig11c_vtot.dat"), r_grid_11, v_tot_11, header="# Distance(A) V_tot/kBT")

# Fig 12: 20.40 wt%, a = 14.9 A, sigma = 29.8 A, K1 = 0.18, Z1 = 0.39, K2 = -1.95, Z2 = 3.00, phi = 0.15
a12 = 14.9
sigma12 = 2.0 * a12
K1_12 = 0.18
Z1_12 = 0.39
K2_12 = -1.95
Z2_12 = 3.00
phi_12 = 0.15

p12 = [P_sphere(q, a12) for q in q_grid_sans]
sq_func_12 = S_SALR_MSA(phi_12, K1_12, Z1_12, K2_12, Z2_12)
s12 = [sq_func_12(q * sigma12) for q in q_grid_sans]
i12 = 2.119 .* p12 .* s12

write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12a_Iq.dat"), q_grid_sans, i12, header="# Q I(Q)")
write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12b_Pq.dat"), q_grid_sans, p12, header="# Q P(Q)")
write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12b_Sq.dat"), q_grid_sans, s12, header="# Q S(Q)")

r_grid_12 = range(sigma12, 100.0, length=150)
v_att_12 = [-K1_12 * exp(-Z1_12 * (r/sigma12 - 1.0)) / (r/sigma12) for r in r_grid_12]
v_rep_12 = [-K2_12 * exp(-Z2_12 * (r/sigma12 - 1.0)) / (r/sigma12) for r in r_grid_12]
v_tot_12 = v_att_12 .+ v_rep_12

write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12c_vatt.dat"), r_grid_12, v_att_12, header="# Distance(A) V_att/kBT")
write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12c_vrep.dat"), r_grid_12, v_rep_12, header="# Distance(A) V_rep/kBT")
write_dat(joinpath(BASE_DIR, "Fig12", "theory_Fig12c_vtot.dat"), r_grid_12, v_tot_12, header="# Distance(A) V_tot/kBT")

println("All theoretical structure factor curves successfully generated!")
