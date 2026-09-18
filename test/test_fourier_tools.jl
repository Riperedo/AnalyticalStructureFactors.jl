# test/test_fourier_tools.jl

using Test
using AnalyticalStructureFactors

@testset "3D Spherical Fourier Transform and Inversion Tests" begin
    # Test with a known analytical system: Hard Sphere PY at low density
    ϕ = 0.1
    ρ = (6.0 * ϕ) / π # number density
    
    q_grid = collect(range(0.01, 30.0, length=600))
    Sq_grid = [S_HS_PY(ϕ, q) for q in q_grid]
    
    r_grid = collect(range(0.1, 5.0, length=200))
    gr_calc = sk_to_gr(q_grid, Sq_grid, r_grid, ρ)
    
    # Core exclusion: for r < 0.9, g(r) should be close to 0
    idx_core = findall(r -> r < 0.85, r_grid)
    mean_core = sum(abs.(gr_calc[idx_core])) / length(idx_core)
    @test mean_core < 0.25
    
    # Long range asymptote: g(r) -> 1 as r -> large
    @test isapprox(gr_calc[end], 1.0; atol=0.1)
    
    # Forward transformation gr_to_sk
    Sq_reconstructed = gr_to_sk(r_grid, gr_calc, [5.0], ρ)
    @test isapprox(Sq_reconstructed[1], S_HS_PY(ϕ, 5.0); atol=0.2)
end
