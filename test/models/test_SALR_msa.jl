# test/models/test_SALR_msa.jl

using Test
using AnalyticalStructureFactors

@testset "Exact Two-Yukawa (SALR) MSA Tests" begin
    ϕ = 0.15
    z1 = 2.0
    K1 = 1.0  # Short-range attraction
    z2 = 0.5
    K2 = -0.5 # Long-range repulsion

    @testset "Hard-sphere and Single-Yukawa limits" begin
        k_test = 4.0
        # When both K1 = 0 and K2 = 0 -> Hard Sphere PY
        S_hs = S_SALR_MSA(ϕ, 0.0, z1, 0.0, z2, k_test)
        @test isapprox(S_hs, S_HS_PY(ϕ, k_test); atol=1e-6)

        # When K2 = 0 -> Single-Yukawa MSA with K1
        S_y1 = S_SALR_MSA(ϕ, 0.5, z1, 0.0, z2, k_test)
        @test isapprox(S_y1, S_Yukawa_MSA(ϕ, 0.5, z1, k_test); atol=1e-6)

        # When K1 = 0 -> Single-Yukawa MSA with K2
        S_y2 = S_SALR_MSA(ϕ, 0.0, z1, 0.3, z2, k_test)
        @test isapprox(S_y2, S_Yukawa_MSA(ϕ, 0.3, z2, k_test); atol=1e-6)
    end

    @testset "General SALR MSA properties" begin
        # Positivity and finiteness
        for k_val in [0.0, 0.5, 1.5, 3.5, 7.0, 15.0]
            S_val = S_SALR_MSA(ϕ, K1, z1, K2, z2, k_val)
            @test S_val > 0.0
            @test isfinite(S_val)
        end

        # Direct correlation consistency
        k_val = 3.0
        S_val = S_SALR_MSA(ϕ, K1, z1, K2, z2, k_val)
        IS_val = IS_SALR_MSA(ϕ, K1, z1, K2, z2, k_val)
        C_val = C_SALR_MSA(ϕ, K1, z1, K2, z2, k_val)

        @test isapprox(IS_val, 1.0 / S_val; atol=1e-12)
        @test isapprox(C_val, (1.0 - IS_val) / (24.0 * ϕ); atol=1e-12)

        # High-k limit: S(k) -> 1
        @test isapprox(S_SALR_MSA(ϕ, K1, z1, K2, z2, 100.0), 1.0; atol=1e-2)

        # Factory closure and vectorized evaluation
        S_func = S_SALR_MSA(ϕ, K1, z1, K2, z2)
        @test isapprox(S_func(2.5), S_SALR_MSA(ϕ, K1, z1, K2, z2, 2.5); atol=1e-12)

        k_vec = [1.0, 2.0, 4.0, 8.0]
        S_vec = S_SALR_MSA(ϕ, K1, z1, K2, z2, k_vec)
        @test length(S_vec) == 4
        @test isapprox(S_vec[2], S_SALR_MSA(ϕ, K1, z1, K2, z2, 2.0); atol=1e-12)

        # Domain errors
        @test_throws DomainError S_SALR_MSA(-0.1, K1, z1, K2, z2, 1.0)
        @test_throws DomainError S_SALR_MSA(1.1, K1, z1, K2, z2, 1.0)
        @test_throws DomainError S_SALR_MSA(ϕ, K1, -z1, K2, z2, 1.0)
        @test_throws DomainError S_SALR_MSA(ϕ, K1, z1, K2, -z2, 1.0)
        @test_throws DomainError S_SALR_MSA(ϕ, K1, z1, K2, z2, -1.0)
    end
end
