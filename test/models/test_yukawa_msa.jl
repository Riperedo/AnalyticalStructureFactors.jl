# test/models/test_yukawa_msa.jl

using Test
using AnalyticalStructureFactors

@testset "Single-Yukawa MSA Tests" begin
    ϕ = 0.2
    z = 2.0
    K = 0.5

    @testset "Hard-sphere limit when K -> 0" begin
        # As interaction amplitude K -> 0, S_Yukawa_MSA should converge to S_HS_PY
        k_test = 5.0
        S_msa_0 = S_Yukawa_MSA(ϕ, 0.0, z, k_test)
        S_py_0 = S_HS_PY(ϕ, k_test)
        @test isapprox(S_msa_0, S_py_0; atol=1e-6)
    end

    @testset "General properties of S_Yukawa_MSA" begin
        # Positivity
        for k_val in [0.0, 0.1, 1.0, 3.0, 7.0, 15.0]
            S_val = S_Yukawa_MSA(ϕ, K, z, k_val)
            @test S_val > 0.0
            @test isfinite(S_val)
        end

        # Direct correlation function consistency C(k) = (1 - 1/S(k)) / (24 phi)
        k_val = 4.0
        S_val = S_Yukawa_MSA(ϕ, K, z, k_val)
        IS_val = IS_Yukawa_MSA(ϕ, K, z, k_val)
        C_val = C_Yukawa_MSA(ϕ, K, z, k_val)

        @test isapprox(IS_val, 1.0 / S_val; atol=1e-12)
        @test isapprox(C_val, (1.0 - IS_val) / (24.0 * ϕ); atol=1e-12)

        # High-k limit: S(k) -> 1
        @test isapprox(S_Yukawa_MSA(ϕ, K, z, 100.0), 1.0; atol=1e-2)

        # Factory closure and vectorized evaluation
        S_func = S_Yukawa_MSA(ϕ, K, z)
        @test isapprox(S_func(3.0), S_Yukawa_MSA(ϕ, K, z, 3.0); atol=1e-12)

        k_vec = [1.0, 2.0, 5.0, 10.0]
        S_vec = S_Yukawa_MSA(ϕ, K, z, k_vec)
        @test length(S_vec) == 4
        @test isapprox(S_vec[3], S_Yukawa_MSA(ϕ, K, z, 5.0); atol=1e-12)

        # Domain errors
        @test_throws DomainError S_Yukawa_MSA(-0.1, K, z, 1.0)
        @test_throws DomainError S_Yukawa_MSA(1.1, K, z, 1.0)
        @test_throws DomainError S_Yukawa_MSA(ϕ, K, -z, 1.0)
        @test_throws DomainError S_Yukawa_MSA(ϕ, K, z, -1.0)
    end
end
