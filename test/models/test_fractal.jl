# test/models/test_fractal.jl

using Test
using AnalyticalStructureFactors
using SpecialFunctions: gamma

@testset "Fractal Structure Factor Tests" begin
    R = 2.5
    D = 2.1
    xi = 100.0

    @testset "Limiting behaviors of S_Fractal_Chen" begin
        # q -> 0 clustering limit: S(0) = 1 + D * Gamma(D) * (xi/R)^D
        expected_S0 = 1.0 + D * gamma(D) * ((xi / R)^D)
        @test isapprox(S_Fractal_Chen(0.0, R, D, xi), expected_S0; rtol=1e-4)
        @test isapprox(S_Fractal_Chen(1e-6, R, D, xi), expected_S0; rtol=1e-3)

        # High-q limit: S(q) -> 1
        @test isapprox(S_Fractal_Chen(50.0, R, D, xi), 1.0; atol=1e-3)

        # Intermediate power-law regime (1/xi << q << 1/R)
        # S(q) ~ q^(-D)
        q1 = 0.05
        q2 = 0.1
        S1 = S_Fractal_Chen(q1, R, D, xi)
        S2 = S_Fractal_Chen(q2, R, D, xi)
        # Ratio of S(q) should scale roughly as (q1/q2)^(-D)
        scaling_ratio = S1 / S2
        expected_ratio = (q1 / q2)^(-D)
        @test isapprox(scaling_ratio, expected_ratio; rtol=0.1)

        # Factory closure and vectorized evaluation
        S_func = S_Fractal_Chen(R, D, xi)
        @test isapprox(S_func(0.1), S_Fractal_Chen(0.1, R, D, xi); atol=1e-12)

        q_vec = [0.01, 0.05, 0.1, 1.0]
        S_vec = S_Fractal_Chen(q_vec, R, D, xi)
        @test length(S_vec) == 4
        @test isapprox(S_vec[3], S_Fractal_Chen(0.1, R, D, xi); atol=1e-12)

        # Domain errors
        @test_throws DomainError S_Fractal_Chen(-0.1, R, D, xi)
        @test_throws DomainError S_Fractal_Chen(0.1, -R, D, xi)
        @test_throws DomainError S_Fractal_Chen(0.1, R, 0.5, xi) # D < 1
        @test_throws DomainError S_Fractal_Chen(0.1, R, 3.5, xi) # D > 3
    end
end
