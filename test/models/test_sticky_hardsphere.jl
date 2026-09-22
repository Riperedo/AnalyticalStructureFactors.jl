# test/models/test_sticky_hardsphere.jl
#
# Comprehensive test suite for Sticky Hard Sphere (SHS) models
# References:
# [1] R. J. Baxter, J. Chem. Phys. 49(6), 2770–2774 (1968).
# [2] S. V. G. Menon, C. Manohar, and K. S. Rao, J. Chem. Phys. 95(12), 9186–9190 (1991).

using Test
using AnalyticalStructureFactors

@testset "Sticky Hard Spheres (Baxter 1968 & Menon et al. 1991)" begin
    η = 0.25
    τ = 0.8
    k_test = 3.5
    k_vec = collect(0.1:0.5:15.0)

    @testset "Baxter Quadratic Solution & Critical Stickiness" begin
        # Critical tau should be positive for moderate densities
        tau_c = baxter_shs_critical_tau(η)
        @test tau_c > 0.0
        @test isfinite(tau_c)

        # Above spinodal, lambda is well-defined and positive
        lambda_val = baxter_shs_lambda(η, τ)
        @test lambda_val > 0.0
        @test isfinite(lambda_val)

        # Below spinodal, should throw an error
        @test_throws ErrorException baxter_shs_lambda(η, tau_c * 0.5)

        # Domain errors
        @test_throws DomainError baxter_shs_lambda(-0.1, τ)
        @test_throws DomainError baxter_shs_lambda(1.0, τ)
        @test_throws DomainError baxter_shs_lambda(η, -0.5)
        @test_throws DomainError baxter_shs_critical_tau(-0.1)
    end

    @testset "Pure Hard-Sphere Asymptotic Limit (τ → ∞)" begin
        # As stickiness vanishes (tau -> infty), S_Baxter_SHS must converge to S_HS_PY
        tau_large = 1e8
        s_shs_hs = S_Baxter_SHS(η, tau_large, k_test)
        s_hs_py = S_HS_PY(η, k_test)
        @test s_shs_hs ≈ s_hs_py atol=1e-6

        # Check across wavevector vector
        s_shs_vec = S_Baxter_SHS(η, tau_large, k_vec)
        s_py_vec = S_HS_PY.(Ref(η), k_vec)
        @test s_shs_vec ≈ s_py_vec atol=1e-5
    end

    @testset "Inverse Structure Factor Consistency" begin
        s_val = S_Baxter_SHS(η, τ, k_test)
        is_val = IS_Baxter_SHS(η, τ, k_test)
        @test s_val * is_val ≈ 1.0 atol=1e-12
        @test s_val > 0.0
    end

    @testset "Small-k (Compressibility) Limit Stability" begin
        k_tiny = 1e-6
        s_zero = S_Baxter_SHS(η, τ, 0.0)
        s_tiny = S_Baxter_SHS(η, τ, k_tiny)
        @test isfinite(s_zero)
        @test isfinite(s_tiny)
        @test s_zero ≈ s_tiny atol=1e-8
        @test s_zero > 0.0
    end

    @testset "Menon et al. (1991) Physical Parameter Mapping" begin
        sigma = 1.0
        delta = 0.05
        u0_kT = -1.2 # Attractive well
        phi = 0.20

        s_menon = S_Menon_SHS(phi, u0_kT, delta, sigma, k_test)
        @test isfinite(s_menon)
        @test s_menon > 0.0

        # Vectorized Menon S(k)
        s_menon_vec = S_Menon_SHS(phi, u0_kT, delta, sigma, k_vec)
        @test length(s_menon_vec) == length(k_vec)
        for (i, k) in enumerate(k_vec)
            @test s_menon_vec[i] ≈ S_Menon_SHS(phi, u0_kT, delta, sigma, k)
        end

        # Domain error check for negative lengths
        @test_throws DomainError S_Menon_SHS(phi, u0_kT, -0.01, sigma, k_test)
        @test_throws DomainError S_Menon_SHS(phi, u0_kT, delta, -1.0, k_test)
    end

    @testset "Vectorization Consistency" begin
        s_vec = S_Baxter_SHS(η, τ, k_vec)
        is_vec = IS_Baxter_SHS(η, τ, k_vec)
        @test length(s_vec) == length(k_vec)
        @test length(is_vec) == length(k_vec)

        for (i, k) in enumerate(k_vec)
            @test s_vec[i] ≈ S_Baxter_SHS(η, τ, k)
            @test is_vec[i] ≈ IS_Baxter_SHS(η, τ, k)
            @test s_vec[i] * is_vec[i] ≈ 1.0 atol=1e-12
        end
    end
end
