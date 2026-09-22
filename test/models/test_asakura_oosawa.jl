# test/models/test_asakura_oosawa.jl
#
# Comprehensive test suite for the Asakura-Oosawa (AO) Colloid-Polymer mixture model
# Reference: M. Schmidt et al., J. Phys.: Condens. Matter 14, 9353–9377 (2002).

using Test
using AnalyticalStructureFactors

@testset "Asakura-Oosawa (AO) Colloid-Polymer Model" begin
    # Test physical parameters
    η_c = 0.2
    η_p = 0.1
    q = 0.6 # Rp / Rc
    k_test = 2.5
    k_vec = collect(0.1:0.5:15.0)

    @testset "Free Volume Fraction and Reservoir Relations" begin
        # Alpha should be in (0, 1) for valid packing fractions
        α = free_volume_fraction_AO(η_c, q)
        @test 0.0 < α < 1.0

        # When η_c -> 0, α -> 1
        @test free_volume_fraction_AO(0.0, q) ≈ 1.0

        # Reservoir conversions should be inverse operations
        η_p_res = polymer_system_to_reservoir_density_AO(η_p, η_c, q)
        @test polymer_reservoir_to_system_density_AO(η_p_res, η_c, q) ≈ η_p

        # Domain error check
        @test_throws DomainError free_volume_fraction_AO(-0.1, q)
        @test_throws DomainError free_volume_fraction_AO(1.0, q)
    end

    @testset "Spinodal Line" begin
        η_p_spin = spinodal_AO(η_c, q)
        @test η_p_spin > 0.0
        @test isfinite(η_p_spin)

        @test_throws DomainError spinodal_AO(0.0, q)
        @test_throws DomainError spinodal_AO(1.0, q)
    end

    @testset "Asymptotic Limit: Zero Polymer (η_p → 0)" begin
        # As η_p -> 0, S_cc should recover pure Hard-Sphere Percus-Yevick
        s_cc_zero = S_AO_cc(η_c, 0.0, q, k_test)
        s_hs_py = S_HS_PY(η_c, k_test)
        @test s_cc_zero ≈ s_hs_py atol=1e-12

        # S_cp should vanish
        s_cp_zero = S_AO_cp(η_c, 0.0, q, k_test)
        @test abs(s_cp_zero) < 1e-12

        # S_pp should be ideal gas value (1.0)
        s_pp_zero = S_AO_pp(η_c, 0.0, q, k_test)
        @test s_pp_zero ≈ 1.0 atol=1e-12

        # S_eff_cc should recover pure HS
        s_eff_zero = S_AO_eff_cc(η_c, 0.0, q, k_test)
        @test s_eff_zero ≈ s_hs_py atol=1e-12
    end

    @testset "Partial Structure Factor Matrix S_AO_matrix" begin
        S_mat = S_AO_matrix(η_c, η_p, q, k_test)
        @test size(S_mat) == (2, 2)
        @test S_mat[1, 1] == S_AO_cc(η_c, η_p, q, k_test)
        @test S_mat[2, 2] == S_AO_pp(η_c, η_p, q, k_test)
        @test S_mat[1, 2] == S_AO_cp(η_c, η_p, q, k_test)
        @test S_mat[2, 1] == S_mat[1, 2] # Matrix symmetry

        # Stability check: diagonal elements positive
        @test S_mat[1, 1] > 0.0
        @test S_mat[2, 2] > 0.0

        # Vectorized matrix
        S_mat_vec = S_AO_matrix(η_c, η_p, q, k_vec)
        @test length(S_mat_vec) == length(k_vec)
        @test S_mat_vec[1] ≈ S_AO_matrix(η_c, η_p, q, k_vec[1])
    end

    @testset "Effective Structure Factors" begin
        s_eff_cc = S_AO_eff_cc(η_c, η_p, q, k_test)
        s_eff_pp = S_AO_eff_pp(η_c, η_p, q, k_test)

        @test isfinite(s_eff_cc)
        @test s_eff_cc > 0.0
        @test isfinite(s_eff_pp)
        @test s_eff_pp > 0.0
    end

    @testset "Verlet-Weis Core Correction Toggle" begin
        s_cc_novw = S_AO_cc(η_c, η_p, q, k_test; VW=false)
        s_cc_vw = S_AO_cc(η_c, η_p, q, k_test; VW=true)
        @test s_cc_novw != s_cc_vw
        @test isfinite(s_cc_vw)
        @test s_cc_vw > 0.0
    end

    @testset "Small k Limit Stability (k → 0)" begin
        k_tiny = 1e-6
        s_cc_tiny = S_AO_cc(η_c, η_p, q, k_tiny)
        s_cp_tiny = S_AO_cp(η_c, η_p, q, k_tiny)
        s_pp_tiny = S_AO_pp(η_c, η_p, q, k_tiny)

        @test isfinite(s_cc_tiny)
        @test isfinite(s_cp_tiny)
        @test isfinite(s_pp_tiny)
        @test !isnan(s_cc_tiny)
        @test !isnan(s_cp_tiny)
        @test !isnan(s_pp_tiny)
    end

    @testset "Vectorization Consistency" begin
        s_cc_vec = S_AO_cc(η_c, η_p, q, k_vec)
        s_cp_vec = S_AO_cp(η_c, η_p, q, k_vec)
        s_pp_vec = S_AO_pp(η_c, η_p, q, k_vec)
        s_eff_vec = S_AO_eff_cc(η_c, η_p, q, k_vec)

        for (i, k) in enumerate(k_vec)
            @test s_cc_vec[i] ≈ S_AO_cc(η_c, η_p, q, k)
            @test s_cp_vec[i] ≈ S_AO_cp(η_c, η_p, q, k)
            @test s_pp_vec[i] ≈ S_AO_pp(η_c, η_p, q, k)
            @test s_eff_vec[i] ≈ S_AO_eff_cc(η_c, η_p, q, k)
        end
    end

    @testset "Direct Correlation Functions in Fourier Space" begin
        C_hs = C_AO_HS(η_c, k_test)
        C_star = C_AO_star(η_c, q, k_test)
        C_cp = C_AO_cp(η_c, q, k_test)
        C_cc = C_AO_cc(η_c, η_p, q, k_test)
        C_eff = C_AO_eff_cc(η_c, η_p, q, k_test)

        @test isfinite(C_hs)
        @test isfinite(C_star)
        @test isfinite(C_cp)
        @test isfinite(C_cc)
        @test isfinite(C_eff)

        # C_cc should be linear combination
        @test C_cc ≈ C_hs + η_p * C_star
    end
end
