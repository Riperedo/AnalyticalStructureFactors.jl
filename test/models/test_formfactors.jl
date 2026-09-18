# test/models/test_formfactors.jl

using Test
using AnalyticalStructureFactors

@testset "Form Factors and Decoupling Factor Tests" begin
    # 1. Homogeneous Sphere
    @testset "P_sphere" begin
        R = 50.0
        # P(0) == 1
        @test isapprox(P_sphere(0.0, R), 1.0; atol=1e-7)
        @test isapprox(P_sphere(1e-6, R), 1.0; atol=1e-7)
        
        # P(q) is strictly between 0 and 1
        q_vals = range(0.01, 1.0, length=50)
        P_vals = P_sphere(q_vals, R)
        @test all(0.0 .<= P_vals .<= 1.0)
        
        # Factory closure consistency
        P_func = P_sphere(R)
        @test isapprox(P_func(0.1), P_sphere(0.1, R); atol=1e-12)

        # Domain error for negative values
        @test_throws DomainError P_sphere(-1.0, R)
        @test_throws DomainError P_sphere(1.0, -R)
    end

    # 2. Core-Shell Sphere
    @testset "P_coreshell" begin
        Rc = 30.0
        Rs = 50.0
        rho_c = 4.0e-6
        rho_s = 2.0e-6
        rho_solv = 1.0e-6

        # P(0) == 1
        @test isapprox(P_coreshell(0.0, Rc, Rs, rho_c, rho_s, rho_solv), 1.0; atol=1e-7)

        # Homogeneous limit when rho_c == rho_s
        P_homo = P_coreshell(0.05, Rc, Rs, rho_s, rho_s, rho_solv)
        P_sph = P_sphere(0.05, Rs)
        @test isapprox(P_homo, P_sph; atol=1e-10)

        # Vectorized consistency
        q_vec = [0.01, 0.05, 0.1]
        P_vec = P_coreshell(q_vec, Rc, Rs, rho_c, rho_s, rho_solv)
        @test length(P_vec) == 3
        @test isapprox(P_vec[2], P_coreshell(0.05, Rc, Rs, rho_c, rho_s, rho_solv); atol=1e-12)
    end

    # 3. Multi-Shell Concentric Spheres
    @testset "P_multishell" begin
        radii = [20.0, 40.0, 60.0]
        slds = [5.0e-6, 3.0e-6, 1.5e-6]
        sld_solv = 0.5e-6

        # P(0) == 1
        @test isapprox(P_multishell(0.0, radii, slds, sld_solv), 1.0; atol=1e-7)

        # Single shell match with P_sphere
        P_1shell = P_multishell(0.05, [50.0], [3.0e-6], 1.0e-6)
        @test isapprox(P_1shell, P_sphere(0.05, 50.0); atol=1e-10)
    end

    # 4. Ellipsoid of Revolution and Kotlarchyk-Chen Decoupling
    @testset "P_ellipsoid & beta_ellipsoid" begin
        Ra = 60.0
        Rb = 30.0

        # Spherical limit Ra == Rb => beta(q) == 1 and P_ellipsoid == P_sphere
        @test isapprox(P_ellipsoid(0.05, 40.0, 40.0), P_sphere(0.05, 40.0); atol=1e-8)
        @test isapprox(beta_ellipsoid(0.05, 40.0, 40.0), 1.0; atol=1e-8)

        # P(0) == 1 and beta(0) == 1
        @test isapprox(P_ellipsoid(0.0, Ra, Rb), 1.0; atol=1e-7)
        @test isapprox(beta_ellipsoid(0.0, Ra, Rb), 1.0; atol=1e-7)

        # General properties
        q_test = 0.05
        P_val = P_ellipsoid(q_test, Ra, Rb)
        beta_val = beta_ellipsoid(q_test, Ra, Rb)
        @test 0.0 < P_val < 1.0
        @test 0.0 < beta_val <= 1.0

        # Effective radius
        R_eff = effective_radius_ellipsoid(Ra, Rb)
        @test isapprox(R_eff, cbrt(Ra * Rb^2); atol=1e-12)

        # Decoupled intensity
        S_val = 1.2
        I_dec = intensity_decoupled(P_val, beta_val, S_val)
        @test isapprox(I_dec, P_val * (1.0 + beta_val * (S_val - 1.0)); atol=1e-12)
    end
end
