# test/test_colloidal_potentials.jl

using Test
using AnalyticalStructureFactors

@testset "Colloidal Potentials & Electrochemistry Tests" begin
    @testset "ionic_strength" begin
        # 1:1 salt (NaCl, 10 mM = 0.01 M)
        c = [0.01, 0.01]
        z = [1.0, -1.0]
        I = ionic_strength(c, z)
        @test isapprox(I, 0.01; atol=1e-12)

        # 2:1 salt (CaCl2, 10 mM = 0.01 M Ca2+, 0.02 M Cl-)
        c2 = [0.01, 0.02]
        z2 = [2.0, -1.0]
        I2 = ionic_strength(c2, z2) # 0.5 * (0.01*4 + 0.02*1) = 0.03 M
        @test isapprox(I2, 0.03; atol=1e-12)
    end

    @testset "bjerrum_length" begin
        # In water at 298.15 K with epsilon_r = 78.4, lambda_B ~ 0.714 nm (7.14e-10 m)
        T_room = 298.15
        eps_water = 78.4
        l_B = bjerrum_length(T_room, eps_water)
        @test isapprox(l_B * 1e9, 0.714; atol=0.01) # in nanometers
    end

    @testset "debye_screening_length and DLVO parameters" begin
        T_room = 298.15
        eps_water = 78.4
        I_molar = 0.01 # 10 mM

        # For 10 mM 1:1 electrolyte, kappa^(-1) ~ 3.04 nm
        inv_kappa = debye_screening_length(T_room, eps_water, I_molar)
        @test isapprox(inv_kappa * 1e9, 3.04; atol=0.1)

        # DLVO parameters for 20 nm diameter sphere with Z = 50 e
        diameter = 20.0e-9
        Z_eff = 50.0
        params = dlvo_yukawa_parameters(Z_eff, diameter, T_room, eps_water, I_molar)
        @test params.K > 0.0
        @test params.z > 0.0
        # z = kappa * sigma = (diameter / inv_kappa)
        @test isapprox(params.z, diameter / inv_kappa; rtol=1e-5)
    end
end
