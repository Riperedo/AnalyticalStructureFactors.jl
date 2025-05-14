# test/test_monodisperse.jl

@testset "Hard Sphere Percus-Yevick (HS-PY) Model" begin

    @testset "Ideal Gas Limit (ϕ -> 0)" begin
        phi_ideal = 0.0
        k_values = [0.0, 0.1, 1.0, 5.0, 10.0] # k = qσ, incluimos k=0

        for k_val in k_values
            # Para phi = 0:
            # IS_HS_PY(0.0, k_val) = 1.0 - 24*0.0*C_HS_PY(0.0, k_val) = 1.0
            # S_HS_PY(0.0, k_val) = 1.0 / IS_HS_PY(0.0, k_val) = 1.0
            @test IS_HS_PY(phi_ideal, k_val) ≈ 1.0 atol=1e-9
            @test S_HS_PY(phi_ideal, k_val) ≈ 1.0 atol=1e-9

            # C_HS_PY(0.0, k) es -_integral1_py(k).
            # Específicamente, C_HS_PY(0.0, 0.0) = -_integral1_py(0.0) = -1/3.
            if k_val == 0.0
                @test C_HS_PY(phi_ideal, 0.0) ≈ -1.0/3.0 atol=1e-7
            else
                # Para k != 0, C_HS_PY(0.0, k_val) no es cero.
                # No necesitamos una prueba específica para su valor exacto aquí
                # si las pruebas para S(k) e IS(k) son suficientes para el gas ideal.
                # Si se quisiera, se podría calcular -_integral1_py(k_val) y comparar.
                # Por ejemplo: C_HS_PY(0.0, 0.1) ≈ -0.333000119...
                # Esta línea es solo para ilustrar, no es una prueba activa:
                # @test C_HS_PY(phi_ideal, k_val) ≈ -((sin(k_val) - k_val*cos(k_val))/k_val^3) atol=1e-7
            end
        end
    end

    @testset "Finite Volume Fraction (e.g., ϕ = 0.1)" begin
        phi_finite = 0.1
        k1 = 1.0 
        k2 = 7.0 # Cerca del primer pico de S(k) para esferas duras

        # 1. Consistencia S(k) = 1 / IS(k)
        @test S_HS_PY(phi_finite, k1) ≈ 1.0 / IS_HS_PY(phi_finite, k1) atol=1e-9
        @test S_HS_PY(phi_finite, k2) ≈ 1.0 / IS_HS_PY(phi_finite, k2) atol=1e-9

        # 2. S(k) debe ser positivo
        @test S_HS_PY(phi_finite, k1) > 0.0
        @test S_HS_PY(phi_finite, k2) > 0.0
        
        # EJEMPLO: Si tuvieras un valor de referencia de un libro/artículo
        # para phi=0.1, k=7.0 (qσ ≈ 7), S(k) ~ 1.52 (ESTE ES UN VALOR INVENTADO PARA EL EJEMPLO)
        # @test S_HS_PY(phi_finite, k2) ≈ 1.52 atol=0.01 
        # ¡NECESITAS REEMPLAZAR ESTO CON VALORES REALES DE LA LITERATURA!
        # Por ejemplo, de Hansen & McDonald (2013), para phi=0.1047 (eta=pi/6 * 0.2), qsigma ~ 6.8, S(qsigma) ~ 1.45
        # Para phi=0.1, qsigma=6.8, S(qsigma) podría ser diferente. Es importante encontrar referencias precisas.
    end

    @testset "Thermodynamic Limit (k -> 0)" begin
        # S(k->0) está relacionado con la compresibilidad isotérmica.
        # Para esferas duras PY, S(0) = (1-ϕ)⁴ / (1+2ϕ)²
        phi_values = [0.05, 0.1, 0.2, 0.3, 0.4]
        for phi_val in phi_values
            expected_S0 = (1.0 - phi_val)^4 / (1.0 + 2.0*phi_val)^2
            @test S_HS_PY(phi_val, 0.0) ≈ expected_S0 atol=1e-7
            @test IS_HS_PY(phi_val, 0.0) ≈ 1.0/expected_S0 atol=1e-7
            # C_HS_PY(phi_val, 0.0) = (_alpha1_py(phi_val)/3 + _alpha2_py(phi_val)/4 + _alpha3_py(phi_val)/6)
            # Y 1/S(0) = 1 - 24*phi*C_HS_PY(phi,0)
            # C_HS_PY(phi_val, 0.0) = (1 - 1/S(0)) / (24*phi_val) si phi_val != 0
            if phi_val > 1e-10 # Evitar división por cero si phi_val es muy pequeño
                 expected_C0 = (1.0 - 1.0/expected_S0) / (24.0*phi_val)
                 @test C_HS_PY(phi_val, 0.0) ≈ expected_C0 atol=1e-7
            else # Caso phi=0 ya cubierto en "Ideal Gas Limit"
                 @test C_HS_PY(phi_val, 0.0) ≈ -1.0/3.0 atol=1e-7
            end
        end
    end
    
    # @testset "Input Validation Warnings" begin
    #     # Esta sección se puede desarrollar más adelante si se desea probar la emisión de advertencias.
    # end

end
