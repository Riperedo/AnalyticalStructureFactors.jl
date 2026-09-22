# test/test_smoke_all_models.jl

# Asumimos que este archivo será incluido por test/runtests.jl,
# donde AnalyticalStructureFactors y Test ya están importados.

@testset "Smoke Tests for All S(k) Models" begin

    # --- Parámetros Comunes para las Pruebas de Humo ---
    phi_smoke_test = 0.25 # Fracción de volumen de ejemplo
    # Un array de k adimensionales (k = qσ)
    k_smoke_array = [0.1, 0.5, 1.0, 2.0, 3.0 * π / 1.0, 7.0, 10.0, 15.0] 

    @testset "Hard Sphere - Percus-Yevick (HS-PY)" begin
        try
            s_values = S_HS_PY.(Ref(phi_smoke_test), k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_values) # Todos los S(k) deben ser positivos
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false # Falla si hay cualquier error durante la ejecución
            println("Error en S_HS_PY smoke test: $e")
        end
    end

    @testset "Hard Sphere - Verlet-Weiss (HS-VW)" begin
        try
            s_values = S_HS_VW.(Ref(phi_smoke_test), k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_values)
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false
            println("Error en S_HS_VW smoke test: $e")
        end
    end

    @testset "WCA via Blip Function" begin
        temp_wca_smoke = 1.0 # Temperatura de ejemplo
        nu_wca_smoke = 6     # Parámetro de suavidad
        try
            s_values = S_WCA_blip.(Ref(phi_smoke_test), Ref(temp_wca_smoke), k_smoke_array; ν = nu_wca_smoke)
            @test all(s_val -> s_val > 0.0, s_values)
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false
            println("Error en S_WCA_blip smoke test: $e")
        end
    end

    @testset "Square-Well - RPA (SW-RPA)" begin
        temp_sw_smoke = 1.5  # Temperatura de ejemplo
        lambda_sw_smoke = 1.5 # Rango del pozo de ejemplo
        try
            s_values = S_SW_RPA.(Ref(phi_smoke_test), Ref(temp_sw_smoke), Ref(lambda_sw_smoke), k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_values)
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false
            println("Error en S_SW_RPA smoke test: $e")
        end
    end

    @testset "Yukawa - RPA" begin
        amp_yukawa_smoke = 1.0 # Amplitud de ejemplo (puede ser positiva o negativa)
        z_yukawa_smoke = 1.0   # Inversa de long. de apantallamiento de ejemplo
        try
            s_values = S_Yukawa_RPA.(Ref(phi_smoke_test), Ref(amp_yukawa_smoke), Ref(z_yukawa_smoke), k_smoke_array)
            #@test all(s_val -> s_val > 0.0, s_values)
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false
            println("Error en S_Yukawa_RPA smoke test: $e")
        end
    end

    @testset "SALR - RPA" begin
        amp_attr_salr_smoke = 1.5  # Amplitud atractiva
        zA_salr_smoke = 1.8        # Inversa long. apant. atractiva
        amp_rep_salr_smoke = 0.75  # Amplitud repulsiva
        zR_salr_smoke = 0.8        # Inversa long. apant. repulsiva
        try
            s_values = S_SALR_RPA.(Ref(phi_smoke_test), 
                                    Ref(amp_attr_salr_smoke), Ref(zA_salr_smoke), 
                                    Ref(amp_rep_salr_smoke), Ref(zR_salr_smoke), 
                                    k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_values)
            @test length(s_values) == length(k_smoke_array)
        catch e
            @test false
            println("Error en S_SALR_RPA smoke test: $e")
        end
    end

    @testset "Asakura-Oosawa Colloid-Polymer Model" begin
        eta_c_smoke = 0.2
        eta_p_smoke = 0.08
        q_smoke = 0.5
        try
            s_cc = S_AO_cc(eta_c_smoke, eta_p_smoke, q_smoke, k_smoke_array)
            s_cp = S_AO_cp(eta_c_smoke, eta_p_smoke, q_smoke, k_smoke_array)
            s_pp = S_AO_pp(eta_c_smoke, eta_p_smoke, q_smoke, k_smoke_array)
            s_eff_cc = S_AO_eff_cc(eta_c_smoke, eta_p_smoke, q_smoke, k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_cc)
            @test all(s_val -> s_val > 0.0, s_pp)
            @test all(s_val -> s_val > 0.0, s_eff_cc)
            @test length(s_cc) == length(k_smoke_array)
        catch e
            @test false
            println("Error en Asakura-Oosawa smoke test: $e")
        end
    end

    @testset "Sticky Hard Spheres (Baxter & Menon)" begin
        eta_shs_smoke = 0.2
        tau_shs_smoke = 1.0
        try
            s_baxter = S_Baxter_SHS(eta_shs_smoke, tau_shs_smoke, k_smoke_array)
            s_menon = S_Menon_SHS(0.2, -1.0, 0.05, 1.0, k_smoke_array)
            @test all(s_val -> s_val > 0.0, s_baxter)
            @test all(s_val -> s_val > 0.0, s_menon)
            @test length(s_baxter) == length(k_smoke_array)
            @test length(s_menon) == length(k_smoke_array)
        catch e
            @test false
            println("Error en Sticky Hard Spheres smoke test: $e")
        end
    end

end
