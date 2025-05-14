# test/models/test_SALR.jl

# Asumimos que este archivo será incluido por test/runtests.jl,
# donde AnalyticalStructureFactors y Test ya están importados.

@testset "SALR RPA Model (S_SALR_RPA)" begin

    # --- Parámetros Comunes para las Pruebas ---
    phi_test = 0.2  # Fracción de volumen de ejemplo
    k_test_val = 1.5 # Un valor k de ejemplo
    k_test_array = [0.5, 1.5, 3.0, 5.0, 7.0] # Un array de k de ejemplo

    # Parámetros de ejemplo para un SALR "típico" (necesitarás ajustarlos según tus necesidades)
    amp_attr_typical = 1.0  # Amplitud atractiva
    zA_typical = 1.0      # Inversa de longitud de apantallamiento atractiva
    amp_rep_typical = 0.5   # Amplitud repulsiva
    zR_typical = 0.5      # Inversa de longitud de apantallamiento repulsiva (rango más largo)

    @testset "Limit: Reference Hard-Sphere System (HS-VW)" begin
        # Si ambas amplitudes (atractiva y repulsiva) son cero,
        # S_SALR_RPA debería reducirse a S_HS_VW.
        amp_zero = 0.0

        s_salr_ref = S_SALR_RPA(phi_test, amp_zero, zA_typical, amp_zero, zR_typical, k_test_val)
        s_hsvw_ref = S_HS_VW(phi_test, k_test_val)
        @test s_salr_ref ≈ s_hsvw_ref atol=1e-7

        # Prueba con el array de k
        s_salr_ref_array = S_SALR_RPA(phi_test, amp_zero, zA_typical, amp_zero, zR_typical, k_test_array)
        s_hsvw_ref_array = S_HS_VW.(Ref(phi_test), k_test_array) # Usamos S_HS_VW directamente
        @test s_salr_ref_array ≈ s_hsvw_ref_array atol=1e-7
    end

    @testset "Limit: Single Yukawa Potential (Attractive part only)" begin
        # Si la amplitud repulsiva es cero, debería ser como un Yukawa atractivo sobre HS-VW.
        amp_rep_zero = 0.0
        
        s_salr_yuk_attr = S_SALR_RPA(phi_test, amp_attr_typical, zA_typical, amp_rep_zero, zR_typical, k_test_val)
        
        # Construimos el S(k) para un Yukawa atractivo de referencia
        inv_s0_k_ref = IS_HS_VW(phi_test, k_test_val)
        beta_u_attr_k_ref = betaU_Yukawa(amp_attr_typical, zA_typical, k_test_val)
        density_factor_ref = (6.0 * phi_test) / π
        s_yuk_attr_ref = 1.0 / (inv_s0_k_ref - density_factor_ref * beta_u_attr_k_ref)
        
        @test s_salr_yuk_attr ≈ s_yuk_attr_ref atol=1e-7
    end

    @testset "Limit: Single Yukawa Potential (Repulsive part only)" begin
        # Si la amplitud atractiva es cero, debería ser como un Yukawa repulsivo (con signo negativo en la perturbación) sobre HS-VW.
        amp_attr_zero = 0.0
        
        s_salr_yuk_rep = S_SALR_RPA(phi_test, amp_attr_zero, zA_typical, amp_rep_typical, zR_typical, k_test_val)
        
        # Construimos el S(k) para un Yukawa repulsivo de referencia
        inv_s0_k_ref = IS_HS_VW(phi_test, k_test_val)
        # La perturbación es -beta_u_rep_k
        beta_u_rep_k_ref = betaU_Yukawa(amp_rep_typical, zR_typical, k_test_val)
        density_factor_ref = (6.0 * phi_test) / π
        s_yuk_rep_ref = 1.0 / (inv_s0_k_ref - density_factor_ref * (-beta_u_rep_k_ref)) # Nótese el -beta_u_rep_k_ref
        
        @test s_salr_yuk_rep ≈ s_yuk_rep_ref atol=1e-7
    end

    @testset "General SALR Case" begin
        # Para un caso general, es difícil tener valores exactos sin referencias.
        # Podemos probar propiedades básicas.
        s_k_general = S_SALR_RPA(phi_test, amp_attr_typical, zA_typical, amp_rep_typical, zR_typical, k_test_val)

        # 1. S(k) debe ser positivo (asumiendo parámetros físicos)
        @test s_k_general > 0.0

        # 2. Consistencia con la función fábrica y vectorizada
        s_k_factory_func = S_SALR_RPA(phi_test, amp_attr_typical, zA_typical, amp_rep_typical, zR_typical)
        @test s_k_factory_func(k_test_val) ≈ s_k_general atol=1e-9

        s_k_vectorized = S_SALR_RPA(phi_test, amp_attr_typical, zA_typical, amp_rep_typical, zR_typical, [k_test_val])
        @test s_k_vectorized[1] ≈ s_k_general atol=1e-9

        # 3. Más pruebas cualitativas o con valores de la literatura si los tienes.
        # Por ejemplo, para ciertos parámetros SALR, se espera un pico a k bajos o intermedios.
        # @test S_SALR_RPA(phi_val, KA, ZA, KR, ZR, k_pico_esperado) > S_SALR_RPA(phi_val, KA, ZA, KR, ZR, k_valle_esperado)
    end
    
    @testset "Physicality Checks (S(k) > 0 for a range of k)" begin
        # Verifica que S(k) sea positivo para un rango de k y parámetros típicos.
        s_values = S_SALR_RPA(phi_test, amp_attr_typical, zA_typical, amp_rep_typical, zR_typical, k_test_array)
        for s_val in s_values
            @test s_val > 0.0
        end
    end

end
