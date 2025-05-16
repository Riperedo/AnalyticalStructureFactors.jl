# test/models/mixtures/test_yukawa_rpa_mixture.jl

# Asumimos que este archivo será incluido por test/runtests.jl,
# donde AnalyticalStructureFactors, Test, y LinearAlgebra ya están importados.

@testset "Yukawa RPA Mixture Model (S_RPA_mixture_Yukawa)" begin

    # --- Parámetros Comunes para las Pruebas ---
    phi_test_yuk = 0.15  # Fracción de volumen de ejemplo
    k_test_val_yuk = 2.0 # Un valor k de ejemplo
    k_test_array_yuk = [0.2, 0.8, 1.5, 3.0, 6.0] # Un array de k de ejemplo

    # Parámetros para una mezcla binaria
    sigma_vec_bin_yuk = [1.0, 0.8]
    phi_vec_bin_yuk = [0.1, 0.05] # Suma = phi_test_yuk = 0.15

    # Parámetros del potencial de Yukawa para la mezcla binaria (ejemplos)
    # Matriz de amplitudes A_ij
    A_matrix_example = [1.0  0.5;  # A11, A12
                        0.5  0.8]  # A21, A22 (A_ij = A_ji)
    
    # Matriz de inversas de longitudes de apantallamiento Z_ij
    Z_matrix_example = [1.8  2.0;  # Z11, Z12
                        2.0  2.2]  # Z21, Z22 (Z_ij = Z_ji)

    @testset "Limit: Reference Hard-Sphere System (Zero Amplitude)" begin
        # Si todas las amplitudes A_ij son cero,
        # S_RPA_mixture_Yukawa debería reducirse a S_HS_VW_mixture.
        
        A_matrix_zero = zeros(Float64, 2, 2)

        # Prueba para un solo k
        s_rpa_ref_scalar = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_test_val_yuk, 
                                                A_matrix_zero, Z_matrix_example)
        s_hsvw_ref_scalar = S_HS_VW_mixture(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_test_val_yuk)
        @test s_rpa_ref_scalar ≈ s_hsvw_ref_scalar atol=1e-7 # Tolerancia puede necesitar ajuste

        # Prueba con el array de k
        s_rpa_ref_array = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, 
                                               A_matrix_zero, Z_matrix_example, 
                                               k_test_array_yuk)
        s_hsvw_ref_array = S_HS_VW_mixture(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_test_array_yuk)
        @test s_rpa_ref_array ≈ s_hsvw_ref_array atol=1e-7
    end

    @testset "Limit: Single Component Yukawa RPA" begin
        # Configurar una mezcla de un solo componente
        sigma_sc_yuk = [1.0]
        phi_sc_val_yuk = 0.2
        phi_sc_vec_yuk = [phi_sc_val_yuk]
        
        A_sc_yuk = 0.5  # Amplitud para el componente único
        Z_sc_yuk = 2.0  # Z para el componente único
        
        A_matrix_sc_yuk = fill(A_sc_yuk, (1,1))
        Z_matrix_sc_yuk = fill(Z_sc_yuk, (1,1))

        for k_val in k_test_array_yuk
            S_matrix_rpa_sc = S_RPA_mixture_Yukawa(sigma_sc_yuk, phi_sc_vec_yuk, k_val, 
                                                   A_matrix_sc_yuk, Z_matrix_sc_yuk)
            @test size(S_matrix_rpa_sc) == (1,1)
            @test eltype(S_matrix_rpa_sc) <: Real
            s11_rpa_mixture = S_matrix_rpa_sc[1,1]

            # S(k) del sistema monodisperso Yukawa-RPA (de models/yukawa.jl)
            # Asumimos que S_Yukawa_RPA toma (ϕ, A, Z, k)
            s_mono_yuk_rpa = S_Yukawa_RPA(phi_sc_val_yuk, A_sc_yuk, Z_sc_yuk, k_val) 

            @test s11_rpa_mixture ≈ s_mono_yuk_rpa atol=1e-7
        end
    end

    @testset "Binary Mixture Properties (Yukawa RPA)" begin
        # Usar los parámetros de ejemplo definidos al principio del testset padre
        for k_val in k_test_array_yuk
            S_matrix_bin_yuk = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_val,
                                                    A_matrix_example, Z_matrix_example)
            
            @test size(S_matrix_bin_yuk) == (2,2)
            @test eltype(S_matrix_bin_yuk) <: Real
            
            # Simetría: S_ij(k) = S_ji(k)
            @test S_matrix_bin_yuk[1,2] ≈ S_matrix_bin_yuk[2,1] atol=1e-7
            
            # Positividad de los elementos diagonales (puede fallar si la RPA predice inestabilidades)
            # Si los parámetros de prueba causan S(k) < 0, esta prueba fallará.
            # Ajusta los parámetros de prueba o la condición si los S(k) < 0 son esperados/útiles.
            @test S_matrix_bin_yuk[1,1] > 0.0
            @test S_matrix_bin_yuk[2,2] > 0.0
        end
    end
    
    @testset "Factory and Vectorized Versions (Yukawa RPA)" begin
        # Usar los parámetros de ejemplo definidos al principio
        s_k_func_mixture_yuk = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, 
                                                    A_matrix_example, Z_matrix_example)
        
        S_matrix_direct = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_test_val_yuk,
                                               A_matrix_example, Z_matrix_example)
        S_matrix_factory = s_k_func_mixture_yuk(k_test_val_yuk)
        @test S_matrix_factory ≈ S_matrix_direct atol=1e-9

        S_matrices_vectorized = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, 
                                                     A_matrix_example, Z_matrix_example, 
                                                     k_test_array_yuk)
        @test length(S_matrices_vectorized) == length(k_test_array_yuk)
        for i_k in 1:length(k_test_array_yuk)
            S_direct_i = S_RPA_mixture_Yukawa(sigma_vec_bin_yuk, phi_vec_bin_yuk, k_test_array_yuk[i_k],
                                              A_matrix_example, Z_matrix_example)
            @test S_matrices_vectorized[i_k] ≈ S_direct_i atol=1e-9
        end
    end

end
