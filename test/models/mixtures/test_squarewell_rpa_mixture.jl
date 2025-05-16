# test/models/mixtures/test_squarewell_rpa_mixture.jl

# Asumimos que este archivo será incluido por test/runtests.jl,
# donde AnalyticalStructureFactors, Test, y LinearAlgebra ya están importados.

@testset "Square-Well RPA Mixture Model (S_RPA_mixture_SquareWell)" begin

    # --- Parámetros Comunes para las Pruebas ---
    phi_test_sw = 0.15  # Fracción de volumen de ejemplo
    k_test_val_sw = 2.0 # Un valor k de ejemplo
    k_test_array_sw = [0.2, 0.8, 1.5, 3.0, 6.0] # Un array de k de ejemplo

    # Parámetros para una mezcla binaria
    sigma_vec_bin = [1.0, 0.8]
    phi_vec_bin = [0.1, 0.05] # Suma = phi_test_sw = 0.15

    # Parámetros del potencial de Pozo Cuadrado para la mezcla binaria (ejemplos)
    # Matriz de temperaturas T_ij (relacionadas con la profundidad del pozo ε_ij, ej. T_ij = kBT/ε_ij)
    # Si T_ij es muy grande, la perturbación es débil.
    temperature_matrix_example = [1.5  1.2;  # T11, T12
                           1.2  1.8]  # T21, T22
    
    # Matriz de rangos del pozo λ_ij (en unidades de σ)
    lambda_matrix_example = [1.5  1.4;  # λ11, λ12
                             1.4  1.6]  # λ21, λ22

    @testset "Limit: Reference Hard-Sphere System (Infinite Temperature / Zero Well Depth)" begin
        # Si todas las temperaturas T_ij son infinitas (o profundidades de pozo son cero),
        # S_RPA_mixture_SquareWell debería reducirse a S_HS_VW_mixture.
        
        large_temp = 1.0e8 
        temp_matrix_inf = [large_temp large_temp; 
                           large_temp large_temp]

        # Prueba para un solo k
        s_rpa_ref_scalar = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, k_test_val_sw, 
                                                    temp_matrix_inf, lambda_matrix_example)
        s_hsvw_ref_scalar = S_HS_VW_mixture(sigma_vec_bin, phi_vec_bin, k_test_val_sw)
        @test s_rpa_ref_scalar ≈ s_hsvw_ref_scalar atol=1e-6

        # Prueba con el array de k, corrigiendo el orden de los argumentos para la versión vectorizada
        s_rpa_ref_array = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, 
                                                   temp_matrix_inf, lambda_matrix_example, 
                                                   k_test_array_sw)
        # La llamada a S_HS_VW_mixture para un array de k ya era correcta si esa es su firma vectorizada
        s_hsvw_ref_array = S_HS_VW_mixture(sigma_vec_bin, phi_vec_bin, k_test_array_sw)
        @test s_rpa_ref_array ≈ s_hsvw_ref_array atol=1e-6
    end

    @testset "Limit: Single Component Square-Well RPA" begin
        sigma_sc = [1.0]
        phi_sc_val = 0.2
        phi_sc_vec = [phi_sc_val]
        
        temp_sc = 1.5
        lambda_sc = 1.5
        
        temp_matrix_sc = fill(temp_sc, (1,1))
        lambda_matrix_sc = fill(lambda_sc, (1,1))

        for k_val in k_test_array_sw
            S_matrix_rpa_sc = S_RPA_mixture_SquareWell(sigma_sc, phi_sc_vec, k_val, 
                                                       temp_matrix_sc, lambda_matrix_sc)
            @test size(S_matrix_rpa_sc) == (1,1)
            @test eltype(S_matrix_rpa_sc) <: Real
            s11_rpa_mixture = S_matrix_rpa_sc[1,1]

            s_mono_sw_rpa = S_SW_RPA(phi_sc_val, temp_sc, lambda_sc, k_val) 
            @test s11_rpa_mixture ≈ s_mono_sw_rpa atol=1e-7
        end
    end

    @testset "Binary Mixture Properties (Square-Well RPA)" begin
        # Estas variables deberían ser accesibles desde el @testset padre
        for k_val in k_test_array_sw
            S_matrix_bin_sw = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, k_val,
                                                       temperature_matrix_example, lambda_matrix_example)
            
            @test size(S_matrix_bin_sw) == (2,2)
            @test eltype(S_matrix_bin_sw) <: Real
            @test S_matrix_bin_sw[1,2] ≈ S_matrix_bin_sw[2,1] atol=1e-7
            @test S_matrix_bin_sw[1,1] > 0.0
            @test S_matrix_bin_sw[2,2] > 0.0
        end
    end
    
    @testset "Factory and Vectorized Versions (Square-Well RPA)" begin
        # Estas variables deberían ser accesibles
        s_k_func_mixture_sw = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, 
                                                       temperature_matrix_example, lambda_matrix_example)
        S_matrix_direct = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, k_test_val_sw,
                                                   temperature_matrix_example, lambda_matrix_example)
        S_matrix_factory = s_k_func_mixture_sw(k_test_val_sw)
        @test S_matrix_factory ≈ S_matrix_direct atol=1e-9

        S_matrices_vectorized = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, 
                                                         temperature_matrix_example, lambda_matrix_example, 
                                                         k_test_array_sw)
        @test length(S_matrices_vectorized) == length(k_test_array_sw)
        for i_k in 1:length(k_test_array_sw)
            S_direct_i = S_RPA_mixture_SquareWell(sigma_vec_bin, phi_vec_bin, k_test_array_sw[i_k],
                                                  temperature_matrix_example, lambda_matrix_example)
            @test S_matrices_vectorized[i_k] ≈ S_direct_i atol=1e-9
        end
    end

end
