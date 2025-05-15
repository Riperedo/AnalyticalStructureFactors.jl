# test/models/mixtures/test_hardsphere_mixture.jl

# Asumimos que este archivo será incluido por test/runtests.jl,
# donde AnalyticalStructureFactors, Test, y LinearAlgebra ya están importados.

@testset "Hard Sphere Mixture (Baxter Solution)" begin

    # --- Parámetros comunes para pruebas de mezclas ---
    k_test_val = 1.5
    k_test_array = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    # --- Pruebas para Qk_mixture ---
    @testset "Qk_mixture Function" begin
        @testset "Binary Mixture (n=2)" begin
            σ_bin = [1.0, 0.8]
            ϕ_bin = [0.1, 0.05]
            ρ_bin = phi_to_rho_mixture(ϕ_bin, σ_bin) # Usa tu función de utils.jl

            Q_matrix = Qk_mixture(σ_bin, ρ_bin, k_test_val)
            @test size(Q_matrix) == (2, 2)
            @test eltype(Q_matrix) <: Complex

            # Prueba para k=0 (debería usar la rama k≈0 en Qk_mixture)
            Q_matrix_k0 = Qk_mixture(σ_bin, ρ_bin, 0.0)
            @test size(Q_matrix_k0) == (2, 2)
            @test eltype(Q_matrix_k0) <: Complex 
            # Aquí se podrían añadir pruebas más específicas para los valores de Q(k=0)
            # si se conocen las expresiones analíticas para ∫Q_ij(r)dr.
        end

        @testset "Single Component (n=1, as mixture)" begin
            σ_mono_vec = [1.0]
            ϕ_mono_vec = [0.2]
            ρ_mono_vec = phi_to_rho_mixture(ϕ_mono_vec, σ_mono_vec)

            Q_matrix_mono = Qk_mixture(σ_mono_vec, ρ_mono_vec, k_test_val)
            @test size(Q_matrix_mono) == (1, 1)
            @test eltype(Q_matrix_mono) <: Complex
        end
    end

    # --- Pruebas para IS_HS_Baxter_mixture y S_HS_Baxter_mixture ---
    @testset "S_HS_Baxter_mixture & IS_HS_Baxter_mixture" begin

        @testset "Single Component Limit (n=1)" begin
            σ_sc = [1.0]
            ϕ_sc_val = 0.2
            ϕ_sc_vec = [ϕ_sc_val]
            ρ_sc_vec = phi_to_rho_mixture(ϕ_sc_vec, σ_sc)

            for k_val in k_test_array
                # S(k) de la mezcla de un componente
                S_matrix_sc = S_HS_Baxter_mixture(σ_sc, ρ_sc_vec, k_val)
                @test size(S_matrix_sc) == (1,1)
                @test eltype(S_matrix_sc) <: Real
                s11_mixture = S_matrix_sc[1,1]

                # S(k) del sistema monodisperso PY (de hardsphere.jl)
                s_mono_py = S_HS_PY(ϕ_sc_val, k_val) # Asume que S_HS_PY toma ϕ y k

                @test s11_mixture ≈ s_mono_py atol=1e-7

                # Consistencia IS y S
                IS_matrix_sc = IS_HS_Baxter_mixture(σ_sc, ρ_sc_vec, k_val)
                @test size(IS_matrix_sc) == (1,1)
                @test eltype(IS_matrix_sc) <: Complex # Q'*Q puede ser complejo antes de inv
                @test S_matrix_sc ≈ real.(inv(IS_matrix_sc)) atol=1e-7
            end
        end

        @testset "Binary Mixture Properties (n=2)" begin
            σ_bin = [1.0, 0.5] # Diámetros diferentes
            ϕ_bin = [0.1, 0.15] # Fracciones de volumen
            ρ_bin = phi_to_rho_mixture(ϕ_bin, σ_bin)

            for k_val in k_test_array
                S_matrix_bin = S_HS_Baxter_mixture(σ_bin, ρ_bin, k_val)
                
                @test size(S_matrix_bin) == (2,2)
                @test eltype(S_matrix_bin) <: Real
                
                # Simetría: S_ij(k) = S_ji(k)
                @test S_matrix_bin[1,2] ≈ S_matrix_bin[2,1] atol=1e-7
                
                # Positividad de los elementos diagonales
                @test S_matrix_bin[1,1] > 0.0
                @test S_matrix_bin[2,2] > 0.0
                
                # Condición de estabilidad termodinámica (más avanzada, opcional por ahora):
                # det(S_matrix_bin) > 0 o S11*S22 - S12^2 > 0
                # @test det(S_matrix_bin) > 0.0 (si S_matrix_bin es definida positiva)

                # Consistencia IS y S
                IS_matrix_bin = IS_HS_Baxter_mixture(σ_bin, ρ_bin, k_val)
                @test size(IS_matrix_bin) == (2,2)
                @test S_matrix_bin ≈ real.(inv(IS_matrix_bin)) atol=1e-7
            end
        end

        @testset "Ideal Gas Limit for Mixtures (ϕ_total -> 0)" begin
            σ_gas = [1.0, 0.8]
            # Fracciones de volumen muy pequeñas
            ϕ_gas = [1e-8, 0.5e-8] 
            ρ_gas = phi_to_rho_mixture(ϕ_gas, σ_gas)

            for k_val in k_test_array
                S_matrix_gas = S_HS_Baxter_mixture(σ_gas, ρ_gas, k_val)
                
                @test S_matrix_gas[1,1] ≈ 1.0 atol=1e-5 # Tolerancia más laxa por ser límite
                @test S_matrix_gas[2,2] ≈ 1.0 atol=1e-5
                @test S_matrix_gas[1,2] ≈ 0.0 atol=1e-5
                @test S_matrix_gas[2,1] ≈ 0.0 atol=1e-5
            end
        end

        @testset "Factory and Vectorized Versions" begin
            σ_fac = [1.0, 0.7]
            ϕ_fac = [0.1, 0.1]
            ρ_fac = phi_to_rho_mixture(ϕ_fac, σ_fac)

            # Prueba de la función fábrica
            s_k_func_mixture = S_HS_Baxter_mixture(σ_fac, ρ_fac)
            S_matrix_direct = S_HS_Baxter_mixture(σ_fac, ρ_fac, k_test_val)
            S_matrix_factory = s_k_func_mixture(k_test_val)
            @test S_matrix_factory ≈ S_matrix_direct atol=1e-9

            # Prueba de la versión vectorizada
            S_matrices_vectorized = S_HS_Baxter_mixture(σ_fac, ρ_fac, k_test_array)
            @test length(S_matrices_vectorized) == length(k_test_array)
            for i in 1:length(k_test_array)
                S_direct_i = S_HS_Baxter_mixture(σ_fac, ρ_fac, k_test_array[i])
                @test S_matrices_vectorized[i] ≈ S_direct_i atol=1e-9
            end
        end
    end
end
