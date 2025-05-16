# test/runtests.jl

using AnalyticalStructureFactors # Importa tu paquete
using Test                   # Importa el módulo de pruebas de Julia

@testset "AnalyticalStructureFactors.jl Tests" begin
    # Aquí puedes tener pruebas generales o incluir otros archivos de prueba.

    @testset "Monodisperse Systems" begin
        # Esta sección incluirá pruebas para los modelos monodispersos.
        # Las pruebas específicas estarán en test_monodisperse.jl
        include("test_monodisperse.jl")
    end

    @testset "SALR Model" begin
        include(joinpath("models", "test_SALR.jl"))
    end
    
    @testset "Smoke Tests for All Models" begin
        include("test_smoke_all_models.jl")
    end
    
    @testset "Mixture Systems" begin
        # Esta sección incluirá pruebas para los modelos de mezclas.
        # Puedes crear un archivo test_mixtures.jl para estas.
        # Por ahora, lo dejaremos vacío o con una prueba simple si no tienes código de mezclas aún.
        #@test true # Una prueba trivial para que el testset no esté vacío
        include("models/mixtures/test_hardsphere_mixture.jl") # Descomenta cuando tengas este archivo
        include("models/mixtures/test_squarewell_rpa_mixture.jl") # Descomenta cuando tengas este archivo
        include("models/mixtures/test_yukawa_rpa_mixture.jl") # Descomenta cuando tengas este archivo
    end

    # Puedes añadir más @testset para otras categorías de pruebas.
end
