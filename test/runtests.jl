# test/runtests.jl

using AnalyticalStructureFactors
using Test

@testset "AnalyticalStructureFactors.jl Tests" begin
    @testset "Monodisperse Systems (PY, VW, WCA, SW, Yukawa RPA)" begin
        include("test_monodisperse.jl")
    end

    @testset "SALR Models (RPA & Exact MSA)" begin
        include(joinpath("models", "test_SALR.jl"))
        include(joinpath("models", "test_SALR_msa.jl"))
    end

    @testset "Exact Single-Yukawa MSA Model" begin
        include(joinpath("models", "test_yukawa_msa.jl"))
    end

    @testset "Fractal Aggregates Structure Factor" begin
        include(joinpath("models", "test_fractal.jl"))
    end

    @testset "Particle Form Factors P(q) and Decoupling Factor beta(q)" begin
        include(joinpath("models", "test_formfactors.jl"))
    end

    @testset "Colloidal Electrostatics & DLVO Parameters" begin
        include("test_colloidal_potentials.jl")
    end

    @testset "3D Spherical Fourier Transform Tools" begin
        include("test_fourier_tools.jl")
    end

    @testset "Multicomponent Mixture Systems" begin
        include(joinpath("models", "mixtures", "test_hardsphere_mixture.jl"))
        include(joinpath("models", "mixtures", "test_squarewell_rpa_mixture.jl"))
        include(joinpath("models", "mixtures", "test_yukawa_rpa_mixture.jl"))
    end

    @testset "Smoke Tests for All Models" begin
        include("test_smoke_all_models.jl")
    end
end
