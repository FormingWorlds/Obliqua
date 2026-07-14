using Test
using Obliqua
using Obliqua.constants


@testset "constants Module Tests" begin

    @testset "Exported Precision Types" begin
        # Verify current type configuration (Float64)
        @test prec === Float64
        @test precc === Complex{Float64}
        
        # Architecture check: Ensure precc is always the complex counterpart of prec,
        # which guarantees type stability if you switch to BigFloat or Double64 later.
        @test precc === Complex{prec}
    end

    @testset "Physical Constants" begin
        # Verify values match expectations
        @test AU ≈ 1.495978707e11
        @test G ≈ 6.6743e-11
        @test M_Earth ≈ 5.9724e24
        
        # Verify type compliance with the module's defined precision
        @test AU isa prec
        @test G isa prec
        @test M_Earth isa prec
    end

    @testset "Directory Paths" begin
        # Verify paths are bound to strings
        @test ROOT_DIR isa String
        @test RES_DIR isa String
        @test OUT_DIR isa String

        # Ensure the script correctly resolved absolute paths
        @test isabspath(ROOT_DIR)
        @test isabspath(RES_DIR)
        @test isabspath(OUT_DIR)

        # Verify that subdirectories are nested under ROOT_DIR as expected
        @test occursin(ROOT_DIR, RES_DIR)
        @test occursin(ROOT_DIR, OUT_DIR)
        
        # Verify names match the target folders
        @test endswith(RES_DIR, "res/") || occursin("res", RES_DIR)
        @test endswith(OUT_DIR, "out/") || occursin("out", OUT_DIR)
    end

    @testset "Grid Resolution" begin
        @test res == 60.0
        @test res isa Float64
    end
end