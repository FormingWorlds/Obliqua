using Test
using Obliqua

# Tests of the helper module (src/helper.jl): argument checks and the tidal work of the energy closure. The column
# diagnostics (face_residuals, boundary_layer_resolution) are tested on solved columns in test_fluid1d_relax.jl.

const HL = Obliqua.fluid1d_relax.helper

@testset "helper" begin

    @testset "argument checks" begin
        @test HL.check_arg(true, "never thrown")
        @test_throws ArgumentError HL.check_arg(false, "f: x must be > 0")
        @test HL.check_option("f", "slip", :noslip, (:freeslip, :noslip))
        err = try HL.check_option("f", "slip", :sticky, (:freeslip, :noslip)) catch e; e end
        @test err isa ArgumentError && err.msg == "f: slip must be one of :freeslip, :noslip, got :sticky"
        @test_throws ArgumentError HL.check_option("f", "n", 3, 1:2)
    end

    @testset "tidal work" begin
        # (2n+1) r ω / (8π G) (-Im k): degree 2, R = 6371 km, ω = 1e-4 rad/s, Im k = -1e-3
        G = Obliqua.constants.G
        P = HL.tidal_work(0.3 - 1e-3im, 2, 6.371e6, 1e-4)
        @test P ≈ 5 * 6.371e6 * 1e-4 / (8π * G) * 1e-3 rtol = 1e-14
        # the degree enters as 2n + 1 (n = 3: 7/5 of n = 2), and the sign follows -Im k: a passive response gives P > 0
        @test HL.tidal_work(0.3 - 1e-3im, 3, 6.371e6, 1e-4) ≈ 1.4P rtol = 1e-14
        @test HL.tidal_work(0.3 + 1e-3im, 2, 6.371e6, 1e-4) ≈ -P rtol = 1e-14
        # the real part does no work
        @test HL.tidal_work(0.9 + 0im, 2, 6.371e6, 1e-4) == 0
    end
end
