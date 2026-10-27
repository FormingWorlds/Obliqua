using Test
using Obliqua

# Melt-fraction laws of the porous skeleton (src/interior.jl): drained (compaction) bulk modulus and the Biot
# coefficient alpha = 1 - K_d/K_s that the driver derives from it (Hay et al. 2025, Eqs. 48-49).
#   K_d = K_s for phi < phi_0 = 0.0769;  K_s (phi_0/phi)^b up to phi_crit;  0 at and above phi_crit (no skeleton).
#   phi_crit = orbit.obliqua.solid.dbulk_phi_crit, default 0.7715 (Rumpf-Gupte to Stokes permeability crossover);
#   0.3 (disaggregation, Hay et al. 2025) is supported as a configured value.

"""Minimal configuration dictionary for the interior helper functions."""
ir_cfg(; b=0.5, phicrit=nothing) = begin
    solid = Dict{String,Any}("dbulk_power" => b, "porosity_thresh" => 0.03)
    phicrit === nothing || (solid["dbulk_phi_crit"] = phicrit)
    Dict{String,Any}("orbit" => Dict{String,Any}("obliqua" => Dict{String,Any}("solid" => solid)),
                     "interior_energetics" => Dict{String,Any}("grain_size" => 1e-3))
end

@testset "interior: drained bulk modulus and skeleton collapse" begin
    Ks = 1.0e11
    phi = [0.0, 0.05, 0.0769 - 1e-9, 0.1, 0.2, 0.3, 0.45, 0.7715 - 1e-9, 0.7715, 0.9, 1.0]
    bulk = fill(Ks, length(phi))

    @testset "three branches with the default threshold (0.7715)" begin
        Kd = Obliqua.interior.get_drained_bulk(bulk, phi, ir_cfg())
        α = 1 .- Kd ./ Ks
        @test all(Kd[1:3] .≈ Ks)                                   # coherent skeleton: K_d = K_s, alpha = 0
        @test all(abs.(α[1:3]) .< 1e-12)
        @test Kd[4] ≈ Ks * (0.0769/0.1)^0.5 rtol=1e-12             # power-law branch (b = 0.5)
        @test Kd[7] ≈ Ks * (0.0769/0.45)^0.5 rtol=1e-12            # mush at 0.45 keeps ~0.41 K_s
        @test Kd[8] ≈ Ks * (0.0769/0.7715)^0.5 rtol=1e-8           # just below the threshold: ~0.32 K_s, not zero
        @test all(iszero, Kd[9:end])                               # at and above 0.7715: no skeleton, alpha = 1
        @test all(α[9:end] .≈ 1)
        @test all(diff(Kd) .<= 0)                                   # K_d non-increasing with melt fraction
        @test all(0 .<= α .<= 1)                                    # Biot coefficient bounded
        # guard: the disaggregation threshold (0.3) removes the skeleton of the same mush, so the two choices are
        # distinguished at phi = 0.3 and 0.45 by a factor > 2 in K_s, far above the tolerance
        Kd_dis = Obliqua.interior.get_drained_bulk(bulk, phi, ir_cfg(phicrit=0.3))
        @test Kd_dis[5] ≈ Kd[5] rtol=1e-12                          # below 0.3 the two thresholds agree
        @test all(iszero, Kd_dis[6:end])
        @test Kd[6] > 0.5 * Ks && Kd[7] > 0.4 * Ks
    end

    @testset "configurable threshold and exponent" begin
        Kd = Obliqua.interior.get_drained_bulk(bulk, phi, ir_cfg(b=1.5, phicrit=0.5))
        @test Kd[7] ≈ Ks * (0.0769/0.45)^1.5 rtol=1e-12             # 0.45 < 0.5: still on the power-law branch
        @test iszero(Kd[8])                                          # 0.77 > 0.5: no skeleton
        @test Kd[4] < Ks * (0.0769/0.1)^0.5                          # larger b: faster loss of compaction stiffness
    end

    @testset "error contract: threshold must lie above phi_0" begin
        @test_throws ArgumentError Obliqua.interior.get_drained_bulk(bulk, phi, ir_cfg(phicrit=0.05))
        @test_throws ArgumentError Obliqua.interior.get_drained_bulk(bulk, phi, ir_cfg(phicrit=1.5))
    end
end
