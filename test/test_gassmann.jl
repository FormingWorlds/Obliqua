using Test
using Obliqua

# Undrained (Gassmann) modulus of layers below the porosity threshold (common.gassmann).
#
#   1. gassmann: K_u = K_s without melt, = K_d for alpha = 0, between them otherwise.
#   2. Crossing the porosity threshold changes k2 (fluid1d_relax) much less with K_u than with the grain modulus K_s.

const GS_CMN = Obliqua.fluid1d_relax.common
const GS_F = Obliqua.fluid1d_relax

@testset "gassmann: undrained modulus below the porosity threshold" begin
    @testset "gassmann" begin
        Ks, Kd, Kl = ComplexF64(1e11), ComplexF64(4e10, 4e9), 2e10
        @test GS_CMN.gassmann(Kd, Ks, Kl, 0.0) === Ks
        @test GS_CMN.gassmann(Ks, Ks, Kl, 0.1) == Ks       # alpha = 0
        Ku = GS_CMN.gassmann(Kd, Ks, Kl, 0.02)
        @test real(Kd) < real(Ku) < real(Ks)
        # closed form check
        α = 1 - Kd / Ks; S = 0.02 / Kl + (α - 0.02) / Ks
        @test Ku ≈ Kd + α^2 / S
    end

    @testset "porosity threshold: Gassmann keeps the bulk stiffness continuous" begin
        # one porous layer just above or just below the threshold in an elastic column; above it the low permeability
        # makes the layer effectively undrained, below it the layer is non-porous with K_u
        thr = 0.03
        function kthr(ϕ1; use_gassmann = true)
            N = 6
            r = collect(range(3000e3, 6000e3, length = N + 1)); ρ = fill(4000.0, N); η = fill(1e21, N)
            ϕ = zeros(N); ϕ[4] = ϕ1; K = zeros(N); K[4] = 1e-14 * (ϕ1 / thr)^3
            μ = fill(ComplexF64(5e10), N); κ = fill(ComplexF64(1e11), N); κd = copy(κ)
            κd[4] = ComplexF64(6e10); κl = fill(2e10, N)
            αB = 1 .- κd ./ κ
            if !use_gassmann && ϕ1 < thr
                κd[4] = κ[4]; αB[4] = 0                    # the former treatment: K_d = K_s below the threshold
                ϕ1 = 0.0; ϕ[4] = 0.0
            end
            cs = GS_F.column_spec(1e-5, r, ρ, η, μ, κ, κd, αB, κl, 1.0, ϕ, K, 4π / 3 * r[1]^3 * 1e4, 1e4,
                                  0.1 + 0im, 5e11 + 0im; σ_R = 0.0, dr_min = 1000, dr_max = 3000, porosity_thresh = thr)
            return GS_F.solve_column(cs).k
        end
        k_above = kthr(thr * 1.001)
        jump_gassmann = abs(kthr(thr * 0.999) - k_above)
        jump_grain = abs(kthr(thr * 0.999; use_gassmann = false) - k_above)
        @test jump_gassmann < 0.2 * jump_grain
        @test jump_gassmann / abs(k_above) < 1e-4
    end
end
