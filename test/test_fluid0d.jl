using Test
using LinearAlgebra
using Obliqua
using Obliqua.fluid0d


@testset "fluid0d Module Tests" begin

    @testset "compute_fluid_lovenumbers" begin
        # Setup mock planetary/magma ocean parameters
        omega = 2.0e-5      # Forcing frequency (rad/s)
        R = 1.737e6         # Radius (meters)
        H_magma = 50.0e3    # Magma ocean thickness (meters)
        g = 1.62            # Gravity (m/s^2)
        ρ_ratio = 0.15      # Density contrast
        n = 2               # Degree
        σ_R = 1.0e-4        # Rayleigh drag coefficient

        # Execute target function
        k2_T, k2_L = compute_fluid_lovenumbers(omega, R, H_magma, g, ρ_ratio, n, σ_R)

        # 1. Verify exact output types conform to specifications
        @test k2_T isa ComplexF64
        @test k2_L isa ComplexF64

        # 2. Analytical validation
        μ_n  = n * (n + 1)
        ξ_n = 3.0 / (2.0 * n + 1.0) * ρ_ratio
        σ_n = sqrt(μ_n * g * H_magma / R^2)
        σ_T = omega - im * σ_R
        k2_T_expected = -ξ_n * σ_n^2 / (omega * σ_T - σ_n^2)

        @test k2_T ≈ k2_T_expected
        
        # 3. Verify currently hardcoded Load Love number
        @test k2_L == 0.0
    end

    @testset "mean_rho (Volume-Weighted Density)" begin
        # Test Case 1: Uniform profile consistency
        # If all shells share an identical density, the aggregate mean must equal that density
        r_uniform = [1000.0, 2000.0, 3000.0, 4000.0] # 3 shell boundaries
        ρ_uniform = [3300.0, 3300.0, 3300.0]        # 3 shell layers
        
        @test mean_rho(ρ_uniform, r_uniform) ≈ 3300.0

        # Test Case 2: Discrete multi-layer analytical validation
        r_two = [1.0, 2.0, 3.0]      # 2 shells (radii boundaries 1->2 and 2->3)
        ρ_two = [2000.0, 4000.0]     # Distinct density layers
        
        # Calculate expected values manually
        V1 = 4/3 * π * (2.0^3 - 1.0^3)
        V2 = 4/3 * π * (3.0^3 - 2.0^3)
        M_manual = (V1 * 2000.0) + (V2 * 4000.0)
        mean_rho_manual = M_manual / (V1 + V2)

        @test mean_rho(ρ_two, r_two) ≈ mean_rho_manual

        # Test Case 3: Dimension Mismatch Assertion
        # Providing mismatched matrix profiles should invoke a broadcast DimensionMismatch error
        r_mismatch = [1.0, 2.0, 3.0]           # yields 2 shell spaces
        ρ_mismatch = [1000.0, 2000.0, 3000.0]  # provides 3 density configurations
        
        @test_throws DimensionMismatch mean_rho(ρ_mismatch, r_mismatch)
    end

    @testset "enclosed_mean_rho (self-attraction density)" begin
        # Earth-like: core of 12274 kg/m^3 to 3355 km under a fully molten
        # mantle of 3247 kg/m^3 to 6338 km (values from a saved PROTEUS run).
        rc, R = 3.355e6, 6.338e6
        ρc, ρm = 12274.0, 3247.0
        M = 4/3 * π * (rc^3 * ρc + (R^3 - rc^3) * ρm)
        ρ_bar = enclosed_mean_rho(M, R)
        @test ρ_bar ≈ M / (4/3 * π * R^3)
        # Consistent with surface gravity: rho_bar = 3 g / (4 pi G R)
        G = 6.6743e-11
        @test ρ_bar ≈ 3 * (G * M / R^2) / (4 * π * G * R)
        # Discriminating: the fluid/enclosed ratio (~0.71) differs from the
        # fluid/core ratio (~0.26) by far more than rounding
        @test isapprox(ρm / ρ_bar, 0.708; atol=0.005)
        @test abs(ρm / ρ_bar - ρm / ρc) > 0.3
        # Bounded by the densest and lightest enclosed material
        @test ρm < ρ_bar < ρc
        # Uniform body: enclosed mean equals the material density
        @test enclosed_mean_rho(4/3 * π * R^3 * 3300.0, R) ≈ 3300.0
        # Error contract: a non-positive radius is rejected
        @test_throws DomainError enclosed_mean_rho(1.0e24, 0.0)
        @test_throws DomainError enclosed_mean_rho(1.0e24, -1.0)
    end
end

@testset "run_fluid0d surface gravity" begin
    # Thin (87 km) magma ocean on an Earth-sized solid interior. The shallow-water
    # mode sigma_n = sqrt(n(n+1) g H)/R must use the gravity at the ocean surface
    # from the full enclosed mass (g_top), not only the ocean shells' own mass.
    G_si  = Float64(Obliqua.constants.G)
    r     = Obliqua.prec.(collect(range(6.25e6, 6.337e6, length=4)))
    ρ     = Obliqua.prec.(fill(3000.0, 3))
    R     = Float64(maximum(r))
    H     = Float64(r[end] - r[1])
    ω     = 4.5e-4
    σ_R   = 1.0e-4
    ρ_rat = Obliqua.prec(0.7)
    g_top = 10.28

    kT, kL = Obliqua.run_fluid0d(ω, ρ, r, ρ_rat; n=2, sigma_R=σ_R, g_top=g_top)
    kT_ref, kL_ref = compute_fluid_lovenumbers(ω, R, H, g_top, Float64(ρ_rat), 2, σ_R)
    @test kT ≈ kT_ref rtol=1e-12
    @test kL ≈ kL_ref rtol=1e-12

    # Legacy fallback (g from the ocean shells only) is ~0.4 m/s^2 here, ~25x too
    # small; it must stay available but give a clearly different response.
    g_shell = G_si * sum(4/3 * π * (Float64.(r[2:end]).^3 .- Float64.(r[1:end-1]).^3) .* 3000.0) / R^2
    @test g_shell < g_top / 10
    kT_legacy, _ = Obliqua.run_fluid0d(ω, ρ, r, ρ_rat; n=2, sigma_R=σ_R)
    kT_legacy_ref, _ = compute_fluid_lovenumbers(ω, R, H, g_shell, Float64(ρ_rat), 2, σ_R)
    @test kT_legacy ≈ kT_legacy_ref rtol=1e-10
    @test abs(kT_legacy - kT) > 0.5 * abs(kT)

    # Passivity for forcing of either sign: sign(ω) * Im(k_T) <= 0.
    kT_neg, _ = Obliqua.run_fluid0d(-ω, ρ, r, ρ_rat; n=2, sigma_R=σ_R, g_top=g_top)
    @test imag(kT) <= 0
    @test imag(kT_neg) >= 0
end
