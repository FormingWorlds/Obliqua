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
end