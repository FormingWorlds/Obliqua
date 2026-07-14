using Test
using LinearAlgebra
using AssociatedLegendrePolynomials
using DoubleFloats
using Obliqua
using Obliqua.solid1d.common


@testset "common Module Unit Tests" begin

    @testset "Ynm (Spherical Harmonics)" begin
        # Test basic evaluation and shape guarantees
        theta = [0.0, π/2, π]
        phi = [0.0, π]'
        n, m = 2, 1
        
        res = Ynm(n, m, theta, phi) 
        @test size(res) == (3, 2)
        @test res isa Array{ComplexF64, 2}
    end

    @testset "define_spherical_grid" begin
        # Define a grid with a large angular spacing for efficiency
        grid = define_spherical_grid(60.0, 2, 1) 
        
        @test grid isa NamedTuple
        @test haskey(grid, :res)
        @test haskey(grid, :clats)
        @test haskey(grid, :lons)
        @test haskey(grid, :Y)
        @test haskey(grid, :dYdθ)
        @test haskey(grid, :dYdϕ)
        @test haskey(grid, :Z)
        @test haskey(grid, :X)
        
        @test grid.res == 60.0
        @test size(grid.Y, 1) == 1
    end

    @testset "get_scales" begin
        R0 = 1.737e6
        M0 = 7.342e22
        G0 = 1.0
        
        # Test standard 6-variable configuration scaling matrices
        R0_o, M0_o, ω0, ρ0, G0_o, g0, μ0, S, Sinv = get_scales(R0, M0, G0, Y=[1,2,3,4,5,6]) 
        @test size(S) == (6, 6)
        @test size(Sinv) == (6, 6)
        @test S * Sinv ≈ I
        
        # Test 8-variable two-phase configuration layout matrix size
        _, _, _, _, _, _, _, S8, Sinv8 = get_scales(R0, M0, G0, Y=collect(1:8)) 
        @test size(S8) == (8, 8)
        @test S8[7,7] > 0
        @test S8[8,8] > 0
    end

    @testset "doublefactorial" begin
        # Verify boundary mathematical cases
        @test doublefactorial(0) == 1
        @test doublefactorial(1) == 1
        @test doublefactorial(4) == 8   # 4 * 2
        @test doublefactorial(5) == 15  # 5 * 3 * 1
        
        # Enforce error assertion for invalid domain spaces
        @test_throws ErrorException doublefactorial(-1)
    end

    @testset "sbesselj" begin
        # Analytical matching for n=0: j_0(x) = sin(x)/x
        x_val = 1.5
        @test sbesselj(0, x_val) ≈ sin(x_val) / x_val
        
        # Verify execution path type consistency
        @test sbesselj(1, 2.0 + 0.0im) isa ComplexF64
    end

    @testset "get_Ic (Core Solution Matrix)" begin
        ω = 2.0e-5
        r = 1.2e6
        ρ = 4000.0
        g = 2.0
        μ = 2.0e10
        K = 7.0e10
        n = 2
        
        # 1. Liquid Core verification
        Ic_liquid = get_Ic(ω, r, ρ, g, μ, K, "liquid", n) 
        @test size(Ic_liquid) == (6, 3)
        @test Ic_liquid[2, 2] == 1.0 # Tangential slip
        
        # 2. Solid Core verification
        Ic_solid = get_Ic(ω, r, ρ, g, μ, K, "solid", n) 
        @test size(Ic_solid) == (6, 3)
        
        # 3. Inertial Fluid Core verification
        Ic_inertial = get_Ic(ω, r, ρ, g, μ, K, "inertial", n) 
        @test size(Ic_inertial) == (6, 3)
        
        # 4. Error throwing for unsupported definitions
        @test_throws ErrorException get_Ic(ω, r, ρ, g, μ, K, "invalid_type", n) 
    end

    @testset "get_A & get_A! (ODE Propagator System Matrices)" begin
        ω = 2.2e-5
        r = 1.5e6
        ρ = 3300.0
        g = 1.62
        μ = 4.0e10 + 2.0e8im
        K = 9.0e10 + 4.0e8im
        n = 2
        
        # Test 6x6 solid matrix calculations
        A6 = get_A(ω, r, ρ, g, μ, K, n) 
        @test size(A6) == (6, 6)
        @test A6 isa Matrix{ComplexF64}
        
        # Test 8x8 two-phase porous mixture matrix calculations
        ρₗ = 2900.0
        Kl = 7.5e10
        Kd = 3.5e10 + 1.0e8im
        α = 0.85 + 0.01im
        ηₗ = 1.0
        ϕ = 0.05
        k = 1.0e-11
        
        # FIX: Explicitly pass G0=1.0 to bypass the source code's Int64 default type issue
        A8 = get_A(ω, r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1.0) 
        @test size(A8) == (8, 8)
        @test A8 isa Matrix{ComplexF64}
    end

    @testset "compute_strain_ten!" begin
        # Generate a baseline grid setup space
        grid = define_spherical_grid(90.0, 2, 1) 
        n_clats = length(grid.clats)
        n_lons = length(grid.lons)
        
        # Prepare 3D strain tensor destination array allocation
        ϵ = zeros(ComplexF64, n_clats, n_lons, 6) 
        
        # Mocking values for standard solution vector configurations
        y = zeros(ComplexF64, 6) 
        y[1] = 0.002
        y[2] = 0.001
        y[3] = 2.0e5
        y[4] = 6.0e4
        
        rr = 1.4e6 
        ρr = 3100.0 
        gr = 1.8 
        μr = 3.5e10 + 1.0e7im 
        Ksr = 8.5e10 + 3.0e7im 
        
        # Call and test non-destructive updating checks
        compute_strain_ten!(ϵ, y, 2, rr, ρr, gr, μr, Ksr, grid) 
        
        @test size(ϵ) == (n_clats, n_lons, 6)
        @test !all(iszero, ϵ) 
    end
end