using Test
using LinearAlgebra
using Obliqua
using Obliqua.solid1d_mush_relax


@testset "solid1d_mush_relax Module Tests" begin

    # Define a realistic physical configuration for a planetary body (SI units)
    # This ensures N > 1 during grid resampling, preventing empty vectors.
    radius = [1000000.0, 1500000.0, 2000000.0] # Layer boundaries in meters
    rho    = [3000.0, 3500.0]                 # Density (kg/m^3)
    visc   = [1e21, 1e22]                     # Viscosity (Pa*s)
    shear  = complex.([1e10, 2e10])           # Complex shear modulus (Pa)
    bulk_s = complex.([5e10, 6e10])           # Solid bulk modulus (Pa)
    bulk_l = [4e10, 4.5e10]                   # Liquid bulk modulus (Pa)
    bulk_d = complex.([4.8e10, 5.5e10])       # Drained bulk modulus (Pa)
    alpha  = complex.([0.8, 0.9])             # Biot coefficient
    visc_l = [1.0, 2.0]                       # Liquid viscosity (Pa*s)
    phi    = [0.1, 0.05]                      # Porosity
    k      = [1e-12, 0.0]                     # Permeability (m^2) - Layer 1 is mushy, Layer 2 is solid
    m_core = 1e23                             # Core mass (kg)
    dr_min = 10000                            # Minimum spacing (10 km)
    dr_max = 50000                            # Maximum spacing (50 km)
    
    # Global physics parameters
    ω = 1e-4                                  # Angular frequency
    n = 2                                     # Tidal degree
    ρ_core = 7000.0                           # Core density
    μ_core = 0.0                              # Liquid core shear modulus
    κ_core = 1e11                             # Core bulk modulus
    M_tot  = 6e24                             # Approximate total mass target

    @testset "1. Grid Generation & Gravity Calculations" begin
        # Test profile resampling
        r_new_b, ρ_new, η_new, μ_new, κs_new, κl_new, κd_new, α_new, ηl_new, φ_new, k_new, g_new, M_tot_calc = 
            solid1d_mush_relax.resample_profiles(
                radius, rho, visc, shear, bulk_s, bulk_l, bulk_d, alpha, visc_l, phi, k, m_core, dr_min, dr_max
            )

        # Basic size and dimension validations
        N = length(r_new_b)
        @test N > 1
        @test length(ρ_new)  == N - 1
        @test length(μ_new)  == N - 1
        @test length(g_new)  == N - 1
        @test length(k_new)  == N - 1

        # Check boundary integrity
        @test r_new_b[1]   ≈ radius[1]
        @test r_new_b[end] ≈ radius[end]

        # Test the underlying gravity compiler routine directly
        g_test, M_enc_test = solid1d_mush_relax.get_g(r_new_b, ρ_new, m_core)
        @test length(g_test) == N - 1
        @test M_enc_test > m_core
    end

    @testset "2. Core Boundary Condition Drivers" begin
        # Generate non-dimensional profiles similar to solve_radial_system execution
        Nr = 5
        rs = collect(range(0.5, 1.0, length=Nr))
        ρs = fill(1.0, Nr)
        gs = fill(1.0, Nr)
        μs = fill(complex(1.0), Nr)
        Ks = fill(complex(2.0), Nr)
        
        # Additional mushy specific arrays
        ρₗs = fill(0.9, Nr)
        Kls = fill(1.5, Nr)
        Kds = fill(complex(1.2), Nr)
        αs  = fill(complex(0.8), Nr)
        ηₗs = fill(1.0, Nr)
        ϕs  = fill(0.1, Nr)
        ks  = fill(1e-4, Nr)

        # Allocate state vector integration arrays
        R_solid = [zeros(precc, 8, 8) for _ in 1:Nr]
        R_mush  = [zeros(precc, 8, 8) for _ in 1:Nr]
        
        idx_6 = [1, 2, 3, 5, 6, 7]
        R6_view = [view(R_solid[i], idx_6, idx_6) for i in 1:Nr]
        R8_view = [view(R_mush[i], 1:8, 1:8) for i in 1:Nr]

        # Test Solid Core Boundary initialization (6x6 context)
        C1l, D2l = solid1d_mush_relax.core_boundary(
            R6_view, (1, 2), rs, ρs, gs, μs, Ks, ω, ρ_core, μ_core, κ_core, "liquid", n; G0=1.0, Y=[1,2,4,5,3,6]
        )
        @test size(C1l) == (3, 6)
        @test size(D2l) == (3, 6)
        @test any(R_solid[1] .!= 0.0)

        # Test Mushy Core Boundary initialization (8x8 context)
        C1l_m, D2l_m = solid1d_mush_relax.core_boundary_mush(
            R8_view, (1, 2), rs, ρs, gs, μs, Ks, ω, ρₗs, Kls, Kds, αs, ηₗs, ϕs, ks, ρ_core, μ_core, κ_core, "liquid", n; G0=1.0, Y=[1,2,5,6,3,7,4,8]
        )
        @test size(C1l_m) == (4, 8)
        @test size(D2l_m) == (4, 8)
        @test any(R_mush[1] .!= 0.0)
    end

    @testset "3. Layer Propagation Routines" begin
        Nr = 5
        rs = collect(range(0.5, 1.0, length=Nr))
        ρs = fill(1.0, Nr)
        gs = fill(1.0, Nr)
        μs = fill(complex(1.0), Nr)
        Ks = fill(complex(2.0), Nr)
        
        # Mushy vectors
        ρₗs = fill(0.9, Nr)
        Kls = fill(1.5, Nr)
        Kds = fill(complex(1.2), Nr)
        αs  = fill(complex(0.8), Nr)
        ηₗs = fill(1.0, Nr)
        ϕs  = fill(0.1, Nr)
        ks  = fill(1e-4, Nr)

        R = [zeros(precc, 8, 8) for _ in 1:Nr]
        idx_6 = [1, 2, 3, 5, 6, 7]
        R6_view = [view(R[i], idx_6, idx_6) for i in 1:Nr]
        R8_view = [view(R[i], 1:8, 1:8) for i in 1:Nr]

        # Test solid layer propagation engine
        Cn_l_6 = zeros(precc, 3, 6)
        Dnp_l_6 = zeros(precc, 3, 6)
        Cn_l_out6, Dnp_l_out6 = solid1d_mush_relax.propagate_solid(
            R6_view, Cn_l_6, Dnp_l_6, (2, 4), rs, ρs, gs, μs, Ks, ω, n; G0=1.0, Y=[1,2,4,5,3,6]
        )
        @test size(Cn_l_out6)  == (3, 6)
        @test size(Dnp_l_out6) == (3, 6)
        @test any(R[2] .!= 0.0)

        # Test mushy layer propagation engine
        Cn_l_8 = zeros(precc, 4, 8)
        Dnp_l_8 = zeros(precc, 4, 8)
        Cn_l_out8, Dnp_l_out8 = solid1d_mush_relax.propagate_mush(
            R8_view, Cn_l_8, Dnp_l_8, (2, 4), rs, ρs, gs, μs, Ks, ω, ρₗs, Kls, Kds, αs, ηₗs, ϕs, ks, n; G0=1.0, Y=[1,2,5,6,3,7,4,8]
        )
        @test size(Cn_l_out8)  == (4, 8)
        @test size(Dnp_l_out8) == (4, 8)
    end

    @testset "4. Interface Transitions Execution Loops" begin
        Nr = 5
        rs = collect(range(0.5, 1.0, length=Nr))
        ρs = fill(1.0, Nr)
        gs = fill(1.0, Nr)
        μs = fill(complex(1.0), Nr)
        Ks = fill(complex(2.0), Nr)
        
        ρₗs = fill(0.9, Nr)
        Kls = fill(1.5, Nr)
        Kds = fill(complex(1.2), Nr)
        αs  = fill(complex(0.8), Nr)
        ηₗs = fill(1.0, Nr)
        ϕs  = fill(0.1, Nr)
        ks  = fill(1e-4, Nr)

        R = [zeros(precc, 8, 8) for _ in 1:Nr]
        for i in 1:Nr
            R[i] .= Matrix{precc}(I, 8, 8)
        end
        R8_view = [view(R[i], 1:8, 1:8) for i in 1:Nr]

        # FIXED: Mushy layer solutions have 4 rows and 8 columns 
        Cn_l_in = zeros(precc, 4, 8)
        Dnp_l_in = zeros(precc, 4, 8)
        
        # Test Interface Transition: Mushy -> Solid
        Cn_l_out, Dnp_l_out = solid1d_mush_relax.interface_mush_solid(
            R8_view, Cn_l_in, Dnp_l_in, (2, 3), rs, ρs, gs, μs, Ks, ω, ρₗs, Kls, Kds, αs, ηₗs, ϕs, ks, n; G0=1.0, Y=[1,2,5,6,3,7,4,8]
        )
        @test size(Cn_l_out)  == (3, 6)
        @test size(Dnp_l_out) == (3, 6)

        # Test Interface Transition: Solid -> Mushy
        Cn_l_out2, Dnp_l_out2 = solid1d_mush_relax.interface_solid_mush(
            R8_view, Cn_l_out, Dnp_l_out, (3, 4), rs, ρs, gs, μs, Ks, ω, ρₗs, Kls, Kds, αs, ηₗs, ϕs, ks, n; G0=1.0, Y=[1,2,5,6,3,7,4,8]
        )
        @test size(Cn_l_out2)  == (4, 8)
        @test size(Dnp_l_out2) == (4, 8)
    end

    @testset "5. Surface Boundaries & System Resolutions" begin
        Nr = 5
        rs = collect(range(0.5, 1.0, length=Nr))
        ρs = fill(1.0, Nr)
        gs = fill(1.0, Nr)
        μs = fill(complex(1.0), Nr)
        Ks = fill(complex(2.0), Nr)

        # Initialize propagation history matrices with Identity properties
        R = [Matrix{precc}(I, 8, 8) for _ in 1:Nr]
        idx_6 = [1, 2, 3, 5, 6, 7]
        R6_view = [view(R[i], idx_6, idx_6) for i in 1:Nr]
        R8_view = [view(R[i], 1:8, 1:8) for i in 1:Nr]

        # FIXED: Populated with rand() instead of zeros to prevent SingularException(1)
        CNm_l6 = rand(precc, 3, 6)
        DN_l6  = rand(precc, 3, 6)
        y_t6, y_l6 = solid1d_mush_relax.surface_boundary(
            R6_view, CNm_l6, DN_l6, (Nr-1, Nr), rs, ρs, gs, μs, Ks, ω, n; G0=1.0, Y=[1,2,4,5,3,6]
        )
        @test length(y_t6) == 6
        @test length(y_l6) == 6

        # FIXED: Random complex population for the mushy surface boundary run
        CNm_l8 = rand(precc, 4, 8)
        DN_l8  = rand(precc, 4, 8)
        y_t8, y_l8 = solid1d_mush_relax.surface_boundary_mush(
            R8_view, CNm_l8, DN_l8, (Nr-1, Nr), rs, ρs, gs, μs, Ks, ω, n; G0=1.0, Y=[1,2,5,6,3,7,4,8]
        )
        @test length(y_t8) == 8
        @test length(y_l8) == 8
    end
end