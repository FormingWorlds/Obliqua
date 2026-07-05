using Test
using LinearAlgebra
using Obliqua
using Obliqua.solid0d

@testset "solid0d Module Tests" begin

    @testset "compute_solid_lovenumbers" begin
        # Setup typical planetary parameters (Earth-like mock variables)
        μc = complex(1.0e11, 1.0e9) # Complex shear modulus (with anelastic dissipation)
        mass_tot = 5.972e24         # Planet mass in kg
        R = 6.371e6                 # Planet radius in meters
        
        # Test Case 1: Standard Tidal Degree n = 2
        n2 = 2
        k2_T, k2_L = compute_solid_lovenumbers(μc, mass_tot, R, n2)
        
        # Verify return types conform to specified outputs
        @test k2_T isa ComplexF64
        @test k2_L isa ComplexF64
        
        # Analytical ground truth check for n = 2 formulas
        An_expected = 4. * (2. * 2^2 + 4. * 2 + 3.) / (3. * 2 * G * mass_tot^2) * π * R^4
        μc_n_expected = An_expected * μc
        factor_expected = 1. / (1. + μc_n_expected)
        
        @test k2_T ≈ (factor_expected * 3. / (2. * (2 - 1.)))
        @test k2_L ≈ (factor_expected * -1.0)

        # Test Case 2: Special Degree n = 1 Boundary Condition
        n1 = 1
        k1_T, k1_L = compute_solid_lovenumbers(μc, mass_tot, R, n1)
        
        @test k1_T == 0.0
        
        # Test Case 3: Purely elastic behavior (imaginary part = 0)
        μc_elastic = complex(8.0e10, 0.0)
        k2_T_el, k2_L_el = compute_solid_lovenumbers(μc_elastic, mass_tot, R, n2)
        @test imag(k2_T_el) == 0.0
        @test imag(k2_L_el) == 0.0
    end

    @testset "mean_cmu (Hill-Averaging)" begin
        # Test Case 1: Uniform profile consistency
        # If all layers have the identical shear modulus, the average must exactly equal that value
        r_uniform = [0.0, 1.0e6, 2.0e6, 3.0e6]
        val = complex(5.0e10, 2.0e8)
        μc_uniform = [val, val, val] # length must match length(r) - 1
        
        @test mean_cmu(μc_uniform, r_uniform) ≈ val

        # Test Case 2: Basic Mathematical verification with two distinctive layers
        r_two = [1.0e6, 2.0e6, 3.0e6]
        μc_two = [complex(2.0e10, 0.0), complex(4.0e10, 0.0)]
        
        # Manually compute volumes and weight constraints
        V1 = (4/3) * π * (2.0e6^3 - 1.0e6^3)
        V2 = (4/3) * π * (3.0e6^3 - 2.0e6^3)
        Vtot = V1 + V2
        f1 = V1 / Vtot
        f2 = V2 / Vtot
        
        # Mathematical Voigt, Reuss, and Hill configurations
        μV_manual = (f1 * 2.0e10) + (f2 * 4.0e10)
        μR_manual = 1.0 / ((f1 / 2.0e10) + (f2 / 4.0e10))
        μHill_manual = 0.5 * (μV_manual + μR_manual)
        
        @test real(mean_cmu(μc_two, r_two)) ≈ μHill_manual
        
        # Test Case 3: Error Handling / Dimension Mismatch Assertions
        # A 4-element radial vector creates 3 shells, which cannot be broadcasted with 2 moduli entries.
        r_mismatch = [0.0, 1.0, 2.0, 3.0]      # produces 3 shells (lengths: 3)
        μc_mismatch = [complex(1.0), complex(2.0)] # length: 2

        @test_throws DimensionMismatch mean_cmu(μc_mismatch, r_mismatch)
    end
end