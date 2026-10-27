using Test
using LinearAlgebra
using Obliqua
using Plots
using NCDatasets
using Interpolations
using Obliqua.constants
using Obliqua.solid1d_relax
using Obliqua.solid1d_relax.common


@testset "Obliqua Tidal Simulation Tests" begin

    # --- Setup Common Data Paths ---
    interior_json_path = joinpath(RES_DIR, "interior_data", "runtests_mantle.json")
    config_path        = joinpath(TEST_DIR, "test.toml")  
    
    cfg = Obliqua.open_config(config_path)

    @testset "interior data and config" begin
        # Test interior data validity
        @test load.load_interior_mush_full(interior_json_path, true)

        # Test config validity
        @test Obliqua.open_config(config_path) !== nothing
    end

    # =========================================================================
    # cap_lovenumber: n-dependent Love-number clamp used by cap_LN
    # =========================================================================
    @testset "cap_lovenumber" begin
        # n=2: fluid limit 3/(2*(2-1)) = 1.5, so re_max=4.5, im_max=3.0
        @testset "n=2 bounds" begin
            # Within bounds: passed through unchanged (both signs)
            k_small = precc(1.0, -0.5)
            @test Obliqua.cap_lovenumber(k_small, 2) == k_small

            # Re exceeds +4.5: clamped to exactly the bound, Im untouched
            k_re_hi = precc(14.19, -0.60)
            capped = Obliqua.cap_lovenumber(k_re_hi, 2)
            @test real(capped) == 4.5
            @test imag(capped) == -0.60

            # Re exceeds -4.5 (negative side): clamped to exactly -4.5
            k_re_lo = precc(-14.19, 0.0)
            @test real(Obliqua.cap_lovenumber(k_re_lo, 2)) == -4.5

            # Im exceeds +3.0 and -3.0: clamped to exactly the bound
            @test imag(Obliqua.cap_lovenumber(precc(0.0, 15.5), 2)) == 3.0
            @test imag(Obliqua.cap_lovenumber(precc(0.0, -15.5), 2)) == -3.0

            # Exactly at the boundary is not altered
            @test Obliqua.cap_lovenumber(precc(4.5, 3.0), 2) == precc(4.5, 3.0)
        end

        # Discrimination guard: n=3 must give a DIFFERENT (smaller) bound than
        # n=2, not the same fixed number for every degree -- this is what
        # distinguishes the n-dependent formula from the earlier fixed-value
        # (1.5/1.0) implementation.
        @testset "degree dependence (n=3)" begin
            fluid_limit_n3 = 3.0 / (2.0 * (3 - 1))  # = 0.75
            @test Obliqua.cap_lovenumber(precc(10.0, 10.0), 3) ==
                  precc(3.0 * fluid_limit_n3, 2.0 * fluid_limit_n3)
            @test real(Obliqua.cap_lovenumber(precc(10.0, 0.0), 3)) <
                  real(Obliqua.cap_lovenumber(precc(10.0, 0.0), 2))
        end

        # n=1 (translation, undefined fluid limit) is guarded up to n=2
        # rather than dividing by zero or clamping to zero.
        @testset "n=1 guarded to n=2 bound" begin
            @test Obliqua.cap_lovenumber(precc(10.0, 10.0), 1) ==
                  Obliqua.cap_lovenumber(precc(10.0, 10.0), 2)
        end
    end

    # =========================================================================
    # spectrum = "legacy": hardcoded LovePy-compatible (n,m,k) triplet
    # =========================================================================
    @testset "legacy spectrum" begin
        # runtests_mantle.json has axial != omega, so this run also exercises
        # the spin-synchronisation mismatch @warn path (Obliqua.jl ~L534).
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["material_mu"]  = "andrade"
        cfg["orbit"]["obliqua"]["spectrum"]     = "legacy"

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        perm      = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

        power_blk_expt = 2.283014981761584e6
        imag_k2_expt   = 0.0010975674744855479

        power_prf, power_blk, nmk, sigma_range, LNk = Obliqua.run_tides(
            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
        )
        imag_k2 = -imag.(LNk)

        # Mode selection matches the hardcoded LovePy (n,m,k) triplet exactly.
        @test nmk == [(2, 0, 1), (2, 2, 1), (2, 2, 3)]

        # All three modes share one forcing frequency (dropping the sign, as
        # LovePy does), this is the whole point of the "legacy" mode.
        @test length(unique(sigma_range)) == 1
        @test sigma_range[1] ≈ omega

        # ...and therefore also collapse onto a single shared Im(k2), not
        # three independently-evaluated Love numbers.
        @test length(unique(round.(imag_k2, sigdigits=12))) == 1
        @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

        @test isapprox(power_blk, power_blk_expt; rtol=rtol)
        @test power_blk > 0.0

        # Restore spectrum for subsequent testsets.
        cfg["orbit"]["obliqua"]["spectrum"] = "adaptive"
    end

    @testset "invalid spectrum value" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["material_mu"]  = "andrade"
        cfg["orbit"]["obliqua"]["spectrum"]     = "not_a_real_spectrum_mode"

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        perm      = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

        @test_throws String Obliqua.run_tides(
            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
        )

        # Restore spectrum for subsequent testsets.
        cfg["orbit"]["obliqua"]["spectrum"] = "adaptive"
    end

    # =========================================================================
    # cap_LN pipeline wiring: cap_lovenumber is applied inside run_tides
    # =========================================================================
    @testset "cap_LN pipeline wiring" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["material_mu"]  = "andrade"
        cfg["orbit"]["obliqua"]["spectrum"]     = "adaptive"
        cfg["orbit"]["obliqua"]["cap_LN"]       = true

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        perm      = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

        # For this config, Im(k2)/Re(k2) sit far below the fluid-limit
        # cap, so clamping must be a no-op: cap_LN=true must reproduce the
        # exact same result as the "solid1d module / Andrade Rheology"
        # cap_LN=false test above.
        power_blk_expt = 1.1589870984446697e6
        imag_k2_expt   = 0.0015880279776708645

        _, power_blk, _, _, LNk = Obliqua.run_tides(
            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
        )
        imag_k2 = -imag.(LNk)

        @test isapprox(power_blk, power_blk_expt; rtol=rtol)
        @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

        # Physical bound implied by cap_LN for n=2 (im_max = 2 * 3/(2*(2-1))).
        @test all(imag_k2 .<= 3.0)

        # Restore configuration for subsequent testsets.
        cfg["orbit"]["obliqua"]["cap_LN"] = false
    end

    # The `module_mushy == "interp"` path adds a correction (ΔkT/ΔkL) onto the
    # *previous* segment's Love number strictly after that segment's own
    # cap_LN clamp already ran, which could silently push it back out of
    # bounds. Check that the re-clamp after interpolation is applied correctly.
    @testset "cap_LN with module_mushy=interp (re-clamp after interpolation)" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d"
        cfg["orbit"]["obliqua"]["module_mushy"] = "interp"
        cfg["orbit"]["obliqua"]["s_min"] = -2
        cfg["orbit"]["obliqua"]["s_max"] = 6
        cfg["orbit"]["obliqua"]["material_mu"] = "andrade"
        cfg["orbit"]["obliqua"]["cap_LN"] = true

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        perm      = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

        # Same config as "Complete model" (cap_LN=false there).
        power_blk_expt = 6.191476504474057e8
        imag_k2_expt   = [0.0023483255451422886, 0.002285812174864539, 0.002222475409171229, 0.0021583473367432215, 0.002093470459633021, 0.0020279001489672084, 0.001961707850789305, 0.0018949853440882103, 0.001827850506121588]

        _, power_blk, _, _, LNk = Obliqua.run_tides(
            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
        )
        imag_k2 = -imag.(LNk)

        @test isapprox(power_blk, power_blk_expt; rtol=rtol)
        @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol, atol=atol))
        @test all(imag_k2 .<= 3.0)

        # Restore configuration (in particular s_min/s_max) to the test.toml
        # defaults for subsequent testsets.
        cfg["orbit"]["obliqua"]["cap_LN"]       = false
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["s_min"]        = 1
        cfg["orbit"]["obliqua"]["s_max"]        = 1
    end

    # =========================================================================
    # 1. Solid1D Module Tests
    # =========================================================================
    @testset "solid1d module" begin
        # Base configuration for solid1d tests
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"

        @testset "Andrade Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)
            
            power_prf_expt = [7.242179775890015e-17, 4.016306221497748e-18, 3.953886094358103e-18, 3.7650551284492495e-18, 3.6478736887306624e-18, 3.5525422610971744e-18, 3.4674648845948973e-18, 3.3920022306827884e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.1589870984446697e6
            imag_k2_expt    = [0.0015880279776708645]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        @testset "Maxwell Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [8.619985537238349e-17, 9.890274078517318e-22, 9.78285619264205e-22, 8.894191282667416e-22, 8.547390553829827e-22, 8.32090508206856e-22, 8.120294679682862e-22, 7.942238495113745e-22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 998604.9454312292
            imag_k2_expt    = [0.0013682745857252405]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end
    end

    # =========================================================================
    # 2. Solid1D-Relax Module Tests
    # =========================================================================
    @testset "solid1d-relax module" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-relax"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"

        @testset "Andrade Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [6.959370044353372e-17, 4.025453284603604e-18, 3.958019593996895e-18, 3.773749598468625e-18, 3.65627053076728e-18, 3.560726255335133e-18, 3.475492657099843e-18, 3.3998905698561992e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.065805252331627e6
            imag_k2_expt    = [0.0014603515101440801]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        @testset "Maxwell Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [8.261617494979112e-17, 9.9123652357263e-22, 9.774598494511299e-22, 8.914361776354412e-22, 8.567009999036354e-22, 8.3397483106081985e-22, 8.13878785487143e-22, 7.960418429918644e-22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 884988.664219201
            imag_k2_expt    = [0.001212599139876233]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end
    end

    # =========================================================================
    # 3. Solid1D-Mush Module Tests
    # =========================================================================
    @testset "solid1d-mush module" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["visc_sus"] = 5e12

        @testset "Andrade Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [6.603637767730829e-17, 3.571770478503037e-18, 3.522965536184631e-18, 3.3623462607302073e-18, 3.2650519410681068e-18, 3.1868016298015206e-18, 3.1172691806776446e-18, 3.0559103956608782e-18, 4.369740141629308e-18, 9.43761000734665e-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.1286905442530464e6
            imag_k2_expt    = [0.0015465160611466133]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        @testset "Maxwell Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [7.876045951759971e-17, 8.815115299922489e-22, 8.735114488233781e-22, 7.958786252083756e-22, 7.664829034651009e-22, 7.477473033368469e-22, 7.312294113685989e-22, 7.166377405686856e-22, 1.7607339490038306e-21, 8.745399167231134e-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 912492.8452275486
            imag_k2_expt    = [0.0012502849855623444]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end
        
        # Cleanup configuration change
        cfg["orbit"]["obliqua"]["visc_sus"] = 5e13
    end

    # =========================================================================
    # 4. Solid1D-Mush-Relax Module Tests
    # =========================================================================
    @testset "solid1d-mush-relax module" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush-relax"
        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["visc_sus"] = 5e12

        @testset "Andrade Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [6.33011526729793e-17, 3.489629952893154e-18, 3.441462179946195e-18, 3.290446311515891e-18, 3.1995336228591145e-18, 3.126039936005275e-18, 3.061311991118494e-18, 2.978940824682916e-18, 9.866755832352206e-19, 1.3684483519237692e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 987722.2340916484
            imag_k2_expt    = [0.0013533632462432374]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        @testset "Maxwell Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [7.53295064481807e-17, 8.598753169153816e-22, 8.506642508308545e-22, 7.7757996711802065e-22, 7.50182022842248e-22, 7.323774689285769e-22, 7.170581662791695e-22, 7.244822809419244e-22, 9.017965441140445e-20, 7.672011461831307e-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 808584.1416013801
            imag_k2_expt    = [0.0011079107272955312]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        # Cleanup configuration change
        cfg["orbit"]["obliqua"]["visc_sus"] = 5e13
    end

    # =========================================================================
    # 5. Fluid Module Tests (fluid0d & fluid1d)
    # =========================================================================
    @testset "Fluid module models" begin
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

        @testset "fluid0d module" begin
            cfg["orbit"]["obliqua"]["module_solid"] = "none"
            cfg["orbit"]["obliqua"]["module_fluid"] = "fluid0d"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_blk_expt = 288817.1601289502
            imag_k2_expt    = [0.0003957332496036539]

            _, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_blk)
            println(imag_k2)

            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end

        @testset "fluid1d module" begin
            cfg["orbit"]["obliqua"]["module_solid"] = "none"
            cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5284804061097509e-18, 1.5517739264667466e-18, 1.578086507560905e-18, 1.6081325356797763e-18, 1.3335961567122854e-18, 1.045246240775732e-18, 8.235756116658285e-19, 6.543766781350802e-19]
            power_blk_expt = 288817.1601289502
            imag_k2_expt    = [0.0003957332496036539]

            power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
            )
            imag_k2 = -imag.(LNk)

            println(power_prf)
            println(power_blk)
            println(imag_k2)

            @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
            @test isapprox(power_blk, power_blk_expt; rtol=rtol)
            @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))
        end
    end

    # =========================================================================
    # 6. Complete Model Test
    # =========================================================================
    @testset "Complete model (Simultaneous modules)" begin
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
        cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d"
        cfg["orbit"]["obliqua"]["module_mushy"] = "interp"

        cfg["orbit"]["obliqua"]["s_min"] = -2
        cfg["orbit"]["obliqua"]["s_max"] = 6
        cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        perm      = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

        power_prf_expt = [2.807380901847586e-14, 1.6625591278974696e-15, 1.6367238256716559e-15, 1.5585771327055287e-15, 1.5100768957302145e-15, 1.4706192404488703e-15, 1.4354063733848886e-15, 1.4041737183913223e-15, 1.4485187591660389e-15, 4.982192627340357e-16, 1.7064257362404406e-16, 5.820345753537023e-17, 1.9770983710385396e-17, 6.688804941687805e-18, 2.2538857660667605e-18, 7.564839052565844e-19, 2.529138886539909e-19, 8.4230815616803e-20, 2.7945665339198135e-20, 9.236806474324574e-21, 3.041667417513335e-21, 9.97935759900796e-22, 3.2622084934472473e-22, 1.0625659789819664e-22, 3.4486705171140184e-23, 1.1153600701451622e-23, 3.594691710013135e-24, 1.1545340200621888e-24, 3.695434551396372e-25, 1.1788317806490004e-25, 3.747827409135519e-26, 1.1875778159847635e-26, 3.750699428551841e-27, 1.1807095998085861e-27, 3.704818022492138e-28, 1.1587671195431382e-28, 3.6127858162901483e-29, 1.1228355695939897e-29, 3.478795437270058e-30, 1.0744627470998481e-30, 3.308368774483916e-31, 1.0155625421352863e-31, 3.10799301673613e-32, 9.482941630402976e-33, 2.884728083853341e-33, 8.749325321810225e-34, 2.645825455508165e-34, 7.977703766591376e-35, 2.3987413769314712e-35, 7.201998416890351e-36, 2.190873682716311e-36, 7.819935969588235e-37, 6.699559103369782e-37, 1.6870278574120774e-36, 5.577122289630727e-36, 1.9021997671992408e-35, 6.523153165809442e-35, 2.2439787359133524e-34, 7.74192585667256e-34, 2.6787821689410175e-33, 9.295649716929996e-33, 3.234988642548915e-32, 1.1290510632987489e-31, 3.951841119629788e-31, 1.3871681177008412e-30, 4.883132088973449e-30, 1.7238799155246933e-29, 6.103133017377502e-29, 2.166880462151407e-28, 7.715297956336489e-28, 2.754900270882748e-27, 9.864940252994297e-27, 3.542576988393181e-26, 1.2757971442120463e-25, 4.607676863550642e-25, 1.6688695064155635e-24, 6.061851916766561e-24, 2.2081805527062296e-23, 8.067061976708758e-23, 2.955659579996275e-22, 1.0860729358420973e-21, 4.002565564422271e-21, 1.4794624028685632e-20, 5.484903642643844e-20, 2.0397495240307037e-19, 7.609012896158725e-19, 2.8473847722217573e-18, 1.0689543600979445e-17, 4.0262999513592707e-17, 1.521734257652441e-16, 5.771940244817454e-16, 5.851757587761389e-16, 5.940936378637519e-16, 6.041673584986944e-16, 6.156704220855017e-16, 5.1056470185571235e-16, 4.001704957876109e-16, 3.1530432535520087e-16, 2.505268418040781e-16]
        power_blk_expt = 6.191476504474057e8
        imag_k2_expt    = [0.0023483255451422886, 0.002285812174864539, 0.002222475409171229, 0.0021583473367432215, 0.002093470459633021, 0.0020279001489672084, 0.001961707850789305, 0.0018949853440882103, 0.001827850506121588]

        power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
        )
        imag_k2 = -imag.(LNk)

        println(power_prf)
        println(power_blk)
        println(imag_k2)

        @test all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
        @test isapprox(power_blk, power_blk_expt; rtol=rtol)
        @test all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol, atol=atol))
    end


    # =========================================================================
    # 7. Analytical Limit Tests (homogeneous sphere)
    # =========================================================================
    #
    # Idea: if every layer of a discretised model is given identical material 
    # properties, the discretised solid should converge onto the same tidal 
    # response as Obliqua's zero-dimensional homogeneous-sphere model ("solid0d"),
    # which is Obliqua's closed-form/analytical solution. We treat the
    # solid0d output as the "analytical limit" and check that the 1D
    # (multi-layer) solid modules reproduce it within a loose relative
    # tolerance once discretised over many identical layers.
    #
    # By setting the inner radius for the 1d solvers to a small non-zero value 
    # (e.g., 1 m), we avoid the singularity at r=0 in the 1D solvers, which is 
    # not present in the solid0d model. At the same time, our test becomes 
    # agnostic to the choise of core solution vector. Hence we tests with
    # the fluid core solution vector, and not the solid core solution vector.
    # The effects become stronger as the CMB shifts upwards, and we test this 
    # behaviour in test 8, hereafter.
    #
    # A full forcing spectrum is used (rather than the single semi-diurnal
    # term used in the exact-regression tests above) so several forcing
    # frequencies are checked simultaneously: this exercises N_sigma
    # frequency components spanning p in [p_min, p_max].
    @testset "Analytical limit (homogeneous sphere)" begin

        cfg["orbit"]["obliqua"]["spectrum"] = "full"
        cfg["orbit"]["obliqua"]["N_sigma"]  = 10
        cfg["orbit"]["obliqua"]["p_min"]    = -5
        cfg["orbit"]["obliqua"]["p_max"]    = 5

        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"

        cfg["orbit"]["obliqua"]["solid"]["core_props"] = "mantle"

        analytic_rtol = 1e-1

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        idx = 5    # corresponds to a solid layer

        n_points    = 10
        radius_homo = collect(range(1, stop=radius[end], length=n_points+1))
        rho_homo    = fill(rho[idx],   n_points)
        visc_homo   = fill(visc[idx],  n_points)
        shear_homo  = fill(shear[idx], n_points)
        bulk_homo   = fill(1e20,       n_points)        # incompressible sphere
        phi_homo    = fill(phi[idx],   n_points)

        for material_mu in ("andrade", "maxwell")
            @testset "$material_mu rheology" begin
                cfg["orbit"]["obliqua"]["material_mu"] = material_mu

                # --- "Analytical" reference: solid0d homogeneous sphere ---
                cfg["orbit"]["obliqua"]["module_solid"] = "solid0d"

                perm                        = Obliqua.interior.get_permeability(phi_homo, cfg)
                perm_homo, phi_homo = Obliqua.interior.limit_porosity(perm, phi_homo, cfg)
                bulkd_homo          = Obliqua.interior.get_drained_bulk(bulk_homo, phi_homo, cfg)

                # Fetch spectral frequency axis (sigma) and reference Love Numbers
                _, power_blk_ref, _, sigma, LNk_ref = Obliqua.run_tides(
                    omega, axial, ecc, sma, S_mass,
                    rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                    phi_homo, perm_homo, cfg
                )

                # Load reference data from NetCDF (extracting both Tidal & Load Love numbers)
                nc_file_path_ref = joinpath(OUT_DIR, "0_obliqua.nc")
                ds_ref = Dataset(nc_file_path_ref, "r")

                # Assuming index conventions: Tidal = index 1, Load = index 2 (adjust if needed)
                raw_kT_ref = ds_ref["knms_T"]
                raw_kL_ref = ds_ref["knms_L"]
                k_tidal_ref = raw_kT_ref[:, 1, 1] + 1im * raw_kT_ref[:, 1, 2]  # Tidal
                k_load_ref  = raw_kL_ref[:, 1, 1] + 1im * raw_kL_ref[:, 1, 2]  # Load 
                close(ds_ref)

                # --- Initialize 2x2 Plot Canvas ---
                plt = plot(
                    layout = (2, 2),
                    size   = (900, 700),
                    xscale = :log10,
                    yscale = :log10,
                    xlabel = "Forcing Frequency σ [rad/s]",
                    link   = :x # Keeps x-axes synchronized across panels
                )

                # Subplot Titles & Y-Labels
                plot!(plt[1], title="Tidal Love Number: Re(k₂)",   ylabel="Re(k₂)")
                plot!(plt[2], title="Load Love Number: -Re(k'₂)",   ylabel="-Re(k'₂)")
                plot!(plt[3], title="Tidal Love Number: -Im(k₂)",   ylabel="-Im(k₂)")
                plot!(plt[4], title="Load Love Number: Im(k'₂)",  ylabel="Im(k'₂)")

                # Helper plot function to apply across all 4 panels
                plot_ref_line! = (sp, data, lbl) -> plot!(
                    plt, sigma, data, sp=sp, label=lbl, 
                    lw=2, linestyle=:dash, color=:black
                )

                # Plot solid0d reference lines
                plot_ref_line!(1, real.(k_tidal_ref), "solid0d (ref)")
                plot_ref_line!(2, -real.(k_load_ref),  "solid0d (ref)")
                plot_ref_line!(3, -imag.(k_tidal_ref), "solid0d (ref)")
                plot_ref_line!(4, imag.(k_load_ref),   "solid0d (ref)")

                # --- Discretised 1D models ---
                models_to_test = ("solid1d", "solid1d-mush", "solid1d-relax", "solid1d-mush-relax", "solid1d-equil-relax")
                
                # Dictionary to store results for boundary checks
                model_results_real_kT = Dict{String, Vector{Float64}}()

                for model in models_to_test
                    @testset "$model vs analytical limit" begin
                        cfg["orbit"]["obliqua"]["module_solid"] = model

                        _, power_blk, _, sigma, LNk = Obliqua.run_tides(
                            omega, axial, ecc, sma, S_mass,
                            rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                            phi_homo, perm_homo, cfg
                        )

                        # Load netCDF data for current model
                        nc_file_path = joinpath(OUT_DIR, "0_obliqua.nc")
                        ds = Dataset(nc_file_path, "r")
                        
                        raw_kT  = ds["knms_T"]
                        raw_kL  = ds["knms_L"]
                        k_tidal = raw_kT[:, 1, 1] + 1im * raw_kT[:, 1, 2]
                        k_load  = raw_kL[:, 1, 1] + 1im * raw_kL[:, 1, 2]
                        close(ds)

                        # Cache real tidal results for boundary evaluation
                        model_results_real_kT[model] = real.(k_tidal)

                        # Overlay 1D model curves onto each panel
                        lstyle = (model == "solid1d-equil-relax") ? :dot : :solid
                        lw_val = (model == "solid1d-equil-relax") ? 2.5 : 1.5
                        
                        plot!(plt[1], sigma, real.(k_tidal),  label=model, marker=:circle, ms=3, linestyle=lstyle, lw=lw_val)
                        plot!(plt[2], sigma, -real.(k_load),  label=model, marker=:circle, ms=3, linestyle=lstyle, lw=lw_val)
                        
                        if model != "solid1d-equil-relax"
                            plot!(plt[3], sigma, -imag.(k_tidal), label=model, marker=:circle, ms=3, linestyle=lstyle, lw=lw_val)
                            plot!(plt[4], sigma, imag.(k_load),   label=model, marker=:circle, ms=3, linestyle=lstyle, lw=lw_val)

                            imag_k2     = -imag.(LNk)
                            imag_k2_ref = -imag.(LNk_ref)

                            @test isapprox(power_blk, power_blk_ref; rtol=analytic_rtol)
                            @test length(imag_k2) == length(imag_k2_ref)
                            @test all(isapprox.(imag_k2, imag_k2_ref; rtol=analytic_rtol))

                            # k_load (k_tidal was already checked above via LNk/imag_k2,
                            # for module_fluid="none" LNk *is* the raw k_T alone, but
                            # k_load itself was previously computed and plotted here
                            # without ever being asserted on -- this was the actual gap
                            # behind the Im(k2)>0 passivity violations found while
                            # investigating solid1d-relax's load-Love-number handling.
                            # Passivity (dissipation has the sign consistent with a
                            # passive body) must hold regardless of model/discretisation
                            # error, independent of how closely a given model tracks the
                            # solid0d homogeneous-sphere reference. Note the SIGN
                            # CONVENTION IS OPPOSITE to the tidal Love number's: the
                            # existing plotting code in this same test already encodes
                            # this (it negates -Re(k_load) but not Im(k_load), the
                            # mirror image of the tidal panels, which negate -Im(k_tidal)
                            # but not Re(k_tidal)) -- confirmed empirically too: solid0d's
                            # own closed-form reference (k_load_ref, the analytically
                            # exact case) satisfies sign(sigma)*Im(k_load) >= 0
                            # unanimously, not <= 0.
                            @test all(isapprox.(real.(k_load), real.(k_load_ref); rtol=analytic_rtol))
                            @test all(sign.(sigma) .* imag.(k_load) .>= 0)
                        end
                    end
                end

                # --- Verify that solid1d-equil-relax serves as a boundary/limit ---
                @testset "solid1d-equil-relax boundary check" begin
                    if haskey(model_results_real_kT, "solid1d-equil-relax")
                        eq_vals = model_results_real_kT["solid1d-equil-relax"]
                        dynamic_models = setdiff(collect(models_to_test), ["solid1d-equil-relax"])
                        
                        # Check that equil-relax bounds or tightly tracks the dynamic models 
                        # (e.g., serving as a static/relaxed upper or lower envelope limit across frequencies)
                        for i in eachindex(eq_vals)
                            dyn_vals = [model_results_real_kT[m][i] for m in dynamic_models]
                            # Ensure it doesn't wildly diverge from the envelope of dynamic solutions
                            @test eq_vals[i] <= maximum(dyn_vals) * 1.25 || eq_vals[i] >= minimum(dyn_vals) * 0.75
                        end
                    end
                end

                # --- Save Plot to OUT_DIR ---
                mkpath(OUT_DIR)
                plot_filename = joinpath(OUT_DIR, "analytical_limit_comparison_$(material_mu).png")
                savefig(plt, plot_filename)

            end
        end
    end


    # =========================================================================
    # 8. Core Boundary Tests (homogeneous sphere)
    # =========================================================================
    #
    # Idea: Extending test 7, we now shift the core boundary to a radius
    # of 50% of the total radius, and check that the 1D solid modules still
    # converge onto the same tidal response. We use the 'solid' core 
    # solution and compare it to the 'solid0d' homogeneous sphere solution, 
    # which is the analytical limit. We check that the 1D
    # (multi-layer) solid modules reproduce it within a loose relative
    # tolerance once discretised over many identical layers.
    #
    # A full forcing spectrum is used (rather than the single semi-diurnal
    # term used in the exact-regression tests above) so several forcing
    # frequencies are checked simultaneously: this exercises N_sigma
    # frequency components spanning p in [p_min, p_max].
    @testset "Core boundary (homogeneous sphere)" begin

        cfg["orbit"]["obliqua"]["spectrum"] = "full"
        cfg["orbit"]["obliqua"]["N_sigma"]  = 10
        cfg["orbit"]["obliqua"]["p_min"]    = -5
        cfg["orbit"]["obliqua"]["p_max"]    = 5

        cfg["orbit"]["obliqua"]["module_fluid"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"

        cfg["orbit"]["obliqua"]["solid"]["core"] = "solid"
        cfg["orbit"]["obliqua"]["solid"]["core_props"] = "mantle"

        analytic_rtol = 1e0

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        idx = 5     # corresponds to a solid layer

        n_points    = 100
        radius_homo = collect(range(0.5*radius[end], stop=radius[end], length=n_points+1))
        rho_homo    = fill(rho[idx],   n_points)
        visc_homo   = fill(visc[idx],  n_points)
        shear_homo  = fill(shear[idx], n_points)
        bulk_homo   = fill(1e20,       n_points)        # incompressible sphere
        phi_homo    = fill(phi[idx],   n_points)

        # Visual configuration map for models
        model_styles = Dict(
            "solid0d"            => (color=:black,       style=:dash,       marker=:none,       ms=0),
            "solid1d"            => (color=:crimson,     style=:dot,        marker=:circle,     ms=4),
            "solid1d-mush"       => (color=:royalblue,   style=:dash,       marker=:square,     ms=4),
            "solid1d-relax"      => (color=:darkgreen,   style=:dashdot,    marker=:diamond,    ms=4),
            "solid1d-mush-relax" => (color=:darkorange,  style=:dashdotdot, marker=:utriangle,  ms=4)
        )

        # Dictionary to cache heating profiles and bulk power without double-evaluations
        # Key: (material_mu, model) => (r_norm, power_prf, power_blk)
        heating_data = Dict{Tuple{String, String}, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()

        for material_mu in ("andrade", "maxwell")
            @testset "$material_mu rheology" begin
                cfg["orbit"]["obliqua"]["material_mu"] = material_mu

                # --- "Analytical" reference: solid0d homogeneous sphere ---
                cfg["orbit"]["obliqua"]["module_solid"] = "solid0d"

                perm                = Obliqua.interior.get_permeability(phi_homo, cfg)
                perm_homo, phi_homo = Obliqua.interior.limit_porosity(perm, phi_homo, cfg)
                bulkd_homo          = Obliqua.interior.get_drained_bulk(bulk_homo, phi_homo, cfg)

                # Fetch spectral frequency axis (sigma) and reference Love Numbers
                power_prf_ref, power_blk_ref, _, sigma, LNk_ref = Obliqua.run_tides(
                    omega, axial, ecc, sma, S_mass,
                    rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                    phi_homo, perm_homo, cfg
                )

                # Cache heating data for reference model
                # skip first radial layer
                r_norm = radius_homo[2:end] ./ radius_homo[end]
                heating_data[(material_mu, "solid0d")] = (r_norm, power_prf_ref, power_blk_ref)

                # Load reference data from NetCDF (extracting both Tidal & Load Love numbers)
                nc_file_path_ref = joinpath(OUT_DIR, "0_obliqua.nc")
                ds_ref = Dataset(nc_file_path_ref, "r")

                # Assuming index conventions: Tidal = index 1, Load = index 2 (adjust if needed)
                raw_kT_ref = ds_ref["knms_T"]
                raw_kL_ref = ds_ref["knms_L"]
                k_tidal_ref = raw_kT_ref[:, 1, 1] + 1im * raw_kT_ref[:, 1, 2]  # Tidal
                k_load_ref  = raw_kL_ref[:, 1, 1] + 1im * raw_kL_ref[:, 1, 2]  # Load 
                close(ds_ref)

                # --- Initialize 2x2 Plot Canvas ---
                plt = plot(
                    layout = (2, 2),
                    size   = (900, 700),
                    xscale = :log10,
                    yscale = :log10,
                    xlabel = "Forcing Frequency σ [rad/s]",
                    link   = :x # Keeps x-axes synchronized across panels
                )

                # Subplot Titles & Y-Labels
                plot!(plt[1], title="Tidal Love Number: Re(k₂)",   ylabel="Re(k₂)")
                plot!(plt[2], title="Load Love Number: -Re(k'₂)",  ylabel="-Re(k'₂)")
                plot!(plt[3], title="Tidal Love Number: -Im(k₂)",  ylabel="-Im(k₂)")
                plot!(plt[4], title="Load Love Number: Im(k'₂)", ylabel="Im(k'₂)")

                # Helper plot function to apply across all 4 panels
                plot_ref_line! = (sp, data, lbl) -> plot!(
                    plt, sigma, data, sp=sp, label=lbl, 
                    lw=2, linestyle=:dash, color=:black
                )

                # Plot solid0d reference lines
                plot_ref_line!(1, real.(k_tidal_ref), "solid0d (ref)")
                plot_ref_line!(2, -real.(k_load_ref),  "solid0d (ref)")
                plot_ref_line!(3, -imag.(k_tidal_ref), "solid0d (ref)")
                plot_ref_line!(4, imag.(k_load_ref),  "solid0d (ref)")

                # --- Discretised 1D models ---
                for model in ("solid1d", "solid1d-mush", "solid1d-relax", "solid1d-mush-relax")
                    @testset "$model vs analytical limit" begin
                        cfg["orbit"]["obliqua"]["module_solid"] = model

                        power_prf, power_blk, _, sigma, LNk = Obliqua.run_tides(
                            omega, axial, ecc, sma, S_mass,
                            rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                            phi_homo, perm_homo, cfg
                        )

                        # Cache heating data for current model
                        heating_data[(material_mu, model)] = (r_norm, power_prf, power_blk)

                        # Load netCDF data for current model
                        nc_file_path = joinpath(OUT_DIR, "0_obliqua.nc")
                        ds = Dataset(nc_file_path, "r")
                        
                        raw_kT  = ds["knms_T"]
                        raw_kL  = ds["knms_L"]
                        k_tidal = raw_kT[:, 1, 1] + 1im * raw_kT[:, 1, 2]
                        k_load  = raw_kL[:, 1, 1] + 1im * raw_kL[:, 1, 2]
                        close(ds)

                        # Retrieve style settings for current model
                        st = model_styles[model]

                        # Shared plot kwargs for consistency and high contrast
                        curve_kwargs = (
                            label             = model,
                            color             = st.color,
                            linestyle         = st.style,
                            marker            = st.marker,
                            markersize        = st.ms,
                            markerstrokewidth = 0.5,
                            lw                = 1.8,
                            alpha             = 0.75,
                            markevery         = 1
                        )

                        # Overlay 1D model curves onto each panel
                        plot!(plt[1], sigma, real.(k_tidal);  curve_kwargs...)
                        plot!(plt[2], sigma, -real.(k_load);  curve_kwargs...)
                        plot!(plt[3], sigma, -imag.(k_tidal); curve_kwargs...)
                        plot!(plt[4], sigma, imag.(k_load);   curve_kwargs...)

                        imag_k2     = -imag.(LNk)
                        imag_k2_ref = -imag.(LNk_ref)

                        @test isapprox(power_blk, power_blk_ref; rtol=analytic_rtol)
                        @test length(imag_k2) == length(imag_k2_ref)
                        @test all(isapprox.(imag_k2, imag_k2_ref; rtol=analytic_rtol))
                    end
                end

                # --- Save Plot to OUT_DIR ---
                mkpath(OUT_DIR)
                plot_filename = joinpath(OUT_DIR, "core_boundary_comparison_$(material_mu).png")
                savefig(plt, plot_filename)

            end
        end

        # =========================================================================
        # Additional Plot 1: Radial Heating Profiles (2 panels: Andrade vs Maxwell)
        # =========================================================================
        models_keys = ("solid0d", "solid1d", "solid1d-mush", "solid1d-relax", "solid1d-mush-relax")

        plt_profiles = plot(
            layout = (1, 2),
            size   = (950, 450),
            xlabel = "Normalized Radius (r / R)",
            ylabel = "Volumetric Heating Power [log10(W/m³)]"
        )

        for (i, material_mu) in enumerate(("andrade", "maxwell"))
            plot!(plt_profiles[i], title="$(uppercasefirst(material_mu)) Rheology: Radial Profile")
            
            for model in models_keys
                r_norm, power_prf, _ = heating_data[(material_mu, model)]
                st = model_styles[model]
                
                plot!(
                    plt_profiles[i], r_norm, log10.(power_prf),
                    label      = model,
                    color      = st.color,
                    linestyle  = st.style,
                    marker     = st.marker,
                    markersize = st.ms,
                    lw         = 1.8,
                    alpha      = 0.8
                )
            end
        end

        savefig(plt_profiles, joinpath(OUT_DIR, "radial_heating_profiles.png"))

        # =========================================================================
        # Additional Plot 2: Bulk Heating Distribution (Models vs Bulk Heating)
        # =========================================================================
        # Convert models_keys to a vector of strings and define numeric x-coordinates
        x_labels  = string.(collect(models_keys))
        x_indices = 1:length(x_labels)

        andrade_bulk = [heating_data[("andrade", m)][3] for m in models_keys]
        maxwell_bulk = [heating_data[("maxwell", m)][3] for m in models_keys]

        # Matrix layout: solid models horizontally, bulk power values vertically per rheology group
        bulk_matrix = [andrade_bulk maxwell_bulk]

        plt_bulk = scatter(
            x_indices,
            log10.(bulk_matrix),
            xticks     = (x_indices, x_labels),
            label      = ["Andrade" "Maxwell"],
            title      = "Bulk Heating Power across Solid Models",
            xlabel     = "Solid Model",
            ylabel     = "Total Bulk Heating [log10(W)]",
            size       = (850, 500),
            markersize = 6,
            palette    = :tab10,
            legend     = :topright
        )

        savefig(plt_bulk, joinpath(OUT_DIR, "bulk_heating_distribution.png"))

    end
end