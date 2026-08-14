using Test
using LinearAlgebra
using Obliqua
using Plots
using NCDatasets


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
            
            power_prf_expt = [1.4223292846297774e-15, 3.8867134531638245e-17, 4.021216844819114e-17, 4.004062618640702e-17, 4.051857756503588e-17, 4.115858556579879e-17, 4.1441808443294514e-18, 1.9802834680439397e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.9744913276554514e7
            imag_k2_expt    = [0.02705420512612409]

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

            power_prf_expt = [1.6921825815955382e-15, 9.359886440781222e-21, 9.731832656808001e-21, 9.255678791828123e-21, 9.292656874682616e-21, 9.438434333854014e-21, 1.168004125142122e-21, 6.447326525864462e-21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 2.009128776085858e7
            imag_k2_expt    = [0.02752880261954255]

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

            power_prf_expt = [7.991612406943655e-16, 2.0112864303716865e-17, 2.1031663579639105e-17, 2.1135660508079922e-17, 2.1565372323119404e-17, 2.2067418219489537e-17, 3.2197613221224515e-18, 1.527510143054395e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 9.947911671051508e6
            imag_k2_expt    = [0.013630489997885423]

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

            power_prf_expt = [9.39995660056893e-16, 4.884160891952851e-21, 5.1319536492709515e-21, 4.925459388280025e-21, 4.985616332965325e-21, 5.100623398713375e-21, 9.078248231422461e-22, 4.977497249012951e-21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 9.801043861173445e6
            imag_k2_expt    = [0.013429253770649935]

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

            power_prf_expt = [3.8579746198494463e-17, 1.3323037849703029e-18, 1.4761161534233225e-18, 1.5420875932937248e-18, 1.6394392061066468e-18, 1.7508362031329686e-18, 1.6593899874877974e-19, 2.0359406953573953e-18, 2.1975605019669924e-17, 5.961261942859147e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.6919082666113111e6
            imag_k2_expt    = [0.0023182291387341502]

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

            power_prf_expt = [4.6788895645919137e-17, 3.2866468512099903e-22, 3.656488963506418e-22, 3.6460970482090325e-22, 3.8435860677237195e-22, 4.1021117428892268e-22, 4.781907108751783e-23, 6.781032915461063e-22, 8.28672492690357e-21, 2.2615657276928835e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 882784.9144282046
            imag_k2_expt    = [0.0012095795926107072]

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

            power_prf_expt = [4.633237959114448e-17, 2.505662685066841e-19, 8.625929205555747e-20, 7.98926167330342e-20, 7.804638802561065e-20, 7.969154317956457e-20, 1.0463857914905028e-19, 1.6915777602897524e-18, 6.388697271255131e-18, 5.741663436242325e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.073767334752778e6
            imag_k2_expt    = [0.0014712610445662296]

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

            power_prf_expt = [4.48320853453346e-17, 6.291039177210475e-23, 2.1317770015202265e-23, 1.8963033923152144e-23, 1.8482132847680073e-23, 1.8962639225331219e-23, 2.910245354007299e-23, 5.454163687103323e-22, 3.425127722599139e-18, 4.9182137516862614e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 527064.8817992565
            imag_k2_expt    = [0.0007221769590603987]

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

            power_blk_expt = 5.198662290572516e6
            imag_k2_expt    = [0.007123134653500875]

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

            power_prf_expt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7512395381499228e-17, 2.793167490862105e-17, 2.840529703140165e-17, 2.894612058526975e-17, 2.4004514792321753e-17, 1.8814267868589255e-17, 1.482423643344309e-17, 1.1778686437017732e-17]
            power_blk_expt = 5.198662290572516e6
            imag_k2_expt    = [0.007123134653500875]

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

        power_prf_expt = [5.551965209735074e-13, 1.6141221914213976e-14, 1.6699147372091e-14, 1.6627465321654235e-14, 1.6825457225022463e-14, 1.7090765227814157e-14, 1.718667877156147e-15, 8.207063650304935e-15, 8.466249937191441e-15, 2.9119738872129418e-15, 9.973655288099364e-16, 3.401854588247475e-16, 1.1555672892537009e-16, 3.9094484665185465e-17, 1.3173429825917841e-17, 4.421469708167416e-18, 1.478221928699407e-18, 4.923092178909089e-19, 1.6333581179091874e-19, 5.398695166234708e-20, 1.7777827466525558e-20, 5.832698755966026e-21, 1.906683795289897e-21, 6.210447117714311e-22, 2.0156664429890025e-22, 6.519016108047458e-23, 2.1010123805122778e-23, 6.74797859053233e-24, 2.1598941912680023e-24, 6.889993260856947e-25, 2.190516578844033e-25, 6.941111771158263e-26, 2.1921952063418952e-26, 6.900968670213089e-27, 2.165378608440543e-27, 6.772720056935811e-28, 2.111587958160958e-28, 6.562708636251623e-29, 2.0332737471240368e-29, 6.279980916772706e-30, 1.9336633890299983e-30, 5.935723133930046e-31, 1.8165484928000425e-31, 5.542555361044136e-32, 1.6860554501387215e-32, 5.113775155626458e-33, 1.546426488076053e-33, 4.662913513154538e-34, 1.4024569677189314e-34, 4.224493489152402e-35, 1.331383802740971e-35, 6.289930923100525e-36, 9.744331189586592e-36, 2.9677590518253196e-35, 1.001748528964325e-34, 3.4229971560591927e-34, 1.1740247737031437e-33, 4.03872573577219e-33, 1.3933978585388317e-32, 4.821293532920577e-32, 1.6730384708915788e-31, 5.822358436981651e-31, 2.0320751356281128e-30, 7.112555260103682e-30, 2.4966362749749025e-29, 8.788700196651074e-29, 3.10265286224518e-28, 1.098446762706548e-27, 3.899968789874166e-27, 1.3886054980816842e-26, 4.9582915455296047e-26, 1.7754998382537867e-25, 6.3759584027734245e-25, 2.2961899059708784e-24, 8.292933678414047e-24, 3.003644687871403e-23, 1.091016957195519e-22, 3.9743010232373613e-22, 1.4519162678640861e-21, 5.319619756058623e-21, 1.9547227580361233e-20, 7.2038495216176e-20, 2.66274826773369e-19, 9.87177345285578e-19, 3.671157510452736e-18, 1.3694762279275003e-17, 5.124745890875036e-17, 1.9239126084706928e-16, 7.246566861091503e-16, 2.7388294901300264e-15, 1.0388384225614868e-14, 1.0532040117949126e-14, 1.069254481915923e-14, 1.0873852449000426e-14, 1.108088550763046e-14, 9.18918563012912e-15, 7.202303249026087e-15, 5.674876263134549e-15, 4.509006923637211e-15]
        power_blk_expt = 1.0342974673648409e10
        imag_k2_expt    = [0.0396102763376021, 0.0385022569474954, 0.037362841788000756, 0.036191520988917605, 0.034987920787144335, 0.03375183877217035, 0.03248328835819582, 0.031182555951041815, 0.029850276065565954]

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

        analytic_rtol = 1e-1

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        idx = 5     # corresponds to a solid layer

        n_points    = 10
        radius_homo = collect(range(1, stop=radius[end], length=n_points+1))
        rho_homo    = fill(rho[idx],   n_points)
        visc_homo   = fill(visc[idx],  n_points)
        shear_homo  = fill(shear[idx], n_points)
        bulk_homo   = fill(1e20,       n_points)        # incompressible sphere
        phi_homo    = fill(phi[idx],   n_points)

        cfg["struct"]["core_density"] = rho[idx]
        cfg["struct"]["core_shear"]   = shear[idx]
        cfg["struct"]["core_bulk"]    = 1e20

        for material_mu in ("andrade", "maxwell")
            @testset "$material_mu rheology" begin
                cfg["orbit"]["obliqua"]["material_mu"] = material_mu

                # --- "Analytical" reference: solid0d homogeneous sphere ---
                cfg["orbit"]["obliqua"]["module_solid"] = "solid0d"

                perm                = Obliqua.interior.get_permeability(phi_homo, cfg)
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
                plot!(plt[2], title="Load Love Number: Re(k'₂)",  ylabel="Re(k'₂)")
                plot!(plt[3], title="Tidal Love Number: -Im(k₂)",  ylabel="-Im(k₂)")
                plot!(plt[4], title="Load Love Number: -Im(k'₂)", ylabel="-Im(k'₂)")

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

                        # Overlay 1D model curves onto each panel
                        plot!(plt[1], sigma, real.(k_tidal),  label=model, marker=:circle, ms=3)
                        plot!(plt[2], sigma, -real.(k_load),  label=model, marker=:circle, ms=3)
                        plot!(plt[3], sigma, -imag.(k_tidal), label=model, marker=:circle, ms=3)
                        plot!(plt[4], sigma, imag.(k_load),   label=model, marker=:circle, ms=3)

                        imag_k2     = -imag.(LNk)
                        imag_k2_ref = -imag.(LNk_ref)

                        @test isapprox(power_blk, power_blk_ref; rtol=analytic_rtol)
                        @test length(imag_k2) == length(imag_k2_ref)
                        @test all(isapprox.(imag_k2, imag_k2_ref; rtol=analytic_rtol))
                    end
                end

                # --- Save Plot to OUT_DIR ---
                mkpath(OUT_DIR)
                plot_filename = joinpath(OUT_DIR, "analytical_limit_comparison_$(material_mu).png")
                savefig(plt, plot_filename)

            end
        end
    end
end