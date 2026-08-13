using Test
using LinearAlgebra
using Obliqua
using Plots

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
            
            power_prf_expt = [0.0, 2.8587303613385215e-17, 2.968521783298822e-17, 2.965109898511842e-17, 3.0088133083343907e-17, 3.063751071681457e-17, 3.506150772139225e-18, 1.8107520734263888e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 2.159878395470904e6
            imag_k2_expt    = [0.002959435290502846]

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

            power_prf_expt = [0.0, 6.836052936946304e-21, 7.134693310200378e-21, 6.80758297455825e-21, 6.854349159327097e-21, 6.979328732637413e-21, 9.845458936642574e-22, 5.891946493834706e-21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 527.5744088946691
            imag_k2_expt    = [7.228751059888524e-7]

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

            power_prf_expt = [0.0, 1.295668025636027e-18, 1.4348106817840115e-18, 1.498477529778287e-18, 1.5933269782322786e-18, 1.7024312107623153e-18, 1.6970017972945087e-19, 2.1904520388500575e-18, 2.3649316137156205e-17, 6.415141327202965e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 1.2934772890254976e6
            imag_k2_expt    = [0.001772304563364745]

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

            power_prf_expt = [0.0, 3.1940736972941545e-22, 3.5511587570600433e-22, 3.539561453098292e-22, 3.731525459766172e-22, 3.9841939181049334e-22, 4.884693770760742e-23, 7.2981134297204295e-22, 8.920908913933924e-21, 2.4325781224855244e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 323062.9235271194
            imag_k2_expt    = [0.00044265631757047625]

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
        cfg["orbit"]["obliqua"]["visc_sus"] = 5e10

        @testset "Andrade Rheology" begin
            cfg["orbit"]["obliqua"]["material_mu"] = "andrade"

            omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
                load.load_interior_mush_full(interior_json_path, false)

            perm      = Obliqua.interior.get_permeability(phi, cfg)
            perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
            bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

            power_prf_expt = [2.893448766976472e-17, 2.0873895247427394e-19, 8.976980979485759e-20, 8.346566445443389e-20, 8.073093789859302e-20, 8.046109918164386e-20, 6.731136442243804e-20, 1.0557967908222806e-18, 9.616823605063416e-18, 8.688680540164449e-17, 8.175831691668493e-17, 9.847317792637044e-17, 1.0815505786481383e-16, 1.069701479490247e-16, 1.0010206272311943e-16, 8.668891445651362e-17, 7.558118583377667e-17, 6.815259446672906e-17, 6.094458393870817e-17, 5.556215133756868e-17, 5.317655500342848e-17, 5.470263217381461e-17, 5.998793303234427e-17, 6.840811958614819e-17, 8.499962317189811e-17, 1.08868075369319e-16, 1.370636578214253e-16, 1.6482395663106023e-16, 1.7999023517875913e-16, 1.7009200696021675e-16, 1.2693604824304955e-16, 6.23427442431497e-17, 1.3560409117722137e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 3.334150314552537e7
            imag_k2_expt    = [0.04568406316493879]

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

            power_prf_expt = [3.5458752212883464e-17, 9.370431858440995e-23, 3.0291617546858735e-23, 2.982319644472709e-23, 3.1546501214811024e-23, 3.4250715498155014e-23, 8.692434753175485e-24, 8.112004881217679e-23, 1.0436805728200542e-17, 5.730792799979111e-17, 9.277217486047095e-18, 6.6677737863741616e-18, 6.292579309081926e-18, 5.273953601588074e-18, 3.989055231854982e-18, 2.659236317918259e-18, 1.6412034195947134e-18, 1.1264169903582682e-18, 8.057809321930727e-19, 7.107672569893748e-19, 7.497586293655229e-19, 8.972004642364159e-19, 1.0387642784727909e-18, 1.0810818737609171e-18, 1.1997570954922017e-18, 1.309788003615809e-18, 1.5654445481830434e-18, 1.9698527908486245e-18, 2.6879943539043226e-18, 3.9664054540061386e-18, 6.317395658011158e-18, 1.0427040703073473e-17, 1.7210473254675018e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            power_blk_expt = 2.7551171681292225e6
            imag_k2_expt    = [0.003775023165160282]

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

        power_prf_expt = [0.0, 1.1839815191662313e-14, 1.2294468997514884e-14, 1.2280369767344691e-14, 1.2461339647784897e-14, 1.2688823660407206e-14, 1.4520148592274163e-15, 7.499709678962462e-15, 7.736557105428425e-15, 2.6609954153341316e-15, 9.114041548345337e-16, 3.1086540654469654e-16, 1.055970753143408e-16, 3.572499221772785e-17, 1.2038032526639168e-17, 4.04039015395603e-18, 1.3508162941942807e-18, 4.498778569021704e-19, 1.4925815380559769e-19, 4.933390079224491e-20, 1.6245584340834563e-20, 5.329987578805413e-21, 1.742349703078666e-21, 5.67517840045946e-22, 1.8419393069386982e-22, 5.957152312449675e-23, 1.9199294116797106e-23, 6.1663808769128325e-24, 1.9737362913209e-24, 6.296155525065903e-25, 2.0017193832381124e-25, 6.342868210968698e-26, 2.0032533324591136e-26, 6.30618498106699e-27, 1.9787480151608687e-27, 6.188989915049269e-28, 1.9295934968424324e-28, 5.997079050030044e-29, 1.858029065087104e-29, 5.738719190202773e-30, 1.7670039678597106e-30, 5.4241324468831985e-31, 1.6599830212731182e-31, 5.064851190027978e-32, 1.5407369731151484e-32, 4.6730269920261086e-33, 1.4131427881952589e-33, 4.261036971421258e-34, 1.2816239550591233e-34, 3.86182468926184e-35, 1.2214798455779766e-35, 5.912152842924285e-36, 9.463661106743774e-36, 2.902860446353964e-35, 9.807856280133717e-35, 3.3526028937195446e-34, 1.1502664624207894e-33, 3.958383353142417e-33, 1.3661841481581326e-32, 4.7289699912733585e-32, 1.6416700266854622e-31, 5.71562477987592e-31, 1.9957082049830574e-30, 6.988483257086165e-30, 2.4542559075191965e-29, 8.643778481077045e-29, 3.053047172945141e-28, 1.0814525129353158e-27, 3.8417073840130144e-27, 1.3686209146553728e-26, 4.8897174318694956e-26, 1.751966770747124e-25, 6.295209490054534e-25, 2.268494915062136e-24, 8.198023481201838e-24, 2.97115975485753e-23, 1.079918077609216e-22, 3.936472120550978e-22, 1.4390643616511384e-21, 5.2761416576805386e-21, 1.9400954410525342e-20, 7.154995451411694e-20, 2.646587508699501e-19, 9.818999594559296e-19, 3.654225766744993e-18, 1.364179186726399e-17, 5.108789157396564e-17, 1.919392341026692e-16, 7.235148326663033e-16, 2.7366590400966954e-15, 1.0388384225614868e-14, 1.0532040117949126e-14, 1.069254481915923e-14, 1.0873852449000426e-14, 1.108088550763046e-14, 9.18918563012912e-15, 7.202303249026087e-15, 5.674876263134549e-15, 4.509006923637211e-15]
        power_blk_expt = 3.440555232391258e9
        imag_k2_expt    = [0.013240088341154182, 0.01283405478099495, 0.012431153843255759, 0.012031733618704925, 0.011636206202843032, 0.011245063588260013, 0.010858898719815609, 0.010478433836753922, 0.01010455930945429]

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

        analytic_rtol = 5e-2

        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)

        idx = 5     # corresponds to a solid layer

        n_points    = 100
        radius_homo = collect(range(0, stop=radius[end], length=n_points+1))
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

                # Fetch spectral frequency axis (sigma) alongside outputs if available, 
                # or compute/extract forcing frequencies:
                _, power_blk_ref, _, sigma, LNk_ref = Obliqua.run_tides(
                    omega, axial, ecc, sma, S_mass,
                    rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                    phi_homo, perm_homo, cfg
                )
                imag_k2_ref = -imag.(LNk_ref)

                # Initialize plot canvas for log10-log10 comparison
                plt = plot(
                    xscale = :log10,
                    yscale = :log10,
                    xlabel = "Forcing Frequency σ [rad/s]",
                    ylabel = "- Im(k₂)",
                    title  = "Analytical Limit Comparison ($material_mu)",
                    legend = :bottomleft
                )

                # Plot reference line (solid0d)
                plot!(plt, sigma, imag_k2_ref, label="solid0d (reference)", lw=2, linestyle=:dash, color=:black)

                # --- Discretised 1D models ---
                for model in ("solid1d", "solid1d-relax", "solid1d-mush-relax")
                    @testset "$model vs analytical limit" begin
                        cfg["orbit"]["obliqua"]["module_solid"] = model

                        _, power_blk, _, sigma, LNk = Obliqua.run_tides(
                            omega, axial, ecc, sma, S_mass,
                            rho_homo, radius_homo, visc_homo, shear_homo, bulk_homo, bulkd_homo,
                            phi_homo, perm_homo, cfg
                        )
                        imag_k2 = -imag.(LNk)

                        # Overlay model curve onto the plot
                        plot!(plt, sigma, imag_k2, label=model, marker=:circle, ms=3)

                        @test isapprox(power_blk, power_blk_ref; rtol=analytic_rtol)
                        @test length(imag_k2) == length(imag_k2_ref)
                        @test all(isapprox.(imag_k2, imag_k2_ref; rtol=analytic_rtol))
                    end
                end

                # --- Save Plot to OUT_DIR ---
                mkpath(OUT_DIR) # Ensure directory exists
                plot_filename = joinpath(OUT_DIR, "analytical_limit_comparison_$(material_mu).png")
                savefig(plt, plot_filename)

            end
        end
    end
end