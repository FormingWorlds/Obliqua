#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get Obliqua root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"

# Include libraries
using LoggingExtras
using Obliqua
using Glob
using Test

@info "Begin Obliqua tests"

# Configure
LOG_LEVEL = Logging.Error

# Prepare
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

total  = 0
failed = 0

# which test suite to run?
# 0 - none
# 10 - fast
# 20 - all
suite::Int64 = 20
if length(ARGS)>0
    if ARGS[1] == "all"
        suite = 20
    elseif ARGS[1] == "fast"
        suite = 10
    elseif ARGS[1] == "none"
        suite = 0
    else
        if !isnothing(tryparse(Int64, ARGS[1]))
            suite = parse(Int64, ARGS[1])
        else
            @warn "Invalid test suite option '$(ARGS[1])'. Running all tests."
            suite = 20
        end
    end
end
suite = min(max(suite, 0), 20)
@info "Using suite $suite"

rm(OUT_DIR,force=true,recursive=true)
if !isdir(OUT_DIR) && !isfile(OUT_DIR)
    mkdir(OUT_DIR)
end

rtol = 1e-3
atol = 1e-18  

# Find test names
test_names = sort([replace(split(basename(f), ".jl")[1], "test_"=>"")
                        for f in glob("test_*.jl", TEST_DIR)])

# Collect tests
test_files = String[]
for test_name in test_names
    push!(test_files, joinpath(TEST_DIR, "test_$test_name.jl"))
end
@info "Collected: $(join(test_names, ", "))"
@info "Running tests..."

# test module imported
if suite >= 0
    @info " "
    @info "Testing if solid tides modules are imported"
    if isdefined(Obliqua, :run_solid0d)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    if isdefined(Obliqua, :run_solid1d)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    if isdefined(Obliqua, :run_solid1d_relax)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    if isdefined(Obliqua, :run_solid1d_mush)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    if isdefined(Obliqua, :run_solid1d_mush_relax)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 3
    @info "--------------------------"
end

if suite >= 0
    @info " "
    @info "Testing if fluid tides modules are imported"
    if isdefined(Obliqua, :run_fluid0d)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    if isdefined(Obliqua, :run_fluid1d)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 2
    @info "--------------------------"
end

if suite >= 0
    @info " "
    @info "Performing unit tests"

    # Run tests
    for test_file in test_files
        @info "Running '$(test_file)'"
        include(test_file)
    end

    @info "--------------------------"
end

if suite > 2
    # test interior data validity
    @info " "
    @info "Testing interior data validity"
    ok = load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", true)
    if ok
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # test config validity
    @info " "
    @info "Testing config validity"
    cfg = Obliqua.open_config("$TEST_DIR/test.toml")
    if cfg !== nothing
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"
end

if suite > 2   
    # test solid1d module with andrade rheology
    @info " "
    @info "Testing solid1d module with andrade rheology"

    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material_mu"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)
    
    power_prf_expt  = [0.0, 0.0, 5.100789685077343e-18, 5.0871195149156065e-18, 5.155039107907815e-18, 5.2428149328626275e-18, 5.701571310249043e-19, 2.895743389398542e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 301993.3820059709
    imag_k2_expt    = [0.0004137871252755984]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    @info " "
    @info "Testing solid1d module with maxwell rheology"

    # update config to use maxwell rheology
    cfg["orbit"]["obliqua"]["material_mu"]     = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 0.0, 1.2290767214931002e-21, 1.1709393963170431e-21, 1.1773812933061318e-21, 1.1974087289679709e-21, 1.6059955605007558e-22, 9.450706691475171e-22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 73.95592715362396
    imag_k2_expt    = [1.0133338118444123e-7]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end

if suite > 2   
    # test solid1d-relax module with andrade rheology
    @info " "
    @info "Testing solid1d-relax module with andrade rheology"

    # update config to use only solid1d-relax
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-relax"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material_mu"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 2.887982804636551e-18, 3.02227075349468e-18, 3.0394288616465934e-18, 3.1032772414702275e-18, 3.177408519534877e-18, 4.855956432156795e-19, 2.3447248162853347e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 229233.153746261
    imag_k2_expt    = [0.00031409207405958716]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    @info " "
    @info "Testing solid1d-relax module with maxwell rheology"

    # update config to use maxwell rheology
    cfg["orbit"]["obliqua"]["material_mu"]     = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 7.038722567353927e-22, 7.401306542343406e-22, 7.108414979921441e-22, 7.199742683540696e-22, 7.36995567723654e-22, 1.371172157469167e-22, 7.649457213312718e-22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 58.02811920965406
    imag_k2_expt    = [7.950932061298529e-8]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end

if suite > 4
    # test solid1d-mush module with andrade rheology
    @info " "
    @info "Testing solid1d-mush module with andrade rheology"

    # update config to use only solid1d-mush
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material_mu"]     = "andrade"

    # lower visc_sus to include mush
    visc_sus = cfg["orbit"]["obliqua"]["visc_sus"] = 5e12

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 0.0, 3.569577803765813e-18, 3.7378606580921296e-18, 3.983399333870779e-18, 4.264242639004841e-18, 4.3142238554531846e-19, 5.600003528864794e-18, 6.078684071116696e-17, 1.627256817757301e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 3.1645221680729627e6
    imag_k2_expt    = [0.004335984193097071]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    @info " "
    @info "Testing solid1d-mush module with maxwell rheology"

    # update config to use maxwell rheology
    cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 0.0, 8.874758543456532e-22, 8.869406161469306e-22, 9.371641594622808e-22, 1.0025409489837838e-21, 1.2478405731561284e-22, 1.875175305328289e-21, 2.3044238096976424e-20, 6.165902070797792e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 804040.11674267
    imag_k2_expt    = [0.0011016845677321097]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    visc_sus = cfg["orbit"]["obliqua"]["visc_sus"] = 5e13

end


if suite > 4
    # test solid1d-mush-relax module with andrade rheology
    @info " "
    @info "Testing solid1d-mush-relax module with andrade rheology"

    # update config to use only solid1d-mush-relax
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush-relax"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material_mu"]     = "andrade"

    # lower visc_sus to include mush
    visc_sus = cfg["orbit"]["obliqua"]["visc_sus"] = 5e10

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 9.486022666253308e-20, 4.0588388169682e-20, 3.766101384948525e-20, 3.642883896772203e-20, 3.639474279676227e-20, 3.154930594195203e-20, 5.133122007486018e-19, 4.1825341116642436e-18, 4.104856630206054e-17, 3.735799820884784e-17, 4.507308473687904e-17, 4.6192445071639857e-17, 4.5352937971520745e-17, 4.056595508613484e-17, 3.53971995001516e-17, 3.1130041711618394e-17, 2.691557389655362e-17, 2.467119574849353e-17, 2.3047408490882098e-17, 2.1259842561787253e-17, 2.1979124595524053e-17, 2.3750005233385358e-17, 2.7008681804655046e-17, 3.3856987874266826e-17, 4.3916864825784067e-17, 5.676266830633033e-17, 7.021200717883206e-17, 7.762928336400809e-17, 7.348633149608228e-17, 5.5059866355714375e-17, 2.7147439648528524e-17, 5.932262422323913e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 1.3996746667039031e7
    imag_k2_expt    = [0.019178147309368464]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    @info " "
    @info "Testing solid1d-mush-relax module with maxwell rheology"

    # update config to use maxwell rheology
    cfg["orbit"]["obliqua"]["material_mu"] = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 4.7859804317301677e-23, 1.6640371609255852e-23, 1.6448054666020436e-23, 1.7445145747654098e-23, 1.8968720803676234e-23, 4.238365559179562e-24, 3.936971824975529e-23, 4.451036298259694e-18, 2.9323097351280397e-17, 4.695322450403615e-18, 3.580160814065616e-18, 3.409651066391215e-18, 3.080974383790756e-18, 2.2339493951071856e-18, 1.4689567805610032e-18, 9.194887725106967e-19, 5.74482324287654e-19, 4.835047889660556e-19, 5.104077873862976e-19, 5.110814770911871e-19, 5.944345890083043e-19, 6.645176059999082e-19, 6.731126204847642e-19, 7.386183209103673e-19, 8.232570027905483e-19, 9.512924161724462e-19, 1.1684125808977266e-18, 1.596610774258235e-18, 2.381721116187253e-18, 3.8287343969745466e-18, 6.353451595189355e-18, 1.0430290600461818e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 1.326600050845652e6
    imag_k2_expt    = [0.0018176889102127136]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    visc_sus = cfg["orbit"]["obliqua"]["visc_sus"] = 5e13

end


if suite > 2
    # test fluid0d model
    @info " "
    @info "Testing fluid0d module with andrade rheology"
    
    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "none"
    cfg["orbit"]["obliqua"]["module_fluid"] = "fluid0d"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material_mu"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_blk_expt  = 5.167672882253204e6
    imag_k2_expt    = [0.007080673397902278]

    _, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    @info " "
    @info "Testing fluid1d module with andrade rheology"
    
    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "none"
    cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7655907378155416e-17, 2.807737397918556e-17, 2.855346664135607e-17, 2.909711127521295e-17, 2.4187219619345477e-17, 1.8967431361232358e-17, 1.4952738276867223e-17, 1.1886968308812462e-17]
    power_blk_expt  = 5.167672882253204e6
    imag_k2_expt    = [0.007080673397902278]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end

if suite > 10
    # test complete model
    @info " "
    @info "Testing all modules simultaneously with andrade rheology and interpolation for mushy layer"

    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d"
    cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d"
    cfg["orbit"]["obliqua"]["module_mushy"] = "interp"

    cfg["orbit"]["obliqua"]["s_min"]        = -2
    cfg["orbit"]["obliqua"]["s_max"]        = 6

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    perm      = Obliqua.interior.get_permeability(phi, cfg)
    perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
    bulkd     = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)

    power_prf_expt  = [0.0, 0.0, 2.1124506193530383e-15, 2.106794301650103e-15, 2.1349166926471606e-15, 2.1712599792776307e-15, 2.3610779298939566e-16, 1.1992826689985647e-15, 1.2371570702642999e-15, 4.2617076026443184e-16, 1.4617855511405637e-16, 4.992863143033294e-17, 1.6982669489719978e-17, 5.752755549704997e-18, 1.9408081209534226e-18, 6.521507258361804e-19, 2.1826972266029966e-19, 7.276791538600777e-20, 2.4166163030860044e-20, 7.994958835522106e-21, 2.635024575927128e-21, 8.652310259531574e-22, 2.8305858447284387e-22, 9.226452597450837e-23, 2.9965700553906005e-23, 9.697546722780387e-24, 3.127254156084172e-24, 1.0049493702780452e-24, 3.218255682063066e-25, 1.0270854568008316e-25, 3.2667555088979715e-26, 1.0355342483249376e-26, 3.271624598632958e-27, 1.0302129438990058e-27, 3.233462276664915e-28, 1.011576651804659e-28, 3.154508191974903e-29, 9.805683976585089e-30, 3.0384266276471903e-30, 9.385459936583588e-31, 2.8900742601550554e-31, 8.871959503937063e-32, 2.7151756380308966e-32, 8.284256691295213e-33, 2.519972879537751e-33, 7.642492484824555e-34, 2.3109376120042628e-34, 6.968767950897512e-35, 2.1007189156825448e-35, 6.496369265317593e-36, 2.6164669804332096e-36, 3.0630675371694664e-36, 8.657057193671852e-36, 2.8926343813974634e-35, 9.851523378154153e-35, 3.369871591274238e-34, 1.1562606178693972e-33, 3.979031338759478e-33, 1.3733111312441776e-32, 4.7536398208657056e-32, 1.6502342054133783e-31, 5.745441755170378e-31, 2.0061193122171595e-30, 7.024940414719192e-30, 2.4670591426709822e-29, 8.688870897134545e-29, 3.068974151368643e-28, 1.0870941784137263e-27, 3.8617486041935594e-27, 1.375760665384033e-26, 4.9152258566080765e-26, 1.761106340290514e-25, 6.328050012994375e-25, 2.2803290818860846e-24, 8.240790505654495e-24, 2.986659547237988e-23, 1.085551738324861e-22, 3.9570077045027497e-22, 1.446571598106951e-21, 5.3036659603136205e-21, 1.950216430503009e-20, 7.192321261240227e-20, 2.660394089390589e-19, 9.870222843275777e-19, 3.6732889425312376e-18, 1.3712957660787654e-17, 5.1354404241706654e-17, 1.929405327616833e-16, 7.272892273861344e-16, 2.7509355012894177e-15, 1.0442577810595348e-14, 1.0586983119551101e-14, 1.0748325133376516e-14, 1.0930578599474054e-14, 1.1138691697445431e-14, 9.259131434648513e-15, 7.260939358124996e-15, 5.7240708933028476e-15, 4.5504605705103284e-15]
    power_blk_expt  = 2.535160929202743e9
    imag_k2_expt    = [0.010364568326949035, 0.009920392683094205, 0.009476788118090197, 0.009033797870420055, 0.008591474242462661, 0.008149880848499443, 0.007709095592658267, 0.007269214676949498, 0.0068303580925931804]

    power_prf, power_blk, _, _, LNk = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg
    )
    imag_k2 = -imag.(LNk)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol, atol=atol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W"
    @info "Modelled total power = $(power_blk) W"
    @info "Expected imag(k2): $(imag_k2_expt)"
    @info "Modelled imag(k2): $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end

