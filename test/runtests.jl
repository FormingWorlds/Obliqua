#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get Obliqua root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"

# Include libraries
using LoggingExtras
using Obliqua

@info "Begin Obliqua tests"

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

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 7.12546350513748e-15, 7.141653473320296e-15, 7.266955578280898e-15, 7.413341032856918e-15, 7.5585214768856e-15, 7.702021368425093e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 107761.54726903763
    imag_k2_expt    = [0.00014765337095640072]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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
    cfg["orbit"]["obliqua"]["material"]     = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 1.722694782808758e-18, 1.6491571129537458e-18, 1.6649092230686594e-18, 1.6982534724631658e-18, 1.7317130709612992e-18, 1.764826259885598e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 24.92921722417735
    imag_k2_expt    = [3.415766617812618e-8]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 7.750115499856373e-15, 8.057499722716283e-15, 8.059575850362661e-15, 8.185436835049663e-15, 8.33553548214869e-15, 8.484799440335237e-15, 8.632786806517406e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 137743.23201529591
    imag_k2_expt    = [0.0001887338577527246]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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
    cfg["orbit"]["obliqua"]["material"]     = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 1.8887574959136746e-18, 1.972957870150714e-18, 1.884584844212067e-18, 1.8986419071321904e-18, 1.9329191192829336e-18, 1.967464038487986e-18, 2.0017743960506084e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 32.41911211996522
    imag_k2_expt    = [4.44202158305658e-8]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    # lower visc_sus to include mush
    visc_sus = cfg["orbit"]["obliqua"]["visc_sus"] = 5e10

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 9.86707030745882e-15, 9.150910190756036e-15, 8.622765164513395e-15, 8.150013792655512e-15, 7.702026907818404e-15, 7.276882015209697e-15, 2.36961014227723e-14, 1.9566536667394827e-14, 1.6121402444527615e-14, 7.538539756442404e-15, 8.508094326483907e-15, 3.568025886298421e-15, 3.582082099457579e-15, 1.8450607517355677e-15, 1.3451995999757726e-15, 8.725598622614597e-16, 4.877653976477049e-16, 3.636754700375336e-16, 1.762407992723698e-16, 1.343339281131328e-16, 6.412353878778215e-17, 4.3760705656795036e-17, 2.3569631734565316e-17, 1.2275542502234671e-17, 8.578724675042095e-18, 2.996216880865113e-18, 2.818774227616443e-18, 7.725592763205111e-19, 7.291475505535899e-19, 4.578348846472072e-19, 1.7918454344299687e-18, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 3.0395423054335793e6
    imag_k2_expt    = [0.004164738526270287]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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
    cfg["orbit"]["obliqua"]["material"] = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 2.811270711855333e-18, 2.4873681738050708e-18, 2.322121690521703e-18, 2.191368048585138e-18, 2.0681585759977523e-18, 1.9514616381743257e-18, 9.401030519482231e-17, 1.3496470773634048e-16, 6.582027976570243e-17, 7.489744673034899e-17, 6.156680255962594e-17, 3.460102640163939e-17, 4.550106850237844e-17, 1.8331157658567176e-17, 2.6909888276869764e-17, 1.138847035948234e-17, 1.405393299838946e-17, 6.943151792071354e-18, 6.977590742356405e-18, 3.80235936859e-18, 3.4214575190463734e-18, 1.8384897659114318e-18, 1.6702466275884794e-18, 7.837030357462777e-19, 7.983327480776783e-19, 3.0102026582678284e-19, 3.561792052557227e-19, 1.152470964370035e-19, 1.3751179251924843e-19, 9.509403483402107e-20, 9.061242335085774e-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 76721.32337867716
    imag_k2_expt    = [0.00010512248857021226]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_blk_expt  = 5.167672882253202e6
    imag_k2_expt    = [0.007080673397902285]

    _, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

    power_prf_expt  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.028715534438163e-14, 9.028715534438163e-14, 9.028715534438163e-14, 9.028715534438163e-14, 7.344385700006137e-14, 5.6140948782225885e-14, 4.2905519887460544e-14, 3.2786815172550926e-14]
    power_blk_expt  = 5.167672882253202e6
    imag_k2_expt    = [0.007080673397902285]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

    power_prf_expt  = [0.0, 0.0, 2.9507745130337324e-12, 2.957492544894895e-12, 3.0093798848852802e-12, 3.0699942838630136e-12, 3.1301087443252187e-12, 3.189526446012955e-12, 3.189526446012955e-12, 1.0963618378318885e-12, 3.7524391601966574e-13, 1.278880490634037e-13, 4.340368867102278e-14, 1.4669889562954858e-14, 4.938022630587577e-15, 1.6554924192323772e-15, 5.528030698211065e-16, 1.8386685971776486e-16, 6.091808976258348e-17, 2.0105661404428437e-17, 6.610565688410192e-18, 2.1653365808715678e-18, 7.066372826281614e-19, 2.2975730444184624e-19, 7.443214658305481e-20, 2.402627776885429e-20, 7.727926641765474e-21, 2.4768823612742404e-21, 7.910948363947936e-22, 2.5179490667963084e-22, 7.986832225390398e-23, 2.524788461375205e-23, 7.95447314090666e-24, 2.497736398632058e-24, 7.81705067304542e-25, 2.43844179976677e-25, 7.581700680138862e-26, 2.3497243537257714e-26, 7.258955978439079e-27, 2.2353675605634917e-27, 6.86201250668302e-28, 2.099866879583562e-28, 6.405887712936706e-29, 1.948154896774264e-29, 5.906543012799335e-30, 1.7853323147817247e-30, 5.380272249605189e-31, 1.6172101556799416e-31, 4.8696574173948217e-32, 1.5395902009428098e-32, 7.449587884676381e-33, 1.1867910509433857e-32, 3.615139994984785e-32, 1.213099710206271e-31, 4.11842994772204e-31, 1.4033105735719724e-30, 4.795701640364832e-30, 1.6435941328164814e-29, 5.649022188495877e-29, 1.947075293104084e-28, 6.730043193775503e-28, 2.3327794465513773e-27, 8.108593064357973e-27, 2.8263685143966506e-26, 9.879140502313731e-26, 3.4626807127680235e-25, 1.217039381153432e-24, 4.289353781166347e-24, 1.5158997325577025e-23, 5.372011529522199e-23, 1.9089206788068747e-22, 6.801741207805623e-22, 2.4301357661387006e-21, 8.705951947519668e-21, 3.1273357857633184e-20, 1.1264259524532833e-19, 4.068161322472987e-19, 1.473189801198765e-18, 5.3491160810343616e-18, 1.9474486452782743e-17, 7.10900292213008e-17, 2.6020054102308763e-16, 9.549107380509521e-16, 3.513749608145839e-15, 1.2963753981674927e-14, 4.795581609850411e-14, 1.7786931534026058e-13, 6.614679959444938e-13, 2.46640058030969e-12, 9.220710510053815e-12, 3.4562908259458124e-11, 3.4562908259458124e-11, 3.4562908259458124e-11, 3.4562908259458124e-11, 3.4562908259458124e-11, 2.8115109373183124e-11, 2.149136679880193e-11, 1.6424699810989694e-11, 1.2551148938850651e-11]
    power_blk_expt  = 2.4310813508705845e9
    imag_k2_expt    = [0.01000287153467124, 0.009561732793922827, 0.00912095768370095, 0.008680565130249904, 0.00824057872897682, 0.007801027901476037, 0.007361949428107609, 0.006923389510482543, 0.0064854065968675805]

    power_prf, power_blk, _, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
    )
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

