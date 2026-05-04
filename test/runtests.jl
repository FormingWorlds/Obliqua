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

    power_prf_expt  = [0.0, 0.0, 3.3836853500831803e-15, 3.199107433705875e-15, 3.1066684894960716e-15, 3.0607994210208172e-15, 3.049518217180047e-15, 3.0708441754457593e-15, 6.484028375243808e-16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 48799.97241357575
    imag_k2_expt    = [6.68650424205084e-5]

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

    power_prf_expt  = [0.0, 0.0, 8.38538695455464e-19, 7.568669610705781e-19, 7.289216195147812e-19, 7.17739811066603e-19, 7.147892586185959e-19, 7.194534263837449e-19, 2.6571851057720503e-25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 11.123362128828322
    imag_k2_expt    = [1.5241075841179643e-8]

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
    # test solid1d module with andrade rheology
    @info " "
    @info "Testing solid1d-relax module with andrade rheology"

    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-relax"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 1.2500871386903165e-14, 1.3490093179591017e-14, 1.379776907343251e-14, 1.4367432888221434e-14, 1.5036981510041056e-14, 1.576511111912771e-14, 1.6552497289716098e-14, 1.8304993632021006e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 248561.9899812635
    imag_k2_expt    = [0.0003405761762193043]

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

    power_prf_expt  = [0.0, 3.061317898445062e-18, 3.31954900407594e-18, 3.2426735038389572e-18, 3.3496980944013153e-18, 3.505053284163295e-18, 3.674843274104151e-18, 3.85852228062806e-18, 7.519018418066546e-25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 57.809566146923785
    imag_k2_expt    = [7.92098622508629e-8]

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
    # test solid1d module with andrade rheology
    @info " "
    @info "Testing solid1d-mush module with andrade rheology"

    # update config to use only solid1d
    cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush"
    cfg["orbit"]["obliqua"]["module_fluid"] = "none"
    cfg["orbit"]["obliqua"]["module_mushy"] = "none"

    cfg["orbit"]["obliqua"]["material"]     = "andrade"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 4.950610209920039e-14, 4.567639539072659e-14, 4.2821578796970674e-14, 4.0271551813907144e-14, 3.787137048538955e-14, 3.560979144308022e-14, 1.2496829309351557e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 4.8052844153493536e8
    imag_k2_expt    = [0.6584133768599382]

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

    power_prf_expt  = [0.0, 0.0, 1.200063876496574e-17, 1.0577506483396584e-17, 9.837824669298708e-18, 9.24980671498098e-18, 8.698457949883253e-18, 8.179102199913946e-18, 1.2481596843673565e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 4.8114338657619107e8
    imag_k2_expt    = [0.6592559659893207]

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

    power_blk_expt  = 5.173115797760875e6
    imag_k2_expt    = [0.007088131204911426]

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

    power_prf_expt  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.038225140970252e-14, 9.038225140970252e-14, 9.038225140970252e-14, 9.038225140970252e-14, 7.352121265265731e-14, 5.620007993230152e-14, 4.295071065802233e-14, 3.282134829196775e-14]
    power_blk_expt  = 5.173115797760875e6
    imag_k2_expt    = [0.007088131204911426]

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

    power_prf_expt  = [0.0, 0.0, 1.4016031926830564e-12, 1.3252327158389514e-12, 1.2870162668501596e-12, 1.2680881334840961e-12, 1.2634868744830177e-12, 1.2723920449064774e-12, 2.700139333088046e-13, 2.700139333088046e-13, 9.135453140439972e-14, 3.077589624238481e-14, 1.032409751431994e-14, 3.448876655072261e-15, 1.1473875813510601e-15, 3.8016457270591504e-16, 1.2545374349729542e-16, 4.1235090092699916e-17, 1.3500251803670166e-17, 4.402790932153665e-18, 1.430357630822549e-18, 4.629244748327109e-19, 1.4925970913368033e-19, 4.794678235714497e-20, 1.534538342612161e-20, 4.893435781788792e-21, 1.5548335355809905e-21, 4.922698894532402e-22, 1.5530558012141815e-22, 4.882585539387812e-23, 1.5296984689356948e-23, 4.776048253436516e-24, 1.4861128774147332e-24, 4.608589397064881e-25, 1.4243931534608167e-25, 4.387827153552935e-26, 1.347220472644357e-26, 4.122956530401987e-27, 1.2576818499209302e-27, 3.8241549088487843e-28, 1.1590793200304743e-28, 3.501981611078859e-29, 1.0547445654520945e-29, 3.166816462678909e-30, 9.478730353074464e-31, 2.8284121594062167e-31, 8.415173815813368e-32, 2.4999722468286247e-32, 7.535872945365564e-33, 2.7103998397967852e-33, 2.461067331488008e-33, 6.40449392644815e-33, 2.1334331756124745e-32, 7.306476076697219e-32, 2.5150336009252767e-31, 8.68411965800176e-31, 3.007289565812955e-30, 1.0444320858282567e-29, 3.6377645206897917e-29, 1.2706705686237368e-28, 4.451129198391182e-28, 1.5636563114412864e-27, 5.508603317277516e-27, 1.9461047005488956e-26, 6.894630582425686e-26, 2.4494674298625543e-25, 8.726581495295039e-25, 3.1176298746233385e-24, 1.116886818577161e-23, 4.0123094730131237e-23, 1.4453629793639127e-22, 5.220996833425348e-22, 1.891127615394336e-21, 6.868726096099303e-21, 2.501593333124611e-20, 9.13565540811195e-20, 3.3453585303030626e-19, 1.2283530149639029e-18, 4.5225047381945086e-18, 1.6695830999220718e-17, 6.180280947734556e-17, 2.29391690004367e-16, 8.537184470396599e-16, 3.1857928013568695e-15, 1.1920232695749766e-14, 4.4721346239449414e-14, 1.6823123189920694e-13, 6.345397076604935e-13, 2.3997738290864037e-12, 9.099974217275072e-12, 3.459931207093036e-11, 3.459931207093036e-11, 3.459931207093036e-11, 3.459931207093036e-11, 3.459931207093036e-11, 2.8144721960568983e-11, 2.1514002847231398e-11, 1.64419993296216e-11, 1.2564368592020458e-11]
    power_blk_expt  = 1.5086896530891235e9
    imag_k2_expt    = [0.006235304368310926, 0.0059553773333793384, 0.005675317002289581, 0.005395132786467927, 0.005114835104873023, 0.004834435546425278, 0.004553947086848024, 0.004273384382550748, 0.00399276417566161]

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

