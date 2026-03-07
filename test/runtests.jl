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
    if isdefined(Obliqua, :run_solid1d_mush)
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
    cfg = Obliqua.open_config("$RES_DIR/config/default.toml")
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

    power_prf_expt  = [0.0, 0.0, 1.2930879671213278e-12, 1.2219771504808403e-12, 1.1736558859289374e-12, 1.1314931658945296e-12, 1.0916026417910732e-12, 1.0538251910982026e-12, 1.5883496284658382e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 1.731021692653176e7
    imag_k2_expt    = [0.0237182181024913]

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

    power_prf_expt  = [0.0, 0.0, 3.3410413059922125e-16, 3.0141311725796533e-16, 2.871071975586501e-16, 2.7665275365727065e-16, 2.6681804901946297e-16, 2.5750291702649744e-16, 7.077614201094658e-25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 4280.179309796564
    imag_k2_expt    = [5.8646420676523305e-6]

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

    power_prf_expt  = [0.0, 0.0, 3.9458586892328636e-15, 3.684375027447792e-15, 3.5332688538838223e-15, 3.4385147729212697e-15, 3.3857777200941584e-15, 3.3721836008236663e-15, 1.2678042929823814e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 61564.69369317412
    imag_k2_expt    = [8.435508570604252e-5]

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
    cfg["orbit"]["obliqua"]["material"] = "maxwell"

    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full("$RES_DIR/interior_data/runtests_mantle.json", false)

    power_prf_expt  = [0.0, 0.0, 9.412009210814184e-19, 8.391526074045509e-19, 7.984510236608038e-19, 7.771471545020786e-19, 7.656366535345375e-19, 7.630908981500375e-19, 1.1998652842288277e-15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    power_blk_expt  = 3139.7068234228814
    imag_k2_expt    = [4.3019825535337385e-6]

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

    power_prf_expt  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9150692840613072e-13, 2.32176469195039e-13, 2.2984543955634825e-13, 2.275218082497258e-13, 1.8312731634158356e-13, 1.3840633350642642e-13, 1.0446190925830448e-13, 7.863674415674511e-14] 
    power_blk_expt  = 1.2647997272570137e7
    imag_k2_expt    = [0.017330109677062184]

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

    power_prf_expt  = [0.0, 0.0, 5.348438930694959e-10, 5.054382700817955e-10, 4.854545590088724e-10, 4.6801736056944e-10, 4.5151996680864546e-10, 4.3589657868479087e-10, 6.567896224392372e-13, 6.567896224392372e-13, 2.2221337785775897e-13, 7.486017119771332e-14, 2.511263039415512e-14, 8.389146324289135e-15, 2.7909384049643753e-15, 9.247231915500194e-16, 3.0515727768368006e-16, 1.0030141378763135e-16, 3.283839902745504e-17, 1.0709474724406045e-17, 3.479242855318201e-18, 1.1260307470708736e-18, 3.6306359011179086e-19, 1.1662712622133663e-19, 3.7326550015850954e-20, 1.1902933304838137e-20, 3.7820216100559424e-21, 1.1974113738139303e-21, 3.777697398081894e-22, 1.187654086455705e-22, 3.720882354278839e-23, 1.1617396519825671e-23, 3.614863513517372e-24, 1.1210064802892951e-24, 3.464734689808167e-25, 1.0673076400893613e-25, 3.2770176513763044e-26, 1.0028797513345896e-26, 3.059221341787089e-27, 9.301983893760305e-28, 2.8193777248554025e-28, 8.518320338741817e-29, 2.5655908770677535e-29, 7.703054921027464e-30, 2.305633125385449e-30, 6.879897808388008e-31, 2.0468875703213448e-31, 6.079540588563313e-32, 1.82809769862545e-32, 6.424949729505258e-33, 5.415085000960695e-33, 1.3628910622963522e-32, 4.522144596189985e-32, 1.5481839048545105e-31, 5.328997698087551e-31, 1.8400365756743336e-30, 6.372000678813409e-30, 2.2129966861902714e-29, 7.707883477528396e-29, 2.692362472722205e-28, 9.431282592662013e-28, 3.313155806917663e-27, 1.167191340904073e-26, 4.1235072197865054e-26, 1.4608699612294494e-25, 5.190058185883491e-25, 1.849033188695096e-24, 6.605795306391037e-24, 2.3665175151102258e-23, 8.501488679062831e-23, 3.062509781173761e-22, 1.1062517926728058e-21, 4.007018930387838e-21, 1.4553811847848123e-20, 5.3005052437310766e-20, 1.935709883553861e-19, 7.088318551714252e-19, 2.602697852906662e-18, 9.582516775281375e-18, 3.5375989609502625e-17, 1.30950986866763e-16, 4.860469684006445e-16, 1.808902768201557e-15, 6.750222438409277e-15, 2.5257205107513347e-14, 9.475789973937655e-14, 3.5645702881982644e-13, 1.34449552742012e-12, 5.084764816251399e-12, 1.928149568428732e-11, 7.331081060740859e-11, 7.331081060740859e-11, 8.88794915973433e-11, 8.798711294233438e-11, 8.709756474695825e-11, 7.010282387306757e-11, 5.298315256901038e-11, 3.998889011221784e-11, 3.010277455399582e-11]
    power_blk_expt  = 1.2866082902183443e10
    imag_k2_expt    = [0.046417892444006494, 0.045639480552623346, 0.04487993352593461, 0.0441413835094393, 0.04342634229137915, 0.04273779394964199, 0.04207931712777366, 0.041455248954265535, 0.040870908664398777]

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

