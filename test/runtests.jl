#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get Love root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras
using .MAIN

@info "Begin Love tests"

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

# -------------
# Test module imported
# -------------
if suite >= 0
    @info " "
    @info "Testing module imported"
    if isdefined(Love, :calc_lovepy_tides)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"
end

if suite > 2
    # -------------
    # Test interior data validity
    # -------------
    @info " "
    @info "Testing interior data validity"
    ok = load.load_interior("$RES_DIR/interior_data/test_mantle.json", true)
    if ok
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Test Love.jl method (andrade)
    # -------------
    @info " "
    @info "Testing Love.jl method (Andrade)"
    omega, ecc, rho, radius, visc, shear, bulk, ncalc = load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

    power_prf_expt  = [0.0, 3.777310932533753e-12, 3.650761379774313e-12, 3.5304847221781085e-12, 3.421375892547284e-12, 3.314966022202826e-12, 3.2185792279633594e-12, 3.1227775210056437e-12, 3.0387623792987696e-12, 2.9503173602715362e-12, 2.8854427562467733e-12, 2.840183803388453e-12, 2.766937542627038e-12, 2.6614024308196217e-12, 2.5882592611305263e-12, 2.5270407684882074e-12, 2.4696571292730092e-12, 2.414940623392922e-12, 2.3625490809018268e-12, 2.312188588780301e-12, 2.263588422530672e-12, 2.2165597967690805e-12, 2.170959666893054e-12, 2.1266604153692343e-12, 2.0835487124355886e-12, 2.0415233385939996e-12, 2.0004922650117834e-12, 1.9603731740913513e-12, 1.921091050223868e-12, 1.8825812613342155e-12, 1.8447801687955216e-12, 1.8076327139644142e-12, 1.771090343999002e-12, 1.73510794435057e-12, 1.6996475025775436e-12, 1.6646747616218843e-12, 1.6301609527692688e-12, 1.5960809957264244e-12, 1.562410030784291e-12, 1.529128970618847e-12, 1.4962208528868597e-12, 1.4636715624790206e-12, 1.4314699344379267e-12, 1.399606576523371e-12, 1.3680752968065698e-12, 1.3368713097274294e-12, 1.3059922477773532e-12, 1.2754390640359429e-12, 1.245212794684355e-12, 1.215315597304841e-12, 1.1857521079276538e-12, 1.1565282363027113e-12, 1.1276510872890982e-12, 1.09912987419273e-12, 1.0709751376920464e-12, 1.0431972650446318e-12, 1.0158101431173686e-12, 9.888284234763632e-13, 9.622670614333102e-13, 9.361420239769226e-13, 9.104710233397912e-13, 8.85273996222416e-13, 8.605702022528852e-13, 8.363799038147769e-13, 8.127249183563586e-13, 7.896303606558283e-13, 7.671229357856957e-13, 7.452278144406235e-13, 7.239703594971455e-13, 7.033797797374361e-13, 6.83485154464037e-13, 6.6431570746088e-13, 6.459045785582216e-13, 6.282835604446689e-13, 6.114879986487898e-13, 5.955527886842414e-13, 5.805150639198377e-13, 5.664147311104618e-13, 5.532905562807219e-13, 5.411779498624159e-13, 5.301191570564058e-13, 5.201543628263393e-13, 5.113322254080463e-13, 5.036997676692369e-13, 4.973082269350431e-13, 4.922088792777169e-13, 4.884482322090376e-13, 4.860823905034184e-13, 4.851463181940493e-13, 4.856425352418459e-13, 4.875534348367834e-13, 4.909133845741023e-13, 4.957926420786741e-13, 5.022726796825203e-13, 5.104530461657479e-13, 5.204340134922985e-13, 5.322983847200073e-13, 5.461978547218397e-13]
    power_blk_expt  = 4.66556078444765e12
    imag_k2_expt    = -0.00464899461127863

    power_prf, power_blk, imag_k2 = Love.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="andrade")
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= isapprox(imag_k2, imag_k2_expt; rtol=rtol)

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W, Modelled total power = $(power_blk) W"
    @info "Expected imag(k2) = $(imag_k2_expt), Modelled imag(k2) = $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Test Love.jl method (andrade)
    # -------------
    @info " "
    @info "Testing Love.jl method (Maxwell)"
    omega, ecc, rho, radius, visc, shear, bulk, ncalc = load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

    power_prf_expt  = [0.0, 2.0653528984529426e-15, 2.0060790864097053e-15, 1.9455449084871923e-15, 1.88919107873945e-15, 1.8318816669547857e-15, 1.780822503763121e-15, 1.728608390258191e-15, 1.6843719042431002e-15, 1.6351431326184818e-15, 1.602588165496073e-15, 1.5740111906069696e-15, 1.495297618306528e-15, 1.388298500401266e-15, 1.337328046865521e-15, 1.3040429771846793e-15, 1.2741497526570228e-15, 1.2458148771697857e-15, 1.2187084469823691e-15, 1.1926556532177095e-15, 1.1675124869852474e-15, 1.143180685425577e-15, 1.1195865042934052e-15, 1.0966643586968862e-15, 1.0743560344968784e-15, 1.0526095330595604e-15, 1.0313775626005167e-15, 1.0106178139796856e-15, 9.902917151200325e-16, 9.703660203447969e-16, 9.508079740039851e-16, 9.315892209489253e-16, 9.126847453382941e-16, 8.940712885196417e-16, 8.757292387987386e-16, 8.576409122281058e-16, 8.397914437438691e-16, 8.22167866668446e-16, 8.04757320761093e-16, 7.8754990844539e-16, 7.705368206938009e-16, 7.537107064711539e-16, 7.370657330152177e-16, 7.205969767324174e-16, 7.04301157076733e-16, 6.881757193609651e-16, 6.722193533200811e-16, 6.564324643239135e-16, 6.408154970766184e-16, 6.25369463766325e-16, 6.100966519303122e-16, 5.950000015799526e-16, 5.800830707003278e-16, 5.653505023140044e-16, 5.508076190637922e-16, 5.364596695648551e-16, 5.223137023250271e-16, 5.083771613673487e-16, 4.946576443258068e-16, 4.811632642609772e-16, 4.679030370789002e-16, 4.548871197521228e-16, 4.421253222706185e-16, 4.296279704137449e-16, 4.174061848707881e-16, 4.0547278632451565e-16, 3.9384142025168055e-16, 3.825249404792695e-16, 3.7153629221952967e-16, 3.6089039475526676e-16, 3.506021069676788e-16, 3.4068636041900957e-16, 3.311600954986214e-16, 3.220395663631613e-16, 3.133428459718256e-16, 3.050877805179145e-16, 2.97293364742345e-16, 2.8998001324622132e-16, 2.8316755956632255e-16, 2.7687411796711393e-16, 2.711213134800712e-16, 2.6592970961588713e-16, 2.6132421052082086e-16, 2.573288923025315e-16, 2.539699825618067e-16, 2.5127373283597636e-16, 2.4926395967002206e-16, 2.4796938759397796e-16, 2.474079562866971e-16, 2.47581059020701e-16, 2.4847993378891895e-16, 2.5012239906355857e-16, 2.5254446170072026e-16, 2.557878445692455e-16, 2.599034824949647e-16, 2.64942663593431e-16, 2.709477727695381e-16, 2.779963783408105e-16]
    power_blk_expt  = 2.430386157481513e9
    imag_k2_expt    = -2.4217564986232166e-6

    power_prf, power_blk, imag_k2 = Love.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="maxwell")
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= isapprox(imag_k2, imag_k2_expt; rtol=rtol)

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W, Modelled total power = $(power_blk) W"
    @info "Expected imag(k2) = $(imag_k2_expt), Modelled imag(k2) = $(imag_k2)"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end