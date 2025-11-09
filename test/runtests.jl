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

rtol   = 1e-3

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

    power_prf_expt  = [0.0, 1.4060891744290011e7, 1.3849717515140837e7, 1.364797510722255e7, 1.3475999993537705e7, 1.3301981632862581e7, 1.3156186861130793e7, 1.3001328722353794e7, 1.2884790059360832e7, 1.273905042057936e7, 1.2685949878706416e7, 1.2713184988180216e7, 1.2608417854282262e7, 1.2344734487517366e7, 1.2219296702963324e7, 1.2141597227024635e7, 1.207491053915965e7, 1.2014247342453409e7, 1.195841386645607e7, 1.1906359581584478e7, 1.1857063933397794e7, 1.1809856860670978e7, 1.1764257158703893e7, 1.1719824570596479e7, 1.1676159155228265e7, 1.1632894324265566e7, 1.1589684549458914e7, 1.1546211826021804e7, 1.1502175035241676e7, 1.1457311446926601e7, 1.1411343030359829e7, 1.1364023970939187e7, 1.1315131119585259e7, 1.1264446847192898e7, 1.121178448482343e7, 1.1156969153101485e7, 1.1099851153132072e7, 1.1040296104827875e7, 1.0978162279651169e7, 1.0913340068922155e7, 1.0845728492311005e7, 1.0775241504805574e7, 1.0701810205575474e7, 1.062537544949735e7, 1.0545899768197304e7, 1.0463355121885581e7, 1.0377731720234249e7, 1.0289046723172579e7, 1.0197319747174637e7, 1.0102581582987234e7, 1.0004886382800227e7, 9.904302960432561e6, 9.800915008918457e6, 9.694830034764985e6, 9.586173867647875e6, 9.475078243924456e6, 9.361714550421234e6, 9.246270789305018e6, 9.128947960539702e6, 9.009967334737463e6, 8.889578578275379e6, 8.76806600569353e6, 8.64572132033498e6, 8.522860594900785e6, 8.399831052729681e6, 8.277030838074323e6, 8.154893712136215e6, 8.03385701634784e6, 7.914379874078655e6, 7.796984312773714e6, 7.682214042134573e6, 7.570638031390361e6, 7.462894764493932e6, 7.359634291378405e6, 7.261575981392939e6, 7.169463873772936e6, 7.084099820068225e6, 7.00635268941591e6, 6.937112396922527e6, 6.877261619846973e6, 6.827804989491275e6, 6.789759695871346e6, 6.764293942572604e6, 6.752599837205235e6, 6.755972465131892e6, 6.775757277197395e6, 6.813284225306799e6, 6.8700675713335695e6, 6.947375783479572e6, 7.046048926297024e6, 7.166649736254762e6, 7.3105182395864455e6, 7.479579539229824e6, 7.676002914951468e6, 7.902320859967681e6, 8.161177027228087e6, 8.455040357920263e6, 8.787579005613863e6]
    power_blk_expt  = 4.66556078444765e12
    imag_k2_expt    = -0.00464899461127863

    power_prf, power_blk, imag_k2 = Love.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="andrade")
    test_pass = true
    if maximum(abs.((power_prf_expt .- power_prf) ./ power_prf_expt)) > rtol
        global test_pass = false
    end
    if maximum(abs.((power_blk_expt .- power_blk) ./ power_blk_expt)) > rtol
        global test_pass = false
    end
    if abs(imag_k2_expt - imag_k2)/imag_k2_expt > rtol  
        global test_pass = false
    end
    @info "Expected values = $(power_blk_expt) W"
    @info "Modelled values = $(power_blk) W"
    @info "Expected values = $(imag_k2_expt)"
    @info "Modelled values = $(imag_k2)"
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

    power_prf_expt  = [0.0, 7688.195578237189, 7610.365904916004, 7520.992573601966, 7441.082281176289, 7350.801305301731, 7279.247143075098, 7196.864191345621, 7141.979291240992, 7060.315119624263, 7045.834617803442, 7045.563258150643, 6813.7918180685265, 6439.528053014826, 6313.589848596476, 6265.495969677391, 6229.708248877672, 6197.886197759412, 6168.684141150872, 6141.448051523617, 6115.629984371139, 6090.879796003774, 6066.949455368788, 6043.613182114618, 6020.666088959543, 5997.920410534611, 5975.199134990148, 5952.339366910054, 5929.186812452183, 5905.606856938667, 5881.456902397063, 5856.610849831656, 5830.9542165993, 5804.375285637657, 5776.778219356933, 5748.073210264595, 5718.1833646330415, 5687.039664055059, 5654.569248551796, 5620.7157173443375, 5585.427055737288, 5548.658860764018, 5510.375527686263, 5470.546422797349, 5429.152009420022, 5386.17758852936, 5341.617821944105, 5295.481263381207, 5247.777672874879, 5198.522448634274, 5147.742971944214, 5095.474103999869, 5041.7583465072985, 4986.650413972144, 4930.21437604529, 4872.5173404660645, 4813.646764912801, 4753.698629274035, 4692.775529526792, 4630.990393238254, 4568.47075266402, 4505.3618810008165, 4441.812811754854, 4377.985059448424, 4314.056054879249, 4250.2294026744075, 4186.727041012242, 4123.772667193398, 4061.6010414892053, 4000.479168678588, 3940.6852061853065, 3882.5102178995476, 3826.280890421742, 3772.329862347791, 3721.025251910288, 3672.7478427931537, 3627.9080450609595, 3586.9505746996165, 3550.330954761336, 3518.5009830357903, 3491.975039038115, 3471.273847751674, 3456.995481711091, 3449.7503760503823, 3450.2016108386933, 3459.0381112696136, 3476.9409624761506, 3504.6855435666994, 3542.9219540999343, 3592.0819284382574, 3652.457493677475, 3724.738398845855, 3809.9113104930275, 3909.087559163853, 4023.563794060813, 4154.693299623679, 4303.740530568818, 4472.582509437505]
    power_blk_expt  = 2.430386157481513e9
    imag_k2_expt    = -2.4217564986232166e-6

    power_prf, power_blk, imag_k2 = Love.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="maxwell")
    test_pass = true
    if maximum(abs.((power_prf_expt .- power_prf) ./ power_prf_expt)) > rtol
        global test_pass = false
    end
    if maximum(abs.((power_blk_expt .- power_blk) ./ power_blk_expt)) > rtol
        global test_pass = false
    end
    if abs(imag_k2_expt - imag_k2)/imag_k2_expt > rtol  
        global test_pass = false
    end
    @info "Expected values = $(power_blk_expt) W"
    @info "Modelled values = $(power_blk) W"
    @info "Expected values = $(imag_k2_expt)"
    @info "Modelled values = $(imag_k2)"
    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end