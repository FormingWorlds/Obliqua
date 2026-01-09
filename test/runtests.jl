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

# -------------
# Test module imported
# -------------
if suite >= 0
    @info " "
    @info "Testing if Love.jl module is imported"
    if isdefined(Love, :calc_lovepy_tides)
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"
end

if suite >= 0
    @info " "
    @info "Testing if Fluid.jl module is imported"
    if isdefined(Fluid, :compute_fluid_lovenumbers)
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

    power_prf_expt  = [0.0, 3.7813635828309455e-12, 3.6546971930744196e-12, 3.534301441625492e-12, 3.4250818342911407e-12, 3.3185594391199533e-12, 3.222072390932556e-12, 3.1261682487295903e-12, 3.042066294779961e-12, 2.953524785677839e-12, 2.8885862844037246e-12, 2.8432714801848047e-12, 2.7698721963170804e-12, 2.664128855429253e-12, 2.5908861220399643e-12, 2.5296024211532913e-12, 2.472160194371673e-12, 2.417388159751837e-12, 2.3649435031821515e-12, 2.3145319679554708e-12, 2.265882545884234e-12, 2.2188062577659135e-12, 2.1731599134276475e-12, 2.128815765885817e-12, 2.085660370496327e-12, 2.0435924051699026e-12, 2.0025197478010034e-12, 1.9623599973725486e-12, 1.9230380622550588e-12, 1.8844892448709153e-12, 1.8466498420949624e-12, 1.8094647394735012e-12, 1.7728853349520439e-12, 1.7368664682522491e-12, 1.701370088399369e-12, 1.6663619036147523e-12, 1.6318131160266646e-12, 1.5976986199348546e-12, 1.5639935304462326e-12, 1.5306787408934798e-12, 1.4977372717420818e-12, 1.465154993569144e-12, 1.432920730098821e-12, 1.4010250795679617e-12, 1.3694618437640921e-12, 1.3382262322764747e-12, 1.3073158751992017e-12, 1.2767317265953066e-12, 1.2464748237094159e-12, 1.216547326306679e-12, 1.1869538751132241e-12, 1.157700385861058e-12, 1.1287939706091908e-12, 1.1002438520014154e-12, 1.0720605813934511e-12, 1.044254556574055e-12, 1.0168396784880594e-12, 9.898306135590174e-13, 9.632423322585955e-13, 9.370908177510246e-13, 9.113938002236352e-13, 8.861712365813649e-13, 8.614424059717804e-13, 8.37227591314081e-13, 8.135486321591494e-13, 7.904306687411948e-13, 7.679004332163209e-13, 7.459831217998293e-13, 7.247041230676753e-13, 7.04092675397919e-13, 6.841778875685489e-13, 6.649890129804959e-13, 6.4655922504855e-13, 6.289203486878291e-13, 6.121077652469177e-13, 5.961564056090266e-13, 5.811034407815793e-13, 5.669888179656782e-13, 5.538513424796234e-13, 5.417264606378834e-13, 5.306564604686435e-13, 5.206815676527048e-13, 5.118504897392851e-13, 5.042102972640054e-13, 4.978122794822511e-13, 4.927077644394243e-13, 4.889433067766877e-13, 4.86575068131097e-13, 4.856380479557968e-13, 4.861347686945705e-13, 4.880476057383465e-13, 4.914109616588214e-13, 4.962951652226687e-13, 5.027817712935375e-13, 5.109704295444634e-13, 5.209615135688573e-13, 5.32837910348075e-13, 5.467514687507863e-13]
    power_blk_expt  = 4.670347977728458e12
    imag_k2_expt    = -0.004653764806501424

    power_prf, power_blk, imag_k2 = Obliqua.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="andrade")
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
    # Test Love.jl method (maxwell)
    # -------------
    @info " "
    @info "Testing Love.jl method (Maxwell)"
    omega, ecc, rho, radius, visc, shear, bulk, ncalc = load.load_interior("$RES_DIR/interior_data/test_mantle.json", false)

    power_prf_expt  = [0.0, 2.0653528984529426e-15, 2.0060790864097053e-15, 1.9455449084871923e-15, 1.88919107873945e-15, 1.8318816669547857e-15, 1.780822503763121e-15, 1.728608390258191e-15, 1.6843719042431002e-15, 1.6351431326184818e-15, 1.602588165496073e-15, 1.5740111906069696e-15, 1.495297618306528e-15, 1.388298500401266e-15, 1.337328046865521e-15, 1.3040429771846793e-15, 1.2741497526570228e-15, 1.2458148771697857e-15, 1.2187084469823691e-15, 1.1926556532177095e-15, 1.1675124869852474e-15, 1.143180685425577e-15, 1.1195865042934052e-15, 1.0966643586968862e-15, 1.0743560344968784e-15, 1.0526095330595604e-15, 1.0313775626005167e-15, 1.0106178139796856e-15, 9.902917151200325e-16, 9.703660203447969e-16, 9.508079740039851e-16, 9.315892209489253e-16, 9.126847453382941e-16, 8.940712885196417e-16, 8.757292387987386e-16, 8.576409122281058e-16, 8.397914437438691e-16, 8.22167866668446e-16, 8.04757320761093e-16, 7.8754990844539e-16, 7.705368206938009e-16, 7.537107064711539e-16, 7.370657330152177e-16, 7.205969767324174e-16, 7.04301157076733e-16, 6.881757193609651e-16, 6.722193533200811e-16, 6.564324643239135e-16, 6.408154970766184e-16, 6.25369463766325e-16, 6.100966519303122e-16, 5.950000015799526e-16, 5.800830707003278e-16, 5.653505023140044e-16, 5.508076190637922e-16, 5.364596695648551e-16, 5.223137023250271e-16, 5.083771613673487e-16, 4.946576443258068e-16, 4.811632642609772e-16, 4.679030370789002e-16, 4.548871197521228e-16, 4.421253222706185e-16, 4.296279704137449e-16, 4.174061848707881e-16, 4.0547278632451565e-16, 3.9384142025168055e-16, 3.825249404792695e-16, 3.7153629221952967e-16, 3.6089039475526676e-16, 3.506021069676788e-16, 3.4068636041900957e-16, 3.311600954986214e-16, 3.220395663631613e-16, 3.133428459718256e-16, 3.050877805179145e-16, 2.97293364742345e-16, 2.8998001324622132e-16, 2.8316755956632255e-16, 2.7687411796711393e-16, 2.711213134800712e-16, 2.6592970961588713e-16, 2.6132421052082086e-16, 2.573288923025315e-16, 2.539699825618067e-16, 2.5127373283597636e-16, 2.4926395967002206e-16, 2.4796938759397796e-16, 2.474079562866971e-16, 2.47581059020701e-16, 2.4847993378891895e-16, 2.5012239906355857e-16, 2.5254446170072026e-16, 2.557878445692455e-16, 2.599034824949647e-16, 2.64942663593431e-16, 2.709477727695381e-16, 2.779963783408105e-16]
    power_blk_expt  = 2.430386157481513e9
    imag_k2_expt    = -2.4217564986232166e-6

    power_prf, power_blk, imag_k2 = Obliqua.calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=ncalc, material="maxwell")
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

if suite > 4
    # -------------
    # Test mush interior data validity
    # -------------
    @info " "
    @info "Testing mush interior data validity"
    ok = load.load_interior("$RES_DIR/interior_data/test_mantle_mush.json", true)
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
    omega, ecc, rho, radius, visc, shear, bulk, phi, ncalc = load.load_interior_mush("$RES_DIR/interior_data/test_mantle_mush.json", false)

    power_prf_expt  = [0.0, 3.710337246759876e-12, 3.577404394688112e-12, 3.4685705186794216e-12, 3.3693962943911173e-12, 3.2716924765504224e-12, 3.182806745136962e-12, 3.0935345162755547e-12, 3.0150753767033445e-12, 2.9314509285268902e-12, 2.8705845432824803e-12, 2.8287021321583345e-12, 2.758423758026706e-12, 2.655442217127958e-12, 2.584421782042107e-12, 2.5249935060740132e-12, 2.469124062722172e-12, 2.4156938816227557e-12, 2.3644038336958e-12, 2.314997930631948e-12, 2.267238510187375e-12, 2.2209654699768102e-12, 2.176060439401945e-12, 2.1324167383534012e-12, 2.089938516103216e-12, 2.048538832067104e-12, 2.00813697381308e-12, 1.9686592106176683e-12, 1.930036592860489e-12, 1.8922082528447515e-12, 1.8551121618852123e-12, 1.8186929297246463e-12, 1.7828998916419356e-12, 1.7476841794369844e-12, 1.7130025572673305e-12, 1.6788141925750806e-12, 1.645082533265391e-12, 1.611773617046761e-12, 1.5788526799441066e-12, 1.5462898649328185e-12, 1.5140566390067991e-12, 1.4821266133861285e-12, 1.4504757340237495e-12, 1.4190811685894463e-12, 1.387922827326786e-12, 1.356981613995432e-12, 1.326240515337058e-12, 1.2956855805083399e-12, 1.2653026884261634e-12, 1.2350786489444604e-12, 1.2050026294846258e-12, 1.175064975220821e-12, 1.1452571686911889e-12, 1.1155727954644753e-12, 1.0860067849474389e-12, 1.0565539376870656e-12, 1.0272126604591608e-12, 9.979822183476823e-13, 9.688622930379236e-13, 9.398537222839842e-13, 9.109592585035203e-13, 8.821840608862061e-13, 8.535328050498925e-13, 8.250113880784356e-13, 7.966274910226316e-13, 7.683922896689487e-13, 7.403187493099858e-13, 7.124185996765918e-13, 6.847040684460189e-13, 6.571913706254313e-13, 6.298969025892628e-13, 6.028375543173759e-13, 5.760341718829457e-13, 5.495067375566012e-13, 5.232787847753632e-13, 4.973738373541958e-13, 4.718177612461197e-13, 4.466391051969356e-13, 4.2186590738485384e-13, 3.9752430473925056e-13, 3.7364644510771696e-13, 3.5026375018812215e-13, 3.274142960340833e-13, 3.051355734081173e-13, 2.834683487024367e-13, 2.6245396097911674e-13, 2.421321259319098e-13, 2.2254878048074735e-13, 2.0374153122975788e-13, 1.8573915721070098e-13, 1.6857180781007828e-13, 1.5229740999820876e-13, 1.3698840530044252e-13, 1.2271978972776673e-13, 1.095690483248367e-13, 9.761068223428319e-14, 8.691375615159178e-14, 7.755705270510724e-14]
    power_blk_expt  = 4.673775086943557e12
    imag_k2_expt    = -0.004657179746957536

    power_prf, power_blk, imag_k2 = Obliqua.calc_lovepy_tides_mush(omega, ecc, rho, radius, visc, shear, bulk, phi; ncalc=ncalc, material="andrade")
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
    # Test Love.jl method (maxwell)
    # -------------
    @info " "
    @info "Testing Love.jl method (Maxwell)"
    omega, ecc, rho, radius, visc, shear, bulk, phi, ncalc = load.load_interior_mush("$RES_DIR/interior_data/test_mantle_mush.json", false)

    power_prf_expt  = [0.0, 2.586681737817711e-15, 1.9648435866065396e-15, 1.9104136884609856e-15, 1.8594038275914604e-15, 1.806824469732757e-15, 1.7598345536654714e-15, 1.7111869057186746e-15, 1.6699715908879005e-15, 1.6233941476134983e-15, 1.5930112216218923e-15, 1.5662999368824204e-15, 1.4894076005862995e-15, 1.3840010288213201e-15, 1.3341775017445436e-15, 1.301816728698303e-15, 1.2727012108486492e-15, 1.2450267329412596e-15, 1.2184857674419958e-15, 1.1929228809315545e-15, 1.1682109297486229e-15, 1.144266261075052e-15, 1.1210276949710421e-15, 1.0984403185403167e-15, 1.0764548382299335e-15, 1.0550265557737052e-15, 1.03411398073964e-15, 1.013679223356684e-15, 9.936868552821639e-16, 9.74105607509137e-16, 9.549036134685393e-16, 9.36052422263822e-16, 9.17526019758457e-16, 8.992993179573792e-16, 8.813501289955468e-16, 8.636575071233003e-16, 8.462027111809957e-16, 8.289683401626569e-16, 8.119365819564001e-16, 7.950921445881629e-16, 7.784204169315239e-16, 7.619078868014965e-16, 7.455422451563377e-16, 7.2931181013360255e-16, 7.132063066308875e-16, 6.97215973013761e-16, 6.813321178897834e-16, 6.655476291862542e-16, 6.49855302680067e-16, 6.34248400783194e-16, 6.187213942454081e-16, 6.03269352675875e-16, 5.878879293211538e-16, 5.725738534582267e-16, 5.573245365217767e-16, 5.421373239701068e-16, 5.270114067509732e-16, 5.119464128870304e-16, 4.969421756270896e-16, 4.819991103661213e-16, 4.671186137983299e-16, 4.523033078628114e-16, 4.3755556051832413e-16, 4.2287835930904805e-16, 4.082755949654496e-16, 3.9375294808261796e-16, 3.7931701636102686e-16, 3.649737525841346e-16, 3.507293494536656e-16, 3.3659203419856035e-16, 3.2257012434341496e-16, 3.0867218027148446e-16, 2.949087797626608e-16, 2.812900537745497e-16, 2.678279400380918e-16, 2.545343659200951e-16, 2.4142244998306053e-16, 2.2850667406805856e-16, 2.1580125223803559e-16, 2.0331941507006968e-16, 1.9107746769139985e-16, 1.7909132854251696e-16, 1.6738031493317398e-16, 1.5596344104026648e-16, 1.4486138468680592e-16, 1.3409510988586057e-16, 1.2368474201468432e-16, 1.136535808484025e-16, 1.0402066406901132e-16, 9.480051650785084e-17, 8.600839988680856e-17, 7.76737697211545e-17, 6.983343750169644e-17, 6.252550855160264e-17, 5.578934378551636e-17, 4.966276256730797e-17, 4.4180799153105354e-17, 3.9383428333889435e-17]
    power_blk_expt  = 2.4432169118834085e9
    imag_k2_expt    = -2.4345416944075893e-6

    power_prf, power_blk, imag_k2 = Obliqua.calc_lovepy_tides_mush(omega, ecc, rho, radius, visc, shear, bulk, phi; ncalc=ncalc, material="maxwell")
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

if suite > 6
    # -------------
    # Test fluid interior data validity
    # -------------
    @info " "
    @info "Testing fluid interior data validity"
    ok = load.load_interior_liquid("$RES_DIR/interior_data/test_mantle_liquid.json", true)
    if ok
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Test Fluid.jl (minimal)
    # -------------
    @info " "
    @info "Testing Fluid.jl"
    omega, axial, ecc, sma, S_mass, rho, radius, visc = load.load_interior_liquid("$RES_DIR/interior_data/test_mantle_liquid.json", false)

    power_prf_expt  = [0.0, 4.8070888926231234e-5]
    power_blk_expt  = 8.501321352156166e20
    imag_k2_expt    = -0.010784662530332856

    power_prf, power_blk, imag_k2 = Obliqua.calc_fluid_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc; N_sigma=301, visc_l=2e2, visc_s=5e21)
    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= isapprox(imag_k2, imag_k2_expt; rtol=rtol)

    @info "Expected profile elements: $(power_prf_expt)"
    @info "Modelled profile elements: $(power_prf)"
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

if suite > 8

    # -------------
    # Test liquid interior data validity
    # -------------
    @info " "
    @info "Testing liquid interior data validity"
    ok = load.load_interior_liquid("$RES_DIR/interior_data/test_mantle_liquid_full.json", true)
    if ok
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Test Fluid.jl (PROTEUS-like)
    # -------------
    @info " "
    @info "Testing Fluid.jl"
    omega, axial, ecc, sma, S_mass, rho, radius, visc = load.load_interior_liquid("$RES_DIR/interior_data/test_mantle_liquid_full.json", false)

    power_prf_expt  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2425507160886101e-11, 1.247057869117265e-11, 1.2516453810781429e-11, 1.2563135220101134e-11, 1.2610676478283761e-11, 1.2659083790612889e-11, 1.2708409541623968e-11, 1.2758673237938062e-11, 1.2809914760255097e-11, 1.2862209173471728e-11, 1.2915579079219025e-11, 1.2969931827930746e-11, 1.302529416461838e-11, 1.3081604494621154e-11, 1.3138973237950605e-11, 1.3197436524617666e-11, 1.3257102950688364e-11, 1.3318044915122327e-11, 1.3380159375468342e-11, 1.3443554451719328e-11, 1.3507697369285256e-11, 1.357106154494241e-11, 1.3631454664018605e-11, 1.3688086768040071e-11, 1.3741180885641205e-11, 1.379128953020106e-11, 1.3839469103782372e-11, 1.3886772628443755e-11, 1.3933736994598853e-11, 1.3982642765511493e-11]
    power_blk_expt  = 2.2257483041633004e13
    imag_k2_expt    = -0.0006584045901003507

    power_prf, power_blk, imag_k2 = Obliqua.calc_fluid_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc; N_sigma=301, visc_l=2e2, visc_s=5e21)
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

if suite > 10

    # -------------
    # Test complete interior data validity
    # -------------
    @info " "
    @info "Testing liquid interior data validity"
    ok = load.load_interior_full("$RES_DIR/interior_data/test_mantle_full.json", true)
    if ok
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

    # -------------
    # Test complete model
    # -------------
    @info " "
    @info "Testing all components simultaneously"

    cfg = Obliqua.open_config("$ROOT_DIR/res/config/runtests.toml")
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, ncalc = load.load_interior_full("$RES_DIR/interior_data/test_mantle_full_test.json", false)

    power_prf_expt  = [0.0, 0.0, 7.426074810661149e-7, 7.141339323280178e-7, 6.880511303104244e-7, 6.626091914965803e-7, 6.393084884211161e-7, 6.162131326367927e-7, 5.95598470435145e-7, 5.741712909931592e-7, 5.57578255630565e-7, 5.450000281121666e-7, 5.270321787210745e-7, 5.029920224873281e-7, 4.853926718049137e-7, 4.7015690367438e-7, 4.5570518109796236e-7, 4.41806776705985e-7, 4.2839891795532446e-7, 4.154425605547741e-7, 4.0290520898949683e-7, 3.9075828339034296e-7, 3.7897643905499967e-7, 3.6753718859421477e-7, 3.564206267092027e-7, 3.4560920184076336e-7, 3.3508751235461023e-7, 3.2484212269440694e-7, 3.1486139498479976e-7, 3.0513533652935257e-7, 2.9855942934464557e-7, 2.8654193184246404e-7, 2.7753387945434864e-7, 2.68754834528521e-7, 2.602011904140257e-7, 2.5187032913422737e-7, 2.437605443654469e-7, 2.358709726657163e-7, 2.282015303083454e-7, 2.20752856110591e-7, 2.1352626051460859e-7, 2.065236789808777e-7, 1.9974763083113782e-7, 1.9320118109163912e-7, 1.8688790604185119e-7, 1.8081186297756487e-7, 1.749775620227741e-7, 1.6938994097483858e-7, 1.6405434131337383e-7, 1.5897648746429468e-7, 1.5416246903772564e-7, 1.496187233903219e-7, 1.453520211130342e-7, 1.4136945156206726e-7, 1.3767840940501013e-7, 1.3428658480378992e-7, 1.3120195176773663e-7, 1.284327574865411e-7, 1.2598751363159958e-7, 1.238749887479976e-7, 1.2210420151490848e-7, 1.2068441098915266e-7, 1.1962511098279694e-7, 1.1893602600971416e-7, 1.1862710506742638e-7, 1.1870851337185436e-7, 1.1919062194703133e-7, 2.713081566103371e-6, 3.266420036202856e-10, 3.370535778986049e-10, 3.4773777112819166e-10, 3.5870528031458646e-10, 3.699685635983859e-10, 3.8154200350601684e-10, 3.9344208721572494e-10, 4.0568760336251063e-10, 4.182998602181457e-10, 4.3130293013677005e-10, 4.4472394417583124e-10, 4.585934965818968e-10, 4.729463732312791e-10, 4.878243703500552e-10, 5.032946848400699e-10, 5.196043062869222e-10, 5.38444101227443e-10, 5.641786650732795e-10, 5.834129313913307e-10, 5.968777024182884e-10, 6.14640254532467e-10, 6.375776727432495e-10, 6.596453055585456e-10, 6.843684771633824e-10, 7.093114694267475e-10, 7.365466168359673e-10, 7.644677970177269e-10, 7.949482294572432e-10, 0.0, 0.0]
    power_blk_expt  = 1.4555219857478472e14
    imag_k2_expt    = [9.338013011939597e-8, 6.869356835427922e-7, 5.052955815066234e-6, 3.714832089304862e-5, 0.0002733212410749532, 0.09171577226656115, 0.015200057976338515, 0.09305722494783335, 0.18000650834364215, 0.0009399615893044621]

    power_prf, power_blk, σ_range, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, cfg
    )

    test_pass = true

    test_pass &= all(isapprox.(power_prf, power_prf_expt; rtol=rtol, atol=atol))
    test_pass &= isapprox(power_blk, power_blk_expt; rtol=rtol)
    test_pass &= all(isapprox.(imag_k2, imag_k2_expt; rtol=rtol, atol=atol))

    @info "First 5 expected profile elements: $(power_prf_expt[1:5])"
    @info "First 5 modelled profile elements: $(power_prf[1:5])"
    @info "Expected total power = $(power_blk_expt) W, Modelled total power = $(power_blk) W"
    @info "First 5 expected imag(k2) elements: $(imag_k2_expt[1:5])"
    @info "First 5 modelled imag(k2) elements: $(imag_k2[1:5])"

    if test_pass
        @info "Pass"
    else
        @warn "Fail"
        failed += 1
    end
    total += 1
    @info "--------------------------"

end

