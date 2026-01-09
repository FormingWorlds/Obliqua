
#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get Love root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
# import Pkg
# Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras
using Obliqua
using Plots
using Printf

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

if suite > 10

    @info " "
    @info "Testing Obliqua.jl"

    # Load configuration
    cfg = Obliqua.open_config("/home/marijn/LovePy/fwlLove.jl/res/config/all_options.toml")

    # Load interior model
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, ncalc =
        load.load_interior_full("$RES_DIR/interior_data/test_mantle_full_test.json", false)

    # Run tidal calculation
    power_prf, power_blk, σ_range, imag_k2 = Obliqua.run_tides(
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk; cfg=cfg
    )

    println(power_prf)
    println(power_blk)
    println(imag_k2)

    # Save k2 spectrum
    plotting.plot_imagk2_spectrum(σ_range, imag_k2; outpath="$OUT_DIR/k2_spectrum.png")

    # Save heating profile
    plotting.save_heat_profile(radius[2:end], power_prf; filename="$OUT_DIR/heat_profile.png")

    @info "k2 spectrum saved to $OUT_DIR/k2_spectrum.png"
    @info "Heating profile saved to $OUT_DIR/heat_profile.png"
    @info "--------------------------"

end

