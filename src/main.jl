# Core file containing functions for running the model

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module MAIN

    # Include system libraries
    using LoggingExtras
    using Printf
    import TOML:parsefile

    # Include local jl files (order matters)
    include("Love.jl")
    include("Fluid.jl")
    include("load.jl")


    # Import submodules
    import .Love
    import .Fluid
    import .load

    # Export submodules (mostly for autodoc purposes)
    export Love
    export Fluid
    export load

    export run_tides

    prec = BigFloat
    precc = Complex{BigFloat}

    # omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk

    function run_tides(omega::prec,
                        axial::prec,
                        ecc::prec,
                        sma::prec,
                        S_mass::prec,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{prec,1},
                        bulk::Array{prec,1};
                        ncalc::Int64
                        )::Tuple{Array{Float64,1},Float64,Float64}
      
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)

        visc_l = 1e2
        visc_s = 5e21

        mask_l = η .< visc_l
        mask_s = η .> visc_s

        r_c = r[1]                           # core radius
        r_s = r[2:end][mask_s]               # solid region
        r_m = r[2:end][.!mask_l .& .!mask_s] # mush region
        r_l = r[2:end][mask_l]               # liquid region

        power_prf_s_p, power_blk_s, imag_k2_s = Love.calc_lovepy_tides(
            omega, ecc, rho[mask_s], 
            vcat(r_c, r_s), visc[mask_s], 
            shear[mask_s], bulk[mask_s]; 
            ncalc=ncalc, 
            material="andrade"
        )

        # Get power profile
        power_prf_s = zeros(prec, length(ρ))
        power_prf_s[mask_s] .= power_prf_s_p

        power_prf_l, power_blk_l, imag_k2_l = Fluid.calc_fluid_tides(
            omega, axial, ecc, 
            sma, S_mass, rho, 
            radius, visc; 
            N_sigma=301, 
            visc_l=2e2, 
            visc_s=5e21
        )

        power_prf = power_prf_s .+ power_prf_l      
        power_blk = power_blk_s + power_blk_l
        imag_k2   = imag_k2_s + imag_k2_l

        save_heat_profile(r[2:end], power_prf)

        return power_prf, power_blk, imag_k2

    end

    using Plots

    """
        save_heat_profile(radius, power_prf; filename="heat_profile.png")

    Create and save a plot of the heat dissipation profile as a function of radius.
    """
    function save_heat_profile(radius::AbstractVector,
                            power_prf::AbstractVector;
                            filename::String = "heat_profile.png")

        # Filter out zero values (and negative ones, if they ever occur)
        mask = power_prf .> 0
        r_nonzero = radius[mask]
        p_nonzero = power_prf[mask]

        plt = plot(
            r_nonzero,
            p_nonzero,
            xlabel = "Radius (m)",
            ylabel = "Power per shell (W/kg)",
            title  = "Tidal Heating Profile",
            lw     = 2,
            legend = false,
            xscale = :log10,
            yscale = :log10,
            minorgrid = true,
        )

        savefig(plt, filename)
        return filename
    end


end