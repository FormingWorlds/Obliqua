# Core file containing functions for running the model

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module main

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
                        ncalc::Int64=1000,
                        N_sigma::Int64=301,
                        material::String="andrade",
                        visc_l::Float64=1e2,
                        visc_l_tol::Float64=5e2,
                        visc_s::Float64=1e22,
                        visc_s_tol::Float64=5e21,
                        sigma_R::Float64=1e-3
                        )::Tuple{Array{Float64,1},Float64,Float64}
      
        # Interior profiles                        
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)

        # Liquidus and Solidus viscosity with tolerance
        η_l = visc_l + visc_l_tol
        η_s = visc_s - visc_s_tol

        # Masks for liquid and solid regions
        mask_l = η .< η_l
        mask_s = η .> η_s

        # Core radius and Arrays with radii that belong to solid/mush/liquid 
        r_c = r[1]                           # core radius
        r_s = r[2:end][mask_s]               # solid region
        r_m = r[2:end][.!mask_l .& .!mask_s] # mush region
        r_l = r[2:end][mask_l]               # liquid region

        # Check/Validate regions
        # ...

        # Calculate solid tides in solid region 
        power_prf_s, power_blk_s, imag_k2_s = Love.calc_lovepy_tides(
            omega, ecc, rho[mask_s], 
            vcat(r_c, r_s), visc[mask_s], # assuming that the bottom of the solid region is r_c
            shear[mask_s], bulk[mask_s]; 
            ncalc=ncalc, 
            material=material
        )
        
        # Calculate liquid tides in liquid region 
        power_prf_l, power_blk_l, imag_k2_l = Fluid.calc_fluid_tides(
            omega, axial, ecc, 
            sma, S_mass, rho, 
            radius, visc; # pass all radii, since density ratio is derived (tbd)
            N_sigma=N_sigma, 
            visc_l=η_l, 
            visc_s=η_s,
            sigma_R=sigma_R # currently = 1e-3, in the future this may turn into an array, as a function of radii (tbd)
        )

        # Get power profile
        power_prf = zeros(prec, length(ρ))

        power_prf[mask_s] .= power_prf_s        # specify solid tidal heating in solid region
        #power_prf[mask_l] .= power_prf_l        # specify liquid tidal heating in liquid region
        
        # Until Liquid tides can be called with only liquid radii array we need to sum as:
        power_prf .+= power_prf_l  # sum with liquid (tbd)

        # Get total (bulk) heating
        power_blk = power_blk_s + power_blk_l

        # Get total k2 Love Number
        imag_k2   = imag_k2_s + imag_k2_l

        # Save heating profile as pdf (image)
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
            yscale = :log10,
            minorgrid = true,
        )

        savefig(plt, filename)
        return filename
    end


end