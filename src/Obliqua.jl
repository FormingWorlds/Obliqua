# Core file containing functions for running the model

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module Obliqua

    # Include system libraries
    using LoggingExtras
    using Printf
    import TOML:parsefile
    using Interpolations
    using LinearAlgebra


    # Include local jl files (order matters)
    include("Love.jl")
    include("Fluid.jl")
    include("load.jl")
    include("plotting.jl")


    # Import submodules
    import .Love
    import .Fluid
    import .load
    import .plotting

    # Export submodules (mostly for autodoc purposes)
    export Love
    export Fluid
    export load
    export plotting

    export run_tides

    prec = BigFloat
    precc = Complex{BigFloat}


    # Constants
    AU = 1.495978707e11  # m
    G = prec(6.6743e-11)  # m^3 kg^-1 s^-2

    res = 20.0


    """
        run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk; kwargs...)

    Compute the tidal heating profile of a planetary interior considering both **solid** and **liquid/mushy layers**.

    # Arguments
    - `omega::prec` : Orbital frequency of the body.
    - `axial::prec` : Axial (spin) frequency of the body.
    - `ecc::prec` : Orbital eccentricity.
    - `sma::prec` : Semi-major axis of the orbit.
    - `S_mass::prec` : Mass of the central body (e.g., star) inducing tides.
    - `rho::Array{prec,1}` : Radial density profile of the planet, from core to surface.
    - `radius::Array{prec,1}` : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}` : Viscosity profile of the planet.
    - `shear::Array{prec,1}` : Shear modulus profile of the solid layers.
    - `bulk::Array{prec,1}` : Bulk modulus profile of the solid layers.

    # Keyword Arguments
    - `ncalc::Int64=1000` : Number of points for tidal calculation in solid layers.
    - `N_sigma::Int64=301` : Number of frequency points for fluid tides.
    - `material::String="andrade"` : Rheology model used for solid tides.
    - `visc_l::Float64=1e2` : Lower threshold viscosity for defining liquid regions.
    - `visc_l_tol::Float64=5e2` : Tolerance added to `visc_l`.
    - `visc_s::Float64=1e22` : Upper threshold viscosity for defining solid regions.
    - `visc_s_tol::Float64=5e21` : Tolerance subtracted from `visc_s`.
    - `sigma_R::Float64=1e-3` : Rayleigh drag coefficient for fluid layers.

    # Returns
    - `power_prf::Array{Float64,1}` : Radial profile of tidal heating (W/m³).
    - `power_blk::Float64` : Total tidal power integrated over the interior (W).
    - `imag_k2::Float64` : Imaginary part of the Love number `k2` for the planet.

    # Notes
    - The function splits the interior into **solid**, **mush**, and **liquid** regions based on viscosity thresholds.
    - Solid tides are computed using `calc_lovepy_tides`, while fluid/mushy tides are computed using `calc_fluid_tides`.
    - The heating profile is automatically saved as a PDF using `save_heat_profile`.
    """
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
        power_prf_s, power_blk_s, imag_k2_s = calc_lovepy_tides(
            omega, ecc, rho[mask_s], 
            vcat(r_c, r_s), visc[mask_s], # assuming that the bottom of the solid region is r_c
            shear[mask_s], bulk[mask_s]; 
            ncalc=ncalc, 
            material=material
        )
        
        # Calculate liquid tides in liquid region 
        power_prf_l, power_blk_l, imag_k2_l = calc_fluid_tides(
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
        plotting.save_heat_profile(r[2:end], power_prf)

        return power_prf, power_blk, imag_k2

    end

    
    # Calculate heating from interior properties
    """
        calc_lovepy_tides(omega, ecc, rho, radius, visc, shear, bulk; ncalc=2000, material="andrade")

    Calculate tidal heating in a solid planetary interior using the Love number approach (`Love.jl`).

    # Arguments
    - `omega::prec` : Forcing angular frequency (rad/s).  
    - `ecc::prec` : Orbital eccentricity.  
    - `rho::Vector{prec}` : Radial density profile (kg/m³), first element = core.  
    - `radius::Vector{prec}` : Radial grid points corresponding to `rho` (m).  
    - `visc::Vector{prec}` : Viscosity profile (Pa·s).  
    - `shear::Vector{prec}` : Shear modulus profile (Pa).  
    - `bulk::Vector{prec}` : Bulk modulus profile (Pa).  

    # Keyword Arguments
    - `ncalc::Int=2000` : Number of subdivisions for radial layers.  
    - `material::String="andrade"` : Rheology model for complex shear modulus. Options: `"maxwell"`, `"andrade"`.

    # Returns
    A tuple `(power_prf, power_blk, k2_im)`:
    - `power_prf::Vector{Float64}` : Tidal heating profile per unit mass (W/kg).  
    - `power_blk::Float64` : Total bulk tidal heating (W).  
    - `k2_im::Float64` : Imaginary part of the complex tidal Love number k2 at the forcing frequency.  

    # Notes
    - Implements Efroimsky (2012) equations for viscoelastic response.  
    - Shear and bulk heating contributions are computed separately and summed.
    """
    function calc_lovepy_tides( omega::prec,
                                    ecc::prec,
                                    rho::Array{prec,1},
                                    radius::Array{prec,1},
                                    visc::Array{prec,1},
                                    shear::Array{prec,1},
                                    bulk::Array{prec,1};
                                    ncalc::Int=2000,
                                    material::String="andrade"
                                    )::Tuple{Array{Float64,1},Float64,Float64}

        # Internal structure arrays.
        # First element is the innermost layer, last element is the outermost layer
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μ = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        # Complex shear modulus (=(modulus of) rigidity)
        if material == "maxwell"
            # Maxwell material. 
            μc = 1im*μ*omega./(1im*omega.+μ./η)
        elseif material == "andrade"
            # Andrade material.
            μc = Love.andrade_mu_complex(omega, μ, η)
        else
            throw("Material type for complex shear modulus not defined, options: 'maxwell', 'andrade'.")
        end
        
        # See Efroimsky, M. 2012, Eqs. 3, 64 
        # To find k2 corresponding to andrade solid rheology insert relevant eqs here
        # Identical implementation exists in Farhat 2025 (Eqs not listed there)

        # Outer radius
        R = r[end]

        # Subdivide input layers such that we have ~ncalc in total
        rr = Love.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # Get gravity at each layer
        g = Love.get_g(rr, ρ);

        # Create grid
        Love.define_spherical_grid(res; n=2)

        # Get y-functions
        tidal_solution = Love.compute_y(rr, ρ, g, μc, κ)

        # Get k2 tidal Love Number (complex-valued)
        k2 = tidal_solution[5, end, end] - 1

        # Get bulk power output in watts
        power_blk = Love.get_total_heating(tidal_solution, omega, R, ecc)

        # Get profile power output (W m-3), converted to W/kg
        (Eμ, Eκ) = Love.get_heating_profile(tidal_solution,
                               rr, ρ, g, μc, κ,
                               omega, ecc)

        Eμ_tot, _ = Eμ   # shear       (W), (W/m3)
        Eκ_tot, _ = Eκ   # compaction  (W), (W/m3)

        power_prf = Eμ_tot .+ Eκ_tot # Compute total volumetric heating (W/m3)

        power_prf = power_prf ./ ρ # Convert to mass heating rate (W/kg)

        # Call Fluid script here
        # ...

        # Sum k2 love numberes
        # ...

        # Sum heating: bulk and profile
        # ...

        #return power_prf[11,:], power_blk, imag(k2)
        return power_prf, power_blk, imag(k2)
    end

    # Calculate heating from interior properties with mush
    """
        calc_lovepy_tides_mush(omega, ecc, rho, radius, visc, shear, bulk, phi;
                                ncalc=2000, material="andrade", visc_l=1e2,
                                bulk_l=1e9, permea=1e-7)

    Calculate tidal heating in a partially molten planetary interior with a mush region.

    # Arguments
    - `omega::prec` : Forcing angular frequency (rad/s).  
    - `ecc::prec` : Orbital eccentricity.  
    - `rho::Vector{prec}` : Radial density profile (kg/m³), first element = core.  
    - `radius::Vector{prec}` : Radial grid points corresponding to `rho` (m).  
    - `visc::Vector{prec}` : Viscosity profile (Pa·s).  
    - `shear::Vector{prec}` : Shear modulus profile (Pa).  
    - `bulk::Vector{prec}` : Bulk modulus profile (Pa).  
    - `phi::Vector{prec}` : Melt fraction profile (0–1) for each layer.

    # Keyword Arguments
    - `ncalc::Int=2000` : Number of subdivisions for radial layers.  
    - `material::String="andrade"` : Rheology model for complex shear modulus. Options: `"maxwell"`, `"andrade"`.  
    - `visc_l::Float64=1e2` : Liquid viscosity in mush region (Pa·s).  
    - `bulk_l::Float64=1e9` : Bulk modulus of liquid in mush region (Pa).  
    - `permea::Float64=1e-7` : Permeability of mush region (m²).

    # Returns
    A tuple `(power_prf, power_blk, k2_im)`:
    - `power_prf::Vector{Float64}` : Tidal heating profile per unit mass (W/kg).  
    - `power_blk::Float64` : Total bulk tidal heating (W).  
    - `k2_im::Float64` : Imaginary part of the complex tidal Love number k2 at the forcing frequency.  

    # Notes
    - Accounts for porosity and partial melt (mush) effects on tidal dissipation.  
    - Uses drained and Biot moduli to handle fluid–solid coupling.  
    - Throws an error if no mush region is identified (`phi` insufficient).
    """
    function calc_lovepy_tides_mush(omega::prec,
                                    ecc::prec,
                                    rho::Array{prec,1},
                                    radius::Array{prec,1},
                                    visc::Array{prec,1},
                                    shear::Array{prec,1},
                                    bulk::Array{prec,1},
                                    phi::Array{prec,1};
                                    ncalc::Int=2000,
                                    material::String="andrade",
                                    visc_l::Float64=1e2,
                                    bulk_l::Float64=1e9,
                                    permea::Float64=1e-7
                                    )::Tuple{Array{Float64,1},Float64,Float64}

        # Internal structure arrays.
        # First element is the innermost layer, last element is the outermost layer
        ρ  = convert(Vector{prec}, rho)      # density --> solid + liquid density
        r  = convert(Vector{prec}, radius)   # radius
        η  = convert(Vector{prec}, visc)     # shear viscosity
        μ  = convert(Vector{precc},shear)    # shear modulus
        κs = convert(Vector{prec}, bulk)     # solid bulk modulus
        ϕ  = convert(Vector{prec}, phi)      # melt fraction
        κd = 0.01.*κs                        # drained bulk modulus

        α = 1.0.-(κd./κs)                    # Biot's modulus

        # allocate zero arrays with same length and precision as r
        κl = zeros(prec, length(r))
        ηl = zeros(prec, length(r))
        k  = zeros(prec, length(r))

        # Find mush index
        ii = Love.find_mush_index(ϕ)
        # If no matches, throw error (because the matrix cannot be resolved, instead use 1 phase model)
        if ii === nothing
            throw("No mush region identified in viscosity profile.")
        end

        # update only the largest index that matches
        κl[ii] = prec(bulk_l)   # liquid bulk modulus
        ηl[ii] = prec(visc_l)   # liquid viscosity
        k[ii]  = prec(permea)   # permeability

        ρs = ρ.*(1.0.-ϕ)        # solid density 
        ρl = ρ.*ϕ               # liquid density

        # set porosity to zero outside mush region (otherwise code cannot solve system)
        ϕ[1:ii-1]   .= 0.0      # zero below ii
        ϕ[ii+1:end] .= 0.0      # zero above ii

        porous = true

        # Complex shear modulus (=(modulus of) rigidity)
        if material == "maxwell"
            # Maxwell material. 
            μc = 1im*μ.*omega./(1im*omega.+μ./η)
        elseif material == "andrade"
            # Andrade material.
            μc = Love.andrade_mu_complex(omega, μ, η)
        else
            throw("Material type for complex shear modulus not defined, options: 'maxwell', 'andrade'.")
        end
        
        # See Efroimsky, M. 2012, Eqs. 3, 64 
        # To find k2 corresponding to andrade solid rheology insert relevant eqs here
        # Identical implementation exists in Farhat 2025 (Eqs not listed there)

        # Outer radius
        R = r[end]

        # Subdivide input layers such that we have ~ncalc in total
        rr = Love.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # Get gravity at each layer
        g = Love.get_g(rr, ρ);

        # Create grid
        Love.define_spherical_grid(res; n=2)

        # Get y-functions        
        tidal_solution = Love.compute_y(rr, ρs, g, μc, κs, omega, ρl, κl, κd, α, ηl, ϕ, k)

        # Get k2 tidal Love Number (complex-valued)
        k2 = tidal_solution[5, end, end] - 1

        # Get bulk power output in watts
        power_blk = Love.get_total_heating(tidal_solution, omega, R, ecc)

        # Get profile power output (W m-3), converted to W/kg
        Eμ, Eκ, El = Love.get_heating_profile(tidal_solution,
                               rr, ρs, g, μc, κs,
                               omega, ρl, κl, κd, 
                               α, ηl, ϕ, k, ecc)

        Eμ_tot, _ = Eμ   # shear       (W), (W/m3)
        Eκ_tot, _ = Eκ   # compaction  (W), (W/m3)
        El_tot, _ = El   # fluid       (W), (W/m3)

        power_prf = Eμ_tot .+ Eκ_tot .+ El_tot# Compute total volumetric heating (W/m3)

        power_prf = power_prf ./ ρ # Convert to mass heating rate (W/kg)

        # Call Fluid script here
        # ...

        # Sum k2 love numberes
        # ...

        # Sum heating: bulk and profile
        # ...

        #return power_prf[11,:], power_blk, imag(k2)
        return power_prf, power_blk, imag(k2)
    end


    # Get fluid tidal heating
    """
        calc_fluid_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc; 
                        N_sigma=301, visc_l=5e2, visc_s=5e21, sigma_R=1e-3)

    Calculate the tidal heating in the fluid layers of a planetary interior.

    This function computes the tidal dissipation power profile and the total tidal heating 
    for liquid layers using Love numbers, Hansen coefficients, and a viscoelastic 
    approximation for the solid interior. It also returns the k2 Love number at the 
    forcing frequency used by the calculation.

    # Arguments
    - `omega::prec` : Orbital frequency of the planet [rad/s].
    - `axial::prec` : Axial rotation frequency of the planet [rad/s].
    - `ecc::prec`   : Orbital eccentricity.
    - `sma::prec`   : Semi-major axis of the orbit [m].
    - `S_mass::prec` : Stellar mass [kg].
    - `rho::Array{prec,1}` : Layered density profile [kg/m³].
    - `radius::Array{prec,1}` : Layer boundaries [m], length = number of layers + 1.
    - `visc::Array{prec,1}` : Layer viscosities [Pa·s].

    # Keyword Arguments
    - `N_sigma::Int=301` : Number of frequency bins for spectral calculations.
    - `visc_l::Float64=5e2` : Threshold viscosity below which a layer is considered liquid [Pa·s].
    - `visc_s::Float64=5e21` : Threshold viscosity above which a layer is considered solid [Pa·s].
    - `sigma_R::Float64=1e-3` : Rayleigh drag coefficient for the fluid interface [1/s].

    # Returns
    A tuple with three elements:

    1. `power_prf::Vector{Float64}` : Tidal heating power profile for each layer [W].
    2. `P_tidal_total::Float64` : Total tidal dissipation power in the fluid layers [W].
    3. `k2_total::Float64` : Complex k2 Love number at the forcing frequency.

    # Errors
    Throws an error if:
    - No liquid layers are found in the interior (`ρ_l` empty), suggesting that `visc_l` needs adjustment.
    - The interior structure is incompatible with computing mean densities.

    # Notes
    - Layers are categorized into liquid, mush, and solid based on `visc_l` and `visc_s`.
    - `r_b` denotes the bottom of the liquid region.
    - Hansen coefficients and Love numbers are computed over a frequency spectrum and interpolated.
    - The total tidal dissipation is computed using both the fluid and solid contributions.

    """
    function calc_fluid_tides( omega::prec,
                        axial::prec,
                        ecc::prec,
                        sma::prec,
                        S_mass::prec,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1};
                        N_sigma::Int=301,
                        visc_l::Float64=5e2,
                        visc_s::Float64=5e21,
                        sigma_R::Float64=1e-3
                        )::Tuple{Array{Float64,1},Float64,Float64}

        # Degree Love number
        n = 2
        m = 2
        k_min = -30
        k_max = 40
        k_range = collect(k_min:k_max)

        # Internal structure arrays.
        # First element is the innermost layer, last element is the outermost layer
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)

        # Convert profiles to scalar quantities (1 dimensional)
        
        R = maximum(r)  # Planet radius (m)
        g = G * sum(4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) .* ρ[1:end]) / R^2  # Surface gravity (m/s^2)

        # g = Love.get_g(r, ρ)
        # g_s = g[end] # Surface gravity (m/s^2)


        # The follwoing code should be 
        # moved to a separate function

        mask_l = η .< visc_l
        mask_s = η .> visc_s

        r_c = r[1]                           # core radius
        r_s = r[2:end][mask_s]               # solid region
        r_m = r[2:end][.!mask_l .& .!mask_s] # mush region
        r_l = r[2:end][mask_l]               # liquid region

        # bottom of liquid region
        r_b = maximum([
            r_c,
            isempty(r_s) ? -Inf : r_s[end],
            isempty(r_m) ? -Inf : r_m[end]
        ])

        ρ_l = ρ[mask_l]
        ρ_s = ρ[mask_s]
        ρ_m = ρ[.!mask_l .& .!mask_s] # mush region

        #######
        #######

        # mean densities
        if length(ρ_l) == 0
            error("No liquid layers found in the interior structure. Adjust visc_l parameter.")
        elseif length(ρ_l) == 1
            ρ_l_mean = ρ_l[1]
        else
            ρ_l_mean = Fluid.mean_rho(ρ_l, r_l, r_b)
        end
        
        if length(ρ_s) == 0
            ρ_s_mean = 5000.0  # arbitrary high density
        elseif length(ρ_s) == 1
            ρ_s_mean = ρ_s[1]
        else
            ρ_s_mean = Fluid.mean_rho(ρ_s, r_s, r_c)
        end
        
        # magma ocean height
        if length(r_l) == 0
            error("No liquid layers found in the interior structure.")
        else
            H_magma = r_l[end] - r_b
        end

        # density ratio
        ρ_ratio = ρ_l_mean / ρ_s_mean

        
        # get hansen coefficients
        k_range2, X_hansen = Fluid.get_hansen(ecc, n, m, k_min, k_max)

        # orbital and axial frequencies
        t_range = 10 .^ range(-15, stop=6, length=N_sigma)   # periods        
        σ_range = 2π ./ (t_range .* 1e3 .* 365.25 .* 24 .* 3600)
        σ_range = reshape(σ_range, :)

        # preallocate (complex for viscoelastic)
        k_T_homo = zeros(precc, n, N_sigma)
        k_L_homo = zeros(precc, n, N_sigma)

        # could include Andrade solid tides here, instead of propagator method
        # ...

        # get 2,2 harmonic
        k_T_22_homo = vec(k_T_homo[2, :])
        k_L_22_homo = vec(k_L_homo[2, :])

        # get Rayleigh drag at interface, needs to be changed
        # employ correlation with mixing length
        σ_R = sigma_R 

        # fluid Love Numbers
        k22_fluid_high_friction, k22_total = Fluid.compute_fluid_lovenumbers(
            n,
            σ_range,
            k_T_22_homo,
            k_L_22_homo,
            ρ_ratio,
            g,
            H_magma,
            σ_R,
            R
        )

        # interpolate k2 love number arrays
        μ_n  = n * (n + 1)
        ξ_n = 3.0 / (2.0*n + 1.0) * ρ_ratio
        σ_n = sqrt(μ_n * g * H_magma / R^2) # characterestic frequency

        for kk in 1:N_sigma
            σ = σ_range[kk]
            σ_T = σ - im*σ_R

            k22_fluid_high_friction[kk] =
                -ξ_n * σ_n^2 / (σ*σ_T - σ_n^2)

            k22_total[kk] =
                k_T_22_homo[kk] + (1 + k_L_22_homo[kk])*k22_fluid_high_friction[kk]
        end

        k22_total = vec(k22_total)    # reshape(-1)

        # Build symmetric full spectrum for interpolation
        full_σ_range   = vcat(-σ_range,     reverse(σ_range))
        full_k22_total = vcat(-k22_total,   reverse(k22_total))
        full_k22_homo  = vcat(-k_T_22_homo, reverse(k_T_22_homo))

        imag_full_k22  = imag.(full_k22_total)
        imag_solid_k22 = imag.(full_k22_homo)

        # interpolation functions for imaginary parts (extrapolate outside)
        interp_full  = extrapolate(interpolate((full_σ_range,), imag_full_k22,
                                            Gridded(Linear())), Flat())
        interp_solid = extrapolate(interpolate((full_σ_range,), imag_solid_k22,
                                            Gridded(Linear())), Flat())

        # calculate tidal heating
        A_22k_e      = zeros(prec,  length(k_range))
        U_22k_e      = zeros(precc, length(k_range))
        P_T_k_total  = zeros(prec,  length(k_range))
        P_T_k_solid  = zeros(prec,  length(k_range))

        for (ikk, kk) in pairs(k_range)
            σ = 2*axial - kk*omega

            # Eq. 33
            A_22k_e[ikk] = sqrt(6π/5) * X_hansen[ikk]

            # Eq. 32
            U_22k_e[ikk] = (G*S_mass/sma) * (R/sma)^2 * A_22k_e[ikk]

            img_full_k22  = interp_full(σ)
            img_solid_k22 = interp_solid(σ)

            prefactor = 5 * R * σ / (8π*G)
            U2 = abs2(U_22k_e[ikk])

            P_T_k_total[ikk] = prefactor * img_full_k22  * U2
            P_T_k_solid[ikk] = prefactor * img_solid_k22 * U2
        end

        # Total tidal heating (negative sum)
        P_tidal_total = -sum(P_T_k_total) # W
        P_tidal_solid = -sum(P_T_k_solid) # W

        # Get power profile
        power_prf = zeros(prec, length(ρ))
        power_prf[mask_l] .= Fluid.heat_profile(P_tidal_total, ρ_l, r_l, r_b)
        
        # Get k2 lovenumber at forcing frequency used by Love.jl
        k2_total = interp_full(2*axial - omega)
        
        return convert(Vector{Float64}, power_prf), convert(Float64, P_tidal_total), convert(Float64, k2_total)

    end

end