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
    using DoubleFloats

    # Include local jl files (order matters)
    include("Love.jl")
    include("Fluid.jl")
    include("Hansen.jl")
    include("load.jl")
    include("plotting.jl")

    # Import submodules
    import .Love
    import .Fluid
    import .Hansen
    import .load
    import .plotting

    # Export submodules (mostly for autodoc purposes)
    export Love
    export Fluid
    export Hansen
    export load
    export plotting

    export run_tides

    const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../")

    prec  = BigFloat
    precc = Complex{BigFloat}

    const AU::prec = prec(1.495978707e11)   # m
    const G::prec  = prec(6.6743e-11)       # m^3 kg^-1 s^-2

    const res::Float64 = 20.0               # angular resolution in degrees


    """
        Open and validate config file.

    Arguments:
    - `cfg_path::String`                : Path to configuration file

    Returns:
    - `cfg_dict::Dict`                  : Dictionary containing the configuration
    """
    function open_config(cfg_path::String)::Dict

        # open file
        cfg_dict = parsefile(cfg_path)

        # check headers
        headers = ["params", "star", "orbit", "struct", "title", "version"]
        for h in headers
            if !haskey(cfg_dict, h)
                error("Key $h is missing from configuration file at '$cfg_path'")
            end
        end

        # check that output dir is named
        if !haskey(cfg_dict["params"]["out"],"path") || (cfg_dict["params"]["out"]=="")
            error("Output directory is missing from configuration file at '$cfg_path'")
        end
        out_path = abspath(cfg_dict["params"]["out"]["path"])

        # check if this is a dangerous path
        if ispath(joinpath(out_path, ".git")) || (joinpath(out_path) == pwd()) || samefile(out_path, ROOT_DIR)
            error("Output directory is unsafe")
        end

        # looks good
        return cfg_dict
    end


    """
        run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, cfg)

    Compute the tidal heating profile of a planetary interior considering solid and fluid layers.

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
    - `cfg::Dict` : Configuration parameters from dictionary.

    # Returns
    - `power_prf::Array{Float64,1}` : Radial profile of tidal heating (W/m³).
    - `power_blk::Float64` : Total tidal power integrated over the interior (W).
    - `σ_range::Array{Float64,1}` : Frequencies at which the Love number `k2` was evaluated.
    - `imag_k2::Array{Float64,1}` : Imaginary part of the Love number `k2` for the planet.
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
                        bulk::Array{prec,1},
                        cfg::Dict
                        )::Tuple{Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}}
      
        # Read configuration options from dict
        @info "Using configuration '$(cfg["title"])'"
        
        # Check that config has these always-required keys
        req_keys = Dict(
            "params.out" => ["path"],
            "star" => ["mass"],
            "orbit" => ["semimajoraxis", "eccentricity"],
            "orbit.obliqua" => [
                "min_frac","visc_l","visc_l_tol","visc_s","visc_s_tol",
                "sigma_R","n","m","N_sigma","ncalc","k_min","k_max",
                "p_min","p_max","material","strain"
            ],
            "struct" => ["mass_tot","core_density"]
        )
        for (section, keys) in req_keys
            path = split(section, ".")
            node = cfg

            # walk down the nested Dict
            for p in path
                if !haskey(node, p)
                    @error "Config: missing required section `$(join(path, "."))`"
                    return false
                end
                node = node[p]
            end

            # check required keys at this level
            for k in keys
                if !haskey(node, k)
                    @error "Config: missing required key `$(join(path, "."))::$k`"
                    return false
                end
            end
        end

        # collection of config params 
        min_frac = cfg["orbit"]["obliqua"]["min_frac"]
        sigma_R  = cfg["orbit"]["obliqua"]["sigma_R"]
        n        = cfg["orbit"]["obliqua"]["n"]
        m        = cfg["orbit"]["obliqua"]["m"]
        N_σ      = cfg["orbit"]["obliqua"]["N_sigma"]
        ncalc    = cfg["orbit"]["obliqua"]["ncalc"]
        k_min    = cfg["orbit"]["obliqua"]["k_min"]
        k_max    = cfg["orbit"]["obliqua"]["k_max"]
        p_min    = cfg["orbit"]["obliqua"]["p_min"]
        p_max    = cfg["orbit"]["obliqua"]["p_max"]
        material = cfg["orbit"]["obliqua"]["material"]
        alpha    = cfg["orbit"]["obliqua"]["alpha"]
        strain   = cfg["orbit"]["obliqua"]["strain"]

        # convert interior profiles to BigFloat                 
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μ = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        # number of layers
        N_layers = length(r)-1

        # liquidus and solidus viscosity with tolerance
        η_l = cfg["orbit"]["obliqua"]["visc_l"] + cfg["orbit"]["obliqua"]["visc_l_tol"]
        η_s = cfg["orbit"]["obliqua"]["visc_s"] - cfg["orbit"]["obliqua"]["visc_s_tol"]

        # get smoothed out region masks and segments
        mask_l, mask_s, mask_c, is_seg, segments = get_layers(r, η, η_l, η_s; min_frac)

        # check if CMB is at the bottom of the mantle
        if any(r[1] .>= r[2:end])
            throw("CMB radius not at bottom of mantle, did you properly order the interior arrays?")
        end

        # find planet radius (m)
        R = maximum(r)

        # tidal mode range (k is the Fourier index in mean anomaly)
        k_range = collect(k_min:k_max)

        # get hansen coefficients
        k_range2, X_hansen = Hansen.get_hansen(ecc, n, m, k_min, k_max)

        # orbital and axial frequencies
        t_range = 10 .^ range(p_min, stop=p_max, length=N_σ)        # periods [1e3 yr]       
        σ_range = 2π ./ (t_range .* 1e3 .* 365.25 .* 24 .* 3600)    # freq    [s-1]
        σ_range = reshape(σ_range, :)

        # get forcing frequency dependent complex shear modulus
        μc = complex_mu(σ_range, μ, η; material=material, α=alpha)

        # initiate forcing frequency dependent k2 love and load number arrays (one spectrum for each segment)
        k22_T = zeros(precc, N_σ, length(segments))
        k22_L = zeros(precc, N_σ, length(segments))
        
        # initiate forcing frequency dependent heating profile 
        prf_total = zeros(prec, N_σ, N_layers)

        # arbitrary high density for bottom boundary
        ρ_mean_lower = Float64(cfg["struct"]["core_density"])

        # loop over segments, starting at CMB
        for (iseg, seg) in pairs(segments)

            # preallocate (complex for viscoelastic)
            #  (T)idal love number
            k22_T_seg = zeros(precc, N_σ)

            #  (L)oad  love number
            k22_L_seg = zeros(precc, N_σ)

            # get start and stop index for segment
            i_start, i_end = is_seg[iseg]

            # perform slices
            r_seg  = r[i_start-1:i_end]
            ρ_seg  = ρ[i_start:i_end]
            η_seg  = η[i_start:i_end]                              
            μc_seg = μc[i_start:i_end, :] 
            κ_seg  = κ[i_start:i_end]

            # preallocate heating profile for segment
            prf_seg = zeros(prec, N_σ, length(r_seg)-1)

            # mean density in current segment
            if length(ρ_seg) == 1
                ρ_mean = ρ_seg[1]
            else
                ρ_mean = Fluid.mean_rho(ρ_seg, r_seg)
            end

            # density ratio
            ρ_ratio = ρ_mean / ρ_mean_lower

            # get k2 spectrum for segment
            for i in 1:N_σ
                # specify forcing frequency
                σ = σ_range[i]
                
                # preallocate k2 for segment
                kT = zero(precc)
                kL = zero(precc)

                # if segment is solid
                if seg == "solid"
                    # if heating profile from strain tensor
                    if strain==true
                        # calculate tides in solid region 
                        prf_seg[i,:], kT, kL = run_solid_strain( 
                            σ, ecc, ρ_seg,  
                            r_seg, η_seg,                               
                            μc_seg[:, i], κ_seg; 
                            ncalc=ncalc
                        )
                    # else no heating profile in segment 
                    #   --> global heating profile from complex shear modulus
                    else
                        # calculate tides in solid region 
                        kT, kL = run_solid( 
                            ρ_seg, r_seg, η_seg,                               
                            μc_seg[:, i], κ_seg; 
                            ncalc=ncalc
                        )

                    end

                # if segment is fluid
                elseif seg == "fluid"
                    # calculate fluid tides in fluid region 
                    kT, kL = run_fluid(
                        σ, ρ_seg, 
                        r_seg, ρ_ratio;
                        n=n, 
                        sigma_R=sigma_R
                    ) 
                
                # if segment is mush
                elseif seg == "mush"
                    # calculate mush tides in mush region 
                    kT, kL = 0., 0. # no expression for this yet

                # if segment is ice
                elseif seg == "ice"
                    # calculate ice tides in ice region 
                    kT, kL = 0., 0. # no expression for this yet
                           
                # if segment is water
                elseif seg == "water"
                    # calculate water tides in water region 
                    kT, kL = 0., 0. # no expression for this yet
                          
                end

                # update k2 spectrum for segment
                k22_T_seg[i] = kT
                k22_L_seg[i] = kL

            # repeat for all probe forcing frequencies
            end
        
            # update previous segment mean density before moving to next segment
            ρ_mean_lower = ρ_mean

            # store k2 spectra (max 1 per segment)
            k22_T[:, iseg] .= k22_T_seg
            k22_L[:, iseg] .= k22_L_seg

            # append segment heating profile to global heating profile
            prf_total[:, i_start:i_end] .= prf_seg[:, :]

            # step to next segment
        end

        # initialize total k2 with the contribution from the top layer
        k22_total = copy(k22_T[:, end])  

        # loop from top (surface) to just above CMB
        for iseg in reverse(1:length(segments)-1)
            for i in 1:N_σ
                k22_total[i] = k22_T[i, iseg] + (1.0 + k22_L[i, iseg]) * k22_total[i]
            end
        end
        
        # extract imaginary part of complex global k2 spectrum
        imag_k2 = .-imag.(k22_total)

        # build symmetric full spectrum for interpolation
        full_σ_range  = vcat(-σ_range, reverse(σ_range))
        imag_full_k22 = vcat(-imag_k2, reverse(imag_k2))

        # interpolation function for imaginary part (extrapolate outside)
        interp_full  = extrapolate(interpolate((full_σ_range,), imag_full_k22,
                                            Gridded(Linear())), Flat())

        # if segment wise heating profiles are calculated, interpolate heating in each layer across forcing frequency domain
        if strain==true
            # build symmetric full spectrum for interpolation
            full_prf_total = vcat(-prf_total,   reverse(prf_total))

            # create an interpolator per radial shell
            prf_itp_shells = Vector{Any}(undef, N_layers)

            # interpolation functions for heating profile (extrapolate outside)
            for j in 1:N_layers
                prf_layer = Float64.(full_prf_total[:, j])
                itp = extrapolate(interpolate((Float64.(full_σ_range),), prf_layer, Gridded(Linear())), Flat())
                prf_itp_shells[j] = itp
            end

            # define a function to get the radial profile at a given σ
            function radial_profile_at_sigma(σ::Float64, prf_itp_shells::Vector)
                return [itp(σ) for itp in prf_itp_shells]
            end
        end

        # calculate tidal heating
        # initialize frequency dependent quentities
        A_22k_e     = zeros(prec,  length(k_range))
        U_22k_e     = zeros(precc, length(k_range))

        # initialize frequency dependent total heating
        P_T_k_total = zeros(prec,  length(k_range))

        # initialize frequency dependent heating profile
        P_T_k_prf = zeros(prec,  length(k_range), length(shear))

        # loop over tidal modes 
        for (ikk, kk) in pairs(k_range)
            # calculate physical forcing frequency
            σ = m*axial - kk*omega

            # calculate coefficients
            A_22k_e[ikk] = sqrt(6π/5) * X_hansen[ikk]
            U_22k_e[ikk] = (G*S_mass/sma) * (R/sma)^2 * A_22k_e[ikk]

            # get imaginary part of complex k2 love number from global spectrum at forcing frequency
            img_full_k22 = interp_full(σ)

            # calculate prefactor and total availible heat
            prefactor = 5 * R * σ / (8π*G)
            U2 = abs2(U_22k_e[ikk])

            # calculate total heat input at forcing frequency
            P_T_k_total[ikk] = prefactor * img_full_k22  * U2

            # obtain heating profile
            # if the forcing frequenccy is zero, then the total heat input is zero.
            if σ == 0.
                continue
            else
                # if segment wise heating profile exists, normalize it to the obtained total heating
                if strain==true
                    # get global heating profile from at forcing frequency
                    unorm_prf = radial_profile_at_sigma(Float64.(σ), prf_itp_shells)
            
                    # determine the unnormalized total heat input
                    shell_volumes = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) 
                    unorm_tot = sum(shell_volumes .* unorm_prf)

                    # determine the normalization
                    norm = P_T_k_total[ikk] / unorm_tot

                    # normalize the heating profile 
                    P_T_k_prf[ikk, :] = unorm_prf .* norm

                # if no segment wise heating profiles exist, assume 
                # disipation propto imaginary part of the complex shear modulus.
                else
                    # get complex shear modulus profile at forcing frequency
                    μc = complex_mu([σ], μ, η; material=material, α=alpha)[:,1]

                    # get heating profile
                    P_T_k_prf[ikk, :] = radial_heating_profile(r, μc, P_T_k_total[ikk])
                end
            end
        end

        # total tidal heating
        power_blk = sum(P_T_k_total) # W

        # get radial heating profile W/m^3
        power_prf = [sum(P_T_k_prf[:,j]) for j in 1:size(P_T_k_prf,2)]
        power_prf ./ ρ # convert to mass heating rate (W/kg)

        # convert everything to Float64
        return Float64.(power_prf), power_blk, Float64.(σ_range), Float64.(imag_k2)

    end


    """
        get_layers(r, η, η_l, η_s; min_frac=0.02)

    Determine the phase profile of a planetary interior considering solid, fluid, and mush layers.

    # Arguments
    - `r::Array{prec,1}`                : Radial positions of layers, from core to surface.
    - `η::Array{prec,1}`                : Viscosity profile of the planet.
    - `η_l::Float64`                    : Liquidus viscosity.
    - `η_s::Float64`                    : Solidus viscosity.
    
    # Keyword Arguments
    - `min_frac::Float64=0.02`          : Minimal segment radius fraction before smoothing.

    # Returns
    - `mask_s::Vector{Bool}`            : Solid region mask.
    - `mask_l::Vector{Bool}`            : Fluid region mask.
    - `mask_c::Vector{Bool}`            : Mush region mask.
    - `is_seg::Vector{Tuple{Int,Int}}`  : Segment [start, stop] index array.
    - `segments::Vector{String}`        : Segment phase array.
    """
    function get_layers(r::Array{prec,1},
                        η::Array{prec,1},
                        η_l::Float64,
                        η_s::Float64;
                        min_frac::Float64=0.02
                        )::Tuple{Vector{Bool},Vector{Bool},Vector{Bool},Vector{Tuple{Int,Int}},Vector{String}}

        # masks for liquid and solid regions
        mask_l = η .< η_l
        mask_s = η .> η_s

        # total mantle thickness
        H = r[end] - r[1]
        N = length(r) - 1 # subtract CMB layer

        # phase encoding:
        #   0 = mush / other
        #   1 = solid
        #   2 = liquid
        phase_prf = zeros(Int, N)

        # build phase profile
        for i in 1:N
            if mask_s[i]
                phase_prf[i] = 1
            elseif mask_l[i]
                phase_prf[i] = 2
            end
        end

        # smooth phase profile
        #   moving up from the CMB
        i = 2
        while i < N
            p = phase_prf[i]
            i_start = i

            # extend segment while in same phase
            while i < N && phase_prf[i] == p
                i += 1
            end

            i_end = i - 1

            # check island layer 
            if i_start > 1 && i_end < N
                # get neighbouring segment phases
                p_below = phase_prf[i_start - 1]
                p_above = phase_prf[i_end + 1]

                # if sandwiched, check island layer thickness
                if p_below == p_above && p_below != p
                    dr = r[i_end] - r[i_start]
                    # if less then threshold fraction smooth phase profile
                    if dr / H < min_frac
                        phase_prf[i_start:i_end] .= p_below
                    end
                end
            end

            i += 1
        end

        # update masks
        mask_s .= phase_prf .== 1
        mask_l .= phase_prf .== 2

        # build segments + bottom mask
        segments = String[]
        mask_c = falses(length(r))
              
        # build boundary indices array
        is_seg = Vector{Tuple{Int,Int}}()

        i = 2
        while i <= N
            p = phase_prf[i]

            # mark bottom of this segment
            mask_c[i-1] = true # "bottom" radius is r[i-1]

            i_start = i

            # step forward through this segment
            while i <= N && phase_prf[i] == p
                i += 1
            end

            i_end = i - 1

            # store boundary indices
            push!(is_seg, (i_start, i_end))

            # store segments
            push!(segments,
                p == 1 ? "solid" :
                p == 2 ? "fluid" :
                        "mush"
            )
        end

        return mask_s, mask_l, mask_c, is_seg, segments
    end


    """
        complex_mu(σ_range, μ_profile, η_profile; material="andrade", α=0.3)

    Return the complex shear modulus μ̃(σ) for Maxwell or Andrade rheology.

    # Arguments
    - `σ_range::AbstractVector`         : Forcing frequency range.
    - `μ_profile::Array{precc,1}`       : Shear profile of the planet (aka unrelaxed rigidity).
    - `η_profile::Array{prec,1}`        : Viscosity profile of the planet.
    
    # Keyword Arguments
    - `material::String="andrade"`      : Material for which to find complex shear modulus.
    - `α::Float64=0.3"`                 : Power-law exponent (free parameter).

    # Returns
    - `μc::Matrix{precc}`               : Complex shear modulus profile at all forcing frequencies.
    """
    function complex_mu(σ_range::AbstractVector,
                            μ_profile::Array{precc,1},
                            η_profile::Array{prec,1};
                            material::String="andrade", 
                            α::Float64=0.3
                            )::Matrix{precc}

        nlayer = length(μ_profile)
        nfreq  = length(σ_range)

        # initialize output array
        μc = zeros(precc, nlayer, nfreq)

        if material == "maxwell"
            @inbounds for i in 1:nlayer
                μ = μ_profile[i]
                η = η_profile[i]
                for j in 1:nfreq
                    σ = σ_range[j]
                    μc[i,j] = 1im*μ*σ / (1im*σ + μ/η)
                end
            end
        elseif material == "andrade"
            @inbounds for i in 1:nlayer
                μ = μ_profile[i]
                η = η_profile[i]
                for j in 1:nfreq
                    σ = σ_range[j]
                    τM = η ./ μ # Maxwell time
                    τA = τM     # Andrade time 
                    term_andrade = gamma(1 + α) .* (1im .* σ .* τA).^(-α)
                    term_maxwell = (1im .* σ .* τM).^(-1)

                    μc[i,j] = μ ./ (1 .+ term_andrade .+ term_maxwell)
                end
            end
        else
            throw(ArgumentError("Material type for complex shear modulus not defined; options: 'maxwell' or 'andrade'"))
        end

        return μc
    end


    """
        run_solid_strain(omega, ecc, rho, radius, visc, shear, bulk; ncalc=2000)

    Calculate k2 Lovenumbers in the solid, and compute 1D heating profile from strain tensor.

    # Arguments
    - `omega::Float64`                  : Forcing frequency range.
    - `ecc::prec`                       : Eccentricity of the orbit.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `μ_profile::Array{precc,1}`       : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    
    # Keyword Arguments
    - `ncalc::Int=1000`                 : Number of sublayers to use for Love.jl

    # Returns
    - `power_prf::Array{prec,1}`        : Heating profile.
    - `k2_T::precc`                     : Complex Tidal k2 Lovenumber.
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.
    """
    function run_solid_strain( omega::Float64,
                        ecc::prec,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1};
                        ncalc::Int=2000
                        )::Tuple{Array{prec,1},precc,precc}

        # internal structure arrays.
        # first element is the innermost layer, last element is the outermost layer
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μc = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        # subdivide input layers such that we have ~ncalc in total
        rr = Love.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # get gravity at each layer
        g = Love.get_g(rr, ρ);

        # create grid
        Love.define_spherical_grid(res; n=2)

        # get y-functions
        M, y1_4 = Love.compute_M(rr, ρ, g, μc, κ)
        #   Tidal
        tidal_solution_T = Love.compute_y(rr, g, M, y1_4; load=false)
        #   Load
        tidal_solution_L = Love.compute_y(rr, g, M, y1_4; load=true)

        # get k2 tidal Love Number (complex-valued)
        k2_T = tidal_solution_T[5, end, end] - 1
        k2_L = tidal_solution_L[5, end, end] - 1
        
        # return zero for now
        k2_L = 0.

        # Get profile power output (W m-3), converted to W/kg
        (Eμ, Eκ) = Love.get_heating_profile(tidal_solution_T,
                               rr, ρ, g, μc, κ,
                               omega, ecc)

        Eμ_tot, _ = Eμ   # shear       (W), (W/m3)
        Eκ_tot, _ = Eκ   # compaction  (W), (W/m3)

        power_prf = Eμ_tot .+ Eκ_tot # Compute total volumetric heating (W/m3)

        return power_prf, k2_T, k2_L
    end


    """
        run_solid(rho, radius, visc, shear, bulk; ncalc=2000)

    Calculate k2 Lovenumbers in the solid.

    # Arguments
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `μ_profile::Array{precc,1}`       : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    
    # Keyword Arguments
    - `ncalc::Int=1000`                 : Number of sublayers to use for Love.jl

    # Returns
    - `k2_T::precc`                     : Complex Tidal k2 Lovenumber.
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.
    """
    function run_solid( rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1};
                        ncalc::Int=1000
                        )::Tuple{precc,precc}

        # internal structure arrays.
        # first element is the innermost layer, last element is the outermost layer
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μc = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        # subdivide input layers such that we have ~ncalc in total
        rr = Love.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # get gravity at each layer
        g = Love.get_g(rr, ρ);

        # create grid
        Love.define_spherical_grid(res; n=2)

        # get y-functions
        M, y1_4 = Love.compute_M(rr, ρ, g, μc, κ)
        #   Tidal
        tidal_solution_T = Love.compute_y(rr, g, M, y1_4; load=false)
        #   Load
        tidal_solution_L = Love.compute_y(rr, g, M, y1_4; load=true)

        # get k2 tidal Love Number (complex-valued)
        k2_T = tidal_solution_T[5, end, end] - 1
        k2_L = tidal_solution_L[5, end, end] - 1
        
        # return zero for now
        k2_L = 0.

        return k2_T, k2_L
    end
    

    """
        run_fluid(omega, rho, radius, ρ_ratio; n=2, sigma_R=1e-3)

    Calculate k2 Lovenumbers in the fluid.

    # Arguments
    - `omega::Float64`                  : Forcing frequency range.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `ρ_ratio::prec`                   : Density contrast between current (fluid) and lower (non-fluid) layer.
    
    # Keyword Arguments
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `sigma_R::Float64=1e-3`           : Rayleigh drag coefficient.

    # Returns
    - `k2_T::precc`                     : Complex Tidal k2 Lovenumber.
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.
    """
    function run_fluid( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        ρ_ratio::prec;
                        n::Int64=2,
                        sigma_R::Float64=1e-3
                        )::Tuple{precc,precc}

        # internal structure arrays
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)

        # surface properties
        R = maximum(r)  # Planet radius (m)
        g = G * sum(4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) .* ρ[1:end]) / R^2  # Surface gravity (m/s^2)

        # fluid magma ocean bottom and height
        r_b = r[1]
        H_magma = r[end] - r_b
             
        # get k2 Lovenumbers
        k2_T, k2_L = Fluid.compute_fluid_lovenumbers(omega, R, H_magma, g, ρ_ratio, n, sigma_R)

        return k2_T, k2_L

    end


    """
        radial_heating_profile(r, mu, Ptot)

    Compute a spherically symmetric volumetric heating profile H(r) [W/m^3]
    
    # Arguments
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.

    - `r::Vector`                       : Radial positions of layers, from core to surface.
    - `mu::Vector`                      : Complex shear modulus profile of the planet.
    - `Ptot::Real`                      : Totally dissipated power (W)

    # Returns
    - `H::Vector`                       : Heating profile.

    Heating is assumed proportional to Im(μ(r)) and normalized so the
    integral of H over the volume equals Ptot.
    """
    function radial_heating_profile(r::Vector, mu::Vector, Ptot::Real)::Vector

        # imaginary part determines dissipation strength
        w = imag.(mu)

        # avoid division by zero if purely elastic somewhere
        w .= w .+ (minimum(w) == 0 ? eps() : 0)

        # compute shell volumes (m^3)
        shell_volumes = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) 

        # total power for unnormalized profile
        P_unnorm = sum(w .* shell_volumes)

        # normalization factor
        norm = Ptot / P_unnorm

        # normalized volumetric profile (W/m^3)
        H = norm .* w

        return H
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
            ρ_l_mean = Fluid.mean_rho(ρ_l, vcat(r_l, r_b))
        end
        
        if length(ρ_s) == 0
            ρ_s_mean = 5000.0  # arbitrary high density
        elseif length(ρ_s) == 1
            ρ_s_mean = ρ_s[1]
        else
            ρ_s_mean = Fluid.mean_rho(ρ_s, vcat(r_s, r_c))
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
        k22_total = zeros(precc, N_sigma)

        for kk in 1:N_sigma
            σ = σ_range[kk]

            k2_T, k2_L = Fluid.compute_fluid_lovenumbers(σ, R, H_magma, g, ρ_ratio, n, sigma_R)

            k22_total[kk] = k2_T
        end

        k22_total = vec(k22_total)    # reshape(-1)

        # Build symmetric full spectrum for interpolation
        full_σ_range   = vcat(-σ_range,     reverse(σ_range))
        full_k22_total = vcat(-k22_total,   reverse(k22_total))

        imag_full_k22  = imag.(full_k22_total)

        # interpolation functions for imaginary parts (extrapolate outside)
        interp_full  = extrapolate(interpolate((full_σ_range,), imag_full_k22,
                                            Gridded(Linear())), Flat())

        # calculate tidal heating
        A_22k_e      = zeros(prec,  length(k_range))
        U_22k_e      = zeros(precc, length(k_range))
        P_T_k_total  = zeros(prec,  length(k_range))

        for (ikk, kk) in pairs(k_range)
            σ = 2*axial - kk*omega

            # Eq. 33
            A_22k_e[ikk] = sqrt(6π/5) * X_hansen[ikk]

            # Eq. 32
            U_22k_e[ikk] = (G*S_mass/sma) * (R/sma)^2 * A_22k_e[ikk]

            img_full_k22  = interp_full(σ)

            prefactor = 5 * R * σ / (8π*G)
            U2 = abs2(U_22k_e[ikk])

            P_T_k_total[ikk] = prefactor * img_full_k22  * U2
        end

        # Total tidal heating (negative sum)
        P_tidal_total = -sum(P_T_k_total) # W

        # Get power profile
        power_prf = zeros(prec, length(ρ))
        power_prf[mask_l] .= Fluid.heat_profile(P_tidal_total, ρ_l, r_l, r_b)
        
        # Get k2 lovenumber at forcing frequency used by Love.jl
        k2_total = interp_full(2*axial - omega)
        
        return convert(Vector{Float64}, power_prf), convert(Float64, P_tidal_total), convert(Float64, k2_total)

    end

end