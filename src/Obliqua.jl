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
    using MultiFloats
    using AssociatedLegendrePolynomials

    using NCDatasets

    # Include local jl files
    include("solid0d.jl")
    include("solid1d.jl")
    include("solid1d_mush.jl")
    include("solid1d_relax.jl")
    include("solid1d_mush_relax.jl")
    include("fluid0d.jl")
    include("Hansen.jl")
    include("load.jl")
    include("plotting.jl")

    # Import submodules
    import .solid0d
    import .solid1d
    import .solid1d_mush
    import .solid1d_relax
    import .solid1d_mush_relax
    import .fluid0d
    import .Hansen
    import .load
    import .plotting

    # Export submodules (mostly for autodoc purposes)
    export solid0d
    export solid1d
    export solid1d_mush
    export solid1d_relax
    export solid1d_mush_relax
    export fluid0d
    export Hansen
    export load
    export plotting

    export run_tides

    const ROOT_DIR::String = abspath(dirname(abspath(@__FILE__)), "../")
    const RES_DIR::String  = joinpath(ROOT_DIR,"res/")
    const OUT_DIR::String  = joinpath(ROOT_DIR,"out/")

    prec  = BigFloat
    precc = Complex{BigFloat}

    # prec  = Float64x4
    # precc = Complex{Float64x4}

    # prec  = Float64
    # precc = Complex{Float64}

    const AU::prec      = prec(1.495978707e11)   # m
    const G::prec       = prec(6.6743e-11)       # m^3 kg^-1 s^-2
    const M_Earth::prec = prec(5.9724e24)        # kg

    const res::Float64 = 20.0                    # angular resolution in degrees


    """
        Create a logger object and return it.

    Arguments:
    - `outpath::String`                 : Output file (empty to disable file logging).

    Optional arguments:
    - `to_term::Bool`                   : Log to terminal?

    Returns:
    - `logger_both`                     : Logger object.
    """
    function make_logger(outpath::String; to_term::Bool=true)

        # Formatting
        color::Int = 39
        level::String = "UNSET"
        term_io::IO = stdout

        # Setup file logger
        to_file::Bool = !isempty(outpath)
        if to_file
            # remove old file
            if isfile(outpath)
                rm(outpath)
            end

            # configure
            logger_file = FormatLogger(outpath; append=false) do io, args
                if args.level == LoggingExtras.Info
                    level = "INFO"
                elseif args.level == LoggingExtras.Warn
                    level = "WARN"
                elseif args.level == LoggingExtras.Debug
                    level = "DEBUG"
                elseif args.level == LoggingExtras.Error
                    level = "ERROR"
                end
                @printf(io, "[ %-5s ] %s \n", level, args.message)
            end;
        end

        # Setup terminal logger
        if to_term
            logger_term = FormatLogger() do io, args
                if args.level == LoggingExtras.Info
                    color = 32
                    level = "INFO"
                elseif args.level == LoggingExtras.Warn
                    color = 93
                    level = "WARN"
                elseif args.level == LoggingExtras.Debug
                    color = 96
                    level = "DEBUG"
                elseif args.level == LoggingExtras.Error
                    color = 91
                    term_io = stderr
                    level = "ERROR"
                end
                # Set color, set bold, print level, unset bold, unset color, message
                @printf(term_io, "[\033[%dm\033[1m %-5s \033[21m\033[0m] %s\n",
                                    color, level, args.message)
            end;
        end

        # Return logger object
        if to_file && to_term
            return TeeLogger(logger_file, logger_term)
        elseif to_term
            return logger_term
        elseif to_file
            return logger_file
        else
            println(stderr, "Warning: using NullLogger to log all messages")
            return NullLogger()
        end
    end


    """
        Setup terminal logging and file logging.

    Arguments:
    - `outpath::String`                 : Output file (empty to disable file logging)
    - `verbosity::Int`                  : Verbosity (0: silent, 1: normal, 2: debug)
    """
    function setup_logging(outpath::String, verbosity::Int)

        # If silent
        if verbosity==0
            global_logger(MinLevelLogger(current_logger(), Logging.Error))
            return nothing
        end

        # Make the logger
        logger_both = make_logger(outpath)
        global_logger(logger_both)

        # Disable debug
        if verbosity == 1
            disable_logging(Logging.Debug)
        end

        return nothing
    end


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
                        phi::Array{prec,1},
                        cfg::Dict
                        )::Tuple{Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}}
      
        # Read configuration options from dict
        @info "Using configuration '$(cfg["title"])'"
        
        # Check that config has these always-required keys
        req_keys = Dict(
            "params.out" => ["path"],
            "orbit.obliqua" => [
                "min_frac","visc_l","visc_lus","visc_sus",
                "n","m","spectrum", "material",
                "module_solid", "module_fluid", "module_mushy"
            ],
            "orbit.obliqua.solid" => [
                "ncalc"
            ],
            "orbit.obliqua.fluid" => [
                "sigma_R"
            ],
            "struct" => ["mass_tot","core_density"]
        )
        for (section, keys) in req_keys
            path = split(section, ".")
            node = cfg

            # walk down the nested Dict
            for p in path
                if !(node isa AbstractDict)
                    @error "Config: expected Dict at `$(join(path, "."))`, got $(typeof(node))"
                    break
                end
                if !haskey(node, p)
                    @error "Config: missing required section `$(join(path, "."))`"
                    break
                end
                node = node[p]
            end

            @debug "Tested config section `$(join(path, "."))`"

            # check required keys at this level
            for k in keys
                @debug "Testing config key `$(join(path, "."))::$k`"
                if !haskey(node, k)
                    @error "Config: missing required key `$(join(path, "."))::$k`"
                end
            end
        end

        @debug "Config file validation complete."

        # collection of config params 
        min_frac = cfg["orbit"]["obliqua"]["min_frac"]

        visc_l   = cfg["orbit"]["obliqua"]["visc_l"]
        visc_lus = cfg["orbit"]["obliqua"]["visc_lus"]
        visc_s   = cfg["orbit"]["obliqua"]["visc_s"]
        visc_sus = cfg["orbit"]["obliqua"]["visc_sus"]

        n        = cfg["orbit"]["obliqua"]["n"]
        m        = cfg["orbit"]["obliqua"]["m"]

        spectrum = cfg["orbit"]["obliqua"]["spectrum"]
        if spectrum == "adaptive"
            s_min    = cfg["orbit"]["obliqua"]["s_min"]
            s_max    = cfg["orbit"]["obliqua"]["s_max"]
            s_min    = nothing_if_none(s_min)
            s_max    = nothing_if_none(s_max)
        elseif spectrum == "full"
            N_σ      = cfg["orbit"]["obliqua"]["N_sigma"]
            p_min    = cfg["orbit"]["obliqua"]["p_min"]
            p_max    = cfg["orbit"]["obliqua"]["p_max"]
        end

        material = cfg["orbit"]["obliqua"]["material"]
        alpha    = cfg["orbit"]["obliqua"]["alpha"]

        module_solid = cfg["orbit"]["obliqua"]["module_solid"]
        ncalc        = cfg["orbit"]["obliqua"]["solid"]["ncalc"]
        dr_min       = cfg["orbit"]["obliqua"]["solid"]["dr_min"]
        dr_max       = cfg["orbit"]["obliqua"]["solid"]["dr_max"]
        core         = cfg["orbit"]["obliqua"]["solid"]["core"]
        bulk_l       = cfg["orbit"]["obliqua"]["solid"]["bulk_l"]
        permea       = cfg["orbit"]["obliqua"]["solid"]["permea"]
        porosity_thresh = cfg["orbit"]["obliqua"]["solid"]["porosity_thresh"]

        module_fluid = cfg["orbit"]["obliqua"]["module_fluid"]
        sigma_R      = cfg["orbit"]["obliqua"]["fluid"]["sigma_R"]
        sigma_R_inf  = cfg["orbit"]["obliqua"]["fluid"]["sigma_R_inf"]
        sigma_R_prf  = cfg["orbit"]["obliqua"]["fluid"]["sigma_R_prf"]
        H_R          = cfg["orbit"]["obliqua"]["fluid"]["H_R"]
        efficiency   = cfg["orbit"]["obliqua"]["fluid"]["efficiency"]

        module_mushy = cfg["orbit"]["obliqua"]["module_mushy"]
        if module_mushy == "interp"
            b_width  = cfg["orbit"]["obliqua"]["mushy"]["b_width"]
            t_width  = cfg["orbit"]["obliqua"]["mushy"]["t_width"]
        end

        mass_tot = cfg["struct"]["mass_tot"]*M_Earth

        # convert "none" to nothing
        module_solid = nothing_if_none(module_solid)
        module_fluid = nothing_if_none(module_fluid)
        module_mushy = nothing_if_none(module_mushy)

        # convert interior profiles to BigFloat                 
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μ = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)
        ϕ = convert(Vector{prec}, phi)

        # number of layers
        N_layers = length(r)-1

        # liquidus and solidus viscosity
        η_l = visc_lus
        η_s = visc_sus

        # get smoothed out region masks and segments
        mask_l, mask_s, mask_c, is_seg, segments = get_layers(r, η, η_l, η_s; min_frac)

        # check if CMB is at the bottom of the mantle
        if any(r[1] .>= r[2:end])
            throw("CMB radius not at bottom of mantle, did you properly order the interior arrays?")
        end

        # get core mass from core density and radius
        μ_core = convert(prec, cfg["struct"]["core_shear"])
        κ_core = convert(prec, cfg["struct"]["core_bulk"])
        ρ_core = convert(prec, cfg["struct"]["core_density"])
        m_core = (4/3)*π*r[1]^3*ρ_core

        # find planet radius (m)
        R = maximum(r)

        # shell masses
        dv = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) 
        dm = dv .* ρ

        # cumulative enclosed mass, including core mass
        M_enc = cumsum(dm) .+ m_core

        # gravity at each layer radius
        g = G .* M_enc ./ r[2:end].^2

        # collect tidal modes (n, m) 
        nm  = [(n_i, m_i) for n_i in n, m_i in m if n_i >= m_i]
        # initiate container for (n, m, k) combinations
        nmk = Vector{Tuple{Int,Int,Int}}()

        # Initialize containers that will hold data per (n, m, k) mode
        σ_range  = Vector{Float64}()
        X_hansen = Vector{Float64}()
        μc       = Matrix{ComplexF64} # layer x frequency

        # Output containers that will hold data per (n, m, k) mode
        knms_T    = Matrix{ComplexF64}
        knms_L    = Matrix{ComplexF64}
        prf_total = Matrix{Float64}
        map_total = Array{Float64, 4}

        # orbital and axial frequencies
        if spectrum == "adaptive"            

            for (n_i, m_i) in nm
                # avoid overwriting inputs by using local scope variables
                s_min_i = s_min === nothing ? nothing : s_min
                s_max_i = s_max === nothing ? nothing : s_max

                # get s range for proper convergence for given eccentricity
                s_min_ecc, s_max_ecc = Hansen.get_k_range(ecc, n_i, m_i)

                if s_min_i === nothing 
                    s_min_i = s_min_ecc
                elseif s_min_i > s_min_ecc
                    @warn "Provided s_min=$s_min_i is larger than estimated s_min=$s_min_ecc for eccentricity $(round(Float64(ecc), digits=2))."
                end
                if s_max_i === nothing
                    s_max_i = s_max_ecc
                elseif s_max_i < s_max_ecc
                    @warn "Provided s_max=$s_max_i is smaller than estimated s_max=$s_max_ecc for eccentricity $(round(Float64(ecc), digits=2))."
                end

                @info "Using adaptive spectrum with s range [$s_min_i, $s_max_i] for eccentricity $(round(Float64(ecc), digits = 2)) and tide (n, m) = ($n_i, $m_i)."

                # build s range for this mode
                s_range = collect(s_min_i:s_max_i)

                # get hansen coefficients for this specific mode
                _, X_hansen_i = Hansen.get_hansen(ecc, n_i, m_i, s_min_i, s_max_i)
                
                # calculate forcing frequencies only for region of interest
                σ_range_i = m_i .* axial .- s_range .* omega
                σ_range_i = Float64.(σ_range_i)

                for i in 1:length(s_range)
                    push!(nmk, (n_i, m_i, s_range[i]))
                    push!(X_hansen, X_hansen_i[i])
                    push!(σ_range, σ_range_i[i])
                end
            end

            N_σ = length(σ_range)

        elseif spectrum == "full"
            # Maintain backward compatibility for the fallback plotting spectrum
            t_range = 10 .^ range(p_min, stop=p_max, length=N_σ)       
            σ_range_i = 2π ./ (t_range .* 1e3 .* 365.25 .* 24 .* 3600)    
            σ_range_i = reshape(σ_range_i, :)

            # set particular nmk for full spectrum
            for i in 1:N_σ
                push!(nmk, (nm[1][1], nm[1][2], 1)) 
                push!(σ_range, σ_range_i[i])
            end
            @info "Using (n, m, k) = ($(nm[1][1]), $(nm[1][2]), 1) for full spectrum."
        end

        # get frequency dependent complex shear modulus per mode
        μc = complex_mu(σ_range, μ, η; material=material, α=alpha)

        # allocate outputs for this specific mode's frequency count
        # initiate forcing frequency dependent k Love numbers (one spectrum for each segment)
        knms_T = zeros(ComplexF64, N_σ, length(segments))
        knms_L = zeros(ComplexF64, N_σ, length(segments))
        # initiate forcing frequency dependent heating profile
        prf_total = zeros(Float64, N_σ, N_layers)
        map_total = zeros(Float64, N_σ, length(segments), length(0:res:180), length(0:res:360-0.001))

        # core density for bottom boundary
        ρ_mean_lower = ρ_core
        # Rayleigh drag efficiency at fluid core
        efficiency_seg = efficiency

        # interpolation activity boolean, used with the interp function for mushy layers
        interp_active = false
        interp_previous = false

        # loop over segments, starting at CMB
        for (iseg, seg) in pairs(segments)
            # check if interpolation is active
            if interp_active
                interp_previous = true
                interp_active = false
            end

            # get start and stop index for segment
            i_start, i_end = is_seg[iseg]

            # perform slices
            r_seg  = r[i_start-1:i_end]
            ρ_seg  = ρ[i_start:i_end]
            η_seg  = η[i_start:i_end]                              
            μc_seg = μc[i_start:i_end, :] 
            κ_seg  = κ[i_start:i_end]
            ϕ_seg  = ϕ[i_start:i_end]
            g_seg  = g[i_start:i_end]

            # mean density in current segment
            if length(ρ_seg) == 1
                ρ_mean = ρ_seg[1]
            else
                ρ_mean = fluid0d.mean_rho(ρ_seg, r_seg)
            end

            # density ratio
            ρ_ratio = ρ_mean / ρ_mean_lower

            # get k2 spectrum for segment
            for iss in 1:N_σ
                # specify mode
                n_i, m_i, s_i = nmk[iss]

                # specify forcing frequency
                σ = σ_range[iss]

                # if forcing frequency is zero, then skip to next frequency (no heating)
                iszero(σ) && continue

                # if segment is solid
                if seg == "solid"
                    # don't model solid tides
                    if module_solid===nothing
                        knms_T[iss, iseg], knms_L[iss, iseg] = 0., 0.
                    # 0D interior but no heating profile in segment 
                    elseif module_solid=="solid0d"
                        knms_T[iss, iseg], knms_L[iss, iseg] = run_solid0d( 
                            μc_seg[:, iss],
                            r_seg,
                            mass_tot;
                            n=n_i
                        )
                    # elseif 1D interior and heating profile from strain tensor
                    elseif module_solid=="solid1d"
                        prf_total[iss, i_start:i_end], knms_T[iss, iseg], knms_L[iss, iseg] = run_solid1d( 
                            σ, ρ_seg,
                            r_seg, η_seg,                               
                            μc_seg[:, iss], 
                            κ_seg, R, 
                            m_core, ρ_core, 
                            μ_core, κ_core;
                            ncalc=ncalc, n=n_i, m=m_i
                        )
                    # elseif 1D interior and heating profile from strain tensor
                    elseif module_solid=="solid1d-relax"
                        prf_total[iss, i_start:i_end], map_total[iss, iseg, :, :], knms_T[iss, iseg], knms_L[iss, iseg] = run_solid1d_relax( 
                            σ, ρ_seg, r_seg,
                            η_seg, μc_seg[:, iss], 
                            κ_seg, R, 
                            m_core, ρ_core, 
                            μ_core, κ_core;
                            dr_min=dr_min, dr_max=dr_max, 
                            n=n_i, m=m_i, core=core
                        )
                    # elseif 1D interior with mush interface and heating profile from strain tensor
                    elseif module_solid=="solid1d-mush"
                        prf_total[iss, i_start:i_end], knms_T[iss, iseg], knms_L[iss, iseg] = run_solid1d_mush( 
                            σ, ρ_seg, r_seg,
                            η_seg, μc_seg[:, iss], 
                            κ_seg, ϕ_seg, R, 
                            m_core, ρ_core, 
                            μ_core, κ_core;
                            ncalc=ncalc, n=n_i, m=m_i, core=core, visc_l=visc_l, bulk_l=bulk_l,
                            permea=permea, porosity_thresh=porosity_thresh
                        )
                    elseif module_solid=="solid1d-mush-relax"
                        # prf_total[iss, i_start:i_end], map_total[iss, iseg, :, :], knms_T[iss, iseg], knms_L[iss, iseg] = run_solid1d_mush_relax( 
                        prf_total[iss, i_start:i_end], knms_T[iss, iseg], knms_L[iss, iseg] = run_solid1d_mush_relax( 
                            σ, ρ_seg, r_seg,
                            η_seg, μc_seg[:, iss], 
                            κ_seg, ϕ_seg, R, 
                            m_core, ρ_core,
                            μ_core, κ_core;
                            dr_min=dr_min, dr_max=dr_max, 
                            n=n_i, m=m_i, core=core, visc_l=visc_l, bulk_l=bulk_l,
                            permea=permea, porosity_thresh=porosity_thresh
                        )
                    else
                        throw("No compatible solid tides module: $module_solid.")
                    end

                # if segment is fluid
                elseif seg == "fluid"
                    # don't model fluid tides
                    if module_fluid===nothing
                        knms_T[iss, iseg], knms_L[iss, iseg] = 0., 0.
                    # 0D interior but no heating profile in segment 
                    elseif module_fluid=="fluid0d"
                        knms_T[iss, iseg], knms_L[iss, iseg] = run_fluid0d(
                            σ, ρ_seg, 
                            r_seg, ρ_ratio;
                            n=n_i, 
                            sigma_R=sigma_R
                        ) 
                    # elseif 1D interior and heating profile from density-contrast/Rayleigh-drag
                    elseif module_fluid=="fluid1d"
                        prf_total[iss, i_start:i_end], knms_T[iss, iseg], knms_L[iss, iseg] = run_fluid1d(
                            σ, ρ_seg, r_seg, 
                            g_seg, ρ_ratio,
                            S_mass, sma, R; n=n_i,
                            sigma_R=sigma_R,
                            sigma_inf=sigma_R_inf,
                            sigma_R_prf=sigma_R_prf,
                            H_R=H_R, efficiency=efficiency_seg
                        )
                    else
                        throw("No compatible fluid tides module: $module_fluid.")
                    end
                
                # if segment is mush
                elseif seg == "mush"
                    # don't model mush tides
                    if module_mushy===nothing
                        knms_T[iss, iseg], knms_L[iss, iseg] = 0., 0.
                    elseif module_mushy=="interp"
                        # turn on interpolation mode
                        interp_active = true

                        # get heating in previous layers
                        if i_start > 1
                            i_sp, i_ep = is_seg[iseg-1]
                            P_b = prf_total[iss, i_ep]

                            # first solve heating spectrum for lower interface
                            prf_total[iss, i_start:i_end], knms_T[iss, iseg], knms_L[iss, iseg] = run_interp(
                                σ, r_seg, R, 0., P_b;
                                t_width=t_width, b_width=b_width
                            )
                        end
                        # Then model next segment to get heating at the upper interface P_t
                    else
                        throw("No compatible mush tides module: $module_mushy.")
                    end

                # if segment is ice
                elseif seg == "ice"
                    # calculate ice tides in ice region 
                    knms_T[iss, iseg], knms_L[iss, iseg] = 0., 0. # no expression for this yet
                    @warn "Ice layers are currently not supported. Skipping this segment..."
                    
                # if segment is water
                elseif seg == "water"
                    # calculate water tides in water region 
                    knms_T[iss, iseg], knms_L[iss, iseg] = 0., 0. # no expression for this yet
                    @warn "Water layers are currently not supported. Skipping this segment..."    
                end

                if interp_previous
                    # Solve heating spectrum for upper interface in previous segment
                    P_t = prf_total[iss, i_start] # get heating in bottom layer of current segment
                    i_sp, i_ep = is_seg[iseg-1]
                    Δprf, ΔkT, ΔkL = run_interp(
                        σ, r[i_sp-1:i_ep], R, P_t, 0.;
                        t_width=t_width, b_width=b_width
                    )
                    prf_total[iss, i_sp:i_ep] .+= Δprf
                    knms_T[iss, iseg-1]        += ΔkT
                    knms_L[iss, iseg-1]        += ΔkL
                end

                # repeat for all probe forcing frequencies
            end
        
            # update previous segment mean density before moving to next segment
            ρ_mean_lower = ρ_mean
            # update Rayleigh drag efficiency away from core
            efficiency_seg = 1.

            # turn off interpolator after completion
            interp_previous = false

            # step to next segment
        end  

        # initialize total k2 with the contribution from the top layer
        knms_total = copy(knms_T[:, end])  

        # loop from top (surface) to just above CMB
        # This implementation is currently only valid for CMB -> Solid -> Fluid -> Surface layering, 
        # and would need to be adapted for more complex layering (e.g., fluid under solid).
        for iseg in reverse(1:length(segments)-1)

            # get the indices for the current segment and the previous segment
            # current segemnt
            i_start, i_end = is_seg[iseg]
            # previous segment
            i_start_ini, i_end_ini = is_seg[iseg+1]

            load = segments[iseg] == "solid" && segments[iseg+1] == "fluid" ? true : false

            if load
                for i in 1:N_σ
                    # calculate the load distribution factor
                    factor = 1 + (imag(knms_L[i, iseg] * knms_total[i])) / (imag(knms_T[i, iseg]) + imag(knms_total[i]) + 1e-40) # add small number to avoid divide by zero

                    # correct heating profile with the load distribution factor
                    prf_total[i, i_start:i_end] .*= factor
                    prf_total[i, i_start_ini:i_end_ini] .*= factor

                    # update the total Lovenumber
                    knms_total[i] = knms_T[i, iseg] + (1.0 + knms_L[i, iseg]) * knms_total[i]
                end
            else
                for i in 1:N_σ
                    # update the total Lovenumber
                    knms_total[i] = knms_T[i, iseg] + knms_total[i]
                end
            end
        end

        # extract imaginary part of complex global k2 spectrum
        imag_k2 = .-imag.(knms_total)

        # if using full spectrum, then calculate heating profile and bulk heating at each 
        # frequency and return the full spectrum of heating profiles and bulk heating for plotting
        # the code assumes s=1 to find a solution for the Hansen coefficients and normalization.
        if spectrum == "full"
            # get (n, m, k) mode
            n, m, k = nmk[1]

            # Hansen coefficient at s=1
            _, X = Hansen.get_hansen(ecc, n, m, 1, 1)

            A = 2 * sqrt(4π * factorial(n-m) / ((2*n+1) * factorial(n+m))) * Plm(n, m, 0.) * X
                        
            U = (G*S_mass/sma) * (R/sma)^n * A

            prefactor = (2*n + 1) * R / (8π*G) .* σ_range

            U2 = abs2.(U)

            # return bulk heating at each frequency
            P_T_1_blk = prefactor .* imag_k2 .* U2

            # return power profile at each frequency
            P_T_1_prf = zeros(Float64, N_σ, length(shear))
            P_T_1_map = zeros(Float64,  N_σ, length(collect(0:res:180)), length(collect(0:res:360-0.001)))

            for iss in 1:N_σ
                unorm_prf = prf_total[iss, :]
                unorm_map = sum(map_total[iss, :, :, :], dims=1) # sum over layers to get surface map

                P_T_1_prf[iss, :]    = Float64.(unorm_prf .* U2)
                P_T_1_map[iss, :, :] = Float64.(unorm_map .* U2)

                # Integrate the radial profile over the volume for this frequency
                # Assuming 'dv' is the differential volume element array
                integrated_profile_power = sum(dv .* P_T_1_prf[iss, :])
                
                # Calculate ratio
                ratio = P_T_1_blk[iss] / integrated_profile_power
                
                # Print formatted results
                @debug("Freq Index: %d | Ratio: %.6f\n", iss, ratio)
            end          

            P_T_blk = Float64(sum(P_T_1_blk)) # W

            # get radial heating profile W/m^3
            P_T_prf = Float64.([sum(P_T_1_prf[:,j]) for j in 1:size(P_T_1_prf,2)])
            P_T_map = dropdims(sum(P_T_1_map, dims=1), dims=1) # sum over frequencies to get total surface map

            # determine the total heat input from heating profile
            P_T_prf_blk = Float64(sum(dv .* P_T_prf))

            # define data file path
            datafile_path = joinpath(OUT_DIR, "obliqua_data.nc")

            # store results in netcdf file
            data_to_nc(
                    nmk, is_seg, segments, knms_total, knms_T, knms_L, 
                    σ_range, P_T_1_prf, P_T_prf, P_T_blk, 
                    P_T_prf_blk, P_T_1_map, P_T_map, datafile_path
                )

            @info("Expected bulk heating: $P_T_blk")
            @info("Obtained bulk heating: $P_T_prf_blk")

            P_T_prf ./ ρ # convert to mass heating rate (W/kg)

            # return full Imk2 spectrum for plotting
            return Float64.(P_T_prf), Float64.(P_T_blk), Float64.(σ_range), Float64.(imag_k2)
        end
            
        # alternatively, if using adaptive spectrum, then calculate heating profile and
        # bulk heating only at the frequencies of interest and return the spectrum of heating 
        # profiles and bulk heating for plotting. Here the s_range is properly chosen to find
        # the true solution for the Hansen coefficients and normalization given the orbit.

        # calculate tidal heating
        # initialize frequency dependent quentities
        A_nms_e   = zeros(Float64,  N_σ)
        U_nms_e   = zeros(ComplexF64, N_σ)

        # initialize frequency dependent total heating
        P_T_s_blk = zeros(Float64,  N_σ)

        # initialize frequency dependent heating profile
        P_T_s_prf = zeros(Float64,  N_σ, length(shear))
        P_T_s_map = zeros(Float64,  N_σ, length(collect(0:res:180)), length(collect(0:res:360-0.001)))

        # loop over tidal modes 
        for iss in 1:N_σ
            # specify mode
            n_i, m_i, s_i = nmk[iss]

            # calculate physical forcing frequency
            σ = m_i*axial - s_i*omega

            # if forcing frequency is zero, then skip to next frequency (no heating)
            iszero(σ) && continue

            # calculate coefficients
            a = (abs(m_i) == 0 && s_i == 0) ? 1.0 : 0.0
            b = (abs(m_i) == 0 && s_i < 0)  ? 1.0 : 0.0
            
            A_nms_e[iss] = (2. - a) * (1. - b) * sqrt(4π * factorial(n_i-m_i) / ((2*n_i+1) * factorial(n_i+m_i))) * Plm.(n_i, m_i, 0.) * X_hansen[iss]
                      
            U_nms_e[iss] = (G*S_mass/sma) * (R/sma)^n_i * A_nms_e[iss]

            # get imaginary part of complex k2 love number from global spectrum at forcing frequency
            img_full_knm = imag_k2[iss] 

            # calculate prefactor and total availible heat
            prefactor = (2*n_i+1) * R * σ / (8π*G)
            U2 = abs2(U_nms_e[iss])

            # calculate total heat input at forcing frequency
            P_T_s_blk[iss] = prefactor * img_full_knm  * U2
            
            # get global heating profile at forcing frequency
            unorm_prf = prf_total[iss, :] 
            unorm_map = sum(map_total[iss, :, :, :], dims=1) # sum over layers to get surface map

            # Hansen renorm
            P_T_s_prf[iss, :] = unorm_prf .* U2
            P_T_s_map[iss, :, :] = unorm_map .* U2

            # For debugging purposes log the ratio between expected and radially integrated heating
            ratio = P_T_s_blk[iss] / sum(dv .* P_T_s_prf[iss, :])
            @debug "Forcing Frequency: $σ, Global heating ratio (blk/prf): $ratio"

        end

        # total tidal heating
        P_T_blk = sum(P_T_s_blk) # W

        # get radial heating profile W/m^3
        P_T_prf = [sum(P_T_s_prf[:,j]) for j in 1:size(P_T_s_prf,2)]
        P_T_map = dropdims(sum(P_T_s_map, dims=1), dims=1) # sum over frequencies to get total surface map

        # determine the total heat input from heating profile
        P_T_prf_blk = Float64.(sum(dv .* P_T_prf))

        # define data file path
        datafile_path = joinpath(OUT_DIR, "obliqua_data.nc")

        # store results in netcdf file
        data_to_nc(
                nmk, is_seg, segments, knms_total, knms_T, knms_L, 
                σ_range, P_T_s_prf, P_T_prf, P_T_blk, 
                P_T_prf_blk, P_T_s_map, P_T_map, datafile_path
            )

        @info("Expected bulk heating: $P_T_blk")
        @info("Obtained bulk heating: $P_T_prf_blk")

        P_T_prf ./ ρ # convert to mass heating rate (W/kg)

        # convert everything to Float64
        return Float64.(P_T_prf), Float64(P_T_blk), Float64.(σ_range), Float64.(imag_k2)

    end 


    """Convert 'none' string into nothing literal."""
    function nothing_if_none(val)
        if val == "none"
            return nothing
        else 
            return val 
        end
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
        # mask_s = (η_s .< η .< 1e13)
        mask_s = η_s .< η 

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
        run_solid0d(μc, radius, mass_tot; n=2)

    Calculate k2 Lovenumbers in the 0D solid.

    # Arguments
    - `μc::Array{precc,1}`              : Forcing frequency range.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `mass_tot::Float64`               : Total mass of planet.
    
    # Keyword Arguments
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).

    # Returns
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_solid0d( μc::Array{precc,1},
                        radius::Array{prec,1},
                        mass_tot::prec;
                        n::Int64=2
                        )::Tuple{ComplexF64,ComplexF64}

        # internal structure arrays
        μc = convert(Vector{precc}, μc)
        r  = convert(Vector{prec}, radius)

        # get mean complex shear modulus in segment
        μc_mean = solid0d.mean_cmu(μc, r)

        # surface properties
        R = maximum(r)  # Segment outer radius (m)

        # get k2 Lovenumbers
        k2_T, k2_L = solid0d.compute_solid_lovenumbers(μc_mean, mass_tot, R, n)

        # return zero for now
        k2_L = 0.

        return ComplexF64(k2_T), ComplexF64(k2_L)

    end


    """
        run_solid1d(omega, rho, radius, visc, shear, bulk, R, m_core, ρ_core, μ_core, κ_core; ncalc=2000, n=2, m=2, core="liquid")

    Use 1D solid tides model to calculate k2 Lovenumbers, and compute 1D heating profile from strain tensor.
    This method ignores inertia effects, since they break the numerical stability.

    # Arguments
    - `omega::prec`                     : Forcing frequency range.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `μ_profile::Array{precc,1}`       : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    - `R::prec`                         : Planet radius.
    - `m_core::prec`                    : Core mass.
    - `ρ_core::prec`                    : Core density.
    - `μ_core::prec`                    : Core shear modulus.
    - `κ_core::prec`                    : Core bulk modulus.
    
    # Keyword Arguments
    - `ncalc::Int=2000`                 : Number of sublayers.
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `m::Int=2`                        : Harmonic of the true anomaly. m=2 corresponds to the semidiurnal tide, m=1 diurnal tide.
    - `core::String="liquid"`           : Core state, either "liquid" or "solid".

    # Returns
    - `power_prf::Array{Float64,1}`     : Heating profile.
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_solid1d( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1},
                        R::prec,
                        m_core::prec,
                        ρ_core::prec,
                        μ_core::prec,
                        κ_core::prec;
                        ncalc::Int=2000,
                        n::Int=2,
                        m::Int=2,
                        core::String="liquid"
                        )::Tuple{Array{Float64,1},ComplexF64,ComplexF64}

        # internal structure arrays.
        # first element is the innermost layer, last element is the outermost layer
        omega = prec(omega)
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μc = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        # subdivide input layers such that we have ~ncalc in total
        rr = solid1d.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # get gravity at each layer
        g = solid1d.get_g(rr, ρ, m_core)

        # create grid
        solid1d.define_spherical_grid(res, n, m)

        # get y-functions
        M, y1_4 = solid1d.compute_M(omega, rr, ρ, g, μc, κ, n, ρ_core, μ_core, κ_core; core=core)
        #   Tidal
        tidal_solution_T = solid1d.compute_y(rr, g, M, y1_4, n; load=false)
        #   Load
        tidal_solution_L = solid1d.compute_y(rr, g, M, y1_4, n; load=true)

        # get k2 tidal Love Number (complex-valued)
        k2_T = tidal_solution_T[5, end, end] - 1
        k2_L = tidal_solution_L[5, end, end] - 1
        
        # Get profile power output (W m-3), converted to W/kg
        Eμ, Eκ = solid1d.get_heating_profile(tidal_solution_T, rr, ρ, g, μc, κ, n, omega; lay=nothing)

        Eμ_tot, _ = Eμ   # shear       (W), (W/m3)
        Eκ_tot, _ = Eκ   # compaction  (W), (W/m3)

        # Renormalization factor
        power_prf = (Eμ_tot .+ Eκ_tot) .* (R ./ maximum(r)) # Compute total volumetric heating (W/m3)

        return Float64.(power_prf), ComplexF64(k2_T), ComplexF64(k2_L)
    end
    

    """
        run_solid1d_relax(omega, rho, radius, visc, shear, bulk, R, m_core, ρ_core, μ_core, κ_core; dr_min=300, dr_max=3000, n=2, m=2, core="liquid")

    Use 1D solid tides model with relaxation method to calculate k2 Lovenumbers, and compute 1D heating profile from strain tensor.
    This method includes inertia effects, but is more computationally expensive. 
    
    # Arguments
    - `omega::prec`                     : Forcing frequency range.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `shear::Array{precc,1}`           : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    - `R::prec`                         : Planet radius.
    - `m_core::prec`                    : Core mass.
    - `ρ_core::prec`                    : Core density.
    - `μ_core::prec`                    : Core shear modulus.
    - `κ_core::prec`                    : Core bulk modulus.

    # Keyword Arguments
    - `dr_min::Int=300`                 : Minimum layer thickness in m.
    - `dr_max::Int=3000`                : Maximum layer thickness in m.
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `m::Int=2`                        : Harmonic of the true anomaly. m=2 corresponds to the semidiurnal tide, m=1 diurnal tide.
    - `core::String="liquid"`           : Core state, either "liquid", "solid", or "inertial".

    # Returns
    - `power_prf::Array{Float64,1}`     : Heating profile.
    - `power_map::Array{Float64,4}`     : Heating map (frequency, colatitude, longitude).
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_solid1d_relax( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1},
                        R::prec,
                        m_core::prec,
                        ρ_core::prec,
                        μ_core::prec,
                        κ_core::prec;
                        dr_min::Int=300,
                        dr_max::Int=3000,
                        n::Int=2,
                        m::Int=2,
                        core::String="liquid"
                        )::Tuple{Array{Float64,1},Matrix{Float64},ComplexF64,ComplexF64}

        # convert inputs
        omega = prec(omega)
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μ = convert(Vector{precc}, shear)
        κ = convert(Vector{prec}, bulk)

        # resample profiles onto new grid
        r_grid, ρ, η, μ, κ, g, M_tot = solid1d_relax.resample_profiles(r, ρ, η, μ, κ, m_core, dr_min, dr_max)

        # use cell centers
        r_centers = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])

        # define angular grid
        SphericalGrid = solid1d_relax.define_spherical_grid(res, n, m)

        # solve y functions across grid
        y_t, y_l = solid1d_relax.compute_y(r_centers, ρ, g, μ, κ, omega, n, ρ_core, μ_core, κ_core, M_tot; core=core)

        # for debugging: plot y-function relaxation solution
        # plotting.plot_relaxation_solution(y_t, r_centers, 
        #         filename="$OUT_DIR/relaxation_solution.png")

        # Love numbers
        k2_T = y_t[5, end] - 1
        k2_L = y_l[5, 1]   - 1

        # heating profile
        Eμ_tot, Eκ_tot = solid1d_relax.get_heating_profile(
            y_t, r_grid, ρ, g, μ, κ, n, omega, SphericalGrid
        )

        Eμ_map, Eκ_map = solid1d_mush_relax.get_heating_map(
            y_t, r_grid, ρ, g, μ, κ, n, omega, SphericalGrid
        )

        power_prf = abs.(Eμ_tot .+ Eκ_tot) 
        power_map = abs.(Eμ_map .+ Eκ_map) 

        # interpolate from grid back to original radius points 
        itp = linear_interpolation(r_centers, power_prf, extrapolation_bc=Line())

        # original centers
        r_orig_centers = 0.5 .* (r[1:end-1] .+ r[2:end])

        power_prf = Float64.(itp.(r_orig_centers) .* (R ./ maximum(r)))

        return power_prf, power_map, k2_T, k2_L
    end
    

    """
        run_solid1d_mush(omega, rho, radius, visc, shear, bulk, phi, R, ρ_core, μ_core, κ_core; ncalc=2000, n=2, m=2, visc_l=1e2, bulk_l=1e9, permea=1e-7, porosity_thresh=1e-5)

    Use 1D solid tides model with mush interface to calculate k2 Lovenumbers, and compute 1D heating profile from strain tensor.

    # Arguments
    - `omega::Float64`                  : Forcing frequency range.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `shear::Array{precc,1}`           : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    - `phi::Array{prec,1}`              : Melt fraction (porosity) profile of the planet.
    - `R::prec`                         : Planet radius.
    - `m_core::prec`                    : Core mass.
    - `ρ_core::prec`                    : Core density.
    - `μ_core::prec`                    : Core shear modulus.
    - `κ_core::prec`                    : Core bulk modulus.
    
    # Keyword Arguments
    - `ncalc::Int=2000`                 : Number of sublayers.
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `m::Int=2`                        : Harmonic of the true anomaly. m=2 corresponds to the semidiurnal tide, m=1 diurnal tide.
    - `core::String="liquid"`           : Core state, either "liquid" or "solid".
    - `visc_l::Float64=1e2`             : Liquid viscosity.
    - `bulk_l::Float64=1e9`             : Liquid bulk modulus.
    - `permea::Float64=1e-7`            : Permeability of mush layer.
    - `porosity_thresh::Float64=1e-5`   : Porosity threshold, below this value no mush.

    # Returns
    - `power_prf::Array{Float64,1}`     : Heating profile.
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_solid1d_mush( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1},
                        phi::Array{prec,1},
                        R::prec,
                        m_core::prec,
                        ρ_core::prec,
                        μ_core::prec,
                        κ_core::prec;
                        ncalc::Int=2000,
                        n::Int=2,
                        m::Int=2,
                        core::String="liquid",
                        visc_l::Float64=1e2,
                        bulk_l::Float64=1e9,
                        permea::Float64=1e-7,
                        porosity_thresh::Float64=1e-5
                        )::Tuple{Array{Float64,1},ComplexF64,ComplexF64}

        # internal structure arrays.
        # first element is the innermost layer, last element is the outermost layer
        omega = prec(omega)
        ρ  = copy(convert(Vector{prec}, rho))
        r  = copy(convert(Vector{prec}, radius))
        η  = copy(convert(Vector{prec}, visc))
        μc = copy(convert(Vector{precc}, shear))
        κs = copy(convert(Vector{prec}, bulk))
        ϕ  = copy(convert(Vector{prec}, phi))
        κd = 0.01.*κs                        # drained bulk modulus

        α  = 1.0.-(κd./κs)                    # Biot's modulus

        # allocate zero arrays with same length and precision as r
        κl = zeros(prec, length(r))
        ηl = zeros(prec, length(r))
        k  = zeros(prec, length(r))

        # implicitely the mush interface occurs at the top of the solid
        # the mush layer index is therefore
        ii = length(ϕ)

        # If the porosity = 0, throw error (because the matrix cannot be resolved, instead use 1 phase model)
        if ϕ[ii] <= prec(porosity_thresh)
            throw("No mush region identified in viscosity profile.")
        end

        # update the liquid arrays
        κl[ii] = prec(bulk_l)   # liquid bulk modulus
        ηl[ii] = prec(visc_l)   # liquid viscosity
        k[ii]  = prec(permea)   # permeability

        # set porosity to zero outside mush region (otherwise code cannot solve system)
        ϕ[1:ii-1]   .= 0.0      # zero below ii
        # limit the mush layer melt fraction to ϕ=<0.1, as the model becomes unstable otherwise
        ϕ[ii] = min(ϕ[ii], 0.3)

        ρs = ρ.*(1.0.-ϕ)        # solid density 
        ρl = ρ.*ϕ               # liquid density

        # subdivide input layers such that we have ~ncalc in total
        rr = solid1d_mush.expand_layers(r, nr=convert(Int,div(ncalc,length(η))))

        # get gravity at each layer
        g = solid1d_mush.get_g(rr, ρ, m_core)

        # create grid
        solid1d_mush.define_spherical_grid(res, n, m)

        # get y-functions
        M, y1_4 = solid1d_mush.compute_M(omega, rr, ρs, g, μc, κs, ρl, κl, κd, α, ηl, ϕ, k, n, ρ_core, μ_core, κ_core; core=core)
        #   Tidal
        tidal_solution_T = solid1d_mush.compute_y(rr, g, M, y1_4, n; load=false)
        #   Load
        tidal_solution_L = solid1d_mush.compute_y(rr, g, M, y1_4, n; load=true)

        # get k2 tidal Love Number (complex-valued)
        k2_T = tidal_solution_T[5, end, end] - 1
        k2_L = tidal_solution_L[5, end, end] - 1
        
        # Get profile power output (W m-3), converted to W/kg
        Eμ, Eκ, El = solid1d_mush.get_heating_profile(tidal_solution_T,
                               rr, ρs, g, μc, κs,
                               omega, ρl, κl, κd, 
                               α, ηl, ϕ, k, n)

        Eμ_tot, _ = Eμ   # shear       (W), (W/m3)
        Eκ_tot, _ = Eκ   # compaction  (W), (W/m3)
        El_tot, _ = El   # fluid       (W), (W/m3)

        # Renormalization factor
        power_prf = (Eμ_tot .+ Eκ_tot .+ El_tot) .* (R ./ maximum(r)) # Compute total volumetric heating (W/m3)

        return Float64.(power_prf), ComplexF64(k2_T), ComplexF64(k2_L)
    end


    """
        run_solid1d_mush_relax(omega, rho, radius, visc, shear, bulk, R, m_core, ρ_core, μ_core, κ_core; dr_min=300, dr_max=3000, n=2, m=2, core="liquid")

    Use 1D solid tides model with relaxation method to calculate k2 Lovenumbers, and compute 1D heating profile from strain tensor.
    This method includes inertia effects, but is more computationally expensive. 
    
    # Arguments
    - `omega::prec`                     : Forcing frequency range.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `visc::Array{prec,1}`             : Viscosity profile of the planet.
    - `shear::Array{precc,1}`           : Complex shear modulus profile of the planet.
    - `bulk::Array{prec,1}`             : Bulk modulus profile of the planet.
    - `phi::Array{prec,1}`              : Melt fraction (porosity) profile of the planet.
    - `R::prec`                         : Planet radius.
    - `m_core::prec`                    : Core mass.
    - `ρ_core::prec`                    : Core density.
    - `μ_core::prec`                    : Core shear modulus.
    - `κ_core::prec`                    : Core bulk modulus.

    # Keyword Arguments
    - `dr_min::Int=300`                 : Minimum layer thickness in m.
    - `dr_max::Int=3000`                : Maximum layer thickness in m.
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `m::Int=2`                        : Harmonic of the true anomaly. m=2 corresponds to the semidiurnal tide, m=1 diurnal tide.
    - `core::String="liquid"`           : Core state, either "liquid", "solid", or "inertial".
    - `visc_l::Float64=1e2`             : Liquid viscosity.
    - `bulk_l::Float64=1e9`             : Liquid bulk modulus.
    - `permea::Float64=1e-7`            : Permeability of mush layer.
    - `porosity_thresh::Float64=1e-5`   : Porosity threshold, below this value no mush.

    # Returns
    - `power_prf::Array{prec,1}`        : Heating profile.
    - `power_map::Array{prec,4}`        : Heating map (frequency, colatitude, longitude).
    - `k2_T::precc`                     : Complex Tidal k2 Lovenumber.
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.
    """
    function run_solid1d_mush_relax( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{precc,1},
                        bulk::Array{prec,1},
                        phi::Array{prec,1},
                        R::prec,
                        m_core::prec,
                        ρ_core::prec,
                        μ_core::prec,
                        κ_core::prec;
                        dr_min::Int=300,
                        dr_max::Int=3000,
                        n::Int=2,
                        m::Int=2,
                        core::String="liquid",
                        visc_l::Float64=1e2,
                        bulk_l::Float64=1e9,
                        permea::Float64=1e-7,
                        porosity_thresh::Float64=1e-5
                        )::Tuple{Array{Float64,1},ComplexF64,ComplexF64}
                        # )::Tuple{Array{Float64,1},Matrix{Float64},ComplexF64,ComplexF64}

        # internal structure arrays.
        # first element is the innermost layer, last element is the outermost layer
        omega = prec(omega)
        ρ  = copy(convert(Vector{prec}, rho))
        r  = copy(convert(Vector{prec}, radius))
        η  = copy(convert(Vector{prec}, visc))
        μc = copy(convert(Vector{precc}, shear))
        κs = copy(convert(Vector{prec}, bulk))
        ϕ  = copy(convert(Vector{prec}, phi))
        κd = 0.01.*κs                        # drained bulk modulus

        α  = 1.0.-(κd./κs)                   # Biot's modulus

        # allocate zero arrays with same length and precision as r
        κl = zeros(prec, length(r))
        ηl = zeros(prec, length(r))
        k  = zeros(prec, length(r))

        # find where the mush interface occurs (excluding the CMB layer)
        # get all indices above the threshold
        ii_all = findall(ϕ .>= porosity_thresh)
        
        # set all porosity values below the threshold to zero (otherwise code cannot solve system)
        ϕ[ϕ .< porosity_thresh] .= 0.0

        # update the liquid arrays
        κl .= prec(bulk_l)   # liquid bulk modulus
        ηl .= prec(visc_l)   # liquid viscosity
        k[ii_all] .= prec(permea)   # permeability

        # resample profiles onto new grid
        r_grid, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, g, M_tot = solid1d_mush_relax.resample_profiles(r, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, m_core, dr_min, dr_max)

        ρs = ρ.*(1.0.-ϕ)        # solid density 
        ρl = ρ.*ϕ               # liquid density

        # use cell centers
        r_centers = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])

        # define angular grid
        SphericalGrid = solid1d_mush_relax.define_spherical_grid(res, n, m)
        
        # solve y functions across grid
        y_t, y_l = solid1d_mush_relax.compute_y(r_centers, ρ, g, μc, κs, omega, ρl, κl, κd, α, ηl, ϕ, k, n, ρ_core, μ_core, κ_core, M_tot; core=core)

        # plotting.plot_relaxation_solution(y_t, r_centers, 
        #         filename="$OUT_DIR/relaxation_solution.png")

        # Love numbers
        k2_T = y_t[5, end] - 1
        k2_L = y_l[5, 1]   - 1

        # heating profile
        Eμ_tot, Eκ_tot, El_tot = solid1d_mush_relax.get_heating_profile(
            y_t, r_grid, ρ, g, μc, κs, omega, ρl, κl, κd, α, ηl, ϕ, k, n, SphericalGrid
        )

        # Eμ_map, Eκ_map, El_map = solid1d_mush_relax.get_heating_map(
        #     y_t, r_grid, ρ, g, μc, κs, omega, ρl, κl, κd, α, ηl, ϕ, k, n, SphericalGrid
        # )

        power_prf = abs.(Eμ_tot .+ Eκ_tot .+ El_tot) 
        # power_map = abs.(Eμ_map .+ Eκ_map .+ El_map) 

        # interpolate from grid back to original radius points 
        itp = linear_interpolation(r_centers, power_prf, extrapolation_bc=Line())

        # original centers
        r_orig_centers = 0.5 .* (r[1:end-1] .+ r[2:end])

        power_prf = Float64.(itp.(r_orig_centers) .* (R ./ maximum(r)))

        return power_prf, k2_T, k2_L
        # return power_prf, power_map, k2_T, k2_L
    end


    """
        run_fluid0d(omega, rho, radius, ρ_ratio; n=2, sigma_R=1e-3)

    Calculate k2 Lovenumbers in the 0D fluid.

    # Arguments
    - `omega::Float64`                  : Forcing frequency.
    - `rho::Array{prec,1}`              : Density profile of the planet.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `ρ_ratio::prec`                   : Density contrast between current (fluid) and lower (non-fluid) layer.
    
    # Keyword Arguments
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).
    - `sigma_R::Float64=1e-3`           : Rayleigh drag coefficient.

    # Returns
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_fluid0d( omega::Float64,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        ρ_ratio::prec;
                        n::Int64=2,
                        sigma_R::Float64=1e-3
                        )::Tuple{ComplexF64,ComplexF64}

        # internal structure arrays
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)

        # surface properties
        R = maximum(r)  # Segment outer radius (m)
        g = G * sum(4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) .* ρ[1:end]) / R^2  # Surface gravity (m/s^2)

        # fluid magma ocean bottom and height
        r_b = r[1]
        H_magma = r[end] - r_b
             
        # get k2 Lovenumbers
        k2_T, k2_L = fluid0d.compute_fluid_lovenumbers(omega, R, H_magma, g, ρ_ratio, n, sigma_R)

        return ComplexF64(k2_T), ComplexF64(k2_L)

    end


    """
        run_fluid1d(omega, rho, radius, gravity, ρ_ratio, S_mass, sma, R; kwargs...)

    Compute tidal heating profile and Love numbers for a 1D fluid model.

    # Arguments
    - `omega::Float64`                  : Forcing frequency
    - `rho::Vector{prec}`               : Density profile
    - `radius::Vector{prec}`            : Radial grid (core → surface)
    - `gravity::Vector{prec}`           : Gravity profile
    - `ρ_ratio::prec`                   : Density ratio of lower layer
    - `S_mass::prec`                    : Stellar mass
    - `sma::prec`                       : Semi-major axis
    - `R::prec`                         : Planet radius

    # Keyword Arguments
    - `n::Int=2`                        : Radial power (dominant term n=2)
    - `sigma_R::Float64=1e-3`           : Rayleigh drag at interface
    - `sigma_inf::Float64=1e-7`         : Drag in fluid interior
    - `sigma_R_prf::String="uniform"`   : Drag profile type
    - `H_R::Float64=1e3`                : Drag scale height
    - `efficiency::Float64=0.3`         : Interface efficiency factor

    # Returns
    - `power_prf::Vector{Float64}`      : Heating profile
    - `k2_T::ComplexF64`                : Tidal Love number
    - `k2_L::ComplexF64`                : Load Love number
    """
    function run_fluid1d(   omega::Float64,
                            rho::Vector{prec},
                            radius::Vector{prec},
                            gravity::Vector{prec},
                            ρ_ratio::prec,
                            S_mass::prec,
                            sma::prec,
                            R::prec;
                            n::Int = 2,
                            sigma_R::Float64 = 1e-3,
                            sigma_inf::Float64 = 1e-7,
                            sigma_R_prf::String = "uniform",
                            H_R::Float64 = 1e3,
                            efficiency::Float64 = 0.3
                        )::Tuple{Vector{Float64}, ComplexF64, ComplexF64}

        # internal structure arrays 
        r = convert(Vector{prec}, radius) 

        # volumes
        Vs  = (4/3) * π * (r[end]^3 - r[1]^3)
        dVs = (4/3) * π .* (r[2:end].^3 .- r[1:end-1].^3)

        # Love numbers (solid-fluid interface and fluid interior)
        k2_T, k2_L  = run_fluid0d(omega, rho, r, ρ_ratio; n=n, sigma_R=efficiency * sigma_R)
        k2_T_inf, _ = run_fluid0d(omega, rho, r, ρ_ratio; n=n, sigma_R=sigma_inf)

        # total heating (bulk)
        prefactor = (2n + 1) * R * omega / (8π * G)

        power_blk     = prefactor * -imag(k2_T)     / Vs
        power_blk_inf = prefactor * -imag(k2_T_inf) / Vs

        D_power_blk = max(power_blk - power_blk_inf, 0)

        # radial positions
        r_mid = 0.5 .* (r[1:end-1] .+ r[2:end])
        z     = abs.(r_mid .- r[1])  # distance from interface

        # profile shape function
        shape = ones(length(z))

        if sigma_R_prf == "exp"
            shape .= exp.(-z ./ H_R)

        elseif sigma_R_prf == "linear"
            shape .= max.(0, 1 .- z ./ H_R)

        elseif sigma_R_prf == "quadratic"
            shape .= max.(0, 1 .- z ./ H_R).^2

        elseif sigma_R_prf == "dynamic"
            l_mix = max.(min.(z, H_R), 1e-12)
            shape .= exp.(-z ./ l_mix)

        elseif sigma_R_prf != "uniform"
            error("Unknown sigma_R_prf: $sigma_R_prf")
        end

        # heating profile
        power_prf = power_blk_inf .+ D_power_blk .* shape

        # normalize to match bulk heating
        unorm = sum(power_prf .* dVs) / Vs
        if unorm > 0
            power_prf .*= power_blk / unorm
        end

        return Float64.(power_prf), ComplexF64(k2_T), ComplexF64(k2_L)
    end


    """
        run_interp(omega, radius, P_t, P_b; t_width=0.1, b_width=0.1)

    Interpolate dissipation and k2 Lovenumbers in a 1D region without active tides.

    # Arguments
    - `omega::Float64`                  : Forcing frequency.
    - `radius::Array{prec,1}`           : Radial positions of layers, from core to surface.
    - `R::prec`                         : Planet Radius.
    - `P_t::prec`                       : Heating at upper interface.
    - `P_b::prec`                       : Heating at lower interface.
    
    # Keyword Arguments
    - `t_width::Float64=0.1`            : Fraction of segment height as standard deviation for upper dissipation peak.
    - `b_width::Float64=0.1`            : Fraction of segment height as standard deviation for lower dissipation peak.

    # Returns
    - `power_prf::Array{Float64,1}`     : Heating profile.
    - `k2_T::ComplexF64`                : Complex Tidal k2 Lovenumber.
    - `k2_L::ComplexF64`                : Complex Load k2 Lovenumber.
    """
    function run_interp(omega::Float64,
                        radius::Vector{prec},
                        R::prec,
                        P_t::Real,
                        P_b::Real;
                        t_width::Real = 0.1,
                        b_width::Real = 0.1
                    )::Tuple{Vector{Float64}, ComplexF64, ComplexF64}

        # radial grid
        r = convert(Vector{prec}, radius)
        r_mid = 0.5 .* (r[1:end-1] .+ r[2:end])

        # shell volumes
        Vs = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3)

        # interface locations
        r_top = r[end]
        r_bot = r[1]

        # decay lengths
        l_t = t_width * (r[end] - r[1])
        l_b = b_width * (r[end] - r[1])

        # exponential decay profiles
        G_top = exp.(-abs.(r_mid .- r_top) ./ l_t)
        G_bot = exp.(-abs.(r_mid .- r_bot) ./ l_b)

        # normalize peaks to unity
        if maximum(G_top) > 0
            G_top ./= maximum(G_top)
        end
        if maximum(G_bot) > 0
            G_bot ./= maximum(G_bot)
        end

        # total power profile
        power_prf = P_t .* G_top .+ P_b .* G_bot

        # allocate Love numbers
        k2_T = zeros(precc, length(r_mid))
        k2_L = zeros(precc, length(r_mid))

        # tidal prefactor
        prefactor = 5 * R * omega / (8π * G)

        # infer complex k2 from dissipation
        @inbounds for is in eachindex(r_mid)
            k2_T[is] = -im * power_prf[is] * Vs[is] / prefactor
        end

        return Float64.(power_prf), ComplexF64(sum(k2_T)), ComplexF64(sum(k2_L))
    end


    """
        data_to_nc(nmk, is_seg, segments, knms_T, knms_L, σ_range, P_T_s_prf, P_T_prf, P_T_blk, P_T_prf_blk, P_T_s_map, P_T_map, datafile_path)
    
    Write model results to a NetCDF file.

    # Arguments
    - `nmk::Vector{Tuple{Int,Int,Int}}`         : Array of (n,m,k) tuples for each segment.
    - `is_seg::Array{Tuple{Int,Int},1}`         : Array of (il,it) tuples indicating segment indices.
    - `segments::Array{String,1}`               : Array of segment labels.
    - `knms_total::Array{ComplexF64}`           : Total k2 Lovenumbers for each (n,m,k).
    - `knms_T::Matrix{ComplexF64}`              : Tidal k2 Lovenumbers for each (n,m,k) and segment.
    - `knms_L::Matrix{ComplexF64}`              : Load k2 Lovenumbers for each (n,m,k) and segment.
    - `σ_range::Array{Float64,1}`               : Array of forcing frequencies.
    - `P_T_s_prf::Array{Float64,2}`             : Tidal heating profiles for each (n,m,k) and depth.
    - `P_T_prf::Array{Float64,1}`               : Tidal heating profile for each depth.
    - `P_T_blk::Float64`                        : Tidal heating in the black hole.
    - `P_T_prf_blk::Float64`                    : Tidal heating profile in the black hole.
    - `P_T_s_map::Array{Float64,3}`             : Tidal heating map for each (n,m,k) and spatial location.
    - `P_T_map::Array{Float64,2}`               : Tidal heating map for each spatial location.
    - `datafile_path::String`                   : Path to the output NetCDF file.
    """
    function data_to_nc(nmk::Vector{Tuple{Int,Int,Int}}, 
                        is_seg::Array{Tuple{Int,Int},1}, 
                        segments::Array{String,1},  
                        knms_total::Array{ComplexF64},
                        knms_T::Matrix{ComplexF64}, 
                        knms_L::Matrix{ComplexF64}, 
                        σ_range::Array{Float64,1}, 
                        P_T_s_prf::Array{Float64,2}, 
                        P_T_prf::Array{Float64,1}, 
                        P_T_blk::Float64, 
                        P_T_prf_blk::Float64, 
                        P_T_s_map::Array{Float64,3}, 
                        P_T_map::Array{Float64,2}, 
                        datafile_path::String
                       )
        
        # helper function to convert Complex arrays to Float64 arrays with an extra dimension
        function pack_complex(A::AbstractArray{Complex{T}}) where T
            # Creates an array with an extra trailing dimension of size 2
            return cat(real.(A), imag.(A), dims=ndims(A)+1)
        end

        @info "Writing results to $datafile_path"

        n_vals = [t[1] for t in nmk]
        m_vals = [t[2] for t in nmk]
        k_vals = [t[3] for t in nmk]

        il_vals = [t[1] for t in is_seg]
        it_vals = [t[2] for t in is_seg]

        NCDataset(datafile_path, "c") do ds
            
            # DEFINE DIMENSIONS
            defDim(ds, "nmk",      length(nmk))
            defDim(ds, "segments", length(segments))
            defDim(ds, "z",        length(P_T_prf))
            defDim(ds, "lon",      length(0:res:360-0.001))
            defDim(ds, "lat",      length(0:res:180))
            defDim(ds, "complex",  2)

            # PART 1
            defVar(ds, "n", n_vals, ("nmk",))
            defVar(ds, "m", m_vals, ("nmk",))
            defVar(ds, "k", k_vals, ("nmk",))

            defVar(ds, "il", il_vals, ("segments",))
            defVar(ds, "it", it_vals, ("segments",))
            
            defVar(ds, "segment_lbl", segments, ("segments",))
            
            # data grids, packed into ("nmk", "segments", "complex")
            defVar(ds, "knms_T", pack_complex(knms_T), ("nmk", "segments", "complex"))
            defVar(ds, "knms_L", pack_complex(knms_L), ("nmk", "segments", "complex"))
            
            defVar(ds, "sigma_range", σ_range, ("nmk",))
            defVar(ds, "knms_total", pack_complex(knms_total), ("nmk", "complex"))

            # PART 2
            defVar(ds, "P_T_s_prf", P_T_s_prf, ("nmk", "z"))

            # PART 3
            defVar(ds, "P_T_prf",     P_T_prf,     ("z",))
            defVar(ds, "P_T_blk",     P_T_blk,     ()) 
            defVar(ds, "P_T_prf_blk", P_T_prf_blk, ()) 

            # PART 4
            defVar(ds, "P_T_s_map", P_T_s_map, ("nmk", "lat", "lon"))

            # PART 5
            defVar(ds, "lon", collect(0:res:360-0.001), ("lon",))
            defVar(ds, "lat", collect(0:res:180), ("lat",))
            defVar(ds, "P_T_map", P_T_map, ("lat", "lon"))
            
            ds.attrib["title"] = "Obliqua Model Run Data"
        end
    end
        

end