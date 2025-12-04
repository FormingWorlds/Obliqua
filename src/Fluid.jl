

module Fluid

    using Interpolations
    using LinearAlgebra
    include("Hansen.jl")
    using .Hansen

    export calc_fluid_tides

    # Precision types
    prec = BigFloat
    precc = Complex{BigFloat}


    # Constants
    AU = 1.495978707e11  # m
    G  = 6.6743e-11  # m^3 kg^-1 s^-2

    # Degree Love number
    n = 2
    m = 2
    k_min = -30
    k_max = 40
    k_range = collect(k_min:k_max)

    # Get fluid tidal heating
    function calc_fluid_tides( omega::prec,
                        axial::prec,
                        ecc::prec,
                        sma::prec,
                        S_mass::prec,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1};
                        N_sigma::Int=301,
                        visc_l::Float64=2e2,
                        visc_s::Float64=5e21
                        )::Tuple{Array{Float64,1},Float64,Float64}

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
            ρ_l_mean = mean_rho(ρ_l, r_l, r_b)
        end
        
        if length(ρ_s) == 0
            ρ_s_mean = 5000.0  # arbitrary high density
        elseif length(ρ_s) == 1
            ρ_s_mean = ρ_s[1]
        else
            ρ_s_mean = mean_rho(ρ_s, r_s, r_c)
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
        k_range2, X_hansen = get_hansen(ecc, n, m, k_min, k_max)

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

        # get friction at interface, needs to be changed
        # employ correlation with mixing length
        σ_R = 10 ^ (-3)

        # fluid Love Numbers
        k22_fluid_high_friction, k22_total = compute_fluid_lovenumbers(
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
        power_prf[mask_l] .= heat_profile(P_tidal_total, ρ_l, r_l, r_b)
        
        # Get k2 lovenumber at forcing frequency used by Love.jl
        k2_total = interp_full(2*axial - omega)
        
        return convert(Vector{Float64}, power_prf), convert(Float64, P_tidal_total), convert(Float64, k2_total)

    end

    function heat_profile( power::prec, 
                            ρ::Array{prec,1},
                            r::Array{prec,1},
                            b::prec
                            )::Array{prec,1}
        # concatenate bottom radius and region radii
        r = vcat(b, r)

        # calculate mean density of region
        V = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) # Volume of each layer (m^3)
        
        # calculate power density
        power_ρ = power / sum(V) # W/m3

        # calculate power profile
        power_prf = power_ρ ./ ρ # W/kg

        return power_prf

    end

    # Identiy liquid region
    function mean_rho( ρ::Array{prec,1},
                        r::Array{prec,1},
                        b::prec
                        )::prec
        """Calculate mean density of a sphere given density profile and radius profile.

        Args:
            ρ (Array{prec,1}): Density profile (kg/m^3).
            r (Array{prec,1}): Radius profile (m).

        Returns:
            prec: Mean density (kg/m^3).
        """
        # concatenate bottom radius and region radii
        r = vcat(b, r)

        # calculate mean density of region
        V = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) # Volume of each layer (m^3)
        M = sum(V .* ρ[1:end])                       # Total mass (kg)
        mean_density = M / sum(V)                    # Mean density (kg/m^3)

        return mean_density

    end

    function compute_fluid_lovenumbers(
            n::Int,
            sigma_range::Array{Float64,1},
            k_T_22_homo::Array{precc,1},
            k_L_22_homo::Array{precc,1},
            rho_ratio::prec,
            P_grav_acc::prec,
            H_magma::prec,
            sigmaR::Float64,
            P_radius::prec
        )::Tuple{Array{precc,1},Array{precc,1}}
        """Calculate k2 lovenumber in liquid.

        Args

        """
        N_sigma = length(sigma_range)

        mu_n  = n * (n + 1)
        ksi_n = 3 / (2*n + 1) * rho_ratio
        sigP_n = sqrt(mu_n * P_grav_acc * H_magma / P_radius^2)

        k22_fluid = zeros(precc, N_sigma)
        k22_total = zeros(precc, N_sigma)

        for i in 1:N_sigma
            sigma = sigma_range[i]
            sigT = sigma - im * sigmaR

            # Eq. (7)
            k22_fluid[i] = -ksi_n * sigP_n^2 / (sigma * sigT - sigP_n^2)

            # Eq. (28), ignoring crust
            k22_total[i] = k_T_22_homo[i] + (1 + k_L_22_homo[i]) * k22_fluid[i]
        end

        return k22_fluid, k22_total
    end

end