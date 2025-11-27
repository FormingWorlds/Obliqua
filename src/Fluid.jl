

module Fluid

    include("Hansen.jl")
    using .Hansen

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
    k_range = np.arange(k_min, k_max + 1)

    # Get fluid tidal heating
    calc_fluid_tides( omega::prec,
                        ecc::prec,
                        sma::prec,
                        rho::Array{prec,1},
                        radius::Array{prec,1},
                        visc::Array{prec,1},
                        shear::Array{prec,1},
                        bulk::Array{prec,1};
                        ncalc::Int=2000,
                        material::String="andrade",
                        N_sigma::Int=301,
                        visc_l::prec=1e2
                        )::Tuple{Array{Float64,1},Float64,Float64}

        # Internal structure arrays.
        # First element is the innermost layer, last element is the outermost layer
        ρ = convert(Vector{prec}, rho)
        r = convert(Vector{prec}, radius)
        η = convert(Vector{prec}, visc)
        μ = convert(Vector{precc},shear)
        κ = convert(Vector{prec}, bulk)

        
        R = findmax(r)  # Planet radius (m)
        g = G * sum(4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) .* ρ[1:end-1]) / R^2  # Surface gravity (m/s^2)

        # g = Love.get_g(r, ρ)
        # g_s = g[end] # Surface gravity (m/s^2)


        # The follwoing code should be 
        # moved to a separate function

        visc_l = 1e2
        visc_s = 5e21

        mask_l = η <= visc_l
        mask_s = η > visc_s

        r_l = r[mask_l]
        r_s = r[mask_s]
        r_m = r[.!mask_l .& .!mask_s] # mush region

        ρ_l = ρ[mask_l]
        ρ_s = ρ[mask_s]
        ρ_m = ρ[.!mask_l .& .!mask_s] # mush region

        #######
        #######

        ρ_l_mean = mean_rho(ρ_l, r_l)
        ρ_s_mean = mean_rho(ρ_s, r_s)
        
        # magma ocean height
        H_magma = r_l[1] - r_s[end] # might be other way around ??

        # density ratio
        ρ_ratio = ρ_l_mean / ρ_s_mean

        
        # get hansen coefficients
        k_range2, X_hansen = get_hansen(ecc, n, m, k_min, k_max)








        r_l, r_s, rho_l, rho_s = find_liquid_region_and_densities(density, radius, visc, visc_thresh)
        H_magma = P_radius - np.max(r_s)
        
        rho_ratio = rho_l / rho_s

        P_n_orb = omega # orbital freq (s^-1)
        Omega   = axial # axial freq (s^-1)


    # Identiy liquid region

    function mean_rho( ρ::Array{prec,1},
                        r::Array{prec,1}
                        )::prec
        """Calculate mean density of a sphere given density profile and radius profile.

        Args:
            ρ (Array{prec,1}): Density profile (kg/m^3).
            r (Array{prec,1}): Radius profile (m).

        Returns:
            prec: Mean density (kg/m^3).
        """
        V = 4/3 * π * (r[2:end].^3 .- r[1:end-1].^3)  # Volume of each layer (m^3)
        M = sum(4/3 * π * (r[2:end].^3 .- r[1:end-1].^3) .* ρ[1:end-1])  # Total mass (kg)
        mean_density = M / V[end]  # Mean density (kg/m^3)
        return mean_density
    end

    function get_hansen( ecc::prec
                        )::Tuple{Array{Int,1}, Array{prec,2}}
        """Calculate Hansen coefficients for given eccentricity.

        Args:
            ecc (prec): Eccentricity.

        Returns:
            Tuple{Array{Int,1}, Array{prec,2}}: k_range and X_hansen matrix.
        """
        k_range = collect(k_min:k_max)
        N_k = length(k_range)
        X_hansen = zeros(prec, N_k, N_k)

        for (i, k) in enumerate(k_range)
            for (j, q) in enumerate(k_range)
                X_hansen[i, j] = hansen_fft(n, k, q, ecc)
            end
        end

        return k_range, X_hansen
    end

    function hansen_fft( n::Int,
                                k::Int,
                                q::Int,
                                ecc::prec
                                )::prec
        """Calculate Hansen coefficient X_n^{k,q}(e).

        Args:
            n (Int): Degree.
            k (Int): Order.
            q (Int): Index.
            ecc (prec): Eccentricity.

        Returns:
            prec: Hansen coefficient.
        """
        # Placeholder implementation
        # Actual implementation would involve series expansion or numerical integration
        return ecc^(abs(k) + abs(q))  # Simplified example
    end

    # Calculate eccentric anomaly

    # Hansen coefficients

    # Calculate Love number

end