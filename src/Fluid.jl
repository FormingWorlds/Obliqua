

module Fluid

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


        visc_s = 5e21

        mask_l = η <= visc_l
        mask_s = η > visc_s

        r_l = r[mask_l]
        r_s = r[mask_s]

        ρ_l = ρ[mask_l]
        ρ_s = ρ[mask_s]

        







        r_l, r_s, rho_l, rho_s = find_liquid_region_and_densities(density, radius, visc, visc_thresh)
        H_magma = P_radius - np.max(r_s)
        
        rho_ratio = rho_l / rho_s

        P_n_orb = omega # orbital freq (s^-1)
        Omega   = axial # axial freq (s^-1)


    # Identiy liquid region

    # Calculate eccentric anomaly

    # Hansen coefficients

    # Calculate Love number

end