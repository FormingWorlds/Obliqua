module solid0d

    include("constants.jl")
    using .constants
    
    export compute_solid_lovenumbers
    
    
    """
        compute_solid_lovenumbers(μc, mass_tot, R, n)    
    
    Calculate k2 lovenumbers in solid.

    # Arguments
    - `μc::precc`                       : Complex shear modulus.
    - `mass_tot::prec`                  : Total mass of the planet.
    - `R::prec`                         : Outer radius of the planet.
    - `n::Int=2`                        : Power of the radial factor (goes with (r/a)^{n}, since r<<a only n=2 contributes significantly).

    # Returns
    - `k2_T::precc`                     : Complex Tidal k2 Lovenumber.
    - `k2_L::precc`                     : Complex Load k2 Lovenumber.
    """
    function compute_solid_lovenumbers(
            μc::precc,
            mass_tot::prec,
            R::prec,
            n::Int
        )::Tuple{precc,precc}
             
        # calculate parameters
        An = 4. * (2. * n^2 + 4. * n + 3.) / (3. * n * G * mass_tot^2) * π * R^4
        μc_n = An * μc
        factor = 1. / (1. + μc_n)

        if n == 1
            k2_T = 0.
            # h2_T = 0.
        else
            k2_T = factor * 3. / (2. * (n - 1.))
            # h2_T = factor * (2. * n + 1.) / (2. * (n - 1.))
        end

        k2_L = factor * -1.0
        # h2_L = factor * (-(2.0 * n + 1.0) / 3.0)

        return k2_T, k2_L

    end


    """
        mean_cmu(μc, r)

    Calculate the Hill-averaged complex shear modulus in a spherical segment (Hill+1952).
    
    # Arguments
    - `μc::Array{precc,1}`              : Complex shear modulus profile of the planet.
    - `r::Array{prec,1}`                : Radial positions of layers, from bottom to top of segment.

    # Returns
    - `mean_μc::precc`                  : Mean complex shear modulus in segment.
    
    # Notes
    - (DOI 10.1088/0370-1298/65/5/307)
    """
    function mean_cmu(
            μc::Vector{precc},
            r::Vector{prec}
        )::precc

        # Volume of each spherical shell
        V = (4/3) * π .* (r[2:end].^3 .- r[1:end-1].^3)

        # Volume fractions
        f = V ./ sum(V)

        # Voigt average (uniform strain)
        μV = sum(f .* μc)

        # Reuss average (uniform stress)
        μR = inv(sum(f ./ μc))

        # Hill average
        mean_μc = 0.5 * (μV + μR)

        return mean_μc
    end

end