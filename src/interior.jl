


module interior

    export get_permeability, limit_porosity, get_drained_bulk


    """
        get_permeability(phi, cfg)

    Calculate permeability profile based on melt fraction (porosity) and configuration parameters.

    # Arguments
    - `phi::Array{prec,1}`              : Melt fraction (porosity) profile of the planet.
    - `cfg::Dict`                       : Configuration dictionary containing interior energetics parameters.

    # Returns
    - `perm::Array{prec,1}`             : Permeability profile of the planet.
    """
    function get_permeability(phi, cfg)

        # extract grain size from config
        Rgrain = cfg["interior_energetics"]["grain_size"] # grain size in m, used to calculate surface area of spheres in the mush layer

        # get native precision type from input arrays
        T  = typeof(phi[1])

        # define mush layer properties --> should move to config file in the future
        ϕ0    = 0.0769 # melt fraction below which the skeleton behaves like a coherent solid
        ϕcrit = 0.7715 # critical porosity above which there is no coherent solid skeleton to resist isotropic stresses
                
        # calculate permeability element-wise, then zero-out elements outside mush zones
        perm = @. ((p) ->
            p > ϕcrit ? (2/9) * Rgrain^2 :
            (p >= ϕ0 ? (5/7) * Rgrain^2 * p^4.5 :
                    (1/1000) * Rgrain^2 * p^2 / (1 - p)^2)
        )(phi)

        return T.(perm)
    end


    """
        limit_porosity(perm, phi, cfg)

    Limit porosity and permeability profiles based on a threshold value from the configuration.

    # Arguments
    - `perm::Array{prec,1}`             : Permeability profile of the planet.
    - `phi::Array{prec,1}`              : Melt fraction (porosity) profile of the planet.
    - `cfg::Dict`                       : Configuration dictionary containing interior energetics parameters.

    # Returns
    - `perm::Array{prec,1}`             : Updated permeability profile of the planet.
    - `phi::Array{prec,1}`              : Updated melt fraction (porosity) profile of the planet.
    """
    function limit_porosity(perm, phi, cfg)

        # extract porosity threshold from config
        porosity_thresh = cfg["orbit"]["obliqua"]["solid"]["porosity_thresh"]

        # get native precision type from input arrays
        T = typeof(perm[1])

        # apply thresholding to porosity array
        phi[phi .< porosity_thresh] .= T(0.0)

        # calculate interior mask
        perm[phi .< porosity_thresh] .= T(0.0)

        return perm, phi
    end


    """
        get_drained_bulk(bulk, phi, cfg)

    Calculate the drained bulk modulus for each layer based on its porosity and bulk modulus.

    Note: The expression given here is a simple power law relationship, a more physically motivated model may be implemented in the future.

    # Arguments
    - `bulk::Array{prec,1}`            : Bulk modulus profile of the planet.
    - `phi::Array{prec,1}`             : Melt fraction (porosity) profile of the planet.
    - `cfg::Dict`                      : Configuration dictionary containing interior energetics parameters.

    # Returns
    - `bulkd::Array{prec,1}`           : Drained bulk modulus profile of the planet.
    """
    function get_drained_bulk(bulk, phi, cfg)

        # extract porosity threshold from config
        b = cfg["orbit"]["obliqua"]["solid"]["dbulk_power"] # power law exponent for drained bulk modulus

        # get native precision type from input arrays
        T  = typeof(bulk[1])
        cT = complex(T)

        # define mush layer properties --> should move to config file in the future
        ϕ0    = 0.0769 # melt fraction below which the skeleton behaves like a coherent solid
        ϕcrit = 0.7715 # critical porosity above which there is no coherent solid skeleton to resist isotropic stresses
        
        # calculate drained bulk modulus
        bulkd = @. ((p, k) -> p < ϕ0 ? cT(k) : (p < ϕcrit ? cT(k * (ϕ0 / p)^b) : cT(0.0)))(phi, bulk)
        
        return T.(bulkd)
    end

end