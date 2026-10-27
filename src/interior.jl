


module interior

    export get_permeability, limit_porosity, get_drained_bulk


    """
        get_permeability(phi, cfg; ϕ0_width=0.02, ϕcrit_width=0.05)

    Calculate permeability profile based on melt fraction (porosity) and configuration parameters.
    Uses continuous tanh smoothing between the Blake-Kozeny-Carman, Rumpf-Gupte, and Stokes regimes.

    # Arguments
    - `phi::Array{prec,1}`              : Melt fraction (porosity) profile of the planet.
    - `cfg::Dict`                       : Configuration dictionary containing interior energetics parameters.

    # Keyword Arguments
    - `ϕ0_width::Float64=0.02`          : Width of the smoothing transition around the Blake-Kozeny-Carman to Rumpf-Gupte regime
    - `ϕcrit_width::Float64=0.05`       : Width of the smoothing transition around the Rumpf-Gupte to Stokes regime

    # Returns
    - `perm::Array{prec,1}`             : Permeability profile of the planet.
    """
    function get_permeability(phi, cfg; ϕ0_width=0.02, ϕcrit_width=0.05)

        # extract grain size from config
        Rgrain = cfg["interior_energetics"]["grain_size"] # grain size in m

        # get native precision type from input arrays
        T = eltype(phi)

        # equality points and smoothing widths
        ϕ0_thresh    = 0.0769452
        ϕcrit_thresh = 0.771462
                
        # calculate permeability element-wise using tanh smoothing
        perm = ((p) -> begin
            # guard against division by zero at limits
            por = max(p, 1e-20)
            one_m_por = max(1.0 - p, 1e-20)

            # Blake-Kozeny-Carman (low porosity)
            F_bkc = (1.0 / 1000.0) * Rgrain^2 * por^2 / one_m_por^2
            
            # Rumpf-Gupte (intermediate porosity)
            F_rg = (5.0 / 7.0) * Rgrain^2 * por^4.5
            
            # Stokes settling (high porosity)
            F_stokes = (2.0 / 9.0) * Rgrain^2

            # smooth transition weights
            w_rg = 0.5 * (1.0 + tanh((por - ϕ0_thresh) / ϕ0_width))
            w_stokes = 0.5 * (1.0 + tanh((por - ϕcrit_thresh) / ϕcrit_width))

            # blend the three regimes
            F = (1.0 - w_rg) * F_bkc + (w_rg - w_stokes) * F_rg + w_stokes * F_stokes
            
            # ensure result is strictly non-negative
            max(F, 0.0)
        end).(phi)

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
        T = eltype(perm)

        # apply thresholding to porosity array
        phi[phi .< porosity_thresh] .= T(0.0)

        # calculate interior mask
        perm[phi .< porosity_thresh] .= T(0.0)

        return perm, phi
    end


    """
        get_drained_bulk(bulk, phi, cfg)

    Calculate the drained (compaction) bulk modulus for each layer from its porosity and grain bulk modulus K_s,
    following the power law of Hay et al. (2025, Eq. 48):

        K_d = K_s                     for phi < phi_0
        K_d = K_s (phi_0 / phi)^b     for phi_0 <= phi < phi_crit
        K_d = 0                       for phi >= phi_crit   (no coherent skeleton; Biot alpha = 1 - K_d/K_s = 1)

    phi_0 = 0.0769 and b = `orbit.obliqua.solid.dbulk_power`. phi_crit = `orbit.obliqua.solid.dbulk_phi_crit` (default
    0.7715) is the melt fraction above which the solid skeleton no longer resists isotropic stress. The values 0.0769 and
    0.7715 are the crossover melt fractions of the Blake-Kozeny-Carman, Rumpf-Gupte and Stokes permeability laws in
    `get_permeability`. Hay et al. (2025) place phi_crit at the disaggregation limit (~0.3), where the Kervazo et al. (2021)
    moduli also lose their shear rigidity. With K_d = 0 and a liquid-like shear modulus the radial stiffness K_d + 4 mu/3
    becomes very small, and the porous solver then loses passivity and energy closure for 0.3 < phi < ~0.7; the default
    therefore keeps a compaction stiffness up to the Stokes crossover.

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
        T  = eltype(bulk)
        cT = complex(T)

        # melt fraction below which the skeleton behaves like a coherent solid
        ϕ0    = 0.0769
        # melt fraction above which there is no coherent solid skeleton to resist isotropic stresses
        ϕcrit = Float64(get(cfg["orbit"]["obliqua"]["solid"], "dbulk_phi_crit", 0.7715))
        if !(ϕ0 < ϕcrit <= 1.0)
            throw(ArgumentError("dbulk_phi_crit = $ϕcrit must lie in ($ϕ0, 1]"))
        end
        
        # calculate drained bulk modulus
        bulkd = @. ((p, k) -> p < ϕ0 ? cT(k) : (p < ϕcrit ? cT(k * (ϕ0 / p)^b) : cT(0.0)))(phi, bulk)
        
        return T.(bulkd)
    end

end