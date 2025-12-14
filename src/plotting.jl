

module plotting

    using Plots

    """
        save_heat_profile(radius, power_prf; filename="heat_profile.png")

    Create and save a plot of the heat dissipation profile as a function of radius.

    # Returns
    - `plt` : The `Plots.Plot` object.
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
        return plt
    end

end