

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

    """
        plot_imag_k2_sigma(filename::String, σ_R_vals::Vector{Float64};
                        height_frac::Float64=0.12,
                        N_sigma::Int=301,
                        visc_l::Float64=2e2,
                        visc_s::Float64=5e21,
                        xmin::Float64=1e-6,
                        xmax::Float64=5.0,
                        ymin::Float64=1e-5,
                        ymax::Float64=2.0,
                        output_path::String="imag_k2_sigma.png")

    Plot `Imag(k₂)` versus forcing frequency `σ` for a range of Rayleigh drag parameters `σ_R`.

    # Arguments
    - `filename::String` : Path to the JSON file containing interior structure data.  
    - `σ_R_vals::Vector{Float64}` : Array of Rayleigh drag parameters to plot.  

    # Keyword Arguments
    - `height_frac::Float64=0.12` : Height fraction (currently informational).  
    - `N_sigma::Int=301` : Number of frequency points.  
    - `visc_l::Float64=2e2` : Liquid viscosity threshold (Pa·s).  
    - `visc_s::Float64=5e21` : Solid viscosity threshold (Pa·s).  
    - `xmin`, `xmax`, `ymin`, `ymax::Float64` : Axis limits for log-log plot.  
    - `output_path::String="imag_k2_sigma.png"` : Filepath to save figure.

    # Returns
    - `plt` : The `Plots.Plot` object.
    """
    function plot_imag_k2_sigma(filename::String, σ_R_vals::Vector{Float64};
                                height_frac::Float64=0.12,
                                N_sigma::Int=301,
                                visc_l::Float64=2e2,
                                visc_s::Float64=5e21,
                                xmin::Float64=1e-6,
                                xmax::Float64=5.0,
                                ymin::Float64=1e-5,
                                ymax::Float64=2.0,
                                output_path::String="imag_k2_sigma.png")

        # Helper for clean log ticks
        function logticks_labeled(min_exp, max_exp)
            major = [10.0^e for e in min_exp:max_exp]
            minor = vcat([m * 10.0^e for e in min_exp:max_exp, m in 2:9]...)
            major_labels = [@sprintf("%.0e", x) for x in major]
            return (major, major_labels), minor
        end

        # Compute log tick positions
        x_exp_min, x_exp_max = floor(Int, log10(xmin)), ceil(Int, log10(xmax))
        y_exp_min, y_exp_max = floor(Int, log10(ymin)), ceil(Int, log10(ymax))
        (x_major, x_major_labels), x_minor = logticks_labeled(x_exp_min, x_exp_max)
        (y_major, y_major_labels), y_minor = logticks_labeled(y_exp_min, y_exp_max)

        # Initialize plot
        plt = plot(
            title="Imag(k₂) vs σ for different σ_R at $height_frac height fraction",
            xlabel="σ",
            ylabel="Imag(k₂)",
            xscale = :log10,
            yscale = :log10,
            xlims = (xmin, xmax),
            ylims = (ymin, ymax),
            xticks = (x_major, x_major_labels),
            yticks = (y_major, y_major_labels),
            minorgrid = true,
            minorticks = :auto,
        )

        # Load interior once
        omega, axial, ecc, sma, S_mass, rho, radius, visc =
            load.load_interior_liquid(filename, false)

        # Loop over σ_R values
        for σ_R in σ_R_vals
            σ_range, imag_k2, _height_frac = Fluid.calc_fluid_tides_full(
                omega, axial, ecc, sma, S_mass, rho, radius, visc;
                N_sigma=N_sigma,
                visc_l=visc_l,
                visc_s=visc_s,
                sigma_R=σ_R,
            )

            mask = (σ_range .< 0) .& (imag_k2 .> 0)

            plot!(
                plt,
                -σ_range[mask],
                imag_k2[mask],
                label = "σ_R = $σ_R"
            )
        end

        savefig(plt, output_path)
        return plt
    end


end