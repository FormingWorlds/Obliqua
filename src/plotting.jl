

module plotting

    using Plots
    using Printf


    """
        plot_imagk2_spectrum(σ_range, imag_k2; outpath)

    Create and save a plot of the k2 Lovenumber as function of forcing frequency.

    # Arguments
    - `σ_range::Array{prec,1}`          : Forcing frequency range.
    - `imag_k2::Array{precc,1}`         : k2 Lovenumbers.
    
    # Keyword Arguments
    - `outpathn::Union{Nothing,String}=nothing` : If provided, the figure is saved to this path using `savefig`.

    # Returns
    - `plt`                             : The `Plots.Plot` object.
    """
    function plot_imagk2_spectrum(σ_range, imag_k2; outpath::Union{Nothing,String}=nothing)

        # helper for scientific-notation log ticks
        function logticks_labeled(min_exp, max_exp)
            major = [10.0^e for e in min_exp:max_exp]
            minor = vcat([m * 10.0^e for e in min_exp:max_exp, m in 2:9]...)
            major_labels = [@sprintf("%.0e", x) for x in major]
            return (major, major_labels), minor
        end

        # keep only positive forcing frequencies
        mask = σ_range .> 0

        σ_range = σ_range[mask]
        imag_k2 = imag_k2[mask]

        # axis limits with padding
        xmin_raw, xmax_raw = minimum(σ_range), maximum(σ_range)
        ymin_raw, ymax_raw = minimum(abs.(imag_k2[imag_k2 .!= 0])), maximum(abs.(imag_k2))

        # multiplicative padding factors 
        pad_low  = 0.8
        pad_high = 1.2

        xmin = xmin_raw * pad_low
        xmax = xmax_raw * pad_high
        ymin = max(ymin_raw * pad_low, 1e-7)
        ymax = ymax_raw * pad_high

        # exponent ranges
        x_exp_min = floor(Int, log10(xmin))
        x_exp_max = ceil(Int,  log10(xmax))
        y_exp_min = floor(Int, log10(ymin))
        y_exp_max = ceil(Int,  log10(ymax))

        (x_major, x_major_labels), x_minor = logticks_labeled(x_exp_min, x_exp_max)
        (y_major, y_major_labels), y_minor = logticks_labeled(y_exp_min, y_exp_max)

        # initialize plot
        plt = plot(
            title = "Imag(k₂) vs σ",
            xlabel = "σ [Hz]",
            ylabel = "Imag(k₂)",
            xscale = :log10,
            yscale = :log10,
            xlims = (xmin, xmax),
            ylims = (ymin, ymax),
            xticks = (x_major, x_major_labels),
            yticks = (y_major, y_major_labels),
            minorgrid = true,
            minorticks = :auto,
        )

        # keep only positive Imag(k2)
        mask = imag_k2 .> 0

        plot!(
            plt,
            σ_range[mask],
            imag_k2[mask],
            label = "Imag(k₂)",
            lw = 2,
            marker = :circle,
            markersize = 3,
        )

        # optional saving
        if outpath !== nothing
            savefig(plt, outpath)
        end

        return plt
    end


    """
        save_heat_profile(radius, power_prf; filename="heat_profile.png")

    Create and save a plot of the heat dissipation profile as a function of radius.

    # Returns
    - `plt` : The `Plots.Plot` object.
    """
    function save_heat_profile(radius::AbstractVector,
                            power_prf::AbstractVector;
                            filename::String = "heat_profile.png")

        # Keep only positive values (required for log scale)
        mask = power_prf .> 0
        r = radius[mask]
        p = power_prf[mask]

        @assert !isempty(p) "No positive heating values to plot."

        # raw limits
        xmin_raw, xmax_raw = minimum(r), maximum(r)
        ymin_raw, ymax_raw = minimum(p), maximum(p)

        # multiplicative padding
        pad_low  = 0.8
        pad_high = 1.2

        xmin = xmin_raw * pad_low
        xmax = xmax_raw * pad_high
        ymin = max(ymin_raw * pad_low, 1e-11)  # safe floor for log scale
        ymax = ymax_raw * pad_high

        plt = plot(
            r,
            p,
            xlabel = "Radius (m)",
            ylabel = "Power per shell (W/kg)",
            title  = "Tidal Heating Profile",
            lw     = 2,
            legend = false,
            yscale = :log10,
            xlims  = (xmin, xmax),
            ylims  = (ymin, ymax),
            minorgrid = true,
        )

        savefig(plt, filename)
        return plt
    end
   
end