

module plotting

    using Plots
    using Printf


    """
        plot_imagk2_spectrum(σ_range, imag_k2; outpath)

    Create and save a plot of the Imk2 Lovenumber as function of forcing frequency.

    # Arguments
    - `σ_range::Array{prec,1}`          : Forcing frequency range.
    - `imag_k2::Array{precc,1}`         : k2 Lovenumbers.
    
    # Keyword Arguments
    - `outpath::Union{Nothing,String}=nothing` : If provided, the figure is saved to this path using `savefig`.

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
        plot_imagk2_spectra(σ_range, imag_k2, segments; outpath)

    Create and save a plot of the segment wise Imk2 Lovenumber as function of forcing frequency.

    # Arguments
    - `σ_range::Array{prec,1}`          : Forcing frequency range.
    - `imag_k2::Array{precc,1}`         : k2 Lovenumbers.
    - `segments::Vector{String}`        : Names of segments corresponding to each column of imag_k2.
    
    # Keyword Arguments
    - `outpath::Union{Nothing,String}=nothing` : If provided, the figure is saved to this path using `savefig`.

    # Returns
    - `plt`                             : The `Plots.Plot` object.
    """
    function plot_imagk2_spectra(σ_range, imag_k2, segments; outpath::Union{Nothing,String}=nothing)

        # helper for scientific-notation log ticks
        function logticks_labeled(min_exp, max_exp)
            major = 10.0 .^ (min_exp:max_exp)
            minor = vec([m * 10.0^e for e in min_exp:max_exp for m in 2:9])
            labels = [@sprintf("%.0e", x) for x in major]
            return (major, labels), minor
        end

        # mask for log-scale safety 
        kabs = abs.(imag_k2)

        σ_mask = isfinite.(σ_range) .& (σ_range .> 0)
        k_mask = isfinite.(kabs) .& (kabs .> 0)

        # require σ valid AND at least one segment valid
        combined_mask = σ_mask .& vec(any(k_mask, dims=2))

        if !any(combined_mask)
            error("No valid data for log-scale plotting.")
        end

        σ = σ_range[combined_mask]
        k = kabs[combined_mask, :]

        # --- limits ---
        xmin_raw, xmax_raw = extrema(σ)
        ymin_raw, ymax_raw = extrema(k[k .> 0])

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

        (x_major, x_labels), _ = logticks_labeled(x_exp_min, x_exp_max)
        (y_major, y_labels), _ = logticks_labeled(y_exp_min, y_exp_max)

        # --- initialize plot ---
        plt = plot(
            title = "Imag(k₂) Spectrum",
            xlabel = "σ [Hz]",
            ylabel = "Imag(k₂)",
            xscale = :log10,
            yscale = :log10,
            xlims = (xmin, xmax),
            ylims = (ymin, ymax),
            xticks = (x_major, x_labels),
            yticks = (y_major, y_labels),
            minorgrid = true,
            minorticks = :auto,
        )

        colors = distinguishable_colors(length(segments))

        # --- plot segments ---
        for (iseg, seg) in pairs(segments)
            plot!(
                plt,
                σ,
                k[:, iseg],
                label = seg,
                lw = 2,
                color = colors[iseg],
                marker = :circle,
                markersize = 3
            )
        end

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

        # Keep only values greater than zero (required for log scale)
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
        ymax = ymax_raw * pad_high
        ymin = max(ymax / 1e8, ymin_raw * pad_low)

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
   
    
    """
        plot_segment_heating(H, k_range, r; mask_floor=0.0, filename="tidal_heating_map.png", title_str="Hansen norm heating")

    Create and save a heatmap of the tidal heating as function of radius and forcing frequency.

    # Arguments
    - `H::AbstractMatrix`               : Tidal heating values for each radius and forcing frequency.
    - `k_range::AbstractVector`         : Forcing frequency range corresponding to columns of H.
    - `r::AbstractVector`               : Radius values corresponding to rows of H.

    # Keyword Arguments
    - `mask_floor::Float64=0.0`         : Minimum heating value to plot (for better color scaling).
    - `filename::String="tidal_heating_map.png"` : Path to save the heatmap figure.
    - `title_str::String="Hansen norm heating"` : Title for the heatmap figure.

    # Returns
    - `plt` : The `Plots.Plot` object.
    """
    function plot_segment_heating(
        H::AbstractMatrix,
        k_range::AbstractVector,
        r::AbstractVector;
        mask_floor = 0.0,
        filename = "tidal_heating_map.png",
        title_str = "Hansen norm heating"
    )

        # radial shell midpoints
        r_mid = 0.5 .* (r[1:end-1] .+ r[2:end])
        R = maximum(r_mid)

        # convert axes to Float64
        r_mid_f64 = Float64.(r_mid)
        k_f64 = Float64.(k_range)

        # ensure H orientation matches (Nr × Nk)
        Nr = length(r_mid)
        Nk = length(k_range)

        if size(H) == (Nr, Nk)
            H_mat = H
        elseif size(H) == (Nk, Nr)
            H_mat = permutedims(H)
        else
            error("H has incompatible dimensions. Expected ($Nr,$Nk) or ($Nk,$Nr), got $(size(H)).")
        end

        H_f64 = Float64.(H_mat)

        # mask low heating
        H_f64 .= max.(H_f64, mask_floor)

        # remove non-finite frequencies
        valid = findall(isfinite, k_f64)
        isempty(valid) && error("No valid forcing frequencies.")

        k_valid = k_f64[valid]
        H_valid = H_f64[:, valid]

        # heatmap
        plt = heatmap(
            k_valid,
            r_mid_f64 ./ R,
            log10.(H_valid);
            xlabel = "Forcing frequency",
            ylabel = "Radius r / R",
            colorbar_title = "log10(tidal heating)",
            title = title_str,
            aspect_ratio = :auto
        )

        savefig(plt, filename)
        @info "Saved heating heatmap to $filename"

        return plt
    end



end