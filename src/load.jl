# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module load

    using JSON
    using LoggingExtras
    using Printf

    export load_interior

    """
        load_interior(fname::String)

    Load interior structure + tidal parameters from a JSON file.

    The JSON file must contain:

    {
        "density": [...],
        "radius": [...],
        "visc": [...],
        "shear": [...],
        "bulk": [...],
        "omega": value,
        "ecc": value,
        "ncalc": value
    }

    Returns:
        (omega, ecc, rho, radius, visc, shear, bulk, ncalc)
    """
    function load_interior(fname::String)

        # Convert to absolute path
        fname = abspath(fname)
        @info "Loading interior from JSON file: $fname"

        # Suppress debug spam
        with_logger(MinLevelLogger(current_logger(), Logging.Info)) do

            # Parse JSON
            params = JSON.parsefile(fname)

            # ----------------------
            # Load vector quantities
            rho    = BigFloat.(params["density"])
            radius = BigFloat.(params["radius"])
            visc   = BigFloat.(params["visc"])
            shear  = BigFloat.(params["shear"])
            bulk   = BigFloat.(params["bulk"])

            # Basic consistency
            N = length(rho)
            if !(length(radius) == N == length(visc) == length(shear) == length(bulk))
                error("JSON input arrays must all have the same length. Got sizes: " *
                      "rho=$(length(rho)), radius=$(length(radius)), visc=$(length(visc)), " *
                      "shear=$(length(shear)), bulk=$(length(bulk))")
            end

            # ----------------------
            # Load scalar quantities
            omega = BigFloat(params["omega"])
            ecc   = BigFloat(params["ecc"])
            ncalc = Int(params["ncalc"])

        end # logger suppressed

        return omega, ecc, rho, radius, visc, shear, bulk, ncalc
    end

end # module
