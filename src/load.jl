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

    prec = BigFloat

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
    function load_interior(fname::String, verify::Bool=false
        )::Union{
            Bool,
            Tuple{
                prec, prec,
                Array{prec,1}, Array{prec,1}, Array{prec,1},
                Array{prec,1}, Array{prec,1}, Int
            }
        }

        # Convert to absolute path
        fname = abspath(fname)
        @info "Loading interior from JSON file: $fname"

        # Suppress debug spam
        omega = ecc = ncalc = nothing
        rho = radius = visc = shear = bulk = nothing

        with_logger(MinLevelLogger(current_logger(), Logging.Info)) do

            # --- Parse JSON ---
            params = JSON.parsefile(fname)

            # --- Load vector quantities ---
            rho    = prec.(params["density"])
            radius = prec.(params["radius"])
            visc   = prec.(params["visc"])
            shear  = prec.(params["shear"])
            bulk   = prec.(params["bulk"])

            # --- Load scalar quantities ---
            omega = prec(params["omega"])
            ecc   = prec(params["ecc"])
            ncalc = Int(params["ncalc"])

        end

        if verify
            ok = true

            # internal consistency
            N = length(rho)
            ok &= (length(radius) == N+1)
            ok &= (length(visc)   == N)
            ok &= (length(shear)  == N)
            ok &= (length(bulk)   == N)

            # # simple physical checks 
            ok &= all(rho .> 0)
            ok &= all(radius .> 0)
            ok &= all(visc .>= 0)
            ok &= all(shear .>= 0)
            ok &= all(bulk .>= 0)

            ok &= (omega ≥ 0)
            ok &= (ecc ≥ 0 && ecc < 1)

            return ok
        end

        return omega, ecc, rho, radius, visc, shear, bulk, ncalc
    end

    function load_interior_mush(fname::String, verify::Bool=false
        )::Union{
            Bool,
            Tuple{
                prec, prec,
                Array{prec,1}, Array{prec,1}, Array{prec,1},
                Array{prec,1}, Array{prec,1}, Array{prec,1}, 
                Int
            }
        }

        # Convert to absolute path
        fname = abspath(fname)
        @info "Loading interior from JSON file: $fname"

        # Suppress debug spam
        omega = ecc = ncalc = nothing
        rho = radius = visc = shear = bulk = phi = nothing

        with_logger(MinLevelLogger(current_logger(), Logging.Info)) do

            # --- Parse JSON ---
            params = JSON.parsefile(fname)

            # --- Load vector quantities ---
            rho    = prec.(params["density"])
            radius = prec.(params["radius"])
            visc   = prec.(params["visc"])
            shear  = prec.(params["shear"])
            bulk   = prec.(params["bulk"])
            phi    = prec.(params["phi"])

            # --- Load scalar quantities ---
            omega = prec(params["omega"])
            ecc   = prec(params["ecc"])
            ncalc = Int(params["ncalc"])

        end

        if verify
            ok = true

            # internal consistency
            N = length(rho)
            ok &= (length(radius) == N+1)
            ok &= (length(visc)   == N)
            ok &= (length(shear)  == N)
            ok &= (length(bulk)   == N)
            ok &= (length(phi)    == N)

            # # simple physical checks 
            ok &= all(rho .> 0)
            ok &= all(radius .> 0)
            ok &= all(visc .>= 0)
            ok &= all(shear .>= 0)
            ok &= all(bulk .>= 0)
            ok &= all(phi .>= 0)

            ok &= (omega ≥ 0)
            ok &= (ecc ≥ 0 && ecc < 1)

            return ok
        end

        return omega, ecc, rho, radius, visc, shear, bulk, phi, ncalc
    end

    function load_interior_liquid(fname::String, verify::Bool=false
        )::Union{
            Bool,
            Tuple{
                prec, prec,
                prec, prec, prec,
                Array{prec,1}, Array{prec,1}, Array{prec,1}
            }
        }

        # Convert to absolute path
        fname = abspath(fname)
        @info "Loading interior from JSON file: $fname"

        # Suppress debug spam
        omega = axial = ecc = sma = S_mass = nothing
        rho = radius = visc = nothing

        with_logger(MinLevelLogger(current_logger(), Logging.Info)) do

            # --- Parse JSON ---
            params = JSON.parsefile(fname)

            # --- Load vector quantities ---
            rho    = prec.(params["density"])
            radius = prec.(params["radius"])
            visc   = prec.(params["visc"])

            # --- Load scalar quantities ---
            omega  = prec(params["omega"])
            axial  = prec(params["axial"])
            ecc    = prec(params["ecc"])
            sma    = prec(params["sma"])
            S_mass = prec(params["S_mass"])

        end

        if verify
            ok = true

            # internal consistency
            N = length(rho)
            ok &= (length(radius) == N+1)
            ok &= (length(visc)   == N)

            # # simple physical checks 
            ok &= all(rho .> 0)
            ok &= all(radius .> 0)
            ok &= all(visc .>= 0)

            ok &= (omega ≥ 0)
            ok &= (ecc ≥ 0 && ecc < 1)
            ok &= (sma > 0)
            ok &= (S_mass > 0)

            return ok
        end

        return omega, axial, ecc, sma, S_mass, rho, radius, visc
    end

    function load_interior_full(fname::String, verify::Bool=false
        )::Union{
            Bool,
            Tuple{
                prec, prec,
                prec, prec, prec,
                Array{prec,1}, Array{prec,1}, Array{prec,1},
                Array{prec,1}, Array{prec,1}, Int
            }
        }

        # Convert to absolute path
        fname = abspath(fname)
        @info "Loading interior from JSON file: $fname"

        # Suppress debug spam
        omega = axial = ecc = sma = S_mass = ncalc = nothing
        rho = radius = visc = shear = bulk = nothing

        with_logger(MinLevelLogger(current_logger(), Logging.Info)) do

            # --- Parse JSON ---
            params = JSON.parsefile(fname)

            # --- Load vector quantities ---
            rho    = prec.(params["density"])
            radius = prec.(params["radius"])
            visc   = prec.(params["visc"])
            shear  = prec.(params["shear"])
            bulk   = prec.(params["bulk"])

            # --- Load scalar quantities ---
            omega  = prec(params["omega"])
            axial  = prec(params["axial"])
            ecc    = prec(params["ecc"])
            sma    = prec(params["sma"])
            S_mass = prec(params["S_mass"])
            ncalc  = Int(params["ncalc"])

        end

        if verify
            ok = true

            # internal consistency
            N = length(rho)
            ok &= (length(radius) == N+1)
            ok &= (length(visc)   == N)
            ok &= (length(shear)  == N)
            ok &= (length(bulk)   == N)

            # # simple physical checks 
            ok &= all(rho .> 0)
            ok &= all(radius .> 0)
            ok &= all(visc .>= 0)
            ok &= all(shear .>= 0)
            ok &= all(bulk .>= 0)

            ok &= (omega ≥ 0)
            ok &= (ecc ≥ 0 && ecc < 1)
            ok &= (sma > 0)
            ok &= (S_mass > 0)

            return ok
        end

        return omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, ncalc

    end

end # module