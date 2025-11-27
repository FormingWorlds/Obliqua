# Core file containing functions for running the model

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module MAIN

    # Include system libraries
    using LoggingExtras
    using Printf
    import TOML:parsefile

    # Include local jl files (order matters)
    include("Love.jl")
    include("Fluid.jl")
    include("load.jl")


    # Import submodules
    import .Love
    import .Liquid
    import .load

    # Export submodules (mostly for autodoc purposes)
    export Love
    export Liquid
    export load
end