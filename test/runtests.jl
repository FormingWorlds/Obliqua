#!/usr/bin/env -S julia --color=yes --startup-file=no

# Get Obliqua root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"

# Include libraries
using LoggingExtras
using Obliqua
using Glob
using Test

@info "Begin Obliqua tests"

# Configure
LOG_LEVEL = Logging.Error
SLOW_TESTS = ["obliqua"]

# Prepare
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

total  = 0
failed = 0

# which test suite to run?
suite::String = "all"
if length(ARGS)>0
    suite = strip(ARGS[1])
    if suite == "0"
        suite = "none"
    end
end
@info "Requested suite '$suite'"

# remove old files
rm(OUT_DIR,force=true,recursive=true)
if !isdir(OUT_DIR) && !isfile(OUT_DIR)
    mkdir(OUT_DIR)
end
# remove stale coverage files recursively across repository
for (dir, _, files) in walkdir(ROOT_DIR)
    for f in files
        if endswith(f, ".cov")
            rm(joinpath(dir, f), force=true)
        end
    end
end
for f in glob("coverage.*", ROOT_DIR)
    rm(f, force=true)
end

# Parameters for test suite
rtol = 1e-3
atol = 1e-18  

# Find test names
test_names = sort([replace(split(basename(f), ".jl")[1], "test_"=>"")
                        for f in glob("test_*.jl", TEST_DIR)])

# Test tidal model imports
if "setup" in test_names
    test_file = joinpath(TEST_DIR, "test_setup.jl")
    @info "Checking tidal model imports '$(test_file)'"
    include(test_file)
    # remove "setup" from test_names to avoid running it again
    filter!(t -> t != "setup", test_names)
end

# Select tests
if suite == "none"
    # no tests
    @info "No tests selected. Exiting."
    exit()

elseif suite == "fast"
    # exclude slow tests
    for slow in SLOW_TESTS
        filter!(t -> !occursin(slow, t), test_names)
    end

elseif suite in test_names
    # run only the specified test suite
    test_names = [suite]

else
    if suite != "all"
        @error "Suite '$suite' not found"
        exit(1)
    end
end

# Collect tests
test_files = String[]
for test_name in test_names
    push!(test_files, joinpath(TEST_DIR, "test_$test_name.jl"))
end
@info "Collected: $(join(test_names, ", "))"
@info "Running tests..."


# Run tests
for test_file in test_files
    # Configure logging to show only error messages during tests
    LoggingExtras.global_logger(Logging.SimpleLogger(LOG_LEVEL))
    @info "Running '$(test_file)'"
    include(test_file)
end

exit(0)

