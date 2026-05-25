#!/usr/bin/env -S julia --color=yes --startup-file=no
# Run this tool from inside the Obliqua root folder
# e.g. as `julia --project test/runprofiler.jl`

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
cd(ROOT_DIR)

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

using Profile
using Obliqua
using Printf

function _print_help()
    println("Obliqua Profiler Tool")
    println("=====================")
    println("Profiles an Obliqua simulation run and generates an interactive flame graph.")
    println("\nUsage:")
    println("  julia --project test/runprofiler.jl [config_path] [report_dir]")
end

function _load_profiler_library()
    try
        @eval using ProfileCanvas
        return true
    catch err
        @error "ProfileCanvas.jl missing. Install with: julia --project -e 'using Pkg; Pkg.add(\"ProfileCanvas\")'"
        return false
    end
end


"""
    print_top_bottlenecks(data, lidict, root_dir; max_lines=10)

Parses raw profile buffers, isolates files living in `src/`, 
and prints a ranked leaderboard of resource utilization.
"""
function print_top_bottlenecks(data, lidict, root_dir; max_lines=10)
    # Count occurrences of instruction pointers
    counts = Dict{UInt64, Int}()
    for ip in data
        counts[ip] = get(counts, ip, 0) + 1
    end
    
    total_samples = length(data)
    total_samples == 0 && return
    
    # Structure data into lookups
    ranked_functions = []
    for (ip, count) in counts
        frame = get(lidict, ip, nothing)
        if frame !== nothing
            file_str = string(frame.file)
            
            # Filter to focus strictly on package business logic
            if contains(file_str, joinpath(root_dir, "src"))
                short_file = relpath(file_str, root_dir)
                raw_func = string(frame.func)
                
                # CLEANUP: Remove internal compiler tokens (e.g., #hansen_fft#3 -> hansen_fft)
                func_name = replace(raw_func, r"^#([^#]+)#\d+$" => s"\1")
                if startswith(func_name, "#")
                    func_name = "anonymous closure (" * func_name * ")"
                end

                push!(ranked_functions, (
                    count = count,
                    pct = (count / total_samples) * 100,
                    func = func_name,
                    line = frame.line,
                    file = short_file
                ))
            end
        end
    end
    
    # Sort by sample counts descending
    sort!(ranked_functions, by = x -> x.count, rev=true)
    
    # Deduplicate matching function/line combinations
    seen = Set{String}()
    unique_ranks = []
    for item in ranked_functions
        key = "$(item.file):$(item.func):$(item.line)"
        if !(key in seen)
            push!(seen, key)
            push!(unique_ranks, item)
        end
    end
    
    lines_to_show = min(max_lines, length(unique_ranks))
    
    println("\n" * "🚨" * " ="^28 * "🚨")
    println("         TOP $lines_to_show OBLIQUA COMPUTE BOTTLENECKS (Pure CPU Time)")
    println("="^65)
    println(@sprintf("| %-7s | %-32s | %-16s |", "Share %", "Function Name", "Source Line"))
    println("-"^65)
    
    for i in 1:lines_to_show
        item = unique_ranks[i]
        
        # Enforce column string widths safely
        func_label = length(item.func) > 30 ? item.func[1:27] * "..." : item.func
        location = "$(basename(item.file)):$(item.line)"
        loc_label = length(location) > 16 ? "..." * location[end-13:end] : location
        
        # FIXED: Explicitly wrap the string macro in a println expression
        println(@sprintf("| %5.1f%%  | %-32s | %-16s |", item.pct, func_label, loc_label))
    end
    println("="^65)
end


function do_profiling(cfg_path, report_dir)
    isfile(cfg_path) || error("Config not found: $cfg_path")
    mkpath(report_dir)

    cfg = Obliqua.open_config(cfg_path)
    json_path = joinpath(ROOT_DIR, "res", "interior_data", "test_mantle_mush_full_test.json")
    
    omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
        load.load_interior_mush_full(json_path, false)

    @info "Executing warm-up run to isolate JIT compilation overhead..."
    try
        Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg)
    catch; end

    @info "Starting Profile Capture..."
    Profile.clear()
    Profile.Allocs.clear()
    
    Profile.Allocs.@profile sample_rate=0.01 begin
        Profile.@profile begin
            power_prf, power_blk, σ_range, imag_k2 = Obliqua.run_tides(
                omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, cfg
            )
        end
    end

    raw_data = Profile.fetch()
    prof_atoms = Profile.getdict(raw_data) # Explicitly look up backtrace symbols
    prof_data, prof_atoms = Profile.flatten(raw_data, prof_atoms) # Pass both to flatten

    # 1. Export HTML Canvas Flame graph
    output_html = joinpath(report_dir, "profile.html")
    ProfileCanvas.html_file(output_html, ProfileCanvas.view())
    
    # 2. Print Direct Clean Analytics Leaderboard
    print_top_bottlenecks(prof_data, prof_atoms, ROOT_DIR, max_lines=10)
    
    @info "Full interactive flame graph exported to: $output_html"
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    any(arg -> arg in ("-h", "--help"), ARGS) && (_print_help(); exit(0))
    _load_profiler_library() || exit(1)

    cfg_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(ROOT_DIR, "test", "test.toml")
    report_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(ROOT_DIR, "profile_report")

    do_profiling(cfg_path, report_dir)
end