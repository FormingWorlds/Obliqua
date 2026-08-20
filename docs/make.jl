using Documenter
using DocumenterPages
using DocumenterTools: Themes
using Obliqua

ASSETS_DIR = joinpath(@__DIR__, "src", "assets")

# following https://github.com/JuliaMusic/JuliaMusic_documentation.jl/blob/master/docs/make.jl
# combine style and defs files into single scss files for compilation
for w in ("light", "dark")
    style = read(joinpath(ASSETS_DIR, "style.scss"), String)
    theme = read(joinpath(ASSETS_DIR, "$(w)defs.scss"), String)
    write(joinpath(ASSETS_DIR, "$(w).scss"), style*"\n"*theme)
end

# compile styles into scss files
Themes.compile(joinpath(ASSETS_DIR, "light.scss"),
                joinpath(ASSETS_DIR, "themes/documenter-light.css"))
Themes.compile(joinpath(ASSETS_DIR, "dark.scss"),
                joinpath(ASSETS_DIR, "themes/documenter-dark.css"))

format = Documenter.HTML(   
    edit_link = nothing,
    collapselevel = 1,
    size_threshold = 300 * 1024,
    size_threshold_warn = 200 * 1024,
    prettyurls = get(ENV, "CI", nothing) == "true",
    assets = [
        # local assets
        "assets/style.css",
        "assets/logo.ico",
    ]
)

# Build the docs
makedocs(
    sitename="Obliqua",
    checkdocs = :none,
    format=format,
    modules = [Obliqua],
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        PageNode("Tutorials" => "tutorials/index.md", [
            "1 - Loading Data" => "tutorials/loading-data.md", 
            "2 - Configuring Model" => "tutorials/configuration-file.md", 
            "3 - Running Model" => "tutorials/running-model.md", 
            "4 - Plotting" => "tutorials/plotting.md"
            ]
        ),
        "How-to guides" => "how-to-guides/index.md",
        PageNode("Reference" => "reference/index.md", [
            "0 - Configuring Model" => "reference/configuration-file.md", 
            "1 - Forcing Frequency" => "reference/forcing-frequency.md",
            "2 - Rheology" => "reference/rheology.md", 
            PageNode("3 -Solid-phase" => "reference/solid-phase.md", [
                "Solid0d" => "reference/solid/solid0d.md",
                "Solid1d" => "reference/solid/solid1d.md", 
                "Solid1d-mush" => "reference/solid/solid1d_mush.md", 
                "Solid1d-relax" => "reference/solid/solid1d_relax.md", 
                "Solid1d-mush-relax" => "reference/solid/solid1d_mush_relax.md",
                "Solid1d-equil-relax" => "reference/solid/solid1d_equil_relax.md", 
                ]
            ),
            "4 - Mush layer" => "reference/mush-layer.md", 
            "5 - Liquid-phase" => "reference/liquid-phase.md",
            "6 - Surface Loading" => "reference/surface-loading.md",
            "7 - Tidal potential" => "reference/tidal-potentials.md"
            ]
        ),
        PageNode("Explanation" => "explanation/index.md", [
            "1 - Technical overview" => "explanation/technical_overview.md",
            "2 - Main loop" => "explanation/main_loop.md",
            "3 - Iterative process" => "explanation/iter_proc.md",
            "4 - Post-processing" => "explanation/post_proc.md",
            "5 - Stabilization techniques" => "explanation/stabilize.md"
            ]
        ),
        "Troubleshooting" => "troubleshooting.md",
        "Development" => "development.md",
        "Related codes" => "ecosystem.md"
    ]
)

deploydocs(
    repo = "https://github.com/FormingWorlds/Obliqua.git",
    push_preview=true
)

