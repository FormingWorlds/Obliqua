using Documenter
using DocumenterPages
using Obliqua

# HTML format configuration
format = Documenter.HTML(
    edit_link = "main",  # branch name for "Edit on GitHub"
    prettyurls = get(ENV, "CI", nothing) == "true",
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
        "Compass" => "compass.md",
        PageNode("Tutorials" => "tutorials/index.md", [
            "1 - Loading Data" => "tutorials/loading-data.md", 
            "2 - Configuring Model" => "tutorials/configuration-file.md", 
            "3 - Running Model" => "tutorials/running-model.md", 
            "4 - Plotting" => "tutorials/plotting.md"
            ]
        ),
        "How-to guides" => "how-to-guides/index.md",
        PageNode("Reference" => "reference/index.md", [
            "1 - Forcing Frequency" => "reference/forcing-frequency.md",
            "2 - Rheology" => "reference/rheology.md", 
            PageNode("3 -Solid-phase" => "reference/solid-phase.md", [
                "Solid0d" => "reference/solid/solid0d.md",
                "Solid1d" => "reference/solid/solid1d.md", 
                "Solid1d-mush" => "reference/solid/solid1d_mush.md", 
                "Solid1d-relax" => "reference/solid/solid1d_relax.md", 
                "Solid1d-mush-relax" => "reference/solid/solid1d_mush_relax.md"
                ]
            ),
            "4 - Mush layer" => "reference/mush-layer.md", 
            "5 - Liquid-phase" => "reference/liquid-phase.md",
            "6 - Surface Loading" => "reference/surface-loading.md",
            "7 - Tidal potential" => "reference/tidal-potential.md"
            ]
        ),
        PageNode("Explanation" => "explanation/index.md", [
            "1 - Technical overview" => "explanation/technical_overview.md",
            "2 - Main loop" => "explanation/main_loop.md",
            "3 - Iterative process" => "explanation/iter_proc.md",
            "4 - Post-processing" => "explanation/post_proc.md"
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

