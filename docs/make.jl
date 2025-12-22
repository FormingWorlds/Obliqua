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
    format=format,
    modules = [Obliqua],
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Compass" => "compass.md",
        PageNode("Tutorials" => "tutorials/index.md", [
            "1 - Loading Data" => "tutorials/loading-data.md", 
            "2 - Running Model" => "tutorials/running-model.md", 
            "3 - Plotting" => "tutorials/plotting.md"
            ]
        ),
        "How-to guides" => "how-to-guides/index.md",
        PageNode("Reference" => "reference/index.md", [
            "1 - Rheology" => "reference/rheology.md", 
            "2 - Solid-phase" => "reference/solid-phase.md", 
            "3 - Mush layer" => "reference/mush-layer.md", 
            "4 - Liquid-phase" => "reference/liquid-phase.md",
            "5 - Forcing Frequency" => "reference/forcing-frequency.md",
            "6 - Surface Loading" => "reference/surface-loading.md",
            "7 - Tidal Potentials" => "reference/tidal-potentials.md"
            ]
        ),
        "Explanation" => "explanation/index.md",
        "Troubleshooting" => "troubleshooting.md",
        "Development" => "development.md",
        "Related codes" => "ecosystem.md"
    ]
)

deploydocs(
    repo = "https://github.com/FormingWorlds/Obliqua.git",
    push_preview=true
)

