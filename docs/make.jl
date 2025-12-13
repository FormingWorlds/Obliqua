using Documenter
using Love

# HTML format configuration
format = Documenter.HTML(
    edit_link = "main",  # branch name for "Edit on GitHub"
    prettyurls = get(ENV, "CI", nothing) == "true",
)

# Build the docs
makedocs(
    sitename="Love.jl",
    format=format,
    modules = [Love],
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Compass" => "compass.md",
        "Tutorials" => "tutorials/index.md",
        "How-to guides" => "how-to-guides/index.md",
        "Reference" => "reference/index.md",
        "Explanation" => "explanation/index.md",
        "Troubleshooting" => "troubleshooting.md",
        "Development" => "development.md",
        "Related codes" => "ecosystem.md"
    ]
)

deploydocs(
    repo = "https://github.com/FormingWorlds/Love.jl.git",
    push_preview=true
)

