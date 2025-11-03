using Documenter
using Love

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
)

makedocs(
    sitename="Love.jl",
    format=format,
    pages = [
        "Home" => "index.md",
        "model/index.md",
        "setup.md" ,
        "usage.md",
        "examples/index.md",
        "troubleshooting.md",
        "Development" => "manual/index.md",
        "Related codes" => "ecosystem.md"
    ]
)

deploydocs(
    repo = "github.com/FormingWorlds/Love.jl.git",
    push_preview=true
)

