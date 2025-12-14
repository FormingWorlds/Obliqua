using Documenter
using fwlLove

# HTML format configuration
format = Documenter.HTML(
    edit_link = "main",  # branch name for "Edit on GitHub"
    prettyurls = get(ENV, "CI", nothing) == "true",
)

# Build the docs
makedocs(
    sitename="Love.jl",
    format=format,
    modules = [fwlLove],
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Compass" => "compass.md",
        PageNode("Tutorials" => "tutorials/index.md", [
            "chapter 1 - Loading Data" => "tutorials/ch1.md", 
            "chapter 2 - Running Model" => "tutorials/ch2.md", 
            "chapter 3 - Plotting" => "tutorials/ch3.md"
            ]
        ),
        "How-to guides" => "how-to-guides/index.md",
        PageNode("Reference" => "reference/index.md", [
            "chapter 1 - Rheology" => "reference/ch1.md", 
            "chapter 2 - Solid-phase" => "reference/ch2.md", 
            "chapter 3 - Mush layer" => "reference/ch3.md", 
            "chapter 4 - Liquid-phase" => "reference/ch4.md",
            "chapter 5 - Forcing Frequency" => "reference/ch5.md",
            "chapter 6 - Surface Loading" => "reference/ch6.md"
            ]
        ),
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

