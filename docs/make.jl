using Documenter
using SelfPropelledVoronoi

makedocs(
    sitename = "SelfPropelledVoronoi",
    format = Documenter.HTML(),
    modules = [SelfPropelledVoronoi],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Overview" => "overview.md",
        "Examples" => "examples.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl.git"
)
