using GeometricSolutions
using Documenter

DocMeta.setdocmeta!(GeometricSolutions, :DocTestSetup, :(using GeometricSolutions); recursive=true)

makedocs(;
    modules=[GeometricSolutions],
    authors="Michael Kraus",
    repo="https://github.com/JuliaGNI/GeometricSolutions.jl/blob/{commit}{path}#{line}",
    sitename="GeometricSolutions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGNI.github.io/GeometricSolutions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/GeometricSolutions.jl",
    devbranch="main",
)
