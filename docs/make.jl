using Documenter
using Pkg

ENV["PYTHON"] = ""
Pkg.build("PyCall")

using IMinuit


makedocs(;
    modules=[IMinuit],
    authors="Feng-Kun Guo",
    repo="https://github.com/fkguo/IMinuit.jl/blob/{commit}{path}#L{line}",
    sitename="IMinuit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fkguo.github.io/IMinuit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "iminuit functions" => "iminuit.md",
        "Example" => "example.md",
        "Data and plot macros" => "data.md",
        "chisq and model_fit" => "modelfit.md",
        "Error analysis" => "contour.md"
    ],
)

deploydocs(;
    repo="github.com/fkguo/IMinuit.jl"
    # target = "build",
    # deps = nothing,
    # make = nothing,
    # branch = "gh-pages"
)
