using ElectronLiquid
using Documenter

DocMeta.setdocmeta!(ElectronLiquid, :DocTestSetup, :(using ElectronLiquid); recursive=true)

makedocs(;
    modules=[ElectronLiquid],
    authors="Kun Chen",
    repo="https://github.com/numericaleft/ElectronLiquid.jl/blob/{commit}{path}#{line}",
    sitename="ElectronLiquid.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/ElectronLiquid.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "SubModules" => [
            "Parameter" => "lib/UEG.md",
            "CounterTerm" => "lib/CounterTerm.md",
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/ElectronLiquid.jl",
    devbranch="master"
)
