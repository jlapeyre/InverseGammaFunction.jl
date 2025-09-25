using InverseGammaFunction
using Documenter

DocMeta.setdocmeta!(InverseGammaFunction, :DocTestSetup, :(using InverseGammaFunction); recursive=true)

makedocs(;
    modules=[InverseGammaFunction],
    authors="John Lapeyre",
    sitename="InverseGammaFunction.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
