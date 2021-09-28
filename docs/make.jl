include("../src/funciones.jl")
include("../src/main.jl")

using Documenter, .Utils

#= makedocs(
    modules = [Utils],
    format = Documenter.HTML(),
    sitename = "Proyecto chido",
    authors  = "Alonso Martinez",
    pages = [
        "Home" => "index.md",
    ]
) =#

makedocs(
    modules = [Utils],
    format = Documenter.LaTeX(platform="none"),
    authors  = "Alonso Martinez",
    pages = [
        "Home" => "index.md",
    ]
)
