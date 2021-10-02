include("../src/Utils.jl")
include("../src/Solvers.jl")

using Documenter, .Utils, .Solvers

makedocs(
    modules = [Utils, Solvers],
    format = Documenter.HTML(),
    sitename = "Optimización Numéria - Proyecto 1",
    authors  = "Alonso Martinez",
    pages = [
        "Home" => "index.md",
    ]
)

# makedocs(
#     modules = [Utils],
#     format = Documenter.LaTeX(platform="none"),
#     authors  = "Alonso Martinez",
#     pages = [
#         "Home" => "index.md",
#     ]
# )
