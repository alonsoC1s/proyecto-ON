include("../src/Solvers.jl")

# using Documenter, .Solvers.Utils, .Solvers

makedocs(
    modules = [Solvers.Utils, Solvers],
    format = Documenter.HTML(),
    sitename = "Optimización Numérica - Proyecto 1",
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
