using CryoGrid
using Documenter

IS_LOCAL = haskey(ENV,"LOCALDOCS") && ENV["LOCALDOCS"] == "true"

const modules = [
       CryoGrid,
       CryoGrid.Utils,
       CryoGrid.Numerics,
       CryoGrid.InputOutput,
       CryoGrid.Tiles,
       CryoGrid.Heat,
       CryoGrid.Hydrology,
       CryoGrid.Soils,
       CryoGrid.Snow,
       CryoGrid.Surface,
       CryoGrid.Presets,
       CryoGrid.Diagnostics,
       # solvers
       CryoGrid.DiffEq,
       CryoGrid.LiteImplicit
]

makedocs(modules=modules,
         sitename="CryoGrid.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         format=Documenter.HTML(prettyurls=!IS_LOCAL),
         pages=["Home" => "index.md",
                "Installation" => "installation.md",
                "Getting Started" => "quickstart.md",
                "Manual" => [
                       "Overview" => "manual/overview.md",
                ],
                "Library" => [
                       "Index" => "api/index.md",
                       "Method interface" => "api/toplevel.md",
                       "Numerics" => "api/numerics.md",
                       "Utilities" => "api/utils.md",
                       "Physics" => [
                            "Heat Conduction" => "api/physics/heat_conduction.md",
                            "Hydrology" => "api/physics/hydrology.md",
                            "Soils" => "api/physics/soils.md",
                            "Snow" => "api/physics/snow.md",
                            "Surface Energy Balance" => "api/physics/seb.md",
                            "Salt" => "api/physics/salt.md"
                       ],
                       "Tiles" => "api/tiles.md",
                       "Solvers" => [
                            "SciML/DiffEq" => "api/solvers/diffeq.md",
                            "CryoGridLite" => "api/solvers/lite_implicit.md",
                       ],
                       "Diagnostics" => "api/diagnostics.md",
                       "Presets" => "api/presets.md",
                ],
                "Contributing" => "contributing.md",
])

deploydocs(repo="github.com/CryoGrid/CryoGrid.jl.git", push_preview=true)
