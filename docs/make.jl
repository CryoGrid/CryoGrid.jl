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
       CryoGrid.Solvers,
       CryoGrid.Snow,
       CryoGrid.Surface,
       CryoGrid.Presets,
       CryoGrid.Diagnostics,
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
                            "Heat Conduction" => "api/heat_conduction.md",
                            "Hydrology" => "api/hydrology.md",
                            "Soils" => "api/soils.md",
                            "Snow" => "api/snow.md",
                            "Surface Energy Balance" => "api/seb.md",
                            "Salt" => "api/salt.md"
                       ],
                       "Tiles" => "api/tiles.md",
                       "Solvers" => "api/solvers.md",
                       "Diagnostics" => "api/diagnostics.md",
                       "Presets" => "api/presets.md",
                ],
                "Contributing" => "contributing.md",
])

deploydocs(repo="github.com/CryoGrid/CryoGrid.jl.git", push_preview=true)
