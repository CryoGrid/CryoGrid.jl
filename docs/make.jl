using CryoGrid
using Documenter

IS_LOCAL = haskey(ENV,"LOCALDOCS") && ENV["LOCALDOCS"] == "true"

const modules = [
       CryoGrid,
       CryoGrid.Utils,
       CryoGrid.Numerics,
       CryoGrid.Physics,
       CryoGrid.Boundaries,
       CryoGrid.HeatConduction,
       CryoGrid.Hydrology,
       CryoGrid.Soils,
       CryoGrid.SEB,
       CryoGrid.Strat,
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
                            "Surface Energy Balance" => "api/seb.md",
                       ],
                       "Stratigraphy" => "api/strat.md",
                       "Diagnostics" => "api/diagnostics.md",
                       "Presets" => "api/presets.md",
                ],
                "Contributing" => "contributing.md",
])

deploydocs(repo="github.com/CryoGrid/CryoGrid.jl.git", push_preview=true)
