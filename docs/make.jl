using CryoGrid
using Documenter

IS_LOCAL = haskey(ENV,"LOCALDOCS") && ENV["LOCALDOCS"] == "true"

const modules = [
       CryoGrid,
       CryoGrid.Utils,
       CryoGrid.Numerics,
       CryoGrid.Forcings,
       CryoGrid.Layers,
       CryoGrid.Processes,
       CryoGrid.HeatConduction,
       CryoGrid.SEB,
       CryoGrid.Setup,
       CryoGrid.Models,
       CryoGrid.Callbacks,
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
                       "Core" => [
                              "Numerics" => "api/numerics.md",
                              "Forcings" => "api/forcings.md",
                              "Utilities" => "api/utils.md",
                       ],
                       "Processes" => [
                            "Heat Conduction" => "api/heat_conduction.md",
                            "Surface Energy Balance" => "api/seb.md",
                       ],
                       "Setup" => "api/setup.md",
                       "Callbacks" => "api/callbacks.md",
                       "Models" => "api/models.md",
                ],
                "Contributing" => "contributing.md",
])

deploydocs(repo="github.com/CryoGrid/CryoGrid.jl.git", push_preview=true)
