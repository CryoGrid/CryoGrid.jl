using CryoGrid
using Documenter

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
         format=Documenter.HTML(prettyurls=false),
         pages=["Home" => "index.md",
                "Installation" => "installation.md",
                "Getting Started" => "quickstart.md",
                "Library" => [
                       "Overview" => "api/overview.md",
                       "Common" => [
                              "Interface" => "api/interface.md",
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
])
