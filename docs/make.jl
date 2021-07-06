using CryoGrid
using Documenter

const modules = [
       CryoGrid,
       Utils,
       Numerics,
       Forcings,
       Layers,
       Processes,
       HeatConduction,
       SEB,
       Setup,
       Models,
       Callbacks,
]

makedocs(modules=modules,
         sitename="CryoGrid.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         format=Documenter.HTML(prettyurls=false),
         pages=["Home" => "index.md",
                "installation" => "installation.md",
                "Getting Started" => "quickstart.md",
                "Library" => [
                       "Overview" => "api/overview.md",
                       "Setup" => "api/setup.md",
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
                       "Callbacks" => "api/callbacks.md",
                       "Models" => "api/models.md",
                ],
])
