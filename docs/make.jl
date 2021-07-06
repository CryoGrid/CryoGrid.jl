using CryoGrid
using CryoGrid: Common, HeatConduction, Layers, Setup
using Documenter

makedocs(modules=[CryoGrid, Common, HeatConduction, Layers, Setup],
         sitename="CryoGrid.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         format=Documenter.HTML(prettyurls=false),
         pages=["Home" => "index.md",
                "installation" => "installation.md",
                "Getting Started" => "quickstart.md",
                "Library" => [
                       "Interface" => "api/interface.md",
                       "Common" => "api/common.md",
                       "Heat Conduction" => "api/heat_conduction.md",
                ],
])
