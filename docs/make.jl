using CryoGrid, Documenter

makedocs(modules=[CryoGrid],
         sitename="CryoGrid.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         pages=["Home" => "index.md",
                "Library" => [
                       "Core" => "api/core.md"
                ],
])
