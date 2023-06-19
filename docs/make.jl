using CryoGrid
using Documenter
using Literate

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
];

examples_dir = joinpath(@__DIR__, "..", "examples")
examples_output_dir = joinpath(@__DIR__, "src", "examples")
# remove existing files
rm(examples_output_dir, recursive=true)
# recreate directory
mkpath(examples_output_dir)
# generate example docs from scripts
example_docs = map(readdir(examples_dir)) do f
       infile = joinpath(examples_dir, f)
       Literate.markdown(infile, examples_output_dir)
       return joinpath("examples", replace(f, "jl" => "md"))
end

makedocs(modules=modules,
         sitename="CryoGrid.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         format=Documenter.HTML(prettyurls=!IS_LOCAL),
         pages=["Home" => "index.md",
                "Installation" => "installation.md",
                "Getting Started" => "quickstart.md",
                "User manual" => [
                       "Overview" => "manual/overview.md",
                       "Architecture" => "manual/architecture.md",
                       "Model interface" => "manual/interface.md",
                       "Coupling layers and processes" => "manual/coupling.md",
                ],
                "Examples" => example_docs,
                "Developer guide" => [
                     "Concepts" => "dev/concepts.md",
                     "Debugging" => "dev/debugging.md",
                     "Contributing" => "dev/contributing.md",
                ],
                "API" => [
                     "Index" => "api/index.md",
                     "CryoGrid" => "api/toplevel.md",
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
                          "Built-in" => "api/solvers/basic_solvers.md",
                          "SciML/DiffEq" => "api/solvers/diffeq.md",
                          "CryoGridLite" => "api/solvers/lite_implicit.md",
                     ],
                     "Diagnostics" => "api/diagnostics.md",
                     "Presets" => "api/presets.md",
              ],
])

deploydocs(repo="github.com/CryoGrid/CryoGrid.jl.git", push_preview=true)
