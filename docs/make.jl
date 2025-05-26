using ArgParse
using CryoGrid
using Documenter
using Literate

s = ArgParseSettings()
@add_arg_table! s begin
    "--local", "-l"
       action = :store_true
       help = "Local docs build mode"
    "--draft", "-d"
       action = :store_true
       help = "Whether to build docs in draft mode, i.e. skipping execution of examples and doctests"
end
parsed_args = parse_args(ARGS, s)

IS_LOCAL = parsed_args["local"] || parse(Bool, get(ENV, "LOCALDOCS", "false"))
IS_DRAFT = parsed_args["draft"] || parse(Bool, get(ENV, "DRAFTDOCS", "false"))
if haskey(ENV, "GITHUB_ACTIONS")
       ENV["JULIA_DEBUG"] = "Documenter"
end

deployconfig = Documenter.auto_detect_deploy_system()

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
       CryoGrid.Diagnostics,
       # solvers
       CryoGrid.DiffEq,
       CryoGrid.LiteImplicit
];

examples_dir = joinpath(@__DIR__, "..", "examples")
examples_output_dir = joinpath(@__DIR__, "src", "examples")
# remove existing files
rm(examples_output_dir, recursive=true, force=true)
# recreate directory
mkpath(examples_output_dir)
# define example files to make pages from
name_lookup = Dict(
       "heat_freeW_samoylov.md" => "Soil heat with free water freeze curve",
       "heat_sfcc_constantbc.md" => "Soil heat with SFCC and constant BCs",
       "heat_sfcc_samoylov.md" => "Soil heat with SFCC",
       "heat_freeW_snow_samoylov.md" => "Soil heat with bulk snow scheme",
       "heat_freeW_bucketW_samoylov.md" => "Soil heat with bucket water scheme",
       "heat_freeW_seb_snow_bucketW_samoylov.md" => "Soil heat w/ SEB, snow cover, and bucket water scheme",
       "heat_freeW_lite_implicit.md" => "Fast heat conduction with CryoGridLite",
       "heat_sfcc_richardseq_samoylov.md" => "Coupled soil heat and water transport",
       "heat_sfcc_salt_constantbc.md" => "Coupled heat and salt diffusion on salty soil column",
       "heat_simple_autodiff_grad.md" => "Computing parameter sensitivities with autodiff",
       "cglite_parameter_ensembles.md" => "Running parameter ensembles",
)
# change file suffixes
example_scripts = map(f -> split(f, ".")[1]*".jl", collect(keys(name_lookup)))
# generate example docs from scripts
example_docfiles = map(filter(∈(example_scripts), readdir(examples_dir))) do f
       infile = joinpath(examples_dir, f)
       @info "Generating docpage for example script $infile and writing to directory $examples_output_dir"
       Literate.markdown(infile, examples_output_dir, execute=false, documenter=true)
       return f
end

example_docpages = map(example_docfiles) do f
       docpage = replace(f, "jl" => "md")
       name_lookup[docpage] => joinpath("examples", docpage)
end

makedocs(
       repo=Remotes.GitHub("CryoGrid", "CryoGrid.jl"),
       modules=modules,
       sitename="CryoGrid.jl",
       authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
       format=Documenter.HTML(
              prettyurls=!IS_LOCAL,
              canonical = "https://cryogrid.github.io/CryoGrid.jl/v0",
       ),
       warnonly=true, # don't fail when there are errors
       doctest=!IS_DRAFT,
       draft=IS_DRAFT,
       pages=["Home" => "index.md",
              "Installation" => "installation.md",
              "Getting Started" => "quickstart.md",
              "User manual" => [
                     "Overview" => "manual/overview.md",
                     "Architecture" => "manual/architecture.md",
                     "Coupling layers and processes" => "manual/coupling.md",
              ],
              "Developer guide" => [
              "Concepts" => "dev/concepts.md",
              "Debugging" => "dev/debugging.md",
              # "Contributing" => "dev/contributing.md",
              ],
              "Examples" => example_docpages,
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
       ],
])

# remove gitignore from build files
rm(joinpath(@__DIR__, "build", ".gitignore"))

deploydocs(
       repo="github.com/CryoGrid/CryoGrid.jl.git",
       push_preview = true,
       versions = ["v0" => "v^", "v#.#", "dev" => "dev"],
       deploy_config = deployconfig,
)
