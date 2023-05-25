using Base: NamedTuple
using CryoGrid
using Test

include("Utils/runtests.jl")
include("Numerics/runtests.jl")
include("IO/runtests.jl")
include("Physics/runtests.jl")
include("Tiles/runtests.jl")
include("Solvers/runtests.jl")
include("Diagnostics/runtests.jl")

# these tests verify that (some of) the example scripts do not error;
# the actual model output is not verified!
@testset "Examples" begin
    examples_dir = joinpath(dirname(Base.current_project()), "examples")
    @test isdir(examples_dir)
    # to run all examples:
    # test_example_scripts = readdir(examples_dir)
    # but this is annoyingly slow, so let's just pick a few
    test_example_scripts = [
        "heat_freeW_lite_implicit.jl",
        "heat_freeW_samoylov.jl",
        "heat_sfcc_constantbc.jl",
        "heat_sfcc_samoylov.jl",
        "readme_examples.jl",
    ]
    for file in test_example_scripts
        @info "Running example script: $file"
        @testset "$file" begin
            @test try
                include(joinpath(examples_dir, file))
                true
            catch ex
                @error ex
                false
            end
        end
    end
end
