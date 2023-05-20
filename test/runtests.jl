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

# @testset "Examples" begin
#     examples_dir = joinpath(dirname(Base.current_project()), "examples")
#     @test isdir(examples_dir)
#     for file in readdir(examples_dir)
#         @testset "$file" begin
#             @test try
#                 include(joinpath(examples_dir, file))
#                 true
#             catch ex
#                 @error ex
#                 false
#             end
#         end
#     end
# end
