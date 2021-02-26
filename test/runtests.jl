using CryoGrid
using Test

@testset "CryoGrid" begin
    include("core/runtests.jl")
    include("processes/runtests.jl")
end
