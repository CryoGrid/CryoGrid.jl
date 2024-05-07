using Test

@testset "Hydrology" begin
    include("water_balance_tests.jl")
    include("water_ET_tests.jl")
end
