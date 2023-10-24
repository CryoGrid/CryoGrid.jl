using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Hydrology
using CryoGrid.Utils

using Test

include("../../testutils.jl")
include("../../types.jl")

@testset "Hydrology" begin
    include("water_balance_tests.jl")
    include("water_ET_tests.jl")
end
