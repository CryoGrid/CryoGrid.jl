using CryoGrid
using CryoGrid.Snow
using Test

@testset "Bulk snowpack, prescribed mass balance" begin
    heat = HeatBalance()
    water = WaterBalance()
    mass = Snow.PrescribedSnowMassBalance(swe)
    snow = Snowpack(Snow.Bulk(); heat, water)
    grid = Grid([-2.0u"m", 0.0u"m"])
    dummystate = Diagnostics.build_dummy_state(grid, snow)
    @test all(map(∈(propertynames(dummystate)), [:swe, :dsn, :θwi]))
end
