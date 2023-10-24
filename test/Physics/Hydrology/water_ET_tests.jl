using CryoGrid
using CryoGrid.Hydrology
using Test

include("../../testutils.jl")
include("../../types.jl")

Hydrology.hydraulicproperties(::TestGroundLayer) = HydraulicProperties()

function evapotranspiration_tests()
    layer = stripparams(TestGroundLayer(WaterBalance(BucketScheme(), DampedET())))
    water = processes(layer)
    grid = Grid((0.0:0.05:1.0)u"m")
    Qe = 10.0u"W/m^2"
    state = Diagnostics.build_dummy_state(grid, layer)
    state.Qe .= Qe
    # test that results are valid when total water content is zero
    Hydrology.evapotranspiration!(layer, water, state)
    Hydrology.evapotranspirative_fluxes!(layer, water, state)
    @test all(isfinite.(state.jw))
    @test all(iszero.(state.jw))
    # test with water present (saturated conditions)
    state.θwi .= state.θw .= 0.5
    Hydrology.evapotranspiration!(layer, water, state)
    Hydrology.evapotranspirative_fluxes!(layer, water, state)
    @test all(isfinite.(state.jw_ET))
    @test !all(iszero.(state.jw_ET))
    # check that damping effect is working (fluxes should decrease in magnitude downward)
    @test all(state.jw_ET[2:end] .< state.jw_ET[1:end-1])
end

@testset "Evapotranspiration" begin
    
end
