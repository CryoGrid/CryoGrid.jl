@testset "Surface Water Balance" begin
    testgrid = CryoGrid.Presets.DefaultGrid_2cm
    @testset "interact!" begin
        sub = TestGroundLayer(WaterBalance(BucketScheme()))
        top = Top(
            SurfaceWaterBalance(
                ConstantForcing(1e-6u"m/s", :rainfall),
                ConstantForcing(0.0u"m/s", :snowfall),
            )
        )
        stop = Diagnostics.build_dummy_state(testgrid, top, with_units=true)
        ssub = Diagnostics.build_dummy_state(testgrid, sub, with_units=true)
        # initialize fully saturated
        ssub.sat .= 1.0
        swb = CryoGrid.processes(top)
        water = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(top, swb, stop)
        CryoGrid.computediagnostic!(sub, water, ssub)
        CryoGrid.interact!(top, sub, stop, ssub)
        CryoGrid.computefluxes!(top, stop)
        @test iszero(ssub.jw[1])
        @test stop.∂runoff∂t[1] > zero(eltype(stop.∂runoff∂t))
        ssub.sat .= 0.5
        swb = CryoGrid.processes(top)
        water = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(top, swb, stop)
        CryoGrid.computediagnostic!(sub, water, ssub)
        CryoGrid.interact!(top, sub, stop, ssub)
        CryoGrid.computefluxes!(top, stop)
        @test ssub.jw[1] > zero(eltype(ssub.jw))
        @test iszero(stop.∂runoff∂t[1])
    end
end