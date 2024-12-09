@testset "Surface Water Balance" begin
    testgrid = CryoGrid.DefaultGrid_2cm
    @testset "interact!" begin
        sub = stripparams(TestGroundLayer(WaterBalance(BucketScheme())))
        top = Top(SurfaceWaterBalance(1e-6u"m/s", 0.0u"m/s"))
        stop = Diagnostics.build_dummy_state(testgrid, top, with_units=true)
        ssub = Diagnostics.build_dummy_state(testgrid, sub, with_units=true)
        # initialize fully saturated
        ssub.sat .= 1.0
        swb = CryoGrid.processes(top)
        water = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(top, swb, stop)
        CryoGrid.computediagnostic!(sub, water, ssub)
        CryoGrid.interact!(top, sub, stop, ssub)
        CryoGrid.computeprognostic!(top, stop)
        @test iszero(ssub.jw[1])
        @test stop.drunoff[1] > zero(eltype(stop.drunoff))
        ssub.sat .= 0.5
        swb = CryoGrid.processes(top)
        water = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(top, swb, stop)
        CryoGrid.computediagnostic!(sub, water, ssub)
        CryoGrid.interact!(top, sub, stop, ssub)
        CryoGrid.computeprognostic!(top, stop)
        @test ssub.jw[1] > zero(eltype(ssub.jw))
        @test iszero(stop.drunoff[1])
    end
end