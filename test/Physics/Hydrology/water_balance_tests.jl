@testset "Bucket scheme" begin
    testgrid = CryoGrid.Presets.DefaultGrid_2cm
    @testset "variables" begin
        sub = TestGroundLayer(WaterBalance(BucketScheme()))
        vars = CryoGrid.variables(sub)
        prognostic_vars = filter(CryoGrid.isprognostic, vars)
        diagnostic_vars = filter(CryoGrid.isdiagnostic, vars)
        # check for core water variables
        @test :sat ∈ map(varname, prognostic_vars)
        @test :jw ∈ map(varname, diagnostic_vars)
        @test :jw_v ∈ map(varname, diagnostic_vars)
        @test :jw_ET ∈ map(varname, diagnostic_vars)
        @test :kw ∈ map(varname, diagnostic_vars)
        # check processes
        proc = CryoGrid.processes(sub)
        @test isa(proc, WaterBalance)
    end
    @testset "computediagnostic!" begin
        sub = TestGroundLayer(WaterBalance(BucketScheme()))
        state = Diagnostics.build_dummy_state(testgrid, sub, with_units=true)
        # initialize fully saturated
        state.sat .= 1.0
        procs = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(sub, procs, state)
        @test allfinite(state.kw)
        @test allfinite(state.θw)
    end
    @testset "computefluxes!" begin
        sub = TestGroundLayer(WaterBalance(BucketScheme()))
        state = Diagnostics.build_dummy_state(testgrid, sub, with_units=true)
        # initialize fully saturated
        state.sat .= 1.0
        procs = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(sub, procs, state)
        CryoGrid.computefluxes!(sub, procs, state)
        @test allfinite(state.kw)
        @test allfinite(state.θw)
        @test all(iszero.(state.jw))
        # initialize at half saturation
        state.sat .= 0.5
        procs = CryoGrid.processes(sub)
        CryoGrid.computediagnostic!(sub, procs, state)
        CryoGrid.computefluxes!(sub, procs, state)
        @test allfinite(state.kw)
        @test allfinite(state.θw)
        @test all(state.jw[2:end-1] .> zero(eltype(state.jw)))
    end
    @testset "interact!" begin
        sub1 = TestGroundLayer(WaterBalance(BucketScheme()))
        sub2 = TestGroundLayer(WaterBalance(BucketScheme()))
        state1 = Diagnostics.build_dummy_state(testgrid[0.0u"m"..10.0u"m"], sub1, with_units=true)
        state2 = Diagnostics.build_dummy_state(testgrid[10.0u"m"..1000.0u"m"], sub2, with_units=true)
        # case 1: both saturated
        state1.sat .= 1.0
        state2.sat .= 1.0
        procs = CryoGrid.processes(sub1)
        CryoGrid.computediagnostic!(sub1, procs, state1)
        CryoGrid.computediagnostic!(sub2, procs, state2)
        CryoGrid.interact!(sub1, sub2, state1, state2)
        # check that layer boundary flux is zero due to full saturation
        @test iszero(state1.jw[end])
        @test iszero(state2.jw[1])
        # case 2: upper saturated
        state1.sat .= 1.0
        state2.sat .= 0.5
        water = CryoGrid.processes(sub1)
        CryoGrid.computediagnostic!(sub1, water, state1)
        CryoGrid.computediagnostic!(sub2, water, state2)
        CryoGrid.interact!(sub1, sub2, state1, state2)
        # check that layer boundary flux is zero due to full saturation
        @test state1.jw[end] > zero(eltype(state1.jw))
        @test state1.jw[end] > zero(eltype(state1.jw))
        # case 2: lower saturated
        state1.sat .= 0.5
        state2.sat .= 1.0
        procs = CryoGrid.processes(sub1)
        CryoGrid.computediagnostic!(sub1, water, state1)
        CryoGrid.computediagnostic!(sub2, water, state2)
        CryoGrid.interact!(sub1, sub2, state1, state2)
        # check that layer boundary flux is zero due to full saturation
        @test iszero(state1.jw[end])
        @test iszero(state2.jw[1])
        # case 3: upper dry
        state1.sat .= minwater(sub1, water)
        state2.sat .= 0.5
        water = CryoGrid.processes(sub1)
        CryoGrid.computediagnostic!(sub1, water, state1)
        CryoGrid.computediagnostic!(sub2, water, state2)
        CryoGrid.interact!(sub1, sub2, state1, state2)
        # check that layer boundary flux is zero due to full saturation
        @test iszero(state1.jw[end])
        @test iszero(state2.jw[1])
    end
end
