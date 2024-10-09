using CryoGrid
using CryoGrid.Salt
using CryoGrid.Utils

using Test

@testset "Salt" begin
    testgrid = CryoGrid.Presets.DefaultGrid_2cm
    @testset "variables" begin
        soil = SalineGround()
        vars = CryoGrid.variables(soil)
        prognostic_vars = filter(CryoGrid.isprognostic, vars)
        diagnostic_vars = filter(CryoGrid.isdiagnostic, vars)
        # check that porosity is defined as a state variable
        @test :por ∈ map(varname, diagnostic_vars)
        # check for salt concentation and its flux
        @test :c ∈ map(varname, prognostic_vars)
        @test :jc ∈ map(varname, diagnostic_vars)
        # check processes
        procs = CryoGrid.processes(soil)
        @test isa(procs[1], WaterBalance)
        @test isa(procs[2], Salt.CoupledHeatSalt)
    end
    @testset "computediagnostic!" begin
        soil = pstrip(SalineGround())
        state = Diagnostics.build_dummy_state(testgrid, soil, with_units=false)
        # initialize variables
        state.T .= -2.0
        state.c .= 800.0
        state.por .= porosity(soil)
        procs = CryoGrid.processes(soil)
        CryoGrid.computediagnostic!(soil, procs, state)
        # check that 
        @test allfinite(state.k)
        @test allfinite(state.θw)
        # TODO: need to set Tmelt for diagnostics
        @test_broken all(state.Tmelt .< zero(eltype(state.T)))
        @test all(state.θw .> zero(eltype(state.θw)))
        @test all(state.θw .< state.por)
    end
    @testset "interact!" begin
        soil1 = pstrip(SalineGround())
        soil2 = pstrip(SalineGround())
        state1 = Diagnostics.build_dummy_state(testgrid[0.0u"m"..10.0u"m"], soil1, with_units=false)
        state2 = Diagnostics.build_dummy_state(testgrid[10.0u"m"..1000.0u"m"], soil2, with_units=false)
        # initialize variables
        state1.T .= -2.0
        state2.T .= -1.0
        state1.c .= 800.0
        state2.c .= 700.0
        state1.por .= porosity(soil1)
        state2.por .= porosity(soil2)
        procs = CryoGrid.processes(soil1)
        CryoGrid.computediagnostic!(soil1, procs, state1)
        CryoGrid.computediagnostic!(soil2, procs, state2)
        CryoGrid.interact!(soil1, soil2, state1, state2)
        # check that 
        @test allfinite(state1.k) && allfinite(state2.k)
        @test allfinite(state1.θw) && allfinite(state2.θw)
        @test allfinite(state1.jc) && allfinite(state2.jc)
        @test allfinite(state1.jH) && allfinite(state2.jH)
        # check that flux is positive since c1 > c2
        @test state1.jc[end] > zero(eltype(state1.jc))
        @test state2.jc[1] == state1.jc[end]
    end
    @testset "computeprognostic!" begin
        soil = pstrip(SalineGround())
        state = Diagnostics.build_dummy_state(testgrid, soil, with_units=false)
        # initialize variables
        state.T .= LinRange(-2.0, -3.0, length(cells(testgrid)))
        state.c .= LinRange(900.0, 500.0, length(cells(testgrid)))
        state.por .= porosity(soil)
        procs = CryoGrid.processes(soil)
        CryoGrid.computediagnostic!(soil, procs, state)
        CryoGrid.computeprognostic!(soil, procs, state)
        # check that all divergences are nonzero for both temperature and salt
        @test all(.!iszero.(state.dc))
        @test all(.!iszero.(state.dT))
    end
end
