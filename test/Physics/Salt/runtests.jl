using CryoGrid
using CryoGrid.Salt
using CryoGrid.Utils

using Test

include("../../testutils.jl")

@testset "Salt" begin
    testgrid = CryoGrid.Presets.DefaultGrid_2cm
    @testset "variables" begin
        sediment = MarineSediment()
        vars = CryoGrid.variables(sediment)
        prognostic_vars = filter(CryoGrid.isprognostic, vars)
        diagnostic_vars = filter(CryoGrid.isdiagnostic, vars)
        # check that porosity is defined as a state variable
        @test :por ∈ map(varname, diagnostic_vars)
        # check for salt concentation and its flux
        @test :c ∈ map(varname, prognostic_vars)
        @test :jc ∈ map(varname, diagnostic_vars)
        # check processes
        procs = CryoGrid.processes(sediment)
        @test isa(procs[1], SaltMassBalance)
        @test isa(procs[2], HeatBalance)
    end
    @testset "updatestate!" begin
        sediment = pstrip(MarineSediment())
        state = Diagnostics.build_dummy_state(testgrid, sediment, with_units=false)
        # initialize variables
        state.T .= -2.0
        state.c .= 800.0
        state.por .= porosity(sediment)
        procs = CryoGrid.processes(sediment)
        CryoGrid.updatestate!(sediment, procs, state)
        # check that 
        @test allfinite(state.k)
        @test allfinite(state.θw)
        # TODO: need to set Tmelt for diagnostics
        @test_broken all(state.Tmelt .< zero(eltype(state.T)))
        @test all(state.θw .> zero(eltype(state.θw)))
        @test all(state.θw .< state.por)
    end
    @testset "interact!" begin
        sediment1 = pstrip(MarineSediment())
        sediment2 = pstrip(MarineSediment())
        state1 = Diagnostics.build_dummy_state(testgrid, sediment1, with_units=false)
        state2 = Diagnostics.build_dummy_state(testgrid, sediment2, with_units=false)
        # initialize variables
        state1.T .= -2.0
        state2.T .= -1.0
        state1.c .= 800.0
        state2.c .= 700.0
        state1.por .= porosity(sediment1)
        state2.por .= porosity(sediment2)
        procs = CryoGrid.processes(sediment1)
        CryoGrid.updatestate!(sediment1, procs, state1)
        CryoGrid.updatestate!(sediment2, procs, state2)
        CryoGrid.interact!(sediment1, sediment2, state1, state2)
        # check that 
        @test allfinite(state1.k) && allfinite(state2.k)
        @test allfinite(state1.θw) && allfinite(state2.θw)
        @test allfinite(state1.jc) && allfinite(state2.jc)
        @test allfinite(state1.jH) && allfinite(state2.jH)
        # check that flux is positive since c1 > c2
        @test state1.jc[end] > zero(eltype(state1.jc))
        @test state2.jc[1] == state1.jc[end]
    end
    @testset "computefluxes!" begin
        sediment = pstrip(MarineSediment())
        state = Diagnostics.build_dummy_state(testgrid, sediment, with_units=false)
        # initialize variables
        state.T .= LinRange(-2.0, -3.0, length(cells(testgrid)))
        state.c .= LinRange(900.0, 500.0, length(cells(testgrid)))
        state.por .= porosity(sediment)
        procs = CryoGrid.processes(sediment)
        CryoGrid.updatestate!(sediment, procs, state)
        CryoGrid.computefluxes!(sediment, procs, state)
        # check that all divergences are nonzero for both temperature and salt
        @test all(.!iszero.(state.∂c∂t))
        @test all(.!iszero.(state.∂T∂t))
    end
end

sediment = pstrip(MarineSediment())
state = Diagnostics.build_dummy_state(testgrid, sediment, with_units=false)

