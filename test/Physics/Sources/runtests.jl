using CryoGrid
using Test

include("../../types.jl")

@testset "Sources" begin
    layer = TestGroundLayer()
    @testset "Constant" begin
        heatsource = parameterize(Source(Heat, Sources.Constant()))
        heatsource = CryoGrid.update(heatsource, [1.0u"W/m^3"])
        state = (dH=zeros(100)u"W/m^3",)
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3")
        state = (dH=ones(100)u"W/m^3",)
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 2.0u"W/m^3")
    end
    @testset "Periodic" begin
        level, amp, freq, shift = 0.0u"W/m^3", 2.0u"W/m^3", 0.5u"Hz", π/2
        heatsource = parameterize(Source(Heat, Sources.Periodic()))
        heatsource = CryoGrid.update(heatsource, [amp, freq, shift, level])
        state = (dH=zeros(100)u"W/m^3", t=1.0u"s")
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ level + amp*sin(2π*freq*1.0u"s" - shift))
        state = (dH=ones(100)u"W/m^3", t=1.5u"s")
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3" + level + amp*sin(2π*freq*1.5u"s" - shift))
    end
end
