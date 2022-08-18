using CryoGrid
using Test

include("../../types.jl")

@testset "Sources" begin
    layer = TestGroundLayer()
    @testset "Constant" begin
        heatsource = Source(Heat, Sources.Constant(1.0u"W/m^3"))
        state = (∂H∂t=zeros(100)u"W/m^3",)
        prognosticstep!(layer, heatsource, state)
        @test all(state.∂H∂t .≈ 1.0u"W/m^3")
        state = (∂H∂t=ones(100)u"W/m^3",)
        prognosticstep!(layer, heatsource, state)
        @test all(state.∂H∂t .≈ 2.0u"W/m^3")
    end
    @testset "Periodic" begin
        level, amp, freq, shift = 0.0u"W/m^3", 2.0u"W/m^3", 0.5u"Hz", π/2
        heatsource = Source(Heat, Sources.Periodic(; level, amp, freq, shift))
        state = (∂H∂t=zeros(100)u"W/m^3", t=1.0u"s")
        prognosticstep!(layer, heatsource, state)
        @test all(state.∂H∂t .≈ level + amp*sin(2π*freq*1.0u"s" - shift))
        state = (∂H∂t=ones(100)u"W/m^3", t=1.5u"s")
        prognosticstep!(layer, heatsource, state)
        @test all(state.∂H∂t .≈ 1.0u"W/m^3" + level + amp*sin(2π*freq*1.5u"s" - shift))
    end
end
