using CryoGrid
using Test

@testset "Sources" begin
    @testset "Constant" begin
        heatsource = Source(HeatBalance, Sources.Constant(1.0u"W/m^3"))
        layer = TestGroundLayer(heatsource)
        state = (dH=zeros(100)u"W/m^3",)
        computeprognostic!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3")
        state = (dH=ones(100)u"W/m^3",)
        computeprognostic!(layer, heatsource, state)
        @test all(state.dH .≈ 2.0u"W/m^3")
    end
    @testset "Periodic" begin
        level, amp, freq, shift = 0.0u"W/m^3", 2.0u"W/m^3", 0.5u"Hz", π/2
        heatsource = Source(HeatBalance, Sources.Periodic(; level, amp, freq, shift))
        layer = TestGroundLayer(heatsource)
        state = (dH=zeros(100)u"W/m^3", t=1.0u"s")
        computeprognostic!(layer, heatsource, state)
        @test all(state.dH .≈ level + amp*sin(2π*freq*1.0u"s" - shift))
        state = (dH=ones(100)u"W/m^3", t=1.5u"s")
        computeprognostic!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3" + level + amp*sin(2π*freq*1.5u"s" - shift))
    end
end
