using CryoGrid
using Test

include("../../types.jl")

@testset "Sources" begin
    layer = TestGroundLayer()
    @testset "Constant" begin
        heatsource = Source(Heat, Sources.Constant(:S₀))
        vars = variables(layer, heatsource)
        @test length(vars) == 1
        @test varname(vars[1]) == :S₀
        @test length(vars[1].default_value) == 1
        state = (dH=zeros(100)u"W/m^3", params=(S₀=1.0u"W/m^3",))
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3")
        state = (dH=ones(100)u"W/m^3", params=(S₀=1.0u"W/m^3",))
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 2.0u"W/m^3")
    end
    @testset "Periodic" begin
        heatsource = Source(Heat, Sources.Periodic(:S₀,0.0u"W/m^3"))
        vars = variables(layer, heatsource)
        @test length(vars) == 3
        @test all([length(var.default_value) == 1 for var in vars])
        amp, freq, shift = 2.0u"W/m^3", 0.5u"Hz", π/2
        params = (S₀_amp=amp, S₀_freq=freq, S₀_shift=shift)
        state = (dH=zeros(100)u"W/m^3", t=1.0u"s", params=params)
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ amp*sin(2π*freq*1.0u"s" - shift))
        state = (dH=ones(100)u"W/m^3", t=1.5u"s", params=params)
        prognosticstep!(layer, heatsource, state)
        @test all(state.dH .≈ 1.0u"W/m^3" + amp*sin(2π*freq*1.5u"s" - shift))
    end
end
