using CryoGrid
using Test
using Plots
@testset "Heat conduction" begin
	x = Grid(exp.(0.0:0.01:1.0)u"m")
	xc = cells(x)
	T = (1.0./(ustrip.(xc).+0.1) .+ 273.15)u"K"
	k = collect(LinRange(0.5,5.0,length(x)))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	∂H = zeros(length(T))u"J/s/m^3"
	println(T)
	@inferred heatconduction(T,ΔT,k,Δk,∂H)
	@test ∂H[1] < 0.0u"J/s/m^3"
	@test ∂H[end] > 0.0u"J/s/m^3"
end
