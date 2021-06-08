using CryoGrid
using LinearAlgebra
using Statistics
using Test

include("../../types.jl")

@testset "Sanity checks" begin
	x = Grid(exp.(0.0:0.01:1.0)u"m")
	xc = cells(x)
	T₀ = uconvert.(u"K",(1.0./(ustrip.(xc).+0.01))u"°C")
	k = collect(LinRange(0.5,5.0,length(x)))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	∂H = zeros(length(T₀))u"J/s/m^3"
	@inferred heatconduction!(T₀,ΔT,k,Δk,∂H)
	# conditions based on initial temperature gradient
	@test ∂H[1] < 0.0u"J/s/m^3"
	@test ∂H[end] > 0.0u"J/s/m^3"
	@test sum(∂H) <= 0.0u"J/s/m^3"
end
@testset "Bounday conditions" begin
	x = Grid(Vector(0.0:0.01:1.0)u"m")
	xc = cells(x)
	k = collect(LinRange(0.5,1.0,length(x)))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	sub = TestGroundLayer()
	heat = Heat{:H}()
	bc = Constant{Heat,Dirichlet}(uconvert(u"K",0.0u"°C"))
	@testset "top: +, bot: -" begin
		T₀ = Vector(LinRange(250,300,length(xc)))u"K"
		∂H = zeros(length(T₀))u"J/s/m^3"
		state = (T=T₀,k=k,dH=∂H,grids=(T=xc,k=x),t=0.0)
		@test boundaryflux(Top(),bc,sub,heat,state,state) > 0.0u"J/s/m^3"
		@test boundaryflux(Bottom(),bc,sub,heat,state,state) < 0.0u"J/s/m^3"
	end
	@testset "top: -, bot: +" begin
		T₀ = Vector(LinRange(300,250,length(xc)))u"K"
		∂H = zeros(length(T₀))u"J/s/m^3"
		state = (T=T₀,k=k,dH=∂H,grids=(T=xc,k=x),t=0.0)
		@test boundaryflux(Top(),bc,sub,heat,state,state) < 0.0u"J/s/m^3"
		@test boundaryflux(Bottom(),bc,sub,heat,state,state) > 0.0u"J/s/m^3"
	end
	@testset "inner edge boundary (positive)" begin
		T₀ = Vector(sin.(ustrip.(xc).*π).+273.15)u"K"
		∂H = zeros(length(T₀))u"J/s/m^3"
		heatconduction!(T₀,ΔT,k,Δk,∂H)
		@test ∂H[1] > 0.0u"J/s/m^3"
		@test ∂H[end] > 0.0u"J/s/m^3"
	end
	@testset "inner edge boundary (negative)" begin
		T₀ = Vector(-sin.(ustrip.(xc).*π).+273.15)u"K"
		∂H = zeros(length(T₀))u"J/s/m^3"
		heatconduction!(T₀,ΔT,k,Δk,∂H)
		@test ∂H[1] < 0.0u"J/s/m^3"
		@test ∂H[end] < 0.0u"J/s/m^3"
	end
end
@testset "Fourier solution" begin
	x = Grid(Vector(0.0:0.01:1.0)u"m")
	xc = cells(x)
	k = ones(length(x))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	# Fourier's solution to heat equation with Dirichlet boundaries
	T₀ = uconvert.(u"K",(sin.(2π.*ustrip.(xc)))u"°C")
	f_analytic(x,t) = exp(-t*4π^2)*sin(2.0*π*x)
	sub = TestGroundLayer()
	heat = Heat{:H}()
	bc = Constant{Heat,Dirichlet}(uconvert(u"K",0.0u"°C"))
	function dTdt(T,p,t)
		dT = similar(T)u"J/s/m^3"
		dT .= zero(eltype(dT))
		T_K = (T)u"K"
		heatconduction!(T_K,ΔT,k,Δk,dT)
		# compute boundary fluxes;
		# while not correct in general, for this test case we can just re-use state for both layers.
		state = (T=T_K,k=k,dH=dT,grids=(T=xc,k=x),t=t)
		dT[1] += boundaryflux(Top(),bc,sub,heat,state,state)
		dT[end] += boundaryflux(Bottom(),bc,sub,heat,state,state)
		# strip units from dT before returning dT to the solver
		return ustrip.(dT)
	end
	tspan = (0.0,0.5)
	prob = ODEProblem(dTdt,ustrip.(T₀),tspan)
	# Forward Euler scheme with small step size
	sol = solve(prob,Euler(),dt=1.0e-5,saveat=0.01)
	# build solution matrices: time x depth
	res_sol = hcat(sol.u...)'u"K"
	res_analytic = (f_analytic.(ustrip.(xc),sol.t')'.+273.15)u"K"
	ϵ = 1.0e-5u"K" # error tolerance
	# verify inifnity norm of abs error < ϵ K (i.e. all predicted values differ < ϵ K)
	@test norm(mean(abs.(res_sol .- res_analytic),dims=1), Inf) < ϵ
	# verify convergence, last values should be almost identical
	@test norm(mean(abs.(res_sol[end,:] .- res_analytic[end,:]),dims=1), Inf) < ϵ
end