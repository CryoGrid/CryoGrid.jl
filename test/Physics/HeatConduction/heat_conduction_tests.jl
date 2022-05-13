using CryoGrid
using Dates
using DiffEqBase
using OrdinaryDiffEq
using LinearAlgebra
using Statistics
using Test

include("../../types.jl")

@testset "Sanity checks" begin
	@testset "Generic" begin
		x = Grid(exp.(0.0:0.01:1.0)u"m")
		xc = cells(x)
		T₀ = (1.0./(ustrip.(xc).+0.01))u"°C"
		k = collect(LinRange(0.5,5.0,length(x)))u"W/m/K"
		ΔT = Δ(xc)
		Δk = Δ(x)
		∂H = zeros(length(T₀))u"W/m^3"
		@inferred heatconduction!(∂H,T₀,ΔT,k,Δk)
		# conditions based on initial temperature gradient
		@test ∂H[1] < 0.0u"W/m^3"
		@test ∂H[end] > 0.0u"W/m^3"
		@test sum(∂H) <= 0.0u"W/m^3"
	end
	# check boundary fluxes
	@testset "Boundary fluxes" begin
		x = Grid(Vector(0.0:0.01:1.0)u"m")
		xc = cells(x)
		k = collect(LinRange(0.5,1.0,length(x)))u"W/m/K"
		ΔT = Δ(xc)
		Δk = Δ(x)
		sub = TestGroundLayer()
		heat = Heat()
		bc = ConstantBC(Heat, CryoGrid.Dirichlet, 0.0u"°C")
		@testset "top: +, bot: -" begin
			T₀ = Vector(LinRange(-23,27,length(xc)))u"°C"
			∂H = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=∂H,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(),heat,sub,state,state) > 0.0u"W/m^2"
			@test boundaryflux(bc,Bottom(),heat,sub,state,state) < 0.0u"W/m^2"
		end
		@testset "top: -, bot: +" begin
			T₀ = Vector(LinRange(27,-23,length(xc)))u"°C"
			∂H = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=∂H,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(),heat,sub,state,state) < 0.0u"W/m^2"
			@test boundaryflux(bc,Bottom(),heat,sub,state,state) > 0.0u"W/m^2"
		end
		@testset "inner edge boundary (positive)" begin
			T₀ = Vector(sin.(ustrip.(xc).*π))u"°C"
			∂H = zeros(length(T₀))u"W/m^3"
			heatconduction!(∂H,T₀,ΔT,k,Δk)
			@test ∂H[1] > 0.0u"W/m^3"
			@test ∂H[end] > 0.0u"W/m^3"
		end
		@testset "inner edge boundary (negative)" begin
			T₀ = Vector(-sin.(ustrip.(xc).*π))u"°C"
			∂H = zeros(length(T₀))u"W/m^3"
			heatconduction!(∂H,T₀,ΔT,k,Δk)
			@test ∂H[1] < 0.0u"W/m^3"
			@test ∂H[end] < 0.0u"W/m^3"
		end
		@testset "Neumann boundary" begin
			bc = ConstantBC(Heat, CryoGrid.Neumann, -1.0u"W/m^2")
			T₀ = Vector(LinRange(-23,27,length(xc)))u"°C"
			∂H = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=∂H,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(),heat,sub,state,state) == -1.0u"W/m^2"
		end
	end
end
@testset "Boundary conditions" begin
	@testset "n-factors" begin
		ts = DateTime(2010,1,1):Hour(1):DateTime(2010,1,1,4)
		forcing = TimeSeriesForcing([1.0,0.5,-0.5,-1.0,0.1], ts, :Tair)
		tgrad = TemperatureGradient(forcing, NFactor(nf=0.5, nt=1.0))
		sub = TestGroundLayer()
		heat = Heat()
		function f1(t)
			state = (T_ub=[Inf], nfactor=[Inf], t=t)
			diagnosticstep!(Top(), tgrad, state)
			return boundaryvalue(tgrad,Top(),heat,sub,state,state)
		end
		Tres = f1.(Dates.datetime2epochms.(ts)./1000.0)
		@test all(Tres .≈ [1.0,0.5,-0.25,-0.5,0.1])
	end
end
@testset "Fourier solution" begin
	x = Grid(Vector(0.0:0.01:1.0)u"m")
	xc = cells(x)
	k = ones(length(x))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	# Fourier's solution to heat equation with Dirichlet boundaries
	T₀ = (sin.(2π.*ustrip.(xc)))u"°C"
	f_analytic(x,t) = exp(-t*4π^2)*sin(2.0*π*x)
	sub = TestGroundLayer()
	heat = Heat()
	bc = ConstantBC(Heat, CryoGrid.Dirichlet, 0.0u"°C")
	function dTdt(T,p,t)
		dH = similar(T)u"W/m^3"
		dH .= zero(eltype(dH))
		T_K = (T)u"°C"
		heatconduction!(dH,T_K,ΔT,k,Δk)
		# compute boundary fluxes;
		# while not correct in general, for this test case we can just re-use state for both layers.
		state = (T=T_K,k=k,dH=dH,grid=x,grids=(T=xc,k=x),t=t)
		dH[1] += boundaryflux(bc,Top(),heat,sub,state,state)/Δk[1]
		dH[end] += boundaryflux(bc,Bottom(),heat,sub,state,state)/Δk[end]
		# strip units from dH before returning it to the solver
		return ustrip.(dH)
	end
	tspan = (0.0,0.5)
	prob = ODEProblem(dTdt,ustrip.(T₀),tspan)
	# Forward Euler scheme with small step size
	sol = solve(prob,Euler(),dt=1.0e-5,saveat=0.01)
	# build solution matrices: time x depth
	res_sol = hcat(sol.u...)'u"°C"
	res_analytic = (f_analytic.(ustrip.(xc),sol.t')')u"°C"
	ϵ = 1.0e-5u"K" # error tolerance
	# verify inifnity norm of abs error < ϵ K (i.e. all predicted values differ < ϵ K)
	@test norm(mean(abs.(res_sol .- res_analytic),dims=1), Inf) < ϵ
	# verify convergence, last values should be almost identical
	@test norm(mean(abs.(res_sol[end,:] .- res_analytic[end,:]),dims=1), Inf) < ϵ
end