using CryoGrid
using Dates
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
		∂H = zeros(length(T₀))u"J/s/m^3"
		@inferred heatconduction!(∂H,T₀,ΔT,k,Δk)
		# conditions based on initial temperature gradient
		@test ∂H[1] < 0.0u"J/s/m^3"
		@test ∂H[end] > 0.0u"J/s/m^3"
		@test sum(∂H) <= 0.0u"J/s/m^3"
	end
	# check boundary fluxes
	@testset "Boundary fluxes" begin
		x = Grid(Vector(0.0:0.01:1.0)u"m")
		xc = cells(x)
		k = collect(LinRange(0.5,1.0,length(x)))u"W/m/K"
		ΔT = Δ(xc)
		Δk = Δ(x)
		sub = TestGroundLayer()
		heat = Heat{:H}()
		bc = Constant{Dirichlet}(0.0u"°C")
		@testset "top: +, bot: -" begin
			T₀ = Vector(LinRange(-23,27,length(xc)))u"°C"
			∂H = zeros(length(T₀))u"J/s/m^3"
			state = (T=T₀,k=k,dH=∂H,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(Top(),bc,sub,heat,state,state) > 0.0u"J/s/m^3"
			@test boundaryflux(Bottom(),bc,sub,heat,state,state) < 0.0u"J/s/m^3"
		end
		@testset "top: -, bot: +" begin
			T₀ = Vector(LinRange(27,-23,length(xc)))u"°C"
			∂H = zeros(length(T₀))u"J/s/m^3"
			state = (T=T₀,k=k,dH=∂H,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(Top(),bc,sub,heat,state,state) < 0.0u"J/s/m^3"
			@test boundaryflux(Bottom(),bc,sub,heat,state,state) > 0.0u"J/s/m^3"
		end
		@testset "inner edge boundary (positive)" begin
			T₀ = Vector(sin.(ustrip.(xc).*π))u"°C"
			∂H = zeros(length(T₀))u"J/s/m^3"
			heatconduction!(∂H,T₀,ΔT,k,Δk)
			@test ∂H[1] > 0.0u"J/s/m^3"
			@test ∂H[end] > 0.0u"J/s/m^3"
		end
		@testset "inner edge boundary (negative)" begin
			T₀ = Vector(-sin.(ustrip.(xc).*π))u"°C"
			∂H = zeros(length(T₀))u"J/s/m^3"
			heatconduction!(∂H,T₀,ΔT,k,Δk)
			@test ∂H[1] < 0.0u"J/s/m^3"
			@test ∂H[end] < 0.0u"J/s/m^3"
		end
	end
end
@testset "Boundary conditions" begin
	@testset "n-factors" begin
		ts = DateTime(2010,1,1):Hour(1):DateTime(2010,1,1,4)
		forcing = TimeSeriesForcing([1.0,0.5,-0.5,-1.0,0.1], ts, :Tair)
		tgrad = TemperatureGradient(forcing, NFactor(:stationary, :nf))
		vars = variables(Top(), tgrad)
		@test length(vars) == 2
		@test getscalar(vars[1].default_value) == 1.0
		@test getscalar(vars[2].default_value) == 1.0
		sub = TestGroundLayer()
		heat = Heat{:H}()
		# stationary case
		f1(t) = tgrad(Top(),sub,heat,(t=t, params=(nfw=0.5, nfs=1.0)), (t=t,))
		Tres = f1.(Dates.datetime2epochms.(ts)./1000.0)
		@test all(Tres .≈ [1.0,0.5,-0.25,-0.5,0.1])
		# two-stage non-stationary case
		tgrad = TemperatureGradient(forcing, NFactor(:twostage, :nf))
		vars = variables(Top(), tgrad)
		@test length(vars) == 5
		@test all([getscalar(vars[i].default_value) == 1.0 for i in 1:4])
		@test getscalar(vars[5].default_value) == 0.0
		tchange = Dates.datetime2epochms(ts[end-1])./1000
		f2(t) = tgrad(Top(),sub,heat,(t=t, params=(nfw1=0.5, nfw2=0.2, nfs1=1.0, nfs2=0.9, nftc=tchange)), (t=t,))
		Tres = f2.(Dates.datetime2epochms.(ts)./1000)
		@test all(Tres .≈ [1.0,0.5,-0.25,-0.2,0.09])
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
	heat = Heat{:H}()
	bc = Constant{Dirichlet}(0.0u"°C")
	function dTdt(T,p,t)
		dT = similar(T)u"J/s/m^3"
		dT .= zero(eltype(dT))
		T_K = (T)u"°C"
		heatconduction!(dT,T_K,ΔT,k,Δk)
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
	res_sol = hcat(sol.u...)'u"°C"
	res_analytic = (f_analytic.(ustrip.(xc),sol.t')')u"°C"
	ϵ = 1.0e-5u"K" # error tolerance
	# verify inifnity norm of abs error < ϵ K (i.e. all predicted values differ < ϵ K)
	@test norm(mean(abs.(res_sol .- res_analytic),dims=1), Inf) < ϵ
	# verify convergence, last values should be almost identical
	@test norm(mean(abs.(res_sol[end,:] .- res_analytic[end,:]),dims=1), Inf) < ϵ
end
