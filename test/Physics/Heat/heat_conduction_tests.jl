using CryoGrid
using CryoGrid.Numerics: nonlineardiffusion!
using Dates
using DiffEqBase
using OrdinaryDiffEq
using LinearAlgebra
using Statistics
using Test

@testset "Sanity checks" begin
	@testset "Generic" begin
		x = Grid(exp.(0.0:0.01:1.0)u"m")
		xc = cells(x)
		T₀ = (1.0./(ustrip.(xc).+0.01))u"°C"
		k = collect(LinRange(0.5,5.0,length(x)))u"W/m/K"
		ΔT = Δ(xc)
		Δk = Δ(x)
		jH = zeros(length(x))u"W/m^2"
		dH = zeros(length(T₀))u"W/m^3"
		@inferred nonlineardiffusion!(dH, jH, T₀, ΔT, k, Δk)
		# conditions based on initial temperature gradient
		@test dH[1] < 0.0u"W/m^3"
		@test dH[end] > 0.0u"W/m^3"
	end
	# check boundary fluxes
	@testset "Boundary fluxes" begin
		x = Grid(Vector(0.0:0.01:1.0)u"m")
		xc = cells(x)
		k = collect(LinRange(0.5,1.0,length(x)))u"W/m/K"
		ΔT = Δ(xc)
		Δk = Δ(x)
		heat = HeatBalance()
		sub = TestGroundLayer(heat)
		bc = ConstantBC(HeatBalance, CryoGrid.Dirichlet, 0.0u"°C")
		# remember that fluxes are positive *downward*!!!
		@testset "top: +, bot: -" begin
			# initial condition: -1°C
			T₀ = -ones(length(xc))u"°C"
			jH = zeros(length(x))u"W/m^2"
			dH = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=dH,jH=jH,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(bc),heat,sub,state,state) > 0.0u"W/m^2"
			@test boundaryflux(bc,Bottom(bc),heat,sub,state,state) < 0.0u"W/m^2"
		end
		@testset "top: -, bot: +" begin
			# initial condition: 1°C
			T₀ = ones(length(xc))u"°C"
			jH = zeros(length(x))u"W/m^2"
			dH = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=dH,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(bc),heat,sub,state,state) < 0.0u"W/m^2"
			@test boundaryflux(bc,Bottom(bc),heat,sub,state,state) > 0.0u"W/m^2"
		end
		@testset "inner edge boundary (positive)" begin
			T₀ = Vector(sin.(ustrip.(xc).*π))u"°C"
			jH = zeros(length(x))u"W/m^2"
			dH = zeros(length(T₀))u"W/m^3"
			nonlineardiffusion!(dH,jH,T₀,ΔT,k,Δk)
			@test dH[1] > 0.0u"W/m^3"
			@test dH[end] > 0.0u"W/m^3"
		end
		@testset "inner edge boundary (negative)" begin
			T₀ = Vector(-sin.(ustrip.(xc).*π))u"°C"
			jH = zeros(length(x))u"W/m^2"
			dH = zeros(length(T₀))u"W/m^3"
			nonlineardiffusion!(dH,jH,T₀,ΔT,k,Δk)
			@test dH[1] < 0.0u"W/m^3"
			@test dH[end] < 0.0u"W/m^3"
		end
		@testset "Neumann boundary" begin
			bc = ConstantBC(HeatBalance, CryoGrid.Neumann, -1.0u"W/m^2")
			T₀ = Vector(LinRange(-23,27,length(xc)))u"°C"
			jH = zeros(length(x))u"W/m^2"
			dH = zeros(length(T₀))u"W/m^3"
			state = (T=T₀,k=k,dH=dH,jH=jH,grid=x,grids=(T=xc,k=x),t=0.0)
			@test boundaryflux(bc,Top(bc),heat,sub,state,state) == -1.0u"W/m^2"
		end
	end
end
@testset "Boundary conditions" begin
	@testset "n-factors" begin
		T_pos = TemperatureBC(10.0, NFactor(nf=0.5, nt=0.75))
		T_neg = TemperatureBC(-10.0, NFactor(nf=0.5, nt=0.75))
		heat = HeatBalance()
		sub = TestGroundLayer(heat)
		zerobc = ConstantBC(HeatBalance, CryoGrid.Dirichlet, 0.0u"°C")
		function f1(t)
			state = (T_ub=[Inf], nfactor=[Inf], t=t)
			computediagnostic!(Top(T_pos), T_pos, state)
			T_ub = boundaryvalue(T_pos, state)
			computediagnostic!(Top(T_neg), T_neg, state)
			T_lb = boundaryvalue(T_neg, state)
			return T_ub, T_lb
		end
		T1, T2 = f1(0.0)
		@test T1 == 7.5
		@test T2 == -5.0
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
	heat = HeatBalance()
	sub = TestGroundLayer(heat)
	bc = ConstantBC(HeatBalance, CryoGrid.Dirichlet, 0.0u"°C")
	function ∂T∂t(u,p,t)
		dH = similar(u)u"W/m^3"
		dH .= zero(eltype(dH))
		jH = similar(u, length(u)+1)u"W/m^2"
		jH .= zero(eltype(jH))
		T = (u)u"°C"
		# compute boundary fluxes;
		# while not correct in general, for this test case we can just re-use state for both layers.
		state = (T=T,dH=dH,jH=jH,k=k,grid=x,grids=(T=xc,k=x),t=t)
		interact!(Top(bc), bc, sub, heat, state, state)
		interact!(sub, heat, Bottom(bc), bc, state, state)
		computeprognostic!(sub, heat, state)
		# strip units from dH before returning it to the solver;
		# note that we do not need to divide by diffusivity since we assume it to be unity
		return ustrip.(dH)
	end
	tspan = (0.0,0.5)
	prob = ODEProblem(∂T∂t,ustrip.(T₀),tspan)
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
