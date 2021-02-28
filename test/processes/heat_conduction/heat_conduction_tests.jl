using CryoGrid
using FastClosures
using LinearAlgebra
using Statistics
using Test

include("../../types.jl")

@testset "Heat conduction" begin
	x = Grid(exp.(0.0:0.01:1.0)u"m")
	xc = cells(x)
	# normalize cell centers (dimensionless 0..1)
	nxc = xc .- minimum(xc)
	nxc /= maximum(nxc)
	T₀ = uconvert.(u"K",(1.0./(nxc.+0.01))u"°C")
	k = collect(LinRange(0.5,5.0,length(x)))u"W/m/K"
	ΔT = Δ(xc)
	Δk = Δ(x)
	@testset "Sanity checks" begin
		∂H = zeros(length(T₀))u"J/s/m^3"
		@inferred heatconduction!(T₀,ΔT,k,Δk,∂H)
		# conditions based on initial temperature gradient
		@test ∂H[1] < 0.0u"J/s/m^3"
		@test ∂H[end] > 0.0u"J/s/m^3"
		@test sum(∂H) <= 0.0u"J/s/m^3"
	end
	@testset "Fourier solution" begin
		∂H = zeros(length(T₀))u"J/s/m^3"
		# Fourier's solution to heat equation with Dirichlet boundaries
		T₀ = uconvert.(u"K",(sin.(2.0.*π.*nxc))u"°C")
		f_analytic(t) = let x=xc; exp(-pi^2*t).*sin.(2.0*π.*nxc) end
		sub = TestGroundLayer()
		heat = Heat{UT"J"}()
		bc = Constant{Heat,Dirichlet}(uconvert(u"K",0.0u"°C"))
		function dTdt(T,p,t)
			dT = similar(∂H)
			dT .= zero(eltype(dT))
			heatconduction!(T,ΔT,k,Δk,dT)
			# compute boundary fluxes;
			# while not correct in general, for this test case we can just re-use state for both layers.
			state = (T=T,k=k,dH=dT,grids=(T=xc,k=x),t=t)
			dT[1] += boundaryflux(Top(),bc,sub,heat,state,state)
			dT[end] += boundaryflux(Bottom(),bc,sub,heat,state,state)
			# assume heat capacity = 1.0 J/(Km^3) and convert to K
			# ODEProblem assumes that du and u have same type, so we do s*K/s = K to make it happy
			return 1.0u"s".*dT./1.0u"J/(K*m^3)"
		end
		tspan = (0.0,2.0)
		prob = ODEProblem(dTdt,T₀,tspan)
		# Forward Euler scheme with small step size
		sol = solve(prob,Euler(),dt=1.0e-5,saveat=0.01)
		# build solution matrices: time x depth
		u_sol = hcat(sol.u...)'
		u_analytic = (hcat(f_analytic.(sol.t)...)'.+273.15)u"K"
		# verify inifnity norm of abs error < 1.0 K (i.e. all predicted values <1.0 K error)
		@test norm(mean(abs.(u_sol .- u_analytic),dims=1), Inf) < 1.0u"K"
		# verify convergence, last values should be almost identical
		@test norm(mean(abs.(u_sol[end,:] .- u_analytic[end,:]),dims=1), Inf) < 1.0e-4u"K"
	end
end
