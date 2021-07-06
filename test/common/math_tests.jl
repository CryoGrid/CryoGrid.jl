using CryoGrid.Numerics: finitediff!, lineardiffusion!, nonlineardiffusion!, ∇
using Test

include("../testutils.jl")

@testset "finitediff!" begin
	f(x) = (1/2)x^2 # function to differntiate
	df(x) = x # analytical 2nd derivative
	x = exp.(0.0:0.001:0.05)
	y = f.(x)
	dy = df.(x[1:end-1])
	δ = x[2:end] .- x[1:end-1]
	∂y = zeros(length(y)-1)
	finitediff!(∂y,y,δ)
	@test allfinite(∂y)
	@test allequal(∂y,dy,atol=0.01)
end

@testset "lineardiffusion!" begin
	f(x) = (1/6)x^3 # function to differntiate
	d2f(x) = x # analytical 2nd derivative
	x = exp.(0.0:0.001:0.05)
	y = f.(x)
	d2y = d2f.(x[2:end-1])
	k = 2.0 # constant diffusion
	δ = x[2:end] .- x[1:end-1]
	∂²y = zeros(length(y)-2)
	lineardiffusion!(∂²y,y,δ,k)
	@test allfinite(∂²y)
	@test allequal(∂²y,k*d2y,atol=0.01)
end

@testset "nonlineardiffusion!" begin
	f(x) = (1/6)x^3 # function to differntiate
	df(x) = (1/2)x^2 # analytical 1st derivative
	d2f(x) = x # analytical 2nd derivative
	k(x) = sin(x) # non-linear diffusion function
	dk(x) = cos(x) # analytical 1st derivative of k
	d2kf(x) = k(x)d2f(x) + dk(x)df(x) # by product rule
	x = exp.(0.0:0.001:0.05)
	xc = (x[1:end-1].+x[2:end])/2
	y = f.(xc)
	kₓ = k.(x)[2:end-1]
	d2y = d2kf.(xc[2:end-1]) # analytical solution
	δx = x[2:end] .- x[1:end-1]
	δxc = xc[2:end] .- xc[1:end-1]
	out = zeros(length(y)-2)
	nonlineardiffusion!(out,y,δxc,kₓ,δx[2:end-1])
	@test allfinite(out)
	@test allequal(out,d2y,atol=0.01)
end

@testset "∇" begin
	# test 1: single arg
	f(x) = x^2 + 2*x
	∇f(x) = 2x + 2
	∇f_gen = ∇(f, :x)
	xs = -10:0.1:10
	true_vals = [∇f(x) for x in xs]
	gen_vals = [∇f_gen(x) for x in xs]
	@test all(true_vals .≈ gen_vals)
	# test 2: with extra args
	f2(x,a,b) = a*x^2 + 2*x*b
	∇f2(x,a,b) = 2x*a + 2*b
	∇f2_gen = ∇(f2, :x)
	xs = -1:0.1:1
	as = -1:0.1:1
	bs = -1:0.1:1
	true_vals = [∇f2(x,a,b) for x in xs for a in as for b in bs]
	gen_vals = [∇f2_gen(x,a,b) for x in xs for a in as for b in bs]
	@test all(true_vals .≈ gen_vals)
	# test 3: with extra args, different order
	f3(a,x,b) = a*x^2 + 2*x*b
	∇f3(a,x,b) = 2x*a + 2*b
	∇f3_gen = ∇(f3, :x)
	xs = -1:0.1:1
	as = -1:0.1:1
	bs = -1:0.1:1
	true_vals = [∇f3(a,x,b) for x in xs for a in as for b in bs]
	gen_vals = [∇f3_gen(a,x,b) for x in xs for a in as for b in bs]
	@test all(true_vals .≈ gen_vals)
end
